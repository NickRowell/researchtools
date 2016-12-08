/**
 * Copyright (c) 2011 University of Dundee
 *
 * Name:
 * PoseDetermination.java
 *
 * Purpose:
 * Factory class that calculates solutions for 6 degree of freedom camera
 * pose given 2D/3D point correspondences.
 * 
 * Language:
 * Java
 *
 * Author:
 * Nicholas Rowell (nickrowell@computing.dundee.ac.uk)
 *
 */
package numeric.vision;


import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import Jama.SingularValueDecomposition;
import numeric.geom.dim3.Quaternion;
import numeric.geom.dim3.Vector3d;

import java.util.List;



public class PoseDetermination {
    
    /** 
     * Frame objects contain solution for 6 degree of freedom camera pose,
     * as well as camera projection matrix and List of 2D/3D point 
     * correspondences, so although we are really interested in the Pose, 
     * using a Frame object allows the extra information to be wrapped up
     * nicely.
     */
    public Frame frame;

    /** 
     * Validity flag, used to indicate e.g. if critical arrangement of points
     * has occurred or if too few Matches to determine solution.
     */
    public boolean ISVALID = false;
    
    /** Goodness-of-fit statistics */
    private double e1=-1,e2=-1;

    /** 
     * Centre of mass of landmarks in world frame. Fixed for a given 
     * correspondence set.
     */
    private Matrix COM_target = new Matrix(4,1);
    /** Body from COM-relative position vectors of Landmarks */
    private Matrix[] targetFrameCOMVectors  = new Matrix[0];
    
    public void setE1(double ee1){this.e1 = ee1;}
    public void setE2(double ee2){this.e2 = ee2;}    
    public double getE1(){ return e1;}
    public double getE2(){ return e2;}
    
    
 
    /**
     * Linear pose estimate, using correspondences between 3D points in 
     * the camera frame and their projected coordinates in the image plane.
     * 
     * This method implements the techniques described in the papers
     * 
     * "Linear N-point Camera Pose Determination" by Quan & Lan 1999
     * 
     * and
     * 
     * "Closed Form Solution of Absolute Orientation Using Unit Quaternions"
     * by K.P.Horn 1987.
     * 
     */
    public PoseDetermination(List<Correspondence_2D_3D> matches, Matrix k){

        // Copy reference to List of Correspondence_2D_3Ds
        frame.correspondences = matches;
        
        // Copy reference to camera matrix
        frame.K = k;
        
        // How many correspondences are there?
        int NCORR = frame.correspondences.size();
        
        // Immediately bail if there are too few matches
        if(NCORR<4)
            return;

        
        /**
         * First, measure intrinsic geometric properties of Landmarks in 
         * body frame that are independent of solution for Pose.
         */
        
        // Set COM of 3D points in target frame. This is an invariant.
        for(Correspondence_2D_3D point: frame.correspondences)      
            COM_target.plusEquals(point.p3d.x);
        COM_target.timesEquals(1.0/(double)NCORR);
        
        // Now set COM-relative position vectors
        targetFrameCOMVectors  = new Matrix[NCORR];

        int f=0;
        for(Correspondence_2D_3D point: frame.correspondences)
            targetFrameCOMVectors[f++] = point.p3d.x.minus(COM_target);
        
        /**
         * Now set same properties in camera frame. These depend on the Pose
         * solution, so are not member fields for this Pose. This is because
         * if the Pose solution is changed later by e.g. non-linear iterative
         * optimisation, then these will be wrong.
         */
        
        // First, estimate position vector of each feature in camera frame
        Matrix[] cameraFrameVectors = new Matrix[NCORR];

        // Watch out for critical configurations of known landmarks
        if(!getCameraFrameVectors(matches, cameraFrameVectors, frame.K)) return;

        // Get 3D centroids of features in camera frame.
        Matrix COM_camera = new Matrix(3,1);

        // Sum camera frame position vectors of all tracked features
        for(Matrix v: cameraFrameVectors)
            COM_camera.plusEquals(v);
        COM_camera.timesEquals(1.0/(double)NCORR);
        
        // Now subtract COM vectors from each feature to get positions relative
        // to centroids.
        for(Matrix v: cameraFrameVectors)
            v.minusEquals(COM_camera);
        

        // Tensor product between corresponding vectors in each frames produces
        // a 3x3 matrix with elements consisting of products of the vector
        // components. There is one such matrix for each pair of vectors.
        // Sum these now.
        Matrix M = new Matrix(3,3);
        for(f=0; f<NCORR; f++)
            // Remember to chop off fourth homogenous component of target body
            // frame vector, otherwise operation is illegal.
            M.plusEquals(cameraFrameVectors[f].times(
                     targetFrameCOMVectors[f].getMatrix(0,2,0,0).transpose()));
        
        // Now, use the elements of M to form N, a 4x4 symmetric matrix.
        // Do the algebra in a separate method for clarity here.        
        Matrix N = getN(M);

        // Now, take eigenvalue decomposition of N
        EigenvalueDecomposition evd = new EigenvalueDecomposition(N);
        
        // Eigenvector corresponding to largest eigenvalue is 
        // least-squares solution for rotation
        Matrix eigenvalues = evd.getD();

        // Get column index of most positive eigenvalue
        int cindex = 0;
        double max = eigenvalues.get(0, 0);
        for(int c=1; c<4; c++) 
            if(eigenvalues.get(c, c) > max){
                max = eigenvalues.get(c, c);
                cindex = c;
            }

        // Now construct the rotation Quaternion from eigenvector
        Matrix rot = evd.getV().getMatrix(0, 3, cindex, cindex);

        // Get least squares rotation cam -> world frame
        Quaternion q = new Quaternion(rot.get(0, 0),
                                      rot.get(1, 0),
                                      rot.get(2, 0),
                                      rot.get(3, 0));

        
        
        // Rotate camera frame COM vector to World frame:
        Vector3d com_cam_world = q.rotate(new Vector3d(COM_camera));
        
        // Subtract from world frame COM vector to get position of camera
        // in world frame
        Vector3d t = (new Vector3d(COM_target.getMatrix(0,2,0,0))).minus(com_cam_world);
        
        // Least squares Pose can now be constructed:
        frame.setAttitude(q);
        frame.setPosition(t);        
        
        // Calculate first pose error statistic e1
        //
        // "Mean mis-alignment error"
        //
        // Rotate all centroid-relative world frame position vectors to camera
        // frame, then compare these to direction of corresponding vectors
        // determined from image features.
        Matrix[] targetFrameCOMVectorsR = new Matrix[NCORR];
        for(f=0; f<NCORR; f++)
            // Remember to chop off homogenous fourth component
            targetFrameCOMVectorsR[f] = q.inverse().toMatrix().times(
                                     targetFrameCOMVectors[f].getMatrix(0,2,0,0));
        
        // Sum cross products of corresponding vectors. These are zero for 
        // perfect rotation.
        e1=0;
        //
        for(f=0; f<NCORR; f++){

            // Take cross product of vectors that are now parallel for a 
            // perfect pose.

            // Normalise prior to cross product so that magnitude is sin(angle)
            Matrix a = cameraFrameVectors[f].times(
                                            1.0/cameraFrameVectors[f].normF());
            Matrix b = targetFrameCOMVectorsR[f].times(
                                            1.0/targetFrameCOMVectorsR[f].normF());

            // Take cross product
            Matrix c = new Matrix(new double[][]{{a.get(1,0)*b.get(2,0) -
                                                        a.get(2,0)*b.get(1,0)},
                                                 {a.get(2,0)*b.get(0,0) -
                                                        a.get(0,0)*b.get(2,0)},
                                                 {a.get(0,0)*b.get(1,0) -
                                                      a.get(1,0)*b.get(0,0)}});
            
            // Sum |sin(angle)|^2
            e1 += c.norm2();
        }
        //
        // Divide by the number of features to get error per feature
        e1 /= (double)NCORR;


        // Calculate second pose error statistic e2
        //
        // "Mean re-projection error"
        //
        this.setE2();
        
        // Finally, set validity flag to true to indicate success
        ISVALID = true;

    }



    /**
     * This method implements the algorithm from "Linear N-Point Camera Pose Determination"
     *
     * @param matches
     * @param cameraFrameVectors
     * @param cam_K
     * @return  true if all vectors have been determined, false if any are NaN or infinite.
     */
    public static boolean getCameraFrameVectors(List<Correspondence_2D_3D> matches,
                                                Matrix[] cameraFrameVectors,
                                                Matrix cam_K){

        // How many matches are there?
        int NCORR = cameraFrameVectors.length;

        // Check there are sufficient matches to determine pose in principle
        if(NCORR < 4) return false;
        
        // How many fourth degree polynomials can be constructed for
        // each landmark?
        int NPOLY = (NCORR - 1) * (NCORR - 2) / 2;

        
        System.out.println(NCORR+" matches in gCFV(), "+NPOLY+" polynomials");        
        
        // Track index of primary feature, the one whose position vector
        // is currently being measured.
        int INDEX_J=0;

        // Get inverse camera matrix for de-projecting image coordinates
        // to camera frame unit vectors.
        Matrix invK = cam_K.inverse();

        
        // Loop over all matches
        for (Correspondence_2D_3D j: matches) {
            
            
            // Initialise design matrix. Set minimum number of rows to
            // 5. This essentially uses zero-padding to fill out a
            // rank-deficient matrix so SVD algorithm can be used
            Matrix A = new Matrix(Math.max(NPOLY,5), 5);

            // Track index of next row to be filled
            int r = 0;

            // Turn homogenous image coordinates j into unit vector in
            // camera frame.
            Matrix jx = invK.times(j.p2d.x); // Vector of unit length in z
            jx.timesEquals(1.0/jx.normF());  // Vector of unit length

            // Get landmark position vector in target body frame.
            Matrix J = j.p3d.x;

            
            // Get second landmark I
            for (Correspondence_2D_3D i: matches) {
                
                // If I is the same landmark as J, skip this loop
                if (i.equals(j)) { continue;}
        
                // Turn homogenous image coordinates i into unit vector in
                // camera frame.
                Matrix ix = invK.times(i.p2d.x); // Vector of unit length in z
                ix.timesEquals(1.0/ix.normF());  // Vector of unit length

                // Get landmark position vector in target body frame
                Matrix I = i.p3d.x;


                
                // Get third landmark K. starting index in List for third Match
                // in set of three required for 4th degree polynomial is 
                // element immediately after i.
                for(int INDEX_K = matches.indexOf(i)+1; 
                        INDEX_K < matches.size(); 
                        INDEX_K ++){

                    Correspondence_2D_3D k = matches.get(INDEX_K);
                    
                    
                    // If k is the same landmark as i or j, skip this loop
                    if (k.equals(i) || k.equals(j)) { continue;}
                    
                    // Turn homogenous image coordinates k into unit vector in
                    // camera frame.
                    Matrix kx = invK.times(k.p2d.x); // Vector of unit length in z
                    kx.timesEquals(1.0/kx.normF());  // Vector of unit length

                    // Get landmark position vector in target body frame
                    Matrix K = k.p3d.x;

                    // ix, jx & kx are unit vectors in the camera frame in
                    // the direction of tracked features i,j,k.

                    // Get cosines of internal angles
                    double cos_jk = jx.transpose().times(kx).get(0,0);
                    double cos_ik = ix.transpose().times(kx).get(0,0);
                    double cos_ij = ix.transpose().times(jx).get(0,0);

                    // Get distances between known landmarks from database
                    double d_jk = K.minus(J).normF();
                    double d_ik = K.minus(I).normF();
                    double d_ij = J.minus(I).normF();

                    // Get the fourth degree polynomial for x^2
                    double[] a = getFourthDegreePolynomial(d_ij, d_jk, d_ik, cos_ij, cos_jk, cos_ik);

                    // Load these coefficients into the design matrix
                    for (int c = 0; c < 5; c++) { A.set(r, c, a[c]);}
                        
                    // Increment row index
                    r++;
                    
                }
                
            }

            //+++ Design matrix is now filled for landmark l1 (J). Solve for squared distance x +++//
            //
            //  This is a homogeneous system of equations of the form
            //
            //  {a_00       a_01       a_02       a_03       a_04        {1}     = {0}
            //   a_10       a_11       a_12       a_13       a_14        {x}
            //   a_20       a_21       a_22       a_23       a_24        {x^2}
            //              ...                                          {x^3}
            //   a_(NPOLY)0 a_(NPOLY)1 a_(NPOLY)2 a_(NPOLY)3 a_(NPOLY)4} {x^4}
            //
            //  Ax = 0
            //
            //  Solve this using singular value decomposition.
            //
            // Note that if the landmark correspondances are wrong, a non-trivial solution does not
            // exist and the singular values will all be large. This could be used to identify
            // wrong matches when doing the random assignments.
            //
            //

            // Initialise distance vector
            Matrix x = new Matrix(5, 1);
                
            // Take singular value decomposition of design matrix A
            SingularValueDecomposition svd = new SingularValueDecomposition(A);
                
            // If NPOLY=3, solution lies in 2D right null space of A. Use 
            // 4-point algorithm to apply non-linear constraints and locate
            // solution.
            if (NPOLY == 3){

                // Right singular vectors corresponding to two smallest
                // singular values are the two right-most columns of V,
                // due to way that singular values and vectors are
                // arranged by algorithm.
                Matrix v4 = svd.getV().getMatrix(0, 4, 3, 3);
                Matrix v5 = svd.getV().getMatrix(0, 4, 4, 4);
                //
                // Solutions are given by linear combinations of these two vectors:
                //
                // x = LAMBDA*v4 + RHO*v5  where LAMBDA and RHO are scalars
                //
                // Not all linear combinations are solutions however, because of the
                // constraints between components of x:
                //
                // x[0]*x[3] = x[1]*x[2]
                // x[i]*x[j] = x[k]*x[l]  generally for i+j = k+l
                //
                // Substituting this into equation for x gives
                //
                // b1*LAMBDA^2 + b2*LAMBDA*RHO + b3*RHO^2 = 0
                //
                // where b123 are functions of the components of v4 & v5
                //
                // There are seven such equations from seven distinct combinations
                // of ijkl. Load the corresponding values of b123 into a matrix.
                //

                Matrix b = new Matrix(7, 3);

                int ii = -1, jj = -1, kk = -1, ll = -1;

                for (int row = 0; row < 7; row++) {

                    //+++ Set values of ijkl +++//
                    switch(row){
                        case 0: ii=4; jj=2; kk=3; ll=3; break;
                        case 1: ii=4; jj=1; kk=3; ll=2; break;
                        case 2: ii=4; jj=0; kk=3; ll=1; break;
                        case 3: ii=4; jj=0; kk=2; ll=2; break;
                        case 4: ii=3; jj=1; kk=2; ll=2; break;
                        case 5: ii=3; jj=0; kk=2; ll=1; break;
                        case 6: ii=2; jj=0; kk=1; ll=1; break;
                    }
                    // set b1
                    b.set(row, 0, v4.get(ii,0)*v4.get(jj,0) - v4.get(kk,0)*v4.get(ll,0));
                    // set b2
                    b.set(row, 1, v4.get(ii,0)*v5.get(jj,0) + v5.get(ii,0)*v4.get(jj,0)
                               - (v4.get(kk,0)*v5.get(ll,0) + v5.get(kk,0)*v4.get(ll,0)));
                    // set b3
                    b.set(row, 2, v5.get(ii,0)*v5.get(jj,0) - v5.get(kk,0)*v5.get(ll,0));
                }

                // Now b*y = 0
                //
                // where y = (LAMDBA^2,LAMBDA*RHO,RHO^2)T
                // This is an overdetermined system.
                //
                // Solve for vector y by singular value decomposition.
                //
                SingularValueDecomposition svd_b = new SingularValueDecomposition(b);

                // Desired vector is right singular vector of smallest singular value
                Matrix rightSingular_b = svd_b.getV();

                // Get right singular vector corresponding to smallest
                // singular value. This is the solution for y.        
                Matrix y = rightSingular_b.getMatrix(0, 2, 2, 2);

                // Now solve for ratio LAMBDA/RHO
                double lambdaOverRho1 = y.get(0, 0) / y.get(1, 0);
                double lambdaOverRho2 = y.get(1, 0) / y.get(2, 0);

                double lambdaOverRho = (lambdaOverRho1 + lambdaOverRho2) / 2.0;

                //+++ Now obtain LAMBDA and RHO explicitly using 1 = LAMBDA*v4[0] + RHO*v5[0]
                double RHO = 1.0 / (lambdaOverRho * v4.get(0, 0) + v5.get(0, 0));
                double LAMBDA = 1.0 / (v4.get(0, 0) + v5.get(0, 0) / lambdaOverRho);

                //+++ Get linear combination of v4 & v5 +++//
                x = v4.times(LAMBDA).plus(v5.times(RHO));
                

            }
            
            // With more than four points, solution is obtained directly
            // from right null space of matrix A.
            else {

                // Get the right singular vectors. There are always 5
                Matrix rightSingular = svd.getV();

                x = rightSingular.getMatrix(0, 4, 4, 4);

            }
            
            // Terms in x are coupled; combine 4 measures of distance
            double d = 0;
            int N = 0;

            for (int a = 0; a < 4; a++) {

                double d2_est = x.get(a + 1, 0) / x.get(a, 0);

                // Watch out for negative values and infinities
                if (d2_est > 0 && !Double.isInfinite(d2_est)) {
                    d += Math.sqrt(d2_est);
                    N++;
                }
            }

            // Store distance to this point
            d /= (double) N;

            if (N == 0) {
                // If individual distance measure is NaN, scrap entire pose determination.
                // This occurs for certain critical configurations.                       
                return false;
            }

            // Multiply unit vector in direction of j by this distance
            // to obtain camera frame position vector.                
            cameraFrameVectors[INDEX_J] = jx.times(d);

            // Increment index of feature
            INDEX_J++;

        }  // Close loop over landmark J

        return true;
    }



    /**
     *
     * This method computes the coefficients of the fourth degree polynomial
     * for landmark J. The algebra used to obtain the coefficients uses an
     * original substitution and so does not implement any specific method from
     * the literature, but the coefficients are labelled using the convention
     * from "Linear N-point Camera Pose Determination".
     *
     *
     *
     * @param d_ij     Distance between landmarks I & J
     * @param d_jk     Distance between landmarks J & K
     * @param d_ik     Distance between landmarks I & K
     * @param cos_ij   Cosine of angle between points i & j in camera frame
     * @param cos_jk   Cosine of angle between points j & k in camera frame
     * @param cos_ik   Cosine of angle between points i & k in camera frame
     * @return
     */
    public static double[] getFourthDegreePolynomial(double d_ij, double d_jk, double d_ik,
                                                     double cos_ij, double cos_jk, double cos_ik){

        double sin_jk = Math.sqrt(1.0 - cos_jk*cos_jk);

        //+++ h(xj,xk) +++//
        //
        // This is the function obtained by eliminating xi between f(xi,xk) and f(xi,xj)
        // where f(ij) = xi^2 + xj^2 - 2.0 * xi * xj * cos(theta_ij) - d_ij^2
        // i.e. the cosine rule.
        //
        //
        // h(xj,xk) = a*xj^4 + b*xk^4 + c*xj^2 + d*xk^2 + e*xj*xk^3 + f*xj^2*xk^2 + g*xj^3*xk + h*xj*xk + i
        //
        double a = 1;
        double b = 1;
        double c = (2 * (d_ik * d_ik - d_ij * d_ij) - 4 * cos_ij * cos_ij * d_ik * d_ik);
        double d = (2 * (d_ij * d_ij - d_ik * d_ik) - 4 * cos_ik * cos_ik * d_ij * d_ij);
        double e = (-4 * cos_ij * cos_ik);
        double f = (4 * (cos_ij * cos_ij + cos_ik * cos_ik) - 2);
        double g = (-4 * cos_ij * cos_ik);
        double h = 4 * cos_ij * cos_ik * (d_ik * d_ik + d_ij * d_ij);
        double i = (d_ik * d_ik - d_ij * d_ij) * (d_ik * d_ik - d_ij * d_ij);
        //
        //+++ Now eliminate xk from h(xj,xk) +++//
        //
        //  Use internal angle alpha between xj and d_jk to get:
        //
        //  xk = xj*cos_jk + d_jk*cos(alpha)
        //
        //  Note that cos(alpha)^2 = 1 - (xj^2 * sin_jk^2)/(d_jk^2)
        //
        //  Substitute this into h(xj,xk) to get a polynomial in terms of xj and
        //  angle alpha:
        //

        double h1 = (a + b * (Math.pow(cos_jk, 4) + Math.pow(sin_jk, 4) - 6 * cos_jk * cos_jk * sin_jk * sin_jk) + e * (Math.pow(cos_jk, 3) - 3 * cos_jk * sin_jk * sin_jk) + f * (cos_jk * cos_jk - sin_jk * sin_jk) + g * cos_jk);
        double h2 = (4 * b * d_jk * (Math.pow(cos_jk, 3) - cos_jk * sin_jk * sin_jk) + e * d_jk * (3 * cos_jk * cos_jk - sin_jk * sin_jk) + 2 * f * d_jk * cos_jk + g * d_jk);
        double h3 = (c + d * (cos_jk * cos_jk - sin_jk * sin_jk) + 2 * b * d_jk * d_jk * (3 * cos_jk * cos_jk - sin_jk * sin_jk) + 3 * e * d_jk * d_jk * cos_jk + f * d_jk * d_jk + h * cos_jk);
        double h4 = (4 * b * d_jk * d_jk * d_jk * cos_jk + 2 * d * d_jk * cos_jk + e * d_jk * d_jk * d_jk + h * d_jk);
        double h5 = (d * d_jk * d_jk + b * d_jk * d_jk * d_jk * d_jk + i);
        
        // This polynomial is of the form h(xj) = h1*xj^4 + cos(alpha)*h2*xj^3 + h3*xj^2 + cos(alpha)*h4*xj + h5 = 0
        //
        // Rearrange like this:
        //
        // h1*xj^4 + h3*xj^2 + h5 = cos(alpha)*(h2*xj^3 + h4*xj)
        //
        // Then square both sides to get rid of cos(alpha) and obtain the desired
        // 8th degree polynomial in xj
        //
        //
        //+++ Get coefficients of 8th degree polynomial +++//
        //
        // Choose syntax to match that in "Linear N-point Camera Pose Determination"
        //
        //
        double a5 = h1 * h1 + h2 * h2 * sin_jk * sin_jk / (d_jk * d_jk);
        double a4 = 2 * h1 * h3 - h2 * h2 + 2 * h2 * h4 * sin_jk * sin_jk / (d_jk * d_jk);
        double a3 = 2 * h1 * h5 + h3 * h3 - 2 * h2 * h4 + h4 * h4 * sin_jk * sin_jk / (d_jk * d_jk);
        double a2 = 2 * h3 * h5 - h4 * h4;
        double a1 = h5 * h5;

        return new double[]{a1,a2,a3,a4,a5};

    }


    public static Matrix getN(Matrix M){

        Matrix N = new Matrix(4,4);

        N.set(0, 0,  M.get(0,0) + M.get(1,1) + M.get(2,2));  // Diagonal elements
        N.set(1, 1,  M.get(0,0) - M.get(1,1) - M.get(2,2));
        N.set(2, 2, -M.get(0,0) + M.get(1,1) - M.get(2,2));
        N.set(3, 3, -M.get(0,0) - M.get(1,1) + M.get(2,2));

        N.set(0, 1,  M.get(1,2) - M.get(2,1));    // Remaining elements of first row
        N.set(0, 2,  M.get(2,0) - M.get(0,2));
        N.set(0, 3,  M.get(0,1) - M.get(1,0));
        N.set(1, 2,  M.get(0,1) + M.get(1,0));    // Remaining elements of second row
        N.set(1, 3,  M.get(2,0) + M.get(0,2));
        N.set(2, 3,  M.get(1,2) + M.get(2,1));    // Remaining elements of third row

        // N is symmetric. Copy upper right off-diagonal elements to lower left
        N.set(1, 0, N.get(0, 1));
        N.set(2, 0, N.get(0, 2));
        N.set(2, 1, N.get(1, 2));
        N.set(3, 0, N.get(0, 3));
        N.set(3, 1, N.get(1, 3));
        N.set(3, 2, N.get(2, 3));

        return N;
    }

    
    /**
     * Use current pose solution and Matches set to determine mean 
     * re-projection error.
     */
    public final void setE2(){
 
        e2=0.0;
        //
        for(Correspondence_2D_3D M: frame.correspondences){
           
            // Project coordinates of 3D point onto image plane, along with
            // position covariance matrix.
            Point2D p2d = frame.project(M.p3d);
            
            // Express distance between projected 3D point and observed 2D
            // point in terms of numbers of sigmas, where distribution is
            // that of projected 3D point.
            e2 += p2d.getSigmasFromMean2(M.p2d.x.get(0,0), M.p2d.x.get(1,0));
           
        }
        //
        // Divide by the number of features to get offset per feature.
        e2 /= (double)frame.correspondences.size();
        
    }

       
    
    
    
}
