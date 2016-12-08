package numeric.vision;

import Jama.Matrix;
import Jama.SingularValueDecomposition;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * This class represents Fundamental Matrices. It contains constructors that
 * can be used to measure F from point correspondences and from known
 * camera matrices K[R|t].
 *
 * To do, to improve estimate of F from point correspondences:
 *
 * 1) Use 7 point algorithm in RANSAC loop to improve outlier rejection.
 * 2) Use non-linear optimisation algorithm to estimate F from all inliers.
 *
 * @author nickrowell
 */
public class FundamentalMatrix{

    // Class members
    //
    /** Consensus level for RANSAC hypothesis acceptance */
    static double THRESHOLD = 0.8;
    /** Outlier threshold (pixels offset perpendicular to epipolar line) */
    static double OUTLIER = 5.0;
    /** Max number of RANSAC loops */
    static int MAXLOOP = 50000;

    // Setters for class members
    public static void setThreshold(double thresh){ THRESHOLD = thresh;}
    public static void setOutlier(double out){ OUTLIER = out;}
    public static void setMaxloop(int loop){ MAXLOOP = loop;}

    // Instance members
    //
    /** Fundamental matrix */
    public Matrix F = new Matrix(3,3);
    /** Epipoles in homogenous coordinates in image */
    public Matrix e_prime = new Matrix(3,1), e = new Matrix(3,1);

    /** Condition number of design matrix used to solve for F */
    public double COND=-1;

    // Data concerning F when measured from point correspondences
    //
    /** Point correspondences */
    public List<Correspondence_2D_2D_3D> points = new ArrayList<Correspondence_2D_2D_3D>();
    /** Inlying point correspondences */
    public List<Correspondence_2D_2D_3D> inliers = new ArrayList<Correspondence_2D_2D_3D>();
    /** Outlying point correspondences */
    public List<Correspondence_2D_2D_3D> outliers = new ArrayList<Correspondence_2D_2D_3D>();
    /** Transformation matrices used to normalise coordinates */
    Matrix T, T_prime;

    /** Default constructor */
    public FundamentalMatrix(){}


    /** Constructor for F estimated from point correspondences */
    public FundamentalMatrix(List<Correspondence_2D_2D_3D> sps)
    {

        // Check if there are sufficient point correspondences
        if(sps.size()<8)
            throw new RuntimeException("Cannot estimate E from "+sps.size()+
                                       " point correspondences");

        // Copy Correspondences by reference to internal data
        for(int p=0; p<sps.size(); p++)
            points.add(sps.get(p));
        

        // Normalise points and record transformation matrices
        this.normalise(points);

        // Use RANSAC to divide all point correspondences between inlier
        // and outlier lists.
        this.RANSAC();

        // Now normalise points again, but this time only inliers. This
        // corrects for any adverse effects that inclusion of outliers
        // might have had on normalisation of entire point list.
        // After call, normalisation matrices T,T_prime represent 
        // normalisation applied to inliers alone.
        this.normalise(inliers);

        /**
         * At this point, algorithm could simply use full list of inliers
         * to calculate F. However, it may be desirable to use a restricted
         * set of the most linearly independent point correspondences. This
         * is possible using the QR decomposition with column pivoting. The
         * permutation matrix is used to arrange the design matrix rows into
         * into order of linear independency, then the permuted design matrix
         * is top-sliced to get a selection of the most linearly independent
         * correspondences.
         */
        
        // Turn inliers list into a Matrix with no fewer than 9 rows,
        // with zero-padding to fill out to 9 rows if necessary.
        Matrix A = getDesignMatrix(inliers,9);

        // Record condition number of design matrix
        COND = A.cond();

        // The following section of code results in a solution for F that
        // differs from Matlab code that I've validated my code against. It's
        // not clear if the differences are an improvement or not - would have
        // to test this by checking effect it has on accuracy of recovered
        // trajectory. I've removed this for now, but should investigate it
        // further when I'm at the point of measuring S/C trajectories.

//        // Get permutation matrix via QR decomposition. Transpose of A has
//        // point correspondences in column vectors.
//        Matrix[] PQR = QR.qr(A.transpose());
//
//        // Multiply transposed A by permutation matrix and transpose back.
//        // This sorts rows into order of linear dependency.
//        A = (A.transpose().times(PQR[0])).transpose();
//
//        // Get first 20% of rows of A, but not less than 9
//        int    NROWS = Math.max(9,A.getRowDimension()/5);
//        Matrix A_SUB = A.getMatrix(0, NROWS-1, 0, 8);
//
//        // Use inliers to calculate F with singularity constraint enforced
//        F = computeF(A_SUB);

        // Use inliers to calculate F with singularity constraint enforced
        F = computeF(A);

        // De-normalise F
        F = T_prime.transpose().times(F.times(T));

        
        // Set the epipoles
        setEpipoles();
        

        // Finished!

    }

    /** 
     * Constructor for F obtained directly from camera matrices.
     */
    public FundamentalMatrix(Matrix p, Matrix p_prime)
    {
        // Use camera matrices to set F and epipoles
        setFundFromCameras(p,p_prime);        
    }


    /**
     * Convert two camera matrices into a FundamentalMatrix.
     * Follows calculation in Multiple View Geometry pp.412 section 17.3
     */
    public final void setFundFromCameras(Matrix P, Matrix P_prime){

        // F obtained from determinants of sub-matrices of P, as in
        // equation 17.3 (pp412) of Multiple View Geometry
        Matrix A1 = P.getMatrix(new int[]{1,2},0,3);
        Matrix A2 = P.getMatrix(new int[]{0,2},0,3);
        Matrix A3 = P.getMatrix(new int[]{0,1},0,3);

        Matrix B1 = P_prime.getMatrix(new int[]{1,2},0,3);
        Matrix B2 = P_prime.getMatrix(new int[]{0,2},0,3);
        Matrix B3 = P_prime.getMatrix(new int[]{0,1},0,3);

        Matrix[] A = new Matrix[]{A1,A2,A3};
        Matrix[] B = new Matrix[]{B1,B2,B3};

        // Matrix to hold temporary values
        Matrix AB = new Matrix(4,4);
        // This will be populated with fundamental matrix elements
        F = new Matrix(3,3);

        for(int r=0; r<3; r++)
            for(int c=0; c<3; c++){

                AB.setMatrix(new int[]{0,1}, 0, 3, A[c]);
                AB.setMatrix(new int[]{2,3}, 0, 3, B[r]);

                F.set(r, c, AB.det()*Math.pow(-1,r+c));
            }
     
        // Set the epipoles for the new Fundamental Matrix.
        setEpipoles();
        

    }

    public final void setEpipoles(){
        // Now set the epipoles:
        SingularValueDecomposition svd = F.svd();
        // right singular vector, i.e. null space of F; F*e1 = (0,0,0)^T
        e = svd.getV().getMatrix(0, 2, 2, 2);
        // left singular vector, i.e. null space of F^T; e2^T*F = (0,0,0)
        e_prime = svd.getU().getMatrix(0, 2, 2, 2);

        // Scale homogenous epipoles to pixel units in first two components
        e.timesEquals(1.0/e.get(2,0));
        e_prime.timesEquals(1.0/e_prime.get(2,0));
        
    }
    
    
    /**
     * Split entire list of point correspondences into inliers and outliers
     * based on RANSAC scheme, using 8 random point correspondences to
     * initialise each F hypothesis.
     */
    private void RANSAC(){

        // Count number of loops
        int N=0;
        // Current consensus level
        double consensus=-1;
        // Active consensus set, reset with each loop
        List<Correspondence_2D_2D_3D> activeSet = new ArrayList<Correspondence_2D_2D_3D>();
        // From all sets tested, record the running best
        List<Correspondence_2D_2D_3D> largestSet = new ArrayList<Correspondence_2D_2D_3D>();

        while(consensus < THRESHOLD && N++ < MAXLOOP){

            // Clear active consensus set
            activeSet.clear();
            
            // Add first 8 points to consensus set
            for(int p=0; p<8; p++)
                activeSet.add(points.get(p));
            
            // Calculate normalised F from 8 point correspondences
            F = computeF(activeSet);
            
            // Now expand consensus set using F to assess each point
            for(int p=8; p<points.size(); p++)
                if(getPointResidual(points.get(p))  < OUTLIER)
                    activeSet.add(points.get(p));
            
            // Check size of consensus set for this F hypothesis
            consensus = (double)activeSet.size()/(double)points.size();
            
            // If active consensus set is larger than running best, reset
            // running best and copy active set members to it
            if(activeSet.size() > largestSet.size() || consensus > THRESHOLD){
                largestSet.clear();
                for(int p=0; p<activeSet.size(); p++)
                    largestSet.add(activeSet.get(p));
            }
            
            // Randomize order of points ready for next loop
            Collections.shuffle(points);
            
        }
        
        System.out.println("Fundamental Matrix: RANSAC used "+N+" random trials.");
        
        // A sufficiently large consensus set has been found, or max number of
        // loops has been reached.
        // Now divide all points correspondences among inlier and outlier lists
        for(int p=0; p<points.size(); p++)
            if(largestSet.contains(points.get(p)))
                inliers.add(points.get(p));
            else
                outliers.add(points.get(p));
        
    }

    /**
     * Computes the fundamental matrix given a set of point
     * correspondences. Uses the 8 point algorithm, and enforces conditions on
     * the singular values of F that the smallest one is zero
     */
    private static Matrix computeF(List<Correspondence_2D_2D_3D> points){
        return computeF(getDesignMatrix(points, 9));
    }


    /**
     * Compute the fundamental matrix for a set of point correspondences.
     * 
     * This is the "8-point algorithm"
     * 
     */
    private static Matrix computeF(Matrix A){

        // Test shape of A
        if(A.getRowDimension() < A.getColumnDimension())
            throw new RuntimeException("Matrix A has too few rows for SVD");

        // Take singular value decomposition of A
        SingularValueDecomposition svd = A.svd();

        // Solution for f is right singular vector corresponding to
        // smallest singular value:
        Matrix f = svd.getV().getMatrix(0, 8, 8, 8);
        //
        // Now build F from entries of f
        Matrix F = new Matrix(new double[][]{{f.get(0, 0), f.get(1, 0), f.get(2, 0)},
                                             {f.get(3, 0), f.get(4, 0), f.get(5, 0)},
                                             {f.get(6, 0), f.get(7, 0), f.get(8, 0)}});
        //
        // Enforce singularity of fundamental matrix
        SingularValueDecomposition svdF = F.svd();
        //
        // Get singular values...
        Matrix S = svdF.getS();
        //
        // Set the smallest one to zero...
        S.set(2, 2, 0.);
        //
        // Rebuild singular F. This stage succesfully reduces the rank
        // of F from 3 to 2.
        Matrix F_SING = svdF.getU().times(S.times(svdF.getV().transpose()));

        return F_SING;

    }


    /**
     * Get matrix used to determine F. If there are fewer points in list
     * than MINSIZE, the returned matrix is zero-padded so that the
     * singular value decomposition algorithm always works.
     */
    public static Matrix getDesignMatrix(List<Correspondence_2D_2D_3D> points,
                                         int MINSIZE){

        // Matrix A should be padded with rows of zeros so that it is
        // always at least square. This way, the singular value decomposition
        // algorithm can be used.
        int NROWS = Math.max(MINSIZE, points.size());

        // Copy data to array
        double[][] a= new double[NROWS][9];
        // insert data and leave any remaining rows at zero
        for (int r = 0; r < points.size(); r++)
            a[r] = points.get(r).getRowEntry();

        // Build design matrix A
        return new Matrix(a);

    }

    /** Get number of inliers found by RANSAC algorithm */
    public int getInliersN(){ return inliers.size();}
    /** Get number of outliers found by RANSAC algorithm */
    public int getOutliersN(){ return outliers.size();} 
    
    
    
    
    /**
     * Calculate residual used to decide whether a point correspondence is
     * and inlier or outlier with respect to a given estimated fundamental
     * matrix. The residual is equal to the perpendicular distance from each
     * point to its corresponding epipolar line, squared then summed across
     * both images.
     *
     * This is used during RANSAC to test for outliers. The calculation uses
     * normalised point coordinates, and F is normalised at the point this
     * method is called.
     */
    public double getPointResidual(Correspondence_2D_2D_3D sp){

        // Get the normalised epipolar lines in each image
        Matrix epl       = F.transpose().times(sp.p2d_prime.x_norm);
        Matrix epl_prime = F.times(sp.p2d.x_norm);

        // Get scaled distance between each point and its epipolar line
        double d       = (sp.p2d.x_norm.transpose().times(epl)).get(0, 0);
        double d_prime = (sp.p2d_prime.x_norm.transpose().times(epl_prime)).get(0, 0);

        // Now scale perpendicular distances to pixel units
        double a2 = epl.get(0, 0);
        double b2 = epl.get(1, 0);

        double a1 = epl_prime.get(0, 0);
        double b1 = epl_prime.get(1, 0);

        d       *= 1.0/Math.sqrt(a2*a2 + b2*b2);
        d_prime *= 1.0/Math.sqrt(a1*a1 + b1*b1);

        // Now inflate d & d_prime by scale factor
        d       /= T.get(0,0);
        d_prime /= T_prime.get(0,0);

        // residual
        return (d*d + d_prime*d_prime);

    }


    /**
     * Use Hartley's isotropic normalisation method to shift and scale
     * point correspondence coordinates. This improves numerical stability
     * of solution.
     * 
     */
    private void normalise(List<Correspondence_2D_2D_3D> pts){

        // Intermediate
        double[][] x   = new double[pts.size()][2];
        double[][] x_p = new double[pts.size()][2];

        // Get mean position of points in each frame
        double x_COM=0, y_COM=0, x_p_COM=0, y_p_COM=0;
        for(int i=0; i<pts.size(); i++){
            x_COM   += pts.get(i).p2d.x.get(0,0);
            y_COM   += pts.get(i).p2d.x.get(1,0);
            x_p_COM += pts.get(i).p2d_prime.x.get(0,0);
            y_p_COM += pts.get(i).p2d_prime.x.get(1,0);
        }
        
        x_COM   /= (double)pts.size();
        y_COM   /= (double)pts.size();
        x_p_COM /= (double)pts.size();
        y_p_COM /= (double)pts.size();

        // Get centroid-relative coordinates for each point
        for(int i=0; i<pts.size(); i++){
            x[i][0]   = pts.get(i).p2d.x.get(0,0) - x_COM;
            x[i][1]   = pts.get(i).p2d.x.get(1,0) - y_COM;
            x_p[i][0] = pts.get(i).p2d_prime.x.get(0,0) - x_p_COM;
            x_p[i][1] = pts.get(i).p2d_prime.x.get(1,0) - y_p_COM;
        }

        // Get mean distance from origin in each frame:
        //
        // 1) Get sum distances from origin in each frame
        double D=0, D_p=0;
        for(int i=0; i<pts.size(); i++){
            D   += Math.sqrt(x[i][0] * x[i][0] + x[i][1] * x[i][1]);
            D_p += Math.sqrt(x_p[i][0] * x_p[i][0] + x_p[i][1] * x_p[i][1]);
        }
        // 2) Take average
        D   /= (double)pts.size();
        D_p /= (double)pts.size();

        // Scale factors
        double s   = Math.sqrt(2.0)/D;
        double s_p = Math.sqrt(2.0)/D_p;

        // Set transformation matrices
        T = new Matrix(new double[][]{{s,0,-s*x_COM},
                                      {0,s,-s*y_COM},
                                      {0,0,1}});

        T_prime = new Matrix(new double[][]{{s_p,0,-s_p*x_p_COM},
                                            {0,s_p,-s_p*y_p_COM},
                                            {0,0,1}});

        // Set normalised point correspondences
        for(int i=0; i<pts.size(); i++){
            pts.get(i).p2d.x_norm       = T.times(pts.get(i).p2d.x);
            pts.get(i).p2d_prime.x_norm = T_prime.times(pts.get(i).p2d_prime.x);
        }

    }




    @Override
    public String toString(){

        String OUT = "";

        // Print outlying points
        OUT += "Outliers\nx = [";
        for(int OL=0; OL<outliers.size(); OL++){
            Matrix x  = outliers.get(OL).p2d.x;
            OUT += "\n" + x.get(0,0) + " " + x.get(1,0) + " " + x.get(2,0);
        }
        OUT += "];\n\nxp = [";
        for(int OL=0; OL<outliers.size(); OL++){
            Matrix xp = outliers.get(OL).p2d_prime.x;
            OUT += "\n" + xp.get(0,0) + " " + xp.get(1,0) + " " + xp.get(2,0);
        }
        OUT += "];";

        // Print inlying points
        OUT += "\n\nInliers\nx = [";
        for(int IL=0; IL<inliers.size(); IL++){
            Matrix x  = inliers.get(IL).p2d.x;
            OUT += "\n" + x.get(0,0) + " " + x.get(1,0) + " " + x.get(2,0);
        }
        OUT += "];\n\nxp = [";
        for(int IL=0; IL<inliers.size(); IL++){
            Matrix xp = inliers.get(IL).p2d_prime.x;
            OUT += "\n" + xp.get(0,0) + " " + xp.get(1,0) + " " + xp.get(2,0);
        }
        OUT += "];";

        DecimalFormat mat = new DecimalFormat("0.000E0");

        // Print fundamental matrix
        OUT += "\n\nFundamental matrix F = \n"+mat.format(F.get(0, 0))+" "+
                mat.format(F.get(0, 1))+" "+" "+mat.format(F.get(0, 2))+"\n"+
                mat.format(F.get(1, 0))+" "+mat.format(F.get(1, 1))+" "+
                mat.format(F.get(1, 2))+"\n"+mat.format(F.get(2, 0))+" "+
                mat.format(F.get(2, 1))+" "+mat.format(F.get(2, 2));

        OUT += "\n\ne = "+e.get(0,0)+" "+e.get(1,0)+" "+e.get(2,0);
        OUT += "\n\nep = "+e_prime.get(0,0)+" "+e_prime.get(1,0)+" "
                +e_prime.get(2,0);


        return OUT;
    }



}

