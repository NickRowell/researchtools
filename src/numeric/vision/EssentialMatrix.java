package numeric.vision;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

import java.util.List;

import numeric.dynamics.Pose;
import numeric.geom.dim3.Quaternion;
import numeric.minimisation.nllsq.algo.LevenbergMarquardt;

/**
 * Class represents the Essential Matrix of projective geometry and provides a
 * range of methods for computing and handling these objects.
 * 
 * To do: Tidy up setCameras() methods. There are way too many. In practise
 *        only setCameras(P,T) and setCameras(P,T,List<Correspondence>) will
 *        be used. P=camera matrix for un-primed frame, T=translation magnitude
 * 
 * @author nickrowell
 */


// Parameters of EssentialMatrix are 3D locations of point correspondences 
// (contained in FundamentalMatrix class), and normalised camera matrices
// that provide estimated camera motion. Data is the image plane coordinates
// of point correspondences.
//
public class EssentialMatrix extends LevenbergMarquardt {

    /** Matrix used in calculation of camera matrices from E */
    static Matrix W = new Matrix(new double[][]{{ 0,-1, 0},
                                                { 1, 0, 0},
                                                { 0, 0, 1}});

    // Instance members
    /**
     * Fundamental matrix is a generalised epipolar geometry entity that cannot
     * be used to recover full camera matrices but is an important first step
     * in computing Euclidean essential matrix.
     */
    public FundamentalMatrix fund;
    
    /** Essential matrix */
    public Matrix E = new Matrix(3,3);
    /** Camera calibration matrix for unprimed image */
    public Matrix K = null;
    /** Camera calibration matrix for primed image */
    public Matrix K_prime = null;
    /** Second (later) normalised camera matrix used in reconstruction */
    public Matrix P_prime = null;
    /** First (earlier) normalised camera matrix used in reconstruction */
    public Matrix P = null;

    // Use quaternion & vector representation of P' to optimize E solution
    // with gold standard algorithm
    
    /** Quaternion to represent orientation of P' */
    private Quaternion ROT = new Quaternion();
    /** Vector to represent translation of P' */
    private Matrix TRANS = new Matrix(3,1);       
    
    /** Default constructor */
    public EssentialMatrix(){}


    /** 
     * Use 8-point algorithm to calculate Fundamental matrix, then transform
     * it to the Essential matrix using the camera calibration matrices
     * provided. The RANSAC algorithm is used to divide the input points
     * into inliers and outliers.
     * 
     * After calling the constructor, call E.setCameras() to obtain a linear
     * estimate of the camera positions and orientations in each frame.
     * Then call E.triangulateInliers() to get linear estimates of the 3D
     * positions of each inlying point.
     * 
     * After that, a call to E.goldStandardAlgorithm() will use Levenberg-
     * Marquardt algorithm to improve solution.
     * 
     */
    public EssentialMatrix(List<Correspondence_2D_2D_3D> sps,
                           Matrix k,
                           Matrix k_prime)
    {
        // Call constructor of superclass FundamentalMatrix. This solves
        // projective geometry of image pair, and divides point correspondences
        // into inliers and outliers.
        fund = new FundamentalMatrix(sps);
        
        // Take copies of camera calibration matrices
        K = k.copy();
        K_prime = k_prime.copy();

        // Convert fundamental matrix to essential matrix using equation 9.12
        // page 257 of multiple view geometry:
        E = K_prime.transpose().times(fund.F).times(K);

        // Now enforce condition on E that its two non-zero singular values
        // are equal:
        SingularValueDecomposition svd = E.svd();
        // Get Matrices U, S, V
        Matrix U = svd.getU();
        Matrix S = svd.getS();
        Matrix V = svd.getV();
        // Get average of singular values:
        double AVE_S = (S.get(0, 0) + S.get(1, 1))/2.0;
        // Set singular values equal to this number
        S.set(0, 0, AVE_S);
        S.set(1, 1, AVE_S);
        S.set(2, 2, 0);
        // Rebuild E:
        E = U.times(S).times(V.transpose());

    }

    /**
     * Construct E directly from (un-normalised) camera
     * matrices, i.e. p = K[R|t]
     */
    public EssentialMatrix(Matrix p,
                           Matrix p_prime,
                           Matrix k,
                           Matrix k_prime)
    {
        // Call constructor of superclass FundamentalMatrix. This solves
        // projective geometry of image pair.
        fund = new FundamentalMatrix(p,p_prime);

        // Take copies of camera calibration matrices
        K = k.copy();
        K_prime = k_prime.copy();

        // Set normalised camera matrices from full camera matrices
        P       = K.inverse().times(p);
        P_prime = K_prime.inverse().times(p_prime);
        
        // Set values of ROT and TRANS from current P_prime matrix
        TRANS = P_prime.getMatrix(new int[]{0,1,2}, new int[]{3});
        ROT = new Quaternion(P_prime.getMatrix(new int[]{0,1,2}, new int[]{0,1,2}));
        
        // Convert fundamental matrix to essential matrix using equation 9.12
        // page 257 of multiple view geometry:
        E = K_prime.transpose().times(fund.F).times(K);

        // Now enforce condition on E that its two non-zero singular values
        // are equal:
        SingularValueDecomposition svd = E.svd();
        // Get atrices U, S, V
        Matrix U = svd.getU();
        Matrix S = svd.getS();
        Matrix V = svd.getV();
        // Get average of singular values:
        double AVE_S = (S.get(0, 0) + S.get(1, 1))/2.0;
        // Set singular values equal to this number
        S.set(0, 0, AVE_S);
        S.set(1, 1, AVE_S);
        S.set(2, 2, 0);
        // Rebuild E:
        E = U.times(S).times(V.transpose());

    }
    
    
    /** Get number of inliers found by RANSAC algorithm */
    public int getInliersN(){ return fund.getInliersN();}
    /** Get number of outliers found by RANSAC algorithm */
    public int getOutliersN(){ return fund.getOutliersN();}    
    

    ///////////////////////////////////////////////////////////////////////
    //
    // There now follows a bunch of methods for manipulating camera
    // matrices given knowledge of calibration matrix K.
    //
    ///////////////////////////////////////////////////////////////////////


    /**
     * Scale-free version of setCameras(Translation magnitude)
     */
    public void setCameras(){ setCameras(1.0, fund.inliers);}
    /**
     * Scale-free version of setCameras(Translation magnitude)
     */
    public void setCameras(double T_MAG){ setCameras(T_MAG, fund.inliers);}
    /**
     * Scale-free version of setCameras(Translation magnitude) with alternative
     * list of points used to locate correct solution. Used in testing.
     */
    public void setCameras(List<Correspondence_2D_2D_3D> points){
        setCameras(1.0, points);
    }
    
    /**
     * Set two consistent camera matrices from an EssentialMatrix.
     *
     * Follows calculation outlined in Multiple View Geometry pp.258-60
     *
     * P is set to [I|0]. P_prime is determined from E and triangulation.
     *
     */
    public void setCameras(double T_MAG, List<Correspondence_2D_2D_3D> points){

        // Initialise un-primed (earlier) camera matrix to [I|0]
        P = Matrix.identity(3, 4);

        // Take SVD of E:
        SingularValueDecomposition svd = E.svd();

        // Translation: get final column of matrix U (left null vectors)
        Matrix u3 = svd.getU().getMatrix(new int[]{0,1,2}, new int[]{2});

        // Scale translation
        u3.timesEquals(T_MAG);

        // Rotation: get two variations on first 3x3 sub-matrix
        Matrix UWVT  = svd.getU().times(W).times(svd.getV().transpose());
        Matrix UWTVT = svd.getU().times(W.transpose()).times(svd.getV().transpose());
        
        // If rotation is in left-handed set (det(R) = -1) this is a rotation
        // and inversion. Multiply by -1 to remove inversion:
        if(UWVT.det()<0){
            UWVT.timesEquals(-1);
        }
        if(UWTVT.det()<0){
            UWTVT.timesEquals(-1);
        }

        // Now build four possible P' solutions
        Matrix[] P_SOLS = {new Matrix(3,4),
                           new Matrix(3,4),
                           new Matrix(3,4),
                           new Matrix(3,4)};

        P_SOLS[0].setMatrix(new int[]{0,1,2}, new int[]{0,1,2}, UWVT);
        P_SOLS[0].setMatrix(new int[]{0,1,2}, new int[]{3}, u3);

        P_SOLS[1].setMatrix(new int[]{0,1,2}, new int[]{0,1,2}, UWVT);
        P_SOLS[1].setMatrix(new int[]{0,1,2}, new int[]{3}, u3.times(-1));

        P_SOLS[2].setMatrix(new int[]{0,1,2}, new int[]{0,1,2}, UWTVT);
        P_SOLS[2].setMatrix(new int[]{0,1,2}, new int[]{3}, u3);

        P_SOLS[3].setMatrix(new int[]{0,1,2}, new int[]{0,1,2}, UWTVT);
        P_SOLS[3].setMatrix(new int[]{0,1,2}, new int[]{3}, u3.times(-1));

        // Get solution from 4 possible camera matrices by a voting
        // technique using all inlying points.
        int SOL = getCorrectSolution(P_SOLS, points);

        P_prime = P_SOLS[SOL].copy();
                
        // Set values of ROT and TRANS from current P_prime matrix
        TRANS = P_prime.getMatrix(new int[]{0,1,2}, new int[]{3});
        ROT = new Quaternion(P_prime.getMatrix(new int[]{0,1,2}, new int[]{0,1,2}));
        
    }


    /**
     * Scale-free version of setCamera(Matrix p, double T_MAG), using this
     * instances own point correspondences to find correct solution.
     */
    public void setCameras(Matrix p){ setCameras(p,1.0, fund.inliers);}

    /**
     * Scale-free version of setCamera(Matrix p, double T_MAG) with alternative
     * list of points used to locate correct solution. Used in testing.
     */
    public void setCameras(Matrix p, List<Correspondence_2D_2D_3D> points){
        setCameras(p,1.0, points);
    }


    public void setCameras(Matrix p, double T_MAG){
        setCameras(p,T_MAG, fund.points);
    }

    /**
     * Set consistent pair of camera matrices, with P (later frame)
     * set equal to the input matrix.
     *
     * @param p Normalised camera matrix. E.P is set equal to this matrix, with
     *          E.P_prime set according to transformation between frames
     *          encoded by Essential Matrix.
     * 
     */
    public void setCameras(Matrix p, double T_MAG, List<Correspondence_2D_2D_3D> points){
        
        // Check that input rotation matrix is orthonormal
        if(Math.abs(p.getMatrix(new int[]{0,1,2}, new int[]{0,1,2}).det()-1)>1E-9){
            System.err.println("setCameras(): Input camera matrix not orthonormal");
            Quaternion Q_ORTH = new Quaternion(p.getMatrix(new int[]{0,1,2}, new int[]{0,1,2}));
            Q_ORTH.normalise();
            p.setMatrix(new int[]{0,1,2}, new int[]{0,1,2}, Q_ORTH.toMatrix());
        }
        
        // p = [R|t] where t is position vector of world origin in P frame,
        // and R is rotation matrix that transforms position vectors from
        // world to P frame according to X_p = R * X_w + t

        // First, set default camera matrices with P = [1|0]
        setCameras(T_MAG, points);

        // Now transform both normalised camera matrices so that the unprimed
        // matrix equals the input matrix

        // Extract rotation and translation from p
        Matrix R0 = p.getMatrix(new int[]{0,1,2}, new int[]{0,1,2});
        Matrix t0 = p.getMatrix(new int[]{0,1,2}, new int[]{3});
        // Set P quantities equal to these
        P.setMatrix(new int[]{0,1,2}, new int[]{0,1,2}, R0);
        P.setMatrix(new int[]{0,1,2}, new int[]{3}, t0);

        // Now extract rotation & translation between default camera frames
        Matrix R = P_prime.getMatrix(new int[]{0,1,2}, new int[]{0,1,2});
        Matrix t = P_prime.getMatrix(new int[]{0,1,2}, new int[]{3});

        // Transform to frame of first camera matrix:
        Matrix R_prime = R.times(R0);
        Matrix T_prime = (R.times(t0)).plus(t);
        
        // Build new camera matrices:
        //
        // P_prime fixed by chaining together inter-frame transformation and
        // input world-P transformation:
        P_prime.setMatrix(new int[]{0,1,2}, new int[]{0,1,2}, R_prime);
        P_prime.setMatrix(new int[]{0,1,2}, new int[]{3}, T_prime);

        // Set values of ROT and TRANS from current P_prime matrix
        TRANS = P_prime.getMatrix(new int[]{0,1,2}, new int[]{3});
        ROT = new Quaternion(P_prime.getMatrix(new int[]{0,1,2}, new int[]{0,1,2}));
        
    }


    /**
     * Get the correct camera matrix solution for this essential matrix, using
     * the given list of point correspondences to figure out which of the four
     * solutions is the correct one.
     * @param P_SOLS
     * @param INLIERS
     * @return
     */
    public int getCorrectSolution(Matrix[] P_SOLS,
                                  List<Correspondence_2D_2D_3D> INLIERS){

        if(P_SOLS.length!=4)
            throw new RuntimeException("Wrong number of possible solutions!");

        if(INLIERS.isEmpty())
            throw new RuntimeException("Require at least one inlying"
                                       +" point correspondence to identify " +
                                      "correct solution for camera matrices!");

        // Array to count votes for correct solution
        int[] VOTES = new int[4];

        // Loop over all inlying points and accumulate votes to find solution
        for(int p=0; p<INLIERS.size(); p++){

            // Get reference to this inlying point
            Correspondence_2D_2D_3D point = INLIERS.get(p);

            // Use each possible solution for P' in conjunction with fixed P
            // to triangulate the 3D location of the observed point.
            //
            // Each solution is a homogenous 4x1 position vector in world
            // coordinates
            Matrix[] X = new Matrix[]{triangulate(P_SOLS[0], P, K_prime, K, point),
                                      triangulate(P_SOLS[1], P, K_prime, K, point),
                                      triangulate(P_SOLS[2], P, K_prime, K, point),
                                      triangulate(P_SOLS[3], P, K_prime, K, point)};

            // record depth to point in each solution pair
            double[][] DEPTH = new double[4][2];

            for (int sol = 0; sol < 4; sol++) {
                
                // Transform position vector of point from world frame to frame
                // of un-primed camera
                Matrix T       = P.times(X[sol]);
                // Transform position vector of point from world frame to frame
                // of primed camera
                Matrix T_prime = P_SOLS[sol].times(X[sol]);

                // Get Z component of position vector relative to each 
                // coordinate frame
                DEPTH[sol][0] = T.get(2, 0);
                DEPTH[sol][1] = T_prime.get(2, 0);

            }

            // Voting based on depth
            for(int s=0; s<4; s++){

                // If both depths are positive, then this is a viable
                // solution.
                if(DEPTH[s][0] > 0 && DEPTH[s][1] > 0)
                    VOTES[s]++;

            }


        }

        // Now find solution with the most votes

        // Record # votes (0) and corresponding index in VOTES (1)
        int[] MAX_INDEX = {-1,-1};

        for(int s=0; s<4; s++){

            // In case of a tie, the first solution encountered wins
            if(VOTES[s] > MAX_INDEX[0])
                MAX_INDEX = new int[]{VOTES[s],s};

        }

        return MAX_INDEX[1];

    }


    /**
     * Method to calculate an estimate of the 3D position p3d, using the given
     * image coordinates and the two camera matrices passed as argument.
     *
     * The camera frames are already normalised, so P = [R|t] etc.
     *
     * Follows method in section 12.2 (pp312) of Multiple View Geometry
     *
     */
    public static Matrix triangulate(Matrix P_prime,
                                     Matrix P,
                                     Matrix K_prime,
                                     Matrix K,
                                     Correspondence_2D_2D_3D point){

        // Triangulation works from basic equations p2d=PX, p2d'=P'p3d

        // Get normalised image coordinates in each frame
        Matrix xhat  = K.inverse().times(point.p2d.x);
        Matrix xphat = K_prime.inverse().times(point.p2d_prime.x);

        // Get rows of camera matrix P
        Matrix P1 = P.getMatrix(0, 0, 0, 3);
        Matrix P2 = P.getMatrix(1, 1, 0, 3);
        Matrix P3 = P.getMatrix(2, 2, 0, 3);

        // Get rows of camera matrix P'
        Matrix Pp1 = P_prime.getMatrix(0, 0, 0, 3);
        Matrix Pp2 = P_prime.getMatrix(1, 1, 0, 3);
        Matrix Pp3 = P_prime.getMatrix(2, 2, 0, 3);

        Matrix A = new Matrix(6,4);

        // Build matrix A. Note that I have included extra rows than
        // is shown in MVG. May as well include as much constraint as poss.
        A.setMatrix(0, 0, 0, 3, P3.times(xhat.get(1, 0)).minus(P2));
        A.setMatrix(1, 1, 0, 3, P3.times(xhat.get(0, 0)).minus(P1));
        A.setMatrix(2, 2, 0, 3, P2.times(xhat.get(0, 0)).minus(P1.times(xhat.get(1, 0))));
        A.setMatrix(3, 3, 0, 3, Pp3.times(xphat.get(1, 0)).minus(Pp2));
        A.setMatrix(4, 4, 0, 3, Pp3.times(xphat.get(0, 0)).minus(Pp1));
        A.setMatrix(5, 5, 0, 3, Pp2.times(xphat.get(0, 0)).minus(Pp1.times(xphat.get(1, 0))));

        // Singular value decomposition of A:
        SingularValueDecomposition svd = A.svd();

        // Solution is right singular vector corresponding to smallest
        // singular value; final column of V.
        Matrix X = svd.getV().getMatrix(0, 3, 3, 3);

        // Normalise:
        X.timesEquals(1.0/X.get(3, 0));

        // p3d is 4x1
        return X;

    }



    /**
     * Triangulate the 3D locations of all inlying points used to estimate
     * essential matrix.
     */
    public final void triangulatePoints(){

        // If Fundamental matrix from which this essential matrix was derived
        // was not estimated from point correspondences, stop now.
        if(fund.points.isEmpty()){
            System.err.println("No inlying points to triangulate!");
            return;
        }

        // If intrinsic or extrinsic parameters of each camera matrix
        // (K, K_prime, P & P_prime) are not fully known, stop now.
        if(K==null || K_prime==null || P==null || P_prime==null){
            System.err.println("Camera properties must be fully defined for"+
                               " triangulation!");
            return;
        }

        // Set homogenous position vector in world frame for each inlying point
        for(Correspondence_2D_2D_3D point : fund.points){

            // Triangulate 3D position in world frame
            point.setX(triangulate(P_prime,P,K_prime,K,point));

            // Set back-projected coordinates of triangulated point
            point.setSubsPoints(K, P, K_prime, P_prime);
            
            
            // At this point, it is possible to exclude points for which the
            // triangulated position is behind the camera. This removes
            // points for which triangulated position is significantly in 
            // error. However, it is likely that the triangulated position will
            // be improved using the gold standard algorithm, so I think it
            // is better to include them in the data.
        }

        // Done

    }




    public static Matrix getConjugateTransform(Matrix p){

        // New matrix to hold transformation
        Matrix p_conj = new Matrix(3,4);

        // Extract rotation and translation matrices
        Matrix r = p.getMatrix(new int[]{0,1,2}, new int[]{0,1,2});
        Matrix t = p.getMatrix(new int[]{0,1,2}, new int[]{3});

        // take inverse rotation...
        r = r.transpose();

        // conjugate translation...
        t = r.times(t).times(-1);

        // Load these back into matrix Pnorm
        p_conj.setMatrix(new int[]{0,1,2}, new int[]{0,1,2}, r);
        p_conj.setMatrix(new int[]{0,1,2}, new int[]{3}, t);

        return p_conj;
    }
    
    
    
    
    
    /**
     * Use the gold standard algorithm to improve the solution for the
     * Essential Matrix. Can only be used for Essential Matrices estimated
     * from point correspondences.
     */
    public void goldStandardAlgorithm(int N_ITER){
        
        if(fund.inliers.size()<8){
            System.err.println("Cannot apply gold standard algorithm with "
                               +fund.inliers.size()+" points");
            return;
        }
         
        // Check that necessary linear estimates of camera matrices and
        // 3D point locations have been obtained
        if(P==null||P_prime==null){
            System.err.println("Missing initial estimate of camera matrices!");
            return;
        }
        for (Correspondence_2D_2D_3D corr : fund.inliers) {
            if(corr.getX()==null){
                System.err.println("Missing initial point triangulation!");
                return;
            }
        }
        
        
        // Configure LM algorithm
        setExitTolerance(1E-92);
        setMaxDamping(1E32);
        
        // Use Levenberg-Marquardt algorithm to improve solution for
        // P' and the triangulated 3D locations of all points.
        fit(N_ITER,true);
        
        // Now use improved value of P' to reset Essential and Fundamental
        // matrix elements.
             
        // Algorithm only adjusts R', t' and triangulated positions of
        // features. Must reset fundamental and essential matrices
        // from new R' and t' values.

        // Reset fundamental matrix (& epipoles)
        fund.setFundFromCameras(K.times(P), K_prime.times(P_prime));
        
        // Reset essential matrix
        E = K_prime.transpose().times(fund.F.times(K));
        
    }
    
    
    /**
     * Get geometry and 2D/3D point correspondences for unprimed frame as
     * a Frame object suitable for use with BundleAdjustment algorithm.
     * @return 
     */
    public Frame getUnprimedFrame(){
    
        // Define Pose object from camera matrix
        Pose conjPose = new Pose(P);
        
        // Set geometry of new Frame
        Frame frame = new Frame(K,conjPose.getConjP());       
                
        // Copy 2D/3D point correspondences for inlying points
        for(Correspondence_2D_2D_3D corr2d2d3d : fund.inliers)
            frame.addCorrespondence_2D_3D(new Correspondence_2D_3D(corr2d2d3d.p2d,corr2d2d3d.p3d)); 
    
        return frame;
    }
    /**
     * Get geometry and 2D/3D point correspondences for primed frame as
     * a Frame object suitable for use with BundleAdjustment algorithm.
     * @return 
     */
    public Frame getPrimedFrame(){
    
        // Define Pose object from camera matrix
        Pose conjPose = new Pose(P_prime);
        
        // Set geometry of new Frame
        Frame frame = new Frame(K_prime,conjPose.getConjP());
        
        // Copy 2D/3D point correspondences for inlying points
        for(Correspondence_2D_2D_3D corr2d2d3d : fund.inliers)
            frame.addCorrespondence_2D_3D(new Correspondence_2D_3D(corr2d2d3d.p2d_prime,corr2d2d3d.p3d)); 
    
        return frame;
    }    
    
    
    
    
    
    
    /**
     * There now follows a bunch of methods to implement the Levenberg-Marquardt
     * class 
     * 
     */
    



    // Data is coordinates of each tracked feature in both frames
    // Parameters are 7 pose parameters plus three for every tracked point
    
    /** Four image ordinates (i,j,i',j') for each inlying point */
    public int getDataN(){ return fund.inliers.size() * 4;}
    
    /** Seven parameters of second camera pose, and 3 per inlying point */
    public int getParametersN(){ return 7 + fund.inliers.size()*3;}    
    
    /** Set reasonable step sizes here. */
    public Matrix getParametersSteps(){
        
        Matrix STEPS = new Matrix(getParametersN(),1);

        // Step size for pose rotation quaternion elements
        double STEP_Q=0.001;
        // Step size for pose translation vector elements
        double STEP_T=5;
        // Step size for triangulated point locations
        double STEP_X=5;
        
        // Enter steps into Matrix
        for(int p=0; p<getParametersN(); p++)
            STEPS.set(p, 0, (p<4) ? STEP_Q : ((p<7) ? STEP_T : STEP_X));
        
        return STEPS;
    }

    
    /** Required by LMA
     * 
     * New parameters are organised in the matrix delta in the following
     * manner:
     * 
     * 0-3: components of Quaternion used to represent rotation of P'
     * 4-6: components of translation vector of P'
     * 7-9: position of first inlying point
     * 10-12: position of second inlying point
     *   :
     * 7+3n - 9+3n: position of nth inlying point
     * 
     *
     */
    public void setParameters(Matrix delta){

        // First four elements are components of Quaternion used to represent
        // rotation. Update Quaternion, normalise it, then reset value of
        // R' in P'

        // First update components...
        ROT.re   = delta.get(0, 0);
        ROT.im.setX(delta.get(1, 0));
        ROT.im.setY(delta.get(2, 0));
        ROT.im.setZ(delta.get(3, 0));
        // Enforce unit quaternion...
        ROT.normalise();
        
        // Convert to rotation matrix and assign to P'
        P_prime.setMatrix(new int[]{0,1,2}, new int[]{0,1,2}, ROT.toMatrix());

        // Elements 4-6 are the new translation vector t' in P'.
        // Update vector t'
        TRANS = delta.getMatrix(new int[]{4,5,6}, new int[]{0}).copy();
        
        // Assign back to P'
        P_prime.setMatrix(new int[]{0,1,2}, new int[]{3}, TRANS);

        // Now loop over all inlying points and update their 3D coordinates
        for(int p=0; p<fund.inliers.size(); p++){
            fund.inliers.get(p).setXi(0,delta.get(7+3*p, 0));
            fund.inliers.get(p).setXi(1,delta.get(8+3*p, 0));
            fund.inliers.get(p).setXi(2,delta.get(9+3*p, 0));
        }
       
    }
        
    
    
    
    /** Required by LMA */
    public Matrix getParameters(){
                
        // Construct a Matrix to hold parameter values
        Matrix PARAMS = new Matrix(getParametersN(),1);
        
        // Write quaternion components to first 4 elements
        PARAMS.set(0, 0, ROT.re);
        PARAMS.set(1, 0, ROT.im.getX());
        PARAMS.set(2, 0, ROT.im.getY());
        PARAMS.set(3, 0, ROT.im.getZ());
        // Write translation components to next 3 elements
        PARAMS.setMatrix(new int[]{4,5,6}, new int[]{0}, TRANS);
        // Write 3D positions of each inlying point to remaining elements.
        // Ignore homogenous component.
        for(int p=0; p<fund.inliers.size(); p++){
            PARAMS.setMatrix(new int[]{7+3*p,8+3*p,9+3*p}, 
                             new int[]{0}, 
                             fund.inliers.get(p).getX().getMatrix(new int[]{0,1,2}, new int[]{0}));
        }

        return PARAMS;
    }      
    
    /**
     * Override method of Levenberg-Marquardt class. Naive calculation of 
     * degrees of freedom by data - parameters doesn't work for Essential
     * Matrix, because system is over-parameterised. Although system has four
     * rotation and three translation parameters, rotation has 3 degrees of
     * freedom and translation 2 (only direction is constrained).
     */
    @Override
    public double getDOF(){
        return getDataN() - getParametersN() + 2;
    }
    
    
    
    /** 
     * Required by LMA
     * 
     * Covariance matrix for full data set is block-diagonal:

     * Entry for      { s1_ii   s1_ij   0   0           - position covariance matrix
     * first          | s1_ij   s1_jj   0   0             in un-primed frame
     * correspondence | 0       0      s1_ii' s1_ij'    - position covariance matrix
     *                { 0       0      s1_ij' s1_jj'      in primed frame
     *                                                 ... ...
     *                                                 ... ...
     *                                                         ... ...
     *                                                         ... ...
     * 
     */
    public Matrix getCovariance(){
    
        Matrix covariance = new Matrix(getDataN(), getDataN());
        
        for(int p=0; p<fund.inliers.size(); p++){
            
            // Insert covariance matrix in unprimed frame
            covariance.setMatrix(new int[]{4*p, 4*p+1}, 
                                 new int[]{4*p, 4*p+1}, 
                                 fund.inliers.get(p).p2d.cov);
            
            // Insert covariance matrix in primed frame
            covariance.setMatrix(new int[]{4*p+2, 4*p+3}, 
                                 new int[]{4*p+2, 4*p+3}, 
                                 fund.inliers.get(p).p2d_prime.cov); 
            
        }
        
        return covariance;
        
    }

    
    /** 
     * Required by LMA.
     * 
     * Data is arranged in the vector as follows:
     * 
     * Element  Quantity
     *  0       i coordinate of feature 1
     *  1       j coordinate of feature 1
     *  2       i' coordinate of feature 1
     *  3       j' coordinate of feature 1
     *  4
     *  .
     *  .
     *  .
     *  N       j' coordinate of feature N
     * 
     * 
     */
    public Matrix getData(){
    
        Matrix data = new Matrix(getDataN(), 1);
        
        for(int p=0; p<fund.inliers.size(); p++){
            data.set(4*p+0, 0, fund.inliers.get(p).p2d.x.get(0, 0));
            data.set(4*p+1, 0, fund.inliers.get(p).p2d.x.get(1, 0));
            data.set(4*p+2, 0, fund.inliers.get(p).p2d_prime.x.get(0, 0));
            data.set(4*p+3, 0, fund.inliers.get(p).p2d_prime.x.get(1, 0));
        }
        
        return data;
    }  

    
    /** 
     * Required by LMA. Used to test sensitivity of parameter solution to 
     * changes in the data.
     */
    public boolean updateData(Matrix delta){
            
        for(int p=0; p<fund.inliers.size(); p++){
            
            // Build matrix representing update to ij coordinates for point p.
            // Data stored in homogenous matrices so must add third component.
            Matrix d_x = new Matrix(new double[][]{{delta.get(4*p+0, 0)},
                                                   {delta.get(4*p+1, 0)},
                                                   {0}});
            
            // Build matrix representing update to ij' coordinates for point p.
            // Data stored in homogenous matrices so must add third component.
            Matrix d_x_prime = new Matrix(new double[][]{{delta.get(4*p+2, 0)},
                                                         {delta.get(4*p+3, 0)},
                                                         {0}});            
            
            // Now add these updates to the original data
            fund.inliers.get(p).p2d.x.plusEquals(d_x);
            fund.inliers.get(p).p2d_prime.x.plusEquals(d_x_prime);
            
        }
        
        return true;
    }      
    
    
    /** Required by LMA. */    
    public Matrix getModel(){
    
        Matrix model = new Matrix(getDataN(), 1);
        
        for(int p=0; p<fund.inliers.size(); p++){
    
            fund.inliers.get(p).setSubsPoints(K,P, K_prime,P_prime);   
    
            model.set(4*p+0, 0, fund.inliers.get(p).p2d.x_sub.get(0, 0));
            model.set(4*p+1, 0, fund.inliers.get(p).p2d.x_sub.get(1, 0));
            model.set(4*p+2, 0, fund.inliers.get(p).p2d_prime.x_sub.get(0, 0));
            model.set(4*p+3, 0, fund.inliers.get(p).p2d_prime.x_sub.get(1, 0));
        }
                
        return model;
    }
  
    
    /** Required by LMA */
    public void printParameters(){
        System.out.println("printParameters() called");
    }
}