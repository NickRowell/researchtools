package misc;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import misc.Quaternion;
import misc.Vector3d;

/**
 * Test operations of Quaternion class.
 * @author nickrowell
 */
public class QuaternionTester {
    
    public static void main(String[] args){
    
        
//        Quaternion test = new Quaternion(0.7071, 0.7071, 0, 0);
//        Matrix test_R = test.toMatrix();
//        System.out.println("Quaternion = "+test.toString());
//        System.out.println("Matrix:");
//        test_R.print(4, 4);
        
        
        Matrix test = new Matrix(new double[][]{{1,0,0},{0,0,-1},{0,1,0}});
        
        Quaternion test_q = new Quaternion(test);
        
        System.out.println("Quaternion = "+test_q.toString());
        
        
        
        
        //testCovariancePropagationMonteCarlo();
        
        
        //testFrustrumPlanes();

        //testCombineRotations();
                
        
    }
    
    /**
     * Method used during PANGU debugging to determine frustrum planes manually.
     * 
     * There are three coordinate frames in use here:
     * 
     * CAM - frame with origin at the camera projection centre and axes aligned
     *       with the image axes.
     * 
     * PANGU - frame in which positions and attitudes of dynamic objects and 
     *         camera are specified via PANGU interface
     * 
     * OPENGL - frame in which positions and attitudes are expressed for 
     *          interface with OpenGL and ModelView matrix.
     * 
     */
    public static void testFrustrumPlanes(){
    
        // Field of view in each axis, in radians.
        double xfov = Math.toRadians(23.097309);
        double yfov = Math.toRadians(23.097309);
        
        // Get accurate camera atitude by chaining two precise rotations
        // First rotation
        Quaternion r1 = new Quaternion(new Vector3d(1,0,0), Math.toRadians(-90));
        // Second rotation
        Quaternion r2 = new Quaternion(new Vector3d(0,0,1), Math.toRadians(52));
        // Combine rotations...
        Quaternion q_CAM_PANGU = r2.multiply(r1);

        // OpenGL eye coordinate frame is rotated 90 degrees wrt PANGU 
        // camera frame, so must correct for this here to get results that
        // can be directly compared to clipping planes obtained from OpenGL
        // modelview matrix derived planes.
        Quaternion q_PANGU_OPENGL = new Quaternion(1.0/Math.sqrt(2.0), -1.0/Math.sqrt(2.0), 0.0, 0.0);
        
        Quaternion q_CAM_OPENGL = q_PANGU_OPENGL.multiply(q_CAM_PANGU);

            
        // Camera origin in PANGU frame for test case
        Vector3d r_CAM_PANGU = new Vector3d( 9213000 + Math.cos(Math.toRadians(38.0)),
                                           -Math.sin(Math.toRadians(38.0)), 
                                            0.0);
        
        // Camera origin in OPENGL frame    
        Vector3d r_CAM_OPENGL = q_PANGU_OPENGL.rotate(r_CAM_PANGU);
        
        
        
        // Origin of bounding sphere for point to be clipped, in PANGU frame
        Vector3d r_PANGU = new Vector3d(9213000.0, 0.0, 0.0);
        // Transform bounding sphere origin to modelview matrix frame:
        Vector3d r_OPENGL = q_PANGU_OPENGL.rotate(r_PANGU);
        

        
        // Note that projection model is also different, but this is not
        // investigated here - we simply want clipping planes in eye coordinates.
        
        
        // Clip plane normals in camera frame - normals point outwards, i.e. 
        // a positive distance from plane to object places object outside 
        // frustrum.
        Vector3d left_clip_plane_normal_CAM   = new Vector3d(-Math.cos(xfov/2.0), 0,-Math.sin(xfov/2.0));
        Vector3d right_clip_plane_normal_CAM  = new Vector3d( Math.cos(xfov/2.0), 0,-Math.sin(xfov/2.0));
        Vector3d top_clip_plane_normal_CAM    = new Vector3d( 0,-Math.cos(yfov/2.0),-Math.sin(yfov/2.0));
        Vector3d bottom_clip_plane_normal_CAM = new Vector3d( 0, Math.cos(yfov/2.0),-Math.sin(yfov/2.0));
        
        // Transform these to world frame
        Vector3d left_clip_plane_normal_world   = q_CAM_OPENGL.rotate(left_clip_plane_normal_CAM);
        Vector3d right_clip_plane_normal_world  = q_CAM_OPENGL.rotate(right_clip_plane_normal_CAM);
        Vector3d top_clip_plane_normal_world    = q_CAM_OPENGL.rotate(top_clip_plane_normal_CAM);
        Vector3d bottom_clip_plane_normal_world = q_CAM_OPENGL.rotate(bottom_clip_plane_normal_CAM);
        
        // Shortest distance to world origin for each clipping plane: simply
        // project any point in the plane onto the surface normal. In this case
        // the known point in the plane is obtained from the camera position,
        // because each clipping plane passes through the camera origin.
        double d_left   = -r_CAM_OPENGL.dot(left_clip_plane_normal_world);
        double d_right  = -r_CAM_OPENGL.dot(right_clip_plane_normal_world);
        double d_top    = -r_CAM_OPENGL.dot(top_clip_plane_normal_world);
        double d_bottom = -r_CAM_OPENGL.dot(bottom_clip_plane_normal_world);
        
        // Print out results:
        System.out.println("Left clip plane    = "+left_clip_plane_normal_world.toString()+", d = "+d_left);
        System.out.println("Right clip plane   = "+right_clip_plane_normal_world.toString()+", d = "+d_right);
        System.out.println("Top clip plane     = "+top_clip_plane_normal_world.toString()+", d = "+d_top);
        System.out.println("Bottom clip plane  = "+bottom_clip_plane_normal_world.toString()+", d = "+d_bottom);
        
        // Distance of bounding sphere origin from each plane:
        double clip_left   = left_clip_plane_normal_world.dot(r_OPENGL) + d_left;
        double clip_right  = right_clip_plane_normal_world.dot(r_OPENGL) + d_right;
        double clip_top    = top_clip_plane_normal_world.dot(r_OPENGL) + d_top;
        double clip_bottom = bottom_clip_plane_normal_world.dot(r_OPENGL) + d_bottom;
        
        System.out.println("Bounding sphere origin in ModelView frame = "+r_OPENGL.toString());
        System.out.println("Distance to left clipping plane   = "+clip_left);
        System.out.println("Distance to right clipping plane  = "+clip_right);
        System.out.println("Distance to top clipping plane    = "+clip_top);
        System.out.println("Distance to bottom clipping plane = "+clip_bottom);
        
        
        System.out.println("N.P for left plane = "  +left_clip_plane_normal_world.dot(r_OPENGL));
        System.out.println("N.P for right plane = " +right_clip_plane_normal_world.dot(r_OPENGL));
        System.out.println("N.P for top plane = "   +top_clip_plane_normal_world.dot(r_OPENGL));
        System.out.println("N.P for bottom plane = "+bottom_clip_plane_normal_world.dot(r_OPENGL));
       
        
        
    }
    
    
    public static void testCameraBodyRotation(){
        
        Quaternion q_CAM_body = new Quaternion(-0.635543,0.635543,0.309975,-0.309975);
        
        // Vectors in body frame
        Vector3d r_BODY = new Vector3d(1,0,0);
        // Express these in camera frame
        Vector3d r_CAM = q_CAM_body.inverse().rotate(r_BODY);
        
        
        // Print result
        System.out.println("r_CAM = "+r_CAM.toString());   
    }
    
    public static void testCombineRotations(){
        
        // First rotation
        Quaternion r1 = new Quaternion(1.0/Math.sqrt(2.0),
                                       0,
                                       -1.0/Math.sqrt(2.0),
                                       0);
        // Second rotation
        Quaternion r2 = new Quaternion(0,
                                       0,
                                       0,
                                       1);
        // Combine rotations...
        Quaternion r3 = r2.multiply(r1);
        // Print combined rotation
        System.out.println("q1 = "+r1.toString());
        System.out.println("q2 = "+r2.toString());
        System.out.println("Combined rotation = "+r3.toString());  
    
    }
    
    
    /**
     * Use a Monte Carlo sampling technique to check that the rotation matrix
     * elements covariance matrix obtained from the quaternion elements
     * covariance matrix by first order propagation is correct. Or at least 
     * agrees with that obtained by numerical simulation under correlated 
     * Gaussian errors on the quaternion elements.
     */
    public static void testCovariancePropagationMonteCarlo(){
    
        // True quaternion
        Quaternion q_true = new Quaternion(0,0,0,1);        
        
        // Get Random instance for drawing Gaussian errors
        Random random = new Random();
        
        // Generate 4x4 matrix of Gaussian pseudo-random numbers
        Matrix tmp_rand = new Matrix(new double[][]{{random.nextGaussian(),random.nextGaussian(),random.nextGaussian(),random.nextGaussian()},
                                                    {random.nextGaussian(),random.nextGaussian(),random.nextGaussian(),random.nextGaussian()},
                                                    {random.nextGaussian(),random.nextGaussian(),random.nextGaussian(),random.nextGaussian()},
                                                    {random.nextGaussian(),random.nextGaussian(),random.nextGaussian(),random.nextGaussian()}});
        
        // Scale down to sigma = 0.01
        tmp_rand.timesEquals(0.01);
        
        // Make a nice positive definite matrix out of this by squaring it.
        // This provides our random covariance matrix on the quaternion parameters.
        Matrix S_q = tmp_rand.transpose().times(tmp_rand);
        
        // Propagate this to covariance matrix on corresponding rotation matrix
        // elements by first orde propagation
        Matrix S_R = q_true.toMatrixCovariance(S_q);
        
        /**
         * Now, use Monte Carlo techniques to estimate numerically the 
         * covariance matrix on the rotation matrix elements.
         */
        
        // Get Eigenvalue decomposition of random quaternion covariance matrix
        EigenvalueDecomposition eig = S_q.eig();
        
        // Get dispersions along the principal axes of the covariance hyper-ellipsoid
        Matrix evals = eig.getD();
        // ... in terms of standard deviations
        double s_q0 = Math.sqrt(evals.get(0, 0));
        double s_q1 = Math.sqrt(evals.get(1, 1));
        double s_q2 = Math.sqrt(evals.get(2, 2));
        double s_q3 = Math.sqrt(evals.get(3, 3));
        
        // Set number of Gaussian random quaternions that will be drawn
        int N = 1000000;
        
        // Store all matrices to List for statistical analysis at end
        List<Matrix> matrices_reshaped = new ArrayList<Matrix>();
        // Also aggregate sum of elements inside loop for calculating mean later
        Matrix mean = new Matrix(9,1);
        
        for(int n=0; n<N; n++){
            
            // Draw random dispersions along principal axes and form these
            // into a column vector
            Matrix q_resid_mat = new Matrix(new double[][]{{s_q0 * random.nextGaussian()},
                                                           {s_q1 * random.nextGaussian()},
                                                           {s_q2 * random.nextGaussian()},
                                                           {s_q3 * random.nextGaussian()}});
            
            // Transform this back to original frame using eigenvectors
            // of covariance matrix
            q_resid_mat = eig.getV().times(q_resid_mat);
            
            // Make a residual Quaternion out of this random error term.
            Quaternion q_resid = new Quaternion(q_resid_mat.get(0, 0),
                                                q_resid_mat.get(1, 0),
                                                q_resid_mat.get(2, 0),
                                                q_resid_mat.get(3, 0));
            
            // Add residual to true quaternion to get error quaternion
            Quaternion q_err = q_true.add(q_resid);
            
            
            // Convert this to rotation Matrix
            Matrix R_err = q_err.toMatrix();
            
            // Reshape elements into a matrix with 6 rows and one column. This
            // will greatly ease calculating the covariance of the elements.
            
            Matrix R_reshaped = new Matrix(new double[][]{{R_err.get(0,0)},
                                                          {R_err.get(0,1)},
                                                          {R_err.get(0,2)},
                                                          {R_err.get(1,0)},
                                                          {R_err.get(1,1)},
                                                          {R_err.get(1,2)},
                                                          {R_err.get(2,0)},
                                                          {R_err.get(2,1)},
                                                          {R_err.get(2,2)}});
            
            // Accumulate results
            mean.plusEquals(R_reshaped);
            matrices_reshaped.add(R_reshaped);
            
        }
        
        // Sample covariance calculation
        mean.timesEquals(1.0/(double)N);
        // now measure sample covariance of parameters
        Matrix sample_covariance = new Matrix(9,9);
        for (Matrix p_i : matrices_reshaped) {
            Matrix diff = p_i.minus(mean);
            sample_covariance.plusEquals(diff.times(diff.transpose()));           
        }
        sample_covariance.timesEquals(1.0/(double)(N-1));         
        
        
        System.out.println("Covariance of quaterion elements:");
        S_q.print(5, 5);
        System.out.println("Propagated covariance of rotation matrix elements:");
        S_R.print(5, 5);
        System.out.println("Sample covariance by Monte Carlo methods:");
        sample_covariance.print(5, 5);
    } 
    
}
