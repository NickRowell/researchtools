package vision;

import images.Rendering;
import infra.gui.IPanel;

import java.util.List;
import java.util.ArrayList;
import java.io.*;
import java.awt.image.*;
import java.awt.Graphics;
import java.awt.Dimension;
import java.text.DecimalFormat;

import javax.swing.JFrame;
import javax.swing.JTextArea;
import javax.swing.JPanel;
import javax.swing.JSeparator;

import stats.Statistics;
import vision.Correspondence_2D_2D_3D;
import vision.EssentialMatrix;
import vision.Frame;
import vision.FundamentalMatrix;
import vision.Point2D;
import vision.Point3D;
import data.Histogram;
import misc.Quaternion;
import misc.Vector3d;
import Jama.*;

import java.awt.BorderLayout;



/**
 * This class examines performance of linear & nonlinear EssentialMatrix 
 * algorithms using synthetic data.
 * 
 * It works like this. A bunch (N) of 3D points are defined in the world frame
 * positioned randomly within a sphere centred on the origin. Then, two camera
 * frames are defined - the position and orientation of each camera is set
 * separately, the intrinsic parameters (field of view, number of pixels) are
 * the same for each. The 3D points are then projected into each image plane,
 * producing N point correspondences.
 * 
 * Each 2D image point is assigned a randomly generated covariance matrix for
 * the image plane coordinates. The magnitude of the deviation along each
 * principal axis is a parameter set in the code below. Each image point then
 * has noise added to it, producing noisy point correspondences.
 * 
 * These point correspondences are then used to calculate the Essential Matrix
 * for the pair of cameras, first by the linear 8-point algorithm and then
 * by the Gold Standard Algorithm.
 * 
 * The performance of each algorithm is assessed by comparing its results to  
 * the known geometry of the camera pair. The magnitude of the translation is
 * not fixed by the Essential Matrix (leading to overall scale ambiguity) and
 * as a result the measured 3D positions of the points relative to the cameras
 * can't be compared to the ground truth. Only the relative orientation of the
 * cameras and the direction of the translation are constrained by the 
 * Essential Matrix, and these are the statistics that are recorded for each
 * realisation of the data.
 * 
 * The accuracy of the triangulated point locations is assessed by setting the
 * translation magnitude manually once the essential matrix has been measured,
 * then triangulating the points at true scale. The linear triangulation method
 * overwrites the values found by the gold standard algorithm, so this is a 
 * somewhat sketchy way to assess the gold standard algorithm and the results
 * are not compiled for analysis at the end.
 * 
 * To do: perhaps it would be interesting to bin the reduced-chi square
 *        statistic?
 * 
 * 
 * @author nickrowell
 */
public class EssentialMatrixTester_Algorithms extends javax.swing.JFrame{
        
    public static void main(String[] args) throws IOException{
        
        // Size of images
        int NPIX = 512;
        
        // Set field of view of camera (degrees)
        double fov = 30;
        
        // Magnitude of random errors on image coordinates (pixels)
        double sigma = 2.0;
        
        // Number of 3D points / 2D point correspondences to use
        int N = 50;
        
        // Set number of Monte Carlo trials
        int N_TRIALS = 120;
        
        // Create and set up the GUI
        EssentialMatrixTester_Algorithms frame = 
                new EssentialMatrixTester_Algorithms("EssentialMatrixTester_MonteCarlo", NPIX);
        
        // Initialise camera frames
        Frame cam_1 = new Frame();
        Frame cam_2 = new Frame();
        
        // Set coordinate frames of each camera, and get true translation
        // magnitude between frames. Note that this is used in debugging and
        // not in assessing the performance of the Essential Matrix algorithms.
        double T = getFrames(cam_1, cam_2, getK(fov, NPIX), getK(fov, NPIX));
        
        // Create a few Histogram objects to store histograms of rotation
        // and translation bearing error.
        Histogram rotation_errors_lin    = new Histogram(0,30,1.0,true);
        Histogram rotation_errors_nonlin = new Histogram(0,30,1.0,true);
        Histogram bearing_errors_lin       = new Histogram(0,10,0.25,true);
        Histogram bearing_errors_nonlin    = new Histogram(0,10,0.25,true);
        
        // Store all rotation and bearing errors for statistical analysis
        List<Matrix> errors_lin = new ArrayList<Matrix>();
        List<Matrix> errors_nonlin = new ArrayList<Matrix>();
        // Take running total of errors
        Matrix errors_lin_sum = new Matrix(2,1);
        Matrix errors_nonlin_sum = new Matrix(2,1);

        // Set RANSAC parameters in FundamentalMatrix class
        FundamentalMatrix.setOutlier(1E9);  // Want all points as inliers
        
        
        
        for(int TRIAL = 0; TRIAL<N_TRIALS; TRIAL++){

            // Get 3D points randomly distributed inside a sphere:
            List<Point3D> p3ds = getPoint3Ds(N, 3.0);            

            // Generate noise-free point correspondences between camera frames
            List<Correspondence_2D_2D_3D> corrs = new ArrayList<Correspondence_2D_2D_3D>();
            for (Point3D p3d : p3ds) {
                // Project 3D point to image planes - perfect noise free coordinates
                Point2D p2d_1 = cam_1.project(p3d);
                Point2D p2d_2 = cam_2.project(p3d);
                // Generate random position covariance matrix for each Point2D
                p2d_1.cov = Statistics.getRandomCovarianceMatrix(sigma, 2);
                p2d_2.cov = Statistics.getRandomCovarianceMatrix(sigma, 2);
                // Make noisy point correspondence (noise to be added later)
                corrs.add(new Correspondence_2D_2D_3D(p2d_1,p2d_2));
            }
        
            // Refresh images (paint over everything in black)
            frame.refreshImage();
            
            // Draw noise-free points onto image plane in white
            frame.drawPoint2Ds(corrs, 0xFFFFFF);
            
            // Now add noise to image points
            for(int p=0; p<N; p++)
                addNoise(corrs.get(p));
            
            // Draw noisy points onto image plane in red
            frame.drawPoint2Ds(corrs, 0xFF0000);
            
            // Calculate EssentialMatrix using linear method
            EssentialMatrix E = new EssentialMatrix(corrs, cam_1.K, cam_2.K);

            // Set unprimed camera frame to cam_1. Set true translation so that
            // triangulated locations of points can be directly compared to 
            // true locations
            E.setCameras(cam_1.getConjP(), T);
            
            // Triangulate inlying points and project onto image plane
            E.triangulatePoints();
                        
            // Draw projected 3D points obtained by linear 8-point algorithm
            frame.drawProjectedPoint3Ds(E, 0x00FF00, 2);
            
            /** Get chi squared after linear algorithm */
            double initialchi2 = E.getChi2();
            /** Calculate orientation error (angle) */
            double rotation_error_lin = cam_2.getRotationDiffMag(E.getPrimedFrame());
            /** Calculate bearing error (angle) */
            double bearing_error_lin = getBearingError(cam_1, cam_2, E.getPrimedFrame());            
            /** Calculate translation error magnitude */
            double trans_error_lin = cam_2.getTranslationDiff(E.getPrimedFrame()).norm();
            /** Calculate mean error in 3D point triangulation (inliers only) */
            Matrix rms_lin = getTriangulationError(E, p3ds, corrs);

            // Get non-linear optimisation for essential matrix
            E.goldStandardAlgorithm(1000);
            
            /** Get chi square after gold standard algorithm */
            double finalchi2 = E.getChi2();             
            
            // Set unprimed camera frame to cam_1. Set true translation so that
            // triangulated locations of points can be directly compared to 
            // true locations.
            E.setCameras(cam_1.getConjP(), T);
            // Triangulate inlying points and project onto image plane
            E.triangulatePoints();           

            // Get errors after non-linear optimisation
            double rotation_error_nonlin = cam_2.getRotationDiffMag(E.getPrimedFrame());
            /** Calculate bearing error (angle) */
            double bearing_error_nonlin = getBearingError(cam_1, cam_2, E.getPrimedFrame());  
            /** Calculate translation error magnitude */
            double trans_error_nonlin = cam_2.getTranslationDiff(E.getPrimedFrame()).norm();
            /** Calculate mean error in 3D point triangulation (inliers only) */
            Matrix rms_nonlin = getTriangulationError(E, p3ds, corrs);
            
            // Draw projected 3D points obtained by gold standard algorithm 
            frame.drawProjectedPoint3Ds(E, 0x0000FF, 2);  
            
            DecimalFormat xpxx = new DecimalFormat("0.00");
            
            frame.details.setText("Initial chi2 = "+xpxx.format(initialchi2)
                                + "\nFinal chi2  = "+xpxx.format(finalchi2)
                                + "\nInliers (outliers) = "+E.fund.inliers.size()+" ("+(E.fund.points.size()-E.fund.inliers.size())+")"
                                + "\nLinear rotation error (deg) = "+Math.rint(rotation_error_lin)
                                + "\nLinear bearing error (deg) = "+Math.rint(bearing_error_lin)
                                + "\nLinear translation error (m) = "+xpxx.format(trans_error_lin)
                                + "\nRMS linear triangulation error (m) = "+printVector(rms_lin)
                                + "\nNon-linear rotation error (deg) = "+Math.rint(rotation_error_nonlin)
                                + "\nNon-linear bearing error (deg) = "+Math.rint(bearing_error_nonlin)
                                + "\nNon-linear translation error (m) = "+xpxx.format(trans_error_nonlin)
                                + "\nRMS non-linear triangulation error (m) = "+printVector(rms_nonlin));
            
            // Repaint frame. Note that because this method runs in a different
            // thread, if the refreshImage() method is called too quickly then
            // there is a risk that the images will be partially erased by the 
            // time they are displayed in the frame.
            frame.repaint();
            
            // Pause to allow frame to be repainted before images are erased
            try{Thread.sleep(1000);}
            catch(InterruptedException ie){}
                        
            // Accumulate stats
            rotation_errors_nonlin.add(rotation_error_nonlin);
            bearing_errors_nonlin.add(bearing_error_nonlin);
            errors_nonlin.add(new Matrix(new double[][]{{rotation_error_nonlin},{trans_error_nonlin}}));
            errors_nonlin_sum.plusEquals(new Matrix(new double[][]{{rotation_error_nonlin},{trans_error_nonlin}}));
            rotation_errors_lin.add(rotation_error_lin);
            bearing_errors_lin.add(bearing_error_lin);
            errors_lin.add(new Matrix(new double[][]{{rotation_error_lin},{trans_error_lin}}));
            errors_lin_sum.plusEquals(new Matrix(new double[][]{{rotation_error_lin},{trans_error_lin}}));

        }
        
        // Sample covariance calculation:
        //
        // Take mean of linear & non-linear errors
        errors_lin_sum.timesEquals(1.0/(double)N_TRIALS);
        errors_nonlin_sum.timesEquals(1.0/(double)N_TRIALS);
        //
        // Measure sample covariance of R & T errors
        Matrix sample_covariance_lin = new Matrix(2,2);
        Matrix sample_covariance_nonlin = new Matrix(2,2);
        for (Matrix rt_lin : errors_lin) {
            Matrix diff = rt_lin.minus(errors_lin_sum);
            sample_covariance_lin.plusEquals(diff.times(diff.transpose()));           
        }
        for (Matrix rt_nonlin : errors_nonlin) {
            Matrix diff = rt_nonlin.minus(errors_nonlin_sum);
            sample_covariance_nonlin.plusEquals(diff.times(diff.transpose()));           
        }
        // Normalise
        sample_covariance_lin.timesEquals(1.0/(double)(N_TRIALS-1));
        sample_covariance_nonlin.timesEquals(1.0/(double)(N_TRIALS-1));        
        
        String stub = "/home/nickrowell/Java_projects/MathsUoD/docs/Vision/EssentialMatrix tests/";
                
        // Write out linear rotation error
        BufferedWriter out = new BufferedWriter(new FileWriter(new File(stub+"rot_lin")));
        out.write(rotation_errors_lin.print(true));
        out.flush();
        // Write out linear bearing error
        out = new BufferedWriter(new FileWriter(new File(stub+"bear_lin")));
        out.write(bearing_errors_lin.print(true));
        out.flush();        
        // Write out nonlinear rotation error
        out = new BufferedWriter(new FileWriter(new File(stub+"rot_nonlin")));
        out.write(rotation_errors_nonlin.print(true));
        out.flush();
        // Write out nonlinear bearing error
        out = new BufferedWriter(new FileWriter(new File(stub+"bear_nonlin")));
        out.write(bearing_errors_nonlin.print(true));
        out.flush();
        
        // Write out covariance on R & T error
        out = new BufferedWriter(new FileWriter(new File(stub+"rt_error_covariance")));
        out.write("Covariance on linear solution errors:\n"+printMatrix(sample_covariance_lin));
        out.write("\nCovariance on nonlinear solution errors:\n"+printMatrix(sample_covariance_nonlin));
        out.flush();
        
        System.exit(0);
        
        
    }
    
    // Instance fields of this class are those elements of GUI that are changed
    // with every Monte Carlo shot.
    //
    // JPanel fro holding IPanels
    JPanel images;
    // IPanels for displaying images
    IPanel im_1, im_2;
    // Labels for displaying fit info    
    JTextArea details;
    
    
    public EssentialMatrixTester_Algorithms(String name, int NPIX) {
        
        super(name);
        setResizable(true);
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);        
        
        // Initialise IPanels for image display
        im_1 = new IPanel(new BufferedImage(NPIX,NPIX, BufferedImage.TYPE_INT_RGB));
        im_2 = new IPanel(new BufferedImage(NPIX,NPIX, BufferedImage.TYPE_INT_RGB));

        // Get JPanel for displaying images in
        images = new JPanel();         
        // Add these to JPanel in locations determined by layout manager
        images.add(im_1);
        images.add(im_2);
        // Set preferred size of JPanel based on image size
        images.setPreferredSize(new Dimension(NPIX*2 + 20, 
                                              NPIX + 20));
        
        // Make a panel for displaying fit statistics then lay it out
        details = new JTextArea("Initial chi2 = "
                + "\nFinal chi2  = "
                + "\nInliers (outliers) = "
                + "\nLinear rotation error (deg) = "
                + "\nLinear bearing error (deg) = "
                + "\nLinear translation error (m) = "
                + "\nRMS linear triangulation error (m) = "
                + "\nNon-linear rotation error (deg) = "
                + "\nNon-linear bearing error (deg) = "
                + "\nNon-linear translation error (m) = "
                + "\nRMS non-linear triangulation error (m) = ");
        
        // Add JPanel with images to JFrame
        getContentPane().add(images, BorderLayout.NORTH);        
        // Add divider line between images and statistics printed below
        getContentPane().add(new JSeparator(), BorderLayout.CENTER);           
        // Add JTextArea to frame
        getContentPane().add(details, BorderLayout.SOUTH);        
        
        //Display the window.
        pack();
        setVisible(true); 
    }
    
    
    public static String printMatrix(Matrix in){
    
        String out = "";
        

        for(int r=0; r<in.getRowDimension(); r++){
            for(int c=0; c<in.getColumnDimension(); c++){  
                
                out += in.get(r,c)+"\t";
            
            }
            out += "\n";
        }
        
        return out;
    }
    
    // Convert a 3x1 Jama Matrix into a string representation
    public static String printVector(Matrix in){
        
        DecimalFormat xpxx = new DecimalFormat("0.00");

        return "("+xpxx.format(in.get(0,0))+", "+xpxx.format(in.get(1,0))+", "+xpxx.format(in.get(2,0))+")";
        
    }    
    
    
    /** Paint all pixels in each image black */
    public void refreshImage(){
    
        // First, wipe images (draw all pixels black)
        for(int i=0; i<im_1.image.getWidth(); i++)
            for(int j=0; j<im_1.image.getHeight(); j++){
                int pixel = 0x000000;
                im_1.image.setRGB(i, j, pixel);
            }
        for(int i=0; i<im_2.image.getWidth(); i++)
            for(int j=0; j<im_2.image.getHeight(); j++){
                int pixel = 0x000000;
                im_2.image.setRGB(i, j, pixel);
            }        
    }
    
    public void drawPoint2Ds(List<Correspondence_2D_2D_3D> corrs, int colour){
    
         
        // Loop over all correspondences.
        for(int p=0; p<corrs.size(); p++){
                   
            // Coordinates of point in unprimed image
            int[] coord_unprimed = new int[]{(int)Math.rint(corrs.get(p).p2d.x.get(0, 0)),
                                             (int)Math.rint(corrs.get(p).p2d.x.get(1, 0))};                
            // Coordinates of point in unprimed image
            int[] coord_primed = new int[]{(int)Math.rint(corrs.get(p).p2d_prime.x.get(0, 0)),
                                           (int)Math.rint(corrs.get(p).p2d_prime.x.get(1, 0))};                
            
            // Draw points in unprimed image
            Rendering.drawBorder(im_1.image, coord_unprimed, colour, 3);
            Rendering.drawBody(im_1.image, coord_unprimed, colour);
            // ...in primed image
            Rendering.drawBorder(im_2.image, coord_primed, colour, 3);
            Rendering.drawBody(im_2.image, coord_primed, colour);
            
        }
    
    }
    
        
  
    public void drawProjectedPoint3Ds(EssentialMatrix E,
                                      int colour, 
                                      int size){
            
       // Draw all triangulated 3D points in projection.
       for(Correspondence_2D_2D_3D corr : E.fund.points){
        
            // Coordinates of noisy image point in unprimed image
            int[] unprimed_coord_nois = new int[]{(int)Math.rint(corr.p2d.x.get(0, 0)),
                                                  (int)Math.rint(corr.p2d.x.get(1, 0))};
            // Coordinates of noisy image point in primed image
            int[] primed_coord_nois = new int[]{(int)Math.rint(corr.p2d_prime.x.get(0, 0)),
                                                (int)Math.rint(corr.p2d_prime.x.get(1, 0))};            

            // Coordinates of measured 3D point projected into unprimed image
            int[] unprimed_coord_proj = new int[]{(int)Math.rint(corr.p2d.x_sub.get(0, 0)),
                                                  (int)Math.rint(corr.p2d.x_sub.get(1, 0))};              
            // Coordinates of measured 3D point projected into primed image
            int[] primed_coord_proj = new int[]{(int)Math.rint(corr.p2d_prime.x_sub.get(0, 0)),
                                                (int)Math.rint(corr.p2d_prime.x_sub.get(1, 0))};              
            

            // Draw line connecting these points in white, separately in each image
            Rendering.drawLine(im_1.image, 
                               unprimed_coord_nois,
                               unprimed_coord_proj,
                               colour);
            Rendering.drawLine(im_2.image, 
                               primed_coord_nois,
                               primed_coord_proj,
                               colour);
            
            // Draw projected 3D points in designated colour - unprimed image
            Rendering.drawBorder(im_1.image, unprimed_coord_proj, colour, size);
            Rendering.drawBody(im_1.image, unprimed_coord_proj, colour);
            // ... primed image
            Rendering.drawBorder(im_2.image, primed_coord_proj, colour, size);
            Rendering.drawBody(im_2.image, primed_coord_proj, colour);

            
        }   
           
       // Now highlight outliers
       for(Correspondence_2D_2D_3D corr : E.fund.outliers){
        
            // Coordinates of measured 3D point projected into unprimed image
            int[] unprimed_coord_proj = new int[]{(int)Math.rint(corr.p2d.x_sub.get(0, 0)),
                                                  (int)Math.rint(corr.p2d.x_sub.get(1, 0))};              
            // Coordinates of measured 3D point projected into primed image
            int[] primed_coord_proj = new int[]{(int)Math.rint(corr.p2d_prime.x_sub.get(0, 0)),
                                                (int)Math.rint(corr.p2d_prime.x_sub.get(1, 0))};              

            // Draw circles into images
            Rendering.drawCircle(unprimed_coord_proj, im_1.image, 5, 0xFFFF00);
            Rendering.drawCircle(primed_coord_proj, im_2.image, 5, 0xFFFF00);   
            
        }
         
    }    
       
    public static double getFrames(Frame cam_1, Frame cam_2, Matrix K1, Matrix K2){
    
                
        // Position first camera on +z axis looking along -z, at distance 10:
        //
        // Columns are basis vectors of camera frame expressed in world frame.
        //
        Matrix r_1 = new Matrix(new double[][]{{-1, 0, 0},
                                               { 0, 1, 0},
                                               { 0, 0,-1}});
        //
        // Vector is world frame position of camera frame origin
        //
        Matrix t_1 = new Matrix(new double[][]{{0},
                                               {0},
                                               {10.}});
        
        cam_1.setK(K1);
        cam_1.setPosition(new Vector3d(t_1));
        cam_1.setAttitude(new Quaternion(r_1));
    
        
        // Position second camera on +x axis looking along -x, at distance 10:
        //
        // Columns are basis vectors of camera frame expressed in world frame.
        //
        Matrix r_2 = new Matrix(new double[][]{{ 0, 0,-1},
                                               { 1, 0, 0},
                                               { 0,-1, 0}});        
        
        //
        // Vector is world frame position of camera frame origin
        //      
        Matrix t_2 = new Matrix(new double[][]{{10},
                                               {0},
                                               {0}});
        
        cam_2.setK(K2);
        cam_2.setPosition(new Vector3d(t_2));
        cam_2.setAttitude(new Quaternion(r_2));
        
        return cam_1.getTranslationDiff(cam_2).norm();
    }
    
    /**
     * Generate N Point3D objects distributed randomly within a sphere of
     * radius r, centred on the origin.
     */
    public static List<Point3D> getPoint3Ds(int N, double r){
    
        // Count total number of Point3Ds created
        int N_P=0;
        
        // List to store objects
        List<Point3D> p3ds = new ArrayList<Point3D>();
        
        while(N_P<N){
        
            // Generate 3D point randomly located within a box of width 2*r
            // centred on the origin
            double x = 2.0*(Math.random()-0.5) * r;
            double y = 2.0*(Math.random()-0.5) * r;
            double z = 2.0*(Math.random()-0.5) * r;
          
            // Check distance from origin
            if(Math.sqrt(x*x + y*y + z*z) <= r){
            
                // Point lies within sphere of radius r
                
                // Homogenous position vector
                Matrix X = new Matrix(new double[][]{{x},{y},{z},{1.0}});
                
                // Set covariance matrix to identity
                Matrix cov = Matrix.identity(3, 3);
                
                // Add point to list
                p3ds.add(new Point3D(X,cov));
            
                // Increment counter
                N_P++;
                
            }
        
        }
    
        return p3ds;
    
    }
    
    
    /**
     * Create a camera projection matrix from a field of view and number of
     * pixels along an edge. Assumes central projection, square image, etc.
     */
    public static Matrix getK(double fov, double pix){

        // Focal length
        double f = (double) (pix / 2.0) * (1.0 / Math.tan(Math.toRadians(fov) / 2.0));
        
        // Coordinates of projection centre
        double px = pix/2.0;
        double py = pix/2.0;
        


        Matrix K = new Matrix(new double[][]{{f,0,px},
                                             {0,f,py},
                                             {0,0,1.}});
  
        return K;
    
    }
   
    public static Matrix getTriangulationError(EssentialMatrix E, 
                                               List<Point3D> p3ds,
                                               List<Correspondence_2D_2D_3D> corrs){
    
        /** Calculate mean error in 3D point triangulation (inliers only) */
        Matrix ms = new Matrix(3,1);
        double N_inliers = 0;
        for(int p=0; p<corrs.size(); p++){
                
            // Find inlying points correspondences
            if(E.fund.inliers.contains(corrs.get(p))){
                    
                N_inliers++;
                    
                // Triangulation error for this point
                Matrix tri_error = corrs.get(p).p3d.x.minus(p3ds.get(p).x);
                // Sum of squares. No element-wise operations in Jama, so this is messy...
                ms.plusEquals(new Matrix(new double[][]{{tri_error.get(0,0)*tri_error.get(0,0)},
                                                        {tri_error.get(1,0)*tri_error.get(1,0)},
                                                        {tri_error.get(2,0)*tri_error.get(2,0)}}));
            }
        }
        ms.timesEquals(1.0/N_inliers);
        // Sqrt of mean-squared error. No element-wise operations in Jama, so this is messy...
        Matrix rms = new Matrix(new double[][]{{Math.sqrt(ms.get(0,0))},
                                               {Math.sqrt(ms.get(1,0))},
                                               {Math.sqrt(ms.get(2,0))}});   
    
        return rms;
            
    }
    
    
    
    
    public static double getBearingError(Frame cam_1, Frame cam_2_true, Frame cam_2_meas){
            
        // Ground truth translation direction:
        //
        // First, get (fixed) position of cam_1 origin in world frame
        Matrix t_1w = cam_1.getP().getMatrix(new int[]{0,1,2}, new int[]{3});    //
        // Next, get ground truth position of cam_2 origin in world frame
        Matrix t_2w_gt = cam_2_true.getP().getMatrix(new int[]{0,1,2}, new int[]{3});    //
        // Difference gives ground truth translation direction in world frame
        Matrix t_gt = t_2w_gt.minus(t_1w);
        
        // Now get calculated position of cam_2 in world frame
        Matrix t_2w_meas = cam_2_meas.getP().getMatrix(new int[]{0,1,2}, new int[]{3});    //
        // Difference gives estimated translation direction in world frame
        Matrix t_meas = t_2w_meas.minus(t_1w);
        
        // Normalise both vectors
        t_gt.timesEquals(1.0/t_gt.normF());
        t_meas.timesEquals(1.0/t_meas.normF());
        // Get dot product
        double dot_prod = t_gt.transpose().times(t_meas).get(0, 0);
        // Convert to angle
        double bearing_error = Math.toDegrees(Math.acos(dot_prod));
        
        return bearing_error;
    }
    
    // Add noise to the image plane coordinates of all Point2Ds in the
    // correspondence. This assumes that when called, the coordinates for the
    // correspondence are noise-free, i.e. represent the true mean when drawing
    // Gaussian errors.
    public static void addNoise(Correspondence_2D_2D_3D corr){

        // Draw new coordinates
        Matrix new_coords       = Statistics.drawRandomVector(corr.p2d.cov, corr.p2d.x.getMatrix(0, 1, 0, 0));
        Matrix new_coords_prime = Statistics.drawRandomVector(corr.p2d_prime.cov, corr.p2d_prime.x.getMatrix(0, 1, 0, 0));

        // Assign noisy coordinates to measured points
        corr.p2d.x.setMatrix(new int[]{0, 1}, new int[]{0}, new_coords);
        corr.p2d_prime.x.setMatrix(new int[]{0, 1}, new int[]{0}, new_coords_prime);
            
    }
    
    
    
}
