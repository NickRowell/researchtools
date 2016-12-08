package vision;

import images.Rendering;
import infra.gui.IPanel;

import java.util.List;
import java.util.ArrayList;
import java.io.*;
import java.awt.image.*;
import java.awt.Dimension;
import java.text.DecimalFormat;
import java.awt.BorderLayout;

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

/**
 * This class examines uncertainty in Essential Matrix measured by gold 
 * standard algorithm. It uses the Levenberg-Marquardt algorithm to implement
 * the gold standard algorithm, and uses propagation of covariance to estimate
 * the uncertainty on the final parameters. It then performs Monte Carlo tests
 * using known covariance to generate new realisations of the data, and 
 * measures sample covariance for results to compare these to those found by
 * propagation of covariance.
 * 
 * It also measured histograms for each parameter to check the Gaussianity.
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
 * for the pair of cameras by the Gold Standard Algorithm.
 * 
 * 
 * 
 * 
 * 
 * @author nickrowell
 */
public class EssentialMatrixTester_MonteCarlo extends javax.swing.JFrame{
        
    public static void main(String[] args) throws IOException{
        
        // Size of images
        int NPIX = 512;
        
        // Set field of view of camera (degrees)
        double fov = 30;
        
        // Scale factor for random errors on image coordinates (pixels).
        // Note NOT equal to standard deviation along principal axes.
        double sigma = 2.0;
        
        // Number of 3D points / 2D point correspondences to use
        int N = 10;
        
        // Set max number of Monte Carlo trials
        int N_TRIALS = 30000;
        // Also, set maximum runtime in seconds so that long simulations
        // still complete if max number of trials will take way too long to
        // compute.
        final long HOURS = 1;
        final long MINS  = 0;
        // Max duration of simulation 
        final long MAX_TIME_MS = HOURS*60*60*1000 + MINS*60*1000;
        // Count how many trials were actaully performed (less than N_TRIALS
        // if program runs out of time).
        int TRIALS=0;
        
        
        // Create and set up the GUI
        EssentialMatrixTester_MonteCarlo frame = 
                new EssentialMatrixTester_MonteCarlo("EssentialMatrixTester_MonteCarlo", NPIX);
        
        // Initialise camera frames
        Frame cam_1 = new Frame();
        Frame cam_2 = new Frame();
        
        // Set coordinate frames of each camera, and get true translation
        // magnitude between frames. Note that this is used in debugging and
        // not in assessing the performance of the Essential Matrix algorithms.
        double T = getFrames(cam_1, cam_2, getK(fov, NPIX), getK(fov, NPIX));
        
        // Store all parameter solutions
        List<Matrix> solutions = new ArrayList<Matrix>();        
        
        // Set RANSAC parameters in FundamentalMatrix class
        FundamentalMatrix.setOutlier(1E9);  // Want all points as inliers
        
        // Get 3D points randomly distributed inside a sphere:
        List<Point3D> p3ds = getPoint3Ds(N, 3.0);   
        
        // Make two lists of Correspondence_2D_2D_3D objects - one with true
        // locations of 3D points and noise-free projections in the image plane,
        // the other with noisy image coordinates that will be used to measure
        // EssentialMatrix.
        List<Correspondence_2D_2D_3D> corr_true  = new ArrayList<Correspondence_2D_2D_3D>();
        List<Correspondence_2D_2D_3D> corr_noisy = new ArrayList<Correspondence_2D_2D_3D>();
        
        // Make all Correspondence_2D_2D_3D objects.
        for (Point3D p3d : p3ds) {
            
            // Project 3D points to image planes - perfect noise free coordinates
            Point2D p2d_1 = cam_1.project(p3d);
            Point2D p2d_2 = cam_2.project(p3d);
            
            // Generate random position covariance matrix for each point
            p2d_1.cov = Statistics.getRandomCovarianceMatrix(sigma, 2);
            p2d_2.cov = Statistics.getRandomCovarianceMatrix(sigma, 2);

            // Ground truth point correspondences (true locations of points in
            // 3D space along with noise-free projections in image).
            corr_true.add(new Correspondence_2D_2D_3D(cam_1.project(p3d),
                                                      cam_2.project(p3d),
                                                      p3d));
            
            // Noisy point correspondences (noise to be added later)
            corr_noisy.add(new Correspondence_2D_2D_3D(p2d_1,
                                                       p2d_2));
        }        
        
        // Add noise to observed image coordinates
        for(int p=0; p<N; p++)
            setNoisyCoordinates(corr_noisy.get(p),corr_true.get(p));
        
        // Get solution for Essential Matrix using noisy data
        EssentialMatrix E = new EssentialMatrix(corr_noisy, cam_1.K, cam_2.K);  // linear 8-pt algorithm    
        E.setCameras(cam_1.getConjP(), T);  // Set frames of both cameras
        E.triangulatePoints();          // Get initial triangulation for points
        E.goldStandardAlgorithm(1000);  // Use gold standard algorithm to improve E solution
        E.setCameras(cam_1.getConjP(), T);  // Set frames of both cameras
        EssentialMatrix.setDELTAD(sigma); // set step size in data for numerical differentiation
        Matrix covariance = E.getFourthOrderCovariance();  // Get propagated parameter covariance
        
        // Bin chi-square statistic
        Histogram chi2 = new Histogram(0,E.getDOF()*3,E.getDOF()*0.1, true);
         
        String stub = "/home/nickrowell/Java_projects/MathsUoD/docs/Vision/"
                     +"EssentialMatrix tests/MonteCarlo/";       
        
        BufferedWriter out = new BufferedWriter(new FileWriter(new File(stub+"parameter_covariance")));
        out.write("Covariance by propagation (fourth order):\n"+printMatrix(covariance.getMatrix(0, 6, 0, 6)));
        out.flush();
     
        
        // Get time of simulation start
        long START_TIME = System.currentTimeMillis();
        
        for(    TRIALS = 0;
                TRIALS<=N_TRIALS && (System.currentTimeMillis() - START_TIME < MAX_TIME_MS);
                TRIALS++){
            
            // Draw perfect noise-free points onto image plane
            frame.refreshImage();             
            
            // Get new realisation of noisy image plane coordinates of points
            for(int p=0; p<N; p++)
                setNoisyCoordinates(corr_noisy.get(p),corr_true.get(p));   
            
            // Draw noise-free points onto image plane in white
            frame.drawPoint2Ds(corr_true, 0xFFFFFF);
            // Draw noisy points onto image plane in white
            frame.drawPoint2Ds(corr_noisy, 0xFF0000);
            
            // Calculate EssentialMatrix using linear method
            E = new EssentialMatrix(corr_noisy, cam_1.K, cam_2.K);

            // Set unprimed camera frame to cam_1. Set true translation so that
            // triangulated locations of points can be directly compared to 
            // true locations
            E.setCameras(cam_1.getConjP(), T);
            
            // Triangulate inlying points and project onto image plane
            E.triangulatePoints();

            // Get non-linear optimisation for essential matrix
            E.goldStandardAlgorithm(10000);
            
            // Store chi-squared
            chi2.add(E.getChi2());            
            
            // Set unprimed camera frame to cam_1. Set true magnitude of 
            // translation so that direction errors can be tested.
            E.setCameras(cam_1.getConjP(), T);
            
            // Store parameters solution
            solutions.add(E.getParameters().getMatrix(0,6,0,0));

            // Draw projected 3D points obtained by gold standard algorithm 
            frame.drawProjectedPoint3Ds(E, 0x0000FF, 2);  
            
            frame.details.setText("Trial "+TRIALS+"/"+N_TRIALS);
            
            // Repaint frame. Note that because this method runs in a different
            // thread, if the refreshImage() method is called too quickly then
            // there is a risk that the images will be partially erased by the 
            // time they are displayed in the frame.
            frame.repaint();
            
            // Pause to allow frame to be repainted before images are erased
            try{Thread.sleep(1000);}
            catch(InterruptedException ie){}

        }
        
        System.out.println("Simulation complete: "+TRIALS+" Monte Carlo shots");
                
        // Histograms to measure distribution of each parameter.
        // 4 Quaternion elements and 3 translation elements
        Histogram q0 = new Histogram(-1.0, 1.0, 0.005, false);
        Histogram q1 = new Histogram(-1.0, 1.0, 0.005, false);
        Histogram q2 = new Histogram(-1.0, 1.0, 0.005, false);
        Histogram q3 = new Histogram(-1.0, 1.0, 0.005, false);
        Histogram t0 = new Histogram(-1.0, 1.0, 0.01, false);
        Histogram t1 = new Histogram(-10.0, 10.0, 0.1, false);
        Histogram t2 = new Histogram(-10.0, 10.0, 0.1, false);
        
        // Calculate statistics:
        //
        // Sample covariance for parameters & distribution of each parameter.
        //
        Matrix mean = new Matrix(7,1);
        
        for (Matrix solution : solutions) {
            
            // Add up solutions
            mean.plusEquals(solution);
            
            // Accumulate histogram for error on each parameter
            q0.add(solution.get(0, 0) - cam_2.getAttitude().getRe());
            q1.add(solution.get(1, 0) - cam_2.getAttitude().getQ1());
            q2.add(solution.get(2, 0) - cam_2.getAttitude().getQ2());
            q3.add(solution.get(3, 0) - cam_2.getAttitude().getQ3());
            t0.add(solution.get(4, 0) - cam_2.getPosition().getX());
            t1.add(solution.get(5, 0) - cam_2.getPosition().getY());
            t2.add(solution.get(6, 0) - cam_2.getPosition().getZ());
            
        }
        
        // Sample mean of parameter solutions
        mean.timesEquals(1.0/(double)TRIALS);
        //
        // Measure sample covariance of Q & T components
        Matrix sample_covariance = new Matrix(7,7);
        for (Matrix solution : solutions) {
            Matrix diff = solution.minus(mean);
            sample_covariance.plusEquals(diff.times(diff.transpose()));           
        }

        // Normalise
        sample_covariance.timesEquals(1.0/(double)(TRIALS-1));       

        // Write out covariance on R & T error
        out.write("\nSample covariance on parameters ("+TRIALS+" shots):\n"
                  +printMatrix(sample_covariance));
        out.flush();        
        
        // Write out parameter distributions
        out = new BufferedWriter(new FileWriter(new File(stub+"q0")));
        out.write(q0.print(true));
        out.flush();
        out = new BufferedWriter(new FileWriter(new File(stub+"q1")));
        out.write(q1.print(true));
        out.flush();        
        out = new BufferedWriter(new FileWriter(new File(stub+"q2")));
        out.write(q2.print(true));
        out.flush();        
        out = new BufferedWriter(new FileWriter(new File(stub+"q3")));
        out.write(q3.print(true));
        out.flush();         
        out = new BufferedWriter(new FileWriter(new File(stub+"t0")));
        out.write(t0.print(true));
        out.flush();         
        out = new BufferedWriter(new FileWriter(new File(stub+"t1")));
        out.write(t1.print(true));
        out.flush();         
        out = new BufferedWriter(new FileWriter(new File(stub+"t2")));
        out.write(t2.print(true));
        out.flush();         

        // Write out reduced chi-2 distribution
        out = new BufferedWriter(new FileWriter(new File(stub+"chi2")));
        out.write("# Degrees of freedom = "+E.getDOF()+"\n");
        out.write(chi2.print(true));
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
    
    
    public EssentialMatrixTester_MonteCarlo(String name, int NPIX) {
        
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
        details = new JTextArea("Trial ");
        
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
    
        DecimalFormat xexxx = new DecimalFormat("0.00E0");
        
        String out = "";
        

        for(int r=0; r<in.getRowDimension(); r++){
            for(int c=0; c<in.getColumnDimension(); c++){  
                
                out += (in.get(r,c)<0 ? "" : " ")+xexxx.format(in.get(r,c))+" ";
            
            }
            out += "\n";
        }
        
        return out;
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
        // Columns are basis vectors of world frame expressed in camera frame.
        //
        Matrix r_1 = new Matrix(new double[][]{{-1, 0, 0},
                                               { 0, 1, 0},
                                               { 0, 0,-1}});
        //
        // Vector is camera frame position of world frame origin
        //
        Matrix t_1 = new Matrix(new double[][]{{0},
                                               {0},
                                               {10.}});
        
        cam_1.setK(K1);
        cam_1.setPosition(new Vector3d(t_1));
        cam_1.setAttitude(new Quaternion(r_1));
    
        
        // Position second camera on +x axis looking along -x, at distance 10:
        //
        // Columns are basis vectors of world frame expressed in camera frame.
        //
        Matrix r_2 = new Matrix(new double[][]{{ 0, 0,-1},
                                               { 1, 0, 0},
                                               { 0,-1, 0}});
        //
        // Vector is camera frame position of world frame origin
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
    
    
    public static void setNoisyCoordinates(Correspondence_2D_2D_3D corr_noisy,
                                           Correspondence_2D_2D_3D corr_true){

        // Draw noisy coordinates using covariance matrix of measured points 
        // and true coordinates
        Matrix new_coords       = Statistics.drawRandomVector(corr_noisy.p2d.cov, corr_true.p2d.x.getMatrix(0, 1, 0, 0));
        Matrix new_coords_prime = Statistics.drawRandomVector(corr_noisy.p2d_prime.cov, corr_true.p2d_prime.x.getMatrix(0, 1, 0, 0));

        // Assign noisy coordinates to measured points
        corr_noisy.p2d.x.setMatrix(new int[]{0, 1}, new int[]{0}, new_coords);
        corr_noisy.p2d_prime.x.setMatrix(new int[]{0, 1}, new int[]{0}, new_coords_prime);
            
    }
    
    
    
}

