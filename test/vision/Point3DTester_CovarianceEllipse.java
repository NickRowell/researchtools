package vision;


import images.Rendering;

import java.io.*;
import java.util.Random;
import java.util.List;
import java.util.LinkedList;
import java.awt.image.BufferedImage;

import javax.imageio.ImageIO;

import stats.Statistics;
import vision.Frame;
import vision.Point2D;
import vision.Point3D;
import dynamics.Pose;
import misc.Quaternion;
import misc.Vector3d;
import Jama.*;

import java.text.DecimalFormat;


/**
 * This class tests that the image plane uncertainty ellipses calculated for
 * projected landmarks are correct. It does this by loading a bunch of trial
 * landmarks from the TestDatabase, then using their mean positions and 
 * covariances to draw many random possible locations in 3D space for each
 * landmark. These are then projected onto the image plane to provide many
 * possible observed locations for each landmark. The covariance matrix for
 * each landmark is then projected to the image plane and drawn onto the image
 * as an ellipse. The idea is then that I can check by eye that the ellipses
 * align with the distribution of observed locations for each landmark.
 * 
 * @author nickrowell
 */
public class Point3DTester_CovarianceEllipse {
    
    static Random rand = new Random();
    
    // Camera projection matrix
    static Matrix K = new Matrix(new double[][]{{512,  0,  512},
                                                {  0, 512, 512},
                                                {  0, 0,   1}});
    // Image boundaries in projection plane
    static int X0=0,XMAX=1023,Y0=0,YMAX=1023;
    
    public static void main(String[] args) throws IOException{
    
        // Load Point3D objects to memory
        List<Point3D> p3ds = new LinkedList<Point3D>();
        File points = new File("test/Vision/Point3DTester_Points");
        BufferedReader in = new BufferedReader(new FileReader(points));
        String data;
        while((data = in.readLine())!=null){
            p3ds.add(parsePoint3D(data));
        }
           
        // Initialise camera frame
        Frame camera = new Frame(K, new Pose(new Vector3d(0,-10,0),
                                             new Quaternion(1.0/Math.sqrt(2.0),
                                                            -1.0/Math.sqrt(2.0),
                                                            0.0,0.0)));
        
        // Initialise grey pixel array for imaging random landmark positions
        int[][] image = new int[XMAX-X0+1][YMAX-Y0+1];
        // Initialise grey pixel array for imaging landmark position pdf
        int[][] pdf = new int[XMAX-X0+1][YMAX-Y0+1];        
       
        
        // Project all Point3D object to image plane
        List<Point2D> p2ds = new LinkedList<Point2D>();
        for (Point3D  p3d: p3ds)
            p2ds.add(camera.project(p3d));
       
        
        // Loop over each Point3D
        for (Point3D p3d : p3ds) {
      
            System.out.print("Point: "+p3ds.indexOf(p3d)+", ");
            
            // Now project Point3D position and covariance to image plane
            Point2D p2d = camera.project(p3d);
            
            // Record fraction of points lying within one sigma confidence
            // boundary.
            double inside_1_sig = 0, outside_1_sig = 0;
            
            
            // Draw many random positions for landmark from it's position
            // distribution.
            for (int t = 0; t < 10000; t++) {

                // Draw a random position for landmark L
                Matrix X_CAM_RAND = drawCameraFramePosition(camera, p3d);

                // Project onto image plane
                Matrix x = K.times(X_CAM_RAND);

                // Normalise homogenous image coordinates
                x.timesEquals(1.0 / x.get(2, 0));

                if(p2d.getSigmasFromMean2(Math.rint(x.get(0, 0) - X0), 
                                          Math.rint(x.get(1, 0) - Y0)) < 1.0){
                    inside_1_sig++;
                }
                else{
                    outside_1_sig++;
                }
                
                // Enter into pixel array if point lies in image
                if ((      x.get(0, 0) >= X0
                        && x.get(0, 0) <= XMAX
                        && x.get(1, 0) >= Y0
                        && x.get(1, 0) <= YMAX)) {

                    image[(int) Math.rint(x.get(0, 0) - X0)]
                         [(int) Math.rint(x.get(1, 0) - Y0)]=255;
                }

            }
            
            DecimalFormat xpxx = new DecimalFormat("0.00");
            
            System.out.println("percentage within 1-sig boundary = "+
                    xpxx.format(100*(inside_1_sig/(inside_1_sig+outside_1_sig)))+
                    "% (theoretical = 39.35%)");
            

            
            // Now loop  over all pixels in image, and draw pdf into frame.
            // To avoid later Point3D overwriting all others, set pixel 
            // intensity equal to maximum of new value and existing value.
            for(int i=0; i<pdf.length; i++)
                for(int j=0; j<pdf[i].length; j++)
                    pdf[i][j] = (int)Math.max(p2d.getPositionPDF(i, j)*300000.0,
                                              pdf[i][j]);
            
        }
        
        // Image is currently black, with single pixels painted white at the
        // points of the randomly drawn landmark positions.
        BufferedImage out_random = renderImage(image, p2ds);

        // Do the same with image of position pdfs for each landmark
        BufferedImage out_pdf = renderImage(pdf, p2ds);
        
        // Save to test directory
        ImageIO.write(out_random, "gif", new File("test/Vision/landmarks_random_positions.gif"));
        ImageIO.write(out_pdf, "gif", new File("test/Vision/landmarks_analytical_pdf.gif"));
        
    }
    
    /**
     * Generate synthetic position of landmark in the camera frame given
     * its mean camera frame position and covariance.
     * 
     * Synthetic position is a sum of mean and residual, where the residual
     * is drawn from a multivariate Gaussian represented by the Landmark's
     * position covariance matrix.
     * 
     * @return 
     */
    public static Matrix drawCameraFramePosition(Frame frame, Point3D p3d){

        // Transform Point3D position covariance matrix to camera frame
        Matrix p3d_cov_cam = frame.toBodyFrame(p3d.cov);
        
        // Transform Point3D position to camera frame
        Matrix p3d_pos_cam = frame.getConjP().times(p3d.x);
        
        // Use method of Statistics class to draw random vector
        return Statistics.drawRandomVector(p3d_cov_cam, p3d_pos_cam);

    }
    
    
    private static Point3D parsePoint3D(String data){
    
        java.util.Scanner scan = new java.util.Scanner(data);
        
        // Read all tokens from record
        int id = scan.nextInt();
        double x = scan.nextDouble();
        double y = scan.nextDouble();
        double z = scan.nextDouble();
        double sxx = scan.nextDouble();
        double syy = scan.nextDouble();
        double szz = scan.nextDouble();
        double sxy = scan.nextDouble();
        double sxz = scan.nextDouble();
        double syz = scan.nextDouble();

        // Position (column) vector
        Matrix R = new Matrix(new double[][]{{x},
                                             {y},
                                             {z},
                                             {1}});


        Matrix sigR = new Matrix(new double[][]{{sxx,sxy,sxz},
                                                {sxy,syy,syz},
                                                {sxz,syz,szz}});

        return new Point3D(R,sigR);
    
    }
    
    
    /**
     * Render an image for display by Java. The current scene is drawn into the
     * image and the feature points and tracks are overlaid in colour.
     *
     * @param image Array of gray-scale pixel values representing scene
     * @return
     */
    public static BufferedImage renderImage(int[][] image_array, 
                                            List<Point2D> p2ds){

        // Draw image
        BufferedImage image = Rendering.arrayToIm(image_array);
        
        // Now draw position and covariance ellipse for all Point2Ds in list
        for (Point2D p2d: p2ds) {

            // Quantize position,
            int[] position = new int[]{(int)Math.rint(p2d.x.get(0,0)),
                                       (int)Math.rint(p2d.x.get(1,0))};

            // Draw feature point
            Rendering.drawDiamondWithBorder(image, position, 0xFF00FFFF, 0xFF000000, 2);

                
            // Confidence ellipse - red
            Rendering.drawConfidenceEllipse(image, position, p2d.cov, 
                                            Math.sqrt(1), 0xFFFF0000);
            
        }
        
        return image;
        
    }  
    
    
}