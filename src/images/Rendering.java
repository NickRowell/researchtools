/**
 * Copyright (c) 2011 University of Dundee
 *
 * Name:
 * Rendering.java
 *
 * Purpose:
 * Various image processing functions and operations used in this project.
 * 
 * Language:
 * Java
 *
 * Author:
 * Nicholas Rowell (nickrowell@computing.dundee.ac.uk)
 *
 */

package images;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import numeric.geom.dim3.Vector3d;
import numeric.vision.Correspondence_2D_2D_3D;
import numeric.vision.FundamentalMatrix;
import numeric.vision.Point2D;

import java.awt.image.BufferedImage;
import java.awt.image.Raster;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;


public abstract class Rendering {

    /** Random number generator used to add Gaussian noise to images */
    static Random noise = new Random();

    /**
     * Write a BufferedImage into an integer array
     * @param im
     * @return
     */
    public static int[][] imToArray(BufferedImage im) {

        Raster in = im.getRaster();

        int WIDTH = im.getWidth();
        int HEIGHT = im.getHeight();

        int[][] image = new int[WIDTH][HEIGHT];

        for (int x = 0; x < WIDTH; x++)
            for (int y = 0; y < HEIGHT; y++)
                image[x][y] = in.getSample(x, y, 0);

        return image;
    }

    /**
     * Write an array to a BufferedImage object. Use colour image but same strength for each colour
     * @param im
     * @return
     */
    public static BufferedImage arrayToIm(int[][] im) {

        BufferedImage image = new BufferedImage(im.length, im[0].length, BufferedImage.TYPE_INT_ARGB);

        // Draw integer array of pixel values into image as graylevels
        for (int x = 0; x < im.length; x++) {
            for (int y = 0; y < im[0].length; y++) {
                int pixel = 0xFF000000 + (im[x][y] << 16) + (im[x][y] << 8) + im[x][y];
                image.setRGB(x, y, pixel);
            }
        }

        return image;
    }



//    public static void renderCentroid(double[][][] rgb, double[] centroid){
//
//        //+++ Set colour to render centroid - red +++//
//        int col = 0xFF0000;
//        // Decode hexadecimal colour to 0->1 score for RGB separately
//        double r = (col >> 16)/255.0;
//        double g = ((col >> 8) & 0x00FF)/255.0;
//        double b = (col & 0x0000FF)/255.0;
//
//        int cross = 5;
//
//        for(int offset=-1; offset <= 1; offset++){
//            for(int j = -cross; j<= cross; j++){
//                try{
//                    // Flip coordinates for Matlab
//                    rgb[(int) Math.rint(centroid[1]) + j][(int) Math.rint(centroid[0]) + j + offset][0] = r;
//                    rgb[(int) Math.rint(centroid[1]) + j][(int) Math.rint(centroid[0]) + j + offset][1] = g;
//                    rgb[(int) Math.rint(centroid[1]) + j][(int) Math.rint(centroid[0]) + j + offset][2] = b;
//                    rgb[(int) Math.rint(centroid[1]) + j][(int) Math.rint(centroid[0]) - j + offset][0] = r;
//                    rgb[(int) Math.rint(centroid[1]) + j][(int) Math.rint(centroid[0]) - j + offset][1] = g;
//                    rgb[(int) Math.rint(centroid[1]) + j][(int) Math.rint(centroid[0]) - j + offset][2] = b;
//
//                }
//                catch(Exception e){
//                    //+++ avoid exceptions if centroid is close to edge and cross overshoots image boundary +++//
//                }
//            }
//        }
//
//        return;
//
//    }


    /**
     * Scale a Harris corner strength map into image range.
     * @param harris
     * @return
     */
    public static int[][] renderHarrisMap(long[][] harris){

        //+++ Get max value from harris map +++//
        long max = harris[0][0];
        long min = harris[0][0];

        for(int i=0; i<harris.length; i++){
            for(int j=0; j<harris[i].length; j++){
                if(harris[i][j] > max) max = harris[i][j];
                if(harris[i][j] < min) min = harris[i][j];
            }
        }

        int[][] image = new int[harris.length][harris[0].length];

        for(int i=0; i<harris.length; i++){
            for(int j=0; j<harris[i].length; j++){
                image[i][j] = (int)(((harris[i][j]-min)*255)/(max-min));
            }
        }

        return image;
    }




    
    /**
     * Turn a correlation surface into an image
     *
     * @return
     */
    public static int[][] renderCorrelationSurface(double[][] corr){

        int[][] image = new int[corr.length][corr[0].length];

        for(int i=0; i<corr.length; i++){
            for(int j=0; j<corr[i].length; j++){
                image[i][j] = (int)(((corr[i][j]+1)*255)/(2));
            }
        }

        return image;
    }
    
    /**
     * Draw circle around Point2d
     * @param point
     * @param image
     * @param r
     * @param col
     */
    public static void drawCircle(Point2D point, BufferedImage image, double r, int col) {
        
        // Get pixel coordinates of Point2D. Assume already homogenised.
        int[] pix = new int[]{(int)Math.rint(point.x.get(0,0)),
                              (int)Math.rint(point.x.get(1,0))};
        
        drawCircle(pix, image,r,col);
        
    }
    
    /**
     * Draw circle around Point2d
     * @param pix
     * @param image
     * @param r
     * @param col
     */
    public static void drawCircle(int[] pix, BufferedImage image, double r, int col) {
        
        // Angular step between points on circumference that lie one pixel apart
        double ang = 2*Math.asin(0.5/r);
        
        // Loop round circumference of circle:
        for(double theta = 0; theta < 2*Math.PI; theta += ang){

            int i = pix[0] + (int)Math.rint(r*Math.sin(theta));
            int j = pix[1] + (int)Math.rint(r*Math.cos(theta));

            try{ image.setRGB(i, j, col);}
            catch(ArrayIndexOutOfBoundsException aioobe){}

        }
        
    }    
    
    /**
     * Draw complete diamond and border onto image.
     */
    public static void drawDiamondWithBorder(BufferedImage image, 
                                             int[] position, 
                                             int borderColour,
                                             int centreColour,
                                             int borderThickness){
        drawBorder(image, position, borderColour, borderThickness);
        drawBody(image, position, centreColour);
    }
    
    /**
     * Borrowed from java_FEIC - draw diamond onto image
     */
    public static void drawBorder(BufferedImage rgb, int[] coords, int colour, int t){

        // Watch out for overshooting image boundary
        for (int i = -(1+t); i < (2+t); i++) {
            for (int j = Math.abs(i) - (1+t); j <= Math.abs(Math.abs(i) - (1+t)); j++) {
                try {
                    rgb.setRGB(coords[0] + i, coords[1] + j, colour);
                }
                catch (ArrayIndexOutOfBoundsException aioobe) {}
            }
        }
    }

    /**
     * Borrowed from java_FEIC - draw diamond onto image
     */
    public static void drawBody(BufferedImage rgb, int[] coords, int colour){


        //+++ Watch out for overshooting image boundary +++//
        for (int i = -1; i < 2; i++) {
            for (int j = Math.abs(i) - 1; j <= Math.abs(Math.abs(i) - 1); j++) {
                try {
                    rgb.setRGB(coords[0] + i, coords[1] + j, colour);
                }
                catch (ArrayIndexOutOfBoundsException aioobe) {}
            }
        }
    }    

    /**
     * Get coordinates of points on circumference of confidence ellipse, for
     * drawing into image plane.
     * 
     * @param position      Mean position in image [pixels]
     * @param covariance    2x2 covariance matrix [pixels]
     * @param sigmas        Number of sigmas at which to draw confidence boundary
     * @return 
     */
    public static float[] drawEllipse(float[] position, Matrix covariance, float sigmas) {
    	
        // Components of covariance matrix in image plane
        float a = (float)covariance.get(0,0);
        float b = (float)covariance.get(0,1);
        float c = (float)covariance.get(1,1);
        
        // Eigenvalues of image covariance matrix:
        float trA  = a+c;
        float detA = a*c - b*b;
        // intermediate quantity
        float temp = (float)Math.sqrt(trA*trA - 4*detA)/2;
        // actual eigenvalues
        float e1 = trA/2.0f + temp;
        float e2 = trA/2.0f - temp;
        
        // Compute normalised eigenvectors from these
        float norm1 = 1.0f/(float)Math.sqrt(b*b/((a-e1)*(a-e1)) + 1);
        float norm2 = 1.0f/(float)Math.sqrt(b*b/((a-e2)*(a-e2)) + 1);
        
        // Components of eigenvector matrix in image plane
        //      [ A   B ]
        //  V = |       |
        //      [ C   D ]
        //
        float A = norm1 * (b/(a-e1));
        float B = norm2 * (b/(a-e2));
        float C = norm1;
        float D = norm2;
        
        // Check determinant: if it is -1, then must reverse one of the 
        // eigenvectors in order to avoid doing a reflection as well as
        // rotation when we apply this to our circular points later.
        double det = A*D - B*C;
        
        if(det<0)
        {
            // Reverse first eigenvector.
            A *= -1;
            C *= -1;
        }
        
        // Length of major axes of ellipse - number of standard deviations
        float s_p = (float)Math.sqrt(e1)*sigmas;
        float s_q = (float)Math.sqrt(e2)*sigmas;
        
        // Number of points to draw around circumference of ellipse
        int N = 36;
        // Corresponding angular step size
        float ang_step = 2.0f*(float)Math.PI/(float)N;
        
        // Set angular step size so that distance travelled round circumference
        // is max one pixel (occurs at the major axis).
        //double ang_step = Math.asin(1.0/Math.max(s_p,s_q));
        
        float c_ang,s_ang,r;
        
        float[] points = new float[N*2];
        
        int index = 0;
        
        for(int n=0; n<N; n++){
            
            // Translate index n to angle
            float ang = n * ang_step;
            
            c_ang = (float)Math.cos(ang);
            s_ang = (float)Math.sin(ang);
            
            // Get radius of ellipse at this angle
            r = s_p*s_q/(float)Math.sqrt(s_p*s_p*s_ang*s_ang + s_q*s_q*c_ang*c_ang);
            
            // Rotate point at (r*c_ang, r*s_ang) back to image frame
            float X0 = A*r*c_ang + C*r*s_ang;
            float X1 = B*r*c_ang + D*r*s_ang;
            
            // Get pixel coordinates by adding mean position
            points[index++] = position[0] + X0;
            points[index++] = position[1] + X1;
            
        }
        
        return points;
        
    }
    
    /**
     * Method to draw confidence ellipse from 2D covariance matrix
     */
    public static void drawConfidenceEllipse(BufferedImage rgb, int[] position, Matrix covariance, double sigmas, int colour) {
        
        // Get principal axes frame (p,q) of covariance ellipse
        EigenvalueDecomposition evd = new EigenvalueDecomposition(covariance);
        
        // Length of major axes - standard deviation along these directions.
        double s_p = Math.sqrt(evd.getD().get(0,0))*sigmas;
        double s_q = Math.sqrt(evd.getD().get(1,1))*sigmas;       
        
        // Eigenvector matrix for principal axes frame
        Matrix V = evd.getV();
                
        // Ensure eigenvector matrix can be used as rotation matrix: check if
        // determinant is minus one (in which case matrix is a roto-reflection)
        // and change sign of one column if so.
        if(V.det()<0){
            
            Matrix correct = new Matrix(new double[][]{{-1.0, 1.0},
                                                       {-1.0, 1.0}});
        
            V.arrayTimesEquals(correct);       
        }
        
        // Set angular step size so that distance travelled round circumference
        // is max one pixel (occurs at the major axis).
        double ang_step = Math.asin(1.0/Math.max(s_p,s_q));
        double c_ang,s_ang,r;
        int i,j;
                
        for(double ang = 0; ang < 2.0*Math.PI; ang += ang_step){
        
            c_ang = Math.cos(ang);
            s_ang = Math.sin(ang);
            
            // Get radius of ellipse at this angle
            r = s_p*s_q/Math.sqrt(s_p*s_p*s_ang*s_ang + s_q*s_q*c_ang*c_ang);

            // Convert to cartesian coordinates and store as column vector.
            Matrix X = new Matrix(new double[][]{{r*c_ang},{r*s_ang}});
        
            // Rotate back to image frame
            X = V.times(X);
                        
            // Get pixel coordinates by adding mean position.
            i = position[0] + (int)Math.rint(X.get(0,0));
            j = position[1] + (int)Math.rint(X.get(1,0));
            
            // Now try drawing this pixel in image.
            try { rgb.setRGB(i, j, colour);}
            catch (ArrayIndexOutOfBoundsException aioobe) {}            
            
        }        
    
    
    }
    
    /**
     * Use a FundamentalMatrix to draw epipolar lines onto BufferedImage
     * image, using the set of point correspondences given.
     *
     * The present frame (passed as argument) and the previous frame have been
     * used to estimate the fundamental matrix. So, the present frame is the
     * primed frame from the pair, and the correct epipolar lines to draw
     * are those projected from the line of sight ray in the un-primed frame.
     * These lines have the equation F * p2d.
     *
     * @param image
     * @param shortpoints
     */
    public static void drawEpipolarLines(BufferedImage image,
                                         FundamentalMatrix FUND,
                                         List<Correspondence_2D_2D_3D> corrs,
                                         int colour){

        int h = image.getHeight();
        int w = image.getWidth();

        // Get lines representing image boundaries
        Vector3d l1 = new Vector3d(0,1,0);       // upper edge
        Vector3d l2 = new Vector3d(1,0,-w);       // right edge
        Vector3d l3 = new Vector3d(0,1,-h);       // bottom edge
        Vector3d l4 = new Vector3d(1,0,0);       // left edge

        // Loop over all points
        for(int p=0; p<corrs.size(); p++){

            // Get epipolar line for this point in primed frame.
            Matrix epl_prime = FUND.F.times(corrs.get(p).p2d.x);

            // Get homogeneous 3-vector form for epipolar line
            Vector3d l = new Vector3d(epl_prime.get(0,0),
                                      epl_prime.get(1,0),
                                      epl_prime.get(2,0));

            // Get coordinates where this line cuts each image boundary
            Vector3d c1 = l.cross(l1);
            Vector3d c2 = l.cross(l2);
            Vector3d c3 = l.cross(l3);
            Vector3d c4 = l.cross(l4);

            // Now find crossing points along line segments that lie
            // within image region
            List<int[]> cross = new ArrayList<int[]>();

            if(c1.getX()/c1.getZ() > 0 && c1.getX()/c1.getZ() < 513) cross.add(new int[]{(int)(c1.getX()/c1.getZ()), 0});
            if(c2.getY()/c2.getZ() > 0 && c2.getY()/c2.getZ() < 513) cross.add(new int[]{512,(int)(c2.getY()/c2.getZ())});
            if(c3.getX()/c3.getZ() > 0 && c3.getX()/c3.getZ() < 513) cross.add(new int[]{(int)(c3.getX()/c3.getZ()), 512});
            if(c4.getY()/c4.getZ() > 0 && c4.getY()/c4.getZ() < 513) cross.add(new int[]{0,(int)(c4.getY()/c4.getZ())});

            //. Skip this line if it has other than two crossings of the
            // image boundary (indicates somethin funny going on)
            if(cross.size()!=2)
                continue;

            // Draw the epipolar line
            drawLine(image, cross.get(0), cross.get(1), colour);

        }

    }
    
    /** 
     * Use this FundamentalMatrix and it's own set of point correspondences
     * to draw epipolar lines into the image passed as argument.
     */
    public static void drawEpipolarLines(BufferedImage image,
                                         FundamentalMatrix FUND, 
                                         int colour){
        drawEpipolarLines(image, FUND, FUND.points, colour);
    }

    /**
     * Draw a straight line between two coordinates in the desired colour
     */
    public static void drawLine(BufferedImage rgb, int[] coords1, int[] coords2, int colour) {

        // If point is at the same coordinates in both frames, r=0 and the
        // trig terms here are NaN. Skip the rest of the drawing method.
        if((coords2[0]-coords1[0] == 0) && (coords2[1]-coords1[1] == 0))
            return;        
        
        //+++ Get shift in each coordinate between frames +++//
        double delta_i = coords2[0] - coords1[0];
        double delta_j = coords2[1] - coords1[1];

        //+++ Absolute shift in position +++//
        double r = Math.sqrt(delta_j * delta_j + delta_i * delta_i);

        //+++ Use polar representation to draw a straight line between these points +++//
        double sin_theta = delta_j / r;
        double cos_theta = delta_i / r;

        for (double r_prime = 0; r_prime <= r; r_prime++) {

            double i = (r_prime * cos_theta) + coords1[0];
            double j = (r_prime * sin_theta) + coords1[1];

            try {
                rgb.setRGB((int) Math.rint(i), (int) Math.rint(j), colour);
            }
            catch (ArrayIndexOutOfBoundsException aioobe) {}

        }

    }

    /**
     * Add Gaussian noise to an image.
     * @param image     Image to be degraded.
     */
    public static int[][] addGaussianNoise(int[][] image, double sigma){

        int[][] noisy = new int[image.length][image[0].length];

        for(int i=0; i<image.length; i++){
            for(int j=0; j<image[i].length; j++){

                //+++ Additive noise drawn from a Gaussian +++//
                noisy[i][j] = image[i][j] + (int)Math.rint(noise.nextGaussian()*sigma);

                //+++ Clamp values to avoid intensities going out of range +++//
                if(noisy[i][j]<0)   noisy[i][j] = 0;
                if(noisy[i][j]>255) noisy[i][j] = 255;

            }
        }

        return noisy;
    }

    /**
     * Add Gaussian blur to an image.
     * @param image     Image to be degraded.
     */
    public static  int[][] addGaussianBlur(int[][] image, int order, double sigma){

        if(order==1 || sigma<1E-9) return image;    // blurring has no effect

        if(order%2==0||order<1) System.err.println("Warning! Kernel width must be odd and positive");

        //+++ Make kernel for convolution +++//
        double[][] kernel = new double[order][order];

        double normal = 0;

        for(int i=0; i<order; i++){
            for(int j=0; j<order; j++){

                //+++ Coordinates relative to centre of mask +++//
                double x = i - (order-1)/2;
                double y = j - (order-1)/2;

                kernel[i][j] = Math.exp(-1.0*(x*x + y*y)/(2*sigma*sigma));
                normal       += kernel[i][j];

            }
        }

        //+++ Normalise mask +++//
        for(int i=0; i<order; i++){
            for(int j=0; j<order; j++){
                kernel[i][j] /= normal;
            }
        }

        //+++ Figure out border region that gradient operator can't be applied to +++//
        int k_border = (kernel.length - 1) / 2;

        //+++ Set up output gradient image. All elements are initialised to zero +++//
        int[][] convolved = new int[image.length][image[0].length];

        //+++ Now loop over pixels that kernel can be applied to. Borders will be left at zero +++//
        for (int x = k_border; x < image.length - k_border; x++) {
            for (int y = k_border; y < image[0].length - k_border; y++) {

                //+++ Loop over kernel elements +++//
                for (int i = 0; i < kernel.length; i++) {
                    for (int j = 0; j < kernel[i].length; j++) {
                        convolved[x][y] += image[x - k_border + i][y - k_border + j] * kernel[i][j];
                    }
                }
            }
        }

        return convolved;
    }
}
