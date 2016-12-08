package numeric.vision;

import Jama.*;



/**
 * Class to represent correspondences between points in two images, and the 
 * point in the scene that they are projected from.
 * 
 */
public class Correspondence_2D_2D_3D {

    /** Image point in un-primed frame */
    public Point2D p2d;
    /** Image point in primed frame */
    public Point2D p2d_prime;
    /** World point */
    public Point3D p3d;
    
    
    /** Basic constructor */
    public Correspondence_2D_2D_3D(){

        p2d = new Point2D();
        p2d_prime = new Point2D();
        p3d = new Point3D();
                        
    }
    
    
    /**
     * Constructor for point correspondences with no link to existing Point3D.
     */
    public Correspondence_2D_2D_3D(Point2D p, Point2D p_prime){

        p2d = p;
        p2d_prime = p_prime;
        p3d = new Point3D();
                        
    }
    
    /**
     * Constructor for point correspondences with Point3D object specified too.
     */
    public Correspondence_2D_2D_3D(Point2D p, Point2D p_prime, Point3D X){
        p2d = p;
        p2d_prime = p_prime;
        p3d = X;    
    }   
    /** Deep copy a Correspondence_2D_2D_3D to a new object */
    public Correspondence_2D_2D_3D(Correspondence_2D_2D_3D copyme){

        p3d = copyme.p3d.copy();
        p2d = copyme.p2d.copy();
        p2d_prime = copyme.p2d_prime.copy();
        
    }
    
    
    
    /** Set the 3D position of this point. */
    public void setX(Matrix XX){ p3d.x = XX.copy();}
    /** Set single component of 3D position of this point. */
    public void setXi(int i, double val){ p3d.x.set(i, 0, val);}
    
    /** Get single component of 3D position of this point. */
    public double getXi(int i){ return p3d.x.get(i, 0);}    
    
    
    /** Get the 3D position of this point. */
    public Matrix getX(){ return p3d.x;}    



    /**
     * Format the normalised coordinates for this point into a row of the
     * matrix used to estimate the fundamental matrix.
     */
    public double[] getRowEntry(){
    
          return new double[]{p2d_prime.x_norm.get(0, 0) * p2d.x_norm.get(0, 0),
                              p2d_prime.x_norm.get(0, 0) * p2d.x_norm.get(1, 0),
                              p2d_prime.x_norm.get(0, 0),
                              p2d_prime.x_norm.get(1, 0) * p2d.x_norm.get(0, 0),
                              p2d_prime.x_norm.get(1, 0) * p2d.x_norm.get(1, 0),
                              p2d_prime.x_norm.get(1, 0),
                              p2d.x_norm.get(0, 0),
                              p2d.x_norm.get(1, 0),
                              1};
    }


    public String getStringX(){
        return p3d.toString();
    }




    // Methods used by GoldStandardAlgorithm to improve triangulation:

    /**
     * Return 're-projection error ' : sum of squared distances between observed
     * and subsidiary points.
     */
    public double getReprojectionError(){
        double d       = p2d.x.minus(p2d.x_sub).normF();
        double d_prime = p2d_prime.x.minus(p2d_prime.x_sub).normF();
        return d*d + d_prime*d_prime;
    }

    /**
     * Project 3D location of point onto image plane. This sets the 'subsidiary'
     * coordinates for the point used in improving camera matrix estimates.
     * 
     * 
     */
    public void setSubsPoints(Matrix K, Matrix P, Matrix K_prime, Matrix P_prime){
        
        // Transform to frame of each camera then project onto image plane
        p2d.x_sub       = K.times(P.times(p3d.x));
        p2d_prime.x_sub = K_prime.times(P_prime.times(p3d.x));
        p2d.x_sub.timesEquals(1.0/p2d.x_sub.get(2,0));
        p2d_prime.x_sub.timesEquals(1.0/p2d_prime.x_sub.get(2,0));

    }




}