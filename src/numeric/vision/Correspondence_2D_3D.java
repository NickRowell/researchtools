package numeric.vision;


import Jama.Matrix;

/**
 * Class to represent point in a single image, and the 3D scene point that it
 * is projected from.
 * 
 * @author nickrowell
 */
public class Correspondence_2D_3D {
    
    public Point2D p2d;
    public Point3D p3d;
    
    public Correspondence_2D_3D(Point2D p2d_in, Point3D p3d_in){
    
        // Copy references
        p2d = p2d_in;
        p3d = p3d_in;
        
    }
    
    /**
     * Project 3D location of point onto image plane. This sets the 'subsidiary'
     * coordinates for the point used in improving camera matrix estimates.
     * 
     * 
     */
    public void setSubsPoint(Matrix K, Matrix P){
                
        // Transform to frame of camera then project onto image plane
        p2d.x_sub       = K.times(P.times(p3d.x));
        p2d.x_sub.timesEquals(1.0/p2d.x_sub.get(2,0));

    }
    
    public String getStringX(){
        return p3d.toString();
    }  
}
