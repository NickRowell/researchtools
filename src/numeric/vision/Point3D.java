package numeric.vision;

import Jama.Matrix;


/**
 *
 * @author nickrowell
 */
public class Point3D {
    
    /** Homogenous position vector (4,1) in world frame */
    Matrix x;

    /** Position covariance matrix (3,3) */
    Matrix cov;
    
    /** Basic constructor */
    public Point3D(){}
    
    /** Main constructor */
    public Point3D(Matrix _x, Matrix _cov){
        x   = _x.copy();
        cov = _cov.copy();
    }
    /** Copy constructor */
    public Point3D copy(){
        return new Point3D(x,cov);
    }
    
    /** Return column vector of xyz coordinates */
    public Matrix getParameters(){
        return x.getMatrix(new int[]{0,1,2}, new int[]{0});
    }
    
    /** Set xyz coordinates */
    public void setParameters(Matrix delta){
        x.setMatrix(new int[]{0,1,2}, new int[]{0}, delta);
    }

    @Override
    public String toString(){
        return ""+x.get(0,0)+" "+x.get(1,0)+" "+x.get(2,0)+" "+x.get(3,0);
    }
    
    
    
    
}
