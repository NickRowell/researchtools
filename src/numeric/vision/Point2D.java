package numeric.vision;

import Jama.Matrix;



/**
 *
 * @author nickrowell
 */
public class Point2D {
        
    /** Homogenous coordinates. */
    public Matrix x = new Matrix(3,1);
    /** Covariance on image coordinates. */
    public Matrix cov = new Matrix(2,2);    
    /**
     * Normalised homogenous coordinates. Used to improve
     * eight point algorithm, and NOT equal to K^{-1}*x_prime.
     */
    public Matrix x_norm = new Matrix(3,1);
    /**
     * Ideal ('subsidiary') image coordinates obtained by back-projection of
     * triangulated 3D location. Set by call to triangulateInliers(), and 
     * updated by Gold Standard Algorithm to constrain solution for E.
     */
    public Matrix x_sub = new Matrix(3,1);
    
    /** Basic constructor */
    public Point2D(){}
    
    /** Constructor used to initialise points */
    public Point2D(Matrix _x, Matrix _cov){
        x   = _x.copy();
        cov = _cov.copy();
    }
    
    /** Constructor used to copy existing points */
    public Point2D(Matrix _x, Matrix _cov, Matrix _x_norm, Matrix _x_sub){
        this(_x, _cov);
        x_norm = _x_norm.copy();
        x_sub  = _x_sub.copy();
    }   
 
    /** Copy constructor */
    public Point2D copy(){
        return new Point2D(x,cov,x_norm,x_sub);
    }

    /** String representation */
    @Override
    public String toString(){
        return x.get(0,0)+" "+x.get(1,0)+" "+x.get(2,0);
    }
    
    /**
     * Get the probability density of observing the Point2D at the point
     * (i,j) in the image, given the current projection of the covariance
     * matrix into the image plane.
     */
    public double getPositionPDF(double i, double j){
        
        // Get sigmas from mean, squared
        double exponent = getSigmasFromMean2(i,j)*(-.5);
        
        double coeff = 1.0/(2*Math.PI*Math.sqrt(Math.abs(cov.det())));
                
        return coeff*Math.exp(exponent);
    
    }
    /**
     * Get offset of point i,j from mean position of Point2D in terms of 
     * number of standard deviations squared.
     */
    public double getSigmasFromMean2(double i, double j){
        
        Matrix x_diff = new Matrix(new double[][]{{i-x.get(0,0)},
                                                  {j-x.get(1,0)}});
        
        return x_diff.transpose().times(cov.solve(x_diff)).get(0, 0);    
    
    }   
    
    
}
