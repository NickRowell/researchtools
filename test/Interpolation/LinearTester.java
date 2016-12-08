package Interpolation;

import functions.Linear;

/**
 * Test class for Linear interpolation.
 * @author nickrowell
 */
public class LinearTester {
        
    public static void main(String[] args){
    
        /** Test data */
        double[] x = new double[]{0,1,2,3};
        double[] y = new double[]{0,1,2,3};    
        
        Linear linear = new Linear(x,y);
        
        // interpolate a bunch of points along x axis
        for(double X = 0; X<11; X+=0.1){
            double Y = linear.interpolateY(X);
            System.out.println(X+" "+Y+" "+linear.integrateWrtX(X));
        }
        
        
        
        
        
    }
    
    
}
