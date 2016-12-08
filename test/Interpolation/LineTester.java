package Interpolation;

import functions.Line;

/**
 * Class intended for testing various methods of Line class.
 * @author nickrowell
 */
public class LineTester {
    
    
    public static void main(String[] args){
    
        Line y = new Line(1,-1,2,-2);
        
        // Get definite integral at x=1.5
        double x0 = 2.0;
        double int_x = y.integrateWrtX(x0);
                
    
    }
    
    
    
}
