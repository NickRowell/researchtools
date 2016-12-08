package fitting;

import fitting.Ellipse;

/**
 * Class intended for testing out various Ellipse operations.
 * 
 * @author nickrowell
 */
public class EllipseTester {
   
    public static void main(String[] args){
        
        Ellipse e = new Ellipse(new double[]{2, 1, 2, 5, 4, 1});
        
        System.out.println(e.toString());
        
        
    }
    
    
}
