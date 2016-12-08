package data;

import data.DiscreteFunction1D;


/**
 * This class tests various methods associated with VaryingHistogram objects.
 * @author nickrowell
 */
public class DiscreteFunction1DTester {
    
    public static void main(String[] args){
        
        testIntegration();
        
    }
    

    
    
    
    
    
    public static void testIntegration(){
        
        // set bin sizes, contents and errors
        double[] bin_centres  = {1};
        double[] bin_widths   = {2};
        double[] bin_contents = {1};
        double[] bin_errors   = {1};
        
        
        // Get Histogram suitable for binning unit Gaussian rvs.
        DiscreteFunction1D func = new DiscreteFunction1D(bin_centres,
                                                         bin_widths,
                                                         bin_contents,
                                                         bin_errors);
        
        // Print histogram
        System.out.println(func.print(true));
        
        // Get integral
        double[] integral = func.integrate(-10, 10);
        
        System.out.println("\nIntegral of histogram = "+integral[0]+" +/- "+integral[1]);

    }
    
    
    
}
