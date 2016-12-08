package data;

import java.util.Random;

import data.Histogram;
import data.VaryingHistogram;


/**
 * This class tests various methods associated with VaryingHistogram objects.
 * @author nickrowell
 */
public class VaryingHistogramTester {
    
    public static void main(String[] args){
        
        testDraw();
        
    }
    
    
    
    public static void testDraw(){
        
        // set parameters of varyinghistogram
        double[] bin_centres = {0.5, 1.5};
        double[] bin_widths  = {1,1};
        
        // Get Histogram suitable for binning unit Gaussian rvs.
        VaryingHistogram hist = new VaryingHistogram(bin_centres,bin_widths);   
    
        hist.add(0.5, 2.0);
        hist.add(1.5, 1.0);
        
        // Print histogram
        System.out.println(hist.print(true));      
        
        // Make a new Histogram for storing values drawn from this varyinghistogram
        Histogram hist2 = new Histogram(-1,4,0.1,true);
        
        for(int i=0; i<100000; i++)
            hist2.add(hist.draw(-10, 10));
        
        // Print histogram
        System.out.println(hist2.print(true));        
        
    }
    
    
    
    
    
    public static void testIntegration(){
        
        // set parameters of varyinghistogram
        double[] bin_centres = {-1,0,1,1.75,2.25};
        double[] bin_widths  = {1,1,1,0.5,0.5};
        
        
        // Get Histogram suitable for binning unit Gaussian rvs.
        VaryingHistogram hist = new VaryingHistogram(bin_centres,bin_widths);
        
        hist.add(-1, 1.0);
        hist.add( 0, 1.0);
        hist.add( 1, 1.0);
        hist.add(1.75, 1.0);
        hist.add(2.25, 1.0);
        
        // Print histogram
        System.out.println(hist.print(true));
        
        System.out.println("\nIntegral of histogram = "+hist.integrate(-10, 10));

    }
    
    
    
}
