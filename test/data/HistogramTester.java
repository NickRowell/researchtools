package data;

import java.util.Random;

import data.Histogram;


/**
 * This class uses native Java methods to draw unit-Gaussian-distributed
 * random variables which are then binned using the Data.Histogram class.
 * The results are printed at the end so that the binned, normalised 
 * histogram can be compared to the known distribution.
 * 
 * Printed distribution should match:
 * 
 * p(x) = 1/(S*sqrt(2*PI)) * exp(-(x-mean)^2/(2*S^2))
 * 
 * where S = 1 and mean = 0.
 * 
 * @author nickrowell
 */
public class HistogramTester {
    
    public static void main(String[] args){
    
        // set sigma
        double sigma = 1.0;
        // set mean
        double mean = 0;
        
        // Get Random & seed with current system time.
        Random random = new Random(System.currentTimeMillis());
    
        // Get Histogram suitable for binning unit Gaussian rvs.
        Histogram hist = new Histogram(-5.0*sigma + mean, 5.0*sigma + mean, 0.1*sigma, false);
        
        for(int N=0; N<1000000; N++)
            hist.add(sigma*random.nextGaussian() + mean);
    
        // Print histogram
        System.out.println(hist.print(true));
        
        System.out.println("\nMean (by histogram integration) = "+hist.getMean());
        // Print gnuplot function for Gaussian
        System.out.println("\nf(x) = (1/("+sigma+" * sqrt(2*pi)))*exp(-(x-"+mean+")**2/(2*("+sigma+"**2)))");
        
    }
    
    
    
}
