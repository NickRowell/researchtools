package data;

import java.util.List;
import java.util.Random;

import data.Binnable;
import data.ObjectVaryingHistogram;




/**
 * 
 * 
 * @author nickrowell
 */
public class ObjectVaryingHistogramTester {
    
    public static void main(String[] args){
    
        //testBinning();
        //System.out.println();
        //testConvolution();
        testSameObjectsInTwoInstances();
    }
    
    /**
     * This method places a million objects into the central histogram bin,
     * then uses a Gaussian error kernel to convolve the histogram, and prints
     * the result so the validity of the convolution operation can be 
     * checked.
     */
    private static void testConvolution(){
        
        // Bin centres
        double[] centres = new double[]{-4.5,-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5,4.5};
        // Bin widths
        //double[] widths = new double[]{1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
        double[] widths = new double[]{0.5,0.5,0.5,0.5,1.0,1.0,0.5,0.5,0.5,0.5};       
        
        // Get ObjectHistogram for storing Doubles
        ObjectVaryingHistogram<DoubleBin> hist = new ObjectVaryingHistogram<DoubleBin>(centres, widths);       
        
        // Add one million objects to centre bin
        for(int N=0; N<100000; N++)
            hist.add(new DoubleBin(0.0));        
        
        // Convolve histogram...
        hist.convolveSym(1.0);
    
        // Print histogram
        System.out.println(hist.print(true));   
        
    }
    
    
    
    
    /**
     * This method simply adds a bunch of Double objects to an ObjectHistogram
     * binned according to a Gaussian random variable. It then prints the 
     * number of objects in each bin, so test that binning has worked correctly.
     */
    private static void testBinning(){
    
        // Bin centres
        double[] centres = new double[]{-4.5,-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5,4.5};
        // Bin widths
        //double[] widths = new double[]{1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
        double[] widths = new double[]{0.5,0.5,0.5,0.5,1.0,1.0,0.5,0.5,0.5,0.5};
        
        // Get ObjectHistogram for storing Doubles
        ObjectVaryingHistogram<DoubleBin> hist = new ObjectVaryingHistogram<DoubleBin>(centres,widths);       
        
        // set sigma
        double sigma = 1.0;
        // set mean
        double mean = 0.0;
        
        // Get Random & seed with current system time.
        Random random = new Random(System.currentTimeMillis());
            
        for(int N=0; N<100000; N++)
            hist.add(new DoubleBin(sigma*random.nextGaussian() + mean));
    
        // Print histogram
        System.out.println(hist.print(true));    
    
    }
    
    
    
    private static void testSameObjectsInTwoInstances(){
        
        double[] bin_centres = {0.0};
        double[] bin_widths  = {1.0};
        
        
        // Get ObjectHistogram for storing TestObjects
        ObjectVaryingHistogram<TestObject> hist1 = new ObjectVaryingHistogram<TestObject>(bin_centres,bin_widths);
        ObjectVaryingHistogram<TestObject> hist2 = new ObjectVaryingHistogram<TestObject>(bin_centres,bin_widths);
    
        // Create a single TestObject and set its boolean field to false
        TestObject test = new TestObject(false);
        
        System.out.println("Initial boolean value = "+test.A);
        
        // Now add it to both ObjectHistograms
        hist1.add(test);
        hist2.add(test);
        
        // Check that it has been added to each histogram:
        System.out.println("hist1: "+hist1.print(false));
        System.out.println("hist2: "+hist2.print(false));
        
        // Switch the boolean flag of object in first histogram:
        switchBoolean(hist1.getBinContents(0));
//        for (TestObject object : hist1.getBinContents(0)) {
//            object.switchA();
//        }
        
        // Now check status of all references to 'test'
        System.out.println("TestObject in hist1: "+hist1.getBinContents(0).get(0).A);
        System.out.println("TestObject in hist2: "+hist2.getBinContents(0).get(0).A);
        System.out.println("Original object: "+test.A);
        
    
    }
    
    private static void switchBoolean(List<TestObject> list){
        
        for (TestObject testObject : list) {
            testObject.switchA();
        }
    
    }
    
    
}
class DoubleBin implements Binnable{

    
    double value;
    
    public DoubleBin(double val){ value = val;}
    
    public double getBinValue() {
        return value;
    }


}

