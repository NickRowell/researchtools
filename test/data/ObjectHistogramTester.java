package data;

import java.util.List;
import java.util.Random;

import stats.Statistics;
import data.Binnable;
import data.ObjectHistogram;


/**
 * 
 * 
 * @author nickrowell
 */
public class ObjectHistogramTester {
    
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
        
         // Get ObjectHistogram for storing Doubles
        ObjectHistogram<Double> hist = new ObjectHistogram<Double>(-5.25, 5.25, 0.5, false);       
        
        // Add one million objects to centre bin
        for(int N=0; N<100000; N++)
            hist.add(new Double(0.0), 0.0);        
        
        // Now get unit Gaussian error kernel
        double[] kernel = Statistics.getSymGausKernel(1.0, 0.5);
        
        // Convolve histogram...
        hist.convolveSym(kernel);
    
        // Print histogram
        System.out.println(hist.print(true));   
        
    }
    
    
    
    
    /**
     * This method simply adds a bunch of Double objects to an ObjectHistogram
     * binned according to a Gaussian random variable. It then prints the 
     * number of objects in each bin, so test that binning has worked correctly.
     */
    private static void testBinning(){
    
        // Get ObjectHistogram for storing Doubles
        ObjectHistogram<Double> hist = new ObjectHistogram<Double>(-5.25, 5.25, 0.5, false);       
        
        // set sigma
        double sigma = 1.0;
        // set mean
        double mean = 0.0;
        
        // Get Random & seed with current system time.
        Random random = new Random(System.currentTimeMillis());
            
        for(int N=0; N<100000; N++)
            hist.add(new Double(0.0), sigma*random.nextGaussian() + mean);
    
        // Print histogram
        System.out.println(hist.print(true));    
    
    }
    
    
    
    
    private static void testSameObjectsInTwoInstances(){
        
        // Get ObjectHistogram for storing TestObjects
        ObjectHistogram<TestObject> hist1 = new ObjectHistogram<TestObject>(-0.5, 0.5, 1.0, false);
        ObjectHistogram<TestObject> hist2 = new ObjectHistogram<TestObject>(-0.5, 0.5, 1.0, false);
    
        // Create a single TestObject and set its boolean field to false
        TestObject test = new TestObject(false);
        
        System.out.println("Initial boolean value = "+test.A);
        
        // Now add it to both ObjectHistograms
        hist1.add(test, 0.0);
        hist2.add(test, 0.0);
        
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
            testObject.switchA();;
        }
    
    }
    
    
    
}


class TestObject implements Binnable{

    boolean A;

    public TestObject(boolean a){ A = a;}
    
    public void switchA(){ A = !A;}

    public double getBinValue(){ return 0.0;}
    
}