/**
 * Copyright (c) 2011 University of Dundee
 *
 * Name:
 * Histogram.java
 *
 * Purpose:
 * This class is useful for binning numerical values drawn from a 
 * continuum of values into a histogram.
 * 
 * Language:
 * Java
 *
 * Author:
 * Nicholas Rowell (nickrowell@computing.dundee.ac.uk)
 *
 */

package numeric.data;

public class Histogram {

    // Data range
    protected double MIN;
    protected double MAX;
    // Bin width
    protected double BIN_WIDTH;
    // Outliers trimmed (true) or placed in extreme bins (false)
    protected boolean TRIM;
    // Number of bins
    protected int N;
    // Array to accumulate histogram results
    protected double[] HIST;

    /** Default constructor. */
    public Histogram(){}
    
    
    /** Main constructor */
    public Histogram(double min, double max, double delta, boolean trim){

        MIN   = min;
        MAX   = max;
        BIN_WIDTH = delta;
        TRIM  = trim;

        N = (int)Math.rint((MAX-MIN)/BIN_WIDTH);
        
        HIST = new double[N];

    }

    /**
     * Add a particular quantity to the histogram in the bin that
     * VALUE falls in.
     */
    public void add(double binCoordinate, double quantityToAddToBin){

        // Get bin number
        int bin = (int)Math.floor((binCoordinate-MIN)/BIN_WIDTH);

        // Ignore outliers if trimming is desired
        if((binCoordinate<MIN || binCoordinate>MAX) && TRIM)
            return;

        add(bin,quantityToAddToBin);

    }
    /** Add a particular quantity to the histogram in the designated bin */
    public void add(int BIN, double QUANTITY){

        // Move outlier to closest bin
        if(BIN<0) BIN=0;
        if(BIN>N-1) BIN=N-1;

        // Accumulate bin
        HIST[BIN]+=QUANTITY;

    }
    
    /** 
     * Sometimes we just want to count the number of objects that fall into
     * each histogram bin. This overridden method provides a simpler way to 
     * do this. Each object contributes precisely 1 to the bin it falls in.
     */
    public void add(double VALUE){
        add(VALUE, 1);
    }
    
    /** Get number of bins. */
    public int getNumberOfBins(){ return N;}
    
    /** Get centred bin coordinate */
    public double getBinCentre(int BIN){
        return (MIN + (double)BIN*BIN_WIDTH + (BIN_WIDTH/2.0));
    }
    /** Get lower boundary of bin */
    public double getBinLowerEdge(int BIN){
        return (MIN + (double)BIN*BIN_WIDTH);
    }
    /** Get upper boundary of bin */
    public double getBinUpperEdge(int BIN){
        return (MIN + (double)BIN*BIN_WIDTH + BIN_WIDTH);
    }    
    /** Get width of bins */
    public double getBinWidth(){
        return BIN_WIDTH;
    }
    
    /** Get list containing objects in a bin */
    public double getBinContents(int BIN){
        return HIST[BIN];
    }
    
    /** Get the max value of the histogram. */
    public double getMax(){
          
        // Initialise maximum value to first bin.
        double max = getBinContents(0);

        for(int x=0; x<HIST.length; x++)
                if(max < getBinContents(x))
                    max = getBinContents(x);

        return max;  
    
    } 
    
    /**
     * Get the mean value of the histogram.
     * @return
     * 		The mean value of the histogram.
     */
    public double getMean() {
          
        // Get normalisation constant for histogram
        double normalisation = integrate();
        
        double mean = 0.0;

        for(int i=0; i<HIST.length; i++)
            
            // <x>  = SUM( x  *  P(x)  *  dx)
            mean += getBinCentre(i) * (getBinContents(i) / normalisation) * BIN_WIDTH;

        return mean;
    }
    
    /**
     * Get the median value of the histogram.
     * @return
     * 		The median value of the histogram.
     */
    public double getMedian() {
          
        // Get normalisation constant for histogram
        double normalisation = integrate();
        
        double integral = 0.0;

        for(int i=0; i<HIST.length; i++) {
        	
        	// PDF in this bin
        	double px = (getBinContents(i) / normalisation);
        	
            // <x>  = SUM(P(x)  *  dx)
            double bin = px * BIN_WIDTH;
            
            if(integral+bin > 0.5)
            {
            	// Median lies in the current bin
            	
            	// At what bin width does the integral reach 0.5?
            	double deltaIntegral = 0.5 - integral;
            	double deltaBin = deltaIntegral / px;
            	
            	// Point at which integral reaches 0.5 = median
            	return getBinLowerEdge(i)+deltaBin;
            }
            else
            {
            	integral += bin;
            }
        }
        throw new RuntimeException("Histogram.getMedian(): couldn't find median!");
    }
    
    /**
     * Multiply (in place) the contents of each bin by the given number.
     */
    public void multiply(double S){
        for(int i=0; i<HIST.length; i++)
            HIST[i] *= S;
    }
    
    /** 
     * Convolve the histogram using the given kernel. The kernel is assumed to
     * be symmetric and is reflected about kernel[0] to get two-sided 
     * distribution.
     */
    public void convolveSym(double[] kernel){
    
        double[] HIST_CONV = new double[N];
    
        // How wide are the kernel borders?
        int border = kernel.length - 1;
        
        // Loop over each bin.
        for(int i=0; i<N; i++){
        
            // Place kernel at this point, and loop over all elements
            for(int k=-border; k<=border; k++){
                
                // Check that point i+k lies within range of histogram before
                // adding contribution from this element.
                if(i+k >= 0 && i+k < N)
                    // Move quantity from bin i in original histogram to 
                    // bin i+k in convolved histogram:
                    HIST_CONV[i+k] += HIST[i] * kernel[Math.abs(k)];
            
            }
            
        }
        
        // Finished calculating convolved histogram. Copy this back to the
        // original array.
        System.arraycopy(HIST_CONV, 0, HIST, 0, HIST.length);        
    
    }
    
    
    
    /**
     * Integrate histogram over full range.
     * 
     * @return
     * 		The sum SIGMA_{n=0}^{n=N-1}(B(n)*dB)
     */
    public double integrate(){
        double integral = 0.0;
        for(double binContents: HIST)
            integral += binContents;
        integral *= BIN_WIDTH;
        return integral;
    }
    
    /**
     * Integrate histogram over restricted range.
     */
    public double integrate(double lower, double upper){
    
        // Sanity check
        if(lower - upper >= 0){
            throw new RuntimeException("Integration limits wrong.");
        }
        
        // Get bin numbers on boundaries
        int BIN_LOWER = (int)Math.floor((lower-MIN)/BIN_WIDTH);
        int BIN_UPPER = (int)Math.floor((upper-MIN)/BIN_WIDTH);
        
        // Further sanity checks on range
        if(BIN_LOWER>N-1) return 0;   // Lower limit above range of histogram
        if(BIN_UPPER<0) return 0;     // Upper limit below range of histogram
        
        // Restrict boundaries to histogram range to avoid trying to access
        // bins outside range.
        if(BIN_LOWER<0) BIN_LOWER=0;
        if(BIN_UPPER>N-1) BIN_UPPER=N-1;
        
        // Now perform integration
        double integral = 0.0;
        
        for(int b=BIN_LOWER; b<=BIN_UPPER; b++){
            integral += getBinContents(b) * BIN_WIDTH;
        }
        
        return integral;
        
    }
        
    /**
     * Get normalised version of histogram.
     * @return 
     */
    public Histogram getNormalised(){
        
        Histogram normalised = new Histogram();
        
        normalised.BIN_WIDTH = BIN_WIDTH;
        normalised.MIN       = MIN;
        normalised.MAX       = MAX;
        normalised.N         = N;
        normalised.TRIM      = TRIM;
        normalised.HIST      = new double[N];
        
        // Get normalisation constant
        double integral = integrate();
        
        for(int i=0; i<HIST.length; i++)
            normalised.HIST[i] = getBinContents(i)/integral;
        
        return normalised;
    }
    
    
    
        
    /** 
     * Get String representation of Histogram.
     * 
     * @param NORM  If true, histogram is normalised so that 
     *              SUM(hist)*BIN_WIDTH = 1.0.
     * 
     */
    public String print(boolean NORM){

        // Set normalisation constant to unity or integral of histogram 
        // depending on whether normalisation is desired.
        double integral = NORM ? integrate() : 1.0;
     
        StringBuilder out = new StringBuilder();

        for(int i=0; i<HIST.length; i++)
            out.append(getBinCentre(i)).append("\t")
                    .append(getBinContents(i)/integral).append("\n");

        return out.toString();
    }

 
    
    
}
