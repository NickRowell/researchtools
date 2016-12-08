/**
 * Copyright (c) 2011 University of Dundee
 *
 * Name:
 * Histogram2D.java
 *
 * Purpose:
 * This class is very similar to Histogram.java but extended to 2 dimensions.
 * 
 * To do:
 * Implement 
 *  - integrate(double, double)
 *  - convolveSym(double[])
 *  - getMean()
 * 
 * Language:
 * Java
 *
 * Author:
 * Nicholas Rowell (nickrowell@computing.dundee.ac.uk)
 *
 */

package numeric.data;

public class Histogram2D {

    // Data range in first coordinate
    protected double MIN_X;
    protected double MAX_X;
    protected double MIN_Y;
    protected double MAX_Y;    
    // Bin width
    protected double BIN_WIDTH_X;
    protected double BIN_WIDTH_Y;    
    // Outliers trimmed (true) or placed in extreme bins (false)
    protected boolean TRIM;
    // Number of bins
    protected int N_X;
    protected int N_Y;
    // 2D array to accumulate histogram results
    protected double[][] HIST;

    /** Main constructor */
    public Histogram2D(double min_x, double max_x, double delta_x, 
                       double min_y, double max_y, double delta_y, 
                       boolean trim){
        
        // Basic sanity checks on bin parameters
        assert(min_x < max_x && delta_x > 0);
        assert(min_y < max_y && delta_y > 0);
        
        // Range & bins in X
        MIN_X   = min_x;
        MAX_X   = max_x;
        BIN_WIDTH_X = delta_x;
        // Range & bins in Y
        MIN_Y   = min_y;
        MAX_Y   = max_y;
        BIN_WIDTH_Y = delta_y;        
        
        TRIM  = trim;

        N_X = (int)Math.rint((MAX_X-MIN_X)/BIN_WIDTH_X);
        N_Y = (int)Math.rint((MAX_Y-MIN_Y)/BIN_WIDTH_Y);
        
        HIST = new double[N_X][N_Y];

    }

    /**
     * Add a particular quantity to the histogram in the bin that
     * VALUE falls in.
     */
    public void add(double X, double Y, double QUANTITY){

        // If point lies outside range and outliers are discarded, return now.
        if((X<MIN_X || X>MAX_X || Y<MIN_Y || Y>MAX_Y) && TRIM)
            return;       
        
        // Get bin number
        int BIN_X = (int)Math.floor((X-MIN_X)/BIN_WIDTH_X);
        int BIN_Y = (int)Math.floor((Y-MIN_Y)/BIN_WIDTH_Y);

        add(BIN_X, BIN_Y, QUANTITY);

    }
    /** Add a particular quantity to the histogram in the designated bin */
    public void add(int BIN_X, int BIN_Y, double QUANTITY){

        // Move outliers to closest bin
        if(BIN_X<0) BIN_X=0;
        if(BIN_X>N_X-1) BIN_X=N_X-1;
        if(BIN_Y<0) BIN_Y=0;
        if(BIN_Y>N_Y-1) BIN_Y=N_Y-1;
        
        // Accumulate bin
        HIST[BIN_X][BIN_Y]+=QUANTITY;

    }
    
    /** 
     * Sometimes we just want to count the number of objects that fall into
     * each histogram bin. This overridden method provides a simpler way to 
     * do this. Each object contributes precisely 1 to the bin it falls in.
     */
    public void add(double X, double Y){
        add(X, Y, 1);
    }
    
    /** Get total number of bins. */
    public int getN(){ return N_X*N_Y;}
    /** Get number of X bins. */
    public int getNX(){ return N_X;}
    /** Get number of Y bins. */
    public int getNY(){ return N_Y;}
    
    /** Get centred bin coordinate */
    public double[] getBinCentre(int BIN_X, int BIN_Y){
        return new double[]{(MIN_X + (double)BIN_X*BIN_WIDTH_X + (BIN_WIDTH_X/2.0)),
                            (MIN_Y + (double)BIN_Y*BIN_WIDTH_Y + (BIN_WIDTH_Y/2.0))};
    }
    /** Get lower boundary of bin in X direction. */
    public double getBinLowerEdgeX(int BIN_X){
        return (MIN_X + (double)BIN_X*BIN_WIDTH_X);
    }
    /** Get lower boundary of bin in Y direction. */
    public double getBinLowerEdgeY(int BIN_Y){
        return (MIN_Y + (double)BIN_Y*BIN_WIDTH_Y);
    }
    /** Get upper boundary of bin in X direction. */
    public double getBinUpperEdgeX(int BIN_X){
        return (MIN_X + (double)BIN_X*BIN_WIDTH_X + BIN_WIDTH_X);
    }
    /** Get upper boundary of bin in Y direction. */
    public double getBinUpperEdgeY(int BIN_Y){
        return (MIN_Y + (double)BIN_Y*BIN_WIDTH_Y + BIN_WIDTH_Y);
    }   

    
    
    /** Get list containing objects in a bin */
    public double getBinContents(int BIN_X, int BIN_Y){
        return HIST[BIN_X][BIN_Y];
    }

    
    /** Get the max value of the histogram. */
    public double getMax(){
          
        // Initialise maximum value to first bin.
        double max = getBinContents(0,0);

        for(int x=0; x<HIST.length; x++)
            for(int y=0; y<HIST[x].length; y++)
                if(max < getBinContents(x,y))
                    max = getBinContents(x,y);

        return max;  
    
    }   
    
    
    
    /** 
     * Get the mean value of the histogram.
     * 
     * Can mean be obtained for bivariate histogram?
     * 
     * Could probably get an expectation value for each coordinate - need to 
     * look this up somewhere...
     */
//    public double getMean(){
//          
//        // Get normalisation constant for histogram
//        double integral = integrate();
//        
//        double mean = 0.0;
//
//        for(int x=0; x<HIST.length; x++)
//            
//            for(int y=0; y<HIST[x].length; y++)
//            
//                // <x>  = SUM( x  *  P(x)  *  dx)
//                mean += getBinCentre(x) * (getBinContents(x,y) / integral) * BIN_WIDTH_X * BIN_WIDTH_Y;
//
//        return mean;  
//    
//    }
    
    /** Multiply (in place) the contents of each bin by the given number */
    public void multiply(double S){
        for(int x=0; x<HIST.length; x++)
            for(int y=0; y<HIST[x].length; y++)
                HIST[x][y] *= S;
    }
    
    /** 
     * Convolve the histogram using the given kernel. The kernel is assumed to
     * be symmetric and is reflected about kernel[0] to get two-sided 
     * distribution.
     * 
     * Need to extend this method to 2 dimensions at some point.
     * 
     */
//    public void convolveSym(double[] kernel){
//    
//        double[] HIST_CONV = new double[N];
//    
//        // How wide are the kernel borders?
//        int border = kernel.length - 1;
//        
//        // Loop over each bin.
//        for(int i=0; i<N; i++){
//        
//            // Place kernel at this point, and loop over all elements
//            for(int k=-border; k<=border; k++){
//                
//                // Check that point i+k lies within range of histogram before
//                // adding contribution from this element.
//                if(i+k >= 0 && i+k < N)
//                    // Move quantity from bin i in original histogram to 
//                    // bin i+k in convolved histogram:
//                    HIST_CONV[i+k] += HIST[i] * kernel[Math.abs(k)];
//            
//            }
//            
//        }
//        
//        // Finished calculating convolved histogram. Copy this back to the
//        // original array.
//        System.arraycopy(HIST_CONV, 0, HIST, 0, HIST.length);        
//    
//    }
    
    
    
    /** Integrate histogram over full range */
    public double integrate(){
        double integral = 0.0;
        for(double[] HIST_X : HIST)
            for(double binContents : HIST_X)
                integral += binContents * BIN_WIDTH_X * BIN_WIDTH_Y;
        return integral;
    }
    
    /** Integrate histogram over restricted range */
//    public double integrate(double lower, double upper){
//    
//        // Sanity check
//        if(lower - upper >= 0){
//            throw new RuntimeException("Integration limits wrong.");
//        }
//        
//        // Get bin numbers on boundaries
//        int BIN_LOWER = (int)Math.floor((lower-MIN)/BIN_WIDTH);
//        int BIN_UPPER = (int)Math.floor((upper-MIN)/BIN_WIDTH);
//        
//        // Further sanity checks on range
//        if(BIN_LOWER>N-1) return 0;   // Lower limit above range of histogram
//        if(BIN_UPPER<0) return 0;     // Upper limit below range of histogram
//        
//        // Restrict boundaries to histogram range to avoid trying to access
//        // bins outside range.
//        if(BIN_LOWER<0) BIN_LOWER=0;
//        if(BIN_UPPER>N-1) BIN_UPPER=N-1;
//        
//        // Now perform integration
//        double integral = 0.0;
//        
//        for(int b=BIN_LOWER; b<=BIN_UPPER; b++){
//            integral += getBinContents(b) * BIN_WIDTH;
//        }
//        
//        return integral;
//        
//    }
        
        
    /** Get String representation of Histogram2D, optionally normalised. */
    public String print(boolean NORM){

        
        // Set normalisation constant to unity or integral of histogram 
        // depending on whether normalisation is desired.
        double integral = NORM ? integrate() : 1.0;
     
        
        StringBuilder out = new StringBuilder();

        for(int x=0; x<HIST.length; x++){
            for(int y=0; y<HIST[x].length; y++){
                
                // (X,Y) coordinates of centre of this bin:
                double[] XY = getBinCentre(x,y);
                
                out.append(String.format("%1.2f",XY[0])).append("\t").append(String.format("%1.2f",XY[1])).append("\t")
                        .append(getBinContents(x,y)/integral).append("\n");
                
            }
            
            // New line between rows
            out.append("\n");
        
        }
        
        return out.toString();
    }

 
    
    
}
