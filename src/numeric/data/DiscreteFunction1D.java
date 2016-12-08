package numeric.data;

import util.ArrayUtil;

/**
 * Class represents a 1D function as a set of discrete bins.
 *
 *
 * @author nrowell
 * @version $Id$
 */
public class DiscreteFunction1D {
	
    /** 
     * Array containing bin contents and error.
     */
    protected ErrorBin[] bins;
    
    /** 
     * Main constructor.
     */
    public DiscreteFunction1D(double[] binCentres, double[] binWidths, double[] binContents,
                              double[] binErrors) {

        // Sanity checks on array sizes
		if(binCentres.length != binWidths.length) {
			throw new IllegalArgumentException("RangeMap: number of bin centres ("+binCentres.length+")"
					+ " does not equal the number of bin widths ("+binWidths.length+")!");
		}
		if(!ArrayUtil.checkIncreasing(binCentres)) {
			throw new IllegalArgumentException("RangeMap: Bin centres not increasing!");
		}
        if(!ArrayUtil.checkNonOverlappingBins(binCentres, binWidths)) {
            throw new IllegalArgumentException("RangeMap: Illegal bin centres/sizes.");
        }
        
        // Create new array of Bins
        bins = new ErrorBin[binCentres.length];
        
        // Loop over all bins and create a new Bin for each
        for(int b=0; b<binCentres.length; b++) {
           bins[b] = new ErrorBin(binCentres[b], binWidths[b], binContents[b], binErrors[b]);
        }
        
    }
    
    /**
     * Copy constructor.
     * 
     * @param copyme
     * 	The instance of {@link DiscreteFunction1D} to copy.
     */
    public DiscreteFunction1D(DiscreteFunction1D copyme) {
        // Create new array of bins
        bins = new ErrorBin[copyme.bins.length];        
        
        // Copy centres, widths and contents to new bins
        for(int b=0; b<bins.length; b++) {
            bins[b] = new ErrorBin(copyme.bins[b]);
        }
    }
    
    /**
     * Set specific bin value and standard deviation.
     */
    public void setBin(int bin, double value, double std) {
        bins[bin].set(value, std);
    }
    
    /**
     * Get centred bin coordinate.
     */
    public double getBinCentre(int bin) {
        return bins[bin].centre;
    }
    
    /**
     * Get all bin centres in an array.
     */
    public double[] getBinCentres() {
        double[] centres = new double[bins.length];
        for(int b=0; b<bins.length; b++)
            centres[b] = getBinCentre(b);
        return centres;
    }
    
    /**
     * Get lower boundary of bin.
     */
    public double getBinLowerEdge(int bin) {
        return bins[bin].centre - bins[bin].width/2.0;
    }
    
    /**
     * Get upper boundary of bin.
     */
    public double getBinUpperEdge(int bin) {
        return bins[bin].centre + bins[bin].width/2.0;
    }
    
    /**
     * Get value of a bin.
     */
    public double getBinContents(int bin) {
        return bins[bin].value;
    }
    
    /**
     * Get all bin values in an array.
     */
    public double[] getBinContents() {
        double[] contents = new double[bins.length];
        for(int b=0; b<bins.length; b++)
            contents[b] = getBinContents(b);
        return contents;
    }
    
    /**
     * Get standard deviation on a bin.
     */
    public double getBinUncertainty(int BIN) {
        return bins[BIN].std;
    }
    
    /**
     * Get all bin standard deviations in an array.
     */
    public double[] getBinUncertainties() {
        double[] errors = new double[bins.length];
        for(int b=0; b<bins.length; b++)
            errors[b] = getBinUncertainty(b);
        return errors;
    }    
    
    /**
     * Get range of independent variable covered by this bin.
     */
    public double getBinWidth(int BIN) {
        return bins[BIN].width;
    }
    
    /**
     * Get all bin widths in an array.
     */
    public double[] getBinWidths() {
        double[] widths = new double[bins.length];
        for(int b=0; b<bins.length; b++)
            widths[b] = getBinWidth(b);
        return widths;
    }
    
    /** 
     * Get number of data points.
     */
    public int size() {
    	return bins.length;
    }

    /**
     * Get the range of the X axis, considering the widths of the bins.
     * @return
     * 	The lower and upper limits on the X axis
     */
    public double[] getRangeX() {
    	return new double[]{getBinLowerEdge(0), getBinUpperEdge(size()-1)};
    }
    
    /**
     * Get the range of the Y axis spanned by the data.
     * @return
     * 	The lower and upper limits on the Y axis
     */
    public double[] getRangeY() {
    	
    	double min =  Double.MAX_VALUE;
    	double max = -Double.MAX_VALUE;
    	
        for(int i=0; i<bins.length; i++)
        {
        	if(getBinContents(i) < min) {
        		min = getBinContents(i);
        	}
        	if(getBinContents(i) > max) {
        		max = getBinContents(i);
        	}
        }
        
        return new double[]{min,max};
    }
    
    /**
     * Print contents each bin. Optionally scale numbers to units
     * of "counts per whole unit" rather than "counts per bin".
     * 
     * @param perUnit
     * 	If true, then the returned quantities are per-unit rather than per-bin.
     * @return
     */
    public String print(boolean perUnit) {
                
        StringBuilder out = new StringBuilder();
        
        for(int i=0; i<bins.length; i++) {
            
            double A = perUnit ? bins[i].width : 1.0;
            
            out.append(getBinCentre(i)).append("\t")
                                       .append(getBinWidth(i))
                                       .append("\t")
                                       .append(getBinContents(i)/A)
                                       .append("\t")
                                       .append(getBinUncertainty(i)/A)
                                       .append("\n");

        }
        return out.toString();
    }   
    
    /** 
     * Convolve the histogram using the given kernel. not currently
     * implemented because I can't see right now how to do convolution
     * when bins are different sizes and non-contiguous.
     * 
     */
    public void convolveSym(double sigma) {
        throw new RuntimeException("Not yet implemented");
    }  
    
    /**
     * Look up contents of bin that contains the given point.
     * @param value
     * @return 
     */
    public double getBinContents(double value) {
        
        for (ErrorBin bin : bins) {
            if(bin.contains(value))
                return bin.value;
        }
        
        return 0.0;
        
    }
    
    /**
     * Integrate over the full width of one bin.
     * @param bin
     * 	Index of the bin to integrate.
     * @return
     * 	The integral [0] and standard error [1]
     */
    public double[] integrate(int bin) {
    	return bins[bin].integral();
    }
    
    /**
     * Integrate histogram between two limits.
     * @param x_min Lower limit on independent variable
     * @param x_max Upper limit on independent variable
     * @return  Integral in first element; standard deviation of integral in
     *          second element.
     */
    public double[] integrate(double x_min, double x_max) {
        
        double I      = 0;
        double sig2_I = 0;
        
        for (ErrorBin bin : bins) {
            double[] integral = bin.integral(x_min, x_max);
            
            // Sum integral of bin area
            I += integral[0];
            // integral method returns standard deviation on integral over
            // single bin, must sum these in quadrature to get correct 
            // uncertainty on total integral
            sig2_I += integral[1]*integral[1];
        }        
        
        return new double[]{I, Math.sqrt(sig2_I)};
    }
    
    
    /**
     * Treat function like a PDF and draw a random x value from it.
     */
    public double draw(double a, double b) {
        
        // Sanity check on range
        assert b >= a;
        
        // First, get integral of function over this range
        double integral = integrate(a,b)[0];
        
        // Integral of random X point that lies in this range
        double int_x = integral*Math.random();
        
        // Now find what value of X this integral corresponds to
        for(ErrorBin bin : bins)
        {
            
            // Integral of this bin over the full range [a:b]
            double int_bin = bin.integral(a, b)[0];
            
            // If this value is less than the integral out to X, then this
            // bin doesn't contain X.
            if(int_x > int_bin)
            {
                // This bin doesn't contain X. Must either keep running sum of
                // int_bin, or subtract int_bin from int_x every loop.
                int_x -= int_bin;
                continue;
            }
            
            // Bin DOES contain X.
            return bin.invertIntegral(int_x);
                
       }
       
       // Shouldn't ever reach this code if the above algorithm works.
       assert false;
       
       return 0;
    }
    
}

/**
 * Class representing a single bin containing a value and error.
 * @author nickrowell
 */
class ErrorBin {

	/**
	 * The centre of the bin.
	 */
    public double centre;
    
    /**
     * The width of the bin.
     */
    public double width;
    
    /**
     * The value in the bin.
     */
    public double value;
    
    /**
     * The standard deviation.
     */
    public double std;
    
    /**
     * Main constructor.
     * @param centre
	 * 	The centre of the bin.
     * @param width
     * 	The width of the bin.
     * @param value
     * 	The value in the bin.
     * @param std
     * 	The standard deviation.
     */
    public ErrorBin(double centre, double width, double value, double std) {
        this.centre = centre;
        this.width = width;
        this.value = value;
        this.std = std;
    }
    
    /**
     * Copy constructor.
     * @param copyme
     * 	The instance of {@link ErrorBin} to copy.
     */
    public ErrorBin(ErrorBin copyme) {
    	this.centre = copyme.centre;
        this.width = copyme.width;
        this.value = copyme.value;
        this.std = copyme.std;
    }

    /**
     * Does this bin contain the given value?
     * 
     * @param value
     * @return 
     */
    public boolean contains(double value) {
        // Absolute difference between value and centre of bin must be smaller
        // than half the width of the bin.
         return value >= centre - width/2.0 && value < centre + width/2.0;
    }
    
    public double getLowerEdge() {
    	return centre - width/2.0;
    }
    
    public double getUpperEdge() {
    	return centre + width/2.0;
    }

    
    /**
     * Get integral over full bin width, and standard deviation.
     * @return
     * 	The integral [0] and standard error [1]
     */
    public double[] integral() {
        return new double[]{width * value, width * std};
    }

    /**
     * Get integral over restricted range 
     * @param a Lower limit on integral
     * @param b Upper limit on integral
     */
    public double[] integral(double a, double b) {
        
        assert a<b;
        
        // Boundaries lie entirely outside range of bin
        if(a >= getUpperEdge() || b <= getLowerEdge())
            return new double[]{0.0,0.0};
    
        // Boundaries straddle range of bin: must clamp to bin edges to 
        // avoid overshooting.
        a = Math.max(a, getLowerEdge());
        b = Math.min(b, getUpperEdge());
        
        // Use clamped values to calculate integral over bin
        return new double[]{(b-a) * value, (b-a) * std};
    }
    
    /** 
     * Get point inside bin where integral (from left side of bin) equals 
     * the given value.
     */
    public double invertIntegral(double int_X) {
        
        if(int_X < 0 || int_X > integral()[0]) {
            throw new RuntimeException("X outside range of bin integral. X = "
                                       +int_X);
        }
        
        return getLowerEdge() + int_X/value;
    }
    
    /**
     * Set contents of bin to specific value.
     */
    public void set(double val, double std) {
        this.value               = val;
        this.std  = std;
    }
    
}