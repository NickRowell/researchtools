package numeric.functions;

import util.ArrayUtil;

/**
 *
 * MonotonicLinear class represents special case of Linear interpolation
 * objects, where linear function is strictly monotonic. This means that
 * the function can be uniquely inverted for X given Y.
 * 
 * This class also automatically extrapolates beyond the range of the data,
 * and includes methods to test whether a given value (x or y) is within the
 * range of the data or lies in extrapolated region.
 * 
 * This class is most useful for interpolating tabular data, of e.g.
 * White Dwarf cooling times or main sequence lifetimes.
 * 
 * @author nickrowell
 */
public class MonotonicLinear extends Linear implements MonotonicFunction1D
{
    
    /** 
     * There are two possible types of monotonic function, those that
     * increase and those that decrease. This enum is used to identify
     * instances of each type.
     */
    public static enum Type{INCREASING,DECREASING};
    

    /** 
     * String identifier (e.g. mass, metallicity).
     */
    public String name;
    
    /**
     * Type of monotonic function (increasing/decreasing). 
     */
    Type type;
    
    /**
     * Main constructor.
     * 
     * @param x
     * 	The X coordinates
     * @param y
     * 	The corresponding Y coordinates
     * @param name
     * 	A suitable name to associate with this instance
     * @param type
     * 	The type (INCREASING/DECREASING) of the monotonicity
     */
    public MonotonicLinear(double[] x, double[] y, String name, Type type)
    {
        // constructor of superclass tests for equal number of x and y points,
        // and that x points show monotonic rise.
        super(x,y);
        
        // Now check for monotonic increase/decrease in y points
        switch(type)
        {
            case INCREASING: 
                if(!ArrayUtil.checkIncreasing(y))
                    throw new IllegalArgumentException("MonotonicLinear: "
                            + "values don't show monotonic increase!");
                break;
                
            case DECREASING: 
                if(!ArrayUtil.checkDecreasing(y))
                    throw new IllegalArgumentException("MonotonicLinear: "
                           + "values don't show monotonic decrease!");
                break;
        }
        
        this.name = name;
        this.type = type;
    }
    
    /**
     * Constructor that infers the sense of the monotonicity from the data.
     * 
     * @param x
     * 	The X coordinates
     * @param y
     * 	The corresponding Y coordinates
     */
    public MonotonicLinear(double[] x, double[] y)
    {
        // constructor of superclass tests for equal number of x and y points,
        // and that x points show monotonic rise.
        super(x,y);
        
        // Check if y increases or decreases monotonically
        boolean isIncreasing = ArrayUtil.checkIncreasing(y);
        boolean isDecreasing = ArrayUtil.checkDecreasing(y);
        
        // If both are false, then function is not monotonic.
        if(!isIncreasing && !isDecreasing)
        {
        	StringBuilder message = new StringBuilder("MonotonicLinear: "
                    + "Y values not monotonic!\n");
        	
        	for(int i=0; i<x.length; i++) {
        		message.append(x[i]).append("\t").append(y[i]).append("\n");
        	}
        	
        	throw new IllegalArgumentException(message.toString());
        }
        
        // Sanity check: both conditions should never be true
        if(isIncreasing && isDecreasing)
        	throw new IllegalArgumentException("MonotonicLinear: "
                    + "Y values evaluated as both increasing and decreasing!");
        
        // Only one of isIncreasing and isDecreasing is true
        type = (isIncreasing) ? Type.INCREASING : Type.DECREASING;
        
    }
    
    
    /**
     * Interpolate the X value at the corresponding Y value. The monotonicity
     * ensures that there is a single unique solution for X
     * 
     * @param y
     * 	The Y value
     * @return 
     * 	The cooresponding interpolated X value and the first derivative wrt y
     */
    public double[] interpolateUniqueX(double y) {
        
    	// Line segment to use for interpolation/extrapolation
    	Line lineToUse = null;
    	
    	if(y < Math.min(Y[0], Y[Y.length-1])) {
    		// Extrapolation
    		if(type==Type.INCREASING) {
    			lineToUse = lines.get(0);
    		}
    		else {
    			lineToUse = lines.get(Y.length-2);
    		}
    	}
    	else if(y > Math.max(Y[0], Y[Y.length-1])) {
    		// Extrapolation
    		if(type==Type.INCREASING) {
    			lineToUse = lines.get(Y.length-2);
    		}
    		else {
    			lineToUse = lines.get(0);
    		}
    	}
    	else {
    		// Interpolation
    		for(Line line : lines) {
                if(line.containsY(y)) {
                	lineToUse = line;
                }
            }
    	}
    	
    	return new double[]{lineToUse.getX(y), 1.0/lineToUse.getGradient()};
    }
    
    /**
     * Check if the given X coordinate requires extrapolation of the input data, i.e. if
     * it lies beyond the range of the X coordinates used to construct this instance.
     * @param x
     * 	The X value
     * @return 
     * 	True, if the given X coordinate lies outside of the range of the input X coordinates
     */
    public boolean isXExtrapolated(double x) {
        return (x < X[0] || x > X[X.length-1]);
    }
    
    /**
     * Check if the given Y coordinate requires extrapolation of the input data, i.e. if
     * it lies beyond the range of the Y coordinates used to construct this instance.
     * @param y
     * 	The Y value
     * @return 
     * 	True, if the given Y coordinate lies outside of the range of the input Y coordinates
     */
    public boolean isYExtrapolated(double y){
        
        // Check if value y lies beyond either end of function. The test 
        // depends on whether the function starts low (INCREASING) or
        // starts high (DECREASING).
        switch(type){
            case INCREASING: return  (y < Y[0] || y > Y[Y.length-1]);
            case DECREASING: return  (y > Y[0] || y < Y[Y.length-1]);
            default: return false;
        }
    }

	@Override
	public double[] getY(double x) {
		return interpolateY(x);
	}

	@Override
	public double[] getX(double y) {
		return interpolateUniqueX(y);
	}    
    
}
