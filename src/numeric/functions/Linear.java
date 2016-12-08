/**
 * Copyright (c) 2011 University of Dundee
 *
 * Name:
 * Linear.java
 *
 * Purpose:
 * This class is used to interpolate values within a 1D array.
 * 
 * Language:
 * Java
 *
 * Author:
 * Nicholas Rowell (nickrowell@computing.dundee.ac.uk)
 *
 */
package numeric.functions;

import java.util.LinkedList;
import java.util.List;

import util.ArrayUtil;

/**
 * Linear interpolate within a set of 1D data.
 * @author nickrowell
 */
public class Linear {
    
    /** 
     * Array of interpolating Line segments 
     */
    protected List<Line> lines = new LinkedList<Line>();
    
    /** 
     * Abscissa
     */
    public double[] X;
    
    /**
     * Ordinate
     */
    public double[] Y;
    
    /**
     * Main constructor for {@link Linear} interpolation object.
     * 
     * @param x
     * 	Ordered array of X points (min to max).
     * @param y
     * 	Corresponding Y points.
     */
    public Linear(double[] x, double[] y) throws RuntimeException
    {
        // Checks on arguments
        if(x.length != y.length)
            throw new RuntimeException("Unequal number of x and y values!");
        if(x.length < 2)
            throw new RuntimeException("Need more than "+x.length+" points for interpolation!");
        if(!ArrayUtil.checkIncreasing(x))
            throw new RuntimeException("Linear: x points must show monotonic rise");
        
        // Now loop over all coordinates and create line segments
        for(int l=0; l<x.length-1; l++) 
            lines.add(new Line(x[l],y[l],x[l+1],y[l+1]));
        
        X = x;
        Y = y;
    }
    
    /**
     * Alternative constructor that takes {@link List}s.
     * 
     * @param x
     * 	Ordered list of X points
     * @param y
     * 	Corresponding Y points
     */
    public Linear(List<Double> x, List<Double> y) {
    	this(ArrayUtil.toArray(x), ArrayUtil.toArray(y));
    }
    
    /**
     * Interpolate (or extrapolate) the Y value corresponding to the given X coordinate,
     * along with the first derivative at the interpolated point.
     * If the X coordinate lies outside the range of the input data, then the Y value is
     * computed by linear extrapolation using the extreme linear segments of the range.
     * 
     * @param x
     * 	The X coordinate
     * @return 
     * 	The interpolated (or extrapolated) Y value and first derivative. To find out if
     * the value was extrapolated, use {@link MonotonicLinear#isXExtrapolated(double)}.
     */
    public double[] interpolateY(double x) {
    	
    	double y=0, dydx=0;
    	
    	// Check for extrapolated values lying before start or after end of range
        if(x<=X[0]) {
        	y = lines.get(0).getY(x);
        	dydx = lines.get(0).getGradient();
        }
        else if(x>=X[X.length-1]) {
        	y = lines.get(X.length-2).getY(x);
        	dydx = lines.get(X.length-2).getGradient();
        }
        else {
	        for(Line line : lines) {
	            if(line.containsX(x)) {
	            	y = line.getY(x);
	            	dydx = line.getGradient();
	                break;
	            }
	        }
        }
        return new double[]{y, dydx};
    }
    
    /**
     * Compute the definite integral of the {@link Linear} function over the
     * range [-INF:x]. Regions where the function is not defined are assumed to
     * be zero.
     * 
     * @param x
     * 	Upper limit on the integration range.
     * @return
     * 	The definite integral of the {@link Linear} function over the
     * range [-INF:x]
     */
    public double integrateWrtX(double x){
    
        // Sum integral of all line segments
        double integral = 0.0;
        
        for(Line line : lines)
            integral += line.integrateWrtX(x);

        return integral;
    }
    
    /**
     * Interpolate the X value(s) corresponding to the given Y coordinate.
     * Extrapolation beyond the ends of the defined range is NOT performed.
     * 
     * As the {@link Linear} function is not necessarily monotonic there are multiple (or zero)
     * possible solutions. An array of all the solutions is returned.
     * 
     * @param y
     * 	The Y coordinate
     * @return 
     * 	The interpolated X value(s) and their first derivatives.
     */
    public double[][] interpolateX(double y)
    {
        List<double[]> solutions = new LinkedList<double[]>();
        for(Line line : lines) {
            if(line.containsY(y)) {
                solutions.add(new double[]{line.getX(y), 1.0/line.getGradient()});
            }
        }
        
        double[][] array = new double[2][solutions.size()];
    	for(int i=0; i<solutions.size(); i++) {
    		array[0][i] = solutions.get(i)[0];
    		array[1][i] = solutions.get(i)[1];
    	}
        return array;
    }
    
}
