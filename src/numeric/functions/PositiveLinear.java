
package numeric.functions;

/**
 * This class is an extension of {@link Linear} interpolation method for the case when
 * all y data values are positive. This means that the function can be 
 * treated like a PDF, and this class simply adds a single method that allows
 * x values to be 'drawn' from the piecewise linear function.
 * 
 * @author nickrowell
 */
public class PositiveLinear extends Linear {
    
    /** 
     * Constructor same as for superclass, but with extra check that all
     * y values are positive or zero.
     */
    public PositiveLinear(double[] x, double[] y) throws RuntimeException{
        super(x,y);
        
        for (double yi : y) {
            if(yi<0) {
                throw new RuntimeException("Negative Y value! Y = "+yi);
            }
        }
    }
    
    /**
     * Draw random value from full range of function.
     * @return
     */
    public double drawX() {
    	double xmin = X[0];
    	double xmax = X[X.length-1];
    	return drawX(xmin, xmax);
    }
    
    /**
     * Treat piecewise linear function like a pdf, and draw a value from it
     * over restricted range.
     * 
     * TODO: record integral of function so that if we repeatedly draw values over the same
     * range we avoid recomputing the normalisation, which will be slow.
     * 
     */
    public double drawX(double a, double b) {
        
        // Sanity check on range
        assert b >= a;
        
        // First, get integral of function between these limits
        double int_a = integrateWrtX(a);
        double int_b = integrateWrtX(b);
        
        // Integral of random X point between these limits
        double int_x = int_a + (int_b - int_a)*Math.random();
        
        // Now find what value of X this integral corresponds to
        double totalInt = 0.0;
        
        for(Line line : lines){
            
            // Add contribution of this line segment to total integral
            totalInt += line.integralWrtX();
            
            // If total integral is now larger than integral out to X,
            // X lies in this line segment.
            if(totalInt >= int_x){
                
                // Get point inside this line segment that has an integral
                // equal to int_x - (totalInt - line.integral())
                
                // Now invert the integral to get the value of x at which
                // this integrated value occurs. Note that the retriction
                // that all y values are positive means there is a unique 
                // solution.
                double integral = int_x - (totalInt - line.integralWrtX());
                
                // 'integral' now contains the contribution from just this line
                // segment.
                
                // Integral defined by polynomial in x:
                // A*x^2 + B*x + C = 0
                double A = line.m/2.0;
                double B = line.c;
                double C = -integral - line.m*line.min_x*line.min_x/2.0 - line.c*line.min_x;
                
                // Check that A is non-zero before using quadratic inversion
                // formula
                if(A==0){
                    return -C/B;
                }
                
                // Solutions for x...
                double disc = B*B - 4*A*C;
                
                double x_p = (-B + Math.sqrt(disc))/(2*A);
                
                // XXX: why do we use x_p rather than x_m?
                // double x_m = (-B - Math.sqrt(disc))/(2*A);
                
                return x_p;
                
            }
            
            // Otherwise, move on to next line segment.
            
        }
        
        // Shouldn't ever reach this code.
        System.err.println("PositiveLinear.drawX(): reached supposedly "
                + "un-reachable code beyond end of loop! int_x = "+int_x);
        
        return -10.0;
    }
}