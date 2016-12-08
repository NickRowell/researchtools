package numeric.integration;

/**
 * Utility methods to compute the definite integral of functions.
 *
 *
 * @author nrowell
 * @version $Id$
 */
public class IntegrationUtils {
	
	/**
	 * Computes the convolution of the two functions using numerical integration, i.e. computed the
	 * numerical integration of a(x)*b(x0-x) over all x.
	 * 
	 * @param a
	 * 	The function a(x)
	 * @param b
	 * 	The function b(x)
	 * @param xmin
	 * 	The lower limit on the integration range
	 * @param xmax
	 * 	The upper limit on the integration range
	 * @param x0
	 * 	The point at which to evaluate the convolution
	 * @param h
	 * 	The step size to use in the numerical integration
	 * @return
	 * 	The value of the convolution of a(x) and b(x) at the point x0, i.e. the integral of
	 * a(x)*b(x0-x) over all x.
	 */
	public static double convolve(final IntegrableFunction a, final IntegrableFunction b, double xmin, double xmax, final double x0, double h) {
		
		// Create an IntegrableFunction to represent the convolution equation
		IntegrableFunction conv = new IntegrableFunction() {
			@Override
			public double evaluate(double x) {
				return a.evaluate(x) * b.evaluate(x0 - x);
			}
		};
		
		return integrate(conv, xmin, xmax, h);
	}
	
	/**
	 * Computes the definite integral of the function using Romberg integration
	 * (trapezium approximation with Richardson's Extrapolation).
	 * 
	 * Implements integration algorithm.
	 * @param y
	 * 	The {@link IntegrableFunction} to integrate.
	 * @param a
	 * 	Lower limit on range.
	 * @param b
	 * 	Upper limit on range.
	 * @param h
	 * 	Step size for trapezium method.
	 * @return
	 * 	Numerical integral of function over range [a:b]
	 */
	public static double integrate(IntegrableFunction y, double a, double b, double h)
	{
		
		// Boundaries on integration steps
		double stepStart = a;
		double stepEnd, stepMid;
		
		// Integrals using step sizes h and h/2
		double T_h  = 0.0;
		double T_h2 = 0.0;
		
		while(stepStart < b)
		{
			stepEnd = stepStart+h;
			// On final step, trim the range down to avoid overshoot
			if(stepEnd > b)
			{
				stepEnd = b;
			}
			// Get midpoint
			stepMid = (stepStart+stepEnd)/2.0;
			
			// Recompute step size (to account for trimmed final step)
			double step = stepEnd - stepStart;
			
			// Values of the function at the start, midpoint and end of this step
			double yStart = y.evaluate(stepStart);
			double yMid   = y.evaluate(stepMid);
			double yEnd   = y.evaluate(stepEnd);
			
			// Integral using full step
			T_h += 0.5 * step * (yStart + yEnd);
			
			// Integral using half step
			T_h2 += 0.5 * (step/2.0) * (yStart + yMid);
			T_h2 += 0.5 * (step/2.0) * (yMid + yEnd);
			
			stepStart += h;
		}
		
		// Apply Richardson's Extrapolation to get solution
		return (4.0*T_h2 - T_h)/3.0;
	}
	
}