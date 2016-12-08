package numeric.functions;

/**
 * Interface represents one-dimensional functions that are monotonic, and
 * thus can be uniquely inverted.
 *
 * @author nrowell
 * @version $Id$
 */
public interface MonotonicFunction1D {
	
	/**
	 * Get the forward value of the function, i.e. y = f(x), and the first
	 * derivative.
	 * @param x
	 * 	The value at which to evaulate the function.
	 * @return
	 * 	f(x) and df(x)/dy
	 */
	public double[] getY(double x);
	
	/**
	 * Get the inverse value of the function, i.e. x = f^{-1}(y),  and the
	 * first derivative. Note that because the function is monotonic, the
	 * inverse is a unique value.
	 * @param y
	 * 	The value at which to evaulate the inverse function.
	 * @return
	 * 	f^{-1}(y) and df^{-1}(x)/dy
	 */
	public double[] getX(double y);
	
}
