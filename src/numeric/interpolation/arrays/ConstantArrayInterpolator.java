package numeric.interpolation.arrays;

/**
 * The simplest type of array interpolation: values are constant across the width of each
 * element.
 *
 *
 * @author nrowell
 * @version $Id$
 */
public class ConstantArrayInterpolator extends ArrayInterpolator {

	/**
	 * Internal copy of the original array.
	 */
	double[] y;
	
	public ConstantArrayInterpolator(double[] y)
	{
		N = y.length;
		this.y = new double[N];
		System.arraycopy(y, 0, this.y, 0, N);
	}

	public double getValue(double x)
	{
		return y[(int)Math.floor(x)];
	}
	
	
}
