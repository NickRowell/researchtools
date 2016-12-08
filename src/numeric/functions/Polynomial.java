package numeric.functions;

/**
 * Simple variable order polynomial implementation.
 *
 *
 * @author nrowell
 * @version $Id$
 */
public class Polynomial {
	
	/**
	 * Coefficients array, where the index equals the power associated with that term in the
	 * polynomial, e.g. for a cubic polynomial with 4 elements in a:
	 * 
	 * y(x) = a[0] + a[1]*x + a[2]*x*x + a[3]*x*x*x
	 * 
	 */
	double[] a;
	
	/**
	 * Order.
	 */
	int N;
	
	
	public Polynomial(double[] a)
	{
		N = a.length;
		this.a = new double[N];
		System.arraycopy(a, 0, this.a, 0, N);
	}
	
	public double getFn(double x)
	{
		double y = 0.0;
		
		double xx = 1.0;
		
		for(int i=0; i<N; i++)
		{
			y += a[i]*xx;
			xx *= xx;
		}
		
		return y;
	}

}
