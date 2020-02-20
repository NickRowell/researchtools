package numeric.functions;

import Jama.Matrix;

/**
 * Simple variable order polynomial implementation, including covariance information on the polynomial coefficients.
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
	public double[] a;
	
	/**
	 * Covariance matrix on the polynomial coefficients, used for error propagation.
	 */
	public double[][] cov;
	
	/**
	 * Order.
	 */
	int N;
	
	/**
	 * Constructor that takes the polynomial coefficients only. The covariance matrix is set to the identity.
	 * 
	 * @param a
	 * 	The polynomial coefficients [1, x, x^2, ...]
	 */
	public Polynomial(double[] a)
	{
		N = a.length;
		this.a = new double[N];
		System.arraycopy(a, 0, this.a, 0, N);
		cov = new double[N][N];
		for(int i=0; i<N; i++) {
			cov[i][i] = 1.0;
		}
	}
	
	/**
	 * Constructor that takes the polynomial coefficients and the associated covariance matrix.
	 * 
	 * @param a
	 * 	The polynomial coefficients [1, x, x^2, ...]
	 * @param cov
	 * 	The covariance matrix on the polynomial coefficients.
	 */
	public Polynomial(double[] a, double[][] cov)
	{
		this(a);
		for(int i=0; i<N; i++) {
			System.arraycopy(cov[i], 0, this.cov[i], 0, N);
		}
	}
	
	/**
	 * Returns a {@link String} expression suitable for plotting the polynomial using Gnuplot.
	 */
	public String toString() {
		StringBuilder polyFn = new StringBuilder();
		polyFn.append("f(x) = ");
		polyFn.append(a[0]);
		
		for(int i=1; i<N; i++)
		{
			polyFn.append((a[i] < 0.0 ? " - " : " + ") + Math.abs(a[i]) + " * x**" + i);
		}
		
		return polyFn.toString();
	}

	/**
	 * Compute the value of the polynomial and the variance on the value considering the covariance on the
	 * polynomial coefficients, at a given set of X coordinates.
	 *  
	 * @param x
	 * 	The X coordinates at which to compute the polynomial value and error.
	 * @return
	 * 	A two dimensional array containing the value of the polynomial [i][0] and the variance on the value [i][1]
	 * considering the covariance on the polynomial coefficients, for each X coordinate [i].
	 */
	public double[][] getFnWithError(double[] x)
	{
		// Number of values to interpolate
		int n = x.length;
		
		double[][] xArr = new double[n][N];
		for(int j=0; j<n; j++) {
			double xx = 1.0;
			for(int i=0; i<N; i++)
			{
				xArr[j][i] = xx;
				xx *= x[j];
			}
		}
		
		double[][] aArr = new double[N][1];
		for(int i=0; i<N; i++)
		{
			aArr[i][0] = a[i];
		}
		
		Matrix xMat = new Matrix(xArr);
		Matrix aMat = new Matrix(aArr);
		Matrix aCov = new Matrix(cov);
		
		// Compute the function
		Matrix yMat = xMat.times(aMat);
		
		// Propagate the covariance on coefficients to the solution
		Matrix yCov = xMat.times(aCov.times(xMat.transpose()));
		
		// Extract the solution and the variance
		double[][] yArr = new double[n][2];
		for(int i=0; i<n; i++)
		{
			yArr[i][0] = yMat.get(i, 0);
			yArr[i][1] = yCov.get(i, i);
		}
		
		return yArr;
	}
	
	/**
	 * Compute the value of the polynomial and the variance on the value considering the covariance on the
	 * polynomial coefficients.
	 *  
	 * @param x
	 * 	The X coordinate at which to compute the polynomial value.
	 * @return
	 * 	A two element array containing the value of the polynomial [0] and the variance on the value [1]
	 * considering the covariance on the polynomial coefficients.
	 */
	public double[] getFnWithError(double x)
	{
		return getFnWithError(new double[]{x})[0];
	}
	
	/**
	 * Compute the value of the polynomial at a set of X coordinates.
	 * 
	 * @param x
	 * 	The X coordinates at which to compute the polynomial value.
	 * @return
	 * 	An array containing the polynomial value at each X coordinate.
	 */
	public double[] getFn(double[] x)
	{
		// Number of values to interpolate
		int n = x.length;
		
		double[][] xArr = new double[n][N];
		for(int j=0; j<n; j++) {
			double xx = 1.0;
			for(int i=0; i<N; i++)
			{
				xArr[j][i] = xx;
				xx *= x[j];
			}
		}
		
		double[][] aArr = new double[N][1];
		for(int i=0; i<N; i++)
		{
			aArr[i][0] = a[i];
		}
		
		Matrix xMat = new Matrix(xArr);
		Matrix aMat = new Matrix(aArr);
		
		Matrix yMat = xMat.times(aMat);
		
		double[] yArr = new double[n];
		for(int i=0; i<n; i++)
		{
			yArr[i] = yMat.get(i, 0);
		}
		
		return yArr;
	}
	
	/**
	 * Compute the value of the polynomial at a given X coordinate.
	 * 
	 * @param x
	 * 	The X coordinate at which to compute the polynomial value.
	 * @return
	 * 	The polynomial value at the X coordinate.
	 */
	public double getFn(double x)
	{
		return getFn(new double[]{x})[0];
	}
	
    /**
     * Get the derivative of the polynomial.
     *
     * @param x
     *            x value at which to evaluate function
     * @return derivative value
     */
    public double getDerivative(double x) {

        double value = 0.0;
        double param = 1;
        for (int i = 1; i < N; i++) {

            value += i * a[i] * param;

            param *= x;
        }

        return value;
    }

}
