package numeric.interpolation.arrays;


import util.ArrayUtil;
import Jama.Matrix;
import numeric.functions.Polynomial;

/**
 * Spline interpolation using third order polynomials to interpolate within each element.
 * The value, first and second derivatives of the polynomials are continuous across neighbouring
 * segments.
 * 
 * @author nrowell
 * @version $Id$
 */
public class CubicSplineArrayInterpolator extends ArrayInterpolator {

	/**
	 * Option to improve numerical conditioning by shifting each polynomial
	 * segment to [0:1] range.
	 */
	boolean shifted = true;
	
	/**
	 * Option to improve numerical conditioning by scaling array elements to 
	 * range [0:1]
	 */
	boolean scaled = true;
	
	/**
	 * Scale factor applied to array elements
	 */
	private double scaleFactor = 1.0;
	
	/**
	 * Array of interpolating polynomials.
	 */
	Polynomial[] polys;
	
	
	public CubicSplineArrayInterpolator(double[] y)
	{
		N = y.length;
		polys = new Polynomial[N];
		
		// Apply optional scaling to the input
		if(scaled)
		{
			scaleFactor = 1.0/ArrayUtil.max(y);
		}
		else
		{
			scaleFactor = 1.0;
		}
		
		// x[n] and x[n+1] form the boundaries of segment n, with sample y[n].
		// This is used if the shifting option is turned off.
		double[] x = new double[N+1];
		for(int i=0; i<x.length; i++)
			x[i] = i;
		
		// Build linear equation set a*x = b
		double[][] a = new double[4*N][4*N];
		double[][] b = new double[4*N][1];
		
		
		// x has size [4*N][1]
		
		// Index of next free row in design and observation arrays
		int row = 0;
		
		// Enter integral constraints into arrays
		for(int i=0; i<N; i++)
		{
			if(shifted)
			{
				a[row][4*i+0] = (1.0)/4.0;
				a[row][4*i+1] = (1.0)/3.0;
				a[row][4*i+2] = (1.0)/2.0;
				a[row][4*i+3] = (1.0);
			}
			else
			{
				a[row][4*i+0] = (x[i+1]*x[i+1]*x[i+1]*x[i+1] - x[i]*x[i]*x[i]*x[i])/4.0;
				a[row][4*i+1] = (x[i+1]*x[i+1]*x[i+1] - x[i]*x[i]*x[i])/3.0;
				a[row][4*i+2] = (x[i+1]*x[i+1] - x[i]*x[i])/2.0;
				a[row][4*i+3] = (x[i+1] - x[i]);
			}
			b[row][0] = y[i]*scaleFactor;
			row++;
		}
		
		// Second derivative continuous across boundaries
		for(int i=0; i<N-1; i++)
		{
			if(shifted)
			{
				a[row][4*i+0] = 6.0*1.0;
				a[row][4*i+1] = 2.0;
				a[row][4*i+2] = 0.0;
				a[row][4*i+3] = 0.0;
				a[row][4*i+4] = -6.0*0.0;
				a[row][4*i+5] = -2.0;
				a[row][4*i+6] = 0.0;
				a[row][4*i+7] = 0.0;
			}
			else
			{
				a[row][4*i+0] = 6.0*x[i+1];
				a[row][4*i+1] = 2.0;
				a[row][4*i+2] = 0.0;
				a[row][4*i+3] = 0.0;
				a[row][4*i+4] = -6.0*x[i+1];
				a[row][4*i+5] = -2.0;
				a[row][4*i+6] = 0.0;
				a[row][4*i+7] = 0.0;
			}
			b[row][0] = 0.0;
			row++;
		}
		
		// First derivative continuous across boundaries
		for(int i=0; i<N-1; i++)
		{
			if(shifted)
			{
				a[row][4*i+0] = 3.0*1.0*1.0;
				a[row][4*i+1] = 2.0*1.0;
				a[row][4*i+2] = 1.0;
				a[row][4*i+3] = 0.0;
				a[row][4*i+4] = -3.0*0.0*0.0;
				a[row][4*i+5] = -2.0*0.0;
				a[row][4*i+6] = -1.0;
				a[row][4*i+7] = 0.0;
			}
			else
			{
				a[row][4*i+0] = 3.0*x[i+1]*x[i+1];
				a[row][4*i+1] = 2.0*x[i+1];
				a[row][4*i+2] = 1.0;
				a[row][4*i+3] = 0.0;
				a[row][4*i+4] = -3.0*x[i+1]*x[i+1];
				a[row][4*i+5] = -2.0*x[i+1];
				a[row][4*i+6] = -1.0;
				a[row][4*i+7] = 0.0;
			}
			b[row][0] = 0.0;
			row++;
		}
		
		// Function itself continuous across boundaries
		for(int i=0; i<N-1; i++)
		{
			if(shifted)
			{
				a[row][4*i+0] = 1.0*1.0*1.0;
				a[row][4*i+1] = 1.0*1.0;
				a[row][4*i+2] = 1.0;
				a[row][4*i+3] = 1.0;
				a[row][4*i+4] = -0.0*0.0*0.0;
				a[row][4*i+5] = -0.0*0.0;
				a[row][4*i+6] = -0.0;
				a[row][4*i+7] = -1.0;
			}
			else
			{
				a[row][4*i+0] = x[i+1]*x[i+1]*x[i+1];
				a[row][4*i+1] = x[i+1]*x[i+1];
				a[row][4*i+2] = x[i+1];
				a[row][4*i+3] = 1.0;
				a[row][4*i+4] = -x[i+1]*x[i+1]*x[i+1];
				a[row][4*i+5] = -x[i+1]*x[i+1];
				a[row][4*i+6] = -x[i+1];
				a[row][4*i+7] = -1.0;
			}
			b[row][0] = 0.0;
			row++;
		}
		
		// Second derivative zero at the start (x[0])
		if(shifted)
		{
			a[row][0] = 6.0*0.0;
			a[row][1] = 2.0;
		}
		else
		{
			a[row][0] = 6.0*x[0];
			a[row][1] = 2.0;
		}
		b[row][0] = 0.0;
		row++;
		
		// Second derivative zero at the end (x[N])
//		if(shifted)
//		{
//			a[row][4*(N-1)+0] = 6.0*1.0;
//			a[row][4*(N-1)+1] = 2.0;
//		}
//		else
//		{
//			a[row][4*(N-1)+0] = 6.0*x[N];
//			a[row][4*(N-1)+1] = 2.0;
//		}
//		b[row][0] = 0.0;
//		row++;
		
		// First derivative zero at the start (x[0])
		if(shifted)
		{
			a[row][0] = 3.0*0.0*0.0;
			a[row][1] = 2.0*0.0;
			a[row][2] = 1.0;
		}
		else
		{
			a[row][0] = 3.0*x[0]*x[0];
			a[row][1] = 2.0*x[0];
			a[row][2] = 1.0;
		}
		b[row][0] = 0.0;
		row++;

		// First derivative zero at the end (x[N])
		if(shifted)
		{
			a[row][4*(N-1)+0] = 3.0*1.0*1.0;
			a[row][4*(N-1)+1] = 2.0*1.0;
			a[row][4*(N-1)+2] = 1.0;
		}
		else
		{
			a[row][4*(N-1)+0] = 3.0*x[N]*x[N];
			a[row][4*(N-1)+1] = 2.0*x[N];
			a[row][4*(N-1)+2] = 1.0;
		}
		b[row][0] = 0.0;
		row++;
		
		// Solve equation set
		Matrix A = new Matrix(a);
		Matrix B = new Matrix(b);
		Matrix X = A.solve(B);
		
//		A.print(5, 5);
//		B.print(5, 5);
//		X.print(5, 5);
		
		// Solution via Householder LSQ
//		X = new HouseHolderLeastSquares(A,B).x;
		
		// Extract polynomial coefficients
		for(int i=0; i<N; i++)
		{
			// Polynomial coefficients are arranged in a*x3 + b*x2 + c*x + d
			double aa = X.get(4*i+0, 0);
			double bb = X.get(4*i+1, 0);
			double cc = X.get(4*i+2, 0);
			double dd = X.get(4*i+3, 0);
			polys[i] = new Polynomial(new double[]{dd,cc,bb,aa});
		}
		
	}
	
	public double getValue(double x)
	{
		int poly = (int)Math.floor(x);
		
		double y = Double.NaN;
		
		if(shifted)
		{
			double x_shifted = x - Math.floor(x);
			y = polys[poly].getFn(x_shifted);
		}
		else
		{
			y = polys[poly].getFn(x);
		}
		
		if(scaled)
			y /= scaleFactor;
		
		return y;
		
	}
}
