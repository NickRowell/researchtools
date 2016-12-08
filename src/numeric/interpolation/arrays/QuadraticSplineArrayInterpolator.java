package numeric.interpolation.arrays;


import Jama.Matrix;
import numeric.functions.Polynomial;

/**
 * Spline interpolation using second order polynomials to interpolate within each element.
 * The value and first derivative of the polynomials are continuous across neighbouring
 * segments.
 *
 * @author nrowell
 * @version $Id$
 */
public class QuadraticSplineArrayInterpolator extends ArrayInterpolator {

	/**
	 * Array of interpolating polynomials.
	 */
	Polynomial[] polys;
	
	
	public QuadraticSplineArrayInterpolator(double[] y)
	{
		N = y.length;
		polys = new Polynomial[N];
		
		// Build linear equation set a*x = b
		double[][] a = new double[3*N][3*N];
		double[][] b = new double[3*N][1];
		
		
		// x has size [3*N][1]
		
		// Index of next free row in design and observation arrays
		int row = 0;
		
		// Enter integral constraints into arrays
		for(int i=0; i<N; i++)
		{
			a[row][3*i+0] = 1.0/3.0;
			a[row][3*i+1] = 1.0/2.0;
			a[row][3*i+2] = 1.0;
				
			b[row][0] = y[i];
			row++;
		}
		
		// First derivative continuous across boundaries
		for(int i=0; i<N-1; i++)
		{
			a[row][3*i+0] = 2.0;
			a[row][3*i+1] = 1.0;
			a[row][3*i+2] = 0.0;
			
			a[row][3*i+3] =  0.0;
			a[row][3*i+4] = -1.0;
			a[row][3*i+5] =  0.0;
			
			b[row][0] = 0.0;
			row++;
		}
		
		// Function itself continuous across boundaries
		for(int i=0; i<N-1; i++)
		{
			a[row][3*i+0] = 1.0;
			a[row][3*i+1] = 1.0;
			a[row][3*i+2] = 1.0;
			
			a[row][3*i+3] =  0.0;
			a[row][3*i+4] =  0.0;
			a[row][3*i+5] = -1.0;
			
			b[row][0] = 0.0;
			row++;
		}
		
		// First derivative zero at the start (x[0])
		a[row][0] = 0.0;
		a[row][1] = 1.0;
		a[row][2] = 0.0;
			
		b[row][0] = 0.0;
		row++;

		// First derivative zero at the end (x[N])
		a[row][3*(N-1)+0] = 2.0;
		a[row][3*(N-1)+1] = 1.0;
		a[row][3*(N-1)+2] = 0.0;
		
		b[row][0] = 0.0;
		row++;
		
		// Solve equation set
		Matrix A = new Matrix(a);
		Matrix B = new Matrix(b);
		Matrix X = A.solve(B);
		
		// Extract polynomial coefficients
		for(int i=0; i<N; i++)
		{
			double aa = X.get(3*i+0, 0);
			double bb = X.get(3*i+1, 0);
			double cc = X.get(3*i+2, 0);
			polys[i] = new Polynomial(new double[]{cc,bb,aa,0});
		}
		
	}
	
	public double getValue(double x)
	{
		int poly = (int)Math.floor(x);
		
		double y = polys[poly].getFn(x - Math.floor(x));
		
		return y;
		
	}
}
