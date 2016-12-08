package numeric.minimisation.llsq;

import Jama.Matrix;



public class HouseHolderLeastSquares
{
	
	/**
	 * The default value that unconstrained parameters are set to, so that a
	 * solution can always be obtained even for underdetermined equation sets.
	 */
	private static double DEFAULT_UNCONSTRAINED_PARAMETER = 1.0;
	
	// Solution vector
	public Matrix x;
	
	// Covariance matrix of solution
	public Matrix C;
	
	// Triangularized form of input information matrix (useful for computing running solutions)
	public Matrix U;
	
	// Upper part of transformed observation vector (useful for computing running solutions)
	public Matrix z1;
	
	// Sum-of-square residuals
	public double r2;
	
	// Number of parameters
	int nPar;
	
	// Number of observations - can be used to determine how many parameters are constrained, or how
	// many rows of the U and z1 matrices have been filled in with default values in order to obtain
	// a solution in the case of underdetermined systems.
	int nObs;
	
	
	
	/**
	 * Solves for x in A*x = b using Householder least squares. Note that for underdetermined
	 * systems (where A has more columns than rows), unconstrained parameters are set to a
	 * default value of {@link HouseHolderLeastSquares#DEFAULT_UNCONSTRAINED_PARAMETER} so
	 * that a solution for x is always returned.
	 * 
	 * @param pA	Design matrix
	 * @param pb	Observation vector
	 */
	public HouseHolderLeastSquares(Matrix pA, Matrix pb)
	{
		// Take copies of matrices
		Matrix QA = pA.copy();
		Matrix Qb = pb.copy();
		
		// Number of rows in A
		int M = QA.getRowDimension();
		// Number of columns in A
		int N = QA.getColumnDimension();
		
		// Do some sanity checks on Matrix sizes
		if(Qb.getRowDimension() != M)
			throw new IllegalArgumentException("Equation set badly formed! A has "+M+" rows; b has "+Qb.getRowDimension());
		
		if(Qb.getColumnDimension() != 1)
			throw new IllegalArgumentException("Equation set badly formed! b has "+Qb.getColumnDimension() + " columns (!=1)");
		
		// Set number of free parameters
		nPar = N;
		
		// Set number of observations
		nObs = M;
		
		// Loop over elementary reflectors - one for each column of A (or row if there are fewer rows, i.e. underconstrained system)
		for(int c=0; c<Math.min(M,N); c++)
		{
			// Set up Householder vector
			Matrix v = new Matrix(M,1);
			
			// The upper c elements are zero.
			// The lower M-c elements are extracted from rows c->M-1, column c
			v.setMatrix(c, M-1, 0, 0, QA.getMatrix(c, M-1, c, c));
			
			// Basis vector for this column
			Matrix ec = new Matrix(M,1);
			ec.set(c, 0, 1.0);
			
			// Get sign of c-th element. This lies on the diagonal of A.
			double sign = Math.signum(QA.get(c,c));
			
			// If diagonal is zero, set the sign arbitrarily to +1
			if(sign==0)
				sign = 1.0;
			
			// Scalar alpha
			double alpha = -sign * Math.sqrt(v.transpose().times(v).get(0,0));

			// Complete Householder vector
			v = v.minus(ec.times(alpha));
			
			// Compute Householder orthogonal matrix for this column
			double beta = 2/(v.transpose().times(v).get(0, 0));
			Matrix H = Matrix.identity(M, M).minus(v.times(v.transpose()).times(beta));
			
			// Apply reflection
			QA = H.times(QA);
			Qb = H.times(Qb);
		}
		
		// Matrix QA is now in upper triangular form.
		// Vector Qb is partitioned into [c1 c2]^T
		
		U = new Matrix(N,N);
		z1 = new Matrix(N,1);
		
		if(M<N) {
			
			// More columns than rows: equation set is underconstrained, there is no unique solution for the 
			// parameters. We can still find a solution however, by setting all the parameters with indices
			// of M=<i<N to a default value. This is achieved by padding out R and b with defaults; we extend
			// R by adding 1s along the main diagonal, and b by adding default values.
			
			U.setMatrix(0, M-1, 0, N-1, QA.getMatrix(0, M-1, 0, N-1));
			z1.setMatrix(0, M-1, 0, 0, Qb.getMatrix(0, M-1, 0, 0));
			
			for(int i=M; i<N; i++) {
				// Pad out R with identity matrix
				U.set(i, i, 1.0);
				// Set default value for parameter i
				z1.set(i, 0, DEFAULT_UNCONSTRAINED_PARAMETER);
			}
			
			// No solution for z2 exists. Set sum-of-square residuals to zero
			r2 = 0.0;
		}
		else {
			
			// Extract upper triangular portion of matrix A
			U.setMatrix(0, N-1, 0, N-1, QA.getMatrix(0, N-1, 0, N-1));
			
			// Extract upper N and lower M-N portions of 
			z1.setMatrix(0, N-1, 0, 0, Qb.getMatrix(0, N-1, 0, 0));
			
			// Compute sum-of-square residual from partitioned z vector
			Matrix z2 = Qb.getMatrix(N, M-1, 0, 0);
			
			// Get sum-of-square residuals from z2^Tz2
			r2 = z2.transpose().times(z2).get(0, 0);
		}
		
		// Solve for x by back-substitution
		x = solveUpperTriangularSystem(U, z1);
		
		// Get parameter covariance matrix from (U^TU)^{-1} 
		C = U.transpose().times(U).inverse();
		
	}
	
	/**
	 * Solves the upper triangular system Rx=b using back-substitution. All
	 * elements of R below the main diagonal must be zero. R must also be square;
	 * note that if the system is not square, then it can be padded out with defaults
	 * in order to achieve a square system and get a solution (although it won't be
	 * unique).
	 * 
	 * @param R
	 * 	The (square) upper triangular matrix
	 * @param b
	 * 	The column vector of observations
	 * @return
	 * 	The parameter solution x
	 */
	private static Matrix solveUpperTriangularSystem(Matrix R, Matrix b)
	{
		// Number of rows in R
		int M = R.getRowDimension();
		
		// Number of columns in R
		int N = R.getColumnDimension();
		
		// Sanity checks on inputs
		if(M!=N)
			throw new IllegalArgumentException("Equation set badly formed! R is not square: M,N = "+M+","+N);
		
		if(b.getRowDimension() != M)
			throw new IllegalArgumentException("Equation set badly formed! R has "+M+" rows; b has "+b.getRowDimension());
		
		if(b.getColumnDimension() != 1)
			throw new IllegalArgumentException("Equation set badly formed! b has "+b.getColumnDimension() + " columns (!=1)");
		
		// Check R is upper triangular (all zero below diagonal)
		for(int r=1; r<M; r++)
			for(int c=0; c<r; c++)
				if(Math.abs(R.get(r, c))>1e-9)
					throw new IllegalArgumentException("R is not upper triangular! Found non-zero element "+R.get(r, c)+" at ("+r+","+c+")");
		
		// Generate empty solution vector x
		Matrix x = new Matrix(N,1);
		
		// Compute first element in back substitution
		x.set(N-1, 0, b.get(N-1, 0)/R.get(N-1, N-1));
		
		// Now apply general rule for remaining elements
		for(int r=N-2; r>=0; r--)
		{
			double sum = 0.0;
			for(int c=r+1; c<N; c++)
			{
				sum += R.get(r, c) * x.get(c, 0);
			}
			
			x.set(r, 0, (b.get(r, 0)-sum)/R.get(r, r) );
		}
		
		return x;
	}
	

	/**
	 * Test method.
	 * @param args
	 */
	public static void main(String[] args)
	{
		
		// Build dummy linear equation set
		
		// Choose some parameters. Can add extra here.
		Matrix x = new Matrix(new double[][]{{2},{3},{4},{5}});
		
		// Specify number of rows in design matrix
		int M = 10;
		
		// Create design matrix and fill with random values
		Matrix A = new Matrix(M,x.getRowDimension());
		for(int r=0; r<A.getRowDimension(); r++)
			for(int c=0; c<A.getColumnDimension(); c++)
				A.set(r, c, Math.random());
		
		// Get observation matrix. No errors added (yet).
		Matrix b = A.times(x);
		
		// Add random error to observations
		for(int r=0; r<M; r++)
			b.set(r, 0, b.get(r, 0) + 0.1*(Math.random() - 0.5));
		
		// Compute solution using householder least squares
		HouseHolderLeastSquares hhlsq = new HouseHolderLeastSquares(A, b);
		
		// Compute solution using normal equations (via JAMA)
		Matrix x_norms = A.solve(b);
		double r2_norms = (A.times(x_norms).minus(b)).normF();
		r2_norms *= r2_norms;
		Matrix C_norms = (A.transpose().times(A)).inverse();
		
		
		System.out.println("Solution via Householder Least Squares = ");
		hhlsq.x.print(5, 5);
		
		System.out.println("Solution via Normal equations = ");
		x_norms.print(5, 5);
		
		System.out.println("Sum-of-square residuals:");
		System.out.println("Householder least squares = "+hhlsq.r2);
		System.out.println("Normal equations method   = "+r2_norms);
		
		System.out.println("Parameter covariance:");
		System.out.println("Householder least squares:");
		hhlsq.C.print(5, 5);
		System.out.println("Normal equations method:");
		C_norms.print(5, 5);
		
		
	}
	
	
}
