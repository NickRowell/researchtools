package numeric.minimisation.nllsq.algo;

import Jama.Matrix;

/**
 * Implementation of the Levenberg-Marquardt nonlinear least squares algorithm.
 * 
 * In order to use this class you must extend it and implement the getModel(Matrix params)
 * method and optionally the getJacobian(double[] params) method. If you would like the
 * fitter to use a numerical approximation to the Jacobian instead then override the
 * {@link #useFiniteDifferencesJacobian()} to return true and also override {@link #finiteDifferencesStepSizePerParam()}
 * to provide an appropriate step size for each parameter.
 * 
 * Normal usage:
 * 
 *	LevenbergMarquardt lma = new LevenbergMarquardt() {
 *
 *		@Override
 *		public double[] getModel(double[] params) {
 *
 *			double[] model = new double[N];
 *
 *			// Implement the parameters -> model function
 *
 *			return model;
 *		}
 *
 *		@Override
 * 		public double[][] getJacobian(double[] params) {
 *	
 * 			double[][] jac = new double[N][M];
 *	
 *			// Compute the Jacobian elements
 *	
 *			return jac;
 *		}
 *	}
 *
 *	lma.setData(electronSamples);
 *	lma.setCovariance(cov);
 *	lma.setInitialGuessParameters(params);
 *	// Perform the optimization
 *	lma.fit(500, true);
 *		        
 *	// Extract the solution
 *	double[] solution = lma.getParametersSolution();
 * 
 * @author nrowell
 * @version $Id$
 */
public abstract class LevenbergMarquardt {

	/**
	 * Status flags to record final state of algorithm.
	 */
	public enum STATUS {FAILED_SINGULAR_MATRIX, FAILED_DAMPING_THRESHOLD_EXCEEDED, FAILED_ITERATION_LIMIT_EXCEEDED, SUCCESS};
	
    /** 
     * Absolute step size used in finite difference Jacobian approximation,
     * for Jacobian of parameter solution with respect to data. This step is
     * applied to the datapoints not the parameters.
     */
    private double h = 1E-2;   

    /** 
     * Exit tolerance 
     */
    private double exitTolerance = 1E-32;

    /**
     * Max damping scale factor. This is multiplied by automatically selected
     * starting value of damping parameter, which is 10^{-3}
     * times the average of the diagonal elements of J^T*J
     */
    private double maxDamping = 1E32;
    
    /**
     * Nx1 column vector of observed values
     * 
     * Y = [y_0, y_1, y_2, ..., y_{N-1}]^T
     * 
     */
    private Matrix data;
    
    /**
     * Number of data points.
     */
    protected int N;
    
    /**
     * NxN covariance matrix for the observed values
     */
    private Matrix dataCovariance;
    
    /**
     * Mx1 column vector of parameters
     * 
     * P = [p_0, p_1, p_2, ..., p_{M-1}]^T
     * 
     */
    private Matrix params;
    
    /**
     * Number of free parameters.
     */
    protected int M;

    /**
     * Set the absolute step size used in the finite difference Jacobian estimation for the
     * calculation of the rate of change of parameters as a function of the data.
     * 
     * @param h
     * 	The finite step size (applied to the data, not the parameters).
     */
    public void setH(double h)  { 
    	this.h = h;
    }
    
    /**
     * Set the exit tolerance - if the (absolute value of the) relative change in the chi-square
     * from one iteration to the next is lower than this, then we're at the minimum and the fit
     * is halted.
     * @param exitTolerance
     * 	The exit tolerance to set
     */
    public void setExitTolerance(double exitTolerance) { 
    	this.exitTolerance = exitTolerance;
    }
    
    /**
     * Set the maximum damping factor. If the damping factor becomes larger than this during the fit,
     * then we're stuck and cannot reach a better solution.
     * 
     * @param maxDamping
     * 	The max damping factor to set
     */
    public void setMaxDamping(double maxDamping) {
    	this.maxDamping = maxDamping;
    }

    /**
     * Set the Nx1 column vector of observed values
     * @param data
     * 	The Nx1 column vector of observed values
     */
    public void setData(double[] data) {
    	N = data.length;
    	this.data = new Matrix(N, 1);
    	for(int r=0; r<N; r++) {
    		this.data.set(r, 0, data[r]);
    	}
    }
    
    /**
     * Set the Mx1 column vector of initial-guess parameters.
     * 
     * @param params
     * 	The Mx1 column vector of initial-guess parameters.
     */
    public void setInitialGuessParameters(double[] params) {
    	M = params.length;
    	this.params = new Matrix(M, 1);
    	for(int r=0; r<M; r++) {
    		this.params.set(r, 0, params[r]);
    	}
    }
    
    /**
     * Set the NxN matrix of covariance of the data points. The format of the input
     * array is such that element [i][j] lies at row i column j in the covariance
     * matrix.
     * 
     * @param dataCovariance
     * 	NxN matrix of covariance of the data points
     */
    public void setCovariance(double[][] dataCovariance) {
    	
    	// Check the input is a square matrix
    	int rows = dataCovariance.length;
    	int cols = dataCovariance[0].length;
    	for(int i=0; i<rows; i++) {
    		if(dataCovariance[i].length!=rows) {
    			throw new RuntimeException("The data covariance matrix is not square!");
    		}
    	}
    	
    	this.dataCovariance = new Matrix(rows, cols);
    	for(int i=0; i<rows; i++) {
    		for(int j=0; j<cols; j++) {
    			this.dataCovariance.set(i, j, dataCovariance[i][j]);
    		}
    	}
    }
    
    /**
     * Get the parameter values
     * @return
     * 	The internal values of the parameters following fitting.
     */
    public double[] getParametersSolution() {
    	double[] params = new double[M];
    	for(int r=0; r<M; r++) {
    		params[r] = this.params.get(r, 0);
    	}
    	return params;
    }
    
    /** 
     * Get the vector of observed values.
     * 
     * @return
     * 	The vector of observed values
     */
    private Matrix getData() {
    	return data;
    }
    
    /** 
     * Get f(X,P): column vector of model values given x points and current 
     * parameters set.
     * 
     * f(X,P) = [f(x_1, P), f(x_2, P), ... , f(x_N, P)]^T
     * 
     */
    public abstract double[] getModel(double[] params);
    
    private Matrix getModel() {
    	
    	double[] params = new double[M];
    	for(int i=0; i<M; i++) {
    		params[i] = this.params.get(i, 0);
    	}
    	
    	double[] modelArray = getModel(params);
    	
    	Matrix model = new Matrix(N,1);
    	
    	for(int i=0; i<N; i++) {
    		model.set(i, 0, modelArray[i]);
    	}
    	
    	return model;
    }
    
    /**
     * Get the Jacobian matrix -> the matrix of partial derivatives of the 
     * model values with respect to the parameters, given the current parameter set.
     * 
     * J = [ df(x_0, P)/dp_0  df(x_0, P)/dp_1  df(x_0, p)/dp_2  ...  df(x_0, P)/dp_{M-1} ]
     *     [ df(x_1, P)/dp_0  df(x_1, P)/dp_1  df(x_1, p)/dp_2  ...  df(x_1, P)/dp_{M-1} ]
     * 
     * 
     * Leading dimension is number of rows, trailing dimension is number of columns, i.e.
     * A[r][c] has r rows and c columns.
     * 
     * @param params
     * @return
     */
    public abstract double[][] getJacobian(double[] params);
    
    /**
     * Implementing classes should override this method to specify whether the fitter should use
     * the {@link #getJacobian(double[])} method to obtain the Jacobian, or use a finite differences
     * method in cases where no analytic Jacobian is possible.
     * 
     * @return
     * 	Boolean specifying whether to use the user-provided Jacobian (false) or an internal finite
     * differences method (true).
     */
    public boolean useFiniteDifferencesJacobian() {
    	return false;
    }
    
    /**
     * Implementing classes should override this to provide appropriate step sizes per parameter for use
     * in the finite differences Jacobian approximation, if they intend to use that.
     * @return
     * 	Appropriate step size per parameter for use in the finite difference Jacobian approximation.
     */
    public double[] finiteDifferencesStepSizePerParam() {
    	throw new RuntimeException("If the finite differences Jacobian approximation is to be used then "
    			+ "the finiteDifferencesStepSizePerParam() must be overridden also!");
    }
    
    /**
     * Callback function that user can override to perform actions when each good step occurs.
     */
    public void goodStepCallback() {
    	// Default action is to do nothing
    }
    
    /**
     * Gets the Jacobian for the current parameters values.
     * 
     * @return
     */
    private Matrix getJacobian() {
    	
    	if(useFiniteDifferencesJacobian()) {
    		return getJacobianByFiniteDifferences();
    	}
    	
    	double[] params = new double[M];
    	for(int i=0; i<M; i++) {
    		params[i] = this.params.get(i, 0);
    	}
    	return new Matrix(getJacobian(params));
    }
    
    /** 
     * First order central difference Jacobian approximation. This is the 
     * derivative of the model values with respect to the parameters.
     * 
     * Jacobian is N rows by M columns.
     * 
     * It's elements are the first derivatives of the model with respect to the
     * parameters, i.e.:
     * 
     * J = [ df(x_0, P)/dp_0  df(x_0, P)/dp_1  df(x_0, p)/dp_2  ...  df(x_0, P)/dp_{M-1} ]
     *     [ df(x_1, P)/dp_0  df(x_1, P)/dp_1  df(x_1, p)/dp_2  ...  df(x_1, P)/dp_{M-1} ]
     * 
     */
    public Matrix getJacobianByFiniteDifferences() {
        
        // Set up Jacobian matrix: N rows by M columns
        Matrix jac = new Matrix(N, M);

        // Matrix used to adjust parameters: M rows X 1 column
        Matrix delta = new Matrix(M, 1);
        
        // Get finite step sizes to use for each parameter
        double[] steps = finiteDifferencesStepSizePerParam();

        // Iterate over each parameter...
        for (int j=0; j<M; j++) {
        	
            // Make 'nudge' vector for this parameter. Set all entries to zero
            // except the parameter that is to be adjusted
            for (int k=0; k<M; k++) {
            	delta.set(k, 0, (k==j) ? steps[j] : 0.0);
            }
              
            // First advance this single parameter
            params.plusEquals(delta);
            
            // Get model values for advanced parameter set: f(x+h)
            Matrix model = getModel();

            // Build Jacobian by finite difference
            for (int i = 0; i < N; i++) {
            	jac.set(i, j, model.get(i,0));
            }

            // Now retard the parameter...
            params.minusEquals(delta.times(2.0));

            // Get model values for retarded parameter set: f(x-h)
            model = getModel();

            // Build Jacobian by finite difference 
            for (int i = 0; i < N; i++) {
            	jac.set(i, j,  (jac.get(i, j) - model.get(i,0)) / (2.0*steps[j]));
            }
            
            // Return the params to original value
            params.plusEquals(delta);
        }
        
        return jac;
    }
    
    /**
     * Perform Levenberg-Marquardt parameter fitting until parameters cannot be improved
     * or maximum number of iterations is reached.
     * @param maxIter
     * 	Maximum allowed number of iterations
     * @param verbose
     * 	Boolean flag controlling the verbose-ness of the logging
     */
    public STATUS fit(int maxIter, boolean verbose) {

        // Covariance weighted chi-square for current parameter set
        double chi2_initial = getChi2();
        
        if(verbose) {
        	System.out.println("LMA: "+N+" data and "+M+" parameters");
        	System.out.println("LMA: Initial chi2 = "+String.format("%3.3f", chi2_initial));
        }
                
        // Get suitable starting value for damping parameter, from 10^{-3}
        // times the average of the diagonal elements of JTWJ:
        Matrix J = getJacobian();
        Matrix JTWJ = J.transpose().times(dataCovariance.solve(J));
        double L=0;
        for(int i=0; i<M; i++) {
            L += JTWJ.get(i, i);
        }
        L /= (M*1000.0);
        
        if(verbose) {
        	System.out.println("LMA: Initial damping parameter L = " + L);
        }
        
        double[] lambda = {L, 10, L*maxDamping, exitTolerance};

        int nIter = 0;

        double chi2 = 0.0;
        double chi2prev = chi2_initial;
        
        // Iterate until fit converges, fails or runs out of iterations
        while((nIter++) < maxIter) {
        	
        	try {
                chi2 = iteration(lambda, chi2prev, verbose);
            }
            catch(RuntimeException re) {
                if(verbose) {
                    System.out.println("LMA: Singular update matrix on exit");
                }
                // Convergence indeterminate; no more iterations though
                return STATUS.FAILED_SINGULAR_MATRIX;
            }
        	
        	if(verbose) {
        		System.out.println("LMA: Iteration "+nIter+ " complete");
        	}
        	
        	// Check that the damping parameter remains in allowed range
        	if(lambda[0]>L*maxDamping) {
	            if(verbose) {
	                System.out.println("LMA: Damping threshold exceeded");
	            }
	
	            // Not converged; exceeded damping parameter so don't continue.
	            return STATUS.FAILED_DAMPING_THRESHOLD_EXCEEDED;
        	}
        	
            // Relative change in chi-squared; negative values indicate improvement in fit
            double rrise = (chi2-chi2prev)/chi2;

            // Residuals dropped by an amount greater than exit tolerance.
            // Succesful LM iteration. Shrink damping parameter and quit loop.
            if (rrise < -exitTolerance) {
            	
            	// Good step!
                
                if(verbose) {
                    System.out.println("LMA: Good step. New chi2 = "+chi2);
                }
                
                // Callback function
                goodStepCallback() ;
                
                // Shrink damping factor, taking care to avoid zeroing it.
                lambda[0] = Math.max(Double.MIN_VALUE, lambda[0] / lambda[1]);

                chi2prev = chi2;
            }

            // Exit tolerance exceeded: residuals changed by a very small
            // amount. We appear to be at the minimum.
            else if (Math.abs(rrise) < exitTolerance) {
            	
                // Cannot improve parameters - no further iterations
                if(verbose) {
                    System.out.println("LMA: Residual threshold exceeded.");
                }
                
                if(verbose) {
                	// Chi-square on exit
                	System.out.println("LMA: Number of iterations = "+nIter);
                	System.out.println("LMA: Final chi2 = "+ String.format("%3.3f", chi2));
                	System.out.println("LMA: Reduced chi2 = "+ String.format("%3.3f", getReducedChi2()));
                	System.out.println("LMA: Reduction factor = "+ String.format("%3.3f", chi2_initial/chi2)+"\n");
                }
                
                // Converged
                return STATUS.SUCCESS;
            }
        }
        
    	return STATUS.FAILED_ITERATION_LIMIT_EXCEEDED;
    }

    /**
     * Each call performs one iteration of parameters.
     *
     * @param   lambda  Damping parameters: element 0 is value of lambda, 1 is boost factor,
     * 					2 is the max permitted damping, 3 is exit tolerance on residual change
     *                  between steps.
     *                  
     * @param chi2prev
     * 	Chi-squared for the current parameters.
     *
     * @return status of the fit following the current iteration
     */
    private double iteration(double[] lambda, double chi2prev, boolean verbose) {
        
        /**
         * Entries of array lambda are as follows:
         * [0] -    Value of the damping parameter. This will be overwritten
         *          by the method to provide the initial value for the next iteration
         * [1] -    Boost/shrink factor on unsuccessful/successful steps
         * [2] -    Max damping
         * [3] -    Exit tolerance
         */

        // Now get Jacobian matrix for current parameters
        Matrix J = getJacobian();
        
        // Compute RHS of LM update equation: 
        // (J^T*W*J + diag(J^T*W*J))*delta = J^T*W*(residuals)
        // Note use of solve() method to avoid having to invert covariance
        // matrix to get W*(residuals).
        Matrix RHS = J.transpose().times(dataCovariance.solve(getResiduals()));
        // Get J^T*W*J
        Matrix JTWJ = J.transpose().times(dataCovariance.solve(J));
        
        // Search for a good step:
        do {

            // Make damping matrix
            Matrix L = new Matrix(M,M);
            for (int i=0; i<M; i++) {
                L.set(i, i, lambda[0] * JTWJ.get(i, i));
            }

            // Add this to Grammian
            Matrix LHS = JTWJ.plus(L);

            Matrix delta = LHS.solve(RHS);
            
            // Adjust parameters...
            params.plusEquals(delta);

            // Get new chi-square statistic
            double chi2 = getChi2();
            
            // XXX
//            System.out.println("Parameters update:");
//            System.out.println("dAC rate = " + params.get(0, 0));
//            System.out.println("dAL rate = " + params.get(1, 0));
//            System.out.println("dAL rate = " + params.get(0, 0));
//            System.out.println("chi2 = " + chi2);
            
            // Relative change in chi-squared; negative values indicate improvement in fit
            double rrise = (chi2-chi2prev)/chi2;
            
            if(rrise > lambda[3]) {
            	
            	// Bad step (residuals increased)! Try again with larger damping.
                // Reset parameters to values before previous nudge.
            	params.minusEquals(delta);

                // Boost damping parameter and try another step.
                lambda[0] *= lambda[1];
                
                if(verbose) {
                    System.out.println("LMA: Bad step. New lambda = "+lambda[0]);
                }
            }
            else {
            	return chi2;
            }
        }
        // Check damping parameter remains within allowed range
        while (lambda[0]<=lambda[2]);
        
        // Damping threshold exceeded.
        return Double.NaN;
    }


    /** 
     * Get x - f(x) 
     */
    private Matrix getResiduals(){
        return getData().minus(getModel());
    }
    
    /** 
     * Chi-square statistic, (x - f(x))^T*C^{-1}*(x - f(x)) 
     */
    public double getChi2(){
        
        // Get residuals vector
        Matrix r = getResiduals();
        
        // Covariance weighted chi-square. Note use of solve() method
        // to avoid having to invert covariance matrix to get weight matrix.
        return r.transpose().times(dataCovariance.solve(r)).get(0,0);        
    
    }
    
    /**
     * Reduced Chi-square statistic.
     */
    public double getReducedChi2(){
        return getChi2()/getDOF();
    }
    
    /** 
     * Degrees of freedom of fit.
     */
    public double getDOF(){ 
        return N - M;
    }
    
    /**
     * This method estimates parameter covariance by propagating data 
     * covariance through the system using the following equation:
     * 
     * S_p = (dp/dx)^T S_x (dp/dx)
     * 
     * It uses a fourth order central difference approximation for the 
     * parameter/data Jacobian.
     * 
     * Note that this method fails for functions that are significantly 
     * non-linear within a STDDEV or two of current solution.
     */
    public Matrix getFourthOrderCovariance(){
        
        Matrix dpdx = getJacobian_dpdx();
                
        // First order propagation.
        return dpdx.transpose().times(dataCovariance.times(dpdx));
    }
    
    /** 
     * Get the covariance matrix for parameters. This method has been tested
     * against Gnuplot 'fit' function and provides the same asymptotic
     * standard error and parameter correlation.
     */
    public Matrix getParameterCovariance(){

        // Get Jacobian matrix for current parameter set
        Matrix J = getJacobian();  
        
        // Get J^T*W*J using new weight matrix
        Matrix JTWJ = J.transpose().times(dataCovariance.solve(J));
        
        // Apply scale factor
        JTWJ.timesEquals(getDOF() / getChi2());
        
        // Invert...
        return JTWJ.inverse();
    }
    
    /**
     * Get the asymptotic standard error for the parameters
     */
    public Matrix getAsymptoticStandardError(){
        
        // Get covariance matrix for current parameter set
        Matrix cov = getParameterCovariance();        
        
        Matrix err = new Matrix(M,1);
    
        for(int p=0; p<M; p++) {
            err.set(p, 0, Math.sqrt(cov.get(p, p)));
        }
        
        return err;
    }
    
    /**
     * Get the correlation matrix for the parameters
     */
    public Matrix getParameterCorrelation() {
    
        Matrix cov = getParameterCovariance();
        
        Matrix err = getAsymptoticStandardError();
        
        Matrix err2 = err.times(err.transpose());
        
        Matrix corr = new Matrix(M,M);
    
        for(int i=0; i<M; i++) {
            for(int j=0; j<M; j++) {
                corr.set(i, j, cov.get(i, j)/err2.get(i, j));
            }
        }
        
        return corr;
    }
    
    /** 
     * Finite difference Jacobian approximation. This is the derivative of the
     * parameters solution with respect to the data, useful in estimating the
     * parameters covariance with respect to the data.
     * 
     * Jacobian is N rows by M columns.
     * 
     * This uses fourth-order central difference approximation for the first
     * derivative.
     * 
     */
    private Matrix getJacobian_dpdx() {

        // Set up Jacobian matrix: getDataN() rows by getParametersN() columns
        Matrix jac = new Matrix(N, M);

        // Matrix used to adjust data: getDataN() rows X 1 column
        Matrix delta = new Matrix(N,1);
        
        // Current parameter set, used to return parameters to original
        // values once method has completed
        Matrix origParams = params.copy();

        // Iterate over each data point...
        for (int j=0; j<N; j++) {
            
            // Get f(x+2h)
            for (int k=0; k<N; k++) {
              delta.set(k, 0, (k==j) ? h : 0.0);    
            }
            data.plusEquals(delta.times(2.0));
            fit(500,false);
            
            // -f(x+2h)
            for (int i = 0; i < M; i++) {
              jac.set(j, i, -params.get(i,0));
            }

            // Get f(x+h) 
            data.plusEquals(delta.times(-1.0));
            //
            fit(500,false);
            // +8f(x+h)
            for (int i = 0; i < M; i++)
              jac.set(j, i,  jac.get(j, i) + 8.0*params.get(i,0));

            
            // Get f(x-h)
            data.plusEquals(delta.times(-2.0));
            //
            fit(500,false);
            // -8f(x-h)
            for (int i = 0; i < M; i++) {
              jac.set(j, i,  jac.get(j, i) - 8.0*params.get(i,0));
            }
            
            // Get f(x-2h)
            data.plusEquals(delta.times(-1.0));
            //
            fit(500,false);
            // +f(x-2h)
            for (int i = 0; i < M; i++)
              jac.set(j, i,  jac.get(j, i) + params.get(i,0));
            
            
            // Divide by step size to complete derivative approximation
            for (int i = 0; i < M; i++)
              jac.set(j, i, jac.get(j, i) / (12.0*h));

            // Return data to its original value
            data.plusEquals(delta.times(2.0));
            
            // Return parameters solution to original value
            params = origParams.copy();

        }
        
        return jac;
    }

}
