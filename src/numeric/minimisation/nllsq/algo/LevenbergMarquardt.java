package numeric.minimisation.nllsq.algo;

import Jama.Matrix;

/**
 * Implementation of the Levenberg-Marquardt nonlinear least squares algorithm.
 * 
 * In order to use this class you must extend it and implement the getModel(Matrix params)
 * method.
 * 
 * Normal usage:
 * 
  	LevenbergMarquardt lma = new LevenbergMarquardt() {
  
  		@Override
		public double[] getModel(double[] params) {
			
			double[] model = new double[N];
			
			// Implement the parameters -> model function
			
			return model;
		}

		@Override
		public double[][] getJacobian(double[] params) {
			
			double[][] jac = new double[N][M];
			
			// Compute the Jacobian elements
			
			return jac;
		}
  	}
  
  	lma.setData(electronSamples);
    lma.setCovariance(cov);
    lma.setInitialGuessParameters(params);
    // Perform the optimization
    lma.fit(500, true);
			        
	// Extract the solution
	double[] solution = lma.getParametersSolution();
 * 
 * 
 * 
 * 
 * 
 *
 *
 * @author nrowell
 * @version $Id$
 */
public abstract class LevenbergMarquardt{

    /** 
     * Absolute step size used in finite difference Jacobian approximation,
     * for Jacobian of parameter solution with respect to data.
     */
    private double h = 1E-2;   

    /** 
     * Exit tolerance 
     */
    private double EXIT_TOL = 1E-32;

    /**
     * Max damping scale factor. This is multiplied by automatically selected
     * starting value of damping parameter, which is 10^{-3}
     * times the average of the diagonal elements of J^T*J
     */
    private double MAX_DAMP = 1E32;
    

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
     * Set the absolute step size used in the finite difference Jacobian estimation.
     * @param D
     * 	The finite step size
     */
    public void setDELTAD(double h)  { 
    	this.h   = h;
    }
    
    /**
     * Set the exit tolerance - if the (absolute value of the) relative change in the chi-square
     * from one iteration to the next is lower than this, then we're at the minimum and the fit
     * is halted.
     * @param E
     * 	The exit tolerance to set
     */
    public void setEXIT_TOL(double E) { 
    	EXIT_TOL = E;
    }
    
    /**
     * Set the maximum damping factor. If the damping factor becomes larger than this during the fit,
     * then we're stuck and cannot reach a better solution.
     * 
     * @param M
     * 	The max damping factor to set
     */
    public void setMAX_DAMP(double M) {
    	MAX_DAMP = M;
    }

    /**
     * Set the Nx1 column vector of observed values
     * @param data
     * 	The Nx1 column vector of observed values
     */
    public void setData(double[] data) {
    	this.N = data.length;
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
    	this.M = params.length;
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
     * Gets the Jacobian for the current parameters values.
     * 
     * @return
     */
    private  Matrix getJacobian() {
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
     * Jacobian is getDataN() rows by getParametersN() columns.
     * 
     */
//    public Matrix getJacobian(){
//        
//        // Set up Jacobian matrix: getDataN() rows by getParametersN() columns
//        Matrix jac = new Matrix(getDataN(), getParametersN());
//
//        // Get parameter values prior to any tinkering
//        Matrix p = getParameters();
//                
//        // Matrix used to adjust parameters: getParametersN() rows X 1 column
//        Matrix delta = new Matrix(getParametersN(),1);
//        
//        // Get finite step sizes to use for each parameter
//        Matrix FINITE_STEP = getParametersSteps();
//
//        // Iterate over each parameter...
//        for (int j=0; j<getParametersN(); j++)
//        {
//            
//            // Make 'nudge' vector for this parameter. Set all entries to zero
//            // except the parameter that is to be adjusted
//            for (int k=0; k<getParametersN(); k++)
//              delta.set(k, 0, (k==j) ? FINITE_STEP.get(j,0) : 0.0);
//              
//            // First advance this single parameter
//            setParameters(p.plus(delta));
//            
//            // Get model values for advanced parameter set: f(x+h)
//            Matrix R = getModel(params);
//
//            // Build Jacobian by finite difference
//            for (int i = 0; i < getDataN(); i++)
//              jac.set(i, j, R.get(i,0));
//
//            // Now retard the parameter...
//            setParameters(p.minus(delta));
//
//            // Get model values for retarded parameter set: f(x-h)
//            R = getModel(params);
//
//            // Build Jacobian by finite difference 
//            for (int i = 0; i < getDataN(); i++)
//              jac.set(i, j,  (jac.get(i, j) - R.get(i,0)) / (2.0*FINITE_STEP.get(j,0)));
//
//        }
//        
//        // Return parameters to original values
//        setParameters(p);
//        
//        return jac;
//
//    }
    
    
    /**
     * Perform LM iteration loop until parameters cannot be improved.
     *
     * @param MAX_ITER  Maximum number of allowed iteration before convergence.
     *
     */
    public void fit(int MAX_ITER, boolean verbose){

        // Covariance weighted chi-square for current parameter set
        double chi2_initial = getChi2();
        
        if(verbose) System.out.println("LMA: "+N+" data and "+M+" parameters");
        if(verbose) System.out.println("LMA: Initial chi2 = "+String.format("%3.3f", chi2_initial));
                
        // get suitable starting value for damping parameter, from 10^{-3}
        // times the average of the diagonal elements of JTWJ:
        
        Matrix J = getJacobian();
        
        // Note use of A.solve(B) method of Jama to get A^{-1}B without
        // having to invert the covariance matrix.
        Matrix JTWJ = J.transpose().times(dataCovariance.solve(J));
        double L=0;
        for(int i=0; i<M; i++)
            L += JTWJ.get(i, i);
        L /= (M*1000.0);

        double[] lambda = {L, 10, L*MAX_DAMP, EXIT_TOL};

        int N_ITER = 0;

        while(!iteration(lambda, verbose) && N_ITER<MAX_ITER){

            // Un-comment the following for very verbose behaviour
//            chi2 = getChi2();
//            System.out.println("LMA: Iteration "+N_ITER+
//                               " complete, " + "residual = "+chi2);
            

            N_ITER++;
        }

        // Chi-square on exit
        double chi2_final = getChi2();

        if(verbose) System.out.println("LMA: Number of iterations = "+N_ITER);
        if(verbose) System.out.println("LMA: Final chi2 = "+ String.format("%3.3f", chi2_final));
        if(verbose) System.out.println("LMA: Reduced chi2 = "+ String.format("%3.3f", getReducedChi2()));
        if(verbose) System.out.println("LMA: Reduction factor = "+ String.format("%3.3f", chi2_initial/chi2_final)+"\n");        

        return;

    }

    /**
     * Each call performs one iteration of parameters.
     *
     * @param   lambda  Damping parameters: element 0 is value of lambda, 1 is boost factor,
     * 					2 is the max permitted damping, 3 is
     *                  shrink factor, 4 is exit tolerance on residual change
     *                  between steps.
     *
     * @return Boolean  States whether another iteration would be appropriate, or
     *                  if change in residuals and/or damping thresholds have
     *                  been reached
     */
    private boolean iteration(double[] lambda, boolean verbose)
    {
        
        /**
         * Entries of array lambda are as follows:
         * [0] -    Value of the damping parameter. This will be overwritten
         *          by the method to provide the initial value for the next iteration
         * [1] -    Boost/shrink factor on unsuccessful/successful steps
         * [2] -    Max damping
         * [3] -    Exit tolerance
         */

        // chi-square prior to parameter update
        double chi2prev = getChi2();
        
        // Now get Jacobian matrix for current parameters
        Matrix J = getJacobian();
        
        // Compute RHS of LM update equation: 
        // (J^T*W*J + diag(J^T*W*J))*delta = J^T*W*(residuals)
        // Note use of solve() method to avoid having to invert covariance
        // matrix to get W*(residuals).
        Matrix RHS = J.transpose().times(dataCovariance.solve(getResiduals()));

        // Get J^T*W*J
        Matrix JTWJ = J.transpose().times(dataCovariance.solve(J));
               
        // Change in chi-square from one iteration to the next
        double rrise = 0;

        // Exit status
        boolean done = true;

        // Search for a good step:
        do
        {

            // Make damping matrix
            Matrix L = new Matrix(M,M);

            for (int i=0; i<M; i++)
                L.set(i, i, lambda[0]);

            // Add this to Grammian
            Matrix LHS = JTWJ.plus(L);

            Matrix delta = new Matrix(M,1);

            try
            {
                delta = LHS.solve(RHS);
            }
            catch(RuntimeException re)
            {
                if(verbose)
                    System.out.println("LMA: Singular update matrix on exit");
                return true;
            }

            // Adjust parameters...
            params.plusEquals(delta);

            // Get new chi-square statistic
            double chi2 = getChi2();
            
            // if rrise is negative, then current residuals are lower than
            // those found on previous step
            rrise = (chi2-chi2prev)/chi2;

            // Residuals dropped by an amount greater than exit tolerance.
            // Succesful LM iteration. Shrink damping parameter and quit loop.
            if (rrise < -lambda[3])
            {
                
                //if(verbose)
                //    System.out.println("LMA: Good step. New chi2 = "+chi2);
                // Good step! Want more iterations.
                done = false;
                lambda[0] /= lambda[1];
                
                // Callback method to notify of successful parameter update.
//                updateCallback();
                
                break;
            }

            // Exit tolerance exceeded: residuals changed by a very small
            // amount. We appear to be at the minimum, so keep previous
            // parameters and quit loop. Algorithm cannot find a better value.
            else if (Math.abs(rrise) < lambda[3])
            {
                params.minusEquals(delta);
                
                //System.out.println("Meh");
                // Cannot improve parameters - no further iterations
                if(verbose)
                    System.out.println("LMA: Residual threshold exceeded.");
                break;

            }
            else
            {
                
                // Bad step (residuals increased)! Try again with larger damping.
                // Reset parameters to values before previous nudge.
            	params.minusEquals(delta);

                // Boost damping parameter and try another step.
                lambda[0] *= lambda[1];
            }
            
        }
        // Check damping parameter remains within allowed range
        while (lambda[0]<=lambda[2]);
        
        if(lambda[0]>lambda[2] && verbose)
            System.out.println("LMA: Damping threshold exceeded");

        return done;
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
        Matrix R = getResiduals();
        
        // Covariance weighted chi-square. Note use of solve() method
        // to avoid having to invert covariance matrix to get weight matrix.
        return R.transpose().times(dataCovariance.solve(R)).get(0,0);        
    
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
        
        // Get inverse of covariance matrix
        // TODO: shouldn't take the inverse here - use dataCovariance.solve(J) further on...
        Matrix W = dataCovariance.inverse();
        
        // Weighted sum of squared residuals
        double WSSR = getChi2();
        
        // Degrees of freedom
        double dof = getDOF();
        
        //System.out.println("RMS of residuals = "+Math.sqrt(WSSR/dof));
        
        // Weight terms in matrix W.
        W.timesEquals(dof / WSSR);
        
        // Get J^T*W*J using new weight matrix
        Matrix JTWJ = J.transpose().times(W).times(J);
        
        // Invert...
        return JTWJ.inverse();
    
    }
    
    
    /** Get the asymptotic standard error for the parameters */
    public Matrix getAsymptoticStandardError(){
        
        // Get covariance matrix for current parameter set
        Matrix COVAR = getParameterCovariance();        
        
        Matrix STD_ERR = new Matrix(M,1);
    
        for(int p=0; p<M; p++)
            STD_ERR.set(p, 0, Math.sqrt(COVAR.get(p, p)));
        
        return STD_ERR;
    }
    
    /** Get the correlation matrix for the parameters */
    public Matrix getParameterCorrelation(){
    
        Matrix COVAR = getParameterCovariance();
        
        Matrix STD_ERR = getAsymptoticStandardError();
        
        Matrix STD_ERR2 = STD_ERR.times(STD_ERR.transpose());
        
        Matrix CORR = new Matrix(M,M);
    
        for(int i=0; i<M; i++)
            for(int j=0; j<M; j++)            
                CORR.set(i, j, COVAR.get(i, j)/STD_ERR2.get(i, j));
        
        return CORR;
    }
    
    /** 
     * Finite difference Jacobian approximation. This is the derivative of the
     * parameters solution with respect to the data, useful in estimating the
     * parameters covariance with respect to the data.
     * 
     * Jacobian is getDataN() rows by getParametersN() columns.
     * 
     * This uses fourth-order central difference approximation for the first
     * derivative.
     * 
     */
    private Matrix getJacobian_dpdx(){

        // Set up Jacobian matrix: getDataN() rows by getParametersN() columns
        Matrix jac = new Matrix(N, M);

        // Matrix used to adjust data: getDataN() rows X 1 column
        Matrix delta = new Matrix(N,1);
              
        // Finite difference in terms of absolute change in parameter
        double FINITE_STEP = h;
            
        // Current parameter set, used to return parameters to original
        // values once method has completed
        Matrix PARAM = params.copy();

        // Iterate over each data point...
        for (int j=0; j<N; j++)
        {
            
            // Get f(x+2h)
            //
            for (int k=0; k<N; k++)
              delta.set(k, 0, (k==j) ? FINITE_STEP : 0.0);            
            //
            data.plusEquals(delta.times(2.0));
            //
            fit(500,false);
            
            // -f(x+2h)
            for (int i = 0; i < M; i++)
              jac.set(j, i, -params.get(i,0));

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
            for (int i = 0; i < M; i++)
              jac.set(j, i,  jac.get(j, i) - 8.0*params.get(i,0));
            
            // Get f(x-2h)
            data.plusEquals(delta.times(-1.0));
            //
            fit(500,false);
            // +f(x-2h)
            for (int i = 0; i < M; i++)
              jac.set(j, i,  jac.get(j, i) + params.get(i,0));
            
            
            // Divide by step size to complete derivative approximation
            for (int i = 0; i < M; i++)
              jac.set(j, i, jac.get(j, i) / (12.0*FINITE_STEP));

            // Return data to its original value
            data.plusEquals(delta.times(2.0));
            
            // Return parameters solution to original value
            params = PARAM.copy();

        }
        
        return jac;

    }

}
