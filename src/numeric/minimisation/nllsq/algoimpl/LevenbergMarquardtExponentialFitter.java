/*
 * Gaia CU5 DU10
 *
 * (c) 2005-2020 Gaia Data Processing and Analysis Consortium
 *
 *
 * CU5 photometric calibration software is free software; you can redistribute
 * it and/or modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1 of the
 * License, or (at your option) any later version.
 *
 * CU5 photometric calibration software is distributed in the hope that it will
 * be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
 * General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this CU5 software; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 *
 *-----------------------------------------------------------------------------
 */

package numeric.minimisation.nllsq.algoimpl;

import numeric.minimisation.nllsq.algo.LevenbergMarquardt;

/**
 * Class fits a function of the form:
 * 
 * f(x) = A * x^{B}
 *
 * to data of the form (x,y), with free parameters A,B.
 * 
 * TODO: determine reasonable starting values for the parameters programmatically.
 * 
 * @author nrowell
 * @version $Id$
 */
public class LevenbergMarquardtExponentialFitter extends LevenbergMarquardt {
	
	/**
	 * Amplitude
	 */
	public double A = 1.0;
	
	/**
	 * Exponent
	 */
	public double B = 1.0;
	
	/**
	 * The abscissa values corresponding to each observed data value.
	 */
	double[] abscissa;
	
	/**
	 * 
	 * @param observations
	 */
	public LevenbergMarquardtExponentialFitter(double[] x, double[] y) {
		
		abscissa = new double[x.length];
		
		// Build data array
		double[] data = new double[x.length];
		double[][] cov = new double[x.length][x.length];
		
		for(int i=0; i<data.length; i++) {
			abscissa[i] = x[i];
			data[i] = y[i];
			cov[i][i] = 1.0;
		}
		
		setData(data);
	    setCovariance(cov);
	}
	
	/**
	 * Perform the fit and extract the solution.
	 */
	public void invoke() {
		
		setInitialGuessParameters(new double[]{A, B});
		
		// Perform the optimization
	    fit(500, true);
	    
	    // Extract the solution
		double[] solution = getParametersSolution();
		A = solution[0];
		B = solution[1];
	}
	
	@Override
	public double[] getModel(double[] params) {
		
		// Extract parameter values
		double A = params[0];
		double B = params[1];
			
		// Build model values array
		double[] model = new double[N];
		
		// Implement the parameters -> model function
		for(int i=0; i<N; i++) {
			
			// Get the data point
			double x = abscissa[i];
			double fx = A * Math.pow(x, B);
			model[i] = fx;
		}
		
		return model;
	}

	@Override
	public double[][] getJacobian(double[] params) {
		
		double[][] jac = new double[N][M];
		
		// Compute the Jacobian elements
		
		double A = params[0];
		double B = params[1];
		
		// Loop over the data points
		for(int i=0; i<N; i++) {
			
			// Get the data point
			double x = abscissa[i];
			
			// Partial derivative wrt A
			jac[i][0] = Math.pow(x, B);
			
			// Partial derivative wrt B
			jac[i][1] = Math.log(x) * A * Math.pow(x, B);
		}
		
		return jac;
	}
}