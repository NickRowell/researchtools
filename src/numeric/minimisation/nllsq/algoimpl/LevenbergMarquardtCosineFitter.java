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

import java.util.List;

import numeric.minimisation.nllsq.algo.LevenbergMarquardt;

/**
 * Class fits a function of the form:
 * 
 * x(t) = A*Cos(w*t + p) + B
 *
 * to data of the form (t,x).
 * 
 * Note that in the current implementation the solution is sensitive to the starting parameter values
 * and can easily get stuck in local minima.
 * 
 * 
 * TODO: determine reasonable starting values for the parameters programmatically.
 * 
 *
 * @author nrowell
 * @version $Id$
 */
public class LevenbergMarquardtCosineFitter extends LevenbergMarquardt {
	
	/**
	 * Cosine amplitude
	 */
	public double A = 0.25;
	
	/**
	 * Bias (offset from zero)
	 */
	public double B = -0.4;
	
	/**
	 * Frequency [rad/s]
	 */
	public double w = 2.0 * Math.PI;
	
	/**
	 * Phase offset [rad]
	 */
	public double p = 0.4 * 2.0 * Math.PI;
	
	/**
	 * The abscissa values corresponding to each observed data value.
	 */
	double[] abscissa;
	
	/**
	 * 
	 * @param observations
	 */
	public LevenbergMarquardtCosineFitter(List<double[]> observations) {
		
		abscissa = new double[observations.size()];
		
		// Build data array
		double[] data = new double[observations.size()];
		double[][] cov = new double[observations.size()][observations.size()];
		
		for(int i=0; i<data.length; i++) {
			abscissa[i] = observations.get(i)[0];
			data[i] = observations.get(i)[1];
			cov[i][i] = 1.0;
		}
		
		setData(data);
	    setCovariance(cov);
	    
	}
	
	/**
	 * Perform the fit and extract the solution.
	 */
	public void invoke() {
		
		setInitialGuessParameters(new double[]{A, B, w, p});
		
		// Perform the optimization
	    fit(500, true);
	    
	    // Extract the solution
		double[] solution = getParametersSolution();
		A = solution[0];
		B = solution[1];
		w = solution[2];
		p = solution[3];
	}
	
	@Override
	public double[] getModel(double[] params) {
		
		// Extract parameter values
		double A = params[0];
		double B = params[1];
		double w = params[2];
		double p = params[3];
			
		// Build model values array
		double[] model = new double[N];
		
		// Implement the parameters -> model function
		for(int i=0; i<N; i++) {
			
			// Get the data point
			double t = abscissa[i];
			double ft = A * Math.cos(w*t + p) + B;
			model[i] = ft;
		}
		
		return model;
	}

	@Override
	public double[][] getJacobian(double[] params) {
		
		double[][] jac = new double[N][M];
		
		// Compute the Jacobian elements
		
		// FIXME: need to use the current parameter values passed into the method!
		
		
		// Loop over the data points
		for(int i=0; i<N; i++) {
			
			// Get the data point
			double t = abscissa[i];
			
			// Partial derivative wrt A
			jac[i][0] = Math.cos(w * t + p);
			
			// Partial derivative wrt B
			jac[i][1] = 1.0;

			// Partial derivative wrt w
			jac[i][2] = -1.0 * A * t * Math.sin(w * t + p);

			// Partial derivative wrt p
			jac[i][3] = -1.0 * A * Math.sin(w * t + p);
		}
		
		return jac;
	}
	
}