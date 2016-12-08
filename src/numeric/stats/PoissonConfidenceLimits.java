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

package numeric.stats;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;

/**
 * This class provides confidence limits for variables described by Poisson processes.
 * See e.g. Gehrels, N., (1986) "Confidence Limits for Small Number of Events in Astrophysical Data".
 *
 * @author nrowell
 * @version $Id$
 */
public class PoissonConfidenceLimits {
	
	/**
	 * Gets the double-sided confidence limits on the mean.
	 * <p>
	 * For example, given <code>n</code> events, the true mean lies in the range [ll:ul]
	 * with <code>dsCl</code>*100% probability, where
	 * 
	 * {ll,ul} = getDoubleSidedConfidenceLimits(n, dsCl)
	 * 
	 * @param n
	 * 	The observed number of events
	 * @param dsCl
	 * 	Double-sided confidence level [0.0:1.0]
	 * @return
	 * 	Array containing the lower and upper double-sided confidence limits on the true mean
	 */
	public static double[] getDoubleSidedConfidenceLimits(int n, double dsCl) {
		
		// Translate double-sided condifence level to equivalent single-sided level
		double ssCl = (1.0 + dsCl)/2.0;
		
		double ll = getLowerSingleSidedConfidenceLimit(n, ssCl);
		double ul = getUpperSingleSidedConfidenceLimit(n, ssCl);
		
		return new double[]{ll, ul};
	}
	
	/**
	 * Gets the upper single-sided confidence limit on the mean. 
	 * <p>
	 * For example, given <code>n</code> events, the true mean lies in the range [0:ul]
	 * with <code>cl</code>*100% probability, where
	 * 
	 * ul = getUpperSingleSidedConfidenceLimit(n, cl)
	 * 
	 * @param n
	 * 	The observed number of events
	 * @param ssCl
	 * 	Single-sided confidence level [0.0:1.0]
	 * @return
	 * 	The upper single-sided confidence limit on the true mean
	 */
	public static double getUpperSingleSidedConfidenceLimit(int n, double ssCl) {
		
		// Degrees of freedom of the corresponding Chi-squared distribution
		double dof = 2*n + 2;
		
		// Solution exploits the relationship between the cumulative Poisson
		// distribution and the chi-square distribution.
		ChiSquaredDistribution chi2 = new ChiSquaredDistribution(dof);
		
		return chi2.inverseCumulativeProbability(ssCl)/2.0;
	}
	
	/**
	 * Gets the lower single-sided confidence limit on the mean. 
	 * <p>
	 * For example, given <code>n</code> events, the true mean lies in the range [ll:inf]
	 * with <code>cl</code>*100% probability, where
	 * 
	 * ll = getLowerSingleSidedConfidenceLimit(n, cl)
	 * 
	 * @param n
	 * 	The observed number of events
	 * @param ssCl
	 * 	Single-sided confidence level [0.0:1.0]
	 * @return
	 * 	The lower single-sided confidence limit on the true mean
	 */
	public static double getLowerSingleSidedConfidenceLimit(int n, double ssCl) {
		
		if(n==0) {
			return 0.0;
		}
		
		// Degrees of freedom of the corresponding Chi-squared distribution
		double dof = 2*n;

		// Solution exploits the relationship between the cumulative Poisson
		// distribution and the chi-square distribution.
		ChiSquaredDistribution chi2 = new ChiSquaredDistribution(dof);
		
		return chi2.inverseCumulativeProbability(1.0 - ssCl)/2.0;
	}
	
	/**
	 * Test application prints out confidence limits for direct comparison with the tables in the
	 * Gehrels, N. (1986) paper.
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		
		// Single-sided confidence limits to match tables in Gehrels paper
		double[] cls = new double[]{0.8413, 0.9, 0.95, 0.975, 0.9772, 0.99};
		
		System.out.println("Poisson single-sided upper confidence limits for n events");
		
		System.out.print("\nn\t");
		for(double cl : cls) {
			System.out.print(String.format("%.4f\t", cl));
		}
		System.out.println("\n");
		
		for(int n=0; n<20; n++) {
			
			System.out.print(n+"\t");
			for(double cl : cls) {
				
				double ul = getUpperSingleSidedConfidenceLimit(n,cl);
				System.out.print(String.format("%.4f\t", ul));
			}
			System.out.println();
		}
		
		System.out.println("\n\nPoisson single-sided lower confidence limits for n events");
		
		System.out.print("\nn\t");
		for(double cl : cls) {
			System.out.print(String.format("%.4f\t", cl));
		}
		System.out.println("\n");
		
		for(int n=0; n<20; n++) {
			
			System.out.print(n+"\t");
			for(double cl : cls) {
				
				double ll = getLowerSingleSidedConfidenceLimit(n,cl);
				System.out.print(String.format("%.4f\t", ll));
			}
			System.out.println();
		}
	}
	
}
