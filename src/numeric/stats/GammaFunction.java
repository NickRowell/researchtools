package numeric.stats;

import java.util.NavigableMap;
import java.util.TreeMap;

import numeric.integration.IntegrableFunction;
import numeric.integration.IntegrationUtils;

/**
 * Class provides utilities to compute values of the Gamma function.
 *
 *
 * @author nrowell
 * @version $Id$
 */
public class GammaFunction {
	
	/**
	 * Current upper limit on k that's been added to the map
	 */
	private static int kMax = 2;
	
	/**
	 * Precomputed values of the Gamma function at half-integer values
	 */
	private static NavigableMap<Integer, Double> gammaHalfInt = new TreeMap<>();
	
	static {
		gammaHalfInt.put(1, Math.sqrt(Math.PI));  // k=1; value of Gamma(k/2) = Gamma(0.5)
		gammaHalfInt.put(2, 1.0);                 // k=2; value of Gamma(k/2) = Gamma(1.0)
	}
	
	
	/**
	 * Compute the Gamma function when the argument is an integer divided by two. This is especially
	 * useful when computing certain probability distributions.
	 * 
	 * @param k
	 * 	The integer that is the numerator.
	 * @return
	 * 	Gamma(k/2).
	 */
	public static double gammaFunctionHalfInt(int k) {
		
		// Lazy initialisation
		if(k > kMax) {
			
			// Compute all absent values up to k
			for(int kk = kMax + 1; kk <= k; kk++) {
				
				// kk/2 - 1
				double kk2m1 = (kk / 2.0) - 1.0;
				
				// Retrieve the cached value for Gamma(kk/2 - 1)
				double gammaKk2m1 = gammaHalfInt.get(kk-2);
				
				// Value of the Gamma function for (kk/2) 
				double gammaKk2 = kk2m1 * gammaKk2m1;
				
				gammaHalfInt.put(kk, gammaKk2);
			}
			
			kMax = k;
		}
		
		return gammaHalfInt.get(k);
	}
	
	/**
	 * Compute the Gamma function for general floating point argument.
	 * 
	 * @param z
	 * 	The argument of the Gamma function.
	 * @return
	 * 	Gamma(x).
	 */
	public static double gammaFunction(double z) {
		
		GammaInt g = new GammaInt(z);
		
		// Parameters of the numerical integration
		double xMax = 1000.0;
		double h = 0.01;
		
		return IntegrationUtils.integrate(g, 0.0, xMax, h);
	}
	
	/**
	 * Class provides an {@link IntegrableFunction} that implements the integrand of the Gamma function.
	 *
	 * @author nrowell
	 * @version $Id$
	 */
	private static class GammaInt implements IntegrableFunction {

		/**
		 * The parameter of the Gamma function.
		 */
		private double z;
		
		/**
		 * Main constructor.
		 * 
		 * @param z
		 * 	The parameter of the Gamma function.
		 */
		public GammaInt(double z) {
			this.z = z;
		}
		
		@Override
		public double evaluate(double x) {
			return Math.pow(x, z-1) * Math.exp(-x);
		}
	}
	
	/**
	 * Test function to demonstrate usage.
	 * 
	 * @param args
	 * 	The command line args (ignored)
	 */
	public static void main(String[] args) {
		
		System.out.println("k\tk/2\tΓ(k/2)\tΓ(k/2)");
		System.out.println("-\t---\t------\t------");
		
		for(int k=1; k<20; k++) {
			System.out.println(k + "\t" + (k/2.0) + "\t" + gammaFunctionHalfInt(k) + "\t" + gammaFunction(k/2.0));
		}
	}
	
	
}
