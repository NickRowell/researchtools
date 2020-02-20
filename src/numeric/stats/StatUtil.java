package numeric.stats;

import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.stat.StatUtils;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import numeric.functions.Function;

/**
 * Utility class for statistical methods.
 * 
 * @author nickrowell
 */
public final class StatUtil {
    
    /**
     * Static Random instance to avoid generating new ones with each call.
     */
    static Random random = new Random();
    
    /**
     * Private constructor to enforce non-instantiation.
     */
    private StatUtil() {
    	
    }
    
    /**
     * Used to map single values from arbitrary objects, for use in creating histograms etc.
     *
     * @author nrowell
     * @version $Id$
     * @param <T>
     */
	@FunctionalInterface
	public interface ValueMapper<T> {
		public double getValue(T object);
	}
    
	/**
	 * Bins the objects according to the defined value and normalises the bin values.
	 * 
	 * @param objects
	 * 	The objects to be used to compute the histogram.
	 * @param mapper
	 * 	The rule for deriving the statistic to be binned from the objects.
	 * @param min
	 * 	The lower limit of the range to be binned.
	 * @param max
	 * 	The upper limit of the range to be binned.
	 * @param binWidth
	 * 	The width of each bin; if there's not a whole number of bins between the limits then the uppermost bin
	 * @param normalise
	 * 	Boolean flag indicating if the histogram should be normalised.
	 * @return
	 * 	An array containing the binned & optionally normalised data.
	 */
	public static <T> double[] computeHistogram(Collection<T> objects, ValueMapper<T> mapper, double min, double max, double binWidth, boolean normalise) {
		double[] data = new double[objects.size()];
		int idx=0;
		for(T object : objects) {
        	data[idx++] = mapper.getValue(object);
		}
		return computeHistogram(data, min, max, binWidth, normalise);
	}
	
	/**
	 * Bins the data and optionally normalises the bin values.
	 * 
	 * @param data
	 * 	An array of the values to be binned.
	 * @param min
	 * 	The lower limit of the range to be binned.
	 * @param max
	 * 	The upper limit of the range to be binned.
	 * @param binWidth
	 * 	The width of each bin; if there's not a whole number of bins between the limits then the uppermost bin
	 * is extended slightly beyond the range.
	 * @param normalise
	 * 	Boolean flag indicating if the histogram should be normalised.
	 * @return
	 * 	An array containing the binned, normalised data.
	 */
	public static double[] computeHistogram(double[] data, double min, double max, double binWidth, boolean normalise) {
		
		// Compute the binned data
        int nBins = (int)Math.ceil((max - min) / binWidth);
        double[] bins = new double[nBins];
        
        // Number of points lying within the histrogram range
        int n = 0;
        
        for(double datum : data) {
	        int bin = (int)Math.floor((datum - min) / binWidth);
	    	if(bin >= 0 && bin < nBins) {
	    		bins[bin]++;
	    		n++;
	    	}
        }
        
        if(normalise) {
	        // Normalise the histogram
	        double norm = 1.0 / (n * binWidth);
	        for(int bin = 0; bin < nBins; bin++) {
	        	bins[bin] *= norm;
	        }
        }
        
        return bins;
	}
	
	/**
	 * Compute the cumulative distribution function for the statistic derived from the objects.
	 * 
	 * @param objects
	 * 	The objects to be processed.
	 * @param mapper
	 * 	Rule for compute the value to be plotted from each object.
	 * @param min
	 * 	Minimum value of the histogram.
	 * @param max
	 * 	Maximum value of the histogram.
	 * @param binWidth
	 * 	The width of each bin; if there's not a whole number of bins between the limits then the uppermost bin
	 * is extended slightly beyond the range.
	 * @param normalise
	 * 	Boolean flag indicating if the histogram should be normalised.
	 * @param ascending
	 * 	Boolean flag indicating if the cumulative distribution should be summed in the ascending or descending direction.
	 * @return
	 * 	An array containing the binned & optionally normalised cumulative distribution, computed in either the ascending or
	 * descending direction.
	 */
	public static <T> double[] computeCdf(Collection<T> objects, ValueMapper<T> mapper, double min, double max, double binWidth, boolean normalise, boolean ascending) {

		double[] histogram = computeHistogram(objects, mapper, min, max, binWidth, normalise);
		
		double[] cdf = new double[histogram.length];
		
		if(ascending) {
			cdf[0] = histogram[0] * binWidth;
			for(int i=1; i<histogram.length; i++) {
				cdf[i] = cdf[i-1] + histogram[i] * binWidth;
			}
		}
		else {
			cdf[histogram.length - 1] = histogram[histogram.length - 1] * binWidth;
			for(int i=histogram.length - 2; i>=0; i--) {
				cdf[i] = cdf[i+1] + histogram[i] * binWidth;
			}
		}
		
        return cdf;
	}
    
    /**
     * Compute the confidence limits of the value of a univariate random variable obtained by transforming
     * one or more Gaussian random variables through an arbitrary non-linear multivariable function.
     * 
     * @param f
     * 	The univariate function.
     * @param mean_x
     * 	The mean values of the input Gaussian random variables.
     * @param sigma_x
     * 	The standard deviations of the input Gaussian random variables.
     * @param lowerPerc
     * 	The lower confidence boundary, expressed as a percentage in the range [0:100)
     * @param upperPerc
     *  The upper confidence boundary, expressed as a percentage in the range (0:100]
     * @return
     * 	The lower and upper values of the transformed random variable at which the lower and upper
     * confidence limits lie.
     */
    public static double[] computeConfidenceLimits(Function f, double[] mean_x, double[] sigma_x, double lowerPerc, double upperPerc) {
    	
    	// Sanity checks on the function
    	if(f.getOutputs() != 1) {
    		throw new IllegalArgumentException("Invalid Function: outputs=1; found outputs="+f.getOutputs());
    	}
    	
    	// Sanity checks on the inputs
    	if(mean_x.length < 1 || mean_x.length != sigma_x.length) {
    		throw new IllegalArgumentException("Invalid inputs: number of means = " + mean_x.length + "; number of sigmas = " + sigma_x.length);
    	}
    	
    	// Sanity check on the confidence limits
    	if(lowerPerc < 0.0 || lowerPerc >= 100.0) {
    		throw new IllegalArgumentException("Invalid lower confidence limit ("+lowerPerc+"); must lie in range [0:100)");
    	}
    	if(upperPerc <= 0.0 || upperPerc > 100.0) {
    		throw new IllegalArgumentException("Invalid upper confidence limit ("+upperPerc+"); must lie in range (0:100]");
    	}
    	
    	// Number of Monte Carlo shots
    	int nTrials = 100000;
    	
    	// Array to buffer all outputs
    	double[] output = new double[nTrials];
    	
    	// Track the minimum and maximum values
    	double min = Double.MAX_VALUE;
    	double max = -Double.MAX_VALUE;
    	
    	for(int nTrial = 0; nTrial < nTrials; nTrial++) {
    		
    		// Random value of the inputs
    		double[] inputs = new double[mean_x.length];
    		for(int i=0; i<mean_x.length; i++) {
    			double x = mean_x[i] + random.nextGaussian() * sigma_x[i];
    			inputs[i] = x;
    		}
    		
    		// Random value of the output
    		output[nTrial] = f.compute(inputs)[0];
    		
    		min = Math.min(min, output[nTrial]);
    		max = Math.max(max, output[nTrial]);
    	}
    	
    	double lower = min;
    	double upper = max;
    	
    	if(lowerPerc > 0.0) {
    		lower = StatUtils.percentile(output, lowerPerc);
    	}
    	if(upperPerc < 100.0) {
    		upper = StatUtils.percentile(output, upperPerc);
    	}
    	
    	return new double[]{lower, upper};
    }
    
    
    /**
     * Compute the binomial coefficient 'N choose M'.
     */
    public static double nChooseM(int n, int k)
    {
        double nk = 1l;
        for (int i = n - k + 1; i <= n; i++) {
            nk *= i;
		}
		for (int i = 1; i <= k; i++) {
	            nk /= i;
		}
	        return nk;
    }
    
    /**
	 * Utility to generate realizations of a random vector described by
	 * a mean and covariance matrix.
	 * 
	 * @param covar		NxN covariance matrix on vector elements.
	 * @param mean		Nx1 column vector containing mean of random vector.
	 * @return			Random vector drawn from input distribution.
	 */
    public static Matrix drawRandomVector(Matrix covar, Matrix mean) {
        
        // How many data to draw?
        int N = covar.getColumnDimension();     
        
        // Basic checks on matrix sizes and shapes etc.
        if(covar.getColumnDimension()!=covar.getRowDimension()) {
            throw new RuntimeException();
        }
        if(mean.getColumnDimension()!=1 || mean.getRowDimension()!=N) {
            throw new RuntimeException();
        }
                
        // Get Eigenvalue decomposition of covariance matrix
        EigenvalueDecomposition eig = covar.eig();
        
        // Get dispersions along the principal axes of the covariance hyper-ellipsoid
        Matrix evals = eig.getD();
        
        // Random vector in principal axes frame
        Matrix rand = new Matrix(N,1);
        
        // Draw elements from unit Gaussian and scale according to axis sigma
        for(int d=0; d<N; d++) {
            rand.set(d, 0, Math.sqrt(evals.get(d, d))*random.nextGaussian());
        }
        
        // Transform this back to original frame using eigenvectors
        // of covariance matrix
        rand = eig.getV().times(rand);    
                
        // Add random error onto mean and return
        return mean.plus(rand);
    }
    
    /**
     * Generate an N-by-N random positive definite matrix that can be used
     * to simulate a covariance matrix. The sigma factor is a scale factor
     * and not necessarily representative of the standard deviation of the
     * matrix.
     * @return 
     */
    public static Matrix getRandomCovarianceMatrix(double sigma, int N){
                
        // Generate random covariance matrix for data
        Matrix covar = new Matrix(N,N);
        for(int i=0; i<N; i++)
            for(int j=0; j<N; j++)
                covar.set(i, j, random.nextGaussian());

        // Scale to sigma
        covar.timesEquals(sigma);
        
        // Make into a nice positive definite matrix by squaring
        covar = covar.transpose().times(covar);
        
        return covar;
    }
    
    /**
     * Compute the chi-square probability density function for k degrees of freedom.
     * 
     * @param k
     * 	Number of degrees of freedom
     * @param x
     * 	The value of chi-square
     * @return
     * 	Probability density for x given the number of degrees of freedom, i.e. P(x;k)
     */
    public static double getChi2(int k, double x) {
        
        double gamma = GammaFunction.gammaFunctionHalfInt(k);
        
        return Math.pow(x, k/2.0 - 1.0) * Math.exp(-x/2) / (Math.pow(2, k/2.0) * gamma);
    }
    
    /**
     * This method provides a 1D Gaussian convolution kernel suitable for 
     * application to binned data. It integrates the Gaussian over the width
     * of the bins (centred on central bin) to obtain the fraction of the
     * contents of the central bin that are scattered into each of the 
     * surrounding bins. The integration proceeds outward from the central bin,
     * and halts when 98% of the probability has been captured (1% lost from 
     * each wing). The returned error kernel will therefore be wider if the 
     * standard deviation is larger. The kernel is re-normalised to correct for
     * the missing probability, i.e. all the elements are scaled so that their
     * sum is one.
     * 
     * Only one half of the error kernel is returned because it is symmetric.
     * The central bin is in element 0.
     * 
     * @param sigma     Standard deviation.
     * @param bin_width Width of the bins.
     * @return 
     */
    public static double[] getSymGausKernel(double sigma, double bin_width){
 
        // Sanity checks
        if(sigma<0)
            throw new IllegalArgumentException("Sigma must be positive");
        
        if(sigma==0) return new double[]{1.0};
        
        // List to store probability in each bin.
        List<Double> values = new LinkedList<Double>();
        
        // Get integral over half-width of first bin and start accumulating 
        // total integrated probability.
        // Value of integral approaches 0.5 as kernel bins carry less and less
        // probability.        
        double integral = Gaussian.Phi(bin_width/2.0, 0, sigma) - 0.5;
        
        // Calculate first bin (centred on mean) using symmetry.
        values.add(integral*2.0);
        
        // Continue calculating bins until integral drops below threshold.
        // Stop adding kernel points once all but 2E-2 of the probability has
        // been captured by kernel.
        int bin = 1;
        
        while(integral < 0.5 - 1E-2){
            
            double bin_contents = Gaussian.Phi(bin* bin_width + bin_width/2.0, 0, sigma) -
                                  Gaussian.Phi((bin-1)* bin_width + bin_width/2.0, 0, sigma);
            
            values.add(bin_contents);
            
            // Sum all probability elements...
            integral += bin_contents;
            
            // On to next bin...
            bin++;
            
        }
        
        // Total integrated probability captured by truncated kernel
        double total = integral*2.0;
        
        // Copy values out to array and re-normalise
        double[] kernel = new double[values.size()];
        for(int b=0; b<values.size(); b++){
            // Renormalisation step
            kernel[b] = values.get(b)/total;
        }
    
        return kernel;
        
    }
    
    /**
     * Compute the sample variance for the set of data points.
     * 
     * @param x
     * 	The datapoints.
     * @return
     * 	The sample variance
     */
    public static double getSampleCovariance(double[] x) {
    	
    	// Load into 2D array in order to use general method
    	double[][] data = new double[x.length][1];
    	for(int i=0; i<x.length; i++) {
    		data[i] = new double[]{x[i]};
    	}
    	
    	return getSampleCovariance(data)[0][0];
    }
    
    /**
     * Compute the sample covariance matrix for the set of data points.
     * 
     * @param x
     * 	Two dimensional array containing the data points for which the sample covariance is to be computed.
     * The leading element indexes the data points, the trailing element indexes the dimensions of each 
     * data point.
     * @return
     * 	The sample covariance matrix; for N-dimensional data this array has size dimensions [N][N].
     */
    public static double[][] getSampleCovariance(double[][] x) {
    	
    	double[] mean = getSampleMean(x);
    	
    	// Number of datapoints
    	int n = x.length;
    	// Number of dimensions
    	int d = x[0].length;
    	
    	double[][] s2 = new double[d][d];
    	for(int i=0; i<x.length; i++) {
    		for(int j=0; j<d; j++) {
    			for(int k=0; k<d; k++) {
    				s2[j][k] += (x[i][j] - mean[j])*(x[i][k] - mean[k]);
    			}
    		}
    	}
    	for(int j=0; j<d; j++) {
			for(int k=0; k<d; k++) {
				s2[j][k] /= (n - 1);
			}
		}
    	
    	return s2;
    }
    
    /**
     * Compute the one-dimensional dispersion about the mean.
     * 
     * @param x
     * @return
     */
    public static double getSampleDispersion(double[][] x) {
    	
    	double[] mean = getSampleMean(x);
    	
    	// Number of datapoints
    	int n = x.length;
    	// Number of dimensions
    	int d = x[0].length;
    	
    	double s2 = 0.0;
    	for(int i=0; i<x.length; i++) {
    		for(int j=0; j<d; j++) {
    			// Compute mean-subtracted scalar product
    			s2 += (x[i][j] - mean[j])*(x[i][j] - mean[j]);
    		}
    	}
    	s2 /= n;
    	
    	return s2;
    }
    
    
    /**
     * Get the mean of the multidimensional data points.
     * 
     * @param x
     * 	Two dimensional array containing the data points for which the sample covariance is to be computed.
     * The leading element indexes the data points, the trailing element indexes the dimensions of each 
     * @return
     * 	The sample mean; for N-dimensional data this array has size dimensions [N].
     */
    public static double[] getSampleMean(double[][] x) {
    	
    	// Number of datapoints
    	int n = x.length;
    	// Number of dimensions
    	int d = x[0].length;
    	
    	// Compute mean of the data points
    	double[] mean = new double[d];
    	for(int i=0; i<n; i++) {
    		for(int j=0; j<d; j++) {
    			mean[j] += x[i][j];
    		}
    	}
    	for(int j=0; j<d; j++) {
    		mean[j] /= n;
    	}
    	
    	return mean;
    }
    
    /**
     * Get an element from x at random.
     * @param x
     * 	The set of points.
     * @return
     * 	An element of x uniformly distributed within the array bounds.
     */
    public static double getRandomElement(double[] x) {
    	return x[random.nextInt(x.length)];
    }
    
    /**
     * Generates a set of points tracing out a 2D confidence ellipse of a given sigma level.
     * 
     * @param position
     * 	The mean position, i.e. the centre of the confidence ellipse.
     * @param covariance
     * 	The 2x2 covariance matrix.
     * @param sigmas
     * 	The size of the confidence ellipse [sigmas].
     * @param nPoints
     * 	The number of points to draw; these are uniformly distributed in angle, so greater distance between points
     * along the ellipse major axes.
     * @return
     */
    public static double[][] getConfidenceEllipsePoints(double[] position, Matrix covariance, double sigmas, int nPoints) {
    	
    	// Get principal axes frame (p,q) of covariance ellipse
        EigenvalueDecomposition evd = new EigenvalueDecomposition(covariance);
        
        // Length of major axes - standard deviation along these directions.
        double s_p = Math.sqrt(evd.getD().get(0, 0)) * sigmas;
        double s_q = Math.sqrt(evd.getD().get(1, 1)) * sigmas;       
        
        // Eigenvector matrix for principal axes frame
        Matrix V = evd.getV();
                
        // Ensure eigenvector matrix can be used as rotation matrix: check if
        // determinant is minus one (in which case matrix is a roto-reflection)
        // and change sign of one column if so.
        if(V.det() < 0){
            Matrix correct = new Matrix(new double[][]{{-1.0, 1.0},
                                                       {-1.0, 1.0}});
            V.arrayTimesEquals(correct);       
        }
        
        // Step in angle between points, such that the first and last point coincide
        double ang_step = 2.0*Math.PI / (nPoints-1);
        
        double[][] points = new double[nPoints][2];
        
        for(int i=0; i<nPoints; i++) {
        	
        	double ang = i * ang_step;
        
        	double c_ang = Math.cos(ang);
        	double s_ang = Math.sin(ang);
            
            // Get radius of ellipse at this angle
            double r = s_p * s_q / Math.sqrt(s_p * s_p * s_ang * s_ang + s_q * s_q * c_ang * c_ang);

            // Convert to cartesian coordinates and store as column vector.
            Matrix X = new Matrix(new double[][]{{r * c_ang},{r * s_ang}});
        
            // Rotate back to image frame
            X = V.times(X);
                        
            // Get pixel coordinates by adding mean position.
            points[i][0] = position[0] + X.get(0,0);
            points[i][1] = position[1] + X.get(1,0);
        }
        return points;
    }
    
}
