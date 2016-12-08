package numeric.interpolation.arrays;

/**
 * Defines the interface that ArrayInterpolators must implement.
 * 
 * ArrayInterpolators provides means to subsample arrays and histograms to obtain higher resolution
 * approximations that have certain mathematical properties. The main property is that if you integrate
 * the subsampled array elements over the extent of each original element then you get the same value,
 * i.e. the subsampled function conserves the contents of the original bins so that if you then
 * supersample it you get the original content back. Additional properties of the spline variants
 * are varying degrees of smoothness in the interpolating function to better approximate the true 
 * underlying function while reducing oscillations. This is different to standard interpolation methods
 * where the function is constrained to pass precisely through each data point, with no constraint on
 * the integral of the function in the vicinity of each point.
 * 
 * This is especially useful for such applications as interpolating line spread functions, luminosity
 * functions, and any other case where the array elements represent binned data at regular intervals,
 * and quantify the value of the underlying function integrated over a finite range rather than at
 * discrete points. A valid subsampling the elements (e.g. to ease deconvolution) must conserve the
 * integral of the function over the range of the original points.
 *
 *
 * @author nrowell
 * @version $Id$
 */
public abstract class ArrayInterpolator {

	/**
	 * Number of elements in the array.
	 */
	int N;
	
	/**
	 * Get the interpolated value at position x in the array. Element zero extends
	 * from [0:1], element 1 from [1:2] etc.
	 * 
	 * @param x
	 * @return
	 */
	public abstract double getValue(double x);
	
	/**
	 * Subsamples the input array elements according to the rules of the interpolation.
	 * 
	 * @param input
	 * @param subSampleFactor
	 * @return
	 */
	public double[] getSubSampledArray(int subSampleFactor) {

		if(subSampleFactor < 1) {
			throw new IllegalArgumentException("Sub-sample factor must be equal to or greater than one!");
		}
		
		// number of elements in sub-sampled array
		final int nSubSamples = N * subSampleFactor;
		
		// Empty sub-sampled array
		final double[] subSamples = new double[nSubSamples];
		
		// Fill out sub-sampled array
		for(int i=0; i<nSubSamples; i++) {
			
			// x coordinate of this subsample
			double x = (double)i/(double)subSampleFactor;
			
			// Interpolated value
			subSamples[i] = getValue(x)/subSampleFactor;
		}
		
		return subSamples;
	}
	
}
