package numeric.transforms;

import numeric.complex.Complex;

/**
 * Standard Cooley-Tukey decimation-in-time FFT algorithm (See e.g. Numerical Recipes 12.2.1).
 * Requires 2^N input elements.
 * 
 * @author nickrowell
 */
public class FFT_1D
{
    
    /** 
     * Alternative interface for real-only input.
     */
    public static Complex[] transformReal(double[] real, Direction type)
    {
    	// Pad out array with zero imaginary parts
    	double[] comp = new double[real.length*2];
    	for(int i=0; i<real.length; i++)
    	{
    		comp[2*i + 0] = real[i];
    		comp[2*i + 1] = 0.0;
    	}
    	
    	// Perform FFT
    	fft1d(comp, type.sign);
    	
    	// Extract results to Complex array
    	Complex[] output = new Complex[real.length];
    	for(int i=0; i<real.length; i++)
    	{
    		output[i] = new Complex(comp[2*i + 0], comp[2*i + 1]);
    	}
    	
    	return output;
    }
    
    
    public static Complex[] transformComplex(Complex[] input, Direction type)
    {
    	// Extract Complex values to double array
    	double[] comp_arr = new double[input.length*2];
    	
    	for(int i=0; i<input.length; i++)
    	{
    		comp_arr[2*i + 0] = input[i].getRe();
    		comp_arr[2*i + 1] = input[i].getIm();
    	}
    	
    	// Perform FFT
    	fft1d(comp_arr, type.sign);
    	
    	// Extract results to Complex array
    	Complex[] output = new Complex[input.length];
    	for(int i=0; i<input.length; i++)
    	{
    		output[i] = new Complex(comp_arr[2*i + 0], comp_arr[2*i + 1]);
    	}
    	
    	return output;
    }
    
    
    
    public static Complex[] transformComplex(double[] comp, Direction type)
    {
    	// Perform FFT
    	fft1d(comp, type.sign);
    	
    	// Extract results to Complex array
    	Complex[] output = new Complex[comp.length/2];
    	for(int i=0; i<comp.length/2; i++)
    	{
    		output[i] = new Complex(comp[2*i + 0], comp[2*i + 1]);
    	}
    	
    	return output;
    }
    
    
    
    
    
    
    
    /**
     * Main FFT interface. Uses a standard Cooley-Tukey decimation-in-time algorithm (See e.g. Numerical Recipes 12.2.1).
     * The transformed function is written back to the input array.
     * 
     * @param data	Untransformed function: array contains pairs of data corresponding to real and imaginary parts.
     * 				The number of data points MUST be equal to a power of 2 (otherwise FFT doesn't apply), so the
     * 				array must have a length equal to twice a power of two.
     * @param type	Direction of transformation.
     * 
     */
    public static void fft1d(double[] data, double isign)
    {
    	// Number of data points
        int n = data.length/2;
        
        // Verify that length is an integer power of two
        if((n&(n-1)) != 0)
        {
            throw new IllegalArgumentException("Error: FFT input must have a length an integer power of two!\n");
        }
        
        
        int nn  = n << 1;
        int mmax, istep, m, i, j=1;

        // Bit-reverse input data by swapping corresponding pairs.
        for(i=1; i<nn; i+=2)
        {
            if(j>i)
            {
                swap(data, j-1, i-1);
                swap(data, j, i);
            }
            m=n;
            while(m >= 2 && j > m)
            {
                j -= m;
                m >>= 1;
            }
            j += m;
        }
        
        // Begin main FFT section

        double theta, wtemp, wpr, wpi, wr, wi, tempr, tempi;

        mmax=2;
        while(nn > mmax)
        {
            istep = mmax << 1;
            theta = isign*(2*Math.PI/mmax);        // Twiddle factor
            wtemp = Math.sin(0.5*theta);
            wpr   = -2.0*wtemp*wtemp;
            wpi   = Math.sin(theta);
            wr    = 1.0;
            wi    = 0.0;

            for(m=1; m<mmax; m+=2)
            {
                for(i=m; i<=nn; i+=istep)
                {
                    j = i+mmax;
                    tempr = wr*data[j-1] - wi*data[j];
                    tempi = wr*data[j]   + wi*data[j-1];
                    data[j-1] = data[i-1]-tempr;
                    data[j]   = data[i] - tempi;
                    data[i-1] += tempr;
                    data[i]   += tempi;
                }
                wr = (wtemp=wr)*wpr - wi*wpi + wr;
                wi = wi*wpr + wtemp*wpi + wi;
            }
            mmax=istep;
        }
        
    }
    
    
    /**
     * Swaps two values in an array.
     * @param in	Array containing values to be swapped.
     * @param a		Index of first value.
     * @param b		Index of second value.
     */
    private static void swap(double[] in, int a, int b)
    {
    	double tmp = in[a];
    	in[a] = in[b];
    	in[b] = tmp;
    }
    
}
