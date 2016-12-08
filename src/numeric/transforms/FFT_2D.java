package numeric.transforms;

import numeric.complex.Complex;



/**
 * Computes the 2d Fourier transform using Cooley-Tukey FFT algorithm.
 * 
 * TODO: update interfaces to operate on 2D arrays of real-only or full complex data.
 * TODO: test against external FFT routines.
 * 
 * 
 * @author nrowell
 *
 */
public class FFT_2D
{
	
	/**
	 * Compute FFT of 2D real-only data stored in 2D array.
	 * 
	 * @param data	Input real-only data, row-packed
	 * @param w		Width of 2D data array.
	 * @param h		Height of 2D array data.
	 * @return
	 */
	public static Complex[][] transformReal(double[][] real, Direction dir)
	{
		
		int w = real.length;
		int h = real[0].length;
		
		// Convert to 1D row-packed array of complex numbers by adding zero imaginary parts
		double[] complex = new double[w*h*2];
	    for(int k1=0; k1<w; k1++)
	        for(int k2=0; k2<h; k2++)
	        {
	            complex[2*(k2*w + k1) + 0] = real[w][h];
	            complex[2*(k2*w + k1) + 1] = 0.0;
	        }
	    
	    // Compute full complex 2D FFT
	    fft2d(complex, w, h, dir.sign);
	    
	    // Create output array of Complex types
	    Complex[][] out = new Complex[w][h];
	    
	    for(int k1=0; k1<w; k1++)
	        for(int k2=0; k2<h; k2++)
	        {
	        	out[k1][k2] = new Complex(complex[2*(k2*w + k1) + 0], complex[2*(k2*w + k1) + 1]);
	        }
	    
	    return out;
	}
	
	
	
	/**
	 * Compute FFT of 2D real-only data stored in 1D row-packed array.
	 * 
	 * @param data	Input real-only data, row-packed
	 * @param w		Width of 2D data array.
	 * @param h		Height of 2D array data.
	 * @return
	 */
	public static Complex[] transformReal(double[] real, int w, int h, Direction dir)
	{
		// Convert to array of complex numbers by adding zero imaginary parts
		double[] complex = new double[w*h*2];
	    for(int k1=0; k1<w; k1++)
	        for(int k2=0; k2<h; k2++)
	        {
	            complex[2*(k2*w + k1) + 0] = real[k2*w + k1];
	            complex[2*(k2*w + k1) + 1] = 0.0;
	        }
	    
	    // Compute full complex 2D FFT
	    fft2d(complex, w, h, dir.sign);
	    
	    // Create output array of Complex types
	    Complex[] out = new Complex[w*h];
	    
	    for(int k1=0; k1<w; k1++)
	        for(int k2=0; k2<h; k2++)
	        {
	        	out[k2*w + k1] = new Complex(complex[2*(k2*w + k1) + 0], complex[2*(k2*w + k1) + 1]);
	        }
	    
	    return out;
	}
	
	
	
	/**
	 * Full complex 2D FFT. Transformed function is written back to the input array.
	 * 
	 * @param input
	 * @param w
	 * @param h
	 * @param isign
	 */
	public static void fft2d(double[] input, int w, int h, double isign)
	{
		// Verify that length is an integer power of two
        if((w&(w-1)) != 0)
        {
            throw new IllegalArgumentException("Error: FFT 2D input width must be an integer power of two!\n");
        }
        if((h&(h-1)) != 0)
        {
            throw new IllegalArgumentException("Error: FFT 2D input height must be an integer power of two!\n");
        }
        
	    double[] column = new double[2*h];     // Array to store temporary column copies

	    for(int k1=0; k1<w; k1++)            // First transform columns
	    {
	        for(int k2=0; k2<h; k2++)        // Copy column to temporary array
	        {
	            column[2*k2+0] = input[2*(k2*w + k1)+0];   // Real part of [k1,k2]
	            column[2*k2+1] = input[2*(k2*w + k1)+1];   // Imaginary part of [k1,k2]
	        }
	        
	        FFT_1D.fft1d(column, isign);
	        
	        // Copy FFT of column back to input
	        for(int k2=0; k2<h; k2++)
	        {
	            input[2*(k2*w + k1)+0] = column[2*k2+0];
	            input[2*(k2*w + k1)+1] = column[2*k2+1];
	        }
	        
	    }
	    
	    // Now loop over rows. Copy out single rows, take 1D FFT then write back to the output array
	    for(int k2=0; k2<h; k2++)
	    {
	    	double[] tmp = new double[2*w];
	    	System.arraycopy(input, 2*w*k2, tmp, 0, 2*w);   	// Extract single rows from array...
	        FFT_1D.fft1d(tmp, isign);							// ...take 1D FFT...
	        System.arraycopy(tmp, 0, input, 2*w*k2, 2*w);		// ...write back to array.
	    }
	    
	    // Now shift elements to move low frequency components to centre.
	    // This amounts to swapping quadrants 1&4 and 2&3 in the array:
	    //
	    //      j=0          j=w/2           j=w-1
	    //        ______________|_____________
	    //       |             :              |  i=0
	    //       |  1      _   :   _       2  |
	    //       |        |\   :   /|         |
	    //       |          \  :  /           |
	    //       |           \ : /            |
	    //       |            \:/             |  i=h/2-1
	    //       |-------------:--------------| 
	    //       |            /:\             |  i=h/2
	    //       |           / : \            |
	    //       |          /  :  \           |
	    //       |        |/_  :  _\/         |
	    //       |  3          :           4  |
	    //       |_____________:______________|  i=h-1
	    //
	    //
	    
	    int q1,q2,q3,q4;
	    
	    double tmpi, tmpr;
	    
	    for(int i=0; i<h/2; i++)
	    {
	        for(int j=0; j<w/2; j++)
	        {
	            q1 = i*w + j;             // Index of element in quadrant 1
	            q2 = i*w + j + w/2;       //  " " 2
	            q3 = (i+h/2)*w + j;       //  " " 3
	            q4 = (i+h/2)*w + j + w/2; //  " " 4
	            
	            // Swap elements in quads 1 & 4
	            tmpr = input[2*q1 + 0];
	            tmpi = input[2*q1 + 1];
	            input[2*q1 + 0] = input[2*q4 + 0];
	            input[2*q1 + 1] = input[2*q4 + 1];
	            input[2*q4 + 0] = tmpr;
	            input[2*q4 + 1] = tmpi;
	            
	            // Swap elements in quads 2 & 3
	            tmpr = input[2*q2 + 0];
	            tmpi = input[2*q2 + 1];
	            input[2*q2 + 0] = input[2*q3 + 0];
	            input[2*q2 + 1] = input[2*q3 + 1];
	            input[2*q3 + 0] = tmpr;
	            input[2*q3 + 1] = tmpi;
	        }
	    }
	}
	
	
	
	
	
}
