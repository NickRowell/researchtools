package numeric.transforms;

import numeric.complex.Complex;

/**
 * Brute-force 1D discrete Fourier Transform that is applicable to any size of input (i.e. does not require 2^N elements).
 * The results have been validated against Matlab built in FFT routine.
 * 
 * @author nickrowell
 */
public class DFT_1D {
    /**
     * Untransformed function. Basic data type is 1D array of complex numbers.
     */
    protected Complex[] xn = null;
    
    /**
     * Transformed function. Basic data type is 1D array of complex numbers.
     */    
    protected Complex[] Xk = null;    
    
    /** Type of transform that this object performs on it's input. */
    protected Direction type;
    
    /** Number of elements. */
    protected int N;
    
    /**
     * Constructor.
     * 
     * @param   comp    2D Complex array specifying untransformed function.
     */
    public DFT_1D(Complex[] comp, Direction atype)
    {
        
        type = atype;
        
        N = comp.length;
        
        xn = comp; 
        
        Xk = new Complex[N];
        
        // Perform transformation
        transform();
    }    
    
    /** 
     * Constructor that takes an array of doubles. These are converted to
     * Complex numbers with zero imaginary part, before transforming.
     */
    public DFT_1D(double[] comp, Direction atype)
    {
        this(Complex.toComplexArray(comp), atype);
    }
    public DFT_1D(int[] comp, Direction atype)
    {
        this(Complex.toComplexArray(comp), atype);
    }
    
    /** Get results of transform. */
    public Complex[] getTransform(){ return Xk;}
    
        
    /** Perform discrete Fourier transform. */
    private void transform()
    {

        // Loop over elements in transformed function
        for(int k=0; k<N; k++)
        {
            Xk[k] = new Complex(0,0);
            
            // Loop over elements in input
            for(int n=0; n<N; n++)
            {
                Xk[k] = Xk[k].add(xn[n].times(type.twiddle(n,k,N)));
            }
            
            // Include normalisation term.
            switch(type)
            {
                case INVERSE: Xk[k].timesEquals(new Complex(1.0/N, 0));
            }       
        
        }

    }
    
}
