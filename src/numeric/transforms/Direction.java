package numeric.transforms;

import numeric.complex.Complex;

public enum Direction {
	
	
	FORWARD (-1),
    INVERSE (1);

    // Sign of kernel
    public final double sign;
    
    Direction(double psign){ sign = psign;}
	
    /** Get twiddle factor. */
    public Complex twiddle(double n, double k, double N){
        double phase = sign * 2.0 * Math.PI * k * n / N;
        Complex twiddle = new Complex(Math.cos(phase), Math.sin(phase));
        return twiddle;    
    }
}
