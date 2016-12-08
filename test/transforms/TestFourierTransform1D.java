package transforms;

import transforms.Direction;
import transforms.FFT_1D;
import transforms.DFT_1D;
import misc.Complex;

/**
 * This class is designed to test FourierTransform1D class.
 * @author nickrowell
 */
public class TestFourierTransform1D {
    
    
    
    public static void main(String[] args){
    
        // Generic values
//        double[] xn = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
        
        
        Complex[] xn = new Complex[]
        		{
        			new Complex(0.0,1.0),
        			new Complex(1.0,4.0),
        			new Complex(3.0,8.0),
        			new Complex(6.0,2.0),
        			new Complex(7.0,6.0),
        			new Complex(4.0,4.0),
        			new Complex(3.0,7.0),
        			new Complex(8.0,9.0),
        			new Complex(0.0,2.0),
        			new Complex(3.0,-5.0),
        			new Complex(8.0,6.0),
        			new Complex(9.0,2.0),
        			new Complex(1.0,-3.0),
        			new Complex(1.0,-1.0),
        			new Complex(2.0,0.0),
        			new Complex(8.0,3.0)
        		};
        
        

//        // Sine wave
//        double[] xn = new double[128];
//        // frequency in units of cycles per array element
//        double freq = 0.1;
//        for(int i=0; i<xn.length; i++)
//            xn[i] = Math.sin(2 * Math.PI * freq * i);
        
        
        
        
        DFT_1D dft = new DFT_1D(xn, Direction.FORWARD);
        Complex[] Xk = dft.getTransform();
        
        
        
        Complex[] Xk_fft = FFT_1D.transformComplex(xn, Direction.FORWARD);
        
        
        System.out.println("Input data:");
        System.out.print("xn = [");
        for(int i=0; i<xn.length; i++)
        {
        	System.out.println(xn[i].toString() + " ; ...");
        }
        System.out.println("]\n");
        
        
        
        
        System.out.println("Output from Fourier Transform algorithms");
        System.out.println("----------------------------------------");
        
        System.out.println("Element			DFT			FFT");
        
        for(int i=0; i<xn.length; i++)
        {
        	System.out.println(i+"\t\t"+Xk[i].toString() + "\t" + Xk_fft[i].toString());
        }
        
        
        System.out.println("\nElement		DFT	(mag)		FFT (mag)");
        
        for(int i=0; i<xn.length; i++)
        {
        	System.out.println(i+"\t\t"+Xk[i].abs() + "\t" + Xk_fft[i].abs());
        }
        
        
        
    }
    
    
    
}
