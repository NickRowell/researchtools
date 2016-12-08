package images;

import images.PPM;
import images.Rendering;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

/**
 * This class is intended for testing PPM reading functions of PPM class.
 * @author nickrowell
 */
public class PPMTester
{
    
    
    public static void main(String args[]) throws IOException
    {
    
    
        // Image file
        File image = new File("/home/nickrowell/Desktop/disk.ppm");
    
        // Read to Java BufferedImage
        BufferedImage im_ppm = PPM.readPPM(image);
        
        // Convert to array of integers
        int[][] im_arr = Rendering.imToArray(im_ppm);
        
        // Sum intensity of all pixels
        double int_flux = 0;
        
        for(int[] row : im_arr)
            for( int pixel : row)
                int_flux += pixel;
        
        System.out.println("Sum of all pixels = "+int_flux);
    
    }
    
    
    
}
