package images;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * This class was hacked up to convert an IMG file representing a colour image
 * from sequential bands (where all red values are listed first, then all green,
 * then all blue) to interleaved (where each consecutive 3 bytes represents
 * red, green and blue intensities in a single pixel).
 * @author nickrowell
 */
public class SequentialToInterleaved
{
    
    public static void main(String[] args) throws IOException
    {
        
        File in  = new File("/home/nickrowell/Projects/PANGU_4/Documents/TN01/Planetary_Science/Available_Data/Shape_Models/Mars/HiRISE/Gale_Crater/mrohr_0001_0001/dtm/psp/orb_010500_010599/psp_010573_1755_psp_010639_1755/pangu_models/correct_PXN/PSP_010573_1755_IRB_A_01_ORTHO.SEQUENTIAL");
        File out = new File("/home/nickrowell/Projects/PANGU_4/Documents/TN01/Planetary_Science/Available_Data/Shape_Models/Mars/HiRISE/Gale_Crater/mrohr_0001_0001/dtm/psp/orb_010500_010599/psp_010573_1755_psp_010639_1755/pangu_models/correct_PXN/PSP_010573_1755_IRB_A_01_ORTHO.INTERLEAVED");
        FileInputStream freader = new FileInputStream(in);
        
        // Size of data and number of bands
        int w = 2049;
        int h = 2049;
        
        int nbytes = w*h*3;
        
        // Store all data temporarily
        int[] data = new int[nbytes];
        
        int bytes_so_far = 0;
        while(bytes_so_far < nbytes)
        {
            data[bytes_so_far++] = freader.read();
        }
        
        // data array is now loaded with sequential data
        
        FileOutputStream fwriter = new FileOutputStream(out);
        
        // Loop over each pixel
        for(int p=0; p<w*h; p++)
        {
            int b1 = data[p+(0*w*h)];  // Band 1  (near-IR)
            int b2 = data[p+(1*w*h)];  // Band 2  (red)
            int b3 = data[p+(2*w*h)];  // Band 3  (blue-green)
            
            // Transform to standard RGB bands
            // See http://www.uahirise.org/pdf/color-products.pdf
            int r = b2;
            int g = b3;
            int b = 2*g - (int)(0.3*r);  // Range doubled; need to halve again
            
            fwriter.write(r);
            fwriter.write(g);
            fwriter.write(b/2);
            
        }
        
        fwriter.flush();
        
    }
    
    
    
    
}