
package images;

import Jama.Matrix;
import java.awt.image.BufferedImage;
import java.awt.image.Raster;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;

/**
 * Utility class containing a bunch of handy functions for reading, converting,
 * manipulating etc images.
 */
public class ImageUtil 
{
    
    /** Read a BufferedImage into an array of integers. */
    public static int[][] getImage(BufferedImage im){
    
        Raster in = im.getRaster();

        int WIDTH = im.getWidth();
        int HEIGHT = im.getHeight();

        int[][] image = new int[WIDTH][HEIGHT];

        for (int x = 0; x < WIDTH; x++)
            for (int y = 0; y < HEIGHT; y++)
                image[x][y] = in.getSample(x, y, 0);        
     
        return image; 
    
    }
       
    /**
     * Read a .ppm image from a URL into a BufferedImage.
     */
    public static BufferedImage readPPM(URL url) throws IOException {
        
        // Byte stream to store byte chunks read from URL target
        ByteArrayOutputStream bais = new ByteArrayOutputStream();

        // Open stream reader on URL
        InputStream in = url.openStream();
        
        // Byte array to store chunks read from URL and transfer to byte stream
        byte[] buf = new byte[512];
        
        int len;
        
        while (true) {
            len = in.read(buf);
            if (len == -1) {
                break;
            }
            bais.write(buf, 0, len);
        }
        
        return readPPM(bais.toByteArray());

    } 

    /**
     * Read a .ppm image from a byte array into a BufferedImage.
     */
    public static BufferedImage readPPM(byte[] ppm) throws IOException {

        String width = "";
        String height = "";
        String maxPixVal = "";

        // Check image signature
        if(ppm[0] != 'P' || ppm[1] != '6') 
            throw new IOException("Image not a RAW ppm!");

        // Position pointer after signature
        int p = 2;

        // Read white space until values in the range 48-57 are found 
        // (indicates numerals in ASCII)
        while ((ppm[p] < 48) || (ppm[p] > 57)) p++;

        // Read in width
        for( ; (ppm[p] > 47) && (ppm[p] < 58); p++) width = width + (char) ppm[p];

        // Read white space separating width and height
        while ((ppm[p] < 48) || (ppm[p] > 57)) p++;

        // Read in height
        for( ; (ppm[p] > 47) && (ppm[p] < 58); p++) height = height + (char) ppm[p];

        // Once this code is reached, width and height of image have been read.
        // Some images have an extra number indicating max allowable pixel 
        // value, some don't. Read remaining line until new line character is 
        // found, and if another numeral is found before then, read it in.

        while (ppm[p] != 10) {
            if ((ppm[p] > 47) && (ppm[p] < 58)) {
                maxPixVal = maxPixVal + (char) ppm[p];
            }
            p++;
        }

        // Move pointer onto first datapoint
        p++;

        // Image header has been read and pointer is positioned at start of 
        // data section.

        int w = Integer.parseInt(width);
        int h = Integer.parseInt(height);

        BufferedImage out = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);

        // Read data section three bytes at a time until EOF character found
        int r, g, b;

        // increment y coordinate (y position or row index)
        for (int j = 0; j < h; j++) {

            // Scan along this row (x position or column index)
            for (int i = 0; i < w; i++) {

                r = ppm[p++];
                g = ppm[p++];
                b = ppm[p++];

                int pixel = 0x000000 + (r << 16) + (g << 8) + b;
                out.setRGB(i, j, pixel);
            }
        }
        return out;
    }
    
    /**
     * Determines index of start of data in array.
     * @param pgm
     * @return 
     */
    public static int getPGMHeaderSize(byte[] pgm)
    {
        // Check image signature
        if(pgm[0] != 'P' || pgm[1] != '5') 
            throw new RuntimeException("Image not a RAW pgm!");

        // Position pointer after signature
        int p = 2;

        // Read white space until values in the range 48-57 are found 
        // (indicates numerals in ASCII)
        while ((pgm[p] < 48) || (pgm[p] > 57)) p++;

        // Read in width
        for( ; (pgm[p] > 47) && (pgm[p] < 58); p++) ;

        // Read white space separating width and height
        while ((pgm[p] < 48) || (pgm[p] > 57)) p++;

        // Read in height
        for( ; (pgm[p] > 47) && (pgm[p] < 58); p++) ;

        // Once this code is reached, width and height of image have been read.
        // Some images have an extra number indicating max allowable pixel 
        // value, some don't. Read remaining line until new line character is 
        // found, and if another numeral is found before then, read it in.

        while (pgm[p] != 10) { p++;}

        // Move pointer onto first datapoint
        p++;
        
        return p;
    }
    
    /**
     * Read PGM header and return image width and height in int array, as well
     * as number of bytes in header.
     * @param pgm
     * @return 
     */
    public static int[] getPGMHeader(byte[] pgm) throws IOException
    {
        int[] out = new int[3];
        
        String width = "";
        String height = "";
        String maxPixVal = "";

        // Check image signature
        if(pgm[0] != 'P' || pgm[1] != '5') 
            throw new IOException("Image not a RAW pgm!");

        // Position pointer after signature
        int p = 2;

        // Read white space until values in the range 48-57 are found 
        // (indicates numerals in ASCII)
        while ((pgm[p] < 48) || (pgm[p] > 57)) p++;

        // Read in width
        for( ; (pgm[p] > 47) && (pgm[p] < 58); p++) width = width + (char) pgm[p];

        // Read white space separating width and height
        while ((pgm[p] < 48) || (pgm[p] > 57)) p++;

        // Read in height
        for( ; (pgm[p] > 47) && (pgm[p] < 58); p++) height = height + (char) pgm[p];

        // Once this code is reached, width and height of image have been read.
        // Some images have an extra number indicating max allowable pixel 
        // value, some don't. Read remaining line until new line character is 
        // found, and if another numeral is found before then, read it in.

        while (pgm[p] != 10) {
            if ((pgm[p] > 47) && (pgm[p] < 58)) {
                maxPixVal = maxPixVal + (char) pgm[p];
            }
            p++;
        }

        // Move pointer onto first datapoint
        p++;

        // Image header has been read and pointer is positioned at start of 
        // data section.
        out[0] = Integer.parseInt(width);
        out[1] = Integer.parseInt(height);
        out[2] = p;
        
        return out;
        
    }
    
    /**
     * Remove header from PPM image to leave a packed array of RGB pixels.
     * @param ppm
     * @return 
     */
    public static byte[] toRGBArray(byte[] ppm)
    {
        String width = "";
        String height = "";
        String maxPixVal = "";

        // Check image signature
        if(ppm[0] != 'P' || ppm[1] != '6') 
            throw new RuntimeException("Image not a RAW ppm!");

        // Position pointer after signature
        int p = 2;

        // Read white space until values in the range 48-57 are found 
        // (indicates numerals in ASCII)
        while ((ppm[p] < 48) || (ppm[p] > 57)) p++;

        // Read in width
        for( ; (ppm[p] > 47) && (ppm[p] < 58); p++) width = width + (char) ppm[p];

        // Read white space separating width and height
        while ((ppm[p] < 48) || (ppm[p] > 57)) p++;

        // Read in height
        for( ; (ppm[p] > 47) && (ppm[p] < 58); p++) height = height + (char) ppm[p];

        // Once this code is reached, width and height of image have been read.
        // Some images have an extra number indicating max allowable pixel 
        // value, some don't. Read remaining line until new line character is 
        // found, and if another numeral is found before then, read it in.

        while (ppm[p] != 10) {
            if ((ppm[p] > 47) && (ppm[p] < 58)) {
                maxPixVal = maxPixVal + (char) ppm[p];
            }
            p++;
        }

        // Move pointer onto first datapoint
        p++;

        // Image header has been read and pointer is positioned at start of 
        // data section.

        int w = Integer.parseInt(width);
        int h = Integer.parseInt(height);

        byte[] out = new byte[w*h*3];

        for(int pp=0; pp<w*h; pp++)
        {
            out[pp*3+0] = ppm[p+pp*3+0]; // R
            out[pp*3+1] = ppm[p+pp*3+1]; // G
            out[pp*3+2] = ppm[p+pp*3+2]; // B
        }

        return out;
    }
    
    
    /**
     * Print the BufferedImage type.
     * 
     * @param image 
     */
    public static void printBufferedImageType(BufferedImage image)
    {
        
        switch(image.getType())
            {
                case BufferedImage.TYPE_3BYTE_BGR: System.out.println("TYPE_3BYTE_BGR"); break;
                case BufferedImage.TYPE_4BYTE_ABGR: System.out.println("TYPE_4BYTE_ABGR");  break;   
                case BufferedImage.TYPE_4BYTE_ABGR_PRE: System.out.println("TYPE_4BYTE_ABGR_PRE");  break;   
                case BufferedImage.TYPE_BYTE_BINARY: System.out.println("TYPE_BYTE_BINARY");     break;
                case BufferedImage.TYPE_BYTE_GRAY: System.out.println("TYPE_BYTE_GRAY");     break;
                case BufferedImage.TYPE_BYTE_INDEXED: System.out.println("TYPE_BYTE_INDEXED");  break;   
                case BufferedImage.TYPE_CUSTOM: System.out.println("TYPE_CUSTOM");     break;
                case BufferedImage.TYPE_INT_ARGB: System.out.println("TYPE_INT_ARGB");     break;
                case BufferedImage.TYPE_INT_ARGB_PRE: System.out.println("TYPE_INT_ARGB_PRE");  break;   
                case BufferedImage.TYPE_INT_BGR: System.out.println("TYPE_INT_BGR");     break;
                case BufferedImage.TYPE_INT_RGB: System.out.println("TYPE_INT_RGB");     break;
                case BufferedImage.TYPE_USHORT_555_RGB: System.out.println("TYPE_USHORT_555_RGB"); break;    
                case BufferedImage.TYPE_USHORT_565_RGB: System.out.println("TYPE_USHORT_565_RGB"); break;
                case BufferedImage.TYPE_USHORT_GRAY: System.out.println("TYPE_USHORT_GRAY"); break;
            
            }
    
    }
    
    
    
}
