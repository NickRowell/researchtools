/**
 * 
 * Name:
 *  ImageGobbler.java
 * 
 * Purpose:
 *  Class used to read images back from Gnuplot output stream.
 * 
 * Language:
 * Java
 *
 * Author:
 * Nicholas Rowell
 * 
 */
package infra.io;

import java.awt.image.BufferedImage;
import java.io.IOException;
import java.io.InputStream;
import javax.imageio.ImageIO;

public class ImageGobbler 
extends Thread
{
    InputStream is;
    public BufferedImage img;

    public ImageGobbler(InputStream is) 
    {
        this.is = is;
    }

    @Override
    public void run() 
    {
        try 
        {
            img = ImageIO.read(is);
        } 
        catch (IOException ex)
        {
            System.err.println("Unable to create image: "+ex.getMessage());
        }
//        if (img==null) 
//            System.out.println("Unable to create image from the gnuplot "
//                    + "output. Null image created.");
        
    }
}
