package images;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import images.Rendering;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

/**
 * Small class that reads an image, calculates the derivatives and writes these
 * to disk as image files.
 * @author nickrowell
 */
public class DifferentiateImage
{
    
    public static void main(String[] args) throws IOException
    {
        // Read image
        BufferedImage image = ImageIO.read(new File("/home/nickrowell/Conferences and Presentations/STAR-Seminar/images/image_representation/cropped.png"));
        
        int w = image.getWidth();
        int h = image.getHeight();
        
        // convert to array of integers
        int[][] I = Rendering.imToArray(image);
        
        // First derivative images
        int[][] Ix = new int[w][h];
        int[][] Iy = new int[w][h];
        
        // Second derivative images
        int[][] Ixx = new int[w][h];
        int[][] Iyy = new int[w][h];
        int[][] Ixy = new int[w][h];
        
        // Calculate first derivatives using Sobel kernels
        
        // Loop over columns/x
        for(int i=1; i<h-1; i++)
        {
            // loop over rows/y
            for(int j=1; j<w-1; j++)
            {
                // X gradient operation
                Ix[i][j] =   I[i-1][j-1] -   I[i+1][j-1]
                          +2*I[i-1][j]   - 2*I[i+1][j]
                            +I[i-1][j+1] -   I[i+1][j+1];
                
                // Y gradient operation
                Iy[i][j] =   I[i-1][j-1] + 2*I[i][j-1] + I[i+1][j-1]
                            -I[i-1][j+1] - 2*I[i][j+1] - I[i+1][j+1];
            }
        }
        
        // Calculate second derivatives by applying sobel kernels to first 
        // derivatives.
        
        // Loop over columns/x
        for(int i=1; i<h-1; i++)
        {
            // loop over rows/y
            for(int j=1; j<w-1; j++)
            {
                // X gradient operation applied to X derivative
                Ixx[i][j] =   Ix[i-1][j-1] -   Ix[i+1][j-1]
                           +2*Ix[i-1][j]   - 2*Ix[i+1][j]
                             +Ix[i-1][j+1] -   Ix[i+1][j+1];
                
                // Y gradient operation applied to X derivative
                Ixy[i][j] =  Ix[i-1][j-1] + 2*Ix[i][j-1] + Ix[i+1][j-1]
                            -Ix[i-1][j+1] - 2*Ix[i][j+1] - Ix[i+1][j+1];
                
                // Y gradient operation applied to Y derivative
                Iyy[i][j] =   Iy[i-1][j-1] + 2*Iy[i][j-1] + Iy[i+1][j-1]
                             -Iy[i-1][j+1] - 2*Iy[i][j+1] - Iy[i+1][j+1];
            }
        }
        
        // Now calculate max eigenvalues of Jacobian and Hessian
        float[][] J_im = new float[w][h];
        float[][] H_im = new float[w][h];
        for(int i=1; i<h-1; i++)
        {
            // loop over rows/y
            for(int j=1; j<w-1; j++)
            {
                
                Matrix J = new Matrix(new double[][]{{Ix[i][j]*Ix[i][j], Ix[i][j]*Iy[i][j]},
                                                     {Ix[i][j]*Iy[i][j], Iy[i][j]*Iy[i][j]}});
                
                Matrix H = new Matrix(new double[][]{{Ixx[i][j], Ixy[i][j]},
                                                     {Ixy[i][j], Iyy[i][j]}});
                
                EigenvalueDecomposition eigJ = new EigenvalueDecomposition(J);
                EigenvalueDecomposition eigH = new EigenvalueDecomposition(H);
                
                // Largest eigenvalues
                J_im[i][j] = Math.max(Math.abs((float)eigJ.getD().get(0,0)), Math.abs((float)eigJ.getD().get(1,1)));
                // Determinant of the Hessian
                H_im[i][j] = (float)H.det(); //   Math.max(Math.abs((float)eigH.getD().get(0,0)), Math.abs((float)eigH.getD().get(1,1)));
                
                
            }
        }
        
        
        
        
        
        ImageIO.write(arrayToIm(Ix), "png", new File("/home/nickrowell/Conferences and Presentations/STAR-Seminar/images/image_representation/cropped_Ix.png"));
        ImageIO.write(arrayToIm(Iy), "png", new File("/home/nickrowell/Conferences and Presentations/STAR-Seminar/images/image_representation/cropped_Iy.png"));
        ImageIO.write(arrayToIm(Ixx), "png", new File("/home/nickrowell/Conferences and Presentations/STAR-Seminar/images/image_representation/cropped_Ixx.png"));
        ImageIO.write(arrayToIm(Ixy), "png", new File("/home/nickrowell/Conferences and Presentations/STAR-Seminar/images/image_representation/cropped_Ixy.png"));
        ImageIO.write(arrayToIm(Iyy), "png", new File("/home/nickrowell/Conferences and Presentations/STAR-Seminar/images/image_representation/cropped_Iyy.png"));
        
        ImageIO.write(arrayToIm(J_im), "png", new File("/home/nickrowell/Conferences and Presentations/STAR-Seminar/images/image_representation/cropped_J.png"));
        ImageIO.write(arrayToIm(H_im), "png", new File("/home/nickrowell/Conferences and Presentations/STAR-Seminar/images/image_representation/cropped_H.png"));
        
        
    }
    
    
    private static void writeImage()
    {
        
    }
    
    public static BufferedImage arrayToIm(float[][] im)
    {

        // Calculate scaling and offset to apply...
        float min = im[0][0];
        float max = im[0][0];
        for (int x = 0; x < im.length; x++)
        {
            for (int y = 0; y < im[0].length; y++)
            {
                if(im[x][y] < min) min = im[x][y];
                if(im[x][y] > max) max = im[x][y];
            }
        }
        
        float range = max - min;
        
        
        BufferedImage image = new BufferedImage(im.length, im[0].length, BufferedImage.TYPE_INT_ARGB);

        // Draw integer array of pixel values into image as graylevels
        for (int x = 0; x < im.length; x++)
        {
            for (int y = 0; y < im[0].length; y++)
            {
                
                int pix = (int)(((im[x][y] - min) * 255) / range);
                
                
                int pixel = 0xFF000000 + (pix << 16) + (pix << 8) + pix;
                image.setRGB(x, y, pixel);
            }
        }

        return image;
    }
    
    public static BufferedImage arrayToIm(int[][] im)
    {

        // Calculate scaling and offset to apply...
        int min = im[0][0];
        int max = im[0][0];
        for (int x = 0; x < im.length; x++)
        {
            for (int y = 0; y < im[0].length; y++)
            {
                if(im[x][y] < min) min = im[x][y];
                if(im[x][y] > max) max = im[x][y];
            }
        }
        
        int range = max - min;
        
        
        BufferedImage image = new BufferedImage(im.length, im[0].length, BufferedImage.TYPE_INT_ARGB);

        // Draw integer array of pixel values into image as graylevels
        for (int x = 0; x < im.length; x++)
        {
            for (int y = 0; y < im[0].length; y++)
            {
                
                int pix = ((im[x][y] - min) * 255) / range;
                
                int pixel = 0xFF000000 + (pix << 16) + (pix << 8) + pix;
                image.setRGB(x, y, pixel);
            }
        }

        return image;
    }
}
