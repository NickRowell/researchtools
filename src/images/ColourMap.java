/**
 * Copyright (c) 2011 University of Dundee
 *
 * Name:
 * ColourMap.java
 *
 * Purpose:
 * ColourMaps map number line to RGB colour. They are used for e.g. mapping
 * age of a tracked feature point to a unique colour for rendering in an image.
 * 
 * Language:
 * Java
 *
 * Author:
 * Nicholas Rowell (nickrowell@computing.dundee.ac.uk)
 *
 */

package images;

public class ColourMap {

    /** Max feature point age. */
    int max;

    /** Map of age/RGB component */
    int[][] rgb;


    /**
     * Colour map that sets the colour of a feature track
     * depending on its age. The operation is functionally the same
     * as the create_colormap() method used in glfeic. The
     * comments are Martin's.
     */
    public ColourMap(int n){

        rgb = new int[n][3];
        max = n;

	// We create a rainbow colour map which fades from pure red
	// through yellow, green, cyan, blue and magenta. That means
	// we have five equal gradients: divide n by 5 and find out
	// the spare slots. Unfortunately that maximum value of n is
	// 1274 which means a 10-bit colour map (1024).
	int n_div_5 = n / 5;
	int spare   = n % 5;

	// Each gradient will have n_div_5 steps to range from values
	// of between 0 and 255. Compute the maximum integer step and
	// the extra bit required to reach the maximum value 255.
	int grad_step = 255 / n_div_5;
	int grad_xtra = 255 % n_div_5;

	// Gradient step counter.
	int i;

        // Element index counter
        int e=0;

	// Temporary RGB values.
	double r, g, b;

	// First gradient: red to yellow.
	for (i = 0, r = 1.0, b = 0.0; i < n_div_5; i++)
	{
		g = (float)(i*grad_step + grad_xtra)/255.0;

                rgb[e][0]   = (int)(r*255);
                rgb[e][1]   = (int)(g*255);
                rgb[e++][2] = (int)(b*255);
	}

	// Second gradient: yellow to green.
	for (i = 0, g = 1.0, b = 0.0; i < n_div_5; i++)
	{
		int j = n_div_5 - i - 1; // Counting down

		r = (float)(j*grad_step + grad_xtra)/255.0;

                rgb[e][0]   = (int)(r*255);
                rgb[e][1]   = (int)(g*255);
                rgb[e++][2] = (int)(b*255);
	}

	// Third gradient: green to cyan.
	for (i = 0, r = 0.0, g = 1.0; i < n_div_5; i++)
	{
		b = (float)(i*grad_step + grad_xtra)/255.0;
		rgb[e][0]   = (int)(r*255);
                rgb[e][1]   = (int)(g*255);
                rgb[e++][2] = (int)(b*255);
	}

	// Fourth gradient: cyan to blue.
	for (i = 0, r = 0.0, b = 1.0; i < n_div_5; i++)
	{
		int j = n_div_5 - i - 1; // Counting down
		g = (float)(j*grad_step + grad_xtra)/255.0;

                rgb[e][0]   = (int)(r*255);
                rgb[e][1]   = (int)(g*255);
                rgb[e++][2] = (int)(b*255);
	}

	// Final gradient: blue to magenta.
	for (i = 0, g = 0.0, b = 1.0; i < n_div_5; i++)
	{
		r = (float)(i*grad_step + grad_xtra)/255.0;

                rgb[e][0]   = (int)(r*255);
                rgb[e][1]   = (int)(g*255);
                rgb[e++][2] = (int)(b*255);
	}

	// Fill in any remaining steps with magenta.
	for (i = 0, r = 1.0, g = 0.0, b = 1.0; i < spare; i++){
                rgb[e][0]   = (int)(r*255);
                rgb[e][1]   = (int)(g*255);
                rgb[e++][2] = (int)(b*255);
        }

    }


    public int getColour(int frames){

        if(frames >= max) frames = max-1;

        return (rgb[frames][0] << 16) + (rgb[frames][1] << 8) + rgb[frames][2];

    }

}
