package numeric.geom.dim2.test;

import java.awt.BorderLayout;
import java.awt.image.BufferedImage;

import javax.swing.JFrame;

import images.Rendering;
import infra.gui.IPanel;
import numeric.geom.dim2.Triangle2d;
import numeric.geom.dim2.Vector2d;

/**
 * Test class for the {@link Triangle2d}.
 *
 *
 * @author nrowell
 * @version $Id$
 */
public class TriangleTest {
	
	/**
	 * Main application entry point.
	 * 
	 * @param args
	 * 	The command line arguments.
	 */
	public static void main(String[] args) {
		
		// Box within which to randomly distribute the points
		double xmin = 0;
		double xmax = 512;
		double ymin = 0;
		double ymax = 512;
		
		// Size of square image to render
		final int width  = 512;
		final int height = 512;
		
		// Randomly generate three points
		Vector2d v0 = new Vector2d(xmin + Math.random()*(xmax - xmin), ymin + Math.random()*(ymax - ymin));
		Vector2d v1 = new Vector2d(xmin + Math.random()*(xmax - xmin), ymin + Math.random()*(ymax - ymin));
		Vector2d v2 = new Vector2d(xmin + Math.random()*(xmax - xmin), ymin + Math.random()*(ymax - ymin));
		
		// Create a Triangle
		Triangle2d tri = new Triangle2d(v0, v1, v2);
		
		// Render an image of the result
		BufferedImage im = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		
		
		// Colour all the pixels according to whether they are inside the triangle or not
		int inside = 0xFF00000;
		int outside = 0x00FF00;
		for(double x=xmin; x<xmax; x++) {
			for(double y=ymin; y<ymax; y++) {
				
				if(tri.contains(new Vector2d(x, y))) {
					im.setRGB((int) Math.rint(x), (int) Math.rint(y), inside);
				}
				else {
					im.setRGB((int) Math.rint(x), (int) Math.rint(y), outside);
				}
				
				
			}
		}
		
		int colour = 0x0000FF;
		
		Rendering.drawLine(im, new int[]{(int)v0.getX(), (int)v0.getY()}, new int[]{(int)v1.getX(), (int)v1.getY()}, colour);
		Rendering.drawLine(im, new int[]{(int)v0.getX(), (int)v0.getY()}, new int[]{(int)v2.getX(), (int)v2.getY()}, colour);
		Rendering.drawLine(im, new int[]{(int)v1.getX(), (int)v1.getY()}, new int[]{(int)v2.getX(), (int)v2.getY()}, colour);

		JFrame frame = new JFrame("TriangleTest");
		frame.setLayout(new BorderLayout());
		
		IPanel ipanel = new IPanel();
		ipanel.setImage(im);
		frame.getContentPane().add(ipanel, BorderLayout.CENTER);
		frame.pack();
		frame.setVisible(true);
		
	}
	
}