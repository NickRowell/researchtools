package numeric.triangulation.dim2.exec;

import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import javax.imageio.ImageIO;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JPanel;

import infra.gui.IPanel;
import numeric.geom.dim2.Vector2d;
import numeric.triangulation.dim2.base.DelaunayTriangulation2D;
import numeric.triangulation.dim2.base.VoronoiTessellation2D;
import numeric.triangulation.dim2.impl.DelaunayIncremental;

/**
 * Class provides an application used to test the 2D Triangulation.
 *
 * @author nrowell
 * @version $Id$
 */
public class StainedGlassDesigner extends JPanel {
	
	/**
	 * The serial version UID.
	 */
	private static final long serialVersionUID = 6627096547262250370L;

	/**
	 * Output directory for images.
	 */
	private static File outputDir = new File("/home/nrowell/Temp/tmp/stained_glass/");
	
	int idx=0;
	Random rnd = new Random();
	
	// Size of square image to save to disk
	static final int width  = 250;
	static final int height = 525;
	
	// Number of points to generate for RANDOM
	static int n = 8;
	
	// Box within which to randomly distribute the points
	static double xmin = 0;
	static double xmax = width;
	
	// Number of vertices to generate, in X direction
	static int xVerts = 4;
	// Step size between point centres in X direction
	static double xstep = xmax/(xVerts-1);
	
	
	static double ymin = 0;
	static double ymax = height;
	
	// Number of interior vertices to generate, in Y direction
	static int yVerts = 6;
	// Step size between point centres in Y direction
	static double ystep = ymax/(yVerts-1);
	
	// Size of random displacement
	static double sigma_x = 15.0;
	static double sigma_y = 30.0;
//	static double sigma_x = 1.0;
//	static double sigma_y = 1.0;
	
	/**
	 * Colours used to paint the Delaunay Triangulation.
	 */
	int[] delTriColours = new int[]{0xFFD333, 0x63eF33, 0xa333FF, 0x883388};

	/**
	 * Colours used to paint the Voronoi Tessellation.
	 */
	int[] vorTessColours = new int[]{0xFF3333, 0x33FF33, 0x3333FF, 0x883388, 0x889388, 0xa83388, 0x883317, 0x813338, 0x58b386, 0x49e6da};

	/**
	 * Options for ways distribute the points.
	 */
	public static enum Style {RANDOM, SEMI_RANDOM};
	
	/**
	 * Currently selected {@link Style}
	 */
	Style style;
	
	/**
	 * Current Voronoi Tessellation coloured image
	 */
	BufferedImage tessColouredIm;
	
	/**
	 * Current Voronoi Tessellation wireframe image
	 */
	BufferedImage tessWireframeIm;
	
	/**
	 * Current Delaunay Triangulation coloured image
	 */
	BufferedImage delColouredIm;

	/**
	 * Current Delaunay Triangulation wireframe image
	 */
	BufferedImage delWireframeIm;
	
	/**
	 * Current Delaunay Triangulation
	 */
	DelaunayTriangulation2D delTri;
	
	/**
	 * Current Voronoi Tessellation
	 */
	VoronoiTessellation2D vorTess;
	
	/**
	 * Main constructor for the {@link StainedGlassDesigner}.
	 */
	public StainedGlassDesigner() {
		
		this.setLayout(new BorderLayout());
		
		style = Style.SEMI_RANDOM;
		
		// Initialise the plots
		redraw();
		
		// Central panel showing the 4 images
		final IPanel tl = new IPanel(delColouredIm);
		final IPanel tr = new IPanel(delWireframeIm);
		final IPanel bl = new IPanel(tessColouredIm);
		final IPanel br = new IPanel(tessWireframeIm);
		
		JPanel plotPanel = new JPanel(new GridLayout(1,4));
		plotPanel.add(tl);
		plotPanel.add(tr);
		plotPanel.add(bl);
		plotPanel.add(br);
		this.add(plotPanel, BorderLayout.CENTER);
		
		JButton saveButton = new JButton("Save Images");
		JButton redrawButton = new JButton("Redraw");
		
		final JComboBox<Style> styleComboBox = new JComboBox<Style>(Style.values());
		styleComboBox.setSelectedItem(style);
		
		JPanel buttonPanel = new JPanel(new GridLayout(1,3));
		buttonPanel.add(saveButton);
		buttonPanel.add(styleComboBox);
		buttonPanel.add(redrawButton);
		this.add(buttonPanel, BorderLayout.SOUTH);
		
		saveButton.addActionListener(new ActionListener() {
    		@Override
            public void actionPerformed(ActionEvent evt) {
                try {
					ImageIO.write(delColouredIm, "png", new File(outputDir, "delaunay_coloured_"+idx+".png"));
					ImageIO.write(delWireframeIm, "png", new File(outputDir, "delaunay_wireframe_"+idx+".png"));
					ImageIO.write(tessColouredIm, "png", new File(outputDir, "voronoi_coloured_"+idx+".png"));
					ImageIO.write(tessWireframeIm, "png", new File(outputDir, "voronoi_wireframe_"+idx+".png"));
					idx++;
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
            }
        });
		
		styleComboBox.addActionListener(new ActionListener() {
    		@Override
            public void actionPerformed(ActionEvent evt) {
    			style = (Style)styleComboBox.getSelectedItem();
            }
        });
		
		redrawButton.addActionListener(new ActionListener() {
    		@Override
            public void actionPerformed(ActionEvent evt) {
    			redraw();
    			tl.setImage(delColouredIm);
    			tr.setImage(delWireframeIm);
    			bl.setImage(tessColouredIm);
    			br.setImage(tessWireframeIm);
            }
        });
	}
	
	/**
	 * Randomly distributes points, then triangulates them, creates images and updates
	 * the internal image fields ready for display in the GUI.
	 */
	private void redraw() {
		
		List<Vector2d> points = new LinkedList<>();
		
		// Put vertices on the boundary
		for(int xx = 0; xx<xVerts; xx++) {
			double x = xx*xstep;
			points.add(new Vector2d(x, 0));
			points.add(new Vector2d(x, height-1));
		}
		for(int yy=1; yy<yVerts-1; yy++) {
			double y = yy*ystep;
			points.add(new Vector2d(0, y));
			points.add(new Vector2d(width-1, y));
		}
		
		double border_offset = 100;
		
		// Put a ring of vertices outside the boundary in order to fill out the voronoi tessellation 
		double xStepBoundary = ((xmax+border_offset) - (xmin-border_offset))/ (xVerts-1);
		double yStepBoundary = ((ymax+border_offset) - (ymin-border_offset))/ (yVerts-1);
		for(int xx = 0; xx<xVerts; xx++) {
			double x = xx*xStepBoundary - border_offset;
			double randx = rnd.nextGaussian()*sigma_x;
			double randy = rnd.nextGaussian()*sigma_y;
			points.add(new Vector2d(x + randx, -border_offset + randy));
			points.add(new Vector2d(x + randx, height-1 + border_offset + randy));
		}
		for(int yy=1; yy<yVerts-1; yy++) {
			double y = yy*yStepBoundary - border_offset;
			double randx = rnd.nextGaussian()*sigma_x;
			double randy = rnd.nextGaussian()*sigma_y;
			points.add(new Vector2d(-border_offset + randx, y + randy));
			points.add(new Vector2d(width-1 + border_offset + randx, y + randy));
		}
		
		switch(style) {
			case RANDOM: 
				{
					// N randomly distributed points
					for(int i=0; i<n; i++) {
						double x = xmin + rnd.nextDouble()*(xmax - xmin);
						double y = ymin + rnd.nextDouble()*(ymax - ymin);
						Vector2d vec = new Vector2d(x, y);
						points.add(vec);
					}
					break;
				}
			case SEMI_RANDOM:
			{
				for(int xIdx=1; xIdx<xVerts-1; xIdx++) {
					for(int yIdx=1; yIdx<yVerts-1; yIdx++) {
						double x = xIdx*xstep;
						double y = yIdx*ystep;
						double randx = rnd.nextGaussian()*sigma_x;
						double randy = rnd.nextGaussian()*sigma_y;
						
						// XXX Offset alternate rows by half a step
						if(yIdx%2==0) {
							x += xstep/4.0;
						}
						else {
							x -= xstep/4.0;
						}
						
						Vector2d vec = new Vector2d(x+randx, y+randy);
						points.add(vec);
					}
				}
				
				// XXX Add point in centre of one of the boxes
				Vector2d vec = new Vector2d(xmin+1.5*xstep,ymin+2.5*ystep);
				points.add(vec);
				
				break;
			}
		}
		
		DelaunayTriangulation2D delTri = new DelaunayIncremental(points);
		delTri.doTriangulation();
		VoronoiTessellation2D vorTess = new VoronoiTessellation2D(delTri);

		delColouredIm = delTri.getColouredImage(width, height, delTriColours);
		delWireframeIm = delTri.getWireframeImage(width, height);
		tessColouredIm = vorTess.getColouredImage(width, height, vorTessColours);
		tessWireframeIm = vorTess.getWireframeImage(width, height);
	}
	
	/**
	 * Main application entry point.
	 * @param args
	 * 	The command line arguments (ignored).
	 */
	public static void main(String[] args) {
		JFrame frame = new JFrame("Stained Glass Designer");
		frame.setLayout(new BorderLayout());
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.getContentPane().add(new StainedGlassDesigner(), BorderLayout.CENTER);
		frame.pack();
		frame.setVisible(true);
	}
	
}