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
	private static File outputDir = new File("/home/nrowell/Temp/");

	// Size of square image to render
	static final int width  = 480;
	static final int height = 480;
	
	// Randomly generate a set of n points to triangulate
	static int n = 8;
	
	// Box within which to randomly distribute the points
	static double xmin = 0;
	static double xmax = width;
	// Step size between point centres in x direction
	static int xVerts = 4;
	static double xstep = width/xVerts;
	
	
	static double ymin = 0;
	static double ymax = height;
	// Step size between point centres in y direction
	static int yVerts = 4;
	static double ystep = height/yVerts;
	
	// Size of random displacement
	static double sigma = 20.0;
	
	
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
		
		style = Style.RANDOM;
		
		// Initialise the plots
		redraw();
		
		// Central panel showing the 4 images
		final IPanel tl = new IPanel(delColouredIm);
		final IPanel tr = new IPanel(delWireframeIm);
		final IPanel bl = new IPanel(tessColouredIm);
		final IPanel br = new IPanel(tessWireframeIm);
		
		JPanel plotPanel = new JPanel(new GridLayout(2,2));
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
					ImageIO.write(delColouredIm, "png", new File(outputDir, "delaunay_coloured.png"));
					ImageIO.write(delWireframeIm, "png", new File(outputDir, "delaunay_wireframe.png"));
					ImageIO.write(tessColouredIm, "png", new File(outputDir, "voronoi_coloured.png"));
					ImageIO.write(tessWireframeIm, "png", new File(outputDir, "voronoi_wireframe.png"));
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
		
		double xStepBoundary = (xmax - xmin)/ (xVerts-1);
		for(int xx = 0; xx<xVerts; xx++) {
			double x = xx*xStepBoundary;
			points.add(new Vector2d(x, 0));
			points.add(new Vector2d(x, height-1));
		}
		double yStepBoundary = (ymax - ymin)/ (yVerts-1);
		for(int yy=1; yy<yVerts-1; yy++) {
			double y = yy*yStepBoundary;
			points.add(new Vector2d(0, y));
			points.add(new Vector2d(width-1, y));
		}
		
		// Put a ring of vertices outside the boundary in order to fill out the voronoi tessellation 
		xStepBoundary = ((xmax+100) - (xmin-100))/ (xVerts-1);
		yStepBoundary = ((ymax+100) - (ymin-100))/ (yVerts-1);
		for(int xx = 0; xx<xVerts; xx++) {
			double x = xx*xStepBoundary;
			points.add(new Vector2d(x, -100));
			points.add(new Vector2d(x, height-1 + 100));
		}
		for(int yy=1; yy<yVerts-1; yy++) {
			double y = yy*yStepBoundary;
			points.add(new Vector2d(-100, y));
			points.add(new Vector2d(width-1 + 100, y));
		}
		
		switch(style) {
			case RANDOM: 
				{
					// N randomly distributed points
					for(int i=0; i<n; i++) {
						double x = xmin + Math.random()*(xmax - xmin);
						double y = ymin + Math.random()*(ymax - ymin);
						Vector2d vec = new Vector2d(x, y);
						points.add(vec);
					}
					break;
				}
			case SEMI_RANDOM:
			{
				Random rnd = new Random();
				for(double x = xmin+xstep; x<xmax-xstep; x+=xstep) {
					for(double y = ymin+ystep; y<ymax-ystep; y+=ystep) {
						double randx = rnd.nextGaussian()*sigma;
						double randy = rnd.nextGaussian()*sigma;
						Vector2d vec = new Vector2d(x+randx, y+randy);
						points.add(vec);
					}
				}
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