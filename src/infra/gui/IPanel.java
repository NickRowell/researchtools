package infra.gui;

import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Insets;
import java.awt.image.BufferedImage;
import java.io.File;

import javax.swing.BorderFactory;
import javax.swing.JPanel;
import javax.swing.border.BevelBorder;

import infra.io.Gnuplot;

/**
 * 
 * Name:
 *  IPanel.java
 * 
 * Purpose:
 *  Utility class that allows an image to be drawn inside a JPanel easily.
 * 
 * Language:
 * Java
 *
 * Author:
 * Nicholas Rowell
 * 
 */
public class IPanel extends JPanel
{
    
    /**
	 * The serial version UID.
	 */
	private static final long serialVersionUID = 1L;
	
	
	public BufferedImage image;
    
	/**
	 * Default constructor, setting up the border.
	 */
    public IPanel() {
//    	this.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));
    	this.setBorder(BorderFactory.createEmptyBorder(0,0,0,0));
    }
    
    /**
     * Construct an IPanel to display the given image
     * @param img
     * 		The BufferedImage to display
     */
    public IPanel(BufferedImage img) { 
    	this();
    	setImage(img);
    }
    
    /**
     * Sets the BufferedImage to display.
     * @param img
     * 		The BufferedImage to display.
     */
    public void setImage(BufferedImage img) {
    	image = img;
    	setPreferredSize();
    	repaint();
    }
    
    @Override
    public void paintComponent(Graphics g)
    {
        Insets insets = this.getBorder().getBorderInsets(this);
    	
        // Watch out for attempting to draw img before it is initialised.
        if(image!=null)
            g.drawImage(image, insets.left, insets.top, null);
    }
    
    public void setPreferredSize()
    {
        
        // Get image size
        int image_x = image.getWidth(this);
        int image_y = image.getHeight(this);
        
        Insets insets = new Insets(0, 0, 0, 0);
        
        // get border
        if(getBorder()!=null)
        	insets = this.getBorder().getBorderInsets(this);
        
        Dimension dims = new Dimension(image_x + insets.left + insets.right,
                                       image_y + insets.top + insets.bottom);
        
        setPreferredSize(dims);
        setMinimumSize(dims);
        setMaximumSize(dims);
    }
    
    /**
     * Run the given Gnuplot script, and display the plot.
     * 
     * TODO: move this outside of the plot panel itself?
     * 
     * @param script 
     */
    public void plotGnuplot(File script)
    {
            // Get plot image from Gnuplot
            BufferedImage img = Gnuplot.executeScript(script);
            
            if(img==null)
            	return;
            setImage(img);
    }
}