package util;

import java.awt.Container;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Hashtable;

import javax.imageio.ImageIO;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPopupMenu;
import javax.swing.JSlider;
import javax.swing.filechooser.FileFilter;

/**
 * Utility methods for GUI building.
 *
 *
 * @author nrowell
 * @version $Id$
 */
public class GuiUtil {

	/**
	 * Make a simple JSlider.
	 * @param start
	 * 	The minimum value.
	 * @param end
	 * 	The maximum value.
	 * @param nLabels
	 * 	The number of labels, positioned uniformly along the slider.
	 * @param format
	 * 	Format string to use in formatting the label text.
	 * @return
	 * 	A JSlider.
	 */
	public static JSlider buildSlider(double start, double end, int nLabels, String format) 
	{
		int min = 0;
		int max = 100;
		
		final JSlider slider = new JSlider(JSlider.HORIZONTAL, min, max, 0);
		slider.setMajorTickSpacing(50);
		slider.setMinorTickSpacing(10);
		slider.setPaintTicks(true);
		slider.setPaintLabels(true);
		
		double delta = (end - start)/(nLabels-1);
		
		Hashtable<Integer, JLabel> labelTable = new Hashtable<>();
		for(int i=0; i<nLabels; i++) {
			double step = start + i*delta;
			labelTable.put(i*(max-min)/(nLabels-1), new JLabel(String.format(format, step)));
		}
		slider.setLabelTable( labelTable );
		
		return slider;
	}
	
	/**
	 * Adds a right-click menu to the {@link JComponent} that provides an option to save the contents of the
	 * component to an image file.
	 * @param parent
	 * 	The parent {@link Container}, used to anchor the menu and any error messages.
	 * @param child
	 * 	The {@link JComponent} which will be rendered to an image.
	 */
	public static void addRightClickMenuSaveJPanelImage(final Container parent, final JComponent child) {
		
	    // Add right-click menu with option to save image to file
	    final JMenuItem saveImageMenuItem = new JMenuItem("Save image to file");
	    
	    saveImageMenuItem.addActionListener(new ActionListener(){
			@Override
			public void actionPerformed(ActionEvent e) {
				
				BufferedImage im = new BufferedImage(child.getWidth(), child.getHeight(), BufferedImage.TYPE_INT_ARGB);
				child.paint(im.getGraphics());
				
				// Use file chooser dialog to select file save location:
				JFileChooser chooser = new JFileChooser();
				chooser.setDialogTitle("Set output file...");
				chooser.addChoosableFileFilter(new FileFilter() {
					public boolean accept(File file) {
			    		String filename = file.getName();
			    		return filename.endsWith(".png");
			    	}
			    	public String getDescription() {
			    		return "*.png";
			    	}
				});
				
				int userSelection = chooser.showSaveDialog(parent);
				 
				if (userSelection == JFileChooser.APPROVE_OPTION) {
					File file = chooser.getSelectedFile();
					// Add extension '.png' if this was not given by the user
					String extension = "";
					int i = file.getName().lastIndexOf('.');
					if (i > 0) {
					    extension = file.getName().substring(i+1);
					}
					if (!extension.equalsIgnoreCase("png")) {
					    file = new File(file.toString() + ".png");
					}
					
					try {
						ImageIO.write(im, "PNG", file);
					} catch (IOException e1) {
						JOptionPane.showMessageDialog(parent, "Error saving to file "+file.getAbsolutePath()+"!",
								"Error", JOptionPane.ERROR_MESSAGE);
					}
				}
			}
	    });
	    
	    JPopupMenu popupMenu = new JPopupMenu();
	    popupMenu.add(saveImageMenuItem);
	    child.setComponentPopupMenu(popupMenu);
	}
	
}