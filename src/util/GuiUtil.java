package util;

import java.util.Hashtable;

import javax.swing.JLabel;
import javax.swing.JSlider;

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
	
	
	
}
