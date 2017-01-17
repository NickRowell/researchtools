package util.test;

import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.lang.reflect.Field;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;

import util.CharUtil;

/**
 * Tests for the {@link CharUtil} class.
 *
 * @author nrowell
 * @version $Id$
 */
public class TestCharUtil {
	
	/**
	 * Main application entry point.
	 * @param args
	 * 	The command line args (ignored).
	 * @throws IllegalAccessException 
	 * @throws IllegalArgumentException 
	 */
	public static void main(String[] args) throws IllegalArgumentException, IllegalAccessException {
		
		// Get all the fields and display the ones associated with Strings
		Field[] fields = CharUtil.class.getFields();
		
		// Map the field names to the field values for all Strings
		Map<String, String> chars = new HashMap<>();
		
		for(Field field : fields) {
			if(field.getType() == String.class) {
				chars.put(CharUtil.class.getName() + "." + field.getName(), (String)field.get(null));
			}
		}
		
		final JPanel charPanel = new JPanel(new GridLayout(chars.size(), 2));
		for(Entry<String, String> entry : chars.entrySet()) {
			charPanel.add(new JLabel(entry.getKey()));
			charPanel.add(new JLabel(entry.getValue(), SwingConstants.CENTER));
		}
		
		java.awt.EventQueue.invokeLater(
                new Runnable() 
                    {
                        @Override
                        public void run() 
                        {
                            JFrame tester = new JFrame("Unicode Characters");
                            tester.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                            tester.setLayout(new BorderLayout());
                            tester.add(charPanel, BorderLayout.CENTER);
                            tester.setLocationRelativeTo(null);
                            tester.pack();
                            tester.setVisible(true);
                        }
                    });
	}
	
}
