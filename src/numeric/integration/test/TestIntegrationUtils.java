package numeric.integration.test;

import java.awt.image.BufferedImage;
import java.io.IOException;

import javax.swing.JFrame;

import infra.gui.IPanel;
import infra.io.Gnuplot;
import numeric.integration.IntegrableFunction;
import numeric.integration.IntegrationUtils;


/**
 * Class tests the numerical integration algorithm by comparing against the analytical
 * integral of a known function.
 *
 *
 * @author nrowell
 * @version $Id$
 */
public class TestIntegrationUtils {
	
	public static void main(String[] args) throws IOException
	{
		// Set up function that we can integrate analytically in order to test
		// numerical integration.
		class TestFunction implements IntegrableFunction
		{
			
			public double evaluate(double x)
			{
				// Has the definite integral ln(x)
				return 1.0/x;
			}
			
			public double getDefiniteIntegral(double a, double b)
			{
				// Corresponds to function y = 1/x
				return Math.log(b) - Math.log(a);
			}
			
		}
		
		TestFunction testFunc = new TestFunction();
		
		StringBuilder script = new StringBuilder();
		script.append("set terminal pngcairo size 640,480\n");
		script.append("set yrange [-0.00001:0.00001]\n");
		script.append("set xlabel 'h'\n");
		script.append("set ylabel 'Error'\n");
		script.append("plot '-' w l t 'Analytic-Numerical'\n");
		
		double a = 0.1;
		double b = 5.0;
		
		for(double h = 1.0; h>0.001; h-=0.001)
		{
			double analytic = testFunc.getDefiniteIntegral(a, b);
			double numerical = IntegrationUtils.integrate(testFunc, a, b, h);
			script.append(h+"\t"+(analytic-numerical)+"\n");
		}
		script.append("e\n");
		
		final BufferedImage plot = Gnuplot.executeScript(script.toString());
		
		// Create and display the form
        java.awt.EventQueue.invokeLater(
                new Runnable() 
                {
                    @Override
                    public void run() 
                    {
                    	JFrame frame = new JFrame("Numerical Integration");
                    	frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                    	
                    	IPanel plotPanel = new IPanel(plot);
                    	
                    	frame.getContentPane().add(plotPanel);
                        frame.pack();
                        frame.setLocationRelativeTo(null);
                        frame.setVisible(true);
                    }
                });
		
	}
	
}
