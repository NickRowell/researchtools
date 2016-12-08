package numeric.minimisation.nllsq.exec;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.jfree.chart.ChartColor;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import numeric.minimisation.nllsq.algoimpl.LevenbergMarquardtCosineFitter;

/**
 * Class to test the {@link LevenbergMarquardtCosineFitter}.
 *
 *
 * @author nrowell
 * @version $Id$
 */
public class TestCosineFitter {

	// True parameter values
	private static final double A = 1.0;
	private static final double B = 0.3876;
	private static final double w = 1.4567;
	private static final double p = 0.234;
	
	/**
	 * Main application entry point.
	 * @param args
	 * 	The command line args - ignored.
	 */
	public static void main(String[] args) {
		
		// Number of random data points to draw
		int N = 100;
		
		// Range of abscissa over which to draw data
		double tMin = 0.0;
		double tMax = 10.0;
		
		// Noise amplitude
		double noise = 0.05;
		
		// Generate noisy realisations of the function
		List<double[]> obs = new ArrayList<>();
		for(int i=0; i<N; i++) {
			
			// Get random abscissa value
			double t = tMin + Math.random() * (tMax - tMin);
			
			// Get the true value of the underlying function
			double f = getTrueValue(t);
			
			// Add observational error
			f += Math.random()*noise;
			
			obs.add(new double[]{t, f});
		}
		
		LevenbergMarquardtCosineFitter cosFit = new LevenbergMarquardtCosineFitter(obs);
		
		// Set initial guess parameters
		cosFit.A = A+Math.random()*0.1;
		cosFit.B = B+Math.random()*0.1;
		cosFit.w = w+Math.random()*0.5;
		cosFit.p = p+Math.random()*0.5;
		
		// Display chart of data points, initial and final fits.
		XYSeries fit = new XYSeries("Fit");
		XYSeries init = new XYSeries("Init");
		XYSeries trueF = new XYSeries("True");
		XYSeries data = new XYSeries("Data");
		
		for(double[] d : obs) {
			data.add(d[0], d[1]);
		}
		// For the fitted function, plot 100 data points within the range
		double step = (tMax - tMin)/100;
		for(int i=0; i<100; i++) {
			double t = tMin + i*step;
			trueF.add(t, evaluate(t, A, B, w, p));
		}
		for(int i=0; i<100; i++) {
			double t = tMin + i*step;
			init.add(t, evaluate(t, cosFit.A, cosFit.B, cosFit.w, cosFit.p));
		}
		
		cosFit.invoke();
		
		for(int i=0; i<100; i++) {
			double t = tMin + i*step;
			fit.add(t, evaluate(t, cosFit.A, cosFit.B, cosFit.w, cosFit.p));
		}
		
		XYSeriesCollection seriesCollection = new XYSeriesCollection();
		seriesCollection.addSeries(data);
		seriesCollection.addSeries(trueF);
		seriesCollection.addSeries(init);
		seriesCollection.addSeries(fit);
		
		XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
		
		renderer.setSeriesLinesVisible(0, false);
        renderer.setSeriesShapesVisible(0, true);
        renderer.setSeriesPaint(0, ChartColor.BLACK);
        renderer.setSeriesStroke(0, new BasicStroke(1.0f));
        
        renderer.setSeriesLinesVisible(1, true);
        renderer.setSeriesShapesVisible(1, false);
        renderer.setSeriesPaint(1, ChartColor.GREEN);
        renderer.setSeriesStroke(1, new BasicStroke(1.0f));
        
        renderer.setSeriesLinesVisible(2, true);
        renderer.setSeriesShapesVisible(2, false);
        renderer.setSeriesPaint(2, ChartColor.RED);
        renderer.setSeriesStroke(2, new BasicStroke(1.0f));
        
        renderer.setSeriesLinesVisible(3, true);
        renderer.setSeriesShapesVisible(3, false);
        renderer.setSeriesPaint(3, ChartColor.BLUE);
        renderer.setSeriesStroke(3, new BasicStroke(1.0f));
		
        NumberAxis xAxis = new NumberAxis("Time");
	    NumberAxis yAxis = new NumberAxis("Amplitude");
        
        // Update the plot dataset
	    XYPlot plot = new XYPlot(null, xAxis, yAxis, renderer);
        plot.setBackgroundPaint(Color.lightGray);
        plot.setDomainGridlinePaint(Color.white);
        plot.setRangeGridlinePaint(Color.white);
        plot.setDataset(seriesCollection);

        final ChartPanel chartPanel = new ChartPanel(new JFreeChart("CosineFitter test", plot));
        
        final JFrame frame = new JFrame("CosineFitter Test Application");
		
		SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
            	frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                frame.setLayout(new BorderLayout());
                frame.add(chartPanel, BorderLayout.CENTER);
                frame.setSize(1500, 750);
                frame.pack();
                frame.setVisible(true);
            }
        });
        
        
	}
	
	private static double evaluate(double t, double A, double B, double w, double p) {
		return B + A * Math.cos(w * t + p);
	}
	
	
	/**
	 * Computes the true value of the function at the abscissa value.
	 * @param t
	 * @return
	 */
	private static double getTrueValue(double t) {
		return evaluate(t, A, B, w, p);
	}
	
}
