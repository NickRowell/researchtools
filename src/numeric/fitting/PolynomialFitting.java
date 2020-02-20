/*
 * Gaia CU5 DU10
 *
 * (c) 2005-2020 Gaia Data Processing and Analysis Consortium
 *
 *
 * CU5 photometric calibration software is free software; you can redistribute
 * it and/or modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1 of the
 * License, or (at your option) any later version.
 *
 * CU5 photometric calibration software is distributed in the hope that it will
 * be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
 * General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this CU5 software; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 *
 *-----------------------------------------------------------------------------
 */

package numeric.fitting;

import java.awt.Color;
import java.awt.FlowLayout;
import java.util.Random;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.jfree.chart.ChartColor;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYErrorRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.xy.YIntervalSeries;
import org.jfree.data.xy.YIntervalSeriesCollection;

import Jama.Matrix;
import numeric.data.FloatList;
import numeric.functions.Polynomial;

/**
 * This utility class provides methods to fit a {@link Polynomial} of desired order to a set of datapoints, optionally
 * with outlier rejection based on sigma-clipping.
 *
 *
 * @author nrowell
 * @version $Id$
 */
public class PolynomialFitting {

	/**
	 * Maximum number of iterations in the outlier rejection.
	 */
	private static int MAX_ITERATIONS = 10;

	/**
	 * Private constructor to prevent instantiation.
	 */
	private PolynomialFitting() {
		
	}
	
	/**
	 * Solves for the polynomial parameters using weighted least squares and an iterative outlier rejection
	 * scheme. Following each iteration, inlying points are checked against the new model to see if they have
	 * become outliers. Also, points that were determined to be outliers in a previous iteration and excluded
	 * from the current fit and rechecked to see if they are now inliers (possible if the early iterations are
	 * plagued by far outliers). The fit iterates until there are no new outliers or inliers from one iteration
	 * to the next, or until a limit on the number of iterations has been reached.
	 * 
	 * @param order
	 * 	Order of the polynomial to fit, i.e. power of the leading term. Must have order < x.length.
	 * @param x
	 * 	The X coordinates of the datapoints
	 * @param y
	 * 	The Y coordinates of the datapoints
	 * @param yerrs
	 * 	The errors (standard deviation) on the Y coordinates of the datapoints, used to weight the fit.
	 * @param clip
	 * 	Outlier threshold, measured in standard deviations from the model.
	 * @param outliers
	 * 	The boolean array that on exit records whether each point was an outlier (true) or inlier (false) with respect to
	 * the final fit.
	 * @return
	 * 	A {@link Polynomial} fitted to the data using iterative outlier rejection based on sigma clipping. If too many
	 * outliers were found such that no fit is possible, then the returned value is null.
	 */
	public static Polynomial fitPolyClipped(int order, double[] x, double[] y, double[] yerrs, double clip, boolean[] outliers) {

		// Sanity checks: this also guarantees that there are enough inliers at the start that
		// the first iteration can get a solution (albeit possibly a bad one)
		checkInputs(order, x, y, yerrs, outliers);
		
		// Indicates when no new inliers/outliers are found from one iteration to the next
		boolean converged = false;
		// Indicates if too many outliers have been detected and the fit is no longer constrained
		boolean constrained = true;
		// Counts number of iterations made
		int iterations = 0;
		
		// Count initial number of inliers
		int nInliers = numberOfInliers(outliers);
		
		Polynomial polyFit = null;
		
		// Loop until we get no new inliers or outliers per iteration, or we hit the iteration limit, or
		// we end up with too few inliers to get a constrained solution.
		while(iterations++ < MAX_ITERATIONS && !converged && constrained) {

			// Fit a polynomial to the inliers
			polyFit = fitPoly(order, x, y, yerrs, outliers);
			
			// Check for new outliers and inliers
			int nNewOutliers = 0;
			int nNewInliers = 0;
			
			if(nInliers > order+1) {
				
				// Only perform outlier rejection if the fit is over-constrained, otherwise if nInliers==order+1,
				// then the regression line passes exactly through every point and we can't measure a dispersion
				// for sigma-clipping.
				
				// Get a robust estimate of the dispersion about the regression line, for outlier detection
				FloatList absRes = new FloatList();
				
				for(int i=0; i<x.length; i++) {
					if(!outliers[i]) {
						absRes.add((float) Math.abs(polyFit.getFn(x[i]) - y[i]));
					}
				}
				float sigma = 1.48f * absRes.getMAD()[1];
			
				// Update outliers set
				for(int i=0; i<x.length; i++) {
					
					// Distance of point from regression line, in units of the robust RMS
					double sigmas = Math.abs((polyFit.getFn(x[i]) - y[i])/sigma);
					
					if(!outliers[i]) {
						// Point was an inlier in the current fit - check to see if it's now an outlier
						if(sigmas > clip) {
							// Found a new outlier
							outliers[i] = true;
							nNewOutliers++;
						}
					}
					else {
						// Point was an outlier in the current fit - check to see if it's now an inlier
						if(sigmas < clip) {
							// Found a new inlier
							outliers[i] = false;
							nNewInliers++;
						}
					}
				}
			
			}
			else {
				// The regression line passes through every point; take no action. This has the effect of
				// indicating that the fit has converged, because the number of outliers and inliers
				// didn't change relative to the previous iteration.
			}
			
			// Check for convergence
			converged = (nNewOutliers==0 && nNewInliers==0);
			
			// Update the number of inliers
			nInliers += nNewInliers;
			nInliers -= nNewOutliers;
			
			// Check we're still constrained
			constrained = (nInliers > order);
		}
		
		// Reached 
		return polyFit;
	}
	
	/**
	 * Solves for the polynomial parameters using weighted least squares, with no outlier rejection.
	 * 
	 * @param order
	 * 	Order of the polynomial to fit, i.e. power of the leading term.
	 * @param x
	 * 	The X coordinates of the datapoints
	 * @param y
	 * 	The Y coordinates of the datapoints
	 * @param yerrs
	 * 	The errors (standard deviation) on the Y coordinates of the datapoints, used to weight the fit.
	 * @return
	 * 	A {@link Polynomial} least-squares fitted to the data points
	 */
	public static Polynomial fitPoly(int order, double[] x, double[] y, double[] yerrs) {
		double[] a = new double[order+1];
		double[][] aCov = new double[order+1][order+1];
		// Initialises to all inliers
		boolean[] outliers = new boolean[x.length];
		
		fitPolyCoeffs(order, x, y, yerrs, a, aCov, outliers);
		
		return new Polynomial(a, aCov);
	}
	
	/**
	 * Solves for the polynomial parameters using weighted least squares, with outliers rejected according to
	 * the given array.
	 * 
	 * @param order
	 * 	Order of the polynomial to fit, i.e. power of the leading term.
	 * @param x
	 * 	The X coordinates of the datapoints.
	 * @param y
	 * 	The Y coordinates of the datapoints.
	 * @param yerrs
	 * 	The errors (standard deviation) on the Y coordinates of the datapoints, used to weight the fit.
	 * @param outliers
	 * 	The boolean array that records whether each point is an outlier (true) or inlier (false).
	 * @return
	 * 	A {@link Polynomial} least-squares fitted to the inlying data points.
	 */
	public static Polynomial fitPoly(int order, double[] x, double[] y, double[] yerrs, boolean[] outliers) {
		double[] a = new double[order+1];
		double[][] aCov = new double[order+1][order+1];
		
		fitPolyCoeffs(order, x, y, yerrs, a, aCov, outliers);
		
		return new Polynomial(a, aCov);
	}
	
	/**
	 * Solves for the polynomial parameters using weighted least squares.
	 * 
	 * @param order
	 * 	Order of the polynomial to fit, i.e. power of the leading term.
	 * @param x
	 * 	The X coordinates of the datapoints
	 * @param y
	 * 	The Y coordinates of the datapoints
	 * @param yerrs
	 * 	The errors (standard deviation) on the Y coordinates of the datapoints, used to weight the fit.
	 * @param a
	 * 	On exit, this array will contain the polynomial coefficients
	 * @param aCov
	 * 	On exit, this array will contain the covariance matrix for the polynomial coefficients
	 * @param outliers
	 * 	The boolean array that records whether each point is an outlier (true) or inlier (false)
	 * 
	 */
	private static void fitPolyCoeffs(int order, double[] x, double[] y, double[] yerrs, double[] a, double [][] aCov, boolean[] outliers) {
		
		// Sanity checks
		checkInputs(order, x, y, yerrs, outliers);
		
		int nInliers = numberOfInliers(outliers);
		
		// Read arrays into Matrices for linear algebra operations
		double[][] xArr = new double[nInliers][order+1];
		double[][] yArr = new double[nInliers][1];
		// Compute the inverse covariance matrix directly as we know it's diagonal
		double[][] yInvCovArr = new double[nInliers][nInliers];
		
		// Index of next free elements in arrays (excluding outliers)
		int idx=0;
		for(int i=0; i<x.length; i++) {
			if(!outliers[i]) {
				
				yArr[idx][0] = y[i];
				yInvCovArr[idx][idx] = 1.0 / (yerrs[i]*yerrs[i]);
				
				double xx = 1.0;
				for(int k=0; k<order+1; k++)
				{
					xArr[idx][k] = xx;
					xx *= x[i];
				}
				idx++;
			}
		}
		
		Matrix xMat = new Matrix(xArr);
		Matrix yMat = new Matrix(yArr);
		Matrix yInvCovMat = new Matrix(yInvCovArr);
		
		Matrix xTEx = xMat.transpose().times(yInvCovMat.times(xMat));
		Matrix xTEy = xMat.transpose().times(yInvCovMat.times(yMat));
		
		// Solve weighted least squares
		Matrix aMat = xTEx.solve(xTEy);
		
		// Compute covariance matrix of parameters
		Matrix aCovMat = xTEx.inverse();
		
		// Read out the solution for the polynomial parameters, and the covariance matrix
		for(int i=0; i<order + 1; i++) {
			a[i] = aMat.get(i, 0);
			for(int j=0; j<order + 1; j++) {
				aCov[i][j] = aCovMat.get(i, j);
			}
		}
		
		return;
	}
	
	/**
	 * Performs basic sanity checks on the inputs.
	 * 
	 * @param order
	 * 	The order of the polynomial to fit.
	 * @param x
	 * 	The x coordinates
	 * @param y
	 * 	The cooresponding y coordinates
	 * @param yerrs
	 * 	The corresponding y coordinate errors
	 * @param outliers
	 * 	The boolean array that records whether each point is an outlier (true) or inlier (false)
	 */
	private static void checkInputs(int order, double[] x, double[] y, double[] yerrs, boolean[] outliers) {
		
		// Sanity checks
		if(x.length != y.length || x.length != yerrs.length || x.length != outliers.length) {
			throw new IllegalArgumentException("Inconsistent number of X and Y values!");
		}
		
		// Check there are sufficient inliers to constrain the fit
		int nInliers = numberOfInliers(outliers);
		
		if(order >= nInliers) {
			throw new IllegalArgumentException("Polynomial order = "+order+"; number of inlying data points = "+
											nInliers+": too few inlying data points to constrain fit!");
		}
	}
	
	/**
	 * Count the number of inliers, as indicated by the outliers array.
	 * 
	 * @param outliers
	 * 	The boolean array that records whether each point is an outlier (true) or inlier (false)
	 * @return
	 * 	The number of inliers.
	 */
	private static int numberOfInliers(boolean[] outliers) {
		int N = 0;
		for(int i=0; i<outliers.length ; i++) {
			if(!outliers[i]) {
				N++;
			}
		}
		return N;
	}
	
	/**
	 * Demo application to test and visually display the results of a polynomial fit to simulated data.
	 * 
	 * @param args
	 * 	The args (not used)
	 */
	public static void main(String[] args) {
		testUnclipped();
		testClipped();
	}
	
	/**
	 * Generates a set of simulated data points, including noise, fits a polynomial and displays the results.
	 * Outliers and outlier rejection is included.
	 */
	private static void testClipped() {

		// Create a LightPolynomialFunction that we'll use to generate some simulated data
		Polynomial poly = new Polynomial(new double[]{1.5, 1.345, 0.4623, 0.1573, 0.256});
		
		// Set the order of the polynomial to fit to the data
		int order = 4;
		
		// Set the range over which to generate simulated data
		double xmin = -2.0;
		double xmax = 2.0;
		
		// Set the number of data points to simulate
		int N = 50;
		
		// Fraction of outliers
		double olFrac = 0.1;
		
		// Sigma-clipping threshold for outier rejection
		double clip = 3.0;
		
		// Errors are drawn for each point from a Gaussian distribution with zero mean and a standard
		// deviation suppiled by this variable:
		double s = 0.5;
		
		Random random = new Random();
		
		double[] x = new double[N];
		double[] y = new double[N];
		double[] yerrs = new double[N];
		boolean[] outliers = new boolean[N];
		
		for(int n=0; n<N; n++) {
			
			// Draw a random x value:
			double xn = xmin + Math.random()*(xmax - xmin);
			
			// Get the noise-free y value
			double yn = poly.getFn(xn);
			
			// Get the error to add to this point
			double yerrn = 0.0;
			
			if(Math.random() < olFrac) {
				// Data point is an outlier
				yerrn = random.nextGaussian() * s * 10;
			}
			else {
				// Data point is not an outlier
				yerrn = random.nextGaussian() * s;
			}
			
			// Assign to the arrays
			x[n] = xn;
			y[n] = yn + yerrn;
			yerrs[n] = s;
		}
		
		// Fit a light polynomial function
		Polynomial polyFit = PolynomialFitting.fitPolyClipped(order, x, y, yerrs, clip, outliers);
		
		// Now make a plot of the result: plot the data points and the fitted function
		
		final YIntervalSeries[] allSeries = new YIntervalSeries[3];
        allSeries[0] = new YIntervalSeries("Inliers");
        allSeries[1] = new YIntervalSeries("Outliers");
        allSeries[2] = new YIntervalSeries("Fit");
        
        // Add the inlying input data
        for(int n=0; n<N; n++) {
        	if(!outliers[n]) {
        		allSeries[0].add(x[n], y[n], y[n]-yerrs[n], y[n]+yerrs[n]);
        	}
        }
        
        // Add the outlying input data
        for(int n=0; n<N; n++) {
        	if(outliers[n]) {
        		allSeries[1].add(x[n], y[n], y[n]-yerrs[n], y[n]+yerrs[n]);
        	}
        }
        
        // Add the fitted function
        for(double p=xmin; p<=xmax; p+=(xmax-xmin)/100.0) {
        	double fit = polyFit.getFn(p);
        	allSeries[2].add(p, fit, fit, fit);
        }
        
        // Make a series collection
        final YIntervalSeriesCollection data = new YIntervalSeriesCollection();
        for(YIntervalSeries series : allSeries) {
        	data.addSeries(series);
        }
        
        // Configure axes
        NumberAxis xAxis = new NumberAxis("X");
        xAxis.setRange(xmin, xmax);
        NumberAxis yAxis = new NumberAxis("Y");
        yAxis.setAutoRange(true);
        
        // Configure box renderer
        XYErrorRenderer xyLineAndShapeRenderer = new XYErrorRenderer();
        xyLineAndShapeRenderer.setDrawXError(false);
        xyLineAndShapeRenderer.setDrawYError(true);
        
        // Plot shapes for the inlier data points
        xyLineAndShapeRenderer.setSeriesLinesVisible(0, false);
        xyLineAndShapeRenderer.setSeriesShapesVisible(0, true);
        xyLineAndShapeRenderer.setSeriesPaint(0, ChartColor.DARK_GREEN);
        // Plot shapes for outliers
        xyLineAndShapeRenderer.setSeriesLinesVisible(1, false);
        xyLineAndShapeRenderer.setSeriesShapesVisible(1, true);
        xyLineAndShapeRenderer.setSeriesPaint(1, ChartColor.RED);
        // Plot lines for the fitted function
        xyLineAndShapeRenderer.setSeriesLinesVisible(2, true);
        xyLineAndShapeRenderer.setSeriesShapesVisible(2, false);
        xyLineAndShapeRenderer.setSeriesPaint(2, ChartColor.BLACK);
        
        // Configure plot
        XYPlot xyplot = new XYPlot(data, xAxis, yAxis, xyLineAndShapeRenderer);
        xyplot.setBackgroundPaint(Color.lightGray);
        xyplot.setDomainGridlinePaint(Color.white);
        xyplot.setDomainGridlinesVisible(true);
        xyplot.setRangeGridlinePaint(Color.white);

        // Configure chart
        final JFreeChart chart = new JFreeChart("Polynomial fitted to simulated data", xyplot);
        chart.setBackgroundPaint(Color.white);
        
        final JFrame frame = new JFrame("Polynomial with outlier rejection");
		
		SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                frame.setLayout(new FlowLayout());
                frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                frame.add(new ChartPanel(chart));
                frame.setSize(1500, 750);
                frame.pack();
                frame.setVisible(true);
            }
        });
	}
	
	/**
	 * Generates a set of simulated data points, including noise, fits a polynomial and displays the results. No
	 * outliers or outlier rejection is included.
	 */
	private static void testUnclipped() {

		// Create a LightPolynomialFunction that we'll use to generate some simulated data
		Polynomial poly = new Polynomial(new double[]{1.5, 1.345, 0.4623, 0.1573, 0.256});
		
		// Set the order of the polynomial to fit to the data
		int order = 4;
		
		// Set the range over which to generate simulated data
		double xmin = -2.0;
		double xmax = 2.0;
		
		// Set the number of data points to simulate
		int N = 50;
		
		// Errors are drawn for each point from a Gaussian distribution with zero mean and a standard
		// deviation suppiled by this variable:
		double s = 0.5;
		
		Random random = new Random();
		
		double[] x = new double[N];
		double[] y = new double[N];
		double[] yerrs = new double[N];
		
		for(int n=0; n<N; n++) {
			
			// Draw a random x value:
			double xn = xmin + Math.random()*(xmax - xmin);
			
			// Get the noise-free y value
			double yn = poly.getFn(xn);
			
			// Get the error to add to this point
			double yerrn = random.nextGaussian() * s;
			
			// Assign to the arrays
			x[n] = xn;
			y[n] = yn + yerrn;
			yerrs[n] = s;
		}
		
		// Fit a light polynomial function
		Polynomial polyFit = PolynomialFitting.fitPoly(order, x, y, yerrs);
		
		// Now make a plot of the result: plot the data points and the fitted function
		
		YIntervalSeries dataPointsSeries = new YIntervalSeries("Data");
        
        // Add the input data
        for(int n=0; n<N; n++) {
        	dataPointsSeries.add(x[n], y[n], y[n]-yerrs[n], y[n]+yerrs[n]);
        }
        
        XYSeries fitSeries = new XYSeries("Fit");
        XYSeries fitUpperSeries = new XYSeries("Fit upper one-sigma confidence limit");
        XYSeries fitLowerSeries = new XYSeries("Fit lower one-sigma confidence limit");
        
        // Add the fitted function
        for(double p=xmin; p<=xmax; p+=(xmax-xmin)/100.0) {
        	double[] fit = polyFit.getFnWithError(p);
        	double std = Math.sqrt(fit[1]);
        	
        	fitSeries.add(p, fit[0]);
        	fitUpperSeries.add(p, fit[0]+std);
        	fitLowerSeries.add(p, fit[0]-std);
        }
        
        // Make a series collection
        final YIntervalSeriesCollection data = new YIntervalSeriesCollection();
        data.addSeries(dataPointsSeries);
        
        XYSeriesCollection lines = new XYSeriesCollection();
        lines.addSeries(fitSeries);
        lines.addSeries(fitLowerSeries);
        lines.addSeries(fitUpperSeries);
        
        // Configure axes
        NumberAxis xAxis = new NumberAxis("X");
        xAxis.setRange(xmin, xmax);
        NumberAxis yAxis = new NumberAxis("Y");
        yAxis.setAutoRange(true);
        
        // Configure box renderer
        XYErrorRenderer dataPointRenderer = new XYErrorRenderer();
        dataPointRenderer.setDrawXError(false);
        dataPointRenderer.setDrawYError(true);
        
        // Plot shapes for the data points
        dataPointRenderer.setSeriesLinesVisible(0, false);
        dataPointRenderer.setSeriesShapesVisible(0, true);
        dataPointRenderer.setSeriesPaint(0, ChartColor.DARK_GREEN);
        
        XYLineAndShapeRenderer lineRenderer = new XYLineAndShapeRenderer();
        
        // Plot lines for the fitted function
        lineRenderer.setSeriesLinesVisible(0, true);
        lineRenderer.setSeriesShapesVisible(0, false);
        lineRenderer.setSeriesPaint(0, ChartColor.BLACK);
        
        lineRenderer.setSeriesLinesVisible(1, true);
        lineRenderer.setSeriesShapesVisible(1, false);
        lineRenderer.setSeriesPaint(1, ChartColor.GRAY);

        lineRenderer.setSeriesLinesVisible(2, true);
        lineRenderer.setSeriesShapesVisible(2, false);
        lineRenderer.setSeriesPaint(2, ChartColor.GRAY);
        
        // Configure plot
        XYPlot xyplot = new XYPlot();
        xyplot.setDomainAxis(xAxis);
        xyplot.setRangeAxis(yAxis);
        xyplot.setBackgroundPaint(Color.lightGray);
        xyplot.setDomainGridlinePaint(Color.white);
        xyplot.setDomainGridlinesVisible(true);
        xyplot.setRangeGridlinePaint(Color.white);
        
        // Add the data to be plotted
        xyplot.setDataset(0, data);
        xyplot.setRenderer(0, dataPointRenderer);
        
        xyplot.setDataset(1, lines);
        xyplot.setRenderer(1, lineRenderer);
        
        // Configure chart
        final JFreeChart chart = new JFreeChart("Polynomial fitted to simulated data", xyplot);
        chart.setBackgroundPaint(Color.white);
        
        final JFrame frame = new JFrame("Polynomial without outlier rejection");
		
		SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                frame.setLayout(new FlowLayout());
                frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                frame.add(new ChartPanel(chart));
                frame.setSize(1500, 750);
                frame.pack();
                frame.setVisible(true);
            }
        });
	}
	
}
