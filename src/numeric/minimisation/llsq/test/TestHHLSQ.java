package numeric.minimisation.llsq.test;

import infra.gui.IPanel;
import infra.io.Gnuplot;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;

import numeric.data.DoubleList;
import numeric.minimisation.llsq.HouseHolderLeastSquares;
import Jama.Matrix;

/**
 * Class tests the solution to a least squares problem by Householder Least Squares, and investigates
 * the merging of solutions across consecutive datasets using a priori information in the HHLSQ problem.
 *
 * TODO: add side plot showing convergence of polynomial coefficient estimates as more datasets are added
 *
 * @author nrowell
 * @version $Id$
 */
public class TestHHLSQ extends JPanel {
	
	/**
	 * The serial version UID
	 */
	private static final long serialVersionUID = 282831465434365334L;
	
	// List of all parameter solutions, for tracking evolution
	List<Matrix> solutions = new LinkedList<>();
	// List of all datapoints used over all sub-fits
	List<double[]> allPoints = new LinkedList<>();
	
	// Previous transformed information matrix and observation vector, carried forward and
	// used to improve current estimate
	Matrix U_prev = null;
	Matrix z1_prev = null;
	
	// Cache parameter solutions
	List<Double> a0s = new LinkedList<>();
	List<Double> a1s = new LinkedList<>();
	List<Double> a2s = new LinkedList<>();
	List<Double> a3s = new LinkedList<>();

	// Set the cubic polynomial coefficients. These are what we're estimating.
	double a0 = 0.5;
	double a1 = -1.3;
	double a2 = 3.1;
	double a3 = 0.003;
	
	int nPar = 4;
	
	// Range of the (uniformly distributed) random error added to observations.
	// Errors added to the observed value of the cubic will be distributed in [-e:e]
	double e = 1.5;
	
	// Number of points to generate per dataset
	int N_OBS = 2;
	
	// Range in which to generate X values
	double x_min = -1.0;
	double x_max =  1.0;
	
	// Use previous solutions to compute a 'running' solution?
	boolean memory=true;
	
	public TestHHLSQ() {
		
		// IPanel to display plot
		final IPanel ipanel = new IPanel();
		ipanel.setPreferredSize(new Dimension(640,480));
		
		// Multiple panels to display evolution of each parameter
		final IPanel param0Panel = new IPanel();
		final IPanel param1Panel = new IPanel();
		final IPanel param2Panel = new IPanel();
		final IPanel param3Panel = new IPanel();
		
		param0Panel.setPreferredSize(new Dimension(200,120));
		param1Panel.setPreferredSize(new Dimension(200,120));
		param2Panel.setPreferredSize(new Dimension(200,120));
		param3Panel.setPreferredSize(new Dimension(200,120));
		
		JPanel paramsPanel = new JPanel();
		paramsPanel.setLayout(new GridLayout(4,1));
		paramsPanel.add(param0Panel);
		paramsPanel.add(param1Panel);
		paramsPanel.add(param2Panel);
		paramsPanel.add(param3Panel);
				
		// JButton to control iterations of solver
		JButton button = new JButton("Add more observations -->");
		
		button.addActionListener(new ActionListener(){

			@Override
			public void actionPerformed(ActionEvent e) {
				try {
					BufferedImage[] plots = nextSolution();
					ipanel.setImage(plots[0]);
					param0Panel.setImage(plots[1]);
					param1Panel.setImage(plots[2]);
					param2Panel.setImage(plots[3]);
					param3Panel.setImage(plots[4]);
						
				} catch (IOException e1) {
				}
			}
			
		});
		
		setLayout(new BorderLayout());
		add(ipanel,BorderLayout.CENTER);
		add(paramsPanel,BorderLayout.EAST);
		add(button, BorderLayout.SOUTH);
	}
	
	/**
	 * Generates a random set of datapoints, optionally combines these with the HouseHolder results
	 * from the previous iterations, then solves for the polynomial coefficients and produces a plot
	 * of the current best fit and all datapoints.
	 * 
	 * @return
	 * 	BufferedImage containing a plot to display.
	 * @throws IOException
	 */
	private BufferedImage[] nextSolution() throws IOException {
		
		// Build information and observation arrays. If we have a previous solution,
		// then leave space for this in the arrays
		int N_ROWS = (U_prev==null ? N_OBS : N_OBS + U_prev.getRowDimension());
		
		double[][] A = new double[N_ROWS][4];
		double[][] Y = new double[N_ROWS][1];
		
		// Loop over rows
		int i=0;
		
		// Enter any previous solution
		if(U_prev!=null) {
			// Loop over rows
			for(; i<U_prev.getRowDimension(); i++) {
				// Loop over columns
				for(int j=0; j<4; j++) {
					A[i][j] = U_prev.get(i, j);
				}
				Y[i][0] = z1_prev.get(i, 0);
			}
		}
		
		// Enter new observations
		for(; i<N_ROWS; i++) {
			
			// Get random x
			double x = x_min + (x_max-x_min)*Math.random();
			
			// Get noise-free cubic at x
			double y = getCubic(a0,a1,a2,a3,x);
			
			// Add noise to get the observed y
			double y_obs = y + e*(2*Math.random()-1.0);
			
			// Load values into system matrices
			A[i][0] = 1.0;
			A[i][1] = x;
			A[i][2] = x*x;
			A[i][3] = x*x*x;
			Y[i][0] = y_obs;
			
			allPoints.add(new double[]{x,y_obs});
		}
		
		
		System.out.println("A = "); new Matrix(A).print(5, 5);
		System.out.println("Y = "); new Matrix(Y).print(5, 5);
		
		
		
		// Solve for coefficients
		HouseHolderLeastSquares hhlsq = new HouseHolderLeastSquares(new Matrix(A), new Matrix(Y));
		
		solutions.add(hhlsq.x);
		
		// Get the coefficient solutions
		double a0_est = hhlsq.x.get(0, 0);
		double a1_est = hhlsq.x.get(1, 0);
		double a2_est = hhlsq.x.get(2, 0);
		double a3_est = hhlsq.x.get(3, 0);

		// Save to lists
		a0s.add(a0_est);
		a1s.add(a1_est);
		a2s.add(a2_est);
		a3s.add(a3_est);
		
		//////////////////////////////////////////////////////////////////
		//                                                              //
		// Plot all points, the current points, the fit and the model   //
		//                                                              //
		//////////////////////////////////////////////////////////////////
		
		// Build GNUplot script
		File script = new File(new File("/home/nrowell/Temp"), "script.p");
		BufferedWriter scriptOut = new BufferedWriter(new FileWriter(script));
		// Write gnuplot header
		scriptOut.write("set terminal pngcairo size 640,480\n");
		scriptOut.write("set xrange ["+x_min+":"+x_max+"]\n");
		scriptOut.write("set yrange [-1:6]\n");
		scriptOut.write("set key off\n");
		scriptOut.write("f(x) = "+a0_est+" + "+a1_est+"*x + "+a2_est+"*x*x + "+a3_est+"*x*x*x\n");
		scriptOut.write("g(x) = "+a0+" + "+a1+"*x + "+a2+"*x*x + "+a3+"*x*x*x\n");
		scriptOut.write("plot g(x) ls 0 lw 2, f(x), '-' w p pt 5 ps 0.25 lc rgbcolor 'red', '-' w p pt 5 lc rgbcolor 'green'\n");
		// Write all points
		for(double[] point : allPoints) {
			scriptOut.write(point[0] + " " + point[1] + "\n");
		}
		scriptOut.write("e\n");
		// Write points for the current dataset
		for(int p=0; p<N_OBS; p++) {
			
			// Get row offset in A & Y
			int off = (U_prev==null ? 0 : U_prev.getRowDimension());
			
			double x = A[p+off][1];
			double y = Y[p+off][0];
			scriptOut.write(x + " " + y + "\n");
		}
		scriptOut.close();
		
		BufferedImage plot = Gnuplot.executeScript(script);
		
		//////////////////////////////////////////////////////////////////
		//                                                              //
		// Plot the evolution of each parameter                         //
		//                                                              //
		//////////////////////////////////////////////////////////////////
		
		BufferedImage plotA0 = getParameterPlot(a0s, a0, "a_0");
		BufferedImage plotA1 = getParameterPlot(a1s, a1, "a_1");
		BufferedImage plotA2 = getParameterPlot(a2s, a2, "a_2");
		BufferedImage plotA3 = getParameterPlot(a3s, a3, "a_3");
		
		// Carry forward the transformed information and observation matrices
		if(memory) {
			U_prev  = hhlsq.U;
			z1_prev = hhlsq.z1;
		}

		return new BufferedImage[]{plot, plotA0, plotA1, plotA2, plotA3};
		
	}
	
	private BufferedImage getParameterPlot(List<Double> estValues, double trueValue, String label) throws IOException {
		
		// Compute Y axis range so that true value is at centre, and range extends to X times the median of
		// absolute deviation
		DoubleList devs = new DoubleList();
		for(Double estValue : estValues) {
			devs.add(Math.abs(estValue - trueValue));
		}
		double yMin = trueValue - 2*devs.getMedian();
		double yMax = trueValue + 2*devs.getMedian();
		
		File script = new File(new File("/home/nrowell/Temp"), "script.p");
		BufferedWriter scriptOut = new BufferedWriter(new FileWriter(script));
		// Write gnuplot header
		scriptOut.write("set terminal pngcairo size 200,128 enhanced\n");
		scriptOut.write("set xrange [0:"+estValues.size()+"]\n");
		scriptOut.write("set yrange ["+yMin+":"+yMax+"]\n");
		scriptOut.write("set key off\n");
		scriptOut.write("set ytics scale 0.25 out font 'Helvetica,10'\n");
		scriptOut.write("set xtics format '' scale 0.25 out nomirror\n");
		scriptOut.write("set ylabel '"+label+"' offset 2,0\n");
		scriptOut.write("g(x) = "+trueValue+"\n");
		scriptOut.write("plot g(x) ls 0 lw 2, '-' w p pt 5 ps 0.25 lc rgbcolor 'red'\n");
		// Write all a0 solutions
		for(Double a0 : estValues) {
			scriptOut.write(a0 + "\n");
		}
		scriptOut.write("e\n");
		
		scriptOut.close();
		return Gnuplot.executeScript(script);
	}
	
	
	public static void main(String[] args) throws IOException {
		
		// Create and display the form
        java.awt.EventQueue.invokeLater(
                new Runnable() 
                {
                    @Override
                    public void run() 
                    {
                    	JFrame frame = new JFrame("Test Householder least squares");
                    	frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                    	frame.getContentPane().add(new TestHHLSQ());
                        frame.pack();
                        frame.setLocationRelativeTo(null);
                        frame.setVisible(true);
                    }
                });
		
	}
	
	
	/**
	 * Computes the cubic polynomial at x.
	 * 
	 * @param A
	 * 	First coefficient.
	 * @param B
	 * 	Second coefficient.
	 * @param C
	 * 	Third coefficient.
	 * @param D
	 * 	Fourth coefficient.
	 * @param x
	 *  The coordinate at which to compute the polynomial
	 * @return
	 * 	A + B*x + C*x*x + D*x*x*x
	 */
	private static double getCubic(double A, double B, double C, double D, double x) {
		return A + B*x + C*x*x + D*x*x*x;
	}
	

}
