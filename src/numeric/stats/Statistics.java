/*
 * CONFIDENTIAL
 * 
 * Copyright (c) 2010 University of Dundee
 *
 * Name:
 * Statistics.java
 *
 * Purpose:
 * This class provides various statistical methods, such as the generation of
 * noisy random correlated data and random covariance matrices.
 *
 * Language:
 * Java
 *
 * Author:
 * Nicholas Rowell (nickrowell@computing.dundee.ac.uk)
 *
 *
 *
 */
package numeric.stats;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

/**
 *
 * @author nickrowell
 */
public class Statistics {
    
    /** Static Random instance to avoid generating new ones with each call */
    static Random random = new Random();
    
    
    
    
    /**
     * Compute the binomial coefficient 'N choose M'. This is used to determine
     * the number of trials required by the Delaunay Triangulation algorithm.
     */
    public static double nChooseM(int n, int k)
    {
        double nk = 1l;
        for (int i = n - k + 1; i <= n; i++)
        {
            nk *= i;
		}
		for (int i = 1; i <= k; i++)
	        {
	            nk /= i;
		}
	        return nk;
    }
    
    
    
    
    /**
     * Draw a random vector given a covariance matrix for the elements
     * and a mean vector.
     * @param covar
     * @return 
     */
    public static Matrix drawRandomVector(Matrix covar, Matrix mean){
        
        // How many data to draw?
        int N = covar.getColumnDimension();     
        
        // Basic checks on matrix sizes and shapes etc.
        if(covar.getColumnDimension()!=covar.getRowDimension())
            throw new RuntimeException();
        if(mean.getColumnDimension()!=1 || mean.getRowDimension()!=N)
            throw new RuntimeException();
                
        // Get Eigenvalue decomposition of covariance matrix
        EigenvalueDecomposition eig = covar.eig();
        
        // Get dispersions along the principal axes of the covariance hyper-ellipsoid
        Matrix evals = eig.getD();
        
        // Random vector in principal axes frame
        Matrix rand = new Matrix(N,1);
        // Draw elements from unit Gaussian and scale according to axis sigma
        for(int d=0; d<N; d++)
            rand.set(d, 0, Math.sqrt(evals.get(d, d))*random.nextGaussian());
        
        // Transform this back to original frame using eigenvectors
        // of covariance matrix
        rand = eig.getV().times(rand);    
                
        // Add random error onto mean and return
        return mean.plus(rand);
        
    }
    
    /**
     * Generate an N-by-N random positive definite matrix that can be used
     * to simulate a covariance matrix. The sigma factor is a scale factor
     * and not necessarily representative of the standard deviation of the
     * matrix.
     * @return 
     */
    public static Matrix getRandomCovarianceMatrix(double sigma, int N){
                
        // Generate random covariance matrix for data
        Matrix covar = new Matrix(N,N);
        for(int i=0; i<N; i++)
            for(int j=0; j<N; j++)
                covar.set(i, j, random.nextGaussian());

        // Scale to sigma
        covar.timesEquals(sigma);
        
        // Make into a nice positive definite matrix by squaring
        covar = covar.transpose().times(covar);
        
        return covar;
    }
    
    
    /**
     * Get probability density for chi-square statistic given number of
     * degrees of freedom.
     */
    public static double getChi2(int dof_int, double x){
    
        double dof = (double)dof_int;
        
        double gamma = getGammaHalfInt(dof_int);
        
        double A = Math.pow(x, dof/2 - 1) / (Math.pow(2, dof/2) * gamma);
        
        return A * Math.exp(-x/2);
    
    }
    
    /**
     * Half-integer values of the Gamma function.
     * 
     * @param numerator Gamma function is evaluated at numerator/2
     * @return Gamma(numerator/2)
     * 
     */
    private static double getGammaHalfInt(int numerator){
    
        if(numerator < 1 || numerator > 7)
            throw new RuntimeException("getGammaHalfInt() invalid argument: "
                                       +numerator);
    
        switch(numerator){
            case 1: return Math.sqrt(Math.PI);            // Gamma(1/2)
            case 2: return 1;                             // Gamma(1)
            case 3: return 0.5*Math.sqrt(Math.PI);        // Gamma(3/2)
            case 4: return 1;                             // Gamma(2)
            case 5: return 0.75*Math.sqrt(Math.PI);       // Gamma(5/2)
            case 6: return 2;                             // Gamma(3)
            case 7: return (15.0/8.0)*Math.sqrt(Math.PI); // Gamma(7/2)
        }
        
        return 0.0;
    }
    
    
    /**
     * This method provides a 1D Gaussian convolution kernel suitable for 
     * application to binned data. It integrates the Gaussian over the width
     * of the bins (centred on central bin) to obtain the fraction of the
     * contents of the central bin that are scattered into each of the 
     * surrounding bins. The integration proceeds outward from the central bin,
     * and halts when 98% of the probability has been captured (1% lost from 
     * each wing). The returned error kernel will therefore be wider if the 
     * standard deviation is larger. The kernel is re-normalised to correct for
     * the missing probability, i.e. all the elements are scaled so that their
     * sum is one.
     * 
     * Only one half of the error kernel is returned because it is symmetric.
     * The central bin is in element 0.
     * 
     * @param sigma     Standard deviation.
     * @param bin_width Width of the bins.
     * @return 
     */
    public static double[] getSymGausKernel(double sigma, double bin_width){
 
        // Sanity checks
        if(sigma<0)
            throw new IllegalArgumentException("Sigma must be positive");
        
        if(sigma==0) return new double[]{1.0};
        
        // List to store probability in each bin.
        List<Double> values = new LinkedList<Double>();
        
        // Get integral over half-width of first bin and start accumulating 
        // total integrated probability.
        // Value of integral approaches 0.5 as kernel bins carry less and less
        // probability.        
        double integral = Gaussian.Phi(bin_width/2.0, 0, sigma) - 0.5;
        
        // Calculate first bin (centred on mean) using symmetry.
        values.add(integral*2.0);
        
        // Continue calculating bins until integral drops below threshold.
        // Stop adding kernel points once all but 2E-2 of the probability has
        // been captured by kernel.
        int bin = 1;
        
        while(integral < 0.5 - 1E-2){
            
            double bin_contents = Gaussian.Phi(bin* bin_width + bin_width/2.0, 0, sigma) -
                                  Gaussian.Phi((bin-1)* bin_width + bin_width/2.0, 0, sigma);
            
            values.add(bin_contents);
            
            // Sum all probability elements...
            integral += bin_contents;
            
            // On to next bin...
            bin++;
            
        }
        
        // Total integrated probability captured by truncated kernel
        double total = integral*2.0;
        
        // Copy values out to array and re-normalise
        double[] kernel = new double[values.size()];
        for(int b=0; b<values.size(); b++){
            // Renormalisation step
            kernel[b] = values.get(b)/total;
        }
    
        return kernel;
        
    }
    
    
    
    
    
}
