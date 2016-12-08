package fitting;

import java.util.List;
import java.util.ArrayList;

import stats.Statistics;
import data.Histogram;
import minimisation.nllsq.LevenbergMarquardt;
import Jama.*;


/**
 * This class tests the accuracy of the estimated parameter covariance that
 * the Levenberg-Marquardt algorithm obtains numerically.
 * 
 * It starts with a randomly generated set of data lying along a cubic function,
 * generates a random covariance matrix and uses this to add noise to the
 * data. It then fits a cubic function to the noisy data, and uses the true
 * data covariance matrix and an analytically determined parameter/data 
 * Jacobian to calculate the parameter covariance by first order propagation.
 * 
 * It then estimates the true parameter covariance by using many realisations
 * of the same data and fitting these with cubic functions. The sample
 * covariance of the resulting parameter sets is calculated in the standard
 * way.
 * 
 * The sample covariance measured from many independent parameter sets is then
 * compared to the parameter covariance estimated from a single set using
 * first order propagation. If the first order propagation technique is
 * correct, then the two will agree to within statistical fluctuations.
 * 
 * It also measures the chi-square distribution for all the fits for
 * comparison to theoretical distribution.
 * 
 * @author nickrowell
 */
public class LevenbergMarquardtTester_MonteCarlo extends LevenbergMarquardt{
    
    // Number of datapoints to generate
    static final int N_POINTS = 7;
    
    // Range in which to generate datapoints
    static final double N_MIN = 0;
    static final double N_MAX = 20;
    
    // True parameters of quadratic used to generate perfect data
    public final static double A_TRUE=1.5;
    public final static double B_TRUE=1.5;
    public final static double C_TRUE=1.5;     
    public final static double D_TRUE=1.5;   
    
    // Perfect data that will be used to generate noisy realisations
    public static Matrix DATA_X_PERFECT = new Matrix(N_POINTS,1);
    public static Matrix DATA_Y_PERFECT = new Matrix(N_POINTS,1);    
    
    // Parameters of cubic that is to be fit to the data. Set these to some
    // appropriate initial guess.
    // f(x) = A*x*x*x + B*x*x + C*x + D
    public double A=1;
    public double B=1;
    public double C=1;  
    public double D=1;
    
    // Finite step sizes in parameters used in Levenberg-Marquardt algorithm
    // to estimate numerical Jacobian approximation.
    public double dA=0.1, dB=0.1, dC=0.1, dD=0.1;
    
    // Data points for the fit: Two matrices of x values and y values
    public Matrix DATA_X;
    public Matrix DATA_Y;
    
    // Covariance matrix of fit, which will be randomly generated
    public Matrix covariance;
    
    // Main constructor - set data and covariance
    public LevenbergMarquardtTester_MonteCarlo(Matrix dat_x,
                                               Matrix dat_y,
                                               Matrix covar){
        // Copy reference to data
        DATA_X = dat_x;
        DATA_Y = dat_y;
        // Copy reference to covariance matrix
        covariance = covar;    
    }
    
    // Get model (predicted y coordinates) given data (x) and current parameters
    public Matrix getModel(){
        
        Matrix MOD = new Matrix(getDataN(),1);
        
        for(int p=0; p<N_POINTS; p++)
            MOD.set(p, 0, A*DATA_X.get(p,0)*DATA_X.get(p,0)*DATA_X.get(p,0) +
                          B*DATA_X.get(p,0)*DATA_X.get(p,0) + 
                          C*DATA_X.get(p,0) + D);

        return MOD;
    }
    
    // Get observed y coordinates
    public Matrix getData(){ return DATA_Y;}
    
    // Update data points (observed y coordinates) temporarily
    public boolean updateData(Matrix delta){
        
        for(int p=0; p<N_POINTS; p++){
            DATA_Y.set(p, 0, DATA_Y.get(p,0) + delta.get(p,0));
        }
        
        return true;
    }
    
    
    public Matrix getCovariance(){ return covariance;}

    
    
    public Matrix getParameters(){return new Matrix(new double[][]{{A},{B},{C},{D}});}
    
    public Matrix getParametersSteps(){
        return new Matrix(new double[][]{{dA},{dB},{dC},{dD}});
    }
    
    
    public void setParameters(Matrix param){
        A = param.get(0,0);
        B = param.get(1,0);
        C = param.get(2,0);
        D = param.get(3,0);
    }
   
    
    
    public int getParametersN(){ return 4;}
    public int getDataN(){ return N_POINTS;}
    public void printParameters(){ System.out.println("A = "+A+"\nB = "+B+"\nC = "+C+"\nD = "+D);}
    
       
    
    
    
    
    
    
    
    
    
    public static void main(String[] args){
        
        // First, generate perfect noise-free data that will be used to
        // test fit of parameters to noisy data.       
      
        for(int d=0; d<N_POINTS; d++){
            
            double RAND_X = Math.random()*(N_MAX - N_MIN);
            double RAND_Y = A_TRUE*RAND_X*RAND_X*RAND_X + B_TRUE*RAND_X*RAND_X + C_TRUE*RAND_X + D_TRUE;
            
            DATA_X_PERFECT.set(d,0,RAND_X);
            DATA_Y_PERFECT.set(d,0,RAND_Y);
                        
        }
        
        // Generate random covariance matrix for data
        Matrix covar = Statistics.getRandomCovarianceMatrix(0.1, N_POINTS);

        // Get data to be fitted by levenberg-marquardt algorithm
        Matrix dat_x = DATA_X_PERFECT;
        Matrix dat_y = Statistics.drawRandomVector(covar, DATA_Y_PERFECT);
        //Matrix dat_y = DATA_Y_PERFECT;
        
        // Initialise LM fit with noisy data
        LevenbergMarquardtTester_MonteCarlo test = new LevenbergMarquardtTester_MonteCarlo(dat_x, dat_y, covar);
        // Perform Levenberg-Marquardt fitting algorithm on data & function
        test.fit(1000, true);
                
        System.out.println("Parameter covariance by fourth order propagation of"
                + " data covariance using numerically determined parameter/"
                + "data Jacobian:");
        Matrix propagated_covariance = test.getFourthOrderCovariance();
        propagated_covariance.print(5,5);
        
        System.out.println("Fitted parameters:");
        test.printParameters();               
        
        /**
         * Now embark on Monte Carlo test campaign.
         */
        
        // Record histogra of chi-square statistic
        Histogram chiSquare = new Histogram(0, test.getDOF()*3.0, test.getDOF()*0.1, false);
        
        
        // Set number of Monte Carlo iterations
        int N = 1000000;
        
        // Store all parameter solutions to List for statistical analysis at end
        List<Matrix> solutions = new ArrayList<Matrix>();       
      
        // Sum them also to make calculating mean easy
        Matrix mean = new Matrix(4,1);
        
        for(int n=0; n<N; n++){
            
            // Get new realisation of y data points with noise and assign
            // it to y data values of LM object here
            test.DATA_Y = Statistics.drawRandomVector(test.covariance, 
                                                      DATA_Y_PERFECT);
            
            // Re-fit new data
            test.fit(1000, false);
            
            // Store chi-square
            chiSquare.add(test.getChi2());
            
            // Store parameters solution to List
            solutions.add(test.getParameters());
            // Sum all parameter solution sets
            mean.plusEquals(test.getParameters());
        }
        
        // Get mean estimated from data
        mean.timesEquals(1.0/(double)N);
        
        // Now measure sample covariance of parameters
        Matrix sample_covariance = new Matrix(test.getParametersN(),
                                              test.getParametersN());
        
        for (Matrix p_i : solutions) {
            Matrix diff = p_i.minus(mean);
            sample_covariance.plusEquals(diff.times(diff.transpose()));
        }
        
        sample_covariance.timesEquals(1.0/(double)(N-1));
        
        System.out.println("Sample covariance:");
        sample_covariance.print(5,5);
        
        System.out.println("Element-by-element division with propagated covariance:");
        for(int i=0; i<4; i++){
            for(int j=0; j<4; j++){
                System.out.print(propagated_covariance.get(i,j)/sample_covariance.get(i,j)+"\t");
            }
            System.out.println();
        }
        
        System.out.println("\nParameter covariance by Gnuplot method:");
        test.getParameterCovariance().print(5, 5);
        
        System.out.println("Measured chi-square distribution "
                           + "("+test.getDOF()+" dof)");
        System.out.println(chiSquare.print(true));
        
    }
    

}
