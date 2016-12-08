package stats;

import java.util.List;
import java.util.ArrayList;

import stats.Statistics;
import Jama.*;



public class StatisticsTester {
    
    
    public static void main(String args[]){

        monteCarloTestRandomVectors();
        
    }
    
    /**
     * This method uses a Monte Carlo technique to test that vectors produced
     * by the Statistics.drawRandomVector(true_covariance, true_mean) method
     * are distributed correctly.
     * 
     * It obtains a random covariance matrix for 5 hardwired data (1,2,3,4,5),
     * then uses the drawRandomVector() method to obtain many realisations
     * of this set of data. It then calculates the sample covariance matrix
     * for these realisations, and prints both this and the true covariance
     * matrix so the user can check they are consistent within statistical
     * uncertainty.
     * 
     */
    public static void monteCarloTestRandomVectors(){
       
        // Select mean vector
        Matrix true_mean = new Matrix(new double[][]{{1},
                                                     {2},
                                                     {3},
                                                     {4},
                                                     {5}});
        
        // Get random covariance matrix for elements
        Matrix covariance = Statistics.getRandomCovarianceMatrix(1.0, 5);
        
        // List to store all random realisations of noisy data
        List<Matrix> realisations = new ArrayList<Matrix>();
    
        // Sum realisations to get mean easily
        Matrix est_mean = new Matrix(5,1);
        
        // Number of random realisations
        int N = 1000000;
        
        for (int i = 0; i < N; i++) {
            
            // Get noisy realisation of data
            Matrix realisation = Statistics.drawRandomVector(covariance, 
                                                             true_mean);
            
            realisations.add(realisation);
            est_mean.plusEquals(realisation);
        }
        
        // Take mean of all random realisations
        est_mean.timesEquals(1.0/(double)N);
        
        // Get sample covariance matrix and sum terms
        Matrix sample_covariance = new Matrix(5,5);
        for (Matrix x : realisations){
            Matrix diff = x.minus(est_mean);
            sample_covariance.plusEquals(diff.times(diff.transpose()));
        }
        
        // Normalise sample covariance
        sample_covariance.timesEquals(1.0/(double)(N-1));
        
        System.out.println("True mean:"); true_mean.print(5,5);
        System.out.println("Estimated mean:"); est_mean.print(5,5);
        
        System.out.println("True covariance:"); covariance.print(5, 5);
        System.out.println("Sample covariance:"); sample_covariance.print(5, 5);
    
    }
    
}
