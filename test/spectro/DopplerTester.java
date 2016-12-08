package spectro;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import spectro.DopplerProfile;

/**
 *
 * @author nickrowell
 */
public class DopplerTester {
    
    
    public static void main(String[] args) throws IOException
    {
    
        // Frequency of line centre [Hz]
        double nu_0 = 800e6;
        
        // Gaussian width [Hz]
        double FWHM_D = 50e6;
        
        DopplerProfile doppler = new DopplerProfile(FWHM_D, nu_0);
    
        // Integrate doppler profile to check normalisation
        BufferedWriter out = new BufferedWriter(new FileWriter(new File("/home/nickrowell/Projects/WBS/WBS_simulations/WBS/C/Noise_analysis/Voigt_profile_tests/doppler_prf")));
       
//        // Sum of power
//        double sum2 = 0.0;
//        // integration step size [Hz]
//        double dnu = 10;
//        // Get line profile about the peak
//        for(double nu = 0; nu < 1.5e9; nu += dnu){
//            sum2 += doppler.getPowerAtFrequency(nu+0.5*dnu) * dnu;
//        }
//        System.out.println("Integrated power spectrum = "+sum2);
        
        // Number of samples in FFT
        double N = 2048;
        
        // Sample period [s]
        double dt = 3.333e-10;
        
        // Frequency step between FFT elements [Hz]
        double df = 1/(N * dt);
        
        for(double nu = 0; nu < 1.5e9; nu += df){
            out.write((nu)+"\t" +doppler.getPowerAtFrequency(nu)+"\n");
        }
        
        
        out.close();
    
    }
    
    
    
    
}
