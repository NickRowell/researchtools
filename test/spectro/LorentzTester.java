package spectro;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import spectro.LorentzProfile;

/**
 *
 * @author nickrowell
 */
public class LorentzTester {
    
    
    public static void main(String[] args) throws IOException
    {
        
        // Frequency of line centre [Hz]
        double nu_0 = 800e6;
        
        // Lorentz width [Hz]
        double FWHM_L = 30e6;
        
        LorentzProfile lorentz = new LorentzProfile(FWHM_L, nu_0);
        
        // Integrate lorentz profile to check normalisation
        BufferedWriter out = new BufferedWriter(new FileWriter(new File("/home/nickrowell/Projects/WBS/WBS_simulations/WBS/C/Noise_analysis/Voigt_profile_tests/lorentz_prf")));
        
        
//        // Sum of power
//        double sum2 = 0.0;
//        // integration step size [Hz]
//        double dnu = 10;
//        for(double nu = 0; nu < 1.5e9; nu += dnu){
//            sum2 += lorentz.getPowerAtFrequency(nu+0.5*dnu) * dnu;
//        }        
//        System.out.println("Integrated power spectrum = "+sum2);
        
        
        // Number of samples in FFT
        double N = 2048;
        
        // Sample period [s]
        double dt = 3.333e-10;
        
        // Frequency step between FFT elements [Hz]
        double df = 1/(N * dt);       
        
        // Get line profile about the peak
        for(double nu = 0; nu < 1.5e9; nu += df){
            out.write((nu+df/2)+"\t" +lorentz.getPowerAtFrequency(nu+df/2)+"\n");
        }
        
        out.close();
        
    }
    
}