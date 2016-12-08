package spectro;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import spectro.Oscillator;
import spectro.VoigtProfile;
import misc.Complex;

/**
 * Class used to test features of VoigtBroadening class.
 * @author nickrowell
 */
public class VoigtProfileTester {
    
    
    
    public static void main(String[] args) throws IOException{

        //printTaps();
        printVoigtPowerSpectrum();
        //generateVoigtTimeSeries();
        
    }
    
    
    private static void printTaps() throws IOException{
        
        // Frequency of line centre [Hz]
        double nu_0 = 800e6;
        // Lorentz width [Hz]
        double FWHM_L = 30e6;
        // Doppler width [Hz]
        double FWHM_D = 50e6;
        
        VoigtProfile voigt = new VoigtProfile(nu_0, FWHM_L, FWHM_D);
        
        System.out.println("Voigt FWHM = "+voigt.getVoigtFWHM() + " Hz");
        
        Complex[] taps = voigt.getFIRTaps(3.333e-10, 100);
        
        
        for(Complex tap : taps)
            System.out.print(tap.getRe() + " ");
        
    
    
    }
    
    
    
    private static void printVoigtPowerSpectrum() throws IOException{
            
        // Frequency of line centre [Hz]
        double nu_0 = 800e6;
        
        // Lorentz width [Hz]
        double FWHM_L = 30e6;
        // Doppler width [Hz]
        double FWHM_D = 50e6;
        
        VoigtProfile voigt = new VoigtProfile(nu_0, FWHM_L, FWHM_D);
        
        System.out.println("Voigt FWHM = "+voigt.getVoigtFWHM());
        
        BufferedWriter out = new BufferedWriter(new FileWriter(new File("/home/nickrowell/Projects/WBS/WBS_simulations/WBS/C/Noise_analysis/Voigt_profile_tests/voigt_prf")));
        
        // Number of samples in FFT
        double N = 2048;
        
        // Sample period [s]
        double dt = 3.333e-10;
        
        // Frequency step between FFT elements [Hz]
        double df = 1/(N * dt);
        
        
        // Total signal power
        
        
        // Get line profiles about the peak
        for(double nu = 0; nu < 1.5e9; nu += df){
            out.write((nu+df/2)+"\t" +voigt.getPowerAtFrequency(nu+df/2)+"\n");
        }
        
        // Integrate...
        double dnu = 1e5;  // step size
        
        double sum2 = 0.0;
        for(double nu = 0; nu < 1.5e9; nu += dnu){
            sum2 += voigt.getPowerAtFrequency(nu+dnu/2) * dnu;
        }

        
        System.out.println("Integrated power spectrum = "+sum2);
        
        out.close();

        
    }
    
    
    private static void generateVoigtTimeSeries() throws IOException{
    
        // Frequency of line centre [Hz]
        double nu_0 = 800e6;
        // Lorentz width [Hz]
        double FWHM_L = 30e6;
        // Doppler width [Hz]
        double FWHM_D = 50e6;
        
        VoigtProfile voigt = new VoigtProfile(nu_0, FWHM_L, FWHM_D);

        // Set frequency array for oscillators. Make sure band limited to
        // resolution of FFT to avoid aliasing (no power at frequencies
        // higher than the Nyquist frequency).
        double nu_min = 0.0e9;
        double nu_max = 1.5e9;   // 1.5Ghz
        double nu_step = 1e5;    // 0.1MHz
        
        int N_OSC = (int)((nu_max - nu_min)/nu_step);
        
        // Array to store magnitude of each oscillator
        Oscillator[] oscillators = new Oscillator[N_OSC];
        
        // Get voigt profile magnitude at each frequency
        for(int i=0; i<N_OSC; i++){
            
            // Translate index i to frequency in Hz
            double nu = (i*nu_step) + nu_step/2.0 + nu_min;
            
            // Get amplitude at this frequency
            double A = Math.sqrt(voigt.getPowerAtFrequency(nu));
            
            // Get random phase angle in range 0->2*pi
            double phi = Math.random() * 2 * Math.PI;
            
            oscillators[i] = new Oscillator(A, nu*2*Math.PI, phi);
            
        }
        
        // Now set sampling rate and frame size
        double t_period = 3.333e-10;  // [seconds]
        int N_SAMPLES = 2048;
        
        // Nyquist frequency. Signal above this frequency will be aliased.
        double f_ny = 1/(2 * t_period);
        
        System.out.println("Nyquist frequency = "+f_ny+" Hz");
        
        // Write samples out to file
        BufferedWriter out = new BufferedWriter(new FileWriter(new File("/home/nickrowell/tmp.txt")));
        
        for(double t=0; t < N_SAMPLES*t_period; t+= t_period){
        
            // Sum amplitude of all oscillators.
            double amplitude = 0;
            
            for(Oscillator oscillator : oscillators){

                amplitude += oscillator.getMagnitude(t);
                
            }
            out.write(amplitude + "\n");

        }
        
        out.close();
    
    
    }
    
    
    
    
    
}

