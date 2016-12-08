package spectro;



import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import spectro.Oscillator;
import spectro.VoigtProfile;

/**
 *
 * @author nickrowell
 */
public class OscillatorTester {
    
    public static void main(String[] args) throws IOException{
    
        generateVoigtTimeSeries();
        generateNotchTimeSeries();
    }
    
    
        
    private static void generateVoigtTimeSeries() throws IOException{
    
        // Frequency of line centre [Hz]
        double nu_0 = 800e6;
        // Lorentz width [Hz]
        double FWHM_L = 30e6;
        // Doppler width [Hz]
        double FWHM_D = 50e6;

        // Set frequency array for oscillators. Make sure band limited to
        // resolution of FFT to avoid aliasing (no power at frequencies
        // higher than the Nyquist frequency).
        double nu_min = 0.0e9;
        double nu_max = 1.5e9;   // 1.5Ghz
        double nu_step = 1e5;    // 0.1MHz
        
        // Initialise oscillators and get sum of amplitudes.
        double sum_osc_amps2 = Oscillator
                .initVoigtOscillators(nu_min, nu_max, nu_step, nu_0, FWHM_L, FWHM_D);
        
        // Now set sampling rate and frame size
        double t_period = 3.333e-10;  // [seconds]
        int N_SAMPLES = 2048;
        
        // Nyquist frequency. Signal above this frequency will be aliased.
        double f_ny = 1/(2 * t_period);
        
        System.out.println("Nyquist frequency = "+f_ny+" Hz");
        System.out.println("Sum of oscillator power = "+sum_osc_amps2);
        System.out.println("Voigt FWHM: "+VoigtProfile.getVoigtFWHM(FWHM_L, FWHM_D));
        
        // Write samples out to file
        BufferedWriter out = new BufferedWriter(new FileWriter(
                new File("/home/nickrowell/voigt_time")));
        
        for(double t=0; t < N_SAMPLES*t_period; t+= t_period)
            out.write(t + "\t" + Oscillator.getOscillatorsAmplitude(t) + "\n");
        
        out.close();
    
    
    }
    
        
    private static void generateNotchTimeSeries() throws IOException{
    
        // Frequency of line centre [Hz]
        double nu_0 = 800e6;
        // Notch full width
        double FWHM = 100e6;

        // Set frequency array for oscillators. Make sure band limited to
        // resolution of FFT to avoid aliasing (no power at frequencies
        // higher than the Nyquist frequency).
        double nu_min = 0.0e9;
        double nu_max = 1.5e9;   // 1.5Ghz
        double nu_step = 1e5;    // 0.1MHz
        
        // Initialise oscillators and get sum of amplitudes.
        double sum_osc_amps2 = Oscillator
                .initNotchOscillators(nu_min, nu_max, nu_step, nu_0, FWHM);
        
        // Now set sampling rate and frame size
        double t_period = 3.333e-10;  // [seconds]
        int N_SAMPLES = 2048;
        
        // Nyquist frequency. Signal above this frequency will be aliased.
        double f_ny = 1/(2 * t_period);
        
        System.out.println("Nyquist frequency = "+f_ny+" Hz");
        System.out.println("Sum of oscillator power = "+sum_osc_amps2);
        
        // Write samples out to file
        BufferedWriter out = new BufferedWriter(new FileWriter(
                new File("/home/nickrowell/notch_time")));
        
        for(double t=0; t < N_SAMPLES*t_period; t+= t_period)
            out.write(t + "\t" + Oscillator.getOscillatorsAmplitude(t) + "\n");
        
        out.close();
    
    
    }
        
}
