package spectro;

/**
 * Class represents a sine wave oscillator.
 * 
 * @author nickrowell
 */
public class Oscillator {
    
    
    /** Static array of oscillators that can be accessed from eg Matlab. */
    public static Oscillator[] oscillators;
    
    
    /** Amplitude. */
    double A;
    
    /** Frequency [rad/s]. */
    double w;
    
    /** Phase [rad]. */
    double phi;
    
    
    public Oscillator(double _A, double _w, double _phi){
    
        A = _A;
        w = _w;
        phi = _phi;
    
    }
    
    /**
     * Get amplitude of wave at time t [seconds].
     * @param t Time [seconds].
     * @return 
     */
    public double getMagnitude(double t){
    
        return A * Math.sin(getPhase(t));
    
    }
    
    /**
     * Get phase of wave at time t.
     * @param t
     * @return 
     */
    public double getPhase(double t){
    
        return w * t + phi;
    
    }
    
        
    public static void clearOscillators(){ oscillators = null;}
    
    public static double getOscillatorsAmplitude(double t){
        
        // Sum amplitude of all oscillators.
        double amplitude = 0;
            
        for(Oscillator oscillator : oscillators)
            amplitude += oscillator.getMagnitude(t);
                
        return amplitude;
    }
    
    /**
     * Populate oscillator array to provide a frequency power spectrum that
     * shows a notch.
     * 
     * @param nu_min    Minimum frequency in oscillator array [Hz]
     * @param nu_max    Maximum frequency in oscillator array [Hz]
     * @param dnu       Frequency step between oscillators [Hz]
     * @param nu_0      Frequency of notch centre [Hz]
     * @param FWHM      Notch full width
     */
    public static double initNotchOscillators(double nu_min, double nu_max, 
                                            double dnu, double nu_0,
                                            double FWHM){
        
        // Number of oscillators
        int N_OSC = (int)((nu_max - nu_min)/dnu);
        
        // Initialise array...
        oscillators = new Oscillator[N_OSC];
        
        // Calculate sum of oscillator power
        double SUM_AMP2 = 0;
        
        // Get magnitude at each frequency
        for(int i=0; i<N_OSC; i++){
            
            // Translate index i to frequency in Hz
            double nu = (i*dnu) + dnu/2.0 + nu_min;
            
            // Initialise oscillator amplitude
            double A = 1.0;
            
            // Is this frequency inside notch?
            if(Math.abs(nu - nu_0) < FWHM/2.0)
                // Yes - set amplitude of this component to zero
                A = 0.0;
            
            // Get random phase angle in range 0->2*pi
            double phi = Math.random() * 2 * Math.PI;
            
            // Initialise oscillator at this frequency
            oscillators[i] = new Oscillator(A, nu*2*Math.PI, phi);
            
            SUM_AMP2 += (A*A);
            
        }              
        
        return SUM_AMP2;
    
    }
    
    
    
    
    /**
     * Method populates oscillator array with a number of random-phase, varying
     * frequency sine waves with amplitudes that have a Voigt profile.
     * 
     * 
     * @param nu_min    Minimum frequency in oscillator array [Hz]
     * @param nu_max    Maximum frequency in oscillator array [Hz]
     * @param dnu       Frequency step between oscillators [Hz]
     * @param nu_0      Frequency of line centre [Hz]
     * @param FWHM_L    FWHM of Lorentz broadening mechanism [Hz]
     * @param FWHM_D    FWHM of Doppler broadening mechanism [Hz] 
     */
    public static double initVoigtOscillators(double nu_min, double nu_max, 
                                            double dnu, double nu_0,
                                            double FWHM_L, double FWHM_D){
        
        // Get Voigt broadening mechanism.
        VoigtProfile voigt = new VoigtProfile(nu_0, FWHM_L, FWHM_D);
        
        // Number of oscillators
        int N_OSC = (int)((nu_max - nu_min)/dnu);
        
        // Initialise array...
        oscillators = new Oscillator[N_OSC];
        
        // Calculate sum of oscillator amplitudes
        double SUM_AMP2 = 0;
        
        // Get voigt profile magnitude at each frequency
        for(int i=0; i<N_OSC; i++){
            
            // Translate index i to frequency in Hz
            double nu = (i*dnu) + dnu/2.0 + nu_min;
            
            // Get amplitude of voigt profile at this frequency
            double A = Math.sqrt(voigt.getPowerAtFrequency(nu));
            
            // Get random phase angle in range 0->2*pi
            double phi = Math.random() * 2 * Math.PI;
            
            // Initialise oscillator at this frequency
            oscillators[i] = new Oscillator(A, nu*2*Math.PI, phi);
            
            SUM_AMP2 += (A);
            
        }       
        
        return SUM_AMP2;
        
    }
    
    

    
    
    
}
