package spectro;

import numeric.complex.Complex;

/**
 * This class represents the Voigt profile line shape.
 * 
 * To do: FIR taps have frequency response that is wider than the corresponding
 *          voigt profile. There appears to be an incorrect term in either the
 *          Lorentz or Gaussian characteristic functions.
 * 
 * @author nickrowell
 */
public class VoigtProfile {
    
    /** Frequency of line centre [Hz]. */
    double nu_0;
    
    /** FWHM of Lorentz broadening component. */
    double FWHM_L;
    
    /** FWHM of Doppler broadening component. */
    double FWHM_D;
    
    /** Doppler broadening component. */
    private DopplerProfile doppler;
    
    /** Lorentz broadening component. */
    private LorentzProfile lorentz;
    
    public VoigtProfile(double _nu_0, double _FWHM_L, double _FWHM_D){
        
        nu_0   = _nu_0;
        FWHM_L = _FWHM_L;
        FWHM_D = _FWHM_D;
        
        doppler = new DopplerProfile(FWHM_D, nu_0);
        lorentz = new LorentzProfile(FWHM_L, nu_0);
        
    }
    
    public LorentzProfile getLorentz(){ return lorentz;}
    public DopplerProfile getDoppler(){ return doppler;}
    
    /**
     * Approximation to the FWHM of the Voigt profile that is accurate to
     * 0.02%. Units are Hz.
     * 
     * @return 
     */
    public double getVoigtFWHM(){
        return getVoigtFWHM(FWHM_L, FWHM_D);
    }
    
    /** Static interface to Voigt FWHM method. */
    public static double getVoigtFWHM(double _FWHM_L, double _FWHM_D){
        return 0.5346 * _FWHM_L + Math.sqrt(0.2166*_FWHM_L*_FWHM_L + _FWHM_D*_FWHM_D);
    }
    
    /**
     * Characteristic function of Voigt profile is product of characteristic
     * functions of each broadening mechanism.
     * 
     * @param t
     * @return 
     */
    public double getCharacteristicFunction(double t){
    
        return doppler.getCharacteristicFunction(t) * 
               lorentz.getCharacteristicFunction(t);
    
    }
    
    /**
     * Normalised absorption, or frequency-domain Voigt profile, is a 
     * convolution of the individual broadening mechanisms.
     * 
     * @param nu    Frequency at which to calculate Voigt profile.
     * @return 
     */
    public double getPowerAtFrequency(double nu){
        
        // Set integration limits based on FWHM of voigt profile.
        double limit = 5 * getVoigtFWHM();
        // Integration limits [Hz] relative to line centre
        double nu_lower = nu-limit;
        double nu_upper = nu+limit;
        
        // Set step size [Hz] as fixed number of steps that is a power of 2
        double h = 2*limit/Math.pow(2,12);
        
        // Use Romberg integration with Richardson's extrapolation to 
        // restrict errors to fourth order.
        double T_h   = 0;    // T(h)
        double T_h_2 = 0;    // T(h/2)
        
        for(double tau = nu_lower; tau < nu_upper; tau += h){
            
            double tau_1   = tau;
            double tau_1p5 = tau+(h/2);
            double tau_2   = tau+h;
            
            // Full step
            T_h += 0.5 * h * (lorentz.getPowerAtFrequency(tau_1)*doppler.getPowerAtFrequency(tau_1-nu+nu_0)
                           +  lorentz.getPowerAtFrequency(tau_2)*doppler.getPowerAtFrequency(tau_2-nu+nu_0));
            
            // Half steps
            T_h_2 += 0.5 * (h/2) * (lorentz.getPowerAtFrequency(tau_1)*doppler.getPowerAtFrequency(tau_1-nu+nu_0)
                                 +  lorentz.getPowerAtFrequency(tau_1p5)*doppler.getPowerAtFrequency(tau_1p5-nu+nu_0));
            
            T_h_2 += 0.5 * (h/2) * (lorentz.getPowerAtFrequency(tau_1p5)*doppler.getPowerAtFrequency(tau_1p5-nu+nu_0)
                                 +  lorentz.getPowerAtFrequency(tau_2)*doppler.getPowerAtFrequency(tau_2-nu+nu_0));
            
        }
        
        // Richardson's extrapolation
        //return (1.0/3.0) * (4.0 * T_h_2  - T_h);
        return T_h;
    }
    
    
    /**
     * Get tap weights for an FIR filter with a frequency response equal to
     * the absorption profile of this Voigt line. It is normal to use only
     * the real parts of these.
     * 
     * @param dt    Sample period of filter.
     * @param N     Number of taps.
     * @return 
     */
    public Complex[] getFIRTaps(double dt, int N){
        
        Complex[] taps = new Complex[N];
        
        for(int n=0; n<N; n++){
        
            // Phasor
            Complex phasor = new Complex(Math.cos(2*Math.PI*nu_0*n*dt),
                                         Math.sin(2*Math.PI*nu_0*n*dt));
            
            // Voigt time domain response
            double voigt = getCharacteristicFunction(n*dt);
            
            
            taps[n] = phasor.times(voigt);
        
        }
        
        return taps;
        
    }
    
    
    
//    public double[] getNormalisedFIRTaps(double dt, int N){
//        
//        double[] taps = getFIRTaps(dt, N);
//        
//        double norm = 0;
//        
//        for(int n=0; n<N; n++) norm += Math.abs(taps[n]);
//        
//        norm *= norm;
//        
//        double[] norm_taps = new double[N];
//        
//        for(int n=0; n<N; n++) norm_taps[n] = taps[n]/norm;
//        
//        return norm_taps;
//        
//    }   
    
    
    
    
    
    
}
