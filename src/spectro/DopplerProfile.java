package spectro;

/**
 * This class represents the Doppler broadening mechanism. The maths is taken
 * from the paper
 * 
 * "An Analytical Derivation of a Popular Approximation of the Voigt Function
 *  for Quantification of NMR Spectra"
 *  Bruce, S., Higinbotham, J., Marshall, I., Beswick, P.
 *  Journal of Magnetic Resonance 142, 57-63 (2000)
 * 
 * @author nickrowell
 */
public class DopplerProfile {
    
    /** Full width at half-maximum [Hz]. */
    double FWHM;
    
    /** Line centre [Hz]. */
    double nu_0;
    
    /**
     * Main constructor.
     * 
     * @param _FWHM
     * @param _nu_0 
     */
    public DopplerProfile(double _FWHM, double _nu_0){
    
        FWHM = _FWHM;
        nu_0 = _nu_0;
    
    }
    
    /**
     * Get the normalised absorption at frequency nu due to Doppler broadening.
     * 
     * @param nu
     * @return 
     */
    public double getPowerAtFrequency(double nu) {
    
        // Exponent of Gaussian...
        double exponent = (nu - nu_0)/(FWHM/(2 * Math.sqrt(Math.log(2))));
        // ...squared
        double exponent2 = exponent * exponent;
        // Get complete function.
        return (2/FWHM) * Math.sqrt(Math.log(2)/Math.PI) * Math.exp(-exponent2);
    }
    
    /**
     * Get characteristic function of Doppler line shape. This is in the time
     * domain. Doppler broadening profile is Gaussian, so characteristic
     * function (which is the Fourier transform of the line profile) is also
     * a Gaussian.
     * 
     * Basic formula for Fourier transform of a Gaussian:
     * 
     * FFT(exp(-a*x^2)) = sqrt(pi/a) * exp(- pi^2 * K^2 / a)
     * 
     * 
     * 
     * @param t
     * @return 
     */
    public double getCharacteristicFunction(double t) {
        
        // Get exponent
        double exponent = Math.PI * Math.PI * FWHM * FWHM * t * t / (4 * Math.log(2));
        
        // Get complete function
        return Math.exp(-exponent);
    }
}