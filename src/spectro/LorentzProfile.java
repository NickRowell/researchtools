package spectro;

/**
 * This class represents the Lorentz broadening mechanism. The maths is taken
 * from the paper
 * 
 * "An Analytical Derivation of a Popular Approximation of the Voigt Function
 *  for Quantification of NMR Spectra"
 *  Bruce, S., Higinbotham, J., Marshall, I., Beswick, P.
 *  Journal of Magnetic Resonance 142, 57-63 (2000)
 * 
 * @author nickrowell
 */
public class LorentzProfile {
    
    
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
    public LorentzProfile(double _FWHM, double _nu_0){
    
        FWHM = _FWHM;
        nu_0 = _nu_0;
    
    }
    
    /**
     * Get the normalised absorption at frequency nu due to Lorentz broadening.
     * 
     * @param nu
     * @return 
     */
    public double getPowerAtFrequency(double nu){
        
        // Intermediate variable...
        double alpha = (nu - nu_0)/(FWHM/2.0);
        // ...squared
        double alpha2 = alpha * alpha;
        // Get complete function
        return (2.0/(Math.PI*FWHM))*(1/(1+alpha2));
        
    }
    
    
    /**
     * Get characteristic function of Doppler line shape. This is in the time
     * domain.
     * 
     * @param t
     * @return 
     */
    public double getCharacteristicFunction(double t){
        
        // Only valid for t>=0
        assert t>=0;
        
        // Normalisation factor
        double norm = Math.PI * FWHM;
        
        return norm * Math.exp(-Math.PI * FWHM * t);
        
        
    }
    
}
