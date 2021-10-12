package shape;

/**
 * This class implements the conversion between image coordinates and Mars areocentric coordinates
 * for a calibrated HiRISE RDR image. Some notes:
 * 
 *  - HiRISE is a pushbroom camera and the image coordinates are line (in the direction of motion) 
 *    and sample (transverse to that).
 * 
 * @author nickrowell
 */
public class Equirectangular extends Projection {
	
	/**
	 * Projection offset parameter for the line coordinate [pixels].
	 */
    public final double LINE_PROJECTION_OFFSET;
    
    /**
	 * Projection offset parameter for the sample coordinate [pixels].
	 */
    public final double SAMPLE_PROJECTION_OFFSET;
    
    /**
     * Map scale parameter [metres/pixel].
     */
    public final double MAP_SCALE;
    
    /**
     * Center of projection, longitude parameter [radians].
     */
    public final double CENTER_LONGITUDE;
    
    /**
     * Center of projection, latitude parameter [radians].
     */
    public final double CENTER_LATITUDE;
    
    /**
     * Mars equatorial radius parameter [km].
     */
    public final double A_AXIS_RADIUS;
    
    /**
     * Mars equatorial radius parameter [km].
     */
    public final double C_AXIS_RADIUS;
    
    /**
     * Main constructor accepting the calibration coefficients necessary to transform between
     * (sample,line) and (longitude,latitude) coordinates.
     * 
     * @param LINE_PROJECTION_OFFSET
	 * 	Projection offset parameter for the line coordinate [pixels].
     * @param SAMPLE_PROJECTION_OFFSET
	 * 	Projection offset parameter for the sample coordinate [pixels].
     * @param MAP_SCALE
     * 	Map scale parameter [metres/pixel].
     * @param CENTER_LONGITUDE
     * 	Center of projection, longitude parameter [degrees].
     * @param CENTER_LATITUDE
     * 	Center of projection, latitude parameter [degrees].
     * @param A_AXIS_RADIUS
     * 	Mars equatorial radius parameter [km].
     * @param C_AXIS_RADIUS
     * 	Mars equatorial radius parameter [km].
     */
    public Equirectangular(double LINE_PROJECTION_OFFSET, double SAMPLE_PROJECTION_OFFSET,
                           double MAP_SCALE, double CENTER_LONGITUDE, double CENTER_LATITUDE,
                           double A_AXIS_RADIUS, double C_AXIS_RADIUS) {
    	
        this.LINE_PROJECTION_OFFSET   = LINE_PROJECTION_OFFSET;
        this.SAMPLE_PROJECTION_OFFSET = SAMPLE_PROJECTION_OFFSET;
        this.MAP_SCALE                = MAP_SCALE;
        this.CENTER_LONGITUDE         = Math.toRadians(CENTER_LONGITUDE);
        this.CENTER_LATITUDE          = Math.toRadians(CENTER_LATITUDE);
        this.A_AXIS_RADIUS            = A_AXIS_RADIUS;
        this.C_AXIS_RADIUS            = C_AXIS_RADIUS;
    }

    /**
     * Calculate indices of sample closest to the specified latitude and longitude coordinates.
     * 
     * @param lon
     * 	Longitude [radians]
     * @param lat
     * 	Latitude [radians]
     * @return 
     * 	A two-element array containing the corresponding sample and line coordinates.
     */
    public int[] getLineSampleIndices(double lon, double lat) {
    	
        // Local radius of Mars
        double clat0 = Math.cos(CENTER_LATITUDE);
        double slat0 = Math.sin(CENTER_LATITUDE);
        double Re    = A_AXIS_RADIUS;
        double Rp    = C_AXIS_RADIUS;
        double R     = Re * Rp / Math.sqrt(Rp*Rp*clat0*clat0 + Re*Re*slat0*slat0);
        
        // Map coordinates of desired centre [units ~metres for small angles]
        //
        // NOTE: map projection info is wrong in its description - the angular
        //       coordinates (lat & lon) should be in radians.
        
        // Map projection coordinates
        double xm = R * (lon - CENTER_LONGITUDE) * clat0;
        double ym = R * (lat);
        
        // Sample indices at floating point precision
        double s =  xm / (MAP_SCALE/1000) + SAMPLE_PROJECTION_OFFSET + 1;
        double l = -(ym / (MAP_SCALE/1000) - LINE_PROJECTION_OFFSET + 1);
        
        // Round off
        int is =  (int)Math.rint(s);
        int il =  (int)Math.rint(l);
    
        return new int[]{is,il};
    }

    /**
     * Calculate the longitude and latitude corresponding to the given pixel coordinate. 
     * 
     * @param s
     * 	The sample coordinate [pixels].
     * @param l
     * 	The line coordinate [pixels].
     * @return
     * 	A two-element array containing the corresponding longitude and latitude coordinates [radians].
     */
    public double[] getLongLat(int s, int l) {
    	
        // Local radius of Mars
        double clat0 = Math.cos(CENTER_LATITUDE);
        double slat0 = Math.sin(CENTER_LATITUDE);
        double Re = A_AXIS_RADIUS;
        double Rp = C_AXIS_RADIUS;
        double R  = Re * Rp / Math.sqrt(Rp*Rp*clat0*clat0 + Re*Re*slat0*slat0);
        
        double xm =  (s - SAMPLE_PROJECTION_OFFSET - 1)*(MAP_SCALE/1000);
        double ym = -(l - LINE_PROJECTION_OFFSET   + 1)*(MAP_SCALE/1000);
        
        double lon = xm/(R*clat0) + CENTER_LONGITUDE;
        double lat = ym/R;
        
        return new double[]{lon,lat};
        
    }
    
    /**
     * {@inheritDoc}
     */
    public String toString() {
    	
    	StringBuilder sb = new StringBuilder();
    	sb.append("[LINE_PROJECTION_OFFSET = " + LINE_PROJECTION_OFFSET + "; ");
    	sb.append("SAMPLE_PROJECTION_OFFSET = " + SAMPLE_PROJECTION_OFFSET + "; ");
    	sb.append("MAP_SCALE = " + MAP_SCALE + "; ");
    	sb.append("CENTER_LONGITUDE = " + CENTER_LONGITUDE + "; ");
    	sb.append("CENTER_LATITUDE = " + CENTER_LATITUDE + "; ");
    	sb.append("A_AXIS_RADIUS = " + A_AXIS_RADIUS + "; ");
    	sb.append("C_AXIS_RADIUS = " + C_AXIS_RADIUS + "]");
    	
    	return sb.toString();
    }
}
