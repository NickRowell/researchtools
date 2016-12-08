package shape;

/**
 *
 * @author nickrowell
 */
public class Sinusoidal extends Projection
{
    // Offset of image array origin from projection centre [pixels]
    double LINE_PROJECTION_OFFSET;
    double SAMPLE_PROJECTION_OFFSET;
    
    // Resolution of map [pixels/degree]
    double MAP_RESOLUTION;
    
    // Centre of projection [degrees]
    double CENTER_LONGITUDE;
    double CENTER_LATITUDE;
    
    public Sinusoidal(double pLINE_PROJECTION_OFFSET, double pSAMPLE_PROJECTION_OFFSET, 
                      double pMAP_RESOLUTION, double pCENTER_LONGITUDE, double pCENTER_LATITUDE)
    {
        LINE_PROJECTION_OFFSET = pLINE_PROJECTION_OFFSET;
        SAMPLE_PROJECTION_OFFSET = pSAMPLE_PROJECTION_OFFSET;
        MAP_RESOLUTION = pMAP_RESOLUTION;
        CENTER_LONGITUDE = pCENTER_LONGITUDE;
        CENTER_LATITUDE  = pCENTER_LATITUDE;
    }
    
    /**
     * Calculate indices of sample closest to the specified latitude & 
     * longitude.
     * @param lon   Longitude [radians]
     * @param lat   Latitude  [radians]
     * @return 
     */
    public int[] getLineSampleIndices(double lon, double lat)
    {
        
        // Sample indices at floating point precision
        double s0 = SAMPLE_PROJECTION_OFFSET + MAP_RESOLUTION*(Math.toDegrees(lon) - CENTER_LONGITUDE)*Math.cos(lat);
        double l0 = LINE_PROJECTION_OFFSET   - MAP_RESOLUTION*Math.toDegrees(lat);
        
        // Round off
        int is0 = (int)Math.rint(s0);
        int il0 = (int)Math.rint(l0);
        
        return new int[]{is0,il0};
    }

    @Override
    public double[] getLongLat(int l, int s) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}
