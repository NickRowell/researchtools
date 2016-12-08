package shape;

/**
 *
 * @author nickrowell
 */
public class Equirectangular extends Projection
{
    // Offset of image array origin from projection centre [pixels]
    double LINE_PROJECTION_OFFSET;
    double SAMPLE_PROJECTION_OFFSET;
    
    // Scale of map [metres/pixel]
    double MAP_SCALE;
    
    // Centre of projection [radians]
    double CENTER_LONGITUDE;
    double CENTER_LATITUDE;
    
    // Equatorial and polar radii [km]
    double A_AXIS_RADIUS;
    double C_AXIS_RADIUS;
    
    
    public Equirectangular(double pLINE_PROJECTION_OFFSET, double pSAMPLE_PROJECTION_OFFSET,
                           double pMAP_SCALE, double pCENTER_LONGITUDE, double pCENTER_LATITUDE,
                           double pA_AXIS_RADIUS, double pC_AXIS_RADIUS)
    {
        LINE_PROJECTION_OFFSET   = pLINE_PROJECTION_OFFSET;
        SAMPLE_PROJECTION_OFFSET = pSAMPLE_PROJECTION_OFFSET;
        MAP_SCALE                = pMAP_SCALE;
        CENTER_LONGITUDE         = Math.toRadians(pCENTER_LONGITUDE);
        CENTER_LATITUDE          = Math.toRadians(pCENTER_LATITUDE);
        A_AXIS_RADIUS            = pA_AXIS_RADIUS;
        C_AXIS_RADIUS            = pC_AXIS_RADIUS;
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

    @Override
    public double[] getLongLat(int s, int l)
    {
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
}
