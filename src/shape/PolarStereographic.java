package shape;

/**
 *
 * @author nickrowell
 */
public class PolarStereographic extends Projection
{
    
    public static enum HEMI{NORTH, SOUTH};
    
    /** Resolution of map [m/pixel]. */
    double MAP_SCALE;
    /** Lines. */
    double LINES;
    /** Radius of body [m]. */
    double RADIUS;
    /** Hemisphere (determines deprojection). */
    HEMI hemi;
    
    
    public PolarStereographic(double pMAP_SCALE,
                              double pRADIUS,
                              double pLINES,
                              HEMI phemi)
    {
        MAP_SCALE = pMAP_SCALE;
        RADIUS    = pRADIUS;
        LINES     = pLINES;
        hemi      = phemi;
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
        // Distance from tangent point to sample, in projection plane
        double R = 2.0 * RADIUS * 1000.0 * Math.tan((Math.PI/2 - Math.abs(lat))/2.0);
        
        // line/samples origin in lower left:
        
        // Components of sample position in projection plane [m]
        double X = R * Math.sin(lon);
        double Y = R * Math.cos(lon);
        // Round off
        int I = (int)Math.rint(X/MAP_SCALE + LINES/2 + 0.5);
        int J = (int)Math.rint(-Y/MAP_SCALE + LINES/2 + 0.5);
        
        
        
        
        
        return new int[]{I,J};
    }

    @Override
    public double[] getLongLat(int I, int J)
    {
        
        // Line/samples origin in lower left:
        double X = (I - LINES/2 - 0.5)*MAP_SCALE;
        double Y = -(J - LINES/2 - 0.5)*MAP_SCALE;
        double LON = Math.atan2(X,Y);
        
        
        
        
        
        
        double R = Math.sqrt(X*X + Y*Y);
        double LAT = 0;
        
        switch(hemi)
        {
            case NORTH: LAT =  Math.PI/2.0 - 2*Math.atan(0.5 * R/(RADIUS*1000)); break;
            case SOUTH: LAT = -Math.PI/2.0 + 2*Math.atan(0.5 * R/(RADIUS*1000)); break;
        }
        
        return new double[]{LON,LAT};
        
    }
}
