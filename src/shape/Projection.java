package shape;

/**
 * Parent class that defines interfaces for all map projection objects.
 * 
 * The class provides functions that map back and forth between latitude
 * and longitude, the projection coordinates, and the array indices of samples
 * in the data.
 * 
 * 
 * 
 * @author nickrowell
 */
public abstract class Projection
{
    public abstract int[] getLineSampleIndices(double lon, double lat);
    
    
    public abstract double[] getLongLat(int l, int s);
    
}
