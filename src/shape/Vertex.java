package shape;

/**
 * Instances of this class represent a single vertex.
 * @author nickrowell
 */
public class Vertex
{
    
    // Vertex number
    public int N;
    
    // Position coordinates.
    public double x,y,z;
    
    // Components of surface normal.
    public double n0,n1,n2;
    
    /** Main constructor. */
    public Vertex(int _N, double _x, double _y, double _z)
    {
        N = _N;
        x = _x;
        y = _y;
        z = _z;
        
        n0=n1=n2=0;
    }
    
    /** Constructor that doesn't use vertex number. */
    public Vertex(double _x, double _y, double _z)
    {
        this(0, _x, _y, _z);
    }    
    
    
    
    
    @Override
    public String toString()
    {
        return String.format("%f %f %f",x,y,z);
    }
    
    /** Get vector representing surface normal. */
    public String getNormalString()
    {
        
        // Normalise components of surface normal
        double norm = Math.sqrt(n0*n0 + n1*n1 + n2*n2);
        
        // Check that at normal vector isn't zero
        if(norm < 1e-9) return " 0 0 0 ";
        
        return (n0/norm) + "\t" + (n1/norm) + "\t" + (n2/norm);        
        
    }
    
    /**
     * Calculates the surface normal for the facet produced by joining vertices
     * r0, r1 and r2 in clockwise winding order.
     * @param r0
     * @param r1
     * @param r2
     * @return 
     */
    public static double[] getClockwiseSurfaceNormal(Vertex r0, Vertex r1, Vertex r2)
    {
        // Vector from v0 to v1
        double[] a = {r1.x-r0.x, r1.y-r0.y, r1.z-r0.z};
        // Vector from v0 to v2
        double[] b = {r2.x-r0.x, r2.y-r0.y, r2.z-r0.z};
        
        // Cross product a x b gives normal direction
        double n0 = (a[1]*b[2] - a[2]*b[1]);
        double n1 = (a[2]*b[0] - a[0]*b[2]);
        double n2 = (a[0]*b[1] - a[1]*b[0]);
        
        // Normalise vector.
        double norm = Math.sqrt(n0*n0 + n1*n1 + n2*n2);
        
        n0 /= norm;
        n1 /= norm;
        n2 /= norm; 
    
        return new double[]{n0,n1,n2};
        
    }

}
