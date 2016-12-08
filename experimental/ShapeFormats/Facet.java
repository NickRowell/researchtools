package ShapeFormats;

import java.util.LinkedList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Scanner;

import shape.Vertex;

/**
 * Class represents triangular polygons expressed as the indices of three
 * vertices.
 * 
 * @author nickrowell
 */
public class Facet 
{
    
    // Facet number
    int N;    
    
    // Numbers of each vertex.
    int v0, v1, v2;
    
    // Components of surface normal
    double[] n;
    
    // Main constructor
    public Facet(int _N, int _v0, int _v1, int _v2,
                 Vertex r0, Vertex r1, Vertex r2)
    {
        N  = _N;
        
        v0 = _v0;
        v1 = _v1;
        v2 = _v2;
        
        n = new double[3];
        
        // Vector from v0 to v2
        double[] a = {r2.x-r0.x, r2.y-r0.y, r2.z-r0.z};
        // Vector from v0 to v1
        double[] b = {r1.x-r0.x, r1.y-r0.y, r1.z-r0.z};
        
        // Cross product a x b gives normal direction
        n[0] = (a[1]*b[2] - a[2]*b[1]);
        n[1] = (a[2]*b[0] - a[0]*b[2]);
        n[2] = (a[0]*b[1] - a[1]*b[0]);
        
        // Normalise vector.
        double norm = Math.sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        
        n[0] /= norm;
        n[1] /= norm;
        n[2] /= norm;
        
        r0.n0 += n[0]; r0.n1 += n[1]; r0.n2 += n[2];
        r1.n0 += n[0]; r1.n1 += n[1]; r1.n2 += n[2];
        r2.n0 += n[0]; r2.n1 += n[1]; r2.n2 += n[2];
        
    }
    
    // Partial constructor
    public Facet(int _v0, int _v1, int _v2)
    {
        v0 = _v0;
        v1 = _v1;
        v2 = _v2;
    }
    
    // Copy constructor
    public Facet(Facet copyme)
    {
        v0 = copyme.v0;
        v1 = copyme.v1;
        v2 = copyme.v2;
    }
    
    
    public Facet(List<Integer> vert12, List<Integer> vert3)
    {
        if(vert12.size()!=2) throw new RuntimeException("Expected 2 vertices, found "+vert12.size());
        if(vert3.size()!=1) throw new RuntimeException("Expected 1 vertex, found "+vert3.size());
        
        v0 = vert12.get(0);
        v1 = vert12.get(1);
        v2 = vert3.get(0);
    }
    
    
    /** Count how many vertices this Facet shares with that facet. */
    public int numberVertsInCommon(Facet compare)
    {
        return getVertsInCommon(compare).size();
    }
    
    /** Check if this Facet shares two vertices with that facet. */
    public boolean hasTwoVertsInCommon(Facet compare)
    {
        return getVertsInCommon(compare).size() == 2;
    }
    
    /** Check if this Facet shares no vertices with that facet. */
    public boolean hasNoVertsInCommon(Facet compare)
    {
        return getVertsInCommon(compare).size() == 0;
    }    
    
    
    /** Check if this Facet shares one vertex with that facet. */
    public boolean hasOneVertInCommon(Facet compare)
    {
        return getVertsInCommon(compare).size() == 1;
    }    
    
    
    
    /** Get a List of the vertices that these facets share. */
    public List<Integer> getVertsInCommon(Facet compare)
    {
        
        List<Integer> verts = new LinkedList<Integer>();
        
        // Compare vertex 1...
        if(v0 == compare.v0 || v0 == compare.v1 || v0 == compare.v2)
            verts.add(v0);
        // Compare vertex 2...
        if(v1 == compare.v0 || v1 == compare.v1 || v1 == compare.v2)
            verts.add(v1);
        // Compare vertex 3...
        if(v2 == compare.v0 || v2 == compare.v1 || v2 == compare.v2)
            verts.add(v2);
        
        return verts;
        
    }
    
    /** Get a List of the vertices in this facet that aren't part of that facet. */
    public List<Integer> getVertsNotInCommon(Facet compare)
    {
        
        List<Integer> verts = new LinkedList<Integer>();
        
        // Compare vertex 1...
        if(v0 != compare.v0 && v0 != compare.v1 && v0 != compare.v2)
            verts.add(v0);
        // Compare vertex 2...
        if(v1 != compare.v0 && v1 != compare.v1 && v1 != compare.v2)
            verts.add(v1);
        // Compare vertex 3...
        if(v2 != compare.v0 && v2 != compare.v1 && v2 != compare.v2)
            verts.add(v2);
        
        return verts;
        
    }
    
    
    public static Facet parseFacet(String str, List<Vertex> verts)
    {
                
        // Open scanner on String
        Scanner scan = new Scanner(str);
        
        try
        {
            int N  = scan.nextInt();
            int v0 = scan.nextInt();
            int v1 = scan.nextInt();
            int v2 = scan.nextInt();
            
            // Create new Facet
            return new Facet(N, v0, v1, v2, verts.get(v0-1), verts.get(v1-1),verts.get(v2-1));
            
        }
        catch(NoSuchElementException nsee)
        {
            throw new RuntimeException("Missing integer on facet "+str);
        }
        catch(NumberFormatException nfe)
        {
            throw new RuntimeException("Number format exception on facet "+str);
        }
    
    }
    
    
    
    public static double[] getSurfaceNormal(Vertex r0, Vertex r1, Vertex r2)
    {
        // Vector from v0 to v2
        double[] a = {r2.x-r0.x, r2.y-r0.y, r2.z-r0.z};
        // Vector from v0 to v1
        double[] b = {r1.x-r0.x, r1.y-r0.y, r1.z-r0.z};
        
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
