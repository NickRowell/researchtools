package numeric.triangulation.dim3.infra;

import java.util.LinkedList;
import java.util.List;
import numeric.geom.dim3.Vector3d;

/**
 * Extension of the base class {@link numeric.geom.dim3.Triangle3d} to add some additional fields
 * in order to support Delaunay triangulation.
 * 
 * This class uses {@link Vertex3d} as the type parameter, as Vertex contains fields supporting
 * Delaunay triangulation.
 *
 * @author nrowell
 * @version $Id$
 */
public class Triangle3d extends numeric.geom.dim3.Triangle3d<Vertex3d>
{
    /**
     * List of tetrahedra this Triangle is part of.
     */
    public List<Tetrahedron> tetras;
    
    /**
     * Main constructor.
     * 
     * @param v0
     * 	The first vertex.
     * @param v1
     * 	The second vertex.
     * @param v2
     * 	The third vertex.
     */
    public Triangle3d(Vertex3d v0, Vertex3d v1, Vertex3d v2)
    {
    	super(v0, v1, v2);
    	tetras = new LinkedList<Tetrahedron>();
    }
    
    /**
     * Check if triangle lies on external hull of triangulation. We can test this simply by checking
     * if the triangle is connected to two tetrahedra, which means it forms the joining face of them
     * and is therefore internal.
     * @return 
     */
    public boolean isExternal()
    {
        if(tetras.size()==1) {
        	return true;
        }
        else if(tetras.size()==2) {
        	return false;
        }
        else {
        	throw new RuntimeException("Queried Triangle part of "+tetras.size()+" Tetrahedrons");
        }
    }
    
    /**
     * Associate this triangle with a Tetrahedron object that it forms
     * part of.
     * 
     * Also determines if this Triangle is part of external hull:
     * 
     * In a Delaunay triangulation, each triangle is connected to at most two
     * tetrahedra. If a triangle is connected to one, then we know that it
     * forms part of the external hull. If it's connected to two, then we know
     * it's an interior triangle because it forms the joining boundary of two
     * closed tetrahedra.
     * 
     * @param tetra 
     */
    public void linkTo(Tetrahedron tetra)
    {
        tetras.add(tetra);
    }
    
    /**
     * Method obtains the surface normal for the triangle, defined so that it
     * points outwards from the tetrahedron that this triangle forms part of.
     * Note that this is only defined for triangles that are only attached to
     * one Tetrahedron.
     */
    public Vector3d getNormal()
    {
        // Should never call this method for triangles that have not previously
        // been identified as EXTERNAL, so we should always find only one
        // adjoining Tetrahedron at this point.
        if(tetras.size()!=1) {
            throw new RuntimeException("Triangle attached to "+tetras.size()
                    +" tetrahedra!");
        }
        
        return tetras.get(0).getOutwardNormalForTriangle(this);
    }
    
    /**
     * Get the distance to the closest point in the Triangle to the point P.
     */
    public double getMinDistanceToP(Vertex3d p)
    {
        return p.minus(getClosestPointToP(p)).norm();
    }

}