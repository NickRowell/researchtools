package numeric.triangulation.dim3.infra;

import java.util.LinkedList;
import java.util.List;

import numeric.geom.dim3.Vector3d;

/**
 * Class representing Vertex objects that are to be triangulated.
 *
 * @author nrowell
 * @version $Id$
 */
public class Vertex3d extends Vector3d {
    
    /**
     * List of {@link Tetrahedron} objects that this vertex is connected to.
     */
    public List<Tetrahedron> tetras;
    
    /**
     * List of {@link Triangle3d} objects that this vertex is connected to.
     */
    public List<Triangle3d> tris;
    
    /**
     * Does this {@link Vertex3d} connect to the external hull?
     */
    public boolean EXTERNAL;
    
    public Vertex3d(Vector3d vec3) {
    	super(vec3.getComponents());
    	tetras = new LinkedList<Tetrahedron>();
    	tris = new LinkedList<Triangle3d>();
    	EXTERNAL = false;
    }
    
    /**
     * Associate this Vertex with a Tetrahedron object that it forms
     * part of.
     * @param tetra 
     */
    public void linkTo(Tetrahedron tetra) {
        tetras.add(tetra);
    }
    
    /**
     * Associate this Vertex with a Triangle object that it forms
     * part of. A check is made to avoid repeat adding of the same triangle,
     * in the case where a Triangle forms part of two neighbouring
     * Tetrahedra that share a joining face. The same triangle may be passed 
     * to this method twice when neighbouring tetrahedra are found.
     * 
     * @param tetra 
     */
    public void linkTo(Triangle3d tri) {
        if(!tris.contains(tri)) {
            tris.add(tri);
        }
    }    
    
}