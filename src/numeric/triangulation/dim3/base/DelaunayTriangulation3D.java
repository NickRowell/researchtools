package numeric.triangulation.dim3.base;

import java.util.LinkedList;
import java.util.List;

import numeric.geom.dim3.Sphere;
import numeric.triangulation.dim3.infra.Tetrahedron;
import numeric.triangulation.dim3.infra.Triangle3d;
import numeric.triangulation.dim3.infra.Vertex3d;

/**
 * Brute-force Delaunay triangulation algorithm implementation.
 * 
 * There are three basic objects that define the triangulation: Vertex,
 * Triangle and Tetrahedron. Only Vertex object contain data: these contain
 * the three coordinates defining the vertex 3D position. Each object contains
 * two lists; these hold references to objects of each of the other two types.
 * For example, Vertex objects hold a list of references to Triangles and 
 * Tetrahedron objects that they form part of in the final triangulation.
 * Triangle objects hold a list of (3) Vertex objects that define their corners
 * and a list of (one or two) Tetrahedra that they form part of.
 * 
 * The objects are thus divided into a geometrical hierarchy:
 * 
 * Vertex -> Triangle -> Tetrahedron
 * 
 * where each object holds references to the objects that define it (backwards)
 * and the objects that it forms part of (forwards). Only references are stored
 * so that the same object in different lists can be identified.
 * 
 * The main class DelaunayTriangulation3D contains the overall lists of each 
 * object in the triangulation, and does the work of identifying the
 * connections between them, starting out with just the Vertex objects.
 *
 * @author nrowell
 * @version $Id$
 */
public abstract class DelaunayTriangulation3D {
    
    /**
     * Handle to list of all Vertices to be triangulated.
     */
    protected List<Vertex3d> vertices;
    
    /**
     * Number of external vertices, i.e. those lying on the outer hull.
     */
    protected int N_VERT_EXTERNAL;
    
    /**
     * List of all Triangles comprising the triangulation.
     */
    protected List<Triangle3d> tris;
    
    /** 
     * Number of external triangles. An external triangle is one which forms
     * part of the external hull. These can be identified because they form
     * part of only one tetrahedron. Internal triangles form part of two 
     * tetrahedra: the face where two tetrahedra meet.
     */
    protected int N_TRI_EXTERNAL;
    
    /**
     * List of all Tetrahedra that form part of triangulation.
     */
    protected List<Tetrahedron> tetras;
       
    /**
     * Number of external tetrahedra, i.e those lying on the outer hull.
     */
    protected int N_TETRA_EXTERNAL;
    
    /**
     * Main constructor.
     * 
     * @param pvertices
     * 	The {@link List} of {@link Vertex3d} to triangulate.
     */
    public DelaunayTriangulation3D(List<Vertex3d> pvertices) {
    	
        // Copy references to existing Vertex list
        vertices         = pvertices;
        tris             = new LinkedList<Triangle3d>();
        tetras           = new LinkedList<Tetrahedron>();
        N_VERT_EXTERNAL  = 0;
        N_TRI_EXTERNAL   = 0;
        N_TETRA_EXTERNAL = 0;
        
        // For all Vertices, clear internal lists of connected Triangle and
        // Tetrahedron objects that may have contents from a previous
        // triangulation; reset all to EXTERNAL.
        for(Vertex3d v : vertices)
        {
            v.EXTERNAL = true;
            v.tetras.clear();
            v.tris.clear();
        }
    }
    
    /**
     * Instance of 3D Delaunay Triangulation algorithm that operates on the
     * internal List of Clusters, and writes status messages etc to the
     * PipedOutputStream passed in as argument.
     * @param _clusters
     * @param os
     */
    public abstract void doTriangulation();
    
    
    /**
     * 
     * @param v1
     * @param v2
     * @param v3
     * @param v4
     * @param sphere
     */
    protected void addTetrahedron(Vertex3d v1, Vertex3d v2, Vertex3d v3, Vertex3d v4, Sphere sphere) {
        
        // Each subset of three 
        // vertices forms a triangle in the final 
        // delaunay triangulation, and the set of four
        // forms a 'simplex', in this case a tetrahedron
        // that doesn't intersect with any others.

        // The connectedness of the points is stored in 
        // terms of the triangles they form. 
        // Some or all of these
        // may already exist however from triangulations
        // with neighbouring points, so check for this before
        // adding these to the main List. The check is
        // done in getUniqueTriangle(Triangle) method, 
        // which either returns a reference to an existing
        // triangle, or adds a new triangle to the list and
        // returns a reference to it.

        // Four triangles formed by connecting the four
        // vertices. Add these to the list held by the
        // class, taking care to avoid adding multiple
        // instances of the same triangle lying on
        // shared face of neighbouring tetrahedra.
        Triangle3d v123 = getUniqueTriangle(new Triangle3d(v1,v2,v3));
        Triangle3d v124 = getUniqueTriangle(new Triangle3d(v1,v2,v4));
        Triangle3d v134 = getUniqueTriangle(new Triangle3d(v1,v3,v4));
        Triangle3d v234 = getUniqueTriangle(new Triangle3d(v2,v3,v4));

        // These triangles and vertices form a unique 
        // Tetrahedron in the complete triangulation.
        // Vertex v1 is opposite Triangle v234
        // Vertex v2 is opposite Triangle v134
        // Vertex v3 is opposite Triangle v124
        // Vertex v4 is opposite Triangle v123
        Tetrahedron tetra = new Tetrahedron(v123,v124,v134,v234,v4,v3,v2,v1, sphere);

        // Link each of these Triangles to this new 
        // Tetrahedron
        v123.linkTo(tetra);
        v124.linkTo(tetra);                           
        v134.linkTo(tetra);                           
        v234.linkTo(tetra);

        // Link this new Tetrahedron to each Vertex
        v1.linkTo(tetra);
        v2.linkTo(tetra);
        v3.linkTo(tetra);
        v4.linkTo(tetra);

        // Link each triangle to each vertex that forms
        // part of it.
        v1.linkTo(v123); v1.linkTo(v124); v1.linkTo(v134);
        v2.linkTo(v123); v2.linkTo(v124); v2.linkTo(v234);
        v3.linkTo(v123); v3.linkTo(v134); v3.linkTo(v234);
        v4.linkTo(v124); v4.linkTo(v134); v4.linkTo(v234);

        // Add the new Tetrahedron to the list of this class.
        tetras.add(tetra);
    }
    
    /**
     * Removes a tetrahedron from the triangulation.
     * @param tet
     */
    public void deleteTetrahedron(Tetrahedron tet)
    {
        
        // remove tet from list
        tetras.remove(tet);
        
        // For each vertex that is connected to tet, remove tet from its
        // list of connected tetrahedra
        for (Vertex3d v: tet.verts)
        {
            v.tetras.remove(tet);
        }
        // Any external triangles that are part of tet should be removed
        // from the list entirely, and references to them in the Vertex
        // list should be removed as well.
        // Any internal triangles must have this Tetrahedron removed
        // from their internal list.
        for (Triangle3d tri: tet.tris)
        {
            
            if(tri.isExternal())
            {
                // Remove Triangle from the triangulation
                tris.remove(tri);
                
                // Remove from the lists belonging to 3 Vertices as well
                tri.v0.tris.remove(tri);
                tri.v1.tris.remove(tri);
                tri.v2.tris.remove(tri);
            }
            else
            {
                // Each internal triangle forms part of another connected
                // tetrahedron. We must purge the present tetrahedron
                // from their internal lists of connected tetrahedra. We 
                // can then reset internal/external status by a call to
                // setExternalTriangles()

                if(!tri.tetras.contains(tet))
                {
//                    throw new RuntimeException("Triangle/Tetrahedron link not reciprocal!");
                    System.err.println("Triangle/Tetrahedron link not reciprocal!");
                }
                
                // Finally, remove the Tetrahedron from this Triangle
                tri.tetras.remove(tet);

            }

        }
        
        // Update external vertices/tetrahedrons/triangles
        
        // Could maybe do this only at end of triangulation processing stages,
        // i.e. user decides when to call it, and not automatically after
        // each addition or removal of a tetrahedron.
        setExternalHull();
    }
    
    /**
     * This method purges tetrahedra from the list in order to satisfy the
     * constraint that all vertices lie on the external hull. This is actually
     * very easy to do: it turns out that any tetrahedron for which the
     * circumscribing sphere centre lies outside of the mesh (i.e. is not
     * contained within any tetrahedron) should be deleted. In some sparsely
     * sampled areas this can cause vertices to become separated from the mesh
     * altogether, so we additionally test that any tetrahedron marked for
     * deletion is not the only one connecting any of it's vertices to the
     * mesh.
     */
    protected void exposeConcaveHull() {
    	
        // Flag indicates whether changes were made to the triangulation on the
        // previous iteration over the tetrahedron list
        boolean deleted_simplex = true;
        
        // Repeat loop until we don't find any more tetrahedra to delete
        while(deleted_simplex) {
            // List of Tetrahedrons to delete
            List<Tetrahedron> deleteus = new LinkedList<Tetrahedron>();
            
            // Loop over list of all tetrahedra
            next_tetrahedron:
            for(Tetrahedron t1 : tetras) {
                // If the centre of this tetrahedron's circumsphere lies inside
                // a single tetrahedron in the triangulation, then it is NOT
                // to be deleted
                for(Tetrahedron t2 : tetras) {
                    if(t2.contains(t1.sphere.centre)) {
                        // Tetrahedron forms part of concave hull
                        continue next_tetrahedron;
                    }
                }
                
                // Candidate for deletion
                deleteus.add(t1);
            }
            
            // Now do the removal
            deleted_simplex = false;
            
            for(Tetrahedron deleteme : deleteus) {
                
                // This tetrahedron has been marked for deletion, because
                // the centre of it's circumscribing sphere does not lie
                // inside any tetrahedron of the mesh including itself.

                // However, we first check each of it's vertices to ensure
                // that it is not the only tetrahedron connecting them to
                // the mesh.
                if( deleteme.verts[0].tetras.size() != 1 &&
                    deleteme.verts[1].tetras.size() != 1 &&
                    deleteme.verts[2].tetras.size() != 1 &&
                    deleteme.verts[3].tetras.size() != 1) {
                    // Set flag: if any tetrahedra are deleted then we need to loop
                    // again, as the removal may mean we can delete additional
                    // tetrahedra.
                    deleted_simplex = true;
                    deleteTetrahedron(deleteme);
                }
            }
            
        }
    }
    
    /**
     * Method iterates over the triangles, clusters and tetrahedra and
     * analyses their connectivity in order to determine the subset of each
     * object that lies in the external hull of the triangulation (which can
     * be convex or concave, or indeed any shape).
     * 
     * First, any triangles that form part of two tetrahedra are identified.
     * These necessarily lie on the inside of the surface mesh. Then, vertices
     * and tetrahedra that are part of the external hull can be identified by
     * examining the triangles that they are connected to or composed of.
     * 
     */
    protected void setExternalHull()
    {
        
        // Count number of external triangles.
        N_TRI_EXTERNAL=0;
        
        // Loop over each Triangle and look for one or two attached Tetrahedrons
        for(Triangle3d t : tris)
        {
            if(t.isExternal()) N_TRI_EXTERNAL++;
        }
        
        // Count external vertices
        N_VERT_EXTERNAL = 0;
        
        // Loop over each Vertex.
        for(Vertex3d c : vertices)
        {
            // Initialise to an internal vertex, then try to prove this wrong
            c.EXTERNAL = false;
            
            // Loop over each Triangle that this Vertex is connected to. Look
            // for a single external one.
            for(Triangle3d t: c.tris)
            {
                
                // Is this an external triangle?
                if(t.isExternal())
                {
                    // Yep - therefore Vertex is also external.
                    c.EXTERNAL = true;
                    N_VERT_EXTERNAL++;
                    // Don't look any further.
                    break;
                }
            }
            
        }
        
        // Count external tetrahedra
        N_TETRA_EXTERNAL = 0;
        
        // Loop over each Tetrahedron.
        for(Tetrahedron tetra : tetras)
        {
                        
            boolean tetraIsExternal = false;
            
            // Loop over each Triangle that forms this Tetrahedron.
            // Look for a single external one.
            for(Triangle3d t: tetra.tris)
            {
                // Is this an external triangle?
                if(t.isExternal())
                {
                    // Yep - therefore Vertex is also external.
                    tetraIsExternal = true;
                    N_TETRA_EXTERNAL++;
                    // Don't look any further.
                    break;
                }
            }
            
            tetra.EXTERNAL = tetraIsExternal;
        }
    }
    
    /**
     * Check if a triangle matching the given one already exists in the
     * internal Triangle list. If so, return a reference to it. Otherwise, add
     * the new triangle to the list and return a reference to it.
     */
    protected Triangle3d getUniqueTriangle(Triangle3d t1)
    {
        
        for (Triangle3d t2 : tris) {
            if(t1.equals(t2)) {
                // Return reference to existing triangle.
                return t2;
            }
        }
        
        // Add new Triangle to List.
        tris.add(t1);
        
        // Return reference to new Triangle.
        return t1;
    }
    
    /** Print some stats about this Delaunay triangulation. */
    public final String print()
    {
        return "Triangulation statistics:"+
               "\nNumber of vertices = "+vertices.size()+
               "\n -of which external = "+N_VERT_EXTERNAL+
               "\nNumber of tetrahedra = "+tetras.size()+
               "\n -of which external = "+N_TETRA_EXTERNAL+
               "\nNumber of triangles = "+tris.size()+
               "\n -of which external = "+N_TRI_EXTERNAL+"\n";
    }
    
}