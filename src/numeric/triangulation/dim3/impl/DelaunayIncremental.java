package numeric.triangulation.dim3.impl;

import java.util.LinkedList;
import java.util.List;

import numeric.geom.dim3.Sphere;
import numeric.geom.dim3.Vector3d;
import numeric.triangulation.dim3.base.DelaunayTriangulation3D;
import numeric.triangulation.dim3.infra.Tetrahedron;
import numeric.triangulation.dim3.infra.Triangle3d;
import numeric.triangulation.dim3.infra.Vertex3d;

/**
 * Efficient Delaunay triangulation algorithm using incremental insertion.
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
public class DelaunayIncremental extends DelaunayTriangulation3D {
	
	/**
     * Main constructor.
     * 
     * @param pvertices
     * 	The {@link List} of {@link Vertex3d} to triangulate.
     */
    public DelaunayIncremental(List<Vertex3d> pvertices) {    
        super(pvertices);
    }
    
    /**
     * Incremental delaunay triangulation.
     * 
     * @param _clusters
     * @param os 
     */
    @Override
    public void doTriangulation() {
    	
        // Get a sphere that encloses all of the points.
        Sphere enclosing_sphere = Sphere.getEnclosingSphere(vertices);
        
        // Now enlarge the sphere - we want one big enough that no circumsphere
        // in the final triangulation contains any of the four artifical points
        // we're going to add that lie at the corners of the enclosing tetrahedron
        Sphere big_sphere = new Sphere(enclosing_sphere, 100);
        
        // Get the corners of the regular tetrahedron that tightly fits this sphere
        Vector3d[] ABCD = big_sphere.getEnclosingTetrahedronVertices();
        
        // Add four artificial points to start of cluster list
        Vertex3d A = new Vertex3d(ABCD[0]);
        Vertex3d B = new Vertex3d(ABCD[1]);
        Vertex3d C = new Vertex3d(ABCD[2]);
        Vertex3d D = new Vertex3d(ABCD[3]);
        
        Sphere sphere = new Sphere(A, B, C, D);
        
        // Make a secondary list containing artificial points and original points.
        // We triangulate on this list.
        List<Vertex3d> temp_vertices = new LinkedList<Vertex3d>();
        
        temp_vertices.add(0, A);
        temp_vertices.add(0, B);
        temp_vertices.add(0, C);
        temp_vertices.add(0, D);
        
        // Add original points
        for(Vertex3d vertex : vertices) temp_vertices.add(vertex);
        
        // Add tetrahedron for these points
        addTetrahedron(A, B, C, D, sphere);
        
        // Loop over all other vertices and add them one after another.
        
        // Because our initial tetrahedron contains all the remaining vertices,
        // we guarantee that every new point will result in some part of the
        // mesh requiring retriangulation.
        for(int v=4; v<temp_vertices.size(); v++)
        {
            // Also check stop flag?
            
            // Get list of tetrahedra that contain this new vertex in their
            // circumsphere. This defines the region of the mesh
            // that requires retriangulation.
            List<Tetrahedron> retriangulation = new LinkedList<Tetrahedron>();
            for(Tetrahedron tetra : tetras)
            {
                if(tetra.circumSphereContains(temp_vertices.get(v)))
                    retriangulation.add(tetra);
            }
            
            retriangulateMesh(temp_vertices.get(v), retriangulation);
            
        }
        
        // Now remove tetrahedra that connect to artifical points. Note that 
        // the artificial points are automatically excluded from the final
        // triangulation.
        int i = A.tetras.size();
        for(int x=0; x<i; x++) deleteTetrahedron(A.tetras.get(0));
        
        i = B.tetras.size();
        for(int x=0; x<i; x++) deleteTetrahedron(B.tetras.get(0));
        
        i = C.tetras.size();
        for(int x=0; x<i; x++) deleteTetrahedron(C.tetras.get(0));
        
        i = D.tetras.size();
        for(int x=0; x<i; x++) deleteTetrahedron(D.tetras.get(0));
        
        // At this point, we attempt to extract the concave hull by deleting
        // tetrahedra for which the centre of the circumscribing sphere does
        // not lie inside any tetrahedron in the triangulation.
        exposeConcaveHull();
        
        setExternalHull();
        
    }
    
    /**
     * Incremental delaunay retriangulation function.
     * 
     * In this method, we remove parts of the mesh to make a convex cavity with
     * the new vertex inside it, then add a bunch of new tetrahedra by
     * connecting the new vertex with the exposed triangular faces around the
     * cavity walls.
     * 
     * @param v1
     * @param removeMe
     */
    private void retriangulateMesh(Vertex3d v1, List<Tetrahedron> removeMe)
    {
        // First, remove all tetrahedra from existing mesh. We need to also
        // maintain a list of the triangles lying in the cavity wall, so this
        // deleteTetrahedron method is overridden for this class.
        
        // List of triangles lying in exposed cavity wall.
        List<Triangle3d> cavity = new LinkedList<Triangle3d>();
        
        for(Tetrahedron tetra : removeMe)
        {
            // Delete tetrahedron and update the list of the triangles that are
            // exposed by it's removal.
            deleteTetrahedronAndBuildCavityWall(tetra, cavity);
        }
        
        // Now fill the cavity with new Tetrahedra by connecting the new Vertex
        // the each exposed face of the cavity wall.
        for(Triangle3d cavity_wall : cavity)
        {
            
            // Circumsphere for new tetrahedron
            Sphere sphere = new Sphere(cavity_wall.v0, cavity_wall.v1, cavity_wall.v2, v1);
            
            addTetrahedron(cavity_wall.v0, cavity_wall.v1, cavity_wall.v2, v1, sphere);
        }
        
    }
    
    /**
     * Removes a tetrahedron from the triangulation
     * @param tet
     */
    private void deleteTetrahedronAndBuildCavityWall(Tetrahedron tet, List<Triangle3d> cavity)
    {
        
        // remove tet from list
        tetras.remove(tet);
        
        // For each vertex that is connected to tet, remove tet from its
        // list of connected tetrahedra
        for (Vertex3d v: tet.verts)
        {
            v.tetras.remove(tet);
        }
        
        // Add the four triangular faces of this tetrahedron to the cavity wall
        // list. If any of the faces are already in the list, then they are
        // removed from the cavity.
        for (Triangle3d tri: tet.tris)
        {
            if(cavity.contains(tri)) cavity.remove(tri);
            else                     cavity.add(tri);
        }
        
        for (Triangle3d tri: tet.tris)
        {
            
            // Any external triangles that are part of this tetrahedron should be
            // removed from the main triangulation.
            // We also need to update the vertices that this triangle was attached
            // to.
            if(tri.isExternal())
            {
                tris.remove(tri);
                
                // Remove from the lists belonging to 3 Vertices as well
                tri.v0.tris.remove(tri);
                tri.v1.tris.remove(tri);
                tri.v2.tris.remove(tri);
            }
            // Any internal triangles must have this Tetrahedron removed
            // from their internal list. They then become part of the exposed
            // cavity wall.
            else
            {
                if(!tri.tetras.contains(tet))
                {
                    System.err.println("Triangle/Tetrahedron link not reciprocal!");
                }
                
                // Remove the Tetrahedron from this Triangle
                tri.tetras.remove(tet);
            }
            
        }
        
    }
    
}