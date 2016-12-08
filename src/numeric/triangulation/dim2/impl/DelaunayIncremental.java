package numeric.triangulation.dim2.impl;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import numeric.geom.dim2.Circle;
import numeric.geom.dim2.Edge2d;
import numeric.geom.dim2.Triangle2d;
import numeric.geom.dim2.Vector2d;
import numeric.triangulation.dim2.base.DelaunayTriangulation2D;

/**
 * 
 * TODO: update comment
 * 
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
public class DelaunayIncremental extends DelaunayTriangulation2D {
	
	/**
     * Main constructor.
     * 
     * @param pvertices
     * 	The {@link List} of {@link Vertex2d} to triangulate.
     */
    public DelaunayIncremental(List<Vector2d> vertices) {    
        super(vertices);
    }
    
    /**
     * Incremental delaunay triangulation.
     */
    @Override
    public void doTriangulation() {
    	
    	// In order to apply the incremental insertion algorithm, we need to start with
    	// a single large triangle that encompasses all of the points. This is obtained by
    	// adding three artificial construction points that form a large enclosing triangle
    	// which itself provides the scaffolding that we need.
    	
        // Get a Circle that encloses all of the points.
        Circle enclosing_circle = Circle.getEnclosingCircle(vertices);
        
        // Now enlarge the Circle - we want one big enough that no circumcircle
        // in the final triangulation contains any of the three artifical points
        // we're going to add that lie at the corners of the enclosing Triangle
        Circle big_circle = new Circle(enclosing_circle, 100);
        
        // Get the corners of the regular Triangle that tightly fits this sphere
        Vector2d[] ABCD = big_circle.getEnclosingTriangleVertices();
        
        // Add three artificial points to triangulation
        Vector2d A = ABCD[0];
        Vector2d B = ABCD[1];
        Vector2d C = ABCD[2];
        
        // Make a secondary list containing artificial points and original points.
        // We triangulate on this list.
        List<Vector2d> temp_vertices = new LinkedList<Vector2d>();
        
        temp_vertices.add(A);
        temp_vertices.add(B);
        temp_vertices.add(C);
        
        // Add original points
        for(Vector2d vertex : vertices) {
        	temp_vertices.add(vertex);
        }
    	vertToTriMap.put(A, new HashSet<Triangle2d>());
    	vertToTriMap.put(B, new HashSet<Triangle2d>());
    	vertToTriMap.put(C, new HashSet<Triangle2d>());
        
        // Add Triangle for these points
        addTriangle(A, B, C);
        
        // Loop over all other vertices and add them one after another.
        
        // Because our initial Triangle contains all the remaining vertices,
        // we guarantee that every new point will result in some part of the
        // mesh requiring retriangulation.
        for(int v=3; v<temp_vertices.size(); v++) {
        	
            // Get list of Triangles that contain this new vertex in their
            // circumcircle. This defines the region of the mesh
            // that requires retriangulation.
            List<Triangle2d> retriangulation = new LinkedList<Triangle2d>();
            for(Triangle2d tri : triToEdgeMap.keySet()) {
                if(tri.circle.contains(temp_vertices.get(v))) {
                    retriangulation.add(tri);
                }
            }
            
            retriangulateMesh(temp_vertices.get(v), retriangulation);
        }
        
        // Remove scaffolding: purge Triangles connecting to artifical points
        Set<Triangle2d> toDelete = new HashSet<>();
        toDelete.addAll(vertToTriMap.get(A));
        toDelete.addAll(vertToTriMap.get(B));
        toDelete.addAll(vertToTriMap.get(C));
        for(Triangle2d tri : toDelete) {
        	deleteTriangle(tri);
        }
        
        // Remove scaffolding: purge artificial points
        vertToTriMap.remove(A);
        vertToTriMap.remove(B);
        vertToTriMap.remove(C);
        
        setExternalHull();
    }
    
    /**
     * TODO: update comment.
     * 
     * Incremental delaunay retriangulation function.
     * 
     * In this method, we remove parts of the mesh to make a convex cavity with
     * the new vertex inside it, then add a bunch of new tetrahedra by
     * connecting the new vertex with the exposed triangular faces around the
     * cavity walls.
     * 
     * @param v1
     * 	The new {@link Vector2d} to add to the triangulation.
     * @param removeMe
     * 	The set of {@link Triangle2d}s that need to be removed from the triangulation
     * to accommodate the new {@link Vector2d}.
     */
    private void retriangulateMesh(Vector2d v1, List<Triangle2d> removeMe) {
    	
        // List of edges bounding the exposed cavity wall
        List<Edge2d> cavity = new LinkedList<Edge2d>();
        
        for(Triangle2d tri : removeMe) {
        	for (Edge2d edge: triToEdgeMap.get(tri)) {
                if(cavity.contains(edge)) {
                	// Remove cavity internal wall
                	cavity.remove(edge);
                }
                else {
                	// New section of cavity exposed
                	cavity.add(edge);
                }
            }
        	deleteTriangle(tri);
        }
        
        // Now fill the cavity with new Triangles by connecting the new vertex
        // to each exposed edge of the cavity wall.
        for(Edge2d cavity_wall : cavity) {
            addTriangle(cavity_wall.v0, cavity_wall.v1, v1);
        }
        
    }
    
}