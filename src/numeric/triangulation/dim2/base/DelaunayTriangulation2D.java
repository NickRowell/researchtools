package numeric.triangulation.dim2.base;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import images.Rendering;

import java.util.Set;

import numeric.geom.dim2.Edge2d;
import numeric.geom.dim2.Triangle2d;
import numeric.geom.dim2.Vector2d;
import numeric.triangulation.dim2.util.TriangulationUtils;

/**
 * 
 * TODO: finish implementing the expose concave hull method
 * TODO: document and tidy up
 * TODO: can the triangulation algorithm be generalized to N dimensions?
 * 
 * 
 * 
 * What kind of data structure do we need? Are Edge2Ds necessary? Do we need to be able to extract the concave hull?
 * 
 * Need to know connectivity of triangles, i.e. which other triangles each triangle is connected to.
 * 
 * In order to determine the outer hull of the triangulation we need to use (N-1)D objects, i.e. the components
 * of the simplex, as distinct entities in the data structure.
 * 
 * Edges need to know which triangles they are part of (so we can detect the outer hull)
 * Triangles need to know which edges they are composed of (so can delete redundant edges when a triangle is deleted)
 * 
 * Can we do this using Map structures within this class, rather than having to add fields to the geometric classes (thus
 * introducing a layer of unnecessary classes)?
 * 
 * Do vertices need to know which edges and triangles they are connected to? Probably not.
 * 
 * 
 * To generalise to N dimensions, we need:
 * 
 * Points (the things being triangulated)
 * Simplexes (N dimensional objects)
 * N-1 dimensional objects (components of the simplex), in order to determine the external hull by the number of simplexes they form part of
 * 
 * Required operations:
 *  1) Merge a new triangle (defined by 3 vertices) into the triangulation structure; the edges connecting the vertices
 *     may already be part of the triangulation so need to check and reuse existing edges.
 * 
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
public abstract class DelaunayTriangulation2D {
    
	// The set of points to be triangulated:
	
	/**
     * Handle to list of all Vertices to be triangulated.
     */
    protected List<Vector2d> vertices;
	
	// Main triangulation data structures:

    /**
     * Maps triangles to the edges they are composed of, so we don't need to keep
     * track of this connection within the triangle & edge classes.
     * 
     * Used to delete edges from the triangulation
     */
    public Map<Triangle2d, Set<Edge2d>> triToEdgeMap;
    
    /**
     * Maps Edges to the (one or two) triangles that they are part of.
     * Used to identify the outer hull of the triangulation
     */
    public Map<Edge2d, Set<Triangle2d>> edgeToTriMap;
    
    /**
     * Maps vertices to the Triangles that they are connected to.
     */
    public Map<Vector2d, Set<Triangle2d>> vertToTriMap;
    
    
	// Derived triangulation data structures, updated whenever the triangulation is revised.
    
    /**
     * List containing the {@link Edge2d}s lying on the outer hull of the triangulation.
     */
    public Set<Edge2d> edgesHull;
    
    /**
     * List containing the {@link Vector2d}s lying on the outer hull of the triangulation.
     */
    public Set<Vector2d> vectorHull;
	
    /**
     * Main constructor.
     * 
     * @param vertices
     * 	The {@link List} of {@link Vertex2d} to triangulate.
     */
    public DelaunayTriangulation2D(List<Vector2d> vertices) {
    	this.vertices = vertices;
    	triToEdgeMap = new HashMap<Triangle2d, Set<Edge2d>>();
    	edgeToTriMap = new HashMap<Edge2d, Set<Triangle2d>>();
    	vertToTriMap = new HashMap<Vector2d, Set<Triangle2d>>();
    	edgesHull = new HashSet<>();
    	vectorHull = new HashSet<>();
    	for(Vector2d vertex : vertices) {
    		vertToTriMap.put(vertex, new HashSet<Triangle2d>());
    	}
    }
    
    /**
     * TODO: rewrite comment.
     */
    public abstract void doTriangulation();
    
    /**
     * Add a triangle to the triangulation.
     * @param v0
     * @param v1
     * @param v2
     */
    protected void addTriangle(Vector2d v0, Vector2d v1, Vector2d v2) {
        
    	// Create a new Triangle from the 3 vertices
        Triangle2d v012 = new Triangle2d(v0, v1, v2);
        
    	// Create Edges (N-1 dimension object) for the triangle, or reuse any existing ones
    	Edge2d e01 = getUniqueEdge(new Edge2d(v0, v1));
    	Edge2d e02 = getUniqueEdge(new Edge2d(v0, v2));
    	Edge2d e12 = getUniqueEdge(new Edge2d(v1, v2));
    	
    	// Map the new Triangle to the Edges that it is composed of
    	Set<Edge2d> edges = new HashSet<>();
    	edges.add(e01);
    	edges.add(e02);
    	edges.add(e12);
    	
    	// Insert the triangle into the triangulation data structure
    	triToEdgeMap.put(v012, edges);
    	
    	// Map each vertex to the Triangle
    	vertToTriMap.get(v0).add(v012);
    	vertToTriMap.get(v1).add(v012);
    	vertToTriMap.get(v2).add(v012);
    	
    	// Map each Edge to the Triangle
    	edgeToTriMap.get(e01).add(v012);
    	edgeToTriMap.get(e02).add(v012);
    	edgeToTriMap.get(e12).add(v012);
    }
    
    /**
     * Removes a Triangle from the triangulation.
     * @param tri
     */
    public void deleteTriangle(Triangle2d tri) {
        
    	// Retrieve the list of Edges that this triangle is composed of
    	Set<Edge2d> edges = triToEdgeMap.get(tri);
    	
    	// Disconnect the triangle from each Edge; if any edge is no longer connected
    	// to a triangle then it should be removed as well
    	for(Edge2d edge : edges) {
    		edgeToTriMap.get(edge).remove(tri);
    		if(edgeToTriMap.get(edge).isEmpty()) {
    			edgeToTriMap.remove(edge);
    		}
    	}
    	
    	// Disconnect the Triangle from each vertex
    	vertToTriMap.get(tri.v0).remove(tri);
    	vertToTriMap.get(tri.v1).remove(tri);
    	vertToTriMap.get(tri.v2).remove(tri);
    	
    	// Remove the triangle from the Map
    	triToEdgeMap.remove(tri);
    }
    
    /**
     * TODO: update comment.
     * 
     * NOTE: this method is not fully implemented currently. It can leave
     * isolated points by deleting triangles a bit too zealously; the way to
     * avoid this is to check that a triangle to be deleted is not the sole
     * triangle connected to any of the vertices that compose it. This requires
     * an additional data structure to track map vertices to the triangles
     * that are connected to them.
     * 
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
    public void exposeConcaveHull() {
    	
        // Flag indicates whether changes were made to the triangulation on the
        // previous iteration over the Triangle2d list
        boolean deleted_simplex = true;
        
        // Repeat loop until we don't find any more Triangle2d to delete
        while(deleted_simplex) {
        	
            // List of Triangle2ds to delete
            List<Triangle2d> deleteus = new LinkedList<>();
            
            // Loop over list of all Triangle2d
            next_triangle:
            for(Triangle2d t1 : triToEdgeMap.keySet()) {
            	
                // If the centre of this Triangle2d's circumsphere lies inside
                // a single Triangle2d in the triangulation, then it is NOT
                // to be deleted
                for(Triangle2d t2 : triToEdgeMap.keySet()) {
                	
                    if(t2.contains(t1.circle.centre)) {
                        // Triangle2d forms part of concave hull
                        continue next_triangle;
                    }
                }
                
                // Candidate for deletion
                deleteus.add(t1);
            }
            
            // Now do the removal
            deleted_simplex = false;
            
            for(Triangle2d deleteme : deleteus) {
                
            	// NOTE: should check if any of the vertices that this
            	// triangle is composed off would be left isolated by the
            	// deletion of this triangle (i.e. if this is the only
            	// triangle to which they are connected).
                deleted_simplex = true;
                deleteTriangle(deleteme);
            }
            
        }
    }
    
    /**
     * TODO: update comment.
     * 
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
    public void setExternalHull() {
    	
    	// Clear existing outer hull data structure
    	edgesHull.clear();
    	vectorHull.clear();
    	
    	// Loop over the Edges: any that are linked to only one Triangle form part of the outer hull
    	for(Entry<Edge2d, Set<Triangle2d>> entry : edgeToTriMap.entrySet()) {
    		
    		Set<Triangle2d> tris = entry.getValue();
    		
    		if(tris.size()==1) {
    			
    			Edge2d edge = entry.getKey();
    			
    			edgesHull.add(edge);
    			
    			// Add the Vectors that form this edge; note we don't need to worry about duplicating vertices
    			// interconnecting two edges because the set is unchanged if it already contains an object.
    			vectorHull.add(edge.v0);
    			vectorHull.add(edge.v1);
    		}
    		else if(tris.size()==2){
    			// Edge is interior to the triangulation
    		}
    		else {
    			// Sanity check: this should not happen!
    			throw new RuntimeException("Found an edge connected to "+tris.size()+" Triangles!");
    		}
    	}
    }
    
    /**
     * Check if an Edge matching the given one already exists in the
     * internal Edge list. If so, return a reference to it. Otherwise, add
     * the new Edge to the list and return a reference to it.
     */
    protected Edge2d getUniqueEdge(Edge2d edge) {
    	if(!edgeToTriMap.containsKey(edge)) {
    		edgeToTriMap.put(edge, new HashSet<Triangle2d>());
    	}
        return edge;
    }
    
    /**
     * Maps each {@link Triangle2d} in the triangulation to the {@link Set} of it's neighbouring
     * {@link Triangle2d}s. Each triangle has precisely three neighbours except for those around the
     * edge of the triangulation which have two.
     * 
     * @return
     * 	A Mapping of each {@link Triangle2d} in the triangulation to the {@link Set} of it's neighbouring
     * {@link Triangle2d}s.
     */
    public Map<Triangle2d, Set<Triangle2d>> getNeighboursMap() {

		Map<Triangle2d, Set<Triangle2d>> neighboursMap = new HashMap<>();
		
		for(Triangle2d tri : triToEdgeMap.keySet()) {
			
			// Add the triangle to the map
			Set<Triangle2d> neighbours = new HashSet<Triangle2d>();
			neighboursMap.put(tri, neighbours);
			
			// Look up it's neighbours by examining the edges it shares with other triangles
			Set<Edge2d> edges = triToEdgeMap.get(tri);
			
			// Constraint: there should only ever be 3 edges
			for(Edge2d edge : edges) {
				
				// Retrieve the triangles connected to this edge
				Set<Triangle2d> tris = edgeToTriMap.get(edge);
				
				// Constraint: there should be precisely one or two triangles connected to this edge,
				// one triangle in the case of edges on the external hull of the triangulation.
				
				// Add all the triangles connected to this edge (this includes the current triangle)
				neighbours.addAll(tris);
			}
			// Remove the current triangle (don't consider it one it's own neighbours)
			neighbours.remove(tri);
		}
		
		return neighboursMap;
    }
    
    /**
     * Print some stats about this Delaunay triangulation.
     */
    public final String toString() {
        return "Triangulation statistics:"+
               "\nNumber of vertices  = "+vertices.size()+
               "\n -of which external = "+vectorHull.size()+
               "\nNumber of edges     = "+edgeToTriMap.size()+
               "\n -of which external = "+edgesHull.size()+
               "\nNumber of triangles = "+triToEdgeMap.size()+"\n";
    }
    
    /**
     * Generate a coloured image of the {@link DelaunayTriangulation2D}.
     * 
     * @param width
     * 	Width of the image to generate
     * @param height
     * 	Height of the image to generate
     * @param colours
     * 	Array of colours to use
     * @return
     * 	A {@link BufferedImage} showing the triangulation.
     */
	public BufferedImage getColouredImage(int width, int height, int[] colours) {

		// Colour the triangulation
		Map<Triangle2d, Integer> colourMap = TriangulationUtils.colourTriangulation(getNeighboursMap(), colours);
		
		// Render an image of the result
		BufferedImage im = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		
		// Fill in each triangle
		for(int i=0; i<width; i++) {
			for(int j=0; j<height; j++) {
				Vector2d point = new Vector2d(i, j);
				for(Triangle2d triangle : triToEdgeMap.keySet()) {
					if(triangle.contains(point)) {
						int colour = colourMap.get(triangle);
						im.setRGB(i, j, colour);
					}
				}
			}
		}
		
		return im;
	}

    /**
     * Generate a wireframe image of the {@link DelaunayTriangulation2D}.
     * 
     * @param width
     * 	Width of the image to generate
     * @param height
     * 	Height of the image to generate
     * @return
     * 	A {@link BufferedImage} showing the triangulation.
     */
	public BufferedImage getWireframeImage(int width, int height) {

		// Render an image of the result
		BufferedImage im = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);

		Graphics2D    graphics = im.createGraphics();
		graphics.setPaint ( new Color ( 255, 255, 255 ) );
		graphics.fillRect ( 0, 0, width, height);

		// Draw border
		for(int i=0; i<width; i++) {
			im.setRGB(i, 0, 0x000000);
			im.setRGB(i, height-1, 0x000000);
		}
		for(int j=0; j<height; j++) {
			im.setRGB(0, j, 0x000000);
			im.setRGB(width-1, j, 0x000000);
		}
		
		// Draw the complete triangulation
		for(Triangle2d triangle : triToEdgeMap.keySet()) {
			
			Vector2d v0 = triangle.v0;
			Vector2d v1 = triangle.v1;
			Vector2d v2 = triangle.v2;
			
			int colour = 0x000000;
			
			Rendering.drawLine(im, new int[]{(int)v0.getX(), (int)v0.getY()}, new int[]{(int)v1.getX(), (int)v1.getY()}, colour);
			Rendering.drawLine(im, new int[]{(int)v0.getX(), (int)v0.getY()}, new int[]{(int)v2.getX(), (int)v2.getY()}, colour);
			Rendering.drawLine(im, new int[]{(int)v1.getX(), (int)v1.getY()}, new int[]{(int)v2.getX(), (int)v2.getY()}, colour);
		}
		
		// Draw vertices
//		for(Vector2d point : tri.vertToTriMap.keySet()) {
//			Rendering.drawCircle(new int[]{(int)point.getX(), (int)point.getY()}, im, 5, 0xf28354);
//		}
		
		// Draw external vertices
//		for(Vector2d point : tri.vectorHull) {
//			Rendering.drawCircle(new int[]{(int)point.getX(), (int)point.getY()}, im, 7, 0x8700f2);
//		}
		
		// Now draw the outer hull
//		for(Edge2d edge : tri.edgesHull) {
//			Vector2d v0 = edge.v0;
//			Vector2d v1 = edge.v1;
//			int colour = 0x930058;
//			Rendering.drawLine(im, new int[]{(int)v0.getX(), (int)v0.getY()}, new int[]{(int)v1.getX(), (int)v1.getY()}, colour);
//		}
		
		return im;
	}
    
}