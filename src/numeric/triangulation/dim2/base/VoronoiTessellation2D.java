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
import numeric.geom.dim2.Polygon2d;
import numeric.geom.dim2.Triangle2d;
import numeric.geom.dim2.Vector2d;
import numeric.triangulation.dim2.util.TriangulationUtils;

/**
 * Class represents the 2D Voronoi Tessellation, the dual graph of the 2D Delaunay Triangulation.
 * 
 * TODO: properly implement tessellation colouring method
 * TODO: fix polygons at edge of tessellation
 * TODO: are all elements of the data structure actually required?
 * 
 * @author nrowell
 * @version $Id$
 */
public class VoronoiTessellation2D {
	
    /**
     * Maps {@link Vector2d} to the {@link Edge2d}s that they are connected to.
     */
    public Map<Vector2d, Set<Edge2d>> vertToEdgeMap;
    
	/**
	 * Maps {@link Polygon2d}s to the {@link Edge2d}s that they are connected to.
	 */
    public Map<Polygon2d, Set<Edge2d>> polyToEdgeMap;
    
    /**
	 * Maps {@link Edge2d}s to the {@link Polygon2d}s that they are connected to.
	 */
    public Map<Edge2d, Set<Polygon2d>> edgeToPolyMap;
    
	/**
	 * Constructs the 2D Voronoi Tessellation for a set of points from the 2D
	 * Delaunay Triangulation of the points.
	 *  
	 * @param delTri
	 * 	The {@link DelaunayTriangulation2D} of the points.
	 */
    public VoronoiTessellation2D(DelaunayTriangulation2D delTri) {
    	
    	vertToEdgeMap = new HashMap<>();
    	polyToEdgeMap = new HashMap<>();
    	edgeToPolyMap = new HashMap<>();
    	
    	// How is the Voronoi Tessellation constructed?
    	//  - vertices lie at the centres of the circumcircles of each triangle
    	//  - each vertex is connected to the vertices associated with it's triangle neighbours
    	//  - each vertex therefore is connected to precisely three other vertices
    	//  - the polygons around the edge of the triangulation are not closed: they stretch off to infinity
    	
    	// Each vertex in the original triangulation corresponds to one voronoi cell
    	for(Entry<Vector2d, Set<Triangle2d>> entry : delTri.vertToTriMap.entrySet()) {
    		
    		// Central vertex of the Voronoi cell
    		Vector2d vertex = entry.getKey();
    		Set<Triangle2d> tris = entry.getValue();
    		
    		// Polygon edges associated with this Voronoi cell
    		Set<Edge2d> edges = new HashSet<>();
    		
    		// Vertices that form the corners of this voronoi cell, in winding order
    		List<Vector2d> verts = new LinkedList<>();
    		
    		// The corners of the Voronoi cell are formed from the centres of the circumcircles
    		// for each attached triangle. We need to first put the triangles into an order list
    		// so we know which pairs of centres to connect.
    		List<Triangle2d> linkedTris = new LinkedList<>();
    		try {
    			linkedTris.addAll(getOrderedList(tris, vertex));
    		}
    		catch(RuntimeException re) {
    			// Occurs around the triangulation edge because the fans aren't closed.
    			// TODO: need to handle the tessellation edges.
    			continue;
    		}

    		for(int i=0; i<linkedTris.size(); i++) {
    			Triangle2d t0 = linkedTris.get(i);
    			Triangle2d t1 = linkedTris.get(i==(linkedTris.size()-1) ? 0 : i+1);
    			
    			Vector2d v0 = t0.circle.centre;
    			Vector2d v1 = t1.circle.centre;
    			Edge2d edge = new Edge2d(v0, v1);
    			edges.add(edge);
    			verts.add(v0);
    			
    			// Insert into the data structure
    			if(!vertToEdgeMap.containsKey(v0)) {
    				vertToEdgeMap.put(v0, new HashSet<Edge2d>());
    			}
    			vertToEdgeMap.get(v0).add(edge);
    			if(!vertToEdgeMap.containsKey(v1)) {
    				vertToEdgeMap.put(v1, new HashSet<Edge2d>());
    			}
    			vertToEdgeMap.get(v1).add(edge);
    		}
    		
    		Polygon2d poly = new Polygon2d(verts);
    		
    		// Insert into the data structure
    		polyToEdgeMap.put(poly, edges);
    		for(Edge2d edge : edges) {
    			if(!edgeToPolyMap.containsKey(edge)) {
    				edgeToPolyMap.put(edge, new HashSet<Polygon2d>());
    			}
    			edgeToPolyMap.get(edge).add(poly);
    		}
    	}
    	
    }
    
    /**
     * Map each {@link Polygon2d} to a set of it's neighbours.
     * @return
     */
    public Map<Polygon2d, Set<Polygon2d>> getNeighboursMap() {
    	
    	Map<Polygon2d, Set<Polygon2d>> neighboursMap = new HashMap<>();
    	
    	for(Entry<Polygon2d, Set<Edge2d>> entry : polyToEdgeMap.entrySet()) {
    		
    		// Set of neighbouring polygons
    		Set<Polygon2d> neighbours = new HashSet<>();
    		
    		// Loop over each edge that composes this polygon
    		for(Edge2d edge : entry.getValue()) {
    			// Get the two polygons that share this edge (one of which will be the central one
    			// which we shouldn't regard as a neighbour of itself).
    			for(Polygon2d poly : edgeToPolyMap.get(edge)) {
    				neighbours.add(poly);
    			}
    		}
    		
    		Polygon2d poly = entry.getKey();
    		
    		// Remove the central polygon from it's neighbours map
    		neighbours.remove(poly);
    		
    		// Add to the neighbours map
    		neighboursMap.put(poly, neighbours);
    	}
    	
    	return neighboursMap;
    }
    
    /**
     * Sorts a set of triangles that are each connected continuously to the same vertex
     * (in a fan pattern) into an ordered list such that neighbouring triangles are
     * consecutive.
     * @param tris
     * 	The set of triangles.
     * @param centralVertex	
     * 	All the triangles contain this vertex, which is at the centre of a star/fan.
     * @return
     * 	A sorted list of the triangles.
     */
    public static List<Triangle2d> getOrderedList(Set<Triangle2d> tris, Vector2d centralVertex) {
    	
    	// Initialise List to store the sorted triangles
    	List<Triangle2d> sortedTris = new LinkedList<>();
    	
    	// Load all the triangles into a new list so we can remove them as we connect the list
    	List<Triangle2d> unsortedTris = new LinkedList<>();
    	unsortedTris.addAll(tris);
    	
    	// Make a set of vertices for each triangle. This is useful for searching later.
    	Map<Triangle2d, Set<Vector2d>> triToVertsMap = new HashMap<>();
    	for(Triangle2d tri : tris) {
    		Set<Vector2d> vertSet = new HashSet<>();
    		vertSet.add(tri.v0);
    		vertSet.add(tri.v1);
    		vertSet.add(tri.v2);
    		triToVertsMap.put(tri, vertSet);
    		
    		// Sanity check: verify that every triangle contains the central vertex
    		if(!vertSet.contains(centralVertex)) {
    			throw new RuntimeException("Triangle "+tri+" does not contain the central vertex "+centralVertex);
    		}
    	}
    	
    	// Pick one triangle to start the list and remove it from the set
    	Triangle2d currentPiece = unsortedTris.remove(0);
    	sortedTris.add(currentPiece);
    	
    	// Pick one vertex of the currentPiece (which is not the centralVertex) and
    	// search for the triangle that also contains this vertex. This gives the next
    	// triangle in the loop.
    	List<Vector2d> verts = new LinkedList<Vector2d>(triToVertsMap.get(currentPiece));
    	verts.remove(centralVertex);
    	Vector2d searchV = verts.remove(0);
    	Vector2d lastV = verts.remove(0);
    	
    	// Loop over the rest of the triangles and look for one connected to this
    	while(!unsortedTris.isEmpty()) {
    		
    		Triangle2d nextPiece = null;
    		
    		for(Triangle2d t : unsortedTris) {
    			// Does this triangle contain the target vertex?
    			if(triToVertsMap.get(t).contains(searchV)) {
    				// Found the next triangle in the list
    				nextPiece = t;
    				break;
    			}
    		}

			// Add it to the sorted list...
			sortedTris.add(nextPiece);
			
			// ...remove it form the unsorted list...
			unsortedTris.remove(nextPiece);
			
			// ...identify the next vertex to continue the search
			for(Vector2d nextVert : triToVertsMap.get(nextPiece)) {
				if(!nextVert.equals(searchV) && !nextVert.equals(centralVertex)) {
					searchV = nextVert;
					break;
				}
			}
    	}
    	
    	// Sanity check: verify that the fan is closed, by checking that the vertex
    	// we're now searching for is the remaining one from the first triangle
    	if(!searchV.equals(lastV)) {
    		throw new RuntimeException("Triangles don't form a continuous loop!");
    	}
    	
    	return sortedTris;
    }
    
    /**
     * Generate a coloured image of the {@link VoronoiTessellation2D}.
     * 
     * @param width
     * 	Width of the image to generate
     * @param height
     * 	Height of the image to generate
     * @param colours
     * 	Array of colours to use
     * @return
     * 	A {@link BufferedImage} showing the tessellation.
     */
	public BufferedImage getColouredImage(int width, int height, int[] colours) {
		
		// Colour the triangulation
		Map<Polygon2d, Integer> colourMap = TriangulationUtils.colourTessellation(getNeighboursMap(), colours);
		
		// Render an image of the result
		BufferedImage im = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		
		// Fill in each polygon
		for(int i=0; i<width; i++) {
			for(int j=0; j<height; j++) {
				Vector2d point = new Vector2d(i, j);
				for(Polygon2d polygon : polyToEdgeMap.keySet()) {
					if(polygon.contains(point)) {
						int colour = colourMap.get(polygon);
						im.setRGB(i, j, colour);
					}
				}
			}
		}
		
		return im;
	}
	
	/**
     * Generate a wireframe image of the {@link VoronoiTessellation2D}.
     * 
     * @param width
     * 	Width of the image to generate
     * @param height
     * 	Height of the image to generate
     * @return
     * 	A {@link BufferedImage} showing the tessellation.
     */
	public BufferedImage getWireframeImage(int width, int height) {
		
		// Render an image of the result
		BufferedImage im = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		
		Graphics2D graphics = im.createGraphics();
		graphics.setPaint(new Color(255, 255, 255));
		graphics.fillRect(0, 0, width, height);

		// Draw border
		for(int i=0; i<width; i++) {
			im.setRGB(i, 0, 0x000000);
			im.setRGB(i, height-1, 0x000000);
		}
		for(int j=0; j<height; j++) {
			im.setRGB(0, j, 0x000000);
			im.setRGB(width-1, j, 0x000000);
		}
		
		// Draw the edges
		for(Edge2d edge : edgeToPolyMap.keySet()) {
			Vector2d v0 = edge.v0;
			Vector2d v1 = edge.v1;
			int colour = 0x000000;
			Rendering.drawLine(im, new int[]{(int)v0.getX(), (int)v0.getY()}, new int[]{(int)v1.getX(), (int)v1.getY()}, colour);
		}
		
		return im;
	}

}