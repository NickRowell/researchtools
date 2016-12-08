package numeric.geom.dim2;

import java.util.LinkedList;
import java.util.List;

import numeric.geom.dim2.EdgeIntersectionTest.IntersectionStatus;

/**
 * Class represents a general 2D polygon composed of three or more points.
 *
 * TODO: compute polygon area properly
 * TODO: implement sanity check for non-intersecting edges
 *
 * @author nrowell
 * @version $Id$
 */
public class Polygon2d {
	
	/**
	 * The polygon vertices, stored in clockwise winding order in an ordered list.
	 */
	List<Vector2d> vertices;
	
	// Derived quantities used to optimise some geometric operations
	
	/**
	 * The polygon edges, stored in clockwise or anticlockwise winding order
	 */
	List<Edge2d> edges;
	
	/**
	 * Are of the {@link Polygon2d}.
	 */
	public final double area;
	
	/**
	 * Edges of bounding box in the X dimension.
	 */
	public final double[] boundx;
	
	/**
	 * Edges of bounding box in the Y dimension.
	 */
	public final double[] boundy;
	
	/**
	 * Constructs a {@link Polygon2d} from an ordered list of {@link Vector2d}.
	 * @param vertices
	 * 	List of the vertices, in clockwise winding order.
	 */
	public Polygon2d(List<Vector2d> vertices) {
		
		// Sanity check: verify that we have enough points
		if(vertices.size() < 3) {
			throw new IllegalArgumentException("Need at least 3 points to define a Polygon2d! Found "+vertices.size());
		}
		
		// Sanity check: verify that no two vertices are the same
		for(int a=0; a<vertices.size()-1; a++) {
			for(int b=a+1; b<vertices.size(); b++) {
				Vector2d v1 = vertices.get(a);
				Vector2d v2 = vertices.get(b);
				if(v1.equals(v2)) {
					throw new IllegalArgumentException("Vertices "+a+" ("+v1+") and "+b+" ("+v2+") are equal!");
				}
			}
		}
		
		// Build edge list
		edges = new LinkedList<>();
		for(int a=0; a<vertices.size(); a++) {
			
			Vector2d v1 = vertices.get(a);
			Vector2d v2 = vertices.get(a == vertices.size()-1 ? 0 : a+1);
			
			edges.add(new Edge2d(v1, v2));
		}
		
		// Sanity check: verify that no two edges cross each other
//		for(int a=0; a<edges.size()-1; a++) {
//			for(int b=a+1; b<vertices.size(); b++) {
//				Edge2d e1 = edges.get(a);
//				Edge2d e2 = edges.get(b);
//				
//				if(e1.intersects(e2)) {
//					// TODO: ignore intersections if they occur at the line end points.
//					System.out.println("Edge1 = "+e1.toString());
//					System.out.println("Edge2 = "+e2.toString());
//					throw new IllegalArgumentException("Two edges intersect each other!");
//				}
//			
//			}
//		}
		
		this.vertices = vertices;
		
		// Determine bounding box edges
		double minx =  Double.MAX_VALUE;
		double maxx = -Double.MAX_VALUE;
		double miny =  Double.MAX_VALUE;
		double maxy = -Double.MAX_VALUE;
		
		for(Vector2d vertex : vertices) {
			minx = Math.min(vertex.getX(), minx);
			maxx = Math.max(vertex.getX(), maxx);
			miny = Math.min(vertex.getY(), miny);
			maxy = Math.max(vertex.getY(), maxy);
		}
		
		boundx = new double[]{minx, maxx};
		boundy = new double[]{miny, maxy};
		
		// TODO: compute this properly
		area = 1.0;
	}
	
	/**
	 * Determines if the point lies inside the polygon.
	 * @param point
	 * 	The {@link Vector2d} to test.
	 * @return
	 * 	True, if the point lies inside the polygon.
	 */
	public boolean contains(Vector2d point) {
		
		// Check for easy rejection based on bounding box
		if(point.getX() < boundx[0] || point.getX() > boundx[1]) {
			// Point lies outside bounding box in (at least) X dimension
			return false;
		}
		if(point.getY() < boundy[0] || point.getY() > boundy[1]) {
			// Point lies outside bounding box in Y dimension
			return false;
		}
		
		// Use RAY CASTING: draw a line from some point outside the polygon to the test point
		// and count the number of times it cuts an edge. If this is odd, then the point lies
		// inside the polygon; if it's an even number, the point lies outside the polygon.
		
		// Get a point outside the polygon (easy s we have a bounding box):
		Vector2d v = new Vector2d(boundx[0] - 1.0, boundy[0] - 1.0);
		
		// Create an Edge2d from this point to the test point
		Edge2d ray = new Edge2d(v, point);
		
		// Count how many of the polygon edges intersect with this edge
		int n = 0;
		for(Edge2d edge : edges) {
			
			EdgeIntersectionTest test = new EdgeIntersectionTest(edge, ray);
			
			if(test.status == IntersectionStatus.INTERSECTING) {
				// Ignore COLLINEAR_INTERSECTING
				n++;
			}
		}
		
		if(n%2==0) {
			// n is even: point is outside the polygon
			return false;
		}
		else {
			// n is odd: point is inside the polygon
			return true;
		}
		
	}
	
}