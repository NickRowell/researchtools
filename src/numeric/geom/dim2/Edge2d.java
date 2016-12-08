package numeric.geom.dim2;

import numeric.geom.dim2.EdgeIntersectionTest.IntersectionStatus;

/**
 * Class represents a basic component of the Delaunay Triangulation in 2D: the 1D edges
 * that connect two {@link Vertex2d} and form the edges of {@link Triangle2d}.
 *
 * TODO: tidy up intersection test, write tests for it.
 *
 *
 * @author nrowell
 * @version $Id$
 */
public class Edge2d {
	
	/**
	 * The first vertex.
	 */
	public Vector2d v0;
	
	/**
	 * The second vertex.
	 */
	public Vector2d v1;
	
	/**
	 * Main constructor for the {@link Edge2d}.
	 * @param v0
	 * 	The first {@link Vertex2d}
	 * @param v1
	 * 	The second {@link Vertex2d}
	 */
	public Edge2d(Vector2d v0, Vector2d v1) {
		this.v0 = v0;
		this.v1 = v1;
	}
	
	/**
	 * Tests if two {@link Edge2d}s cross each other.
	 * @param that
	 * 	The second {@link Edge2d}
	 * @return
	 * 	True, if the given {@link Edge2d} crosses this one.
	 */
	public boolean intersects(Edge2d that) {
		
		EdgeIntersectionTest t = new EdgeIntersectionTest(this, that);
		
		if(t.status == IntersectionStatus.INTERSECTING || t.status == IntersectionStatus.COLLINEAR_INTERSECTING) {
			return true;
		}
		else {
			return false;
		}
	}
	
	@Override
	public String toString() {
		return "Vertex 1 = "+v0.toString()+"; Vertex 2 = "+v1.toString();
	}
	
    @Override
    public boolean equals(Object o) {
    	
        // Self check
        if (this == o) {
            return true;
        }
        // Null check
        if (o == null) {
            return false;
        }
        // type check and cast
        if (getClass() != o.getClass()) {
            return false;
        }
        
        Edge2d edge = (Edge2d) o;
        
        // Check if the two vertices are equal
        if(v0.equals(edge.v0) && v1.equals(edge.v1)) {
        	return true;
        }
        else if (v1.equals(edge.v0) && v0.equals(edge.v1)) {
        	return true;
        }
        
        // Fields don't match
        return false;
    }
    
    @Override
    public int hashCode() {
        return v0.hashCode()+v1.hashCode();
    }
    
}