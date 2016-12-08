package numeric.geom.dim2;

/**
 * TODO: implement contains method.
 * 
 * 
 * Generic class representing triangles in 2D.
 *
 * @author nrowell
 * @version $Id$
 * @param <T>
 */
public class Triangle2d {

	/**
	 * The first vertex.
	 */
	public final Vector2d v0;

	/**
	 * The second vertex.
	 */
	public final Vector2d v1;

	/**
	 * The third vertex.
	 */
	public final Vector2d v2;
	
	
	// Derived quantities used to optimise some geometric operations
	
	/**
	 * Are of the {@link Triangle2d}.
	 */
	public final double area;
	
	/**
	 * Circumscribing {@link Circle}.
	 */
	public final Circle circle;
	
	/**
	 * Edges of bounding box in the X dimension.
	 */
	public final double[] boundx;
	
	/**
	 * Edges of bounding box in the Y dimension.
	 */
	public final double[] boundy;
	
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
	public Triangle2d(Vector2d v0, Vector2d v1, Vector2d v2) {
		this.v0 = v0;
		this.v1 = v1;
		this.v2 = v2;
		circle = Circle.fitCircle(v0, v1, v2);
		
		double minx = Math.min(v0.getX(), Math.min(v1.getX(), v2.getX()));
		double maxx = Math.max(v0.getX(), Math.max(v1.getX(), v2.getX()));
		double miny = Math.min(v0.getY(), Math.min(v1.getY(), v2.getY()));
		double maxy = Math.max(v0.getY(), Math.max(v1.getY(), v2.getY()));
		
		boundx = new double[]{minx, maxx};
		boundy = new double[]{miny, maxy};
		
		area = 0.5 * (v0.getX()*(v1.getY() - v2.getY())
				    + v1.getX()*(v2.getY() - v0.getY())
				    + v2.getX()*(v0.getY() - v1.getY()));
	}
	
	/**
	 * Determines if the {@link Triangle2d} contains the point.
	 * @param point
	 * 	The point to test
	 * @return
	 * 	True, if the {@link Triangle2d} contains the point.
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
		
		// Check against bounding circle
		double d2 = point.minus(this.circle.centre).norm2();
		if(d2 > this.circle.r2) {
			return false;
		}
		
		// Use more sophisticated algorithm to determine if it's inside the Triangle
		double v0x = v0.getX();
		double v0y = v0.getY();
		double v1x = v1.getX();
		double v1y = v1.getY();
		double v2x = v2.getX();
		double v2y = v2.getY();
		
	    double sign = Math.signum(area);
	    double s = (v0y * v2x - v0x * v2y + (v2y - v0y) * point.getX() + (v0x - v2x) * point.getY()) * sign;
	    double t = (v0x * v1y - v0y * v1x + (v0y - v1y) * point.getX() + (v1x - v0x) * point.getY()) * sign;
	    
	    return s > 0 && t > 0 && (s + t) < 2 * area * sign;
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
        
		Triangle2d t = (Triangle2d) o;
        
        // Check if the three vertices are equal
        if(v0 != t.v0 && v0 != t.v1 && v0 != t.v2) {
        	return false;
        }
        if(v1 != t.v0 && v1 != t.v1 && v1 != t.v2) {
        	return false;
        }
        if(v2 != t.v0 && v2 != t.v1 && v2 != t.v2) {
        	return false;
        }
        
        // Fields match
        return true;
    }
    

    @Override
    public int hashCode() {
        return v0.hashCode() + v1.hashCode() + v2.hashCode();
    }
    
    public String toString() {
    	return v0 + " " + v1 + " " + v2;
    }
    
    
	
    /**
     * Use line-of-sight intersection test to find the triangle lying along
     * the given line of sight direction.
     * @return
     */
//    public static <U extends Vector3d> Triangle2d<U> lookupTriangle(Vector3d line_of_sight, List<Triangle2d<U>> triangles)
//    {
//        // We will intersect two or more triangles if the LOS vector pierces
//        // the database. Need to record the minimum distance in order to find
//        // the triangle on the nearside.
//        double min_distance = Double.MAX_VALUE;
//        Triangle2d<U> target = null;
//        
//        // List of Triangles forming database...
//        for(Triangle2d<U> triangle : triangles)
//        {
//            // Check if the LOS intersects the given triangle
//            IntersectionTest intersection = new IntersectionTest(line_of_sight, triangle);
//
//            if(intersection.occurs && intersection.distance < min_distance)
//            {
//                min_distance = intersection.distance;
//                target = triangle;
//            }
//        }
//        return target;
//    }
	
}