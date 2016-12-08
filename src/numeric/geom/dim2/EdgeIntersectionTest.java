package numeric.geom.dim2;

/**
 * Utility class to examine whether two 2D edges intersect.
 *
 *
 * @author nrowell
 * @version $Id$
 */
public class EdgeIntersectionTest {

	/**
	 * Tolerance for double equals.
	 */
	private static final double epsilon = 1e-16;
	
	/**
	 * Enumerates the different configurations of edges.
	 */
	public static enum IntersectionStatus {INTERSECTING, COLLINEAR_INTERSECTING, NON_INTERSECTING, PARALLEL_NON_INTERSECTING, COLLINEAR_NON_INTERSECTING};
	
	/**
	 * The {@link IntersectionStatus}
	 */
	public final IntersectionStatus status;
	
	/**
	 * 
	 */
	public final Vector2d intersection_point;
	
	/**
	 * Main constructor.
	 * 
	 * @param e0
	 * @param e1
	 */
	public EdgeIntersectionTest(Edge2d e0, Edge2d e1) {
		
		// Aliases for vectors to make notation easier
		Vector2d p = e0.v0;
		Vector2d q = e1.v0;
		
		// End of this edge = p + r
		Vector2d r = e0.v1.minus(p);
		
		// End of that edge = q + s
		Vector2d s = e1.v1.minus(q);
		
		// Parametric line representation of each edge:
		//
		//   p + tr
		//   q + us
		//
		// where t,u are scalars in the range [0:1].
		//
		// If the two edges cross, then we can find values of t and u such that p + tr = q + us
		// The intersection occurs between the two end points (i.e. within an edge) if 0<=t<=1 and 0<=u<=1
		
		// Solve for t,u
		
		// A couple of quantities...
		double rxs = r.getX()*s.getY() - r.getY()*s.getX();
		double qmpxr = (q.getX() - p.getX())*r.getY() - (q.getY() - p.getY())*r.getX();
		
		if(Math.abs(rxs) < epsilon) {
			
			// Lines are parallel
			
			if(Math.abs(qmpxr) < epsilon) {
				// Lines are parallel and collinear: do they intersect within the range of the lines?
				
				// Express end points of the second line in terms of the first line equation
				double t0 = (q.minus(p)).dot(r) / r.norm2();
				double t1 = (q.add(s).minus(p)).dot(r) / r.norm2();
				
				// If either of these are within [0:1] range, then the lines overlap
				if((t0 >= 0.0 && t0 <= 1.0) || (t1 >= 0.0 && t1 <= 1.0)) {
					// Lines overlap
					status = IntersectionStatus.COLLINEAR_INTERSECTING;
					intersection_point = null;
				}
				else {
					// Lines don't overlap
					status = IntersectionStatus.COLLINEAR_NON_INTERSECTING;
					intersection_point = null;
				}
			}
			else {
				// Lines are parallel but not collinear
				status = IntersectionStatus.PARALLEL_NON_INTERSECTING;
				intersection_point = null;
			}
		}
		else {
			
			// Lines do intersect - but is the intersection point within the range of each line?
			double qmpxs = (q.getX() - p.getX())*s.getY() - (q.getY() - p.getY())*s.getX();
			
			double u = qmpxr / rxs;
			double t = qmpxs / rxs;
			
			// Compute the intersection point
			Vector2d intersection = p.add(r.mult(t));
			
			if(u >= 0.0 && u <= 1.0 && t >= 0.0 && t <= 1.0) {
				// Lines intersect between the end points of each
				status = IntersectionStatus.INTERSECTING;
				intersection_point = intersection;
			}
			else {
				// Lines don't intersect between the end points of each
				status = IntersectionStatus.NON_INTERSECTING;
				intersection_point = intersection;
			}
		}
	}
	
	public IntersectionStatus getIntersectionStatus() {
		return this.status;
	}
	
	public Vector2d getIntersectionPoint() {
		return this.intersection_point;
	}
	
	
	
}
