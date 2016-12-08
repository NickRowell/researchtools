package numeric.geom.dim2.test;

import numeric.geom.dim2.Edge2d;
import numeric.geom.dim2.EdgeIntersectionTest;
import numeric.geom.dim2.Vector2d;

/**
 * Tests for the {@link Edge2d} class.
 *
 *
 * @author nrowell
 * @version $Id$
 */
public class Edge2dTest {
	
	/**
	 * Main constructor.
	 * 
	 * @param args
	 * 	The command line args (ignored)
	 */
	public static void main(String[] args) {
		
		Vector2d v0 = new Vector2d(-1, 0);
		Vector2d v1 = new Vector2d( 1, 0);
		Vector2d v2 = new Vector2d( 0, 1);
		Vector2d v3 = new Vector2d( 0,-1);
//		Vector2d v4 = new Vector2d( 1, 1);
//		Vector2d v5 = new Vector2d( 2, 2);
		Vector2d v6 = new Vector2d( 2, 1);
		Vector2d v7 = new Vector2d( 2,-1);
		Vector2d v8 = new Vector2d( 0,-2);
		Vector2d v9 = new Vector2d( 0,-3);
		
		Edge2d e0 = new Edge2d(v0, v1);
		Edge2d e1 = new Edge2d(v2, v3);
//		Edge2d e2 = new Edge2d(v4, v5);
		Edge2d e3 = new Edge2d(v6, v7);
		Edge2d e4 = new Edge2d(v2, v8);
		Edge2d e5 = new Edge2d(v8, v9);
		
		
		EdgeIntersectionTest t1 = new EdgeIntersectionTest(e0, e1);
		System.out.println("Intersection of e0 & e1:");
		System.out.println("status             = "+t1.status);
		System.out.println("intersection point = "+t1.intersection_point);
		
		EdgeIntersectionTest t2 = new EdgeIntersectionTest(e0, e3);
		System.out.println("Intersection of e0 & e3:");
		System.out.println("status             = "+t2.status);
		System.out.println("intersection point = "+t2.intersection_point);

		EdgeIntersectionTest t3 = new EdgeIntersectionTest(e1, e3);
		System.out.println("Intersection of e1 & e3:");
		System.out.println("status             = "+t3.status);
		System.out.println("intersection point = "+t3.intersection_point);

		EdgeIntersectionTest t4 = new EdgeIntersectionTest(e1, e4);
		System.out.println("Intersection of e1 & e4:");
		System.out.println("status             = "+t4.status);
		System.out.println("intersection point = "+t4.intersection_point);

		EdgeIntersectionTest t5 = new EdgeIntersectionTest(e1, e5);
		System.out.println("Intersection of e1 & e5:");
		System.out.println("status             = "+t5.status);
		System.out.println("intersection point = "+t5.intersection_point);
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
	}
	
}
