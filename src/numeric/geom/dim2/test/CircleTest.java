package numeric.geom.dim2.test;

import numeric.geom.dim2.Circle;
import numeric.geom.dim2.Vector2d;

/**
 * Tests for the {@link Circle} class.
 *
 *
 * @author nrowell
 * @version $Id$
 */
public class CircleTest {
	
	/**
	 * Main application entry point.
	 * 
	 * @param args
	 * 	The command line arguments (ignored).
	 */
	public static void main(String[] args) {
		
		Vector2d a = new Vector2d(-Math.sqrt(3.0), 1.5);
		Vector2d b = new Vector2d(Math.sqrt(3.0), 1.5);
		Vector2d c = new Vector2d(0, 4.5);
		
		Circle circ = new Circle(a, b, c);
		
		System.out.println("Circle centre = "+circ.centre.toString());
		System.out.println("Circle radius = " + Math.sqrt(circ.r2));
		
		Circle circ2 = new Circle(0, 1, 1);
		
		Vector2d[] vs = circ2.getEnclosingTriangleVertices();
		
		System.out.println("Enclosing vertices:");
		for(Vector2d v : vs) {
			System.out.println(v.toString());
		}
		
		
	}
	
}