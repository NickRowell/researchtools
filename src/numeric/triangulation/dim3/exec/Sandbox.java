package numeric.triangulation.dim3.exec;

import numeric.geom.dim3.Sphere;
import numeric.geom.dim3.Vector3d;

/**
 * Sandbox application.
 *
 * @author nrowell
 * @version $Id$
 */
public class Sandbox {
	
	/**
	 * Main entry point
	 * @param args
	 * 	The command line args.
	 */
	public static void main(String[] args) {
		
		
		Sphere sphere = new Sphere(0.0, 0.0, 0.0, 2.0);
		
		Vector3d[] vs = sphere.getEnclosingTetrahedronVertices();
		
		for(Vector3d v : vs) {
			
			System.out.println(v.toString());
			
		}
		
		
	}
	
}