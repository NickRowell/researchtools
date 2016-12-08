package numeric.geom.test;

import numeric.geom.dimN.VectorNd;

/**
 * Class tests the {@link VectorNd} class.
 *
 *
 * @author nrowell
 * @version $Id$
 */
public class TestVectorNd {
	
	/**
	 * Main application entry point.
	 * @param args
	 * 	The command line args.
	 */
	@SuppressWarnings("rawtypes")
	public static void main(String[] args) {
		
		VectorNd a = new VectorNd(1.0, 0.0);
		VectorNd b = new VectorNd(0.0, 1.0);
		VectorNd c = new VectorNd(0.0, 1.0);
		
		System.out.println("a = " + a.toString());
		System.out.println("b = " + b.toString());
		System.out.println("c = " + c.toString());
		
		System.out.println("a.equals(b) = "+a.equals(b));
		
		System.out.println("a.equals(c) = "+a.equals(c));
		System.out.println("b.equals(c) = "+b.equals(c));
		
		System.out.println("a.isParallelTo(b) = " +a.isParallelTo(b));
		System.out.println("a.isParallelTo(c) = " +a.isParallelTo(c));
		System.out.println("b.isParallelTo(c) = " +b.isParallelTo(c));
		
		
	}
	
	
}
