package numeric.triangulation.dim2.exec;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import numeric.geom.dim2.Triangle2d;
import numeric.geom.dim2.Vector2d;
import numeric.triangulation.dim2.base.VoronoiTessellation2D;

/**
 * Tests for the voronoi tessellation methods.
 *
 *
 * @author nrowell
 * @version $Id$
 */
public class TestVoronoi {
	
	/**
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		
		// Central vertex
		Vector2d v0 = new Vector2d(0,0);
		// Outer vertices
		Vector2d v1 = new Vector2d(-1,-1);
		Vector2d v2 = new Vector2d(-1,1);
		Vector2d v3 = new Vector2d(1,1);
		Vector2d v4 = new Vector2d(1,-1);
		
		Triangle2d t0 = new Triangle2d(v0, v1, v2);
		Triangle2d t1 = new Triangle2d(v0, v2, v3);
		Triangle2d t2 = new Triangle2d(v0, v3, v4);
		Triangle2d t3 = new Triangle2d(v0, v4, v1);
		
		Set<Triangle2d> tris = new HashSet<>();
		
		tris.add(t0);
		tris.add(t1);
		tris.add(t2);
		tris.add(t3);
		
		List<Triangle2d> trisList = VoronoiTessellation2D.getOrderedList(tris, v0);
		
		for(Triangle2d tri : trisList) {
			System.out.println(tri);
		}
		
		
	}
	
	
}
