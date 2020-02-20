package numeric.triangulation.dim2.impl;

import java.util.List;

import numeric.geom.dim2.Circle;
import numeric.geom.dim2.Vector2d;
import numeric.stats.StatUtil;
import numeric.triangulation.dim2.base.DelaunayTriangulation2D;

/**
 * Brute-force Delaunay triangulation algorithm implementation.
 * 
 * There are three basic objects that define the triangulation: Vertex,
 * Triangle and Tetrahedron. Only Vertex object contain data: these contain
 * the three coordinates defining the vertex 3D position. Each object contains
 * two lists; these hold references to objects of each of the other two types.
 * For example, Vertex objects hold a list of references to Triangles and 
 * Tetrahedron objects that they form part of in the final triangulation.
 * Triangle objects hold a list of (3) Vertex objects that define their corners
 * and a list of (one or two) Tetrahedra that they form part of.
 * 
 * The objects are thus divided into a geometrical hierarchy:
 * 
 * Vertex -> Triangle -> Tetrahedron
 * 
 * where each object holds references to the objects that define it (backwards)
 * and the objects that it forms part of (forwards). Only references are stored
 * so that the same object in different lists can be identified.
 * 
 * The main class DelaunayTriangulation3D contains the overall lists of each 
 * object in the triangulation, and does the work of identifying the
 * connections between them, starting out with just the Vertex objects.
 * 
 * @author nrowell
 * @version $Id$
 */
public class DelaunayBruteForce extends DelaunayTriangulation2D {
    
	/**
	 * Controls logging from main triangulation loop.
	 */
	private static boolean withLogging = false;
	
    /**
     * Main constructor.
     * 
     * @param vertices
     * 	The {@link List} of {@link Vertex2d} to triangulate.
     */
    public DelaunayBruteForce(List<Vector2d> vertices) {
        super(vertices);
    }
    
    /**
     * Instance of 2D Delaunay Triangulation algorithm that operates on the internal List of Vector2d.
     */
    @Override
    public void doTriangulation() {
    	
    	double N_TRIALS = 0.0, N = 0.0, N_TRIALS_100 = 0.0, N_100 = 0.0;
    	long T_0 = 0L;
    	
    	if(withLogging) {
	        // Number of trials in brute force Delaunay triangulation = number of
	        // ways to select 3 objects from N.
	        N_TRIALS = StatUtil.nChooseM(vertices.size(), 3);
	        // Counts all trials
	        N = 0;
	        // Number of trials in one percent
	        N_TRIALS_100 = N_TRIALS/100;
	        // Counts trials to one percent then resets to zero
	        N_100 = 0;
	        // Calculate some stats to enable time-to-completion estimate
	        T_0 = System.currentTimeMillis();
    	}
        
        // Proceed through all possible combinations of three Vector2d objects

        for(int i=0; i<vertices.size()-2; i++) {
            for(int j=i+1; j<vertices.size()-1; j++) {
                for(int k=j+1; k<vertices.size(); k++) {
                    
                    triangulateSimplex(i, j, k);
                    
                    if(withLogging) {
	                    N_100++;
	                    N++;
	                    
	                    if(N_100>N_TRIALS_100) {
	                        // Percentage complete
	                        int complete = (int)Math.floor(N*100/N_TRIALS);
	                        
	                        // Approx. time left [ms] based on number of trials
	                        // multiplied by average time per trial to this point
	                        double T_REMAINING_MS = ((N_TRIALS - N) * (System.currentTimeMillis() - T_0)) / N;
	                        
	                        System.out.println(complete+"% complete; "+(T_REMAINING_MS/1000.0)+" [sec] remaining...");
	                        
	                        N_100=0;
	                    }
                    }
                    
                }
            }
        }
    
        // Triangulation is complete. Now do some housekeeping to identify
        // triangles and vertices that do not lie on the outer hull.
        setExternalHull();
    }
    
    
    /**
     * Check if a set of three vertices form a simplex (Triangle) in the triangulation,
     * and add it to the triangulation if so.
     * 
     * @param i
     * 	Index of the first {@link Vector2d} in the internal list.
     * @param j
     * 	Index of the second {@link Vector2d} in the internal list.
     * @param k
     * 	Index of the third {@link Vector2d} in the internal list.
     */
    private void triangulateSimplex(int i, int j, int k) {
        
    	Vector2d v1 = vertices.get(i);
        Vector2d v2 = vertices.get(j);
        Vector2d v3 = vertices.get(k);
        
        // Create new Circle using these vertices
        Circle circle = null;
        try {
        	circle = new Circle(v1, v2, v3);
        }
        catch(RuntimeException e) {
        	// Points are in a critical configuration that cannot be
        	// fitted by a circle; triangulation is potentially not
        	// unique.
        	return;
        }

        // Now check for other vertices that lie within this 
        // circle, which would indicate that these three
        // vertices don't form a simplex in the triangulation.
        for(int m=0; m<vertices.size(); m++) {
        	
            // Ignore vertices that are part of the set of three
            if(m==i||m==j||m==k) {
            	continue;
            }
            
            // Does vertex m lie inside the Circle?
            else if(circle.contains(vertices.get(m))) {
            	return;
            }
        }
        
        // If this point is reached, then the set of three
        // vertices currently selected do form part of the
        // Delaunay triangulation. 
        addTriangle(v1, v2, v3);
    }    
    
}