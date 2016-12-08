package numeric.triangulation.dim3.impl;

import java.util.List;

import numeric.geom.dim3.Sphere;
import numeric.triangulation.dim3.base.DelaunayTriangulation3D;
import numeric.triangulation.dim3.infra.Vertex3d;

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
public class DelaunayBruteForce extends DelaunayTriangulation3D {
    
    /**
     * Main constructor.
     * 
     * @param pvertices
     * 	The {@link List} of {@link Vertex3d} to triangulate.
     */
    public DelaunayBruteForce(List<Vertex3d> pvertices) {   
        super(pvertices);
    }
    
    /**
     * Instance of 3D Delaunay Triangulation algorithm that operates on the
     * internal List of Clusters, and writes status messages etc to the
     * PipedOutputStream passed in as argument.
     * @param _clusters
     * @param os
     */
    @Override
    public void doTriangulation() {
    	
        // Number of trials in brute force Delaunay triangulation = number of
        // ways to select 4 objects from N.
//        double N_TRIALS = Statistics.nChooseM(clusters.size(), 4);
        // Counts all trials
//        double N = 0;
        // Number of trials in one percent
//        double N_TRIALS_100 = N_TRIALS/100;
        // Counts trials to one percent then resets to zero
//        double N_100 = 0;
        // Calculate some stats to enable time-to-completion estimate
//        long T_0 = System.currentTimeMillis();
        
        // Proceed through all possible combinations of four Vertex objects
        for(int i=0; i<vertices.size()-3; i++)
        {
            for(int j=i+1; j<vertices.size()-2; j++)
            {
                for(int k=j+1; k<vertices.size()-1; k++)
                {
                    for(int l=k+1; l<vertices.size(); l++)
                    {
                        
                        triangulateSimplex(i,j,k,l);
                        
//                        N_100++;
//                        N++;
//                        
//                        if(N_100>N_TRIALS_100)
//                        {
//                            // Percentage complete
//                            int complete = (int)Math.floor(N*100/N_TRIALS);
//                            
//                            // Approx. time left [ms] based on number of trials
//                            // multiplied by average time per trial to this point
//                            double T_REMAINING_MS = ((N_TRIALS - N) * (System.currentTimeMillis() - T_0)) / N;
//                            
//                            N_100=0;
//                            
//                        }
                    }
                }
            }
        }
        
        // Triangulation is complete. Now do some housekeeping to identify
        // triangles and vertices that do not lie on the outer hull.
        setExternalHull();
        
    }
    
    
    /**
     * Check if a set of four vertices form a simplex (tetrahedron) in the
     * triangulation.
     */
    private void triangulateSimplex(int i, int j, int k, int l)
    {
        
    	Vertex3d v1 = vertices.get(i);
        Vertex3d v2 = vertices.get(j);
        Vertex3d v3 = vertices.get(k);
        Vertex3d v4 = vertices.get(l);
        
        // Create new Sphere using these vertices
        Sphere sphere = new Sphere(v1,v2,v3,v4);

        // If null Sphere is returned, then vertices are in
        // critical configuration, Sphere fit is impossible
        // and triangulation is potentially not unique.
        // Skip this trial.
        if(sphere==null) return;

        // Now check for other vertices that lie within this 
        // sphere, which would indicate that these four
        // vertices don't form a simplex in the triangulation.
        for(int m=0; m<vertices.size(); m++)
        {
            // Ignore vertices that are part of the set of four
            if(m==i||m==j||m==k||m==l) continue;
            
            // Does vertex m lie inside the Sphere?
            else if(sphere.contains(vertices.get(m))) return;
        }
        
        // If this point is reached, then the set of four
        // vertices currenty selected do form part of the
        // Delaunay triangulation. 
        addTetrahedron(v1, v2, v3, v4, sphere);
    
    }    
    
}