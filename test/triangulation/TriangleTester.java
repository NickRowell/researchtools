package triangulation;

/**
 * This class tests the Triangle class, in particular the geometrical functions
 * that calculate the closest point on a Triangle to the point P.
 * @author nickrowell
 */
public class TriangleTester {
    
    
    public static void main(String[] args){
        
        // Make a test Triangle
        Vertex v0 = new Vertex(1,0,0);
        Vertex v1 = new Vertex(-1,0,0);
        Vertex v2 = new Vertex(0,1,0);
        
        Triangle T = new Triangle(v0,v1,v2);
    
        // Make a test point
        Vertex P0 = new Vertex(0.0, 0.5, 1.0);
        
        // Make a test ray
        double[] R = new double[]{0, 0, -1};
        
        // Does ray cut triangle:
        boolean cuts = T.intersects(P0, R);
        System.out.println("\nRay "+(cuts ? "cuts " : "does not cut ")+"Triangle\n");
                

        System.out.println("Closest point to P0 lies at "+ T.getClosestPointToP(P0).toString());
        
        System.out.println("at distance "+T.getMinDistanceToP(P0));
   
    }
    
    
}
