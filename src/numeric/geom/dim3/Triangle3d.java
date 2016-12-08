package numeric.geom.dim3;

import java.util.List;

import Jama.Matrix;
import numeric.triangulation.dim3.infra.IntersectionTest;
import numeric.triangulation.dim3.infra.Vertex3d;


/**
 * Generic class representing triangles, for applications such as 3D graphics and Delaunay triangulations,
 * depending on what type is used.
 *
 * TODO: the intersects method is probably duplicated in the IntersectionTest class.
 *
 * @author nrowell
 * @version $Id$
 * @param <T>
 */
public class Triangle3d<T extends Vector3d> {

	/**
	 * The first vertex.
	 */
	public T v0;

	/**
	 * The second vertex.
	 */
	public T v1;

	/**
	 * The third vertex.
	 */
	public T v2;
	
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
	public Triangle3d(T v0, T v1, T v2) {
		this.v0 = v0;
		this.v1 = v1;
		this.v2 = v2;
	}
	
	/**
	 * Gets the {@link Vector3d} that is normal to the Triangle.
	 * 
	 * @return
	 */
	public Vector3d getClockwiseNormal() {
		return Vector3d.getClockwiseSurfaceNormal(v0, v1, v2);
	}

    /**
     * Two triangles are equal if they are defined by the same vertices.
     */
    public boolean equals(Triangle3d<T> t)
    {
        // Look for a single vertex not shared between the two
        if(v0 != t.v0 && v0 != t.v1 && v0 != t.v2) return false;
        if(v1 != t.v0 && v1 != t.v1 && v1 != t.v2) return false;
        if(v2 != t.v0 && v2 != t.v1 && v2 != t.v2) return false;
        
        // All vertices are common: triangles are equivalent
        return true;
    }
    
    /**
     * Use line-of-sight intersection test to find the triangle lying along
     * the given line of sight direction.
     * @return
     */
    public static <U extends Vector3d> Triangle3d<U> lookupTriangle(Vector3d line_of_sight, List<Triangle3d<U>> triangles)
    {
        // We will intersect two or more triangles if the LOS vector pierces
        // the database. Need to record the minimum distance in order to find
        // the triangle on the nearside.
        double min_distance = Double.MAX_VALUE;
        Triangle3d<U> target = null;
        
        // List of Triangles forming database...
        for(Triangle3d<U> triangle : triangles)
        {
            // Check if the LOS intersects the given triangle
            IntersectionTest intersection = new IntersectionTest(line_of_sight, triangle);

            if(intersection.occurs && intersection.distance < min_distance)
            {
                min_distance = intersection.distance;
                target = triangle;
            }
        }
        return target;
    }
	
    

    /**
     * 
     * TODO: this is probably duplicated in the IntersectionTest class.
     * 
     * Method tests if the ray from point P along vector R intersects with
     * this Triangle. Note that intersecting back-projected rays don't
     * qualify.
     * 
     */
    public boolean intersects(Vertex3d p, double[] R){
    
        // Check that R is unit vector
        double magR2 = R[0]*R[0] + R[1]*R[1] + R[2]*R[2];
        
        if(Math.abs(magR2 - 1) > 9E-9){
            System.err.println("R is not a unit vector! |R| != 1");
            System.exit(1);
        }
        
        double[] E1 = new double[]{v1.getX() - v0.getX(), v1.getY() - v0.getY(), v1.getZ() - v0.getZ()};
        double[] E2 = new double[]{v2.getX() - v0.getX(), v2.getY() - v0.getY(), v2.getZ() - v0.getZ()};        
        
        //System.out.println("E1 = ["+E1[0]+" "+E1[1]+" "+E1[2]+"]");
        //System.out.println("E2 = ["+E2[0]+" "+E2[1]+" "+E2[2]+"]");
        
        // Check that R is not coplanar with Triangle
        double[] E1crossE2 = new double[]{E1[1]*E2[2] - E1[2]*E2[1],
                                          E1[2]*E2[0] - E1[0]*E2[2],
                                          E1[0]*E2[1] - E1[1]*E2[0]};
        
        //System.out.println("E1crossE2 = ["+E1crossE2[0]+" "+E1crossE2[1]+" "+E1crossE2[2]+"]");        
        
        double E1crossE2dotR = E1crossE2[0]*R[0] + E1crossE2[1]*R[1] + E1crossE2[2]*R[2];
        
        if(Math.abs(E1crossE2dotR) < 9E-9){
            // R is coplanar with triangle
            System.out.println("R coplanar!");
            return false;
        }
        
        // R definitely cuts plane of triangle. Find out where.
        
        Matrix E0E1R = new Matrix(new double[][]{{E1[0], E2[0], -R[0]},
                                                 {E1[1], E2[1], -R[1]},
                                                 {E1[2], E2[2], -R[2]}});
        
        Matrix D = new Matrix(new double[][]{{p.getX()-v0.getX()}, {p.getY()-v0.getY()}, {p.getZ()-v0.getZ()}});
        
        Matrix sta = E0E1R.solve(D);
        
        double s = sta.get(0, 0);
        double t = sta.get(1, 0);
        double a = sta.get(2, 0);

        System.out.println("s = "+s+"\nt = "+t+"\na = "+a);
        
        // Check value of s,t
        if(s >=0 && s <=1 && t >=0 && t <=1){
            
            // Ray cuts plane within Triangle area. Check for forward
            // projection
            if(a >= 0){
            
                // Forward projected ray cuts Triangle
                return true;
            }
        
        }
        
        return false;
        
    }
    
    
    
    
    

    /**
     * Get the closest point in the Triangle to the point P. Method uses
     * Lagrange multipliers to find minimum of distance function subject to
     * constraint that point lies in triangle domain.
     */
    public Vector3d getClosestPointToP(Vector3d p) {
        
        // In triangle equation, vertex 0 provides origin B, vertex
        // 1 defines direction of E1 and vertex 2 defines direction of E2.
        Vector3d D = p.minus(v0);
        Vector3d E1 = v1.minus(v0);
        Vector3d E2 = v2.minus(v0);
        
//        double a = D.dot(D);
        double b = D.dot(E1);
        double c = D.dot(E2);
        double d = E1.dot(E2);
        double e = E1.dot(E1);
        double f = E2.dot(E2);
        
        // Solution for global minimum s',t'
        double s = (f*b-d*c)/(f*e-d*d);
        double t = (e*c-d*b)/(e*f-d*d);
                
        // Which region of s,t plane does this fall in?
        
        // Region 6
        if(s < 0 && t < 0)
        {
            // Region 6
            
            // Solution for s',t' must satisfy at least one of the two
            // constraints s=0 or t=0.
            
            // Find solution for t' on applying s=0 constraint:
            t = c/f;
            
            // Find solution for s' on applying t=0 constraint:
            s = b/e;           
            
            // If t' is greater than zero and s' less than zero, then solution 
            // lies on s=0 boundary alone.
            if(t >= 0 && s < 0)
            {
                
                // Clamp t_prime to allowed range
                if(t > 1) t = 1;
            
                s = 0;
            }
            
            // If s' is greater than zero, then solution lies on t=0
            // boundary alone.
            else if(s >= 0 && t < 0)
            {
               
                // Clamp s_prime to allowed range
                if(s > 1) s = 1;
            
                t = 0;

            }
            // If both are negative, then solution lies on corner.
            else if(s <= 0 && t <= 0)
            {
                s = 0;
                t = 0;
            }
            else
            {
                // Both positive. This should not happen so indicates and
                // error if this code is reached.
                System.err.println("s' and t' both positive!");
                System.exit(1);
            }
                
                
                
            
        }
        // Region 2        
        else if(s < 0 && t < 1 - s)
        {
            // Region 2
            
            // Solution for minimum t, constrained to lie along the 
            // boundary s=0;
            t = c/f;
            
            // Clamp t_prime to allowed range
            if(t < 0) t = 0;
            if(t > 1) t = 1;
            
            s = 0;
            
        
        }
        // Region 4        
        else if(s < 0 && t >= 1 - s)
        {
            // Region 4
            
            // Solution for s',t' must satisfy at least one of the two
            // constraints s=0 or s+t=1.
            
            // Find solution for t' on applying s=0 constraint:
            t = c/f;
            
            // Find solution for s' on applying s+t=1 constraint:
            s = (c - b + d - f)/(2*d - e - f);
            
            // If t' is greater than 1, then solution constrained to lie
            // along s+t=1 boundary (possibly at either end), determined by
            // solution for s'.
            if(t >= 1)
            {
                // If s' is < 0, solution lies at top right corner
                if(s < 0)
                {
                    s = 0.0;
                    t = 1.0;
                }
                // If 0 < s' < 1 then s' is a valid solution. t' constrained to
                // equal 1-s in this case.
                else if(s <=1)
                {
                    t = 1.0 - s;
                }
                // Solution lies at bottom right corner of triangle
                else if(s > 1)
                {
                    s = 1.0;
                    t = 0.0;
                }
                else
                {
                    // This should not happen so indicates an
                    // error if this code is reached.
                    System.err.println("Error!");
                    System.exit(1);                
                }
                    
                
            }
            // If t' is in [0,1] range, then it is a valid solution. s'
            // is constrained to equal zero in this case.
            else if(t >= 0 && t < 1)
            {
                s = 0;
            }
            // If t' is less than zero, then solution lies at corner of
            // triangle. Both s' and t' are zero.
            else if(t < 0)
            {
                s = 0;
                t = 0;
            }  
            else
            {
                // This should not happen so indicates an
                // error if this code is reached.
                System.err.println("Error!");
                System.exit(1);
            }
            
        }
        // Region 3        
        else if(t < 0 && s < 1 - t)
        {
            // Region 3
            
            // Solution for minimum s, constrained to lie along the 
            // boundary t=0;
            s = b/e;
            
            // Clamp s_prime to allowed range
            if(s < 0) s = 0;
            if(s > 1) s = 1;
            
            t = 0;
            
        }
        // Region 5        
        else if(t < 0 && s >= 1 - t)
        {
            // Region 5
            
            // Solution for s',t' must satisfy at least one of the two
            // constraints t=0 or s+t=1.
            
            // Find solution for s' on applying t=0 constraint:
            s = b/e;
            
            // Find solution for t' on applying s+t=1 constraint:
            t = 1.0 - (c - b + d - f)/(2*d - e - f);
            
            // If s' is greater than 1, then solution constrained to lie
            // along s+t=1 boundary (possibly at either end), determined by
            // solution for t'.
            if(s >= 1)
            {
                // If t' is < 0, solution lies at bottom right corner
                if(t < 0)
                {
                    s = 1.0;
                    t = 0.0;
                }
                // If 0 < t' < 1 then t' is a valid solution. s' constrained to
                // equal 1-t in this case.
                else if(t <=1)
                {
                    s = 1.0 - t;
                }
                // Solution lies at top left corner of triangle
                else if(t > 1)
                {
                    s = 0.0;
                    t = 1.0;
                }
                else
                {
                    // This should not happen so indicates an
                    // error if this code is reached.
                    System.err.println("Error!");
                    System.exit(1);                
                }
                    
                
            }
            // If s' is in [0,1] range, then it is a valid solution. t'
            // is constrained to equal zero in this case.
            else if(s >= 0 && s < 1)
            {
                t = 0;
            }
            // If s' is less than zero, then solution lies at corner of
            // triangle. Both s' and t' are zero.
            else if(s < 0)
            {
                s = 0;
                t = 0;
            }
            else
            {
                // This should not happen so indicates an
                // error if this code is reached.
                System.err.println("Error!");
                System.exit(1);
            }
            
        }
        // Region 0       
        else if(t < 1 - s)
        {
            // Region 0.
            
            // Global minimum of distance lies inside Triangle, so solution
            // for s_prime, t_prime is correct as it stands.
        
        }
        // Region 1        
        else if(t >= 1 - s)
        {
            // Region 1.
            
            // s-coordinate of minimum, constrained to lie on line
            // s + t = 1
            s = (c - b + d - f)/(2*d - e - f);
            
            // Clamp s_prime to permitted range
            if(s < 0) s = 0;
            if(s > 1) s = 1;
            
            t = 1 - s;
        
        }
        else
        {
            // Sanity check. The above logic should cover all possible
            // s_prime,t_prime cases, so if program ends up here then
            // something has screwed up.
            System.err.println("Illegal region in plane of Triangle!");
            System.exit(1);
            
        }
        
        // Now have location of minimum distance in terms of s', t'.
        // Get coordinates of this point.
        //
        // T = v0 + s*E1 + t*E2
        //
        return v0.add(E1.mult(s)).add(E2.mult(t));
        
    }
}
