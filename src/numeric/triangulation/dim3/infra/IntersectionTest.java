package numeric.triangulation.dim3.infra;

import numeric.geom.dim3.Triangle3d;
import numeric.geom.dim3.Vector3d;

/**
 * Class represents the result of a test to check if a particular Vector3d
 * intersects a Triangle, and at what distance along the vector.
 * 
 * TODO: could use enum to indicate the result of the intersection test.
 * 
 * @author nickrowell
 */
public class IntersectionTest {
	
    /**
     * Does the vector intersect the triangle?.
     */
    public boolean occurs;
    
    /**
     * At what distance does the intersection occur?.
     */
    public double distance;
    
    /**
     * Indicate pathological case where LOS ray lies in plane of triangle.
     */
    public boolean COPLANAR;
    
    /**
     * Main constructor.
     * 
     * @param p
     * 	{@link Vector3d} to be tested for intersection with the triangle. It passes through the origin.
     * The ray equation is r = (0,0,0) + d*p
     * @param t
     * 	Triangle. We use the vertex positions expressed in the camera
     * frame.
     */
    public <T extends Vector3d> IntersectionTest(Vector3d p, Triangle3d<T> t) {
                
        // Check that p is unit vector
        if(Math.abs(p.norm2() - 1) > 9E-9)
        {
            System.err.println("p is not a unit vector! |p| != 1");
            return;
        }
        
        // Construct normalised basis vectors for non-orthogonal triangle
        // frame, with vertex v0 at the origin.
        
        Vector3d E1 = t.v1.minus(t.v0);
        // Length of this vector
        double e1 = E1.norm();
        // Convert to unit basis vector for further calculations
        E1.multEquals(1/e1);
        
        Vector3d E2 = t.v2.minus(t.v0);
        // Length of this vector
        double e2 = E2.norm();
        // Convert to unit basis vector for further calculations
        E2.multEquals(1/e2);
        
        // Complete the set: this vector forms the surface normal for the triangle
        Vector3d E3 = E1.cross(E2);
        
        if(Math.abs(E3.dot(p)) < 9E-9)
        {
            // R is coplanar with triangle
            COPLANAR = true;
            occurs = false;
            distance = 0.0;
            return;
        }
        else
        {
            COPLANAR = false;
        }
        
        // Now solve for plane equation of triangle, in camera frame. The plane
        // equation is
        //
        //  N.r - d = 0
        //
        // where N is the plane normal, d is the minimum distance to the plane
        // and r is any vector lying in the plane. We know N currently (E3), 
        // and have three examples of vectors lying in the plane (vertices of 
        // the triangle), so we can solve for d like this:
        double d = E3.dot(t.v0);
        
        // Equivalent values:
        //double d = E3.dot(t.v1.camP);
        //double d = E3.dot(t.v2.camP);
        
        // Now use the plane equation to find the coordinates of the point
        // where the LOS ray intersects the triangle plane.
        
        // Distance along p where the intersection occurs. Note that this
        // may be negative, in which case the intersection occurs in the
        // backwards direction along the LOS vector.
        distance = d / E3.dot(p);
        
        // Coordinates of intersection point, in camera frame
        Vector3d I_cam = p.mult(distance);
        
        // Now get the position vector of the intersection point relative to 
        // the origin of the triangle frame, expressed in the camera frame.
        Vector3d I = I_cam.minus(t.v0);
        
        // In order to determine whether this lies within the boundaries of the
        // triangle, we first need to express it as a linear combination of the
        // triangle basis vectors E1 and E2. We parameterise I as:
        //
        // I = alpha*E1 + beta*E2
        //
        
        double E1E2 = E1.dot(E2);    // intermediate quantity
        
        double alpha = I.dot(E1.minus(E2.mult(E1E2)))/(1 - E1E2*E1E2);
        double beta  = I.dot(E2.minus(E1.mult(E1E2)))/(1 - E1E2*E1E2);
        
        // If either alpha or beta are negative, then the ray intersects the
        // triangle plane outside the boundaries of the triangle.
        if(alpha < 0 || beta < 0)
        {
            occurs = false;
        }
        
        // Now, the point I lies between the legs of the triangle, so we need
        // to test if it is sufficiently close to the apex that it is within
        // the boundary opposite the apex.
        // 
        // Points on the line have the property that
        //
        //  alpha/e1 + beta/e2 = 1
        // 
        // If alpha/e1 + beta/e2 > 1 then we have crossed the line, and I lies
        // outside of the triangle.
        //
        else if(alpha/e1 + beta/e2 > 1)
        {
            occurs = false;
        }
        
        // If we have reached this point, then all tests have been passed and
        // the LOS vector does intersect the triangle inside the boundaries.
        else
        {
            occurs = true;
        }
    }
    
}