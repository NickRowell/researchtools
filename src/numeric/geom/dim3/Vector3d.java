package numeric.geom.dim3;

import numeric.geom.base.BaseVectorNd;

/**
 * Class represents vectors in three dimensions, and provides some methods unique to them.
 *
 * @author nrowell
 * @version $Id$
 */
public class Vector3d extends BaseVectorNd<Vector3d> {
	
	/**
	 * Main constructor.
	 * 
	 * @param components
	 * 	The vector components.
	 */
	public Vector3d(double... components) {
    	super(components);
		if(components.length!=3) {
			throw new IllegalArgumentException("Wrong number of components! Expected 3, found "+components.length);
		}
    }
	
	/**
	 * {@inheritDoc}
	 */
	@Override
    public Vector3d newInstance(double... components) {
    	return new Vector3d(components);
    }

    /**
     * Set x component
     * @param x
     * 	x component
     */
    public final void setX(double x) {
    	this.components[0] = x;
    }
    
    /**
     * Set y component
     * @param y 
     * 	y component
     */
    public final void setY(double y) {
    	this.components[1] = y;
    }
    
    /**
     * Set z component
     * @param z 
     * 	z component
     */
    public final void setZ(double z) {
    	this.components[2] = z;
    }
    
    /**
     * Get x component
     * @return
     * 	The x component
     */
    public double getX() {
    	return this.components[0];
    }
    
    /**
     * Get y component
     * @return
     * 	The y component
     */
    public double getY() {
    	return this.components[1];
    }
    
    /**
     * Get z component
     * @return
     * 	The z component
     */
    public double getZ() {
    	return this.components[2];
    }

    /**
     * Cross product of two vectors - unique operation for 3D vectors.
     * @param b
     * 	Second vector
     * @return
     * 	this cross b
     */
    public Vector3d cross(Vector3d b) {
    	return new Vector3d(getY() * b.getZ() - getZ() * b.getY(),
                            getZ() * b.getX() - getX() * b.getZ(),
                            getX() * b.getY() - getY() * b.getX());
    }

    /**
     * Get a random unit {@link Vector3d} uniformly distributed on the unit sphere.
     * 
     * @return
     * 	A random unit {@link Vector3d} uniformly distributed on the unit sphere.
     */
    public static Vector3d getRandVecOnUnitSphere() {
    	
        // Random longitude
        double a = 2 * Math.PI * Math.random();
        // Random latitude
        double d = Math.asin(2*Math.random() - 1);
        
        double x = Math.cos(d) * Math.cos(a);
        double y = Math.cos(d) * Math.sin(a);
        double z = Math.sin(d);
        
        return new Vector3d(x,y,z);
    }
    
    /**
     * For the given input {@link Vector3d}, this method returns a unit {@link Vector3d} with
     * the property that it points within a cone centred on the input vector and with an opening
     * angle specified by the second parameter.
     * 
     * Tip of random vector is uniformly distributed on the surface of the unit sphere (within
     * the cone opening angle).
     * 
     * @param in
     * 	Input {@link Vector3d} (defines axis of cone)
     * @param opening
     * 	Cone opening half-angle [radians]
     * @return
     * 	A unit vector randomly distributed on the cap of a cone aligned with the
     * input {@link Vector3d} and with an opening angle equal to the given angle.
     */
    public static Vector3d getRandVecOnUnitSphereWithinCone(Vector3d in, double opening) {
    	
        // First, need to define the rotation matrix that aligns the
        // +z axis with the direction of the input vector.
        //
        // The basis vectors for the rotated frame are e1, e2, e3.
        //
        // e3 - points along direction of input vector
        // e2 - lies in XY plane of non-rotated frame
        // e1 - completes the set
        //
        Vector3d e3 = in.normalise();
        
        Vector3d e2 = new Vector3d(0,1,0);
        Vector3d e1 = new Vector3d(1,0,0);
        
        if(e3.isParallelTo(new Vector3d(0,0,1))) {
            // No action
        }
        else if(e3.isParallelTo(new Vector3d(0,0,-1))) {
            // Reflect basis vector to ensure a right hand set
            e2 = new Vector3d(0,-1,0);
        }
        else {
            // General case
            
            // Get projection of e3 axis in xy plane, rotated by -90 deg.
            e2 = new Vector3d(-e3.getY(), e3.getX(), 0);
            // Normalise to get basis vector
            e2.normaliseInPlace();
            
            // e1 completes set
            e1 = e2.cross(e3);
        }
        
        // Now draw a random unit vector in the rotated frame, lying within
        // the given angle of the +z_prime axis
        
        // Random longitude
        double a = 2 * Math.PI * Math.random();
        
        // Random latitude within cone opening angle
        double d = Math.acos(1 - (1-Math.cos(opening))*Math.random());
        
        // Components of vector in rotated frame
        double x_prime = Math.sin(d) * Math.cos(a);
        double y_prime = Math.sin(d) * Math.sin(a);
        double z_prime = Math.cos(d);
        
        // Components of random vector in original frame
        double x = e1.getX() * x_prime + e2.getX() * y_prime + e3.getX() * z_prime;
        double y = e1.getY() * x_prime + e2.getY() * y_prime + e3.getY() * z_prime;
        double z = e1.getZ() * x_prime + e2.getZ() * y_prime + e3.getZ() * z_prime;
        
        return new Vector3d(x, y, z);
    }
    
    /**
     * Calculates the surface normal for the facet produced by joining vertices
     * r0, r1 and r2 in clockwise winding order.
     * 
     * @param r0
     * 	The first {@link Vector3d} defining a vertex of the triangle.
     * @param r1
     * 	The second {@link Vector3d} defining a vertex of the triangle.
     * @param r2
     * 	The third {@link Vector3d} defining a vertex of the triangle.
     * @return
     * 	The surface normal for the facet
     */
    public static Vector3d getClockwiseSurfaceNormal(Vector3d r0, Vector3d r1, Vector3d r2) {
    	
        // Vector from v0 to v1
        double[] a = {r1.getX()-r0.getX(), r1.getY()-r0.getY(), r1.getZ()-r0.getZ()};
        // Vector from v0 to v2
        double[] b = {r2.getX()-r0.getX(), r2.getY()-r0.getY(), r2.getZ()-r0.getZ()};
        
        // Cross product a x b gives normal direction
        double n0 = (a[1]*b[2] - a[2]*b[1]);
        double n1 = (a[2]*b[0] - a[0]*b[2]);
        double n2 = (a[0]*b[1] - a[1]*b[0]);
        
        // Normalise vector.
        double norm = Math.sqrt(n0*n0 + n1*n1 + n2*n2);
        
        n0 /= norm;
        n1 /= norm;
        n2 /= norm; 
    
        return new Vector3d(n0, n1, n2);
    }
}