package numeric.geom.dim3;

import java.util.List;

import Jama.Matrix;
import Jama.SingularValueDecomposition;
import numeric.triangulation.dim3.infra.Vertex3d;

/**
 * This class represents sphere objects and provides a routine to fit a sphere
 * to a set of points by a linear method.
 * 
 * @author nrowell
 * @version $Id$
 */
public class Sphere {
    
    /**
     * Coordinates of centre.
     */
    public Vector3d centre = new Vector3d(0,0,0);
    
    /**
     * Radius squared.
     */
    public double r2 = -1;
    
    /**
     * Main constructor.
     * 
     * @param X
     * 	The X coordinate of the centre of the {@link Sphere}.
     * @param Y
     * 	The Y coordinate of the centre of the {@link Sphere}.
     * @param Z
     * 	The Z coordinate of the centre of the {@link Sphere}.
     * @param R2
     * 	The radius of the {@link Sphere} squared.
     */
    public Sphere(double X, double Y, double Z, double R2) {
        centre = new Vector3d(X,Y,Z);
        r2 = R2;
    }
    
    /**
     * New sphere that is a scaled version of the given sphere.
     * 
     * @param scale_me
     * 	{@link Sphere} to be scaled
     * @param scale_factor
     * 	Radius scale factor
     * 
     */
    public Sphere(Sphere scale_me, double scale_factor) {
        centre = new Vector3d(scale_me.centre.getComponents());
        r2 = scale_me.r2 * scale_factor * scale_factor;
    }
    
    /**
     * Constructor creates the smallest sphere with given two points lying
     * on it's surface.
     * @param v1
     * 	The first {@link Vector3d}.
     * @param v2
     * 	The second {@link Vector3d}.
     */
    public Sphere(Vector3d v1, Vector3d v2) {
        // Centre lies half way along line joining points
        Vector3d a = v2.minus(v1).mult(0.5);
        r2     = a.norm2();
        centre = v1.add(a);
    }
    
    /**
     * Constructor creates the smallest sphere with given three points lying
     * on it's surface.
     * @param v1
     * 	The first {@link Vector3d}.
     * @param v2
     * 	The second {@link Vector3d}.
     * @param v3
     * 	The third {@link Vector3d}.
     */
    public Sphere(Vector3d v1, Vector3d v2, Vector3d v3) {
        Vector3d a = v2.minus(v1);
        Vector3d b = v3.minus(v1);
        
        double Denominator = 2.0f * (a.cross(b)).norm2();
        
        Vector3d o = a.cross(b).cross(a).mult(b.norm2()).add(
                     b.cross(a.cross(b)).mult(a.norm2())).mult(1.0/Denominator);
        
        r2     = o.norm2();
        centre = v1.add(o);
    }
    
    /**
     * Constructs the unique {@link Sphere} whose surface contains the four
     * {@link Vector3d}.
     * 
     * Current best implementation of sphere fitting to four points algorithm.
     * This method is about ten times faster than the deprecated method in 
     * tests using random data, and gives the same results (to within machine
     * precision errors).
     */
    public Sphere(Vector3d v1, Vector3d v2, Vector3d v3, Vector3d v4) {
        
        Vector3d a = v2.minus(v1);
        Vector3d b = v3.minus(v1);
        Vector3d c = v4.minus(v1);
        
        Matrix A = new Matrix(new double[][]{{a.getX(), a.getY(), a.getZ()},
                                             {b.getX(), b.getY(), b.getZ()},
                                             {c.getX(), c.getY(), c.getZ()}});
        
        double denominator = 2.0f * A.det();
        
        // If points are collinear or coplanar, we can't fit a sphere.
        if(denominator==0) {
        	return;
        }
        
        Vector3d o = (a.cross(b).mult(c.norm2())).add(
                     (c.cross(a).mult(b.norm2())).add(
                     (b.cross(c).mult(a.norm2())))).mult(1.0/denominator);
        
        r2     = o.norm2();
        centre = v1.add(o);
    }
    
    /**
     * Get a Sphere that encloses all of the given points. This is non-optimal
     * in the sense that the enclosing sphere is not the smallest possible.
     * Finding the smallest enclosing sphere is a much harder challenge. Here,
     * we simply find the centroid of the point set and the maximum distance to
     * any point.
     * 
     * @param points
     * 	The {@link List} of {@link Vector3d} to enclose.
     * @return
     * 	A (non-optimal) enclosing {@link Sphere} for the given set of points.
     */
    public static <T extends Vector3d> Sphere getEnclosingSphere(List<T> points) {
        // Get points centroid
        Vector3d centre = new Vector3d(0,0,0);
        for(Vector3d point : points) centre.addEquals(point);
        centre.multEquals(1.0/points.size());
        
        // Now find maximum distance of any point from centroid
        double max_rad_2 = 0;
        for(Vector3d point : points) {
            if(point.minus(centre).norm2() > max_rad_2) {
                max_rad_2 = point.minus(centre).norm2();
            }
        }
        
        return new Sphere(centre.getX(), centre.getY(), centre.getZ(), max_rad_2);
    }
    
    /**
     * Get the corners of the (regular) tetrahedron that tightly encloses this
     * {@link Sphere}, in the sense that the {@link Sphere} lies fully inside
     * the tetrahedron (e.g. the vertices of the tetrahedron don't lie on the
     * surface of the sphere).
     * Note that there is no single solution as we can rotate the
     * tetrahedron around 3 axes. The base lies in the XY plane with the
     * edge connecting the first two vertices lies along the XY direction.
     * 
     * @return
     * 	Array of four {@link Vector3d} lying at the corners of a regular tetrahedron
     * that exactly fits this {@link Sphere}.
     */
    public Vector3d[] getEnclosingTetrahedronVertices() {
    	
        // ABC form the base. Edge AB parallel to XY axis.
        double r = Math.sqrt(r2);
        
        Vector3d a = centre.add(new Vector3d(-Math.sqrt(6)*r, -Math.sqrt(2)*r, -r));
        Vector3d b = centre.add(new Vector3d( Math.sqrt(6)*r, -Math.sqrt(2)*r, -r));
        Vector3d c = centre.add(new Vector3d( 0, 2*Math.sqrt(2)*r, -r));
        
        // D lies at the peak
        Vector3d d = centre.add(new Vector3d( 0, 0, 3*r));
        
        return new Vector3d[]{a,b,c,d};
    }
    
    /**
     * This algorithm fits a sphere to a set of four points by linear least
     * squares. It is now deprecated as I found a much faster algorithm - see
     * fitSphere() method. This method may be slightly more robust however as
     * it's possible the checks on the singular values catch more special
     * cases and critical configurations of points.
     */
    public static Sphere fitSphereDeprecated(Vector3d v1, Vector3d v2, Vector3d v3, Vector3d v4) {
           
        // Fit sphere by linear least squares
        double[][] x = new double[][]{{v1.getX()*v1.getX() + v1.getY()*v1.getY() + v1.getZ()*v1.getZ(), v1.getX(), v1.getY(), v1.getZ(), 1.0},
                                      {v2.getX()*v2.getX() + v2.getY()*v2.getY() + v2.getZ()*v2.getZ(), v2.getX(), v2.getY(), v2.getZ(), 1.0},
                                      {v3.getX()*v3.getX() + v3.getY()*v3.getY() + v3.getZ()*v3.getZ(), v3.getX(), v3.getY(), v3.getZ(), 1.0},
                                      {v4.getX()*v4.getX() + v4.getY()*v4.getY() + v4.getZ()*v4.getZ(), v4.getX(), v4.getY(), v4.getZ(), 1.0}};
        
        Matrix X = new Matrix(x);
        
        Matrix D = X.transpose().times(X);
        
        // Get SVD
        SingularValueDecomposition svd = D.svd();
                        
        // Get right singular vector corresponding to lowest singular value
        Matrix rsv = svd.getV().getMatrix(new int[]{0,1,2,3,4}, new int[]{4});

        // Test for two or more coincident points, or four coplanar points.
        if(svd.getS().get(3, 3)<9E-9) {
            return null;
        }
        
        // Test for three collinear points. All other critical configurations
        // of points are detected by check on second lowest singular value,
        // however three collinear points is a funny case for which the
        // second lowest singular value is non-zero, even though no sphere
        // solution is possible.
        if(rsv.get(0, 0)==0)
            return null;
        
        // Recover canonical sphere parameters from conic coefficients
        double a = rsv.get(0, 0);
        double b = rsv.get(1, 0);
        double c = rsv.get(2, 0);
        double d = rsv.get(3, 0);
        double e = rsv.get(4, 0);
                
        double r2 = (b*b + c*c + d*d)/(4*a*a) - e/a;
        
        return new Sphere(-b/(2*a), -c/(2*a), -d/(2*a), r2);
    
    }
    
    /**
     * Determines if the given {@link Vector3d} lies within this {@link Sphere}.
     * @param vm
     * 	The {@link Vector3d} to check.
     * @return
     * 	True if the given {@link Vertex3d} lies within this {@link Sphere}.
     */
    public boolean contains(Vector3d vm)
    {
        // Check if this is smaller than the radius of the Sphere
        return (vm.minus(centre)).norm2()<r2;
    }
    
    @Override
    public String toString()
    {
        return "(x,y,z) = "+centre.toString()+" r = "+Math.sqrt(r2);
    }
    
}