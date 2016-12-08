package numeric.geom.dim2;

import java.util.List;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

/**
 * This class represents circles and provides fitting routines.
 * 
 * TODO: move remaining Circle methods from the other class.
 * TODO: detect critical configurations of points in fitting method.
 * 
 * 
 * @author nrowell
 * @version $Id$
 */
public class Circle {
    
	/**
	 * Precision level for comparing floating point numbers.
	 */
	private static final double epsilon = 1e-9;
	
	
    /**
     * Coordinates of centre.
     */
    public Vector2d centre;
    
    /**
     * Radius squared.
     */
    public double r2;
    
    /**
     * Main constructor for the {@link Circle}.
     * 
     * @param x
     * 	The X coordinate of the centre of the {@link Circle}.
     * @param y
     * 	The Y coordinate of the centre of the {@link Circle}.
     * @param r2
     * 	The radius of the {@link Circle} squared.
     */
    public Circle(double x, double y, double r2) {
    	this.centre = new Vector2d(x, y);
        this.r2 = r2;
    }
    
    /**
     * New {@link Circle} that is a scaled version of the given {@link Circle}.
     * 
     * @param scale_me
     * 	{@link Circle} to be scaled
     * @param scale_factor
     * 	Radius scale factor
     * 
     */
    public Circle(Circle scale_me, double scale_factor) {
    	this.centre = new Vector2d(scale_me.centre.getComponents());
    	this.r2 = scale_me.r2 * scale_factor * scale_factor;
    }
    
    /**
     * Constructor creates the smallest {@link Circle} with given two points lying
     * on it's circumference.
     * @param v1
     * 	The first {@link Vector2d}.
     * @param v2
     * 	The second {@link Vector2d}.
     */
    public Circle(Vector2d v1, Vector2d v2) {
        // Centre lies half way along line joining points
        Vector2d a = v2.minus(v1).mult(0.5);
        r2     = a.norm2();
        centre = v1.add(a);
    }
    
    /**
     * Constructs the unique {@link Circle} whose circumference
     * contains the three {@link Vector2d}s. If these are collinear
     * then there is no solution.
     * 
     * @param v1
     * 	The first {@link Vector2d}.
     * @param v2
     * 	The second {@link Vector2d}.
     * @param v3
     * 	The third {@link Vector2d}.
     * 
     */
    public Circle(Vector2d v1, Vector2d v2, Vector2d v3) {
    	
    	Circle c = fitCircle(v1, v2, v3);
    	
    	if(c==null) {
    		throw new RuntimeException("Cannot fit a Circle to points "+v1+", "+v2+", "+v3);
    	}
    	
    	centre = c.centre;
    	r2 = c.r2;
    }
    
    /**
     * Get a {@link Circle} that encloses all of the given points. This is non-optimal
     * in the sense that the enclosing {@link Circle} is not the smallest possible.
     * Finding the smallest enclosing {@link Circle} is a much harder challenge. Here,
     * we simply find the centroid of the point set and the maximum distance to
     * any point.
     * 
     * @param points
     * 	The {@link List} of {@link Vector2d} to enclose.
     * @return
     * 	A (non-optimal) enclosing {@link Circle} for the given set of points.
     */
    public static <T extends Vector2d> Circle getEnclosingCircle(List<T> points) {
    	
        // Get centroid of points
        Vector2d centre = new Vector2d(0,0);
        for(Vector2d point : points) {
        	centre.addEquals(point);
        }
        centre.multEquals(1.0/points.size());
        
        // Now find maximum distance of any point from centroid
        double max_rad_2 = 0;
        for(Vector2d point : points) {
        	max_rad_2 = Math.max(max_rad_2, point.minus(centre).norm2());
        }
        
        return new Circle(centre.getX(), centre.getY(), max_rad_2);
    }
    
    /**
     * Fit a {@link Circle} to the given set of {@link Vector2d} by minimizing, in the
     * least squares sense, the algebraic distance between the points and the {@link Circle}.
     * The algebraic distance is computed from the conic representation of the circle,
     * which is:
     * 
     * 		a(x*x + y*y) + dx + ey + f = 0 
     * 
     * where the parameters a,d,e,f represent the circle.
     * 
     * i.e. the algebraic distance of the point (x_i, y_i) is:
     * 
     * 		a(x_i*x_i + y_i*y_i) + dx_i + ey_i + f
     * 
     * Note: minimizing the geometric distance is not possible with least squares.
     * 
     * 
     * @return boolean
     * 	Specifies whether inliers set could be used to calculate {@link Circle}
     * parameters (true) or not (false).
     * @param points
     * @return
     */
    public static Circle fitCircle(Vector2d... points) {

        // Build design matrix
        double[][] x = new double[points.length][4];

        int N=0;

        for(Vector2d xy : points){
            x[N][0] = xy.getX()*xy.getX() + xy.getY()*xy.getY();
            x[N][1] = xy.getX();
            x[N][2] = xy.getY();
            x[N][3] = 1.0;
            N++;
        }

        // Convert to Matrix. Vector of algebraic distances of each point from 
        // circle is given by X*(a,d,e,f)^T. The solution is found by squaring
        // this and setting derivative equal to zero, i.e.
        //
        // X^T * X * (a,d,e,f)^T = 0
        //
        Matrix X = new Matrix(x);

        // Make design matrix X^T * X
        Matrix D = X.transpose().times(X);
        
        // Get SVD
        SingularValueDecomposition svd = D.svd();
                
        // Get right singular vector corresponding to lowest singular value
        Matrix rsv = svd.getV().getMatrix(new int[]{0,1,2,3}, new int[]{3});

        // Check that second lowest singular value is not zero. This would
        // indicate that no solution has been found, i.e. multiple points lie
        // at the same coordinates.
        if(svd.getS().get(2, 2)<epsilon) {
            // Critical configuration of points detected.
            return null;
        }

        // Set new Circle parameters
        double a = rsv.get(0,0);
        double d = rsv.get(1,0);
        double e = rsv.get(2,0);
        double f = rsv.get(3,0);

        double x0 = -d/(2*a);
        double y0 = -e/(2*a);
        double r2  = (d*d + e*e)/(4*a*a) - f/a;
        
        return new Circle(x0, y0, r2);
    }
    
    /**
     * Get the corners of the (equilateral) {@link Triangle2d} that tightly encloses this
     * {@link Circle}, in the sense that the {@link Circle} lies fully inside
     * the {@link Triangle2d} (e.g. the vertices of the {@link Triangle2d} don't lie on the
     * circumference of the circle). The base lies parallel to the X axis.
     * Note that there is no single solution as we can rotate the {@link Triangle2d}.
     * 
     * @return
     * 	Array of three {@link Vector2d} lying at the corners of a regular {@link Triangle2d}
     * that exactly fits this {@link Circle}.
     */
    public Vector2d[] getEnclosingTriangleVertices() {
    	
    	// Centre of circle
    	double x0 = centre.getX();
    	double y0 = centre.getY();
    	double r = Math.sqrt(r2);
    	
    	// Side lengths of right angle triangle with base parallel to the X axis
    	// and one corner at the circle centre.
    	double x = Math.sqrt(3) * r;
    	double h = 2 * r;
    	
    	// Two vertices lying along the X axis
    	Vector2d v1 = new Vector2d(x0 - x, y0 - r);
    	Vector2d v2 = new Vector2d(x0 + x, y0 - r);
    	// Vertex at the peak of the triangle
    	Vector2d v3 = new Vector2d(x0, y0 + h);
    	
        return new Vector2d[]{v1, v2, v3};
    }
    
    /**
     * Determines if the given {@link Vertex2d} lies within this {@link Circle}.
     * @param vm
     * 	The {@link Vector2d} to check.
     * @return
     * 	True if the given {@link Vertex2d} lies within this {@link Circle}.
     */
    public boolean contains(Vector2d vm) {
        // Check if this is smaller than the radius of the Circle
        return (vm.minus(centre)).norm2()<r2;
    }
    
    @Override
    public String toString() {
        return "(x,y) = "+centre.toString()+" r = "+Math.sqrt(r2);
    }
    
}