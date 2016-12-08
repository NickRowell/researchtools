package numeric.fitting;

import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import numeric.geom.dim3.Vector3d;



/**
 * Planes are fitted to data according to two representations of a plane:
 * 
 * 1) N.r - d = 0
 *    ...where N is the surface normal and d is the orthogonal distance from plane to origin.
 *    In this case, the plane is fitted so as to minimise the sum-of-squared orthogonal distance
 *    to the plane. This representation is suited to finding the principal plane of a set of points
 *    in 3D, as it can fit any orientation of plane including vertical.
 * 
 * 2) z = a*x + b*y + c
 *    In this case, the plane is fitted so as to minimise the sum-of-square vertical distance
 *    to the plane. This representation is suited to fitting a plane to a set of data measured over
 *    a regular grid. It cannot fit vertical planes and is not suited to non-gridded data.
 * 
 * @author nrowell
 *
 */
public class PlaneFitting3D {
	
	/**
	 * Uses same plane fitting technique as main fitPlane method, but uses 
	 * outlier rejection via RANSAC to robustify fit.
	 * @param points		List of Vector3d to fit plane to.
	 * @param N				On exit, contain parameters of best fitting plane.
	 * @param clip			Rejection threshold [distance] for RANSAC outlier rejection.
	 * @param REALIZATIONS	Number of RANSAC iterations to perform.
	 * @param debug			Flag indicating whether to print debug info to console.
	 * @return
	 */
	public static double fitPrincipalPlaneRANSAC(List<Vector3d> points, double[] N, double clip, int REALIZATIONS, boolean debug)
	{
		// RANSAC approach is to select a large number of minimal subsets of data, and
		// for each set fit a unique plane and see how many of the remaining data
		// are consistent with it. The plane with the largest 'consensus' among the
		// random sets tried becomes the solution.
		
		// We then fit a plane to the remaining inliers via least squares, and
		// this gives us our solution.
		// We therefore need to keep track of which solution has found the largest
		// consensus thus far.
		
		double[] best_N = new double[4];
		double best_consensus;
		
		double[] current_N = new double[4];
		double current_consensus;
		
		List<Vector3d> inliers  = new LinkedList<>();
		List<Vector3d> outliers = new LinkedList<>();
		
		// Initialise RANSAC solution with first three points
		fitPlane(points.get(0), points.get(1), points.get(2), best_N, false);
		best_consensus = getOrthogonalInliers(points, inliers, outliers, best_N, clip);
		
		// Enter main RANSAC loop
		for(int i=0; i<REALIZATIONS; i++)
		{
			// Shuffle list
			Collections.shuffle(points);
			
			fitPlane(points.get(0), points.get(1), points.get(2), current_N, false);
			current_consensus = getOrthogonalInliers(points, inliers, outliers, current_N, clip);
			
			if(current_consensus > best_consensus)
			{
				// Found new highest consensus
				best_consensus = current_consensus;
				System.arraycopy(current_N, 0, best_N, 0, 4);
			}
		}
		
		// Find set out inliers corresponding to solution with highest consensus
		best_consensus = getOrthogonalInliers(points, inliers, outliers, best_N, clip);
		
		if(debug)
		{
			System.out.println(String.format("RANSAC: consensus set contains %d/%d", inliers.size(), points.size()));
		}
		
		// Fit a plane to these points by least squares
		return fitPrincipalPlane(inliers, N, false);
		
	}
	
	/**
	 * Divides a List of Vector3d into two lists according to whether each Vector3d is
	 * an inlier or outlier with respect to the plane defined by parameters N and
	 * clipping threshold clip. Outliers are identified according to their orthogonal
	 * distance to the plane.
	 * 
	 * @param points	Set of all points.
	 * @param inliers	On exit, contains inliers.
	 * @param outliers	On exit, contains outliers.
	 * @param N			Parameters of plane to test points against.
	 * @param clip		Points at a greater distance from the plane than this are identified as outliers.
	 * @return			The number of inliers.
	 */
	private static int getOrthogonalInliers(List<Vector3d> points, List<Vector3d> inliers, List<Vector3d> outliers,
								   double[] N, double clip)
	{
		inliers.clear();
		outliers.clear();
		
		for(Vector3d p : points)
        {
        	// Orthogonal distance from point to plane:
        	double d = N[0]*p.getX() + N[1]*p.getY() + N[2]*p.getZ() - N[3];
        	
        	if(Math.abs(d) > clip)
        	{
        		// Found outlier
        		outliers.add(p);
        	}
        	else
        	{
        		// Found inlier
        		inliers.add(p);
        	}
        }
		return inliers.size();
	}
	
	/**
	 * Uses same plane fitting technique as main fitPlane method, but uses iterative
	 * outlier rejection via sigma clipping to robustify fit.
	 * @param points	List of Vector3d to fit plane to.
	 * @param N			On exit, contain parameters of best fitting plane.
	 * @param clip		Sigma-clipping threshold: points lying at distances greater than clip*rms from plane are rejected as outliers.
	 * @param debug		Flag indicating whether to print debug info to console.
	 * @return
	 */
	public static double fitPrincipalPlaneSigmaClipping(List<Vector3d> points, double[] N, double clip, boolean debug)
	{
		// Fit plane to all points and get RMS deviation
		double rms = fitPrincipalPlane(points, N, debug);
		
		// Find inliers and outliers
		List<Vector3d> inliers  = new LinkedList<>();
		List<Vector3d> outliers = new LinkedList<>();
		// Distance threshold separating inliers & outliers
		double d = rms*clip;
		// Get inliers/outliers and number of outliers
		int n_outlier = points.size() - getOrthogonalInliers(points, inliers, outliers, N, d);
		
        if(debug)
        {
        	System.out.println(String.format("Clipping %d points...", n_outlier));
        }
        
        if(n_outlier==0)
        {
        	// There were no outliers: the current plane solution is robust.
        	return rms;
        }
        else
        {
        	// Check for pathological case where so many outliers were removed that we've
        	// ended up with too few points to fit a plane. May need to use more robust
        	// method such as RANSAC in these cases.
        	if(points.size()<3)
        	{
        		System.out.println("fitRobustPlane: clipped too many points!");
        		return rms;
        	}
        	
        	return fitPrincipalPlaneSigmaClipping(points, N, clip, debug);
        }
		
		
	}
	
	
	/**
	 * Fits a plane to a set of points in three dimensions. Places no constraints on the
	 * arrangement of the points (e.g. grid is not necessary) and will fit any orientation
	 * of plane including vertical. It represents the plane in terms of the normal vector
	 * N=(nx, ny, nz) and minimum distance to origin d, where points r lying in the plane
	 * satisfy:
	 * 
	 * 		N.r - d = 0
	 * 
	 * For points lying out of the plane, N.r-d gives the orthogonal distance to the plane.
	 * 
	 * The algorithm minimises the sum of squared orthogonal distance to the plane among
	 * the set of points.
	 * 
	 * @param points	List of double arrays each containing (x,y,z) coordinates of one point.
	 * @param N			4-element array that on exit contains the plane parameters (nx,ny,nz,d).
	 * @param debug		Flag indicating whether to print debug info to console.
	 * @return			RMS deviation from the plane for all points.
	 */
	public static double fitPrincipalPlane(List<Vector3d> points, double[] N, boolean debug)
	{
		// First things first, perform some sanity checks:
		assert(N.length==4) : "Require a 4-element array to return plane parameters!";
		
		int M=points.size();
		
		assert(M>2) : "Need at least three points to fit a plane!";
		
		if(M==3)
		{
			System.out.println("Warning: only 3 points provided to fitPlane, unique plane will be returned.");
			fitPlane(points.get(0), points.get(1), points.get(2), N, debug);
			// All three points lie perfectly in the plane: RMS distance from plane is zero.
			return 0.0;
		}
		
		// Compute centroid of points
		double x0=0,y0=0,z0=0;
		
		for(Vector3d point : points)
		{
			x0 += point.getX();
			y0 += point.getY();
			z0 += point.getZ();
			
			if(debug)
			{
				System.out.println("Point "+(points.indexOf(point)+1)+"/"+M+": ("+point.getX()+","+point.getY()+","+point.getZ()+")");
			}
			
		}
		
		x0 /= (double)M;
		y0 /= (double)M;
		z0 /= (double)M;
		
		if(debug)
		{
			System.out.println("Centroid: ("+x0+","+y0+","+z0+")");
		}
		
		// Compute the covariance matrix of points
		double[][] C = new double[][]{{0,0,0},{0,0,0},{0,0,0}};
		
		double xi,yi,zi;	// Centroid-relative coordinates
		
		for(Vector3d point : points)
		{
			xi = point.getX() - x0;
			yi = point.getY() - y0;
			zi = point.getZ() - z0;
			
			C[0][0]+=xi*xi;  C[0][1]+=xi*yi;  C[0][2]+=xi*zi;
			C[1][0]+=yi*xi;  C[1][1]+=yi*yi;  C[1][2]+=yi*zi;
			C[2][0]+=zi*xi;  C[2][1]+=zi*yi;  C[2][2]+=zi*zi;
		}
		
		Matrix c = new Matrix(C);
		
		if(debug)
		{
			System.out.println("Covariance matrix:");
			c.print(5, 5);
		}
		
		// Take eigenvalue decomposition of c
		EigenvalueDecomposition evd = new EigenvalueDecomposition(c);
		
		// Find the eigenvector corresponding to the smallest eigenvalue.
		int smallest_eval_index = 0;
		double smallest_eval = evd.getD().get(0, 0);
		
		for(int i=1; i<3; i++)
		{
			if(evd.getD().get(i,i)<smallest_eval)
			{
				smallest_eval       = evd.getD().get(i,i);
				smallest_eval_index = i;
			}
		}
		
		// Extract the corresponding eigenvector from column i:
		double nx = evd.getV().get(0, smallest_eval_index);
		double ny = evd.getV().get(1, smallest_eval_index);
		double nz = evd.getV().get(2, smallest_eval_index);
		
		// Compute the plane parameter d.
		
		// A point r lying precisely in the plane satisfies the plane equation
		//
		// N.r-d = 0
		//
		// so we can use any point known to lie in the plane to solve for d. Note that
		// likely none of the input points satisfy this; however, we know that the centre of mass
		// must lie in the plane so we can use this:
		double d = nx*x0 + ny*y0 + nz*z0;
		
		N[0] = nx;
		N[1] = ny;
		N[2] = nz;
		N[3] = d;
		
		// Return value is the RMS distance from the plane for all points. Note that each
		// eigenvalue is equal to the sum-of-squared-distances for all points to the plane
		// defined by the corresponding eigenvector. So, the smallest eigenvalue (corresponding
		// to the least-squares solution) actually contains the sum-of-square distances
		// for points to the plane solution. In order to get the RMS, we just divide by the
		// number of points and square root.
		
		if(debug)
		{
			System.out.println(String.format("Plane parameters = %f %f %f %f", N[0], N[1], N[2], N[3]));
			System.out.println(String.format("RMS deviation    = %f", Math.sqrt(smallest_eval/(double)M)));
		}
		
		return Math.sqrt(smallest_eval/(double)M);
		
	}
	
	
	
	/**
	 * This method fits a plane to the set of points by minimising the sum of squared
	 * vertical distances between the points and the plane. It uses the following parameterisation
	 * for the plane when performing the fit:
	 * 
	 * f(x,y) = a*x + b*y + c
	 * 
	 * The fit is performed by solving the following system of linear equations for a,b,c:
	 * 
	 * 	|x_1	y_1	  1| * |a| = |z_1|
	 * 	|x_2	y_2	  1|   |b|   |z_2|
	 * 	|x_3	y_3	  1|   |c|   |z_3|
	 * 	|   ...        |         |...|
	 * 	|x_N	y_N	  1|         |z_N|
	 * 
	 * This method CANNOT SOLVE vertical planes.
	 * 
	 * @param points	
	 * @param N			Plane parameters a,b,c
	 * @param debug		
	 * @return
	 */
	public static double fitPlane(List<Vector3d> points, double[] N, boolean debug)
	{
		// First things first, perform some sanity checks:
		assert(N.length==3) : "Require a 3-element array to return plane parameters for vertical residuals case!";
		
		int M=points.size();
		
		assert(M>2) : "Need at least three points to fit a plane!";
		
		// Build the design and observation matrices
		double[][] Xarr = new double[M][3];
		double[][] Zarr = new double[M][1];
		
		for(int i=0; i<M; i++)
		{
			Vector3d p = points.get(i);
			Xarr[i][0] = p.getX();
			Xarr[i][1] = p.getY();
			Xarr[i][2] = 1.0;
			Zarr[i][0] = p.getZ();
		}
		
		Matrix X = new Matrix(Xarr);
		Matrix Z = new Matrix(Zarr);
		
		// Solve the least squares problem
		Matrix a = X.solve(Z);
		
		// Read plane parameters from solution vector
		N[0] = a.get(0, 0);
		N[1] = a.get(1, 0);
		N[2] = a.get(2, 0);
		
		// Get the residuals vector at the solution
		Matrix e = Z.minus(X.times(a));
		
		// Compute sum-of-square residuals (e is a column vector)
		double SSQ = e.transpose().times(e).get(0, 0);
		
		// Compute RMS
		double RMS = Math.sqrt(SSQ/M);
		
		return RMS;
		
	}
	
	
	public static void transformAbcToNd(double[] abc, double[] N)
	{
		double norm = Math.sqrt(abc[0]*abc[0] + abc[1]*abc[1] + 1);
		
		N[0] =  abc[0]/norm;
		N[1] =  abc[1]/norm;
		N[2] = -1.0/norm;
		N[3] = -abc[2]/norm;
	}
	
	public static void transformNdToAbc(double[] abc, double[] N)
	{
		abc[0] = -N[0]/N[2];
		abc[1] = -N[1]/N[2];
		abc[2] =  N[3]/N[2];
	}
	
	/**
	 * Finds the unique plane passing through 3 points r0, r1, r2.
	 * @param r0
	 * @param r1
	 * @param r2
	 * @param N
	 * @param debug
	 */
	public static void fitPlane(Vector3d r0, Vector3d r1, Vector3d r2, double[] N, boolean debug)
	{
		// Get the surface normal
		Vector3d normal = Vector3d.getClockwiseSurfaceNormal(r0, r1, r2);
		
		// Get the minimum distance to the origin. Can use any of r0, r1, r2 here.
		double d = normal.dot(r0);
		
		N[0] = normal.getX();
		N[1] = normal.getY();
		N[2] = normal.getZ();
		N[3] = d;
		
	}
	
	
}
