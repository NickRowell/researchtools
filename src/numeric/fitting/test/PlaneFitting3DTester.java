package numeric.fitting.test;

import java.util.LinkedList;
import java.util.List;

import numeric.fitting.PlaneFitting3D;
import numeric.geom.dim3.Vector3d;

public class PlaneFitting3DTester
{
	
	
	public static void main(String[] args)
	{
		
//		testPlaneMinimiseVerticalResiduals();
//		testPlaneMinimiseOrthogonalResiduals();
		testPlaneParameterTransform();
		
	}
	
	
	
	
	private static void testPlaneParameterTransform()
	{
		
		// Create a set of points lying in a plane parallel to the XY plane by
		// offset in Z direction.
		List<Vector3d> r = new LinkedList<>();
		
		for(int i=0; i<5; i++)
		{
			r.add(new Vector3d(Math.random(), Math.random(), Math.random()));
		}
		
		double[] abc = new double[3];
		double rms = PlaneFitting3D.fitPlane(r, abc, false);
		
		
		System.out.println("PlaneFitting3DTester");
		System.out.println("++++++++++++++++++++\n");
		System.out.println("Parameters of fitted plane:");
		System.out.println(String.format("(a,b,c) = (%f, %f, %f)",abc[0], abc[1], abc[2]));
		
		double[] N = new double[4];
		
		PlaneFitting3D.transformAbcToNd(abc, N);
		
		System.out.println(String.format("N = (%f, %f, %f)",N[0], N[1], N[2]));
		System.out.println(String.format("D = %f",N[3]));
		
		PlaneFitting3D.transformNdToAbc(abc, N);
		
		System.out.println(String.format("(a,b,c) = (%f, %f, %f)",abc[0], abc[1], abc[2]));
		
	}
	
	
	private static void testPlaneMinimiseVerticalResiduals()
	{
		// Create a set of points lying in a plane parallel to the XY plane by
		// offset in Z direction.
		List<Vector3d> r = new LinkedList<>();
		
		for(int i=0; i<100000; i++)
		{
			r.add(new Vector3d(Math.random(), Math.random(), -5.32+0.1*Math.random()));
		}
		
		// Add one point at very large distance
		r.add(new Vector3d(0, 0, 1000000000));
		
		double[] N = new double[3];
		double rms = PlaneFitting3D.fitPlane(r, N, false);
		
		System.out.println("PlaneFitting3DTester");
		System.out.println("++++++++++++++++++++\n");
		System.out.println("Parameters of fitted plane:");
		System.out.println(String.format("(a,b,c) = (%f, %f, %f)",N[0], N[1], N[2]));
		System.out.println(String.format("RMS = %f",rms));
		
	}
	
	
	
	
	private static void testPlaneMinimiseOrthogonalResiduals()
	{
		// Create a set of points lying in a plane parallel to the XY plane by
		// offset in Z direction.
		List<Vector3d> r = new LinkedList<>();
		
		for(int i=0; i<100000; i++)
		{
			r.add(new Vector3d(Math.random(), Math.random(), -5.32+0.1*Math.random()));
		}
		
		// Add one point at very large distance
		r.add(new Vector3d(0, 0, 1000000000));
		
		double[] N = new double[4];
		double rms = PlaneFitting3D.fitPrincipalPlane(r, N, false);
		
		System.out.println("PlaneFitting3DTester");
		System.out.println("++++++++++++++++++++\n");
		System.out.println("Parameters of fitted plane:");
		System.out.println(String.format("N = (%f, %f, %f)",N[0], N[1], N[2]));
		System.out.println(String.format("D = %f",N[3]));
		System.out.println(String.format("RMS = %f",rms));
	}
	
	
	
	
	
}
