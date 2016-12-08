package gl.java.util;

import Jama.Matrix;
import numeric.geom.dim3.Vector3d;


/**
 * This class represents the projection model of a pinhole camera.
 * 
 * @author nickrowell
 */
public class PerspectiveCamera 
{
    
    /** Camera projection matrix K. */
    public Matrix K;
    
    /** Inverse of camera matrix. Very common operation so do it once. */
    public Matrix invK;
    
    /** Extent of image in focal plane (inclusive) [pixels] */
    public int X0;
    public int XMAX;
    public int Y0;
    public int YMAX;
    
    /** Field of view [degrees]. */
    public double fovx;
    public double fovy;
    
    /** Pixel aspect ratio [dimensionless] */
    public double aspect;
    
    /**  Focal length [pixels]. */
    public double fx;
    public double fy;
    
    /** Viewport size. */
    public int width;
    public int height;
    
    /**
     * Constructor from field of view [degrees] and image size [pixels].
     * Assumes that the principal point is in the centre of the image.
     * 
     * @param fovx
     * 	Field of view in the X (horizontal) direction [degrees]
     * @param fovy
     * 	Field of view in the Y (vertical) direction [degrees]
     * @param width
     * 	Width of the image (horizontal/X direction) [pixels]
     * @param height
     * 	Width of the image (vertical/Y direction) [pixels]
     */
    public PerspectiveCamera(double fovx, double fovy, int width, int height)
    {
        // Field of view in degrees
        this.fovx = fovx;
        this.fovy = fovy;
        
        this.width  = width;
        this.height = height;
        
        // Focal length in pixels
        this.fx = ((double) (this.width / 2))  * (1.0 / Math.tan(Math.toRadians(this.fovx) / 2.0));
        this.fy = ((double) (this.height / 2)) * (1.0 / Math.tan(Math.toRadians(this.fovy) / 2.0));
        
        this.aspect = (width*fy)/(height*fx);
        
        // Boundaries of image in projection plane
        this.X0 = 0;
        this.XMAX = this.width-1;
        this.Y0 = 0;
        this.YMAX = this.height-1;
        
        // Image coordinates of projection centre (assumed in image centre)
        double px = this.width/2;
        double py = this.height/2;
        
        this.K = new Matrix(new double[][]{{fx,0,px},{0,fy,py},{0,0,1}});
        this.invK = this.K.inverse();
    }
    
    /**
     * Converts the camera matrix from the pinhole format to the GL format; this requires the
     * near and far plane positions.
     * @param n
     * 	The near plane distance.
     * @param f
     * 	The far plane distance.
     * @return
     * 	The camera projection matrix expressed in the GL format.
     */
    public float[] getGlProjectionMatrix(double n, double f) {

		// Extract focal lengths and principal point from camera matrix
		float fi = (float)K.get(0,0);
		float fj = (float)K.get(1,1);
		float pi = (float)K.get(0,2);
		float pj = (float)K.get(1,2);

		// Compute clip plane positions
		float l = (float)-n * pi/fi;
		float r =  (float)n * (width-pi)/fi;
		// OpenGL projection has image coordinate origin at bottom left. Must
		// flip Y here in order to produce images with the correct pinhole 
		// camera geometry.
		float t =  (float)n * (height-pj)/fj;
		float b = (float)-n * pj/fj;
		
		// Transposed version (OpenGL wants column-major format)
		float[] projection = {(float)(2*n/(r-l)),          0,            0,            0,
			                  	       0,  (float)(2*n/(t-b)),            0,            0,
			                 (r+l)/(r-l),(t+b)/(t-b), (float)(-(f+n)/(f-n)),           -1,
			                           0,          0, (float)(-2*f*n/(f-n)),            0};
		
		return projection;
    }
    
    /** Get the image width. */
    public int getImageWidth()
    {
        return XMAX - X0 + 1;
    }
    
    /** Get the image height. */
    public int getImageHeight()
    {
        return YMAX - Y0 + 1;
    }
    
    /**
     * Get the unit vector towards a given pixel coordinate.
     * @param i
     * 	Pixel coordinate in the X (horizontal) direction [pixels]
     * @param j
     * 	Pixel coordinate in the Y (vertical) direction [pixels]
     * @return
     * 	Unit vector in the camera frame pointing in the direction defined by the pixel
     */
    public Vector3d getUnitVector(double i, double j)
    {
        // Deproject the pixel
    	Matrix a = invK.times(new Matrix(new double[][]{{i},{j},{1}}));
    	
    	// Vector of length of 1 along optical axis (Z)
        Vector3d unitZ = new Vector3d(a.get(0,0), a.get(1,0), a.get(2,0));
        
        // Unit vector
        return unitZ.normalise();
    }
    
    /**
     * Project the given camera frame position vector into the image plane.
     * @param X_vec
     * 	Camera frame position vector.
     * @return
     * 	Array containing the coordinates of the point in the image plane.
     */
    public float[] projectVector(Vector3d X_vec)
    {
        // Convert to Matrix object for projection
        Matrix X = new Matrix(new double[][]{{X_vec.x},{X_vec.y},{X_vec.z}});
        // Projection
        Matrix x = K.times(X);
        // Perspective division
        x.timesEquals(1.0/x.get(2,0));
        
        return new float[]{(float)x.get(0,0), (float)x.get(1,0)};
    }
    
    /**
     * Project a covariance matrix for a position in the camera frame, 
     * expressed in the camera frame, into the image plane. Image plane
     * covariance matrix is returned in row-major format.
     * 
     * @param X_CAM
     * 	Mean position in camera frame
     * @param cov_CAM
     * 	The covariance matrix in camera frame 3D coordinates
     * @return
     * 	The covariance matrix for the image plane coordinates
     */
    public Matrix projectMatrix(Vector3d X_CAM, Matrix cov_CAM)
    {

        // Now project the covariance matrix to the image plane. This is done
        // by first order propagation of uncertainty, i.e.
        //
        // cov_image_plane = J * cov_camera_frame * J^T
        //
        // where J is the Jacobian for the transformation from camera frame
        // coordinates to image plane coordinates.
        
        // Get some numbers
        double x = X_CAM.x;
        double y = X_CAM.y;
        double z = X_CAM.z;
        
        // Jacobian
        Matrix J = new Matrix(new double[][]{{fx/z,   0, -fx*x/(z*z)},
                                             {0,   fy/z, -fy*y/(z*z)}});
        
        // Now calculate (J * cov_camera_frame) * J^T
        Matrix cov_IM = J.times(cov_CAM).times(J.transpose());
        
        // Sanity check on image plane covariance. All eigenvalues should be
        // positive, for a symmetric positive-definite matrix.
        
        double det = cov_IM.get(0,0)*cov_IM.get(1,1) - cov_IM.get(0,1)*cov_IM.get(1,0);
        double tr  = cov_IM.get(0,0)+cov_IM.get(1,1);
        
        // Watch out for complex eigenvalues
        if(tr*tr<4*(det))
        {
            System.err.println("Illegal image plane covariance matrix!"
                    + " Complex eigenvalues:");
            System.err.println(cov_IM.get(0,0)+"\t"+cov_IM.get(1,0));
            System.err.println(cov_IM.get(0,1)+"\t"+cov_IM.get(1,1));
            
            System.err.println("X_CAM: \n"+X_CAM.toString());
            System.err.println("cov_CAM: ");
            cov_CAM.print(5,5);
            System.exit(1);
        }
        
        // Watch out for negative eigenvalues
        if(!(det>=0 && tr>=0))
        {
            System.err.println("Illegal image plane covariance matrix!"
                    + " Negative eigenvalue(s):");
            System.err.println(cov_IM.get(0,0)+"\t"+cov_IM.get(1,0));
            System.err.println(cov_IM.get(0,1)+"\t"+cov_IM.get(1,1));
            System.err.println("X_CAM: \n"+X_CAM.toString());
            System.err.println("cov_CAM: ");
            cov_CAM.print(5,5);
            System.exit(1);
        }     
        
        return cov_IM;
    }    
    
    
}