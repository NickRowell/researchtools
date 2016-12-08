/**
 * Copyright (c) 2011 University of Dundee
 *
 * Name:
 * Frame.java
 *
 * Purpose:
 * This class extends the Pose class to represent the projective frame of
 * a pinhole camera. It complements the position and attitude parameters of
 * the Pose class (the 'extrinsic' camera parameters) with a perspective
 * projection matrix K (the 'intrinsic' camera parameters).
 * 
 * It also has a List of correspondences between 3D points and their 2D images.
 * This is useful in determining the 6 degree of freedom Pose of a camera,
 * which is done in the PoseDetermination class.
 * 
 * Language:
 * Java
 *
 * Author:
 * Nicholas Rowell (nickrowell@computing.dundee.ac.uk)
 *
 */
package numeric.vision;

import Jama.Matrix;

import java.util.ArrayList;
import java.util.List;

import numeric.dynamics.Pose;
import numeric.geom.dim3.Quaternion;
import numeric.geom.dim3.Vector3d;

/**
 * This 
 *
 * @author nickrowell
 */
public class Frame extends Pose {

    /** Camera calibration matrix */
    Matrix K;
    
    /** Map between Point2D objects and Point3D objects */
    List<Correspondence_2D_3D> correspondences = new ArrayList<Correspondence_2D_3D>();

    /** Basic constructor */
    public Frame(){}
    
    /** Constructor. */
    public Frame(Matrix k, Matrix P){
        // Construct Pose from camera matrix
        super(P);
        K = k;
    }
    
    /** Constructor. */
    public Frame(Matrix k, Pose pose){
        super(pose);
        K = k;
    }    

    public void setK(Matrix k){ K = k;}
    
    /** Project a Point3D onto the image plane */
    public Point2D project(Point3D p3d){
        
        // Transform 3D position from world to camera frame
        Matrix X_cam = super.getConjP().times(p3d.x);
        
        // Transform covariance matrix to camera frame
        Matrix X_cov_cam = super.toBodyFrame(p3d.cov);
        
        // Project point onto image plane
        Matrix x = K.times(X_cam);
        // Homogenise...
        x.timesEquals(1.0/x.get(2,0));
        
        // Project covariance matrix onto image plane. This is done
        // by first order propagation of uncertainty, i.e.
        //
        // cov_image_plane = J * cov_camera_frame * J^T
        //
        // where J is the Jacobian for the transformation from camera frame
        // coordinates to image plane coordinates.
        
        // Get some numbers
        double rx = X_cam.get(0, 0);
        double ry = X_cam.get(1, 0);
        double rz = X_cam.get(2, 0);
        double f  = K.get(0, 0); 
        
        // Jacobian
        Matrix J = new Matrix(new double[][]{{f/rz,   0, -f*rx/(rz*rz)},
                                             {0,   f/rz, -f*ry/(rz*rz)}});
        
        // Complete the transformation
        Matrix x_cov = J.times(X_cov_cam).times(J.transpose());
        
        // Create & return new Point2D
        return new Point2D(x,x_cov);
    
    }
    
    
    /** Add a Correspondence_2D_2D_3D to List */
    public void addCorrespondence_2D_3D(Correspondence_2D_3D corr){
        correspondences.add(corr);
    }
    

    /**
     * Get geometric parameters for this frame, formatted like this:
     * 
     * element  quantity
     * 0,0      Real part of quaternion
     * 1,0      i imaginary part of quaternion
     * 2,0      j   "        "         "
     * 3,0      k   "        "         "
     * 4,0      x translation vector component
     * 5,0      y   "        "         "  
     * 6,0      z   "        "         "
     * 
     */
    public Matrix getParameters(){
                  
        // Construct a Matrix to hold parameter values
        Matrix params = new Matrix(7,1);
        
        // Write quaternion components to first 4 elements
        params.set(0, 0, attitude.re);
        params.set(1, 0, attitude.im.getX());
        params.set(2, 0, attitude.im.getY());
        params.set(3, 0, attitude.im.getZ());
        // Write translation components to next 3 elements
        params.set(4, 0, position.getX());
        params.set(5, 0, position.getY());
        params.set(6, 0, position.getZ());

        return params;  
    }
    
    /** 
     * New parameters are organised in the matrix delta in the following
     * manner:
     * 
     * 0-3: components of Quaternion
     * 4-6: components of translation vector
     */
    public void setParameters(Matrix delta){

        // First update components...
        attitude.setRe(delta.get(0, 0));
        attitude.setQ1(delta.get(1, 0));
        attitude.setQ2(delta.get(2, 0));
        attitude.setQ3(delta.get(3, 0));
        // Enforce unit quaternion...
        attitude.normalise();
        
        // Elements 4-6 are the new translation vector
        position.setX(delta.get(4, 0));
        position.setX(delta.get(5, 0));
        position.setX(delta.get(6, 0));
        
        return;
    }   
    
    /**
     * Get the angle of rotation required to align this Frame with the rotation
     * matrix given as argument.
     * 
     */
    public double getRotationDiffMag(Frame f2){
                
        // Quaternion required to align this with current Frame
        Quaternion Q_error = getRotationDiff(f2);    
        
        // Magnitude of angular displacement
        double rotation_error = Math.toDegrees(2.0*Math.acos(Q_error.re));
        // Shift to 0:+PI range
        if(rotation_error > 180.0) rotation_error = 360-rotation_error;
    
        return rotation_error;
        
    }
    
    public Quaternion getRotationDiff(Frame f2){
        // Quaternion required to align this with current Frame
        return f2.attitude.multiply(attitude.inverse());
    }
    
    /**
     * Get the translation vector in the world frame between two frames.
     */
    public Vector3d getTranslationDiff(Frame f2){
        
        // Difference gives translational offset in world frame
        return position.minus(f2.position);
    
    }
    
    
    
    
    
}
