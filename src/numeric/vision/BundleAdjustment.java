package numeric.vision;


import java.util.ArrayList;
import java.util.List;

import Jama.*;
import numeric.minimisation.nllsq.algo.LevenbergMarquardt;

/**
 * Class designed for sketching out BundleAdjustment algorithm architecture.
 *
 * To do: Add printStats() method that prints out statistics of fit, i.e. 
 *        number of frames, number of point correspondences in each, number
 *        of unique scene points, number constrained by two or more frames.
 * 
 * 
 * 
 * @author nickrowell
 */
public class BundleAdjustment extends LevenbergMarquardt{

    /** Recovered motion of camera in terms of Frame objects */
    public List<Frame> frames;
    /** Recovered structure of scene in terms of Point3D objects */
    public List<Point3D> scene;

    /** Verbose behaviour */
    private boolean verbose = false;
    
    /** Basic constructor */
    public BundleAdjustment(boolean v){
        verbose = v;
        frames = new ArrayList<Frame>();
        scene = new ArrayList<Point3D>();
    }

    
    /** Add a Frame and bunch of Point3D objects to BA */
    public void addFrame(Frame input){
    
        if(verbose)
            System.out.println("BA: Adding frame "+frames.size()+" with "
                               +input.correspondences.size()+" scene points");
        
        // Add Frame to List
        frames.add(input);

        // Count number of new unique scene points
        int new_points=0, old_points=0;
        
        // Frames maintain maps between the Point2D objects in their image
        // space and the Point3D objects in the scene space. Because multiple
        // Frames may observe the same Point3D, must only add Point3D objects
        // to List if they are not already contained in it.
        for(Correspondence_2D_3D corr: input.correspondences){
            if(!scene.contains(corr.p3d)){
                scene.add(corr.p3d);
                new_points++;
            }
            else{ old_points++;}
        }
        
        if(verbose)
            System.out.println("BA: Added "+new_points+" new points, reused "
                               +old_points+" existing points, total points = "
                                +scene.size());
        
    }

    public int getFramesN(){return frames.size();}
    public int getSceneN(){return scene.size();}
    
    /** Required by Levenberg-Marquardt algorithm */
    public int getParametersN(){ return 7*frames.size() + 3*scene.size();}
    
    /** Required by Levenberg-Marquardt algorithm */
    public Matrix getParametersSteps(){
           
        Matrix steps = new Matrix(getParametersN(),1);
        
        // Finite steps for frame geometry parameters
        Matrix frame_steps = new Matrix(new double[][]{{0.001},
                                                       {0.001},
                                                       {0.001},
                                                       {0.001},
                                                       {1.0},
                                                       {1.0},
                                                       {1.0}});
        
        // Finite steps for scene point parameters (position coordinates)
        Matrix scene_steps = new Matrix(new double[][]{{1.0},
                                                       {1.0},
                                                       {1.0}});       
        
        // Number of frames
        int N = frames.size();
        
        // Enter geometry parameters for each frame
        for (int f=0; f<N; f++) 
            // Enter 7 frame geometry parameters into matrix
            steps.setMatrix(7*f, 7*f+6, 0, 0, frame_steps);
        
        // Enter positions for each scene point
        for (int p=0; p<scene.size(); p++)
            steps.setMatrix(7*N + 3*p, 7*N + 3*p + 2, 0, 0, scene_steps);
        
        return steps; 
    }
    
    
    /**
     * BundleAdjustment parameters are arranged as follows:
     * 
     * 0-6:       Geometry of frame 0 (inclusive)
     * 7-13:      Geometry of frame 1
     * :
     * :
     * 7N-(7N+6): Geometry of frame N   (there are N+1 frames)
     * 
     * 7(N+1) + 0:      X component of scene point 0
     * 7(N+1) + 1:      Y component of scene point 0
     * 7(N+1) + 2:      Z component of scene point 0
     * :
     * :
     * 7(N+1) + 3*P + 0:  X component of scene point P
     * 7(N+1) + 3*P + 1:  Y component of scene point P
     * 7(N+1) + 3*P + 2:  Z component of scene point P
     * 
     * @return 
     */
    public Matrix getParameters(){
    
        Matrix params = new Matrix(getParametersN(),1);
        
        // Number of frames
        int N = frames.size();
        
        // Enter geometry parameters for each frame
        for (int f=0; f<N; f++) 
            // Enter 7 frame geometry parameters into matrix
            params.setMatrix(7*f, 7*f+6, 0, 0, frames.get(f).getParameters());
        
        // Enter positions for each scene point
        for (int p=0; p<scene.size(); p++)
            params.setMatrix(7*N + 3*p, 7*N + 3*p + 2, 0, 0, scene.get(p).getParameters());
        
        return params;
    
    }
    public void setParameters(Matrix delta){
    
        // Number of frames
        int N = frames.size();
        
        // Set geometry parameters for each frame
        for (int f=0; f<frames.size(); f++) 
            // Enter 7 frame geometry parameters into matrix
            frames.get(f).setParameters(delta.getMatrix(7*f, 7*f+6, 0, 0));
            
        
        // Set positions for each scene point
        for (int p=0; p<scene.size(); p++)
            scene.get(p).setParameters(delta.getMatrix(7*N + 3*p, 7*N + 3*p + 2, 0, 0));
        
        return;
    
    }    
    public void printParameters(){
        System.out.println("BundleAdjustment.printParameters() called.");
    }

    
    /** Required by Levenberg-Marquardt algorithm. */
    public int getDataN(){
    
        // Sum Point2D objects across all Frames
        int sum_points = 0;
        
        for (Frame frame : frames) 
            sum_points += frame.correspondences.size();
        
        // Two data for each Point2D
        return sum_points*2;
    
    }
    
    /** Required by Levenberg-Marquardt algorithm. */
    public Matrix getData(){
   
        // Generate data vector
        Matrix data = new Matrix(getDataN(),1);

        // Track index of next data to enter
        int index=0;
        
        // Loop over each Frame
        for (Frame frame : frames) {
            // Loop over each point in this frame
            for (Correspondence_2D_3D p2d3d : frame.correspondences){
                // Set i coordinate
                data.set(index++, 0, p2d3d.p2d.x.get(0, 0));
                // Set j coordinate
                data.set(index++, 0, p2d3d.p2d.x.get(1, 0));                
            }
        }
    
        return data;
    
    }
    
    /** Required by Levenberg-Marquardt algorithm. */
    public Matrix getModel(){
   
        // Generate model vector
        Matrix model = new Matrix(getDataN(),1);

        // Track index of next data to enter
        int index=0;
        
        // Loop over each Frame
        for (Frame frame : frames) {
            
            // Loop over each point in this frame
            for (Correspondence_2D_3D p2d3d : frame.correspondences){
                
                // Set subsidiary (model) image coordinates by projecting
                // Point3D onto image plane              
                
                p2d3d.setSubsPoint(frame.K, frame.getP());
                // Set model i coordinate
                model.set(index++, 0, p2d3d.p2d.x_sub.get(0, 0));
                // Set model j coordinate
                model.set(index++, 0, p2d3d.p2d.x_sub.get(1, 0));                
            }
        }
    
        return model;
    
    }    
    
    
    
    /** Required by Levenberg-Marquardt algorithm. */
    public Matrix getCovariance(){
   
        // Generate covariance matrix for all data
        Matrix cov = new Matrix(getDataN(), getDataN());

        // Track index of next point to enter
        int index=0;
        
        // Loop over each Frame
        for (Frame frame : frames) {
            // Loop over each point in this frame
            for (Correspondence_2D_3D p2d3d : frame.correspondences){
                
                // Enter covariance matrix for this point into diagonal
                // block of covariance matrix.
                cov.setMatrix(new int[]{2*index, 2*index+1}, 
                              new int[]{2*index, 2*index+1}, 
                              p2d3d.p2d.cov);
                
                // Increment point counter
                index++;
                              
            }
        }
    
        return cov;
        
    }
    
    /** 
     * Required by LMA. Used to test sensitivity of parameter solution to 
     * changes in the data.
     */
    public boolean updateData(Matrix delta){
        
        // Track index of next point to update
        int index=0;
        
        // Loop over each Frame
        for (Frame frame : frames) {
            // Loop over each Point2D in this Frame
            for (Correspondence_2D_3D p2d3d : frame.correspondences){
                
                // Make update vector (need to add extra zero element because
                // Point2D coordinates are homogenous).
                Matrix update = new Matrix(3,1);
                update.setMatrix(0,1,0,0, delta.getMatrix(2*index, 2*index+1, 0, 0));
                
                p2d3d.p2d.x.plusEquals(update);
                
                // Increment point counter
                index++;              
               
            }
        }
        return true;
    }       
    
    
    
    
    

}