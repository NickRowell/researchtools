package misc;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import dynamics.Pose;
import misc.Quaternion;
import misc.Vector3d;

/**
 *
 * @author nickrowell
 */
public class VidForIrene
{
    
    // Main body - positions & orientations of other bodies referred to this frame
    Pose main      = new Pose();
    
    // Satellite - initial position & orientation
    double sat_r = 1575;  // Satellite orbital radius
    Pose satellite = new Pose(new Vector3d(0,sat_r,0), new Quaternion(1,0,0,0));
    
    // Camera frame
    double cam_az = Math.toRadians(60);  // angle between camera position and +y axis
    double cam_r  = 8000;                // camera distance from origin
    Pose camera = new Pose(new Vector3d(Math.sin(cam_az)*cam_r, Math.cos(cam_az)*cam_r, 0), new Quaternion(1,0,0,0));
    
    // Sun direction in cartesian coordinates in main body frame
    Vector3d sun = new Vector3d(1000000000,0,0);
    
    
    public VidForIrene() throws IOException
    {
        // Turn camera to point boresight towards origin
        camera.pointVectorAt(new Vector3d(0,0,1), new Vector3d(0,0,0));
        
        // Rotation applied to main body each frame
        Quaternion rot_main = new Quaternion(new Vector3d(0,1,1), Math.toRadians(-0.125));
        
        // Angular displacement of satellite through orbit each frame
        double sat_disp = Math.toRadians(0.25);
        
        File fli = new File("/opt/pangu/movies/for_irene/flight.fli");
        BufferedWriter out = new BufferedWriter(new FileWriter(fli));
        
        
        // Iteration loop
        for(int i=0; i<2880; i++)
        {
            // Update model positions etc.
            main.rotate(rot_main);    // Rotate main body
            
            satellite.setPosition(new Vector3d( 0, 
                                                sat_r*Math.cos(i*sat_disp),
                                               -sat_r*Math.sin(i*sat_disp)));
            
            // Satellite tidally locked
            satellite.pointVectorAt(new Vector3d(0,0,1), new Vector3d(0,0,0));
        
            // Camera moves through ten degrees over course of flight
            cam_az = Math.toRadians(60-((i/2880.0)*10));
            camera.setPosition(new Vector3d(Math.sin(cam_az)*cam_r, Math.cos(cam_az)*cam_r, 0));
            camera.pointVectorAt(new Vector3d(0,0,1), new Vector3d(0,0,0));
            
            // Write flight file entry
            out.write("dynamic_object 3 "+main.getPosition().x+" "+
                                          main.getPosition().y+" "+
                                          main.getPosition().z+" "+
                                          main.getAttitude().re+" "+
                                          main.getAttitude().im.x+" "+
                                          main.getAttitude().im.y+" "+
                                          main.getAttitude().im.z+"\n");
            
            out.write("dynamic_object 4 "+satellite.getPosition().x+" "+
                                          satellite.getPosition().y+" "+
                                          satellite.getPosition().z+" "+
                                          satellite.getAttitude().re+" "+
                                          satellite.getAttitude().im.x+" "+
                                          satellite.getAttitude().im.y+" "+
                                          satellite.getAttitude().im.z+"\n");
            
            out.write("sun_position 10000000000 120 0\n");
            
            out.write("quaternion "+camera.getPosition().x+" "+
                                    camera.getPosition().y+" "+
                                    camera.getPosition().z+" "+
                                    camera.getAttitude().re+" "+
                                    camera.getAttitude().im.x+" "+
                                    camera.getAttitude().im.y+" "+
                                    camera.getAttitude().im.z+"\n");
            
        }
        
        out.flush();
        out.close();
        
        
        
    }
    
    public static void main(String[] args) throws IOException
    {
        new VidForIrene();
    }
    
}
