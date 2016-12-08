package dynamics;


import dynamics.Pose;
import misc.Quaternion;
import misc.Vector3d;


/**
 *
 * @author nickrowell
 */
public class PoseTester {
    
    
    
    /**
     * Test code. 
     */
    public static void main(String[] args){
    
        
        double c45 = Math.cos(Math.toRadians(45));
        double s45 = Math.sin(Math.toRadians(45));
        
        // Rotate 90 degrees about world -y axis
        Quaternion R1 = new Quaternion(c45,0,-s45,0);
        
        // Rotate 90 degrees about world -x axis
        Quaternion R2 = new Quaternion(c45,-s45,0,0);

        // and chain rotations
        Quaternion initialAtt = (R2.multiply(R1));
        
        Pose body = new Pose(new Vector3d(0,1,0),initialAtt);
                                                   
        // Test orientation of basis vectors of each frame
        System.out.println("Body x axis in world frame = "+body.getBodyXAxis().toString());
        System.out.println("Body y axis in world frame = "+body.getBodyYAxis().toString());
        System.out.println("Body z axis in world frame = "+body.getBodyZAxis().toString());
        System.out.println("World x axis in body frame = "+body.getWorldXAxis().toString());
        System.out.println("World y axis in body frame = "+body.getWorldYAxis().toString());
        System.out.println("World z axis in body frame = "+body.getWorldZAxis().toString());
        
        // Check transformations of position vectors between frames
        Vector3d worldFramePoint = new Vector3d(1,1,0);
        System.out.println("Point (1,1,0) in World frame has Body frame position "+body.toBodyFrame(worldFramePoint));
        Vector3d bodyFramePoint = new Vector3d(1,1,0);
        System.out.println("Point (1,1,0) in Body frame has World frame position "+body.toWorldFrame(bodyFramePoint));
        
        // Check aligning of body frame with certain point
        System.out.println("Pointing x axis at origin");
        body.pointVectorAt(new Vector3d(1,0,0), new Vector3d(0,0,0));
        System.out.println("New body x axis in world frame = "+body.getBodyXAxis().toString());
        System.out.println("New body y axis in world frame = "+body.getBodyYAxis().toString());
        System.out.println("New body z axis in world frame = "+body.getBodyZAxis().toString());
        
        
        
        // Test relative transformation between two Pose objects
        
        // Body 1 is placed at world origin with body x axis pointing along world
        // z axis, body z axis pointing along world -x axis, and y axes aligned
        Pose body1 = new Pose(new Vector3d(0,0,0),new Quaternion(c45,0,-s45,0));
        // Body 2 is placed at world origin with body y axis pointing along world
        // z axis, body z axis pointing along world -y axis, and x axes aligned        
        Pose body2 = new Pose(new Vector3d(0,0,0),new Quaternion(c45,s45,0,0)); 
        
        Pose b12 = body1.getRelativePose(body2);
        
        // Results should show that:
        // body 2 x axis points along body 1 -z axis
        // body 2 y axis points along body 1 x axis
        // body 2 z axis points along body 1 -y axis
        // body 1 x axis points along body 2 y axis
        // body 1 y axis points along body 2 -z axis
        // body 1 z axis points along body 2 -x axis        
        
        System.out.println("body2 x axis in body1 frame = "+b12.getBodyXAxis().toString());
        System.out.println("body2 y axis in body1 frame = "+b12.getBodyYAxis().toString());
        System.out.println("body2 z axis in body1 frame = "+b12.getBodyZAxis().toString());      
        System.out.println("body1 x axis in body2 frame = "+b12.getWorldXAxis().toString());
        System.out.println("body1 y axis in body2 frame = "+b12.getWorldYAxis().toString());
        System.out.println("body1 z axis in body2 frame = "+b12.getWorldZAxis().toString());       
        
    }
    
    
}
