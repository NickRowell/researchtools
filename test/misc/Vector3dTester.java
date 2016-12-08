package misc;

import misc.Vector3d;

/**
 *
 * @author nickrowell
 */
public class Vector3dTester 
{
    
    public static void main(String[] args)
    {
        
        testRandomVectorWithinCone();
        
    }
    
    
    private static void testRandomVector()
    {
        
        // Number of random vectors to draw
        int N = 1000;
        
        for(int n=0; n<N; n++)
        {
        
            Vector3d rand = Vector3d.getRandVecOnUnitSphere();
            
            System.out.println(rand.x+"\t"+rand.y+"\t"+rand.z);
        
        }
        
    }
    
    private static void testRandomVectorWithinCone()
    {
        
        // Number of random vectors to draw
        int N = 1000;
        
        Vector3d in = new Vector3d(1,1,1);
        double opening = Math.toRadians(10);
                
        for(int n=0; n<N; n++)
        {
        
            Vector3d rand = Vector3d.getRandVecOnUnitSphereWithinCone(in, opening);
            
            System.out.println(rand.x+"\t"+rand.y+"\t"+rand.z);
        
        }
        
    }
    
    private static void testCrossProduct()
    {
        Vector3d z = new Vector3d(1,0,0);
        
        Vector3d r = new Vector3d(0,1,0);
        
        r = r.normalise();
        
        Vector3d cross = z.cross(r);
        
        System.out.println("cross = "+cross.toString());
    }
    
    
    
}
