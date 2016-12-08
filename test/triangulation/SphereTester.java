package triangulation;

import java.io.*;

import triangulation.Sphere;
import triangulation.Vertex;

/**
 * Unit testing of Sphere class.
 * 
 * @author nickrowell
 */
public class SphereTester {

    public static void main(String[] args) throws IOException{
        
        // Test Sphere fitting
        Vertex v1 = new Vertex( 2, 0, 0);
        Vertex v2 = new Vertex(-2, 0, 0);
        Vertex v3 = new Vertex( 0, 2, 0);
        Vertex v4 = new Vertex( 0,-1, 0);
        
        Sphere sphere = Sphere.fitSphere(v1, v2, v3, v4);
        
        System.out.println(sphere.toString());
        
        System.exit(0);
    }
}
