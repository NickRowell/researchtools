package triangulation;

import java.util.List;
import java.util.ArrayList;
import java.io.*;

import triangulation.DelaunayTriangulation3D;

/**
 *
 * @author nickrowell
 */
public class DelaunayTester {
    
    public static void main(String[] args) throws IOException{
        
        // Generate a bunch of Vertex objects distributed uniformly within
        // a box of side length 5 units.
        List<Vertex> verts = new ArrayList<Vertex>();
        
        for(int i=0; i<20; i++)
            verts.add(new Vertex(5*(Math.random()-0.5), 
                                 5*(Math.random()-0.5), 
                                 5*(Math.random()-0.5)));
        
        // Triangulate these points
        DelaunayTriangulation3D del = new DelaunayIncremental(verts);
    
        // Write out delaunay triangulation vertices and triangles as a 
        // ply shape model.
//        del.writePLY(new File("/home/nickrowell/Java_projects/java_AbsNav/"
//                            + "development/point_cloud_to_mesh/"
//                            + "delaunay/model_external.ply"), 1);
        
//        del.writePLY(new File("/home/nickrowell/Java_projects/java_AbsNav/"
//                            + "development/point_cloud_to_mesh/"
//                            + "delaunay/model_internal.ply"), 2);        
        
        
    
    }
    
    
}
