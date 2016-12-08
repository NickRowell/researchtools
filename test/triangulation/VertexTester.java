package triangulation;

import java.util.List;
import java.util.ArrayList;

import triangulation.Vertex;

/**
 *
 * @author nickrowell
 */
public class VertexTester {
    
    static List<Vertex> verts = new ArrayList<Vertex>();
    
    public static void main(String[] args){
    
        verts.add(new Vertex( 1, 1, 1));
        verts.add(new Vertex( 2, 2, 2));
        verts.add(new Vertex( 3, 3, 3));
        verts.add(new Vertex( 4, 4, 4));       
    
        
        for (Vertex v : verts)
            // Is v a reference to the Vertex in the List or a completely
            // new object? Try incrementing the x coordinates if objects in the
            // List.
            v.x+=1.0;
        
        for (Vertex v : verts)
            // Has x coordinate been increased? If so, v is a reference, if not
            // then v is a new object.
            System.out.println(v.toString());       
    
    }
    
    
}
