package ShapeFormats;

import java.util.LinkedList;
import java.util.List;

public class TriStrip 
{
    
    List<Integer> vertex = new LinkedList<Integer>();

    // Current working Facet
    Facet current;
    
    // Preceeding facet: recorded so that we can figure out the correct 
    // order to write out the vertices of the current Facet when the strip
    // is terminated.
    Facet previous;
    
    // Basic constructor
    public TriStrip(){}
    
    
    public TriStrip(Facet _current)
    {
        current = _current;
    }
    
    /** Write remaining vertices from current Facet to vertex list. */
    public void finalise(){
        
        // Note that the order of these vertices is carefully calculated so 
        // that winding order is preserved.
        vertex.add(current.v0);
        vertex.add(current.v1);
        vertex.add(current.v2);
    }
    
    
    /**
     * The given Facet has two vertices in common with the current Facet. 
     * Record the single isolated vertex from the current Facet, then update
     * the current Facet with the vertices of the given Facet. Note that these
     * have to be carefully arranged into the correct order so that if the 
     * strip terminates when the next Facet is checked, we can write out
     * the vertices v1,v2,v3 and know that the winding order is preserved.
     * 
     * We couldn't simply write out the vertices for the
     * current Facet in any order because it would likely not link into the
     * preceeding strip correctly.
     * 
     * @param newFacet 
     */
    public void add(Facet newFacet){
        
        // Add the single isolated vertex from first facet to the
        // vertex list.
        List<Integer> vert = current.getVertsNotInCommon(newFacet);

        // Sanity check...
        if (vert.size() != 1) {
            //System.out.println("verts not in common = "+vert.size());
            //System.out.println("Vertices in current: "+current.v1+", "+current.v2+", "+current.v3);
            //System.out.println("Vertices in newFacet: "+newFacet.v1+", "+newFacet.v2+", "+newFacet.v3);
            throw new RuntimeException("More than one isolated vertex");
        }

        vertex.add(vert.get(0));
        
        if(previous != null)
        {
            // Update current facet: arrange order of vertices so when they are
            // written out in natural order (1,2,3) they have the correct
            // winding order.
            int v1 = newFacet.getVertsInCommon(previous).get(0);
            int v2 = current.getVertsNotInCommon(previous).get(0);
            int v3 = newFacet.getVertsNotInCommon(current).get(0);
            
            previous = new Facet(current);
            current = new Facet(v1,v2,v3);
                
        }
        
        // The given facet is only the second in the strip. We don't need to be
        // as careful with the vertex order in this case.
        else
        {
            previous = new Facet(current);
            // Order of vertices 1 & 2 doesn't matter in this case.
            current = new Facet(newFacet.getVertsInCommon(current), newFacet.getVertsNotInCommon(current));
        }
        
    }
    
}
