/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package shape.vrml;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.StreamTokenizer;
import java.util.LinkedList;
import java.util.List;

/**
 *
 * @author nickrowell
 */
public class VrmlCoordinateNode extends VrmlNode
{
    List<Triple<Float>> points = new LinkedList<Triple<Float>>();
    
    
    
    public VrmlCoordinateNode(String _name)
    {
        super(_name);
    }
    
    
    void writePXN(BufferedWriter out, int indent) throws IOException
    {
        String indentation = "";
        for(int i=0; i<indent; i++) 
            indentation += indentChar;
        
        out.write(indentation+"<positions>\n");
        for(Triple<Float> point : points)
            out.write(indentation+"\t<vec3> "+point.x+" "+point.y+" "+point.z+" </vec3>\n");
        out.write(indentation+"</positions>\n");
        
    }
    
    
    // Stream tokenizer is positioned on opening curly brace on entry.
    // Stream tokenizer is positioned just after closing curly brace on exit.
    static VrmlCoordinateNode parse(StreamTokenizer st, String name, List<VrmlNode> DEFs) throws IOException
    {
        VrmlCoordinateNode coordinate = new VrmlCoordinateNode(name);
        
        // move tokenizer to next token: either first field, or closing brace
        // for empty objects.
        st.nextToken();
        
        // Loop over tokens until closing curly brace is reached
        while(!foundClosingCurlyBrace(st))
        {
            // Read field:
            if(st.ttype == StreamTokenizer.TT_WORD)
            {
                String field = st.sval;
                
                if(field.equalsIgnoreCase("point"))
                {
                    // Move to next token [
                    st.nextToken();
                    parseTriples(st, coordinate.points);
                }
                else
                    throw new RuntimeException("Unrecognized Coordinate field: "+field);
                
            }
            
        }
        
        // Position stream tokenizer just after closing curly brace
        st.nextToken();
        
        return coordinate;
        
    }
    
    
    
}
