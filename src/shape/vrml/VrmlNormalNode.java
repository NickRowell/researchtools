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
public class VrmlNormalNode extends VrmlNode
{
    List<Triple<Float>> normals = new LinkedList<Triple<Float>>();
    
    
    
    public VrmlNormalNode(String _name)
    {
        super(_name);
    }
    
    
    void writePXN(BufferedWriter out, int indent) throws IOException
    {
        String indentation = "";
        for(int i=0; i<indent; i++)
            indentation += indentChar;
        
        out.write(indentation+"<normals>\n");
        for(Triple<Float> normal : normals)
            out.write(indentation+"\t<vec3> "+normal.x+" "+normal.y+" "+normal.z+" </vec3>\n");
        out.write(indentation+"</normals>\n");
        
    }
    
    
    // Stream tokenizer is positioned on opening curly brace on entry.
    // Stream tokenizer is positioned just after closing curly brace on exit.
    static VrmlNormalNode parse(StreamTokenizer st, String name, List<VrmlNode> DEFs) throws IOException
    {
        VrmlNormalNode normal = new VrmlNormalNode(name);
        
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
                
                if(field.equalsIgnoreCase("vector"))
                {
                    // Move to next token [
                    st.nextToken();
                    parseTriples(st, normal.normals);
                }
                else
                    throw new RuntimeException("Unrecognized Coordinate field: "+field);
                
            }
            
        }
        
        // Position stream tokenizer just after closing curly brace
        st.nextToken();
        
        return normal;
        
    }
    
    
    
}
