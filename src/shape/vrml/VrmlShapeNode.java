/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package shape.vrml;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.StreamTokenizer;
import java.util.List;

/**
 *
 * @author nickrowell
 */
public class VrmlShapeNode extends VrmlNode
{
    // Fields
    VrmlAppearanceNode appearance;
    AbstractGeometryNode geometry;
    
    public VrmlShapeNode(String _name)
    {
        super(_name);
    }
    
    
    public void writePXN(BufferedWriter out, int indent) throws IOException
    {
        String indentation = "";
        for(int i=0; i<indent; i++) 
            indentation += indentChar;
        
        // Don't use any of the default materials applied by VRML
//        out.write(indentation+"<appearance>\n");
//        appearance.material.writePXN(out, indent+1);
//        out.write(indentation+indentChar+"<child>\n");
        
        geometry.writePXN(out, indent+2);
        
//        out.write(indentation+indentChar+"</child>\n");
//        out.write(indentation+"</appearance>\n");
        
    }
    
    
    
        // Stream tokenizer is positioned on opening curly brace on entry.
    // Stream tokenizer is positioned just after closing curly brace on exit.
    static VrmlShapeNode parse(StreamTokenizer st, String name, List<VrmlNode> DEFs) throws IOException
    {
        VrmlShapeNode shape = new VrmlShapeNode(name);
        
        // move tokenizer to next token: either first field, or closing brace
        // for empty objects.
        st.nextToken();
        
        // Loop over tokens until closing curly brace is reached
        while(!VrmlNode.foundClosingCurlyBrace(st))
        {
            // Read field:
            if(st.ttype == StreamTokenizer.TT_WORD)
            {
                String field = st.sval;
                
                if(field.equalsIgnoreCase("appearance"))
                {   
                    st.nextToken();
                    shape.appearance = (VrmlAppearanceNode)parseNode(st, DEFs);
                }
                else if(field.equalsIgnoreCase("geometry"))
                {
                    st.nextToken();
                    shape.geometry = (VrmlIndexedFaceSetNode)parseNode(st, DEFs);
                }
                else
                    throw new RuntimeException("Unrecognized Shape field: "+field+", "+st.toString());
                
            }
            
        }
        
        // Position stream tokenizer just after closing curly brace
        st.nextToken();
        
        return shape;
    }
    
    
    
    
    
}
