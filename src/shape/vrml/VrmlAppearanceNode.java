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
public class VrmlAppearanceNode extends VrmlNode
{
    // Fields
    VrmlMaterialNode material;
    // Texture
    // TextureTransform
    
    public VrmlAppearanceNode(String _name)
    {
        super(_name);
    }
    
    
        // Stream tokenizer is positioned on opening curly brace on entry.
    // Stream tokenizer is positioned just after closing curly brace on exit.
    static VrmlAppearanceNode parse(StreamTokenizer st, String name, List<VrmlNode> DEFs) throws IOException
    {
        VrmlAppearanceNode appearance = new VrmlAppearanceNode(name);
        
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
                
                if(field.equalsIgnoreCase("material"))
                {
                    st.nextToken();
                    appearance.material = (VrmlMaterialNode)parseNode(st, DEFs);
                }
                else if(field.equalsIgnoreCase("texture"))
                    throw new RuntimeException("Unsupported Appearance field: "+field);
                else if(field.equalsIgnoreCase("textureTransform"))
                    throw new RuntimeException("Unsupported Appearance field: "+field);
                else
                    throw new RuntimeException("Unrecognized Appearance field: "+field+
                            ", "+st.toString());
                
            }
            
        }
        
        // Position stream tokenizer just after closing curly brace
        st.nextToken();
        
        return appearance;
    }
    
    // Appearance node can have children in PXN format.
//    public void writePXN(BufferedWriter out, int indent) throws IOException
//    {
//        
//        String indentation = "";
//        for(int i=0; i<indent; i++) 
//            indentation += indentChar;
//        
//        out.write(indentation+"<appearance>\n");
//        material.writePXN(out, indent+1);
//        out.write(indentation+"</appearance>\n");
//        
//        // Texture/TextureTransform nodes not supported
//        
//        
//    }
    
    
}
