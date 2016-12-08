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
public class VrmlGroupNode extends AbstractGroupingNode
{
    public VrmlGroupNode(String _name)
    {
        super(_name);
    }
    
    int writeNodeHeader(BufferedWriter out, int indent) throws IOException
    {
        // Avoid additional <group> element for single children
        if(children.size() > 1)
        {
            for(int i=0; i<indent; i++) out.write(indentChar);
            out.write("<group>\n");
            out.flush();
            return 1;
        }
        out.flush();
        return 0;
    }
    
    void writeNodeFooter(BufferedWriter out, int indent) throws IOException
    {
        // Avoid additional <group> element for single children
        if(children.size() > 1)
        {
            for(int i=0; i<indent; i++) out.write(indentChar);
            out.write("</group>\n");
        }
        out.flush();
    }
    
    // Group nodes always apply no operations to children
    boolean isIdentity()
    {
        return true;
    }
    
    
    
    
        // Stream tokenizer is positioned on opening curly brace on entry.
    // Stream tokenizer is positioned just after closing curly brace on exit.
    static VrmlGroupNode parse(StreamTokenizer st, String name, List<VrmlNode> DEFs) throws IOException
    {
        VrmlGroupNode group = new VrmlGroupNode(name);
        
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
                
                if(field.equalsIgnoreCase("children"))
                {
                    // Move to next token [
                    st.nextToken();
                    parseChildren(st, group.children, DEFs);
                }
                else
                {
                    throw new RuntimeException("Unrecognized Group field: "+field);
                }
                
            }
            
        }
        
        // Position stream tokenizer just after closing curly brace
        st.nextToken();
        
        return group;
    }
    
}
