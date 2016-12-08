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
public abstract class AbstractGroupingNode extends VrmlNode
{
    List<VrmlNode> children = new LinkedList<VrmlNode>();
    
    public AbstractGroupingNode(String _name)
    {
        super(_name);
    }
    
    // Returns number of extra indentations to use for sub-nodes
    abstract int writeNodeHeader(BufferedWriter out, int indent) throws IOException;
    abstract void writeNodeFooter(BufferedWriter out, int indent) throws IOException;
    // identifies groups that apply no operations to their children
    abstract boolean isIdentity();
    
    
    static void parseChildren(StreamTokenizer st, List<VrmlNode> children, List<VrmlNode> DEFs) throws IOException
    {
        // Move to start of first child, or closing ] for empty children field
        st.nextToken();
        
        while(!VrmlNode.foundClosingSquareBrace(st))
        {
            children.add(parseNode(st, DEFs));
        }
        // Position stream tokenizer just after closing square brace
        st.nextToken();
        
    }
}
