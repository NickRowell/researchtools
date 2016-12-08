/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package shape.vrml;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.StreamTokenizer;
import java.util.List;

import numeric.geom.dim3.Quaternion;
import numeric.geom.dim3.Vector3d;

/**
 *
 * @author nickrowell
 */
public class VrmlTransformNode extends AbstractGroupingNode
{
    // fields - default values
    Vector3d scale       = new Vector3d(1f,1f,1f);
    Quaternion rotation  = new Quaternion(1,0,0,0);
    Vector3d translation = new Vector3d(0f,0f,0f);
    
    // Fields with no counterpart in PANGU
    Quaternion scaleOrientation = new Quaternion(1,0,0,0);
    Vector3d center             = new Vector3d(0f,0f,0f);
    
    public VrmlTransformNode(String _name)
    {
        super(_name);
    }
    
    // Detect identity transforms, i.e. transform nodes that don't apply any 
    // transformation operations to their children. For identity nodes with
    // only one child, we can simply skip the node and write the child in
    // place of the parent.
    public boolean isIdentity()
    {
        return scale.equals(new Vector3d(1,1,1)) &&
               rotation.equals(Quaternion.getIdentity()) &&
               translation.equals(new Vector3d(0,0,0)) &&
               scaleOrientation.equals(Quaternion.getIdentity()) &&
               center.equals(new Vector3d(0,0,0));
    }
    
    int writeNodeHeader(BufferedWriter out, int indent) throws IOException
    {
        
        if(!scaleOrientation.equals(Quaternion.getIdentity()))
            throw new RuntimeException("ScaleOrientation field defined for vrml"
                    + " Transform node: "+name+". This field is not supported "
                    + "by PANGU.");
        if(!center.equals(new Vector3d(0,0,0)))
            throw new RuntimeException("center field defined for vrml Transform"
                    + " node: "+name+". This field is not supported by PANGU.");
        
        String indentation = "";
        for(int i=0; i<indent; i++) 
            indentation += indentChar;
        
        
        out.write(indentation+"<transform>\n");
        out.write(indentation+indentChar+"<scaling>" + scale.getX() + " " + scale.getY() + " " + scale.getZ() + "</scaling>\n");
        
        // We use a <rotate> node to perform rotation
        out.write(indentation+indentChar+"<rotation> 0 0 0 </rotation>\n");
        
        out.write(indentation+indentChar+"<translation>" + translation.getX() + " " + translation.getY() + " " + translation.getZ() + "</translation>\n");
        out.write(indentation+indentChar+"<child>\n");
        
        // Use rotate node instead
        out.write(indentation+indentChar+indentChar+"<rotate>\n");
        out.write(indentation+indentChar+indentChar+"<quat> "+rotation.re+" "+rotation.im.getX()+" "+rotation.im.getY()+" "+rotation.im.getZ()+" </quat>\n");
        out.write(indentation+indentChar+indentChar+"<child>\n");
        
        // Avoid additional <group> element for single children
        if(children.size() > 1)
        {
            out.write(indentation+indentChar+indentChar+indentChar+"<group>\n");
            out.flush();
            return 4;
        }
        out.flush();
        return 3;
    }
    
    void writeNodeFooter(BufferedWriter out, int indent) throws IOException
    {
        String indentation = "";
        for(int i=0; i<indent; i++) 
            indentation += indentChar;
        
        // Avoid additional <group> element for single children
        if(children.size() > 1) out.write(indentation+indentChar+indentChar+indentChar+"</group>\n");
        
        // Use <rotate> node
        out.write(indentation+indentChar+indentChar+"</child>\n");
        out.write(indentation+indentChar+indentChar+"</rotate>\n");
        
        out.write(indentation+indentChar+"</child>\n");
        out.write(indentation+"</transform>\n");
        
        
        
        out.flush();
    }
    
    
    // Stream tokenizer is positioned on opening curly brace on entry.
    // Stream tokenizer is positioned just after closing curly brace on exit.
    public static VrmlTransformNode parse(StreamTokenizer st, String name, List<VrmlNode> DEFs) throws IOException
    {
        VrmlTransformNode transform = new VrmlTransformNode(name);
        
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
                
                if(field.equalsIgnoreCase("translation"))
                {
                    st.nextToken();
                    transform.translation = parseVector3d(st);
                }
                else if(field.equalsIgnoreCase("rotation"))
                {
                    st.nextToken();
                    transform.rotation = parseQuaternion(st);
                }
                else if(field.equalsIgnoreCase("scale"))
                {
                    st.nextToken();
                    transform.scale = parseVector3d(st);
                }
                else if(field.equalsIgnoreCase("scaleOrientation"))
                {
                    st.nextToken();
                    transform.scaleOrientation = parseQuaternion(st);
                }
                else if(field.equalsIgnoreCase("center"))
                {
                    st.nextToken();
                    transform.center = parseVector3d(st);
                }
                else if(field.equalsIgnoreCase("children"))
                {
                    // Move to next token [
                    st.nextToken();
                    parseChildren(st, transform.children, DEFs);
                    
                }
                else
                {
                    throw new RuntimeException("Unrecognized Transform field: "+field);
                }
                
            }
            
        }
        
        // Position stream tokenizer just after closing curly brace
        st.nextToken();
        
        return transform;
    }
}
