/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package shape.vrml;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.StreamTokenizer;

/**
 *
 * @author nickrowell
 */
public class VrmlMaterialNode extends VrmlNode
{
    Triple diffuseColor     = new Triple<Float>(1.0f, 1.0f, 1.0f);
    Triple specularColor    = new Triple<Float>(1.0f, 1.0f, 1.0f);
    Triple emissiveColor    = new Triple<Float>(0.0f, 0.0f, 0.0f);
    float  shininess        = 0.2f;
    
    // Not supported by PANGU PXN format
    float  ambientIntensity = 0.0f;
    float  transparency     = 0.0f;
    
    public VrmlMaterialNode(String _name)
    {
        super(_name);
    }
    
    public void writePXN(BufferedWriter out, int indent) throws IOException
    {
        
        String indentation = "";
        for(int i=0; i<indent; i++) 
            indentation += indentChar;
        
        out.write(indentation+"<material>\n");
        
        out.write(indentation+indentChar+"<diffuse> " + diffuseColor.x + " " + diffuseColor.y + " " + diffuseColor.z + " </diffuse>\n");
        out.write(indentation+indentChar+"<specular shine=\""+(int)(shininess*128)+"\"> " + specularColor.x + " " + specularColor.y + " " + specularColor.z + " </specular>\n");
        out.write(indentation+indentChar+"<emissive> " + emissiveColor.x + " " + emissiveColor.y + " " + emissiveColor.z + " </emissive>\n");
        
        // Other material properties available in PANGU. Should set these
        // to defaults.
        out.write(indentation+indentChar+"<!-- Other material properties supported by PANGU -->\n");
        out.write(indentation+indentChar+"<!--reflect> 0.0 1.0 0.0 </reflect-->\n");
        out.write(indentation+indentChar+"<!--refract> 4.2 2.77 2.62 </refract-->\n");
        out.write(indentation+indentChar+"<!--oren> 0.1 </oren-->\n");
        out.write(indentation+indentChar+"<!--hapke> 0.33 0.1 0.95 15 </hapke-->\n");
        out.write(indentation+indentChar+"<!--cook> 1.0 0.5 1.0 1.0 </cook-->\n");
        
        out.write(indentation+"</material>\n");
    }
    
    
    // Stream tokenizer is positioned on opening curly brace on entry.
    // Stream tokenizer is positioned just after closing curly brace on exit.
    // Material nodes contain no other nodes as fields, so we don't need to pass
    // in the list of node definitions.
    static VrmlMaterialNode parse(StreamTokenizer st, String name) throws IOException
    {
        VrmlMaterialNode material = new VrmlMaterialNode(name);
        
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
                
                if(field.equalsIgnoreCase("ambientIntensity"))
                {
                    st.nextToken();
                    material.ambientIntensity = VrmlNode.parseFloat(st);
                }
                else if(field.equalsIgnoreCase("diffuseColor"))
                {
                    st.nextToken();
                    material.diffuseColor = VrmlNode.parseTripleFloat(st);
                }
                else if(field.equalsIgnoreCase("specularColor"))
                {
                    st.nextToken();
                    material.specularColor = VrmlNode.parseTripleFloat(st);
                }
                else if(field.equalsIgnoreCase("emissiveColor"))
                {
                    st.nextToken();
                    material.emissiveColor = VrmlNode.parseTripleFloat(st);
                }
                else if(field.equalsIgnoreCase("shininess"))
                {
                    st.nextToken();
                    material.shininess = VrmlNode.parseFloat(st);
                }
                else if(field.equalsIgnoreCase("transparency"))
                {
                    st.nextToken();
                    material.transparency = VrmlNode.parseFloat(st);
                }
                else
                    throw new RuntimeException("Unrecognized Material field: "+field);
                
            }
            
        }
        
        // Position stream tokenizer just after closing curly brace
        st.nextToken();
        
        return material;
    }
    
    
    
    
}
