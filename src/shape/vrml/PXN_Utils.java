/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package shape.vrml;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

import util.RelativePath;

/**
 *
 * @author nickrowell
 */
public class PXN_Utils
{
    
    // Default material
    static VrmlMaterialNode material = new VrmlMaterialNode("default");
    
        
    public static void writePXNHeader(BufferedWriter out) throws IOException
    {
        out.write("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
        out.write("<pangu_model ver=\"0\">\n");
        out.flush();
    }
    
    public static void writePXNFooter(BufferedWriter out) throws IOException
    {
        out.write("</pangu_model>\n");
        out.flush();
    }
    public static void writeInclude(BufferedWriter out, File top, File pxn, int indent, String indentChar, boolean addAppearance) throws IOException
    {
        String indentation = "";
        for(int i=0; i<indent; i++) 
            indentation += indentChar;
        
        if(addAppearance)
        {
            // Wrap in default material
            out.write(indentation+"<appearance>\n");
            material.writePXN(out, indent+1);
            out.write(indentation+indentChar+"<child>\n");
            out.write(indentation+indentChar+indentChar+"<include type=\"pxn\" file=\""+RelativePath.getRelativePath(top,pxn)+"\"/>\n");

            // Wrap in default material
            out.write(indentation+indentChar+"</child>\n");
            out.write(indentation+"</appearance>\n");
        }
        else
        {
            out.write(indentation+"<include type=\"pxn\" file=\""+RelativePath.getRelativePath(top,pxn)+"\"/>\n");
        }
        
        
        out.flush();
    }
    
    
}
