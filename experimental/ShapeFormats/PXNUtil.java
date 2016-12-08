package ShapeFormats;



import infra.gui.StreamGobbler;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;

import shape.Vertex;




/**
 * Utility class for PXN file writing operations.
 * 
 * @author nickrowell
 */
public class PXNUtil {
    
    public static void writePXNHeader(BufferedWriter out) throws IOException
    {
        out.write("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
        out.write("<pangu_model ver=\"0\">\n");
        out.flush();
    }   
    
    public static void writePXNVertices(BufferedWriter out, 
            Vertex[] verts) throws IOException
    {
        out.write("\t\t<positions>\n");
        // Write vertices out to position node
        for(Vertex vert : verts)
            out.write("\t\t\t<vec3> "+vert.toString()+" </vec3>\n");
        out.write("\t\t</positions>\n");
        out.flush();
    }
     
    public static void writePXNNormals(BufferedWriter out, 
            Vertex[] verts) throws IOException
    {
        out.write("\t\t<normals>\n");
        // Write vertices out to position node
        for(Vertex vert : verts)
        {
            out.write("\t\t\t<vec3> "+vert.getNormalString()+" </vec3>\n");
        }
        out.write("\t\t</normals>\n");          
        
        // Set up ready to write triangular strips
        out.write("\t\t<tristrips>\n");       
        out.flush(); 
        
    }    
    
    
    public static void writeTriangularStrip(BufferedWriter out, int[] tristrip) throws IOException
    {
        out.write("\t\t\t<tristrip count=\"" + tristrip.length + "\">\t");
        for (int v = 0; v < tristrip.length; v++) {
            // Subtract one because PANGU vertex indices are zero-based.
            out.write((tristrip[v]) + " ");
        }
        out.write("\t</tristrip>\n");
        out.flush();
    }    
    
    public static void writePXNFooter(BufferedWriter out) throws IOException
    {
        out.write("</pangu_model>\n");
        out.flush();
    }
    
    

    public static void writePanFile(File pan_file, File pxn_file) throws IOException
    {
        
        String[] commands = new String[]{"pxn2pan","-lic","/opt/pangu/license/nrowell.txt","-o",pan_file.toString(), pxn_file.toString()};
        
        for(String command : commands)
            System.out.print(command+" ");
        
        // Convert PXN file to PAN format
        Process proc = Runtime.getRuntime().exec(commands);
        
        // any error message?
        StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");
        // any output?
        StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

        // kick them off
        errorGobbler.start();
        outputGobbler.start();      
        
        try
        {
            proc.waitFor();
        }
        catch(InterruptedException ie)
        {
            System.err.println("Exception thrown while running pxn2pan");
        }
        
    }
    
    
}
