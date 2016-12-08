package shape.pxn;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * This class provides a program that generates a 3D arrow model in PXN format
 * for inclusion in PANGU models.
 * 
 * @author nickrowell
 */
public class Generate3DArrow
{
    
    /** Output PXN file path. */
    public File pxn_file = new File("/opt/pangu/Models/Artificial/Tests/arrow/bluearrow.pxn");
    
    /** Arrow stem length. */
    public double stem_length = 200;
    /** Stem radius. */
    public double stem_radius = 2;
    /** Head length (base to tip). */
    public double head_length = 30;
    /** Head radius (base). */
    public double head_radius = 5;
    
    /** RGB colour. */
    public double r = 0.0;
    public double g = 0.0;
    public double b = 1.0;
    
    public Generate3DArrow() throws IOException
    {
        
        // Open writer on output file
        BufferedWriter out = new BufferedWriter(new FileWriter(pxn_file));
        
        // Write PXN header
        out.write("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
        out.write("<pangu_model ver=\"0\">\n");      
        out.write("<group count=\"2\">\n");
        
        
        writeStem(out);
        
        writeHead(out);
        
        
        // Write PXN footer
        out.write("</group>\n");
        out.write("</pangu_model>\n");
        out.flush();
        out.close();
        
    }
    
    private void writeStem(BufferedWriter out) throws IOException
    {
        
        // Cease recursion and write out mesh for this node
        out.write("<mesh>\n");
        
        // Number of points around ring at each end of stem
        int N = 8;
        double ang_step = 2 * Math.PI / N;
        
        // Write out vertex positions
        out.write("<positions>\n");
        for(int i=0; i<N; i++)
        {
            double angle = i * ang_step;
            double x = stem_radius * Math.cos(angle);
            double y = stem_radius * Math.sin(angle);
            double z = -stem_length;
            out.write("<vec3> "+x+" "+y+" "+z+" </vec3>\n");
            z = stem_length;
            out.write("<vec3> "+x+" "+y+" "+z+" </vec3>\n");
        }
        out.write("</positions>\n");
        
        // Write out surface normals
        out.write("<normals>\n");
        for(int i=0; i<N; i++)
        {
            double angle = i * ang_step;
            // Outward facing normals
            double x = Math.cos(angle);
            double y = Math.sin(angle);
            out.write("<vec3> "+x+" "+y+" 0 </vec3>\n");
            out.write("<vec3> "+x+" "+y+" 0 </vec3>\n");
        }
        out.write("</normals>\n");

        
        // Write out colours
        out.write("<colours>\n");
        for(int i=0; i<N; i++)
        {
            out.write("<rgb> "+r+" "+g+" "+b+" </rgb>\n");
            out.write("<rgb> "+r+" "+g+" "+b+" </rgb>\n");
        }
        out.write("</colours>\n");
        
        // Write out triangular strip
        out.write("<tristrips>\n");
        out.write("<tristrip>\t");
        for(int i=0; i<N; i++)
        {
            out.write( (2*i+0) + " " + (2*i+1) + " ");
        }
        // Connect back to first two vertices
        out.write("0 1 ");
        out.write("</tristrip>\n");
        out.write("</tristrips>\n");
        out.write("</mesh>\n");
    }
    
    
    private void writeHead(BufferedWriter out) throws IOException
    {
        
        // Cease recursion and write out mesh for this node
        out.write("<mesh>\n");
        
        // Number of points around ring at base of head
        int N = 8;
        double ang_step = 2 * Math.PI / N;
        
        // Write out vertex positions
        out.write("<positions>\n");
        // Single point at tip of head
        out.write("<vec3> 0 0 "+(stem_length + head_length)+" </vec3>\n");
        // Points around base of head
        for(int i=0; i<N; i++)
        {
            double angle = i * ang_step;
            double x = head_radius * Math.cos(angle);
            double y = head_radius * Math.sin(angle);
            double z = stem_length;
            out.write("<vec3> "+x+" "+y+" "+z+" </vec3>\n");
        }
        out.write("</positions>\n");
        
        // Write out surface normals
        out.write("<normals>\n");
        // Single point at tip of head
        out.write("<vec3> 0 0 1 </vec3>\n");
        // Points around base of head
        for(int i=0; i<N; i++)
        {
            double angle = i * ang_step;
            // Outward facing normals
            double x = Math.cos(angle);
            double y = Math.sin(angle);
            out.write("<vec3> "+x+" "+y+" 0 </vec3>\n");
        }
        out.write("</normals>\n");

        
        // Write out colours
        out.write("<colours>\n");
        for(int i=0; i<N+1; i++)
        {
            out.write("<rgb> "+r+" "+g+" "+b+" </rgb>\n");
        }
        out.write("</colours>\n");
        
        // Write out triangular strips
        out.write("<tristrips>\n");
        
        // One strip for each face on side of head
        for(int i=1; i<N; i++)
        {
            out.write("<tristrip>\t");
            out.write(" 0 " + i + " " + (i+1));
            out.write("</tristrip>\n");
        }
        // Final strip joining back to start
        out.write("<tristrip>\t");
        out.write(" 0 " + N + " 1");
        out.write("</tristrip>\n");
        
        // One strip covering base of head
        out.write("<tristrip>\t");
        for(int i=0; i<N; i++)
        {
            out.write((i+1) + " " + (N-i) + " ");
        }
        out.write("</tristrip>\n");
        
        
        out.write("</tristrips>\n");
        out.write("</mesh>\n");
    }
    
            
    
    public static void main(String[] args) throws IOException
    {
        
        new Generate3DArrow();
        
    }
    
    
    
    
}
