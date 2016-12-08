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
public class GenerateChessboardPXN
{
    
    /** Output PXN file path. Filename will be derived from settings. */
    public File pxn_file = new File("/opt/pangu/Models/Artificial/Tests/chessboard");
    
    /** Width of board [squares]. */
    public int board_width = 8;
    /** Size of each square [metres]. */
    public double square_width = 2;
    
    /** Colour of white squares. */
    public double[] white = new double[]{1.0, 1.0, 1.0};
    /** Colour of black squares. */
    public double[] black = new double[]{0.0, 0.0, 0.0};
    
    public GenerateChessboardPXN() throws IOException
    {
        
        // Build filename from chessboard parameters
        String name = "chessboard_"+board_width+"x"+board_width+
                      "_"+square_width;
        
        // Open writer on output file
        BufferedWriter out = new BufferedWriter(new FileWriter(new File(pxn_file, name+".pxn")));
        
        // Write PXN header
        out.write("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
        out.write("<pangu_model ver=\"0\">\n");      
        out.write("<group count=\"1\">\n");
        
        
        writeChessboard(out);
        
        // Write PXN footer
        out.write("</group>\n");
        out.write("</pangu_model>\n");
        out.flush();
        out.close();
        
        System.out.println("pxn2pan -o "+name+".pan "+name+".pxn");
        
    }
    
    private void writeChessboard(BufferedWriter out) throws IOException
    {
        
        // Cease recursion and write out mesh for this node
        out.write("<mesh>\n");
        
        
        double wm2 = (board_width*square_width)/2.0; // half total width of board  [metres]
        double hm2 = (board_width*square_width)/2.0; // half total height of board [metres]

        out.write("<positions>\n");
        // Loop over all squares
        for(int i=0; i<board_width; i++)
        {
            for(int j=0; j<board_width; j++)
            {
                // Coordinates in the XY plane of the bottom left corner of this square
                double x = i*square_width - wm2;
                double y = j*square_width - hm2;
                double z = 0.0;                   // Lies in XY plane at z=0

                double s = square_width;   // use concise alias for square_width
                
                out.write("<vec3> "+(x+0)+" "+(y+0)+" "+z+" </vec3>\n");  // bottom left [looking along -z]
                out.write("<vec3> "+(x+s)+" "+(y+0)+" "+z+" </vec3>\n");  // bottom right [looking along -z]
                out.write("<vec3> "+(x+s)+" "+(y+s)+" "+z+" </vec3>\n");  // top right [looking along -z]
                out.write("<vec3> "+(x+0)+" "+(y+s)+" "+z+" </vec3>\n");  // top left [looking along -z]
            }
        }
        out.write("</positions>\n");
        
        // Write out surface normals (all along positive z)
        out.write("<normals>\n");
        for(int i=0; i<board_width; i++)
        {
            for(int j=0; j<board_width; j++)
            {
                out.write("<vec3> 0 0 1 </vec3>\n");
                out.write("<vec3> 0 0 1 </vec3>\n");
                out.write("<vec3> 0 0 1 </vec3>\n");
                out.write("<vec3> 0 0 1 </vec3>\n");
            }
        }
        out.write("</normals>\n");

        
        // Write out colours
        out.write("<colours>\n");
        for(int i=0; i<board_width; i++)
        {
            for(int j=0; j<board_width; j++)
            {
                // Detect square colour
                if((i+j)%2==0)
                {
                    // White square
                    out.write("<rgb> "+white[0]+" "+white[1]+" "+white[2]+" </rgb>\n");
                    out.write("<rgb> "+white[0]+" "+white[1]+" "+white[2]+" </rgb>\n");
                    out.write("<rgb> "+white[0]+" "+white[1]+" "+white[2]+" </rgb>\n");
                    out.write("<rgb> "+white[0]+" "+white[1]+" "+white[2]+" </rgb>\n");
                }
                else
                {
                    // Black square
                    out.write("<rgb> "+black[0]+" "+black[1]+" "+black[2]+" </rgb>\n");
                    out.write("<rgb> "+black[0]+" "+black[1]+" "+black[2]+" </rgb>\n");
                    out.write("<rgb> "+black[0]+" "+black[1]+" "+black[2]+" </rgb>\n");
                    out.write("<rgb> "+black[0]+" "+black[1]+" "+black[2]+" </rgb>\n");
                }
            }
        }
        out.write("</colours>\n");
        
        // Write out triangular strips
        out.write("<tristrips>\n");
        // Each square forms one very short strip consisting of two triangles
        for(int i=0; i<board_width; i++)
        {
            for(int j=0; j<board_width; j++)
            {
                out.write("<tristrip>\t");
                // Index of first vertex in this square. The other 3 corners 
                // are k+1, k+2, k+3
                int k = 4*(i*board_width+j);
                out.write( k + " " + (k+1) + " "+(k+3)+" "+(k+2));
                out.write("</tristrip>\n");
            }
        }
        out.write("</tristrips>\n");
        out.write("</mesh>\n");
    }
    
    public static void main(String[] args) throws IOException
    {
        
        new GenerateChessboardPXN();
        
    }
    
    
    
    
}
