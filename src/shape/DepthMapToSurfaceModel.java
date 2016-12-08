package shape;

import Jama.Matrix;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;

/**
 * This class is used to transform a depth map (image where each pixel contains
 * a distance) and camera matrix (used to convert depth values to a 3D point in
 * the camera frame) into a triangulated 3D shape model.
 * 
 * @author nickrowell
 */
public class DepthMapToSurfaceModel
{
    
    
    public static void main(String[] args) throws IOException
    {
        // Camera matrix
        Matrix K_pangu = new Matrix(new double[][]{{955.4,  0.0, 256},
                                                   {0.0,  955.4, 256},
                                                   {0.0,    0.0,   1}});
        
        Matrix K_mer = new Matrix(new double[][]{{1236.077, 0.0,      512},
                                                 {0.0,      1236.077, 512},
                                                 {0.0,      0.0,      1}});
        
        
        // Input depth map
        File map = new File("/home/nickrowell/Conferences and "
                + "Presentations/STAR-Seminar/images/dense_stereo_correlation/PANGU images/depth_map.matrix");
        
        // Output PXN file
        File outf = new File("/home/nickrowell/test.pxn");
        
        DepthMapToSurfaceModel model = new DepthMapToSurfaceModel(map, 512, 512, K_pangu, outf);
        
        
        
    }
    
    public DepthMapToSurfaceModel(File map, int w, int h, Matrix K, File pxn_file) throws IOException
    {
        
        // Inverse of the camera matrix is used to get line-of-sight vectors
        Matrix invK = K.inverse();
        
        // Each vertex corresponds to a pixel in depth map
        Vertex[][] verts = new Vertex[h][w];
        
        // Open reader on output file
        BufferedWriter out = new BufferedWriter(new FileWriter(pxn_file));
        
        // Write PXN header
        out.write("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
        out.write("<pangu_model ver=\"0\">\n");  
        
        
        // Read depth map
        BufferedReader in = new BufferedReader(new FileReader(map));
        
        String row;
        
        for(int r=0; r<h; r++)
        {
            row = in.readLine();
            
            // Use a Scanner to parse depth values
            Scanner scan = new Scanner(row);
            
            for(int c=0; c<w; c++)
            {
                double range = scan.nextDouble();
                
                if(range > 0)
                {
                    // Get unit vector towards feature
                    Matrix v = invK.times(new Matrix(new double[][]{{h-r},{c},{1}}));
                    v.timesEquals(v.normF());
                    // Scale according to depth value
                    v.timesEquals(range);
                    verts[r][c] = new Vertex(v.get(0,0), v.get(1,0), v.get(2,0));
                }
                // Invalid depth value
                else
                    verts[r][c] = null;
            }
        }
        
        // Now loop along rows and link 3D points into triangular strips
        
        // Write out vertices
        out.write("<mesh>\n");
        out.write("<positions>\n");
        // Number of vertex in list
        int vert_number = 0;
        for(int r=0; r<h; r++)
        {
            for(int c=0; c<w; c++)
            {
                if(verts[r][c]!=null)
                {
                    out.write("<vec3> "+verts[r][c].toString()+" </vec3>\n");

                    // So that valid vertices know their position in the list
                    verts[r][c].N = vert_number;
                    vert_number++;
                }
            }
        }
        out.write("</positions>\n");
        
        
        // Write out surface normals
        out.write("<normals>\n");
        for(int r=0; r<h; r++)
        {
            for(int c=0; c<w; c++)
            {
                if(verts[r][c]!=null)
                {
                    
                    // Check among it's eight nearest neighbours, and average the
                    // outward-facing surface normals for any facets that can be
                    // constructed.
                    if(verts[r-1][c-1]!=null && verts[r-1][c]!=null)
                    {
                        double[] norm = Vertex.getClockwiseSurfaceNormal(verts[r][c], verts[r-1][c-1], verts[r-1][c]);
                        verts[r][c].n0 += norm[0];
                        verts[r][c].n1 += norm[1];
                        verts[r][c].n2 += norm[2];
                    }
                    if(verts[r-1][c]!=null && verts[r-1][c+1]!=null)
                    {
                        double[] norm = Vertex.getClockwiseSurfaceNormal(verts[r][c], verts[r-1][c], verts[r-1][c+1]);
                        verts[r][c].n0 += norm[0];
                        verts[r][c].n1 += norm[1];
                        verts[r][c].n2 += norm[2];
                    }
                    if(verts[r-1][c+1]!=null && verts[r][c+1]!=null)
                    {
                        double[] norm = Vertex.getClockwiseSurfaceNormal(verts[r][c], verts[r-1][c+1], verts[r][c+1]);
                        verts[r][c].n0 += norm[0];
                        verts[r][c].n1 += norm[1];
                        verts[r][c].n2 += norm[2];
                    }
                    if(verts[r][c+1]!=null && verts[r+1][c+1]!=null)
                    {
                        double[] norm = Vertex.getClockwiseSurfaceNormal(verts[r][c], verts[r][c+1], verts[r+1][c+1]);
                        verts[r][c].n0 += norm[0];
                        verts[r][c].n1 += norm[1];
                        verts[r][c].n2 += norm[2];
                    }
                    if(verts[r+1][c+1]!=null && verts[r+1][c]!=null)
                    {
                        double[] norm = Vertex.getClockwiseSurfaceNormal(verts[r][c], verts[r+1][c+1], verts[r+1][c]);
                        verts[r][c].n0 += norm[0];
                        verts[r][c].n1 += norm[1];
                        verts[r][c].n2 += norm[2];
                    }
                    if(verts[r+1][c]!=null && verts[r+1][c-1]!=null)
                    {
                        double[] norm = Vertex.getClockwiseSurfaceNormal(verts[r][c], verts[r+1][c], verts[r+1][c-1]);
                        verts[r][c].n0 += norm[0];
                        verts[r][c].n1 += norm[1];
                        verts[r][c].n2 += norm[2];
                    }
                    if(verts[r+1][c-1]!=null && verts[r][c-1]!=null)
                    {
                        double[] norm = Vertex.getClockwiseSurfaceNormal(verts[r][c], verts[r+1][c-1], verts[r][c-1]);
                        verts[r][c].n0 += norm[0];
                        verts[r][c].n1 += norm[1];
                        verts[r][c].n2 += norm[2];
                    }
                    if(verts[r][c-1]!=null && verts[r-1][c-1]!=null)
                    {
                        double[] norm = Vertex.getClockwiseSurfaceNormal(verts[r][c], verts[r][c-1], verts[r-1][c-1]);
                        verts[r][c].n0 += norm[0];
                        verts[r][c].n1 += norm[1];
                        verts[r][c].n2 += norm[2];
                    }
                    
                    out.write("<vec3> "+verts[r][c].getNormalString()+" </vec3>\n");
                }
            }
        }
        out.write("</normals>\n");
        
        out.write("<tristrips>\n");
        
        
        // We proceed along strips like this:
        //
        //  0    0
        //  |   /|  _
        //  |  / |  /|
        //  V /  V /
        //  0    0
        //
        // We always start at a top left point and check that it has valid vertices
        // below it and to the right. If so, we write these vertices out and proceed
        // to the right. If not, we terminate the triangular strip, move to the
        // right and start a new triangular strip.
        
        // Indicates if we are currently in an active triangular strip
        boolean active_strip = false;
        
        for(int c=0; c<w-1; c++)
        {
            
            for(int r=0; r<h; r++)
            {
                
                // If this vertex is not null
                if(verts[r][c]!=null)
                {
                    
                    if(verts[r][c+1]!=null)
                    {
                        // If we are already writing a triangular strip, we
                        // add this vertex to it because it will form a triangle
                        // with previous vertices
                        if(active_strip)
                        {
                            out.write(verts[r][c].N + " "+verts[r][c+1].N + " ");
                        }
                        
                        // Otherwise, we must look ahead to ensure that there
                        // are at least 3 valid points to form a minimal tristrip
                        else
                        {
                            if(verts[r+1][c]!=null)
                            {
                                // Yes - start a new triangular strip and write
                                // the first two points to it.
                                active_strip = true;
                                out.write("<tristrip>\t");
                                out.write(verts[r][c].N + " "+verts[r][c+1].N + " ");
                            }
                            
                        }
                        
                    }
                    else
                    {
                        // Current point is valid but next in strip is not.
                        // Terminate any active strip here.
                        if(active_strip)
                        {
                            active_strip=false;
                            out.write(verts[r][c].N + " ");
                            out.write("\t</tristrip>\n");
                        }
                    }
                    
                }
                // Point we have moved to is not part of the surface model
                else
                {
                    if(active_strip)
                    {
                        active_strip=false;
                        out.write("\t</tristrip>\n");
                    }
                }
                
            }
            
        }
        out.write("</tristrips>\n");
        
        
        
        
        
        
        out.write("</mesh>\n");
        out.write("</pangu_model>\n");  // Closes pangu model
        out.flush();
        out.close();
    }
    
}
