package shape;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;

import numeric.geom.dim3.Vector3d;

/**
 *
 * @author nickrowell
 */
public class LatLonToPXN
{
    
    
    public static void main(String[] args) throws IOException
    {
        
        File input = new File("/home/nickrowell/Projects/PANGU_4/Documents/TN01/"
                + "Planetary_Science/Available_Data/Shape_Models/Moon/SELENE/"
                + "LALT_South_Pole/LALT_GT_SP_NUM.TAB");
        
        // Number of header lines.
        int NHEAD = 71;
        
        BufferedReader in = new BufferedReader(new FileReader(input));
        
        // Skip past header
        for(int i=0; i<NHEAD; i++) in.readLine();
        
        // Number of longitude samples
        int N_LONG = 360*32;
        // Number of latitude samples
        int N_LAT  = 11*128;
        
        // Number to skip
        int N_LAT_SKIP = 10*128;
        
        Vector3d[][] grid = new Vector3d[N_LAT-N_LAT_SKIP][N_LONG];
        
        String line;
        
        double lon=0, lat=0, r=0;
        double x,y,z;
        
        // Altitude reference point
        double ref_alt = 1737.4;  //km
        
        for(int i=0; i<N_LAT; i++)
        {
            for(int j=0; j<N_LONG; j++)
            {
                line=in.readLine();
                
                if(i>=N_LAT_SKIP)
                {
                    // Parse longitude, latitude and height...
                    Scanner scan = new Scanner(line);

                    lon = Math.toRadians(scan.nextDouble());
                    lat = Math.toRadians(scan.nextDouble());
                    r   = scan.nextDouble() + ref_alt;

                    // Convert to 3D point in Lunar frame
                    x = r*Math.cos(lat)*Math.cos(lon);
                    y = r*Math.cos(lat)*Math.sin(lon);
                    z = r*Math.sin(lat);

                    grid[i-N_LAT_SKIP][j] = new Vector3d(x,y,z);
                }
            }
            
            if(i%128==0) System.out.println("Latitude "+Math.toDegrees(lat));
            
        }
        
        System.out.println("Parsed TAB; writing PXN...");
        
        
        // Open writer on output file
        BufferedWriter out = new BufferedWriter(new FileWriter("/home/nickrowell/Projects/PANGU_4/Documents/TN01/"
                + "Planetary_Science/Available_Data/Shape_Models/Moon/SELENE/"
                + "LALT_South_Pole/LALT_GT_SP_NUM.PXN"));
        
        // Write PXN header
        out.write("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
        out.write("<pangu_model ver=\"0\">\n"); 
        
        // Cease recursion and write out mesh for this node
        out.write("<mesh>\n");
        
        // Write out vertex positions
        out.write("<positions>\n");
        for(int l=0; l<N_LAT-N_LAT_SKIP; l++)
        {
            for(int s=0; s<N_LONG; s++)
            {
                Vector3d pixel = grid[l][s];
                out.write("<vec3> "+pixel.getX()+" "+pixel.getY()+" "+pixel.getZ()+" </vec3>\n"); 
            }
        }
        out.write("</positions>\n");
        // Write out surface normals
        out.write("<normals>\n");
        for(int l=0; l<N_LAT-N_LAT_SKIP; l++)
        {
            for(int s=0; s<N_LONG; s++)
            {
                // Dummy normals along edges currently (can set these to normal
                // of nearest triangular facet)
                if(l==0 || l==(N_LAT-N_LAT_SKIP-1))
                {
                    out.write("<vec3> 0 0 1 </vec3>\n");
                    continue;
                }
                if(s==0 || s==(N_LONG-1))
                { 
                    out.write("<vec3> 0 0 1 </vec3>\n");
                    continue;
                }
                
                
                // Everywhere else: sample the surrounding vertices to estimate
                // the normal at the present point.
                //
                // Use the 4 nearest neighbours. These vertices are laid out:
                //
                //           l-->
                //
                //   s        E
                //   |    D   A   B
                //   V        C
                //
                Vector3d A = grid[l][s];
                Vector3d B = grid[l+1][s];
                Vector3d C = grid[l][s+1];
                Vector3d D = grid[l-1][s];
                Vector3d E = grid[l][s-1];
                
                // Calculate normals of each facet
                Vector3d ABC = Vector3d.getClockwiseSurfaceNormal(A, B, C);
                Vector3d ACD = Vector3d.getClockwiseSurfaceNormal(A, C, D);
                Vector3d ADE = Vector3d.getClockwiseSurfaceNormal(A, D, E);
                Vector3d AEB = Vector3d.getClockwiseSurfaceNormal(A, E, B);
                
                // Take average
                Vector3d norm = (ABC.add(ACD).add(ADE).add(AEB)).mult(0.25);
                
                out.write("<vec3> "+norm.getX()+" "+norm.getY()+" "+norm.getZ()+" </vec3>\n");
            }
        }
        out.write("</normals>\n");

        // Write out triangular strips
        out.write("<tristrips>\n");
        for(int l=0; l<N_LAT-N_LAT_SKIP; l++)
        {
            out.write("<tristrip>\t");
            for(int s=0; s<=N_LONG; s++)
            {
                out.write((l*(N_LONG+1) + s) + " ");
                out.write(((l+1)*(N_LONG+1) + s) + " ");
            }
            out.write("\t</tristrip>\n");
        }
        out.write("</tristrips>\n");
        
        out.write("</mesh>\n");
        
        out.write("</pangu_model>\n");  // Closes pangu model
        out.flush();
        out.close();
        
    }

}