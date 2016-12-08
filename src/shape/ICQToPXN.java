package shape;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.NoSuchElementException;
import java.util.Scanner;

import numeric.data.Histogram;


/**
 * This program is used to convert asteroid shape models from Implicitly
 * Connected Quadrilateral (ICQ) format to PXN format.
 * 
 * It takes care to purge duplicated vertices from the vertex list, and 
 * corrects the triangular strip indices to account for this. Vertex normals
 * are derived by averaging the surface normals for all triangular facets that
 * each is connected to.
 * 
 * Each cube face is broken down recursively into quarters so that the longest
 * triangular strip is less than 33 triangles. Each is given a unique node in
 * the PXN file. This greatly improves efficiency of back projections (e.g. for
 * RADAR simulation) because it allows far more appropriate bounding volumes to
 * be calculated.
 * 
 * 
 * 
 * @author nickrowell
 */
public class ICQToPXN 
{
    
    /** Magnitude parameter q. */
    public int q;
    
    /** Array of all 6*(q+1)^2 vertices. */
    public Vertex[] allVerts;
    
    /** Array of 6*q^2 unique vertices. */
    public Vertex[] uniqueVerts;
    
    /** 6*(q+1)^2 long index of vertex number to position in uniqueVerts array. */
    public int[] vertNumbers;
    
    /** Handle to input ICQ file. */
    File icq_file;
    /** Final output PXN file. */
    File pxn_file;
    
    
    /** Flag indicating whether to scale vertices by factor 1000 (km->m units). */
    boolean KM_TO_M = true;
    
    /** Flag indicating whether texture coordinates should be generated for each vertex. */
    boolean TEXTURE = false;
    
    /** Texture file name (must be ppm file of size 2^N-by-2^N pixels). */
    String tex_name = "/home/nickrowell/Projects/NEOGNC_2/Software/PANGU/models/itokawa/v2_textured/2048x2048_8sided.ppm";
    
    /** 
     * Texture scale factor (number of times texture is repeated over each face
     *  of the ICQ cube).
     */
    double tex_scale = 5.0;
    
    /**
     * We calculate some statistics of the model while it is being processed
     * and output these to the console.
     * 
     */
    // Max distance of any vertex from the origin
    double max_radius2 = 0.0;
    // Extent of model along each axis
    double max_x = 0.0;
    double min_x = 0.0;
    double max_y = 0.0;
    double min_y = 0.0;
    double max_z = 0.0;
    double min_z = 0.0;
    // Histogram of distances between neighbouring vertices
    Histogram VERTEX_PITCH = new Histogram(0.0, 2.0, 0.01, true);
    
    
    public static void main(String[] args) throws IOException
    {
        // Parent directory containing triangular facet model
        String parent = "/home/nickrowell/Shape_Models/Tethys/CO_SA_ISSNA_5_TETHYSSHAPE_V1_0/data";
        // Name of input file (triangular facet format)
        File input_file = new File(parent, "tethys_quad512q.tab");
        // Name out output PXN file (will be written to same directory)
        File pxn_file = new File(parent, "tethys_q512.pxn");
        
        ICQToPXN icqToPXN = new ICQToPXN(input_file, pxn_file, true);
        
    }
    
    public ICQToPXN(File _icq_file, File _pxn_file, boolean verbose) throws IOException
    {
        
        icq_file = _icq_file;
        pxn_file = _pxn_file;
        
        // Parse ICQ model from file
        try 
        {
            parseIcq(icq_file, verbose);
        }
        catch (IOException ex) 
        {
            System.err.println("Error reading file "+icq_file.toString());
            System.exit(1);
        }       
        
        System.out.println("Parsed ICQ model: "+allVerts.length+" vertices.");
        
     	// Write PXn model to file
        try 
        {
            writePxn(pxn_file, verbose);
        }
        catch (IOException ex) 
        {
            System.err.println("Error reading file "+icq_file.toString());
            System.exit(1);
        } 
        
    }
    
    public final void parseIcq(File input, boolean verbose) throws IOException
    {
        
        // Open reader on input file
        BufferedReader in = new BufferedReader(new FileReader(input));

        // Read first line
        String line = in.readLine();
        
        // Open scanner on first line
        Scanner scan = new Scanner(line);
        
        // Read parameter q
        q = scan.nextInt();
        
        allVerts    = new Vertex[6*(q+1)*(q+1)];
        uniqueVerts = new Vertex[6*q*q + 2];
        vertNumbers = new int[6*(q+1)*(q+1)];
        
        // Read all vertices
        if(verbose) System.out.println("ICQToPXN: reading vertex list...");
        
        // Loop over six faces
        for(int k = 1; k<=6; k++)
            // Loop over j coordinate
            for(int j=0; j<=q; j++)
                // Loop over i coordinate
                for(int i = 0; i<=q; i++)
                {
                    // Read in the corresponding line
                    line = in.readLine();
                    // Parse vertex and add to array
                    allVerts[getVertexNumber(i,j,k)] = parseVertex(getVertexNumber(i,j,k), line);
                    
                    if(KM_TO_M)
                    {
                        allVerts[getVertexNumber(i,j,k)].x *= 1000;
                        allVerts[getVertexNumber(i,j,k)].y *= 1000;
                        allVerts[getVertexNumber(i,j,k)].z *= 1000;
                    }
                    
                    // Initialise stats on first vertex
                    if(i==0 && j==0 && k==1)
                    {
                        Vertex tmp = allVerts[getVertexNumber(i,j,k)];
                        
                        max_radius2 = tmp.x*tmp.x + tmp.y*tmp.y + tmp.z*tmp.z;
                        // Extent of model along each axis
                        max_x = tmp.x;
                        min_x = tmp.x;
                        max_y = tmp.y;
                        min_y = tmp.y;
                        max_z = tmp.z;
                        min_z = tmp.z;
                    }
                    // Track stats on remaining vertices
                    else
                    {
                        Vertex tmp = allVerts[getVertexNumber(i,j,k)];
                        
                        double r2 = tmp.x*tmp.x + tmp.y*tmp.y + tmp.z*tmp.z;
                        if(r2 > max_radius2)
                            max_radius2 = r2;
                        if(tmp.x > max_x) max_x = tmp.x;
                        if(tmp.x < min_x) min_x = tmp.x;
                        if(tmp.y > max_y) max_y = tmp.y;
                        if(tmp.y < min_y) min_y = tmp.y;
                        if(tmp.z > max_z) max_z = tmp.z;
                        if(tmp.z < min_z) min_z = tmp.z;
                    }
                    
                }
        
        if(verbose)
        {
            System.out.println("Statistics of model:");
            System.out.println("Extent along X axis = ["+min_x+":"+max_x+"]");
            System.out.println("Extent along Y axis = ["+min_y+":"+max_y+"]");
            System.out.println("Extent along Z axis = ["+min_z+":"+max_z+"]");
            System.out.println("Maximum distance of any vertex from the origin"
                    + " = "+Math.sqrt(max_radius2));
            System.out.println("Number of vertices = "+allVerts.length);
            System.out.println("Number of unique vertices = "+uniqueVerts.length);
            
            // Now calculate histogram of vertex pitch (distance between
            // neighbouring vertices)
            // Loop over six faces
            for(int k = 1; k<=6; k++)
                // Loop over j coordinate
                for(int j=0; j<q; j++)
                    // Loop over i coordinate
                    for(int i = 0; i<q; i++)
                    {
                        Vertex r11 = allVerts[getVertexNumber(i,j,k)];
                        Vertex r12 = allVerts[getVertexNumber(i+1,j,k)];
                        Vertex r21 = allVerts[getVertexNumber(i,j+1,k)];
                        
                        // Distance the current vertex and the next one in the i direction...
                        double di = Math.sqrt((r11.x-r12.x)*(r11.x-r12.x) + 
                                              (r11.y-r12.y)*(r11.y-r12.y) +
                                              (r11.z-r12.z)*(r11.z-r12.z));
                        // ...and in the j direction...
                        double dj = Math.sqrt((r11.x-r21.x)*(r11.x-r21.x) + 
                                              (r11.y-r21.y)*(r11.y-r21.y) +
                                              (r11.z-r21.z)*(r11.z-r21.z));
                        
                        VERTEX_PITCH.add(di);
                        VERTEX_PITCH.add(dj);
                        
                    }
            
            // Now print vertex pitch histogram
            System.out.println("Vertex pitch frequency:");
            System.out.println(VERTEX_PITCH.print(true));
            
        }
        

        // Close file reader.
        in.close();
        
        // Copy non-edge vertices to uniqueVerts array:
        int index = 0;
        
        for(int k = 1; k<=6; k++)
            for(int j=1; j<q; j++)
                for(int i = 1; i<q; i++)
                {
                    uniqueVerts[index++] = allVerts[getVertexNumber(i,j,k)];
                    vertNumbers[getVertexNumber(i,j,k)] = index-1;
                }
  
        
        if(verbose) System.out.println("ICQToPAN: purging duplicate vertices...");
        
        // Now copy single instance of each duplicated vertex to uniqueVerts array
        for(int a=1; a<q; a++)
        {
            uniqueVerts[index++] = allVerts[getVertexNumber(0,a,1)];
            vertNumbers[getVertexNumber(0,a,1)] = index-1;
            vertNumbers[getVertexNumber(a, 0, 3)] = index-1;
            
            uniqueVerts[index++] = allVerts[getVertexNumber(0,a,3)];
            vertNumbers[getVertexNumber(0,a,3)] = index-1;
            vertNumbers[getVertexNumber(q, a, 4)] = index-1;            
            
            uniqueVerts[index++] = allVerts[getVertexNumber(0, a, 2)];
            vertNumbers[getVertexNumber(0, a, 2)] = index-1;
            vertNumbers[getVertexNumber(q, a, 3)] = index-1;             
            
            uniqueVerts[index++] = allVerts[getVertexNumber(a, q, 1)];
            vertNumbers[getVertexNumber(a, q, 1)] = index-1;
            vertNumbers[getVertexNumber(a, 0, 2)] = index-1;              
            
            uniqueVerts[index++] = allVerts[getVertexNumber(q, q - a, 1)];
            vertNumbers[getVertexNumber(q, q - a, 1)] = index-1;
            vertNumbers[getVertexNumber(a, 0, 5)] = index-1;              
            
            uniqueVerts[index++] = allVerts[getVertexNumber(q - a, 0, 1)];
            vertNumbers[getVertexNumber(q - a, 0, 1)] = index-1;
            vertNumbers[getVertexNumber(a, 0, 4)] = index-1;              
            
            uniqueVerts[index++] = allVerts[getVertexNumber(q-a,q,4)];
            vertNumbers[getVertexNumber(q-a,q,4)] = index-1;
            vertNumbers[getVertexNumber(a, q, 6)] = index-1;              
            
            uniqueVerts[index++] = allVerts[getVertexNumber(a, q, 2)];
            vertNumbers[getVertexNumber(a, q, 2)] = index-1;
            vertNumbers[getVertexNumber(a, 0, 6)] = index-1;              
            
            uniqueVerts[index++] = allVerts[getVertexNumber(a, q, 5)];
            vertNumbers[getVertexNumber(a, q, 5)] = index-1;
            vertNumbers[getVertexNumber(q, a, 6)] = index-1;              
            
            uniqueVerts[index++] = allVerts[getVertexNumber(0, a, 4)];
            vertNumbers[getVertexNumber(0, a, 4)] = index-1;
            vertNumbers[getVertexNumber(q, a, 5)] = index-1;              
            
            uniqueVerts[index++] = allVerts[getVertexNumber(q - a, q, 3)];
            vertNumbers[getVertexNumber(q - a, q, 3)] = index-1;
            vertNumbers[getVertexNumber(0, a, 6)] = index-1;              
            
            uniqueVerts[index++] = allVerts[getVertexNumber(q, a, 2)];
            vertNumbers[getVertexNumber(q, a, 2)] = index-1;
            vertNumbers[getVertexNumber(0, a, 5)] = index-1;              
            
        }
        
        // Purge 8 triplicated corner vertices.
        uniqueVerts[index++] = allVerts[getVertexNumber(0, 0, 1)];
        vertNumbers[getVertexNumber(0, 0, 1)] = index-1;
        vertNumbers[getVertexNumber(q, 0, 4)] = index-1;
        vertNumbers[getVertexNumber(0, 0, 3)] = index-1;
        
        uniqueVerts[index++] = allVerts[getVertexNumber(q, 0, 3)];
        vertNumbers[getVertexNumber(q, 0, 3)] = index-1;
        vertNumbers[getVertexNumber(0, q, 1)] = index-1;
        vertNumbers[getVertexNumber(0, 0, 2)] = index-1;
        
        uniqueVerts[index++] = allVerts[getVertexNumber(q, 0, 5)];
        vertNumbers[getVertexNumber(q, 0, 5)] = index-1;
        vertNumbers[getVertexNumber(q, 0, 1)] = index-1;
        vertNumbers[getVertexNumber(0, 0, 4)] = index-1;
        
        uniqueVerts[index++] = allVerts[getVertexNumber(q, 0, 2)];
        vertNumbers[getVertexNumber(q, 0, 2)] = index-1;
        vertNumbers[getVertexNumber(q, q, 1)] = index-1;
        vertNumbers[getVertexNumber(0, 0, 5)] = index-1;
        
        uniqueVerts[index++] = allVerts[getVertexNumber(q, q, 3)];
        vertNumbers[getVertexNumber(q, q, 3)] = index-1;
        vertNumbers[getVertexNumber(0, 0, 6)] = index-1;
        vertNumbers[getVertexNumber(0, q, 2)] = index-1;
        
        uniqueVerts[index++] = allVerts[getVertexNumber(q, q, 4)];
        vertNumbers[getVertexNumber(q, q, 4)] = index-1;
        vertNumbers[getVertexNumber(0, q, 6)] = index-1;
        vertNumbers[getVertexNumber(0, q, 3)] = index-1;
        
        uniqueVerts[index++] = allVerts[getVertexNumber(q, q, 2)];
        vertNumbers[getVertexNumber(q, q, 2)] = index-1;
        vertNumbers[getVertexNumber(q, 0, 6)] = index-1;
        vertNumbers[getVertexNumber(0, q, 5)] = index-1;
        
        uniqueVerts[index++] = allVerts[getVertexNumber(q, q, 5)];
        vertNumbers[getVertexNumber(q, q, 5)] = index-1;
        vertNumbers[getVertexNumber(q, q, 6)] = index-1;
        vertNumbers[getVertexNumber(0, q, 4)] = index-1;
        
        // Now find all triangular Facets that each vertex is connected to
        if(verbose) System.out.println("ICQToPAN: calculating normals for vertices...");
        
        // Loop over six faces
        for(int k = 1; k<=6; k++)
        {   
            // Loop over j coordinate
            for(int j=0; j<q; j++)
            {
                // Loop over i coordinate
                for(int i = 0; i<q; i++)
                {
                    // There are two triangular facets. The first connects
                    // vertices (i,j),(i,j+1) and (i+1,j). The second connects
                    // vertices (i,j+1),(i+1,j) and (i+1,j+1).
                    
                    // Get vertices
                    Vertex A = uniqueVerts[vertNumbers[getVertexNumber(i, j, k)]];
                    Vertex B = uniqueVerts[vertNumbers[getVertexNumber(i+1, j, k)]];
                    Vertex C = uniqueVerts[vertNumbers[getVertexNumber(i, j+1, k)]];
                    Vertex D = uniqueVerts[vertNumbers[getVertexNumber(i+1, j+1, k)]];                    
                    
                    // Surface normal for facet formed by ACB (clockwise winding order)
                    double[] ACB = Vertex.getClockwiseSurfaceNormal(A, C, B);
                    A.n0 += ACB[0]; A.n1 += ACB[1]; A.n2 += ACB[2];
                    C.n0 += ACB[0]; C.n1 += ACB[1]; C.n2 += ACB[2];
                    B.n0 += ACB[0]; B.n1 += ACB[1]; B.n2 += ACB[2];
                    
                    // Surface normal for facet formed by BCD (clockwise winding order)
                    double[] BCD = Vertex.getClockwiseSurfaceNormal(B, C, D);
                    B.n0 += BCD[0]; B.n1 += BCD[1]; B.n2 += BCD[2];
                    C.n0 += BCD[0]; C.n1 += BCD[1]; C.n2 += BCD[2];
                    D.n0 += BCD[0]; D.n1 += BCD[1]; D.n2 += BCD[2];
                    
                }
            }
        }
        
        
    }
    
    /** Calculate vertex number L_v from face number and i,j coordinate.  */
    public int getVertexNumber(int i, int j, int k)
    {
        //return 1 + i + (q+1)*j + (q+1)*(q+1)*(k-1);
        
        // Zero-based indexing, more suited to using arrays to store data.
        return i + (q+1)*j + (q+1)*(q+1)*(k-1);
    }
    
    
    /**
     * Method for parsing vertices from a string, applicable to ICQ models
     * where the string doesn't include a vertex number.
     * @param str
     * @return 
     */
    private Vertex parseVertex(int N, String str)
    {
        // Open scanner on String
        Scanner scan = new Scanner(str);

        try
        {
            // Create new Vertex
            return new Vertex(N, scan.nextDouble(), scan.nextDouble(), scan.nextDouble());
        } 
        catch (NoSuchElementException nsee) 
        {
            throw new RuntimeException("Missing coordinate on vertex " + str);
        } 
        catch (NumberFormatException nfe) 
        {
            throw new RuntimeException("Number format exception on vertex " + str);
        }
    
    }
    
    
    private void writePxn(File output, boolean verbose) throws IOException
    {
    	// Open writer on output file
        BufferedWriter out = new BufferedWriter(new FileWriter(output));
        
        // Write PXN header
        out.write("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
        out.write("<pangu_model ver=\"0\">\n");        
        
        if(TEXTURE)
        {
            out.write("<texture type=\"sphere\">\n");
            out.write("<data><include type=\"raw\" file=\""+tex_name+"\"/></data>\n");
            out.write("<child>\n");
        }
        
        // Top level node - contains 6 faces
        out.write("<group count=\"6\">\n");
        // Loop over 6 faces of ICQ cube
        for(int f=1; f<=6; f++)
        {
            // Recursively divide face into sub-nodes and write out to file.
            //
            // On my graphics card, there seems to be a sweet spot in
            // image + radar generation time when a limit of 32 triangles is
            // placed on the length of any strip. This means halting the 
            // recursion when sub-nodes are 16*16 vertices in size.
            mkMesh(out, f, q, 0, 0, 16);
        }
        out.write("</group>\n");        // Closes top level node (6 faces)
        
        if(TEXTURE)
        {
            out.write("</child>\n");
            out.write("</texture>\n");
        }
        
        out.write("</pangu_model>\n");  // Closes pangu model
        out.flush();
        out.close();
    }
    
    
    
    /**
     * 
     * @param out       BufferedWriter to PXN file
     * @param f         Face index
     * @param Q         Current resolution index (size of sub-mesh)
     * @param stride_i  } Coordinates in face (i,j) sub-mesh origin
     * @param stride_j  }
     * @param max_tri   Maximum number of triangles in the smallest strip. This
     *                  determines point at which recursion stops and we write
     *                  out a mesh.
     */
    private void mkMesh(BufferedWriter out, int f, int Q, int stride_i, int stride_j, int max_tri) 
            throws IOException
    {
        
        // Has recursion limit been reached?
        if(Q<=max_tri)
        {
            // Cease recursion and write out mesh for this node
            out.write("<mesh>\n");
            // Write out vertex positions
            out.write("<positions>\n");
            for(int i=stride_i; i<=stride_i + Q; i++)
            {
                for(int j=stride_j; j<=stride_j + Q; j++)
                {
                    out.write("<vec3> "+uniqueVerts[vertNumbers[getVertexNumber(i, j, f)]].toString()+" </vec3>\n"); 
                }
            }
            out.write("</positions>\n");
            // Write out surface normals
            out.write("<normals>\n");
            for(int i=stride_i; i<=stride_i + Q; i++)
            {
                for(int j=stride_j; j<=stride_j + Q; j++)
                {
                    out.write("<vec3> "+uniqueVerts[vertNumbers[getVertexNumber(i, j, f)]].getNormalString()+" </vec3>\n"); 
                }
            }
            out.write("</normals>\n");
            
            // Optionally write a texcoord for each vertex based on projection
            // of vertex onto corresponding face of ICQ model
            if(TEXTURE)
            {
                out.write("<texcoords>\n");
                for(int i=stride_i; i<=stride_i + Q; i++)
                {
                    for(int j=stride_j; j<=stride_j + Q; j++)
                    {
                        // Coordinate of this vertex in the ICQ face is (i,j)
                        // i and j range from 0 to q-1: map this to range 0:1
                        // for texture coordinate
                        double u = (double)i / (double)(q-1);
                        double v = (double)j / (double)(q-1);
                        
                        // Apply texture scale
                        u *= tex_scale;
                        v *= tex_scale;
                        
                        // Third texture coordinate always zero (only 2D
                        // textures are supported).
                        
                        out.write("<vec3> "+u+" "+v+" 0 </vec3>\n"); 
                    }
                }
                out.write("</texcoords>\n");
            }
            
            // Write out triangular strips
            out.write("<tristrips>\n");
            for(int i=0; i<Q; i++)
            {
                out.write("<tristrip>\t");
                for(int j=0; j<=Q; j++)
                {
                    out.write((i*(Q+1) + j) + " ");
                    out.write(((i+1)*(Q+1) + j) + " ");
                }
                out.write("\t</tristrip>\n");
            }
            out.write("</tristrips>\n");
            out.write("</mesh>\n");
            
        }
        else
        {   
            out.write("<group count=\"4\">\n");
            // Quarter sub-mesh and continue recursion on each part
            mkMesh(out, f, Q/2, stride_i, stride_j, max_tri);
            mkMesh(out, f, Q/2, stride_i + Q/2, stride_j, max_tri);
            mkMesh(out, f, Q/2, stride_i, stride_j + Q/2, max_tri);
            mkMesh(out, f, Q/2, stride_i + Q/2, stride_j + Q/2, max_tri);
            out.write("</group>\n");
        }
    
    
    }
    
    /**
     * Run pxn2pan to convert new pxn file to pan file.
     * @param pan_file
     * @param pxn_file
     * @throws IOException 
     */
//    public final void writePanFile(File pan_file, File pxn_file) throws IOException
//    {
//        
//        String[] commands = new String[]{"pxn2pan","-lic","/opt/pangu/license/nrowell.txt","-o",pan_file.toString(), pxn_file.toString()};
//        
//        for(String command : commands)
//            System.out.print(command+" ");
//        
//        // Passed basic checks: attempt to launch PANGU
//        ProcessBuilder pangu = new ProcessBuilder(commands);
//        
//        // Get working environment for pangu process
//        Map<String, String> env = pangu.environment();
//        
//        // Add PANGU /bin directory to LD_LIBRARY_PATH envvar for PANGU process
//        env.put("LD_LIBRARY_PATH", env.get("LD_LIBRARY_PATH") + ":/opt/pangu/pangu_3.30/bin");       
//
//        // Run pxn2pan
//        Process proc = pangu.start();
//        
//    }
    
}