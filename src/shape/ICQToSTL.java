package shape;

import java.io.*;
import java.util.LinkedList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Scanner;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import numeric.data.Histogram;


/**
 * This program is used to convert asteroid shape models from Implicitly
 * Connected Quadrilateral (ICQ) format to binary or ascii STL format.
 * 
 * It takes care to purge duplicated vertices from the vertex list, and 
 * corrects the triangular strip indices to account for this. Vertex normals
 * are derived by averaging the surface normals for all triangular facets that
 * each is connected to.
 * 
 * @author nickrowell
 */
public class ICQToSTL 
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
    
    /** Final output STL file. */
    File stl_file;
    
    /** Flag indicates whether to write ascii or binary file */
    boolean ASCII = true;
    
    /** Flag indicating whether to scale vertices by factor 1000 (km->m units). */
    boolean KM_TO_M = false;
    
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
    	
    	// Here we convert all ICQ models in directory to STL format using multiple threads
    	
    	final List<Future<String>> futures = new LinkedList<Future<String>>();
    	
        // setup the executor
        ExecutorService executor = Executors.newFixedThreadPool(8);
    	
    	// Parent directory of all ICQ models
    	String top_level = "/home/nrowell/Shape_Models/";
    	
//    	for(String body : new String[]{"Dione","Eros","Itokawa","Mimas","Phobos","Phoebe","Tethys"})
    	for(String body : new String[]{"Dione"})	
    	{
    	
    		// Directories for ICQ and STL models for this body
    		String icq = top_level + body+ "/icq/" ;
    		String stl = top_level + body+ "/stl/" ;
    		
    		// Loop over resolution levels
//	    	for(String q : new String[]{"64","128","256","512"})
	    	for(String q : new String[]{"256"})
	    	{
		        // Name of input ICQ file
		        final File icq_file = new File(icq, "quad"+q+"q.tab");
		        // Name out output STL file
		        final File stl_file = new File(stl, body+"_q"+q+".stl");
		        
		        // Create Callable to perform conversion of one file
		        final Callable<String> worker = new Callable<String>() 
                {
                    @Override
                    public String call()
                    {
                    	try
                    	{
                    		// Perform the conversion
                    		new ICQToSTL(icq_file, stl_file, false);
                    	}
                    	catch(IOException ioe)
                    	{
                    		System.out.println("Caught exception processing "+icq_file.getName()+": "+ioe.getMessage());
                    	}
                    	
                    	// Returns name of file it just finished processing, so we can keep track of execution.
                    	return icq_file.getPath();
                    }
                };
                
                // Add to list
                futures.add(executor.submit(worker));
                
	    	}
    	}
    	
		// shutdown the execution
        executor.shutdown();
    	
        // Wait for jobs to complete
        for (final Future<String> future : futures)
        {
        	try 
        	{
				System.out.println("Finished converting "+future.get());
			} 
        	catch (InterruptedException e) 
        	{
				// TODO Auto-generated catch block
				e.printStackTrace();
			} 
        	catch (ExecutionException e) 
        	{
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
        }
    	
    	
    }
    
    public ICQToSTL(File picq_file, File pstl_file, boolean verbose) throws IOException
    {
        
        icq_file = picq_file;
        stl_file = pstl_file;
        
        // Parse ICQ model from file
        try 
        {
            parseICQModel(icq_file, verbose);
        }
        catch (IOException ex) 
        {
            System.err.println("Error reading file "+icq_file.toString());
            System.exit(1);
        }       
        
        System.out.println("Parsed ICQ model "+icq_file.getPath()+": "+allVerts.length+" vertices.");
        
        
        
        // Write STL model to file
        try 
        {
            if(ASCII) writeAsciiStl(stl_file, verbose);
            else      writeBinaryStl(stl_file, verbose);
        }
        catch (IOException ex) 
        {
            System.err.println("Error writing file "+stl_file.toString());
            System.exit(1);
        } 
        
    }
    
    public final void parseICQModel(File input, boolean verbose) throws IOException
    {
        
        // Open reader on input file
        BufferedReader in = new BufferedReader(new FileReader(input));

        // Read first line
        String line = in.readLine();
        
        // Open scanner on first line
        Scanner scan = new Scanner(line);
        
        // Read parameter q
        q = scan.nextInt();
        
        scan.close();
        
        allVerts    = new Vertex[6*(q+1)*(q+1)];
        uniqueVerts = new Vertex[6*q*q + 2];
        vertNumbers = new int[6*(q+1)*(q+1)];
        
        // Read all vertices
        if(verbose) System.out.println("ICQToSTL: reading vertex list for "+input.getName()+" ...");
        
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
  
        
        if(verbose) System.out.println("ICQToSTL: purging duplicate vertices...");
        
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
        if(verbose) System.out.println("ICQToSTL: calculating normals for vertices...");
        
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
        try
        {
        	// Open scanner on String
        	Scanner scan = new Scanner(str);
        	
        	double x = scan.nextDouble();
        	double y = scan.nextDouble();
        	double z = scan.nextDouble();
        	
        	scan.close();
        	
            // Create & return new Vertex
            return new Vertex(N, x, y, z);
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
    
    
    private void writeAsciiStl(File output, boolean verbose) throws IOException
    {
    	// Open writer on output file
        BufferedWriter out = new BufferedWriter(new FileWriter(output));
        
        // Write STL header
        out.write("solid "+output.getName()+"\n");
        
        // Loop over 6 faces of ICQ cube
        for(int f=1; f<=6; f++)
        {
            // Each face is consists of a grid of vertices of size
        	// q*q, with coordinates of vertices labelled by i and j.
        	// Conceptually, the origin is at the top left corner of a
        	// face with the i coordinate increasing right and the j coordinate
        	// increasing down. This matches the convention used for pixel 
        	// coordinates. We compose the grid into a sequence of triangular
        	// facets as follows:
        	//  - Each grid cell is divided into two facets
        	//  - The cell with top left corner at (i,j) is divided into
        	//    two facets with vertices at (i,j),(i+1,j+1),(i,j+1) and
        	//    (i,j),(i+1,j),(i+1,j+1). This is the clockwise winding order.
        	//  - The outward-pointing surface normal is computed from the
        	//    vertices using standard geometry methods.
        	
        	
        	// Loop over grid points. For each value of (i,j) we compute the triangular
        	// facets for the grid cell extending to (i+1,j+1), so we exclude the final
        	// coordinate (i=q,j=q) from the loop
        	
        	// Handles to vertices of current facet
        	Vertex r0,r1,r2;
        	double[] n;
        	
        	// Loop over i coordinate
	        for(int i=0; i<q; i++)
	        {
	            // Loop over j coordinate
	            for(int j = 0; j<q; j++)
	            {
	            	// Do first facet with vertices (i,j),(i,j+1),(i+1,j+1)
	            	r0 = uniqueVerts[vertNumbers[getVertexNumber(i,   j,   f)]];
	            	r1 = uniqueVerts[vertNumbers[getVertexNumber(i,   j+1, f)]];
	            	r2 = uniqueVerts[vertNumbers[getVertexNumber(i+1, j+1, f)]];
	            	n = Vertex.getClockwiseSurfaceNormal(r0, r1, r2);
	            	
	            	out.write(String.format("facet normal %f %f %f \n", n[0], n[1], n[2]));
	            	out.write("  outer loop\n");
	            	out.write(String.format("    vertex %s\n",r0.toString()));
	            	out.write(String.format("    vertex %s\n",r1.toString()));
	            	out.write(String.format("    vertex %s\n",r2.toString()));
	            	out.write("  endloop\n");
	            	out.write("endfacet\n");
	            	
	            	// Do second facet with vertices (i,j),(i+1,j+1),(i+1,j)
	            	r0 = uniqueVerts[vertNumbers[getVertexNumber(i,   j,   f)]];
	            	r1 = uniqueVerts[vertNumbers[getVertexNumber(i+1, j+1, f)]];
	            	r2 = uniqueVerts[vertNumbers[getVertexNumber(i+1, j,   f)]];
	            	n = Vertex.getClockwiseSurfaceNormal(r0, r1, r2);
	            	
	            	out.write(String.format("facet normal %f %f %f \n", n[0], n[1], n[2]));
	            	out.write("  outer loop\n");
	            	out.write(String.format("    vertex %s\n",r0.toString()));
	            	out.write(String.format("    vertex %s\n",r1.toString()));
	            	out.write(String.format("    vertex %s\n",r2.toString()));
	            	out.write("  endloop\n");
	            	out.write("endfacet\n");
	            	
	            }
	        }
        	
        }
        
        // Write STL footer
        out.write("endsolid "+output.getName()+"\n");
        
        out.flush();
        out.close();
    }
    
    
    /**
     * NOTE: neither Blender nor Meshlab are able to load binary STL files written
     * using this method.
     * @param output
     * @param verbose
     * @throws IOException
     */
    private void writeBinaryStl(File output, boolean verbose) throws IOException
    {
    	// Use DataOutputStream to write binary data to file
        DataOutputStream out = new DataOutputStream(new FileOutputStream(output));
        
        // Write STL header (80 characters)
        for(int c=0; c<80; c++)
        	out.writeByte(0);
        
        // Write number of triangles:
        //  - 6 faces
        //  - (q-1)*(q-1) grid cells on each face
        //  - 2 triangles per cell
        int n_tri = 6 * (q-1)*(q-1) * 2;
        
        // NEED TO WRITE AN UNSIGNED INTEGER HERE.
        // For q=512; n_tri = 3133452, which overflows range of java signed int.
        out.writeInt(n_tri);
        
        // Now loop over each triangle and write to file. See
        // writeAsciiStl(File, boolean) for the logic.
        
        // Loop over 6 faces of ICQ cube
        for(int f=1; f<=6; f++)
        {
            
        	// Handles to vertices of current facet
        	Vertex r0,r1,r2;
        	double[] n;
        	
        	// Loop over i coordinate
	        for(int i=0; i<q; i++)
	        {
	            // Loop over j coordinate
	            for(int j = 0; j<q; j++)
	            {
	            	// Do first facet with vertices (i,j),(i,j+1),(i+1,j+1)
	            	r0 = uniqueVerts[vertNumbers[getVertexNumber(i,   j,   f)]];
	            	r1 = uniqueVerts[vertNumbers[getVertexNumber(i,   j+1, f)]];
	            	r2 = uniqueVerts[vertNumbers[getVertexNumber(i+1, j+1, f)]];
	            	n = Vertex.getClockwiseSurfaceNormal(r0, r1, r2);
	            	
	            	// Write normal
	            	out.writeFloat((float)n[0]);
	            	out.writeFloat((float)n[1]);
	            	out.writeFloat((float)n[2]);
	            	
	            	// Write each vertex
	            	out.writeFloat((float)r0.x);
	            	out.writeFloat((float)r0.y);
	            	out.writeFloat((float)r0.z);
	            	out.writeFloat((float)r1.x);
	            	out.writeFloat((float)r1.y);
	            	out.writeFloat((float)r1.z);
	            	out.writeFloat((float)r2.x);
	            	out.writeFloat((float)r2.y);
	            	out.writeFloat((float)r2.z);
	            	
	            	// Write 'attribute byte count': no attributes in this case
	            	out.writeShort(0);
	            	
	            	// Do second facet with vertices (i,j),(i+1,j+1),(i+1,j)
	            	r0 = uniqueVerts[vertNumbers[getVertexNumber(i,   j,   f)]];
	            	r1 = uniqueVerts[vertNumbers[getVertexNumber(i+1, j+1, f)]];
	            	r2 = uniqueVerts[vertNumbers[getVertexNumber(i+1, j,   f)]];
	            	n = Vertex.getClockwiseSurfaceNormal(r0, r1, r2);
	            	
	            	// Write normal
	            	out.writeFloat((float)n[0]);
	            	out.writeFloat((float)n[1]);
	            	out.writeFloat((float)n[2]);
	            	
	            	// Write each vertex
	            	out.writeFloat((float)r0.x);
	            	out.writeFloat((float)r0.y);
	            	out.writeFloat((float)r0.z);
	            	out.writeFloat((float)r1.x);
	            	out.writeFloat((float)r1.y);
	            	out.writeFloat((float)r1.z);
	            	out.writeFloat((float)r2.x);
	            	out.writeFloat((float)r2.y);
	            	out.writeFloat((float)r2.z);

	            	// Write 'attribute byte count': no attributes in this case
	            	out.writeShort(0);
	            	
	            }
	        }
        	
        }
        
        out.flush();
        out.close();
    }
    
}