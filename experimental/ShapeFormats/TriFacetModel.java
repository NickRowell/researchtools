package ShapeFormats;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Scanner;

import shape.Vertex;

/**
 * Supports TriFacetToPAN converter.
 * 
 * Could be sped up considerably by using arrays rather than lists. The number
 * of vertices & facets is included in file header usually.
 * 
 * @author nickrowell
 */
public class TriFacetModel {
    
    /**
     * List of Vertex objects in model. Each Vertex maintains a list of the
     * Facets that it is connected to, so that a reasonable surface normal
     * can be calculated by averaging the normals for each connected Facet.
     */
    public List<Vertex> verts = new LinkedList<Vertex>();
    
    /**
     * List of Facet objects in the model.
     */
    public List<Facet> facets = new LinkedList<Facet>();   
    
    /** Main constructor. */
    public TriFacetModel(File input)
    {
        try 
        {
            parseTriFacetModel(input);
        } 
        catch (IOException ex) 
        {
            System.err.println("Error reading file "+input.toString());
        }
    }
    
    
    public void parseTriFacetModel(File input) throws IOException
    {
        
        // Open reader on input file
        BufferedReader in = new BufferedReader(new FileReader(input));

        
         // Read first line
        String line = in.readLine();
        
        // Open scanner on first line
        Scanner scan = new Scanner(line);
        
        // Read number of vertices
        int n_verts  = scan.nextInt();
        // Read number of facets
        int n_facets = scan.nextInt();
        
        // List to store vertices. By construction, vertex N appears at 
        // position N-1 in the List. Vertices must be in consecutive order
        // in the file.
        verts = new LinkedList<Vertex>();
        
        // Now read all vertices into a List
        for(int n_verts_read = 1; n_verts_read <= n_verts; )
        {
            // Read next line
            String str = in.readLine();
            scan = new Scanner(str);
            // Skip blank lines
            if(!scan.hasNext()) continue;
            
            // Parse Vertex from this line
            Vertex vertex = parseVertex(str);
            
            // Check that vertices are in order. Could handle out-of-order
            // vertices by sorting list. But it's unlikely that vertices will
            // be out of order so not much point implementing this.
            if(vertex.N != n_verts_read)
                throw new RuntimeException("Vertex order broken at vertex "+n_verts_read);
            
            verts.add(vertex);
            n_verts_read++;
        }
        
        // List to store Facets. By construction, Facet N appears at 
        // position N-1 in the List. Facets must be in consecutive order
        // in the file.
        facets = new LinkedList<Facet>();
        
        // Now read all facets into a List
        for(int n_facets_read = 1; n_facets_read <= n_facets;)
        {
            // Read next line
            String str = in.readLine();
            scan = new Scanner(str);
            // Skip blank lines
            if(!scan.hasNext()) continue;
            // Parse Facet from this line
            Facet facet = Facet.parseFacet(str, verts);
            
            // Check that facets are in order. Could handle out-of-order
            // facets by sorting list. But it's unlikely that facets will
            // be out of order so not much point implementing this.
            if(facet.N != n_facets_read)
                throw new RuntimeException("Facet order broken at facet "+n_facets_read);
            
            facets.add(facet);
            n_facets_read++;
        }
        
    }
        
    /**
     * Method for parsing vertices from a string, applicable to triangular
     * facet models where the string includes a vertex number.
     * @param str
     * @return 
     */
    private Vertex parseVertex(String str)
    {
        // Open scanner on String
        Scanner scan = new Scanner(str);

        try
        {
            // Create new Vertex
            return new Vertex(scan.nextInt(), scan.nextDouble(), scan.nextDouble(), scan.nextDouble());
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
    
}
