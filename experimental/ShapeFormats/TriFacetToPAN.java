package ShapeFormats;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import shape.Vertex;


/**
 * Convert triangular facet model to PAN format, via intermediate PXN format.
 * This program attempts to generate long triangular strips from the list of 
 * facets, rather than simply specifying a separate triangular strip for each 
 * facet.
 * 
 * Currently in working state, but could do with some work to make it more
 * efficient. Best to use ICQToPAN if ICQ version is available, because this
 * is much higher fidelity and a lot faster.
 * 
 * To do:
 * 
 * Currently there is a problem concerning duplicated vertices: the Facets
 * connecting to duplicated vertices aren't added to all the vertex copies,
 * so the duplicate vertices don't get the surface normals set correctly.
 * 
 * This code could also be sped up considerably by using arrays rather than 
 * lists.
 * 
 * 
 * @author nickrowell
 */
public class TriFacetToPAN 
{
    
    File tri_file;
    File pxn_file;
    File pan_file;
    
    public static void main(String[] args) throws IOException
    {
        // Parent directory containing triangular facet model
        String parent = "/home/nickrowell/Shape_Models/Eros/NEAR_A_MSI_5_EROSSHAPE_V1_0/data/vertex";

        // Name of input file (triangular facet format)
        File input_file = new File(parent, "ver64q.tab");
        // Name out output PXN file (will be written to same directory)
        File pxn_file = new File(parent, "ver64q.pxn");
        // Name out output PAN file (will be written to same directory)
        File pan_file = new File(parent, "ver64q.pan");
        
        TriFacetToPAN triFacetToPAN = new TriFacetToPAN(input_file, pxn_file, pan_file);
        
    }
    
    public TriFacetToPAN(File _tri_file, File _pxn_file, File _pan_file) throws IOException
    {
        
        tri_file = _tri_file;
        pxn_file = _pxn_file;
        pan_file = _pan_file;
        
        // Parse TriFacet model from file
        TriFacetModel triFacetModel = new TriFacetModel(tri_file);
        
        System.out.println("Parsed Triangular Facet model: "+triFacetModel.verts.size()
                +" vertices and "+triFacetModel.facets.size()+" facets.");
        
        // Open writer on output file
        BufferedWriter out = new BufferedWriter(new FileWriter(pxn_file));
        
        // Write PXN header and vertex data
        PXNUtil.writePXNHeader(out);
        
        Vertex[] verts = new Vertex[triFacetModel.verts.size()];
        for(int i=0; i<verts.length; i++) verts[i] = triFacetModel.verts.get(i);
        
        PXNUtil.writePXNVertices(out, verts);
        PXNUtil.writePXNNormals(out, verts);

        System.out.println("Computing triangular strips...");
        
        // Count number of triangular strips...
        int n_strips = 0;
        
        // Now initialise tristrip searching algorithm by using first facet
        // to create a new TriStrip.
        TriStrip tristrip = new TriStrip(triFacetModel.facets.get(0));
        
        // Now attempt to add more facets to grow the strip
        for(int n_facet = 0; n_facet < triFacetModel.facets.size()-3; n_facet++)
        {
            
            // We keep 4 facets in memory as we slide along triangle strip. This
            // is necessary in order to detect reversals in the triangle winding
            // order (see comments below).
            Facet A = triFacetModel.facets.get(n_facet+0);   // Current 'active' Facet
            Facet B = triFacetModel.facets.get(n_facet+1);
            Facet C = triFacetModel.facets.get(n_facet+2);
            Facet D = triFacetModel.facets.get(n_facet+3);
            
            
            // A necessary but not sufficient condition for the triangular
            // strip to continue from facet A to B is that A and B must 
            // share precisely two vertices. If this condition is not met, then
            // strip ends on A and we start a new one from B.            
            if(!A.hasTwoVertsInCommon(B))
            {
                    
                // Triangular strip splits between A and B.

                // Finalise existing TriStrip (write remaining vertices
                // from current facet)
                tristrip.finalise();
                writeTriangularStrip(out, tristrip);
                //System.out.println("Completed "+(n_facet*100)/triFacetModel.facets.size()+" %");
                n_strips++;
                // Create new TriStrip intialised with new Facet
                tristrip = new TriStrip(B);

            }

            // If A and B share two vertices, then it is still possible that the
            // strip splits between A and B if the winding order of the 
            // triangles reverses. This can only be detected by looking ahead
            // a further two facets. The problem is that we need to 
            // distinguish this situation:
            //
            //       O------O------O
            //       |\     |\     |
            //       | \ B  | \ D  |
            //       |  \   |  \   |   (correct winding order)
            //       |   \  |   \  |
            //       | A  \ | C  \ |
            //       |     \|     \|
            //       O------O------O
            //
            //      from this situation:
            //
            //       O------O------O
            //       |\     |     /|
            //       | \ B  | C  / |
            //       |  \   |   /  |   (winding order reverses between A and D)
            //       |   \  |  /   |
            //       | A  \ | / D  |
            //       |     \|/     |
            //       O------O------O
            //
            // Notice that each facet is connected to the preceeding one by
            // two vertices, and the only way to detect the reversal is that A
            // and D share one vertex. Also, we can split the strip anywhere 
            // between A and D and still make two valid triangular strips. In 
            // this code, we split between A and B so that we keep the same 
            // onward processing as in the case where A and B don't share two 
            // vertices.            
            else if(A.hasOneVertInCommon(C)  &&
                    !A.hasNoVertsInCommon(D) &&     // <-- this is the condition that indicates winding order reversal.
                                                    // If A & D share one vertex, winding order has reversed; if two vertices,
                                                    // then strip has closed back on itself and we should split anyway.
                    B.hasTwoVertsInCommon(C) &&
                    B.hasOneVertInCommon(D) &&
                    C.hasTwoVertsInCommon(D) )
            {
                // Winding order reverses (or strip loops back on itself).
                tristrip.finalise();
                writeTriangularStrip(out, tristrip);
                //System.out.println("Completed "+(n_facet*100)/triFacetModel.facets.size()+" %");
                n_strips++;
                tristrip = new TriStrip(B);
            }
            // A and B form continuous triangular strip
            else
            {
                tristrip.add(B);
            }
            
            
            // When end of Facet list is reached, don't forget to add the 
            // remaining two Facets to the triangular strip. Don't need to 
            // check for winding order reversal here, because at least three
            // more facets would be required for winding order to reverse.
            if(n_facet == triFacetModel.facets.size()-4)
            {
                
                if(B.hasTwoVertsInCommon(C))
                    tristrip.add(C);
                else
                {
                    tristrip.finalise();
                    writeTriangularStrip(out, tristrip);
                    n_strips++;
                    tristrip = new TriStrip(C);              
                }
                
                if(C.hasTwoVertsInCommon(D) && C.hasOneVertInCommon(B))
                    tristrip.add(D);
                else
                {
                    tristrip.finalise();
                    writeTriangularStrip(out, tristrip);
                    n_strips++;
                    tristrip = new TriStrip(D);
                }                
                
                tristrip.finalise();
                writeTriangularStrip(out, tristrip);
                n_strips++;
                break;    
            }
            
            
        }
        
        // Print some stats...
        System.out.println("Found "+n_strips+" triangular strips with "+
         (triFacetModel.facets.size()/(double)n_strips)+" facets on average.");
        
        
        // Finish PXN file
        PXNUtil.writePXNFooter(out);
        out.close();
        
        // Now convert PXN file to PAN file by calling PANGU.
        PXNUtil.writePanFile(pan_file, pxn_file);

    }
        
    public static void writeTriangularStrip(BufferedWriter out, TriStrip tristrip) throws IOException
    {
        out.write("\t\t\t<tristrip count=\"" + tristrip.vertex.size() + "\">\t");
        for (int v = 0; v < tristrip.vertex.size(); v++) {
            // Subtract one because PANGU vertex indices are zero-based.
            out.write((tristrip.vertex.get(v)-1) + " ");
        }
        out.write("\t</tristrip>\n");
        out.flush();
    } 

}
