/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package shape.vrml;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.StreamTokenizer;
import java.util.LinkedList;
import java.util.List;

/**
 *
 * @author nickrowell
 */
public class VrmlIndexedFaceSetNode extends AbstractGeometryNode
{
    boolean solid;
    
    // We support only triangle faces. -1 is a marker between consecutive
    // triangular faces in list.
    List <Triple<Integer>> coordIndices = new LinkedList<Triple<Integer>>();
    VrmlCoordinateNode coord;
    
    VrmlNormalNode normals;
    
    boolean normalsDefined = false;
    
    public VrmlIndexedFaceSetNode(String _name)
    {
        super(_name);
    }
    
    void writePXN(BufferedWriter out, int indent) throws IOException
    {
        
        String indentation = "";
        for(int i=0; i<indent; i++) 
            indentation += indentChar;
        
        out.write(indentation+"<mesh>\n");
        
        if(coord.points.size() != normals.normals.size())
            System.err.println("IndexedFaceSet "+name+" has "+coord.points.size()+
                 " vertices and "+normals.normals.size()+" normals!");
        
        // Writes <positions> element
        coord.writePXN(out, indent+1);
        
        // Write <normals> element
        normals.writePXN(out, indent+1);
        
//        // Optional
//        out.write(indentation+"\t<colours>\n");
//        out.write(indentation+"\t</colours>\n");
//        // Optional
//        out.write(indentation+"\t<texcoords>\n");
//        out.write(indentation+"\t</texcoords>\n");
//        // Optional
//        out.write(indentation+"\t<texgen>\n");
//        out.write(indentation+"\t</texgen>\n");
//        // Optional
//        out.write(indentation+"\t<delete>\n");
//        out.write(indentation+"\t</delete>\n");
        
        // Could do something smarter here to connect individual triangles into
        // longer strips?
        out.write(indentation+indentChar+"<tristrips>\n");
        for(Triple<Integer> strip : coordIndices)
            out.write(indentation+indentChar+indentChar+"<tristrip>" +strip.x+ " " +strip.y+ " " +strip.z+ "</tristrip>\n");
        out.write(indentation+indentChar+"</tristrips>\n");
        
        out.write(indentation+"</mesh>\n");
        
    }
    
    
    // Stream tokenizer is positioned on opening curly brace on entry.
    // Stream tokenizer is positioned just after closing curly brace on exit.
    static VrmlIndexedFaceSetNode parse(StreamTokenizer st, String name, List<VrmlNode> DEFs) throws IOException
    {
        VrmlIndexedFaceSetNode indexedFaceSet = new VrmlIndexedFaceSetNode(name); 
        
        // move tokenizer to next token: either first field, or closing brace
        // for empty objects.
        st.nextToken();
        
        // Loop over tokens until closing curly brace is reached
        while(!VrmlNode.foundClosingCurlyBrace(st))
        {
            // Read field:
            if(st.ttype == StreamTokenizer.TT_WORD)
            {
                String field = st.sval;
                
                if(field.equalsIgnoreCase("solid"))
                {
                    // Move on to boolean value
                    st.nextToken();
                    indexedFaceSet.solid = parseBoolean(st);
                }
                else if(field.equalsIgnoreCase("coordIndex"))
                {
                    // Move on to opening [
                    st.nextToken();
                    parseCoordIndex(st, indexedFaceSet.coordIndices);
                }
                else if(field.equalsIgnoreCase("coord"))
                {
                    // Move on to opening [
                    st.nextToken();
                    indexedFaceSet.coord = (VrmlCoordinateNode)parseNode(st, DEFs);
                }
                else if(field.equalsIgnoreCase("normal"))
                {
                    indexedFaceSet.normalsDefined = true;
                    
                    // Move on to opening [
                    st.nextToken();
                    indexedFaceSet.normals = (VrmlNormalNode)parseNode(st, DEFs);
                }
                
                else
                    throw new RuntimeException("Unrecognized IndexedFaceSet field: "+field
                            +", "+st.toString());
                
            }
            
        }
        
        // Position stream tokenizer just after closing curly brace
        st.nextToken();
        
        return indexedFaceSet;
    }
    
    
    // On entry, stream tokenizer is positioned on opening [
    // Coordinate indices are triple seperated by -1 and ending with ]
    // On exit, stream tokenizer is position after closing ]
    static void parseCoordIndex(StreamTokenizer st, List<Triple<Integer>> coordIndices) throws IOException
    {
        // Move to start of first point, or closing ] for empty coordIndex field
        st.nextToken();
        
        while(!VrmlNode.foundClosingSquareBrace(st))
        {
            coordIndices.add(VrmlNode.parseTripleInt(st));
            
            // Sets of coordIndex values that define single faces in set are
            // seperated by the value -1.
            switch(st.ttype)
            {
                case StreamTokenizer.TT_NUMBER:
                {
                    // Check for -1
                    if(Math.abs(st.nval+1)<1e-5)
                    {
                        // Move to next token (start of next face)
                        st.nextToken();
                        break;
                    }
                    else
                        // We've encountered a face with more than 3 or less
                        // than 3 vertices
                        throw new RuntimeException("Badly formed coordIndex set"
                                + " at "+st.toString()+": number of vertices"
                                + " != 3");
                }
                case StreamTokenizer.TT_WORD:
                {
                    // Most likely a closing square brace. Do nothing and let
                    // next iteration of loop detect it and exit.
                }
                    
            
            }
            
            
        }
        // Position stream tokenizer just after closing square brace
        st.nextToken();
    } 
    
    
    
    
}
