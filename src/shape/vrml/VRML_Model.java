/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package shape.vrml;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StreamTokenizer;
import java.util.LinkedList;
import java.util.List;

/**
 *
 * @author nickrowell
 */
class VRML_Model
{
    
    // Collection of 'defined' nodes: they are stored here whenever a DEF tag
    // is found in the VRML, and can be looked up later when a corresponding
    // USE tag is found.
    
    // The idea is we avoid reproducing nodes and have one single PXN file
    // for repeated structures.
    public List<VrmlNode> DEFs = new LinkedList<VrmlNode>();
    
    // Collection of top level nodes in the model. Often there will be just one.
    public List<VrmlNode> topLevelNodes = new LinkedList<VrmlNode>();
    
    // Name of model, extracted from vrml filename
    String name;
    

    
    
    public VRML_Model(File wrl) throws FileNotFoundException, IOException
    {
        // Extract name
        name = wrl.getName();
        // Remove file extension
        name = name.substring(0, name.lastIndexOf('.'));
        
        System.out.println("Model name = "+name);
        
        // Use StreamTokenizer to parse VRML
        FileReader file = new FileReader(wrl);
        StreamTokenizer streamTokenizer = new StreamTokenizer(file);
        
        // Indicate any special characters
        streamTokenizer.commentChar('#');
        
        // These would otherwise be unknown token types:
        streamTokenizer.wordChars('_', '_');   // Underscore is used in DEF names
        streamTokenizer.wordChars('{', '{');   // Open & closing curly braces
        streamTokenizer.wordChars('}', '}');
        streamTokenizer.wordChars('[', '[');   // Open & closing square brackets
        streamTokenizer.wordChars(']', ']');
        streamTokenizer.wordChars('(', '(');   // Sometimes used in node names
        streamTokenizer.wordChars(')', ')');
        
        // Treat commas as white space rather than tokens
        streamTokenizer.whitespaceChars(',', ',');
        
        // Can use command sed 's/(/_/g;s/)/_/g' Rover_Assembly_vrml.wrl  > tmp
        // to remove other symbols from e.g. badly chosen node names.
        
        // Line ends are not significant (only opening/closing braces)
        streamTokenizer.eolIsSignificant(false);
        
        // Loop over all tokens
        while(streamTokenizer.nextToken() != StreamTokenizer.TT_EOF)
        {
            // Construct the next top level node
            VrmlNode top = VrmlNode.parseNode(streamTokenizer, DEFs);
            topLevelNodes.add(top);
        }
        
        // Now apply transformation to all top level nodes so that the first
        // one is positioned at the origin with no rotation.
        
        
        
        
    }
    
    /**
     * Writes VRML file out to disk as a PXN by matching corresponding node
     * types. Starts a new folder for each grouping node, and places each new
     * shape node into a new file.
     * @param topLevel 
     */
    void writePXN(File topLevel) throws IOException
    {
        // Ideas for writing out to PXN file:
        //  - Very top level contains the topmost .pxn file.
        //  - Iterate over top level nodes. Each defines a new object, contained
        //    in a new folder at the top level.
        //  - Grouping nodes create a new folder to hold each child
        //  - They are represented by a .pxn file that <include>'s each child node
        //  - In the case of Transform node, we include the <transform> node at the top level
        //  - Shape nodes always include an <appearance> node with defaults set
        //  - Need to write IndexedFaceSet out to <mesh> node.
        //  - Should try and skip grouping nodes with single child - just 
        //    write the child.
        //  - Keep track of level in heirarchy so that we can properly indent files?
        //  - Check use of Euler angles for orientation of child node
        
        
        // Top level PXN file contains no geometry itself, only <include> nodes
        // pointing to pxn files at lower levels in model tree.
        File topFile = new File(topLevel, name+".pxn");
        BufferedWriter out = new BufferedWriter(new FileWriter(topFile));
        PXN_Utils.writePXNHeader(out);

        // Now write out nodes
        for(VrmlNode node : topLevelNodes)
        {
            VrmlNode.writeNode(out, topLevel, topLevel, node, 0);
        }
        
        PXN_Utils.writePXNFooter(out);
        
    }
    


    

    
    
}