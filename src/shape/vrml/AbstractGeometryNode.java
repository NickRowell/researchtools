/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package shape.vrml;

import java.io.BufferedWriter;
import java.io.IOException;

/**
 *
 * @author nickrowell
 */
/**
 * Not a true Node: several different nodes can exist in the geometry field
 * of a Shape node, but there is no node named Geometry. We supply an abstract
 * class because we want to ensure that any geometry nodes can write out their
 * geometry information as a PXN mesh object.
 * 
 * @author nickrowell
 */
public abstract class AbstractGeometryNode extends VrmlNode
{
    Triple ambientColor;
    Triple diffuseColor;
    Triple specularColor;
    Triple emissiveColor;
    float shininess = 0.2f;
    float transparency = 0;

    // Objects of type Geometry must support writing a PXN <mesh> node
    abstract void writePXN(BufferedWriter out, int indent) throws IOException;
    
    public AbstractGeometryNode(String _name)
    {
        super(_name);
    }
    
    
    
    
    
    
    
}
