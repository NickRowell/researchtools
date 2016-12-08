package numeric.geom.dim2;

import numeric.geom.base.BaseVectorNd;

/**
 * Class represents vectors in two dimensions, and provides some methods unique to them.
 *
 * @author nrowell
 * @version $Id$
 */
public class Vector2d extends BaseVectorNd<Vector2d> {
	
	/**
	 * Main constructor.
	 * 
	 * @param components
	 * 	The vector components.
	 */
	public Vector2d(double... components) {
    	super(components);
		if(components.length!=2) {
			throw new IllegalArgumentException("Wrong number of components! Expected 2, found "+components.length);
		}
    }
	
	/**
	 * {@inheritDoc}
	 */
	@Override
    public Vector2d newInstance(double... components) {
    	return new Vector2d(components);
    }

    /**
     * Set x component
     * 
     * @param x 
     * 	x component
     */
    public final void setX(double x) {
    	this.components[0] = x;
    }
    
    /**
     * Set y component
     * 
     * @param y 
     * 	y component
     */
    public final void setY(double y) {
    	this.components[1] = y;
    }
    
    /**
     * Get x component
     * 
     * @return
     * 	The x component
     */
    public double getX() {
    	return this.components[0];
    }
    
    /**
     * Get y component
     * 
     * @return
     * 	The y component
     */
    public double getY() {
    	return this.components[1];
    }
    
}