package numeric.geom.dimN;

import numeric.geom.base.BaseVectorNd;

/**
 * Class represents general multidimensional vectors.
 *
 * @author nrowell
 * @version $Id$
 * @param <T>
 * 	The type of the class that extends this one.
 */
public class VectorNd extends BaseVectorNd<VectorNd> {

	/**
	 * Main constructor.
	 * 
	 * @param components
	 * 	The vector components.
	 */
	public VectorNd(double... components) {
		super(components);
	}
	
	/**
	 * {@inheritDoc}
	 */
	@Override
	public VectorNd newInstance(double... components) {
		return new VectorNd(components);
	}
	
}