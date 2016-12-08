package numeric.integration;

/**
 * Defines the interface of an integrable function.
 *
 *
 * @author nrowell
 * @version $Id$
 */
public interface IntegrableFunction {

	/**
	 * This method implements the function that is to be integrated.
	 * @param x	The independent variable x.
	 * @return	The value of the function at x, i.e. y(x).
	 */
	public double evaluate(double x);
}
