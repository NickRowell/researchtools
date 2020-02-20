package numeric.functions;

/**
 * This interface represents a general mathematical function that takes an arbitrary number of
 * inputs and maps them to an arbitrary number of outputs.
 *
 *
 * @author nrowell
 * @version $Id$
 */
public abstract class Function {
	
	/**
	 * Dimensionality of function input.
	 */
	int nInputs;
	
	/**
	 * Dimensionality of function output.
	 */
	int nOutputs;
	
	
	/**
	 * Main constructor for a {@link Function}.
	 * @param nInputs
	 * 	Dimensionality of function input.
	 * @param nOutputs
	 * 	Dimensionality of function output.
	 */
	public Function(int nInputs, int nOutputs) {
		this.nInputs = nInputs;
		this.nOutputs = nOutputs;
	}
	
	/**
	 * Get the dimensionality of function input.
	 * @return
	 * 	The dimensionality of function input.
	 */
	public int getInputs() {
		return nInputs;
	}
	
	/**
	 * Get the dimensionality of function output.
	 * @return
	 *  The dimensionality of function output.
	 */
	public int getOutputs() {
		return nOutputs;
	}
	
	/**
	 * Compute the value of the function at the given coordinate.
	 * 
	 * @param inputs
	 * 	The function input.
	 * @return
	 * 	The function output.
	 */
	public abstract double[] compute(double[] inputs);
	
}
