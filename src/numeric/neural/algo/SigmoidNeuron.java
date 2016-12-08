package numeric.neural.algo;

/**
 * Class represents a single sigmoid neuron in a neural network.
 * 
 * @author nrowell
 * @version $Id$
 */
public class SigmoidNeuron {
	
	/**
	 * The weights applied to the neuron inputs.
	 */
	public double[] w;
	
	/**
	 * The bias parameter for the neuron.
	 */
	public double b;
	
	/**
	 * Main constructor for the {@link SigmoidNeuron}.
	 * @param p
	 * 	The parameters of the {@link SigmoidNeuron}; they are interpreted as follows:
	 * 
	 *  p[0]   -> weight 0
	 *  p[1]   -> weight 1
	 *   .
	 *   .
	 *  p[N-1] -> weight N-1
	 *  p[N]   -> bias
	 */
	public SigmoidNeuron(double[] p) {
		
		// Sanity check on parameters:
		if(p.length < 2) {
			throw new IllegalArgumentException("Need at least two parameters for a "
					+ "neuron! Found "+p.length);
		}
		
		w = new double[p.length-1];
		
		// Load the parameters
		setParameters(p);
	}
	
	/**
	 * Constructor for the {@link SigmoidNeuron} that sets the number of inputs
	 * but not specific values for their weights or the bias. All weights and the
	 * bias parameter will be drawn randomly from [-1:1].
	 * @param N
	 * 	The number of inputs.
	 */
	public SigmoidNeuron(int N) {
		this.w = new double[N];
		initialise();
	}
	
	/**
	 * Sets the values of the neuron weights and bias.
	 * 
	 * @param p
	 * 	The parameters of the {@link SigmoidNeuron}; they are interpreted as follows:
	 * 
	 *  p[0]   -> weight 0
	 *  p[1]   -> weight 1
	 *   .
	 *   .
	 *  p[N-1] -> weight N-1
	 *  p[N]   -> bias
	 */
	public void setParameters(double[] p) {

		// Sanity check on parameters:
		if(p.length != w.length + 1) {
			throw new IllegalArgumentException("Expected "+(w.length + 1)+" parameters,"
					+ "found "+p.length+"!");
		}
		
		// Load the neuron weights
		for(int i=0; i<p.length-1; i++) {
			w[i] = p[i];
		}
		
		// Load the neuron bias
		b = p[p.length-1];
	}
	
	/**
	 * Get the parameters of the neuron.
	 * @return
	 * The parameters of the {@link SigmoidNeuron}; they are interpreted as follows:
	 * 
	 *  p[0]   -> weight 0
	 *  p[1]   -> weight 1
	 *   .
	 *   .
	 *  p[N-1] -> weight N-1
	 *  p[N]   -> bias
	 */
	public double[] getParameters() {
		double[] p = new double[w.length + 1];
		for(int i=0; i<w.length; i++) {
			p[i] = w[i];
		}
		p[w.length] = b;
		return p;
	}
	
	/**
	 * Assigns random values to the input weights and threshold. This is done
	 * prior to training the network.
	 */
	public void initialise() {
		for(int i=0; i<w.length; i++) {
			w[i] = 2.0 * Math.random() - 1.0;
		}
		b = 2.0 * Math.random() - 1.0;
	}

	/**
	 * Computes the weighted inputs of the neuron.
	 * 
	 * @param x
	 * 	The neuron inputs.
	 * @return
	 * 	The weighted input z.
	 */
	public double getZ(double[] x) {

		// Sanity check on inputs...
		if(w.length != x.length) {
			throw new IllegalArgumentException("Number of inputs ("+x.length+") "
					+ "doesn't match the number of weights ("+w.length+")!");
		}
		
		// Dot product of inputs and neuron weights
		double z = 0;
		for(int i=0; i<w.length; i++) {
			z += w[i] * x[i];
		}
		
		// Apply the neuron bias
		z += b;
		
		return z;
	}
	
	/**
	 * Computes the output of the {@link SigmoidNeuron} from the given inputs.
	 * 
	 * @param x
	 * 	The neuron inputs
	 * @return
	 * 	The output of the {@link SigmoidNeuron}.
	 */
	public double computeOutput(double[] x) {
		return computeOutput(getZ(x));
	}
	
	/**
	 * Compute the output of the neuron given the weighted input z.
	 * 
	 * @param z
	 * 	The weighted input.
	 * @return
	 * 	The neuron output.
	 */
	public double computeOutput(double z) {
		// Compute the activation - the sigmoid function
		return 1.0 / (1.0 + Math.exp(-z));
	}

	/**
	 * Compute the derivative of the neuron output with respect to the weighted input.
	 * 
	 * @param x
	 * 	The neuron inputs.
	 * @return
	 * 	The derivative of the neuron output with respect to the weighted input.
	 */
	public double computeDerivative(double[] x) {
		return computeDerivative(getZ(x));
	}
	
	/**
	 * Compute the derivative of the neuron output with respect to the weighted input.
	 * 
	 * @param z
	 * 	The weighted input.
	 * @return
	 * 	The derivative of the neuron output with respect to the weighted input.
	 */
	public double computeDerivative(double z) {
		double sz = computeOutput(z);
		return sz * (1.0 - sz);
	}
	
}