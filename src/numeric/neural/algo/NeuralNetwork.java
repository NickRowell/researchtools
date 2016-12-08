package numeric.neural.algo;

/**
 * Class represents a network of interconnected neurons.
 * 
 * Usage:
 * 1) The application will determine the number of inputs and outputs of the network.
 *    The main design factor to choose is the number of neuron layers and the number of
 *    neurons in each layer. Try out different values.
 * 2) Compile a training set that consists of example inputs and desired outputs.
 * 3) Create a network: the number of inputs is specified in the constructor, the number
 *    of outputs is determined by the number of neurons in the final layer, so make sure
 *    that this matches the number of outputs required by the application.
 * 4) Train the network: feed it with the training data using the 
 *    {@link NeuralNetwork#trainNetwork(double[][], double[][], int, int)} method.
 *    On exit, the network is trained and ready to go.
 * 5) In classification problems, the network can easily be put in a wrapper class that
 *    converts the outputs to a particular class.
 * 
 * @author nrowell
 * @version $Id$
 */
public class NeuralNetwork {
	
	/**
	 * Two-dimensional array containing all the neurons in the network split into layers.
	 * The leading index specifies the layer; the trailing index specifies the neuron in
	 * the layer, i.e.
	 * 
	 *  - neurons.length       = the number of layers
	 *  - neurons[0].length    = the number of neurons in the first layer
	 *  - neurons[1][2]        = the third neuron in the second layer etc
	 *  - neurons[numLayers-1] = array of neurons in the final layer
	 *  
	 */
	SigmoidNeuron[][] neurons;
	
	/**
	 * For a given inputs to the network, this array stores the outputs of each layer of
	 * neurons. This is useful when applying the back propagation algorithm, i.e. for a
	 * given input to the N-layer network;
	 * 
	 *  - a[0]   -> contains the outputs of the neurons in the first layer
	 *  - a[1]   -> contains the outputs of the neurons in the second layer
	 *    ...
	 *  - a[N-1] -> contains the outputs of the neurons in the final layer (the total network outputs)
	 * 
	 */
	double[][] a;
	
	/**
	 * The number of layers in the network.
	 */
	int numLayers;
	
	/**
	 * The number of inputs to the network (equals the number of inputs to the neurons in the
	 * first layer).
	 */
	int numInputs;
	
	/**
	 * The number of outputs from the network (equals the number of neurons in the final layer
	 * of the network).
	 */
	int numOutputs;
	
	/**
	 * Total number of parameters of the network (equals the total number of parameters for all
	 * the neurons).
	 */
	int numParams;
	
	/**
	 * Constructor for a {@link NeuralNetwork}. The network architecture is fully connected, i.e.
	 * each neuron receives as inputs the outputs from all the neurons in the preceeding layer,
	 * and in the first layer each neuron receives all the inputs. The number of outputs from
	 * the network is equal to the number of neurons in he final layer.
	 * 
	 * The parameters of the individual {@link SigmoidNeuron}s are set to random values.
	 * 
	 * @param neuronsPerLayer
	 * 	Number of neurons in each layer (can be different in each layer)
	 * @param numInputs
	 * 	Number of inputs to the first layer
	 */
	public NeuralNetwork(int[] neuronsPerLayer, int numInputs) {

		this.numLayers = neuronsPerLayer.length;
		this.numInputs = numInputs;
		this.numOutputs = neuronsPerLayer[this.numLayers-1];
		this.numParams = 0;
		
		// Create the neurons array; number of layers in leading dimension...
		neurons = new SigmoidNeuron[this.numLayers][];
		a = new double[this.numLayers][];
		
		// ...then neurons per layer in the trailing dimension
		for(int i=0; i<this.numLayers; i++) {
			neurons[i] = new SigmoidNeuron[neuronsPerLayer[i]];
			a[i] = new double[neuronsPerLayer[i]];
		}
		
		// Now we create the individual neurons and configure them according to the network layout
		for(int i=0; i<this.numLayers; i++) {
			
			// Number of inputs to the neurons in this layer: in the first layer this is equal
			// to the number specified in the argument, and in the following layers it's equal
			// to the number of neurons in the preceeding layer (assuming that each neuron is
			// connected to all the neurons in the preceeding layer).
			int n = (i==0) ? this.numInputs : neuronsPerLayer[i-1];
			
			for(int j=0; j<neuronsPerLayer[i]; j++) {
				// Create neuron j in layer i
				neurons[i][j] = new SigmoidNeuron(n);
				
				// Add up the number of parameters for the entire network. For each neuron this is
				// equal to the number of inputs plus one.
				this.numParams += n+1;
			}
		}
	}
	
	/**
	 * Constructor for a {@link NeuralNetwork}. The network architecture is fully connected, i.e.
	 * each neuron receives as inputs the outputs from all the neurons in the preceeding layer,
	 * and in the first layer each neuron receives all the inputs. The parameters of the individual
	 * {@link SigmoidNeuron}s are read from the given array.
	 * 
	 * @param neuronsPerLayer
	 * 	Number of neurons in each layer (can be different in each layer)
	 * @param numInputs
	 * 	Number of inputs to the first layer
	 * @param p
	 * 	The parameters of each constituent neuron
	 */
	public NeuralNetwork(int[] neuronsPerLayer, int numInputs, double[] p) {
		
		// Initialise the network
		this(neuronsPerLayer, numInputs);
		
		setNetworkParameters(p);
	}
	
	/**
	 * Resets the parameters of all the {@link SigmoidNeuron}s in the network to random
	 * starting values.
	 */
	public void reset() {
		for(int i=0; i<this.numLayers; i++) {
			for(int j=0; j<neurons[i].length; j++) {
				neurons[i][j].initialise();
			}
		}
	}
	
	/**
	 * Extract the parameters of the {@link NeuralNetwork}, i.e. the weights and biases of all
	 * the constituent neurons. The ordering of the elements is such that the array can be
	 * passed into the constructor {@link NeuralNetwork#NeuralNetwork(int, int[], int, double[])}
	 * to recreate it.
	 * 
	 * @return
	 * 	The parameters of the {@link NeuralNetwork}, i.e. the weights and biases of all
	 * the constituent neurons.
	 */
	public double[] getNetworkParameters() {
		
		double[] networkParams = new double[numParams];
		
		// Index into the parameters array of the parameters for the current neuron
		int paramIdx = 0;
		
		// Extract the parameters from the neurons in each layer
		for(int l=0; l<neurons.length; l++) {
			
			// Number of parameters per neuron in this layer.
			int numParamsPerNeuron = (l==0) ? this.numInputs + 1 : neurons[l-1].length + 1;
			
			// Loop over each neuron in this layer
			for(int j=0; j<neurons[l].length; j++) {
				
				// Extract the parameters of this neuron
				double[] neuronParams = neurons[l][j].getParameters();
				
				// The relevant elements are [paramIdx] to [paramIdx+numParamsPerNeuron-1]
				for(int p=0; p<numParamsPerNeuron; p++) {
					networkParams[paramIdx++] = neuronParams[p];
				}
			}
		}
		return networkParams;
	}
	
	/**
	 * Set the parameters of the {@link NeuralNetwork}, i.e. the weights and biases of all
	 * the constituent neurons.
	 * 
	 * @param params
	 * 	The parameters of the {@link NeuralNetwork}, i.e. the weights and biases of all
	 * the constituent neurons.
	 */
	public void setNetworkParameters(double[] params) {

		// Check we've got the right number of parameters for the network architecture
		if(params.length != numParams) {
			throw new IllegalArgumentException("Expected "+numParams+" parameters,"
					+ "found "+params.length+"!");
		}
		
		// Index into the parameters array of the parameters for the current neuron
		int paramIdx = 0;
		
		// Extract the parameters from the neurons in each layer
		for(int l=0; l<neurons.length; l++) {
			
			// Number of parameters per neuron in this layer.
			int numParamsPerNeuron = (l==0) ? this.numInputs + 1 : neurons[l-1].length + 1;
			
			double[] newParams = new double[numParamsPerNeuron];
			
			// Loop over each neuron in this layer
			for(int j=0; j<neurons[l].length; j++) {
				
				// Extract the new parameters for this neuron
				for(int p=0; p<numParamsPerNeuron; p++) {
					newParams[p] = params[paramIdx++];
				}
				// Set the new parameters
				neurons[l][j].setParameters(newParams);
			}
		}
	}
	
	/**
	 * Computes the output values for the {@link NeuralNetwork} for the given inputs. On exit,
	 * the elements of the array {@link NeuralNetwork#a} contain the outputs of each layer of
	 * neurons in the network.
	 * 
	 * @param inputs
	 * 	The inputs to the network.
	 * @return
	 * 	The output from each neuron in the final layer.
	 */
	public double[] getNetworkOutput(double[] inputs) {
		
		// Sanity check on inputs
		if(inputs.length != numInputs) {
			throw new IllegalArgumentException("Number of provided inputs ("+inputs.length+") "
					+ "doesn't match the number set on construction of the network ("+numInputs+")!");
		}
		
		// Compute the outputs of each layer of neurons in turn. This variable will be used
		// to cache the outputs of each successive layer for input to the next layer; after
		// the loop is complete it will contain the outputs from the final layer.
		for(int i=0; i<numLayers; i++) {
			
			// Get the inputs to this layer - either the original inputs (for layer one)
			// or the outputs from the previous layer.
			double[] layerInputs = (i==0) ? inputs : a[i-1];
			
			// Compute the output of each neuron and store it in the array
			for(int j=0; j<neurons[i].length; j++) {
				// Retrieve neuron j in layer i
				SigmoidNeuron neuron = neurons[i][j];
				// Compute it's output and store in the array
				a[i][j] = neuron.computeOutput(layerInputs);
			}
			
		}
		
		// Return the outputs from the final layer in the network
		return a[numLayers-1];
	}
	
	/**
	 * Single iteration implementation of the Back Propagation algorithm.
	 * 
	 * @param trainingInputs
	 * 	The training set inputs, stored in a two-dimension array. The leading dimension
	 * loops over each object, the trailing dimension loops over the inputs for each object.
	 * 
	 * @param desiredOutputs
	 * 	The desired outputs for each object in the training set, stored in a two-dimension
	 * array. The leading dimension loops over each object, the trailing dimension loops
	 * over the desired output for each object.
	 */
	public void backPropagation(double[][] trainingInputs, double[][] desiredOutputs) {
		
		// The Back Propagation algorithm
		
		// The partial derivatives of the cost function with respect to each parameter, averaged
		// over all of the training examples.
		double[] dC_by_dp = new double[numParams];
		
		// Accumulate the deltas for each parameter in the network across all the training examples
		for(int s=0; s<trainingInputs.length; s++) {

			// Get the network outputs for this training example
			getNetworkOutput(trainingInputs[s]);
			
			// Compute the 'error' vector delta for each layer of the network. This is the
			// gradient of the cost function wrt the weighted input of each neuron in the layer.
			double[][] delta = new double[numLayers][];
			
			// Index of the final network layer
			int l = numLayers-1;
			
			// First, compute the error in the output layer
			delta[l] = new double[neurons[l].length];
			
			// Retrieve the inputs to the final layer; handle single layer cases
			double[] layerInputs = (l==0) ? trainingInputs[s] : a[l-1];
			
			for(int j=0; j<neurons[l].length; j++) {
				
				// Derivative of output of neuron j wrt the weighted input zj
				double s_prime_zj = neurons[l][j].computeDerivative(layerInputs);
				
				// Derivative of the cost function wrt the output of neuron j in the final layer
				double dC = 2 * (a[l][j] - desiredOutputs[s][j]);
				
				// Component of output error in final layer
				delta[l][j] = dC * s_prime_zj;
			}
			
			// Now backpropagate the error through the remaining network layers
			for(l=numLayers-2; l>=0; l--) {
				
				delta[l] = new double[neurons[l].length];
						
				// Retrieve the inputs to this layer of neurons; handle the first layer
				layerInputs = (l==0) ? trainingInputs[s] : a[l-1];
				
				// Retrieve the error vector for the l+1 layer
				double[] d = delta[l+1];
				
				// Loop over the neurons in layer l
				for(int j=0; j<neurons[l].length; j++) {

					// Derivative of output of neuron j wrt the weighted input zj
					double s_prime_zj = neurons[l][j].computeDerivative(layerInputs);
					
					// Compute the matrix multiplication of the weights and errors for layer l+1
					double wd = 0;
					for(int i=0; i<d.length; i++) {
						wd += neurons[l+1][i].w[j] * d[i];
					}
					
					// Component of output error in layer l
					delta[l][j] = wd * s_prime_zj;
				}
			}
			
			// Delta now contains all the information to compute the gradient of the cost function
			// with respect to any parameter in the network. Use this to apply updates to the
			// parameters and train the network.
			int paramIdx = 0;
			
			for(l=0; l<numLayers; l++) {
				
				// Retrieve the inputs to this layer of neurons; handle the first layer
				layerInputs = (l==0) ? trainingInputs[s] : a[l-1];
				
				for(int j=0; j<neurons[l].length; j++) {
					
					// Compute partial derivatives for the weights of this neuron
					for(int k=0; k<neurons[l][j].w.length; k++) {
						dC_by_dp[paramIdx++] += delta[l][j] * layerInputs[k];
					}
					
					// Compute partial derivatives for the bias of this neuron
					dC_by_dp[paramIdx++] += delta[l][j];
				}
			}
			
		}
		
		// Now apply parameter updates
		double[] networkParams = this.getNetworkParameters();
		for(int p=0; p<numParams; p++) {
			networkParams[p] +=  -dC_by_dp[p] / trainingInputs.length;
		}
		this.setNetworkParameters(networkParams);
	}
	
	/**
	 * Computes the cost function value for the training examples.
	 * 
	 * @param inputs
	 * 	The training set inputs, stored in a two-dimension array. The leading dimension
	 * loops over each object, the trailing dimension loops over the inputs for each object.
	 * @param outputs
	 * 	The desired outputs for each object in the training set, stored in a two-dimension
	 * array. The leading dimension loops over each object, the trailing dimension loops
	 * over the desired output for each object.
	 * 
	 * @return
	 * 	The cost function.
	 */
	public double getTotalNetworkError(double[][] inputs, double[][] outputs) {

		// Some consistency checks
		if(inputs.length != outputs.length) {
			throw new IllegalArgumentException("Number of training inputs ("+inputs.length+") "
					+ "doesn't match the number of outputs ("+outputs.length+")!");
		}
		for(int i=0; i<inputs.length; i++) {
			if(inputs[i].length != numInputs) {
				throw new IllegalArgumentException("Training input "+i+": expected "+numInputs
						+ " inputs, found "+inputs[i].length+"!");
			}
			if(outputs[i].length != numOutputs) {
				throw new IllegalArgumentException("Training output "+i+": expected "+numOutputs
						+ " outputs, found "+outputs[i].length+"!");
			}
		}

		double e = 0.0;
		
		// Sum the squared error for each training example
		for(int s=0; s<inputs.length; s++) {
			
			// Get the network outputs for this training example
			double[] o = getNetworkOutput(inputs[s]);
			
			// Sum the error in quadrature
			for(int i=0; i<numOutputs; i++) {
				e += (o[i] - outputs[s][i])*(o[i] - outputs[s][i]);
			}
		}
		
		// Average
		e /= (2*inputs.length);
		
		return e;
	}
	
	/**
	 * Trains the network using the given training inputs and desired outputs. Uses many
	 * random starting points to improve the chances of finding the global minimum, and
	 * makes many iterations within each realization in order to converge on the best
	 * parameter values.
	 * 
	 * @param inputs
	 * 	The training set inputs, stored in a two-dimension array. The leading dimension
	 * loops over each object, the trailing dimension loops over the inputs for each object.
	 * @param outputs
	 * 	The desired outputs for each object in the training set, stored in a two-dimension
	 * array. The leading dimension loops over each object, the trailing dimension loops
	 * over the desired output for each object.
	 * @param nTrials
	 * 	Number of times to train the network from a new random starting point; use a larger
	 * value to have a better chance of finding the global minimum.
	 * @param nIter
	 * 	Number of iterations to make from each random starting point; use a larger value to
	 * improve convergence for each trial.
	 */
	public void trainNetwork(double[][] inputs, double[][] outputs, int nTrials, int nIter) {

		// Variable to store the network parameters that give the best fit
		double[] bestParams = null;
		
		// Lowest cost function value
		double minCost = Double.MAX_VALUE;
		
		for(int t=0; t<nTrials; t++) {
		
			// Reset the network parameters to random starting values
			this.reset();
			
			// Train the network using back propagation
			for(int i=0; i<nIter; i++) {
				this.backPropagation(inputs, outputs);
			}
			
			double finalError = this.getTotalNetworkError(inputs, outputs);
			
			if(finalError < minCost) {
				minCost = finalError;
				bestParams = this.getNetworkParameters();
			}
		}
		
		// Set the network parameters to the best found in all the trials
		this.setNetworkParameters(bestParams);
	}
	
}