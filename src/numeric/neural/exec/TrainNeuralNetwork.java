package numeric.neural.exec;

import java.io.File;
import java.io.FileFilter;
import java.util.LinkedList;
import java.util.List;

import numeric.neural.algo.NeuralNetwork;
import util.FileUtil;

/**
 * Processing some {@link Source}s to classify them using the {@link NeuralNetwork} classifier.
 */
public class TrainNeuralNetwork {
	
	/**
	 * Main application entry point.
	 * @param args
	 * 	The command line arguments (ignored)
	 */
	public static void main(String[] args) {
		
		// The directory containing all the files of {@link Source}s to process
		File sourceDirectory = new File("/home/nrowell/Projects/SummerProjects/GaiaSourceClassification/data/Source/TrainingSet");
		
		// Use a {@link FileFilter} to locate files containing {@link Source}s (they
		// have names starting 'Source_')
		FileFilter directoryFilter = new FileFilter() {
			public boolean accept(File file) {
				if(!file.isFile()) {
					return false;
				}
				return file.getName().startsWith("Source_");
			}
		};
		// Array of all files containing {@link Source}s
		File[] files = sourceDirectory.listFiles(directoryFilter);
		
		// Load the training set
		List<Source> sources = new LinkedList<>();
		
		// Process each file in turn
		for(File file : files) {
			// Load all the {@link Source}s from the file
			sources.addAll((List<Source>) FileUtil.deserialize(file));
		}
		
		System.out.println("Loaded "+sources.size()+" Sources from "+files.length+" Files");
		
		// Transform the training set into the form required for feeding to network
		
		// Trailing dimension has size equal to the number of network inputs
		double[][] inputs = new double[sources.size()][3];
		
		// Trailing dimension has size equal to the number of classes
		double[][] outputs = new double[sources.size()][6];
		
		// Load each {@link Source} into the training data arrays
		for(int s=0; s<sources.size(); s++) {
			
			Source source = sources.get(s);
			
			// Construct the network inputs
			inputs[s][0] = source.getFluxRatio();
			inputs[s][1] = source.getEigenvalues()[0];
			inputs[s][2] = source.getEigenvalues()[1];
			
			// Construct the desired network outputs, i.e. [0, 0, 1, 0, 0, 0]
			outputs[s][source.getType().ordinal()] = 1.0;
		}

		// Create a {@link NeuralNetwork} with random starting values
		NeuralNetwork nn = new NeuralNetwork(new int[]{4,4,6}, 3);
		
		// Number of random starting points to try for the network
		int nTrials = 1000;
		// Number of iterations to perform for each trial
		int nIter = 1000;
		
		nn.trainNetwork(inputs, outputs, nTrials, nIter);
		
		System.out.println("Lowest cost function = "+nn.getTotalNetworkError(inputs, outputs));
		
		double[] bestParams = nn.getNetworkParameters();
		
		// Print out the network parameters so they can be pasted into the source classifier
		System.out.print("double[] params = new double[]{"+bestParams[0]);
		
		for(int i=1; i<bestParams.length; i++) {
			System.out.print(", "+bestParams[i]);
		}
		System.out.println("};");
		
	}
}