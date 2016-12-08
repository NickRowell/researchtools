package util;

import java.util.LinkedList;
import java.util.List;

/**
 * Array utilities.
 *
 * @author nrowell
 * @version $Id$
 */
public class ArrayUtil {
	
	/**
	 * Tolerance used for testing equality of doubles.
	 */
	private static final double EPSILON = 1e-9;
	
	/**
	 * Gets the largest value in the array.
	 * 
	 * @param input
	 * 		The array.
	 * @return
	 * 		The (valid) element closest to positive infinity.
	 */
	public static double max(double[] input)
	{
		double max = -Double.MAX_VALUE;
		for(double value : input)
		{
			if(Double.isNaN(value) || Double.isInfinite(value))
				continue;
			
			if(value > max)
				max = value;
		}
		return max;
	}
	
    /**
     * Checks that all the entries in the array are positive or zero. If any NaN and
     * Infinity elements are present then the returned value is false.
     * @param array
     * 	The array to check
     * @return
     * 	True if all the entries are positive or zero, false otherwise. If any NaN and
     * Infinity elements are present then the returned value is false.
     */
    public static boolean checkNonNegative(double[] array) {
    	
        for(double value : array) {
        	// Check for NaN or Infinity
        	if(Double.isNaN(value) || Double.isInfinite(value)) {
        		return false;
        	}
        	// Check for negative value
            if(value<0) {
                return false;
            }
        }
        // Test passed
        return true;
    }
	
	/**
	 * Converts a {@link List} of {@link Double}s to an array.
	 * @param list
	 * 	The {@link List} to convert.
	 * @return
	 * 	An array containing the elements of the list.
	 */
	public static double[] toArray(List<Double> list) {
		double[] array = new double[list.size()];
    	for(int i=0; i<list.size(); i++) {
    		array[i] = list.get(i);
    	}
    	return array;
	}
	
	/**
	 * Convert the array of generic types into a List.
	 * @param array
	 * 	Array of generic types
	 * @return
	 * 	List containing all the elements in the array
	 */
	public static <T> List<T> toList(T[] array) {
		
		List<T> list = new LinkedList<>();
		
		for(T el : array) {
			list.add(el);
		}
		
		return list;
	}
	
	/**
	 * Extracts sequences of bytes from an array and returns them as a String
	 * @param bytes
	 * 	The byte array to extract bytes from
	 * @param start
	 * 	The starting index of the elements to extract (inclusive).
	 * @param end
	 * 	The ending index of the elements to extract (exclusive). If start=end then no bytes
	 * are extracted and the String has length zero.
	 * @return
	 * 	A String constructed from the elements of the byte array from bytes[start] to bytes[end-1].
	 */
	public static String getStringFromBytes(byte[] bytes, int start, int end) {
		byte[] sub = new byte[end - start];
		for(int i=start; i<end; i++) {
			sub[i-start] = bytes[i];
		}
		return new String(sub);
	}

    /**
     * Check that elements in the given array show a monotonic increase.
     * Equal values are not allowed.
     * @param x
     * 	input array to be tested.
     * @return  
     * 	true if x shows monotonic increase, false otherwise.
     */
    public static final boolean checkIncreasing(double[] x){
        for(int l=0; l<x.length-1; l++) {
            if(x[l] >= x[l+1]) {
                return false;
            }
        }
        return true;
    }

    /**
     * Check that elements in the given array show a monotonic decrease.
     * Equal values are not allowed.
     * @param x 
     * 	input array to be tested.
     * @return  
     * 	true if x shows monotonic decrease, false otherwise.
     */
    public static final boolean checkDecreasing(double[] x) {
        for(int l=0; l<x.length-1; l++) {
            if(x[l] <= x[l+1]) {
                return false;
            }
        }
        return true;
    }
    
    /**
     * Check that elements in the given array show a monotonic increase or
     * decrease.
     * @param x
     * 	input array to be tested.
     * @return
     * 	true if x is monotonic, false otherwise.
     */
    public static final boolean checkMonotonic(double[] x) {
    	return checkIncreasing(x) || checkDecreasing(x);
    }
    
    /**
     * Checks if any pair of consecutive bins are overlapping. Note that the bin
     * centres must be in order - increasing or decreasing is fine.
     * 
     * @param binCentres
     * 	The bin centres
     * @param binWidths
     * 	The bin widths
     * @return
     * 	True if no consecutive pair of bins overlaps, false otherwise.
     */
    public static boolean checkNonOverlappingBins(double[] binCentres, double[] binWidths) {
        
        // Note start at second bin, to
        // compare bins recursively against previous bin.
        for(int bin=1; bin<binCentres.length; bin++) {
        
            // Check non-overlap. Difference between centres of consecutive 
            // bins must be >= half the combined bin widths.
            double overlap = ((binWidths[bin] + binWidths[bin-1])/2.0) - 
                              (binCentres[bin] - binCentres[bin-1]);
            
            // Allow for small floating point errors when testing 
            // contiguous bins.
            if (overlap > EPSILON * binWidths[bin]) {
                return false;
            }
            
        }
        // Test passed
        return true;
    }
    
}