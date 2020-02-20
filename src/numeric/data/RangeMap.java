package numeric.data;

import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map.Entry;
import java.util.TreeMap;

import util.ArrayUtil;

/**
 * Class provides a useful data structure that maps numbers to a corresponding range, and then to a value. This is
 * useful for binning objects according to a particular property.
 * 
 * Ranges are non-overlapping.
 *
 * @author nrowell
 * @version $Id$
 */
public class RangeMap<T> {
	
	/**
	 * Array of the {@link Range}s used to populate the maps; useful for indexing.
	 */
	Range[] ranges;
	
	/**
	 * Used to map continuous values to the corresponding {@link Range} that they fall in, or
	 * null if they aren't contained by any {@link Range} (i.e. lower or higher than the full set
	 * of {@link Range}s, or in a hole between consecutive non-contiguous {@link Range}s).
	 */
	TreeMap<Double, Range> binMap;
	
	/**
	 * Maps {@link Range} to a corresponding {@link Collection<T>}.
	 */
	HashMap<Range, Collection<T>> objectsMap;
	
	/**
	 * Main constructor.
	 * @param binCentres
	 * 		Centre of each {@link Range}
	 * @param binWidths
	 * 		Width of each {@link Range}
	 * @throws IllegalArgumentException
	 * 		If the arrays are different lengths, non-ascending or if they define
	 * overlapping ranges (which are not allowed).
	 */
	public RangeMap(double[] binCentres, double[] binWidths) {
		
        // Sanity checks on array sizes
		if(binCentres.length != binWidths.length) {
			throw new IllegalArgumentException("RangeMap: number of bin centres ("+binCentres.length+")"
					+ " does not equal the number of bin widths ("+binWidths.length+")!");
		}
		if(!ArrayUtil.checkIncreasing(binCentres)) {
			throw new IllegalArgumentException("RangeMap: Bin centres not increasing!");
		}
        if(!ArrayUtil.checkNonOverlappingBins(binCentres, binWidths)) {
            throw new IllegalArgumentException("RangeMap: Illegal bin centres/sizes.");
        }
        
        init(binCentres, binWidths);
	}
	
	/**
	 * Alternative constructor that creates a set of uniformly distributed bins based on the
	 * range and bin size.
	 * @param min
	 * 	Lower boundary
	 * @param max
	 * 	Upper boundary
	 * @param step
	 * 	Step size (bin width). If there's not a whole number of bins within the range, then the
	 * final bin is trimmed to fit.
	 */
	public RangeMap(double min, double max, double step) {
		
		// Sanity checks
		if(min >= max) {
			throw new IllegalArgumentException("RangeMap: lower limit ("+min+") must be less than upper limit ("+max+")!");
		}
		if(step <= 0.0) {
			throw new IllegalArgumentException("RangeMap: step size ("+step+") must be positive!");
		}
		
		int nBins = (int)Math.ceil((max - min)/step);
		
		double[] binCentres = new double[nBins];
		double[] binWidths = new double[nBins];
		
		for(int b=0; b<nBins; b++) {
			
			// Bin lower edge
			double lower = min + b*step;
			// Bin upper edge; fix upper edge of final bin to upper boundary
			double upper = (b==nBins-1) ? max : min + (b+1)*step;
			
			binCentres[b] = (lower + upper)/2.0;
			binWidths[b] = (upper - lower);
		}

        init(binCentres, binWidths);
	}
	
	/**
	 * Initialise the data structure using the (checked) bin centres and widths.
	 * 
	 * @param binCentres
	 * 		Centre of each {@link Range}
	 * @param binWidths
	 * 		Width of each {@link Range}
	 */
	private void init(double[] binCentres, double[] binWidths) {
		
        ranges = new Range[binCentres.length];
        binMap = new TreeMap<Double, Range>();
        objectsMap = new HashMap<Range, Collection<T>>();
        
        // Create the {@link Range}s and enter them into the map.
        for(int i=0; i<binCentres.length; i++) {
        	double lower = binCentres[i] - binWidths[i]/2.0;
        	double upper = binCentres[i] + binWidths[i]/2.0;
        	
        	Range range = new Range(lower, upper);
        	
        	ranges[i] = range;
        	binMap.put(lower, range);
        	objectsMap.put(range, new LinkedList<T>());
        }
	}
	
	/**
	 * Retrieve the {@link Collection<T>} corresponding to the {@link Range} containing
	 * the given key, or null if there is no such {@link Range}.
	 * @param key
	 * 	The key.
	 * @return
	 * 	The {@link Collection<T>} corresponding to the {@link Range} containing
	 * the given key, or null if there is no such {@link Range}.
	 */
	public Collection<T> get(double key) {
		
		// Find the {@link Range}
		Entry<Double, Range> entry = binMap.floorEntry(key);
		
		if(entry==null) {
			// The key is lower than the lowest {@link Range} in the map
			return null;
		}
		
		Range range = entry.getValue();
		
		if(!range.contains(key)) {
			// The key is not contained in this range (falls in a gap between this range
			// and the next higher range), or is higher than the highest range in the map.
			return null;
		}
		else {
			// We found the {@link Range} containing the key; retrieve the corresponding {@link Collection}.
			return objectsMap.get(range);
		}
	}
	
	/**
	 * Retrieve the {@link Collection<T>} corresponding to the {@link Range} of the given
	 * index.
	 * @param index
	 * 	The index of the {@link Range} that we're interested in.
	 * @return
	 * 	The {@link Collection} corresponding to the {@link Range} of the given index, or null if
	 * there is no such {@link Range}.
	 */
	public Collection<T> get(int index) {
		if(index < 0 || index >= size()) {
			return null;
		}
		else {
			return objectsMap.get(ranges[index]);
		}
	}
	
	/**
	 * Retrieve the {@link Range} of the given index.
	 * @param index
	 * 	The index of the {@link Range} that we're interested in.
	 * @return
	 * 	The {@link Range} of the given index, or null if there is no such {@link Range}.
	 */
	public Range getRange(int index) {
		if(index < 0 || index >= size()) {
			return null;
		}
		else {
			return ranges[index];
		}
	}
	
	/**
	 * Retrieve the full array of {@link Range}s.
	 * @return
	 * 	The full array of {@link Range}s.
	 */
	public Range[] getRanges() {
		return ranges;
	}
	
	/**
	 * Add the object to the {@link Collection} corresponding to the {@link Range} containing
	 * the given key.
	 * @param key
	 * 	The key.
	 * @param object
	 * 	The object to add to the collection; if the key is not contained by any {@link Range} then
	 * the object is not added to the map.
	 * @return
	 * 	Boolean stating whether the object was added to the map (true) or not (false).
	 */
	public boolean add(double key, T object) {
		
		Collection<T> objects = get(key);
		
		if(objects==null) {
			return false;
		}
		else {
			objects.add(object);
			return true;
		}
	}
	
	/**
	 * Get the number of {@link Range} to {@link Collection} entries contained in the map
	 * @return
	 * 	The number of {@link Range} to {@link Collection} entries contained in the map
	 */
	public int size() {
		return objectsMap.size();
	}
	
	/**
	 * Get the number of objects stored in the {@link RangeMap}.
	 * @return
	 * 	The number of objects stored in the {@link RangeMap}.
	 */
	public int numberOfObjects() {
		int n=0;
		for(Collection<T> objects : objectsMap.values()) {
			n += objects.size();
		}
		return n;
	}
    
}