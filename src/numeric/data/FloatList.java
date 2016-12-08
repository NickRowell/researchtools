package numeric.data;

import java.util.ArrayList;
import java.util.Collections;

/**
 * Class provides a concrete List type useful for accumulating float values and computing
 * various statistical quantities.
 *
 * @author nrowell
 * @version $Id$
 */
public class FloatList extends ArrayList<Float>
{
	private static final long serialVersionUID = 1L;

	/** 
	 * The modification count of the underlying List at the time of the last sorting operation.
	 * Used to minimise the amount of unnecessary sorting.
	 */
    private int oldModCount;

	/**
	 * Computes the median value of the elements in the List. Note that on exit, the list
	 * elements will have been sorted into ascending order.
	 * @return	The median value in the List.
	 */
	public float getMedian()
	{
		if(size()==0) 
			return Float.NaN;
		
		// Sort list elements into ascending order.
		sortIfRequired();
		
		// Get median value
        int N = size();
        if(N==1) return get(0);
        
        // Even number of elements: return average of middle two
        // Odd number of elements: return central object
        float median = (N%2==0) ? (get(N/2-1) + get(N/2))/2.0f :
        							get((N-1)/2);
		
		return median;
	}
	
	/**
	 * Computes percentiles for the FLoatList, i.e. for a given argument in the range 0->100%, we compute the
	 * value below which the given percentage of elements fall.
	 * @param percentile	The desired percentile level [0.0->100.0].
	 * @return				The value below which the fraction of data specified fall.
	 */
	public float getPercentile(float percentile)
	{
		// Sanity check
		if(percentile < 0.0 || percentile > 100.0)
			throw new IllegalArgumentException("Percentile must lie in range 0:100! Found "+percentile);
		
		if(size()==0) 
			return Float.NaN;
		
		// Sort list elements into ascending order.
		sortIfRequired();
		
		int N = size();
		
		// Check for extreme percentiles, or lists with only 1 element, both of which will fail the algorithm below
		if(percentile==0.0)
			return get(0);
		if(percentile==100.0)
			return get(N-1);
		if(N==1)
			return get(0);
		
		// Main algorithm for general percentile computation
		
		// Get number of elements in given percentile
		float percN = percentile * ((float)N/100.0f);
		
		// Get the indices of the two elements closest to this point in the list
		int lower = (int)Math.floor(percN);
		int upper = (int)Math.ceil(percN);
		
		// Take the average of these two points as the percentile
		return (get(lower) + get(upper))/2.0f;
		
	}
	
	/**
	 * Computes the mean value of the elements in the list.
	 * @return	The mean value in the list.
	 */
	public float getMean()
	{
		if(size()==0) 
			return Float.NaN;
		
		return getSum()/this.size();
	}
	
	/**
	 * Computes the sum of the array elements
	 * @return
	 * 	The sum of the array elements
	 */
	public float getSum() {
		float sum = 0.0f;
		for(Float value : this)
		{
			sum += value;
		}
		return sum;
	}
	
	
	/**
	 * Gets the minimum and maximum values stored in the List.
	 * @return	2-element double array containing the minimum and maximum values in the list.
	 */
	public float[] getMinMax()
	{
		if(size()==0) 
			return new float[]{Float.NaN,Float.NaN};
		
		// Sort list elements into ascending order.
		sortIfRequired();
		return new float[]{get(0), get(size()-1)};
	}

	/**
	 * Computes the median of absolute deviations.
	 * @return	A 2-element float array containing the median value in element 0 and
	 * 			the median of absolute deviations in element 1.
	 */
	public float[] getMAD()
	{
		if(size()==0) 
			return new float[]{Float.NaN,Float.NaN};
		
		// Get the median value
		float median = getMedian();
		
		// Create a new list to store the absolute deviations
		FloatList absDev = new FloatList();
		
		// Compute & store absolute deviation in the new list
		for(Float value : this)
		{
			absDev.add(Math.abs(value - median));
		}
		
		// Get the median of the absolute deviations
		return new float[]{median, absDev.getMedian()};
	}
	
	/**
	 * Copies the list elements to an array of floats, using Apache Commons ArrayUtils
	 * class.
	 * @return	Array of floats representing the List elements in proper sequence.
	 */
//	public float[] toPrimitiveArray()
//	{
//		return ArrayUtils.toPrimitive(toArray(new Float[0]));
//	}
	
	/**
	 * Sorts the list elements into ascending order, but only if the list has been modified
	 * since the last sorting operation.
	 * @return	Whether a sort occurred.
	 */
	protected boolean sortIfRequired()
	{
		if(isSortRequired())
		{
			// List has been structurally modified since last sort: resort
			Collections.sort(this);
			this.oldModCount = this.modCount;
			return true;
		}
		return false;
	}
	
	/**
	 * Determines if the list is already sorted into numerical order by checking the number
	 * of structural modifications that occurred since the last sort.
	 * @return	True if the list requires sorting.
	 */
	public boolean isSortRequired()
	{
		return this.oldModCount != this.modCount;
	}
	
	/**
	 * {@inheritDoc}
	 * 
	 * Overrides the set method. This is necessary because set operations aren't considered a
	 * structural modification although they can alter the order of elements. So, if set is
	 * called then we want to indicate that a sort is necessary the next time sortIfRequired()
	 * is called. This is done by setting the oldModCount to one less than the current modCount,
	 * which spoofs the sortIfRequired() method into thinking that a structural modification has
	 * occurred since the last sort.
	 */
	public Float set(int index, Float element)
	{
		this.oldModCount = this.modCount-1;
		return super.set(index, element);
	}

	/**
	 * Thread safe add method, for use when handling a FloatList within a multi-threaded
	 * environment.
	 * @param addMe	The Float to add to the list.
	 * @return		true (as specified by Collection.add(E)).
	 */
	public synchronized boolean synchronizedAdd(Float addMe)
	{
		return super.add(addMe);
	}
	
	/**
	 * Thread safe addAll method, for use when handling a FloatList within a multi-threaded
	 * environment.
	 * @param addMe	The ArrayList<Float> to add to the list.
	 * @return		true (as specified by Collection.add(E)).
	 */
	public synchronized boolean synchronizedAddAll(ArrayList<Float> addMe)
	{
		return super.addAll(addMe);
	}
	
}
