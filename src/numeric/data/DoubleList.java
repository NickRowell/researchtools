package numeric.data;

import java.util.ArrayList;
import java.util.Collections;

/**
 * Class provides a concrete List type useful for accumulating double values and computing
 * various statistical quantities.
 * 
 * @author nrowell
 * @version $Id$
 */
public class DoubleList extends ArrayList<Double> {
	
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
	public double getMedian() {
		if(size()==0) 
			throw new IllegalStateException("getMedian() called on empty DoubleList!");

		// Sort list elements into ascending order.
		sortIfRequired();
		
		// Get median value
        int N = size();
        if(N==1) return get(0);
        double median = (N%2==0) ? (get(N/2-1) + get(N/2))/2.0 :
        							get((N-1)/2);
		
		return median;
	}

	/**
	 * Computes the mean value of the elements in the list.
	 * @return	The mean value in the list.
	 */
	public double getMean()
	{
		if(size()==0) 
			throw new IllegalStateException("getMean() called on empty DoubleList!");
		double sum = 0.0f;
		for(Double value : this)
		{
			sum += value;
		}
		return sum/this.size();
	}
	
	/**
	 * Gets the minimum and maximum values stored in the List.
	 * @return	2-element double array containing the minimum and maximum values in the list.
	 */
	public double[] getMinMax()
	{
		// Sort list elements into ascending order.
		sortIfRequired();
		return new double[]{get(0), get(size()-1)};
	}

	/**
	 * Computes the median of absolute deviations.
	 * @return	A 2-element double array containing the median value in element 0 and
	 * 			the median of absolute deviations in element 1.
	 */
	public double[] getMAD()
	{
		// Get the median value
		double median = getMedian();
		
		// Create a new list to store the absolute deviations
		DoubleList absDev = new DoubleList();
		
		// Compute & store absolute deviation in the new list
		for(Double value : this)
		{
			absDev.add(Math.abs(value - median));
		}
		
		// Get the median of the absolute deviations
		return new double[]{median, absDev.getMedian()};
	}

	/**
	 * Copies the list elements to an array of floats, using Apache Commons ArrayUtils
	 * class.
	 * @return	Array of floats representing the List elements in proper sequence.
	 */
//	public double[] toPrimitiveArray()
//	{
//		return ArrayUtils.toPrimitive(toArray(new Double[0]));
//	}

	/**
	 * Sorts the list elements into ascending order, but only if the list has been modified
	 * since the last sorting operation.
	 * @return	Whether a sort occurred.
	 */
	private boolean sortIfRequired()
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
	public Double set(int index, Double element)
	{
		this.oldModCount = this.modCount-1;
		return super.set(index, element);
	}
	
	/**
	 * Thread safe add method, for use when handling a DoubleList within a multi-threaded
	 * environment.
	 * @param addMe	The Double to add to the list.
	 * @return		true (as specified by Collection.add(E)).
	 */
	public synchronized boolean synchronizedAdd(Double addMe)
	{
		return super.add(addMe);
	}

	/**
	 * Thread safe addAll method, for use when handling a DoubleList within a multi-threaded
	 * environment.
	 * @param addMe	The ArrayList<Double> to add to the list.
	 * @return		true (as specified by Collection.add(E)).
	 */
	public synchronized boolean synchronizedAddAll(ArrayList<Double> addMe)
	{
		return super.addAll(addMe);
	}
}
