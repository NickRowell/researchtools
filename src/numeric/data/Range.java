package numeric.data;

/**
 * Class represents a range of values; inclusive of the lower boundary.
 *
 * @author nrowell
 * @version $Id$
 */
public class Range {
	
	/**
	 * The (inclusive) lower limit.
	 */
	public final double lower;
	
	/**
	 * The (exclusive) upper limit.
	 */
	public final double upper;
	
	/**
	 * The main contructor for a {@link Range}. The limits cannot be NaN or infinities,
	 * cannot be equal and must be in the right order (upper > lower).
	 * @param lower
	 * 	The (inclusive) lower limit.
	 * @param upper
	 * 	The (exclusive) upper limit.
	 */
	public Range(double lower, double upper) {
		
		// Check health of the limits
		if(Double.isNaN(lower) || Double.isNaN(upper) || Double.isInfinite(lower) || Double.isInfinite(upper)) {
			throw new IllegalArgumentException("Limits cannot be NaN or Infinity! lower = "+lower+"; upper = "+upper);
		}
		
		// Check validity of the limits
		if(lower >= upper) {
			throw new IllegalArgumentException("Lower limit must be less than upper limit! lower = "+lower+"; upper = "+upper);
		}
		
		this.lower = lower;
		this.upper = upper;
	}
	
	/**
	 * Determines if the given value falls within the boundaries of the {@link Range}.
	 * @param value
	 * 	The value to test
	 * @return
	 * 	True if the value lies within the boundaries of the {@link Range}; false otherwise.
	 */
	public boolean contains(double value) {
		return (value >= lower && value < upper);
	}
	
	/**
	 * Get the midpoint of the range.
	 * @return
	 * 	The midpoint of the range.
	 */
	public double mid() {
		return (lower + upper)/2.0;
	}
	
	/**
	 * Get the width of the range.
	 * @return
	 * 	The width of the range.
	 */
	public double width() {
		return upper - lower;
	}
	
	/**
	 * {@inheritDoc}
	 */
	@Override
	public String toString() {
		return "["+lower+":"+upper+"]";
	}
}