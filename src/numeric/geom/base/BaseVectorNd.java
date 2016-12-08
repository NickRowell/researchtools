package numeric.geom.base;

import Jama.Matrix;

/**
 * Base class for vector objects in variable number of dimensions. Extending classes should implement
 * any operations specific to their number of dimensions. Note that all extending classes must provide
 * the constructor:
 * 
 * {@code public Vector3d(double... components)}
 * 
 * and implement the factory method:
 * 
 * {@code public T newInstance(double... components)}
 * 
 * @author nrowell
 * @version $Id$
 * @param <T>
 * 	The type of the class that extends this one.
 */
public abstract class BaseVectorNd<T extends BaseVectorNd<T>> {

    /**
     * Angular threshold that defines parallel vectors.
     */
    static double PARALLEL_THRESHOLD = Math.toRadians(0.0001);
    
    /**
     * Magnitude threshold that defines unit vectors.
     */
    static double UNIT_THRESHOLD = 1E-9;
    
    /**
     * Magnitude threshold that defines equality of floating point numbers.
     */
    static double ZERO_THRESHOLD = 1E-9; 
    
	/**
	 * The vector components.
	 */
	protected double[] components;
	
	/**
	 * The number of components.
	 */
	protected int n;
	
	/**
	 * Default constructor.
	 */
	public BaseVectorNd() {
	}
	
	/**
	 * Main constructor.
	 * 
	 * @param components
	 * 	The vector components.
	 */
	public BaseVectorNd(double... components) {
		n = components.length;
		this.components = new double[n];
		System.arraycopy(components, 0, this.components, 0, components.length);
	}
	
	/**
	 * Copy constructor.
	 * @param copyme
	 * 	The {@link BaseVectorNd} to copy.
	 */
	public BaseVectorNd(BaseVectorNd<T> copyme) {
		this(copyme.components);
	}

	/**
	 * Constructs a {@link BaseVectorNd} from a row- or column-matrix.
	 * 
	 * @param M
	 * 	The {@link Matrix}.
	 */
	public BaseVectorNd(Matrix M) {
        if(M.getColumnDimension()==1) {
        	components = new double[M.getRowDimension()];
        	for(int c=0; c<components.length; c++) {
        		components[c] = M.get(c, 0);
        	}
        }
        else if(M.getRowDimension()==1) {
        	components = new double[M.getColumnDimension()];
        	for(int c=0; c<components.length; c++) {
        		components[c] = M.get(0, c);
        	}
        }    
        else {
            throw new RuntimeException("Matrix must be have row or column dimension of one!");
        }
    }

	/**
	 * Factory method: create a new instance of the generic type.
	 * 
	 * @param components
	 * 	The vector components.
	 * @return
	 * 	A new instance of the generic type.
	 */
	public abstract T newInstance(double... components);
	
	/**
	 * Get the components of this {@link BaseVectorNd}.
	 * @return
	 * 	The components of this {@link BaseVectorNd}.
	 */
	public double[] getComponents() {
		return components;
	}
	
    /**
     * Checks if two {@link BaseVectorNd} are of equal length.
     * @param a
     * 	The first {@link BaseVectorNd}.
     * @param b
     * 	The second {@link BaseVectorNd}.
     */
	private static <T extends BaseVectorNd<T>> void checkEqualLength(BaseVectorNd<T> a, BaseVectorNd<T> b) {
		if(a.n != b.n) {
    		throw new IllegalArgumentException("Vectors are of different lengths ("+a.n+" & "+b.n+")!");
    	}
	}
	
    /**
     * Converts the {@link BaseVectorNd} to a JAMA row Matrix containing the vector components.
     * 
     * @return
     * 	A (1xN) Jama Matrix.
     */
    public Matrix toRowMatrix() {
    	Matrix r = new Matrix(1, components.length);
    	for(int c=0; c<components.length; c++) {
    		r.set(0, c, components[c]);
    	}
    	return r;
    }
    
    /**
     * Converts the {@link BaseVectorNd} to a JAMA column Matrix containing the vector components.
     * 
     * @return
     * 	A (Nx1) Jama Matrix.
     */
    public Matrix toColumnMatrix() {
    	Matrix r = new Matrix(components.length, 1);
    	for(int c=0; c<components.length; c++) {
    		r.set(c, 0, components[c]);
    	}
    	return r;
    }

    /**
     * Vector subtraction.
     * @param b
     *  {@link BaseVectorNd} to subtract from this.
     * @return
     *  C = this - b
     */
    public T minus(T b) {
    	checkEqualLength(this, b);
    	double[] components = new double[this.n];
    	for(int c=0; c<components.length; c++) {
    		components[c] = this.components[c] - b.components[c];
    	}
    	return newInstance(components);
    }

    /**
     * Vector subtraction in place
     * @param b
     *  Vector to subtract
     */
    public void minusEquals(BaseVectorNd<T> b) {
    	checkEqualLength(this, b);
    	for(int c=0; c<components.length; c++) {
    		this.components[c] -= b.components[c];
    	}
    }
    
    /**
     * Vector addition.
     * @param b
     *  {@link BaseVectorNd} to add from this.
     * @return
     *  C = this - b
     */
    public T add(T b) {
    	checkEqualLength(this, b);
    	double[] components = new double[this.n];
    	for(int c=0; c<components.length; c++) {
    		components[c] = this.components[c] + b.components[c];
    	}
    	return newInstance(components);
    }

    /**
     * Vector addition in place
     * @param b
     *  Vector to add
     */
    public void addEquals(BaseVectorNd<T> b) {
    	checkEqualLength(this, b);
    	for(int c=0; c<components.length; c++) {
    		this.components[c] += b.components[c];
    	}
    }
    
    /**
     * Scalar multiplication
     * @param b
     * 	Scalar to multiply vector by.
     * @return
     * 	(bx,by,bz)
     */
    public T mult(double b) {
    	double[] components = new double[this.n];
    	for(int c=0; c<components.length; c++) {
    		components[c] = this.components[c] * b;
    	}
    	return newInstance(components);
    }

    /**
     * Scalar multiplication in place
     * @param b  Scalar to multiply vector by.
     */
    public void multEquals(double b) {
    	for(int c=0; c<components.length; c++) {
    		components[c] *= b;
    	}
    }

    /**
     * Get the outer, or tensor, product between one vector and another.
     * @param b
     * @return
     */
    public double[][] outerProduct(BaseVectorNd<T> b) {
    	double[][] product = new double[this.n][b.n];
    	for(int c=0; c<this.n; c++) {
    		for(int d=0; d<b.n; d++) {
    			product[c][d] = this.components[c] * b.components[d];
    		}
    	}
        return product;
    }
    
    /**
     * Dot product of two vectors.
     * @param b
     * 	Second vector
     * @return
     *  this dot b
     */
    public double dot(BaseVectorNd<T> b) {
    	checkEqualLength(this, b);
    	double acc = 0.0;
    	for(int c=0; c<components.length; c++) {
    		acc += this.components[c] * b.components[c];
    	}
    	return acc;
    }

    /**
     * Get a unit vector
     * @return unit vector in direction of this vector
     */
    public T normalise() { 
    	return mult(1.0/norm());
    }
    
    /**
     * In-place normalisation.
     */
    public void normaliseInPlace() { 
        multEquals(1.0/norm());
    }   
    
    /**
     * Magnitude squared
     * @return |this|^2
     */
    public double norm2() {
    	double norm2 = 0.0;
    	for(int c=0; c<components.length; c++) {
    		norm2 += this.components[c] * this.components[c];
    	}
    	return norm2;
    }

    /**
     * Magnitude of a vector
     * @return  |this|
     */
    public double norm() {
    	return Math.sqrt(norm2());
    }

    /**
     * Check if this is a unit vector.
     */
    public boolean isUnitVector() {
        return Math.abs(1.0 - norm()) < UNIT_THRESHOLD;
    }
    

    /**
     * Check if this vector is parallel to another, to within a fixed angular
     * threshold.
     * 
     * @param b
     * @return 
     */
    public boolean isParallelTo(BaseVectorNd<T> b) throws IllegalArgumentException {
        
        // Check that each vector is of finite length
        if(this.norm2()==0) {
        	throw new IllegalArgumentException("Vector is of zero length! "+this.toString());
        }
        if(b.norm2()==0) {
        	throw new IllegalArgumentException("Vector is of zero length! "+b.toString());
        }
        
        // Get unit vectors in direction of each vector
        BaseVectorNd<T> a = this.normalise();
        b.normaliseInPlace();
        
        // Get cosine of internal angle
        double cosAng = a.dot(b);
        
        // Check if it is greater than set threshold
        if(cosAng > Math.cos(PARALLEL_THRESHOLD)) {
            // vectors are parallel
            return true;
        }
        
        // vectors aren't parallel
        return false;
        
    }

    /**
     * TODO: could move this to superclass.
     * Overrides toString() method in parent class
     * @return String representation of Vector3d of the form (x,y,z)
     */
    @Override
    public String toString() {
    	
    	StringBuilder format = new StringBuilder();
    	format.append("(");
    	for(int c=0; c<(n-1); c++) {
    		format.append(String.format("%.5f, ", components[c]));
    	}
    	format.append(String.format("%.5f)", components[n-1]));
    	
        return format.toString();
    }
    
    /**
     * {@inheritDoc}
     */
    @SuppressWarnings("unchecked")
	@Override
    public boolean equals(Object obj)
    {
        if(obj == null) {
        	return false;
        }
        if(obj == this) {
        	return true;
        }
        if (!(obj instanceof BaseVectorNd<?>)) {
            return false;
        }
        
        return this.vectorEquals((BaseVectorNd<T>)obj);
    }
    
    /**
     * Assess the equality of two {@link BaseVectorNd}s based on the equality of their
     * floating point components, to within a threshold set by {@link BaseVectorNd#ZERO_THRESHOLD}.
     * 
     * @param that
     * 	The {@link BaseVectorNd} to be compared to this.
     * @return
     * 	True if the {@link BaseVectorNd}s are component-wise equal.
     */
    public boolean vectorEquals(BaseVectorNd<T> that) {
    	
    	if(this.n != that.n) {
    		// Different numbers of components
    		return false;
    	}
    	
    	for(int c=0; c<components.length; c++) {
    		if(Math.abs(this.components[c] - that.components[c]) > ZERO_THRESHOLD) {
    			// Difference in components
    			return false;
    		}
    	}
        return true;
    }
    
}