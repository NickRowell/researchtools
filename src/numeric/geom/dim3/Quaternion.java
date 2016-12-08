package numeric.geom.dim3;

import Jama.Matrix;

/**
 * This class provides objects representing Quaternions and operations on them.
 *
 * @author nrowell
 * @version $Id$
 */
public class Quaternion {

    
    /** Magnitude threshold that defines unit Quaternions. */
    static double UNIT_THRESHOLD = 1E-9;
    
    /** Magnitude threshold that defines equality of floating point numbers. */
    static double ZERO_THRESHOLD = 1E-9;
    
    
    
    /**
     * Quaternion elements; real part and imaginary vector
     */
    public double   re;   // q0
    public Vector3d im;   // q1 q2 q3

    /**
     * Main constructor taking quaternion elements as argument.
     * @param Q0   Real component
     * @param Q1   Imaginary i component
     * @param Q2   Imaginary j component
     * @param Q3   Imaginary k component
     */
    public Quaternion(double Re, double Q1, double Q2, double Q3)
    {
    	setRe(Re);
        setIm(Q1,Q2,Q3);
    }

    /**
     * Constructor using axis-angle definition of a rotation to create quaternion
     * @param axis   {@link #Vector3d Vector3d} object defining axis of rotation. Must be normalised.
     * @param angle  Angle through which to rotate, in radians.
     */
    public Quaternion(Vector3d axis, double angle)
    {

        if(!axis.isUnitVector()){
            System.err.println("Quaternion(Vector3d axis, double angle): "+
                    "normalising rotation axis vector.");
            axis.normaliseInPlace();
        }

        double cosAng2 = Math.cos(angle/2.0);
        double sinAng2 = Math.sin(angle/2.0);

        setRe(cosAng2);
        setIm(axis.mult(sinAng2));
    }

    /** Constructor for Quaternion from rotation matrix. */
    public Quaternion(Matrix M)
    {

        // Check dimensions
        if(M.getRowDimension()!=3 || M.getColumnDimension()!=3)
            throw new RuntimeException("Illegal conversion of matrix to "+
                                       "quaternion.");
        
        // Check determinant
        if(Math.abs(M.det() - 1)>0.0001)
            throw new RuntimeException("Illegal conversion of matrix to "+
                                       "quaternion. Matrix does not have a " +
                                       "determinant of 1: "+M.det());
        
        double qw=0,qx=0,qy=0,qz=0;

        double tr = M.trace();

        if (tr > 0) {
            double S = Math.sqrt(tr + 1.0) * 2; // S=4*qw
            qw = 0.25 * S;
            qx = (M.get(2,1) - M.get(1,2)) / S;
            qy = (M.get(0,2) - M.get(2,0)) / S;
            qz = (M.get(1,0) - M.get(0,1)) / S;
        }
        else if ((M.get(0,0) > M.get(1, 1)) & (M.get(0,0) > M.get(2,2))) {
            double S = Math.sqrt(1.0 + M.get(0,0) - M.get(1,1) - M.get(2,2)) * 2; // S=4*qx
            qw = (M.get(2,1) - M.get(1,2)) / S;
            qx = 0.25 * S;
            qy = (M.get(0,1) + M.get(1,0)) / S;
            qz = (M.get(0,2) + M.get(2,0)) / S;
        }
        else if (M.get(1,1) > M.get(2,2)) {
            double S = Math.sqrt(1.0 + M.get(1,1) - M.get(0,0) - M.get(2,2)) * 2; // S=4*qy
            qw = (M.get(0,2) - M.get(2,0)) / S;
            qx = (M.get(0,1) + M.get(1,0)) / S;
            qy = 0.25 * S;
            qz = (M.get(1,2) + M.get(2,1)) / S;
        }
        else {
            double S = Math.sqrt(1.0 + M.get(2,2) - M.get(0,0) - M.get(1,1)) * 2; // S=4*qz
            qw = (M.get(1,0) - M.get(0,1)) / S;
            qx = (M.get(0,2) + M.get(2,0)) / S;
            qy = (M.get(1,2) + M.get(2,1)) / S;
            qz = 0.25 * S;
        }

        // Now set internal quaternion parameters
        re = qw;
        im = new Vector3d(qx,qy,qz);

    }

    /**
     * Constructor setting real and imaginary parts explicitly.
     * @param real
     * @param imag
     */
    public Quaternion(double real, Vector3d imag)
    {

        setRe(real);
        setIm(imag);

    }

    /**
     * Copy constructor.
     * @param copyme  Quaternion object to copy
     */
    public Quaternion(Quaternion copyme)
    {

	setRe(copyme.re);
        setIm(copyme.im);

    }

    /**
     * Default constructor with all elements set to zero.
     */
    public Quaternion(){
	setRe(0.0);
	setIm(0.0, 0.0, 0.0);
        
    }

    public boolean equals(Object obj)
    {
        if(obj == null) return false;
        if(obj == this) return true;
        if (!(obj instanceof Quaternion))
            return false;
        Quaternion that = (Quaternion)obj;
        
        return this.equals(that);
        
    }
    
    public boolean equals(Quaternion that)
    {
        return (Math.abs(re - that.re) < ZERO_THRESHOLD && im.equals(that.im));
    }
    
    /**
     * Set real (q0) component.
     */
    public final void setRe(double Re){re = Re;}
    /**
     * Set imaginary vector by specifying each component.
     * @param x x component of imaginary vector
     * @param y y component of imaginary vector
     * @param z z component of imaginary vector
     */
    public final void setIm(double x, double y, double z){im = new Vector3d(x,y,z);}
    /**
     * Set imaginary vector by copying an existing vector
     * @param u     Vector3d object
     */
    public final void setIm(Vector3d u){im = new Vector3d(u.getComponents());}

    /**
     * Set i (q1) component
     */
    public final void setQ1(double Q1){im.setX(Q1);}
    /**
     * Set j (q2) component
     */
    public final void setQ2(double Q2){im.setY(Q2);}
    /**
     * Set k (q3) component
     */
    public final void setQ3(double Q3){im.setZ(Q3);}

    /**
     * Get real (q0) component
     */
    public double getRe(){return this.re;}
    /**
     * Get imaginary vector
     */
    public Vector3d getIm(){return this.im;}

    /**
     * Get i (q1) component
     */
    public double getQ1(){return this.im.getX();}
    /**
     * Get j (q2) component
     */
    public double getQ2(){return this.im.getY();}
    /**
     * Get k (q3) component
     */
    public double getQ3(){return this.im.getZ();}

    /**
     * Overrides toString() method of parent class.
     * @return String representation in the form (q0,q1,q2,q3) with components rounded off to third decimal place
     */
    @Override
    public String toString()
    {
    	return String.format("(%.5f, %.5f, %.5f, %.5f)", re, im.getX(), im.getY(), im.getZ());
    }

    /**
     * Returns the identity quaternion corresponding to zero rotation
     * @return Identity quaternion (1,0,0,0)
     */
    public static Quaternion getIdentity(){ return new Quaternion(1,0,0,0);}

    /**
     * Static quaternion multiplication. Invokes methods of the Vector3d class
     * to do certain operations.
     *
     * @param A    First quaternion
     * @param B    Second quaternion
     * @return     C = A*B
     */
    public static Quaternion multiply(Quaternion A, Quaternion B)
    {

        Quaternion mult = new Quaternion(A.re*B.re - A.im.dot(B.im),
                                         new Vector3d(B.im.mult(A.re).add(
                                                      A.im.mult(B.re).add(
                                                      A.im.cross(B.im))).getComponents()));

        return mult;

    }

    /**
     * Non-static quaternion multiplication
     * @param A   Second quaternion
     * @return    C = this * A
     */
    public Quaternion multiply(Quaternion A)
    {

        Quaternion mult = new Quaternion(this.re*A.re - this.im.dot(A.im),
                                          A.im.mult(this.re).add(
                                          this.im.mult(A.re).add(
                                          this.im.cross(A.im))));

        return mult;

    }

    /**
     * Quaternion addition.
     */
    public static Quaternion add(Quaternion A, Quaternion B)
    {
        return new Quaternion(A.re   + B.re,
                              A.im.getX() + B.im.getX(),
                              A.im.getY() + B.im.getY(),
                              A.im.getZ() + B.im.getZ());
    }
    
    /**
     * Quaternion addition.
     */
    public Quaternion add(Quaternion A)
    {
         return new Quaternion(A.re   + re,
                               A.im.getX() + im.getX(),
                               A.im.getY() + im.getY(),
                               A.im.getZ() + im.getZ());   
    }

    /**
     * Get quaternion inverse, i.e. quaternion corresponding to reverse rotation
     * @return   (q0,-q1,-q2,-q3)
     */
    public Quaternion inverse()
    {
        return new Quaternion(this.re, this.im.mult(-1.0));
    }

    /**
     * Reverses sign of all Quaternion elements. If this Quaternion represents
     * a rotation, then the resulting Quaternion is mathematically identical to 
     * original.
     */
    public final void changeSign(){ this.re *= -1.0; this.im.multEquals(-1.0);}

    /**
     * Get quaternion magnitude. Normal usage should not change the 
     * normalisation of quaternions, but accumulation of floating point errors
     * over many operations may result in non-unity magnitude.
     * This can be checked intermittently if desired with this method.
     * @return sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3)
     */
    public double norm(){ return Math.sqrt(re*re + im.norm2());}


    /** Normalise the Quaternion in place */
    public void normalise()
    {

        // Get current normalisation
        double norm = norm();

        re /= norm;
        im.multEquals(1.0/norm);
      
    }

    /**
     * Rotate a vector using a quaternion
     * @param v  Vector to rotate
     * @return   Rotated vector
     */
    public Vector3d rotate(Vector3d v)
    {

        // Represent vector as a quaternion with zero real component
        Quaternion V = new Quaternion(0,v);
        
        // Now do rotation operation. Quaternion multiplication is 
        // associative so difference where brackets are placed
        Quaternion V_prime = this.multiply(V).multiply(this.inverse());

        // Convert rotated quaternion back to vector representation by 
        // discarding real component
        return V_prime.im;

    }

    public double[] getEulerAngles()
    {

        double yaw, pitch, roll;

        // Check for gimbal lock at +/- 90 degrees pitch
        double test = this.getQ2()*this.getQ3() + this.getQ1()*this.getRe();

        if (test > 0.499) { // singularity at north pole
            //System.out.println("gimbal lock north");
	    yaw   = 2 * Math.atan2(this.getQ2(),this.getRe());
            pitch = Math.PI/2;
            roll  = 0;
        }
        else if (test < -0.499) { // singularity at south pole
            //System.out.println("gimbal lock south");
            yaw   = -2 * Math.atan2(this.getQ2(),this.getRe());
            pitch = - Math.PI/2;
            roll  = 0;
        }
        else{
            yaw   = Math.atan2(2*this.getQ3()*this.getRe() - 2*this.getQ2()*this.getQ1() , 1 - 2*this.getQ2()*this.getQ2() - 2*this.getQ3()*this.getQ3());
            pitch = Math.asin(2*test);
            roll  = Math.atan2(2*this.getQ2()*this.getRe() - 2*this.getQ3()*this.getQ1() , 1 - 2*this.getQ1()*this.getQ1() - 2*this.getQ3()*this.getQ3());
        }

        return new double[]{yaw,pitch,roll};

    }


    /**
     * Convert Euler angles to a rotation quaternion.
     * @param yaw
     * @param pitch
     * @param roll
     * @return
     */
    public static Quaternion getQuaternion(double yaw, double pitch, double roll)
    {

            double c1 = Math.cos(yaw/2);
            double s1 = Math.sin(yaw/2);
            double c2 = Math.cos(pitch/2);
            double s2 = Math.sin(pitch/2);
            double c3 = Math.cos(roll/2);
            double s3 = Math.sin(roll/2);
            double c1c2 = c1*c2;
            double s1s2 = s1*s2;

            return new Quaternion(c1c2*c3 - s1s2*s3, c1*s2*c3 - s1*c2*s3, c1c2*s3 + s1s2*c3, s1*c2*c3 + c1*s2*s3);

    }

    /**
     * Convert a normalised quaternion to an orthogonal rotation matrix.
     * 
     * @return
     */
    public Matrix toMatrix(){

          Matrix R = new Matrix(3,3);

          R.set(0, 0, 1 - 2*im.getY()*im.getY() - 2*im.getZ()*im.getZ());
          R.set(0, 1, 2*im.getX()*im.getY() - 2*im.getZ()*re);
          R.set(0, 2, 2*im.getX()*im.getZ() + 2*im.getY()*re);

          R.set(1, 0, 2*im.getX()*im.getY() + 2*im.getZ()*re);
          R.set(1, 1, 1 - 2*im.getX()*im.getX() - 2*im.getZ()*im.getZ());
          R.set(1, 2, 2*im.getY()*im.getZ() - 2*im.getX()*re);

          R.set(2, 0, 2*im.getX()*im.getZ() - 2*im.getY()*re);
          R.set(2, 1, 2*im.getY()*im.getZ() + 2*im.getX()*re);
          R.set(2, 2, 1 - 2*im.getX()*im.getX() - 2*im.getY()*im.getY());

          return R;

      }

    /**
     * 
     * 
     * This method uses first order error propagation to transform the 
     * covariance matrix for the quaternion elements to the covariance matrix
     * for the elements of the rotation matrix that this quaternion can be
     * expressed as.
     * 
     * The method has been tested and found to work using a Monte Carlo 
     * technique in QuaternionTester class. Essentially, a quaternion is 
     * specified and an appropriate covariance matrix is randomly created.
     * This is then used to draw many random realisations of the quaternion,
     * and in each case convert it to a rotation matrix. These are all stored
     * and used to estimate the sample variance/covariance of the matrix
     * elements. This is then compared to the covariance matrix obtained by
     * passing the quaternion covariance matrix to this method. The two agree
     * to within statistical uncertainty associated with finite sample size.
     * 
     * @param S_q Covariance matrix for quaternion elements.
     */
    public Matrix toMatrixCovariance(Matrix S_q)
    {
    
        // Get some handy notation
        double q0 = re;
        double q1 = im.getX();
        double q2 = im.getY();
        double q3 = im.getZ();
        
        // Define Jacobian of rotation matrix elements wrt quaternion elements.
        double[][] drdq = new double[][]{{    0,    0,-4*q2,-4*q3},
                                         {-2*q3, 2*q2, 2*q1,-2*q0},
                                         { 2*q2, 2*q3, 2*q0, 2*q1},
                                         { 2*q3, 2*q2, 2*q1, 2*q0},
                                         {    0,-4*q1,    0,-4*q3},
                                         {-2*q1,-2*q0, 2*q3, 2*q2},
                                         {-2*q2, 2*q3,-2*q0, 2*q1},
                                         { 2*q1, 2*q0, 2*q3, 2*q2},
                                         {    0,-4*q1,-4*q2,    0}};
     
        // Convert to Matrix
        Matrix DRDQ = new Matrix(drdq);
        
        // First order uncertainty propagation
        Matrix S_R = DRDQ.times(S_q).times(DRDQ.transpose());
        
        return S_R;
        
    }
    
    
    /** Check if this is a unit quaternion. */
    public boolean isUnitQuaternion()
    {
        return Math.abs(1.0 - norm()) < UNIT_THRESHOLD;
    }    
    
    
    
    

}
