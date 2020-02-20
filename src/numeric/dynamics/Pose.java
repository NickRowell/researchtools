package numeric.dynamics;

import Jama.Matrix;
import numeric.geom.dim3.Quaternion;
import numeric.geom.dim3.Vector3d;

public class Pose
{
    
    /** Linear position in world frame. */
    protected Vector3d position;
    
    /** Orientation in world frame parameterised by unit quaternion. */
    protected Quaternion attitude;    
    
    /** Default constructor - no transformation. */
    public Pose()
    {
        setPosition(new Vector3d(0,0,0));
        setAttitude(new Quaternion(1,0,0,0));
    }   
    
    /** Main constructor. */
    public Pose(Vector3d pos, Quaternion att)
    {
        setPosition(pos);
        setAttitude(att);
    }
    
    /** Copy constructor. */
    public Pose(Pose copyme)
    {
        setPosition(copyme.position);
        setAttitude(copyme.attitude);
    }
        
    
    /** Constructor from 3x4 extrinsic camera matrix. */
    public Pose(Matrix cam)
    {
        // Orthogonality check made in Quaternion constructor
        setAttitude(new Quaternion(cam.getMatrix(new int[]{0,1,2}, new int[]{0,1,2})));
        setPosition(new Vector3d(cam.get(0,3),cam.get(1,3),cam.get(2,3)));
    }
    
    
    /** Set position in world frame. */
    public final void setPosition(Vector3d pos){this.position = new Vector3d(pos.getComponents());}
    
    public final void setAttitude(Quaternion att){this.attitude = new Quaternion(att);}
    
    /** Get position in world frame. */
    public Vector3d   getPosition(){ return this.position;}
    
    public Quaternion getAttitude(){ return this.attitude;}
    
    /** Get direction of body frame X axis in world frame. */
    public Vector3d getBodyXAxis(){ return attitude.rotate(new Vector3d(1,0,0));}
    /** Get direction of body frame Y axis in world frame. */
    public Vector3d getBodyYAxis(){ return attitude.rotate(new Vector3d(0,1,0));}
    /** Get direction of body frame Z axis in world frame. */
    public Vector3d getBodyZAxis(){ return attitude.rotate(new Vector3d(0,0,1));}
    
    /** Get direction of world frame X axis in body frame. */
    public Vector3d getWorldXAxis(){ return attitude.inverse().rotate(new Vector3d(1,0,0));}
    /** Get direction of world frame Y axis in body frame. */
    public Vector3d getWorldYAxis(){ return attitude.inverse().rotate(new Vector3d(0,1,0));}
    /** Get direction of world frame Z axis in body frame. */
    public Vector3d getWorldZAxis(){ return attitude.inverse().rotate(new Vector3d(0,0,1));}    
    
    /** Convert position vector of point in world frame to body frame. */
    public Vector3d toBodyFrame(Vector3d world){ return attitude.inverse().rotate(world.minus(position));}
    /** Convert position vector of point in body frame to world frame. */
    public Vector3d toWorldFrame(Vector3d body){ return attitude.rotate(body).add(position);}
    
    /** Convert vector form body frame to world frame. */
    public Vector3d getWorldVector(Vector3d bodyVector){ return attitude.rotate(bodyVector);}    
    /** Convert vector form world frame to body frame. */
    public Vector3d getBodyVector(Vector3d worldVector){ return attitude.inverse().rotate(worldVector);}    
    
    
    /** Rotate rank-2 tensor to body frame. */
    public Matrix toBodyFrame(Matrix world)
    {
        // Get orthonormal rotation matrix
        Matrix R    = attitude.toMatrix();
        Matrix invR = attitude.inverse().toMatrix();
        // Perform transformation
        return invR.times(world).times(R);
    }
    
    /** Rotate rank-2 tensor to world frame. */
    public Matrix toWorldFrame(Matrix body)
    {
        // Get orthonormal rotation matrix
        Matrix R    = attitude.toMatrix();
        Matrix invR = attitude.inverse().toMatrix();
        // Perform transformation
        return R.times(body).times(invR);
    }    
    
    /** Get the conjugate Pose. */
    public Pose conjugate()
    {
        // Conjugate rotation is simply inverse
        Quaternion conj_att = attitude.inverse();
        // Conjugate position is position of World origin in body frame
        Vector3d conj_pos = toBodyFrame(new Vector3d(0,0,0));    
        return new Pose(conj_pos, conj_att);
    }
     
    /**
     * This method rotates the body frame so that the given body-frame direction
     * points towards the given world-frame point. This is useful for orientating
     * body so that e.g. camera boresight (which typically points along z axis)
     * points at some other known point such as the centre of an asteroid
     * model or something.
     * @param bodyVector    Unit vector in the body frame.
     * @param worldPoint    Point in the world frame.
     */
    public void pointVectorAt(Vector3d bodyVector, Vector3d worldPoint)
    {
        
        // Get body frame direction vector in world frame. Normalise.
        Vector3d worldVector = attitude.rotate(bodyVector).normalise();

        // Get unit vector from Body origin to point, in world frame
        Vector3d pos = worldPoint.minus(position).normalise();

        // Check if vectors are parallel, to within 0.01 degrees, If this is the 
        // case then worldPoint already lies on bodyVector, so take no action.
        if(pos.isParallelTo(worldVector))
            return;
        
        // Vector perpendicular to these from cross product
        Vector3d rotAxis = worldVector.cross(pos);
        
        // Normalise to get unit vector defining rotation axis.
        rotAxis.normaliseInPlace();
        
        // Angle between vectors
        double angle = Math.acos(pos.dot(worldVector));
        
        // Make a quaternion that represents rotation through this angle about
        // this axis.
        Quaternion rotation = new Quaternion(rotAxis, angle);
        
        // Apply this rotation to the Body
        rotate(rotation);

    }
    
    /**
     * Express the Pose attitude and position as a 3x4 matrix, where the
     * first 3 columns contain the rotation matrix and the fourth contains
     * the translation vector. This form is used widely in computer vision.
     * 
     * In this form, multiplication of P by a homogenous 4-vector expressing
     * position of a point in body/camera frame, i.e. P*X, will result in a
     * 3-vector expressing position of point in world frame.
     * 
     */
    public Matrix getP()
    {
        
        Matrix P = new Matrix(3,4);
        
        // Set rotation matrix elements
        P.setMatrix(new int[]{0,1,2}, new int[]{0,1,2}, attitude.toMatrix());

        // Set translation vector elements
        P.set(0, 3, position.getX());
        P.set(1, 3, position.getY());
        P.set(2, 3, position.getZ());
    
        return P;
    }
    
    /**
     * Conjugate Pose expressed as 3-by-4 matrix.
     * 
     * In this form, multiplication of P by a homogenous 4-vector expressing
     * position of a point in world frame, i.e. P*X, will result in a
     * 3-vector expressing position of point in body/camera frame.
     * 
     */
    public Matrix getConjP()
    {
        
        Pose conjP = this.conjugate();
    
        return conjP.getP();
    }    
    
    
    
    /** 
     * Get relative Pose between two Pose objects. The transformation rules
     * using the relative Pose operate with 'this' Pose corresponding to 
     * the World frame, and 'that' Pose corresponding to a body in the frame.
     */
    public Pose getRelativePose(Pose that)
    {
    
        // Position of 'that' Pose frame origin in 'this' frame
        Vector3d r12 = attitude.inverse().rotate(that.position.minus(this.position));
        
        Quaternion R12 = this.attitude.inverse().multiply(that.attitude);
        
        return new Pose(r12,R12);
    
    }
    
    /** Rotate the Pose by a given world frame rotation quaternion. */
    public void rotate(Quaternion r)
    {
        
        // Check that input quaternion has unit magnitude.
        if(!r.isUnitQuaternion()){
            System.err.println("Rotation quaternion is not normalised!");
            r.normalise();
        }        
        
        attitude = r.multiply(attitude);
        
        // Re-normalise the quaternion to correct for build up of floating
        // point errors.
        if(!attitude.isUnitQuaternion())
            attitude.normalise();
    }
    
    
    
    
    /** Translate the Pose by a given world frame Vector3d */
    public void translate(Vector3d t)
    {
        position.addEquals(t);
    }
    
    /** 
     * Get the Pose object obtained when the current Pose is transformed by
     * the input. The input pose is defined in the frame of this Pose.
     * 
     * This method is useful for example if the current Pose represents one
     * half of a stereo pair of cameras. Call this method and pass as argument
     * the transformation between the cameras, which is defined in the frame
     * of the current camera rather than the World frame because that way the
     * transformation is constant as the cameras move around.
     * 
     */
    public Pose transform(Pose trans)
    {
        
        // Transform position of input Pose to the World frame
        Vector3d position_trans = toWorldFrame(trans.position);
        Quaternion quat_trans = attitude.multiply(trans.attitude);
        
        return new Pose(position_trans, quat_trans);
    }
    
    
}
