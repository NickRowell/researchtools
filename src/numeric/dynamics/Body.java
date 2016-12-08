/**
 * Copyright (c) 2011 University of Dundee
 *
 * Name:
 * Body.java
 *
 * Purpose:
 * This class represents a rigid body moving with six degrees of freedom (3
 * translation & 3 rotation). It extends the basic Pose operations with
 * methods to integrate the motion of the body under the action
 * of inertial forces.
 * 
 * Language:
 * Java
 *
 * Author:
 * Nicholas Rowell (nickrowell@computing.dundee.ac.uk)
 *
 */

package numeric.dynamics;

import numeric.geom.dim3.Quaternion;
import numeric.geom.dim3.Vector3d;

public class Body extends Pose{

    /** Gravitational constant. */
    public static final double G = 6.673E-11;
    
    /** Linear velocity in world frame. */
    Vector3d velocity;

    /** Angular velocity in world frame. */
    Vector3d angVelocity;

    /** Mass of body. */
    double mass;

    /**
     * Main constructor.
     * @param pos      Initial position vector
     * @param vel      Initial velocity vector
     * @param att      Initial attitude quaternion
     * @param angVel   Angular velocity vector that is constant throughout flight. Units are radians per second.
     * @param m        Mass of body
     */
    public Body(Vector3d pos, Vector3d vel, Quaternion att, Vector3d angVel, double m){

        // Call Pose constructor to set position & attitude
        super(pos,att);
        
        setLinearVelocity(vel);
        setAngularVelocity(angVel);
        setMass(m);
      
    }

    /** Setters. */
    public final void setMass(double m){this.mass = m;}
    public final void setLinearVelocity(Vector3d vel){this.velocity = new Vector3d(vel);}
    public final void setAngularVelocity(Vector3d angVel){this.angVelocity = new Vector3d(angVel);}

    /** Getters. */
    public double     getMass(){return this.mass;}
    public Vector3d   getAngularVelocity(){ return this.angVelocity;}
    public Vector3d   getLinearVelocity(){ return this.velocity;}
    
    /** Integrate linear acceleration. */
    public void leapVelocity(Vector3d force, double dT){
	this.velocity = this.velocity.add(force.mult(dT/this.mass));
    }

    /** Integrate linear velocity. */
    public void leapPosition(double dT){
	this.position = this.position.add(this.velocity.mult(dT));
    }

    /** Integrate angular velocity. */
    public void leapAttitude(double dT){
        
        Quaternion rotation = Quaternion.getIdentity();
    
        if(angVelocity.norm() > 9E-9){

            // Get angle through which spacecraft rotates 
            // each integration step, in radians
            double angle = angVelocity.norm() * dT;

            // Get unit vector in direction of rotation axis
            Vector3d u = angVelocity.normalise();

            rotation = new Quaternion(u, angle);
            
        }
        
        // Update attitude by quaternion multiplication
        attitude = rotation.multiply(attitude);
        
    }


    /** Apply friction to linear motion. */
    public void friction(double mag, double dT){

        // Get vector in direction of velocity
        Vector3d fric = new Vector3d(velocity);

        // Multiply by negative 1 to make it a braking force
        fric.multEquals(-1);

        // Multiply by magnitude
        fric.multEquals(mag);

        this.leapVelocity(fric, dT);
    
    }
    
    /**
     * Get gravitational force acting on this body due to mass of other body.
     * @param that  Other body.
     * @return 
     */
    public Vector3d getForce(Body that){

        // Vector joining centres
        Vector3d relpos = this.position.minus(that.position);
        
        // Magnitude of the force. Watch out for divide-by-zero
        double constant = (relpos.norm2()>0.0) ? -1.0 * (G * mass * that.mass)/relpos.norm2() : 0.0;

        return relpos.normalise().mult(constant);
    }   
}
