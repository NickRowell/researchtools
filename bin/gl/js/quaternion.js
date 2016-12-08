/**
 * Javascript object representing a quaternion.
 */
function Quat(q0, q1, q2, q3) 
{
	this.q0 = q0;    // Scalar part
	this.q1 = q1;    // Vector part (q1,q2,q3)
	this.q2 = q2;
	this.q3 = q3;
}

// Constructs a unit quaternion using the axis/angle form.
// axis  - Vec3 unit vector
// angle - rotation angle in radians
function createAxisAngleQuaternion(axis, angle)
{
	q0 = Math.cos(angle/2.0);
        q1 = axis.x * Math.sin(angle/2.0);
        q2 = axis.y * Math.sin(angle/2.0);
        q3 = axis.z * Math.sin(angle/2.0);

	return new Quat(q0, q1, q2, q3);
}



Quat.prototype = {
    constructor: Quat,
    toMat4x4:function()
    {
	// shorthand
	q0 = this.q0;
	q1 = this.q1;	
	q2 = this.q2;
	q3 = this.q3;

	// 4x4 matrix with no translation part, only rotation
	return new Mat4x4([1 - 2*q2*q2 - 2*q3*q3, 2*q1*q2 - 2*q3*q0, 2*q1*q3 + 2*q2*q0, 0,
                           2*q1*q2 + 2*q3*q0, 1 - 2*q1*q1 - 2*q3*q3, 2*q2*q3 - 2*q1*q0, 0,
                           2*q1*q3 - 2*q2*q0, 2*q2*q3 + 2*q1*q0, 1 - 2*q1*q1 - 2*q2*q2, 0,
                           0, 0, 0, 1]);
    },
    toMat3x3:function()
    {
	// shorthand
	q0 = this.q0;
	q1 = this.q1;	
	q2 = this.q2;
	q3 = this.q3;

	// 3x3 rotation matrix
	return new Mat3x3([1 - 2*q2*q2 - 2*q3*q3, 2*q1*q2 - 2*q3*q0, 2*q1*q3 + 2*q2*q0,
                           2*q1*q2 + 2*q3*q0, 1 - 2*q1*q1 - 2*q3*q3, 2*q2*q3 - 2*q1*q0,
                           2*q1*q3 - 2*q2*q0, 2*q2*q3 + 2*q1*q0, 1 - 2*q1*q1 - 2*q2*q2]);
    }
}
