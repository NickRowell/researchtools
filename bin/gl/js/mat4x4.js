/**
 * JavaScript object representing a 4x4 matrix.
 */
function Mat4x4(els) 
{
	this.els = els;
}

Mat4x4.prototype = {
    constructor: Mat4x4,
    // Note transposition: assuming the Mat4x4 is in column-major format (for GL) then we're
    // extracting the 3x3 submatrix in row-major format in order to use it for geometry operations
    // in the conventional representation.
    getSub3x3:function()
    {
	return new Mat3x3([ this.els[0],  this.els[4], this.els[8], 
			            this.els[1],  this.els[5], this.els[9],
			            this.els[2],  this.els[6], this.els[10]]);
    },
    // Multiply one Mat4x4 by another
    multiply:function(B)
    {
        return new Mat4x4(
             [this.els[0]*B.els[0] + this.els[1]*B.els[4] + this.els[2]*B.els[8]  + this.els[3]*B.els[12],
              this.els[0]*B.els[1] + this.els[1]*B.els[5] + this.els[2]*B.els[9]  + this.els[3]*B.els[13],
              this.els[0]*B.els[2] + this.els[1]*B.els[6] + this.els[2]*B.els[10] + this.els[3]*B.els[14],
              this.els[0]*B.els[3] + this.els[1]*B.els[7] + this.els[2]*B.els[11] + this.els[3]*B.els[15],
             
              this.els[4]*B.els[0] + this.els[5]*B.els[4] + this.els[6]*B.els[8]  + this.els[7]*B.els[12],
              this.els[4]*B.els[1] + this.els[5]*B.els[5] + this.els[6]*B.els[9]  + this.els[7]*B.els[13],
              this.els[4]*B.els[2] + this.els[5]*B.els[6] + this.els[6]*B.els[10] + this.els[7]*B.els[14],
              this.els[4]*B.els[3] + this.els[5]*B.els[7] + this.els[6]*B.els[11] + this.els[7]*B.els[15],
             
              this.els[8]*B.els[0] + this.els[9]*B.els[4] + this.els[10]*B.els[8]  + this.els[11]*B.els[12],
              this.els[8]*B.els[1] + this.els[9]*B.els[5] + this.els[10]*B.els[9]  + this.els[11]*B.els[13],
              this.els[8]*B.els[2] + this.els[9]*B.els[6] + this.els[10]*B.els[10] + this.els[11]*B.els[14],
              this.els[8]*B.els[3] + this.els[9]*B.els[7] + this.els[10]*B.els[11] + this.els[11]*B.els[15],
                
              this.els[12]*B.els[0] + this.els[13]*B.els[4] + this.els[14]*B.els[8]  + this.els[15]*B.els[12],
              this.els[12]*B.els[1] + this.els[13]*B.els[5] + this.els[14]*B.els[9]  + this.els[15]*B.els[13],
              this.els[12]*B.els[2] + this.els[13]*B.els[6] + this.els[14]*B.els[10] + this.els[15]*B.els[14],
              this.els[12]*B.els[3] + this.els[13]*B.els[7] + this.els[14]*B.els[11] + this.els[15]*B.els[15]]);
    },
    transpose:function()
    {
	return new Mat4x4([this.els[0],  this.els[4], this.els[8],  this.els[12],
			           this.els[1],  this.els[5], this.els[9],  this.els[13],
			           this.els[2],  this.els[6], this.els[10], this.els[14],
			           this.els[3],  this.els[7], this.els[11], this.els[15]]);
    },
    getNormalsMatrix:function()
    {
    	// This function applies to Mat4x4 instances that represent GL ModelView matrices.
	// The 'normals matrix' in graphics is the inverse transpose of the model view matrix.
	// The calculation here exploits the special structure of the model view matrix, with
	// the upper left 3x3 submatrix containing the rotation part and the (3,0) (3,1) (3,2)
	// elements containing the translation part.

	// Extract rotation part of matrix and take transpose, which is equivalent
	// to the inverse because the rotation matrix is orthogonal.
        var Rt = this.getSub3x3().transpose();
        
        // Extract translation part of matrix
        var t = new Vec3(this.els[12], this.els[13], this.els[14]);
        
        // Inverse translation
        var RT = Rt.multiply(t);
        
        // Build Mat4x4
	return new Mat4x4([Rt.els[0], Rt.els[1], Rt.els[2], -RT.x,
                           Rt.els[3], Rt.els[4], Rt.els[4], -RT.y,
                           Rt.els[6], Rt.els[7], Rt.els[8], -RT.z,
                                 0.0,       0.0,       0.0,   1.0]);

   }
}








