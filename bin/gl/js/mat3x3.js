/**
 * JavaScript object representing a 3x3 matrix.
 */
function Mat3x3(els) 
{
	this.els = els;
}

Mat3x3.prototype = {
    constructor: Mat3x3,
    transpose:function ()  
    {
        return new Mat3x3([ this.els[0],  this.els[3], this.els[6], 
			                this.els[1],  this.els[4], this.els[7],
			                this.els[2],  this.els[5], this.els[8]]);
    },
    get:function(i,j)
    {
    	return this.els[i*3+j]
    },
    multiply:function(vec3)
    {
	var x = this.els[0]*vec3.x + this.els[1]*vec3.y + this.els[2]*vec3.z;
	var y = this.els[3]*vec3.x + this.els[4]*vec3.y + this.els[5]*vec3.z;
	var z = this.els[6]*vec3.x + this.els[7]*vec3.y + this.els[8]*vec3.z;

	return new Vec3(x, y, z);
    }

}
