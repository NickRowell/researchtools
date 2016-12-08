/**
 * Javascript object representing a 3-vector.
 */
function Vec3(x, y, z) 
{
	this.x = x;
	this.y = y;
	this.z = z;
}

    Vec3.prototype = {
    constructor: Vec3,
    minus:function(b)
    {
    	return new Vec3(this.x-b.x, this.y-b.y, this.z-b.z);
    },
    plusEquals:function(b)
    {
    	this.x += b.x;
    	this.y += b.y;
    	this.z += b.z;
    },
    cross:function(b)
    {
    	return new Vec3((this.y * b.z - this.z * b.y),
                        (this.z * b.x - this.x * b.z),
                        (this.x * b.y - this.y * b.x));
    },
    dot:function(b)
    {
    	return this.x*b.x + this.y*b.y + this.z*b.z;
    },
    norm2:function()
    {
    	return this.x*this.x + this.y*this.y + this.z*this.z;
    },
    norm:function()
    {
    	return Math.sqrt(this.norm2());
    },
    normalise:function()
    {
    	var A = this.norm();
    	return new Vec3(this.x/A, this.y/A, this.z/A);
    },
    get:function(i)
    {
    	if(i==0) {
    		return this.x;
    	}
    	if(i==1) {
    		return this.y;
    	}
    	if(i==2) {
    		return this.z;
    	}
    }

}

// Computes the surface normal for the triangle formed by the three
// vertices r0, r1, r2. Clockwise winding order is assumed.
function getClockwiseSurfaceNormal(r0, r1, r2)
{
	// Vector from v0 to v1
	a = r1.minus(r0);

	// Vector from v0 to v2
	b = r2.minus(r0);

	// Cross product a x b gives normal direction
	n = a.cross(b);

	return n.normalise();
}
