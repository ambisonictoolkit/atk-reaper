//---------------------------------------------------------------------
// Joseph Anderson 2011
//
// expansion of Point to support spherical coordinates
//
// NOTE: include license, etc before public release
//---------------------------------------------------------------------

+ Point {

	// spherical support
	phi { ^0.0 }

//	// complex support
//	real { ^this }
//	imag { ^0.0 }

	// conversion
	asCartesian { ^Cartesian.new(this.x, this.y, 0) }
	@@ { arg aNumber; ^Cartesian.new(this.x, this.y, aNumber) }
}
