//---------------------------------------------------------------------
// Joseph Anderson 2011
//
// expansion of Number to support spherical coordinates
//
// NOTE: include license, etc before public release
//---------------------------------------------------------------------

+ Number {

//	performBinaryOpOnPoint { arg op, aPoint, adverb;
//		^Point.new(this.perform(op, aPoint.x, adverb), this.perform(op, aPoint.y, adverb));
//	}

	// spherical support
	phi { ^0.0 }

//	// complex support
//	real { ^this }
//	imag { ^0.0 }

	// conversion
	@@ { arg aPoint; ^Cartesian.new(this, aPoint.x, aPoint.y) }
//	complex { arg imaginaryPart; ^Complex.new(this, imaginaryPart) }
//	polar { arg angle; ^Polar.new(this, angle) }

}
