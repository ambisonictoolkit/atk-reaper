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
//	@@ { arg aPoint; ^Cartesian.new(this, aPoint.x, aPoint.y) }
	@ { arg aValue;								// overload default method
		aValue.isKindOf(SimpleNumber).if(
			{ ^Point.new(this, aValue) },			// default SC
			{ ^Cartesian.new(this, aValue.x, aValue.y) }
		)
		}
	asCartesian { ^Cartesian.new(this, this, this) }
//	complex { arg imaginaryPart; ^Complex.new(this, imaginaryPart) }
//	polar { arg angle; ^Polar.new(this, angle) }

}
