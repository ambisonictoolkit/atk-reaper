//---------------------------------------------------------------------
// Joseph Anderson 2011
//
// expansion of Polar to support spherical coordinates
//
// NOTE: include license, etc before public release
//---------------------------------------------------------------------

+ Polar {

	// spherical support
	phi { ^0.0 }
	
	// Point, Cartesian
	x { this.real }
	y { this.imag }
	z { ^0 }

//	// complex support
//	real { ^this }
//	imag { ^0.0 }

//	// conversion
//	asCartesian { ^Cartesian.new(this.x, this.y, 0) }
////	@@ { arg aNumber; ^Cartesian.new(this.x, this.y, aNumber) }
//	@ { arg aValue;								// overload default method
//		aValue.isNumber.if(
//			{ ^Cartesian.new(this.x, this.y, aValue) },
//			{ ^Rect.fromPoints(this, aValue) }		// default SC
//		)
//		}
		
//	// mirror
//	mirrorX { ^x.neg @ y }
//	mirrorY { ^x @ y.neg }
//	mirrorO { ^x.neg @ y.neg }
	
}
