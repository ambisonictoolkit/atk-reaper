//---------------------------------------------------------------------
// Joseph Anderson 2011
//
// extension of Polar to support spherical and cartesian (3d) coordinates
//
// NOTE: include license, etc before public release
//---------------------------------------------------------------------

+ Polar {

	// spherical support
	phi { ^0.0 }
	angles { ^[ this.theta, 0.0 ] }
	phases { ^[ this.theta, 0.0 ] }
	
	// Point, Cartesian
	x { this.real }
	y { this.imag }
	z { ^0 }

//	// conversion
//	asCartesian { ^Cartesian.new(this.x, this.y, 0) }
////	@@ { arg aNumber; ^Cartesian.new(this.x, this.y, aNumber) }
//	@ { arg aValue;								// overload default method
//		aValue.isNumber.if(
//			{ ^Cartesian.new(this.x, this.y, aValue) },
//			{ ^Rect.fromPoints(this, aValue) }		// default SC
//		)
//		}
		
	// mirror
	mirrorX { ^this.asPoint.mirrorX.asPolar }
	mirrorY { ^this.asPoint.mirrorY.asPolar }
	mirrorZ { ^this }
	mirrorO { ^this.neg }
	
}
