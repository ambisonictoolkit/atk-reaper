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
//	@@ { arg aNumber; ^Cartesian.new(this.x, this.y, aNumber) }
	@ { arg aValue;								// overload default method
		aValue.isKindOf(SimpleNumber).if(
			{ ^Cartesian.new(this.x, this.y, aValue) },
			{ ^Rect.fromPoints(this, aValue) }		// default SC
		)
		}
		
	// mirror
	mirrorX { ^x.neg @ y }
	mirrorY { ^x @ y.neg }
	mirrorO { ^x.neg @ y.neg }
	
}
