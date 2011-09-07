//---------------------------------------------------------------------
// Joseph Anderson 2011
//
// extension of Number to support spherical and cartesion (3d) coordinates
//
// NOTE: include license, etc before public release
//---------------------------------------------------------------------

+ Number {

	performBinaryOpOnCartesian { arg op, aCartesian, adverb;
		^Cartesian.new(
			this.perform(op, aCartesian.x, adverb),
			this.perform(op, aCartesian.y, adverb),
			this.perform(op, aCartesian.z, adverb)
		);
	}

	// spherical support
	phi { ^0.0 }

	// conversion
	@ { arg aValue;								// overload default method
		aValue.isKindOf(SimpleNumber).if(
			{ ^Point.new(this, aValue) },			// default SC
			{ ^Cartesian.new(this, aValue.x, aValue.y) }
		)
		}
	asCartesian { ^Cartesian.new(this, this, this) }
	spherical { arg theta, phi; ^Spherical.new(this, theta, phi) }
}
