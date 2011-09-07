//---------------------------------------------------------------------
// Joseph Anderson 2011
//
// extension of Number to support spherical and cartesion (3d) coordinates
//
// NOTE: include license, etc before public release
//---------------------------------------------------------------------

+ SimpleNumber {

	// conversion
	asCartesian { ^Cartesian.new(this, this, this) }
}
