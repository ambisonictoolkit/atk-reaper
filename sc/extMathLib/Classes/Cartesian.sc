//---------------------------------------------------------------------
// Joseph Anderson 2011
//
// Cartesian class modelled on Point.sc
//
// NOTE: include license, etc before public release
//---------------------------------------------------------------------

Cartesian {
	var <>x = 0, <>y = 0, <>z = 0;

	*new { arg x=0, y=0, z=0;
		^super.newCopyArgs(x, y, z);
	}

	set { arg argX=0, argY=0, argZ=0; x = argX; y = argY; z = argZ;}

	asCartesian { ^this }
	asPoint { ^Point.new(x,y) }					// implemented as a projection
	asComplex { ^Complex.new(x,y) }				// implemented as a projection
											// add tricomplex
	asPolar { ^Polar.new(this.rho, this.theta) }	// implemented as a projection
											// add spherical
	asRect { ^Rect.new(0,0,x,y) }				// implemented as a projection
	asArray { ^[this.x, this.y, this.z] }

	== { arg aCartesian;
		^aCartesian respondsTo: #[\x, \y, \z] and:
			{ x == aCartesian.x and:
				{ y == aCartesian.y and:
					{ z == aCartesian.z }
				}
			}
	}
	hash { ^ ((x.hash << 1) bitXor: y.hash) bitXor: z.hash }

//	+ { arg delta;
//		var deltaPoint;
//		deltaPoint = delta.asPoint;
//		^(this.x + deltaPoint.x) @ (this.y + deltaPoint.y)
//	}
//	- { arg delta;
//		var deltaPoint;
//		deltaPoint = delta.asPoint;
//		^(this.x - deltaPoint.x) @ (this.y - deltaPoint.y)
//	}
//
//	* { arg scale;
//		var scalePoint;
//		scalePoint = scale.asPoint;
//		^(this.x * scalePoint.x) @ (this.y * scalePoint.y)
//	}
//	/ { arg scale;
//		var scalePoint;
//		scalePoint = scale.asPoint;
//		^(this.x / scalePoint.x) @ (this.y / scalePoint.y)
//	}
//	div { arg scale;
//		var scalePoint;
//		scalePoint = scale.asPoint;
//		^(this.x div: scalePoint.x) @ (this.y div: scalePoint.y)
//	}
//	translate { arg delta;
//		^(this.x + delta.x) @ (this.y + delta.y)
//	}
//	scale { arg scale;
//		^(this.x * scale.x) @ (this.y * scale.y)
//	}
//	rotate { arg angle; // in radians
//		var sinr, cosr;
//		sinr = angle.sin;
//		cosr = angle.cos;
//		^((x * cosr) - (y * sinr)) @ ((y * cosr) + (x * sinr))
//	}
//
	abs { ^(x.abs @ y.abs) @@ z.abs }

	rho { ^(x.squared + y.squared + z.squared).sqrt }
	theta { ^atan2(y, x) }
	phi { ^atan2(z, (x.squared + y.squared).sqrt) }

//	dist { arg aPoint;
//		aPoint = aPoint.asPoint;
//		^hypot(x - aPoint.x, y - aPoint.y)
//	}
	transpose { ^(y @ x) @@ z }
	transposeXY { ^(y @ x) @@ z }
	transposeYZ { ^(x @ z) @@ y }
	transposeXZ { ^(z @ y) @@ x }

//	round { arg quant;
//		quant = quant.asPoint;
//		^x.round(quant.x) @ y.round(quant.y)
//	}
//	trunc { arg quant;
//		quant = quant.asPoint;
//		^x.trunc(quant.x) @ y.trunc(quant.y)
//	}
//
//	mod {|that|
//		var thatPoint;
//		thatPoint = that.asPoint;
//		^(this.x mod: thatPoint.x) @ (this.y mod: thatPoint.y)
//	}
//
	printOn { arg stream;
		stream << this.class.name << "( " << x << ", " << y << ", " << z << " )";
	}
	storeArgs { ^[x,y,z] }
}


//PointArray : Point
//{
//	*new { arg n;
//		^super.new(Signal.new(n), Signal.new(n))
//	}
//	add { arg point;
//		x = x.add(point.x);
//		y = y.add(point.y);
//	}
//}