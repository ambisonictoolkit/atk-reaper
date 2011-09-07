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
	asPolar { ^Polar.new(this.rho, this.theta) }	// implemented as a projection
	asSpherical { ^Spherical.new(this.rho, this.theta, this.phi) }	asRect { ^Rect.new(0,0,x,y) }				// implemented as a projection
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

	+ { arg delta;
		var deltaCart;
		deltaCart = delta.asCartesian;
		^(this.x + deltaCart.x) @ (this.y + deltaCart.y) @ (this.z + deltaCart.z)
	}
	- { arg delta;
		var deltaCart;
		deltaCart = delta.asCartesian;
		^(this.x - deltaCart.x) @ (this.y - deltaCart.y) @ (this.z - deltaCart.z)
	}

	* { arg scale;
		var scaleCart;
		scaleCart = scale.asCartesian;
		^(this.x * scaleCart.x) @ (this.y * scaleCart.y) @ (this.z * scaleCart.z)
	}
	/ { arg scale;
		var scaleCart;
		scaleCart = scale.asCartesian;
		^(this.x / scaleCart.x) @ (this.y / scaleCart.y) @ (this.z / scaleCart.z)
	}
	div { arg scale;
		var scaleCart;
		scaleCart = scale.asCartesian;
		^(this.x div: scaleCart.x) @ (this.y div: scaleCart.y) @ (this.z div: scaleCart.z)
	}
	translate { arg delta;
		^(this.x + delta.x) @ (this.y + delta.y) @ (this.z + delta.z)
	}
	scale { arg scale;
		^(this.x * scale.x) @ (this.y * scale.y) @ (this.z * scale.z)
	}
	rotate { arg angle; // XY-plane, in radians
		var sinr, cosr;
		sinr = angle.sin;
		cosr = angle.cos;
		^((x * cosr) - (y * sinr)) @ ((x * sinr) + (y * cosr)) @ z
	}
	tilt { arg angle; // YZ-plane, in radians
		var sinr, cosr;
		sinr = angle.sin;
		cosr = angle.cos;
		^x @ ((y * cosr) - (z * sinr)) @ ((y * sinr) + (z * cosr))
	}
	tumble { arg angle; // XZ-plane, in radians
		var sinr, cosr;
		sinr = angle.sin;
		cosr = angle.cos;
		^((x * cosr) - (z * sinr)) @ y @ ((x * sinr) + (z * cosr))
	}
	rotateXY { arg angle; // XY-plane, in radians
		var sinr, cosr;
		sinr = angle.sin;
		cosr = angle.cos;
		^((x * cosr) - (y * sinr)) @ ((x * sinr) + (y * cosr)) @ z
	}
	rotateYZ { arg angle; // YZ-plane, in radians
		var sinr, cosr;
		sinr = angle.sin;
		cosr = angle.cos;
		^x @ ((y * cosr) - (z * sinr)) @ ((y * sinr) + (z * cosr))
	}
	rotateXZ { arg angle; // XZ-plane, in radians
		var sinr, cosr;
		sinr = angle.sin;
		cosr = angle.cos;
		^((x * cosr) - (z * sinr)) @ y @ ((x * sinr) + (z * cosr))
	}

	mirrorX { ^x.neg @ y @ z }
	mirrorY { ^x @ y.neg @ z }
	mirrorZ { ^x @ y @ z.neg }
	mirrorO { ^x.neg @ y.neg @ z.neg }

	abs { ^x.abs @ y.abs @ z.abs }

	rho { ^(x.squared + y.squared + z.squared).sqrt }
	theta { ^atan2(y, x) }
	phi { ^atan2(z, (x.squared + y.squared).sqrt) }

	dist { arg aCart;
		aCart = aCart.asCartesian;
		^((x - aCart.x).squared + (y - aCart.y).squared + (z - aCart.z).squared).sqrt
	}
	transpose { ^y @ x @ z }				// same as Point
	transposeXY { ^y @ x @ z }
	transposeYZ { ^x @ z @ y }
	transposeXZ { ^z @ y @ x }

	round { arg quant;
		quant = quant.asCartesian;
		^x.round(quant.x) @ y.round(quant.y) @ z.round(quant.z)
	}
	trunc { arg quant;
		quant = quant.asCartesian;
		^x.trunc(quant.x) @ y.trunc(quant.y) @ z.trunc(quant.z)
	}

	mod {|that|
		var thatCart;
		thatCart = that.asCartesian;
		^(this.x mod: thatCart.x) @ (this.y mod: thatCart.y) @ (this.z mod: thatCart.z)
	}

	printOn { arg stream;
		stream << this.class.name << "( " << x << ", " << y << ", " << z << " )";
	}
	storeArgs { ^[x,y,z] }
}


CartesianArray : Cartesian
{
	*new { arg n;
		^super.new(Signal.new(n), Signal.new(n), Signal.new(n))
	}
	add { arg cart;
		x = x.add(cart.x);
		y = y.add(cart.y);
		z = z.add(cart.z);
	}
}