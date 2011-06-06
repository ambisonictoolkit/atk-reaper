// -------------------------------------------------------
// ATK (Sandbox)
//
// Coded by Joseph Anderson 2011
//
// -------------------------------------------------------


// Implemented below as global functions, interface:
//
// ~atkDominate.value(gain, azimuth, elevation)
// ~atkFocus.value(angle, azimuth, elevation)
//
// NOTE: ~atkSquishX is found in file 'directivity.scd'

(
	~atkDominate = { arg bSig, gain = 0, azimuth = 0, elevation = 0;
		var w, x, y, z;
		
		#w, x, y, z = bSig;
		
		#w, x, y, z = AtkRotate.ar(w, x, y, z, -1 * azimuth);
		#w, x, y, z = AtkTumble.ar(w, x, y, z, -1 * elevation);
		#w, x, y, z = AtkDominateX.ar(w, x, y, z, gain);
		#w, x, y, z = AtkTumble.ar(w, x, y, z, elevation);
		AtkRotate.ar(w, x, y, z, azimuth);
	};

	~atkFocus = { arg bSig, angle = 0, azimuth = 0, elevation = 0;
		var w, x, y, z;
		
		#w, x, y, z = bSig;
		
		#w, x, y, z = AtkRotate.ar(w, x, y, z, -1 * azimuth);
		#w, x, y, z = AtkTumble.ar(w, x, y, z, -1 * elevation);
		#w, x, y, z = AtkFocusX.ar(w, x, y, z, angle);
		#w, x, y, z = AtkTumble.ar(w, x, y, z, elevation);
		AtkRotate.ar(w, x, y, z, azimuth);
	};

	~atkPress = { arg bSig, angle = 0, azimuth = 0, elevation = 0;
		var w, x, y, z;
		
		#w, x, y, z = bSig;
		
		#w, x, y, z = AtkRotate.ar(w, x, y, z, -1 * azimuth);
		#w, x, y, z = AtkTumble.ar(w, x, y, z, -1 * elevation);
		#w, x, y, z = AtkPressX.ar(w, x, y, z, angle);
		#w, x, y, z = AtkTumble.ar(w, x, y, z, elevation);
		AtkRotate.ar(w, x, y, z, azimuth);
	};

	~atkPush = { arg bSig, angle = 0, azimuth = 0, elevation = 0;
		var w, x, y, z;
		
		#w, x, y, z = bSig;
		
		#w, x, y, z = AtkRotate.ar(w, x, y, z, -1 * azimuth);
		#w, x, y, z = AtkTumble.ar(w, x, y, z, -1 * elevation);
		#w, x, y, z = AtkPushX.ar(w, x, y, z, angle);
		#w, x, y, z = AtkTumble.ar(w, x, y, z, elevation);
		AtkRotate.ar(w, x, y, z, azimuth);
	};

	~atkZoom = { arg bSig, angle = 0, azimuth = 0, elevation = 0;
		var w, x, y, z;
		
		#w, x, y, z = bSig;
		
		#w, x, y, z = AtkRotate.ar(w, x, y, z, -1 * azimuth);
		#w, x, y, z = AtkTumble.ar(w, x, y, z, -1 * elevation);
		#w, x, y, z = AtkZoomX.ar(w, x, y, z, angle);
		#w, x, y, z = AtkTumble.ar(w, x, y, z, elevation);
		AtkRotate.ar(w, x, y, z, azimuth);
	}
	
	~atkMirror = { arg bSig, azimuth = 0, elevation = 0;
		var w, x, y, z;
		
		#w, x, y, z = bSig;
		
		#w, x, y, z = AtkRotate.ar(w, x, y, z, -1 * azimuth);
		#w, x, y, z = AtkTumble.ar(w, x, y, z, -1 * elevation);
		#w, x, y, z = AtkTumble.ar(w, -1 * x, y, z, elevation);
		AtkRotate.ar(w, x, y, z, azimuth);
	}

	~atkDirectD = { arg bSig, azimuth = 0, elevation = 0;
		var w, x, y, z;
		
		#w, x, y, z = bSig;
		
		#w, x, y, z = AtkRotate.ar(w, x, y, z, -1 * azimuth);
		#w, x, y, z = AtkTumble.ar(w, x, y, z, -1 * elevation);
		#w, x, y, z = ~atkSquishX.value([w, x, y, z], angle);
		#w, x, y, z = AtkTumble.ar(w, x, y, z, elevation);
		AtkRotate.ar(w, x, y, z, azimuth);
	}
	
)
