// -------------------------------------------------------
// ATK (Sandbox)
//
// Coded by Joseph Anderson 2011
//
// -------------------------------------------------------


// RTT (Rotate-Tilt-Tumble)
//
// Implemented below as global functions, interface:
//
// ~atkRtt.value(bSig, rotate, tilt, tumble)
(
	~atkRtt = { arg bSig, rotate = 0, tilt = 0, tumble = 0;
		var w, x, y, z;
		
		#w, x, y, z = bSig;
		
		#w, x, y, z = AtkRotate.ar(w, x, y, z, rotate);
		#w, x, y, z = AtkTilt.ar(w, x, y, z, tilt);
		AtkTumble.ar(w, x, y, z, tumble);
	}
)
