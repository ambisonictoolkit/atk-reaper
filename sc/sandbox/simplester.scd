// -------------------------------------------------------
// ATK (Sandbox)
//
// Coded by Joseph Anderson 2011
//
// -------------------------------------------------------


// simple stereo
//
// Implemented below as global functions, interface:
//
// ~atkSimpleStereo.value(bSig, angle)
// NOTE: angle is framed as the 'distortion' angle rather than width,
//		so that behaviour matches that of AtkFocusX.
//		Angle = 0 gives left and right at 'hard' left and right
(
	~atkSimpleStereo = { arg sterSig, angle = pi/4;
		var l, r, w, x, y, z;
		
		#l, r = sterSig;
		
		[ 0.5 * (l + r), 1/2.sqrt * sin(angle) * (l + r),
			1/2.sqrt * cos(angle) * (l - r), Silent.ar(1) ]
	}
)
