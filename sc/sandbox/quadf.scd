// -------------------------------------------------------
// ATK (Sandbox)
//
// Coded by Joseph Anderson 2011
//
// -------------------------------------------------------


// QuadF (quad decoder)
//
// Implemented below as global functions, interface:
//
// ~atkQuadF.value(bSig, angle = pi/4, directivity = 1)
//
// NOTE: angle is the 1/2 angle of the front (and back) pair
//		The default, pi/4, gives a square.
//		Also, unlike ATKPantoF, output order is in 'stereo pairs'
//		rather than counter-clockwise.
(
	~atkQuadF = { arg bSig, angle = pi/4, directivity = 1;
		var w, x, y, z, x2, y2, g, fl, fr, bl, br;
		
		#w, x, y, z = bSig;

		g = 2**((1 - directivity) / 2);
		
		x2 = g * sin(angle) * x;
		y2 = g * cos(angle) * y;

		fl = w + x2 + y2;
		fr = w + x2 - y2;
		bl = w - x2 + y2;
		br = w - x2 - y2;

		[ fl, fr, bl, br ]		
	}
)
