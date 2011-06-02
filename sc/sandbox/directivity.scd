// -------------------------------------------------------
// ATK (Sandbox)
//
// Coded by Joseph Anderson 2011
//
// -------------------------------------------------------


// directivity transforms
//
// Implemented below as global functions, interface:
//
// ~atkDirect.value(bSig, angle)
(
	~atkDirect = { arg bSig, angle;
		var w, x, y, z, g0, g1;
		
		#w, x, y, z = bSig;
		
		g0 = (1 + sin(angle)).sqrt;
		g1 = (1 - sin(angle)).sqrt;
		
		[g0 *w, g1 * x, g1 * y, g1 * z]
	};

	~atkSquishX = { arg bSig, angle;
		var w, x, y, z, g0, g1;
		
		#w, x, y, z = bSig;
		
		g0 = (1 + sin(angle)).sqrt;
		g1 = (1 - sin(angle)).sqrt;
		
		[g0 *w, g1 * x, g0 * y, g0 * z]
	};

	~atkSquishY = { arg bSig, angle;
		var w, x, y, z, g0, g1;
		
		#w, x, y, z = bSig;
		
		g0 = (1 + sin(angle)).sqrt;
		g1 = (1 - sin(angle)).sqrt;
		
		[g0 *w, g0 * x, g1 * y, g0 * z]
	};

	~atkSquishZ = { arg bSig, angle;
		var w, x, y, z, g0, g1;
		
		#w, x, y, z = bSig;
		
		g0 = (1 + sin(angle)).sqrt;
		g1 = (1 - sin(angle)).sqrt;
		
		[g0 *w, g0 * x, g0 * y, g1 * z]
	};
)
