// -------------------------------------------------------
// ATK (Sandbox)
//
// Coded by Joseph Anderson 2011
//
// -------------------------------------------------------


// mirror transforms (reflection across axes)
//
// Implemented below as global functions, interface:
//
// AtkMirror is 'directable' with interface:
// ~atkMirror.value(bSig, azimuth, elevation)
// Mirror
(
	~atkMirrorX = { arg bSig;
		var w, x, y, z;
		
		#w, x, y, z = bSig;
		
		[w, -1 * x, y, z]
	};

	~atkMirrorY = { arg bSig;
		var w, x, y, z;
		
		#w, x, y, z = bSig;
		
		[w, x, -1 * y, z]
	};

	~atkMirrorZ = { arg bSig;
		var w, x, y, z;
		
		#w, x, y, z = bSig;
		
		[w, x, y, -1 * z]
	};

	~atkMirrorO = { arg bSig;
		var w, x, y, z;
		
		#w, x, y, z = bSig;
		
		[w, -1 * x, -1 * y, -1 * z]
	};
)
