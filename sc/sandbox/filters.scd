// -------------------------------------------------------
// ATK (Sandbox)
//
// Coded by Joseph Anderson 2011
//
// -------------------------------------------------------

// The following filters are needed for the (FOA) ATK:
//
//	NFC					- Nearfield Compensation (note, currently named as Distance)
//	Proximity				- Proximity simulation
//	Psychoacoustic shelf	- built on RMShelf2
//
// This file implements a prototype for the psychoacoustic shelf filter


(
	~atkPsychoShelf = { arg bSig, freq, k;
		var k2;
		
		// expand k from degree (order) gains to channel gains
		k2 = k.collect({ arg item, i;
			Array.fill(2 * i + 1, {item})}).flat;

		RMShelf2.ar(bSig, freq, k2)
	};
)
