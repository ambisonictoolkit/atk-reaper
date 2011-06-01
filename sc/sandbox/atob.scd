// -------------------------------------------------------
// ATK (Sandbox)
//
// Coded by Joseph Anderson 2011
//
// -------------------------------------------------------


// AtoB, BtoA
//
// Implemented below as global functions, interface:
//
// ~a2b.value(aSig, orientation, weight)
// ~b2a.value(bSig, orientation, weight)
(
	var a2b0, a2b1, a2b2, a2b3, a2b4, a2b5, a2b6, a2b7;
	var b2a0, b2a1, b2a2, b2a3, b2a4, b2a5, b2a6, b2a7;

	var w0, w1, w2, w3;

	// some constants
	var recipSqrt2 = 2.sqrt.reciprocal;
	var negRecipSqrt2 = -1 * recipSqrt2;
	var sqrt3Div2 = 3.sqrt/2;
	var negSqrt3Div2 = -1 * sqrt3Div2;
	var sqrt3Div6 = 3.sqrt/6;
	var negSqrt3Div6 = -1 * sqrt3Div6;
	var sqrt6Div3 = 6.sqrt/3;
	var negSqrt6Div3 = -1 * sqrt6Div3;
	var recipSqrt6 = 6.sqrt.reciprocal;
	var negRecipSqrt6 = -1 * recipSqrt6;
	var sqrt6Div4 = 6.sqrt/4;
	
	var sqrt2divSqrt3 = (2/3).sqrt;
	var sqrt6 = 6.sqrt;
	var sqrt2 = 2.sqrt;
	
	// transform matricies
	// NOTE: The broadcast rules of SC3 are different from Numpy,
	// 		so that the below matricies appear to be inverted.
	//		Additionally, these matricies are inverted as to those
	//		appearing in "A2B matrix tests.sc"

	// 0 - orthogonal (front left up)
	// [ FLU, FRD, BLD, BRU ]
	b2a0 = [
		[ 0.5, 0.5, 0.5, 0.5 ],
		[ 0.5, 0.5, -0.5, -0.5 ],
		[ 0.5, -0.5, 0.5, -0.5 ],
		[ 0.5, -0.5, -0.5, 0.5 ]
	];
	a2b0 = b2a0;

	// 1 - front left down
	// [ FLD, FRU, BLU, BRD ]
	b2a1 = [
		[ 0.5, 0.5, 0.5, 0.5 ],
		[ 0.5, 0.5, -0.5, -0.5 ],
		[ 0.5, -0.5, 0.5, -0.5 ],
		[ -0.5, 0.5, 0.5, -0.5 ]
	];
	a2b1 = [
		[ 0.5, 0.5, 0.5, -0.5 ],
		[ 0.5, 0.5, -0.5, 0.5 ],
		[ 0.5, -0.5, 0.5, 0.5 ],
		[ 0.5, -0.5, -0.5, -0.5 ]
	];

	// 2 - front left-right
	// [ FL, FR, BU, BD ]
	b2a2 = [
		[ 0.5, 0.5, 0.5, 0.5 ],
		[ 0.5, 0.5, -0.5, -0.5 ],
		[ recipSqrt2, negRecipSqrt2, 0, 0 ],
		[ 0, 0, recipSqrt2, negRecipSqrt2 ]
	];
	a2b2 = [
		[ 0.5, 0.5, recipSqrt2, 0 ],
		[ 0.5, 0.5, negRecipSqrt2, 0 ],
		[ 0.5, -0.5, 0, recipSqrt2 ],
		[ 0.5, -0.5, 0, negRecipSqrt2 ]
	];

	// 3 - front up-down
	// [ FU, FD, BL, BR ]
	b2a3 = [
		[ 0.5, 0.5, 0.5, 0.5 ],
		[ 0.5, 0.5, -0.5, -0.5 ],
		[ 0, 0, recipSqrt2, negRecipSqrt2 ],
		[ recipSqrt2, negRecipSqrt2, 0, 0 ]
	];
	a2b3 = [
		[ 0.5, 0.5, 0, recipSqrt2 ],
		[ 0.5, 0.5, 0, negRecipSqrt2 ],
		[ 0.5, -0.5, recipSqrt2, 0 ],
		[ 0.5, -0.5, negRecipSqrt2, 0 ]
		
	];

	// 4 - front & back down
	// [ F, BD, BLU, BRU ]
	b2a4 = [
		[ 0.5, 0.5, 0.5, 0.5 ],
		[ sqrt3Div2, negSqrt3Div6, negSqrt3Div6, negSqrt3Div6 ],
		[ 0, 0, recipSqrt2, negRecipSqrt2 ],
		[ 0, negSqrt6Div3, recipSqrt6, recipSqrt6 ]
	];
	a2b4 = [
		[ 0.5, sqrt3Div2, 0, 0 ],
		[ 0.5, negSqrt3Div6, 0, negSqrt6Div3 ],
		[ 0.5, negSqrt3Div6, recipSqrt2, recipSqrt6 ],
		[ 0.5, negSqrt3Div6, negRecipSqrt2, recipSqrt6 ]
	];

	// 5 - front & back up
	// [ F, BU, BLD, BRD ]
	b2a5 = [
		[ 0.5, 0.5, 0.5, 0.5 ],
		[ sqrt3Div2, negSqrt3Div6, negSqrt3Div6, negSqrt3Div6 ],
		[ 0, 0, recipSqrt2, negRecipSqrt2 ],
		[ 0, sqrt6Div3, negRecipSqrt6, negRecipSqrt6 ]
	];
	a2b5 = [
		[ 0.5, sqrt3Div2, 0, 0 ],
		[ 0.5, negSqrt3Div6, 0, sqrt6Div3 ],
		[ 0.5, negSqrt3Div6, recipSqrt2, negRecipSqrt6 ],
		[ 0.5, negSqrt3Div6, negRecipSqrt2, negRecipSqrt6 ]
	];

	// 6 - front left-right up
	// [ FLU, FRU, FD, B ]
	b2a6 = [
		[ 0.5, 0.5, 0.5, 0.5 ],
		[ sqrt3Div6, sqrt3Div6, sqrt3Div6, negSqrt3Div2 ],
		[ recipSqrt2, negRecipSqrt2, 0, 0 ],
		[ recipSqrt6, recipSqrt6, negSqrt6Div3, 0 ]
	];
	a2b6 = [
		[ 0.5, sqrt3Div6, recipSqrt2, recipSqrt6 ],
		[ 0.5, sqrt3Div6, negRecipSqrt2, recipSqrt6 ],
		[ 0.5, sqrt3Div6, 0, negSqrt6Div3 ],
		[ 0.5, negSqrt3Div2, 0, 0 ]
	];

	// 7 - front left-right down
	// [ FLD, FRD, FU, B ]
	b2a7 = [
		[ 0.5, 0.5, 0.5, 0.5 ],
		[ sqrt3Div6, sqrt3Div6, sqrt3Div6, negSqrt3Div2 ],
		[ recipSqrt2, negRecipSqrt2, 0, 0 ],
		[ negRecipSqrt6, negRecipSqrt6, sqrt6Div3, 0 ]
	];
	a2b7 = [
		[ 0.5, sqrt3Div6, recipSqrt2, negRecipSqrt6 ],
		[ 0.5, sqrt3Div6, negRecipSqrt2, negRecipSqrt6 ],
		[ 0.5, sqrt3Div6, 0, sqrt6Div3 ],
		[ 0.5, negSqrt3Div2, 0, 0 ]
	];
	
	// weights
	w0 = 1; 				// canonical
	w1 = sqrt2divSqrt3;	// decorrelated (on the sphere)
	w2 = sqrt6;			// cardioid
	w3 = sqrt2;			// unscaled
	
	// functions
	~a2b = { arg aSig, orientation = 0, weight = 0;
		var bSig;
		
		bSig = switch(orientation,
			0, { a2b0 * aSig },
			1, { a2b1 * aSig },
			2, { a2b2 * aSig },
			3, { a2b3 * aSig },
			4, { a2b4 * aSig },
			5, { a2b5 * aSig },
			6, { a2b6 * aSig },
			7, { a2b7 * aSig }
		).sum;
		
		bSig = switch(weight,
			0, { bSig.put(0, w0.reciprocal * bSig.at(0)) },
			1, { bSig.put(0, w1.reciprocal * bSig.at(0)) },
			2, { bSig.put(0, w2.reciprocal * bSig.at(0)) },
			3, { bSig.put(0, w3.reciprocal * bSig.at(0)) }
		)
		
		
	};

	~b2a = { arg bSig, orientation = 0, weight = 0;
		var aSig;
		
		bSig = switch(weight,
			0, { bSig.put(0, w0 * bSig.at(0)) },
			1, { bSig.put(0, w1 * bSig.at(0)) },
			2, { bSig.put(0, w2 * bSig.at(0)) },
			3, { bSig.put(0, w3 * bSig.at(0)) }
		);

		aSig = switch(orientation,
			0, { b2a0 * bSig },
			1, { b2a1 * bSig },
			2, { b2a2 * bSig },
			3, { b2a3 * bSig },
			4, { b2a4 * bSig },
			5, { b2a5 * bSig },
			6, { b2a6 * bSig },
			7, { b2a7 * bSig }
		).sum
	}
)

