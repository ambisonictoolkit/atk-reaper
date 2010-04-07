// A2B, B2A matrix tests

// This first set does not include W scaling
// e.g, scale for canonical ('can') 2.sqrt.reciprocal
(
	var a2b0, a2b1, a2b2, a2b3, a2b4, a2b5, a2b6, a2b7;
	var b2a0, b2a1, b2a2, b2a3, b2a4, b2a5, b2a6, b2a7;

	var w0;

	var transform, btest, atest;

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
	
	// matrix transform
	transform = { arg testInMatrix, transformMatrix;
		var testOutMatrix;
		testOutMatrix = [];

		transformMatrix.do({ arg item, i;
			testOutMatrix = testOutMatrix.add((item * testInMatrix).sum)
		});
		
		testOutMatrix
	};


	// test matrix input
//	btest = [ 2.sqrt.reciprocal, 1, 0, 0 ];		// front
//	btest = [ 2.sqrt.reciprocal, 1.neg, 0, 0 ];		// back
//	btest = [ 2.sqrt.reciprocal, 2.sqrt.reciprocal, 2.sqrt.reciprocal, 0 ]; // front left
//	btest = [ 2.sqrt.reciprocal, 2.sqrt.reciprocal, 2.sqrt.reciprocal.neg, 0 ]; // front right
//	btest = [ 2.sqrt.reciprocal, 2.sqrt.reciprocal.neg, 2.sqrt.reciprocal, 0 ]; // back left
//	btest = [ 2.sqrt.reciprocal, 2.sqrt.reciprocal.neg, 2.sqrt.reciprocal.neg, 0 ]; // back right
//	btest = [ 2.sqrt.reciprocal, 0, 1, 0 ];		// left
//	btest = [ 2.sqrt.reciprocal, 0, 0 , 1 ];		// up
//	btest = [ 2.sqrt.reciprocal, 2.sqrt.reciprocal.neg, 0, 2.sqrt.reciprocal.neg ]; // back down
//	btest = [ 1, 0, 0, 0 ];
//	btest = [ 0, 1, 0, 0 ];
	btest = [ 0, 0, 1, 0 ];
//	btest = [ 0, 0, 0, 1 ];
//	btest = [ 1, 1, 1, 1 ];
//	btest = [ 1, 1, -1, 1 ];
	atest = [];


	//transform matricies

	// 0 - orthogonal (front left up)
	// [ FLU, FRD, BLD, BRU ]
	a2b0 = [
		[ 0.5, 0.5, 0.5, 0.5 ],
		[ 0.5, 0.5, -0.5, -0.5 ],
		[ 0.5, -0.5, 0.5, -0.5 ],
		[ 0.5, -0.5, -0.5, 0.5 ]
	];
	b2a0 = a2b0;

	// 1 - front left down
	// [ FLD, FRU, BLU, BRD ]
	a2b1 = [
		[ 0.5, 0.5, 0.5, 0.5 ],
		[ 0.5, 0.5, -0.5, -0.5 ],
		[ 0.5, -0.5, 0.5, -0.5 ],
		[ -0.5, 0.5, 0.5, -0.5 ]
	];
	b2a1 = [
		[ 0.5, 0.5, 0.5, -0.5 ],
		[ 0.5, 0.5, -0.5, 0.5 ],
		[ 0.5, -0.5, 0.5, 0.5 ],
		[ 0.5, -0.5, -0.5, -0.5 ]
	];

	// 2 - front left-right
	// [ FL, FR, BU, BD ]
	a2b2 = [
		[ 0.5, 0.5, 0.5, 0.5 ],
		[ 0.5, 0.5, -0.5, -0.5 ],
		[ recipSqrt2, negRecipSqrt2, 0, 0 ],
		[ 0, 0, recipSqrt2, negRecipSqrt2 ]
	];
	b2a2 = [
		[ 0.5, 0.5, recipSqrt2, 0 ],
		[ 0.5, 0.5, negRecipSqrt2, 0 ],
		[ 0.5, -0.5, 0, recipSqrt2 ],
		[ 0.5, -0.5, 0, negRecipSqrt2 ]
	];

	// 3 - front up-down
	// [ FU, FD, BL, BR ]
	a2b3 = [
		[ 0.5, 0.5, 0.5, 0.5 ],
		[ 0.5, 0.5, -0.5, -0.5 ],
		[ 0, 0, recipSqrt2, negRecipSqrt2 ],
		[ recipSqrt2, negRecipSqrt2, 0, 0 ]
	];
	b2a3 = [
		[ 0.5, 0.5, 0, recipSqrt2 ],
		[ 0.5, 0.5, 0, negRecipSqrt2 ],
		[ 0.5, -0.5, recipSqrt2, 0 ],
		[ 0.5, -0.5, negRecipSqrt2, 0 ]
		
	];

	// 4 - front & back down
	// [ F, BD, BLU, BRU ]
	a2b4 = [
		[ 0.5, 0.5, 0.5, 0.5 ],
		[ sqrt3Div2, negSqrt3Div6, negSqrt3Div6, negSqrt3Div6 ],
		[ 0, 0, recipSqrt2, negRecipSqrt2 ],
		[ 0, negSqrt6Div3, recipSqrt6, recipSqrt6 ]
	];
	b2a4 = [
		[ 0.5, sqrt3Div2, 0, 0 ],
		[ 0.5, negSqrt3Div6, 0, negSqrt6Div3 ],
		[ 0.5, negSqrt3Div6, recipSqrt2, recipSqrt6 ],
		[ 0.5, negSqrt3Div6, negRecipSqrt2, recipSqrt6 ]
	];

	// 5 - front & back up
	// [ F, BU, BLD, BRD ]
	a2b5 = [
		[ 0.5, 0.5, 0.5, 0.5 ],
		[ sqrt3Div2, negSqrt3Div6, negSqrt3Div6, negSqrt3Div6 ],
		[ 0, 0, recipSqrt2, negRecipSqrt2 ],
		[ 0, sqrt6Div3, negRecipSqrt6, negRecipSqrt6 ]
	];
	b2a5 = [
		[ 0.5, sqrt3Div2, 0, 0 ],
		[ 0.5, negSqrt3Div6, 0, sqrt6Div3 ],
		[ 0.5, negSqrt3Div6, recipSqrt2, negRecipSqrt6 ],
		[ 0.5, negSqrt3Div6, negRecipSqrt2, negRecipSqrt6 ]
	];

	// 6 - front left-right up
	// [ FLU, FRU, FD, B ]
	a2b6 = [
		[ 0.5, 0.5, 0.5, 0.5 ],
		[ sqrt3Div6, sqrt3Div6, sqrt3Div6, negSqrt3Div2 ],
		[ recipSqrt2, negRecipSqrt2, 0, 0 ],
		[ recipSqrt6, recipSqrt6, negSqrt6Div3, 0 ]
	];
	b2a6 = [
		[ 0.5, sqrt3Div6, recipSqrt2, recipSqrt6 ],
		[ 0.5, sqrt3Div6, negRecipSqrt2, recipSqrt6 ],
		[ 0.5, sqrt3Div6, 0, negSqrt6Div3 ],
		[ 0.5, negSqrt3Div2, 0, 0 ]
	];

	// 7 - front left-right down
	// [ FLD, FRD, FU, B ]
	a2b7 = [
		[ 0.5, 0.5, 0.5, 0.5 ],
		[ sqrt3Div6, sqrt3Div6, sqrt3Div6, negSqrt3Div2 ],
		[ recipSqrt2, negRecipSqrt2, 0, 0 ],
		[ negRecipSqrt6, negRecipSqrt6, sqrt6Div3, 0 ]
	];
	b2a7 = [
		[ 0.5, sqrt3Div6, recipSqrt2, negRecipSqrt6 ],
		[ 0.5, sqrt3Div6, negRecipSqrt2, negRecipSqrt6 ],
		[ 0.5, sqrt3Div6, 0, sqrt6Div3 ],
		[ 0.5, negSqrt3Div2, 0, 0 ]
	];

	// display input btest
	btest.postln;

	// b2a transform
	atest = transform.value(btest, b2a6);
	atest.postln;

	// a2b transform
	btest = transform.value(atest, a2b6);
)



// Scale for decorrelated ('dec') 3.sqrt.reciprocal
(
	var a2b0, a2b1, a2b2, a2b3, a2b4, a2b5, a2b6, a2b7;
	var b2a0, b2a1, b2a2, b2a3, b2a4, b2a5, b2a6, b2a7;

	var w0;

	var transform, btest, atest;

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
	
	// matrix transform
	transform = { arg testInMatrix, transformMatrix;
		var testOutMatrix;
		testOutMatrix = [];

		transformMatrix.do({ arg item, i;
			testOutMatrix = testOutMatrix.add((item * testInMatrix).sum)
		});
		
		testOutMatrix
	};


	// test matrix input
//	btest = [ 2.sqrt.reciprocal, 1, 0, 0 ];		// front
//	btest = [ 2.sqrt.reciprocal, 1.neg, 0, 0 ];		// back
//	btest = [ 2.sqrt.reciprocal, 2.sqrt.reciprocal, 2.sqrt.reciprocal, 0 ]; // front left
//	btest = [ 2.sqrt.reciprocal, 2.sqrt.reciprocal, 2.sqrt.reciprocal.neg, 0 ]; // front right
//	btest = [ 2.sqrt.reciprocal, 2.sqrt.reciprocal.neg, 2.sqrt.reciprocal, 0 ]; // back left
//	btest = [ 2.sqrt.reciprocal, 2.sqrt.reciprocal.neg, 2.sqrt.reciprocal.neg, 0 ]; // back right
//	btest = [ 2.sqrt.reciprocal, 0, 1, 0 ];		// left
//	btest = [ 2.sqrt.reciprocal, 0, 0 , 1 ];		// up
//	btest = [ 2.sqrt.reciprocal, 2.sqrt.reciprocal.neg, 0, 2.sqrt.reciprocal.neg ]; // back down
//	btest = [ 1, 0, 0, 0 ];
//	btest = [ 0, 1, 0, 0 ];
//	btest = [ 0, 0, 1, 0 ];
//	btest = [ 0, 0, 0, 1 ];
//	btest = [ 1, 1, 1, 1 ];
	btest = [ 1, 1, -1, 1 ];
	atest = [];


	//transform matricies

	// 0 - orthogonal (front left up)
	// [ FLU, FRD, BLD, BRU ]
	a2b0 = [
		[ sqrt6Div4, sqrt6Div4, sqrt6Div4, sqrt6Div4 ],
		[ 0.5, 0.5, -0.5, -0.5 ],
		[ 0.5, -0.5, 0.5, -0.5 ],
		[ 0.5, -0.5, -0.5, 0.5 ]
	];
	b2a0 = [
		[ recipSqrt6, 0.5, 0.5, 0.5 ],
		[ recipSqrt6, 0.5, -0.5, -0.5 ],
		[ recipSqrt6, -0.5, 0.5, -0.5 ],
		[ recipSqrt6, -0.5, -0.5, 0.5 ]
	];

	// 1 - front left down
	// [ FLD, FRU, BLU, BRD ]
	a2b1 = [
		[ sqrt6Div4, sqrt6Div4, sqrt6Div4, sqrt6Div4 ],
		[ 0.5, 0.5, -0.5, -0.5 ],
		[ 0.5, -0.5, 0.5, -0.5 ],
		[ -0.5, 0.5, 0.5, -0.5 ]
	];
	b2a1 = [
		[ recipSqrt6, 0.5, 0.5, -0.5 ],
		[ recipSqrt6, 0.5, -0.5, 0.5 ],
		[ recipSqrt6, -0.5, 0.5, 0.5 ],
		[ recipSqrt6, -0.5, -0.5, -0.5 ]
	];

	// 2 - front left-right
	// [ FL, FR, BU, BD ]
	a2b2 = [
		[ sqrt6Div4, sqrt6Div4, sqrt6Div4, sqrt6Div4 ],
		[ 0.5, 0.5, -0.5, -0.5 ],
		[ recipSqrt2, negRecipSqrt2, 0, 0 ],
		[ 0, 0, recipSqrt2, negRecipSqrt2 ]
	];
	b2a2 = [
		[ recipSqrt6, 0.5, recipSqrt2, 0 ],
		[ recipSqrt6, 0.5, negRecipSqrt2, 0 ],
		[ recipSqrt6, -0.5, 0, recipSqrt2 ],
		[ recipSqrt6, -0.5, 0, negRecipSqrt2 ]
	];

	// 3 - front up-down
	// [ FU, FD, BL, BR ]
	a2b3 = [
		[ sqrt6Div4, sqrt6Div4, sqrt6Div4, sqrt6Div4 ],
		[ 0.5, 0.5, -0.5, -0.5 ],
		[ 0, 0, recipSqrt2, negRecipSqrt2 ],
		[ recipSqrt2, negRecipSqrt2, 0, 0 ]
	];
	b2a3 = [
		[ recipSqrt6, 0.5, 0, recipSqrt2 ],
		[ recipSqrt6, 0.5, 0, negRecipSqrt2 ],
		[ recipSqrt6, -0.5, recipSqrt2, 0 ],
		[ recipSqrt6, -0.5, negRecipSqrt2, 0 ]
		
	];

	// 4 - front & back down
	// [ F, BD, BLU, BRU ]
	a2b4 = [
		[ sqrt6Div4, sqrt6Div4, sqrt6Div4, sqrt6Div4 ],
		[ sqrt3Div2, negSqrt3Div6, negSqrt3Div6, negSqrt3Div6 ],
		[ 0, 0, recipSqrt2, negRecipSqrt2 ],
		[ 0, negSqrt6Div3, recipSqrt6, recipSqrt6 ]
	];
	b2a4 = [
		[ recipSqrt6, sqrt3Div2, 0, 0 ],
		[ recipSqrt6, negSqrt3Div6, 0, negSqrt6Div3 ],
		[ recipSqrt6, negSqrt3Div6, recipSqrt2, recipSqrt6 ],
		[ recipSqrt6, negSqrt3Div6, negRecipSqrt2, recipSqrt6 ]
	];

	// 5 - front & back up
	// [ F, BU, BLD, BRD ]
	a2b5 = [
		[ sqrt6Div4, sqrt6Div4, sqrt6Div4, sqrt6Div4 ],
		[ sqrt3Div2, negSqrt3Div6, negSqrt3Div6, negSqrt3Div6 ],
		[ 0, 0, recipSqrt2, negRecipSqrt2 ],
		[ 0, sqrt6Div3, negRecipSqrt6, negRecipSqrt6 ]
	];
	b2a5 = [
		[ recipSqrt6, sqrt3Div2, 0, 0 ],
		[ recipSqrt6, negSqrt3Div6, 0, sqrt6Div3 ],
		[ recipSqrt6, negSqrt3Div6, recipSqrt2, negRecipSqrt6 ],
		[ recipSqrt6, negSqrt3Div6, negRecipSqrt2, negRecipSqrt6 ]
	];

	// 6 - front left-right up
	// [ FLU, FRU, FD, B ]
	a2b6 = [
		[ sqrt6Div4, sqrt6Div4, sqrt6Div4, sqrt6Div4 ],
		[ sqrt3Div6, sqrt3Div6, sqrt3Div6, negSqrt3Div2 ],
		[ recipSqrt2, negRecipSqrt2, 0, 0 ],
		[ recipSqrt6, recipSqrt6, negSqrt6Div3, 0 ]
	];
	b2a6 = [
		[ recipSqrt6, sqrt3Div6, recipSqrt2, recipSqrt6 ],
		[ recipSqrt6, sqrt3Div6, negRecipSqrt2, recipSqrt6 ],
		[ recipSqrt6, sqrt3Div6, 0, negSqrt6Div3 ],
		[ recipSqrt6, negSqrt3Div2, 0, 0 ]
	];

	// 7 - front left-right down
	// [ FLD, FRD, FU, B ]
	a2b7 = [
		[ sqrt6Div4, sqrt6Div4, sqrt6Div4, sqrt6Div4 ],
		[ sqrt3Div6, sqrt3Div6, sqrt3Div6, negSqrt3Div2 ],
		[ recipSqrt2, negRecipSqrt2, 0, 0 ],
		[ negRecipSqrt6, negRecipSqrt6, sqrt6Div3, 0 ]
	];
	b2a7 = [
		[ recipSqrt6, sqrt3Div6, recipSqrt2, negRecipSqrt6 ],
		[ recipSqrt6, sqrt3Div6, negRecipSqrt2, negRecipSqrt6 ],
		[ recipSqrt6, sqrt3Div6, 0, sqrt6Div3 ],
		[ recipSqrt6, negSqrt3Div2, 0, 0 ]
	];

	// display input btest
	btest.postln;

	// b2a transform
	atest = transform.value(btest, b2a0);
	atest.postln;

	// a2b transform
	btest = transform.value(atest, a2b0);

)



// Scale for unscaled ('uns')
(
	var a2b0, a2b1, a2b2, a2b3, a2b4, a2b5, a2b6, a2b7;
	var b2a0, b2a1, b2a2, b2a3, b2a4, b2a5, b2a6, b2a7;

	var w0;

	var transform, btest, atest;

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
	var sqrt2Div4 = 2.sqrt/4;
	
	// matrix transform
	transform = { arg testInMatrix, transformMatrix;
		var testOutMatrix;
		testOutMatrix = [];

		transformMatrix.do({ arg item, i;
			testOutMatrix = testOutMatrix.add((item * testInMatrix).sum)
		});
		
		testOutMatrix
	};


	// test matrix input
//	btest = [ 2.sqrt.reciprocal, 1, 0, 0 ];		// front
//	btest = [ 2.sqrt.reciprocal, 1.neg, 0, 0 ];		// back
	btest = [ 2.sqrt.reciprocal, 2.sqrt.reciprocal, 2.sqrt.reciprocal, 0 ]; // front left
//	btest = [ 2.sqrt.reciprocal, 2.sqrt.reciprocal, 2.sqrt.reciprocal.neg, 0 ]; // front right
//	btest = [ 2.sqrt.reciprocal, 2.sqrt.reciprocal.neg, 2.sqrt.reciprocal, 0 ]; // back left
//	btest = [ 2.sqrt.reciprocal, 2.sqrt.reciprocal.neg, 2.sqrt.reciprocal.neg, 0 ]; // back right
//	btest = [ 2.sqrt.reciprocal, 0, 1, 0 ];		// left
//	btest = [ 2.sqrt.reciprocal, 0, 0 , 1 ];		// up
//	btest = [ 2.sqrt.reciprocal, 2.sqrt.reciprocal.neg, 0, 2.sqrt.reciprocal.neg ]; // back down
//	btest = [ 1, 0, 0, 0 ];
//	btest = [ 0, 1, 0, 0 ];
//	btest = [ 0, 0, 1, 0 ];
//	btest = [ 0, 0, 0, 1 ];
//	btest = [ 1, 1, 1, 1 ];
//	btest = [ 1, 1, -1, 1 ];
	atest = [];


	//transform matricies

	// 0 - orthogonal (front left up)
	// [ FLU, FRD, BLD, BRU ]
	a2b0 = [
		[ sqrt2Div4, sqrt2Div4, sqrt2Div4, sqrt2Div4 ],
		[ 0.5, 0.5, -0.5, -0.5 ],
		[ 0.5, -0.5, 0.5, -0.5 ],
		[ 0.5, -0.5, -0.5, 0.5 ]
	];
	b2a0 = [
		[ recipSqrt2, 0.5, 0.5, 0.5 ],
		[ recipSqrt2, 0.5, -0.5, -0.5 ],
		[ recipSqrt2, -0.5, 0.5, -0.5 ],
		[ recipSqrt2, -0.5, -0.5, 0.5 ]
	];

	// 1 - front left down
	// [ FLD, FRU, BLU, BRD ]
	a2b1 = [
		[ sqrt2Div4, sqrt2Div4, sqrt2Div4, sqrt2Div4 ],
		[ 0.5, 0.5, -0.5, -0.5 ],
		[ 0.5, -0.5, 0.5, -0.5 ],
		[ -0.5, 0.5, 0.5, -0.5 ]
	];
	b2a1 = [
		[ recipSqrt2, 0.5, 0.5, -0.5 ],
		[ recipSqrt2, 0.5, -0.5, 0.5 ],
		[ recipSqrt2, -0.5, 0.5, 0.5 ],
		[ recipSqrt2, -0.5, -0.5, -0.5 ]
	];

	// 2 - front left-right
	// [ FL, FR, BU, BD ]
	a2b2 = [
		[ sqrt2Div4, sqrt2Div4, sqrt2Div4, sqrt2Div4 ],
		[ 0.5, 0.5, -0.5, -0.5 ],
		[ recipSqrt2, negRecipSqrt2, 0, 0 ],
		[ 0, 0, recipSqrt2, negRecipSqrt2 ]
	];
	b2a2 = [
		[ recipSqrt2, 0.5, recipSqrt2, 0 ],
		[ recipSqrt2, 0.5, negRecipSqrt2, 0 ],
		[ recipSqrt2, -0.5, 0, recipSqrt2 ],
		[ recipSqrt2, -0.5, 0, negRecipSqrt2 ]
	];

	// 3 - front up-down
	// [ FU, FD, BL, BR ]
	a2b3 = [
		[ sqrt2Div4, sqrt2Div4, sqrt2Div4, sqrt2Div4 ],
		[ 0.5, 0.5, -0.5, -0.5 ],
		[ 0, 0, recipSqrt2, negRecipSqrt2 ],
		[ recipSqrt2, negRecipSqrt2, 0, 0 ]
	];
	b2a3 = [
		[ recipSqrt2, 0.5, 0, recipSqrt2 ],
		[ recipSqrt2, 0.5, 0, negRecipSqrt2 ],
		[ recipSqrt2, -0.5, recipSqrt2, 0 ],
		[ recipSqrt2, -0.5, negRecipSqrt2, 0 ]
		
	];

	// 4 - front & back down
	// [ F, BD, BLU, BRU ]
	a2b4 = [
		[ sqrt2Div4, sqrt2Div4, sqrt2Div4, sqrt2Div4 ],
		[ sqrt3Div2, negSqrt3Div6, negSqrt3Div6, negSqrt3Div6 ],
		[ 0, 0, recipSqrt2, negRecipSqrt2 ],
		[ 0, negSqrt6Div3, recipSqrt6, recipSqrt6 ]
	];
	b2a4 = [
		[ recipSqrt2, sqrt3Div2, 0, 0 ],
		[ recipSqrt2, negSqrt3Div6, 0, negSqrt6Div3 ],
		[ recipSqrt2, negSqrt3Div6, recipSqrt2, recipSqrt6 ],
		[ recipSqrt2, negSqrt3Div6, negRecipSqrt2, recipSqrt6 ]
	];

	// 5 - front & back up
	// [ F, BU, BLD, BRD ]
	a2b5 = [
		[ sqrt2Div4, sqrt2Div4, sqrt2Div4, sqrt2Div4 ],
		[ sqrt3Div2, negSqrt3Div6, negSqrt3Div6, negSqrt3Div6 ],
		[ 0, 0, recipSqrt2, negRecipSqrt2 ],
		[ 0, sqrt6Div3, negRecipSqrt6, negRecipSqrt6 ]
	];
	b2a5 = [
		[ recipSqrt2, sqrt3Div2, 0, 0 ],
		[ recipSqrt2, negSqrt3Div6, 0, sqrt6Div3 ],
		[ recipSqrt2, negSqrt3Div6, recipSqrt2, negRecipSqrt6 ],
		[ recipSqrt2, negSqrt3Div6, negRecipSqrt2, negRecipSqrt6 ]
	];

	// 6 - front left-right up
	// [ FLU, FRU, FD, B ]
	a2b6 = [
		[ sqrt2Div4, sqrt2Div4, sqrt2Div4, sqrt2Div4 ],
		[ sqrt3Div6, sqrt3Div6, sqrt3Div6, negSqrt3Div2 ],
		[ recipSqrt2, negRecipSqrt2, 0, 0 ],
		[ recipSqrt6, recipSqrt6, negSqrt6Div3, 0 ]
	];
	b2a6 = [
		[ recipSqrt2, sqrt3Div6, recipSqrt2, recipSqrt6 ],
		[ recipSqrt2, sqrt3Div6, negRecipSqrt2, recipSqrt6 ],
		[ recipSqrt2, sqrt3Div6, 0, negSqrt6Div3 ],
		[ recipSqrt2, negSqrt3Div2, 0, 0 ]
	];

	// 7 - front left-right down
	// [ FLD, FRD, FU, B ]
	a2b7 = [
		[ sqrt2Div4, sqrt2Div4, sqrt2Div4, sqrt2Div4 ],
		[ sqrt3Div6, sqrt3Div6, sqrt3Div6, negSqrt3Div2 ],
		[ recipSqrt2, negRecipSqrt2, 0, 0 ],
		[ negRecipSqrt6, negRecipSqrt6, sqrt6Div3, 0 ]
	];
	b2a7 = [
		[ recipSqrt2, sqrt3Div6, recipSqrt2, negRecipSqrt6 ],
		[ recipSqrt2, sqrt3Div6, negRecipSqrt2, negRecipSqrt6 ],
		[ recipSqrt2, sqrt3Div6, 0, sqrt6Div3 ],
		[ recipSqrt2, negSqrt3Div2, 0, 0 ]
	];

	// display input btest
	btest.postln;

	// b2a transform
	atest = transform.value(btest, b2a0);
	atest.postln;

	// a2b transform
	btest = transform.value(atest, a2b0);
)

