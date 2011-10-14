// -------------------------------------------------------
// ATK Matrix
//
// Coded by Joseph Anderson 2011
//
// -------------------------------------------------------

// NOTE: Really need to clean up these comments!!!

//------------------------------------------------------------------------
// (Gerzon's) Diametric Decoder Theorem (DDT)
//------------------------------------------------------------------------
//
// Much of the code below is a transcoding of Aaron Heller's Octave
// code available at: http://www.ai.sri.com/ajh/ambisonics/
//
// Benjamin, et al., "Localization in Horizontal-Only Ambisonic Systems"
// Preprint from AES-121, 10/2006, San Francisco
//
// Heller's original functions are noted through comments in each
// functions help field.
//
// Transcoding to Python/Numpy for use in muse/ATK/SC3 by
// Joseph Anderson <josephlloydanderson@mac.com>
//
// aes_paper.m (expanded version of speaker_matrix.m) contains the
// following functions:
//
//   velocity_gain_matrix()**            : compute alpha, beta, and gamma
//   speaker_matrix()                    : compute alpha, beta, and gamma
//   decoder_gain_matrix()               : compute decoder matrix
//
// ----------------------------------------
// the following functions are not included in SC3
// as are included in muse/ATK and not immediately
// useful in SC3 implementations
//
//   rV()                                : compute the Makita direction and rV
//   rE()                                : compute rE (and direction?)
//
//   _virtual_mic()                      : virtual mic angle and directivity
//   decoder_matrix_to_virtual_mic()     : computes loudspeaker 'virtual mics' 
//
// ----------------------------------------
// the following functions are not included
// as they duplicate muse/ATK functionality
//
//   az2dir()                            : convert azimuth to directon cosines
//   degrees()                           : convert radians to degrees
//   radians()                           : convert degrees to radians
//   gain_to_db()
//
//   rectangular_speaker_arrays()        : example decodes
//   hexagonal_speaker_arrays()          : example decodes
//
//
//
// NOTE: speaker_matrix() and velocity_gain_matrix() are the same code.
//       It appears these two separate names are used (in error) in Heller's
//       code. The expanded version, aes_paper.m defines velocity_gain_matrix(),
//       but calls speaker_matrix().
//
//------------------------------------------------------------------------
//
//
//------------------------------------------------------------------------
// DDT and related decoder matrix gains
//
//   NOTE:   These are the functions that compute gains to generate
//           loudspeaker feeds, and are not the functions which return
//           decoded B-format. See decoders, below.
//
//
//   speaker_matrix                  Heller's DDT (helper function)
//   decoder_gain_matrix             Heller's DDT (returns decoder gains)
//   panto_reg_decoder_gain_matrix   pantophonic
//   peri_reg_decoder_gain_matrix    periphonic
//   quad_decoder_gain_matrix        quad
//
//------------------------------------------------------------------------

//FOA {
//	var <w, <x, <y, <z;
//	
//	*new {arg w, x, y, z;
//		^super.newCopyArgs(w, x, y, z);
//	}
//	
//	*ar {arg w, x, y, z;
//		^this.new(w, x, y, z);
//	}
//	
//	madd {arg mul = 1, add = 0;
//		^MulAdd.ar([w, x, y, z], mul, add);	
//	}
//	
//	sig {arg mul = 1, add = 0;
//		^[w, x, y, z] * mul + add;
//	}
//	
//	asUGenInput {
//		^[w, x, y, z];
//	}
//	
//	asAudioRateInput {
//		^[w, x, y, z];
//	}
//}


//-----------------------------------------------------------------------
// matrix decoders

//   speaker_matrix                  Heller's DDT (helper function)
AtkSpeakerMatrix {
	var <positions, <k, m, n;

	*new { arg directions, k;
		var positions;
		
		switch (directions.rank,					// 2D or 3D?
			1, { positions = Matrix.with(			// 2D
					directions.collect({ arg item;
						Polar.new(1, item).asPoint.asArray
					})
				)
			},
			2, { positions = Matrix.with(			// 3D
					directions.collect({ arg item;
						Spherical.new(1, item.at(0), item.at(1)).asCartesian.asArray
					})
				)
			}
		);
		
		^super.newCopyArgs(positions, k).initDiametric;
	}

	*newPositions { arg positions, k;
		^super.newCopyArgs(positions, k).initDiametric;
	}

	initDiametric {
		
    		// n = number of output channel (speaker) pairs
    		// m = number of dimensions,
    		//        2=horizontal, 3=periphonic 
		m = this.positions.cols;
		n = this.positions.rows;
	}
	
	dim { ^m }
	numChans { ^n * 2 }
	
	matrix {
		var s, directions, pos, dir;
		
	    	// scatter matrix accumulator
	    	s = Matrix.newClear(m, m);

		// output channel (speaker) directions matrix
        	// NOTE: this isn't the user supplied directions arg
	    	directions = Matrix.newClear(m, n);
	
		n.do({ arg i;

			// allow entry of positions as
	    		// transpose for convenience
	    		// e.g., output channel (speaker) positions are now in columns
	    		// rather than rows, then
	        	// get the i'th output channel (speaker) position
	        	// e.g., select the i'th column
        		pos = positions.flop.getCol(i);

        		// normalize to get direction cosines
        		dir = pos /  pos.squared.sum.sqrt;
        		
        		// form scatter matrix and accumulate
        		s = s + Matrix.with(dir * dir.flop);

        		// form matrix of output channel (speaker) directions
        		directions.putCol(i, dir)

			});
			
		// return resulting matrix
	 	^sqrt(1/2) * n * k * ( s.inverse * directions);
	}
}


FOADecoderMatrix {
	var <kind;
	var <matrix;
	var <dirChans;
	var <>shelfFreq, <shelfK;


	*newDiametric { arg directions = [ pi/4, 3*pi/4 ], k = 'single';
		^super.newCopyArgs('diametric').initDiametric(directions, k);
	}
	
	*newPanto { arg numChans = 4, orientation = 'flat', k = 'single';
		^super.newCopyArgs('panto').initPanto(numChans, orientation, k);
	}
	
	*newPeri { arg numChanPairs = 4, elevation = 0.61547970867039,
				orientation = 'flat', k = 'single';
		^super.newCopyArgs('peri').initPeri(numChanPairs, elevation,
			orientation, k);
	}
	
	*newQuad { arg angle = pi/4, k = 'single';
		^super.newCopyArgs('quad').initQuad(angle, k);
	}
	
	*newStereo { arg angle = pi/2, pattern = 0.5;
		^super.newCopyArgs('stereo').initStereo(angle, pattern);
	}

	*newMono { arg theta = 0, phi = 0, pattern = 0;
		^super.newCopyArgs('mono').initMono(theta, phi, pattern);
	}

	*new5_0 { arg irregKind = 'focused';
		^super.newCopyArgs('5.0').init5_0(irregKind);
	}

	*newBtoA { arg orientation = 'flu', weight = 'dec';
		^super.newCopyArgs('BtoA').initBtoA(orientation, weight);
	}

	initK2D { arg k;

		if ( k.isNumber, {
				^k
			}, {
				switch ( k,
					'velocity', 	{ ^1 },
					'energy', 	{ ^2.reciprocal.sqrt },
					'controlled', { ^2.reciprocal },
					'single', 	{ ^2.reciprocal.sqrt },
					'dual', 		{
						shelfFreq = 400.0;
						shelfK = [(3/2).sqrt, 3.sqrt/2];
						^1;
					}
				)
			}
		)
	}

	initK3D { arg k;

		if ( k.isNumber, {
				^k
			}, {
				switch ( k,
					'velocity', 	{ ^1 },
					'energy', 	{ ^3.reciprocal.sqrt },
					'controlled', { ^3.reciprocal },
					'single', 	{ ^3.reciprocal.sqrt },
					'dual', 		{
						shelfFreq = 400.0;
						shelfK = [2.sqrt, (2/3).sqrt];
						^1;
					}
				)
			}
		)
	}

	initDiametric { arg directions, k;
		
		var positions, positions2;
		var speakerMatrix, n;

		switch (directions.rank,					// 2D or 3D?
			1, {									// 2D

				// find positions
				positions = Matrix.with(
					directions.collect({ arg item;
						Polar.new(1, item).asPoint.asArray
					})
				);

				// list all of the output channels (speakers)
				// i.e., expand to actual pairs
				positions2 = positions ++ (positions.neg);
		

			    // set output channel (speaker) directions for instance
				dirChans = positions2.asArray.collect({ arg item;
					item.asPoint.asPolar.angle
				});
			
				// initialise k
				k = this.initK2D(k);
			},
			2, {									// 3D

				// find positions
				positions = Matrix.with(
					directions.collect({ arg item;
						Spherical.new(1, item.at(0), item.at(1)).asCartesian.asArray
					})
				);

				// list all of the output channels (speakers)
				// i.e., expand to actual pairs
				positions2 = positions ++ (positions.neg);
		

			    // set output channel (speaker) directions for instance
				dirChans = positions2.asArray.collect({ arg item;
					item.asCartesian.asSpherical.angles
				});
				
				// initialise k
				k = this.initK3D(k);
			}
		);


	    	// get velocity gains
	    	// NOTE: this comment from Heller seems to be slightly
	    	//       misleading, in that the gains returned will be
	    	//       scaled by k, which may not request a velocity
	    	//       gain. I.e., k = 1 isn't necessarily true, as it
	    	//       is assigned as an argument to this function.
	    	speakerMatrix = AtkSpeakerMatrix.newPositions(positions2, k).matrix;
	    
	    	// n = number of output channels (speakers)
		n = speakerMatrix.cols;

		// build decoder matrix 
		// resulting rows (after flop) are W, X, Y, Z gains
		matrix = speakerMatrix.insertRow(0, Array.fill(n, {1}));

		// return resulting matrix
		// ALSO: the below code calls for the complex conjugate
		//       of decoder_matrix. As we are expecting real vaules,
		//       we may regard this call as redundant.
		// res = sqrt(2)/n * decoder_matrix.conj().transpose()
		matrix = 2.sqrt/n * matrix.flop;
	}
	
	initPanto { arg numChans, orientation, k;

		var g0, g1, theta;

	    	g0 = 1.0;								// decoder gains
	    	g1 = 2.sqrt;							// 0, 1st order


		// return theta from output channel (speaker) number
		theta = numChans.collect({ arg channel;
			switch (orientation,
				'flat',	{ ((1.0 + (2.0 * channel))/numChans) * pi },
				'point',	{ ((2.0 * channel)/numChans) * pi }
			)
		});
		theta = (theta + pi).mod(2pi) - pi;

	    // set output channel (speaker) directions for instance
		dirChans = theta;

		// initialise k
		k = this.initK2D(k);


		// build decoder matrix
		matrix = Matrix.newClear(numChans, 3); // start w/ empty matrix
	
		numChans.do({ arg i;
			matrix.putRow(i, [
				g0,
	              k * g1 * theta.at(i).cos,
	              k * g1 * theta.at(i).sin
			])
			});
		matrix = 2.sqrt/numChans * matrix
	}
	
	initPeri { arg numChanPairs, elevation, orientation, k;

		var theta, directions, upDirs, downDirs, upMatrix, downMatrix;

		// generate output channel (speaker) pair positions
		// start with polar positions. . .
		theta = [];
		numChanPairs.do({arg i;
			theta = theta ++ [2 * pi * i / numChanPairs]}
		);
		if ( orientation == 'flat',
			{ theta = theta + (pi / numChanPairs) });       // 'flat' case

		// collect directions [ [theta, phi], ... ]
		// upper ring only
		directions = [
			theta,
			Array.newClear(numChanPairs).fill(elevation)
		].flop;


	    // prepare output channel (speaker) directions for instance
		upDirs = (directions + pi).mod(2pi) - pi;

		downDirs = upDirs.collect({ arg angles;
			Spherical.new(1, angles.at(0), angles.at(1)).neg.angles
		});
		
		// initialise k
		k = this.initK2D(k);


		// build decoder matrix
		matrix = FOADecoderMatrix.newDiametric(directions, k).matrix;

		// reorder the lower polygon
		upMatrix = matrix[..(numChanPairs-1)];
		downMatrix = matrix[(numChanPairs)..];

		if ( (orientation == 'flat') && (numChanPairs.mod(2) == 1),
			{									 // odd, 'flat'

				downDirs = downDirs.rotate((numChanPairs/2 + 1).asInteger);
				downMatrix = downMatrix.rotate((numChanPairs/2 + 1).asInteger)

			}, {     								// 'flat' case, default

				downDirs = downDirs.rotate((numChanPairs/2).asInteger);
				downMatrix = downMatrix.rotate((numChanPairs/2).asInteger)
			}
		);
		
		dirChans = upDirs ++ downDirs;		// set output channel (speaker) directions
		matrix = upMatrix ++ downMatrix;			// set matrix

	}
	
	initQuad { arg angle, k;

		var g0, g1, g2;

	    // set output channel (speaker) directions for instance
	    dirChans = [ angle, pi - angle, (pi - angle).neg, angle.neg ];


		// initialise k
		k = this.initK2D(k);

		// calculate g1, g2 (scaled by k)
		g0	= 1;
		g1	= k / (2.sqrt * angle.cos);
		g2	= k / (2.sqrt * angle.sin);

		// build decoder matrix 
	    matrix = 2.sqrt/4 * Matrix.with([
	    		[ g0, g1, 	g2 		],
	        	[ g0, g1.neg, g2 		],
	        	[ g0, g1.neg, g2.neg	],
	        	[ g0, g1, 	g2.neg	]
	    ])
	}
	
	initStereo { arg angle, pattern;

		var g0, g1, g2;
	    
	    // set output channel (speaker) directions for instance
	    dirChans = [ pi/6, pi.neg/6 ];

		// calculate g0, g1, g2 (scaled by pattern)
		g0	= (1.0 - pattern) * 2.sqrt;
		g1	= pattern * angle.cos;
		g2	= pattern * angle.sin;

		// build decoder matrix, and set for instance
	    matrix = Matrix.with([
	    		[ g0, g1, g2		],
	        	[ g0, g1, g2.neg	]
	    ])
	}
	
	initMono { arg theta, phi, pattern;

	    // set output channel (speaker) directions for instance
	    dirChans = [ 0 ];

		// build decoder matrix, and set for instance
	    matrix = Matrix.with([
	    		[
	    			(1.0 - pattern) * 2.sqrt,
	    			pattern * theta.cos * phi.cos,
	    			pattern * theta.sin * phi.cos,
	    			pattern * phi.sin
	    		]
	    ])
	}

	init5_0 { arg irregKind;

	    // set output channel (speaker) directions for instance
	    dirChans = [ 0, pi/6, 110/180 * pi, 110/180 * pi.neg, pi.neg/6 ];

		// build decoder matrix
		// Wigging's Matricies (credit contribution/copyright at top)
		matrix = switch (irregKind,
			'focused', {[
		    		[ 0.2000,  0.1600,  0.0000 ],
		        	[ 0.4250,  0.3600,  0.4050 ],
		        	[ 0.4700, -0.3300,  0.4150 ],
		        	[ 0.4700, -0.3300, -0.4150 ],
		        	[ 0.4250,  0.3600, -0.4050 ]
			]},
			'equal', {[
		    		[ 0.0000,  0.0850,  0.0000 ],
		        	[ 0.3650,  0.4350,  0.3400 ],
		        	[ 0.5550, -0.2850,  0.4050 ],
		        	[ 0.5550, -0.2850, -0.4050 ],
		        	[ 0.3650,  0.4350, -0.3400 ]
		    ]},
			'four', {[
		    		[ 0.0000,  0.0000,  0.0000 ],
		        	[ 0.4250,  0.3850,  0.3300 ],
		        	[ 0.6300, -0.2750,  0.2850 ],
		        	[ 0.6300, -0.2750, -0.2850 ],
		        	[ 0.4250,  0.3850, -0.3300 ]
		    ]}
		);
		matrix = Matrix.with(matrix)
	}

	initBtoA { arg orientation, weight;

		var recSqrt2 = 2.sqrt.reciprocal;
		var sqrt3Div2 = 3.sqrt/2;
		var sqrt3Div6 = 3.sqrt/6;
		var sqrt6Div3 = 6.sqrt/3;
		var recSqrt6 = 6.sqrt.reciprocal;
		var g0;
		
		
		// build decoder matrix, and set for instance
		g0 = switch ( weight,
			'dec', { (2/3).sqrt },	// decorrelated (on the sphere)
			'can', { 1 },			// canonical
			'uns', { 2.sqrt },		// unscaled, W_gain = 1
			'car', { 6.sqrt }		// cardioid
		);

		matrix = switch ( orientation,

			// 0 - orthogonal (front left up)
			// [ FLU, FRD, BLD, BRU ]
			'flu', {[
				[ 0.5, 0.5, 0.5, 0.5 ],
				[ 0.5, 0.5, -0.5, -0.5 ],
				[ 0.5, -0.5, 0.5, -0.5 ],
				[ 0.5, -0.5, -0.5, 0.5 ]
			]},
			// 1 - front left down
			// [ FLD, FRU, BLU, BRD ]
			'fld', {[
				[ 0.5, 0.5, 0.5, -0.5 ],
				[ 0.5, 0.5, -0.5, 0.5 ],
				[ 0.5, -0.5, 0.5, 0.5 ],
				[ 0.5, -0.5, -0.5, -0.5 ]
			]},
			// 2 - front left-right
			// [ FL, FR, BU, BD ]
			'flr', {[
				[ 0.5, 0.5, recSqrt2, 0 ],
				[ 0.5, 0.5, recSqrt2.neg, 0 ],
				[ 0.5, -0.5, 0, recSqrt2 ],
				[ 0.5, -0.5, 0, recSqrt2.neg ]
			]},
			// 3 - front up-down
			// [ FU, FD, BL, BR ]
			'fud', {[
				[ 0.5, 0.5, 0, recSqrt2 ],
				[ 0.5, 0.5, 0, recSqrt2.neg ],
				[ 0.5, -0.5, recSqrt2, 0 ],
				[ 0.5, -0.5, recSqrt2.neg, 0 ]
			]},
			// 4 - front & back down
			// [ F, BD, BLU, BRU ]
			'fbd', {[
				[ 0.5, sqrt3Div2, 0, 0 ],
				[ 0.5, sqrt3Div6.neg, 0, sqrt6Div3.neg ],
				[ 0.5, sqrt3Div6.neg, recSqrt2, recSqrt6 ],
				[ 0.5, sqrt3Div6.neg, recSqrt2.neg, recSqrt6 ]
			]},
			// 5 - front & back up
			// [ F, BU, BLD, BRD ]
			'fbu', {[
				[ 0.5, sqrt3Div2, 0, 0 ],
				[ 0.5, sqrt3Div6.neg, 0, sqrt6Div3 ],
				[ 0.5, sqrt3Div6.neg, recSqrt2, recSqrt6.neg ],
				[ 0.5, sqrt3Div6.neg, recSqrt2.neg, recSqrt6.neg ]
			]},
			// 6 - front left-right up
			// [ FLU, FRU, FD, B ]
			'flru', {[
				[ 0.5, sqrt3Div6, recSqrt2, recSqrt6 ],
				[ 0.5, sqrt3Div6, recSqrt2.neg, recSqrt6 ],
				[ 0.5, sqrt3Div6, 0, sqrt6Div3.neg ],
				[ 0.5, sqrt3Div2.neg, 0, 0 ]
			]},
			// 7 - front left-right down
			// [ FLD, FRD, FU, B ]
			'flrd', {[
				[ 0.5, sqrt3Div6, recSqrt2, recSqrt6.neg ],
				[ 0.5, sqrt3Div6, recSqrt2.neg, recSqrt6.neg ],
				[ 0.5, sqrt3Div6, 0, sqrt6Div3 ],
				[ 0.5, sqrt3Div2.neg, 0, 0 ]
			]}
		);
		matrix = Matrix.with(matrix);
		matrix = matrix.putCol(0, g0 * matrix.getCol(0));
		
			    
	    // set output channel (speaker) directions for instance
	    dirChans = matrix.removeCol(0).asArray.collect({arg item;
			item.asCartesian.asSpherical.angles
		})
	}
	
	dim { ^matrix.cols - 1}

	numChans { ^matrix.rows }

	printOn { arg stream;
		stream << this.class.name << "(" <<* [kind, this.dim, this.numChans] <<")";
	}

}


//-----------------------------------------------------------------------
// martrix encoders

FOAEncoderMatrix {
	var <kind;
	var <matrix;
	var <dirChans;


	*newB {
		^super.newCopyArgs('b').initB;
	}

	*newAtoB { arg orientation = 'flu', weight = 'dec';
		^super.newCopyArgs('AtoB').initAtoB(orientation, weight);
	}

	*newOmni {
		^super.newCopyArgs('omni').initOmni;
	}

	*newDirection { arg theta = 0, phi = 0;
		^super.newCopyArgs('dir').initDirection(theta, phi);
	}

	*newStereo { arg angle = 0;
		^super.newCopyArgs('stereo').initStereo(angle);
	}

	*newQuad {
		^super.newCopyArgs('quad').initQuad;
	}

	*new5_0 {
		^super.newCopyArgs('5.0').init5_0;
	}

	*new7_0 {
		^super.newCopyArgs('7.0').init7_0;
	}

	*newDirections { arg directions, pattern = nil;
		^super.newCopyArgs('dirs').initDirections(directions, pattern);
	}

	*newPanto { arg numChans = 4, orientation = 'flat';
		^super.newCopyArgs('panto').initPanto(numChans, orientation);
	}

	*newPeri { arg numChanPairs = 4, elevation = 0.61547970867039,
				orientation = 'flat';
		^super.newCopyArgs('peri').initPeri(numChanPairs, elevation,
			orientation);
	}
	
	*newZoomH2 { arg angles = [pi/3, 3/4*pi], pattern = 0.5857, k = 1;
		^super.newCopyArgs('zoomH2').initZoomH2(angles, pattern, k);
	}

	init2D {

		var g0 = 2.sqrt.reciprocal;
	    
		// build encoder matrix, and set for instance
		matrix = Matrix.newClear(3, dirChans.size); // start w/ empty matrix
	
		dirChans.do({ arg theta, i;
			matrix.putCol(i, [
				g0,
	              theta.cos,
	              theta.sin
			])
		})
	}

	init3D {

		var g0 = 2.sqrt.reciprocal;
	    
		// build encoder matrix, and set for instance
		matrix = Matrix.newClear(4, dirChans.size); // start w/ empty matrix
	
		dirChans.do({ arg thetaPhi, i;
			matrix.putCol(i, [
				g0,
	              thetaPhi.at(1).cos * thetaPhi.at(0).cos,
	              thetaPhi.at(1).cos * thetaPhi.at(0).sin,
	              thetaPhi.at(1).sin
			])
		})
	}

	initInv2D { arg pattern;

		var g0 = 2.sqrt.reciprocal;
	    
		// build 'decoder' matrix, and set for instance
		matrix = Matrix.newClear(dirChans.size, 3); 	// start w/ empty matrix

		if ( pattern.isArray,
			{
				dirChans.do({ arg theta, i;			// mic positions, indivd patterns
					matrix.putRow(i, [
						(1.0 - pattern.at(i)),
			              pattern.at(i) * theta.cos,
			              pattern.at(i) * theta.sin
					])
				})
			}, {
				dirChans.do({ arg theta, i;			// mic positions
					matrix.putRow(i, [
						(1.0 - pattern),
			              pattern * theta.cos,
			              pattern * theta.sin
					])
				})
			}
		);

		// invert to encoder matrix
		matrix = matrix.pseudoInverse;

		matrix = matrix.putRow(0, matrix.getRow(0) * g0); // scale W
	}

	initInv3D { arg pattern;

		var g0 = 2.sqrt.reciprocal;
	    
		// build 'decoder' matrix, and set for instance
		matrix = Matrix.newClear(dirChans.size, 4); 	// start w/ empty matrix

		if ( pattern.isArray,
			{
				dirChans.do({ arg thetaPhi, i;		// mic positions, indivd patterns
					matrix.putRow(i, [
						(1.0 - pattern.at(i)),
			              pattern.at(i) * thetaPhi.at(1).cos * thetaPhi.at(0).cos,
			              pattern.at(i) * thetaPhi.at(1).cos * thetaPhi.at(0).sin,
			              pattern.at(i) * thetaPhi.at(1).sin
					])
				})
			}, {
				dirChans.do({ arg thetaPhi, i;		// mic positions
					matrix.putRow(i, [
						(1.0 - pattern),
			              pattern * thetaPhi.at(1).cos * thetaPhi.at(0).cos,
			              pattern * thetaPhi.at(1).cos * thetaPhi.at(0).sin,
			              pattern * thetaPhi.at(1).sin
					])
				})
			}
		);

		// invert to encoder matrix
		matrix = matrix.pseudoInverse;

		matrix = matrix.putRow(0, matrix.getRow(0) * g0); // scale W
	}

	initB {

	    // set input channel directions for instance
	    dirChans = [ inf ];

		// build encoder matrix, and set for instance
	    matrix = Matrix.newIdentity(4)
	}

	initAtoB { arg orientation, weight;
		var bToAMatrix;

		// retrieve corresponding A-format decoder
		bToAMatrix = FOADecoderMatrix.newBtoA(orientation, weight);

	    // set input channel directions for instance
	    dirChans = bToAMatrix.dirChans;

		// build encoder matrix, and set for instance
	    matrix = bToAMatrix.matrix.inverse
	}

	initOmni {

	    // set input channel directions for instance
	    dirChans = [ inf ];

		// build encoder matrix, and set for instance
	    matrix = Matrix.with([
		    	[ 2.sqrt.reciprocal ]
		])
	}

	initDirection { arg theta, phi;

	    // set input channel directions for instance
	    (phi == 0).if (
		    {
				dirChans = [ theta ];
    				this.init2D
			}, {
	    			dirChans = [ [theta, phi] ];
    				this.init3D
			}
		)
	}

	initStereo { arg angle;

	    // set input channel directions for instance
	    dirChans = [ pi/2 - angle, (pi/2 - angle).neg ];

	    this.init2D
	}

	initQuad {
		
	    // set input channel directions for instance
	    dirChans = [ pi/4, pi * 3/4, pi.neg * 3/4, pi.neg/4 ];

	    this.init2D
	}
	
	init5_0 {

	    // set input channel directions for instance
	    dirChans = [ 0, pi/6, 110/180 * pi, 110/180 * pi.neg, pi.neg/6 ];
	    
	    this.init2D
	}

	init7_0 {

	    // set input channel directions for instance
	    dirChans = [ 0, pi/6, pi/2, 135/180 * pi, 135/180 * pi.neg, pi.neg/2, pi.neg/6 ];
	    
	    this.init2D
	}

	initDirections { arg directions, pattern;

	    // set input channel directions for instance
	    dirChans = directions;

		switch (directions.rank,					// 2D or 3D?
			1, {									// 2D
				if ( pattern == nil, {
					this.init2D					// plane wave
				}, {
					this.initInv2D(pattern)			// mic inversion
				})
			},
			2, {									// 3D
				if ( pattern == nil, {
					this.init3D					// plane wave
				}, {
					this.initInv3D(pattern)			// mic inversion
				})
			}
		)
	}

	initPanto { arg numChans, orientation;

		var theta;

		// return theta from output channel (speaker) number
		theta = numChans.collect({ arg channel;
			switch (orientation,
				'flat',	{ ((1.0 + (2.0 * channel))/numChans) * pi },
				'point',	{ ((2.0 * channel)/numChans) * pi }
			)
		});
		theta = (theta + pi).mod(2pi) - pi;

	    // set input channel directions for instance
		dirChans = theta;

		this.init2D
	}
	
	initPeri { arg numChanPairs, elevation, orientation;

		var theta, directions, upDirs, downDirs, upMatrix, downMatrix;

		// generate input channel pair positions
		// start with polar positions. . .
		theta = [];
		numChanPairs.do({arg i;
			theta = theta ++ [2 * pi * i / numChanPairs]}
		);
		if ( orientation == 'flat',
			{ theta = theta + (pi / numChanPairs) });       // 'flat' case

		// collect directions [ [theta, phi], ... ]
		// upper ring only
		directions = [
			theta,
			Array.newClear(numChanPairs).fill(elevation)
		].flop;


	    // prepare output channel (speaker) directions for instance
		upDirs = (directions + pi).mod(2pi) - pi;

		downDirs = upDirs.collect({ arg angles;
			Spherical.new(1, angles.at(0), angles.at(1)).neg.angles
		});
		
		// reorder the lower polygon
		if ( (orientation == 'flat') && (numChanPairs.mod(2) == 1),
			{									 // odd, 'flat'
				downDirs = downDirs.rotate((numChanPairs/2 + 1).asInteger);
			}, {     								// 'flat' case, default
				downDirs = downDirs.rotate((numChanPairs/2).asInteger);
			}
		);
		
	    // set input channel directions for instance
		dirChans = upDirs ++ downDirs;

		this.init3D
	}

	initZoomH2 { arg angles, pattern, k;

	    // set input channel directions for instance
	    dirChans = [ angles.at(0), angles.at(0).neg, angles.at(1), angles.at(1).neg ];

		this.initInv2D(pattern);

		matrix = matrix.putRow(2, matrix.getRow(2) * k); // scale Y
	}
	
	dim { ^matrix.rows - 1}	

	numChans { ^matrix.cols }

	printOn { arg stream;
		stream << this.class.name << "(" <<* [kind, this.dim, this.numChans] <<")";
	}
}


////-----------------------------------------------------------------------
//// martrix transforms
//
//AtkTransMatrix {
//	var <kind;
//	var <matrix;
//
//
//	*newMirrorX {
//		^super.newCopyArgs('mirrorX').initMirrorX;
//	}
//
//	*newMirrorY {
//		^super.newCopyArgs('mirrorY').initMirrorY;
//	}
//
//	*newMirrorZ {
//		^super.newCopyArgs('mirrorZ').initMirrorZ;
//	}
//
//	*newMirrorO {
//		^super.newCopyArgs('mirrorO').initMirrorO;
//	}
//
//	*newRotate { arg angle = 0;
//		^super.newCopyArgs('rotate').initRotate(angle);
//	}
//
//	*newTilt { arg angle = 0;
//		^super.newCopyArgs('tilt').initTilt(angle);
//	}
//
//	*newTumble { arg angle = 0;
//		^super.newCopyArgs('tumble').initTumble(angle);
//	}
//
//	*newDirect { arg angle = 0;
//		^super.newCopyArgs('direct').initDirect(angle);
//	}
//
//	*newDirectX { arg angle = 0;
//		^super.newCopyArgs('directX').initDirectX(angle);
//	}
//
//	*newDirectY { arg angle = 0;
//		^super.newCopyArgs('directY').initDirectY(angle);
//	}
//
//	*newDirectZ { arg angle = 0;
//		^super.newCopyArgs('directZ').initDirectZ(angle);
//	}
//
//	initMirrorChan { arg chan;
//		matrix = matrix.put(chan, chan, matrix.get(chan, chan).neg)
//	}
//
//	initMirrorX {
//		var chan;
//		
//		// ambisonic channel index
//		chan = 1;
//
//		// build identity matrix 
//		matrix = Matrix.newIdentity(4);
//
//		// mirror it 
//		this.initMirrorChan(chan)
//	}
//
//	initMirrorY {
//		var chan;
//		
//		// ambisonic channel index
//		chan = 2;
//
//		// build identity matrix 
//		matrix = Matrix.newIdentity(4);
//
//		// mirror it 
//		this.initMirrorChan(chan)
//	}
//
//	initMirrorZ {
//		var chan;
//		
//		// ambisonic channel index
//		chan = 3;
//
//		// build identity matrix 
//		matrix = Matrix.newIdentity(4);
//
//		// mirror it 
//		this.initMirrorChan(chan)
//	}
//
//	initMirrorO {
//		var chans;
//		
//		// ambisonic channel index
//		chans = [1, 2, 3];
//
//		// build identity matrix 
//		matrix = Matrix.newIdentity(4);
//
//		// mirror it
//		chans.do({arg chan;
//			this.initMirrorChan(chan)
//		})
//	}
//
//	initRotate { arg angle;
//		var cosAngle, sinAngle;
//
//		// build transform matrix, and set for instance
//		// calculate cos, sin
//		cosAngle	= angle.cos;
//		sinAngle	= angle.sin;
//
//	    matrix = Matrix.with([
//	    		[ 1, 0, 			0,			0 ],
//	    		[ 0, cosAngle,	sinAngle.neg,	0 ],
//	    		[ 0, sinAngle, 	cosAngle,		0 ],
//	    		[ 0, 0, 			0,			1 ]
//	    ])
//	}
//	
//	initTilt { arg angle;
//		var cosAngle, sinAngle;
//
//		// build transform matrix, and set for instance
//		// calculate cos, sin
//		cosAngle	= angle.cos;
//		sinAngle	= angle.sin;
//
//	    matrix = Matrix.with([
//	    		[ 1, 0, 0,		0 			],
//	    		[ 0, 1, 0,		0 			],
//	    		[ 0,	0, cosAngle,	sinAngle.neg 	],
//	    		[ 0,	0, sinAngle, 	cosAngle 		]
//	    ])
//	}
//
//	initTumble { arg angle;
//		var cosAngle, sinAngle;
//
//		// build transform matrix, and set for instance
//		// calculate cos, sin
//		cosAngle	= angle.cos;
//		sinAngle	= angle.sin;
//
//	    matrix = Matrix.with([
//	    		[ 1, 0, 			0,	0 			],
//	    		[ 0, cosAngle,	0,	sinAngle.neg	],
//	    		[ 0, 0,			1, 	0 			],
//	    		[ 0, sinAngle,	0, 	cosAngle 		]
//	    ])
//	}
//
//	initDirect { arg angle;
//		var g0, g1;
//
//		// build transform matrix, and set for instance
//		g0 = (1 + angle.sin).sqrt;
//		g1 = (1 - angle.sin).sqrt;
//
//	    matrix = Matrix.with([
//	    		[ g0,	0,	0,	0 	],
//	    		[ 0, 	g1,	0,	0	],
//	    		[ 0, 	0,	g1, 	0 	],
//	    		[ 0, 	0,	0, 	g1 	]
//	    ])
//	}
//
//	initDirectX { arg angle;
//		var g0, g1;
//
//		// build transform matrix, and set for instance
//		g0 = (1 + angle.sin).sqrt;
//		g1 = (1 - angle.sin).sqrt;
//
//	    matrix = Matrix.with([
//	    		[ g0,	0,	0,	0 	],
//	    		[ 0, 	g1,	0,	0	],
//	    		[ 0, 	0,	g0, 	0 	],
//	    		[ 0, 	0,	0, 	g0 	]
//	    ])
//	}
//
//	initDirectY { arg angle;
//		var g0, g1;
//
//		// build transform matrix, and set for instance
//		g0 = (1 + angle.sin).sqrt;
//		g1 = (1 - angle.sin).sqrt;
//
//	    matrix = Matrix.with([
//	    		[ g0,	0,	0,	0 	],
//	    		[ 0, 	g0,	0,	0	],
//	    		[ 0, 	0,	g1, 	0 	],
//	    		[ 0, 	0,	0, 	g0 	]
//	    ])
//	}
//
//	initDirectZ { arg angle;
//		var g0, g1;
//
//		// build transform matrix, and set for instance
//		g0 = (1 + angle.sin).sqrt;
//		g1 = (1 - angle.sin).sqrt;
//
//	    matrix = Matrix.with([
//	    		[ g0,	0,	0,	0 	],
//	    		[ 0, 	g0,	0,	0	],
//	    		[ 0, 	0,	g0, 	0 	],
//	    		[ 0, 	0,	0, 	g1 	]
//	    ])
//	}
//
//	dim { ^matrix.rows - 1}	
//
//	numChans { ^matrix.cols }
//
//	printOn { arg stream;
//		stream << this.class.name << "(" <<* [kind, this.dim, this.numChans] <<")";
//	}
//}
//

//------------------------------------------------------------------------
// kernel decoders

FOADecoderKernel {
	var <kind, <subjectID;
	var <kernel;
	var <dirChans;
	

	*newSpherical { arg subjectID = 0004, kernelSize = 512, server = Server.default;
		^super.newCopyArgs('spherical', subjectID).initKernel(kernelSize, server);
	}
	
	*newListen { arg subjectID = 1002, server = Server.default;
		^super.newCopyArgs('listen', subjectID).initKernel(512, server);
	}
	
	*newCIPIC { arg subjectID = 0021, server = Server.default;
		^super.newCopyArgs('cipic', subjectID).initKernel(256, server);
	}
	
	*newUHJ { arg kernelSize = 512, server = Server.default;
		^super.newCopyArgs('uhj', 0).initKernel(kernelSize, server);
	}
	
	initPath {
		
		var kernelLibPath;
		var decodersPath;
		
		kernelLibPath = PathName.new("/Library/Application Support/ATK/kernels");

		if ( kernelLibPath.isFolder.not, {	// is kernel lib installed for all users?
			kernelLibPath = PathName.new("~") +/+ kernelLibPath // no? set for single user
		});

		decodersPath	= PathName.new("/FOA/decoders");

		^kernelLibPath +/+ decodersPath +/+ PathName.new(kind.asString)
	}

	initKernel { arg kernelSize, server;
		
		var databasePath, subjectPath;
		var chans;
		var sampleRate;
		var errorMsg;
		
		// constants
		chans = 2;			// stereo kernel
		
		// init dirChans (output channel (speaker) directions) and kernel sr
		if ( kind == 'uhj', {
		    dirChans = [ pi/6, pi.neg/6 ];
			sampleRate = "None";
		}, {
			dirChans = [ 5/9 * pi, 5/9 * pi.neg ];
			sampleRate = server.sampleRate.asString;
		});
		

		// init kernel root, generate subjectPath and kernelFiles
		databasePath = this.initPath;

		subjectPath = databasePath +/+ PathName.new(
			sampleRate ++ "/" ++ 
			kernelSize ++ "/" ++
			subjectID.asString.padLeft(4, "0")
		);
		
		
		// attempt to load kernel
		if ( server.serverRunning.not, {		// is server running?
			
			// throw server error!
			Error(
				"Please boot server: %. Decoder kernel failed to load.".format(
					server.name.asString
				)
			).throw
		}, {
			if ( subjectPath.isFolder.not, {	// does kernel path exist?

				case
				// --> missing kernel database
					{ databasePath.isFolder.not }
					{
						errorMsg = "ATK kernel database missing!" +
							"Please install % database.".format(kind)
					}

				// --> unsupported SR
					{ PathName.new(subjectPath.parentLevelPath(2)).isFolder.not }
					{
						"Supported samplerates:".warn;
						PathName.new(subjectPath.parentLevelPath(3)).folders.do({
							arg folder;
							("\t" + folder.folderName).postln;
					});

						errorMsg = "Samplerate = % is not available for".format(sampleRate)
							+
							"% kernel decoder.".format(kind)
					}

				// --> unsupported kernelSize
					{ PathName.new(subjectPath.parentLevelPath(1)).isFolder.not }
					{
						"Supported kernel sizes:".warn;
						PathName.new(subjectPath.parentLevelPath(2)).folders.do({
							arg folder;
							("\t" + folder.folderName).postln;
					});

						errorMsg = "Kernel size = % is not available for".format(kernelSize)
						+
						"% kernel decoder.".format(kind)
					}

				// --> unsupported subject
					{ subjectPath.isFolder.not }
					{
						"Supported subjects:".warn;
						PathName.new(subjectPath.parentLevelPath(1)).folders.do({
							arg folder;
							("\t" + folder.folderName).postln;
					});

						errorMsg = "Subject % is not available for".format(subjectID)
						+
						"% kernel decoder.".format(kind)
					};

				// throw error!
				"\n".post;
				Error(errorMsg).throw
			}, {
				// Else... everything is fine! Load kernel.
				kernel = subjectPath.files.collect({ arg kernelPath;
					chans.collect({ arg chan;
						Buffer.readChannel(server, kernelPath.fullPath, channels: [chan],
							action: { arg buf;
								(
									"Kernel %, channel % loaded.".format(
										kernelPath.fileName, chan
									)
								).postln
							}
						)
					})
				})
			})
		})
	}

	free {
		kernel.shape.at(0).do({ arg i;
			kernel.shape.at(1).do({ arg j;
				kernel.at(i).at(j).free;
				(
					"Kernel %, channel % freed.".format(
						PathName.new(kernel.at(i).at(j).path).fileName, j
					)
				).postln
			})
		})
	}

	dim { ^kernel.shape.at(0) - 1}

	numChans { ^kernel.shape.at(1) }

	kernelSize { ^kernel.at(0).at(0).numFrames }

	printOn { arg stream;
		stream << this.class.name << "(" <<*
			[kind, this.dim, this.numChans, subjectID, this.kernelSize] <<")";
	}
}


//------------------------------------------------------------------------
// kernel encoders

FOAEncoderKernel {
	var <kind, <subjectID;
	var <kernel;
	var <dirChans;
	

	*newUHJ { arg kernelSize = 512, server = Server.default;
		^super.newCopyArgs('uhj', 0).initKernel(kernelSize, server);
	}

	*newSuper { arg kernelSize = 512, server = Server.default;
		^super.newCopyArgs('super', 0).initKernel(kernelSize, server);
	}

	initPath {
		
		var kernelLibPath;
		var encodersPath;
		
		kernelLibPath = PathName.new("/Library/Application Support/ATK/kernels");

		if ( kernelLibPath.isFolder.not, {	// is kernel lib installed for all users?
			kernelLibPath = PathName.new("~") +/+ kernelLibPath // no? set for single user
		});

		encodersPath	= PathName.new("/FOA/encoders");

		^kernelLibPath +/+ encodersPath +/+ PathName.new(kind.asString)
	}

	initKernel { arg kernelSize, server;
		
		var databasePath, subjectPath;
		var chans;
		var sampleRate;
		var errorMsg;
		
		// constants
		chans = 3;			// horizontal only kernel (at the moment), [w, x, y]
		
		// init dirChans (output channel (speaker) directions) and kernel sr
		switch ( kind,
			'super', {
				dirChans = [ pi/4, pi.neg/4 ];	 // approx, doesn't include phasiness
				sampleRate = "None"
			},
			'uhj', {
				dirChans = [ inf, inf ];
				sampleRate = server.sampleRate.asString;
			}
		);


		// init kernel root, generate subjectPath and kernelFiles
		databasePath = this.initPath;

		subjectPath = databasePath +/+ PathName.new(
			sampleRate ++ "/" ++ 
			kernelSize ++ "/" ++
			subjectID.asString.padLeft(4, "0")
		);

		// attempt to load kernel
		if ( server.serverRunning.not, {		// is server running?
			
			// throw server error!
			Error(
				"Please boot server: %. Encoder kernel failed to load.".format(
					server.name.asString
				)
			).throw
		}, {
			if ( subjectPath.isFolder.not, {	// does kernel path exist?

				case
				// --> missing kernel database
					{ databasePath.isFolder.not }
					{
						errorMsg = "ATK kernel database missing!" +
							"Please install % database.".format(kind)
					}

				// --> unsupported SR
					{ PathName.new(subjectPath.parentLevelPath(2)).isFolder.not }
					{
						"Supported samplerates:".warn;
						PathName.new(subjectPath.parentLevelPath(3)).folders.do({
							arg folder;
							("\t" + folder.folderName).postln;
					});

						errorMsg = "Samplerate = % is not available for".format(sampleRate)
							+
							"% kernel encoder.".format(kind)
					}

				// --> unsupported kernelSize
					{ PathName.new(subjectPath.parentLevelPath(1)).isFolder.not }
					{
						"Supported kernel sizes:".warn;
						PathName.new(subjectPath.parentLevelPath(2)).folders.do({
							arg folder;
							("\t" + folder.folderName).postln;
					});

						errorMsg = "Kernel size = % is not available for".format(kernelSize)
						+
						"% kernel encoder.".format(kind)
					}

				// --> unsupported subject
					{ subjectPath.isFolder.not }
					{
						"Supported subjects:".warn;
						PathName.new(subjectPath.parentLevelPath(1)).folders.do({
							arg folder;
							("\t" + folder.folderName).postln;
					});

						errorMsg = "Subject % is not available for".format(subjectID)
						+
						"% kernel encoder.".format(kind)
					};

				// throw error!
				"\n".post;
				Error(errorMsg).throw
			}, {
				// Else... everything is fine! Load kernel.
				kernel = subjectPath.files.collect({ arg kernelPath;
					chans.collect({ arg chan;
						Buffer.readChannel(server, kernelPath.fullPath, channels: [chan],
							action: { arg buf;
								(
									"Kernel %, channel % loaded.".format(
										kernelPath.fileName, chan
									)
								).postln
							}
						)
					})
				})
			})
		})
	}

	free {
		kernel.shape.at(0).do({ arg i;
			kernel.shape.at(1).do({ arg j;
				kernel.at(i).at(j).free;
				(
					"Kernel %, channel % freed.".format(
						PathName.new(kernel.at(i).at(j).path).fileName, j
					)
				).postln
			})
		})
	}

	dim { ^kernel.shape.at(1) - 1}

	numChans { ^kernel.shape.at(0) }

	kernelSize { ^kernel.at(0).at(0).numFrames }

	printOn { arg stream;
		stream << this.class.name << "(" <<*
			[kind, this.dim, this.numChans, subjectID, this.kernelSize] <<")";
	}
}


//------------------------------------------------------------------------
// Extension to PathName... here for the time being
+ PathName {

	parentLevelPath { arg index;
		
		var ci = this.colonIndices;

		^if( index == 0, {
			fullPath
		}, {		
			if((fullPath.last.isPathSeparator) && (ci.size > 1), {
				fullPath.copyRange(0, ci[ci.size - (1 + index)])
			}, {
				fullPath.copyRange(0, ci[ci.size - index])
			})
		})
	}
}