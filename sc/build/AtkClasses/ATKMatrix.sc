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
		
    		// n = number of speaker pairs
    		// m = number of dimensions,
    		//        2=horizontal, 3=periphonic 
		m = this.positions.cols;
		n = this.positions.rows;
	}
	
	dim { ^m }
	numSpeakers { ^n * 2 }
	
	matrix {
		var s, directions, pos, dir;
		
	    	// scatter matrix accumulator
	    	s = Matrix.newClear(m, m);

		// speaker directions matrix
        	// NOTE: this isn't the user supplied directions arg
	    	directions = Matrix.newClear(m, n);
	
		n.do({ arg i;

			// allow entry of positions as
	    		// transpose for convenience
	    		// e.g., speaker positions are now in columns
	    		// rather than rows, then
	        	// get the i'th speaker position
	        	// e.g., select the i'th column
        		pos = positions.flop.getCol(i);

        		// normalize to get direction cosines
        		dir = pos /  pos.squared.sum.sqrt;
        		
        		// form scatter matrix and accumulate
        		s = s + Matrix.with(dir * dir.flop);

        		// form matrix of speaker directions
        		directions.putCol(i, dir)

			});
			
		// return resulting matrix
	 	^sqrt(1/2) * n * k * ( s.inverse * directions);
	}
}


//        - orientation    : Orientation of the A-format channel tetrahedron//            flu = front left up:          FLU, FRD, BLD, BRU//            fld = front left down:        FLD, FRU, BLU, BRD//            flr = front left-right:       FL, FR, BU, BD//            fud = front up-down:          FU, FD, BL, BR//            fbd = front and back down:    F, BD, BLU, BRU //            fbu = front and back up:      F, BU, BLD, BRD//            flru = front left-right up:   FLU, FRU, FD, B//            flrd = front left-right down: FLD, FRD, FU, B////        - weight : The scaling on W in the output a-format signal://            can = canonical scaling of 1/sqrt(2)//            dec = scaling for decorrelated soundfields, weight of 1/sqrt(3) on W//            car = scaling for cardioid response, weight of sqrt(3) on W//            uns = unscaled, weight of 1 on W////            dec is the usual choice for use in reverberators

//   bToA_matrix
AtkBtoAMatrix {
	var <orientation, <weight;
	var decoderMatrix;

	*new { arg orientation = 'flu', weight = 'can';
		^super.newCopyArgs(orientation).initWeight(weight);
	}
	
	initWeight { arg weightArg;
		switch ( weightArg,
			'can', { weight = 1 },			// canonical
			'dec', { weight = (2/3).sqrt },	// decorrelated (on the sphere)
			'car', { weight = 6.sqrt },		// cardioid
			'uns', { weight = 2.sqrt }		// unscaled, W_gain = 1
		)
	}

	dim { ^3 }

	numSpeakers { ^4 }			// hardwired for FOA

	dirSpeakers {
		^this.matrix.removeCol(0).asArray.collect({arg item;
			item.asCartesian.asSpherical.angles
		})
	}

	matrix {
		var recSqrt2 = 2.sqrt.reciprocal;
		var sqrt3Div2 = 3.sqrt/2;
		var sqrt3Div6 = 3.sqrt/6;
		var sqrt6Div3 = 6.sqrt/3;
		var recSqrt6 = 6.sqrt.reciprocal;

		switch ( orientation,

			// 0 - orthogonal (front left up)
			// [ FLU, FRD, BLD, BRU ]
			'flu', { decoderMatrix = [
					[ 0.5, 0.5, 0.5, 0.5 ],					[ 0.5, 0.5, -0.5, -0.5 ],					[ 0.5, -0.5, 0.5, -0.5 ],					[ 0.5, -0.5, -0.5, 0.5 ]
				]
			},
			// 1 - front left down
			// [ FLD, FRU, BLU, BRD ]
			'fld', { decoderMatrix = [
					[ 0.5, 0.5, 0.5, -0.5 ],					[ 0.5, 0.5, -0.5, 0.5 ],					[ 0.5, -0.5, 0.5, 0.5 ],					[ 0.5, -0.5, -0.5, -0.5 ]				]
			},
			// 2 - front left-right
			// [ FL, FR, BU, BD ]
			'flr', { decoderMatrix = [
					[ 0.5, 0.5, recSqrt2, 0 ],					[ 0.5, 0.5, recSqrt2.neg, 0 ],					[ 0.5, -0.5, 0, recSqrt2 ],					[ 0.5, -0.5, 0, recSqrt2.neg ]				]
			},
			// 3 - front up-down
			// [ FU, FD, BL, BR ]
			'fud', { decoderMatrix = [
					[ 0.5, 0.5, 0, recSqrt2 ],					[ 0.5, 0.5, 0, recSqrt2.neg ],					[ 0.5, -0.5, recSqrt2, 0 ],					[ 0.5, -0.5, recSqrt2.neg, 0 ]				]
			},
			// 4 - front & back down
			// [ F, BD, BLU, BRU ]
			'fbd', { decoderMatrix = [
					[ 0.5, sqrt3Div2, 0, 0 ],					[ 0.5, sqrt3Div6.neg, 0, sqrt6Div3.neg ],					[ 0.5, sqrt3Div6.neg, recSqrt2, recSqrt6 ],					[ 0.5, sqrt3Div6.neg, recSqrt2.neg, recSqrt6 ]				]
			},
			// 5 - front & back up
			// [ F, BU, BLD, BRD ]
			'fbu', { decoderMatrix = [
					[ 0.5, sqrt3Div2, 0, 0 ],					[ 0.5, sqrt3Div6.neg, 0, sqrt6Div3 ],					[ 0.5, sqrt3Div6.neg, recSqrt2, recSqrt6.neg ],					[ 0.5, sqrt3Div6.neg, recSqrt2.neg, recSqrt6.neg ]				]
			},
			// 6 - front left-right up
			// [ FLU, FRU, FD, B ]
			'flru', { decoderMatrix = [
					[ 0.5, sqrt3Div6, recSqrt2, recSqrt6 ],					[ 0.5, sqrt3Div6, recSqrt2.neg, recSqrt6 ],					[ 0.5, sqrt3Div6, 0, sqrt6Div3.neg ],					[ 0.5, sqrt3Div2.neg, 0, 0 ]				]
			},
			// 7 - front left-right down
			// [ FLD, FRD, FU, B ]
			'flrd', { decoderMatrix = [
					[ 0.5, sqrt3Div6, recSqrt2, recSqrt6.neg ],					[ 0.5, sqrt3Div6, recSqrt2.neg, recSqrt6.neg ],					[ 0.5, sqrt3Div6, 0, sqrt6Div3 ],					[ 0.5, sqrt3Div2.neg, 0, 0 ]				]
			}
		);
		decoderMatrix = Matrix.with(decoderMatrix);
		^decoderMatrix.putCol(0, weight * decoderMatrix.getCol(0))
	}
}


//   speaker_kernel                  (helper function)
//	*put together a kernel helper function for kernel decoders?
//AtkSpeakerKernel {



//   decoder_gain_matrix             Heller's DDT (returns decoder gains)
//   panto_reg_decoder_gain_matrix   pantophonic
//   peri_reg_decoder_gain_matrix    periphonic
//   quad_decoder_gain_matrix        quad

//   b_to_ITU5					Wiggins coefficients (Single Band Decoders)
//   b_to_7           				Wiggins coefficients (Single Band Decoders)

//	b_to_stereo					virtual stereo microphone decoding
//	b_to_mono						virtual mono microphone decoding

//	NOTE: b_to_a conversion can go in here, too, and there is an
//			attraction to do so--in terms of consistency
//			and that AtkDecoderMatrix reports back angles!

AtkDecoderMatrix {
	var <kind;
	var <matrix;
	var <dirSpeakers;
	var <>shelfFreq, <shelfK;


	*newDiametric { arg directions = [ pi/4, 3*pi/4 ], k = 'single';
		^super.newCopyArgs('diametric').initDiametric(directions, k);
	}
	
	*newPanto { arg numSpeakers = 4, orientation = 'flat', k = 'single';
		^super.newCopyArgs('panto').initPanto(numSpeakers, orientation, k);
	}
	
	*newPeri { arg numSpeakerPairs = 4, elevation = 0.61547970867039,
				orientation = 'flat', k = 'single';
		^super.newCopyArgs('peri').initPeri(numSpeakerPairs, elevation,
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
		var speakerMatrix, m, n;

		switch (directions.rank,					// 2D or 3D?
			1, {									// 2D

				// find positions
				positions = Matrix.with(
					directions.collect({ arg item;
						Polar.new(1, item).asPoint.asArray
					})
				);

				// list all of the speakers
				// i.e., expand to actual pairs
				positions2 = positions ++ (positions.neg);
		

			    // set speaker directions for instance
				dirSpeakers = positions2.asArray.collect({ arg item;
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

				// list all of the speakers
				// i.e., expand to actual pairs
				positions2 = positions ++ (positions.neg);
		

			    // set speaker directions for instance
				dirSpeakers = positions2.asArray.collect({ arg item;
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
	    
	    	// n = number of speakers
	    	// m = number of dimensions,
		//        2=horizontal, 3=periphonic 
		m = speakerMatrix.rows;
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
	
	initPanto { arg numSpeakers, orientation, k;

		var g0, g1, theta;

	    	g0 = 1.0;								// decoder gains
	    	g1 = 2.sqrt;							// 0, 1st order


		// return theta from speaker number
		theta = numSpeakers.collect({ arg speaker;
			switch (orientation,
				'flat',	{ ((1.0 + (2.0 * speaker))/numSpeakers) * pi },
				'point',	{ ((2.0 * speaker)/numSpeakers) * pi }
			)
		});
		theta = (theta + pi).mod(2pi) - pi;

	    // set speaker directions for instance
		dirSpeakers = theta;

		// initialise k
		k = this.initK2D(k);


		// build decoder matrix 
		matrix = Matrix.newClear(numSpeakers, 3); // start w/ empty matrix
	
		numSpeakers.do({ arg i;
			matrix.putRow(i, [
				g0,
	              k * g1 * theta.at(i).cos,
	              k * g1 * theta.at(i).sin
			])
			});
		matrix = 2.sqrt/numSpeakers * matrix
	}
	
	initPeri { arg numSpeakerPairs, elevation, orientation, k;

		var theta, directions, upDirs, downDirs, upMatrix, downMatrix;

		// generate speaker pair positions
		// start with polar positions. . .
		theta = [];
		numSpeakerPairs.do({arg i;
			theta = theta ++ [2 * pi * i / numSpeakerPairs]}
		);
		if ( orientation == 'flat',
			{ theta = theta + (pi / numSpeakerPairs) });       // 'flat' case

		// collect directions [ [theta, phi], ... ]
		// upper ring only
		directions = [
			theta,
			Array.newClear(numSpeakerPairs).fill(elevation)
		].flop;


	    // prepare speaker directions for instance
		upDirs = (directions + pi).mod(2pi) - pi;

		downDirs = upDirs.collect({ arg angles;
			Spherical.new(1, angles.at(0), angles.at(1)).neg.angles
		});
		
		// initialise k
		k = this.initK2D(k);


		// build decoder matrix
		matrix = AtkDecoderMatrix.newDiametric(directions, k).matrix;

		// reorder the lower polygon
		upMatrix = matrix[..(numSpeakerPairs-1)];
		downMatrix = matrix[(numSpeakerPairs)..];

		if ( (orientation == 'flat') && (numSpeakerPairs.mod(2) == 1),
			{									 // odd, 'flat'

				downDirs = downDirs.rotate((numSpeakerPairs/2 + 1).asInteger);
				downMatrix = downMatrix.rotate((numSpeakerPairs/2 + 1).asInteger)

			}, {     								// 'flat' case, default

				downDirs = downDirs.rotate((numSpeakerPairs/2).asInteger);
				downMatrix = downMatrix.rotate((numSpeakerPairs/2).asInteger)
			}
		);
		
		dirSpeakers = upDirs ++ downDirs;		// set speaker directions
		matrix = upMatrix ++ downMatrix;			// set matrix

	}
	
	initQuad { arg angle, k;

		var alpha, beta;

	    // set speaker directions for instance
	    dirSpeakers = [ angle, pi - angle, (pi - angle).neg, angle.neg ];


		// initialise k
		k = this.initK2D(k);

		// calculate alpha, beta (scaled by k)
		alpha	= k / (2.sqrt * angle.cos);
		beta		= k / (2.sqrt * angle. sin);

		// build decoder matrix 
	    matrix = 2.sqrt/4 * Matrix.with([
	    		[ 1, alpha, 		beta 	],
	        	[ 1, alpha.neg, 	beta 	],
	        	[ 1, alpha.neg, 	beta.neg	],
	        	[ 1, alpha, 		beta.neg	]
	    ])
	}
	
	initStereo { arg angle, pattern;

		var alpha, beta, gamma;
	    
	    // set speaker directions for instance
	    dirSpeakers = [ pi/6, pi.neg/6 ];

		// calculate alpha, beta, gamma (scaled by pattern)
		alpha	= pattern * angle.cos;
		beta		= pattern * angle.sin;
		gamma	= (1.0 - pattern) * 2.sqrt;

		// build decoder matrix, and set for instance
	    matrix = Matrix.with([
	    		[ gamma, alpha, beta		],
	        	[ gamma, alpha, beta.neg	]
	    ])
	}
	
	initMono { arg theta, phi, pattern;

	    // set speaker directions for instance
	    dirSpeakers = [ 0 ];

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

	    // set speaker directions for instance
	    dirSpeakers = [ 0, pi/6, 110/180 * pi, 110/180 * pi.neg, pi.neg/6 ];

		// build decoder matrix
		// Wigging's Matricies (credit contribution/copyright at top)
		switch (irregKind,
			'focused', {
			    matrix = Matrix.with([
			    		[ 0.2000,  0.1600,  0.0000 ],
			        	[ 0.4250,  0.3600,  0.4050 ],
			        	[ 0.4700, -0.3300,  0.4150 ],
			        	[ 0.4700, -0.3300, -0.4150 ],
			        	[ 0.4250,  0.3600, -0.4050 ]
			    ])
			},
			'equal', {
			    matrix = Matrix.with([
			    		[ 0.0000,  0.0850,  0.0000 ],
			        	[ 0.3650,  0.4350,  0.3400 ],
			        	[ 0.5550, -0.2850,  0.4050 ],
			        	[ 0.5550, -0.2850, -0.4050 ],
			        	[ 0.3650,  0.4350, -0.3400 ]
			    ])
			}
		)
	}

	dim { ^matrix.cols - 1}

	numSpeakers { ^matrix.rows }

}


//	b_to_uhj            			"Ambisonic Decoders for HDTV" (1992)
//	b_to_binaural       			HRTF decoding

// kernel matrix for kernel decoders
//AtkDecoderKernel {
