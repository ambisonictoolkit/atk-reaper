// -------------------------------------------------------
// ATK (Sandbox)
//
// Transcoding Aaron Heller's Octave code
// into SC3 language--using the MatLib quark extension.
//
//
// Coded by Joseph Anderson 2011
//
// -------------------------------------------------------

//------------------------------------------------------------------------// (Gerzon's) Diametric Decoder Theorem (DDT)//------------------------------------------------------------------------//// Much of the code below is a transcoding of Aaron Heller's Octave// code available at: http://www.ai.sri.com/ajh/ambisonics///// Benjamin, et al., "Localization in Horizontal-Only Ambisonic Systems"// Preprint from AES-121, 10/2006, San Francisco//// Heller's original functions are noted through comments in each// functions help field.//// Transcoding to Python/Numpy for use in muse/ATK/SC3 by// Joseph Anderson <josephlloydanderson@mac.com>//// aes_paper.m (expanded version of speaker_matrix.m) contains the// following functions:////   velocity_gain_matrix()**            : compute alpha, beta, and gamma//   speaker_matrix()                    : compute alpha, beta, and gamma//   decoder_gain_matrix()               : compute decoder matrix//// ----------------------------------------
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
// ----------------------------------------// the following functions are not included// as they duplicate muse/ATK functionality////   az2dir()                            : convert azimuth to directon cosines//   degrees()                           : convert radians to degrees//   radians()                           : convert degrees to radians//   gain_to_db()////   rectangular_speaker_arrays()        : example decodes//   hexagonal_speaker_arrays()          : example decodes////
//// NOTE: speaker_matrix() and velocity_gain_matrix() are the same code.//       It appears these two separate names are used (in error) in Heller's//       code. The expanded version, aes_paper.m defines velocity_gain_matrix(),//       but calls speaker_matrix().////------------------------------------------------------------------------//////------------------------------------------------------------------------// DDT and related decoder matrix gains////   NOTE:   These are the functions that compute gains to generate//           loudspeaker feeds, and are not the functions which return//           decoded B-format. See decoders, below.//////   speaker_matrix                  Heller's DDT (helper function)//   decoder_gain_matrix             Heller's DDT (returns decoder gains)//   panto_reg_decoder_gain_matrix   pantophonic//   peri_reg_decoder_gain_matrix    periphonic//   quad_decoder_gain_matrix        quad////------------------------------------------------------------------------

// NOTE: the addition of a spherical coordinates class (like Polar.sc)
//		to convert between spherical and cartesian coordinates
//		will replace the following global functions, and improve
//		usability -------> TBD!!
(
	~spherToCart = { arg spher;		// [ r, theta, phi ]
			var r, theta, phi;
			
			r 		= spher.at(0);
			theta	= spher.at(1);
			phi 		= spher.at(2);
		
			[
				r * phi.cos * theta.cos,
				r * phi.cos * theta.sin,
				r * phi.sin
			]
		};

	~cartToSpher = { arg cart;		// [ x, y, z ]
			var x, y, z;
			
			x = cart.at(0);
			y = cart.at(1);
			z = cart.at(2);
			
			r = (x.squared + y.squared + z.squared).sqrt;
		
			[
				r,
				y.atan2(x),
				(z/r).asin
			]
		};

)



//   speaker_matrix                  Heller's DDT (helper function)
(
	~atkSpeakerMatrix = { arg positions, k;
		var m, n, s;
		var directions, pos, dir;

		// allow entry of positions as    		// transpose for convenience
    		// e.g., speaker positions are now in columns
    		// rather than rows
    		positions = positions.flop;

    		// n = number of speaker pairs    		// m = number of dimensions,    		//        2=horizontal, 3=periphonic 		m = positions.rows;
		n = positions.cols;
		
	    	// scatter matrix accumulator	    	s = Matrix.newClear(m, m);

		// speaker directions matrix	    	directions = Matrix.newClear(m, n);	
		n.do({ arg i;

	        	// get the i'th speaker position	        	// e.g., select the i'th column        		pos = positions.getCol(i);

        		// normalize to get direction cosines        		dir = pos /  pos.squared.sum.sqrt;
        		
        		// form scatter matrix and accumulate        		s = s + Matrix.with(dir * dir.flop);
        		// form matrix of speaker directions        		directions.putCol(i, dir)

			});
			
		// return resulting matrix
	 	sqrt(1/2) * n * k * ( s.inverse * directions);
	}
)


//   decoder_gain_matrix             Heller's DDT (returns decoder gains)
(
	~atkDecoderGainMatrix = { arg positions, k;
		var positions2;
		var sm, m, n;
		var decoderMatrix;

		// list all of the speakers		// i.e., expand to actual pairs		positions2 = positions ++ (-1 * positions);		
	    	// get velocity gains	    	// NOTE: this comment from Heller seems to be slightly	    	//       misleading, in that the gains returned will be	    	//       scaled by k, which may not request a velocity	    	//       gain. I.e., k = 1 isn't necessarily true, as it	    	//       is assigned as an argument to this function.
	    	sm = ~atkSpeakerMatrix.value(positions2, k);
	    
	    	// n = number of speakers
	    	// m = number of dimensions,		//        2=horizontal, 3=periphonic 		m = sm.rows;
		n = sm.cols;
		
		// build decoder matrix 		// rows are W, X, and Y gains		// NOTE: this matrix construction can be simplified		//       with a concatenation (hstack) of a column		//       of ones and sm	    	decoderMatrix = Matrix.newClear(m + 1, n) + 1;
	    	n.do({ arg i;
			m.do({ arg j;
				decoderMatrix.put(j + 1, i, sm.at(j, i))
				});
		    });

		// return resulting matrix
		// ALSO: the below code calls for the complex conjugate
		//       of decoder_matrix. As we are expecting real vaules,
		//       we may regard this call as redundant.
		// res = sqrt(2)/n * decoder_matrix.conj().transpose()
		2.sqrt/n * decoderMatrix.flop;

	}
)


//   panto_reg_decoder_gain_matrix   pantophonic
(
	~atkPantoRegDecoderGainMatrix = { arg numSpeakers, orientation, k;
		var theta;
		var g0, g1;
		var decoderMatrix;

		// define (local) function to return theta from speaker number		theta = { arg speaker;
			if ( orientation == 'point',
				{ ((2.0 * speaker)/numSpeakers) * pi },        // 'point' case
				{ ((1.0 + (2.0 * speaker))/numSpeakers) * pi } // 'flat' case, default
			)
		};

	    // define constants	    g0 = 1.0;	    g1 = 2.sqrt;

		// calculate decoding matrix	
		decoderMatrix = Matrix.newClear(numSpeakers, 3); // start w/ empty matrix
	
		numSpeakers.do({ arg i;
			decoderMatrix.putRow(i, [
				g0,
	              k * g1 * theta.value(i).cos,
	              k * g1 * theta.value(i).sin
			])
			});
		
		// return resulting matrix
		2.sqrt/numSpeakers * decoderMatrix
	}
)


//   peri_reg_decoder_gain_matrix    periphonic
(
	~atkPeriRegDecoderGainMatrix = { arg numSpeakerPairs, elevation, orientation, k;
		var theta, spher, positions;
		var top, bottom;
		var decoderMatrix;

		// generate speaker pair positions		// start with polar positions. . .
		theta = [];
		numSpeakerPairs.do({arg i;
			theta = theta ++ [2 * pi * i / numSpeakerPairs]}
		);
		if ( orientation == 'flat',
			{ theta = theta + (pi / numSpeakerPairs) });       // 'flat' case

		// convert to spherical coordinates [ [r], [theta], [phi] ]...
		spher = [
			Array.newClear(numSpeakerPairs).fill(1),
			theta,
			Array.newClear(numSpeakerPairs).fill(elevation)
		];
		
		// . . . then convert from spherical to cartesian [ [x, y, z] ]		positions = Matrix.with(~spherToCart.value(spher)).flop;

		// compute the decoder
		decoderMatrix = ~atkDecoderGainMatrix.value(positions, k);

		// reorder the bottom polygon
		top = decoderMatrix[..(numSpeakerPairs-1)];
		bottom = decoderMatrix[(numSpeakerPairs)..];

		bottom = if ( (orientation == 'flat') && (numSpeakerPairs.mod(2) == 1),
			{ bottom.rotate((numSpeakerPairs/2 + 1).asInteger) }, // odd, 'flat'
			{ bottom.rotate((numSpeakerPairs/2).asInteger) }     // 'flat' case, default
		);
		
		decoderMatrix = top ++ bottom;

		decoderMatrix
	}
)


//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
// test
(
	var positions, k;
	
	positions = Matrix.with([[1, 1, 1], [1, -1, 1], [1, 1, -1], [1, -1, -1]]);
//	positions = Matrix.with([[1, 1], [1, -1]]);
	k = 1;
//	k = 2.sqrt.reciprocal;

//	~atkSpeakerMatrix.value(positions, k);
//	~atkDecoderGainMatrix.value(positions, k);

//	~atkPantoRegDecoderGainMatrix.value(4, 'flat', k);
//	~atkPantoRegDecoderGainMatrix.value(4, 'point', k);
//	~atkPantoRegDecoderGainMatrix.value(5, 'flat', k);
//	~atkPantoRegDecoderGainMatrix.value(5, 'point', k);
//	~atkPantoRegDecoderGainMatrix.value(6, 'flat', k);
//	~atkPantoRegDecoderGainMatrix.value(6, 'point', k);

//	~atkPeriRegDecoderGainMatrix.value(3, pi/3, 'flat', k);
//	~atkPeriRegDecoderGainMatrix.value(3, pi/3, 'point', k);
//	~atkPeriRegDecoderGainMatrix.value(4, pi/3, 'flat', k);
//	~atkPeriRegDecoderGainMatrix.value(4, pi/3, 'point', k);
//	~atkPeriRegDecoderGainMatrix.value(5, pi/3, 'flat', k);
	~atkPeriRegDecoderGainMatrix.value(5, pi/3, 'point', k);

)