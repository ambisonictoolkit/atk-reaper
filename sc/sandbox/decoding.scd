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

//#------------------------------------------------------------------------//# (Gerzon's) Diametric Decoder Theorem (DDT)//#------------------------------------------------------------------------//#//# Much of the code below is a transcoding of Aaron Heller's Octave//# code available at: http://www.ai.sri.com/ajh/ambisonics///#//# Benjamin, et al., "Localization in Horizontal-Only Ambisonic Systems"//# Preprint from AES-121, 10/2006, San Francisco//#//# Heller's original functions are noted through comments in each//# functions help field.//#//# Transcoding to Python/Numpy for use in muse/ATK/SC3 by//# Joseph Anderson <josephlloydanderson@mac.com>//#//# aes_paper.m (expanded version of speaker_matrix.m) contains the//# following functions://#//#   velocity_gain_matrix()**            : compute alpha, beta, and gamma//#   speaker_matrix()                    : compute alpha, beta, and gamma//#   decoder_gain_matrix()               : compute decoder matrix//#//#   rV()                                : compute the Makita direction and rV//#   rE()                                : compute rE (and direction?)//#//#   _virtual_mic()                      : virtual mic angle and directivity//#   decoder_matrix_to_virtual_mic()     : computes loudspeaker 'virtual mics' //#//# ----------------------------------------//# the following functions are not included//# as they duplicate muse/ATK functionality//#//#   az2dir()                            : convert azimuth to directon cosines//#   degrees()                           : convert radians to degrees//#   radians()                           : convert degrees to radians//#   gain_to_db()//#//#   rectangular_speaker_arrays()        : example decodes//#   hexagonal_speaker_arrays()          : example decodes//#//#//# NOTE: speaker_matrix() and velocity_gain_matrix() are the same code.//#       It appears these two separate names are used (in error) in Heller's//#       code. The expanded version, aes_paper.m defines velocity_gain_matrix(),//#       but calls speaker_matrix().//#//#------------------------------------------------------------------------//////#------------------------------------------------------------------------//# DDT and related decoder matrix gains//#//#   NOTE:   These are the functions that compute gains to generate//#           loudspeaker feeds, and are not the functions which return//#           decoded B-format. See decoders, below.//#//#//#   speaker_matrix                  Heller's DDT (helper function)//#   decoder_gain_matrix             Heller's DDT (returns decoder gains)//#   panto_reg_decoder_gain_matrix   pantophonic//#   peri_reg_decoder_gain_matrix    periphonic//#   quad_decoder_gain_matrix        quad//#//#------------------------------------------------------------------------

// speaker_matrix
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


// decoder_gain_matrix
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


//--------------------------------------------------------------------------------
// test
(
	var positions, k;
	
	positions = Matrix.with([[1, 1, 1], [1, -1, 1], [1, 1, -1], [1, -1, -1]]);
//	positions = Matrix.with([[1, 1], [1, -1]]);
	k = 1;
//	k = 2.sqrt.reciprocal;

//	~atkSpeakerMatrix.value(positions, k);
	~atkDecoderGainMatrix.value(positions, k);

)