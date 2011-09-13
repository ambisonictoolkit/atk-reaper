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

// NOTE: The Spherical and Cartesian classes (sclang extensions) are required
//		for peri_reg_decoder_gain_matrix(). Additionally, these classes are
//		necessary from a usability perspective to specify speaker positions
//		(with correct radii) for decoder_gain_matrix(), the DDT decoder


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

		// list all of the speakers		// i.e., expand to actual pairs		positions2 = positions ++ (positions.neg);		
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

		// generate speaker pair positions
		// start with polar positions. . .
		theta = [];
		numSpeakerPairs.do({arg i;
			theta = theta ++ [2 * pi * i / numSpeakerPairs]}
		);
		if ( orientation == 'flat',
			{ theta = theta + (pi / numSpeakerPairs) });       // 'flat' case

		// convert to spherical coordinates [ [r, theta, phi], ... ]
		spher = [
			Array.newClear(numSpeakerPairs).fill(1),
			theta,
			Array.newClear(numSpeakerPairs).fill(elevation)
		].flop;
		
		// . . . then convert from spherical to cartesian [ [x, y, z], ... ]
		positions = Matrix.with(
			spher.collect({ arg item;
				Spherical.new(item.at(0), item.at(1), item.at(2)).asCartesian.asArray
			})
		);

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


//   quad_decoder_gain_matrix        quad
(
	~atkQuadDecoderGainMatrix = { arg angle, k;
		var alpha, beta;
		var decoderMatrix;

		// calculate alpha, beta (scaled by k)
		alpha   = k / (2.sqrt * angle.cos);
		beta    = k / (2.sqrt * angle. sin);


		// fill decoding matrix
	    decoderMatrix = Matrix.with([
	    		[1, alpha, beta],
	        	[1, alpha.neg, beta],
	        	[1, alpha.neg, beta.neg],
	        	[1, alpha, beta.neg]
	    ]);
	    
	    2.sqrt/4 * decoderMatrix
	}
)


//------------------------------------------------------------------------// Muse/ATK decoders built using the above decoder gain matrices//   Developed and coded by J Anderson//////   ++++ Single Band Regular Decoders////   panto_sbr//   peri_sbr////   ++++ Single Band Diametric Decoders////   decode_sbd//   quad_sbd//////   ++++ Dual Band Regular Decoders (Shelf-filtered)////   panto_dbr//   peri_dbr////   ++++ Dual Band Diametric Decoders (Shelf-filtered)////   decode_dbd//   quad_dbd        (See special decoders)////------------------------------------------------------------------------

//------------------------------------------------------------------------// Single Band Regular Decoders////   panto_sbr           pantophonic//   peri_sbr            periphonic//------------------------------------------------------------------------


//   panto_sbr           pantophonic
(
	~atkPantoSBR = { arg bSig, numSpeakers = 4, orientation = 'flat', k = 0.7071;
		var decoder;

		//return decoder gain matrix		decoder = ~atkPantoRegDecoderGainMatrix.value(numSpeakers, orientation, k);

		// decode here!		decoder.asArray.collect({arg channel;
			(channel * bSig[..2]).sum			// panto, need 1st 3 bSig channels
			})


	}
)


//   peri_sbr            periphonic
(
	~atkPeriSBR = { arg bSig, numSpeakerPairs = 4, elevation = 0.6155, orientation = 'flat',
						k = 0.5774;
		var decoder;

		//return decoder gain matrix
		decoder = ~atkPeriRegDecoderGainMatrix.value(numSpeakerPairs, elevation, orientation, k);

		// decode here!
		decoder.asArray.collect({arg channel;
			(channel * bSig).sum
			})


	}
)


//------------------------------------------------------------------------// Single Band Diametric Decoder////   decode_sbd//   quad_sbd//------------------------------------------------------------------------


//   decode_sbd
(
	~atkDecodeSBD = { arg bSig, positions, k = 0.7071;
		var decoder;

		//return decoder gain matrix
		decoder = ~atkDecoderGainMatrix.value(positions, k);

		// pantophonic or periphonic?		//    if panto, discard z		if ( positions.rows == 2, { bSig = bSig[..2]} );


		// decode here!
		decoder.asArray.collect({arg channel;
			(channel * bSig).sum
			})


	}
)


//   quad_sbd
(
	~atkQuadSBD = { arg bSig, angle = 0.7854, k = 0.7071;
		var decoder;

		//return decoder gain matrix
		decoder = ~atkQuadDecoderGainMatrix.value(angle, k);

		// decode here!
		decoder.asArray.collect({arg channel;
			(channel * bSig[..2]).sum			// panto, need 1st 3 bSig channels
			})


	}
)


//------------------------------------------------------------------------// Dual Band Regular Decoders////   panto_dbr           pantophonic//   peri_dbr            periphonic//------------------------------------------------------------------------


//   panto_dbr           pantophonic
(
	~atkPantoDBR = { arg bSig, numSpeakers = 4, orientation = 'flat', freq = 400.0;
		var decoder;
		var k, psychoK;
		
		// constants
		k = 1.0;							// k for gain matrix (velocity)
		psychoK = [(3/2).sqrt, 3.sqrt/2];	// panto k (energy at HF)

		
		// return decoder gain matrix
		decoder = ~atkPantoRegDecoderGainMatrix.value(numSpeakers, orientation, k);
		
		// shelf filter
		bSig = ~atkPsychoShelf.value(bSig, freq, psychoK);

		// decode here!
		decoder.asArray.collect({arg channel;
			(channel * bSig[..2]).sum			// panto, need 1st 3 bSig channels
			})


	}
)


//   peri_dbr            periphonic
(
	~atkPeriDBR = { arg bSig, numSpeakerPairs = 4, elevation = 0.6155, orientation = 'flat',
						freq = 400.0;
		var decoder;
		var k, psychoK;

		// constants
		k = 1.0;							// k for gain matrix (velocity)
		psychoK = [2.sqrt, (2/3).sqrt];		// peri k (energy at HF)

		
		//return decoder gain matrix
		decoder = ~atkPeriRegDecoderGainMatrix.value(numSpeakerPairs, elevation, orientation, k);

		// shelf filter
		bSig = ~atkPsychoShelf.value(bSig, freq, psychoK);

		// decode here!
		decoder.asArray.collect({arg channel;
			(channel * bSig).sum
			})


	}
)


//------------------------------------------------------------------------// Dual Band Diametric Decoder////   decode_dbd//   quad_dbd//------------------------------------------------------------------------

//   decode_dbd
(
	~atkDecodeDBD = { arg bSig, positions, freq = 400.0;
		var decoder;
		var k, psychoK;

		// constants
		k = 1.0;									// k for gain matrix (velocity)
		if ( positions.rows == 2,
			{ psychoK = [(3/2).sqrt, 3.sqrt/2] },		// panto k (energy at HF)
			{ psychoK = [2.sqrt, (2/3).sqrt] }		// peri k (energy at HF)
		);


		//return decoder gain matrix
		decoder = ~atkDecoderGainMatrix.value(positions, k);

		// shelf filter
		bSig = ~atkPsychoShelf.value(bSig, freq, psychoK);

		// pantophonic or periphonic?
		//    if panto, discard z
		if ( positions.rows == 2, { bSig = bSig[..2]} );


		// decode here!
		decoder.asArray.collect({arg channel;
			(channel * bSig).sum
			})


	}
)


//   quad_dbd
(
	~atkQuadDBD = { arg bSig, angle = 0.7854, freq = 400.0;
		var decoder;
		var k, psychoK;
		
		// constants
		k = 1.0;							// k for gain matrix (velocity)
		psychoK = [(3/2).sqrt, 3.sqrt/2];	// panto k (energy at HF)


		//return decoder gain matrix
		decoder = ~atkQuadDecoderGainMatrix.value(angle, k);

		// shelf filter
		bSig = ~atkPsychoShelf.value(bSig, freq, psychoK);

		// decode here!
		decoder.asArray.collect({arg channel;
			(channel * bSig[..2]).sum			// panto, need 1st 3 bSig channels
			})


	}
)


//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
// test (gain matricies)
(
	var positions, k;
	
//	positions = Matrix.with([[1, 1, 1], [1, -1, 1], [1, 1, -1], [1, -1, -1]]);
	positions = Matrix.with([[1, 1], [1, -1]]);
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

	~atkPeriRegDecoderGainMatrix.value(3, pi/3, 'flat', k);
//	~atkPeriRegDecoderGainMatrix.value(3, pi/3, 'point', k);
//	~atkPeriRegDecoderGainMatrix.value(4, pi/3, 'flat', k);
//	~atkPeriRegDecoderGainMatrix.value(4, pi/3, 'point', k);
//	~atkPeriRegDecoderGainMatrix.value(5, pi/3, 'flat', k);
//	~atkPeriRegDecoderGainMatrix.value(5, pi/3, 'point', k);

//	~atkQuadDecoderGainMatrix.value(0.5, k);
//	~atkQuadDecoderGainMatrix.value(pi/3, k);

)


// -------------------------------------------------------
// test (decoders)

// -------------------------------------------------------
// pantoSBR
// mono pink noise, sound synth function
(
{
	var fl, bl, br, fr,					// for quad output
		azim, bSig;
	var numSpeakers;

	// constants
	numSpeakers = 4;

	// kr....
	azim = MouseX.kr(pi, -pi);


	// synth and process
	bSig = AtkMonoToB.ar(PinkNoise.ar, azim); // mono pink noise

	// decode to output
	#fl, bl, br, fr = ~atkPantoSBR.value(bSig, numSpeakers);
	
	[fl, fr, bl, br]		// reorder output to match speaker arrangement

}.scope;
)

// periSBR
// mono pink noise, sound synth function
(
{
	var flu, blu, bru, fru,					// for 3D cube output
		fld, bld, brd, frd,
		azim, bSig;
	var numSpeakerPairs;
	
	// constants
	numSpeakerPairs = 4;

	// kr....
	azim = MouseX.kr(pi, -pi);


	// synth and process
	bSig = AtkMonoToB.ar(PinkNoise.ar, azim); // mono pink noise

	// decode to output
	#flu, blu, bru, fru, fld, bld, brd, frd = ~atkPeriSBR.value(bSig, numSpeakerPairs);
	
	[flu, fru, blu, bru, fld, frd, bld, brd]	// reorder output to match speaker arrangement

}.scope;
)


// -------------------------------------------------------
// decodeSBD (2D, narrow-quad)
// mono pink noise, sound synth function
(
{
	var fl, bl, br, fr,					// for quad output
		azim, bSig;
	var angle;
	var positions;
	

	// constants
	angle = (30/180) * pi;
	positions = Matrix.with([[angle.cos, angle.sin], [angle.cos, angle.sin.neg]]);


	// kr....
	azim = MouseX.kr(pi, -pi);


	// synth and process
	bSig = AtkMonoToB.ar(PinkNoise.ar, azim); // mono pink noise

	// decode to output
	#fl, fr, br, bl = ~atkDecodeSBD.value(bSig, positions);
	
	[fl, fr, bl, br]		// reorder output to match speaker arrangement

}.scope;
)


// -------------------------------------------------------
// decodeSBD (3D, cube)
// mono pink noise, sound synth function
(
{
	var flu, blu, bru, fru,					// for 3D cube output
		fld, bld, brd, frd,
		azim, bSig;
	var positions;
	var k;
	

	// constants
	positions = Matrix.with([[1, 1, 1], [-1, 1, 1], [-1, -1, 1], [1, -1, 1]]);
	k = 3.sqrt.reciprocal;					// energy for 3D


	// kr....
	azim = MouseX.kr(pi, -pi);


	// synth and process
	bSig = AtkMonoToB.ar(PinkNoise.ar, azim); // mono pink noise

	// decode to output
	#flu, blu, bru, fru, brd, frd, fld, bld = ~atkDecodeSBD.value(bSig, positions);
	
	[flu, fru, blu, bru, fld, frd, bld, brd]	// reorder output to match speaker arrangement

}.scope;
)


// -------------------------------------------------------
// quadSBD
// mono pink noise, sound synth function
(
{
	var fl, bl, br, fr,					// for quad output
		azim, bSig;
	var angle;


	// constants
	angle = (30 / 180) * pi;


	// kr....
	azim = MouseX.kr(pi, -pi);


	// synth and process
	bSig = AtkMonoToB.ar(PinkNoise.ar, azim); // mono pink noise

	// decode to output
	#fl, bl, br, fr = ~atkQuadSBD.value(bSig, angle);
	
	[fl, fr, bl, br]		// reorder output to match speaker arrangement

}.scope;
)


// -------------------------------------------------------
// -------------------------------------------------------
// -------------------------------------------------------
// pantoDBR
// mono pink noise, sound synth function
// including NFC (called distance!)
(
{
	var fl, bl, br, fr,					// for quad output
		azim, bSig;
	var numSpeakers, distance;

	// constants
	numSpeakers = 4;
	distance = 1.1;				// speaker distance (radius) in meters
	

	// kr....
	azim = MouseX.kr(pi, -pi);


	// synth and process
	bSig = AtkMonoToB.ar(PinkNoise.ar, azim); // mono pink noise


	// split...
	#w, x, y, z = bSig;
	
	// nearfield compensate
	bSig = AtkDistance.ar(w, x, y, z, distance);

	// decode to output
	#fl, bl, br, fr = ~atkPantoDBR.value(bSig, numSpeakers);
	
	[fl, fr, bl, br]		// reorder output to match speaker arrangement

}.scope;
)


// periDBR
// mono pink noise, sound synth function
// including NFC (called distance!)
(
{
	var flu, blu, bru, fru,					// for 3D cube output
		fld, bld, brd, frd,
		azim, bSig;
	var numSpeakerPairs, distance;
	
	// constants
	numSpeakerPairs = 4;
	distance = 1.1;				// speaker distance (radius) in meters

	// kr....
	azim = MouseX.kr(pi, -pi);


	// synth and process
	bSig = AtkMonoToB.ar(PinkNoise.ar, azim); // mono pink noise


	// split...
	#w, x, y, z = bSig;
	
	// nearfield compensate
	bSig = AtkDistance.ar(w, x, y, z, distance);

	// decode to output
	#flu, blu, bru, fru, fld, bld, brd, frd = ~atkPeriDBR.value(bSig, numSpeakerPairs);
	
	[flu, fru, blu, bru, fld, frd, bld, brd]	// reorder output to match speaker arrangement

}.scope;
)


// -------------------------------------------------------
// decodeDBD (2D, narrow-quad)
// mono pink noise, sound synth function
// including NFC (called distance!)
(
{
	var fl, bl, br, fr,					// for quad output
		azim, bSig;
	var angle;
	var positions, distance;
	

	// constants
	angle = (30/180) * pi;
	positions = Matrix.with([[angle.cos, angle.sin], [angle.cos, angle.sin.neg]]);
	distance = 1.1;				// speaker distance (radius) in meters


	// kr....
	azim = MouseX.kr(pi, -pi);


	// synth and process
	bSig = AtkMonoToB.ar(PinkNoise.ar, azim); // mono pink noise


	// split...
	#w, x, y, z = bSig;
	
	// nearfield compensate
	bSig = AtkDistance.ar(w, x, y, z, distance);

	// decode to output
	#fl, fr, br, bl = ~atkDecodeDBD.value(bSig, positions);
	
	[fl, fr, bl, br]		// reorder output to match speaker arrangement

}.scope;
)


// -------------------------------------------------------
// quadDBD
// mono pink noise, sound synth function
// including NFC (called distance!)
(
{
	var fl, bl, br, fr,					// for quad output
		azim, bSig;
	var angle, distance;


	// constants
	angle = (30 / 180) * pi;
	distance = 1.1;				// speaker distance (radius) in meters


	// kr....
	azim = MouseX.kr(pi, -pi);


	// synth and process
	bSig = AtkMonoToB.ar(PinkNoise.ar, azim); // mono pink noise


	// split...
	#w, x, y, z = bSig;
	
	// nearfield compensate
	bSig = AtkDistance.ar(w, x, y, z, distance);

	// decode to output
	#fl, bl, br, fr = ~atkQuadDBD.value(bSig, angle);
	
	[fl, fr, bl, br]		// reorder output to match speaker arrangement

}.scope;
)
