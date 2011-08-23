// -------------------------------------------------------
// ATK (Sandbox)
//
// Coded by Joseph Anderson 2011
//
// Stereo encoding and decoding
// -------------------------------------------------------



// -------------------------------------------------------
// Stereo encoders
//
//   uhj_to_b            "Ambisonic Decoders for HDTV" (1992)
//   superstereo
//   simplestereo
//
//------------------------------------------------------------------------

// simple stereo

// NOTE: angle is framed as the 'distortion' angle rather than width,
//		so that behaviour matches that of AtkFocusX.
//		Angle = 0 gives left and right at 'hard' left and right
(
	~atkSimpleStereo = { arg sterSig, angle = pi/4;
		var l, r, w, x, y, z;
		
		#l, r = sterSig;
		
		[ 0.5 * (l + r), 1/2.sqrt * sin(angle) * (l + r),
			1/2.sqrt * cos(angle) * (l - r), Silent.ar(1) ]
	}
)


//------------------------------------------------------------------------
// Stereo decoders
//
//   b_to_uhj            "Ambisonic Decoders for HDTV" (1992)
//   b_to_stereo         virtual stereo microphone decoding
//   b_to_binaural       HRTF decoding
//
//------------------------------------------------------------------------

// b_to_uhj(b, decoder_kernels)// SEE: UHU.scd


//   b_to_stereo(a, angle = 0.7854, k = 1.0)
(
	~atkBToStereo = { arg bSig, angle = 0.7854, k = 1.0;
		var decoder;

		//return decoder gain matrix
		decoder = [	        [(1-k) * 2.sqrt, k * angle.cos, k * angle.sin],	        [(1-k) * 2.sqrt, k * angle.cos, k.neg * angle.sin]        ];

		// decode here!
		decoder.asArray.collect({arg channel;
			(channel * bSig[..2]).sum			// panto, need 1st 3 bSig channels
			})


	}
)


// b_to_binaural(b, decoder_kernels)
// SEE: binaural.scd




// -------------------------------------------------------
// test (decoders)


// -------------------------------------------------------
// b_to_stereo
// mono pink noise, sound synth function
(
{
	var azim, bSig;
	var angle, pattern;

	// constants
	angle	= 0.7854;				// +-45 deg (Blumlein)	pattern	= 1.0;				// fig-8//	angle	= 110/2 * pi/180;		// +-55 deg (Hyper-cardioid coincident)//	pattern	= 0.75;				// hyper-cardioid//	angle	= 135/2 * pi/180;		// +-67.5 deg (Cardioid coincident)//	pattern	= 0.5;				// cardioid//	angle	= pi/3;				// +-60 deg (Uniform reverb)//	pattern	= 0.5;				// cardioid//	angle	= pi/2;				// +-90 deg (Sphere)//	pattern	= 0.5;				// cardioid	


	// kr....
	azim = MouseX.kr(pi, -pi);


	// synth and process
	bSig = AtkMonoToB.ar(PinkNoise.ar, azim); // mono pink noise

	// decode to output
	~atkBToStereo.value(bSig, angle, pattern);
	
}.scope;
)

