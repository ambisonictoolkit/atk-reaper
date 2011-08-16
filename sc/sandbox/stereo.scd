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

// b_to_uhj(b, decoder_kernels)//
