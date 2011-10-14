// examples....


// define encoder / decoder matrices

~encoder = FOAEncoderMatrix.newStereo(pi/4)
~encoder = FOAEncoderKernel.newSuper				// remember to free!!

~decoder = FOADecoderMatrix.newQuad(pi/6, 'dual')
~decoder = FOADecoderKernel.newUHJ				// remember to free!!

// inspect
~encoder.kind
~decoder.kind


~encoder.free				// free the kerne when finishedl!
~decoder.free				// free the kerne when finishedl!


//// ------------------------------------------------------------
//// encode
//sig = FOAEncode.ar(Pan2.ar(sig, azim.neg/pi), ~encodeMatrix);
//
//
//// ------------------------------------------------------------
//// transform
//sig = FOATransform.ar('tras', sig, arg0, arg1, arg2, arg3, arg4);
//	
//
//
//// decode to output
//outSig = FOADecode.ar(sig, ~decoderMatrix);





// ------------------------------------------------------------
// quad (spec'd as narrow quad)
// mono pink noise, sound synth function
// stereo encoder
(
{
	var fl, bl, br, fr;					// for quad output
	var azim, sig, foaSig;

	// kr.... for SC3's internal encoder
	azim = MouseX.kr(pi, -pi);


	// ------------------------------------------------------------
	// test sig
	sig = PinkNoise.ar; 					// mono pink noise


	// ------------------------------------------------------------
	// encode
	foaSig = FOAEncode.ar(Pan2.ar(sig, azim.neg/pi), ~encoder);
	
	
	// test... out as array?
//	foaSig = foaSig.asAudioRateInput;
//	foaSig = foaSig.x;

	// ------------------------------------------------------------
	// distance filter
//	sig = AtkDistance.ar(sig.at(0), sig.at(1), sig.at(2), sig.at(3), 1);

//	// decode to output (quad)
//	#fl, bl, br, fr = FOADecode.ar(foaSig, ~decoder);
//
//	[fl, fr, bl, br]		// reorder output to match speaker arrangement

	// decode to output (UHJ)
	FOADecode.ar(foaSig, ~decoder);


}.scope;
)
