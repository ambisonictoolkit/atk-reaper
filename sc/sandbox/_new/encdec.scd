// examples....


// define encoder / decoder matrices

~encoder = FoaEncoderMatrix.newStereo(pi/4)
~encoder = FoaEncoderKernel.newSuper				// remember to free!!
~encoder = FoaEncoderKernel.newGreatHall				// remember to free!!

~decoder = FoaDecoderMatrix.newQuad(pi/6, 'dual')
~decoder = FoaDecoderKernel.newUHJ				// remember to free!!

// inspect
~encoder.kind
~decoder.kind


~encoder.free				// free the kerne when finishedl!
~decoder.free				// free the kerne when finishedl!


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
//	foaSig = FoaEncode.ar(Pan2.ar(sig, azim.neg/pi), ~encoder);
	foaSig = FoaEncode.ar(sig, ~encoder);
	
	
	// test... out as array?
//	foaSig = foaSig.asAudioRateInput;
//	foaSig = foaSig.x;

	// ------------------------------------------------------------
	// distance filter
//	sig = AtkDistance.ar(sig.at(0), sig.at(1), sig.at(2), sig.at(3), 1);

	// decode to output (quad)
//	#fl, bl, br, fr = FoaDecode.ar(foaSig, ~decoder);
//
//	[fl, fr, bl, br]		// reorder output to match speaker arrangement

	// decode to output (UHJ)
//	FoaDecode.ar(foaSig, ~decoder);


}.scope;
)
