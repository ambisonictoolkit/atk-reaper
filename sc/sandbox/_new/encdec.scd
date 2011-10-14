// examples....


// define encoder / decoder matrices

~encodeMatrix = FOAEncoderMatrix.newStereo(pi/4)		// kernel options available
~decoderMatrix = FOADecoderMatrix.newQuad(pi/6, 'dual')	// kernel options available


// ------------------------------------------------------------
// encode
sig = FOAEncode.ar(Pan2.ar(sig, azim.neg/pi), ~encodeMatrix);


// ------------------------------------------------------------
// transform
sig = FOATransform.ar('tras', sig, arg0, arg1, arg2, arg3, arg4);
	


// decode to output
outSig = FOADecode.ar(sig, ~decoderMatrix);





// ------------------------------------------------------------
// quad (spec'd as narrow quad)
// mono pink noise, sound synth function
// stereo encoder
(
{
	var fl, bl, br, fr;					// for quad output
	var azim, sig;

	// kr.... for SC3's internal encoder
	azim = MouseX.kr(pi, -pi);


	// ------------------------------------------------------------
	// test sig
	sig = PinkNoise.ar; 					// mono pink noise


	// ------------------------------------------------------------
	// encode
	sig = FOAEncode.ar(Pan2.ar(sig, azim.neg/pi), ~encodeMatrix);
	

	// ------------------------------------------------------------
	// distance filter
	sig = AtkDistance.ar(sig.at(0), sig.at(1), sig.at(2), sig.at(3), 1);

	// decode to output
	#fl, bl, br, fr = FOADecode.ar(sig, ~decoderMatrix);

	[fl, fr, bl, br]		// reorder output to match speaker arrangement

}.scope;
)
