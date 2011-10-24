// examples....


// define encoder / decoder matrices

~encoder = FOAEncoderMatrix.newStereo(pi/4)		// kernel options available
~decoder = FOADecoderMatrix.newQuad(pi/6, 'single')	// kernel options available
~decoder = FOADecoderMatrix.newQuad(pi/6, 'dual')	// kernel options available

~decoder = FOADecoderKernel.newUHJ	// kernel options available
~decoder.kind
~decoder.free


~encoder = FOAEncoderKernel.newSuper
~encoder.kernel
~encoder.free

~transMatrix = FOATransMatrix.newRotate(pi/2)	// kernel options available (TBD)
~transMatrix = FOATransMatrix.newRotate		// kernel options available (TBD)

//	NAME					SPARSE?		UGEN MODE?	NOTES
// ------------------------------------------------------------
//
//	mirrorO				x
//	mirrorX				x
//	mirrorY				x
//	mirrorZ				x
//	mirror(theta, phi)					x			name?
//
//	rotate(angle)						x
//	tilt(angle)						x
//	tumble(angle)						x
//	rtt(rotAngle, tilAngle, tumAngle)	x
//
//	directO(angle)		x			x
//	directX(angle)		x			x
//	directY(angle)		x			x
//	directZ(angle)		x			x
//	direct(angle, theta, phi)			x			include?
//
//	dominateX(gain)					x
//	dominateY(gain)					x
//	dominateZ(gain)					x
//	dominate(gain, theta, phi)			x
//
//	zoomX(angle)						x
//	zoomY(angle)						x
//	zoomZ(angle)						x
//	zoom(angle, theta, phi)				x
//
//	focusX(angle)						x
//	focusY(angle)						x
//	focusZ(angle)						x
//	focus(angle, theta, phi)			x
//
//	pushX(angle)						x
//	pushY(angle)						x
//	pushZ(angle)						x
//	push(angle, theta, phi)				x
//
//	pressX(angle)						x
//	pressY(angle)						x
//	pressZ(angle)						x
//	press(angle, theta, phi)			x
//
//	asymmetry(angle)					x			include?
//	balance(angle)					x			a synonym for...?
//	scatter(parms...)								scatter matrix, best in a-format


// what about filters? these are transforms too!!
//
//	NFC								x
//	proximity							x
//	psycho							x


// ------------------------------------------------------------
// encode
sig = AtkEncode.ar(Pan2.ar(sig, azim.neg/pi), ~encodeMatrix);


// ------------------------------------------------------------
// transform
sig = AtkTransform.ar(sig, ~transMatrix);		// tiv
sig = AtkTransform.ar(sig, ~transMatrix, tvArg0, tvArg1, tvArg2); // tv

	// or 
sig = AtkTransform.ar('tras', sig, arg0, arg1, arg2, arg3, arg4);
	


// decode to output
outSig = AtkDecode.ar(sig, ~decoder);





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


//	// ------------------------------------------------------------
//	// encode
//	sig = FOAEncode.ar(Pan2.ar(sig, azim.neg/pi), ~encodeMatrix);
//	
//
//	// ------------------------------------------------------------
//	// distance filter
////	sig = AtkDistance.ar(sig.at(0), sig.at(1), sig.at(2), sig.at(3), 1);
//
//	// decode to output
//	#fl, bl, br, fr = FOADecode.ar(sig, ~decoder);
//
//	[fl, fr, bl, br]		// reorder output to match speaker arrangement

}.scope;
)


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
	sig = FOAEncode.ar(Pan2.ar(sig, azim.neg/pi), ~encoder);
//	sig = FOAEncode.ar(Pan2.ar(sig, azim.neg/pi), ~encoder);
//	sig = ATKMatrixMix.ar(Pan2.ar(sig, azim.neg/pi), ~encoder.matrix);
//	sig = ATKKernelConv.ar(Pan2.ar(sig, azim.neg/pi), ~encoder.kernel);
	

	// ------------------------------------------------------------
	// distance filter
	sig = AtkDistance.ar(sig.at(0), sig.at(1), sig.at(2), sig.at(3), 1);

	// decode to output
//	#fl, bl, br, fr = FOADecode.ar(sig, ~decoder);
////	#fl, bl, br, fr = ATKMatrixMix.ar(sig, ~decoder.matrix);
////	#fl, bl, br, fr = ATKKernelConv.ar(sig, ~decoder.kernel);
//
//	[fl, fr, bl, br]		// reorder output to match speaker arrangement

	// uhj
	FOADecode.ar(sig, ~decoder);


}.scope;
)
