// examples....


// define encoder / decoder matrices

~encodeMatrix = AtkEncoderMatrix.newStereo(pi/4)		// kernel options available
~decoderMatrix = AtkDecoderMatrix.newQuad(pi/6, 'dual')	// kernel options available

~transMatrix = AtkTransMatrix.newRotate(pi/2)	// kernel options available (TBD)
~transMatrix = AtkTransMatrix.newRotate		// kernel options available (TBD)

//	NAME					SPARSE?		UGEN MODE?	NOTES
// ------------------------------------------------------------
//
//	mirrorX				x
//	mirrorY				x
//	mirrorZ				x
//	mirrorO				x
//	mirror(theta, phi)					x			name?
//
//	rotate(angle)						x
//	tilt(angle)						x
//	tumble(angle)						x
//	rtt(rotAngle, tilAngle, tumAngle)	x
//
//	direct(angle)			x			x
//	directX(angle)		x			x
//	directY(angle)		x			x
//	directZ(angle)		x			x
//	directD(theta, phi)				x			include?
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
outSig = AtkDecode.ar(sig, ~decoderMatrix);





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
	sig = AtkEncode.ar(Pan2.ar(sig, azim.neg/pi), ~encodeMatrix);
	

	// ------------------------------------------------------------
	// distance filter
	sig = AtkDistance.ar(sig.at(0), sig.at(1), sig.at(2), sig.at(3), 1);

	// decode to output
	#fl, bl, br, fr = AtkDecode.ar(sig, ~decoderMatrix);

	[fl, fr, bl, br]		// reorder output to match speaker arrangement

}.scope;
)
