// examples....


// define encoder / decoder matrices

~encoder = FOAEncoderMatrix.newDirection
~decoder = FOADecoderMatrix.newQuad(pi/6, 'single')
~decoder = FOADecoderMatrix.newQuad(pi/6, 'dual')

~xformer = FOAXformerMatrix.newRotate(pi/2)
~xformer = FOAXformerMatrix.newRotate

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
// quad (spec'd as narrow quad)
// mono pink noise, sound synth function
// mono encoder
~encoder = FOAEncoderMatrix.newDirection
~decoder = FOADecoderMatrix.newQuad(pi/6, 'single')
~decoder = FOADecoderMatrix.newQuad(pi/6, 'dual')

(
{
	var fl, bl, br, fr;					// for quad output
	var azim, mouseY, sig;

	// kr.... for SC3's internal encoder
	azim = MouseX.kr(pi, -pi);
	mouseY = MouseY.kr(1, -1);			// seems to be inverted from documentation
//	mouseY = MouseY.kr(pi/2, -pi/2);		// seems to be inverted from documentation


	// ------------------------------------------------------------
	// test sig
	sig = PinkNoise.ar; 					// mono pink noise


	// ------------------------------------------------------------
	// encode
	sig = FOAEncode.ar(sig, ~encoder);
//	sig = FOAPanB.ar(sig);


	// ------------------------------------------------------------
	// transform

//	sig = AtkDirectO.ar(sig.at(0), sig.at(1), sig.at(2), sig.at(3),
//		mouseY);
//	sig = AtkDirectX.ar(sig.at(0), sig.at(1), sig.at(2), sig.at(3),
//		mouseY);
//	sig = AtkDirectY.ar(sig.at(0), sig.at(1), sig.at(2), sig.at(3),
//		mouseY);
//	sig = AtkDirectZ.ar(sig.at(0), sig.at(1), sig.at(2), sig.at(3),
//		mouseY);


//	sig = FOAXform.ar(sig, ~xformer);
	sig = FOATransform.ar(sig, 'rotate', azim);
//	sig = FOATransform.ar(sig, 'rotate', [azim, mul: 6.dbamp]);
//	sig = FOATransform.ar(sig, 'rotate', azim, 0.dbamp);
//	sig = FOATransform.ar(sig, 'tilt', azim);
//	sig = FOATransform.ar(sig, 'tumble', azim);

//	sig = FOATransform.ar(sig, 'dominateX', 6.neg);
//	sig = FOATransform.ar(sig, 'dominateX', mouseY * 6);
//	sig = FOATransform.ar(sig, 'dominateY', mouseY * 6);

//	sig = FOATransform.ar(sig, 'directO', mouseY * pi/2);
//	sig = FOATransform.ar(sig, 'directX', mouseY * pi/2);
//	sig = FOATransform.ar(sig, 'directY', mouseY * pi/2);
//	sig = FOATransform.ar(sig, 'directZ', mouseY * pi/2);


//	sig = FOATransform.ar(sig, 'zoomX', pi/2);
//	sig = FOATransform.ar(sig, 'zoomX', mouseY * pi/2);
//	sig = FOATransform.ar(sig, 'zoomY', mouseY * pi/2);
//	sig = FOATransform.ar(sig, 'zoomZ', mouseY * pi/2);

//	sig = FOATransform.ar(sig, 'focusX', mouseY * pi/2);
//	sig = FOATransform.ar(sig, 'focusY', mouseY * pi/2);
//	sig = FOATransform.ar(sig, 'focusZ', mouseY * pi/2);

//	sig = FOATransform.ar(sig, 'pushX', mouseY * pi/2);
//	sig = FOATransform.ar(sig, 'pushY', mouseY * pi/2);
//	sig = FOATransform.ar(sig, 'pushZ', mouseY * pi/2);

//	sig = FOATransform.ar(sig, 'pressX', mouseY * pi/2);
//	sig = FOATransform.ar(sig, 'pressY', mouseY * pi/2);
//	sig = FOATransform.ar(sig, 'pressZ', mouseY * pi/2);

//	sig = FOATransform.ar(sig, 'asymmetry', mouseY * pi/2);

//	sig = FOATransform.ar(sig, 'balance', mouseY * pi/2);

//	sig = FOATransform.ar(sig, 'rtt', azim, mouseY);
//	sig = FOATransform.ar(sig, 'rtt', [azim, tumAngle: mouseY]);

//	sig = FOATransform.ar(sig, 'dominate', mouseY * 6, azim);
//	sig = FOATransform.ar(sig, 'zoom', mouseY, azim);
//	sig = FOATransform.ar(sig, 'focus', mouseY, azim);
//	sig = FOATransform.ar(sig, 'push', mouseY, azim);
//	sig = FOATransform.ar(sig, 'press', mouseY, azim);
//	sig = FOATransform.ar(sig, 'direct', mouseY, azim);
//	sig = FOATransform.ar(sig, 'mirror', mouseY, azim);


	// ------------------------------------------------------------
	// distance filter
//	sig = AtkDistance.ar(sig.at(0), sig.at(1), sig.at(2), sig.at(3), 1);
	sig = AtkNFC.ar(sig.at(0), sig.at(1), sig.at(2), sig.at(3), 1);

	// decode to output
	#fl, bl, br, fr = FOADecode.ar(sig, ~decoder);
	[fl, fr, bl, br]		// reorder output to match speaker arrangement

}.scope;
)
