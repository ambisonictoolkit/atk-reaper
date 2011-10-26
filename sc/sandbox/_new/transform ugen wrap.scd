// examples....


// define encoder / decoder matrices

~encoder = FoaEncoderMatrix.newDirection
~decoder = FoaDecoderMatrix.newQuad(pi/6, 'single')
~decoder = FoaDecoderMatrix.newQuad(pi/6, 'dual')

~xformer = FoaXformerMatrix.newRotate(pi/2)
~xformer = FoaXformerMatrix.newRotate

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
~encoder = FoaEncoderMatrix.newDirection
~encoder = FoaEncoderMatrix.newDirection(pi.neg/2)
~decoder = FoaDecoderMatrix.newQuad(pi/6, 'single')
~decoder = FoaDecoderMatrix.newQuad(pi/6, 'dual')

(
{
	var fl, bl, br, fr;					// for quad output
	var azim, mouseX, mouseY, sig;

	// kr.... for SC3's internal encoder
	azim = MouseX.kr(pi, -pi);
//	mouseX = MouseX.kr(pi/2, -pi/2);
	mouseY = MouseY.kr(1, -1);			// seems to be inverted from documentation
//	mouseY = MouseY.kr(pi/2, -pi/2);		// seems to be inverted from documentation


	// ------------------------------------------------------------
	// test sig
	sig = PinkNoise.ar; 					// mono pink noise


	// ------------------------------------------------------------
	// encode
	sig = FoaEncode.ar(sig, ~encoder);
//	sig = FoaPanB.ar(sig);


	// ------------------------------------------------------------
	// transform

//	sig = FoaDirectO.ar(sig.at(0), sig.at(1), sig.at(2), sig.at(3),
//		mouseY);
//	sig = FoaDirectX.ar(sig.at(0), sig.at(1), sig.at(2), sig.at(3),
//		mouseY);
//	sig = FoaDirectY.ar(sig.at(0), sig.at(1), sig.at(2), sig.at(3),
//		mouseY);
//	sig = FoaDirectZ.ar(sig.at(0), sig.at(1), sig.at(2), sig.at(3),
//		mouseY);
//	sig = FoaAsymmetry.ar(sig.at(0), sig.at(1), sig.at(2), sig.at(3),
//		mouseX);


//	sig = FoaXform.ar(sig, ~xformer);
	sig = FoaTransform.ar(sig, 'rotate', azim);
//	sig = FoaTransform.ar(sig, 'rotate', [azim, mul: 6.dbamp]);
//	sig = FoaTransform.ar(sig, 'rotate', azim, 0.dbamp);
//	sig = FoaTransform.ar(sig, 'tilt', azim);
//	sig = FoaTransform.ar(sig, 'tumble', azim);

//	sig = FoaTransform.ar(sig, 'dominateX', 6.neg);
//	sig = FoaTransform.ar(sig, 'dominateX', mouseY * 6);
//	sig = FoaTransform.ar(sig, 'dominateY', mouseY * 6);

//	sig = FoaTransform.ar(sig, 'directO', mouseY * pi/2);
//	sig = FoaTransform.ar(sig, 'directX', mouseY * pi/2);
//	sig = FoaTransform.ar(sig, 'directY', mouseY * pi/2);
//	sig = FoaTransform.ar(sig, 'directZ', mouseY * pi/2);


//	sig = FoaTransform.ar(sig, 'zoomX', pi/2);
//	sig = FoaTransform.ar(sig, 'zoomX', mouseY * pi/2);
//	sig = FoaTransform.ar(sig, 'zoomY', mouseY * pi/2);
//	sig = FoaTransform.ar(sig, 'zoomZ', mouseY * pi/2);

//	sig = FoaTransform.ar(sig, 'focusX', mouseY * pi/2);
//	sig = FoaTransform.ar(sig, 'focusY', mouseY * pi/2);
//	sig = FoaTransform.ar(sig, 'focusZ', mouseY * pi/2);

//	sig = FoaTransform.ar(sig, 'pushX', mouseY * pi/2);
//	sig = FoaTransform.ar(sig, 'pushY', mouseY * pi/2);
//	sig = FoaTransform.ar(sig, 'pushZ', mouseY * pi/2);

//	sig = FoaTransform.ar(sig, 'pressX', mouseY * pi/2);
//	sig = FoaTransform.ar(sig, 'pressY', mouseY * pi/2);
//	sig = FoaTransform.ar(sig, 'pressZ', mouseY * pi/2);

//	sig = FoaTransform.ar(sig, 'asymmetry', mouseY * pi/2);

//	sig = FoaTransform.ar(sig, 'balance', mouseY * pi/2);

//	sig = FoaTransform.ar(sig, 'rtt', azim, mouseY);
//	sig = FoaTransform.ar(sig, 'rtt', [azim, tumAngle: mouseY]);

//	sig = FoaTransform.ar(sig, 'dominate', mouseY * 6, azim);
//	sig = FoaTransform.ar(sig, 'zoom', mouseY, azim);
//	sig = FoaTransform.ar(sig, 'focus', mouseY, azim);
//	sig = FoaTransform.ar(sig, 'push', mouseY, azim);
//	sig = FoaTransform.ar(sig, 'press', mouseY, azim);
//	sig = FoaTransform.ar(sig, 'direct', mouseY, azim);
//	sig = FoaTransform.ar(sig, 'mirror', mouseY, azim);


	// ------------------------------------------------------------
	// distance filter
//	sig = FoaDistance.ar(sig.at(0), sig.at(1), sig.at(2), sig.at(3), 1);
	sig = FoaNFC.ar(sig.at(0), sig.at(1), sig.at(2), sig.at(3), 1);

	// decode to output
	#fl, bl, br, fr = FoaDecode.ar(sig, ~decoder);
	[fl, fr, bl, br]		// reorder output to match speaker arrangement

}.scope;
)
