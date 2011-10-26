// add comments here....

FOAPanB : Panner {
	
	*ar { arg in, azimuth=0, elevation=0;
		^this.multiNew('audio', in, azimuth, elevation )
	}
	
	init { arg ... theInputs;
		inputs = theInputs;		
		channels = [ OutputProxy(\audio,this,0), OutputProxy(\audio,this,1),
					OutputProxy(\audio,this,2), OutputProxy(\audio,this,3) ];
		^channels
	}
}


Atk : Panner {
		
	init { arg ... theInputs;
		inputs = theInputs;		
		channels = [ OutputProxy(\audio,this,0), OutputProxy(\audio,this,1),
					OutputProxy(\audio,this,2), OutputProxy(\audio,this,3) ];
		^channels
	}
	
 	checkInputs { ^this.checkNInputs(4) }
 	
 	}

	
AtkDirectO : Atk {			// check matrix!!, appears to be the old matrix
	*ar { arg w, x, y, z, angle = pi/2, mul = 1, add = 0;
		^this.multiNew('audio', w, x, y, z, angle).madd(mul, add);
	}
}


AtkDirectX : Atk {
	*ar { arg w, x, y, z, angle = pi/2, mul = 1, add = 0;
		^this.multiNew('audio', w, x, y, z, angle).madd(mul, add);
	}
}

AtkDirectY : AtkDirectX { }
AtkDirectZ : AtkDirectX { }

AtkRotate : Atk { 
	*ar { arg w, x, y, z, angle = 0, mul = 1, add = 0;
		^this.multiNew('audio', w, x, y, z, angle).madd(mul, add);
	}
} 
AtkTilt : AtkRotate { }
AtkTumble : AtkRotate { }

AtkFocusX : AtkRotate { }
AtkFocusY : AtkRotate { }
AtkFocusZ : AtkRotate { }

AtkPushX : AtkRotate { }
AtkPushY : AtkRotate { }
AtkPushZ : AtkRotate { }

AtkPressX : AtkRotate { }
AtkPressY : AtkRotate { }
AtkPressZ : AtkRotate { }

AtkZoomX : AtkRotate { }
AtkZoomY : AtkRotate { }
AtkZoomZ : AtkRotate { }


AtkDominateX : Atk {	
	*ar { arg w, x, y, z, gain = 0;
		^this.multiNew('audio', w, x, y, z, gain);
	}	
}

AtkDominateY : AtkDominateX { }
AtkDominateZ : AtkDominateX { }

AtkAsymmetry : AtkRotate { }


FOARTT { 
	*ar { arg w, x, y, z, rotAngle = 0, tilAngle = 0, tumAngle = 0, mul = 1, add = 0;

		#w, x, y, z = AtkRotate.ar(w, x, y, z, rotAngle);
		#w, x, y, z = AtkTilt.ar(w, x, y, z, tilAngle);
		^AtkTumble.ar(w, x, y, z, tumAngle, mul, add);
	}
} 

FOAMirror { 
	*ar { arg w, x, y, z, theta = 0, phi = 0, mul = 1, add = 0;

		#w, x, y, z = AtkRotate.ar(w, x, y, z, theta.neg);
		#w, x, y, z = AtkTumble.ar(w, x, y, z, phi.neg);
		#w, x, y, z = FOAXform.ar([w, x, y, z], FOAXformerMatrix.newMirrorX);
		#w, x, y, z = AtkTumble.ar(w, x, y, z, phi);
		^AtkRotate.ar(w, x, y, z, theta, mul, add);
	}
} 

FOADirect { 
	*ar { arg w, x, y, z, angle = 0, theta = 0, phi = 0, mul = 1, add = 0;

		#w, x, y, z = AtkRotate.ar(w, x, y, z, theta.neg);
		#w, x, y, z = AtkTumble.ar(w, x, y, z, phi.neg);
		#w, x, y, z = AtkSquishX.ar(w, x, y, z, angle); // rename to DirectX
		#w, x, y, z = AtkTumble.ar(w, x, y, z, phi);
		^AtkRotate.ar(w, x, y, z, theta, mul, add);
	}
} 

FOADominate { 
	*ar { arg w, x, y, z, gain = 0, theta = 0, phi = 0, mul = 1, add = 0;

		#w, x, y, z = AtkRotate.ar(w, x, y, z, theta.neg);
		#w, x, y, z = AtkTumble.ar(w, x, y, z, phi.neg);
		#w, x, y, z = AtkDominateX.ar(w, x, y, z, gain);
		#w, x, y, z = AtkTumble.ar(w, x, y, z, phi);
		^AtkRotate.ar(w, x, y, z, theta, mul, add);
	}
} 

FOAZoom { 
	*ar { arg w, x, y, z, angle = 0, theta = 0, phi = 0, mul = 1, add = 0;

		#w, x, y, z = AtkRotate.ar(w, x, y, z, theta.neg);
		#w, x, y, z = AtkTumble.ar(w, x, y, z, phi.neg);
		#w, x, y, z = AtkZoomX.ar(w, x, y, z, angle);
		#w, x, y, z = AtkTumble.ar(w, x, y, z, phi);
		^AtkRotate.ar(w, x, y, z, theta, mul, add);
	}
} 

FOAFocus { 
	*ar { arg w, x, y, z, angle = 0, theta = 0, phi = 0, mul = 1, add = 0;

		#w, x, y, z = AtkRotate.ar(w, x, y, z, theta.neg);
		#w, x, y, z = AtkTumble.ar(w, x, y, z, phi.neg);
		#w, x, y, z = AtkFocusX.ar(w, x, y, z, angle);
		#w, x, y, z = AtkTumble.ar(w, x, y, z, phi);
		^AtkRotate.ar(w, x, y, z, theta, mul, add);
	}
} 

FOAPush { 
	*ar { arg w, x, y, z, angle = 0, theta = 0, phi = 0, mul = 1, add = 0;

		#w, x, y, z = AtkRotate.ar(w, x, y, z, theta.neg);
		#w, x, y, z = AtkTumble.ar(w, x, y, z, phi.neg);
		#w, x, y, z = AtkPushX.ar(w, x, y, z, angle);
		#w, x, y, z = AtkTumble.ar(w, x, y, z, phi);
		^AtkRotate.ar(w, x, y, z, theta, mul, add);
	}
} 

FOAPress { 
	*ar { arg w, x, y, z, angle = 0, theta = 0, phi = 0, mul = 1, add = 0;

		#w, x, y, z = AtkRotate.ar(w, x, y, z, theta.neg);
		#w, x, y, z = AtkTumble.ar(w, x, y, z, phi.neg);
		#w, x, y, z = AtkPressX.ar(w, x, y, z, angle);
		#w, x, y, z = AtkTumble.ar(w, x, y, z, phi);
		^AtkRotate.ar(w, x, y, z, theta, mul, add);
	}
} 


//------------------------------------------------------------------------
// Filters

AtkProximity : Atk { 
	*ar { arg w, x, y, z, distance = 1, mul = 1, add = 0;
		^this.multiNew('audio', w, x, y, z, distance).madd(mul, add);
	}

}

AtkNFC : Atk { 
	*ar { arg w, x, y, z, distance = 1, mul = 1, add = 0;
		^this.multiNew('audio', w, x, y, z, distance).madd(mul, add);
	}
		
}

AtkPsychoShelf : Atk { 
	*ar { arg w, x, y, z, freq = 400, k0 = (3/2).sqrt, k1 = 3.sqrt/2, mul = 1, add = 0;
		^this.multiNew('audio', w, x, y, z, freq, k0, k1).madd(mul, add);
	}
		
}


//------------------------------------------------------------------------
// ATKMatrixMix & ATKKernelConv

ATKMatrixMix {
	*ar { arg in, matrix, mul = 1, add = 0;
		
		var out;

		// wrap input as array if needed, for mono inputs
		in.isArray.not.if({ in = [in] });
		
		out = Mix.fill( matrix.cols, { arg i; // fill input
			UGen.replaceZeroesWithSilence(
				matrix.flop.asArray.at(i) * in.at(i)
			)
		});

		^out.madd(mul, add)
	}
}

ATKKernelConv {
	*ar { arg in, kernel, mul = 1, add = 0;
		
		var out;

		// wrap input as array if needed, for mono inputs
		in.isArray.not.if({ in = [in] });
		
		out = Mix.ar(
			kernel.shape.at(0).collect({ arg i;
				kernel.shape.at(1).collect({ arg j;
					Convolution2.ar(
						in.at(i),
						kernel.at(i).at(j),
						framesize: kernel.at(i).at(j).numFrames
					)
				})
			})
		);

		^out.madd(mul, add)
	}
}


//------------------------------------------------------------------------
// Decoder built using ATKMatrixMix & ATKKernelConv

FOADecode {
	*ar { arg in, decoder, mul = 1, add = 0;

		switch ( decoder.class, 

			FOADecoderMatrix, {

				if ( decoder.shelfFreq.isNumber, { // shelf filter?
					in = AtkPsychoShelf.ar(in.at(0), in.at(1), in.at(2), in.at(3),
						decoder.shelfFreq, decoder.shelfK.at(0), decoder.shelfK.at(1))
				});

				^ATKMatrixMix.ar(in, decoder.matrix, mul, add)
			},
			
			FOADecoderKernel, {
				^ATKKernelConv.ar(in, decoder.kernel, mul, add)
			}
		)
	}
}


//------------------------------------------------------------------------
// Encoder built using ATKMatrixMix & ATKKernelConv

FOAEncode {
	*ar { arg in, encoder, mul = 1, add = 0;
		
		var out;

		switch ( encoder.class, 

			FOAEncoderMatrix, {
				out = ATKMatrixMix.ar(in, encoder.matrix, mul, add)
			},
			
			FOAEncoderKernel, {
				out = ATKKernelConv.ar(in, encoder.kernel, mul, add)
			}
		);

		if ( out.size < 4, {			// 1st order, fill missing harms with zeros
			out = out ++ Silent.ar(4 - out.size)
		});
		
		^out
	}
}


//------------------------------------------------------------------------
// Transformer built using ATKMatrixMix & ATKKernelConv

FOAXform {
	*ar { arg in, xformer, mul = 1, add = 0;
		
		var out;

//		switch ( xformer.class,
//
//			FOAXformerMatrix, {
//				out = ATKMatrixMix.ar(in, xformer.matrix, mul, add)
//			},
//			
//			FOAXformerKernel, {
//				out = ATKKernelConv.ar(in, xformer.kernel, mul, add)
//			}
//		);
//
//		^out

		// for now...
		^ATKMatrixMix.ar(in, xformer.matrix, mul, add)
	}
}


//------------------------------------------------------------------------
// Transformer: UGen wrapper

FOATransform {
	*ar { arg in, kind ... args;
		
		var argDict, argDefaults;
		var ugen;
		
		argDict = { arg ugen, args, argDefaults;
			var index, userDict;
			var ugenKeys;
			var ugenDict;
			
			// find index dividing ordered and named args
			index = args.detectIndex({arg item; item.isKindOf(Symbol)});
		
			// find ugen args, drop [ 'this', w, x, y, z ]
			ugenKeys = ugen.class.findRespondingMethodFor(\ar).argNames.drop(5);
		
			ugenDict = Dictionary.new;
			ugenKeys.do({arg key, i; ugenDict.put(key, argDefaults.at(i))});
			
			// build user dictionary
			userDict = Dictionary.new(ugenKeys.size);
			(index == nil).not.if({
				userDict = userDict.putAll(Dictionary.newFrom(args[index..]));
			}, {
				index = args.size;
			});
			userDict = userDict.putAll(Dictionary.newFrom((index).collect({arg i;
				[ugenKeys.at(i), args.at(i)]}).flat));
				
			// merge
			ugenDict.merge(userDict, {
				arg ugenArg, userArg; (userArg != nil).if({userArg})
			})
		};
		

		switch ( kind, 

			'rotate', {

				ugen = AtkRotate;
				argDefaults = [0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\angle), argDict.at(\mul), argDict.at(\add)
				)
			},
			
			'tilt', {

				ugen = AtkTilt;
				argDefaults = [0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\angle), argDict.at(\mul), argDict.at(\add)
				)
			},
			
			'tumble', {

				ugen = AtkTumble;
				argDefaults = [0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\angle), argDict.at(\mul), argDict.at(\add)
				)
			},
				
			'directO', {

				ugen = AtkDirectO;
				argDefaults = [0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\angle), argDict.at(\mul), argDict.at(\add)
				)
			},

			'directX', {

				ugen = AtkDirectX;
				argDefaults = [0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\angle), argDict.at(\mul), argDict.at(\add)
				)
			},

			'directY', {

				ugen = AtkDirectY;
				argDefaults = [0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\angle), argDict.at(\mul), argDict.at(\add)
				)
			},

			'directZ', {

				ugen = AtkDirectZ;
				argDefaults = [0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\angle), argDict.at(\mul), argDict.at(\add)
				)
			},

			'dominateX', {

				ugen = AtkDominateX;
				argDefaults = [0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\gain), argDict.at(\mul), argDict.at(\add)
				)
			},
			
			'dominateY', {

				ugen = AtkDominateY;
				argDefaults = [0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\gain), argDict.at(\mul), argDict.at(\add)
				)
			},

			'dominateZ', {

				ugen = AtkDominateZ;
				argDefaults = [0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\gain), argDict.at(\mul), argDict.at(\add)
				)
			},

			'zoomX', {

				ugen = AtkZoomX;
				argDefaults = [0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\angle), argDict.at(\mul), argDict.at(\add)
				)
			},
			
			'zoomY', {

				ugen = AtkZoomY;
				argDefaults = [0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\angle), argDict.at(\mul), argDict.at(\add)
				)
			},

			'zoomZ', {

				ugen = AtkZoomZ;
				argDefaults = [0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\angle), argDict.at(\mul), argDict.at(\add)
				)
			},
			
			'focusX', {

				ugen = AtkFocusX;
				argDefaults = [0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\angle), argDict.at(\mul), argDict.at(\add)
				)
			},

			'focusY', {

				ugen = AtkFocusY;
				argDefaults = [0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\angle), argDict.at(\mul), argDict.at(\add)
				)
			},

			'focusZ', {

				ugen = AtkFocusZ;
				argDefaults = [0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\angle), argDict.at(\mul), argDict.at(\add)
				)
			},

			'pushX', {

				ugen = AtkPushX;
				argDefaults = [0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\angle), argDict.at(\mul), argDict.at(\add)
				)
			},

			'pushY', {

				ugen = AtkPushY;
				argDefaults = [0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\angle), argDict.at(\mul), argDict.at(\add)
				)
			},

			'pushZ', {

				ugen = AtkPushZ;
				argDefaults = [0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\angle), argDict.at(\mul), argDict.at(\add)
				)
			},

			'pressX', {

				ugen = AtkPressX;
				argDefaults = [0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\angle), argDict.at(\mul), argDict.at(\add)
				)
			},

			'pressY', {

				ugen = AtkPressY;
				argDefaults = [0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\angle), argDict.at(\mul), argDict.at(\add)
				)
			},

			'pressZ', {

				ugen = AtkPressZ;
				argDefaults = [0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\angle), argDict.at(\mul), argDict.at(\add)
				)
			},

			'asymmetry', {

				ugen = AtkAsymmetry;
				argDefaults = [0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\angle), argDict.at(\mul), argDict.at(\add)
				)
			},

			'balance', {

				ugen = AtkZoomY;
				argDefaults = [0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\angle), argDict.at(\mul), argDict.at(\add)
				)
			},

			'rtt', {

				ugen = FOARTT;
				argDefaults = [0, 0, 0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\rotAngle), argDict.at(\tilAngle), argDict.at(\tumAngle),
					argDict.at(\mul), argDict.at(\add)
				)
			},
			
			'mirror', {

				ugen = FOAMirror;
				argDefaults = [0, 0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\theta), argDict.at(\phi),
					argDict.at(\mul), argDict.at(\add)
				)
			},

			'direct', {

				ugen = FOADirect;
				argDefaults = [0, 0, 0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\angle), argDict.at(\theta), argDict.at(\phi),
					argDict.at(\mul), argDict.at(\add)
				)
			},

			'dominate', {

				ugen = FOADominate;
				argDefaults = [0, 0, 0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\gain), argDict.at(\theta), argDict.at(\phi),
					argDict.at(\mul), argDict.at(\add)
				)
			},

			'zoom', {

				ugen = FOAZoom;
				argDefaults = [0, 0, 0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\angle), argDict.at(\theta), argDict.at(\phi),
					argDict.at(\mul), argDict.at(\add)
				)
			},

			'focus', {

				ugen = FOAFocus;
				argDefaults = [0, 0, 0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\angle), argDict.at(\theta), argDict.at(\phi),
					argDict.at(\mul), argDict.at(\add)
				)
			},

			'push', {

				ugen = FOAPush;
				argDefaults = [0, 0, 0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\angle), argDict.at(\theta), argDict.at(\phi),
					argDict.at(\mul), argDict.at(\add)
				)
			},

			'press', {

				ugen = FOAPress;
				argDefaults = [0, 0, 0, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\angle), argDict.at(\theta), argDict.at(\phi),
					argDict.at(\mul), argDict.at(\add)
				)
			},
			
			'nfc', {

				ugen = AtkNFC;
				argDefaults = [1, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\distance),
					argDict.at(\mul), argDict.at(\add)
				)
			},

			'proximity', {

				ugen = AtkProximity;
				argDefaults = [1, 1, 0];
				
				argDict = argDict.value(ugen, args, argDefaults);
				
				^ugen.ar(
					in.at(0), in.at(1), in.at(2), in.at(3),
					argDict.at(\distance),
					argDict.at(\mul), argDict.at(\add)
				)
			}
		)
	}
}
