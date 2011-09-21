/* FOA wrappers */
// NOTE: we may wish to rename all ATK as FOA
// 		for clarity, and to assist with further development
//		of HOA implementations


FOA {
	var <w, <x, <y, <z;
	
	*new {arg w, x, y, z;
		^super.newCopyArgs(w, x, y, z);
	}
	
	*ar {arg w, x, y, z;
		^this.new(w, x, y, z);
	}
	
	madd {arg mul = 1, add = 0;
		^MulAdd.ar([w, x, y, z], mul, add);	
	}
	
	sig {arg mul = 1, add = 0;
		^[w, x, y, z] * mul + add;
	}
	
	asUGenInput {
		^[w, x, y, z];
	}
	
	asAudioRateInput {
		^[w, x, y, z];
	}
}

ATKMonoToFOA {
	*ar {arg in, azimuth = 0, elevation = 0;	
		var w, x, y, z;
		#w, x, y, z = AtkMonoToB.ar(in, azimuth, elevation);
		^FOA.ar(w, x, y, z);
	}	
}

ATKRotate {
	*ar {arg in, angle = 0, mul = 1, add = 0;
		var w, x, y, z;
		(in.isKindOf(FOA)).if({
			#w, x, y, z = AtkRotate.ar(in.w, in.x, in.y, in.z, angle, mul, add);
			^FOA.ar(w, x, y, z);
		})
	}
}

// ATKPantoF is now redundant, and replaced by matrix style decoder
// See AtkDecode built using Mix and ATKMatrix below
/*
ATKPantoF {
	*ar {arg numChannels, in, orientation = 1, directivity = 1, mul = 1, add = 1;
		(in.isKindOf(FOA)).if({
			^AtkPantoF.ar(numChannels, in.w, in.x, in.y, orientation, directivity);
		})
	}
}
*/
		
AtkMonoToB : Panner {
	
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

AtkSterToB : Panner {

	*ar { arg l, r, azimuth=0, width = 0.5pi, elevation=0;
		^this.multiNew('audio', l, r, azimuth, width, elevation)
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

// See AtkDecode built using Mix and ATKMatrix below
/*
AtkDecode : UGen {
	*ar {arg w, x, y, z, azimuth, elevation, mul = 1, add = 0;
		^this.multiNew('audio', w, x, y, z, azimuth, elevation).madd(mul, add);
		}
 	
 	checkInputs {
 		inputs[0..3].do({arg input, i; 	
 			if (rate !== input.rate) { 
 				^("input " + i + "is not" + rate + "rate: " + input + input.rate);
 			};
 			})
 		^this.checkValidInputs 
 		}
	}
*/
// ATKPantoF is now redundant, and replaced by matrix style decoder
// See AtkDecode built using Mix and ATKMatrix below
/*
AtkPantoF : Panner {
	*ar {arg numChans, w, x, y, orientation = 1, directivity = 1, mul = 1, add = 0;
		^this.multiNew('audio', numChans, w, x, y, orientation, directivity).madd(mul, add)
		}
	
	init { arg numChans ... theInputs;
		inputs = theInputs;
		channels = Array.fill(numChans, { arg i; OutputProxy(rate,this, i) });
		^channels
	}
 	checkInputs { ^this.checkNInputs(3) }
 	
	}
*/
	
AtkDirect : Atk {
	*ar { arg w, x, y, z, angle = pi/2, mul = 1, add = 0;
		^this.multiNew('audio', w, x, y, z, angle).madd(mul, add);
	}
}


AtkSquishX : Atk {
	*ar { arg w, x, y, z, angle = pi/2, mul = 1, add = 0;
		^this.multiNew('audio', w, x, y, z, angle).madd(mul, add);
	}
}

AtkSquishY : AtkSquishX { }
AtkSquishZ : AtkSquishX { }

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
	*ar { arg w, x, y, z, dom = 0;
		^this.multiNew('audio', w, x, y, z, dom);
	}	
}

AtkDominateY : AtkDominateX { }
AtkDominateZ : AtkDominateX { }


//------------------------------------------------------------------------
// Filters

AtkProximity : Atk { 
	*ar { arg w, x, y, z, distance = 0, mul = 1, add = 0;
		^this.multiNew('audio', w, x, y, z, distance).madd(mul, add);
	}

}

// RENAME to NFC
AtkDistance : Atk { 
	*ar { arg w, x, y, z, distance = 0, mul = 1, add = 0;
		^this.multiNew('audio', w, x, y, z, distance).madd(mul, add);
	}
		
}

// realised directly, rather than using RMShelf
AtkPsychoShelf { 
	*ar { arg w, x, y, z, frequency = 400, k = [(3/2).sqrt, 3.sqrt/2], mul = 1, add = 0;
		
		var k2;
		var wc, c;
		var a0, a1, a2, b1, b2;

		// expand k from degree (order) gains to channel gains
		k2 = k.collect({ arg item, i;
			Array.fill(2 * i + 1, {item})}).flat;

		// calculate coefficients
//		wc = (pi * frequency / SampleRate.ir(0)).tan;
//		wc = (pi * frequency / 44100).tan;
		wc = (pi * frequency / Server.default.sampleRate).tan; // I'm sure there's a better way!!
		c = (wc - 1) / (wc + 1);

		a0 = (((1 - k2)/4) * (1 + (c**2))) + (((1 + k2)/2) * c);
		a1 = ((1 - k2) * c) + (((1 + k2)/2) * (1 + (c**2)));
		a2 = a0;

		b1 = Array.fill( k2.size, { (2*c).neg } );
		b2 = Array.fill( k2.size, { (c**2).neg } );

		^SOS.ar([w, x, y, z], a0, a1, a2, b1, b2, mul, add);
	}
		
}



//------------------------------------------------------------------------
// Decoder built using Mix and ATKMatrix
//
// Likely, we'll want to integrate this much better with FOA, checking for
// valid inputs etc.

AtkDecode {
	*ar { arg in, decoderMatrix, mul = 1, add = 0;

		^Mix.fill( decoderMatrix.matrix.cols, { arg speaker;
			decoderMatrix.matrix.flop.asArray.at(speaker) * in.at(speaker)
		}).madd(mul, add);
	}
}
