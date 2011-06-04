/* FOA wrappers */

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

ATKPantoF {
	*ar {arg numChannels, in, orientation = 1, directivity = 1, mul = 1, add = 1;
		(in.isKindOf(FOA)).if({
			^AtkPantoF.ar(numChannels, in.w, in.x, in.y, orientation, directivity);
		})
	}
}
		
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



AtkProximity : Atk { 
	*ar { arg w, x, y, z, distance = 0, mul = 1, add = 0;
		^this.multiNew('audio', w, x, y, z, distance).madd(mul, add);
	}

}

AtkDistance : Atk { 
	*ar { arg w, x, y, z, distance = 0, mul = 1, add = 0;
		^this.multiNew('audio', w, x, y, z, distance).madd(mul, add);
	}
		
}
