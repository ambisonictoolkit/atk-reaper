// ATK SynthDef / Synth examples
Server.default = s = Server.local.boot;


/* ENCODERS */
(
/* omnidirection soundfield, push to a plane wave: */
SynthDef(\foaMonoEnc, {arg outBus, rho = 1, theta = 0, phi = 0;
	var angle, src, foa, encoderMatrix;
	src = PinkNoise.ar(-6.dbamp);
	angle = pi/2 * rho;

	// Encode into our foa signal
	// Note: as the encoder is a matrix encoder, it can be
	//			included in the SynthDef
	encoderMatrix = FoaEncoderMatrix.newOmni;
	foa = FoaEncode.ar(src, encoderMatrix); // encode as omni
	foa = FoaTransform.ar(foa, 'push', angle, theta, phi);   // image / pan
	Out.ar(outBus, foa);
}).add;

/* plane wave, with proximity effect: */
SynthDef(\foaMonoProx, {arg outBus, rho = 1, theta = 0, phi = 0;
	var src, foa, freq;
	freq = 20.0;					// HPF, to control proximity boost
	src = PinkNoise.ar(-6.dbamp);
	src = HPF.ar(src, freq);		// precondition signal for proximity

	// Encode into our foa signal
	foa = FoaPanB.ar(src, theta, phi); // encode as plane wave
	foa = FoaTransform.ar(foa, 'proximity', rho);   // image
	Out.ar(outBus, foa);
}).add;

/* decorrelated soundfield, push to a plane wave: */
SynthDef(\foaDecorEnc, {arg outBus, rho = 1, theta = 0, phi = 0;
	var angle, src, foa, encoderMatrix;
	src = PinkNoise.ar(-9.dbamp * [1, 1, 1, 1]);
	angle = pi/2 * rho;

	// Encode into our foa signal
	// Note: as the encoder is a matrix encoder, it can be
	//			included in the SynthDef
	encoderMatrix = FoaEncoderMatrix.newAtoB;
	foa = FoaEncode.ar(src, encoderMatrix); // encode as AtoB (decorrelated)
	foa = FoaTransform.ar(foa, 'push', angle, theta, phi);   // image / pan
	Out.ar(outBus, foa);
}).add;

/* decorrelated soundfield, push to a plane wave: */
SynthDef(\foaDecProxEnc, {arg outBus, rho = 1, theta = 0, phi = 0;
	var angle, src, foa, encoderMatrix, freq;
	freq = 20.0;					// HPF, to control proximity boost
	src = PinkNoise.ar(-9.dbamp * [1, 1, 1, 1]);
	src = HPF.ar(src, freq);		// precondition signal for proximity
	angle = pi/2 * rho;

	// Encode into our foa signal
	// Note: as the encoder is a matrix encoder, it can be
	//			included in the SynthDef
	encoderMatrix = FoaEncoderMatrix.newAtoB;
	foa = FoaEncode.ar(src, encoderMatrix); // encode as AtoB (decorrelated)
	foa = FoaTransform.ar(foa, 'push', angle, theta, phi);   // image / pan
	foa = FoaTransform.ar(foa, 'proximity', rho);   // image
	Out.ar(outBus, foa);
}).add;
)

/* DECODERS */
/* simple stereo matrix decoder: */
~decoder = FoaDecoderMatrix.newStereo(131/2 * pi/180, 0.5);

//OR

/* HRTF (kernel) decoder: */
~decoder = FoaDecoderKernel.newListen(1013);

//OR

/* dual band narrow quad decoder: */
~decoder = FoaDecoderMatrix.newQuad(pi/6, 'dual');

SynthDef(\foaDecode, {arg inBus;
        var foa, out;
        foa = In.ar(inBus, 4);
        out = FoaDecode.ar(foa, ~decoder);
        Out.ar(0, out);
}).add;

// use this for quad at Audio Room 4
SynthDef(\foaDecode, {arg inBus;
	var foa, fl, fr, bl, br;
	foa = In.ar(inBus, 4);
	foa = FoaNFC.ar(foa, 1.27); // nfc for r = 1.27 m
	#fl, bl, br, fr = FoaDecode.ar(foa, ~decoder);
	Out.ar(0, [fl, fr, bl, br]);
}).add;


// allocate four channels for routing
a = Bus.audio(s, 4);

// start the decoder, reading bus 'a' at the \tail
b = Synth(\foaDecode, [\inBus, a], 1, \addToTail);



// start encoder (listed in order of sophistication)
c = Synth(\foaMonoEnc, [\outBus, a]);	// simple 'in head' localisation
// OR
c = Synth(\foaMonoProx, [\outBus, a]);	// proximity effect
// OR
c = Synth(\foaDecorEnc, [\outBus, a]);	// diffuse field
// OR
c = Synth(\foaDecProxEnc, [\outBus, a]);	// diffuse field with proximity effect


// GUI with unwrapping
// NOTE: theta is unwrapped here
//		FoaPush.kr to be updated to do this internally
(
	var w, xyslider, zslider, azislider;
	var cart;
	var thetau, delta;

	w = Window("Ambisonics Controller", Rect(200, 200, 640, 480)).front;
	w.addFlowLayout;

	cart = Cartesian.new;
	thetau = 0.0;				// init wrapped theta val

	xyslider = Slider2D(w, Rect(0, 0, 400, 400));
	zslider = Slider(w, Rect(0, 0, 50, 400));

	// init sliders to front centre
	xyslider.setXY(0.5, 1.0);
	zslider.value = 0.5;

	xyslider.action_({arg hori;
		cart.x_(hori.y.linlin(0,1,-1,1)).y_(hori.x.linlin(0,1,1,-1));

		// unwrap theta
		// (TBA: FoaPush.kr to be updated to do this internally)
		delta = cart.asSpherical.theta - thetau;
		delta = delta.wrap(-pi, pi);
		case
			{ delta > 0 } { thetau = thetau + (delta % 2pi) }
			{ delta < 0 } { thetau = thetau - (-1 * delta % 2pi) };

		c.set(\theta, thetau);
		c.set(\rho, cart.asSpherical.rho.clip(0.05, 1));	// clip range
	});
	zslider.action_{ arg elev;
		cart.z_(elev.value.linlin(0,1,-1,1));

		c.set(\phi, cart.asSpherical.phi);
		c.set(\rho, cart.asSpherical.rho.clip(0.05, 1));	// clip range
	}
)





// GUI without wrapping
// clicks across the rear!! (set rho = 1; e.g., pan across the bottom of the screen)
// no clicks!!! now fixed!!
(
	var w, xyslider, zslider, azislider;
	var cart;

	w = Window("Ambisonics Controller", Rect(200, 200, 640, 480)).front;
	w.addFlowLayout;

	cart = Cartesian.new;

	xyslider = Slider2D(w, Rect(0, 0, 400, 400));
	zslider = Slider(w, Rect(0, 0, 50, 400));

	// init sliders to front centre
	xyslider.setXY(0.5, 1.0);
	zslider.value = 0.5;

	xyslider.action_({arg hori;
		cart.x_(hori.y.linlin(0,1,-1,1)).y_(hori.x.linlin(0,1,1,-1));

		c.set(\theta, cart.asSpherical.theta);
		c.set(\rho, cart.asSpherical.rho.clip(0.05, 1));	// clip range
	});
	zslider.action_{ arg elev;
		cart.z_(elev.value.linlin(0,1,-1,1));

		c.set(\phi, cart.asSpherical.phi);
		c.set(\rho, cart.asSpherical.rho.clip(0.05, 1));	// clip range
	}
)




c.free; // free the encoder
b.free; // free the decoder
a.free; // free the audio bus
~kernelDecoder.free; // free kernel decoder, if used