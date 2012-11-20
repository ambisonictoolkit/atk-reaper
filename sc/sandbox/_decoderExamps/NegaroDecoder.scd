/*
encode:
    - a plane wave
    - a diffuse field

decode: by routing to a single decoder
*/

s.boot; // boot server

~ambiBus = Bus.audio(s, 4); // allocate four channels for routing

// define planewave synthDef
SynthDef.new(\foaPW, {arg outBus, duration = 0.05, theta, phi;
	var src, foa;

	src = WhiteNoise.ar(-12.dbamp) *
	EnvGen.kr(
		Env.new([0, 1, 0], [0.5, 0.5], \sin),
		timeScale: duration,
		doneAction: 2);

	// Encode into our foa signal
	foa = FoaPanB.ar(src, theta, phi);

	// send out
	Out.ar(outBus, foa);
}).add;

// define diffuse synthDef
SynthDef.new(\foaDiff, {arg outBus, duration = 3.0, theta = 0, phi = 0;
	var src, foa, angle;

	src = PinkNoise.ar(-12.dbamp * [1, 1, 1, 1]) * EnvGen.kr(
		Env.new([0, 1, 0], [0.5, 0.5], \sin),
		timeScale: duration,
		doneAction: 2);

	angle = Line.kr(pi/2, -pi/2,
		duration,
		doneAction: 2);

	// Encode into our foa signal
	foa = FoaTransform.ar(src, 'push', angle, theta, phi);

	// send out
	Out.ar(outBus, foa);
}).add;


// define decoder and decoder SynthDef, choose one:
~decoder = FoaDecoderMatrix.newStereo(131/2 * pi/180, 0.5); // simple stereo
~decoder = FoaDecoderKernel.newUHJ; // UHJ, a kernel decoder
~decoder = FoaDecoderKernel.newListen; // Listen HRTF, a kernel decoder

SynthDef.new(\foaDecode, {arg inBus;
	// var foa, decoderMatrix, out;
	var foa, out;

	foa = In.ar(inBus, 4);
	out = FoaDecode.ar(foa, ~decoder);

	// send out
	Out.ar(0, out);
}).add;

// start the decoder, reading from bus ~ambiBus, at the \tail
~decoderSynth = Synth.new(\foaDecode, [\inBus, ~ambiBus], 1, \addToTail);

// play some stuff
(
// diffuse synth
Synth.new(\foaDiff, [\outBus, ~ambiBus, \theta, pi/2, ]);

// planewave synth
Routine.run({
	20.do({
		Synth.new(\foaPW, [\outBus, ~ambiBus, \theta, pi.rand2, \phi, 0]);
		0.1.wait;
	})
});
)

// clean up
(
~decoderSynth.free; // free the decoder
~ambiBus.free; // free the audio bus
(~decoder.class == FoaDecoderKernel).if({~decoder.free});
)

s.quit;