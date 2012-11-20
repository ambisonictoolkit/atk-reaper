Server.default = s = Server.local.boot;


// test cases to illustrate clicking w/ supplying wrapped angles to Atk UGens
//
// when supplied with wrapped angle values generated at .kr, a number of the transform
// UGens click. this isn't the case when wrapped .ar values are supplied.
//
// one possible solution would be to make sure all of the resulting matrix coeffs
// are smoothed w/in the UGens, to avoid discontinuities
//
// another possibility could be to unwrap supplied angles, so that the resulting
// matrix coeffs don't have large discontinuities. the examples below include
// unwrapped .kr vals which don't click… see code below


// ---------------------------------------------------------
// UGen: FoaPanB

// case: .ar, wrapped angle
// outcome: no clicks
(
{
	var freq = 800, gain = -6.0;
	var rotRate = 1.0;
	var angle;

	angle = LFSaw.ar(rotRate, mul: pi);

	FoaDecode.ar(
		FoaPanB.ar(
			SinOsc.ar(freq, mul: gain.dbamp),
			angle
		),
		FoaDecoderMatrix.newStereo
	)
}.scope;
)

// case: .kr, wrapped angle
// outcome: no clicks
(
{
	var freq = 800, gain = -6.0;
	var rotRate = 1.0;
	var angle;

	angle = LFSaw.kr(rotRate, mul: pi);

	FoaDecode.ar(
		FoaPanB.ar(
			SinOsc.ar(freq, mul: gain.dbamp),
			angle
		),
		FoaDecoderMatrix.newStereo
	)
}.scope;
)

// ---------------------------------------------------------
// UGen: FoaRotate

// case: .ar, wrapped angle
// outcome: no clicks
(
{
	var freq = 800, gain = -6.0;
	var rotRate = 1.0;
	var angle;

	angle = LFSaw.ar(rotRate, mul: pi);

	FoaDecode.ar(
		FoaRotate.ar(
			FoaEncode.ar(
				SinOsc.ar(freq, mul: gain.dbamp),
				FoaEncoderMatrix.newDirection	// init position [0, 0]
			),
			angle
		),
		FoaDecoderMatrix.newStereo
	)
}.scope;
)

// case: .kr, wrapped angle
// outcome: no clicks w/ fix
(
{
	var freq = 800, gain = -6.0;
	var rotRate = 1.0;
	var angle;

	angle = LFSaw.kr(rotRate, mul: pi);

	FoaDecode.ar(
		FoaRotate.ar(
			FoaEncode.ar(
				SinOsc.ar(freq, mul: gain.dbamp),
				FoaEncoderMatrix.newDirection	// init position [0, 0]
			),
			angle
		),
		FoaDecoderMatrix.newStereo
	)
}.scope;
)

// case: .kr, unwrapped angle
// outcome: no clicks previously, now this clicks as an overflow!
(
{
	var freq = 800, gain = -6.0;
	var rotRate = 1.0;
	var angle;

	angle = Integrator.kr(DC.kr(pi * rotRate/(s.sampleRate/2))) * s.options.blockSize;

	FoaDecode.ar(
		FoaRotate.ar(
			FoaEncode.ar(
				SinOsc.ar(freq, mul: gain.dbamp),
				FoaEncoderMatrix.newDirection	// init position [0, 0]
			),
			angle
		),
		FoaDecoderMatrix.newStereo
	)
}.scope;
)

// ---------------------------------------------------------
// UGen: FoaRTT

// case: .ar, wrapped angle
// outcome: no clicks
(
{
	var freq = 800, gain = -6.0;
	var rotRate = 1.0;
	var angle;

	angle = LFSaw.ar(rotRate, mul: pi);

	FoaDecode.ar(
		FoaRTT.ar(
			FoaEncode.ar(
				SinOsc.ar(freq, mul: gain.dbamp),
				FoaEncoderMatrix.newDirection	// init position [0, 0]
			),
			angle
		),
		FoaDecoderMatrix.newStereo
	)
}.scope;
)

// case: .kr, wrapped angle
// outcome: no clicks
(
{
	var freq = 800, gain = -6.0;
	var rotRate = 1.0;
	var angle;

	angle = LFSaw.kr(rotRate, mul: pi);

	FoaDecode.ar(
		FoaRTT.ar(
			FoaEncode.ar(
				SinOsc.ar(freq, mul: gain.dbamp),
				FoaEncoderMatrix.newDirection	// init position [0, 0]
			),
			angle
		),
		FoaDecoderMatrix.newStereo
	)
}.scope;
)

// case: .kr, unwrapped angle
// outcome: no clicks previously, now this clicks as an overflow!
(
{
	var freq = 800, gain = -6.0;
	var rotRate = 1.0;
	var angle;

	angle = Integrator.kr(DC.kr(pi * rotRate/(s.sampleRate/2))) * s.options.blockSize;

	FoaDecode.ar(
		FoaRTT.ar(
			FoaEncode.ar(
				SinOsc.ar(freq, mul: gain.dbamp),
				FoaEncoderMatrix.newDirection	// init position [0, 0]
			),
			angle
		),
		FoaDecoderMatrix.newStereo
	)
}.scope;
)

// ---------------------------------------------------------
// UGen: FoaPush

// case: .ar, wrapped angle
// outcome: no clicks
(
{
	var freq = 800, gain = -6.0;
	var rotRate = 1.0;
	var angle;

	angle = LFSaw.ar(rotRate, mul: pi);

	FoaDecode.ar(
		FoaPush.ar(
			FoaEncode.ar(
				SinOsc.ar(freq, mul: gain.dbamp),
				FoaEncoderMatrix.newOmni
			),
			pi/2,
			angle
		),
		FoaDecoderMatrix.newStereo
	)
}.scope;
)

// case: .kr, wrapped angle
// outcome: no clicks w/ fix
(
{
	var freq = 800, gain = -6.0;
	var rotRate = 1.0;
	var angle;

	angle = LFSaw.kr(rotRate, mul: pi);

	FoaDecode.ar(
		FoaPush.ar(
			FoaEncode.ar(
				SinOsc.ar(freq, mul: gain.dbamp),
				FoaEncoderMatrix.newOmni
			),
			pi/2,
			angle
		),
		FoaDecoderMatrix.newStereo
	)
}.scope;
)

// case: .kr, unwrapped angle
// outcome: no clicks previously, now this clicks as an overflow!
(
{
	var freq = 800, gain = -6.0;
	var rotRate = 1.0;
	var angle;

	angle = Integrator.kr(DC.kr(pi * rotRate/(s.sampleRate/2))) * s.options.blockSize;

	FoaDecode.ar(
		FoaPush.ar(
			FoaEncode.ar(
				SinOsc.ar(freq, mul: gain.dbamp),
				FoaEncoderMatrix.newOmni
			),
			pi/2,
			angle
		),
		FoaDecoderMatrix.newStereo
	)
}.scope;
)

// ---------------------------------------------------------
// UGen: FoaFocus

// case: .ar, wrapped angle
// outcome: no clicks
(
{
	var freq = 800, gain = 0.0;
	var rotRate = 1.0;
	var angle;

	angle = LFSaw.ar(rotRate, mul: pi);

	FoaDecode.ar(
		FoaFocus.ar(
			FoaEncode.ar(
				SinOsc.ar(freq, mul: gain.dbamp),
				FoaEncoderMatrix.newOmni
			),
			pi/2,
			angle
		),
		FoaDecoderMatrix.newStereo
	)
}.scope;
)

// case: .kr, wrapped angle
// outcome: no clicks with fix
(
{
	var freq = 800, gain = 0.0;
	var rotRate = 1.0;
	var angle;

	angle = LFSaw.kr(rotRate, mul: pi);

	FoaDecode.ar(
		FoaFocus.ar(
			FoaEncode.ar(
				SinOsc.ar(freq, mul: gain.dbamp),
				FoaEncoderMatrix.newOmni
			),
			pi/2,
			angle
		),
		FoaDecoderMatrix.newStereo
	)
}.scope;
)

// case: .kr, unwrapped angle
// outcome: no clicks previously, now clicks as overflow
(
{
	var freq = 800, gain = 0.0;
	var rotRate = 1.0;
	var angle;

	angle = Integrator.kr(DC.kr(pi * rotRate/(s.sampleRate/2))) * s.options.blockSize;

	FoaDecode.ar(
		FoaFocus.ar(
			FoaEncode.ar(
				SinOsc.ar(freq, mul: gain.dbamp),
				FoaEncoderMatrix.newOmni
			),
			pi/2,
			angle
		),
		FoaDecoderMatrix.newStereo
	)
}.scope;
)

// ---------------------------------------------------------
// UGen: FoaMirror

// case: .ar, unwrapped angle
// outcome: no clicks
(
{
	var freq = 800, gain = -6.0;
	var rotRate = 1.0;
	var angle;

	angle = LFSaw.ar(rotRate, mul: pi);

	FoaDecode.ar(
		FoaMirror.ar(
			FoaEncode.ar(
				SinOsc.ar(freq, mul: gain.dbamp),
				FoaEncoderMatrix.newDirection	// init position [0, 0]
			),
			angle
		),
		FoaDecoderMatrix.newStereo
	)
}.scope;
)

// case: .kr, wrapped angle
// outcome: no clicks with fix
(
{
	var freq = 800, gain = -6.0;
	var rotRate = 1.0;
	var angle;

	angle = LFSaw.kr(rotRate, mul: pi);

	FoaDecode.ar(
		FoaMirror.ar(
			FoaEncode.ar(
				SinOsc.ar(freq, mul: gain.dbamp),
				FoaEncoderMatrix.newDirection	// init position [0, 0]
			),
			angle
		),
		FoaDecoderMatrix.newStereo
	)
}.scope;
)

// case: .kr, unwrapped angle
// outcome: no clicks previously, now clicks as overflow
(
{
	var freq = 800, gain = -6.0;
	var rotRate = 1.0;
	var angle;

	angle = Integrator.kr(DC.kr(pi * rotRate/(s.sampleRate/2))) * s.options.blockSize;

	FoaDecode.ar(
		FoaMirror.ar(
			FoaEncode.ar(
				SinOsc.ar(freq, mul: gain.dbamp),
				FoaEncoderMatrix.newDirection	// init position [0, 0]
			),
			angle
		),
		FoaDecoderMatrix.newStereo
	)
}.scope;
)