// ------------------------------------------------------------
// switch to local server, and boot
(
        s = Server.local;
        Server.default = s;
        s.boot;
)
s.quit



// ------------------------------------------------------------
// omni encoder
// mono pink noise source


// define encoder matrix
~encoder = FoaEncoderMatrix.newOmni

// inspect
~encoder.kind
~encoder.numChannels
~encoder.dirChannels



// ------------------------------------------------------------
// spreader encoder
// mono pink noise source


// define encoder matrix
~encoder = FoaEncoderKernel.newDiffuse
~encoder = FoaEncoderKernel.newDiffuse(0000, 2048)
~encoder = FoaEncoderKernel.newDiffuse(0001, 2048)
~encoder = FoaEncoderKernel.newDiffuse(0002, 2048)
~encoder = FoaEncoderKernel.newDiffuse(0003, 2048)
~encoder = FoaEncoderKernel.newDiffuse(0004, 2048)
~encoder = FoaEncoderKernel.newDiffuse(0005, 2048)

~encoder = FoaEncoderKernel.newDiffuse(0005, 8192)

~encoder = FoaEncoderKernel.newSpread(0006, 2048)


// free kernel
~encoder.free

// inspect
~encoder
~encoder.kind
~encoder.numChannels
~encoder.dirChannels * 180/pi



(
{
    var sig;                            // audio signal
    var angle, azim;                    // angle and azimuth control


    // display encoder and decoder
    "Ambisonic encoding via % encoder".format(~encoder.kind).postln;
    "Ambisonic decoding via % decoder".format(~decoder.kind).postln;

    // angle ---> top           = push to plane wave
    //            bottom        = omni-directional
    angle = MouseY.kr(pi/2, 0);

    // azimuth -> hard left     = back
    //            centre        = centre
    //            hard right    = back
    azim = MouseX.kr(pi, -pi);

    // ------------------------------------------------------------
    // test sig
	// sig = PinkNoise.ar;                         // mono pink noise
	// sig = Impulse.ar(2.0);                         // impulse
	sig = Saw.ar(110.0, Decay.ar(Impulse.ar(2.0), 1.0/2.0));     // gated saw


    // ------------------------------------------------------------
    // encode
    sig = FoaEncode.ar(sig, ~encoder);

	// // test for A-format swapping
	// sig = sig * [sqrt(3)/sqrt(2), 1, 1, 1];
	// sig = FoaEncode.ar(sig, FoaEncoderMatrix.newAtoB);


    // ------------------------------------------------------------
    // transform
	sig = FoaTransform.ar(sig, 'push', angle, azim);
	// sig = FoaTransform.ar(sig, 'focus', angle, azim);



    // ------------------------------------------------------------
    // decode (via ~renderDecode)
    ~renderDecode.value(sig, ~decoder)

}.scope;
)

// free kernel
~encoder.free
// ------------------------------------------------------------




// ------------------------------------------------------------
// spreader encoder
// stereo soundfile read from disk


// define encoder matrix
~encoder = FoaEncoderKernel.newDiffuse
~encoder = FoaEncoderKernel.newDiffuse(0000, 2048)
~encoder = FoaEncoderKernel.newDiffuse(0001, 2048)
~encoder = FoaEncoderKernel.newDiffuse(0002, 2048)
~encoder = FoaEncoderKernel.newDiffuse(0003, 2048)
~encoder = FoaEncoderKernel.newDiffuse(0004, 2048)
~encoder = FoaEncoderKernel.newDiffuse(0005, 2048)

~encoder = FoaEncoderKernel.newDiffuse(0005, 8192)

// free kernel
~encoder.free

// inspect
~encoder
~encoder.kind
~encoder.numChannels
~encoder.dirChannels * 180/pi


// read a whole sound into memory
// remember to free the buffer later!
// (boot the server, if you haven't!)
~sndbuf = Buffer.read(s, Atk.userSoundsDir ++ "/stereo/The_City_Waites-The_Downfall.wav")
~sndbuf = Buffer.read(s, Atk.userSoundsDir ++ "/stereo/The_City_Waites-An_Old.wav")
~sndbuf = Buffer.read(s, Atk.userSoundsDir ++ "/stereo/Aurora_Surgit-Lux_Aeterna.wav")
~sndbuf = Buffer.read(s, Atk.userSoundsDir ++ "/stereo/Aurora_Surgit-Dies_Irae.wav")

// free buffer
~sndbuf.free


(
{
    var sig;                            // audio signal
    var angle, azim;                    // angle and azimuth control


    // display encoder and decoder
    "Ambisonic encoding via % encoder".format(~encoder.kind).postln;
    "Ambisonic decoding via % decoder".format(~decoder.kind).postln;

    // angle ---> top           = push to plane wave
    //            bottom        = omni-directional
    angle = MouseY.kr(pi/2, 0);

    // azimuth -> hard left     = back
    //            centre        = centre
    //            hard right    = back
    azim = MouseX.kr(pi, -pi);

    // ------------------------------------------------------------
    // test sig
    sig = PlayBuf.ar(~sndbuf.numChannels, ~sndbuf, doneAction:2);   // soundfile
	sig = 0.5 * sig.sum; // to mono


    // ------------------------------------------------------------
    // encode
    sig = FoaEncode.ar(sig, ~encoder);

    // ------------------------------------------------------------
    // transform
    sig = FoaTransform.ar(sig, 'push', angle, azim);



    // ------------------------------------------------------------
    // decode (via ~renderDecode)
    ~renderDecode.value(sig, ~decoder)

}.scope;
)

// free buffer
~sndbuf.free
// free kernel
~encoder.free
// ------------------------------------------------------------


// ------------------------------------------------------------
// push and rtt transformer
// spreader encoder
// stereo soundfile read from disk


// define encoder matrix
~encoder = FoaEncoderKernel.newDiffuse
~encoder = FoaEncoderKernel.newDiffuse(0000, 2048)
~encoder = FoaEncoderKernel.newDiffuse(0001, 2048)
~encoder = FoaEncoderKernel.newDiffuse(0002, 2048)
~encoder = FoaEncoderKernel.newDiffuse(0003, 2048)
~encoder = FoaEncoderKernel.newDiffuse(0004, 2048)
~encoder = FoaEncoderKernel.newDiffuse(0005, 2048)

~encoder = FoaEncoderKernel.newDiffuse(0003, 4096)

~encoder = FoaEncoderKernel.newDiffuse(0003, 8192)
~encoder = FoaEncoderKernel.newDiffuse(0005, 8192)

// free kernel
~encoder.free

// inspect
~encoder
~encoder.kind
~encoder.numChannels
~encoder.dirChannels * 180/pi


// read a whole sound into memory
// remember to free the buffer later!
// (boot the server, if you haven't!)
~sndbuf = Buffer.read(s, Atk.userSoundsDir ++ "/stereo/The_City_Waites-The_Downfall.wav")
~sndbuf = Buffer.read(s, Atk.userSoundsDir ++ "/stereo/The_City_Waites-An_Old.wav")
~sndbuf = Buffer.read(s, Atk.userSoundsDir ++ "/stereo/Aurora_Surgit-Lux_Aeterna.wav")
~sndbuf = Buffer.read(s, Atk.userSoundsDir ++ "/stereo/Aurora_Surgit-Dies_Irae.wav")

// free buffer
~sndbuf.free


(
{
    var sig;                            // audio signal
    var angle, azim;                    // angle and azimuth control


    // display encoder and decoder
    "Ambisonic encoding via % encoder".format(~encoder.kind).postln;
    "Ambisonic decoding via % decoder".format(~decoder.kind).postln;

    // angle ---> top           = push to plane wave
    //            bottom        = omni-directional
    angle = MouseY.kr(pi/2, 0);

    // azimuth -> hard left     = back
    //            centre        = centre
    //            hard right    = back
    azim = MouseX.kr(pi, -pi);

    // ------------------------------------------------------------
    // test sig
    sig = PlayBuf.ar(~sndbuf.numChannels, ~sndbuf, doneAction:2);   // soundfile
	sig = 0.5 * sig.sum; // to mono


    // ------------------------------------------------------------
    // encode
    sig = FoaEncode.ar(sig, ~encoder);

    // ------------------------------------------------------------
    // transform
	sig = FoaRTT.ar(sig,
		LFNoise2.kr(1.0/5.0, pi),
		LFNoise2.kr(1.0/5.0, pi),
		LFNoise2.kr(1.0/5.0, pi)
	);
    sig = FoaTransform.ar(sig, 'push', angle, azim);



    // ------------------------------------------------------------
    // decode (via ~renderDecode)
    ~renderDecode.value(sig, ~decoder)

}.scope;
)

// free buffer
~sndbuf.free
// free kernel
~encoder.free

// quit
s.quit

// ------------------------------------------------------------