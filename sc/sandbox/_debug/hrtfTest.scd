// ------------------------------------------------------------
// switch to local server, and boot
(
        s = Server.local;
        Server.default = s;
        s.boot;
)
s.quit

// ------------------------------------------------------------
// choose a decoder

// stereophonic / binaural
~decoder = FoaDecoderMatrix.newStereo(131/2 * pi/180, 0.5) // Cardioids at 131 deg
~decoder = FoaDecoderMatrix.newStereo(5/9 * pi, 0.5) // Cardioids at Duda angle
~decoder = FoaDecoderKernel.newUHJ                         // UHJ (kernel)

~decoder = FoaDecoderKernel.newSpherical                   // synthetic binaural

~decoder = FoaDecoderKernel.newCIPIC                       // CIPIC binaural
~decoder = FoaDecoderKernel.newCIPIC(0003)
~decoder = FoaDecoderKernel.newCIPIC(0012)
~decoder = FoaDecoderKernel.newCIPIC(0033)
~decoder = FoaDecoderKernel.newCIPIC(0040)

~decoder = FoaDecoderKernel.newListen                      // Listen binaural
~decoder = FoaDecoderKernel.newListen(1002)
~decoder = FoaDecoderKernel.newListen(1003)
~decoder = FoaDecoderKernel.newListen(1004)
~decoder = FoaDecoderKernel.newListen(1005)



// inspect
~decoder.kind
~decoder.kernelSize

// free kernel
~decoder.free


// ------------------------------------------------------------
// define ~renderDecode
(
~renderDecode = { arg in, decoder;
	// decode to stereo (or binaural)
	FoaDecode.ar(in, decoder)
}
)

// ------------------------------------------------------------
// now we're ready to try the examples below!
// ------------------------------------------------------------

// ------------------------------------------------------------
// B-format examples, produced via the ATK
// B-format soundfile read from disk


// read a whole sound into memory
// remember to free the buffer later!
// (boot the server, if you haven't!)
~sndbuf = Buffer.read(s, Atk.userSoundsDir ++ "/b-format/Anderson-Pacific_Slope.wav")
~sndbuf = Buffer.read(s, Atk.userSoundsDir ++ "/b-format/Howle-Calling_Tunes.wav")
~sndbuf = Buffer.read(s, Atk.userSoundsDir ++ "/b-format/Pampin-On_Space.wav")

~sndbuf = Buffer.read(s, Atk.userSoundsDir ++ "/b-format/Hodges-Purcell.wav")
~sndbuf = Buffer.read(s, Atk.userSoundsDir ++ "/b-format/Leonard-Orfeo_Trio.wav")
~sndbuf = Buffer.read(s, Atk.userSoundsDir ++ "/b-format/Courville-Dialogue.wav")
~sndbuf = Buffer.read(s, Atk.userSoundsDir ++ "/b-format/Leonard-Chinook.wav")
~sndbuf = Buffer.read(s, Atk.userSoundsDir ++ "/b-format/Leonard-Fireworks.wav")
~sndbuf = Buffer.read(s, Atk.userSoundsDir ++ "/b-format/Anderson-Nearfield.wav")



(
{
    var sig;                            // audio signal


    // display encoder and decoder
    "Ambisonic decoding via % decoder".format(~decoder.kind).postln;

    // ------------------------------------------------------------
    // test sig
	sig = PlayBuf.ar(~sndbuf.numChannels, ~sndbuf, doneAction:2);   // soundfile
	// sig = sig * 3.dbamp;

    // ------------------------------------------------------------
    // decode (via ~renderDecode)
    ~renderDecode.value(sig, ~decoder)

}.scope;
)

// free buffer
~sndbuf.free
// ------------------------------------------------------------