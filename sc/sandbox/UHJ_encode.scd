// -------------------------------------------------------
// ATK (Sandbox)
//
// Coded by Joseph Anderson 2011
//
// -------------------------------------------------------

//------------------------------------------------------------------------
// Stereo encoders
//
//   uhj_to_b            "Ambisonic Decoders for HDTV" (1992)
//
//
//	NOTE: UHJ kernels generated externally via Muse/ATK

// NOTE: to run this example, the UHJ decoder (and kernels) found in UHJ.scd must
//		be loaded


//   uhj_to_b
(
~atkUHJtoB = { arg sig, encoderKernels;
	
	var res;

	// encode here!
	// convolve and mix to W, X, Y
	res = encoderKernels.shape.at(0).collect({ arg i;	// channel [L, R]
		encoderKernels.shape.at(1).collect({ arg j;	// harmonic [W, X, Y]
			Convolution2.ar(
			sig.at(i),
			encoderKernels.at(i).at(j),
			framesize: encoderKernels.at(i).at(j).numFrames)
		})
	});
	Mix.ar(res) ++ Silent.ar(1);					// mix and append Z (zeros)
}
)


// UHJ kernel reading
(
var kernelFiles;
var sourceDir, kernelDir;
var kernelSR;
var kernelSize;
var subjectID;
var kernelNames;
var kernelExt;
var chans;


// params
sourceDir		= "/Volumes/Audio/test";      			// corresponds to 'target_dir'
kernelDir		= "/ATK_kernels/FOA/encoders/UHJ";		// corresponds to 'file_dir'

kernelNames	= [ "UHJ_L", "UHJ_R" ];					// kernel names
kernelExt		= ".wav";								// kernel file extension

kernelSR 		= s.sampleRate.asString;				// shelf filtered UHJ encoding
												// presume decoding will correctly												// employ psychoacoustic shelf filter

//kernelSize	= 512;
kernelSize	= 1024;
//kernelSize	= 2048;
//kernelSize	= 4096;
//kernelSize	= 8192;
subjectID		= 0;									// only one 'subject' for UHJ

chans 		= 3;									// kernels are 3 channel, [ W, X, Y ]


// collect kernel file names
kernelFiles = kernelNames.collect({ arg name;
	sourceDir ++ kernelDir ++ "/SR_" ++ kernelSR.padLeft(6, "0") ++ "/N_" ++ 
				kernelSize.asString.padLeft(4, "0") ++ "/" ++
				subjectID.asString.padLeft(4, "0") ++ "/" ++ name ++ kernelExt
	});


// open kernels -- load into buffers
// load each W, X, Y buffer separately
c = kernelFiles.collect({ arg kernel;
	chans.collect({ arg chan;
		Buffer.readChannel(s, kernel, channels: [chan])
		})
	})
)

// encoding / decoding test... b-format, UHJ, b-format, quad
(
{
var fl, bl, br, fr;					// for quad output
var sig;
var bSig, bSig2;
var azim, angle;
var uhj;


// constants
angle = (30 / 180) * pi;


// kr....
azim = MouseX.kr(pi, -pi);

// test signal
bSig = AtkMonoToB.ar(PinkNoise.ar, azim);

// decode to UHJ
uhj = ~atkBtoUHJ.value(bSig, b);

// (re-)encode to b-format
bSig2 = ~atkUHJtoB.value(uhj, c);


// decode to output (single band decoder, quad)
// NOTE: as the UHJ encoder is shelf filtered, it is best
//		to decode via a dual band decoder.
#fl, bl, br, fr = ~atkQuadDBD.value(bSig2, angle);
//#fl, bl, br, fr = ~atkQuadDBD.value(bSig, angle);	// compare with original

[fl, fr, bl, br]		// reorder output to match speaker arrangement


}.scope;
)


// free the buffers!!
(
b.shape.at(0).do({ arg i;
	b.shape.at(1).do({ arg j;
		b.at(i).at(j).free
		})
	})

c.shape.at(0).do({ arg i;
	c.shape.at(1).do({ arg j;
		c.at(i).at(j).free
		})
	})
)

