// -------------------------------------------------------
// ATK (Sandbox)
//
// Coded by Joseph Anderson 2011
//
// -------------------------------------------------------

//------------------------------------------------------------------------
// Stereo encoders
//
//   superstereo
//
//
//	NOTE: superstereo kernels generated externally via Muse/ATK


//   superstereo
(
~atkSuperStereo = { arg sig, encoderKernels;
	
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


// superstereo kernel reading
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
kernelDir		= "/ATK_kernels/FOA/encoders/Super";		// corresponds to 'file_dir'

kernelNames	= [ "Super_L", "Super_R" ];				// kernel names
kernelExt		= ".wav";								// kernel file extension

kernelSR		= "None";								// for superstereo encoding, kernelSR = None
//kernelSR 		= s.sampleRate.asString;

//kernelSize	= 512;
kernelSize	= 1024;
//kernelSize	= 2048;
//kernelSize	= 4096;
//kernelSize	= 8192;
subjectID		= 0;									// only one 'subject' for superstereo

chans 		= 3;									// kernels are 3 channel, [ W, X, Y ]


// collect kernel file names
kernelFiles = kernelNames.collect({ arg name;
	sourceDir ++ kernelDir ++ "/SR_" ++ kernelSR.padLeft(6, "0") ++ "/N_" ++ 
				kernelSize.asString.padLeft(4, "0") ++ "/" ++
				subjectID.asString.padLeft(4, "0") ++ "/" ++ name ++ kernelExt
	});


// open kernels -- load into buffers
// load each W, X, Y buffer separately
b = kernelFiles.collect({ arg kernel;
	chans.collect({ arg chan;
		Buffer.readChannel(s, kernel, channels: [chan])
		})
	})
)

// encoding / decoding test... (stereo --> superstereo --> quad)
(
{
var fl, bl, br, fr;					// for quad output
var sig;
var bSig;
var azim;
var angle;


// constants
angle = (30 / 180) * pi;


// kr....
azim = MouseX.kr(pi/4, -pi/4);			

// test signal (Blumlein stereo)
sig = PinkNoise.ar * [(azim + (pi/4)).sin, (azim + (pi/4)).cos];

// encode to b-format
bSig = ~atkSuperStereo.value(sig, b);
//bSig = ~atkSimpleStereo.value(sig);


// decode to output (single band decoder, quad)
#fl, bl, br, fr = ~atkQuadSBD.value(bSig, angle);

[fl, fr, bl, br]		// reorder output to match speaker arrangement


//sig
//bSig

}.scope;
)


// free the buffers!!
(
b.shape.at(0).do({ arg i;
	b.shape.at(1).do({ arg j;
		b.at(i).at(j).free
		})
	})
)


////------------------------------------------------------------------------
//// Stereo encoders
////
////   b_to_uhj            "Ambisonic Decoders for HDTV" (1992)
////
////
////	NOTE: UHJ kernels generated externally via Muse/ATK
//
//
////   b_to_uhj
//(
//~atkBtoUHJ = { arg bSig, decoderKernels;
//	
//	var res;
//
//	// decode here!
//	// convolve and mix to L/R
//	res = decoderKernels.shape.at(0).collect({ arg i;	// harmonic [W, X, Y]
//		decoderKernels.shape.at(1).collect({ arg j;	// channel [L, R]
//			Convolution2.ar(
//			bSig.at(i),
//			decoderKernels.at(i).at(j),
//			framesize: decoderKernels.at(i).at(j).numFrames)
//		})
//	});
//	Mix.ar(res);								// mix
//}
//)
//
//
//// UHJ kernel reading
//(
//var kernelFiles;
//var sourceDir, kernelDir;
//var kernelSR;
//var kernelSize;
//var subjectID;
//var kernelNames;
//var kernelExt;
//var chans;
//
//
//// params
//sourceDir		= "/Volumes/Audio/test";      			// corresponds to 'target_dir'
//kernelDir		= "/ATK_kernels/FOA/decoders/UHJ";		// corresponds to 'file_dir'
//
//kernelNames	= [ "UHJ_W", "UHJ_X", "UHJ_Y" ];			// kernel names
//kernelExt		= ".wav";								// kernel file extension
//
//kernelSR		= "None";								// for UHJ decoding, kernelSR = None
////kernelSR 		= s.sampleRate.asString;
//
////kernelSize	= 512;
//kernelSize	= 1024;
////kernelSize	= 2048;
////kernelSize	= 4096;
////kernelSize	= 8192;
//subjectID		= 0;									// only one 'subject' for UHJ
//
//chans 		= 2;									// kernels are 2 channel, L/R
//
//
//// collect kernel file names
//kernelFiles = kernelNames.collect({ arg name;
//	sourceDir ++ kernelDir ++ "/SR_" ++ kernelSR.padLeft(6, "0") ++ "/N_" ++ 
//				kernelSize.asString.padLeft(4, "0") ++ "/" ++
//				subjectID.asString.padLeft(4, "0") ++ "/" ++ name ++ kernelExt
//	});
//
//
//// open kernels -- load into buffers
//// load each L/R buffer separately
//b = kernelFiles.collect({ arg kernel;
//	chans.collect({ arg chan;
//		Buffer.readChannel(s, kernel, channels: [chan])
//		})
//	})
//)
//
//// decoding test... (to stereo UHJ)
//(
//{
//var sig;
//var bSig;
//var azim;
//var kernel;
//
//
//// kr....
//azim = MouseX.kr(pi, -pi);
//
//// test signal
//bSig = AtkMonoToB.ar(PinkNoise.ar, azim);
//
//~atkBtoUHJ.value(bSig, b)
//
//}.scope;
//)
//
//
//// free the buffers!!
//(
//b.shape.at(0).do({ arg i;
//	b.shape.at(1).do({ arg j;
//		b.at(i).at(j).free
//		})
//	})
//)
