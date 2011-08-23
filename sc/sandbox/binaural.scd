// -------------------------------------------------------
// ATK (Sandbox)
//
// Coded by Joseph Anderson 2011
//
// -------------------------------------------------------

//------------------------------------------------------------------------
// Stereo decoders
//
//   b_to_binaural       HRTF decoding
//
//
//	NOTE: binaural kernels generated externally via Muse/ATK


//   b_to_binaural
(
~atkBtoBinaural = { arg bSig, decoderKernels;
	
	var res;

	// decode here!
	// convolve and mix to L/R
	res = decoderKernels.shape.at(0).collect({ arg i;	// harmonic [W, X, Y]
		decoderKernels.shape.at(1).collect({ arg j;	// channel [L, R]
			Convolution2.ar(
			bSig.at(i),
			decoderKernels.at(i).at(j),
			framesize: decoderKernels.at(i).at(j).numFrames)
		})
	});
	Mix.ar(res);								// mix
}
)


// binaural kernel reading
(
var kernelFiles;
var sourceDir, kernelDir;
var hrirDatabase;
var kernelSR;
var kernelSize;
var subjectID;
var kernelNames;
var kernelExt;
var chans;


// params
sourceDir		= "/Volumes/Audio/test";      			// corresponds to 'target_dir'
kernelDir		= "/ATK_kernels/FOA/decoders";			// corresponds to 'file_dir'

//hrirDatabase	= "/Spherical_HRIR";
//hrirDatabase	= "/Listen_HRIR";
//hrirDatabase	= "/CIPIC_HRIR";
hrirDatabase	= "/CIPIC_EQ_HRIR";

kernelNames	= [ "HRIR_W", "HRIR_X", "HRIR_Y", "HRIR_Z" ]; // kernel names (Listen, CIPIC)
//kernelNames	= [ "HRIR_W", "HRIR_X", "HRIR_Y"]; // kernel names (Spherical)
kernelExt		= ".wav";								// kernel file extension

kernelSR 		= s.sampleRate.asString;

kernelSize	= 256;								// for CIPIC database
//kernelSize	= 512;								// for Listen database

//kernelSize	= 512;								// spherical (Duda)
//kernelSize	= 1024;
//kernelSize	= 2048;
//kernelSize	= 4096;
//kernelSize	= 8192;

//subjectID		= 0004;							// Spherical subject 4

//subjectID		= 1002;							// Listen subject 0
//subjectID		= 1012;							// Listen subject 8

//subjectID		= 0003;							// CIPIC subject 0
//subjectID		= 0021;							// CIPIC KEMAR 0
subjectID		= 0165;							// CIPIC KEMAR 1



chans 		= 2;									// kernels are 2 channel, L/R


// collect kernel file names
kernelFiles = kernelNames.collect({ arg name;
	sourceDir ++ kernelDir ++ hrirDatabase ++"/SR_" ++ kernelSR.padLeft(6, "0") ++ "/N_" ++ 
				kernelSize.asString.padLeft(4, "0") ++ "/" ++
				subjectID.asString.padLeft(4, "0") ++ "/" ++ name ++ kernelExt
	});


// open kernels -- load into buffers
// load each L/R buffer separately
b = kernelFiles.collect({ arg kernel;
	chans.collect({ arg chan;
		Buffer.readChannel(s, kernel, channels: [chan])
		})
	})
)

// decoding test... (to stereo binaural)
(
{
var sig;
var bSig;
var azim;
var kernel;
var scale;

// params
scale = 12.neg.dbamp;


// kr....
azim = MouseX.kr(pi, -pi);

// test signal
bSig = scale * AtkMonoToB.ar(PinkNoise.ar, azim);

~atkBtoBinaural.value(bSig, b)

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
