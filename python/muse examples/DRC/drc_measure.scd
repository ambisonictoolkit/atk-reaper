/*

NOTE: edit this!!!

NOTE: use the Ctk version that Mike has prepared / edited as a model

This code implements is used to playback and record:

1) DUT measurement of the PreSonus, et al
2) RIR measurement of RM117


The sine sweeps to play back are:

1) 'dut_lss_signal.wav'
2) 'rm117_lss_signal.wav'


These signals are generate via Muse (Python), from te following files:

1) 'RM117_DUT_SS_generate.py'
2) 'RM117_RIR_SS_generate.py'


NOTE: set SR to 48000

*/


// ------------------------------------------------------------
// switch to local server, and boot
(
s = Server.local;
Server.default = s;

s.boot;
)

s.quit;

// have a look!
s.meter;



// ------------------------------------------------------------
// ------------------------------------------------------------
// DUT measurement
//
// Connect the 1st output of the device under test (DUT) into the 1st input.
// Rm 117: M32-DA>patchbay > Fireface800 Analog in > FF800 ADAT out > ADI 648 ADAT in > ADI Madi out > RME PCI card in > SC

// paths to the DUT log sine sweep, and resulting measurement
~dutPath = "/home/audioadmin/Desktop/RM117/SS filters/dut_lss_signal_96k.wav";
~dutMeasurePath = "/home/audioadmin/Desktop/RM117/SS filters/dut_lss_measure_96k.wav";


// file_name = 'dut'           # file name & dir
// working_dir = '/Users/josephla/Sound/Seattle Feb-Mar 2013/Turnkey DRC'
// filter_dir = 'filters'
// signal_dir = 'signals'



// read the DUT sweep into memory
~sweepBuf = Buffer.read(s, ~dutPath); // remember to free the buffer later.
~sweepBuf.plot; // have a look

// setup record buffer - same size as sweep
~recordBuf = Buffer.alloc(s, ~sweepBuf.numFrames, 1);

~sweepRecSynth = SynthDef(\sweepRecSynth, {arg out = 0, in = 0,
	playBufnum = 0, recBufnum = 1;
	RecordBuf.ar(SoundIn.ar(in), recBufnum, doneAction: 2, loop: 0);
	Out.ar(
		out,
        PlayBuf.ar(1, playBufnum, doneAction:2)
    )
}).load(s);

// ~sweepRecSynth.play(s,
// 	[
// 	\out, 0, \in, 0, \playBufnum, ~sweepBuf, \recBufnum, ~recordBuf]
// );

s.sendBundle(1.0, [\s_new, \sweepRecSynth, x = s.nextNodeID, 0, 1, \out, 0, \in, 16, \playBufnum, ~sweepBuf.bufnum, \recBufnum, ~recordBuf.bufnum]);

// have a look!
~recordBuf.plot;

// write to disk!
~recordBuf.write(~dutMeasurePath, "wav", "int24");


// free everything when done!!
~sweepBuf.free;
~recordBuf.free;




// ------------------------------------------------------------
// ------------------------------------------------------------
// RIR measurement
//
// Connect a measurement microphone (Earthworks M50) to the input.


// paths to the log sine sweep, and resulting measurement
(
~filterPath = "/home/audioadmin/Desktop/RM117/SS filters";
~signalName = "rm117_lss_signal_96k";
~measurementName = "rm117_lss_measure_96k_micFacingStraightUp_accuratePosition";

~fileExt = "wav";

// read the log sine sweep test signal into memory
~sweepBuf = Buffer.read(s, ~filterPath ++ "/" ++ ~signalName ++ "." ++ ~fileExt, action: {
	{~sweepBuf.plot}.defer;
	~recordBuf = Buffer.alloc(s, ~sweepBuf.numFrames, 1);
}); // remember to free the buffer later.
// ~sweepBuf.plot; // have a look
// ~sweepBuf.duration;

// setup record buffer - same size as sweep
// ~recordBuf = Buffer.alloc(s, ~sweepBuf.numFrames, 1); //moved to the sweepBuf action

~sweepRecSynth = SynthDef(\sweepRecSynth, {arg out = 0, in = 0,
	playBufnum = 0, recBufnum = 1, amp= 1;
	RecordBuf.ar(SoundIn.ar(in), recBufnum, doneAction: 2, loop: 0);
	Out.ar(
		out,
		PlayBuf.ar(1, playBufnum, doneAction:2) * amp
	)
}).load(s);
)

// s.sendBundle(1.0, [\s_new, \sweepRecSynth, x = s.nextNodeID, 0, 1, \out, 0, \in, 32, \playBufnum, ~sweepBuf.bufnum, \recBufnum, ~recordBuf.bufnum, \amp, 1]);

// gain = -12.dbamp gives ~-3db gain with Fireface 800 at gain of 3 o'clock [channel 0]

// measure here!
(
var measureRoutine, channels, sweepDur, micChannel = 16, startDelay = 20.0, additionalWaitTime = 0;
var recBufPath;

sweepDur = ~sweepBuf.numFrames / s.sampleRate;

// output channels to measure
// channels = [0, 1, 2, 3];
// channels = [27];
// channels = [0, 4, 9, 15, 23];
channels = (0..23);
//then
// channels = (24..27); //!! change gain for SUBs!

// post info!
("The measurement sequence will take " ++ (channels.size * (sweepDur + 1.0)).asString ++ " seconds.").postln;
// s.meter;

measureRoutine = Routine.new({

	// wait first
	startDelay.wait;

	channels.do({arg thisChannel;

		// post
		("Measuring channel:" ++ thisChannel.asString.padLeft(2, "0")).postln;

		//set path
		recBufPath = ~filterPath ++ "/" ++ ~measurementName ++ "_" ++ thisChannel.asString.padLeft(2, "0") ++ "." ++ ~fileExt;

		//check if file exists - check added by Marcin, who overwritten his previous measurement
		if(File.exists(recBufPath), {
			("File at " ++ recBufPath ++ " exists! Aborting channel number " ++ thisChannel).error;
			}, {

				// measure!
				s.sendBundle(
					0.5, [
						\s_new, \sweepRecSynth, x = s.nextNodeID, 0, 1,
						\out, thisChannel, \in, micChannel,
						\playBufnum, ~sweepBuf.bufnum, \recBufnum, ~recordBuf.bufnum,
						\amp, -12.dbamp // for mains
						// \amp, -18.dbamp // for sub
						// \amp, -24.dbamp // for sub - even lower gain - test
				]);

				// wait for sweep to complete + 0.5 sec
				(sweepDur + 0.5).wait;

				// write to disk!
				~recordBuf.write(recBufPath, ~fileExt, "int24"
				);
			});

		// wait for 0.5 sec + additional time
		(0.5 + additionalWaitTime).wait;
	});
});

measureRoutine.play;

);


// have a look!
~recordBuf.plot;

// free everything when done!!
(
~sweepBuf.free;
~recordBuf.free;
)