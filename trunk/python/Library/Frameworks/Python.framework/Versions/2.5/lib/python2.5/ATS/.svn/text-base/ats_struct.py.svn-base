# python 2.5

# ATS
# by Juan Pampin
# pampin@u.washington.edu
#
# ported to python by
# Johnathan Lyon
# jgl@u.washington.edu
# 2007-05-23

#
# FILE: ats_struct.py
#
# This file contains the basic class definitions for 
# ATS and ATS helper objects



# Class: ats-sound
# ====================
# main data abstraction
# amp, frq, and pha contain sinusoidal
# modeling information as arrays of
# arrays of data arranged by partial
# par-energy and band-energy hold
# noise modeling information (experimental format)
class ats_sound:

    def __init__ (self, name="new-sound", 
		  # global sound info
		  sampling_rate=0, frame_size=0, window_size=0,  partials=0, frames=0,  bands=[], 
		  # Info deduced from analysis
		  optimized=None, ampmax=0.0,  frqmax=0.0,  frq_av=[], amp_av=[],  dur=0.0,
		  # Sinusoidal Data
		  time=[], frq=[],  amp=[], pha=[],
		  # Noise Data
		  energy=[], band_energy=[]):

	self.name = name 
	self.sampling_rate = sampling_rate
	self.frame_size = frame_size 
	self.window_size = window_size 
	self.partials = partials
	self.frames = frames 
	self.bands = bands
	self.optimized = optimized
	self.ampmax = ampmax 
	self.frqmax = frqmax 
	self.frq_av = frq_av
	self.amp_av = amp_av 
	self.dur = dur
	self.time = time
	self.frq = frq 
	self.amp = amp
	self.pha = pha
	self.energy = energy
	self.band_energy = band_energy

# Class: ats-fft
# ==================
# abstraction used to handle all 
# fft data in a single variable	
class ats_fft:
    
    def __init__ (self, size=0, rate=0.0, fd=[]):
	self.size = size
	self.rate = rate
	self.fd = fd


# Class: ats-peak
# ===================
# abstraction used to keep peak data used
# for peak detection and tracking
class ats_peak:

    def __init__ (self, amp=0.0, frq=0.0, pha=0.0, smr=0.0, track=0):
	self.amp = amp
	self.frq = frq
	self.pha = pha
	self.smr = smr
	self.track = track

# Class: ats-sieve
# ====================
# abstraction used for peak fitting
# by the sieve algorithm
class ats_sieve:

    def __init__ (self, ctrfrq=[], limits=[], tracks=[]):
	self.ctrfrq = ctrfrq
	self.limits = limits
	self.tracks = tracks
