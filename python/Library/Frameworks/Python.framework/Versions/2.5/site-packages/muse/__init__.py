"""\
Muse
==========

This is the start of the Muse project, a Music V type language.

Built on pyaudiolab.
"""

# Likely need to sort out all the various imports.
# Of particular interest, name collisions
# between numpy array functions and math / python functions

# __all__ var, list muse modules here
# __all__ = ["muse", "filters", "decoders", "encoders", "transforms"]

# import numpy and pyaudiolab into muse namespace
from numpy import *
# from pyaudiolab import *
from scikits.audiolab import *
# from scikits.samplerate import *
# import scikits.samplerate
from scikits.samplerate import resample as src


# now import the muse modules into muse namespace
from sndfiles import *
from muse import *
# from generators import *
from filters import *
from decoders import *
from encoders import *
from transforms import *
from generators import *
from nonlinear import *
from analysis import *

# import muse
# import filters
# import decoders
# import encoders
# import transforms

print "Welcome to Muse!!"
