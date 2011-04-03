#! /usr/bin/env python

# vim:syntax=python

# Joseph Anderson 2011

# NOTE: This module is designed to emulate the behaviour
#       of David Cournapeau's scikits.audiolab
#       [http://cournape.github.com/audiolab/]
#       by wrapping soundfile packages found
#       in the Python Standard Library.
#
#       The intention is to give 'minimal' soundfile
#       read/write capabilities without the need to
#       install Erik de Castro Lopo's libsndfile
#       [http://www.mega-nerd.com/libsndfile/].
#
#       Compatability is retained through the use of
#       a common user interface.


# Muse
# ==========

# This is the start of the Muse project, a Music V type language.

# Built on numpy
# """

# Privacy...
__all__ = ["available_file_formats", "available_encodings", "Format", "Sndfile"]

import copy
import numpy


# **************************
##_COMPAT_MODES = {"read": 'r', "write": 'w', "rwrite": 'rw'}
# do we need this??

_aformats = ['aifc', 'aiff', 'wav']

_aformats_descriptions_dic = {
    'aifc'  :   'AIFC (Apple/SGI)',
    'aiff'  :   'AIFF (Apple/SGI)',
    'wav'   :   'WAV (Microsoft)'
}

# complete set of encoding supported by aifc, wav:
##['pcmu8', 'pcms8', 'pcm16', 'pcm24', 'pcm32',\
##               'float32', 'float64', 'ulaw', 'alaw', 'ima_adpcm']
# NOTE: compression not currently supported
_aencodings = ['pcmu8', 'pcms8', 'pcm16', 'pcm24', 'pcm32',\
               'float32', 'float64']

_aencodings_descriptions_dic = {
    'pcmu8'     :   'Unsigned 8 bit PCM',
    'pcms8'     :   'Signed 8 bit PCM',
    'pcm16'     :   'Signed 16 bit PCM',
    'pcm24'     :   'Signed 24 bit PCM',
    'pcm32'     :   'Signed 32 bit PCM',
    'float32'   :   '32 bit float',
    'float64'   :   '64 bit float',
    'ulaw'      :   'U-Law',
    'alaw'      :   'A-Law',
    'ima_adpcm' :   'IMA ADPCM'
}

_aencodings_aifc    = _aencodings
_aencodings_aiff    = _aencodings

_aencodings_wav     = copy.copy(_aencodings)
_aencodings_wav.remove('pcms8')


# dictionary to associate available encodings
# for each file format
_aencodings_dict = {
    'aifc'  :   _aencodings_aifc,
    'aiff'  :   _aencodings_aiff,
    'wav'   :   _aencodings_wav
}

# complete set of endianness supported by scikits.pyaudiolab:
# ['file', 'little', 'big', 'cpu']
_aendianness = ['file', 'little', 'big']


# NOTE: AIFC files do support little endian, but this is resolved
#       in the AIFC file spec as 'sowt' compression. Support for
#       'sowt' could be added (for reading) through the replacement
#       of the ._read_comm_chunk method in class Aifc_read from the
#       aifc.py module from the Python Standard Library.
#
#       AIFC files with 'sowt' compression are 'standard' on OSX.

_aendianness_aifc   = copy.copy(_aendianness)
_aendianness_aifc.remove('little')

_aendianness_aiff   = _aendianness_aifc

_aendianness_wav    = _aendianness


# dictionary to associate available endianness
# for each file format
_aendianness_dict = {
    'aifc'  :   _aendianness_aifc,
    'aiff'  :   _aendianness_aiff,
    'wav'   :   _aendianness_wav
}


# **************************

def available_file_formats():
    return _aformats

def available_encodings(major):
    return _aencodings_dict[major]

def sndfile_version():  #e.g., libsndfile is not installed
    return (0, 0, 0, 0)

# **************************

#NOTE:  scikits.audiolab.Format is a 'type'
#       Format is defined as a 'new-style' class (which is also a type)
class Format(object):
    """Format(type=wav, encoding=pcm16, endianness=file)
    This class represents an audio file format. It knows about audio file
    format (wav, aiff, etc...), encoding (pcm, etc...) and endianness.
    
    Parameters
    ----------
    type : str
        the major file format (wav, etc...).
    encoding : str
        the encoding (pcm16, etc..).
    endianness : str
        the endianess.

    Notes
    -----
    The possible values for type, and encoding depend on those supported
    by the Python Standard Library's Multimedia Services (presently with
    the exception of data compression). You can query the possible values
    with the functions available_file_formats() and available_encodings().

    See also
    --------
    Sndfile class."""

    def __init__(self, type = 'wav', encoding = 'pcm16', endianness = 'file'):

        # Check if values are kosher
        if type not in _aformats:
            raise ValueError, 'file format %s not recognized' % type

        if encoding not in _aencodings:
            raise ValueError, 'encoding %s not recognized' % encoding

        if endianness not in _aendianness:
            raise ValueError, 'encoding %s not recognized' % endianness

        if encoding not in _aencodings_dict[type] or \
           endianness not in _aendianness_dict[type]:
            raise ValueError, 'The combination ' +\
                  '(type=wav|encoding=pcms8|endianness=file) ' +\
                  'you requested is not supported. You can use ' +\
                  'available_formats and available_encodings functions ' +\
                  'to query which formats and encodings are available.'

        # Store params...
        self._file_format   = type
        self._encoding      = encoding
        self._endianness    = endianness

        # ... and pretty names
        self._file_format_description   = _aformats_descriptions_dic[type]
        self._encoding_description      = _aencodings_descriptions_dic[encoding]


    # methods found in scikits.audiolab.pysndfile.compat
    # check copy and deepcopy... do we really need these?
    def __copy__(self):
        return copy.copy(self)

    def __deepcopy__(self):
        return self.__copy__()

    def __repr__(self):
        return 'Major Format: %s\n' % self._file_format_description +\
               'Encoding Format: %s\n' % self._encoding_description +\
               'Endianness: %s' % self._endianness

    def __str__(self):
        return 'Major Format: %s\n' % self._file_format_description +\
               'Encoding Format: %s\n' % self._encoding_description +\
               'Endianness: %s' % self._endianness

    # New-class properties are equivalent to CPython getset.
    # Use properties to emulate scikits.audiolab behaviours.
    # The following methods (to be hidden?) are used to return the
    # instance values (properties, below), which can only be set by
    # creating a new instance.
    # (OR calling an as yet unwritten set method.)

    def _get_file_format(self):
            return self._file_format

    def _get_encoding(self):
            return self._encoding

    def _get_endianness(self):
            return self._endianness

    def _get_file_format_description(self):
            return self._file_format_description

    def _get_encoding_description(self):
            return self._encoding_description

    # Properties
    file_format             = property(_get_file_format, \
                                       doc = 'File format (wav, etc...).')
    encoding                = property(_get_encoding, \
                                       doc = 'File encoding (pcm16, etc...).')
    endianness              = property(_get_endianness, \
                                       doc = 'File endianness (file, ' +\
                                       'little, etc...).')
    file_format_description = property(_get_file_format_description, \
                                       doc = 'File format description: the ' +\
                                       'full description from sndfile.')
    encoding_description    = property(_get_encoding_description, \
                                       doc = 'File encoding description: the' +\
                                       'full description from sndfile.')


#NOTE:  scikits.audiolab.Sndfile is a 'type'
#       Format is defined as a 'new-style' class (which is also a type)

##<attribute 'encoding'>
##<attribute 'endianness'>
##<attribute 'file_format'>

##<attribute 'format'>
##<attribute 'channels'>
##<attribute 'nframes'>
##<attribute 'samplerate'>
##
##<method 'read_frames'>
##<method 'seek'>
##<method 'sync'>
##<method 'write_frames'>
##<method 'close'>
class Sndfile(object):

    pass
