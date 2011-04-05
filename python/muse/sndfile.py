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
import numpy as np
import sndhdr
import aifc
import wave


# **************************
# NOTE: '_a' appended to var names = '_available_'

# At present, no aifc compressions are supported, not even the
# 'sowt' flagged little endian (byte reverse, uncompressed) format.
# It could be the thing to do is to remove 'aifc' from the list below.
#
# The danger is, that scikits.audiolab will read/write 'sowt' files when
# using the 'aiff' encoding flag, where this wrapper for PSL modules will
# not.
#
# For the time being, we'll keep both 'aifc' and 'aiff' here.
#
# Also note, a filename ending with '.aiff' forces the use of PSL aifc.aiff().
# This may or may not cause difficulties.
_aformats = ['aifc', 'aiff', 'wav']

_aformats_descriptions_dic = {
    'aifc'  :   'AIFC (Apple/SGI)',
    'aiff'  :   'AIFF (Apple/SGI)',
    'wav'   :   'WAV (Microsoft)'
}

# complete set of encoding supported by aifc, wav:
# ['pcmu8', 'pcms8', 'pcm16', 'pcm24', 'pcm32',\
#               'float32', 'float64', 'ulaw', 'alaw', 'ima_adpcm']
#
# NOTE: compression not currently supported
#
#       PSL's aifc regards the following encodings as compressed:
#           ['pcmu8', 'pcm24'] = 'raw'
#           ['float32'] = 'FL32'
#
#       PSL's wave regards the following encodings as compressed:
#           ['float32'] = 'WAVE_FORMAT_IEEE_FLOAT'
#
#       In other words, no support for float with PSL
_aencodings = ['pcmu8', 'pcms8', 'pcm16', 'pcm24', 'pcm32',\
               'float32']

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

_aencodings_dtype_dic = {
    'pcmu8'     :   'u',
    'pcms8'     :   'i',
    'pcm16'     :   'i',
    'pcm24'     :   'i',
    'pcm32'     :   'i',
    'float32'   :   'f',
    'float64'   :   'f'
}

# set-up individual type support for encodings
# 'aifc'
_aencodings_aifc    = copy.copy(_aencodings)
_aencodings_aifc.remove('pcmu8')
_aencodings_aifc.remove('pcm24')
_aencodings_aifc.remove('float32')

# 'aiff'
_aencodings_aiff    = _aencodings_aifc

# 'wav'
_aencodings_wav     = copy.copy(_aencodings)
_aencodings_wav.remove('pcms8')
_aencodings_wav.remove('float32')


# dictionary to associate available encodings
# for each file format
# NOTE: We may like to restrict this to 'aiff',
#       and force the returned type string from PSL aifc
#       to 'aiff'.
_aencodings_dict = {
    'aifc'  :   _aencodings_aifc,
    'aiff'  :   _aencodings_aiff,
    'wav'   :   _aencodings_wav
}

# complete set of endianness supported by scikits.pyaudiolab:
# ['file', 'little', 'big', 'cpu']
#
# PSL sndfile modules automatically set endianness:
#       'aifc'/'aiff'   = 'big'
#       'wav'           = 'cpu'
#
# The most equivalent support seems to be to set the default to 'file'
# which will set endianness:
#       'aifc'/'aiff'   = 'big'
#       'wav' (intel)   = 'little'
#       'wav' (ppc)     = 'big'
#
# This results in an inconsistency with scikits.pyaudiolab
_aendianness = ['file']


# NOTE: AIFC files do support little endian, but this is resolved
#       in the AIFC file spec as 'sowt' compression. Support for
#       'sowt' could be added (for reading) through the replacement
#       of the ._read_comm_chunk method in class Aifc_read from the
#       aifc.py module from the Python Standard Library.
#
#       AIFC files with 'sowt' compression are 'standard' on OSX.

##_aendianness_aifc   = copy.copy(_aendianness)
##_aendianness_aifc.remove('little')
##
##_aendianness_aiff   = _aendianness_aifc

_aendianness_aifc   = _aendianness
_aendianness_aiff   = _aendianness
_aendianness_wav    = _aendianness


# dictionary to associate available endianness
# for each file format
_aendianness_dict = {
    'aifc'  :   _aendianness_aifc,
    'aiff'  :   _aendianness_aiff,
    'wav'   :   _aendianness_wav
}


# dictionary to associate available PSL soundfile modules
# for each file format
_apsl_sndfile_dict = {
    'aifc'  :   aifc,
    'aiff'  :   aifc,
    'wav'   :   wave
}

# List of numpy floats
# NOTE: Not all of these are supported on all systems
##_numpy_floats = [ np.half, np.single, np.double, np.float_, np.longfloat, \
##                 np.float16, np.float32, np.float64, np.float96, np.float128 ]
_numpy_floats = [ np.single, np.double, np.float_, np.longfloat, \
                 np.float32, np.float64, np.float128 ]

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
class Sndfile(object):
    """Sndfile(filename, mode=r, Format format=None, int channels=0, int samplerate=0)

    Sndfile is the core class to read/write audio files. Once an instance is
    created, it can be used to read and/or writes data from numpy arrays, query
    the audio file meta-data, etc...

    Parameters
    ----------
    filename : string or int
       name of the file to open (string), or file descriptor (integer)
    mode : string
       'r' for read or 'w' for write. (Mode 'rw' for read and
       write is not supported.)
    format : Format
       Required when opening a new file for writing, or to read raw audio
       files (without header).
    channels : int
       number of channels.
    samplerate : int
       sampling rate.

    Returns
    -------
       sndfile: as Sndfile instance.

    Notes
    -----
    format, channels and samplerate need to be given only in the write modes
    and for raw files."""

    def __init__(self, filename, mode = 'r', \
                format = None, channels = 0, samplerate = 0):

        # Check if values are kosher
        if type(filename) is not str:
            raise ValueError, 'filename should be a string'

        if mode not in ['r', 'w']:
            raise ValueError, 'mode %s not recognized\nnote: rw not supported' \
                  % mode

        if type(channels) is not int or type(samplerate) is not int:
            raise ValueError, 'an integer is required'

        if mode is 'w':
            if type(format) is not Format:
                raise ValueError, 'For write mode, you should provide a ' +\
                  'format argument!'

            if channels <= 0 or samplerate <= 0:
                raise ValueError, 'Bad value of samplerate (%s)' % samplerate +\
                    ' or channels (%s)' % channels

        # Set-up for reading... (including opening files!)
        if mode is 'r':
            # 1) retrieve header: sndhdr determines type, etc.
            #    (type, sampling_rate, channels, frames, bits_per_sample)
            _sndhdr = sndhdr.what(filename)

            # 2) populate params
            #    NOTE: reading nframes is deferred until the file is
            #          opened, as sndhdr returns (-1) for wav files.
            #          In other words, _sndhdr[3] is not set correctly for wav.

            #    NOTE: sndhdr doesn't give a robust answer for encoding.
            #          Instead, we need to intuit it from file_type and
            #          bits_per_sample. This only works because the encodings
            #          supported by PSL are at present very limited!
            file_type   = _sndhdr[0]
            bps         = _sndhdr[4]
            
            if bps is not 8:
                encoding    = 'pcm' + str(bps)
            elif file_type is 'wav':
                encoding    = 'pcmu' + str(bps)
            elif file_type is 'aiff' or file_type is 'aifc':
                encoding    = 'pcms' + str(bps)
            else:
                raise ValueError, "Couldn't determine correct encoding."

            self._format        = Format(file_type, encoding, 'file')
            self._samplerate    = _sndhdr[1]
            self._channels      = _sndhdr[2]
            self._bps           = bps

            # 3) assign appropriate PSL soundfile read/writer,
            #    e.g., psl_sndfile = wave
            _apsl_sndfile = _apsl_sndfile_dict[file_type]
           
            # 4) open w/ PSL... and read nframes (see note above)
            self._sndfile = _apsl_sndfile.open(filename, 'rb')
            self._nframes = long(self._sndfile.getnframes())


        # Set-up for writing, e.g., 'w'... (including opening files!)
        else:
            # 1) populate params,
            self._format        = format
            self._channels      = channels
            self._samplerate    = samplerate
            self._nframes       = long(-1)
            self._bps           = None  #Need to populate via a dict
            
            # 2) assign appropriate PSL soundfile read/writer,
            #    e.g., psl_sndfile = wave
            _apsl_sndfile = _apsl_sndfile_dict[file_type]

            # 3) open w/ PSL
            self._sndfile = _apsl_sndfile.open(filename, 'wb')



    # New-class properties are equivalent to CPython getset.
    # Use properties to emulate scikits.audiolab behaviours.
    # The following methods (to be hidden?) are used to return the
    # instance values (properties, below), which can only be set by
    # creating a new instance.
    # (OR calling an as yet unwritten set method.)

    def _get_channels(self):
            return self._channels

    def _get_nframes(self):
            return self._nframes

    def _get_samplerate(self):
            return self._samplerate

    def _get_format(self):
            return self._format

    def _get_file_format(self):
            return self._format.file_format

    def _get_encoding(self):
            return self._format.encoding

    def _get_endianness(self):
            return self._format.endianness

    # Properties
    channels                = property(_get_channels, \
                                       doc = 'Number of channels.')
    nframes                 = property(_get_nframes, \
                                       doc = 'Number of frames of the file.')
    samplerate              = property(_get_samplerate, \
                                       doc = 'Sampling rate (in Hz).')
    format                  = property(_get_format, \
                                       doc = 'Format instance attached to ' +\
                                       'the Sndfile instance.')
    file_format             = property(_get_file_format, \
                                       doc = 'File format (wav, etc...).')
    encoding                = property(_get_encoding, \
                                       doc = 'File encoding (pcm16, etc...).')
    endianness              = property(_get_endianness, \
                                       doc = 'File endianness (file, ' +\
                                       'little, etc...).')

    def close(self):
        """Close the file."""
        self._sndfile.close()

    #------------------
    # Functions to read
    #------------------
    def read_frames(self, nframes, dtype = np.float64):
        """Read the given number of frames and put the data into a numpy array of
        the requested dtype.
        
        Parameters
        ----------
        nframes : int
            number of frames to read.
        dtype : numpy dtype
            dtype of the returned array containing read data (see note).
        
        Notes
        -----
        One column per channel.
        
        Updates the read pointer.
        
        Notes
        -----
        if float are requested when the file contains integer data, you will
        get normalized data (that is the max possible integer will be 1.0, and
        the minimal possible value -1.0).
        
        if integers are requested when the file contains floating point data,
        it may give wrong results because there is an ambiguity: if the
        floating data are normalized, you can get a file with only 0 ! Getting
        integer data from files encoded in normalized floating point is not
        supported (this is an audiolab limitation: sndfile supports it)."""

        # First, set up...
        
        # Read sample width (in bytes)
        # NOTE: this could be determined from self._bps
        sampwidth = self._sndfile.getsampwidth()

        # Determine actual endianness...
        #           ... from defaults for standard file spec
        # NOTE: this may result in a problem for 'wav' files on ppc
        if self._format.endianness is 'file':
            if self._format.file_format is 'wav':
                endian_flag = '<'       # 'wav' is little endian

            elif self._format.file_format is 'aiff' \
                 or self._format.file_format is 'aifc':
                endian_flag = '>'       # 'aiff' is big endian

        # file data type (as a numpy dtype)
        file_dtype = endian_flag +\
                     _aencodings_dtype_dic[self._format.encoding] +\
                     str(sampwidth)

        # Read frames, returned 'raw' as a string containing for each frame
        # the samples of all channels.
        raw_frames = self._sndfile.readframes(nframes)
        
        # Catch 'pcm24', and 'upsample' to 'pcm32'...
        #               ...re-setting file_dtype correctly
        # NOTE: for some reason 'is' doesn't work, use '==' instead
        if self._format.encoding == 'pcm24':
            tmp_frames = ''             #temp holder to build upsamp str

            if endian_flag is '<':      #little endian
                for x in range(len(raw_frames) / sampwidth):

                    #extract sample and append a padding byte
                    samp = raw_frames[x * sampwidth:(x + 1) * sampwidth]

                    #negative or positive?
                    if ord(samp[2]) & 0x80:
                        tmp_frames += (samp + '\xff')   #neg val
                    else:
                        tmp_frames += (samp + '\x00')   #pos val     

            else:                       #big endian
                for x in range(len(raw_frames) / sampwidth):

                    #extract sample and append a padding byte
                    samp = raw_frames[x * sampwidth:(x + 1) * sampwidth]

                    #negative or positive?
                    if ord(samp[0]) & 0x80:
                        tmp_frames += ('\xff' + samp)   #neg val
                    else:
                        tmp_frames += ('\x00' + samp)   #pos val     

            #assign to raw_frames
            raw_frames = tmp_frames

            #assign new file_dtype (to 'pcm32')
            file_dtype = endian_flag +\
                         _aencodings_dtype_dic[self._format.encoding] +\
                         '4'
            

        # Convert to numpy array
        frames = np.fromstring(
            raw_frames,
            dtype = file_dtype
            ).astype(dtype)

        # Catch 'pcmu8', to center around 0
        # NOTE: for some reason 'is' doesn't work, use '==' instead
        if self._format.encoding == 'pcmu8':
            frames -= pow(2, self._bps - 1)

        # If dtype is float, scale to +/-1
        if dtype in _numpy_floats:
            frames *= pow(2, 1 - self._bps)

        return frames


    def seek(self):
        """Seek into audio file: similar to python seek function, taking only in
        account audio data.
        
        Parameters
        ----------
        offset : int
            the number of frames (eg two samples for stereo files) to move
            relatively to position set by whence.
        whence : int
            only 0 (beginning), 1 (current) and 2 (end of the file) are
            valid.
        mode : string
            If set to 'rw', both read and write pointers are updated. If
            'r' is given, only read pointer is updated, if 'w', only the
            write one is (this may of course make sense only if you open
            the file in a certain mode).
        
        Returns
        -------
        offset : int
            the number of frames from the beginning of the file
        
        Notes
        -----
        
        Offset relative to audio data: meta-data are ignored.
        
        if an invalid seek is given (beyond or before the file), an IOError is
        launched; note that this is different from the seek method of a File
        object."""

        pass

    def sync(self):
        """call the operating system's function to force the writing of all
        file cache buffers to disk the file.
        
        No effect if file is open as read"""

        pass

    def write_frames(self):
        """write given number frames into file.
        
        Parameters
        ----------
        input : ndarray
            array containing data to write.
        
        Notes
        -----
        One column per channel.
        
        updates the write pointer.
        
        if the input type is float, and the file encoding is an integer type,
        you should make sure the input data are normalized normalized data
        (that is in the range [-1..1] - which will corresponds to the maximum
        range allowed by the integer bitwidth)."""

        pass


##    read_frames(...)
##        Sndfile.read_frames(self, scikits.audiolab.pysndfile.sndfile.sf_count_t nframes, dtype=<???>)
##        Read the given number of frames and put the data into a numpy array of
##        the requested dtype.
##    
##    seek(...)
##        Sndfile.seek(self, scikits.audiolab.pysndfile.sndfile.sf_count_t offset, int whence=0, mode=rw)
##        Seek into audio file: similar to python seek function, taking only in
##        account audio data.
##    
##    sync(...)
##        Sndfile.sync(self)
##        call the operating system's function to force the writing of all
##        file cache buffers to disk the file.
##        
##        No effect if file is open as read
##    
##    write_frames(...)
##        Sndfile.write_frames(self, ndarray input)
##        write given number frames into file.

#['__class__', '__delattr__', '__doc__', '__format__', '__getattribute__', '__hash__', '__init__', '__new__', '__pyx_vtable__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', 'channels', 'close', 'encoding', 'endianness', 'file_format', 'format', 'nframes', 'read_frames', 'samplerate', 'seek', 'sync', 'write_frames']
