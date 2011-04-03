#! /usr/bin/env python

# vim:syntax=python

# Joseph Anderson 2007

# NOTE: this module is to be deprecated
#       the functionality to read/write multiple
#       soundfiles is to be replace through the use
#       of 'external', command line driven (de)interleaving
#       by sndfile-(de)interleave [http://www.mega-nerd.com/libsndfile/]


# Muse
# ==========

# This is the start of the Muse project, a Music V type language.

# Built on numpy an pyaudiolab.
# """

from numpy import *
# from pyaudiolab import *
from scikits.audiolab import *
from muse import *


# **************************


# need to rewrite sndfiles. . . ?
# and create a wrapper for sndfile??

# class sndfiles(sndfile):
class sndfiles:
    """
    Class to open, read and write lists of audio files
    
    A loose wrapper for pyaudiolab's sndfile()
    """

    def __init__(self, filenames, mode = 'read', format = None, fchannels = 0, \
            samplerate = 0):
        """sndfiles(filenames, mode = 'read', format = None, 
        fchannels = 0, samplerate = 0):
        Args:
            - filenames  : list of names of the files to open (string)
            - mode      : 'read' for read, 'write' for write, 
            or 'rwrite' for read and write.
            - format    : (write modes and raw files only) when opening a new file for writing, 
            give the format to write in. Should be an instance of the formatinfo class.
            - fchannels : (write mode and raw files only) number of channels for split output files
            - samplerate: (write mode and raw files only) sampling rate

        Returns on success a valid sndfiles object
        """

        self.filenames = filenames

        # Open the files. . . and do error checking
        self.sfiles = [ sndfile(filename, mode, format, fchannels, samplerate) for filename in filenames ]

        # check that attributes are the same
        if mode == 'read' or mode == 'rwrite':

            # check for same number of channels
            if cmp([ self.sfiles[0].get_channels() ] * len(self.sfiles),
                   [ sfile.get_channels() for sfile in self.sfiles ]):

                raise Exception("Number of channels is not consistent across input files!")
            else:
                self.fchannels = self.sfiles[0].get_channels()

            # check for same sample rate
            if cmp([ self.sfiles[0].get_samplerate() ] * len(self.sfiles),
                   [ sfile.get_samplerate() for sfile in self.sfiles ]):

                raise Exception("Sample rate is not consistent across input files!")
            else:
                self.samplerate = self.sfiles[0].get_samplerate()

            # check for same number of frames
            if cmp([ self.sfiles[0].get_nframes() ] * len(self.sfiles),
                   [ sfile.get_nframes() for sfile in self.sfiles ]):

                raise Exception("Number of frames is not equal across input files!")

            # check for same file format
            if cmp([ self.sfiles[0].get_file_format() ] * len(self.sfiles),
                   [ sfile.get_file_format() for sfile in self.sfiles ]):

                raise Exception("File format is not consistent across input files!")

            # check for same encoding
            if cmp([ self.sfiles[0].get_encoding() ] * len(self.sfiles),
                   [ sfile.get_encoding() for sfile in self.sfiles ]):

                raise Exception("Encoding is not consistent across input files!")

            # check for same endianness
            if cmp([ self.sfiles[0].get_endianness() ] * len(self.sfiles),
                   [ sfile.get_endianness() for sfile in self.sfiles ]):

                raise Exception("Endianness is not consistent across input files!")

        # shouldn't need to do anything else for write mode
        # as it should be handled by sndfile
#         elif mode == 'write':
#             pass


    # not sure if we need the __del__ method
    def __del__(self):
        for sfile in self.sfiles:
            sfile.close()

    def close(self):
        for sfile in self.sfiles:
            sfile.close()

    def sync(self):
        """call the operating system's function to force the writing of 
        all file cache buffers to disk the file. No effect if file is open as read"""
        for sfile in self.sfiles:
            sfile.sync()

    def seek(self, offset, whence = 0):
        return [ sfile.seek(offset, whence) for sfile in self.sfiles ][0]

    # Functions to get information about the file
    def get_nframes(self):
        """ Return the number of frames of the files"""
        return self.sfiles[0].get_nframes()
    
    def get_samplerate(self):
        """ Return the samplerate in Hz of the files"""
        return self.sfiles[0].get_samplerate()
    
    def get_fchannels(self):
        """ Return the number of channels for output files"""
        return self.sfiles[0].get_channels()
    
    def get_file_format(self):
        """return user friendly file format string"""
        return self.sfiles[0].get_file_format()

    def get_encoding(self):
        """return user friendly encoding string"""
        return self.sfiles[0].get_encoding()

    def get_endianness(self):
        """return user friendly file format string"""
        return self.sfiles[0].get_endianness()

    #------------------
    # Functions to read
    #------------------
    def read_frames(self, nframes, dtype = float64):
        """Read nframes frames, and returns an array of the required
        type.  
        
        Note: - read_frames updates the read pointer.
              - One column is one channel."""

        a = asarray([ sfile.read_frames(nframes, dtype) for sfile in self.sfiles ])

        if a.ndim > 2:
            return column_stack(tuple(a))

        else:
            return interleave(a)


    #-------------------
    # Functions to write
    #-------------------
    def write_frames(self, input, nframes):
        """ write nframes frames of input array into the files. 
        
        Notes: - one channel is one column
               - updates the write pointer."""

        # First, get the number of channels and frames from input
        if not(input.ndim == 2):
            raise Exception("Expect array of rank = 2, got %d" % input.ndim)

        else:
            (foo, nc)   = input.shape

        # Number of channels should be the one expected
        ec = len(self.sfiles) * self.get_fchannels()

        if not(nc == ec):
            raise Exception("Expected %d channels, got %d" % (ec, nc))

        # Massage input into correct shape for output
        if self.get_fchannels() > 1:
            x = hsplit(input, len(self.sfiles))

        else:
            x = deinterleave(input)

        # Writing to the file
        for (sfile, y) in zip(self.sfiles, x):
            sfile.write_frames(y, nframes)


    # Syntactic sugar
    def __repr__(self):
        return self.__str__()

    def __str__(self):
        repstr = "----------------------------------------\n"
        for (i, filename) in zip(range(len(self.filenames)), self.filenames):
            repstr  += "File %d        : %s\n" % ((i+1), filename)
        repstr  += "File Channels : %d\n" % self.fchannels
        repstr  += "Sample rate   : %d\n" % self.samplerate
        repstr  += "Frames        : %d\n" % self.get_nframes()
        repstr  += "Raw Format    : %#010x -> %s\n" % \
                (self.sfiles[0]._format.get_format_raw(), self.sfiles[0]._format.get_major_str())
        repstr  += "File format   : %s\n" % self.get_file_format()
        repstr  += "Encoding      : %s\n" % self.get_encoding()
        repstr  += "Endianness    : %s\n" % self.get_endianness()
        repstr  += "Sections      : %d\n" % self.sfiles[0]._sfinfo.sections
        if self.sfiles[0]._sfinfo.seekable:
            seek    = 'True'
        else:
            seek    = 'False'
        repstr  += "Seekable      : %s\n" % seek
        repstr  += "Duration      : %s\n" % self._generate_duration_str()
        return repstr

    def _generate_duration_str(self):
        if self.samplerate < 1:
            return None
        tsec    = self.get_nframes() / self.samplerate
        hrs     = tsec / 60 / 60
        tsec    = tsec % (60 ** 2)
        mins    = tsec / 60
        tsec    = tsec % 60
        secs    = tsec
        ms      = 1000 * self.get_nframes() / self.samplerate % 1000

        return "%02d:%02d:%02d.%3d" % (hrs, mins, secs, ms)
