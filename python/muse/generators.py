#! /usr/bin/env python

# vim:syntax=python

# Joseph Anderson 2007


# Muse
# ==========

# This is the start of the Muse project, a Music V type language.



# seem to need to have these imports here. . . to make sure names are defined
from muse import *
import itertools as it 
# from numpy import *
from scipy.signal import *
from numpy.random import *


# #=========================
# # Definition of constants
# #=========================


# *********************************************************************
# functional extensions to numpy. . .

# envelopes. . .

def lin(vals = [0., 1.], nframes = 9):
    """lin(vals = [0., 1.], nframes = 9)

    Return evenly spaced numbers.
    
    Return nframes of evenly spaced samples between vals.

    """
    return linspace(vals[0], vals[1], nframes)


# need to adjust import order. . . exp is being replaced
# by numpy.exp
def exp(vals = [1., 2.], nframes = 9):
    """exp(vals = [1., 2.], nframes = 9)

    Return evenly spaced numbers on an exponential scale.
    
    Return nframes of evenly spaced samples between vals.

    """
    vals = log2(vals)
    return logspace(vals[0], vals[1], nframes, base = 2)


def linseg(vals = [0., 1., 0.], nframes = [5, 9]):
    """linseg(vals = [0., 1., 0.], nframes = [5, 9])

    Return evenly spaced numbers.
    
    Return nframes of evenly spaced samples between vals.     """
    nframes = append([1], nframes)

    res = array([])
    for i in range(len(nframes) - 1):
        res = append(res[:-1], lin([vals[i], vals[i + 1]], 1 + nframes[i + 1] - nframes[i]))
    return res


def expseg(vals = [1., 2., 1.], nframes = [5, 9]):
    """expseg(vals = [1., 2., 1.], nframes = [5, 9])

    Return evenly spaced numbers on an exponential scale.
    
    Return nframes of evenly spaced samples between vals.

    """
    nframes = append([1], nframes)

    res = array([])
    for i in range(len(nframes) - 1):
        res = append(res[:-1], exp([vals[i], vals[i + 1]], 1 + nframes[i + 1] - nframes[i]))
    return res


def linen(nframes = [3, 7, 9], shape = 'line'):
    """linen(nframes = [3, 7, 9], shape = 'line')

    Linen shaped envelope of nframes.

    Shapes: line, sine, hann
    
    """
    incr = linseg([0., 1., 1., 0.], nframes)

    try:
        res = {
            'line': (lambda: incr), 
            'sine': (lambda: sin(.5 * pi * incr)),
            'hann': (lambda: sin(.5 * pi * incr)**2)
            }[shape]()

        return res

    except KeyError:
        print "Invalid shape: %s" % shape


# def win(nframes = 9, warp = 0., shape = 'hann'):
#     """win(nframes = 9, warp = 0., shape = 'hann')

#     Return a window array of nframes.
#     Shapes: sine, hann, bart
#     Warp: warps increment by incr**(2**warp)

#     Endpoints are zeros. (Or nearly!)

#     """
#     incr = (linspace(0., 1., nframes))**(2**warp)

#     try:
#         res = {
#             'sine': (lambda: sin(pi * incr)),
#             'hann': (lambda: sin(pi * incr)**2),
#             'bart': (lambda: (1. - 2. / (nframes  - 1.) * abs(.5 * (nframes  - 1.) * (2. * incr - 1.)))), 
#             }[shape]()

#         return res

#     except KeyError:
#         print "Invalid shape: %s" % shape
def win(nframes = 9, warp = .5, shape = 'hann'):
    """win(nframes = 9, warp = .5, shape = 'hann')

    Return a window array of nframes.
    Shapes: sine, hann, bart
    Warp: warps increment by Wn_warp

    Endpoints are zeros. (Or nearly!)

    """
    Wn = lin([0., 1.], nframes)
    incr = Wn_warp(Wn, 1.-warp, Wn_w = True)

    try:
        res = {
            'sine': (lambda: sin(pi * incr)),
            'hann': (lambda: sin(pi * incr)**2),
            'bart': (lambda: (1. - 2. / (nframes  - 1.) * abs(.5 * (nframes  - 1.) * (2. * incr - 1.)))), 
            }[shape]()

        return res

    except KeyError:
        print "Invalid shape: %s" % shape


# random lines. . .
def rnd_linseg(vals = [0., 1.], npoints = 3, nframes = 9):
    """rnd_linseg(vals = [0., 1.], npoints = 3, nframes = 9)

    Return evenly spaced numbers.
    
    Return nframes of evenly spaced samples between vals [low, high],
    with npoints breakpoints.     """

    # generate breakpoints
    rnframes = unique(random_integers(1, nframes - 1, npoints - 2))

    if rnframes[-1] != nframes:
        rnframes = append(rnframes, nframes)


    # generate values at breakpoints
    rvals = uniform(vals[0], vals[1], len(rnframes) + 1)

    return linseg(rvals, rnframes)

# noise. . .

def white(nframes = 9):
    """white(nframes = 9)

    Return nframes of uniformly distributed random numbers.

    (A convenience function for random.uniform)
    
    nframes may be a shape tuple for multichannel.

    """
    return uniform(-1., 1., nframes)


# fixed osc. . .
# variable osc is found in filters.py
# as phase incrementor needs state

def fphasor(Wn = .5, phase = 0., nframes = 9):
    """fphasor(Wn = .5, phase = 0., nframes = 9)

    Return nframes of a fixed frequency un-wrapped phasor.
    Wn is the normalized frequency (may be an array, for multichannel).
    Phase in radians (may be an array, for multichannel).

    """
    if isscalar(Wn) and isscalar(phase):
        n = arange(0, nframes)
    else:
        n = reshape(arange(0, nframes), (nframes, 1))
    
    return (n * pi * Wn + phase)


def fsinosc(Wn = .5, phase = 0., nframes = 9):
    """fsinosc(Wn = .5, phase = 0., nframes = 9)

    Return nframes of a fixed frequency sine oscillator.
    Wn is the normalized frequency (may be an array, for multichannel).
    Phase in radians (may be an array, for multichannel).

    """
    return sin(fphasor(Wn, phase, nframes))


def fcososc(Wn = .5, phase = 0., nframes = 9):
    """fcososc(Wn = .5, phase = 0., nframes = 9)

    Return nframes of a fixed frequency cosine oscillator.
    Wn is the normalized frequency (may be an array, for multichannel).
    Phase in radians (may be an array, for multichannel).

    """
    return cos(fphasor(Wn, phase, nframes))


# *********************************************************************
# generators

def xlinspace(start, stop, num = 50):
    """xlinspace(start, stop, num = 50)

    Return evenly spaced numbers.

    Like linspace(), but instead of returning a list, returns an object that
    generates the numbers in the range on demand.  For looping, this is 
    slightly faster than linspace() and more memory efficient.
    
    Return num evenly spaced samples from start to stop.  The last sample is stop.

    """
    for x in xrange(num):
        yield float(start) + (float(stop) - float(start)) / (num - 1) * x


def xlogspace(start, stop, num = 50, base = 10.):
    """xlogspace(start, stop, num = 50)

    Return evenly spaced numbers on a logarithmic scale.

    Like logspace(), but instead of returning a list, returns an object that
    generates the numbers in the range on demand.  For looping, this is 
    slightly faster than logspace() and more memory efficient.
    
    Return num evenly spaced exponents from base**start to
    base**stop. The last number is base**stop

    """
    for x in xrange(num):
        yield pow(base, float(start) + (float(stop) - float(start)) / (num - 1) * x)


def lin_gen(vals = [0., 1.], nframes = 9):
    """lin_gen(vals = [0., 1.], nframes = 9)

    Like lin, but instead of returning an array, returns an object that
    generates the numbers in the range on demand.  For looping, this is 
    slightly faster than lin and more memory efficient.
    
    Return evenly spaced numbers.
    
    Return nframes of evenly spaced samples between vals.

    """
    for x in xlinspace(vals[0], vals[1], nframes):
        yield x


def exp_gen(vals = [1., 2.], nframes = 9):
    """exp_gen(vals = [1., 2.], nframes = 9)

    Like exp(), but instead of returning an array, returns an object that
    generates the numbers in the range on demand.  For looping, this is 
    slightly faster than exp() and more memory efficient.
    
    Return evenly spaced numbers on an exponential scale.
    
    Return nframes of evenly spaced samples between vals.

    """
    vals = log2(vals)
    for x in xlogspace(vals[0], vals[1], nframes, base = 2):
        yield x


def linseg_gen(vals = [0., 1., 0.], nframes = [5, 9]):
    """linseg_gen(vals = [0., 1., 0.], nframes = [5, 9])

    Like linseg(), but instead of returning an array, returns an object that
    generates the numbers in the range on demand.  For looping, this is 
    slightly faster than linseg and more memory efficient.
    
    Return evenly spaced numbers.
    
    Return nframes of evenly spaced samples between vals.     """
    nframes = append([1], nframes)

    res = iter([])
    for i in xrange(len(nframes) - 1):
        res = it.chain(it.islice(res, nframes[i] - 1), lin_gen([vals[i], vals[i + 1]], 1 + nframes[i + 1] - nframes[i]))
    return res


def expseg_gen(vals = [1., 2., 1.], nframes = [5, 9]):
    """expseg(vals = [1., 2., 1.], nframes = [5, 9])

    Like expseg(), but instead of returning an array, returns an object that
    generates the numbers in the range on demand.  For looping, this is 
    slightly faster than expseg() and more memory efficient.
    
    Return evenly spaced numbers on an exponential scale.
    
    Return nframes of evenly spaced samples between vals.

    """
    nframes = append([1], nframes)

    res = iter([])
    for i in xrange(len(nframes) - 1):
        res = it.chain(it.islice(res, nframes[i] - 1), exp_gen([vals[i], vals[i + 1]], 1 + nframes[i + 1] - nframes[i]))
    return res


def linen_gen(nframes = [3, 7, 9], shape = 'line'):
    """linen_gen(nframes = [3, 7, 9], shape = 'line')

    Like linen(), but instead of returning an array, returns an object that
    generates the numbers in the range on demand.  For looping, this is 
    slightly faster than linen() and more memory efficient.
    
    Linen shaped envelope of nframes.

    Shapes: line, sine, hann
    
    """

    for incr in linseg_gen([0., 1., 1., 0.], nframes):
        try:
            res = {
                'line': (lambda: incr), 
                'sine': (lambda: sin(.5 * pi * incr)),
                'hann': (lambda: sin(.5 * pi * incr)**2)
                }[shape]()

            yield res

        except KeyError:
            print "Invalid shape: %s" % shape
            break


# def win_gen(nframes = 9, warp = 0., shape = 'hann'):
#     """win_gen(nframes = 9, warp = 0., shape = 'hann')

#     Like win(), but instead of returning an array, returns an object that
#     generates the numbers in the range on demand.  For looping, this is 
#     slightly faster than win() and more memory efficient.
    
#     Return a window array of nframes.
#     Shapes: sine, hann, bart
#     Warp: warps increment by incr**(2**warp)

#     Endpoints are zeros. (Or nearly!)
    
#     """

#     for incr_prewarped in xlinspace(0., 1., nframes):
#         incr = incr_prewarped**(2**warp)
#         try:
#             res = {
#                 'sine': (lambda: sin(pi * incr)),
#                 'hann': (lambda: sin(pi * incr)**2),
#                 'bart': (lambda: (1. - 2. / (nframes  - 1.) * abs(.5 * (nframes  - 1.) * (2. * incr - 1.)))), 
#                 }[shape]()

#             yield res

#         except KeyError:
#             print "Invalid shape: %s" % shape
#             break
def win_gen(nframes = 9, warp =.5, shape = 'hann'):
    """win_gen(nframes = 9, warp = .5, shape = 'hann')

    Like win(), but instead of returning an array, returns an object that
    generates the numbers in the range on demand.  For looping, this is 
    slightly faster than win() and more memory efficient.
    
    Return a window array of nframes.
    Shapes: sine, hann, bart
    Warp: warps increment by Wn_warp

    Endpoints are zeros. (Or nearly!)
    
    """

    for incr_prewarped in xlinspace(0., 1., nframes):
        incr = Wn_warp(incr_prewarped, 1.-warp, Wn_w = True)

        try:
            res = {
                'sine': (lambda: sin(pi * incr)),
                'hann': (lambda: sin(pi * incr)**2),
                'bart': (lambda: (1. - 2. / (nframes  - 1.) * abs(.5 * (nframes  - 1.) * (2. * incr - 1.)))), 
                }[shape]()

            yield res

        except KeyError:
            print "Invalid shape: %s" % shape
            break


def rnd_linseg_gen(vals = [0., 1.], npoints = 3, nframes = 9):
    """rnd_linseg_gen(vals = [0., 1.], npoints = 3, nframes = 9)

    Like rnd_linseg(), but instead of returning an array, returns an object that
    generates the numbers in the range on demand.  For looping, this is 
    slightly faster than linseg and more memory efficient.    Return evenly spaced numbers.
    
    Return nframes of evenly spaced samples between vals [low, high],
    with npoints breakpoints.     """

    # generate breakpoints
    rnframes = unique(random_integers(1, nframes - 1, npoints - 2))

    if rnframes[-1] != nframes:
        rnframes = append(rnframes, nframes)


    # generate values at breakpoints
    rvals = uniform(vals[0], vals[1], len(rnframes) + 1)

    return linseg_gen(rvals, rnframes)


# # noise. . .

# def white_gen(nframes = 9):
#     """white_gen(nframes = 9)

#     Like white_ar(), but instead of returning an array, returns an object that
#     generates samples on demand.  For looping, this is 
#     slightly faster than white_ar() and more memory efficient.
    
#     Return nframes of uniformly distributed random numbers.

#     (A convenience function for random.uniform)
    
#     nframes may be a tuple.

#     """

#     if isscalar(nframes):
#         nframes = (nframes, 1)

#     for x in xrange(nframes[0]):
#         yield uniform(-1., 1., nframes[1:])


# fixed osc. . .
# variable osc is found in filters.py
# as phase incrementor needs state

def fphasor_gen(Wn = .5, phase = 0., nframes = 9):
    """fphasor_gen(Wn = .5, phase = 0., nframes = 9)

    Like fphasor(), but instead of returning an array, returns an object that
    generates the numbers in the range on demand.  For looping, this is 
    slightly faster than fphasor() and more memory efficient.
    """
    if not (isscalar(Wn) and isscalar(phase)):
        Wn = asarray(Wn)
        phase = asarray(phase)

    for n in xrange(nframes):
        yield (n * pi * Wn + phase)


def fsinosc_gen(Wn = .5, phase = 0., nframes = 9):
    """fsinosc_gen(Wn = .5, phase = 0., nframes = 9)

    Like fsinosc(), but instead of returning an array, returns an object that
    generates the numbers in the range on demand.  For looping, this is 
    slightly faster than fsinosc() and more memory efficient.
    """
    for x in fphasor_gen(Wn, phase, nframes):
        yield sin(x)


# ******************************************************************************
# 
# Classes. . .
# 
# **********************

# superclass MusObj
# this needs to handle the various environment
# set-ups along with multichannel expansion
class MusObj:

    sr = 44100                  # default sample rate
    T = 1./sr                   # default sample period

    def __init__(self, sr = None):

        if sr is None:          # assign sampling rate and period
            self.sr = MusObj.sr
            self.T = MusObj.T
        else:
            self.sr = sr
            self.T = 1./sr
        
        self.nchannels = None


# ******************************************************************************
# generate objects: generate signals from control values
# three ways to return array
#    .ar() returns complete array
#    .ar(nframes) returns nframes from array (incrementing)
#    .fr() returns one frame from array (incrementing)
class GenObj(MusObj):

    def __init__(self, vals, durs, sr, gen = True, asnframes = False):
        MusObj.__init__(self, sr) # Run superclass init

        self.vals = asarray(vals) # store vals

        if asnframes:           # convert and set nframes
            nframes = durs
        else:
            nframes = dur_to_nframes(durs, self.sr)

        if isscalar(nframes):
            self.nframes = nframes
        else:
            self.nframes = asarray(nframes)

        self.nchannels = nchannels(asarray(vals)) # store nchannels

        if gen is True:
            self._gen_set()            # initialize / reset generator


    def _gen_set(self):      # set up generator!

        if self.nchannels is 1:
            self.gen = self.it(self.vals, self.nframes)

        else:
            res = array([])

            for nchan_vals in deinterleave(self.vals):
                res = append(res, self.it(nchan_vals, self.nframes))
                self.gen = it.izip(*tuple(res))


    def reset(self):            # reset the generator!
        self._gen_set()


    def ar(self, nframes = None): # array return

        if nframes is None:       # return the complete array

            if self.nchannels is 1:
                res = self.fun(self.vals, self.nframes)

            else:               # multichannel
                res = array([])

                for nchan_vals in deinterleave(self.vals):
                    res = append(res, self.fun(nchan_vals, self.nframes))
                res.shape = (self.nchannels, -1)

            return interleave(res)

        else:                     # or return nframes, by using the generator
            return self._gr(nframes)


    def _gr(self, nframes = None): # array return
        if nframes is None:       # return the complete array
            return self.ar()

        else:                     # or return nframes, by using the generator
            res = array([])

            for i in range(nframes):
                try:
                    res = append(res, self.gen.next())

                except StopIteration:
                    pass

            if len(res) > 0:
                return reshape(res, (-1, self.nchannels))


    def fr(self): # frame return
        try:
            res = asarray(self.gen.next())

            if self.nchannels is 1:
                res.shape = (1,)
            return res

        except StopIteration:
            pass


# envelopes. . .

class Lin(GenObj):

    def __init__(self, vals = [0., 1.], dur = 1., sr = None, gen = True, asnframes = False):
        self.fun = lin          # define function
        self.it = lin_gen       # define iterator
        GenObj.__init__(self, vals, dur, sr, gen, asnframes) # . . . then run superclass init


class LinSeg(GenObj):

    def __init__(self, vals = [0., 1., 0.], durs = [.5, 1.], sr = None, gen = True, asnframes = False):
        self.fun = linseg       # define function
        self.it = linseg_gen    # define iterator
        GenObj.__init__(self, vals, durs, sr, gen, asnframes) # . . . then run superclass init


class Exp(GenObj):

    def __init__(self, vals = [1., 2.], dur = 1., sr = None, gen = True, asnframes = False):
        self.fun = exp          # define function
        self.it = exp_gen       # define iterator
        GenObj.__init__(self, vals, dur, sr, gen, asnframes) # . . . then run superclass init


class ExpSeg(GenObj):

    def __init__(self, vals = [1., 2., 1.], durs = [.5, 1.], sr = None, gen = True, asnframes = False):
        self.fun = expseg       # define function
        self.it = expseg_gen    # define iterator
        GenObj.__init__(self, vals, durs, sr, gen, asnframes) # . . . then run superclass init


class Linen(GenObj):

    def __init__(self, durs = [.25, .75, 1.], shape = 'line', sr = None, gen = True, asnframes = False):
        self.fun = linen        # define function
        self.it = linen_gen     # define iterator
        self.shape = shape
        GenObj.__init__(self, None, durs, sr, gen, asnframes) # . . . then run superclass init

    def _gen_set(self):      # set up generator!
        self.gen = self.it(self.nframes, self.shape)


    def ar(self, nframes = None): # array return
        if nframes is None:       # return the complete array
            return interleave(self.fun(self.nframes, self.shape))
        else:
            return self._gr(nframes)


class Win(GenObj):

    def __init__(self, dur = 1., warp = 0., shape = 'hann', sr = None, gen = True, asnframes = False):
        self.fun = win          # define function
        self.it = win_gen       # define iterator
        self.warp = warp
        self.shape = shape
        GenObj.__init__(self, None, dur, sr, gen, asnframes) # . . . then run superclass init

    def _gen_set(self):         # set up generator!
        self.gen = self.it(self.nframes, self.warp, self.shape)

    def ar(self, nframes = None): # array return
        if nframes is None:       # return the complete array
            return interleave(self.fun(self.nframes, self.warp, self.shape))
        else:
            return self._gr(nframes)


# possible GenObj additions:
#    other forms:  rand envelope generation


# white noise is a special case. . . no need for generator
class WhiteNoise(GenObj):       # if dur is array, generates multichannel output [dur, nchannels]

    def __init__(self, dur = 1., sr = None, gen = False, asnframes = False):
        self.fun = white     # define function
        GenObj.__init__(self, None, dur, sr, gen, asnframes) # . . . then run superclass init
        self.nchannels = (isscalar(dur) and 1) or dur[-1]    # then reset nchannels

    def ar(self, nframes = None): # array return
        if nframes is None:       # return the complete array
            nframes = (isscalar(self.nframes) and self.nframes) or self.nframes[0]
        return self.fun([nframes, self.nchannels])

    def fr(self):               # array return
        return self.fun(self.nchannels)


# fixed sin oscil. . .

class FSinOsc(GenObj):

    def __init__(self, freq = 440., phase = 0., dur = 1., sr = None, gen = True, asnframes = False):
        self.fun = fsinosc      # define function
        self.it = fsinosc_gen   # define iterator
        self.freq = asarray(freq)
        self.phase = asarray(phase)
        GenObj.__init__(self, None, dur, sr, gen, asnframes) # . . . then run superclass init
        self.nchannels = max(                                # then reset nchannels
            (isscalar(freq) and 1) or len(freq),
            (isscalar(phase) and 1) or len(phase))

    def _gen_set(self):         # set up generator!
        self.gen = self.it(freq_to_Wn(self.freq, self.T), self.phase, self.nframes)

    def ar(self, nframes = None): # array return
        if nframes is None:       # return the complete array

            res = self.fun(freq_to_Wn(self.freq, self.T), self.phase, self.nframes)

            if self.nchannels is 1:
                res = interleave(res)
            return res
        else:
            return self._gr(nframes)
