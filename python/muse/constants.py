#! /usr/bin/env python

# vim:syntax=python

# Joseph Anderson 2007


# Muse
# ==========

# This is the start of the Muse project, a Music V type language.

# Built on pyaudiolab.
# """


# import any needed names. . .
from numpy import *
from math import sqrt


# #=========================
# # Definition of constants
# #=========================

speed_of_sound = 343.           # m/sec @ 20 degC, standard atmosphere

twoPi = 2 * pi
halfPi = .5 * pi
recipTwoPi = 1. / twoPi


# for ATK
sqrt2 = sqrt(2.)

rec_sqrt2 = 1./sqrt(2.)
rec_sqrt6 = 1./sqrt(6.)

sqrt3div2 = sqrt(3.)/2.
sqrt3div6 = sqrt(3.)/6.
sqrt6div3 = sqrt(6.)/3.

sqrt2div_sqrt3 = sqrt(2.)/sqrt(3.)
sqrt3div_sqrt2 = sqrt(3.)/sqrt(2.)

rec_sqrt2sum2 = 1./(sqrt(2.) + 2.)
sqrt2div_sqrt2sum2 = sqrt(2.)/(sqrt(2.) + 2.)
twodiv_sqrt2sum1 = 2./(sqrt(2.) + 1)
twodiv_sqrt2sum_sqrt3 = 2./(sqrt(2.) + sqrt(3.))

sqrt6 = sqrt(6.)

# decoder high shelf gains
k_2D = array([sqrt3div_sqrt2, sqrt3div2])       # pantophonic
k_3D = array([sqrt2, sqrt2div_sqrt3])           # periphonic


# b to a conversion dictionary
b_to_a_dict = {
    'flu' : [[.5, .5, .5, .5],
             [.5, .5, -.5, -.5],
             [.5, -.5, .5, -.5],
             [.5, -.5, -.5, .5]],

    'fld' : [[.5, .5, .5, -.5],
             [.5, .5, -.5, .5],
             [.5, -.5, .5, .5],
             [.5, -.5, -.5, -.5]],

    'flr' : [[.5, .5, rec_sqrt2, 0.],
             [.5, .5, -rec_sqrt2, 0.],
             [.5, -.5, 0., rec_sqrt2],
             [.5, -.5, 0., -rec_sqrt2]],

    'fud' : [[.5, .5, 0., rec_sqrt2],
             [.5, .5, 0., -rec_sqrt2],
             [.5, -.5, rec_sqrt2, 0.],
             [.5, -.5, -rec_sqrt2, 0.]],

    'fbd' : [[.5, sqrt3div2, 0., 0.],
             [.5, -sqrt3div6, 0., -sqrt6div3],
             [.5, -sqrt3div6, rec_sqrt2, rec_sqrt6],
             [.5, -sqrt3div6, -rec_sqrt2, rec_sqrt6]],

    'fbu' : [[.5, sqrt3div2, 0., 0.],
             [.5, -sqrt3div6, 0., sqrt6div3],
             [.5, -sqrt3div6, rec_sqrt2, -rec_sqrt6],
             [.5, -sqrt3div6, -rec_sqrt2, -rec_sqrt6]],

    'flru' : [[.5, sqrt3div6, rec_sqrt2, rec_sqrt6],
              [.5, sqrt3div6, -rec_sqrt2, rec_sqrt6],
              [.5, sqrt3div6, 0., -sqrt6div3],
              [.5, -sqrt3div2, 0., 0.]],

    'flrd' : [[.5, sqrt3div6, rec_sqrt2, -rec_sqrt6],
              [.5, sqrt3div6, -rec_sqrt2, -rec_sqrt6],
              [.5, sqrt3div6, 0., sqrt6div3],
              [.5, -sqrt3div2, 0., 0.]]
}

b_to_a_weight_dict = {
    'can' : [1., 1., 1., 1.],

    'dec' : [sqrt2div_sqrt3, 1., 1., 1.],

    'car' : [sqrt6, 1., 1., 1.],

    'uns' : [sqrt2, 1., 1., 1.]
}

