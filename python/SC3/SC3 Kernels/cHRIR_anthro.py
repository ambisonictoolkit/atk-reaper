# CIPIC Database HRIR model
#
# http://interface.cipic.ucdavis.edu/sound/hrtf.html

from muse import *
import pylab as pl
import scipy.io as sio


# consts
duda_r = 0.0875     # duda model standard radius

#a = .5              # model slope (good!)
a = .66              # model slope (good!)
#a = .75              # model slope (good!)
#a = .85              # model slope

# cHRIR dir
database_dir = '/Users/josephla/Documents/Developer/CIPIC_hrtf_database/anthropometry/'

#-----------------------------------------------------------
# generate lHRIR path
anthro_file = database_dir + 'anthro.mat'

# Read matlab file:
anthro_matfile = sio.loadmat(anthro_file)

whd = zeros((len(anthro_matfile['X']), 3))     # head width, height, depth

for i in range(len(anthro_matfile['X'])):
    whd[i] = anthro_matfile['X'][i][:3]
    whd[i] *= .005                             # scale to meters, radius

whd = whd[logical_not(isnan(whd)[:,0])]         # discard nans
wd = delete(whd, 1, 1)
r = mean(wd, 1)                                 # radius = mean(w, d)

var = (amax(r) - amin(r)) / ( 2 * median(r))    # +/- % variation from
                                                # the median
med_r = median(r)                               # median of r

std_r = repeat(duda_r, len(r))                  # duda's standard radius

mod_r = (var * (                                # model to match duda's rs
    (1-a) * lin([-1, 1], len(r)) + \
    a * (lin([-1, 1], len(r)))**3) + 1) \
    * med_r

# print out values                                USE THESE for model!!
print 'a = ', a                                 # 0.66
print 'var = ', var                             # 0.123494002733
print 'med_r = ', med_r                         # 0.0864243630344


# display
pl.plot(lin([0, 1], len(r)), sort(whd, 0))  # sorted, w, h, d

pl.plot(lin([0, 1], len(r)), sort(r))       # sorted r

pl.plot(lin([0, 1], len(r)), std_r)         # duda model r
pl.plot(lin([0, 1], len(r)), lin([med_r, med_r], len(r)))

pl.plot(lin([0, 1], len(r)), mod_r)         # 'new' model r


pl.show()
pl.close()
