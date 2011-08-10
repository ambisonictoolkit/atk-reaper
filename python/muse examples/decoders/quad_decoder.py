# **************************************
# 
# Quad example decode: 'Swarm of Echoes', Dan Peterson
#
# Decode one 'long' interleaved file example to quad.
#
# Decoded to a dual-band compressed 'quad' array
# (+/-30deg) with nearfield (distance) compensation.
# (Alternative decoders can also be using this model.) This example
# uses 'block' decoding--and keeps track of the states of the varying
# decoder filters (nearfield distance compensation, psychoacoustic
# shelf), and is the method is suitable for decoding large files.
#
# NOTE: 'Swarm of Echoes' is encoded as FL32 (float32)
#       'aifc' format, which is unreadable by the Python Standard Library
#       file i/o. This means Muse running with audiolab installed is
#       required to correctly run this file.
#
# J Anderson, 2011
#
# **************************************

from muse import *
import os


# params
freq_ps     = 400           #psychoacoustic shelf crossover
nf_r        = 1.2           # NFC distance, in meters (for compensation)

num_chans   = 4             # decode to 4 speakers == quad
angle = deg_to_rad(60/2)    # 1/2 angle for 'narrow' quad decode
                            # front speakers are at +- 30 deg

gain        = -3            # in dB, gain applied after decode


sndfile_dir = '/Volumes/Parking/Dan Peterson/'   # read/write directory
read_file   = 'SwarmofEchoes.aif'

blocksize = 88200       # number frames to read/write at once


# --------------------------------------------------
# paths and files

# decoded, output file
write_dir = sndfile_dir + 'decode' + '/'

root, ext = os.path.splitext(read_file)
write_file = root + '_quad' + ext

# create decode dir if it doesn't exist
if not os.path.exists(write_dir):
    os.makedirs(write_dir)


# --------------------------------------------------
# set up sndfiles for reading and writing:
read_sndfile = Sndfile(sndfile_dir + read_file)  # read file

sr = read_sndfile.samplerate                # set sr

write_sndfile = Sndfile(                    # write file
    write_dir + write_file,
    'w',
    read_sndfile.format,                    # same format as input
    num_chans,
    sr
    )


# --------------------------------------------------
# calcs
scale = db_to_amp(gain)

Wn_ps = freq_to_Wn(400, 1./sr)  #psychoacoustic shelf crossover

max_db = -inf                   # max db for file, initialise




# --------------------------------------------------
# set up list with blocksize blocks to iterate over
numframes = read_sndfile.nframes
blocks = [ blocksize ] * (numframes / blocksize)

if (numframes % blocksize) != 0:
    blocks = blocks + [ numframes % blocksize ]

print '\nDecoding: to Quad!'

# --------------------------------------------------
# set up filter states for decoding (psycho_shelf, nfc):
zi_ps = zeros((4,2))
zi_nf = zeros((4,1))

# --------------------------------------------------
# decode loop
for (block_num, block) in zip(range(len(blocks)), blocks):

    # read in signal
    b_sig = read_sndfile.read_frames(block)

    # nearfield compensate, decode to quad (and scale)
    b_sig, zi_nf = nfc(b_sig, nf_r, 1./sr, zi_nf)
    b_dec, zi_ps = quad_dbd(b_sig, angle, Wn_ps, zi_ps)
    b_dec *= scale

    peak_db = amp_to_db(peak(b_dec))

    print "Decoding %d of %d: block size = %d frames. Peak dB = %f"\
        % ((block_num + 1), len(blocks), block, peak_db)

    if peak_db > max_db:
        max_db = peak_db

    # write to output
    # re-order channels from FL, BL, BR, FR
    #           ... to FL, FR, BL, BR
    b_dec = interleave(
        array([
            b_dec[:, 0],
            b_dec[:, 3],
            b_dec[:, 1],
            b_dec[:, 2]
            ])
        )
    write_sndfile.write_frames(b_dec)


# ----- close files
read_sndfile.close()
write_sndfile.close()


# ----- final report

print "\nMax db = ", max_db
print "Done!!"
