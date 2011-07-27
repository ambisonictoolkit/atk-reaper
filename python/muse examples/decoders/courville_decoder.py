# **************************************
# 
# Courville example decodes: ST350, AKT Blue Line, ZoomH2
#
# Decode the three example recordings made by
# Courville [1] using the ST350, 3 AKG Blue Line
# microphones, and the ZoomH2 recorder.
#
# These are decoded to a dual-band compressed 'quad' array
# (+/-30deg) with nearfield (distance) compensation.
# (Alternative decoders can also be using this model.) This example
# uses 'block' decoding--and keeps track of the states of the varying
# decoder filters (nearfield distance compensation, psychoacoustic
# shelf), and is the method is suitable for decoding large files.
#
# Additionally, the use of python to call external applications
# from the command line is also illustrated. These are:
#
#   Flac [2]                : to de-compress soundfiles
#   libsamplerate [3]       : to resample soundfiles
#
#
# To run these examples, download:
#
#   http://www.radio.uqam.ca/ambisonic/audio/comparative_quad.zip
#   http://www.radio.uqam.ca/ambisonic/audio/comparative_b-format.zip
#
# Decompress these files and place in a folder named 'courville'.
# You'll need to point the var named 'courville_dir' to this directory.
#
# NOTE: sndfile-resample writes out files as WAVE_FORMAT_EXTENSIBLE
#       'wav' format, which is unreadable by the Python Standard Library
#       file i/o. This means Muse running with audiolab installed is
#       required to correctly run this file.
# 
#
# [1] D. Courville, "Comparative Surround Recording," Ambisonic Studio |
# Comparative Surround Recording, 2007. [Online].
# Available: http://www.radio.uqam.ca/ambisonic/comparative_recording.html.
# [Accessed: 26-Jul-2011].
#
# [2] J. Coalson, "FLAC - Free Lossless Audio Codec." [Online].
# Available: http://flac.sourceforge.net/. [Accessed: 26-Jul-2011].
#
# [3] E. de Castro Lopo, "Secret Rabbit Code (aka libsamplerate)."
# [Online]. Available: http://www.mega-nerd.com/SRC/. [Accessed: 26-Jul-2011].
#
# J Anderson, 2011
#
# **************************************

from muse import *
import os


# params
exec_dir = '/usr/local/bin'  # directory where flac and src are installed

sr = 44100          # sample rate (output)

Wn_ps       = freq_to_Wn(400, 1./sr) #psychoacoustic shelf crossover
nf_r        = 1.2           # NFC distance, in meters (for compensation)

num_chans = 4               # decode to 4 speakers == quad
angle = deg_to_rad(60/2)    # 1/2 angle for 'narrow' quad decode
                            # front speakers are at +- 30 deg

gains       = [1, 1, -1.3]  # in dB, ['st350, 'akg', 'h2']

a_h2        = 0.58          # mic patern for H2
k_h2        = db_to_amp(4.8) # Y scale for H2
                            # NOTE: Y scaling is equivalent
                            #       to adjusting the front and/or back
                            #       mic angles of the H2. As the H2 patterns
                            #       are not consistent across freq, it is
                            #       easiest to just scale Y 'to taste' after
                            #       b-format encoding

file_type   = 'wav'         # write file...
encoding    = 'pcm24'
endianness  = 'file'

courville_dir =   '/Volumes/Audio/courville/'   # point this to the
                                        # dir where you've decompressed
                                        # the source folders
mic_names = ['st350', 'akg', 'h2']

blocksize = 88200       # number frames to read/write at once


# --------------------------------------------------
# constants, dicts
mic_type = {
    'st350' :   'b-format',     # associated microphone type
    'akg'   :   'b-format',
    'h2'    :   'quad'
    }

# --------------------------------------------------
for i in range(len(mic_names)):

    max_db = -inf                   # max db for file, initialise

    mic_name = mic_names[i]
    scale = db_to_amp(gains[i])

    # --------------------------------------------------
    # paths and files

    # generate flac (source) file name
    flac_file = courville_dir + 'comparative_' + mic_type[mic_name] + '/' + \
                mic_name + '_' + mic_type[mic_name] + '.flac'

    # create temp dir if it doesn't exist
    temp_dir = courville_dir + 'temp/'
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    # decompressed and resampled wav (source) file names
    wav_file = temp_dir + mic_name + '_' + mic_type[mic_name] + '.' + file_type
    rss_file = temp_dir + mic_name + '_' + mic_type[mic_name] + \
              '_SR' + str(sr) + '.' + file_type

    # decoded, output file
    dec_file = courville_dir + 'decode' + '/' + mic_name + '_decode' + \
               '.' + file_type

    # create decode dir if it doesn't exist
    if not os.path.exists(os.path.dirname(dec_file)):
        os.makedirs(os.path.dirname(dec_file))


    # --------------------------------------------------
    # decompress source (flac)
    if not os.path.exists(wav_file):
        print '\nDecompressing flac...'
        os.system(exec_dir + '/flac -d ' + flac_file + ' -o ' + wav_file)
        print '     ... decompressed!'

    # --------------------------------------------------
    # resample decompressed source (wav)
    if not os.path.exists(rss_file):
        print '\nResampling to ' + str(sr) + '...'
        os.system(exec_dir + '/sndfile-resample -to ' + str(sr) + \
                  ' ' + wav_file + ' ' + rss_file)
        print '     ... resampled!'


    # --------------------------------------------------
    # set up sndfiles for reading and writing:
    rss_sndfile = Sndfile(rss_file)             # read resampled file

    dec_sndfile = Sndfile(                      # write to decode file
        dec_file,
        'w',
        Format(file_type, encoding, endianness),
        num_chans,
        sr
        )


    # set up list with blocksize blocks to iterate over
    numframes = rss_sndfile.nframes
    blocks = [ blocksize ] * (numframes / blocksize)

    if (numframes % blocksize) != 0:
        blocks = blocks + [ numframes % blocksize ]

    print '\nDecoding:', mic_name

    # --------------------------------------------------
    # set up filter states for decoding (psycho_shelf, nfc):
    zi_ps = zeros((4,2))
    zi_nf = zeros((4,1))

    # --------------------------------------------------
    # decode loop
    for (block_num, block) in zip(range(len(blocks)), blocks):

        # read in signal
        in_sig = rss_sndfile.read_frames(block)

        # is it proper b-format?
        if mic_name is 'st350':             # yes!
            b_sig = in_sig              
        elif mic_name is 'akg':             # no... add zeros as Z
            b_sig = hstack((
                in_sig,
                zeros((nframes(in_sig), 1))
                ))
        elif mic_name is 'h2':              # ZoomH2 to b-format
            in_sig = hstack((               # re-order channels as required
                in_sig[:, 2:],              # by zoomH2_to_b()
                in_sig[:, :2]
                ))
            b_sig = zoomH2_to_b(            # encode to b-format
                in_sig,
                a = a_h2,                   # a to match observed st350
                k = k_h2                    # k to match observed st350
                )
            b_sig = mirror_x(b_sig)         # mirror across x-axis, as
                                            # Courville's H2 is reversed
                                            # from zoomH2_to_b() convention

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
        dec_sndfile.write_frames(b_dec)


    # ----- close files
    rss_sndfile.close()
    dec_sndfile.close()


    # ----- final report

    print "\nMax db = ", max_db
    print "Done!!"
