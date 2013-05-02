"""\
Muse
==========

This is the start of the Muse project, a Music V type language.

"""

# ***** Need to sort out...
#
#       Namespace collisions. At the moment all numpy is imported into muse,
#       which can cause some problems: particularly with numpy, muse, math,
#       and python functions.


# ***** Nice starting messages
print 'Starting up Muse!\n'

try:                                            # try numpy...
    print 'Importing numpy...'
    from numpy import *
    print '                         ... success!'

except ImportError:
    print '     ... oh, no! \nNumpy not found! \Muse depends on numpy!'

else:
    print 'Importing muse...'
    from muse import *                          # import main muse...
    from filters    import *
    from decoders   import *
    from encoders   import *
    from transforms import *
    from generators import *
    from nonlinear  import *
    from analysis   import *
    print '                         ... success!'

    
    try:
        print 'Importing audiolab...'           # try audiolab...
        from scikits.audiolab import *
        print '                         ... success!'

#        print 'Importing sndfiles...'
#        from sndfiles import *        # old soundfiles module to be deprecated
#        print '                         ... success'
#
#        print '                             WARNING: Sndfiles uses the old ' +\
#              'audiolab API, \n' +\
#              '                                      and is soon to be ' +\
#              'deprecated.'

        try:                                    # try samplerate...
            print 'Importing samplerate...'
            from scikits.samplerate import resample as src
            print '                         ... success!'

        except ImportError:
            print '                         ... samplerate unavailable.'

    except ImportError:
        print '                         ... audiolab unavailable!'
        print '                         ... using muse sndfile (PSL), instead.'
        from sndfile import *        # the new PSL aifc/wave interface

        print 'Samplerate is unavailable...'     

    print "\nWelcome to Muse!!"



# ***** Old stuff below...

# from scikits.audiolab import *
# from scikits.samplerate import *
# import scikits.samplerate
# from scikits.samplerate import resample as src

# from sndfiles import *        # the old soundfiles module to be deprecated
#from sndfile import *        # the new PSL aifc/wave interface

#from muse import *
#from filters    import *
#from decoders   import *
#from encoders   import *
#from transforms import *
#from generators import *
#from nonlinear  import *
#from analysis   import *
