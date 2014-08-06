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
#print 'Importing HOA!\n'

print 'Importing HOA functions...'
from muse.hoa import *
from muse.hoa.encoders import *
from muse.hoa.decoders import *
#from muse.hoa.constants import *

#from muse.hoa import *
#
#from encoders import *
#from decoders import *
#from constants import *


print '                         ... success!'
