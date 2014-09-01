This is the Ambisonic Toolkit (ATK) as a JSFX plugin suite for Reaper.
It can be used with Reaper on OSX and Windows, and can also be used with
other VST hosts on the Windows platform using the ReaJS plugin. ReaJS is
part of the [ReaPlugs](http://reaper.fm/reaplugs/index.php) plugin suite.

The Ambisonic Toolkit (ATK) is intended to bring together a number of
tools and methods for working with Ambisonic surround sound. The
intention is for the toolset to be both ergonomic and comprehensive,
providing both classic and novel algorithms to creatively manipulate and
synthesise complex Ambisonic soundfields.

The tools are framed for the user to think in terms of the soundfield
kernel. By this, it is meant the ATK addresses the holistic problem of
creatively controlling a complete soundfield, allowing and encouraging
the composer to think beyond the placement of sounds in a sound-space
and instead attend to the impression and image of a soundfield. This
approach takes advantage of the model the Ambisonic technology presents,
and is viewed to be the idiomatic mode for working with the Ambisonic
technique.

We hope you enjoy the ATK!

For more information please visit the [Ambisonic Toolkit
website](http:ambisonictoolkit.net/) or send us an
[e-mail](mailto:info[at]ambisonictoolkit.net).

&nbsp;

&nbsp;

Installing
==========

&nbsp;

Requirements
------------

* ATK for Reaper requires [Reaper 4.7 or above](http://reaper.fm).

&nbsp;

Windows
-------

1. Start Reaper.
2. From the Options menu choose "Show REAPER resource path in explorer/finder".
3. Unzip the ATK for Reaper archive.
4. Follow the instructions suggested by folder names.

&nbsp;

OSX: What Gets Installed Where?
-------------------------------

All files are installed into the following two folders in your home library folder:

    ~/Library/Application Support/ATK
    ~/Library/Application Support/Reaper/Effects/ATK

&nbsp;

Additionally an alias (or rather a symlink) is created at:

    ~/Library/Application Support/Reaper/Data/ATK

&nbsp;

that points to:

    ~/Library/Application Support/ATK

&nbsp;

If you are using Ambisonic Toolkit with SuperCollider as well,
the convolution kernels are installed in the same place and have 
the exact same content. We do not expect this to cause any conflicts.

If you want to take a look at the installed files and do not see the 
Library folder in Finder, please press the ALT button while clicking
the "Go" menu in Finder. The Library folder will show up as an 
additional option.

&nbsp;

Need Some Sound Files to Play Around With?
------------------------------------------

You can find a collection of sound files here:

* [http://www.ambisonictoolkit.net/wiki/tiki-index.php?page=Downloads](http://www.ambisonictoolkit.net/wiki/tiki-index.php?page=Downloads)

&nbsp;

Additional sound files can be grabbed from these fine sources:

* [http://ambisonia.com/](http://ambisonia.com/).
* [http://www.freesound.org/browse/tags/B-format/](http://www.freesound.org/browse/tags/B-format/).
* [http://www.surround-library.com/](http://www.surround-library.com/) (commercial library ambisonic sound effects).
* [http://www.spheric-collection.com/](http://www.spheric-collection.com/) (commercial library ambisonic sound effects).

&nbsp;

And most of the catalogue of Nimbus Records are UHJ recordings:

* [http://www.wyastone.co.uk/](http://www.wyastone.co.uk/).

&nbsp;

&nbsp;

Feedback and Bug Reports
========================

&nbsp;


Known Issues and Limitatons:
----------------------------

**Encoders:**

* Diffuser encoder: Only works if the Reaper project is set to one of 
  the following sample rates: 44100, 48000, 88200, 96000 or 192000 Hz.
* Spreader encoder: Only works if the Reaper project is set to one of
  the following sample rates: 44100, 48000, 88200, 96000 or 192000 Hz.
* SuperStereo encoder: Only works if the Reaper project is set to one
  of the following sample rates: 44100, 48000, 88200, 96000 or 192000 Hz.
* UHJ encoder: Only works if the Reaper project is set to one of the
  following sample rates: 44100, 48000, 88200, 96000 or 192000 Hz.
* ZoomH2: Not yet implemented.

&nbsp;

**Transformers:**

(no known issues)

&nbsp;

**Decoders:**

* Binaural decoder: The Cipic and Listen HRTFs only work if the Reaper
  project is set to 44100 Hz sample rate.
* Diametric: Not yet implemented.
* Periphonic: Not yet implemented.
* UHJ decoder: Only works if the Reaper project is set to one of the
  following sample rates: 44100, 48000, 88200, 96000 or 192000 Hz.
* Quad: K and shelf filter remains to be implemented.
* Pantophonic: K and shelf filter remains to be implemented.

&nbsp;

Reporting issues
----------------

For issues pertaining to the ATK for Reaper plugins, please e-mail [trond.lossius@bek.no](mailto:trond.lossius@bek.no).

&nbsp;

Alternatively you can use the JS plugins issue tracker:

* [http://www.ambisonictoolkit.net/wiki/tiki-view_tracker.php?trackerId=6](http://www.ambisonictoolkit.net/wiki/tiki-view_tracker.php?trackerId=6)

&nbsp;

If you end up using the plugins for some project, we would be more
than interested to know! You can either mail us information, or add it
to the wiki at:

* [http://www.ambisonictoolkit.net/wiki/tiki-index.php?page=Examples](http://www.ambisonictoolkit.net/wiki/tiki-index.php?page=Examples)

&nbsp;


List of Changes
---------------

Version 1.0.b2

* Fixed issue where matrix-based transform plugins could blow up when used on 2-channel (stereo) tracks.
* Fixed issue where multi-channel matrix-based decoders could blow up when used on 2-channel (stereo) tracks.
* Fixed issue that prevented Omni encoder from producing sound.
* OSX installer is now distributed as disk image.
* Info for beta testers has been merged into this readme document.
* This Readme file is now versioned and maintained as markdown document, and converted to html by instaler script using [Pandoc](http://johnmacfarlane.net/pandoc/).

&nbsp;

Version 1.0.b1:

* First beta release.

&nbsp;

&nbsp;

Credits
=======

&nbsp;

Copyright the ATK Community, Joseph Anderson, Joshua Parmenter and
Trond Lossius, 2014.

* J Anderson : [[e-mail]](mailto:j.anderson[at]ambisonictoolkit.net)
* J Parmenter : [[e-mail]](mailto:j.parmenter[at]ambisonictoolkit.net)
* T Lossius : [[e-mail]](trond.lossius[at]bek.no)

&nbsp;

The port of ATK as a set of Reaper JS plugins by Trond Lossius is
supported by [BEK, Bergen Centre for Electronic Arts](www.bek.no).

The filter kernels distributed with the Ambisonic Toolkit are licensed
under a Creative Commons Attribution-Share Alike 3.0 Unported License and
are copyright the Ambisonic Toolkit Community and Joseph Anderson,
2011.\
 <http://creativecommons.org/licenses/by-sa/3.0/>

&nbsp;

&nbsp;

Third Party Notices
===================

&nbsp;

Diametric Decoder Theorem (DDT) decoding
----------------------------------------

Support for Gerzon's Diametric Decoder Theorem (DDT) decoding algorithm
is derived from Aaron Heller's Octave code available at:
http:www.ai.sri.com/ajh/ambisonics/

Benjamin, et al., "Localization in Horizontal-Only Ambisonic Systems"
Preprint from AES-121, 10/2006, San Francisco

Implementation in the SuperCollider3 version of the ATK is by [Joseph
Anderson](mailto:j.anderson[at]ambisonictoolkit.net).

&nbsp;

Irregular array decoding
------------------------

Irregular array decoding coefficients (5.0, 7.0) are kindly provided by
Bruce Wiggins: http:www.brucewiggins.co.uk/

B. Wiggins, "An Investigation into the Real-time Manipulation and
Control of Three-dimensional Sound Fields," PhD Thesis, University of
Derby, Derby, 2004.

&nbsp;

CIPIC HRTF Database (University of California)
----------------------------------------------

V. R. Algazi, R. O. Duda, D. M. Thompson, and C. Avendano, "The CIPIC
HRTF Database," in Proceedings of the 2001 IEEE ASSP Workshop on
Applications of Signal Processing to Audio and Acoustics, New Paltz, NY,
2001.

"The CIPIC HRTF Database - CIPIC International Laboratory." [Online].
Available: <http://interface.cipic.ucdavis.edu/sound/hrtf.html>.
[Accessed: 07-Jul-2011].

**CIPIC Notices:**

Copyright (c) 2001 The Regents of the University of California. All
Rights Reserved

Disclaimer

THE REGENTS OF THE UNIVERSITY OF CALIFORNIA MAKE NO REPRESENTATION OR
WARRANTIES WITH RESPECT TO THE CONTENTS HEREOF AND SPECIFICALLY DISCLAIM
ANY IMPLIED WARRANTIES OR MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR
PURPOSE.

Further, the Regents of the University of California reserve the right
to revise this software and/or documentation and to make changes from
time to time in the content hereof without obligation of the Regents of
the University of California to notify any person of such revision or
change.

Use of Materials

The Regents of the University of California hereby grant users
permission to reproduce and/or use materials available therein for any
purpose- educational, research or commercial. However, each reproduction
of any part of the materials must include the copyright notice, if it is
present. In addition, as a courtesy, if these materials are used in
published research, this use should be acknowledged in the publication.
If these materials are used in the development of commercial products,
the Regents of the University of California request that written
acknowledgment of such use be sent to:

CIPIC- Center for Image Processing and Integrated Computing University
of California 1 Shields Avenue Davis, CA 95616-8553

&nbsp;

Listen HRTF Database (IRCAM)
----------------------------

"LISTEN HRTF DATABASE." [Online]. Available:
<http://recherche.ircam.fr/equipes/salles/listen/>. [Accessed:
07-Jul-2011].

**IRCAM Notices:**

Copyright (c) 2002 IRCAM (Institut de Recherche et Coordination
Acoustique/Musique). All Rights Reserved

Use of Materials

The Listen database is public and available for any use. We would
however appreciate an acknowledgment of the database somewhere in the
description of your work (e.g. paper) or in your development.

Contacts:

Olivier Warusfel, Room Acoustics Team, IRCAM 1, place Igor Stravinsky
75004 PARIS, France

&nbsp;

MESA GLU Library
----------------

Code for calculating the inverse of a 4x4 matrix is based on the MESA
implementation of the GLU library: <http://www.mesa3d.org/>

The default Mesa license is as follows:

Copyright (C) 1999-2007 Brian Paul All Rights Reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL BRIAN PAUL BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.

