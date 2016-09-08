<div id="container">
  <div id="content">

# ATK for Reaper : Read Me

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
website](http:www.ambisonictoolkit.net/) or send us an
[e-mail](mailto:info[at]ambisonictoolkit.net).

&nbsp;

&nbsp;

# Installing

&nbsp;

## Requirements

* ATK for Reaper requires [Reaper 5.0 or above](http://reaper.fm).

&nbsp;

## Windows

1. Start Reaper.
2. From the Options menu choose "Show REAPER resource path in explorer/finder".
3. Unzip the ATK for Reaper archive.
4. Follow the instructions suggested by the folder names in the unzipped archive. If the Data and Effects folders already contains ATK subfolder from a previous install, these can safely be deleted before copying.

&nbsp;

## Mac OSX

When you run the installer, the required files are installed into the following folders in your home library folder:

    ~/Library/Application Support/ATK
    ~/Library/Application Support/Reaper/ColorThemes
    ~/Library/Application Support/Reaper/Effects/ATK
    ~/Library/Application Support/Reaper/presets

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

## Need Some Sound Files to Play Around With?

You can find a collection of sound files here:

* [http://www.ambisonictoolkit.net/download/recordings/](http://www.ambisonictoolkit.net/download/recordings/)

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

# Feedback and Bug Reports

Known issues are logged at [GitHub](https://github.com/ambisonictoolkit/atk-reaper/issues).

If you experience problems or have questions pertaining to the ATK for Reaper plugins, please create an issue in the [ATK-Reaper issue tracker](https://github.com/ambisonictoolkit/atk-reaper/issues).

If you use the plugins for some project, please [let us know](mailto:info[at]ambisonictoolkit.net). We [plan on](https://github.com/ambisonictoolkit/ambisonictoolkit.github.io/issues/9) adding a gallery of example artistic and creative projects that make use of the Ambisonic Toolkit.

&nbsp;


# List of Changes

## Version 1.0.0.b10

**Changes in this version**

* Added pre-install script to OSX installer that erase earlier versions of the `REAPER/Effects/ATK` folder, preventing possible conflicts when a newer version is installed on top of an older one (issue #47).

## Version 1.0.0.b9

*Released 2016-07-24*

**IMPORTANT: This version breaks backwards compatibility regarding how azimuth angles are described in encoders, tranformsers and decoders.**

**It also breaks backwards compatibility with respect to channel configurations for the Quad Decoder.**

*Quad decoder*

The Quad decoder now returns 4 channels only: L - R - Ls - Rs. Refer to the discussion in issue [#26](https://github.com/ambisonictoolkit/atk-reaper/issues/26) for further details.

*Azimuth angles*

Azimuth angles are now defined to be positive in the counter-clockwise direction. The new azimuth convention aligns with standard ambisonic conventions, and description of azimuths is now consistent between the SuperCollider and Reaper implementations of ATK. 

This is a change from prior beta-versions up to and including v.1.0.0.b7, and unfortunately it breaks backwards compatibility. As an example a planewave encoding done using previous beta verison now will be mirrored between left and right.

When we initially set out to make the Reaper port of ATK, we wanted this new design from the start. However, due to limitations in Reaper version 4 we had to implement azimuth as positive in clockwise rather than counter-clockwise direction. This was disussed in the paper [ATK Reaper: The Ambisonic Toolkit as JSFX plugins](http://www.ambisonictoolkit.net/assets/files/2014-ICMC-ATK-Reaper.pdf) by T. Lossius and J. Anderson, presented at the joint ICMC/SMC conference in Athens in 2014:

> Due to current limitations in the Reaper JSFX API, there are some considerations to be made with respect to description of azimuths. The desired behaviour is that azimuths increase anti-clockwise. At the same time it is also preferable that for sources coming from the front half circle, moving the azimuth slider to the right results in the direction of the sound also moving clockwise to the right, similar to what happens when moving a regular stereo pan pot. The Harpex and Blue Ripple plugins both functions this way. In order to achieve this, horizontal azimuth sliders will need to be implemented with increasing values to the left. This is however currently not supported in Reaper JSFX, and for the time being ATK for Reaper uses azimuth values described according to a navigational coordinate system. If Reaper in the future is updated to support sliders with increasing values to the left, this design decision is likely to be reconsidered, and changed to use acoustic conventions.

REAPER v.5 added support for JSFX sliders with increasing values to the left, and hence this change is now being implemented. This implied changes to the following plugins:

* Changed (breaks backwards compatibility)
    * Encoders
        * Planewave encoder now has positive azimuth in counter-clockwise direction (issue #34)
    * Transforms
        * Direct transform now has positive azimuth in counter-clockwise direction (issue #35)
        * Dominate transform now has positive azimuth in counter-clockwise direction (issue #36)
        * FocusPressPullZoom transform now has positive azimuth in counter-clockwise direction (issue #37)
        * Mirror transform now has positive azimuth in counter-clockwise direction (issue #38)
    * Decoders:
        * Mono decoder now has positive azimuth in counter-clockwise direction (issue #40)

If you need to maintain Reaper projects created using ATK for Reaper beta versions prior to 1.0.0.b8, you can insert a mirror transform 

**Other changes in this version**

* New features
    * New Reaper-theme, provided by Michael McCrea - When selecting the MultiChanMix theme in Reaper, two new Track Mixer Panel Layouts are available with larger multi-channel level meters (issue #23).
    * New encoder: AmbiXtoB - encodes a first order AmbiX format signal to FuMa B-format for further processing in Ambisonic Tooit. Supports 3DN and S3DN.
    * New decoder: BtoAmbiX - decodes a FuMa B-format (used for processing in Ambisonic Tooit) to first order AmbiX format. Supports 3DN and S3DN. S3DN is used to author ambisonic sound files for Google VR 360° videos.
    * Spreader encoder: Now support all sample rates equal to or higher than 44.1 kHz (issues #20, #30).
    * UHJ encoders: Now support all sample rates equal to or higher than 44.1 kHz (issues #20, #31).
    * Binaural decoder: CIPIC and Listen HRTFs now work woth all samle rates equal to or higher than 44.1 kHz (issues #18, #29, #32, #33).
* Refactoring (these changes should not affect the behaviour of plugins):
    * Updated to install atk-kernels v. 1.2.1.
    * Installs atk-matrices, so far used for AmbiX conversions.
    * ATK for Reaper now uses v. 0.3 of the CookDSP library by Olivier Bélanger: https://github.com/belangeo/cookdsp
    * Introducing ATK library header file, simplifies maintenance of libraries and library dependencies (issue #39).
    * Have started to implement matrix library as object-oriented classes, adding new feature to read matrix file from disk. This is used for the AmbiX encoders and decoders.
* On a side note:
    * A new repository has been set up providing [a set of example Ambisonic Toolkit Reaper projects](https://github.com/ambisonictoolkit/atk-reaper-examples). This will soon be available from the ATK web site as a separate download.

## Version 1.0.0.b7
*Released 2016-05-23*

* Refactoring (these changes should not affect the behaviour of plugins):
    * The repository has been moved to [Github](https://github.com/ambisonictoolkit/atk-reaper). If you experience problems with the installers in this version, [please let us know](https://github.com/ambisonictoolkit/atk-reaper/issues).
    * Removed redundant arguments from `drawBackgroundImage()` function - #16.
* New features:
    * Stereo, Quadrophonic, 5_0 and Pantohponic 2D decoders: Added GUI for level metering and display of speaker positions/channel numbering
    * New decoder: Periphonic 3D - 3D dual ring, regular cylindrical decoder, with GUI for level metering and display of speaker positions/channel numbering

## Version 1.0.0.b6
*Released 2015-11-04*

* In Reaper 5 the name of the plugins are retrieved from the desc field rather than filename. The desc field of these plugins was previously used to provide documentation on the plugin. Now it has been changed to have more useful names of the plugins.

## Version 1.0.b5
*Released 2015-08-07*

* Added new utils plugin: MuteSoloChannels - can be used to mute or solo individula channels of a four-channel track. This can be useful when monitoring effect processing of signals within a Bformat => BtoA => SoundFXs => AtoB => B-format signal processing chain.
* The 4channels plugins is now properly installed on Windows.
* Set up the infrastructure required to provide Factory Presets for plugins. A number of presets are provided for the MuteSoloChannels effect.
* Fixed bug in binaural decoder that would prevent the use of Spherical HRTFs at 88.2kHz and 96khz. Thanks to Juan Pampin for pointing this out!
* Resolved: JS plugin automation data generated by GUI is not always properly recorded. This was a bug in Reaper, and was fixed in Reaper 4.75. Issue on the Reaper issue tracker:  http://forum.cockos.com/showthread.php?t=147151

## Version 1.0.b4
*Released 2014-10-04*

* New utility plugin "4channels": Can be used to extract 4-channel A- or B-format recording from sound file with additional channels, e.g. recorded in the field using a portable recorder such as Tascam DR680 or Sound Devices 788T.
* Fixed an issue where matrix-based encoders could cause very loud signals if used on a track with insufficient number of channels.


## Version 1.0.b3
*Released 2014-09-09*

* Fixed issue where multi-channel matrix-based decoders could blow up when used on 2-channel (stereo) tracks.
* Encode/UHJ: Corrected description of plugin in GUI
* Decode/Binaural: Improving text describing in GUI of this plugin
* Decode/5_0: Text description now provides information on output channels
* Decode/Quad: Text description now provides information on output channels

## Version 1.0.b2
*Released 2014-09-3*

* Fixed issue where matrix-based transform plugins could blow up when used on 2-channel (stereo) tracks.
* Fixed issue that prevented Omni encoder from producing sound.
* OSX installer is now distributed as disk image.
* Creation of Zip archive for Windows distribution is automoated using a Terminal script.
* Info for beta testers has been merged into this readme document.
* This Readme file is now versioned and maintained as markdown document, and converted to html by installer script using [Pandoc](http://johnmacfarlane.net/pandoc/). A minimum of CSS is used for it to look OK as a stand-alone HTML document as well as in the OSX installer.

&nbsp;

## Version 1.0.b1:
*Released 2014-08-29*

* First beta release.

&nbsp;

&nbsp;

# Credits

&nbsp;

Copyright the ATK Community, Joseph Anderson, Joshua Parmenter and
Trond Lossius, 2014.

* J Anderson : [[e-mail]](mailto:j.anderson[at]ambisonictoolkit.net)
* J Parmenter : [[e-mail]](mailto:j.parmenter[at]ambisonictoolkit.net)
* T Lossius : [[e-mail]](mailto:trond.lossius[at]bek.no)

&nbsp;

The port of ATK as a set of Reaper JS plugins by Trond Lossius is
supported by [BEK, Bergen Centre for Electronic Arts](www.bek.no).

The filter kernels distributed with the Ambisonic Toolkit are licensed
under a Creative Commons Attribution-Share Alike 3.0 Unported [(CC BY-SA 3.0)](http://creativecommons.org/licenses/by-sa/3.0/) License and
are copyright the Ambisonic Toolkit Community and Joseph Anderson,
2011.

&nbsp;

&nbsp;

# Third Party Notices

&nbsp;

## CookDSP

ATK for Reaper makes use of the [CookDSP library](http://ajaxsoundstudio.com/software/cookdsp/) by Olivier Bélanger.

CookDSP is copyright (c) - 2014 - Olivier Belanger.
Permission to use, copy, modify, and/or distribute this
software for any purpose with or without fee is hereby
granted, provided that the above copyright notice and this
permission notice appear in all copies. CookDSP IS PROVIDED
"AS IS", WITHOUT WARRANTY OF ANY KIND FROM THE AUTHOR.



&nbsp;

## Diametric Decoder Theorem (DDT) decoding

Support for Gerzon's Diametric Decoder Theorem (DDT) decoding algorithm
is derived from Aaron Heller's Octave code available at:
http://www.ai.sri.com/ajh/ambisonics/

Benjamin, et al., "Localization in Horizontal-Only Ambisonic Systems"
Preprint from AES-121, 10/2006, San Francisco

Implementation in the SuperCollider3 version of the ATK is by [Joseph
Anderson](mailto:j.anderson[at]ambisonictoolkit.net).

&nbsp;

## Irregular array decoding

Irregular array decoding coefficients (5.0) are kindly provided by
Bruce Wiggins: http://www.brucewiggins.co.uk/

B. Wiggins, "An Investigation into the Real-time Manipulation and
Control of Three-dimensional Sound Fields," PhD Thesis, University of
Derby, Derby, 2004.

&nbsp;

## CIPIC HRTF Database (University of California)

V. R. Algazi, R. O. Duda, D. M. Thompson, and C. Avendano, "The CIPIC
HRTF Database," in Proceedings of the 2001 IEEE ASSP Workshop on
Applications of Signal Processing to Audio and Acoustics, New Paltz, NY,
2001.

"The CIPIC HRTF Database - CIPIC International Laboratory." [Online].
Available: <http://interface.cipic.ucdavis.edu/sound/hrtf.html>.
[Accessed: 07-Jul-2011].

### CIPIC Notices:

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

## Listen HRTF Database (IRCAM)

"LISTEN HRTF DATABASE." [Online]. Available:
<http://recherche.ircam.fr/equipes/salles/listen/>. [Accessed:
07-Jul-2011].

### IRCAM Notices:

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

## MESA GLU Library

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

  </div> <!-- id="content" -->
</div> <!-- id="container" -->
