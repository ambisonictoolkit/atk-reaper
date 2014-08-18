# OSX installer

The OSX installer is created using the [Packages app](http://s.sudre.free.fr/Software/Packages/about.html).

* All files are installed into various subfolders of the current users home directory, and as such this is a per-user install.
* [These instructions](http://s.sudre.free.fr/Software/Packages/Q&A_3.html) were followed in order to ensure that items are installed into the current users home directory.
* At the end of the install process, a script creates the symbolic link from "~/Application Support/Reaper/Data/ATK" to "~/Application Support/ATK". This is _not_ done using `bristow` as suggested [here](http://s.sudre.free.fr/Software/Packages/Q&A_1.html), but rather using the standard `ln` Terminal command.