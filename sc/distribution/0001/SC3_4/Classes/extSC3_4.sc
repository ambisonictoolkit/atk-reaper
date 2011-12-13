/*
SuperCollider Licensing

SuperCollider is copyright © James McCartney and many other contributors.

SuperCollider is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
*/
/*
	This file is distributed as part of SuperCollider3 (3.4) version of
	the Ambisonic Toolkit (ATK).
	
	The SuperCollider3 version of the Ambisonic Toolkit (ATK) is free software:
	you can redistribute it and/or modify it under the terms of the GNU General
	Public License as published by the Free Software Foundation, either version 3
	of the License, or (at your option) any later version.
	
	The SuperCollider3 version of the Ambisonic Toolkit (ATK) is distributed in
	the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
	implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
	the GNU General Public License for more details.
	
	You should have received a copy of the GNU General Public License along with the
	SuperCollider3 version of the Ambisonic Toolkit (ATK). If not, see
	<http://www.gnu.org/licenses/>.
*/


//---------------------------------------------------------------------
//	The Ambisonic Toolkit (ATK) is a soundfield kernel support library.
//
// 	Extension: Dictionary, String
//
//	The Ambisonic Toolkit (ATK) is intended to bring together a number of tools and
//	methods for working with Ambisonic surround sound. The intention is for the toolset
//	to be both ergonomic and comprehensive, providing both classic and novel algorithms
//	to creatively manipulate and synthesise complex Ambisonic soundfields.
//	
//	The tools are framed for the user to think in terms of the soundfield kernel. By
//	this, it is meant the ATK addresses the holistic problem of creatively controlling a
//	complete soundfield, allowing and encouraging the composer to think beyond the placement
//	of sounds in a sound-space and instead attend to the impression and image of a soundfield.
//	This approach takes advantage of the model the Ambisonic technology presents, and is
//	viewed to be the idiomatic mode for working with the Ambisonic technique.
//	
//	
//	We hope you enjoy the ATK!
//	
//	For more information visit http://ambisonictoolkit.net/ or
//	email info[at]ambisonictoolkit.net
//
//---------------------------------------------------------------------

// NOTE: this extension to Dictionary is drawn from the SuperCollider 3.5 code base
// NOTE: these extensions to String are drawn from the SuperCollider 3.5 code base


+ Dictionary {
	
	merge {|that, func, fill = true|
		var commonKeys, myKeys = this.keys, otherKeys = that.keys;
		var res = ();

		if (myKeys == otherKeys) {
			commonKeys = myKeys
		} {
			commonKeys = myKeys.sect(otherKeys);
		};

		commonKeys.do { |key|
			res[key] = func.value(this[key], that[key], key)
		};

		if (fill) {
			myKeys.difference(otherKeys).do { |key| res[key] = this[key] };
			otherKeys.difference(myKeys).do { |key| res[key] = that[key] };
		};
		^res
	}
}

+ String {

	wrapExtend { arg size;
		^this.dup(size div: this.size).join
	}

	padLeft { arg size, string = " ";
		^string.wrapExtend(max(0, size - this.size)) ++ this
	}
}