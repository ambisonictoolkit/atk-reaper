//---------------------------------------------------------------------
// Joseph Anderson 2011
//
// extension of SequenceableCollection to support spherical and cartesian (3d) coordinates
//
// NOTE: include license, etc before public release
//---------------------------------------------------------------------

+ SequenceableCollection {

	// unary math ops
	phi { ^this.performUnaryOp('phi') }

	// math op dispatch support
//	performBinaryOpOnComplex { arg aSelector, aComplex, adverb;
//		^this.collect({ arg item;
//			aComplex.perform(aSelector, item, adverb)
//		})
//	}

	// conversion
	asCartesian { ^Cartesian(this[0] ? 0, this[1] ? 0, this[2] ? 0) }


}
