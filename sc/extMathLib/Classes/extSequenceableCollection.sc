//---------------------------------------------------------------------
// Joseph Anderson 2011
//
// expansion of SequenceableCollection to support spherical coordinates
//
// NOTE: include license, etc before public release
//---------------------------------------------------------------------

+ SequenceableCollection {
	// unary math ops
//	real { ^this.performUnaryOp('real') }
//	imag { ^this.performUnaryOp('imag') }

//	magnitude { ^this.performUnaryOp('magnitude') }
//	magnitudeApx { ^this.performUnaryOp('magnitudeApx') }
//	phase { ^this.performUnaryOp('phase') }
//	angle { ^this.performUnaryOp('angle') }
//
//	rho { ^this.performUnaryOp('rho') }
//	theta { ^this.performUnaryOp('theta') }

	// math op dispatch support
//	performBinaryOpOnComplex { arg aSelector, aComplex, adverb;
//		^this.collect({ arg item;
//			aComplex.perform(aSelector, item, adverb)
//		})
//	}

	// conversion
//	asPoint { ^Point(this[0] ? 0, this[1] ? 0) }
//	asRect { ^Rect(this[0] ? 0, this[1] ? 0, this[2] ? 0, this[3] ? 0) }
	asCartesian { ^Cartesian(this[0] ? 0, this[1] ? 0, this[2] ? 0) }


}
