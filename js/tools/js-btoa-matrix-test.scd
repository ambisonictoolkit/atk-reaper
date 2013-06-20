// simple test to verify js b-to-a matricies

(
var orientations = ['flu', 'fld', 'flr', 'fud', 'fbd', 'fbu', 'flru', 'flrd'];
var weights = ['dec', 'can', 'uns', 'car'];
var b_format, a_format;

b_format = Matrix.with([[0.13], [-0.17], [0.29], [-0.37]]); // Trond's DC

orientations.do({arg orientation;
	('-- ' ++ orientation).postln;

	weights.do({arg weight;
		// calculate coefficients
		a_format = (FoaDecoderMatrix.newBtoA(orientation, weight).matrix * b_format);

		('---- ' ++ weight ++ ' ').post;
		a_format.asArray.flat.postln;
	});
});
)