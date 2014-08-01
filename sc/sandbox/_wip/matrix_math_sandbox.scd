x = Matrix.newIdentity(4)

y = x + (Matrix.newIdentity(4) * -2)

z = y * Matrix.newClear(4,1).putCol(0, [1,1,1,1])

(
var processPoints, aeds;
// soundfield
x = Matrix.newIdentity(4);
// subtract an xform from it
y = x + (FoaXformerMatrix.newDominateX(12).matrix * -1);

// point at 0,0 through the filter matrix
z = y * FoaEncoderMatrix.newDirection(0, 0).matrix.insertRow(3,[0]);
// point at pi,0 through the filter matrix
q = y * FoaEncoderMatrix.newDirection(pi, 0).matrix.insertRow(3,[0]);

// aeds.collect({|aed|aed[1]}).postln.do(_.postln);
// aeds = aeds.collect({|aed|[aed[0], aed[1].abs.ampdb]});


// based on Muse AED function plugins
processPoints = {|xFormedMatrix|
	var w,x,y,z, az, el, g0square, g1square, term, omni, gain;
	var b, pv_mean, b_sqrd_mean, p_sqrd, v_sqrd, a, e, d, d_norm;

	b = xFormedMatrix.getCol(0);
	b[0] = b[0] * sqrt(2);			// scale W

	b_sqrd_mean = b**2;

	p_sqrd = b_sqrd_mean[0]; 		// W component, pressure?
	v_sqrd = b_sqrd_mean[1..].sum;	// summed XYZ, velocity?

	// directivity
	d = (pi/2) - (2 * atan2(sqrt(v_sqrd), sqrt(p_sqrd)));

	pv_mean = b * b[0];				// mul all WXYZ by W ( [w,x,y,z] * (W*sqrt(2)) )
	// azimuth
	a = atan2(pv_mean[2], pv_mean[1]); // atan2(y,x) - see pv_mean above
	// elevation - omitted as it's not used currently
	// e = atan2(pv_mean[3], sqrt((pv_mean[1]**2) + (pv_mean[2]**2))); // atan2(z,sqrt(x^2 + y^2))

	d_norm = 1 - (d.abs / 0.5pi);
	// d_norm = (d.abs / 0.5pi);

	postf("w,x,y,z: %\naz: %\nd_norm: %\n", b.round(0.01), a.round(0.01), d_norm);

	/* -------- old code ----*/
	#w,x,y,z = xFormedMatrix.getCol(0);
	az  = atan2(y, x);
	el  = atan2(z, sqrt(x.squared + y.squared));
	// Omni transform angle
	g0square	= w.squared;
	g1square	= x.squared + y.squared + z.squared;
	term 		= ((2*g0square) + g1square);
	(term == 0).if{term = 0.0000000001}; // protect from NaN
	omni		= asin(( (2*g0square) - g1square) / term);

	// Gain
	term		= sin(omni);
	(term <= -1).if{term = -0.99999999}; // protect from NaN
	gain		= w * sqrt(2 / (1 + term));

	// DOES A NEGATIVE GAIN MEAN direction is maintained but phase is inverted?
	// or direction is mirrored across the origin?

	//fold the gain if it's negative ?????
	// if(gain<0, {gain = 1+gain.abs});

	// // to retain the sign of the gain
	// a = [-1,2,-3,4];
	// b = a**2;
	// c = b * (a.abs/a); // how to protect from /0?


	// Normalise omni to [0,1] range
	omni = 2 * omni / pi;
	omni = 1 - omni;
	omni = omni.clip(0.0,1.0);

	// ignore elev for now, implement with Spherical to include elev
	// omni is normalized and scales rho (distance from center)
	// postf("w: %\nx: %\ny: %\nz: %\naz: %\nel: %\nomni: %\ngain:: %\n", b[0],b[1],b[2],b[3], a, e, d, gain); // debug
	// [
	// 	Polar(omni, az), gain.abs.ampdb
	// 	// [gain.abs, 0.00000001].maxItem.ampdb // protect from -inf
	// ].postln; // debug

	[
		Polar(d_norm, a), gain //.abs.ampdb
		// [gain.abs, 0.00000001].maxItem.ampdb // protect from -inf
	].postln; // debug

};

processPoints.(z);
processPoints.(q);
nil
)