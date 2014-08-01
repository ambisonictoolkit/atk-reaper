(
// Todo:
// --add gain labels on every other point
// add makeup gain slider
// add matrix subtraction
// consider annular wedge representations
// --add circle on original 0,0 point
// --post min/max gain
// implement units for tranform (pi, db, receiver index)
// add azimuth controls to the transforms
// add option to display directivity as expanding/transparent circle or just proximity to origin

// constants
var kInvSqrt2 = 1/(2.sqrt);
var numPoints = 24, initPointsMatrices, transformedPoints, aeds;
var controls;
// var xFormDuples = [];
// functions
var processPoints, updateMatrix, getXForm, chainXForms, createXForm;
var createNewXFormView, removeXForm;

// drawing vars
var scrnB, winW, winH, uv, ctlv, codev, pnts, arcH;
var degreeSpec, degreeSl, degreeNB, degreePiNB, xFormMenu, xForms;
var rotSpec, rotSl, rotNB, rotPiNB;
var selDex, adjustRhos, getColCount, layoutFunc, colors;

var tv, maxGainTxt, minGainTxt;

// // based on ATK-JS plugins
// processPoints = {|xFormedMatrix|
// 	var w,x,y,z, az, el, g0square, g1square, term, omni, gain;
//
// 	#w,x,y,z = xFormedMatrix.getCol(0);
// 	az  = atan2(y, x);
// 	el  = atan2(z, sqrt(x.squared + y.squared));
//
// 	// Omni transform angle
// 	g0square	= w.squared;
// 	g1square	= x.squared + y.squared + z.squared;
// 	term 		= ((2*g0square) + g1square);
// 	(term == 0).if{term = 0.0000000001}; // protect from NaN
// 	omni		= asin(( (2*g0square) - g1square) / term);
//
// 	// Gain
// 	term		= sin(omni);
// 	(term <= -1).if{term = -0.99999999}; // protect from NaN
// 	gain		= w * sqrt(2 / (1 + term));
//
//
//
// 	// Normalise omni to [0,1] range
// 	omni = 2 * omni / pi;
// 	omni = 1 - omni;
// 	omni = omni.clip(0.0,1.0);
//
// 	// ignore elev for now, implement with Spherical to include elev
// 	// omni is normalized and scales rho (distance from center)
// 	postf("w: %\nx: %\ny: %\nz: %\naz: %\nel: %\nomni: %\ngain:: %\n", w,x,y,z, az, el, omni, gain); // debug
// 	[
// 		Polar(omni, az), gain.abs.ampdb
// 		// [gain.abs, 0.00000001].maxItem.ampdb // protect from -inf
// 	].postln; // debug
// };

// based on Muse AED function plugins
processPoints = {|xFormedMatrix|
	var w,x,y,z, az, el, g0square, g1square, term, omni, gain;
	var b, pv_mean, b_sqrd_mean, p_sqrd, v_sqrd, a, e, d, d_norm;

	b = xFormedMatrix.getCol(0);
	b[0] = b[0] * sqrt(2); // scale W

	b_sqrd_mean = b**2;

	p_sqrd = b_sqrd_mean[0]; // W component, pressure?
	v_sqrd = b_sqrd_mean[1..].sum; // summed XYZ, velocity?

	d = (pi/2) - (2 * atan2(sqrt(v_sqrd), sqrt(p_sqrd)));

	pv_mean = b * b[0]; // mul all WXYZ by W ( [w,x,y,z] * (W*sqrt(2)) )
	a = atan2(pv_mean[2], pv_mean[1]); // atan2(y,x) - see pv_mean above
	e = atan2(pv_mean[3], sqrt((pv_mean[1]**2) + (pv_mean[2]**2))); // atan2(z,sqrt(x^2 + y^2))

	d_norm = 1 - (d.abs / 0.5pi);
	// d_norm = (d.abs / 0.5pi);

	postf("w,x,y,z: %\naz, el: %\nd_norm: %\n", b.round(0.01), [a, e].round(0.01), d_norm);

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


initPointsMatrices = numPoints.collect{|i|
	var ang;
	ang = (2pi/numPoints) * i;
	// Matrix.with([[kInvSqrt2, cos(ang), sin(ang), 0]]);
	FoaEncoderMatrix.newDirection(ang).matrix.addRow([0]); // need to add z manually
};

chainXForms = {|...xFormDegPrs|
	var mtx;
	"in chainXForms".postln;
	xFormDegPrs.postln;
	xFormDegPrs.do{ |xNameDeg|
		var name, deg;
		#name, deg = xNameDeg;

		(name != 'subtract').if({
			// tranform via matrix multiplication
			mtx.isNil.if(
				{ mtx = xForms[name]['getMatrix'].(deg) }, // first mtx
				{ mtx = xForms[name]['getMatrix'].(deg) * mtx; }
			);
			},{ // subtract
				mtx = Matrix.newIdentity(4) + (mtx * -1);
				// mtx = mtx/2;
		});
	};
	// mtx * inputMtx;
	// mtx * (Matrix.newIdentity(4) * -1) * inputMtx
	mtx
};

updateMatrix = {
	var xFormDuples = [];

	"in updateMatrix".postln;
	controls.do{|ctlDict|
		ctlDict.xform.notNil.if{
			ctlDict.degree.isNil.if{"trasform degree is nil!".error};
			xFormDuples = xFormDuples.add( [ctlDict.xform, ctlDict.degree] );
		};
	};

	"found these duples: ".post; xFormDuples.postln;

	transformedPoints = (xFormDuples.size > 0).if({
		initPointsMatrices.collect{ |pointMtx|
			var xformMtx;
			// xformMtx = chainXForms.(pointMtx, *xFormDuples);
			// xformMtx * pointMtx;
			xformMtx = chainXForms.(*xFormDuples);
			xformMtx * pointMtx
		};
		},{ initPointsMatrices } // no transform
	);

	// transformedPoints = initPointsMatrices.collect{ |pointMtx|
	// 	(Matrix.newIdentity(4) * -2) * pointMtx;
	// };

	// calculate and set aeds var for gui update
	aeds = transformedPoints.collect{|ptMtx|
		"\n******-------------".postln;
		processPoints.(ptMtx);
	};
	aeds.collect({|aed|aed[1]}).postln.do(_.postln);
	aeds = aeds.collect({|aed|[aed[0], aed[1].abs.ampdb]});
};

updateMatrix.(); // initialize

/***************************************
Drawing
***************************************/
scrnB = Window.screenBounds;
winW= 600;
winH= 600;
w = Window("soundfield transform", Rect(scrnB.center.x - (winW/2), scrnB.center.y - (winH/2), winW, winH), resizable: true).front;
// w.layout_(VLayout(
uv = UserView( w, Rect(0,0, w.view.bounds.width, w.view.bounds.height/2) ).resize_(5).background_(Color.gray.alpha_(0.99));
ctlv = View( w, Rect(0,w.view.bounds.height/2, w.view.bounds.width, w.view.bounds.height/2) ).resize_(5);
codev = View();
// ));
// uv = UserView();
// ctlv = View();

// controls
// xForms = ['push', 'press', 'focus', 'zoom'];
xForms = IdentityDictionary(know: true).putPairs([
	'push',		IdentityDictionary(know: true).putPairs(
		[ 'min', -pi/2, 'max', pi/2, 'default', 0,
			'getMatrix', {|deg| FoaXformerMatrix.newPushX(deg).matrix} ]
	),
	'press',	IdentityDictionary(know: true).putPairs(
		[ 'min', -pi/2, 'max', pi/2, 'default', 0,
			'getMatrix', {|deg| FoaXformerMatrix.newPressX(deg).matrix} ]),
	'focus',	IdentityDictionary(know: true).putPairs(
		[ 'min', -pi/2, 'max', pi/2, 'default', 0,
			'getMatrix', {|deg| FoaXformerMatrix.newFocusX(deg).matrix} ]),
	'zoom',		IdentityDictionary(know: true).putPairs(
		[ 'min', -pi/2, 'max', pi/2, 'default', 0,
			'getMatrix', {|deg| FoaXformerMatrix.newZoomX(deg).matrix} ]),
	'dominate',		IdentityDictionary(know: true).putPairs(
		[ 'min', -24, 'max', 24, 'default', 0,
			'getMatrix', {|gain| FoaXformerMatrix.newDominateX(gain).matrix} ]),
	'direct',		IdentityDictionary(know: true).putPairs(
		[ 'min', 0, 'max', pi/2, 'default', 0,
			'getMatrix', {|deg| FoaXformerMatrix.newDirectX(deg).matrix} ]),
	'rotate',	IdentityDictionary(know: true).putPairs(
		[ 'min', 2pi,	'max', -2pi, 'default', 0,
			'getMatrix', {|deg| FoaXformerMatrix.newRotate(deg).matrix} ]),
	'asymmetry',	IdentityDictionary(know: true).putPairs(
		[ 'min', -pi/2,	'max', pi/2, 'default', 0,
			'getMatrix', {|deg| FoaXformerMatrix.newAsymmetry(deg).matrix} ]),
	'balance',	IdentityDictionary(know: true).putPairs(
		[ 'min', -pi/2,	'max', pi/2, 'default', 0,
			'getMatrix', {|deg| FoaXformerMatrix.newBalance(deg).matrix} ]),
	'gain',	IdentityDictionary(know: true).putPairs(
		[ 'min', -24,	'max', 24, 'default', 0,
			'getMatrix', {|gainDB| gainDB.dbamp } ]),
	'subtract', IdentityDictionary(know: true).putPairs(
		[ 'min', 0,	'max', 0, 'default', 0,
			'getMatrix', {|deg| FoaXformerMatrix.newBalance(deg)} ]),
]);

// order determines sf processing chain
// list can grow with inserted items and shrink when removed
controls = List();

createXForm = {
	var updateSpecs, xFormMenu, degSpec, degSl, degNB, degPiNB, removeBut, addBut;
	var dict, view, layout, deg;

	"in createXForm".postln;

	view = View();
	// store xform, degree, view, layout here
	dict = IdentityDictionary(know: true);

	updateSpecs = { |xf| // xf is a dict from xForms
		var min, max;
		"in updatespecs".postln;
		degSpec = ControlSpec(xf.min, xf.max, default: xf.default);
		// .min/maxItem to handle negative max on rotate
		min = [xf.min, xf.max].minItem;
		max = [xf.min, xf.max].maxItem;
		degNB.clipHi_(max).clipLo_(min);
		degPiNB.clipHi_(max).clipLo_(min);
		// initialize widgets
		deg.notNil.if({
			if( (deg<min) or: (deg>max), { deg = degSpec.default });
			},{ deg = degSpec.default }
		);
		degSl.value   = degSpec.unmap(deg);
		degNB.value   = deg;
		degPiNB.value = deg / pi;
	};

	degSpec = ControlSpec();

	xFormMenu = PopUpMenu().items_(['-']++xForms.keys)
	.action_({ |mn|
		if( mn.item != '-', {
			var xfName, xfAttributes;
			xfName = mn.item;
			xfAttributes = xForms[xfName];
			updateSpecs.(xfAttributes);
			// update control dict
			dict.xform_(xfName).degree_(deg);
			updateMatrix.();
			w.refresh;
			},{ // nil mutes the transform
				dict.xform_(nil).degree_(nil);
				updateMatrix.();
				w.refresh;
		});
	}).maxWidth_(75).value_(0);

	degSl = Slider()
	.action_(
		{ |sl| var val;
			val = degSpec.map(sl.value);
			deg = val;
			degNB.value_(val.round(0.001));
			degPiNB.value_(val.round(0.001)/pi);
			dict.degree_(deg);
			updateMatrix.(deg);
			w.refresh;
		}
	)
	.step_(0.001).orientation_('horizontal');

	degNB = NumberBox()
	.action_(
		{ |nb| var val;
			val = nb.value;
			deg = val;
			degSl.value_(degSpec.unmap(val));
			degPiNB.value_(val / pi);
			dict.degree_(deg);
			updateMatrix.(deg);
			w.refresh;
		}
	)
	.step_(0.01).maxWidth_(60);

	degPiNB = NumberBox()
	.action_(
		{ |nb| var val;
			val = nb.value * pi;
			deg = val;
			degSl.value_(degSpec.unmap(val));
			degNB.value_(val);
			dict.degree_(deg);
			updateMatrix.(deg);
			w.refresh;
		}
	)
	.step_(0.01/pi).maxWidth_(60);

	addBut = Button().states_([["+"]])
	.action_({ |but|
		createNewXFormView.(view);
	})
	.maxWidth_(20);

	removeBut = Button().states_([["X"]])
	.action_({ |but|
		removeXForm.(dict);
		updateMatrix.(deg);
		w.refresh;
	})
	.maxWidth_(20);

	layout = HLayout(
		xFormMenu, degNB,
		degSl,
		degPiNB, StaticText().string_("pi").align_('left'),
		removeBut, addBut
	).margins_([10,0,10,3]);

	view.layout_(layout);

	dict.view = view;
	dict.layout = layout;
	dict; // return
};

createNewXFormView = { |behindThisView|
	var behindDex, dict;
	dict = createXForm.();
	behindThisView.notNil.if({
		controls.do{|ctldict, i|
			(ctldict.view === behindThisView).if{behindDex = i};
		};
		behindDex.isNil.if{"preceeding view not found!".error};
		"behindDex: ".post; behindDex.postln;
		controls.insert(behindDex+1, dict);
		ctlv.layout.insert(dict.view, behindDex+2); // +2 to account for label row in ctlv[0]
		},{
			controls.insert(0, dict);
			ctlv.layout.insert(dict.view, 0);
	});
	// return
	dict.view;
};

removeXForm = {|dict| var rmvDex;
	controls.do{|ctldict, i|
		(ctldict === dict).if{rmvDex = i};
	};
	rmvDex.isNil.if{"view not found!".error};
	dict.view.remove;
	dict.layout.destroy;
	controls.removeAt(rmvDex);
};


// ctlv.layout_(
// 	VLayout(
// 		createNewXFormView.(),
// 		codev.layout_(VLayout(TextView().enterInterpretsSelection_(true)))
// 	)
// );

ctlv.layout_( VLayout().margins_(0) );
ctlv.layout.add(
	HLayout(
		StaticText().string_("Transform").maxWidth_(60), 20,
		StaticText().string_("Degree of Tranform").minWidth_(150), 280,
		HLayout(
			StaticText().string_("remove").minWidth_(50),
			StaticText().string_("add").maxWidth_(30) //, a: \right],
		)
	)
	// .setAlignment(1, \right).setAlignment(2, \right).setAlignment(3, \right)
	.margins_([10, 5,5,0])
);
ctlv.layout.add(createNewXFormView.());
ctlv.layout.add( codev.layout_(VLayout(TextView().enterInterpretsSelection_(true))) );

// // tv = View(uv, Rect(3*uv.bounds.width/4, 0, uv.bounds.width/4, uv.bounds.height/2)).layout_(
// tv = View(uv, Rect(0, 0, uv.bounds.width/4, uv.bounds.height/2)).layout_(
// 	VLayout(
// 		Slider().maxWidth_(25),
// 		Button().maxWidth_(25),
// 		Slider().maxWidth_(25)
// 	).setAlignment(1, 'right')
// );
// // w.view.onResize_({|winview| tv.bounds_(Rect(3*uv.bounds.width/4, 0, uv.bounds.width/4, uv.bounds.height/2))});

// tv = View(w, uv.bounds);
tv = View(uv,Rect(0, uv.bounds.height/2, uv.bounds.width,uv.bounds.height/2)).layout_(
	HLayout(
		VLayout(
			nil,
			maxGainTxt = StaticText().string_("one").font_(Font('Helvetica', 16)).maxWidth_(160).minHeight_(21),
			minGainTxt = StaticText().string_("two").font_(Font('Helvetica', 16)).maxWidth_(160).minHeight_(21),
		),
		VLayout(
			[ Button()
				.states_([
					["Pan Test", Color.black, Color.grey],
					["Stop test", Color.black, Color.red],
			]).maxWidth_(65), a: \right],
			Button().states_([
				["Continuous", Color.black, Color.grey],["Pulsed", Color.black, Color.green]]),
			HLayout(
				StaticText().string_("Period"),
				NumberBox().maxWidth_(25).value_(5)
			).setAlignment(0,\right).setAlignment(1,\right),
			VLayout(
				HLayout(
					StaticText().string_("Gain").maxHeight_(25),
					NumberBox().maxWidth_(25)
				).setAlignment(0, \right).setAlignment(1, \right),
				// Slider().maxWidth_(150).orientation_(\horizontal)
			).setAlignment(0,\right) //.setAlignment(1,\right) //.setAlignment(0,\right)
		)
		// .setAlignment(0, 'right')
		.setAlignment(1, 'right') // continuous/pulses buttons
		.setAlignment(2, 'right')
		.setAlignment(3, 'right')
	)
);

w.view.onResize_({ |winview|
	uv.bounds_(Rect(0,0, w.view.bounds.width, w.view.bounds.height/2));
	tv.bounds_(Rect(0, uv.bounds.height/2, uv.bounds.width,uv.bounds.height/2));
});

/*******************
Draw the soundfield
*******************/

uv.drawFunc_({|view|
	var cen, getColor, circleViewRatio;
	var alphaSpec, gainSpec, colorSpec, maxMinStr;
	var testSigPnt;

	ctlv.bounds_(Rect(0,w.view.bounds.height/2, w.view.bounds.width, w.view.bounds.height/2));
	// uv.bounds_(Rect(0,0, w.view.bounds.width, w.view.bounds.height/2));
	// tv.bounds_(uv.bounds);
	r = uv.bounds.height * 0.02;
	d = r*2;
	circleViewRatio = 0.8;
	arcH = uv.bounds.height * circleViewRatio / 2;	// height of the "fan" arc

	// center drawing origin
	cen = view.bounds.center;
	Pen.translate(cen.x, cen.y);
	// draw background "fan"
	Pen.strokeColor_(Color.red).fillColor_(Color.blue);
	Pen.fillOval(Rect(r.neg, r.neg, d,d));
	Pen.strokeOval(Rect(r.neg, r.neg, d,d));
	Pen.addAnnularWedge( 0@0, 5, arcH, 0, 2pi );
	Pen.fillColor_(Color.gray(0.9)).fill;

	alphaSpec = ControlSpec(0.1, 1, warp:3);
	// gainSpec = ControlSpec(-100, 6, 'db');
		gainSpec = ControlSpec(-90, 6, -2);
	colorSpec = ControlSpec(768,0);

	// postf("min gain: %\n", aeds.collect({|me| me[1]}).minItem);
	// postf("max gain: %\n", aeds.collect({|me| me[1]}).maxItem);
	// aeds.do(_.postln);

		// draw test signal planewave circle
	testSigPnt = Polar(1.1,rrand(0, 2pi)).asPoint
		.rotate(pi/2)	// convert ambi to screen coords
		* Point(1,-1)	// flip Y for drawing
		* arcH;			// scale normalized points to arcH
	Pen.fillColor_(Color.magenta);
	Pen.fillOval( Rect(testSigPnt.x-r, testSigPnt.y-r, d, d) );

	getColor = { |gain|
		var i;
		i = colorSpec.map(gainSpec.unmap(gain));
		case
		{i < 256} {Color.new255(255, i, 0)}
		{i < 384} {Color.new255(255 - (i-256), 255, 0)}
		{i < 512} {Color.new255(0, 255, (i-384)*2)}
		{i < 768} {Color.new255(0, 255 - (i-512), 255)}
		{i >= 768} {Color.new255(0, 0, 255)}; // catch all
	};

	aeds.do{|pntGainArr, i|
		var drawPnt, polar, gain, omniRad, omniDiam, fullOmni, gainColor, gainLabelColor, gainPnt;
		#polar, gain = pntGainArr;
		fullOmni = 2 * arcH;
		omniDiam = 1-polar.rho * fullOmni;
		omniDiam = omniDiam.clip(d, fullOmni);
		omniRad= omniDiam/2;

		gainColor = getColor.(gain);
		// gainLabelColor = getColor.(gain + 3);
		// cartesian point in view coordinates
		drawPnt = polar.asPoint
		.rotate(pi/2)	// convert ambi to screen coords
		* Point(1,-1)	// flip Y for drawing
		* arcH;			// scale normalized points to arcH
		// original 0deg azimuth point circle
		if( i==0, {
			Pen.strokeColor_(Color.fromHexString("#CC0000"));
			Pen.strokeOval( Rect(drawPnt.x-r, drawPnt.y-r, d, d) );
		});
		// line fron center to point
		Pen.strokeColor_(Color.gray.alpha_(0.5));
		Pen.line(drawPnt, 0@0).stroke;
		// directivity circle
		Pen.fillColor_(gainColor.alpha_(alphaSpec.map(polar.rho)));
		Pen.fillOval( Rect(drawPnt.x-omniRad, drawPnt.y-omniRad, omniDiam, omniDiam) );
		// gain labels
		gainPnt = polar.rho = 1;
		gainPnt = gainPnt.asPoint.rotate(pi/2)
		// flip y for screen bounds and
		// scale in/out toward/away from origin, >1 outside circle
		// * Point(0.4,-0.4)
		* Point(1.15,-1.15)
		* arcH;
		Pen.fillColor_(gainColor.alpha_(1));
		QPen.stringCenteredIn(gain.round(0.1).asString, Rect(gainPnt.x-(r*10), gainPnt.y-(r*10), d*10, d*10));
	};

	// // draw min/max gain
	// maxMinStr = format("Max Gain: %\nMin Gain: %",
	// 	aeds.collect({|me| me[1]}).maxItem.round(0.1),
	// 	aeds.collect({|me| me[1]}).minItem.round(0.1)
	// );
	// Pen.translate(cen.x.neg, cen.y.neg);
	// Pen.fillColor_(Color.black);
	// Pen.font_(Font('Helvetica', 20));
	// QPen.stringCenteredIn(
	// 	maxMinStr,
	// 	Rect(0,0,uv.bounds.width/3, ctlv.bounds.height)
	// );

	maxGainTxt.string_(format("Max Gain:  %", aeds.collect({|me| me[1]}).maxItem.round(0.1)));
	minGainTxt.string_(format("Min Gain:   %", aeds.collect({|me| me[1]}).minItem.round(0.1)));
});

w.refresh;
)



/********************************************************************
Matrix: Interprete processed points

In the resulting matrixOut:
- Row 1 is x position of transformed point
- Row 2 is y position of transformed point
- Row 3 is z position of transformed point
- Row 4 is gain value in dB for transformed point
********************************************************************/
function interpreteProcessedPoints(aMatrixOut, aMatrixIn, aNumPoints) local(i, lG0square, lG1square, lW, lX, lY, lZ, lAzi, lEle, lOmni, lGgain)
(
  i = 0;
  loop(aNumPoints,
    lW = aMatrixIn[i];
    lX = aMatrixIn[aNumPoints + i];
    lY = aMatrixIn[aNumPoints*2 + i];
    lZ = aMatrixIn[aNumPoints*3 + i];

    lAzi  = atan2(lY, lX);
    lEle  = atan2(lZ, sqrt(lX*lX + lY*lY));

    // Omni transform angle
    lG0square        = lW * lW;
    lG1square        = lX*lX + lY*lY + lZ*lZ;
    lOmni = asin((2.*lG0square - lG1square) / (2.*lG0square + lG1square));

    // Gain
    lGain = lW * sqrt(2. / (1. + sin(lOmni)));

    // Normalise omni to [0,1] range
    lOmni = 2 * lOmni / $pi;
    lOmni = 1. - lOmni;

    // From this we calculate relative xyz position of the transformed point
    aMatrixOut[i]                = lOmni * cos(lAzi) * cos(lEle);
    aMatrixOut[aNumPoints   + i] = lOmni * sin(lAzi) * cos(lEle);
    aMatrixOut[aNumPoints*2 + i] = lOmni * sin(lEle);
    aMatrixOut[aNumPoints*3 + i] = amp2db(lGain);

    i += 1;
  );
);


/********************************************************************
Matrix: Generate 8x4 matrix with equally distributed source signals for display purposes
********************************************************************/
function generateDisplaySignalMatrix(aMatrix, aNumPoints) local(i, lAngle, lInc)
(
  // Encoding equally distributed points in the horisontal plane
  i = 0;
  lAngle = 0.;
  lInc = 2. * $pi / aNumPoints;
  loop(aNumPoints,
    aMatrix[i]                 = kInvSqrt2;
    aMatrix[i + aNumPoints]    = cos(lAngle);
    aMatrix[i + 2*aNumPoints]  = sin(lAngle);
    aMatrix[i + 3*aNumPoints]  = 0.;

    lAngle += lInc;
    i += 1;
  );
);

/**********************************************************
Focus/Press/Push/Zoom
**********************************************************/
azimuth   = 0;
elevation = 0.;
mDrawCenterPoint = 0;   // This depends on whether transform mode affects directioness or not

// Setting memory locations for a number of 4x4 matrixes as 1D arrays
matrixRotate1   = 1000;
matrixTumble1   = 1020;
matrixTransform = 1040;
matrixTumble2   = 1060;
matrixRotate2   = 1080;

// These are used for matrix multiplications
matrixTemp1     = 1100;
matrixTemp2     = 1120;

// Resulting matrixes used for DSP processing
matrixDSP       = 1200;
matrixNewDSP    = 1220;
matrixInc       = 1240;

// Matrixes for graphic display of transformations
matrixDisplayPoints      = 1300;
matrixDisplayTransformed = 1400;
matrixDisplayInterpreted = 1500;

// Matrixes for testing and debugging of matrixMultiplication()
matrixCounter = 2000;
matrixTest    = 2100;
matrixUnity   = 2200;

// Generate test matrixes
i = 0;
loop(
  32,
  matrixCounter[i] = i * 0.01;
  i += 1;
);

i = 0;
loop(
  16,
  matrixUnity[i] = 0.;
  i += 1;
);

matrixUnity[0] = 1.;
matrixUnity[5] = 1.;
matrixUnity[10] = 1;
matrixUnity[15] = 1.;

matrixMultiplication  (matrixTest,    4, 4, matrixUnity,     4, 8, matrixCounter);

// This only needs to be generated once:
numDisplayedPoints = 24;
generateDisplaySignalMatrix(matrixDisplayPoints, numDisplayedPoints);

function calculateMatrix()
(
  // Required for graphics
  cosAzimuth   = cos(azimuth);
  sinAzimuth   = sin(azimuth);

  cosElevation = cos(elevation);
  sinElevation = sin(elevation);

  // Generate matrixes
  generateRotationMatrix(matrixRotate1, (-azimuth));
  generateTumblelMatrix(matrixTumble1, -elevation);

  // Determine what transform to do
  mode == 0 ? (         // Focus
    generateFocusMatrix(matrixTransform, transformAngle);
    mDrawCenterPoint = 0;
  ) : (
    mode == 1 ? (       // Press
      generatePressMatrix(matrixTransform, transformAngle);
      mDrawCenterPoint = 1;
    ) : (
      mode == 2? (      // Push
        generatePushMatrix(matrixTransform, transformAngle);
        mDrawCenterPoint = 1;
      ) : (             // Zoom
        generateZoomMatrix(matrixTransform, transformAngle);
        mDrawCenterPoint = 0;
      );
    );
  );
  generateRotationMatrix(matrixRotate2, azimuth);
  generateTumblelMatrix(matrixTumble2, elevation);

  /* Build processing matrix by performing a series of matrix multiplications
                        result =         left X                   right             */
  matrixMultiplication  (matrixTemp1,    4, 4, matrixTumble1,     4, 4, matrixRotate1);
  matrixMultiplication  (matrixTemp2,    4, 4, matrixTransform,   4, 4, matrixTemp1);
  matrixMultiplication  (matrixTemp1,    4, 4, matrixTumble2,     4, 4, matrixTemp2);
  matrixMultiplication  (matrixNewDSP,   4, 4, matrixRotate2,     4, 4, matrixTemp1);

  // Transformation and interpretation of display points
// matrixNewDSP is the matrix resulting from the transform
  matrixMultiplication  (matrixDisplayTransformed, 4, 4, matrixNewDSP, 4, numDisplayedPoints, matrixDisplayPoints);
  interpreteProcessedPoints(matrixDisplayInterpreted, matrixDisplayTransformed, numDisplayedPoints);

);


/********************************************************************
Graphics: Draw points illustrating the effect of various transforms
********************************************************************/
function displayTransformedPoints(aMatrix, aNumPoints, aDrawCenterPoint) local(i, lPosX, lPosY, lPosZ, lDirectedness, lSquareSize)
(
  // Draw transparent circles indicating the location and directedness of each source after the transform
  // TODO: Also take z value into account in this visualisation?
  i = 0;
  loop(
    aNumPoints,
    lPosX = gCenterX - gRadius * aMatrix[aNumPoints+i];
    lPosY = gCenterY - gRadius * aMatrix[i];
    lPosZ = aMatrix[2 * aNumPoints + i];

    // Calculate degree of directedness
    lDirectedness = sqrt(aMatrix[i]*aMatrix[i] + aMatrix[aNumPoints+i]*aMatrix[aNumPoints+i] + aMatrix[2*aNumPoints+i]*aMatrix[2*aNumPoints+i]);

	// Size, saturation and luminocity depends on posZ (mtm - i.e. directness)
    lSquareSize = 2. * lPosZ + 8.;
    gainToColor(aMatrix[3 * aNumPoints + i], -24., 6., lPosZ);

    //drawCenteredSquare(lPosX, lPosY, lSquareSize);
    drawSourceLocation(lPosX, lPosY, lPosZ, lDirectedness);
    i += 1;
  );

  // Draw darker circles at the azi/ele/directedness center location of each source after the transform
  ( aDrawCenterPoint ?
    (
      i = 0;
      loop(
        aNumPoints,
        lPosX = gCenterX - gRadius * aMatrix[aNumPoints+i];
        lPosY = gCenterY - gRadius * aMatrix[i];
        gxColorSet(0.25, 0., 0., 0.);
        drawSourceLocationCenter(lPosX, lPosY);
        i += 1;
      );
    );
  );
);

/**********************************************************
Graphics
**********************************************************/


  // This is the regular interactive GUI mode:

  // Set colors
  gxColorSet(1, 1, 1, 1);
  determineGeometry();
  drawBackgroundImage(800, 800, gCenterX - gRadius, gCenterY - gRadius, 2 * gRadius, 2 * gRadius);

  // Reset azimuth and transformAngle on mouse click
  mouse_cap == 1 ?
  (
    mouseX = mouse_x - gCenterX;
    mouseY = mouse_y - gCenterY;

    // Update azimuth and slider2 and notify GUI of the change
    azimuth = atan2(mouseX, -mouseY);
    slider2 = azimuth * kRad2Deg;
    slider_automate(slider2);

    // Update transformAngle and slider4 and notify GUI of the change
    // COMMENT: There is a risk of dividing by zero here (for elevation = 90), but Reaper seems to handle this grazefully.
    // At 90 deg elevation the GUI has a dead lock and can't be used for changing the degree of transformation, but the slider is still available.
    mNormTransformAngle = sqrt(mouseX*mouseX + mouseY*mouseY) * gInverseRadius / cos(elevation);
    // Restrict the range to [0, 1]
    mNormTransformAngle > 1. ? mNormTransformAngle = 1.;

    transformAngle = mNormTransformAngle * $pi * 0.5;
    slider4 = mNormTransformAngle * 90.;
    slider_automate(slider4);

    calculateMatrix();
  );

  // Draw points illustrating the effect of the transform
  displayTransformedPoints(matrixDisplayInterpreted, numDisplayedPoints, mDrawCenterPoint);

  // Draw circle for interaction, map vertical position to saturation, lightness and radius
  (sinElevation < 0) ? (
    mSaturation = 1. - 0.25*sinElevation;
    mLightness = 0.4;
  ) : (
    mSaturation = 1.;
    mLightness = 0.4 + 0.25*sinElevation;
  );
  gxColorSetHsl(1.0, 195., mSaturation, mLightness);

  posX = gCenterX + gRadius * mNormTransformAngle * (sinAzimuth * cosElevation);
  posY = gCenterY - gRadius * mNormTransformAngle * (cosAzimuth * cosElevation);
  mInteractionCircleRadius = 4. * sinElevation + 11.;
  drawInteractionCircle(posX, posY, mInteractionCircleRadius);

/**********************************************************
--Graphics
**********************************************************/

// From Muse
# aed (azimuth, elevation, directivity)
# return constants
def aed(a):
    """aed(a)

    Analyze an ambisonic B-format sound field, returning azimuth,
    elevation and directivity.

    Inputs:
        - b         : Input b-format signal

    Outputs: ([a, e, d])

      [a, e, d] -- Azimuth, elevation and directivity in radians.
                   (See direct for details on directivity.)

    """

    // normalise W
b = copy(a)
b[:, 0] *= sqrt(2) // normalize all Ws by sqrt(2)

// pv & b**2 mean

// grab w, multiply input array (W,X,Y,Z) all by w
pv_mean = mean(b * interleave(b[:, 0]), 0) // mul all WXYZ by (W*sqrt(2)) then take average (if input is an array of soundfields)?
b_sqrd_mean = mean((b**2), 0) // square all components and average them (if input is an array of soundfields)

// # p**2 and v**2
p_sqrd = b_sqrd_mean[0] // W component, pressure?
v_sqrd = sum(b_sqrd_mean[1:]) // summed XYZ, velocity?

// # calculate azimuth, elevation
a = arctan2(pv_mean[2], pv_mean[1]) // atan2(y,x) - see pv_mean above
e = arctan2(pv_mean[3], sqrt((pv_mean[1])**2 + (pv_mean[2])**2)) // atan2(z,sqrt(x^2 + y^2))

// # calculate directivity
// # pi/2 - 2 * arctan(v/p)
d = pi/2 - 2 * arctan2(sqrt(v_sqrd), sqrt(p_sqrd))

// # return [azimuth, elevation, directivity]
    res = array([a, e, d])

    return res
