(
// Todo:
// --add gain labels on every other point
// --add makeup gain slider
// --add matrix subtraction
// consider annular wedge representations
// --add circle on original 0,0 point
// --post min/max gain
// --implement units for tranform (pi, db, receiver index)
// --add azimuth controls to the transforms
// --add option to display directivity as expanding/transparent circle or just proximity to origin
// consider a view down positive X axis
// --filter out points <-90dB
// --draw most transparent circles first
// --implement matrix input section as well as
// "post trasnform matrix"
// add test signal feature
// add addition
// add indeces in chain for selecting input and receivers for addition/subtraction
// --make matrix text view a pop-out window
// double check how often window needs to be refreshed... maybe just views
// make controls view grow and shrink the window dynamically
// play soundfield through the matrix
// crossfade matrices as parameters are changed

// constants
var kInvSqrt2 = 1/(2.sqrt);
var numPoints = 24, initPointsMatrices, transformedPoints, aeds;
var gControls;
// var xFormDuples = [];
// functions
var processPoints, updateMatrix, getXForm, chainXForms, createXForm;
var createNewXFormView, removeXForm, showMatrixWindow;

// drawing vars
var scrnB, win, winW, winH, uv, ctlv, codev, pnts, arcH, codeWin;
var degreeSpec, degreeSl, degreeNB, degreePiNB, xFormMenu, xForms;
var rotSpec, rotSl, rotNB, rotPiNB;
var selDex, adjustRhos, getColCount, layoutFunc, colors;

var tv, maxGainTxt, minGainTxt;
var getColor, alphaSpec, gainSpec, colorSpec, gainThresh;
var dirDisplay = 'size';
var amp, curXformMatrix, evalTxtView;

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

	// amp
	amp = sqrt((b**2).sum / 2);

	// postf("w,x,y,z: %\naz, el: %\nd_norm: %\n", b.round(0.01), [a, e].round(0.01), d_norm);

	// // TODO: get gain from above code, not the old code
	// /* -------- old code ----*/
	// #w,x,y,z = xFormedMatrix.getCol(0);
	// az  = atan2(y, x);
	// el  = atan2(z, sqrt(x.squared + y.squared));
	// // Omni transform angle
	// g0square	= w.squared;
	// g1square	= x.squared + y.squared + z.squared;
	// term 		= ((2*g0square) + g1square);
	// (term == 0).if{term = 0.0000000001}; // protect from NaN
	// omni		= asin(( (2*g0square) - g1square) / term);
	//
	// // Gain
	// term		= sin(omni);
	// (term <= -1).if{term = -0.99999999}; // protect from NaN
	// gain		= w * sqrt(2 / (1 + term));
	//
	// //debug
	// "AMP AND GAIN:".postln;
	// [amp, gain].postln;

	// // Normalise omni to [0,1] range
	// omni = 2 * omni / pi;
	// omni = 1 - omni;
	// omni = omni.clip(0.0,1.0);

	// ignore elev for now, implement with Spherical to include elev
	[
		// 	Polar(omni, az), gain.abs.ampdb
		// Polar(d_norm, a), gain //.abs.ampdb
		a, e, d_norm, amp
		// [gain.abs, 0.00000001].maxItem.ampdb // protect from -inf
	].postln; // debug

};


initPointsMatrices = numPoints.collect{|i|
	var ang;
	ang = (2pi/numPoints) * i;
	// Matrix.with([[kInvSqrt2, cos(ang), sin(ang), 0]]);
	FoaEncoderMatrix.newDirection(ang).matrix.addRow([0]); // need to add z manually
};

chainXForms = {|...xFormCtlPrs|
	var mtx;
	"in chainXForms".postln;
	xFormCtlPrs.postln;
	xFormCtlPrs.do{ |xNameDeg|
		var name, ctlVals;
		#name, ctlVals = xNameDeg;

		(name != 'subtract').if({
			// tranform via matrix multiplication
			mtx.isNil.if(
				{ mtx = xForms[name]['getMatrix'].(*ctlVals) }, // first mtx
				{ mtx = xForms[name]['getMatrix'].(*ctlVals) * mtx; }
			);
			},{ // subtract
				mtx = Matrix.newIdentity(4) + (mtx * -1);
				// mtx = mtx/2;
		});
	};
	// mtx * inputMtx;
	// mtx * (Matrix.newIdentity(4) * -1) * inputMtx
	curXformMatrix = mtx;
	mtx
};

~updateMatrix = { |xfMatrix|
	"updating from external call".postln;
	transformedPoints = initPointsMatrices.collect{ |pointMtx|
		xfMatrix * pointMtx
	};

	"\n******-------------".postln;
	// calculate and set aeds var for gui update
	aeds = transformedPoints.collect{|ptMtx|
		processPoints.(ptMtx);
	};
	// aeds.collect({|aed|aed[1]}).postln.do(_.postln); // debug: post gains
	// aeds = aeds.collect({|aed|[aed[0], aed[3].abs.ampdb]});
	aeds = aeds.collect({|aedg| aedg[3] = aedg[3].abs.ampdb }); // debug: convert to db in processPoints
	win !? {win.refresh};
};

updateMatrix = {
	var xFormDuples = [];

	"in updateMatrix".postln;

	gControls.do{|ctlDict|
		ctlDict.xform.notNil.if{
			var ctlStates;
			ctlStates = ctlDict.controlStates;
			(ctlStates.size == 0).if{"trasform control states aren't assigned!".error};
			// xFormDuples = xFormDuples.add( [ctlDict.xform, ctlDict.controls] );
			xFormDuples = xFormDuples.add([
				ctlDict.xform,
				ctlStates.indices.collect{|dex| ctlStates[dex] }
			]);
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

	"\n******-------------".postln;
	// calculate and set aeds var for gui update
	aeds = transformedPoints.collect{|ptMtx|
		processPoints.(ptMtx);
	};

	aeds = aeds.collect({|aedg| aedg[3] = aedg[3].abs.ampdb }); // debug: convert to db in processPoints

	win !? {win.refresh};
};

updateMatrix.(); // initialize

/***************************************
Drawing
***************************************/
scrnB = Window.screenBounds;
winW= 600;
winH= 600;
win = Window("soundfield transform",
	Rect(scrnB.center.x - (winW/2), scrnB.height-winH-45, winW, winH),
	resizable: true)
.onClose_({codeWin !? {codeWin.close}}).front;

win.view.palette_(
	QPalette.dark
	.setColor(Color.gray, 'base')
	.setColor(Color.gray, 'alternateBase')
);

uv = UserView( win, Rect(0,0, win.view.bounds.width, win.view.bounds.height/2) ).resize_(5); //.background_(Color.fromHexString("#003366"));

ctlv = View( win, Rect(0,win.view.bounds.height/2, win.view.bounds.width, win.view.bounds.height/2) ).resize_(5);
x = ctlv;

// controls
// xForms = ['push', 'press', 'focus', 'zoom'];
xForms = IdentityDictionary(know: true).putPairs([

	'push', IdentityDictionary(know: true)
	.put( 'controls', [
		// this order defines the order of arguments passed to the transform matrix
		'Degree',	ControlSpec(-pi/2, pi/2, default: 0, units: "pi"),
		'Azimuth',	ControlSpec(pi, -pi, default: 0, units: "pi")
		]
	).put( 'getMatrix',
		{|deg=0, az=0, el=0| FoaXformerMatrix.newPush(deg, az, el).matrix}
	),

	'press', IdentityDictionary(know: true)
	.put( 'controls', [
		// this order defines the order of arguments passed to the transform matrix
		'Degree',	ControlSpec(-pi/2, pi/2, default: 0, units: "pi"),
		'Azimuth',	ControlSpec(pi, -pi, default: 0, units: "pi")
		]
	).put( 'getMatrix',
		{|deg=0, az=0, el=0| FoaXformerMatrix.newPress(deg, az, el).matrix}
	),

	'focus', IdentityDictionary(know: true)
	.put( 'controls', [
		// this order defines the order of arguments passed to the transform matrix
		'Degree',	ControlSpec(-pi/2, pi/2, default: 0, units: "pi"),
		'Azimuth',	ControlSpec(pi, -pi, default: 0, units: "pi")
		]
	).put( 'getMatrix',
		{|deg=0, az=0, el=0| FoaXformerMatrix.newFocus(deg, az, el).matrix}
	),

	'zoom',	IdentityDictionary(know: true)
	.put( 'controls', [
		// this order defines the order of arguments passed to the transform matrix
		'Degree',	ControlSpec(-pi/2, pi/2, default: 0, units: "pi"),
		'Azimuth',	ControlSpec(pi, -pi, default: 0, units: "pi")
		]
	).put( 'getMatrix',
		{|deg=0, az=0, el=0| FoaXformerMatrix.newZoom(deg, az, el).matrix}
	),

	'direct',		IdentityDictionary(know: true)
	.put( 'controls', [
		// this order defines the order of arguments passed to the transform matrix
		'Degree',	ControlSpec(-pi/2, pi/2, default: 0, units: "pi"),
		'Azimuth',	ControlSpec(pi, -pi, default: 0, units: "pi")
		]
	).put( 'getMatrix',
		{|deg=0, az=0, el=0| FoaXformerMatrix.newDirect(deg, az, el).matrix}
	),

	'rotate',	IdentityDictionary(know: true)
	.put( 'controls', [
		// this order defines the order of arguments passed to the transform matrix
		'Azimuth',	ControlSpec(2pi, -2pi, default: 0, units: "pi")
		]
	).put( 'getMatrix',
		{|az=0| FoaXformerMatrix.newRotate(az).matrix}
	),

	'asymmetry',	IdentityDictionary(know: true)
	.put( 'controls', [
		// this order defines the order of arguments passed to the transform matrix
		'Degree',	ControlSpec(-pi/2, pi/2, default: 0, units: "pi")
		]
	).put( 'getMatrix',
		{|deg=0| FoaXformerMatrix.newAsymmetry(deg).matrix}
	),

	'balance',	IdentityDictionary(know: true)
	.put( 'controls', [
		// this order defines the order of arguments passed to the transform matrix
		'Degree',	ControlSpec(-pi/2, pi/2, default: 0, units: "pi")
		]
	).put( 'getMatrix',
		{|deg=0| FoaXformerMatrix.newBalance(deg).matrix}
	),

	'dominate',		IdentityDictionary(know: true)
	.put( 'controls', [
		// this order defines the order of arguments passed to the transform matrix
		'Gain',		ControlSpec(-24, 24, default: 0, units: "dB"),
		'Azimuth',	ControlSpec(pi, -pi, default: 0, units: "pi")
		]
	).put( 'getMatrix',
		{|gain=0, az=0, el=0| FoaXformerMatrix.newDominate(gain, az, el).matrix}
	),

	'gain',	IdentityDictionary(know: true).put( 'controls', [
		// this order defines the order of arguments passed to the transform matrix
			'Gain',		ControlSpec(-24, 24, default: 0, units: "dB")
		]
	).put( 'getMatrix',
		{|gain=0| gain.dbamp }
	),

	'subtract', IdentityDictionary(know: true)
	.put( 'controls', [
		// this order defines the order of arguments passed to the transform matrix
			'This index',	nil,
			'from this index',	nil
		]
	).put( 'getMatrix',
		{|subtractIndex, fromIndex| }
	),
]);

// order determines sf processing chain
// list can grow with inserted items and shrink when removed
gControls = List();

createXForm = {
	var updateControls, addMenu, addSliderCtl, addAddRmvButs;
	var xFormMenu, degSpec, degSl, degNB, degPiNB, rmvBut, addBut;
	var xfStuff, view, layout, deg;
	var curXform; // the currently selected xform

	"in createXForm".postln;
	xfStuff = IdentityDictionary(know: true); // store xform, degree, view, layout here
	xfStuff.put('controlStates', Order()); // to store this xform's controls' states
	degSpec = ControlSpec(); // initialize var

	view = View();
	layout = HLayout().margins_([10,0,10,3]);
	view.layout_(layout);
	xfStuff.view = view;
	xfStuff.layout = layout;

	updateControls = { |xfName| // xf is a dict from xForms
		var xfAttributes, min, max, controls;
		"in updateControls".postln;
		if(xfName != curXform, {

			xfAttributes = xForms[xfName];
			controls = xfAttributes.controls;
			controls = controls.clump(2);

			// rebuild the view with the new controls
			view.children.do(_.remove); // clear the view's elements
			// add menu back, new xform selected
			addMenu.(xfName);
			controls.do{|pair, i|
				var name, ctl;
				#name, ctl = pair;
				ctl.isKindOf(ControlSpec).if({
					addSliderCtl.(name, ctl, i);
					xfStuff.controlStates.put(i, ctl.default);
				});
			};
			addAddRmvButs.();
		})
	};

	addSliderCtl = { |name, spec, ctlOrder|
		var min, max, nb, txt, sl;
		// .min/maxItem to handle negative max on rotate
		min = [spec.minval, spec.maxval].minItem;
		max = [spec.minval, spec.maxval].maxItem;

		sl = Slider()
		.action_(
			{ |sldr| var val;
				val = spec.map(sldr.value);
				if(spec.units == "pi",
					{ nb.value_(val.round(0.001) / pi) },
					{ nb.value_(val.round(0.001)) }
				);
				xfStuff.controlStates[ctlOrder] = val;
				updateMatrix.();
			}
		).orientation_('horizontal').maxHeight_(20)
		.value_(spec.unmap(spec.default));

		nb = NumberBox()
		.action_(
			{ |nb| var val;
				if(spec.units == "pi",
					{ val = nb.value * pi; },
					{ val = nb.value; }
				);
				sl.value_(spec.unmap(val));
				xfStuff.controlStates[ctlOrder] = val;
				updateMatrix.();
			}
		).clipHi_(max).clipLo_(min)
		.step_(0.01).maxWidth_(60);

		txt = StaticText().string_(spec.units).align_('left');

		[nb, txt, sl].do{ |me| layout.add(me) };
	};

	addMenu = { |selectedName|
		xFormMenu = PopUpMenu().items_(['-']++xForms.asSortedArray.collect{|me|me[0]})
		.action_({ |mn|
			if( mn.item != '-', {
				xfStuff.xform_(mn.item);
				updateControls.(mn.item);
				},
				{ xfStuff.xform_(nil) } // nil mutes the transform
			);
			updateMatrix.();
		}).maxWidth_(95).value_(0);

		layout.add(xFormMenu);
		selectedName !? { xFormMenu.value_(xFormMenu.items.indexOf(selectedName)) };
	};


	addAddRmvButs = {
		rmvBut = Button().states_([["X"]])
		.action_({ |but|
			removeXForm.(xfStuff);
			updateMatrix.();
		}).maxWidth_(20);

		addBut = Button().states_([["+"]])
		.action_({ |but|
			createNewXFormView.(view);
		}).maxWidth_(20);

		layout.add(rmvBut);
		layout.add(addBut);
	};

	addMenu.(); // add the first dropdown

	xfStuff; // return this transform's dictionary
};

createNewXFormView = { |behindThisView|
	var behindDex, dict;
	dict = createXForm.();
	behindThisView.notNil.if({
		gControls.do{|ctldict, i|
			(ctldict.view === behindThisView).if{behindDex = i};
		};
		behindDex.isNil.if{"preceeding view not found!".error};
		"behindDex: ".post; behindDex.postln;
		gControls.insert(behindDex+1, dict);
		ctlv.layout.insert(dict.view, behindDex+2); // +2 to account for label row in ctlv[0]
		},{
			gControls.insert(0, dict);
			ctlv.layout.insert(dict.view, 0);
	});
	// return
	dict.view;
};

removeXForm = {|dict| var rmvDex;
	gControls.do{|ctldict, i|
		(ctldict === dict).if{rmvDex = i};
	};
	rmvDex.isNil.if{"view not found!".error};
	dict.view.remove;
	dict.layout.destroy;
	gControls.removeAt(rmvDex);
};


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

showMatrixWindow = { |postCurrentMatrix=false|
	var mwinH = 200;
	var evalMtxString, postMtxString;

	postMtxString =
	"The current transform matrix:\n\n" ++
	"matrix = " ++ curXformMatrix.asString ++
	";\n\nout = AtkMatrixMix.ar(foaSig, matrix);";

	evalMtxString = "var matrix; \n\n/*set you Matrix here, select all and evaluate, or click Evaluate*/\nmatrix = ;\n\n~updateMatrix.(matrix)";

	codeWin.isNil.if( {
		codeWin = Window("Evaluate Matrix",
			Rect(win.bounds.left,
				scrnB.height - win.bounds.height - mwinH
				- (2*20) - 25, // (title bar: 20, OS menu bar: 25)
				win.bounds.width, mwinH))
		.onClose_({codeWin = nil})
		.front;

		View(codeWin, Rect(0,0, codeWin.bounds.width, codeWin.bounds.height))
		.layout_(
			VLayout(
				evalTxtView = TextView().enterInterpretsSelection_(true)
				.string_(postCurrentMatrix.if({postMtxString},{evalMtxString})),

				HLayout(
					Button().states_([["Evaluate"]])
					.action_({evalTxtView.string.interpret}),
					Button().states_([["Reset"]])
					.action_({evalTxtView.string = evalMtxString}),
				),
			)
		);
		},{
			postCurrentMatrix.if({
				evalTxtView.string = postMtxString;
			},{
					evalTxtView.string = evalMtxString;
			});
		}
	);
};

tv = View(uv,Rect(0, 0, uv.bounds.width,uv.bounds.height)).layout_(
	HLayout(
		VLayout(
			VLayout(
				maxGainTxt = StaticText().string_("one")
				.font_(Font('Helvetica', 16)).maxWidth_(180).minHeight_(21),
				minGainTxt = StaticText().string_("two")
				.font_(Font('Helvetica', 16)).maxWidth_(160).minHeight_(21),
			),
			nil,
			Button().states_(
				[["Enter Custom Matrix"]]).action_({|but| showMatrixWindow.() }).maxWidth_(130),
			Button().states_(
				[["Post Matrix"]]).action_({|but| showMatrixWindow.(true) }).maxWidth_(130),
		),
		VLayout(
			VLayout(
				StaticText().string_("Directivity Display"),
				PopUpMenu().items_(['Size and Radius', 'Radius Only'])
				.action_({|mn|
					dirDisplay = switch(mn.value, 0,{'size'},1,{'radius'});
					updateMatrix.();
				}).maxWidth_(130),
				nil
			)
			.setAlignment(0, \topRight)
			.setAlignment(1, \topRight),
			VLayout(
				Button()
					.states_([
						["Pan Test"],
						["Stop test"],
				]).maxWidth_(65),
				Button().states_([
					["Continuous"],["Pulsed"]]).maxWidth_(70),
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
			.setAlignment(0, \bottomRight)
			.setAlignment(1, \bottomRight)
			.setAlignment(2, \bottomRight)
			.setAlignment(3, \bottomRight)
		)
	)
);

win.view.onResize_({ |winview|
	uv.bounds_(Rect(0,0, win.view.bounds.width, win.view.bounds.height/2));
	tv.bounds_(Rect(0, 0, uv.bounds.width,uv.bounds.height));
});

gainThresh = -90;

alphaSpec = ControlSpec(0.1, 1, warp:5);
// gainSpec = ControlSpec(-100, 6, 'db');
gainSpec = ControlSpec(-90, 6, -2);
colorSpec = ControlSpec(768,0);

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


/*******************
Draw the soundfield
*******************/

uv.drawFunc_({|view|
	var cen, circleViewRatio, maxMinStr, testSigPnt, azLineClr;
	var azPnt, drawPnt, omniRad, omniDiam, fullOmni, gainColor, gainPnt;
	var az, el, dir, gain, aeds_sortOrder;

	ctlv.bounds_(Rect(0,win.view.bounds.height/2, win.view.bounds.width, win.view.bounds.height/2));

	r = uv.bounds.height * 0.02;
	d = r*2;
	circleViewRatio = 0.8;
	arcH = uv.bounds.height * circleViewRatio / 2;	// height of the "fan" arc
	fullOmni = 2 * arcH;

	cen = view.bounds.center; // center drawing origin
	Pen.translate(cen.x, cen.y);
	// Pen.strokeColor_(Color.red).fillColor_(Color.blue); // draw background "fan"
	// Pen.fillOval(Rect(r.neg, r.neg, d,d));
	// Pen.strokeOval(Rect(r.neg, r.neg, d,d));
	Pen.addAnnularWedge( 0@0, 5, arcH, 0, 2pi );
	Pen.fillColor_(Color.gray(0.9)).fill;

	// background circles
	Pen.strokeColor_(Color.gray.alpha_(0.2));
	3.do{|i|
		var val;
		val = (i+1 / 3);
		Pen.strokeOval( Rect(
			(arcH * val).neg, (arcH * val).neg,
			arcH*2*val, arcH*2*val
		));
	};

	// draw test signal planewave circle
	testSigPnt = Polar(1.1,rrand(0, 2pi)).asPoint
	.rotate(pi/2)	// convert ambi to screen coords
	* Point(1,-1)	// flip Y for drawing
	* arcH;			// scale normalized points to arcH
	Pen.fillColor_(Color.magenta);
	Pen.fillOval( Rect(testSigPnt.x-r, testSigPnt.y-r, d, d) );

	azLineClr = Color.gray.alpha_(0.2);// line fron center to point

	// get sort order by directivity to draw most transparent first
	aeds_sortOrder = aeds.collect({|arr| arr[2]}).order;
	aeds_sortOrder.do{|sortDex, i|

		#az, el, dir, gain = aeds[sortDex];

		(gain > gainThresh).if{
			omniDiam =	1 - dir * fullOmni;
			omniDiam =	omniDiam.clip(d, fullOmni);
			omniRad	=	omniDiam/2;

			gainColor = getColor.(gain);

			// cartesian point in view coordinates
			azPnt = Point(cos(az), sin(az)) // = Polar(dir, az).asPoint
			.rotate(pi/2)	// convert ambi to screen coords
			* Point(1,-1)	// flip Y for drawing
			* arcH;			// scale normalized points to arcH
			drawPnt = azPnt * dir;

			// line fron center to point
			Pen.strokeColor_(azLineClr);
			Pen.line(drawPnt, 0@0).stroke;

			// directivity circle

			switch(dirDisplay,
				'size',{
					Pen.fillColor_(gainColor.alpha_(alphaSpec.map(dir)));
					Pen.fillOval( Rect(drawPnt.x-omniRad, drawPnt.y-omniRad, omniDiam, omniDiam) );
					gainPnt = azPnt * 1.15; // scale in/out toward/away from origin
				},
				'radius', {
					Pen.fillColor_(gainColor);
					Pen.fillOval( Rect(drawPnt.x-r, drawPnt.y-r, d, d) );
					gainPnt = drawPnt * dir.linlin(0,1,1.75,1.15); // scale in/out toward/away from origin
				}
			);

			// gain labels
			Pen.fillColor_(gainColor.alpha_(1));
			QPen.stringCenteredIn(
				gain.round(0.1).asString,
				Rect(gainPnt.x-(r*10), gainPnt.y-(r*10), d*10, d*10)
			);
		}
	};

	// original 0 deg azimuth point circle ---
	#az, el, dir, gain = aeds[0];
	omniDiam = 1 - dir * fullOmni;
	omniDiam = omniDiam.clip(d, fullOmni);
	omniRad= omniDiam/2;
	azPnt = Point(cos(az), sin(az)) // = Polar(dir, az).asPoint
	.rotate(pi/2)	// convert ambi to screen coords
	* Point(1,-1)	// flip Y for drawing
	* arcH;			// scale normalized points to arcH
	drawPnt = azPnt * dir;

	Pen.strokeColor = Color.fromHexString("#CC0000");
	Pen.width = 3;
	Pen.strokeOval( Rect(drawPnt.x-r, drawPnt.y-r, d, d) );
	// --- original 0 deg azimuth point circle

	maxGainTxt.string_(format("Max Gain:  %", aeds.collect({|me| me[3]}).maxItem.round(0.1)));
	minGainTxt.string_(format("Min Gain:   %", aeds.collect({|me| me[3]}).minItem.round(0.1)));
});

win.refresh;
)


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
