FoaXformChain {
	classvar <xForms;

	// copyArgs
	var <xfAmtPrs;
	var <transformedMatrix, <testDef, <testSynth, <view;
	var <server;

	*initClass {
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
	}

	*getMatrixFromChain {|...xFormAmtPairs|
		var mtx;
		xFormAmtPairs.do{ |xNameDeg|
			var name, deg;
			#name, deg = xNameDeg;
			(name != 'subtract').if({
				// tranform via matrix multiplication
				mtx.isNil.if(
					{ mtx = FoaXformChain.xForms[name]['getMatrix'].(deg) }, // first mtx
					{ mtx = FoaXformChain.xForms[name]['getMatrix'].(deg) * mtx; }
				);
				},{ // subtract
					mtx = Matrix.newIdentity(4) - mtx;
			});
		};
		^mtx;
	}

	*new { |...xFormAmtPairs|
		^super.newCopyArgs(xFormAmtPairs).init;
	}

	init {
		server = Server.local;
		xfAmtPrs  !? {this.updateMatrix(*xfAmtPrs)};
	}

	updateMatrix { |...newxfAmtPrs|
		newxfAmtPrs !? {
			transformedMatrix = FoaXformChain.getMatrixFromChain(*xfAmtPrs);
		}
	}

	loadSynthWithMatrix { |foaXfMatrix, loadCond|
		fork {
			testDef = SynthDef(\foaPanThruXform, {
				arg outbus=0, panRate=0.2, amp=0.25, whichSrc=0, fadeTime=0.2, gate=1;
				var env, sig, panCtl, bf;
				env =EnvGen.kr(Env([0,1,0],fadeTime, \sin, 1), gate, doneAction: 2);
				sig = SelectX.ar(Lag.kr(whichSrc,1), [
					PinkNoise.ar(amp),
					Decay.ar(Impulse.ar(2.5), 0.3, PinkNoise.ar(amp))
				]);
				panCtl = LFSaw.kr(panRate, 1).range(-pi,pi);
				bf = FoaPanB.ar(sig * env, panCtl);
				bf = AtkMatrixMix.ar(bf, foaXfMatrix);
				Out.ar(outbus, bf)
			}).send(server);

			server.sync;
			loadCond.test_(true).signal;
		}
	}

	playTest { |outbus=0, panRate=0.2, amp=0.25, whichSrc=0, fadeTime=0.2|
		transformedMatrix.notNil.if({
			fork {
				var cond;
				cond = Condition();
				this.loadSynthWithMatrix(transformedMatrix, cond);
				cond.wait;

				testSynth = Synth(\foaPanThruXform,
					[\outbus, outbus, \panRate, panRate,
						\amp, amp, \whichSrc, whichSrc, \fadeTime, fadeTime]
				);
			}
		},{ "No transform matrix defined yet".warn })
	}

	stopTest {
		testSynth.set(\gate, 0); // release
	}

	*view {
		^super.new.display;
	}

	display {
		view !? {view = FoaXformChainView(this)};
	}

}

FoaXformChainView {
	var xFormChain; // copyArgs
	var <numPoints = 24, initPointsMatrices, transformedPoints, <aeds, processPoints;

	var <win, scrnB,  winW, winH, uv, ctlv, codev, arcH;
	var xFormMenu, xForms, controls;
	var alphaSpec, colorSpec, gainSpec;

	// doesn't need to have an xForm chain
	*new { |xFormChain|
		^super.newCopyArgs(xFormChain).init;
	}

	init {
		this.numPoints_(numPoints);
		controls = List();

		// init drawing stuff
		scrnB = Window.screenBounds;
		winW= 600;
		winH= 600;
		win = Window("soundfield transform",
			Rect(scrnB.center.x - (winW/2), scrnB.center.y - (winH/2), winW, winH),
			resizable: true).front;

		uv = UserView( win, Rect(0,0, win.view.bounds.width, win.view.bounds.height/2) )
		.resize_(5).background_(Color.gray.alpha_(0.99));

		ctlv = View( win,
			Rect(0,win.view.bounds.height/2, win.view.bounds.width, win.view.bounds.height/2) )
		.resize_(5);

		codev = View();

		this.defineDrawFunc;
		ctlv.layout_( VLayout() );
		ctlv.layout.add(this.getNewXFormView().view); // creat first xform widget, no xform selected
		ctlv.layout.add( codev.layout_(VLayout(TextView().enterInterpretsSelection_(true))) );
		this.updateMatrix; // initialize, refresh window
	}

	defineDrawFunc {

		alphaSpec = ControlSpec(0.1, 1, warp:3);
		// gainSpec = ControlSpec(-100, 6, 'db');
		gainSpec = ControlSpec(-90, 6, -2);
		colorSpec = ControlSpec(768,0);
		/*******************
		Draw the soundfield
		*******************/

		uv.drawFunc_({ |view|
			var r, d, cen, getColor, circleViewRatio;

			ctlv.bounds_(Rect(0,win.view.bounds.height/2, win.view.bounds.width, win.view.bounds.height/2));
			uv.bounds_(Rect(0,0, win.view.bounds.width, win.view.bounds.height/2));
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



			// postf("min gain: %\n", aeds.collect({|me| me[1]}).minItem);
			// postf("max gain: %\n", aeds.collect({|me| me[1]}).maxItem);
			// aeds.do(_.postln);

			aeds.do{|pntGainArr, i|
				var drawPnt, polar, gain, omniRad, omniDiam, fullOmni, gainColor, gainLabelColor, gainPnt;
				#polar, gain = pntGainArr;
				fullOmni = 2 * arcH;
				omniDiam = 1-polar.rho * fullOmni;
				omniDiam = omniDiam.clip(d, fullOmni);
				omniRad= omniDiam/2;

				gainColor = this.getColor(gain);
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
			// draw min/max gain
			Pen.fillColor_(Color.black);
			QPen.stringCenteredIn(
				"max gain\n" ++ aeds.collect({|me| me[1]}).maxItem.round(0.1).asString,
				Rect(ctlv.bounds.width/3.neg, ctlv.bounds.height/4, d*7.5, d*3) );
			QPen.stringCenteredIn(
				"min gain\n" ++ aeds.collect({|me| me[1]}).minItem.round(0.1).asString,
				Rect(ctlv.bounds.width/3.neg, ctlv.bounds.height/2.9, d*7.5, d*3) );
		});
	}

	getPointAED { |xFormedMatrix|
		var w,x,y,z, az, el, g0square, g1square, term, omni, gain;
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
		(term == -1).if{term = -0.99999999}; // protect from NaN
		gain		= w * sqrt(2 / (1 + term));

		// Normalise omni to [0,1] range
		omni = 2 * omni / pi;
		omni = 1 - omni;
		omni = omni.clip(0.0,1.0);

		// ignore elev for now, implement with Spherical to include elev
		// omni is normalized and scales rho (distance from center)
		^[
			Polar(omni, az), // gain.ampdb
			[gain, 0.00000001].maxItem.ampdb // protect from -inf
		];
	}

	updateMatrix {
		var xFormDuples = [];

		"in updateMatrix".postln;
		controls.do{|ctlDict|
			ctlDict.xform.notNil.if{
				ctlDict.degree.isNil.if{"trasform degree is nil!".error};
				xFormDuples = xFormDuples.add( [ctlDict.xform, ctlDict.degree] );
			};
		};

		transformedPoints = (xFormDuples.size > 0).if({
			initPointsMatrices.collect{ |pointMtx|
				var xformMtx;
				// this.xFormPointThruChain(*xFormDuples);
				xformMtx = FoaXformChain.getMatrixFromChain(*xFormDuples);
				xformMtx * pointMtx;
			};
			},{ initPointsMatrices } // no transform
		);

		// calculate and set aeds (az, el, directivity) var for gui update
		aeds = transformedPoints.collect{|ptMtx|
			this.getPointAED(ptMtx);
		};
		win.refresh;
	}

	getNewXFormView { |behindThisView|
		var xFormView, behindDex, dict;
		xFormView = FoaXformView(this);
		dict = xFormView.dict;
		behindThisView.notNil.if({
			controls.do{|ctldict, i|
				(ctldict.view === behindThisView).if{behindDex = i};
			};
			behindDex.isNil.if{ "preceeding view not found!".error };
			controls.insert(behindDex+1, dict);
			ctlv.layout.insert(dict.view, behindDex+1);
			},{
				controls.insert(0, dict);
				ctlv.layout.insert(dict.view, 0);
		});
		^xFormView
	}

	removeXForm { |dict|
		var rmvDex;
		controls.do{|ctldict, i|
			(ctldict === dict).if{rmvDex = i};
		};
		rmvDex.isNil.if{"view not found!".error};
		dict.view.remove;
		dict.layout.destroy;
		controls.removeAt(rmvDex);
	}

	numPoints_ { |numPts|
		numPoints = numPts;
		initPointsMatrices = numPoints.collect{|i|
			var ang;
			ang = (2pi/numPoints) * i;
			// Matrix.with([[kInvSqrt2, cos(ang), sin(ang), 0]]);
			FoaEncoderMatrix.newDirection(ang).matrix.addRow([0]); // need to add z manually
		};
	}

	getColor { |gain|
		var i;
		i = colorSpec.map(gainSpec.unmap(gain));
		^case
		{i < 256} {Color.new255(255, i, 0)}
		{i < 384} {Color.new255(255 - (i-256), 255, 0)}
		{i < 512} {Color.new255(0, 255, (i-384)*2)}
		{i < 768} {Color.new255(0, 255 - (i-512), 255)}
		{i >= 768} {Color.new255(0, 0, 255)}; // catch all
	}
}


FoaXformView {
	var chainView; // copyArgs
	var updateSpecs, xFormMenu, ctlSpec, ctlSl, ctlNB, ctlPiNB, removeBut, addBut;
	var <dict, <view, <layout, <deg;

	*new { |chainView| ^super.newCopyArgs(chainView).init; }

	init {
		view = View();
		// this xform's dictionary: stores degree, view, layout here
		dict = IdentityDictionary(know: true);
		ctlSpec = ControlSpec();
		this.createWidgets;
	}

	updateSpecs { |xf| // xf is a dict from FoaXformChain.xForms
		var min, max;
		"in updatespecs".postln;
		ctlSpec = ControlSpec(xf.min, xf.max, default: xf.default);
		// .min/maxItem to handle negative max on rotate
		min = [xf.min, xf.max].minItem;
		max = [xf.min, xf.max].maxItem;
		ctlNB.clipHi_(max).clipLo_(min);
		ctlPiNB.clipHi_(max).clipLo_(min);
		// initialize widgets
		deg.notNil.if({
			if( (deg<min) or: (deg>max), { deg = ctlSpec.default });
			},{ deg = ctlSpec.default }
		);
		ctlSl.value   = ctlSpec.unmap(deg);
		ctlNB.value   = deg;
		ctlPiNB.value = deg / pi;
	}

	createWidgets {

		xFormMenu = PopUpMenu().items_(['-']++FoaXformChain.xForms.keys)
		.action_({ |mn|
			if( mn.item != '-', {
				var xfName, xfAttributes;
				xfName = mn.item;
				xfAttributes = FoaXformChain.xForms[xfName];
				this.updateSpecs(xfAttributes);
				// update control dict
				dict.xform_(xfName).degree_(deg);
				chainView.updateMatrix;
				},{ // nil mutes the transform
					dict.xform_(nil).degree_(nil);
					chainView.updateMatrix;
			});
		}).maxWidth_(75).value_(0);

		ctlSl = Slider()
		.action_(
			{ |sl| var val;
				val = ctlSpec.map(sl.value);
				deg = val;
				ctlNB.value_(val.round(0.001));
				ctlPiNB.value_(val.round(0.001)/pi);
				dict.degree_(deg);
				chainView.postln;
				chainView.updateMatrix;
			}
		)
		.step_(0.001).orientation_('horizontal');

		ctlNB = NumberBox()
		.action_(
			{ |nb| var val;
				val = nb.value;
				deg = val;
				ctlSl.value_(ctlSpec.unmap(val));
				ctlPiNB.value_(val / pi);
				dict.degree_(deg);
				chainView.updateMatrix;
			}
		)
		.step_(0.01).maxWidth_(60);

		ctlPiNB = NumberBox()
		.action_(
			{ |nb| var val;
				val = nb.value * pi;
				deg = val;
				ctlSl.value_(ctlSpec.unmap(val));
				ctlNB.value_(val);
				dict.degree_(deg);
				chainView.updateMatrix;
			}
		)
		.step_(0.01/pi).maxWidth_(60);

		addBut = Button().states_([["+"]])
		.action_({ |but|
			chainView.getNewXFormView(view);
		})
		.maxWidth_(20);

		removeBut = Button().states_([["X"]])
		.action_({ |but|
			chainView.removeXForm(dict);
			chainView.updateMatrix;
		})
		.maxWidth_(20);

		layout = HLayout(
			xFormMenu, ctlNB,
			VLayout( ctlSl, StaticText().string_("Degree of Tranform") ),
			ctlPiNB, StaticText().string_("pi").align_('left'),
			removeBut, addBut
		).margins_(0);

		view.layout_(layout);

		dict.view = view;
		dict.layout = layout;
	}
}




/*(
// Todo:
// --add gain labels on every other point
// add makeup gain slider
// add matrix subtraction
// consider annular wedge representations
// --add circle on original 0,0 point
// --post min/max gain
// implement units for tranform (pi, db, receiver index)
// add azimuth controls to the transforms


// constants
var kInvSqrt2 = 1/(2.sqrt);
var numPoints = 24, initPointsMatrices, transformedPoints, aeds;
var controls;
// var xFormDuples = [];
// functions
var processPoints, updateMatrix, getXForm, chainXForms, createXForm;
var createNewXFormView, removeXForm;

// drawing vars
var scrnB, winW, winH, uv, ctlv, codev, arcH;
var degreePiNB, xFormMenu, xForms;

processPoints = {|xFormedMatrix|
	var w,x,y,z, az, el, g0square, g1square, term, omni, gain;
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
	(term == -1).if{term = -0.99999999}; // protect from NaN
    gain		= w * sqrt(2 / (1 + term));

	// Normalise omni to [0,1] range
	omni = 2 * omni / pi;
	omni = 1 - omni;
	omni = omni.clip(0.0,1.0);

	// debug ---
	[ Polar(omni, az), // gain.ampdb
		[gain, 0.00000001].maxItem.ampdb // protect from -inf
	].postln;
	[w,x,y,z, az, el].postln;
	omni.postln;
	// --- debug

	// ignore elev for now, implement with Spherical to include elev
	// omni is normalized and scales rho (distance from center)
	[
		Polar(omni, az), // gain.ampdb
		[gain, 0.00000001].maxItem.ampdb // protect from -inf
	];
};

initPointsMatrices = numPoints.collect{|i|
	var ang;
	ang = (2pi/numPoints) * i;
	// Matrix.with([[kInvSqrt2, cos(ang), sin(ang), 0]]);
	FoaEncoderMatrix.newDirection(ang).matrix.addRow([0]); // need to add z manually
};

/*chainXForms = {|...xFormDegPrs|
	var mtx;
	"in chainXForms".postln;
	xFormDegPrs.postln;
	xFormDegPrs.do{ |xNameDeg|
		var name, deg, xMtx;
		#name, deg = xNameDeg;
		mtx.isNil.if(
			{ mtx = xForms[name]['getMatrix'].(deg).matrix; }, // make last xform the first mtx
			{
				xMtx = xForms[name]['getMatrix'].(deg).matrix;
				mtx = xMtx * mtx;
			}
		);
	};
	mtx;
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
	// "Duples: ".post; xFormDuples.postln; // debug
	transformedPoints = (xFormDuples.size > 0).if({
		initPointsMatrices.collect{ |pointMtx|
			var xformMtx;
			xformMtx = chainXForms.(*xFormDuples);
			xformMtx * pointMtx;
		};
		},{ initPointsMatrices } // no transform
	);
	// set aeds var for gui update
	aeds = transformedPoints.collect{|ptMtx|
		processPoints.(ptMtx);
	};

};*/
chainXForms = {|inputMtx...xFormDegPrs|
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
				mtx = Matrix.newIdentity(4) - mtx;
		});
	};
	mtx * inputMtx;
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

	transformedPoints = (xFormDuples.size > 0).if({
		initPointsMatrices.collect{ |pointMtx|
			// var xformMtx;
			// xformMtx = chainXForms.(pointMtx, *xFormDuples);
			// xformMtx * pointMtx;
			chainXForms.(pointMtx, *xFormDuples);
		};
		},{ initPointsMatrices } // no transform
	);

	// calculate and set aeds var for gui update
	aeds = transformedPoints.collect{|ptMtx|
		processPoints.(ptMtx);
	};

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
	.action_({ |but| removeXForm.(dict) })
	.maxWidth_(20);

	layout = HLayout(
		xFormMenu, degNB,
		VLayout( degSl, StaticText().string_("Degree of Tranform") ),
		degPiNB, StaticText().string_("pi").align_('left'),
		removeBut, addBut
	);

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
		controls.insert(behindDex+1, dict);
		ctlv.layout.insert(dict.view, behindDex+1);
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

ctlv.layout_( VLayout() );
ctlv.layout.add(createNewXFormView.());
ctlv.layout.add( codev.layout_(VLayout(TextView().enterInterpretsSelection_(true))) );


/*******************
Draw the soundfield
*******************/

uv.drawFunc_({|view|
	var cen, getColor, circleViewRatio;
	var alphaSpec, gainSpec, colorSpec;

	ctlv.bounds_(Rect(0,w.view.bounds.height/2, w.view.bounds.width, w.view.bounds.height/2));
	uv.bounds_(Rect(0,0, w.view.bounds.width, w.view.bounds.height/2));
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

	postf("min gain: %\n", aeds.collect({|me| me[1]}).minItem);
	postf("max gain: %\n", aeds.collect({|me| me[1]}).maxItem);
	aeds.do(_.postln);

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
	// draw min/max gain
	Pen.fillColor_(Color.black);
	QPen.stringCenteredIn(
		"max gain\n" ++ aeds.collect({|me| me[1]}).maxItem.round(0.1).asString,
		Rect(ctlv.bounds.width/3.neg, ctlv.bounds.height/4, d*7.5, d*3) );
	QPen.stringCenteredIn(
		"min gain\n" ++ aeds.collect({|me| me[1]}).minItem.round(0.1).asString,
		Rect(ctlv.bounds.width/3.neg, ctlv.bounds.height/2.9, d*7.5, d*3) );

	// // draw transducer layout
	// Pen.push;
	// Pen.translate((uv.bounds.width/4).neg, (uv.bounds.height/4));
	// tw = uv.bounds.width/2/maxCols-tgap;
	// th = 10;
	// colgroups.do{|colDict, i|
	// 	Pen.fillColor_(colDict.color);
	// 	colDict.cols.do{
	// 		Pen.fillRect(Rect(0,0, tw, th));
	// 		Pen.translate(tw+tgap,0);
	// 	};
	// };
	// Pen.pop;
	//
	// // draw group angles
	// Pen.push;
	// Pen.translate((uv.bounds.width/4).neg, (uv.bounds.height/4 - th));
	// colgroups.do{|colDict, i|
	// 	Pen.translate((tw+tgap)*(colDict.cols/2),0);
	// 	QPen.stringCenteredIn(colDict.theta.raddeg.round.asString, Rect(tw.neg, th.neg, tw*3, th*2), color: colDict.color);
	// 	Pen.translate((tw+tgap)*(colDict.cols/2),0);
	// };
	// Pen.pop;
});

w.refresh;
)*/

/*
v = FoaXformChainView()
*/