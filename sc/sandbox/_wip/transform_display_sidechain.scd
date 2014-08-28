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
var chainControls, displayControl;
// var xFormDuples = [];
// functions
var processPoints, updateMatrix, getXForm, chainXForms, initNewXForm;
var createNewXFormView, removeXForm, showMatrixWindow;

// drawing vars
var scrnB, sfWin, winW, winH, uv, ctlv, codev, pnts, arcH, codeWin;
var degreeSpec, degreeSl, degreeNB, degreePiNB, xFormMenu, xForms;
var rotSpec, rotSl, rotNB, rotPiNB;
var selDex, adjustRhos, getColCount, layoutFunc, colors;

var tv, maxGainTxt, minGainTxt;
var getColor, alphaSpec, gainSpec, colorSpec, gainThresh;
var dirDisplay = 'size';
var amp, curXformMatrix, evalTxtView;

var xfWin, chainv, abcs;


initPointsMatrices = numPoints.collect{|i|
	var ang;
	ang = (2pi/numPoints) * i;
	FoaEncoderMatrix.newDirection(ang).matrix.addRow([0]); // need to add z manually
};

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

	// ignore elev for now, implement with Spherical to include elev
	[a, e, d_norm, amp].postln; // debug
	[a, e, d_norm, amp]; // return
};

chainXForms = {|...xFormCtlPrs|
	var mtx;
	//debug
	"in chainXForms".postln;
	xFormCtlPrs.postln;

	mtx = Matrix.newIdentity(4); // init original input soundfield

	xFormCtlPrs.do{ |xNameDeg, i|
		var name, ctlVals;
		#name, ctlVals = xNameDeg;

		// (name != 'subtract').if({
			// tranform via matrix multiplication
			//mtx = xForms[name]['getMatrix'].(*ctlVals) * mtx;

		// pass in preceding soundfield followed by the control values for the transform operation, in order
		mtx = xForms[name]['getMatrix'].(mtx, *ctlVals);

		// store resulting matrix by index
		chainControls[chainDex][i].put('mtx', mtx);

	// },{ // subtract
	// 	mtx = Matrix.newIdentity(4) + (mtx * -1);
	// 	// mtx = mtx/2;
	// });
	};

	// mtx * inputMtx;
	// mtx * (Matrix.newIdentity(4) * -1) * inputMtx
	curXformMatrix = mtx;
	mtx
};

// for use in custom matrix window
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

	aeds = aeds.collect({|aedg| aedg[3] = aedg[3].abs.ampdb }); // debug: convert to db in processPoints
	sfWin !? {sfWin.refresh};
};

updateMatrix = { |whichMatrix = 'display'|
	var xFormDuples = [], xformMtx;;

	"in updateMatrix".postln;

	switch( whichMatrix,
		// update the matrix chain
		'chain', {
			chainControls.do{ |chain, i|
				chain.do{ |ctlDict|
					ctlDict.xform.notNil.if{
						var ctlStates;
						ctlStates = ctlDict.controlStates;
						(ctlStates.size == 0).if{"trasform control states aren't assigned!".error};
						xFormDuples = xFormDuples.add([
							ctlDict.xform,
							ctlStates.indices.collect{|dex| ctlStates[dex] }
						]);
					};
				}
			};
		},
		// update the single matrix in "Display Transform"
		'display', {
			displayControl.xform.notNil.if{
				var ctlStates;
				ctlStates = displayControl.controlStates;
				(ctlStates.size == 0).if{"trasform control states aren't assigned!".error};
				xFormDuples = xFormDuples.add([
					displayControl.xform,
					ctlStates.indices.collect{|dex| ctlStates[dex] }
				]);
			};
		}
	);

	"found these duples: ".post; xFormDuples.postln;
	(xFormDuples.size > 0).if({
		// calculate the transform matrix from the chain
		xformMtx = chainXForms.(*xFormDuples);
		transformedPoints = initPointsMatrices.collect{ |pointMtx| xformMtx * pointMtx };
		},{
			transformedPoints = initPointsMatrices } // no transform
	);

	"\n******-------------".postln;
	// calculate and set aeds var for gui update
	aeds = transformedPoints.collect{|ptMtx|
		processPoints.(ptMtx);
	};

	aeds = aeds.collect({|aedg| aedg[3] = aedg[3].abs.ampdb }); // debug: convert to db in processPoints

	sfWin !? {sfWin.refresh};
};

getXMtxByKey = { |key|
	var str, chainDex, xfDex;
	str = key.asString;
	chainDex = str.keep(1).asSymbol;
	xfDex = str.keep(1).asInt;

	chainControls[chainDex][xfDex].mtx;
};

/***************************************
-----  Drawing  -----------------------
***************************************/

scrnB = Window.screenBounds;
winW= 600;
winH= 450;
sfWin = Window( "Soundfield Transform",
	Rect(scrnB.center.x - (winW/2), scrnB.height-winH-45, winW, winH),
	resizable: true
).onClose_(
	{codeWin !? {codeWin.close}}
).front;

sfWin.view.palette_(
	QPalette.dark
	.setColor(Color.gray, 'base')
	.setColor(Color.gray, 'alternateBase')
);

// view containing the soundfield representation
uv = UserView( sfWin, Rect(0,0, sfWin.view.bounds.width, sfWin.view.bounds.height*2/3) ).resize_(5);
// view containing the soundfield display controls and stats
ctlv = View( sfWin, Rect(0,sfWin.view.bounds.height/3, sfWin.view.bounds.width, sfWin.view.bounds.height/3) ).resize_(5);
ctlv.layout_( VLayout().margins_(0) );

// window for xform chain
xfWin = Window( "Transform Chain",
	Rect(scrnB.center.x - (winW/2), scrnB.height-winH-45, 300, 200),
	resizable: true
)
.onClose_().front;

// window pallette
xfWin.view.palette_(
	QPalette.dark
	.setColor(Color.gray, 'base')
	.setColor(Color.gray, 'alternateBase')
);

// view containing the xform chain
chainv = View( xfWin, xfWin.view.bounds ).resize_(5);
chainv.layout_( VLayout().margins_(0) );


// controls
// 'controls' order defines the order of arguments passed to the transform matrix
xForms = IdentityDictionary(know: true).putPairs([

	'push', IdentityDictionary(know: true)
	.put( 'controls', [
		// this order defines the order of arguments passed to the transform matrix
		'Degree',	ControlSpec(-pi/2, pi/2, default: 0, units: "pi"),
		'Azimuth',	ControlSpec(pi, -pi, default: 0, units: "pi")
		]
	).put( 'getMatrix',
		{|mtx, deg=0, az=0, el=0| FoaXformerMatrix.newPush(deg, az, el).matrix * mtx}
	),

	'press', IdentityDictionary(know: true)
	.put( 'controls', [
		'Degree',	ControlSpec(-pi/2, pi/2, default: 0, units: "pi"),
		'Azimuth',	ControlSpec(pi, -pi, default: 0, units: "pi")
		]
	).put( 'getMatrix',
		{|mtx, deg=0, az=0, el=0| FoaXformerMatrix.newPress(deg, az, el).matrix * mtx}
	),

	'focus', IdentityDictionary(know: true)
	.put( 'controls', [
		'Degree',	ControlSpec(-pi/2, pi/2, default: 0, units: "pi"),
		'Azimuth',	ControlSpec(pi, -pi, default: 0, units: "pi")
		]
	).put( 'getMatrix',
		{|mtx, deg=0, az=0, el=0| FoaXformerMatrix.newFocus(deg, az, el).matrix * mtx}
	),

	'zoom',	IdentityDictionary(know: true)
	.put( 'controls', [
		'Degree',	ControlSpec(-pi/2, pi/2, default: 0, units: "pi"),
		'Azimuth',	ControlSpec(pi, -pi, default: 0, units: "pi")
		]
	).put( 'getMatrix',
		{|mtx, deg=0, az=0, el=0| FoaXformerMatrix.newZoom(deg, az, el).matrix * mtx}
	),

	'direct',		IdentityDictionary(know: true)
	.put( 'controls', [
		'Degree',	ControlSpec(-pi/2, pi/2, default: 0, units: "pi"),
		'Azimuth',	ControlSpec(pi, -pi, default: 0, units: "pi")
		]
	).put( 'getMatrix',
		{|mtx, deg=0, az=0, el=0| FoaXformerMatrix.newDirect(deg, az, el).matrix * mtx}
	),

	'rotate',	IdentityDictionary(know: true)
	.put( 'controls', [
		'Azimuth',	ControlSpec(2pi, -2pi, default: 0, units: "pi")
		]
	).put( 'getMatrix',
		{|mtx, az=0| FoaXformerMatrix.newRotate(az).matrix * mtx}
	),

	'asymmetry',	IdentityDictionary(know: true)
	.put( 'controls', [
		'Degree',	ControlSpec(-pi/2, pi/2, default: 0, units: "pi")
		]
	).put( 'getMatrix',
		{|mtx, deg=0| FoaXformerMatrix.newAsymmetry(deg).matrix * mtx}
	),

	'balance',	IdentityDictionary(know: true)
	.put( 'controls', [
		'Degree',	ControlSpec(-pi/2, pi/2, default: 0, units: "pi")
		]
	).put( 'getMatrix',
		{|mtx, deg=0| FoaXformerMatrix.newBalance(deg).matrix * mtx}
	),

	'dominate',		IdentityDictionary(know: true)
	.put( 'controls', [
		'Gain',		ControlSpec(-24, 24, warp: \db, default: 0, units: "dB"),
		'Azimuth',	ControlSpec(pi, -pi, default: 0, units: "pi")
		]
	).put( 'getMatrix',
		{|mtx, gain=0, az=0, el=0| FoaXformerMatrix.newDominate(gain, az, el).matrix * mtx}
	),

	'gain',	IdentityDictionary(know: true)
	.put( 'controls', [
		'Gain',		ControlSpec(-24, 24, warp: \db, default: 0, units: "dB")
		]
	).put( 'getMatrix',
		{|mtx, gainDB=0| gainDB.dbamp * mtx}
	),

	//  ADDITION AND SUBTRACTION TRANSFORMS  //

	// subtract another soundfield from me
	'subtract', IdentityDictionary(know: true)
	.put( 'controls', [
		'this index',	nil,
		'amount',		ControlSpec(-inf, 0, warp: \db, default: 0, units: "dB")
		]
	).put( 'getMatrix',{ |receiverMtx, subtrIndex, amountDB|
		var subMtx;
		subMtx = chainControls[subtrIndex].mtx;
		subMtx ?? {"matrix to 'subtract' not found".error}
		receiverMtx + (subMtx * -1 * amountDB.dbamp);
		}
	),

	// subtract me from another soundfield
	'subtract from', IdentityDictionary(know: true)
	.put( 'controls', [
		'from this index',	nil,
		'amount',		ControlSpec(-inf, 0, warp: \db, default: 0, units: "dB")
		]
	).put( 'getMatrix', { | subtrMtx, receiverDex, amountDB|
		var receiver;
		receiver = chainControls[receiverDex].mtx;
		receiver ?? {"receiver of 'subtract from' not found".error};
		receiver + (subtrMtx * -1 * amountDB.dbamp);
		}
	),

	// add another soundfield to me
	'add', IdentityDictionary(know: true)
	.put( 'controls', [
		'this index',	nil,
		'amount',		ControlSpec(-inf, 0, warp: \db, default: 0, units: "dB")
		]
	).put( 'getMatrix', { |receiverMtx, addIndex, amountDB|
		var addMtx;
		addMtx = chainControls[addIndex].mtx;
		addMtx ?? {"matrix to 'add' not found".error}
		receiverMtx + (addMtx * amountDB.dbamp);
		}
	),

	// add me to another soundfield
	'add to', IdentityDictionary(know: true)
	.put( 'controls', [
		'this index',	nil,
		'amount',		ControlSpec(-inf, 0, warp: \db, default: 0, units: "dB")
		]
	).put( 'getMatrix', { | addMtx, receiverDex, amountDB|
		var receiver;
		receiver = chainControls[receiverDex].mtx;
		receiver ?? {"receiver of 'add to' not found".error};
		receiver + (addMtx * amountDB.dbamp);
		}
	),

]);

initNewXForm = { |target = 'chain'|
	var rebuildControls, addMenu, addSliderCtl, addAddRmvButs;
	var xFormMenu, degSpec, degSl, degNB, degPiNB, rmvBut, addBut;
	var xfStuff, view, layout, deg;

	xfStuff = IdentityDictionary(know: true); // store xform, degree, view, layout here
	xfStuff.put('controlStates', Order()); // to store this xform's controls' states
	degSpec = ControlSpec(); // initialize var

	view = View();
	layout = HLayout().margins_([10,0,10,3]);
	view.layout_(layout);
	xfStuff.view = view;
	xfStuff.layout = layout;

	// when selecting a new control from the dropdown menu
	rebuildControls = { |xfName|
		var xfAttributes, controls;

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
		if(target == 'chain', {addAddRmvButs.()});
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
				updateMatrix.(target);
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
				updateMatrix.(target);
			}
		).clipHi_(max).clipLo_(min)
		.step_(0.01).maxWidth_(60);

		txt = StaticText().string_(spec.units).align_('left');

		[nb, txt, sl].do{ |me| layout.add(me) };
	};

	addMenu = { |selectedName|
		xFormMenu = PopUpMenu().items_(['-']++xForms.asSortedArray.collect{|me|me[0]})
		.action_({ |mn|
			if( mn.item != '-',
				{
					xfStuff.xform_(mn.item);
					rebuildControls.(mn.item);
				},
				{ xfStuff.xform_(nil) } // nil mutes the transform
			);
			updateMatrix.(target);
		}).maxWidth_(95).value_(0);

		layout.add(xFormMenu);
		selectedName !? { xFormMenu.value_(xFormMenu.items.indexOf(selectedName)) };
	};


	addAddRmvButs = {
		rmvBut = Button().states_([["X"]])
		.action_({ |but|
			removeXForm.(xfStuff);
			updateMatrix.(target);
		}).maxWidth_(20);

		addBut = Button().states_([["+"]])
		.action_({ |but|
			createNewXFormView.(view);
		}).maxWidth_(20);

		layout.add(rmvBut);
		layout.add(addBut);
	};

	// add the first dropdown
	addMenu.();

	xfStuff; // return this transform's dictionary
};

createNewXFormView = { |behindThisView|
	var whichChain, behindDex, dict;
	dict = initNewXForm.('chain');
	behindThisView.notNil.if({
		chainControls.do{|chain, j|
			chain.do{|ctldict, i|
				(ctldict.view === behindThisView).if{behindDex = i; whichChain = j};
			}
		};
		behindDex.isNil.if{"preceeding view not found!".error};
		"behindDex: ".post; behindDex.postln; // debug
		chainControls[whichChain].insert(behindDex+1, dict);
		chainv.layout.insert(dict.view, behindDex+2); // +2 to account for label row in ctlv[0]
		},{
			chainControls[0].insert(0, dict); // TODO, how to detect which chain? to insert to?
			chainv.layout.insert(dict.view, 0);
	});
	// return
	dict.view;
};

removeXForm = {|dict|
	var whichChain, rmvDex;
	chainControls.do{ |chain, j|
		chain.do{|ctldict, i|
			(ctldict === dict).if{ rmvDex = i; whichChain = j};
		}
	};
	rmvDex.isNil.if{"view not found!".error};
	dict.view.remove;
	dict.layout.destroy;
	chainControls[whichChain].removeAt(rmvDex);
};

// order determines sf processing chain
// list can grow with inserted items and shrink when removed
// First list added in for the first chain
chainControls = List().add(List());
displayControl = initNewXForm.('display');
updateMatrix.('display'); // initialize
abcs = ['A','B','C','D','E','F','G','H','I','J'];

// add label row to ctlv layout
ctlv.layout.add(
	StaticText().string_("Display Transform").align_('center')
);
// init the first xform
ctlv.layout.add(displayControl.view);

// add label row to chainv layout
chainv.layout.add(
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

chainv.layout.add(createNewXFormView.());


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
			Rect(sfWin.bounds.left,
				scrnB.height - sfWin.bounds.height - mwinH
				- (2*20) - 25, // (title bar: 20, OS menu bar: 25)
				sfWin.bounds.width, mwinH))
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

// Transform View
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

sfWin.view.onResize_({ |winview|
	uv.bounds_(Rect(0,0, sfWin.view.bounds.width, sfWin.view.bounds.height/2));
	tv.bounds_(Rect(0, 0, uv.bounds.width,uv.bounds.height));
});

gainThresh = -90;  // below this thresh, the point isn't displayed

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

	ctlv.bounds_(Rect(0,sfWin.view.bounds.height/2, sfWin.view.bounds.width, sfWin.view.bounds.height/2));

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

sfWin.refresh;
)


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
/*

w = Window().front
l = HLayout( );
w.layout = l;
l.insert(View().minWidth_(200).maxWidth_(200).background_(Color.rand), 0);
(
{
	var newSize = 0;
	w.view.children.postln;
	3.do{|i| w.view.children.reverse.at(i).remove};
	0.1.wait;
	w.view.children.postln;
	w.view.children.do{|me| newSize = newSize + me.bounds.width};
	newSize.postln;
	w.bounds_(w.bounds.width_(newSize));
}.fork(AppClock)
)
*/