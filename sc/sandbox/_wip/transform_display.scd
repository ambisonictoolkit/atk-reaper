(
// constants
var kInvSqrt2 = 1/(2.sqrt);
var numPoints = 8, initPointsMatrices, transformedPoints, aeds;
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

// xFormDuples = [['push', pi/4],['focus', pi/4],['rotate', pi/3]];
// xFormDuples = [];

processPoints = {|xFormedMatrix|
	var w,x,y,z, az, el, g0square, g1square, omni, gain;
	#w,x,y,z = xFormedMatrix.getCol(0);
    az  = atan2(y, x);
	el  = atan2(z, sqrt(x.squared + y.squared));

    // Omni transform angle
    g0square	= w.squared;
    g1square	= x.squared + y.squared + z.squared;
	omni		= asin(( (2*g0square) - g1square) / ((2*g0square) + g1square));

    // Gain
    gain = w * sqrt(2 / (1 + sin(omni)));

	// Normalise omni to [0,1] range
	omni = 2 * omni / pi;
	omni = 1 - omni;

	// ignore el for now, implement with Spherical to include el
	// omni is normalized and scales rho (distance from center)
	[Polar(omni, az), gain.ampdb];
};

initPointsMatrices = numPoints.collect{|i|
	var ang;
	ang = (2pi/numPoints) * i;
	// Matrix.with([[kInvSqrt2, cos(ang), sin(ang), 0]]);
	FoaEncoderMatrix.newDirection(ang).matrix.addRow([0]); // need to add z manually
};

getXForm = {|nameDegPr|
	var name, deg;
	#name, deg = nameDegPr;
	switch( name,
		'push',		{FoaXformerMatrix.newPushX(deg)},
		'press',	{FoaXformerMatrix.newPressX(deg)},
		'focus',	{FoaXformerMatrix.newFocusX(deg)},
		'zoom',		{FoaXformerMatrix.newZoomX(deg)},
		'rotate',	{FoaXformerMatrix.newRotate(deg)}
	).matrix
};

chainXForms = {|...xFormDegPrs|
	var mtx;
	xFormDegPrs.do{ |xNameDeg|
		var xMtx;
		mtx.isNil.if(
			{ mtx = getXForm.(xNameDeg) }, // make last xform the first mtx
			{
				xMtx = getXForm.(xNameDeg);
				mtx = xMtx * mtx;
			}
		);
	};
	mtx;
};

updateMatrix = {
	var xFormDuples = [];
	controls.do{|ctlDict|
		ctlDict.xForm.notNil.if{
			ctlDict.degree.isNil.if{"trasform degree is nil!".error};
			xFormDuples = xFormDuples.add( [ctlDict.xform, ctlDict.degree] );
		};
	};

	transformedPoints = (xFormDuples.size > 0).if({
		initPointsMatrices.collect{|pointMtx|
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
uv = UserView( w, Rect(0,0, w.view.bounds.width, w.view.bounds.height/2) ).resize_(5);
ctlv = View( w, Rect(0,w.view.bounds.height/2, w.view.bounds.width, w.view.bounds.height/2) ).resize_(5);
codev = View();
// ));
// uv = UserView();
// ctlv = View();

// controls
// xForms = ['push', 'press', 'focus', 'zoom'];
xForms = IdentityDictionary(know: true).putPairs([
	'push',		IdentityDictionary(know: true).putPairs([ 'min', -pi/2, 'max', pi/2, 'default', 0 ]),
	'press',	IdentityDictionary(know: true).putPairs([ 'min', -pi/2, 'max', pi/2, 'default', 0 ]),
	'focus',	IdentityDictionary(know: true).putPairs([ 'min', -pi/2, 'max', pi/2, 'default', 0 ]),
	'zoom',		IdentityDictionary(know: true).putPairs([ 'min', -pi/2, 'max', pi/2, 'default', 0 ]),
	'rotate',	IdentityDictionary(know: true).putPairs([ 'min', -2pi,  'max', 2pi,	 'default', 0 ])
]);

// order determines sf processing chain
// list can grow with inserted items and shrink when removed
controls = List();

createXForm = {
	var updateSpecs, xFormMenu, degSpec, degSl, degNB, degPiNB, removeBut, addBut;
	var dict, view, layout, deg;

	"createXForm".postln;

	view = View();
	// store xform, degree, view, layout here
	dict = IdentityDictionary(know: true);

	updateSpecs = { |xf| // xf is a dict from xForms
		var min, max;
		min = xf.min;
		max = xf.max;
		degSpec = ControlSpec(min, max, default: xf.default);
		degNB.clipHi_(max).clipLo_(min);
		degPiNB.clipHi_(max).clipLo_(min);
		// initialize widgets
		deg.notNil.if({
			if( (deg<min) or: (deg>max), { deg = degSpec.default.postln; });
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
			var xf, xfd;
			xf = mn.item;
			xfd = xForms[xf];
			updateSpecs.(xfd);
			// update control dict
			dict.xform_(xf).degree_(deg);
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
	behindThisView !? {
		controls.do{|ctldict, i|
			behindDex = (ctldict.view === behindThisView).if{i};
		};
		behindDex.isNil.if{"preceeding view not found!".error};
		controls.insert(behindDex+1, dict);
		ctlv.layout.insert(dict.view, behindDex+1);
	};
	controls.insert(0, dict);
	ctlv.layout.insert(dict.view, 0);
	// return
	dict.view;
};

removeXForm = {|dict| var rmvDex;
	controls.do{|ctldict, i|
		rmvDex = (ctldict === dict).if{i};
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
	var cen, getColor;
	var alphaSpec, gainSpec, colorSpec;

	ctlv.bounds_(Rect(0,w.view.bounds.height/2, w.view.bounds.width, w.view.bounds.height/2));
	uv.bounds_(Rect(0,0, w.view.bounds.width, w.view.bounds.height/2));
	r = uv.bounds.height * 0.02;
	d = r*2;
	arcH = uv.bounds.height/2.1;	// height of the "fan" arc

	// center drawing origin
	cen = view.bounds.center;
	Pen.translate(cen.x, cen.y);
	// draw background "fan"
	Pen.strokeColor_(Color.red).fillColor_(Color.blue);
	Pen.fillOval(Rect(r.neg, r.neg, d,d));
	Pen.strokeOval(Rect(r.neg, r.neg, d,d));
	Pen.addAnnularWedge( 0@0, 15, arcH, 0, 2pi );
	Pen.fillColor_(Color.gray(0.9)).fill;

	alphaSpec = ControlSpec(0.1, 1, warp:3);
	// gainSpec = ControlSpec(-100, 6, 'db');
		gainSpec = ControlSpec(-90, 6, -2);
	colorSpec = ControlSpec(768,0);

	postf("min gain: %\n", aeds.collect({|me| me[1]}).minItem);
	postf("max gain: %\n", aeds.collect({|me| me[1]}).maxItem);

	getColor = { |gain|
		var i;
		i = colorSpec.map(gainSpec.unmap(gain));
		case
		{i < 256} {Color.new255(255, i, 0)}
		{i < 384} {Color.new255(255 - (i-256), 255, 0)}
		{i < 512} {Color.new255(0, 255, (i-384)*2)}
		{i < 768} {Color.new255(0, 255 - (i-512), 255)};
	};
	/*767.do{|i|
	var color;
	color = case
	{i < 256} {Color.new255(255, i, 0)}
	{i < 384} {Color.new255(255 - (i-256), 255, 0)}
	{i < 512} {Color.new255(0, 255, (i-384)*2)}
	{i < 768} {Color.new255(0, 255 - (i-512), 255)}
	;
	Button(w.view, Rect(0,0, 10,10)).states_([["", color, color]])

};*/


	aeds.do{|pntGainArr, i|
		var drawPnt, polar, gain, omniRad, omniDiam, fullOmni;
		#polar, gain = pntGainArr;
		fullOmni = 2 * arcH;
		omniDiam = 1-polar.rho * fullOmni;
		omniDiam = omniDiam.clip(d, fullOmni);
		omniRad= omniDiam/2;

		// cartesian point in view coordinates
		drawPnt = polar.asPoint
		.rotate(pi/2)	// convert ambi to screen coords
		* Point(1,-1)	// flip Y for drawing
		* arcH;			// scale normalized points to arcH

		// line fron center to point
		Pen.strokeColor_(Color.gray.alpha_(0.5));
		Pen.line(drawPnt, 0@0).stroke;
		// directivity circle
		Pen.fillColor_(getColor.(gain).alpha_(alphaSpec.map(polar.rho)));
		Pen.fillOval( Rect(drawPnt.x-omniRad, drawPnt.y-omniRad, omniDiam, omniDiam) );
		// point circle
		// Pen.fillColor_(Color.fromHexString("#CC0000"));
		// Pen.fillOval( Rect(drawPnt.x-r, drawPnt.y-r, d, d) );
	};

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
