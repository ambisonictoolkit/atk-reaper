// simple test to find nearest angular displacement

(
var mag, a, b, modA, modB, modDiff, delta, unwrapA;

// test magnitude
mag = 10.0 * pi;

// assign vals
a = mag.rand2;
b = mag.rand2;

// let's see
"a = ".post; a.postln; "b = ".post; b.postln;

// wrap vals, (0 - 2pi)
modA = a % 2pi;
modB = b % 2pi;

// let's see
"modA = ".post; modA.postln; "modB = ".post; modB.postln;

// difference needed to increment from B to A
modDiff = modA - modB;

// let's see
"modDiff = ".post; modDiff.postln;

// delta
if ( abs(modDiff) <= pi,
	{
		true.postln;
		delta = modDiff;
	},{
		false.postln;
		if ( modDiff < 0,
			{
				delta = modDiff + 2pi;
			},{
				delta = modDiff - 2pi;
			}
		)
	}
);

// let's see
"delta = ".post; delta.postln;

// unwrapped A
unwrapA = modB + delta;

// compare
"unwrapA = ".post; unwrapA.postln;
"mod(unwrapA) = ".post; (unwrapA % 2pi).postln;
)