DXARTS Room 207A (Studio 1) "Narrow-Quad" room correction protocol
(http://www.dxarts.washington.edu)

see:

1. Angelo Farina, "Simultaneous Measurement of Impulse Response and Distortion with a Swept-Sine Technique," (presented at the Audio Engineering Society Convention no. 108, Audio Engineering Society, 2000).



Loudspeakers:	4 Genelec 1029As
DA/AD:		RME Fireface 800 (828)
Microphone:	Earthworks QTC50


----------------------------------------------------------------

The Python files in this directory are used to correct loudspeaker responses measured in DXARTS Room 207A on xx February 2013.

1) Generate DUT measurement sweep and deconvolution filter via:

	RM207A_DUT_SS_generate.py

This will generate a measurement sweep to be played and recorded via loop-back (plug output into input of Fireface 800 using a TRS cable) and the required deconvolution filter. These are placed in folder 'SS Filters' in a directory defined by the working_dir argument.



----------------------------------------------------------------
[below this line, tbe]


2) Record measurements in the hall. Place the results in a folder titled 'SS Source'.

Notate skiptimes in files for measurements:

[Filename] 
[mic position_spkr position_speaker direction] [skiptime]

Z006004
1_C_F		505002
1_C_L		8177179
1_C_B		15348195
1_C_R		29784837

Z006005
1_L_F		1430801
1_L_L		12552683
1_L_B		28123371
1_L_R		37065097

Z006006
1_R_F		9832901
1_R_L		17412607
1_R_B		32857528
1_R_R		42542653

Z006007
2_C_F		2460564
2_C_L		10212023

Z007001
2_C_B		1932204
2_C_R		11515703

Z007002
2_L_F		4271878
2_L_L		11524418
2_L_B		18571014
2_L_R		26509028

Z007003
2_R_F		1114324
2_R_L		8757690
2_R_B		17022161
2_R_R		24644640


3) Trim measured sweeps via:

	meany_RIR_SS_trim.py

This will trim measured sweeps recorded in the hall. These are placed in folder 'SS Trim' in a directory defined by the working_dir argument.

4) Generate DUT filter via:

	meany_RIR_SS_DUT.py


5) Deconvolve sweeps via:

	meany_RIR_SS_deconvolve.py


xxxxxxx
This will trim measured sweeps recorded in the hall. These are placed in folder 'SS Trim' in a directory defined by the working_dir argument.
