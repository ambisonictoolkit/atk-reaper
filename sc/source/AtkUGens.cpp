/*
 *  AmbisonicUGens.cpp
 *  xSC3plugins
 *
 *  Created by Josh Parmenter on 2/4/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */
/*
 SuperCollider real time audio synthesis system
 Copyright (c) 2002 James McCartney. All rights reserved.
 http://www.audiosynth.com
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */


#include "SC_PlugIn.h"

const double sqrt3div6 = sqrt(3.) * 0.1666666667;
const double sqrt3div2 = sqrt(3.) * 0.5;
const double rsqrt6 = 1. / sqrt(6.);
const double sqrt6div3 = sqrt(6.) * 0.3333333333;

static InterfaceTable *ft;

typedef struct  
{
    float coefs[4][4];
} AtkMatrix;
/*
 
 struct AtoB: public Unit
 {
 };
 */

struct FOAPanB : public Unit 
{
    float m_azimuth, m_elevation, m_W_amp, m_X_amp, m_Y_amp, m_Z_amp;
};

/*
struct FOAPanBWFocus : public Unit 
{
    float m_azimuth, m_elevation, m_focus, m_W_amp, m_X_amp, m_Y_amp, m_Z_amp;
};
 
struct AtkSterToB : public Unit 
{
    float m_azimuth, m_azimuthL, m_azimuthR, m_elevation, m_width;
    float m_W_amp, m_X_ampL, m_Y_ampL, m_X_ampR, m_Y_ampR, m_Z_amp;
};

struct AtkPantoF : public Unit
{
    float m_orientation, m_directivity, m_angle, m_offset;
    float *m_Xamps;
    float *m_Yamps;
};
 */

struct AtkDirectO : public Unit
{
    AtkMatrix matrix;
    float m_angle;
};

struct AtkDirectX : public AtkDirectO { };
struct AtkDirectY : public AtkDirectO { };
struct AtkDirectZ : public AtkDirectO { };

struct AtkRotate : public AtkDirectO { };
struct AtkTilt : public AtkDirectO { };
struct AtkTumble : public AtkDirectO { };

struct AtkFocusX : public AtkDirectO { };
struct AtkFocusY : public AtkDirectO { };
struct AtkFocusZ : public AtkDirectO { };

struct AtkPushX : public AtkDirectO { };
struct AtkPushY : public AtkDirectO { };
struct AtkPushZ : public AtkDirectO { };

struct AtkPressX : public AtkDirectO { };
struct AtkPressY : public AtkDirectO { };
struct AtkPressZ : public AtkDirectO { };

struct AtkZoomX : public AtkDirectO { };
struct AtkZoomY : public AtkDirectO { };
struct AtkZoomZ : public AtkDirectO { };

struct AtkDominateX : public Unit
{
    float m_gain;
    AtkMatrix matrix;	
};

struct AtkDominateY : AtkDominateX { };
struct AtkDominateZ : AtkDominateX { };

struct AtkNFC : public Unit
{
    float m_distanceStart, m_y1x, m_y1y, m_y1z;
};

struct AtkProximity : public Unit
{
    float m_distanceStart, m_y1x, m_y1y, m_y1z;
};

struct AtkPsychoShelf : public Unit
{
	float m_freq, m_k0, m_k1;
	float m_y1w, m_y2w, m_y1x, m_y2x, m_y1y, m_y2y, m_y1z, m_y2z;
};

extern "C"
{
    void load(InterfaceTable *inTable);
    
    void FOAPanB_next_aa(FOAPanB *unit, int inNumSamples);
    void FOAPanB_next_kk(FOAPanB *unit, int inNumSamples);
    void FOAPanB_Ctor(FOAPanB *unit);

    /*
    void FOAPanBWFocus_next_aa(FOAPanBWFocus *unit, int inNumSamples);
    void FOAPanBWFocus_next_kk(FOAPanBWFocus *unit, int inNumSamples);
    void FOAPanBWFocus_Ctor(FOAPanBWFocus *unit);
	 
    void AtkSterToB_next_aaa(AtkSterToB *unit, int inNumSamples);
    void AtkSterToB_next_kkk(AtkSterToB *unit, int inNumSamples);
    void AtkSterToB_Ctor(AtkSterToB *unit);
    
    void AtkPantoF_next(AtkPantoF *unit, int inNumSamples);
    void AtkPantoF_Ctor(AtkPantoF *unit);
    void AtkPantoF_Dtor(AtkPantoF *unit);
	 */
    
    void AtkDirectO_next_a(AtkDirectO *unit, int inNumSamples);
    void AtkDirectO_next_k(AtkDirectO *unit, int inNumSamples);
    void AtkDirectO_Ctor(AtkDirectO* unit);
    
    void AtkDirectX_next_a(AtkDirectX *unit, int inNumSamples);
    void AtkDirectX_next_k(AtkDirectX *unit, int inNumSamples);
    void AtkDirectX_Ctor(AtkDirectX* unit);
    
    void AtkDirectY_next_a(AtkDirectY *unit, int inNumSamples);
    void AtkDirectY_next_k(AtkDirectY *unit, int inNumSamples);
    void AtkDirectY_Ctor(AtkDirectY* unit);
    
    void AtkDirectZ_next_a(AtkDirectZ *unit, int inNumSamples);
    void AtkDirectZ_next_k(AtkDirectZ *unit, int inNumSamples);
    void AtkDirectZ_Ctor(AtkDirectZ* unit);
    
    void AtkRotate_next_a(AtkRotate *unit, int inNumSamples);
    void AtkRotate_next_k(AtkRotate *unit, int inNumSamples);
    void AtkRotate_Ctor(AtkRotate* unit);
    
    void AtkTilt_next_a(AtkTilt *unit, int inNumSamples);
    void AtkTilt_next_k(AtkTilt *unit, int inNumSamples);
    void AtkTilt_Ctor(AtkTilt* unit);
    
    void AtkTumble_next_a(AtkTumble *unit, int inNumSamples);
    void AtkTumble_next_k(AtkTumble *unit, int inNumSamples);
    void AtkTumble_Ctor(AtkTumble* unit);
    
    void AtkFocusX_next_a(AtkFocusX *unit, int inNumSamples);
    void AtkFocusX_next_k(AtkFocusX *unit, int inNumSamples);
    void AtkFocusX_Ctor(AtkFocusX* unit);
    
    void AtkFocusY_next_a(AtkFocusY *unit, int inNumSamples);
    void AtkFocusY_next_k(AtkFocusY *unit, int inNumSamples);
    void AtkFocusY_Ctor(AtkFocusY* unit);
    
    void AtkFocusZ_next_a(AtkFocusZ *unit, int inNumSamples);
    void AtkFocusZ_next_k(AtkFocusZ *unit, int inNumSamples);
    void AtkFocusZ_Ctor(AtkFocusZ* unit);	    
    
    void AtkPushX_next_a(AtkPushX *unit, int inNumSamples);
    void AtkPushX_next_k(AtkPushX *unit, int inNumSamples);
    void AtkPushX_Ctor(AtkPushX* unit);	
    
    void AtkPushY_next_a(AtkPushY *unit, int inNumSamples);
    void AtkPushY_next_k(AtkPushY *unit, int inNumSamples);
    void AtkPushY_Ctor(AtkPushY* unit);	
    
    void AtkPushZ_next_a(AtkPushZ *unit, int inNumSamples);
    void AtkPushZ_next_k(AtkPushZ *unit, int inNumSamples);
    void AtkPushZ_Ctor(AtkPushZ* unit);	
    
    void AtkPressX_next_a(AtkPressX *unit, int inNumSamples);
    void AtkPressX_next_k(AtkPressX *unit, int inNumSamples);
    void AtkPressX_Ctor(AtkPressX* unit);
    
    void AtkPressY_next_a(AtkPressY *unit, int inNumSamples);
    void AtkPressY_next_k(AtkPressY *unit, int inNumSamples);
    void AtkPressY_Ctor(AtkPressY* unit);
    
    void AtkPressZ_next_a(AtkPressZ *unit, int inNumSamples);
    void AtkPressZ_next_k(AtkPressZ *unit, int inNumSamples);
    void AtkPressZ_Ctor(AtkPressZ* unit);
    
    void AtkZoomX_next_a(AtkZoomX *unit, int inNumSamples);
    void AtkZoomX_next_k(AtkZoomX *unit, int inNumSamples);
    void AtkZoomX_Ctor(AtkZoomX* unit);
    
    void AtkZoomY_next_a(AtkZoomY *unit, int inNumSamples);
    void AtkZoomY_next_k(AtkZoomY *unit, int inNumSamples);
    void AtkZoomY_Ctor(AtkZoomY* unit);
    
    void AtkZoomZ_next_a(AtkZoomZ *unit, int inNumSamples);
    void AtkZoomZ_next_k(AtkZoomZ *unit, int inNumSamples);
    void AtkZoomZ_Ctor(AtkZoomZ* unit);
    
    void AtkDominateX_next_a(AtkDominateX *unit, int inNumSamples);
    void AtkDominateX_next_k(AtkDominateX *unit, int inNumSamples);
    void AtkDominateX_Ctor(AtkDominateX* unit);						
    
    void AtkDominateY_next_a(AtkDominateY *unit, int inNumSamples);
    void AtkDominateY_next_k(AtkDominateY *unit, int inNumSamples);
    void AtkDominateY_Ctor(AtkDominateY* unit);
    
    void AtkDominateZ_next_a(AtkDominateZ *unit, int inNumSamples);
    void AtkDominateZ_next_k(AtkDominateZ *unit, int inNumSamples);
    void AtkDominateZ_Ctor(AtkDominateZ* unit);
    
    void AtkNFC_next_k(AtkNFC *unit, int inNumSamples);
    void AtkNFC_next_a(AtkNFC *unit, int inNumSamples);
    void AtkNFC_Ctor(AtkNFC* unit);
    
    void AtkProximity_next_k(AtkProximity *unit, int inNumSamples);
    void AtkProximity_next_a(AtkProximity *unit, int inNumSamples);
    void AtkProximity_Ctor(AtkProximity* unit);
    
    void AtkPsychoShelf_next_k(AtkPsychoShelf *unit, int inNumSamples);
    void AtkPsychoShelf_next_a(AtkPsychoShelf *unit, int inNumSamples);
    void AtkPsychoShelf_Ctor(AtkPsychoShelf* unit);
    
    
    //	void AtoB_Ctor(AtoB* unit);
    //	void AtoB_next_00(AtoB *unit, int inNumSamples);
    //	void AtoB_next_01(AtoB *unit, int inNumSamples);	
    //	void AtoB_next_02(AtoB *unit, int inNumSamples);		
    //	void AtoB_next_10(AtoB *unit, int inNumSamples);
    //	void AtoB_next_11(AtoB *unit, int inNumSamples);
    //	void AtoB_next_12(AtoB *unit, int inNumSamples);
    //	void AtoB_next_20(AtoB *unit, int inNumSamples);
    //	void AtoB_next_21(AtoB *unit, int inNumSamples);
    //	void AtoB_next_22(AtoB *unit, int inNumSamples);
    //	void AtoB_next_30(AtoB *unit, int inNumSamples);
    //	void AtoB_next_31(AtoB *unit, int inNumSamples);
    //	void AtoB_next_32(AtoB *unit, int inNumSamples);
    //	void AtoB_next_40(AtoB *unit, int inNumSamples);
    //	void AtoB_next_41(AtoB *unit, int inNumSamples);
    //	void AtoB_next_42(AtoB *unit, int inNumSamples);
    //	void AtoB_next_50(AtoB *unit, int inNumSamples);
    //	void AtoB_next_51(AtoB *unit, int inNumSamples);
    //	void AtoB_next_52(AtoB *unit, int inNumSamples);
    //	void AtoB_next_60(AtoB *unit, int inNumSamples);
    //	void AtoB_next_61(AtoB *unit, int inNumSamples);
    //	void AtoB_next_62(AtoB *unit, int inNumSamples);	
    //	void AtoB_next_70(AtoB *unit, int inNumSamples);
    //	void AtoB_next_71(AtoB *unit, int inNumSamples);
    //	void AtoB_next_72(AtoB *unit, int inNumSamples);	
}

inline float calcmatrixval(float coef, float curval){
    float val;
    if(coef == 0.){
	val = 0.;
    } else {
	if(coef == 1.){
	    val = curval;
	} else {
	    val = coef * curval;
	}
    }
    return val;
}

// can perhaps optimize a bit here ... check for 0s???
#define CALC_MATRIX \
float curvals[4] = {Win[i], Xin[i], Yin[i], Zin[i]}; \
for(int j = 0; j < 4; j++){ \
Wout[i] += calcmatrixval(matrix.coefs[0][j], curvals[j]); \
Xout[i] += calcmatrixval(matrix.coefs[1][j], curvals[j]); \
Yout[i] += calcmatrixval(matrix.coefs[2][j], curvals[j]); \
Zout[i] += calcmatrixval(matrix.coefs[3][j], curvals[j]); \
} 

#define SETUP_TRANSFORMS \
float *Win = IN(0); \
float *Xin = IN(1); \
float *Yin = IN(2); \
float *Zin = IN(3); \
float *Wout = OUT(0); \
float *Xout = OUT(1); \
float *Yout = OUT(2); \
float *Zout = OUT(3); \
ClearUnitOutputs(unit, inNumSamples); \
AtkMatrix matrix = unit->matrix; \

#define ZERO_MATRIX \
for(int i = 0; i < 4; i++){ \
for(int j = 0; j < 4; j++){ \
unit->matrix.coefs[i][j] = 0.f; \
} \
} \

#define SIN_COS \
sina = sin(azimuth); \
sinb = sin(elevation); \
\
cosa = cos(azimuth); \
cosb = cos(elevation); \

/* FOAPanB - basic encoder (places sound on the sphere) */

void FOAPanB_Ctor(FOAPanB *unit)
{
    if((INRATE(1) == calc_FullRate) && (INRATE(2) == calc_FullRate)){
	SETCALC(FOAPanB_next_aa);//aa
    } else {
	SETCALC(FOAPanB_next_kk);//ak
    }
    
    float azimuth = unit->m_azimuth = IN0(1);
    float elevation = unit->m_elevation = IN0(2);
    float sina, sinb, cosa, cosb;
    
    SIN_COS
    
    unit->m_W_amp = rsqrt2;
    unit->m_X_amp = cosa * cosb;
    unit->m_Y_amp = sina * cosb;
    unit->m_Z_amp = sinb;
    
    FOAPanB_next_kk(unit, 1);
}

void FOAPanB_next_kk(FOAPanB *unit, int inNumSamples)
{
    float azimuth = IN0(1);
    float elevation = IN0(2);
    float *in = IN(0);
    float *Wout = OUT(0); 
    float *Xout = OUT(1); 
    float *Yout = OUT(2); 
    float *Zout = OUT(3);
    float Wamp = unit->m_W_amp;
    float Xamp = unit->m_X_amp;
    float Yamp = unit->m_Y_amp;
    float Zamp = unit->m_Z_amp;
    float sina, sinb, cosa, cosb;
    
    if((unit->m_azimuth == azimuth) && (unit->m_elevation == elevation)){
	for(int i = 0; i < inNumSamples; i++){
	    Wout[i] = in[i] * Wamp;
	    Xout[i] = in[i] * Xamp;
	    Yout[i] = in[i] * Yamp;
	    Zout[i] = in[i] * Zamp;
	}
    } else {
	
	SIN_COS
	
	float nextXamp = cosa * cosb;
	float nextYamp = sina * cosb;
	float nextZamp = sinb;
	
	float xSlope = CALCSLOPE(nextXamp, Xamp);
	float ySlope = CALCSLOPE(nextYamp, Yamp);
	float zSlope = CALCSLOPE(nextZamp, Zamp);
	
	for(int i = 0; i < inNumSamples; i++){
	    Wout[i] = in[i] * Wamp;
	    Xout[i] = in[i] * Xamp;
	    Yout[i] = in[i] * Yamp;
	    Zout[i] = in[i] * Zamp;
	    
	    Xamp += xSlope;
	    Yamp += ySlope;  
	    Zamp += zSlope;
	}	
	
	unit->m_X_amp = Xamp;
	unit->m_Y_amp = Yamp;
	unit->m_Z_amp = Zamp;
	unit->m_azimuth = azimuth;
	unit->m_elevation = elevation;
    }
}

void FOAPanB_next_aa(FOAPanB *unit, int inNumSamples)
{
    float *pazimuth = IN(1);
    float *pelevation = IN(2);
    float *in = IN(0);
    float *Wout = OUT(0); 
    float *Xout = OUT(1); 
    float *Yout = OUT(2); 
    float *Zout = OUT(3);
    float Wamp = unit->m_W_amp;
    float Xamp = unit->m_X_amp;
    float Yamp = unit->m_Y_amp;
    float Zamp = unit->m_Z_amp;
    float sina, sinb, cosa, cosb, azimuth, elevation;
    
    for(int i = 0; i < inNumSamples; i++){
	if((unit->m_azimuth == pazimuth[i]) && (unit->m_elevation == pelevation[i])){
	    Wout[i] = in[i] * Wamp;
	    Xout[i] = in[i] * Xamp;
	    Yout[i] = in[i] * Yamp;
	    Zout[i] = in[i] * Zamp;
	} else {
	    
	    azimuth = pazimuth[i];
	    elevation = pelevation[i];
	    
	    SIN_COS
	    
	    Xamp = cosa * cosb;
	    Yamp = sina * cosb;
	    Zamp = sinb;
	    
	    Wout[i] = in[i] * Wamp;
	    Xout[i] = in[i] * Xamp;
	    Yout[i] = in[i] * Yamp;
	    Zout[i] = in[i] * Zamp;
	    
	    unit->m_azimuth = azimuth;
	    unit->m_elevation = elevation;	    
	}
    }
    unit->m_X_amp = Xamp;
    unit->m_Y_amp = Yamp;
    unit->m_Z_amp = Zamp;
}


/* FOAPanB - basic encoder (places sound on the sphere) */

/*
void FOAPanBWFocus_Ctor(FOAPanBWFocus *unit)
{
    if((INRATE(1) == calc_FullRate) && (INRATE(2) == calc_FullRate) && (INRATE(3) == calc_FullRate)){
	SETCALC(FOAPanBWFocus_next_aa);//aa
    } else {
	SETCALC(FOAPanBWFocus_next_kk);//ak
    }
    
    float azimuth = unit->m_azimuth = IN0(1);
    float elevation = unit->m_elevation = IN0(2);
    float focus = unit->m_focus = IN0(3);
    float sina, sinb, cosa, cosb;
    float sinf;
    
    SIN_COS
    
    sinf = sin(focus);
    unit->m_W_amp = rsqrt2;
    unit->m_X_amp = cosa * cosb * sinf;
    unit->m_Y_amp = sina * cosb* sinf;
    unit->m_Z_amp = sinb * sinf;
    
    FOAPanBWFocus_next_kk(unit, 1);
}

void FOAPanBWFocus_next_kk(FOAPanBWFocus *unit, int inNumSamples)
{
    float azimuth = IN0(1);
    float elevation = IN0(2);
    float focus = IN0(3);
    float *in = IN(0);
    float *Wout = OUT(0); 
    float *Xout = OUT(1); 
    float *Yout = OUT(2); 
    float *Zout = OUT(3);
    float Wamp = unit->m_W_amp;
    float Xamp = unit->m_X_amp;
    float Yamp = unit->m_Y_amp;
    float Zamp = unit->m_Z_amp;
    float sina, sinb, cosa, cosb, sinf;
    
    if((unit->m_azimuth == azimuth) && (unit->m_elevation == elevation) && (unit->m_focus == focus)){
	for(int i = 0; i < inNumSamples; i++){
	    Wout[i] = in[i] * Wamp;
	    Xout[i] = in[i] * Xamp;
	    Yout[i] = in[i] * Yamp;
	    Zout[i] = in[i] * Zamp;
	}
    } else {
	
	SIN_COS
	
	sinf = sin(focus);
	float nextXamp = cosa * cosb * sinf;
	float nextYamp = sina * cosb * sinf;
	float nextZamp = sinb * sinf;
	
	float xSlope = CALCSLOPE(nextXamp, Xamp);
	float ySlope = CALCSLOPE(nextYamp, Yamp);
	float zSlope = CALCSLOPE(nextZamp, Zamp);
	
	for(int i = 0; i < inNumSamples; i++){
	    Wout[i] = in[i] * Wamp;
	    Xout[i] = in[i] * Xamp;
	    Yout[i] = in[i] * Yamp;
	    Zout[i] = in[i] * Zamp;
	    
	    Xamp += xSlope;
	    Yamp += ySlope;  
	    Zamp += zSlope;
	}	
	
	unit->m_X_amp = Xamp;
	unit->m_Y_amp = Yamp;
	unit->m_Z_amp = Zamp;
	unit->m_azimuth = azimuth;
	unit->m_elevation = elevation;
	unit->m_focus = focus;
    }
}

void FOAPanBWFocus_next_aa(FOAPanBWFocus *unit, int inNumSamples)
{
    float *pazimuth = IN(1);
    float *pelevation = IN(2);
    float *pfocus = IN(3);
    float *in = IN(0);
    float *Wout = OUT(0); 
    float *Xout = OUT(1); 
    float *Yout = OUT(2); 
    float *Zout = OUT(3);
    float Wamp = unit->m_W_amp;
    float Xamp = unit->m_X_amp;
    float Yamp = unit->m_Y_amp;
    float Zamp = unit->m_Z_amp;
    float sina, sinb, cosa, cosb, sinf, azimuth, elevation, focus;
    
    for(int i = 0; i < inNumSamples; i++){
	if((unit->m_azimuth == pazimuth[i]) && (unit->m_elevation == pelevation[i]) && (unit->m_focus == pfocus[i])){
	    Wout[i] = in[i] * Wamp;
	    Xout[i] = in[i] * Xamp;
	    Yout[i] = in[i] * Yamp;
	    Zout[i] = in[i] * Zamp;
	} else {
	    
	    azimuth = pazimuth[i];
	    elevation = pelevation[i];
	    focus = pfocus[i];
	    
	    SIN_COS
	    
	    sinf = sin(focus);
	    Xamp = cosa * cosb * sinf;
	    Yamp = sina * cosb * sinf;
	    Zamp = sinb * sinf;
	    
	    Wout[i] = in[i] * Wamp;
	    Xout[i] = in[i] * Xamp;
	    Yout[i] = in[i] * Yamp;
	    Zout[i] = in[i] * Zamp;
	    
	    unit->m_azimuth = azimuth;
	    unit->m_elevation = elevation;
	    unit->m_focus = focus;
	}
    }
    unit->m_X_amp = Xamp;
    unit->m_Y_amp = Yamp;
    unit->m_Z_amp = Zamp;
}
*/
/* AtkSterToB - basic encoder (places stereo sound on the sphere) */
/*
void AtkSterToB_Ctor(AtkSterToB *unit)
{
    if((INRATE(2) == calc_FullRate) && (INRATE(3) == calc_FullRate) && (INRATE(4) == calc_FullRate)){
	SETCALC(AtkSterToB_next_aaa);//aaa
    } else {
	SETCALC(AtkSterToB_next_kkk);//kkk
    }
    
    float azimuth = unit->m_azimuth = IN0(2);
    float width = unit->m_width = IN0(3);
    float elevation = unit->m_elevation = IN0(4);
    float width2 = width * 0.5;
    float azimuthL = unit->m_azimuthL = azimuth + width2;
    float azimuthR = unit->m_azimuthR = azimuth - width2;
    
    float sinaL, sinaR, sinb, cosaL, cosaR, cosb;
    
    sinaL = sin(azimuthL); 
    sinaR = sin(azimuthR); 
    sinb = sin(elevation); 
    
    cosaL = cos(azimuthL); 
    cosaR = cos(azimuthR);
    cosb = cos(elevation);     
    
    unit->m_W_amp = rsqrt2;
    unit->m_X_ampL = cosaL * cosb;
    unit->m_X_ampR = cosaR * cosb;
    unit->m_Y_ampL = sinaL * cosb;
    unit->m_Y_ampR = sinaR * cosb;
    unit->m_Z_amp = sinb;
    
    AtkSterToB_next_kkk(unit, 1);
}

void AtkSterToB_next_kkk(AtkSterToB *unit, int inNumSamples)
{
    float azimuth = IN0(2);
    float width = IN0(3);
    float elevation = IN0(4);
    float *inL = IN(0);
    float *inR = IN(1);
    float *Wout = OUT(0); 
    float *Xout = OUT(1); 
    float *Yout = OUT(2); 
    float *Zout = OUT(3);
    float Wamp = unit->m_W_amp;
    float XampL = unit->m_X_ampL;
    float XampR = unit->m_X_ampR;
    float YampL = unit->m_Y_ampL;
    float YampR = unit->m_Y_ampR;
    float Zamp = unit->m_Z_amp;
    float sinaL, sinaR, sinb, cosaL, cosaR, cosb;
    
    if((unit->m_azimuth == azimuth) && (unit->m_elevation == elevation) && (unit->m_width == width)){
	for(int i = 0; i < inNumSamples; i++){
	    Wout[i] = (inL[i] + inR[i]) * Wamp;
	    Xout[i] = (inL[i] * XampL) + (inR[i] * XampR);
	    Yout[i] = (inL[i] * YampL) + (inR[i] * YampR);
	    Zout[i] = (inL[i] + inR[i]) * Zamp;
	}
    } else {
	
	float width2 = width * 0.5;
	float azimuthL = unit->m_azimuthL = azimuth + width2;
	float azimuthR = unit->m_azimuthR = azimuth - width2;
	
	sinaL = sin(azimuthL); 
	sinaR = sin(azimuthR); 
	sinb = sin(elevation); 
	
	cosaL = cos(azimuthL); 
	cosaR = cos(azimuthR); 
	cosb = cos(elevation); 
	
	float nextXampL = cosaL * cosb;
	float nextXampR = cosaR * cosb;
	float nextYampL = sinaL * cosb;
	float nextYampR = sinaR * cosb;
	float nextZamp = sinb;
	
	float xSlopeL = CALCSLOPE(nextXampL, XampL);
	float xSlopeR = CALCSLOPE(nextXampR, XampR);
	float ySlopeL = CALCSLOPE(nextYampL, YampL);
	float ySlopeR = CALCSLOPE(nextYampR, YampR);
	float zSlope = CALCSLOPE(nextZamp, Zamp);
	
	for(int i = 0; i < inNumSamples; i++){
	    Wout[i] = (inL[i] + inR[i]) * Wamp;
	    Xout[i] = (inL[i] * XampL) + (inR[i] * XampR);
	    Yout[i] = (inL[i] * YampL) + (inR[i] * YampR);
	    Zout[i] = (inL[i] + inR[i]) * Zamp;
	    
	    XampL += xSlopeL;
	    XampR += xSlopeR;
	    YampL += ySlopeL;  
	    YampR += ySlopeR; 
	    Zamp += zSlope;
	}	
	
	unit->m_X_ampL = XampL;
	unit->m_X_ampR = XampR;
	unit->m_Y_ampL = YampL;
	unit->m_Y_ampR = YampR;
	unit->m_Z_amp = Zamp;
	unit->m_azimuth = azimuth;
	unit->m_elevation = elevation;
    }
}

void AtkSterToB_next_aaa(AtkSterToB *unit, int inNumSamples)
{
    float *pazimuth = IN(2);
    float *pwidth = IN(3);
    float *pelevation = IN(4);
    float *inL = IN(0);
    float *inR = IN(1);
    float *Wout = OUT(0); 
    float *Xout = OUT(1); 
    float *Yout = OUT(2); 
    float *Zout = OUT(3);
    float Wamp = unit->m_W_amp;
    float XampL = unit->m_X_ampL;
    float XampR = unit->m_X_ampR;
    float YampL = unit->m_Y_ampL;
    float YampR = unit->m_Y_ampR;
    float Zamp = unit->m_Z_amp;
    float sinaL, sinaR, sinb, cosaL, cosaR, cosb;
    
    for(int i = 0; i < inNumSamples; i++){
	if((unit->m_azimuth == pazimuth[i]) && (unit->m_elevation == pelevation[i]) && (unit->m_width == pwidth[i])){
	    Wout[i] = (inL[i] + inR[i]) * Wamp;
	    Xout[i] = (inL[i] * XampL) + (inR[i] * XampR);
	    Yout[i] = (inL[i] * YampL) + (inR[i] * YampR);
	    Zout[i] = (inL[i] + inR[i]) * Zamp;
	} else {
	    
	    float azimuth = pazimuth[i];
	    float width = pwidth[i];
	    float elevation = pelevation[i];
	    float width2 = width * 0.5;
	    float azimuthL = unit->m_azimuthL = azimuth + width2;
	    float azimuthR = unit->m_azimuthR = azimuth - width2;
	    
	    sinaL = sin(azimuthL); 
	    sinaR = sin(azimuthR); 
	    sinb = sin(elevation); 
	    
	    cosaL = cos(azimuthL); 
	    cosaR = cos(azimuthR); 
	    cosb = cos(elevation); 
	    
	    XampL = cosaL * cosb;
	    XampR = cosaR * cosb;
	    YampL = sinaL * cosb;
	    YampR = sinaR * cosb;
	    Zamp = sinb;
	    
	    Wout[i] = (inL[i] + inR[i]) * Wamp;
	    Xout[i] = (inL[i] * XampL) + (inR[i] * XampR);
	    Yout[i] = (inL[i] * YampL) + (inR[i] * YampR);
	    Zout[i] = (inL[i] + inR[i]) * Zamp;
	    
	    unit->m_X_ampL = XampL;
	    unit->m_X_ampR = XampR;
	    unit->m_Y_ampL = YampL;
	    unit->m_Y_ampR = YampR;
	    unit->m_Z_amp = Zamp;
	    unit->m_azimuth = azimuth;
	    unit->m_elevation = elevation;
	    unit->m_width = width;
	}
    }
}
*/

/* AtkPantoF 2D decoder */
/*
void AtkPantoF_Ctor(AtkPantoF* unit)
{
    // should be 0.0 or 1.0
    unit->m_orientation = IN0(3);
    // should be -1.0, 0.0 or 1.0
    unit->m_directivity = IN0(4);
    float g1 = powf(2., unit->m_directivity * -0.5);	
    int numOutputs = unit->mNumOutputs;
    float iNumOutputsPi = pi / (float)numOutputs;
    unit->m_Xamps = (float*)RTAlloc(unit->mWorld, numOutputs * sizeof(float));
    unit->m_Yamps = (float*)RTAlloc(unit->mWorld, numOutputs * sizeof(float));
    float angle;
    if(unit->m_orientation > 0){
	for(int i = 0; i < numOutputs; i++){
	    angle = (float)(1 + (i * 2)) * iNumOutputsPi;
	    unit->m_Xamps[i] = g1 * cos(angle);
	    unit->m_Yamps[i] = g1 * sin(angle);
	}
    } else {
	for(int i = 0; i < numOutputs; i++){
	    angle = (float)i * (iNumOutputsPi * 2);
	    unit->m_Xamps[i] = g1 * cos(angle);
	    unit->m_Yamps[i] = g1 * sin(angle);
	}	    
    }
    AtkPantoF_next(unit, 1);
    SETCALC(AtkPantoF_next);
    
}

void AtkPantoF_Dtor(AtkPantoF *unit){
    RTFree(unit->mWorld, unit->m_Xamps);
    RTFree(unit->mWorld, unit->m_Yamps);
}

void AtkPantoF_next(AtkPantoF *unit, int inNumSamples)
{
    float *Win = IN(0);
    float *Xin = IN(1);
    float *Yin = IN(2);
    int numOutputs = unit->mNumOutputs;
    for(int i = 0; i < numOutputs; i++){
	float *out = OUT(i);
	for(int j = 0; j < inNumSamples; j++){
	    out[j] = Win[j] + (Xin[j] * unit->m_Xamps[i]) + (Yin[j] * unit->m_Yamps[i]);
	}
    }
}
*/

////////////////////// AtkRotate ///////////////////////
// uses 'angle' var. 

#define FILL_ROTATE_MATRIX \
matrix.coefs[0][0] = matrix.coefs[3][3] = 1.; \
double cosa = cos(unit->m_angle); \
double sina = sin(unit->m_angle); \
matrix.coefs[1][1] = matrix.coefs[2][2] = cosa; \
matrix.coefs[1][2] = -sina; \
matrix.coefs[2][1] = sina;

void AtkRotate_Ctor(AtkRotate* unit)
{
    ZERO_MATRIX
    unit->m_angle = IN0(4);
    AtkMatrix matrix = unit->matrix;
    FILL_ROTATE_MATRIX;
    unit->matrix = matrix;
    if(INRATE(4) == calc_FullRate)
	SETCALC(AtkRotate_next_a);
    else
	SETCALC(AtkRotate_next_k);
    AtkRotate_next_k(unit, 1); 
}

void AtkRotate_next_a(AtkRotate *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float *angle = IN(4);
    for(int i = 0; i < inNumSamples; i++){
	if(angle[i] != unit->m_angle){
	    unit->m_angle = angle[i];
	    FILL_ROTATE_MATRIX
	}
	CALC_MATRIX
    }
    unit->matrix = matrix;
    
}

void AtkRotate_next_k(AtkRotate *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float angle = IN0(4);
    float angleslope;
    if(angle != unit->m_angle){
	angleslope = CALCSLOPE(angle, unit->m_angle);
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	    unit->m_angle += angleslope;
	    FILL_ROTATE_MATRIX
	}	
    } else {
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	}
    }
    unit->matrix = matrix;
    unit->m_angle = angle;	
}

////////////////////// AtkTilt ///////////////////////
// uses 'angle' var. 

#define FILL_TILT_MATRIX \
matrix.coefs[0][0] = matrix.coefs[1][1] = 1.; \
double cosa = cos(unit->m_angle); \
double sina = sin(unit->m_angle); \
matrix.coefs[2][2] = matrix.coefs[3][3] = cosa; \
matrix.coefs[2][3] = -sina; \
matrix.coefs[3][2] = sina;

void AtkTilt_Ctor(AtkTilt* unit)
{
    ZERO_MATRIX
    unit->m_angle = IN0(4);
    AtkMatrix matrix = unit->matrix;
    FILL_TILT_MATRIX;
    unit->matrix = matrix;
    if(INRATE(4) == calc_FullRate)
	SETCALC(AtkTilt_next_a);
    else
	SETCALC(AtkTilt_next_k);
    AtkTilt_next_k(unit, 1); 
}

void AtkTilt_next_a(AtkTilt *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float *angle = IN(4);
    for(int i = 0; i < inNumSamples; i++){
	if(angle[i] != unit->m_angle){
	    unit->m_angle = angle[i];
	    FILL_TILT_MATRIX
	}
	CALC_MATRIX
    }
    unit->matrix = matrix;
    
}

void AtkTilt_next_k(AtkTilt *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float angle = IN0(4);
    float angleslope;
    if(angle != unit->m_angle){
	angleslope = CALCSLOPE(angle, unit->m_angle);
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	    unit->m_angle += angleslope;
	    FILL_TILT_MATRIX
	}	
    } else {
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	}
    }
    unit->matrix = matrix;
    unit->m_angle = angle;	
}

////////////////////// AtkTumble ///////////////////////
// uses 'angle' var. 

#define FILL_TUMBLE_MATRIX \
matrix.coefs[0][0] = matrix.coefs[2][2] = 1.; \
double cosa = cos(unit->m_angle); \
double sina = sin(unit->m_angle); \
matrix.coefs[1][1] = matrix.coefs[3][3] = cosa; \
matrix.coefs[1][3] = -sina; \
matrix.coefs[3][1] = sina;

void AtkTumble_Ctor(AtkTumble* unit)
{
    ZERO_MATRIX
    unit->m_angle = IN0(4);
    AtkMatrix matrix = unit->matrix;
    FILL_TUMBLE_MATRIX;
    unit->matrix = matrix;
    if(INRATE(4) == calc_FullRate)
	SETCALC(AtkTumble_next_a);
    else
	SETCALC(AtkTumble_next_k);
    AtkTumble_next_k(unit, 1); 
}

void AtkTumble_next_a(AtkTumble *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float *angle = IN(4);
    for(int i = 0; i < inNumSamples; i++){
	if(angle[i] != unit->m_angle){
	    unit->m_angle = angle[i];
	    FILL_TUMBLE_MATRIX
	}
	CALC_MATRIX
    }
    unit->matrix = matrix;
    
}

void AtkTumble_next_k(AtkTumble *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float angle = IN0(4);
    float angleslope;
    if(angle != unit->m_angle){
	angleslope = CALCSLOPE(angle, unit->m_angle);
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	    unit->m_angle += angleslope;
	    FILL_TUMBLE_MATRIX
	}	
    } else {
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	}
    }
    unit->matrix = matrix;
    unit->m_angle = angle;	
}

////////////////////// AtkFocusX ///////////////////////
// uses 'angle' var. 

#define FILL_FOCUSX_MATRIX \
double sina = sin(unit->m_angle); \
double roneplussinaa = 1. / (1. + sin(fabs(unit->m_angle))); \
double cosa = cos(unit->m_angle); \
double sina1sina = sina * roneplussinaa; \
matrix.coefs[0][0] = matrix.coefs[1][1] = roneplussinaa; \
matrix.coefs[2][2] = matrix.coefs[3][3] = cosa * roneplussinaa; \
matrix.coefs[0][1] = rsqrt2 * sina1sina; \
matrix.coefs[1][0] = sqrt2 * sina1sina; \

void AtkFocusX_Ctor(AtkFocusX* unit)
{
    ZERO_MATRIX
    unit->m_angle = IN0(4);
    AtkMatrix matrix = unit->matrix;
    
    FILL_FOCUSX_MATRIX
    
    unit->matrix = matrix;
    if(INRATE(4) == calc_FullRate)
	SETCALC(AtkFocusX_next_a);
    else
	SETCALC(AtkFocusX_next_k);
    AtkFocusX_next_k(unit, 1); 
}

void AtkFocusX_next_a(AtkFocusX *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float *angle = IN(4);
    for(int i = 0; i < inNumSamples; i++){
	if(angle[i] != unit->m_angle){
	    unit->m_angle = angle[i];	    
	    
	    FILL_FOCUSX_MATRIX
	    
	}
	CALC_MATRIX
    }
    unit->matrix = matrix;
    
}

void AtkFocusX_next_k(AtkFocusX *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float angle = IN0(4);
    float angleslope;
    if(angle != unit->m_angle){
	angleslope = CALCSLOPE(angle, unit->m_angle);
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	    unit->m_angle += angleslope;
	    FILL_FOCUSX_MATRIX
	}	
    } else {
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	}
    }
    unit->matrix = matrix;
    unit->m_angle = angle;	
}

////////////////////// AtkFocusY ///////////////////////
// uses 'angle' var. 

#define FILL_FOCUSY_MATRIX \
double sina = sin(unit->m_angle); \
double roneplussinaa = 1. / (1. + sin(fabs(unit->m_angle))); \
double cosa = cos(unit->m_angle); \
double sina1sina = sina * roneplussinaa; \
matrix.coefs[0][0] = matrix.coefs[2][2] = roneplussinaa; \
matrix.coefs[1][1] = matrix.coefs[3][3] = cosa * roneplussinaa; \
matrix.coefs[0][2] = rsqrt2 * sina1sina; \
matrix.coefs[2][0] = sqrt2 * sina1sina; \


void AtkFocusY_Ctor(AtkFocusY* unit)
{
    ZERO_MATRIX
    unit->m_angle = IN0(4);
    AtkMatrix matrix = unit->matrix;
    
    FILL_FOCUSY_MATRIX
    
    unit->matrix = matrix;
    if(INRATE(4) == calc_FullRate)
	SETCALC(AtkFocusY_next_a);
    else
	SETCALC(AtkFocusY_next_k);
    AtkFocusY_next_k(unit, 1); 
}

void AtkFocusY_next_a(AtkFocusY *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float *angle = IN(4);
    for(int i = 0; i < inNumSamples; i++){
	if(angle[i] != unit->m_angle){
	    unit->m_angle = angle[i];
	    FILL_FOCUSY_MATRIX
	}
	CALC_MATRIX
    }
    unit->matrix = matrix;
    
}

void AtkFocusY_next_k(AtkFocusY *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float angle = IN0(4);
    float angleslope;
    if(angle != unit->m_angle){
	angleslope = CALCSLOPE(angle, unit->m_angle);
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	    unit->m_angle += angleslope;
	    FILL_FOCUSY_MATRIX
	}	
    } else {
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	}
    }
    unit->matrix = matrix;
    unit->m_angle = angle;	
}

////////////////////// AtkFocusZ ///////////////////////
// uses 'angle' var. 

#define FILL_FOCUSZ_MATRIX \
double sina = sin(unit->m_angle); \
double roneplussinaa = 1. / (1. + sin(fabs(unit->m_angle))); \
double cosa = cos(unit->m_angle); \
double sina1sina = sina * roneplussinaa; \
matrix.coefs[0][0] = matrix.coefs[3][3] = roneplussinaa; \
matrix.coefs[1][1] = matrix.coefs[2][2] = cosa * roneplussinaa; \
matrix.coefs[0][3] = rsqrt2 * sina1sina; \
matrix.coefs[3][0] = sqrt2 * sina1sina; \

void AtkFocusZ_Ctor(AtkFocusZ* unit)
{
    ZERO_MATRIX
    unit->m_angle = IN0(4);
    AtkMatrix matrix = unit->matrix;
    FILL_FOCUSZ_MATRIX;
    unit->matrix = matrix;
    if(INRATE(4) == calc_FullRate)
	SETCALC(AtkFocusZ_next_a);
    else
	SETCALC(AtkFocusZ_next_k);
    AtkFocusZ_next_k(unit, 1); 
}

void AtkFocusZ_next_a(AtkFocusZ *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float *angle = IN(4);
    for(int i = 0; i < inNumSamples; i++){
	if(angle[i] != unit->m_angle){
	    unit->m_angle = angle[i];
	    FILL_FOCUSZ_MATRIX
	}
	CALC_MATRIX
    }
    unit->matrix = matrix;
    
}

void AtkFocusZ_next_k(AtkFocusZ *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float angle = IN0(4);
    float angleslope;
    if(angle != unit->m_angle){
	angleslope = CALCSLOPE(angle, unit->m_angle);
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	    unit->m_angle += angleslope;
	    FILL_FOCUSZ_MATRIX
	}	
    } else {
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	}
    }
    unit->matrix = matrix;
    unit->m_angle = angle;	
}

////////////////////// AtkDirectX /////////////////////
// uses 'angle' var. 

/*
#define FILL_DIRECTX_MATRIX \
double sq2cosa2 = sqrt2 * cos(unit->m_angle * 0.5); \
double sq2sina2 = sqrt2 * sin(unit->m_angle * 0.5); \
matrix.coefs[0][0] = matrix.coefs[2][2] = matrix.coefs[3][3] =  sq2cosa2; \
matrix.coefs[1][1] = sq2sina2; 
*/
#define FILL_DIRECTX_MATRIX \
matrix.coefs[0][0]= matrix.coefs[2][2] = matrix.coefs[3][3]  = sqrt(1.0 + sin(unit->m_angle)); \
matrix.coefs[1][1] = sqrt(1.0 - sin(unit->m_angle)); \


void AtkDirectX_Ctor(AtkDirectX* unit)
{
    ZERO_MATRIX
    unit->m_angle = IN0(4);
    AtkMatrix matrix = unit->matrix;
    FILL_DIRECTX_MATRIX;
    unit->matrix = matrix;
    if(INRATE(4) == calc_FullRate)
	SETCALC(AtkDirectX_next_a);
    else
	SETCALC(AtkDirectX_next_k);
    AtkDirectX_next_k(unit, 1); 
}

void AtkDirectX_next_a(AtkDirectX *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float *angle = IN(4);
    for(int i = 0; i < inNumSamples; i++){
	if(angle[i] != unit->m_angle){
	    unit->m_angle = angle[i];
	    FILL_DIRECTX_MATRIX
	}
	CALC_MATRIX
    }
    unit->matrix = matrix;
    
}

void AtkDirectX_next_k(AtkDirectX *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float angle = IN0(4);
    float angleslope;
    if(angle != unit->m_angle){
	angleslope = CALCSLOPE(angle, unit->m_angle);
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	    unit->m_angle += angleslope;
	    FILL_DIRECTX_MATRIX
	}	
    } else {
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	}
    }
    unit->matrix = matrix;
    unit->m_angle = angle;	
}

////////////////////// AtkDirectY /////////////////////
// uses 'angle' var. 

/*
#define FILL_DIRECTY_MATRIX \
double sq2cosa2 = sqrt2 * cos(unit->m_angle * 0.5); \
double sq2sina2 = sqrt2 * sin(unit->m_angle * 0.5); \
matrix.coefs[0][0] = matrix.coefs[1][1] = matrix.coefs[3][3] =  sq2cosa2; \
matrix.coefs[2][2] = sq2sina2; 
*/
#define FILL_DIRECTY_MATRIX \
matrix.coefs[0][0]= matrix.coefs[1][1] = matrix.coefs[3][3]  = sqrt(1.0 + sin(unit->m_angle)); \
matrix.coefs[2][2] = sqrt(1.0 - sin(unit->m_angle)); \

void AtkDirectY_Ctor(AtkDirectY* unit)
{
    ZERO_MATRIX
    unit->m_angle = IN0(4);
    AtkMatrix matrix = unit->matrix;
    FILL_DIRECTY_MATRIX;
    unit->matrix = matrix;
    if(INRATE(4) == calc_FullRate)
	SETCALC(AtkDirectY_next_a);
    else
	SETCALC(AtkDirectY_next_k);
    AtkDirectY_next_k(unit, 1); 
}

void AtkDirectY_next_a(AtkDirectY *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float *angle = IN(4);
    for(int i = 0; i < inNumSamples; i++){
	if(angle[i] != unit->m_angle){
	    unit->m_angle = angle[i];
	    FILL_DIRECTY_MATRIX
	}
	CALC_MATRIX
    }
    unit->matrix = matrix;
    
}

void AtkDirectY_next_k(AtkDirectY *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float angle = IN0(4);
    float angleslope;
    if(angle != unit->m_angle){
	angleslope = CALCSLOPE(angle, unit->m_angle);
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	    unit->m_angle += angleslope;
	    FILL_DIRECTY_MATRIX
	}	
    } else {
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	}
    }
    unit->matrix = matrix;
    unit->m_angle = angle;	
}

////////////////////// AtkDirectZ /////////////////////
// uses 'angle' var. 

/*
#define FILL_DIRECTZ_MATRIX \
double sq2cosa2 = sqrt2 * cos(unit->m_angle * 0.5); \
double sq2sina2 = sqrt2 * sin(unit->m_angle * 0.5); \
matrix.coefs[0][0] = matrix.coefs[1][1] = matrix.coefs[2][2] =  sq2cosa2; \
matrix.coefs[3][3] = sq2sina2; 
*/
#define FILL_DIRECTZ_MATRIX \
matrix.coefs[0][0]= matrix.coefs[1][1] = matrix.coefs[2][2]  = sqrt(1.0 + sin(unit->m_angle)); \
matrix.coefs[3][3] = sqrt(1.0 - sin(unit->m_angle)); \

void AtkDirectZ_Ctor(AtkDirectZ* unit)
{
    ZERO_MATRIX
    unit->m_angle = IN0(4);
    AtkMatrix matrix = unit->matrix;
    FILL_DIRECTZ_MATRIX;
    unit->matrix = matrix;
    if(INRATE(4) == calc_FullRate)
	SETCALC(AtkDirectZ_next_a);
    else
	SETCALC(AtkDirectZ_next_k);
    AtkDirectZ_next_k(unit, 1); 
}

void AtkDirectZ_next_a(AtkDirectZ *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float *angle = IN(4);
    for(int i = 0; i < inNumSamples; i++){
	if(angle[i] != unit->m_angle){
	    unit->m_angle = angle[i];
	    FILL_DIRECTZ_MATRIX
	}
	CALC_MATRIX
    }
    unit->matrix = matrix;
    
}

void AtkDirectZ_next_k(AtkDirectZ *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float angle = IN0(4);
    float angleslope;
    if(angle != unit->m_angle){
	angleslope = CALCSLOPE(angle, unit->m_angle);
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	    unit->m_angle += angleslope;
	    FILL_DIRECTZ_MATRIX
	}	
    } else {
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	}
    }
    unit->matrix = matrix;
    unit->m_angle = angle;	
}

////////////////////// AtkPushX ///////////////////////
// uses 'angle' var. 

#define FILL_PUSHX_MATRIX \
double cosa = cos(unit->m_angle); \
double cosa2 = cosa * cosa; \
double sina = sin(unit->m_angle); \
double asina = sin(fabs(unit->m_angle)); \
double sqrt2asinasina = sqrt2 * (asina * sina); \
matrix.coefs[0][0] = 1.f; \
matrix.coefs[1][1] = matrix.coefs[2][2] = matrix.coefs[3][3] = cosa2; \
matrix.coefs[1][0] = sqrt2asinasina; \

void AtkPushX_Ctor(AtkPushX* unit)
{
    ZERO_MATRIX
    unit->m_angle = IN0(4);
    AtkMatrix matrix = unit->matrix;
    FILL_PUSHX_MATRIX;
    unit->matrix = matrix;
    if(INRATE(4) == calc_FullRate)
	SETCALC(AtkPushX_next_a);
    else
	SETCALC(AtkPushX_next_k);
    AtkPushX_next_k(unit, 1); 
}

void AtkPushX_next_a(AtkPushX *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float *angle = IN(4);
    for(int i = 0; i < inNumSamples; i++){
	if(angle[i] != unit->m_angle){
	    unit->m_angle = angle[i];
	    FILL_PUSHX_MATRIX
	}
	CALC_MATRIX
    }
    unit->matrix = matrix;
    
}

void AtkPushX_next_k(AtkPushX *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float angle = IN0(4);
    float angleslope;
    if(angle != unit->m_angle){
	angleslope = CALCSLOPE(angle, unit->m_angle);
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	    unit->m_angle += angleslope;
	    FILL_PUSHX_MATRIX
	}	
    } else {
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	}
    }
    unit->matrix = matrix;
    unit->m_angle = angle;	
}

////////////////////// AtkPushY ///////////////////////
// uses 'angle' var.

#define FILL_PUSHY_MATRIX \
double cosa = cos(unit->m_angle); \
double cosa2 = cosa * cosa; \
double sina = sin(unit->m_angle); \
double asina = sin(fabs(unit->m_angle)); \
double sqrt2asinasina = sqrt2 * (asina * sina); \
matrix.coefs[0][0] = 1.f; \
matrix.coefs[1][1] = matrix.coefs[2][2] = matrix.coefs[3][3] = cosa2; \
matrix.coefs[2][0] = sqrt2asinasina; \



void AtkPushY_Ctor(AtkPushY* unit)
{
    ZERO_MATRIX
    unit->m_angle = IN0(4);
    AtkMatrix matrix = unit->matrix;
    FILL_PUSHY_MATRIX;
    unit->matrix = matrix;
    if(INRATE(4) == calc_FullRate)
	SETCALC(AtkPushY_next_a);
    else
	SETCALC(AtkPushY_next_k);
    AtkPushY_next_k(unit, 1); 
}

void AtkPushY_next_a(AtkPushY *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float *angle = IN(4);
    for(int i = 0; i < inNumSamples; i++){
	if(angle[i] != unit->m_angle){
	    unit->m_angle = angle[i];
	    FILL_PUSHY_MATRIX
	}
	CALC_MATRIX
    }
    unit->matrix = matrix;
    
}

void AtkPushY_next_k(AtkPushY *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float angle = IN0(4);
    float angleslope;
    if(angle != unit->m_angle){
	angleslope = CALCSLOPE(angle, unit->m_angle);
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	    unit->m_angle += angleslope;
	    FILL_PUSHY_MATRIX
	}	
    } else {
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	}
    }
    unit->matrix = matrix;
    unit->m_angle = angle;	
}

////////////////////// AtkPushZ ///////////////////////
// uses 'angle' var. 

#define FILL_PUSHZ_MATRIX \
double cosa = cos(unit->m_angle); \
double cosa2 = cosa * cosa; \
double sina = sin(unit->m_angle); \
double asina = sin(fabs(unit->m_angle)); \
double sqrt2asinasina = sqrt2 * (asina * sina); \
matrix.coefs[0][0] = 1.f; \
matrix.coefs[1][1] = matrix.coefs[2][2] = matrix.coefs[3][3] = cosa2; \
matrix.coefs[3][0] = sqrt2asinasina; \

void AtkPushZ_Ctor(AtkPushZ* unit)
{
    ZERO_MATRIX
    unit->m_angle = IN0(4);
    AtkMatrix matrix = unit->matrix;
    FILL_PUSHZ_MATRIX;
    unit->matrix = matrix;
    if(INRATE(4) == calc_FullRate)
	SETCALC(AtkPushZ_next_a);
    else
	SETCALC(AtkPushZ_next_k);
    AtkPushZ_next_k(unit, 1); 
}

void AtkPushZ_next_a(AtkPushZ *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float *angle = IN(4);
    for(int i = 0; i < inNumSamples; i++){
	if(angle[i] != unit->m_angle){
	    unit->m_angle = angle[i];
	    FILL_PUSHZ_MATRIX
	}
	CALC_MATRIX
    }
    unit->matrix = matrix;
    
}

void AtkPushZ_next_k(AtkPushZ *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float angle = IN0(4);
    float angleslope;
    if(angle != unit->m_angle){
	angleslope = CALCSLOPE(angle, unit->m_angle);
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	    unit->m_angle += angleslope;
	    FILL_PUSHZ_MATRIX
	}	
    } else {
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	}
    }
    unit->matrix = matrix;
    unit->m_angle = angle;	
}

////////////////////// AtkPressX ///////////////////////
// uses 'angle' var. 

#define FILL_PRESSX_MATRIX \
double cosa = cos(unit->m_angle); \
double cosa2 = cosa * cosa; \
double sina = sin(unit->m_angle); \
double asina = sin(fabs(unit->m_angle)); \
double sqrt2asinasina = sqrt2 * (asina * sina); \
matrix.coefs[0][0] = 1.f; \
matrix.coefs[1][0] = sqrt2asinasina; \
matrix.coefs[1][1] = cosa2; \
matrix.coefs[2][2] = matrix.coefs[3][3] = cosa; \


void AtkPressX_Ctor(AtkPressX* unit)
{
    ZERO_MATRIX
    unit->m_angle = IN0(4);
    AtkMatrix matrix = unit->matrix;
    FILL_PRESSX_MATRIX;
    unit->matrix = matrix;
    if(INRATE(4) == calc_FullRate)
	SETCALC(AtkPressX_next_a);
    else
	SETCALC(AtkPressX_next_k);
    AtkPressX_next_k(unit, 1); 
}

void AtkPressX_next_a(AtkPressX *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float *angle = IN(4);
    for(int i = 0; i < inNumSamples; i++){
	if(angle[i] != unit->m_angle){
	    unit->m_angle = angle[i];
	    FILL_PRESSX_MATRIX
	}
	CALC_MATRIX
    }
    unit->matrix = matrix;
    
}

void AtkPressX_next_k(AtkPressX *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float angle = IN0(4);
    float angleslope;
    if(angle != unit->m_angle){
	angleslope = CALCSLOPE(angle, unit->m_angle);
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	    unit->m_angle += angleslope;
	    FILL_PRESSX_MATRIX
	}	
    } else {
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	}
    }
    unit->matrix = matrix;
    unit->m_angle = angle;	
}

////////////////////// AtkPressY ///////////////////////
// uses 'angle' var. 


#define FILL_PRESSY_MATRIX \
double cosa = cos(unit->m_angle); \
double cosa2 = cosa * cosa; \
double sina = sin(unit->m_angle); \
double asina = sin(fabs(unit->m_angle)); \
double sqrt2asinasina = sqrt2 * (asina * sina); \
matrix.coefs[0][0] = 1.f; \
matrix.coefs[2][0] = sqrt2asinasina; \
matrix.coefs[2][2] = cosa2; \
matrix.coefs[1][1] = matrix.coefs[3][3] = cosa; \

void AtkPressY_Ctor(AtkPressY* unit)
{
    ZERO_MATRIX
    unit->m_angle = IN0(4);
    AtkMatrix matrix = unit->matrix;
    FILL_PRESSY_MATRIX;
    unit->matrix = matrix;
    if(INRATE(4) == calc_FullRate)
	SETCALC(AtkPressY_next_a);
    else
	SETCALC(AtkPressY_next_k);
    AtkPressY_next_k(unit, 1); 
}

void AtkPressY_next_a(AtkPressY *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float *angle = IN(4);
    for(int i = 0; i < inNumSamples; i++){
	if(angle[i] != unit->m_angle){
	    unit->m_angle = angle[i];
	    FILL_PRESSY_MATRIX
	}
	CALC_MATRIX
    }
    unit->matrix = matrix;
    
}

void AtkPressY_next_k(AtkPressY *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float angle = IN0(4);
    float angleslope;
    if(angle != unit->m_angle){
	angleslope = CALCSLOPE(angle, unit->m_angle);
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	    unit->m_angle += angleslope;
	    FILL_PRESSY_MATRIX
	}	
    } else {
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	}
    }
    unit->matrix = matrix;
    unit->m_angle = angle;	
}

////////////////////// AtkPressZ ///////////////////////
// uses 'angle' var. 

#define FILL_PRESSZ_MATRIX \
double cosa = cos(unit->m_angle); \
double cosa2 = cosa * cosa; \
double sina = sin(unit->m_angle); \
double asina = sin(fabs(unit->m_angle)); \
double sqrt2asinasina = sqrt2 * (asina * sina); \
matrix.coefs[0][0] = 1.f; \
matrix.coefs[3][0] = sqrt2asinasina; \
matrix.coefs[3][3] = cosa2; \
matrix.coefs[1][1] = matrix.coefs[2][2] = cosa; \

void AtkPressZ_Ctor(AtkPressZ* unit)
{
    ZERO_MATRIX
    unit->m_angle = IN0(4);
    AtkMatrix matrix = unit->matrix;
    FILL_PRESSZ_MATRIX;
    unit->matrix = matrix;
    if(INRATE(4) == calc_FullRate)
	SETCALC(AtkPressZ_next_a);
    else
	SETCALC(AtkPressZ_next_k);
    AtkPressZ_next_k(unit, 1); 
}

void AtkPressZ_next_a(AtkPressZ *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float *angle = IN(4);
    for(int i = 0; i < inNumSamples; i++){
	if(angle[i] != unit->m_angle){
	    unit->m_angle = angle[i];
	    FILL_PRESSZ_MATRIX
	}
	CALC_MATRIX
    }
    unit->matrix = matrix;
    
}

void AtkPressZ_next_k(AtkPressZ *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float angle = IN0(4);
    float angleslope;
    if(angle != unit->m_angle){
	angleslope = CALCSLOPE(angle, unit->m_angle);
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	    unit->m_angle += angleslope;
	    FILL_PRESSZ_MATRIX
	}	
    } else {
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	}
    }
    unit->matrix = matrix;
    unit->m_angle = angle;	
}

////////////////////// AtkZoomX ///////////////////////
// uses 'angle' var. 

#define FILL_ZOOMX_MATRIX \
double cosa = cos(unit->m_angle); \
double sina = sin(unit->m_angle); \
matrix.coefs[0][0] = matrix.coefs[1][1] = 1.; \
matrix.coefs[2][2] = matrix.coefs[3][3] = cosa; \
matrix.coefs[0][1] = rsqrt2 * sina; \
matrix.coefs[1][0] = sqrt2 * sina;

void AtkZoomX_Ctor(AtkZoomX* unit)
{
    ZERO_MATRIX
    unit->m_angle = IN0(4);
    AtkMatrix matrix = unit->matrix;
    FILL_ZOOMX_MATRIX;
    unit->matrix = matrix;
    if(INRATE(4) == calc_FullRate)
	SETCALC(AtkZoomX_next_a);
    else
	SETCALC(AtkZoomX_next_k);
    AtkZoomX_next_k(unit, 1); 
}

void AtkZoomX_next_a(AtkZoomX *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float *angle = IN(4);
    for(int i = 0; i < inNumSamples; i++){
	if(angle[i] != unit->m_angle){
	    unit->m_angle = angle[i];
	    FILL_ZOOMX_MATRIX
	}
	CALC_MATRIX
    }
    unit->matrix = matrix;
    
}

void AtkZoomX_next_k(AtkZoomX *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float angle = IN0(4);
    float angleslope;
    if(angle != unit->m_angle){
	angleslope = CALCSLOPE(angle, unit->m_angle);
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	    unit->m_angle += angleslope;
	    FILL_ZOOMX_MATRIX
	}	
    } else {
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	}
    }
    unit->matrix = matrix;
    unit->m_angle = angle;	
}

////////////////////// AtkZoomY ///////////////////////
// uses 'angle' var. 

#define FILL_ZOOMY_MATRIX \
double cosa = cos(unit->m_angle); \
double sina = sin(unit->m_angle); \
matrix.coefs[0][0] = matrix.coefs[2][2] = 1.; \
matrix.coefs[1][1] = matrix.coefs[3][3] = cosa; \
matrix.coefs[0][2] = rsqrt2 * sina; \
matrix.coefs[2][0] = sqrt2 * sina; \


void AtkZoomY_Ctor(AtkZoomY* unit)
{
    ZERO_MATRIX
    unit->m_angle = IN0(4);
    AtkMatrix matrix = unit->matrix;
    FILL_ZOOMY_MATRIX;
    unit->matrix = matrix;
    if(INRATE(4) == calc_FullRate)
	SETCALC(AtkZoomY_next_a);
    else
	SETCALC(AtkZoomY_next_k);
    AtkZoomY_next_k(unit, 1); 
}

void AtkZoomY_next_a(AtkZoomY *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float *angle = IN(4);
    for(int i = 0; i < inNumSamples; i++){
	if(angle[i] != unit->m_angle){
	    unit->m_angle = angle[i];
	    FILL_ZOOMY_MATRIX
	}
	CALC_MATRIX
    }
    unit->matrix = matrix;
    
}

void AtkZoomY_next_k(AtkZoomY *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float angle = IN0(4);
    float angleslope;
    if(angle != unit->m_angle){
	angleslope = CALCSLOPE(angle, unit->m_angle);
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	    unit->m_angle += angleslope;
	    FILL_ZOOMY_MATRIX
	}	
    } else {
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	}
    }
    unit->matrix = matrix;
    unit->m_angle = angle;	
}

////////////////////// AtkZoomZ ///////////////////////
// uses 'angle' var. 

#define FILL_ZOOMZ_MATRIX \
double cosa = cos(unit->m_angle); \
double sina = sin(unit->m_angle); \
matrix.coefs[0][0] = matrix.coefs[3][3] = 1.; \
matrix.coefs[1][1] = matrix.coefs[2][2] = cosa; \
matrix.coefs[0][3] = rsqrt2 * sina; \
matrix.coefs[3][0] = sqrt2 * sina; \

void AtkZoomZ_Ctor(AtkZoomZ* unit)
{
    ZERO_MATRIX
    unit->m_angle = IN0(4);
    AtkMatrix matrix = unit->matrix;
    FILL_ZOOMZ_MATRIX;
    unit->matrix = matrix;
    if(INRATE(4) == calc_FullRate)
	SETCALC(AtkZoomZ_next_a);
    else
	SETCALC(AtkZoomZ_next_k);
    AtkZoomZ_next_k(unit, 1); 
}

void AtkZoomZ_next_a(AtkZoomZ *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float *angle = IN(4);
    for(int i = 0; i < inNumSamples; i++){
	if(angle[i] != unit->m_angle){
	    unit->m_angle = angle[i];
	    FILL_ZOOMZ_MATRIX
	}
	CALC_MATRIX
    }
    unit->matrix = matrix;
    
}

void AtkZoomZ_next_k(AtkZoomZ *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float angle = IN0(4);
    float angleslope;
    if(angle != unit->m_angle){
	angleslope = CALCSLOPE(angle, unit->m_angle);
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	    unit->m_angle += angleslope;
	    FILL_ZOOMZ_MATRIX
	}	
    } else {
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	}
    }
    unit->matrix = matrix;
    unit->m_angle = angle;	
}

////////////////////// AtkDominateX ///////////////////////
// uses 'angle' var. 

#define FILL_DOMINATEX_MATRIX \
double dominate = pow(10, unit->m_gain * 0.05); \
double rsqrt8 = 1./sqrt(8.); \
double rdom = 1./dominate; \
double hdomprdom = 0.5 * (dominate + rdom); \
double dommrdom = dominate - rdom; \
matrix.coefs[0][0] = matrix.coefs[1][1] = hdomprdom; \
matrix.coefs[0][1] = rsqrt8 * dommrdom; \
matrix.coefs[1][0] = rsqrt2 * dommrdom; \
matrix.coefs[2][2] = matrix.coefs[3][3] = 1.; \

void AtkDominateX_Ctor(AtkDominateX* unit)
{
    ZERO_MATRIX
    unit->m_gain = IN0(4);
    AtkMatrix matrix = unit->matrix;
    FILL_DOMINATEX_MATRIX;
    unit->matrix = matrix;
    if(INRATE(4) == calc_FullRate)
	SETCALC(AtkDominateX_next_a);
    else
	SETCALC(AtkDominateX_next_k);
    AtkDominateX_next_k(unit, 1); 
}

void AtkDominateX_next_a(AtkDominateX *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float *gain = IN(4);
    for(int i = 0; i < inNumSamples; i++){
	if(gain[i] != unit->m_gain){
	    unit->m_gain = gain[i];
	    FILL_DOMINATEX_MATRIX
	}
	CALC_MATRIX
    }
    unit->matrix = matrix;
    
}

void AtkDominateX_next_k(AtkDominateX *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float gain = IN0(4);
    float gainslope;
    if(gain != unit->m_gain){
	gainslope = CALCSLOPE(gain, unit->m_gain);
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	    unit->m_gain += gainslope;
	    FILL_DOMINATEX_MATRIX
	}	
    } else {
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	}
    }
    unit->matrix = matrix;
    unit->m_gain = gain;	
}

////////////////////// AtkDominateY ///////////////////////
// uses 'angle' var. 

#define FILL_DOMINATEY_MATRIX \
double dominate = pow(10, unit->m_gain * 0.05); \
double rsqrt8 = 1./sqrt(8.); \
double rdom = 1./dominate; \
double hdomprdom = 0.5 * (dominate + rdom); \
double dommrdom = dominate - rdom; \
matrix.coefs[0][0] = matrix.coefs[2][2] = hdomprdom; \
matrix.coefs[0][2] = rsqrt8 * dommrdom; \
matrix.coefs[2][0] = rsqrt2 * dommrdom; \
matrix.coefs[1][1] = matrix.coefs[3][3] = 1.; \

void AtkDominateY_Ctor(AtkDominateY* unit)
{
    ZERO_MATRIX
    unit->m_gain = IN0(4);
    AtkMatrix matrix = unit->matrix;
    FILL_DOMINATEY_MATRIX;
    unit->matrix = matrix;
    if(INRATE(4) == calc_FullRate)
	SETCALC(AtkDominateY_next_a);
    else
	SETCALC(AtkDominateY_next_k);
    AtkDominateY_next_k(unit, 1); 
}

void AtkDominateY_next_a(AtkDominateY *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float *gain = IN(4);
    for(int i = 0; i < inNumSamples; i++){
	if(gain[i] != unit->m_gain){
	    unit->m_gain = gain[i];
	    FILL_DOMINATEY_MATRIX
	}
	CALC_MATRIX
    }
    unit->matrix = matrix;
    
}

void AtkDominateY_next_k(AtkDominateY *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float gain = IN0(4);
    float gainslope;
    if(gain != unit->m_gain){
	gainslope = CALCSLOPE(gain, unit->m_gain);
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	    unit->m_gain += gainslope;
	    FILL_DOMINATEY_MATRIX
	}	
    } else {
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	}
    }
    unit->matrix = matrix;
    unit->m_gain = gain;	
}

////////////////////// AtkDominateZ ///////////////////////
// uses 'angle' var. 

#define FILL_DOMINATEZ_MATRIX \
double dominate = pow(10, unit->m_gain * 0.05); \
double rsqrt8 = 1./sqrt(8.); \
double rdom = 1./dominate; \
double hdomprdom = 0.5 * (dominate + rdom); \
double dommrdom = dominate - rdom; \
matrix.coefs[0][0] = matrix.coefs[3][3] = hdomprdom; \
matrix.coefs[0][3] = rsqrt8 * dommrdom; \
matrix.coefs[3][0] = rsqrt2 * dommrdom; \
matrix.coefs[1][1] = matrix.coefs[2][2] = 1.; \

void AtkDominateZ_Ctor(AtkDominateZ* unit)
{
    ZERO_MATRIX
    unit->m_gain = IN0(4);
    AtkMatrix matrix = unit->matrix;
    FILL_DOMINATEZ_MATRIX;
    unit->matrix = matrix;
    if(INRATE(4) == calc_FullRate)
	SETCALC(AtkDominateZ_next_a);
    else
	SETCALC(AtkDominateZ_next_k);
    AtkDominateZ_next_k(unit, 1); 
}

void AtkDominateZ_next_a(AtkDominateZ *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float *gain = IN(4);
    for(int i = 0; i < inNumSamples; i++){
	if(gain[i] != unit->m_gain){
	    unit->m_gain = gain[i];
	    FILL_DOMINATEZ_MATRIX
	}
	CALC_MATRIX
    }
    unit->matrix = matrix;
    
}

void AtkDominateZ_next_k(AtkDominateZ *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float gain = IN0(4);
    float gainslope;
    if(gain != unit->m_gain){
	gainslope = CALCSLOPE(gain, unit->m_gain);
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	    unit->m_gain += gainslope;
	    FILL_DOMINATEZ_MATRIX
	}	
    } else {
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	}
    }
    unit->matrix = matrix;
    unit->m_gain = gain;	
}



////////////////////// AtkDirectO ///////////////////////
// uses 'angle' var

/*
#define FILL_DIRECT_MATRIX \
matrix.coefs[0][0] = sqrt2 * cos(unit->m_angle * 0.5); \
matrix.coefs[1][1] = matrix.coefs[2][2] = matrix.coefs[3][3] = sqrt2 * sin(unit->m_angle * 0.5); \
*/
#define FILL_DIRECT_MATRIX \
matrix.coefs[0][0] = sqrt(1.0 + sin(unit->m_angle)); \
matrix.coefs[1][1] = matrix.coefs[2][2] = matrix.coefs[3][3] = sqrt(1.0 - sin(unit->m_angle)); \

void AtkDirectO_Ctor(AtkDirectO* unit)
{
    ZERO_MATRIX
    unit->m_angle = IN0(4);
    AtkMatrix matrix = unit->matrix;
    FILL_DIRECT_MATRIX;
    unit->matrix = matrix;
    if(INRATE(4) == calc_FullRate)
	SETCALC(AtkDirectO_next_a);
    else
	SETCALC(AtkDirectO_next_k);
    AtkDirectO_next_k(unit, 1); 
}

void AtkDirectO_next_a(AtkDirectO *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float *angle = IN(4);
    
    for(int i = 0; i < inNumSamples; i++){
	if(angle[i] != unit->m_angle){
	    unit->m_angle = angle[i];
	    FILL_DIRECT_MATRIX
	}
	CALC_MATRIX
    }
    unit->matrix = matrix;
}


void AtkDirectO_next_k(AtkDirectO *unit, int inNumSamples)
{
    SETUP_TRANSFORMS
    float angle = IN0(4);
    float angleslope;
    
    if(angle != unit->m_angle){
	angleslope = CALCSLOPE(angle, unit->m_angle);
	FILL_DIRECT_MATRIX	
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	    unit->m_angle += angleslope;
	    FILL_DIRECT_MATRIX
	}	
    } else {
	for(int i = 0; i < inNumSamples; i++){
	    CALC_MATRIX
	}
    }
    
    unit->matrix = matrix;
    unit->m_angle = angle;	
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
// AtkNFC - 
void AtkNFC_Ctor(AtkNFC* unit)
{
    unit->m_y1x = 0.0f;
    unit->m_y1y = 0.0f;
    unit->m_y1z = 0.0f;
    unit->m_distanceStart = IN0(4);
    if (INRATE(4) == calc_FullRate) {
	SETCALC(AtkNFC_next_a);
    } else {
	SETCALC(AtkNFC_next_k);
    };
    ClearUnitOutputs(unit, 1);
}

void AtkNFC_next_k(AtkNFC *unit, int inNumSamples)
{       
    float *Wout = OUT(0);
    float *Xout = OUT(1);
    float *Yout = OUT(2);
    float *Zout = OUT(3);
    
    float *Win = IN(0);
    float *Xin = IN(1);
    float *Yin = IN(2);
    float *Zin = IN(3);
    float distanceEnd = IN0(4);
    float distanceStart = unit->m_distanceStart;
    
    float distanceInc = CALCSLOPE(distanceEnd, distanceStart);
    
    float y1x = unit->m_y1x;
    float y1y = unit->m_y1y;
    float y1z = unit->m_y1z;
    
    for(int i = 0; i < inNumSamples; i++){
	float freq = 53.0 / distanceStart;
	float wc = (twopi * freq) * SAMPLEDUR;
	//	a0 = (1 + (wc.cos.neg * 2 + 2).sqrt).reciprocal;
	float a0 = 1 / (sqrt((cos(wc) * -2) + 2) + 1);
	float y0x = Xin[i] + a0 * y1x;
	Xout[i] = a0 * y0x + -a0 * y1x;
	y1x = y0x;
	float y0y = Yin[i] + a0 * y1y;
	Yout[i] = a0 * y0y + -a0 * y1y;
	y1y = y0y;
	float y0z = Zin[i] + a0 * y1z;
	Zout[i] = a0 * y0z + -a0 * y1z;
	y1z = y0z;
	// W is passed straight out
	Wout[i] = Win[i];
	distanceStart += distanceInc;
    }
    
    unit->m_y1x = zapgremlins(y1x);
    unit->m_y1y = zapgremlins(y1y);
    unit->m_y1z = zapgremlins(y1z);
    unit->m_distanceStart = distanceEnd;
    
}

void AtkNFC_next_a(AtkNFC *unit, int inNumSamples)
{       
    float *Wout = OUT(0);
    float *Xout = OUT(1);
    float *Yout = OUT(2);
    float *Zout = OUT(3);
    
    float *Win = IN(0);
    float *Xin = IN(1);
    float *Yin = IN(2);
    float *Zin = IN(3);
    float *distance = IN(4);
    
    float y1x = unit->m_y1x;
    float y1y = unit->m_y1y;
    float y1z = unit->m_y1z;
    
    for(int i = 0; i < inNumSamples; i++){
	float freq = 53.0 / distance[i];
	float wc = (twopi * freq) * SAMPLEDUR;
	float a0 = 1 / (sqrt((cos(wc) * -2) + 2) + 1);
	float y0x = Xin[i] + a0 * y1x;
	Xout[i] = a0 * y0x + -a0 * y1x;
	y1x = y0x;
	float y0y = Yin[i] + a0 * y1y;
	Yout[i] = a0 * y0y + -a0 * y1y;
	y1y = y0y;
	float y0z = Zin[i] + a0 * y1z;
	Zout[i] = a0 * y0z + -a0 * y1z;
	y1z = y0z;
	// W is passed straight out
	Wout[i] = Win[i];
    }
    
    unit->m_y1x = zapgremlins(y1x);
    unit->m_y1y = zapgremlins(y1y);
    unit->m_y1z = zapgremlins(y1z);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// AtkProximity

void AtkProximity_Ctor(AtkProximity* unit)
{
    unit->m_y1x = 0.0f;
    unit->m_y1y = 0.0f;
    unit->m_y1z = 0.0f;
    unit->m_distanceStart = IN0(4);
    if (INRATE(4) == calc_FullRate) {
	SETCALC(AtkProximity_next_a);
    } else {
	SETCALC(AtkProximity_next_k);
    };
    ClearUnitOutputs(unit, 1);
}

void AtkProximity_next_k(AtkProximity *unit, int inNumSamples)
{       
    float *Wout = OUT(0);
    float *Xout = OUT(1);
    float *Yout = OUT(2);
    float *Zout = OUT(3);
    
    float *Win = IN(0);
    float *Xin = IN(1);
    float *Yin = IN(2);
    float *Zin = IN(3);
    float distanceEnd = IN0(4);
    float distanceStart = unit->m_distanceStart;
    
    float distanceInc = CALCSLOPE(distanceEnd, distanceStart);
    
    float y1x = unit->m_y1x;
    float y1y = unit->m_y1y;
    float y1z = unit->m_y1z;
    
    for(int i=0; i<inNumSamples;i++){
	float freq = 53.0 / distanceStart;
	float wc = (twopi * freq) * SAMPLEDUR;
	//	a0 = 1 + (wc.cos.neg * 2 + 2).sqrt;
	float a0 = 1 + sqrt((cos(wc) * -2) + 2);
	float y0x = Xin[i] + y1x;
	Xout[i] = a0 * y0x - y1x;
	y1x = y0x;
	float y0y = Yin[i] + y1y;
	Yout[i] = a0 * y0y - y1y;
	y1y = y0y;
	float y0z = Zin[i] + y1z;
	Zout[i] = a0 * y0z - y1z;
	y1z = y0z;
	// W is passed straight out
	Wout[i] = Win[i];
	distanceStart += distanceInc;
    }
    
    unit->m_y1x = zapgremlins(y1x);
    unit->m_y1y = zapgremlins(y1y);
    unit->m_y1z = zapgremlins(y1z);
    unit->m_distanceStart = distanceEnd;
    
}

void AtkProximity_next_a(AtkProximity *unit, int inNumSamples)
{       
    float *Wout = OUT(0);
    float *Xout = OUT(1);
    float *Yout = OUT(2);
    float *Zout = OUT(3);
    
    float *Win = IN(0);
    float *Xin = IN(1);
    float *Yin = IN(2);
    float *Zin = IN(3);
    float *distance = IN(4);
    
    float y1x = unit->m_y1x;
    float y1y = unit->m_y1y;
    float y1z = unit->m_y1z;
    
    for(int i = 0; i<inNumSamples; i++){
	float freq = 53.0 / distance[i];
	float wc = (twopi * freq) * SAMPLEDUR;
	float a0 = 1 + sqrt((cos(wc) * -2) + 2);
	float y0x = Xin[i] + y1x;
	Xout[i] = a0 * y0x - y1x;
	y1x = y0x;
	float y0y = Yin[i] + y1y;
	Yout[i] = a0 * y0y - y1y;
	y1y = y0y;
	float y0z = Zin[i] + y1z;
	Zout[i] = a0 * y0z - y1z;
	y1z = y0z;
	// W is passed straight out
	Wout[i] = Win[i];
    }
    
    unit->m_y1x = zapgremlins(y1x);
    unit->m_y1y = zapgremlins(y1y);
    unit->m_y1z = zapgremlins(y1z);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
// AtkPsychoShelf

void AtkPsychoShelf_Ctor(AtkPsychoShelf* unit)
{
    unit->m_y1w = 0.0f;
    unit->m_y2w = 0.0f;
    unit->m_y1x = 0.0f;
    unit->m_y2x = 0.0f;
    unit->m_y1y = 0.0f;
    unit->m_y2y = 0.0f;
    unit->m_y1z = 0.0f;
    unit->m_y2z = 0.0f;
    unit->m_freq = IN0(4);
    unit->m_k0 = IN0(5);
    unit->m_k1 = IN0(6);
    if (INRATE(4) == calc_FullRate || INRATE(5) == calc_FullRate || INRATE(6) == calc_FullRate) {
		SETCALC(AtkPsychoShelf_next_a);
    } else {
		SETCALC(AtkPsychoShelf_next_k);
    };
    ClearUnitOutputs(unit, 1);
}

void AtkPsychoShelf_next_k(AtkPsychoShelf *unit, int inNumSamples)
{       
    float *Wout = OUT(0);
    float *Xout = OUT(1);
    float *Yout = OUT(2);
    float *Zout = OUT(3);
    
    float *Win = IN(0);
    float *Xin = IN(1);
    float *Yin = IN(2);
    float *Zin = IN(3);
	float freqEnd = IN0(4);
	float k0End = IN0(5);
	float k1End = IN0(6);

	float freqStart = unit->m_freq;
	float k0Start = unit->m_k0;
	float k1Start = unit->m_k1;

	float freqInc = CALCSLOPE(freqEnd, freqStart);
	float k0Inc = CALCSLOPE(k0End, k0Start);
	float k1Inc = CALCSLOPE(k1End, k1Start);
    
	float y0w;
    float y1w = unit->m_y1w;
    float y2w = unit->m_y2w;
	float y0x;
    float y1x = unit->m_y1x;
    float y2x = unit->m_y2x;
	float y0y;
    float y1y = unit->m_y1y;
    float y2y = unit->m_y2y;
	float y0z;
    float y1z = unit->m_y1z;
    float y2z = unit->m_y2z;
    
    for(int i=0; i<inNumSamples;i++){

		// wc = (pi * frequency / Server.default.sampleRate).tan
		// c = (wc - 1) / (wc + 1)
		float wc = tan(pi * freqStart * SAMPLEDUR);
		float c = (wc - 1.0) / (wc + 1.0);

		// a0 = (((1 - k2)/4) * (1 + (c**2))) + (((1 + k2)/2) * c);
		// a1 = ((1 - k2) * c) + (((1 + k2)/2) * (1 + (c**2)));
		// a2 = a0;
		//
		// b1 = Array.fill( k2.size, { (2*c).neg } );
		// b2 = Array.fill( k2.size, { (c**2).neg } );
		
		// 0 order coefficients
		float a0w = (((1.0 - k0Start) * 0.25) * (1.0 + (c*c))) + (((1.0 + k0Start) * 0.5) * c);
		float a1w = ((1.0 - k0Start) * c) + (((1.0 + k0Start) * 0.5) * (1.0 + (c*c)));
		float a2w = a0w;
		float b1w = -2.0 * c;
		float b2w = -c * c;

		// 1st order coefficients
		float a0xyz = (((1.0 - k1Start) * 0.25) * (1.0 + (c*c))) + (((1.0 + k1Start) * 0.5) * c);
		float a1xyz = ((1.0 - k1Start) * c) + (((1.0 + k1Start) * 0.5) * (1.0 + (c*c)));
		float a2xyz = a0xyz;
		float b1xyz = -2.0 * c;
		float b2xyz = -c * c;
		

		// filter 0 order
		y0w = Win[i] + b1w * y1w + b2w * y2w;
		Wout[i] = a0w * y0w + a1w * y1w + a2w * y2w;
		y2w = y1w;
		y1w = y0w;
		
		
		// filter 1st order
		y0x = Xin[i] + b1xyz * y1x + b2xyz * y2x;
		Xout[i] = a0xyz * y0x + a1xyz * y1x + a2xyz * y2x;
		y2x = y1x;
		y1x = y0x;

		y0y = Yin[i] + b1xyz * y1y + b2xyz * y2y;
		Yout[i] = a0xyz * y0y + a1xyz * y1y + a2xyz * y2y;
		y2y = y1y;
		y1y = y0y;
		
		y0z = Zin[i] + b1xyz * y1z + b2xyz * y2z;
		Zout[i] = a0xyz * y0z + a1xyz * y1z + a2xyz * y2z;
		y2z = y1z;
		y1z = y0z;
		
		
		freqStart += freqInc;
		k0Start += k0Inc;
		k1Start += k1Inc;
    }
    
    unit->m_y1w = zapgremlins(y1w);
    unit->m_y2w = zapgremlins(y2w);
    unit->m_y1x = zapgremlins(y1x);
    unit->m_y2x = zapgremlins(y2x);
    unit->m_y1y = zapgremlins(y1y);
    unit->m_y2y = zapgremlins(y2y);
    unit->m_y1z = zapgremlins(y1z);
    unit->m_y2z = zapgremlins(y2z);

    unit->m_freq = freqEnd;
    unit->m_k0 = k0End;
    unit->m_k1 = k1End;
    
}

void AtkPsychoShelf_next_a(AtkPsychoShelf *unit, int inNumSamples)
{       
    float *Wout = OUT(0);
    float *Xout = OUT(1);
    float *Yout = OUT(2);
    float *Zout = OUT(3);
    
    float *Win = IN(0);
    float *Xin = IN(1);
    float *Yin = IN(2);
    float *Zin = IN(3);
	float *freq = IN(4);
	float *k0 = IN(5);
	float *k1 = IN(6);
    
	float y0w;
    float y1w = unit->m_y1w;
    float y2w = unit->m_y2w;
	float y0x;
    float y1x = unit->m_y1x;
    float y2x = unit->m_y2x;
	float y0y;
    float y1y = unit->m_y1y;
    float y2y = unit->m_y2y;
	float y0z;
    float y1z = unit->m_y1z;
    float y2z = unit->m_y2z;
    
    for(int i = 0; i<inNumSamples; i++){
		// wc = (pi * frequency / Server.default.sampleRate).tan
		// c = (wc - 1) / (wc + 1)
		float wc = tan(pi * freq[i] * SAMPLEDUR);
		float c = (wc - 1.0) / (wc + 1.0);
		
		// a0 = (((1 - k2)/4) * (1 + (c**2))) + (((1 + k2)/2) * c);
		// a1 = ((1 - k2) * c) + (((1 + k2)/2) * (1 + (c**2)));
		// a2 = a0;
		//
		// b1 = Array.fill( k2.size, { (2*c).neg } );
		// b2 = Array.fill( k2.size, { (c**2).neg } );
		
		// 0 order coefficients
		float a0w = (((1.0 - k0[i]) * 0.25) * (1.0 + (c*c))) + (((1.0 + k0[i]) * 0.5) * c);
		float a1w = ((1.0 - k0[i]) * c) + (((1.0 + k0[i]) * 0.5) * (1.0 + (c*c)));
		float a2w = a0w;
		float b1w = -2.0 * c;
		float b2w = -c * c;
		
		// 1st order coefficients
		float a0xyz = (((1.0 - k1[i]) * 0.25) * (1.0 + (c*c))) + (((1.0 + k1[i]) * 0.5) * c);
		float a1xyz = ((1.0 - k1[i]) * c) + (((1.0 + k1[i]) * 0.5) * (1.0 + (c*c)));
		float a2xyz = a0xyz;
		float b1xyz = -2.0 * c;
		float b2xyz = -c * c;
		
		
		// filter 0 order
		y0w = Win[i] + b1w * y1w + b2w * y2w;
		Wout[i] = a0w * y0w + a1w * y1w + a2w * y2w;
		y2w = y1w;
		y1w = y0w;
		
		
		// filter 1st order
		y0x = Xin[i] + b1xyz * y1x + b2xyz * y2x;
		Xout[i] = a0xyz * y0x + a1xyz * y1x + a2xyz * y2x;
		y2x = y1x;
		y1x = y0x;
		
		y0y = Yin[i] + b1xyz * y1y + b2xyz * y2y;
		Yout[i] = a0xyz * y0y + a1xyz * y1y + a2xyz * y2y;
		y2y = y1y;
		y1y = y0y;
		
		y0z = Zin[i] + b1xyz * y1z + b2xyz * y2z;
		Zout[i] = a0xyz * y0z + a1xyz * y1z + a2xyz * y2z;
		y2z = y1z;
		y1z = y0z;
    }
    
    unit->m_y1w = zapgremlins(y1w);
    unit->m_y2w = zapgremlins(y2w);
    unit->m_y1x = zapgremlins(y1x);
    unit->m_y2x = zapgremlins(y2x);
    unit->m_y1y = zapgremlins(y1y);
    unit->m_y2y = zapgremlins(y2y);
    unit->m_y1z = zapgremlins(y1z);
    unit->m_y2z = zapgremlins(y2z);
}


/*
 ///////////////////////////////////////////////////////////////////////////////////////////////////////
 // AtoB - a, b, c, d in, w, x, y, z out
 
 void AtoB_Ctor(AtoB* unit)
 {
 float orientation = ZIN0(4);
 float mode = ZIN0(5);
 // set-up the proper encoding function
 switch ((int)orientation) {
 case 0 : 
 switch ((int)mode) {
 case 0: SETCALC(AtoB_next_00); break;
 //		    case 1: SETCALC(AtoB_next_01); break;
 //		    case 2: SETCALC(AtoB_next_02); break;
 default : break;
 }
 break;
 case 1 : 
 switch ((int)mode) {
 case 0: SETCALC(AtoB_next_10); break;
 //		    case 1: SETCALC(AtoB_next_11); break;
 //		    case 2: SETCALC(AtoB_next_12); break;
 default : break;
 }
 break;
 case 2 : 
 switch ((int)mode) {
 case 0: SETCALC(AtoB_next_20); break;
 //		    case 1: SETCALC(AtoB_next_21); break;
 //		    case 2: SETCALC(AtoB_next_22); break;
 default : break;
 }
 break;
 case 3 : 
 switch ((int)mode) {
 case 0: SETCALC(AtoB_next_30); break;
 //		    case 1: SETCALC(AtoB_next_31); break;
 //		    case 2: SETCALC(AtoB_next_32); break;
 default : break;
 }
 break;
 case 4 : 
 switch ((int)mode) {
 case 0: SETCALC(AtoB_next_40); break;
 //		    case 1: SETCALC(AtoB_next_41); break;
 //		    case 2: SETCALC(AtoB_next_42); break;
 default : break;
 }
 break;
 case 5 : 
 switch ((int)mode) {
 case 0: SETCALC(AtoB_next_50); break;
 //		    case 1: SETCALC(AtoB_next_51); break;
 //		    case 2: SETCALC(AtoB_next_52); break;
 default : break;
 }
 break;
 case 6 :
 switch ((int)mode) {
 case 0: SETCALC(AtoB_next_60); break;
 //		    case 1: SETCALC(AtoB_next_61); break;
 //		    case 2: SETCALC(AtoB_next_62); break;
 default : break;
 }
 break;
 case 7 : 
 switch ((int)mode) {
 case 0: SETCALC(AtoB_next_70); break;
 //		    case 1: SETCALC(AtoB_next_71); break;
 //		    case 2: SETCALC(AtoB_next_72); break;
 default : break;
 }
 break;
 default : break;
 }
 }
 
 #define SETUP_AtoB \
 float* a_in = ZIN(0); \
 float* b_in = ZIN(1); \
 float* c_in = ZIN(2); \
 float* d_in = ZIN(3); \
 float* w_out = ZOUT(0); \
 float* x_out = ZOUT(1); \
 float* y_out = ZOUT(2); \
 float* z_out = ZOUT(3); \
 float a, b, c, d;
 
 #define ASSIGN_ABCD \
 a = ZXP(a_in); \
 b = ZXP(b_in); \
 c = ZXP(c_in); \
 d = ZXP(d_in); 
 
 //some useful SC constants
 #ifndef __FP__
 const double pi     = acos(-1.);
 #endif
 const double pi2    = pi * .5;
 const double pi32   = pi * 1.5;
 const double twopi  = pi * 2.;
 const double rtwopi = 1. / twopi;
 const double log001 = log(0.001);
 const double log01  = log(0.01);
 const double log1   = log(0.1);
 const double rlog2  = 1./log(2.);
 const double sqrt2  = sqrt(2.);
 const double rsqrt2 = 1. / sqrt2;
 
 
 void AtoB_next_00(AtoB *unit, int inNumSamples)
 {
 SETUP_AtoB
 LOOP(inNumSamples,
 ASSIGN_ABCD
 
 ZXP(w_out) = (a + b + c + d) * 0.5;
 ZXP(x_out) = (a + b - c - d) * 0.5;
 ZXP(y_out) = (a - b - c + d) * 0.5;
 ZXP(z_out) = (a - b + c - d) * 0.5;
 );
 }
 
 void AtoB_next_10(AtoB *unit, int inNumSamples)
 {
 SETUP_AtoB
 LOOP(inNumSamples,
 ASSIGN_ABCD
 
 ZXP(w_out) = (a + b + c + d) * 0.5;
 ZXP(x_out) = (a + b - c - d) * 0.5;
 ZXP(y_out) = (a - b - c + d) * 0.5;
 ZXP(z_out) = (-a + b - c + d) * 0.5;
 )
 }
 
 void AtoB_next_20(AtoB *unit, int inNumSamples)
 { 
 SETUP_AtoB
 LOOP(inNumSamples,
 ASSIGN_ABCD
 
 ZXP(w_out) = (a + b + c + d) * 0.5;
 ZXP(x_out) = (a + b - c - d) * 0.5;
 ZXP(y_out) = (a - b) * rsqrt2;
 ZXP(z_out) = (c - d) * rsqrt2;
 )
 }
 
 void AtoB_next_30(AtoB *unit, int inNumSamples)
 {
 SETUP_AtoB
 LOOP(inNumSamples,
 ASSIGN_ABCD
 
 ZXP(w_out) = (a + b + c + d) * 0.5;
 ZXP(x_out) = (a + b - c - d) * 0.5;
 ZXP(y_out) = (-c + d) * rsqrt2;
 ZXP(z_out) = (a - b) * rsqrt2;
 )
 }
 
 
 void AtoB_next_40(AtoB *unit, int inNumSamples)
 {
 SETUP_AtoB
 LOOP(inNumSamples, 
 ASSIGN_ABCD
 
 ZXP(w_out) = (a + b + c + d) * 0.5;
 ZXP(x_out) = (a * sqrt3div2) + ((-b - c - d) * sqrt3div6);
 ZXP(y_out) = (-b + c) * rsqrt2;
 ZXP(z_out) = ((b + c) * rsqrt6) - (d * sqrt6div3);
 )
 }
 
 void AtoB_next_50(AtoB *unit, int inNumSamples)
 {
 SETUP_AtoB
 
 LOOP(inNumSamples,
 ASSIGN_ABCD
 
 ZXP(w_out) = (a + b + c + d) * 0.5;
 ZXP(x_out) = ((a + b + c) * sqrt3div6) - (d * sqrt3div2);
 ZXP(y_out) = (a - b) * rsqrt2;
 ZXP(z_out) = ((-a - b) * rsqrt6) - (c * sqrt6div3);
 )
 }
 
 void AtoB_next_60(AtoB *unit, int inNumSamples)
 {
 SETUP_AtoB
 
 LOOP(inNumSamples,
 ASSIGN_ABCD
 
 ZXP(w_out) = (a + b + c + d) * 0.5;
 ZXP(x_out) = (a * sqrt3div2) + ((-b - c - d) * sqrt3div6);
 ZXP(y_out) = (-c + d) * rsqrt2;
 ZXP(z_out) = (b * sqrt6div3) + ((-c - d) * rsqrt6);
 )
 }
 
 void AtoB_next_70(AtoB *unit, int inNumSamples)
 {
 SETUP_AtoB
 
 LOOP(inNumSamples,
 ASSIGN_ABCD
 
 ZXP(w_out) = (a + b + c + d) * 0.5;
 ZXP(x_out) = ((a + b + c) * sqrt3div6) - (d * sqrt3div2);
 ZXP(y_out) = (a - b) * rsqrt2;
 ZXP(z_out) = ((a + b) * rsqrt6) + (c * sqrt6div3);
 )
 }
 ///////////////////////////////////////////////////////////////////////////////////////////////////////
 // Dominante - 4 in 4 out
 
 void Dominate_Ctor(Dominate* unit)
 {
 SETCALC(Dominate_next);
 Dominate_next(unit, 1);
 unit->m_dominance = ZIN0(4) * pi2;
 }
 
 // four channels in, four channels out... 
 // regardless... call in w, x, y, z and out a, b, c, d
 void Dominate_next(Dominate *unit, int inNumSamples)
 {       
 float *Wout = ZOUT(0);
 float *Xout = ZOUT(1);
 float *Yout = ZOUT(2);
 float *Zout = ZOUT(3);
 
 float *Win = ZIN(0);
 float *Xin = ZIN(1);
 float *Yin = ZIN(2);
 float *Zin = ZIN(3);
 float dominanceEnd = ZIN0(4); // dom * 0.5pi
 float dominance = unit->m_dominance;
 
 float cosDom, sinDom;
 float w, x, y, z;
 
 dominanceEnd = dominanceEnd * (float)pi2;
 
 float dominance_slope = CALCSLOPE(dominanceEnd, dominance);
 
 LOOP(inNumSamples,
 
 cosDom = cos(dominance);
 sinDom = sin(dominance);
 
 w = ZXP(Win);
 x = ZXP(Xin);
 y = ZXP(Yin);
 z = ZXP(Zin);
 
 ZXP(Wout) = w + (sqrt2 * sinDom * w);
 ZXP(Xout) = (rsqrt2 * sinDom * x) + x;
 ZXP(Yout) = y * cosDom;
 ZXP(Zout) = z * cosDom;
 
 dominance += dominance_slope;
 );
 
 unit->m_dominance = dominanceEnd;
 }
 
 ///////////////////////////////////////////////////////////////////////////////////////////////////////
 // Push - 4 in 4 out
 
 void Push_Ctor(Push* unit)
 {
 SETCALC(Push_next);
 Push_next(unit, 1);
 unit->m_push = ZIN0(4) * pi2;
 }
 
 // four channels in, four channels out... 
 // regardless... call in w, x, y, z and out a, b, c, d
 void Push_next(Push *unit, int inNumSamples)
 {       
 float *Wout = ZOUT(0);
 float *Xout = ZOUT(1);
 float *Yout = ZOUT(2);
 float *Zout = ZOUT(3);
 
 float *Win = ZIN(0);
 float *Xin = ZIN(1);
 float *Yin = ZIN(2);
 float *Zin = ZIN(3);
 float pushEnd = ZIN0(4); // dom * 0.5pi
 float push = unit->m_push;
 
 float piover4 = pi * 0.25;
 float cosPushFac, sinPushFac, cosPush;
 float w, x, y, z;
 
 pushEnd = pushEnd * (float)pi2;
 
 float push_slope = CALCSLOPE(pushEnd, push);
 
 LOOP(inNumSamples,
 
 cosPushFac = cos(piover4 - push);
 sinPushFac = sin(piover4 - push);
 cosPush = cos(push);
 
 w = ZXP(Win);
 x = ZXP(Xin);
 y = ZXP(Yin);
 z = ZXP(Zin);
 
 ZXP(Wout) = ((1 + (sqrt2 * cosPushFac)) * w) + ((sqrt2 - (2 * sinPushFac)) * w);
 ZXP(Xout) = ((rsqrt2 - cosPushFac) * x) + ((1 + (sqrt2 * sinPushFac)) * x);
 ZXP(Yout) = y * cosPush;
 ZXP(Zout) = z * cosPush;
 
 push += push_slope;
 );
 
 unit->m_push = pushEnd;
 }
 */
////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////

void load(InterfaceTable *inTable)
{
    ft = inTable;
    
    DefineSimpleCantAliasUnit(FOAPanB);
/*
    DefineSimpleCantAliasUnit(FOAPanBWFocus);
    DefineSimpleCantAliasUnit(AtkSterToB);
    DefineDtorCantAliasUnit(AtkPantoF);
 */
    
    DefineSimpleCantAliasUnit(AtkDirectO);
    DefineSimpleCantAliasUnit(AtkDirectX);
    DefineSimpleCantAliasUnit(AtkDirectY);
    DefineSimpleCantAliasUnit(AtkDirectZ);
    DefineSimpleCantAliasUnit(AtkRotate);
    DefineSimpleCantAliasUnit(AtkTilt);
    DefineSimpleCantAliasUnit(AtkTumble);
    DefineSimpleCantAliasUnit(AtkFocusX);	
    DefineSimpleCantAliasUnit(AtkFocusY);
    DefineSimpleCantAliasUnit(AtkFocusZ);
    DefineSimpleCantAliasUnit(AtkPushX);
    DefineSimpleCantAliasUnit(AtkPushY);
    DefineSimpleCantAliasUnit(AtkPushZ);
    DefineSimpleCantAliasUnit(AtkPressX);
    DefineSimpleCantAliasUnit(AtkPressY);
    DefineSimpleCantAliasUnit(AtkPressZ);
    DefineSimpleCantAliasUnit(AtkZoomX);
    DefineSimpleCantAliasUnit(AtkZoomY);
    DefineSimpleCantAliasUnit(AtkZoomZ);
    DefineSimpleCantAliasUnit(AtkDominateX);
    DefineSimpleCantAliasUnit(AtkDominateY);
    DefineSimpleCantAliasUnit(AtkDominateZ);
    
    
    DefineSimpleCantAliasUnit(AtkNFC);
    DefineSimpleCantAliasUnit(AtkProximity);
    DefineSimpleCantAliasUnit(AtkPsychoShelf);
    //DefineSimpleUnit(Dominate);
    //DefineSimpleUnit(Push);
    //DefineSimpleUnit(AtoB);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

