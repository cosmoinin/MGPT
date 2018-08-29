/*==============================================================================
 HEADER: globaldefs.h		[mglpt]
 Written by: Mario A. Rodriguez-Meza
 Starting date: January 2018
 Purpose: Definitions of global variables and parameters
 Language: C
 Use: '#include "global_defs.h"
 Use in routines and functions:
 External headers: stdinc.h, data_struc_defs.h
 Comments and notes:
 Info: Mario A. Rodriguez-Meza
 Depto. de Fisica, ININ
 Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
 e-mail: marioalberto.rodriguez@inin.gob.mx
 http://www.inin.gob.mx/
 
 Major revisions:
 Copyright: (c) 2005-2018 Mar.  All Rights Reserved
 ================================================================================
 Legal matters:
 The author does not warrant that the program and routines it contains
 listed below are free from error or suitable for particular applications,
 and he disclaims all liability from any consequences arising from their	use.
 ==============================================================================*/

#ifndef _globaldefs_h
#define _globaldefs_h

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
//

#include "general_libs/stdinc.h"
#include "general_libs/numrec.h"
#include "general_libs/diffeqs.h"
#include "general_libs/quads.h"
#include "general_libs/mathfns.h"
#include "general_libs/mathutil.h"
#include "general_libs/inout.h"
#include "general_libs/vectmath.h"
#include "general_libs/getparam.h"
#include "general_libs/machines.h"
#include "general_libs/strings.h"


#include "data_struc_defs.h"
#include "models.h"


#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#define H02     2997.92458
#define FOURPI2   39.4784176043574

typedef struct {
	real x;
	string dxstr;
	real xstop;
    int maxnsteps;
	string integration_method;
    int ngausslegpoints;

// Power spectrum table
    string fnamePS;
    real kmin;
    real kmax;
    int Nk;
//
// Post processing parameters:
    bool postprocessing;
    string options;
//
    string paramfile;

// Modified gravity model parameters:
    string mgmodel;
//    string model_paramfile;
    int nHS;
    real fR0;
//    string beta2str;
    real omegaBD;
    real screening;
//
    real om;
    real h;
//
    real dxmin;
    real eps;
//
} cmdline_data, *cmdline_data_ptr;

typedef struct {
	real cpuinit;
	real dx;
    int method_int;

// Modified gravity model parameters:
    real beta2;
//
    char integration_method_comment[100];

	string headline0;
	string headline1;
	string headline2;
	string headline3;

    char model_comment[100];

	FILE *outlog;

	int stopflag;
    
    real xnow;
    real xout;
    real xoutinfo;

	char mode[2];

// 
    char fnamePSPath[100];
    char logfilePath[100];

    real kf;
    real k1;
    real k2;

    real x;
    real k;
    real p;
} global_data, *global_data_ptr;


global global_data gd;
global cmdline_data cmd;

global real *yout;
#define NEQS            2
#define NEQS3Order      20
#define NEQS3Orderv2    10
#define NEQS2Order      12
#define NEQS2Orderv2    8
#define NEQS1Order      2

typedef struct _pointPSTable {
    real k;
    real ps;
} pointPSTable, *pointPSTableptr;

global int nPSTable;
global pointPSTableptr PSLCDMtab;
global int nPSLogT;
global pointPSTableptr PSLCDMLogtab;

global int nPSLT;
global pointPSTableptr PSLT;
global int nPSLTLog;
global pointPSTableptr PSLTLog;

global real *kPS;
global real *pPS;
global real *pPS2;


#define kPos(x)    (((pointPSTableptr) (x))->k)
#define PS(x)    (((pointPSTableptr) (x))->ps)

//
typedef struct {
    real eta;
    real y1;
    real y2;
} global_D1, *global_D1_ptr;

#define etaD1(x)    (((global_D1_ptr) (x))->eta)
#define DpD1(x)    (((global_D1_ptr) (x))->y1)
#define DppD1(x)    (((global_D1_ptr) (x))->y2)
//

//
typedef struct {
    real eta;
    real y1;
    real y2;
    real y3;
    real y4;
    real y5;
    real y6;
    real y7;
    real y8;
    real y9;
    real y10;
    real y11;
    real y12;
} global_D2, *global_D2_ptr;

#define etaD2(x)    (((global_D2_ptr) (x))->eta)
#define Dpk1D2(x)    (((global_D2_ptr) (x))->y1)
#define Dpk2D2(x)    (((global_D2_ptr) (x))->y3)
#define Dak1k2D2(x)    (((global_D2_ptr) (x))->y5)
#define Dbk1k2D2(x)    (((global_D2_ptr) (x))->y7)
#define DFLk1k2D2(x)    (((global_D2_ptr) (x))->y9)
#define DdIk1k2D2(x)    (((global_D2_ptr) (x))->y11)
//

//
typedef struct {
    real eta;
    real y1;
    real y2;
    real y3;
    real y4;
    real y5;
    real y6;
    real y7;
    real y8;
} global_D2v2, *global_D2v2_ptr;

#define etaD2v2(x)    (((global_D2v2_ptr) (x))->eta)
#define Dpk1D2v2(x)    (((global_D2v2_ptr) (x))->y1)
#define Dpk2D2v2(x)    (((global_D2v2_ptr) (x))->y3)
#define DA2D2(x)    (((global_D2v2_ptr) (x))->y5)
#define DB2D2(x)    (((global_D2v2_ptr) (x))->y7)
//

//
typedef struct {
    real eta;
    real y1;
    real y2;
    real y3;
    real y4;
    real y5;
    real y6;
    real y7;
    real y8;
    real y9;
    real y10;
    real y11;
    real y12;
    real y13;
    real y14;
    real y15;
    real y16;
    real y17;
    real y18;
    real y19;
    real y20;
} global_D3, *global_D3_ptr;

#define etaD3(x)    (((global_D3_ptr) (x))->eta)
#define DpkD3(x)    (((global_D3_ptr) (x))->y1)
#define DppD3(x)    (((global_D3_ptr) (x))->y3)
#define Da2D3(x)    (((global_D3_ptr) (x))->y5)
#define Db2D3(x)    (((global_D3_ptr) (x))->y7)
#define DFL2D3(x)    (((global_D3_ptr) (x))->y9)
#define DdI2D3(x)    (((global_D3_ptr) (x))->y11)
#define DIA3D3(x)    (((global_D3_ptr) (x))->y13)
#define DIB3D3(x)    (((global_D3_ptr) (x))->y15)
#define DIFL3D3(x)    (((global_D3_ptr) (x))->y17)
#define DIdI3D3(x)    (((global_D3_ptr) (x))->y19)
//

//
typedef struct {
    real eta;
    real y1;
    real y2;
    real y3;
    real y4;
    real y5;
    real y6;
    real y7;
    real y8;
    real y9;
    real y10;
} global_D3v2, *global_D3v2_ptr;

#define etaD3v2(x)    (((global_D3v2_ptr) (x))->eta)
#define DpkD3v2(x)    (((global_D3v2_ptr) (x))->y1)
#define DppD3v2(x)    (((global_D3v2_ptr) (x))->y3)
#define D2fD3v2(x)    (((global_D3v2_ptr) (x))->y5)
#define D2mfD3v2(x)    (((global_D3v2_ptr) (x))->y7)
#define D3symmD3v2(x)    (((global_D3v2_ptr) (x))->y9)
//

// GL structure
typedef struct {
    int npts;
    real x1;
    real x2;
    real *xgl;
    real *wgl;
} global_GL, *global_GL_ptr;

global_GL_ptr pGL;

#define nGL(x)    (((global_GL_ptr) (x))->npts)
#define x1GL(x)    (((global_GL_ptr) (x))->x1)
#define x2GL(x)    (((global_GL_ptr) (x))->x2)
#define xGL(x)    (((global_GL_ptr) (x))->xgl)
#define wGL(x)    (((global_GL_ptr) (x))->wgl)

//
// QRs structure
typedef struct {
    int eta;
    real k;
    real Q1;
    real Q2;
    real Q3;
    real Q8;
    real Q9;
    real Q13;
    real QI;
    real Q5;
    real Q7;
    real Q11;
    real Q12;
    real RI;
    real R1p2;
    real R1;
    real R2;
} global_QRs, *global_QRs_ptr;

#define etaQRs(x)    (((global_QRs_ptr) (x))->eta)
#define kQRs(x)    (((global_QRs_ptr) (x))->k)
#define Q1(x)    (((global_QRs_ptr) (x))->Q1)
#define Q2(x)    (((global_QRs_ptr) (x))->Q2)
#define Q3(x)    (((global_QRs_ptr) (x))->Q3)
#define Q8(x)    (((global_QRs_ptr) (x))->Q8)
#define Q9(x)    (((global_QRs_ptr) (x))->Q9)
#define Q13(x)    (((global_QRs_ptr) (x))->Q13)
#define QI(x)    (((global_QRs_ptr) (x))->QI)
#define Q5(x)    (((global_QRs_ptr) (x))->Q5)
#define Q7(x)    (((global_QRs_ptr) (x))->Q7)
#define Q11(x)    (((global_QRs_ptr) (x))->Q11)
#define Q12(x)    (((global_QRs_ptr) (x))->Q12)
#define RI(x)    (((global_QRs_ptr) (x))->RI)
#define R1p2(x)    (((global_QRs_ptr) (x))->R1p2)
#define R1(x)    (((global_QRs_ptr) (x))->R1)
#define R2(x)    (((global_QRs_ptr) (x))->R2)
//

#endif // ! _globaldefs_h

