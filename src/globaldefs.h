/*==============================================================================
 HEADER: globaldefs.h		[mgpt]
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

#if !defined(global)
#  define global extern
#endif

#define IPName(param,paramtext)    \
{strcpy(tag[nt],paramtext);    \
addr[nt]=&(param);    \
id[nt++]=INT;}

#define RPName(param,paramtext)    \
{strcpy(tag[nt],paramtext);    \
addr[nt]=&param;    \
id[nt++]=DOUBLE;}

#define BPName(param,paramtext)    \
{strcpy(tag[nt],paramtext);    \
addr[nt]=&param;    \
id[nt++]=BOOLEAN;}

#define SPName(param,paramtext,n)    \
{strcpy(tag[nt],paramtext);    \
param=(string) malloc(n);    \
addr[nt]=param;    \
id[nt++]=STRING;}
//

#include "models.h"


#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#define invH0     2997.92458   //This is H_0^{-1} in Mpc/h units
#define FOURPI2   39.4784176043574    //but see 0903.5321
#define PI2     9.8696044010893586188
#define TWOPI2     19.739208802178716
#define SIXPI2  59.21762640653615
#define INVSQRTDTWOPI 0.39894228040143267794

typedef struct {
// Differential equations evolution parameters:
	real x;
	string dxstr;
	real xstop;
    int maxnsteps;
	string integration_method;
    real dxmin;
    real eps;

// Power spectrum table
    string fnamePS;
    real kmin;
    real kmax;
    int Nk;
//
// CLPT correlation functions table:
    real rmin;
    real rmax;
    int Nr;
// Post processing parameters:
    bool postprocessing;
    string options;
//
    string paramfile;

// Modified gravity model parameters:
    string mgmodel;
    string suffixModel;
    string model_paramfile;
    int nHS;
    real fR0;
    real omegaBD;
    real screening;
// DGP:
    real eps_DGP;
    real rc_DGP;
//
// Background cosmology:
    real om;
    string olstr;
    real h;
//
// Quadrature parameters:
    string quadratureMethod;
    int nquadSteps;
    int ngausslegpoints;
    real epsquad;
//
} cmdline_data, *cmdline_data_ptr;

typedef struct {
	real cpuinit;
	real dx;
    int method_int;
    int quadmethod_int;

// Modified gravity model parameters:
    real beta2;
//
    char integration_method_comment[100];
    char quadraturemethod_comment[100];

	string headline0;
	string headline1;
	string headline2;
	string headline3;

    char model_comment[100];

	FILE *outlog;
    
    real ol;

    real xnow;
    real xout;
    real xoutinfo;
    real xstop;

	char mode[2];

// 
    char fnamePSPath[100];
    char logfilePath[100];
    char clptDir[100];
    char tmpDir[100];
    char fpfnamekfun[100];
    char fpfnameSPTPowerSpectrum[100];
    char fpfnameqfunctions[100];
    char fpfnameclptfunctions[100];

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
#define NEQS3Order    10
#define NEQS2Order    8
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

/*
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
*/

/*
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
*/

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

//#define etaD2v2(x)    (((global_D2v2_ptr) (x))->eta)
//#define Dpk1D2v2(x)    (((global_D2v2_ptr) (x))->y1)
//#define Dpk2D2v2(x)    (((global_D2v2_ptr) (x))->y3)
#define etaD2(x)    (((global_D2v2_ptr) (x))->eta)
#define Dpk1D2(x)    (((global_D2v2_ptr) (x))->y1)
#define Dpk2D2(x)    (((global_D2v2_ptr) (x))->y3)
#define DA2D2(x)    (((global_D2v2_ptr) (x))->y5)
#define DB2D2(x)    (((global_D2v2_ptr) (x))->y7)
//

/*
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
*/

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

//#define etaD3v2(x)    (((global_D3v2_ptr) (x))->eta)
//#define DpkD3v2(x)    (((global_D3v2_ptr) (x))->y1)
//#define DppD3v2(x)    (((global_D3v2_ptr) (x))->y3)
//#define D2fD3v2(x)    (((global_D3v2_ptr) (x))->y5)
//#define D2mfD3v2(x)    (((global_D3v2_ptr) (x))->y7)
//#define D3symmD3v2(x)    (((global_D3v2_ptr) (x))->y9)
//
#define etaD3(x)    (((global_D3v2_ptr) (x))->eta)
#define DpkD3(x)    (((global_D3v2_ptr) (x))->y1)
#define DppD3(x)    (((global_D3v2_ptr) (x))->y3)
#define D2fD3(x)    (((global_D3v2_ptr) (x))->y5)
#define D2mfD3(x)    (((global_D3v2_ptr) (x))->y7)
#define D3symmD3(x)    (((global_D3v2_ptr) (x))->y9)
//

// GL structure
typedef struct {
    int npts;
    real x1;
    real x2;
    real *xgl;
    real *wgl;
} global_GL, *global_GL_ptr;

// STATIC problem: gcc version 11
//global_GL_ptr pGL;
global global_GL_ptr pGL;

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

//
// qfunctions structure
typedef struct {
    real q;
    real U10L;
    real U10loop;
    real U11;
    real U20;
    real XL;
    real Xloop;
    real X10;
    real YL;
    real Yloop;
    real Y10;
    real VT;
    real TT;
} global_qfunctions, *global_qfunctions_ptr;

#define qqfun(x)    (((global_qfunctions_ptr) (x))->q)
#define U10Lqfun(x)    (((global_qfunctions_ptr) (x))->U10L)
#define U10loopqfun(x)    (((global_qfunctions_ptr) (x))->U10loop)
#define U11qfun(x)    (((global_qfunctions_ptr) (x))->U11)
#define U20qfun(x)    (((global_qfunctions_ptr) (x))->U20)
#define XLqfun(x)    (((global_qfunctions_ptr) (x))->XL)
#define Xloopqfun(x)    (((global_qfunctions_ptr) (x))->Xloop)
#define X10qfun(x)    (((global_qfunctions_ptr) (x))->X10)
#define YLqfun(x)    (((global_qfunctions_ptr) (x))->YL)
#define Yloopqfun(x)    (((global_qfunctions_ptr) (x))->Yloop)
#define Y10qfun(x)    (((global_qfunctions_ptr) (x))->Y10)
#define VTqfun(x)    (((global_qfunctions_ptr) (x))->VT)
#define TTqfun(x)    (((global_qfunctions_ptr) (x))->TT)
//

//
// correlation functions structure
typedef struct {
    real q;
    real xi;
    real Lapxi;
    real nabla4xi;
} global_corrfunctions, *global_corrfunctions_ptr;

#define qcorrfun(x)    (((global_corrfunctions_ptr) (x))->q)
#define xicorrfun(x)    (((global_corrfunctions_ptr) (x))->xi)
#define Lapxicorrfun(x)    (((global_corrfunctions_ptr) (x))->Lapxi)
#define nabla4xicorrfun(x)    (((global_corrfunctions_ptr) (x))->nabla4xi)
//

// BEGIN :: CLPT correlation auxiliary functions and structures
typedef struct {
    real r;
    real xi;
} global_zacorrfunctions, *global_zacorrfunctions_ptr;

#define rzacorrfun(x)    (((global_zacorrfunctions_ptr) (x))->r)
#define xizacorrfun(x)    (((global_zacorrfunctions_ptr) (x))->xi)

typedef struct {
    real r;
    real xiA;
    real xiW;
    real xi10L;
    real xi10loop;
    real xi20L;
    real xi20loop;
    real xi01;
    real xi02;
    real xi11;
    real Lapxi;
    real nabla4xi;
} global_clptcorrfunctions, *global_clptcorrfunctions_ptr;

#define rclptcorrfun(x)    (((global_clptcorrfunctions_ptr) (x))->r)
#define xiAclptcorrfun(x)    (((global_clptcorrfunctions_ptr) (x))->xiA)
#define xiWclptcorrfun(x)    (((global_clptcorrfunctions_ptr) (x))->xiW)
#define xi10Lclptcorrfun(x)    (((global_clptcorrfunctions_ptr) (x))->xi10L)
#define xi10loopclptcorrfun(x)    (((global_clptcorrfunctions_ptr) (x))->xi10loop)
#define xi20Lclptcorrfun(x)    (((global_clptcorrfunctions_ptr) (x))->xi20L)
#define xi20loopclptcorrfun(x)    (((global_clptcorrfunctions_ptr) (x))->xi20loop)
#define xi01clptcorrfun(x)    (((global_clptcorrfunctions_ptr) (x))->xi01)
#define xi02clptcorrfun(x)    (((global_clptcorrfunctions_ptr) (x))->xi02)
#define xi11clptcorrfun(x)    (((global_clptcorrfunctions_ptr) (x))->xi11)
#define Lapxiclptcorrfun(x)    (((global_clptcorrfunctions_ptr) (x))->Lapxi)
#define nabla4xiclptcorrfun(x)    (((global_clptcorrfunctions_ptr) (x))->nabla4xi)
// END :: CLPT correlation auxiliary functions and structures

// STATIC problem: gcc version 11
// From inout.h
global real *inout_xval;
global real *inout_yval;
global real *inout_zval;
global real *inout_wval;

// STATIC problem: gcc version 11
// From diffeqs.h
global double dxsav,*xp,**yp;
global int kmax,kount;
global int nrhs;

// STATIC problem: gcc version 11
// From stdinc.h
global long idum;                // seed for random generators

#endif // ! _globaldefs_h

