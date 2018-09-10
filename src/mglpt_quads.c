/*==============================================================================
 MODULE: mglpt_quads.c				[mglpt]
 ==============================================================================*/

#include "globaldefs.h"
#include "protodefs.h"
#include "models.h"


//global_QRs QsRs_functions(real eta, real ki);
//global_QRs QsRs_functions_LCDM(real eta, real ki);
global_QRs QsRs_functions_trapezoid1(real eta, real ki);
global_QRs QsRs_functions_trapezoid2(real eta, real ki);
global_QRs QsRs_functions_trapezoidAA(real eta, real ki);
global_QRs QsRs_functions_romo(real eta, real ki);

local real KQ8_function(real k, real r, real x);
local real KQ9_function(real k, real r, real x);
local real KQ13_function(real k, real r, real x);
local real KQI_function(real k, real r, real x);

local real KQ5_function(real k, real r, real x);
local real KQ7_function(real k, real r, real x);
local real KQ11_function(real k, real r, real x);
local real KQ12_function(real k, real r, real x);

local real KRI_function(real k, real r, real x);
local real KR1p2_function(real k, real r, real x);

local real GaussLegendreQ1_func(real y);
local real GaussLegendreQ2_func(real y);
local real GaussLegendreQ3_func(real y);
local real GaussLegendreQ8_func(real y);
local real GaussLegendreQ9_func(real y);
local real GaussLegendreQ13_func(real y);
local real GaussLegendreQI_func(real y);
local real GaussLegendreQ5_func(real y);
local real GaussLegendreQ7_func(real y);
local real GaussLegendreQ11_func(real y);
local real GaussLegendreQ12_func(real y);
local real GaussLegendreRI_func(real y);
local real GaussLegendreR1p2_func(real y);
local real GaussLegendreR1_func(real y);
local real GaussLegendreR2_func(real y);
local real funcR1int(real y);

local real Q1_function(real eta, real ki);
local real Q2_function(real eta, real ki);
local real Q3_function(real eta, real ki);
local real Q8_function(real eta, real ki);
local real Q9_function(real eta, real ki);
local real Q13_function(real eta, real ki);
local real QI_function(real eta, real ki);
local real Q5_function(real eta, real ki);
local real Q7_function(real eta, real ki);
local real Q11_function(real eta, real ki);
local real Q12_function(real eta, real ki);
local real RI_function(real eta, real ki);
local real R1p2_function(real eta, real ki);
local real R1_function(real eta, real ki);
local real R2_function(real eta, real ki);

local void quadrature(real ki);

#define abskmq      (1.0+rsqr(rr)-2.0*rr*xv)

// BEGIN Q1

local real GaussLegendreQ1_func(real y)
{
    global_D2v2_ptr ptmp;
    int j;
    real *xGL, *wGL;
    real kmin, kmax, ki;
    real Q1p, Q1aA, Q1aB, KQ1;
    real PSLA, PSLB;
    real rmin, rmax;
    real kk, rr, deltar;
    real mumin, mumax;
    real xv, w, k2, psl;
    real KA, KB;
    int Nx;
    
    gd.p = rpow(10.0,y);
    
    kk = gd.p;
    ki = gd.k;
    rr = kk/ki;
    PSLB = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
    Q1aB = 0.0;
    kmin = kPos(PSLT+1);
    kmax = kPos(PSLT+nPSLT-1);
    rmax = kmax/ki;
    rmin = kmin/ki;
    
    mumin = MAX(-1.0, (1.0 + rsqr(rr) - rsqr(rmax))/(2.0*rr));
    mumax = MIN( 1.0, (1.0 + rsqr(rr) - rsqr(rmin))/(2.0*rr));
    if (rr>=0.5)
        mumax = 0.5/rr;
    Nx=10;
    xGL=dvector(1,Nx);
    wGL=dvector(1,Nx);
    gauleg(mumin,mumax,xGL,wGL,Nx);
    for (j=1; j<=Nx; j++) {
        xv = xGL[j];
        w = wGL[j];
        k2 = ki * rsqrt(1.0 + rsqr(rr) - 2.0*rr*xv);
        ptmp = DsSecondOrder_func(ki, ki*rr, k2);
        if (model_int_flag==LCDM) {
            KQ1 = rsqr(
                       KA_LCDM + KB_LCDM*(-1.0+(1.0-rsqr(xv))/abskmq)
                       );
        } else {
            KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
            KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
            KQ1 = rsqr(
                   KA + KB*(-1.0+(1.0-rsqr(xv))/abskmq)
                   );
        }
        psl = psInterpolation_nr(ki*rsqrt(abskmq), kPS, pPS, nPSLT);
        Q1aB += w*KQ1*psl;
    }

    free_dvector(wGL,1,Nx);
    free_dvector(xGL,1,Nx);

    return 2.0*rpow(gd.p,3.0)*PSLB*Q1aB;
}

local real Q1_function(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;

    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    result=(rlog(10.0)/FOURPI2)
            *qromo(GaussLegendreQ1_func,ymin,ymax,midpnt,cmd.epsquad);

    return result;
}

// END DE Q1


// BEGIN Q2

local real GaussLegendreQ2_func(real y)
{
    global_D2v2_ptr ptmp;
    int j;
    real *xGL, *wGL;
    real kmin, kmax, ki;
    real Q2p, Q2aA, Q2aB, KQ2;
    real PSLA, PSLB;
    real rmin, rmax;
    real kk, rr, deltar;
    real mumin, mumax;
    real xv, w, k2, psl;
    real KA, KB;
    int Nx;
    
    gd.p = rpow(10.0,y);
    
    kk = gd.p;
    ki = gd.k;
    rr = kk/ki;
    PSLB = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
    Q2aB = 0.0;
    kmin = kPos(PSLT+1);
    kmax = kPos(PSLT+nPSLT-1);
    rmax = kmax/ki;
    rmin = kmin/ki;
    
    mumin = MAX(-1.0, (1.0 + rsqr(rr) - rsqr(rmax))/(2.0*rr));
    mumax = MIN( 1.0, (1.0 + rsqr(rr) - rsqr(rmin))/(2.0*rr));
    if (rr>=0.5)
        mumax = 0.5/rr;
    Nx=10;
    xGL=dvector(1,Nx);
    wGL=dvector(1,Nx);
    gauleg(mumin,mumax,xGL,wGL,Nx);
    for (j=1; j<=Nx; j++) {
        xv = xGL[j];
        w = wGL[j];
        k2 = ki * rsqrt(1.0 + rsqr(rr) - 2.0*rr*xv);
        ptmp = DsSecondOrder_func(ki, ki*rr, k2);
        if (model_int_flag==LCDM) {
            KQ2 = (rr*xv*(1.0-rr*xv)/abskmq)
            *(
              KA_LCDM - KB_LCDM*((rsqr(xv)+rsqr(rr)-2.0*rr*xv)/abskmq)
              );
        } else {
            KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
            KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
            KQ2 = (rr*xv*(1.0-rr*xv)/abskmq)
                *(
                  KA - KB*((rsqr(xv)+rsqr(rr)-2.0*rr*xv)/abskmq)
                  );
        }
        psl = psInterpolation_nr(ki*rsqrt(abskmq), kPS, pPS, nPSLT);
        Q2aB += w*KQ2*psl;
    }

    free_dvector(wGL,1,Nx);
    free_dvector(xGL,1,Nx);
    
    return 2.0*rpow(ki,2.0)*gd.p*PSLB*Q2aB;
}

local real Q2_function(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    result=(rlog(10.0)/FOURPI2)
        *qromo(GaussLegendreQ2_func,ymin,ymax,midpnt,cmd.epsquad);

    return result;
}

// END DE Q2


// BEGIN Q3

local real GaussLegendreQ3_func(real y)
{
    int j;
    real *xGL, *wGL;
    real kmin, kmax, ki;
    real Q3p, Q3aA, Q3aB, KQ3;
    real PSLA, PSLB;
    real rmin, rmax;
    real kk, rr, deltar;
    real mumin, mumax;
    real xv, w, k2, psl;
    int Nx;
    
    gd.p = rpow(10.0,y);
    
    kk = gd.p;
    ki = gd.k;
    rr = kk/ki;
    PSLB = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
    Q3aB = 0.0;
    kmin = kPos(PSLT+1);
    kmax = kPos(PSLT+nPSLT-1);
    rmax = kmax/ki;
    rmin = kmin/ki;
    
    mumin = MAX(-1.0, (1.0 + rsqr(rr) - rsqr(rmax))/(2.0*rr));
    mumax = MIN( 1.0, (1.0 + rsqr(rr) - rsqr(rmin))/(2.0*rr));
    if (rr>=0.5)
        mumax = 0.5/rr;
    Nx=10;
    xGL=dvector(1,Nx);
    wGL=dvector(1,Nx);
    gauleg(mumin,mumax,xGL,wGL,Nx);
    for (j=1; j<=Nx; j++) {
        xv = xGL[j];
        w = wGL[j];
        k2 = ki * rsqrt(1.0 + rsqr(rr) - 2.0*rr*xv);
        KQ3 = rsqr(xv)*rsqr(1.0-rr*xv)/rsqr(abskmq);
        psl = psInterpolation_nr(ki*rsqrt(abskmq), kPS, pPS, nPSLT);
        Q3aB += w*KQ3*psl;
    }
    
    free_dvector(wGL,1,Nx);
    free_dvector(xGL,1,Nx);
    
    return 2.0*rpow(ki,2.0)*gd.p*PSLB*Q3aB;
}

local real Q3_function(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    result=(rlog(10.0)/FOURPI2)
    *qromo(GaussLegendreQ3_func,ymin,ymax,midpnt,cmd.epsquad);
    
    return result;
}

// END Q3


// BEGIN Q8

local real KQ8_function(real ki, real rr, real xv)
{
    real k2, KA, KB, KQ8;
    global_D2v2_ptr ptmp;

    k2 = ki * rsqrt(abskmq);
    ptmp = DsSecondOrder_func(ki, ki*rr, k2);
    if (model_int_flag==LCDM) {
        KQ8 = rsqr(rr)*(
                        KA_LCDM - KB_LCDM*rsqr(-rr+xv)/abskmq
                        );
    } else {
        KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KQ8 = rsqr(rr)*(
                    KA - KB*rsqr(-rr+xv)/abskmq
                    );
    }
    return KQ8;
}

local real GaussLegendreQ8_func(real y)
{
    int j;
    real *xGL, *wGL;
    real kmin, kmax, ki;
    real Q8, KQ8;
    real PSLA, PSLB;
    real rmin, rmax;
    real kk, rr;
    real xmin, xmax;
    real xv, w, psl;
    int Nx;
    
    gd.p = rpow(10.0,y);
    
    kk = gd.p;
    ki = gd.k;
    rr = kk/ki;
    PSLB = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
    Q8 = 0.0;
    kmin = kPos(PSLT+1);
    kmax = kPos(PSLT+nPSLT-1);
    rmax = kmax/ki;
    rmin = kmin/ki;
    
    xmin = MAX(-1.0, (1.0 + rsqr(rr) - rsqr(rmax))/(2.0*rr));
    xmax = MIN( 1.0, (1.0 + rsqr(rr) - rsqr(rmin))/(2.0*rr));
    if (rr>=0.5)
        xmax = 0.5/rr;
    Nx=10;
    xGL=dvector(1,Nx);
    wGL=dvector(1,Nx);
    gauleg(xmin,xmax,xGL,wGL,Nx);
    for (j=1; j<=Nx; j++) {
        xv = xGL[j];
        w = wGL[j];
        KQ8 = KQ8_function(ki, rr, xv);
        psl = psInterpolation_nr(ki*rsqrt(abskmq), kPS, pPS, nPSLT);
        Q8 += w*KQ8*psl;
    }
    
    free_dvector(wGL,1,Nx);
    free_dvector(xGL,1,Nx);
    
    return gd.p*PSLB*Q8;
}

local real Q8_function(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    result=2.0*rlog(10.0)*(rsqr(ki)/FOURPI2)
    *qromo(GaussLegendreQ8_func,ymin,ymax,midpnt,cmd.epsquad);
    
    return result;
}

// END DE Q8


// BEGIN Q9

local real KQ9_function(real ki, real rr, real xv)
{
    real KQ9;
    
    KQ9 = rr*xv*(1-rr*xv)/abskmq;
    return KQ9;
}

local real GaussLegendreQ9_func(real y)
{
    int j;
    real *xGL, *wGL;
    real kmin, kmax, ki;
    real Q9, KQ9;
    real PSLA, PSLB;
    real rmin, rmax;
    real kk, rr;
    real xmin, xmax;
    real xv, w, psl;
    int Nx;

    gd.p = rpow(10.0,y);
    
    kk = gd.p;
    ki = gd.k;
    rr = kk/ki;
    PSLB = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
    Q9 = 0.0;
    kmin = kPos(PSLT+1);
    kmax = kPos(PSLT+nPSLT-1);
    rmax = kmax/ki;
    rmin = kmin/ki;
    
    xmin = MAX(-1.0, (1.0 + rsqr(rr) - rsqr(rmax))/(2.0*rr));
    xmax = MIN( 1.0, (1.0 + rsqr(rr) - rsqr(rmin))/(2.0*rr));
    if (rr>=0.5)
        xmax = 0.5/rr;
    Nx=10;
    xGL=dvector(1,Nx);
    wGL=dvector(1,Nx);
    gauleg(xmin,xmax,xGL,wGL,Nx);
    for (j=1; j<=Nx; j++) {
        xv = xGL[j];
        w = wGL[j];
        KQ9 = KQ9_function(ki, rr, xv);
        psl = psInterpolation_nr(ki*rsqrt(abskmq), kPS, pPS, nPSLT);
        Q9 += w*KQ9*psl;
    }
    
    free_dvector(wGL,1,Nx);
    free_dvector(xGL,1,Nx);
    
    return gd.p*PSLB*Q9;
}

local real Q9_function(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    result=2.0*rlog(10.0)*(rsqr(ki)/FOURPI2)
    *qromo(GaussLegendreQ9_func,ymin,ymax,midpnt,cmd.epsquad);
    
    return result;
}

// END DE Q9

// BEGIN Q13

local real KQ13_function(real ki, real rr, real xv)
{
    real KQ13;
    
    KQ13 = rsqr(rr);
    return KQ13;
}

local real GaussLegendreQ13_func(real y)
{
    int j;
    real *xGL, *wGL;
    real kmin, kmax, ki;
    real Q13, KQ13;
    real PSLA, PSLB;
    real rmin, rmax;
    real kk, rr;
    real xmin, xmax;
    real xv, w, psl;
    int Nx;
    
    gd.p = rpow(10.0,y);
    
    kk = gd.p;
    ki = gd.k;
    rr = kk/ki;
    PSLB = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
    Q13 = 0.0;
    kmin = kPos(PSLT+1);
    kmax = kPos(PSLT+nPSLT-1);
    rmax = kmax/ki;
    rmin = kmin/ki;
    
    xmin = MAX(-1.0, (1.0 + rsqr(rr) - rsqr(rmax))/(2.0*rr));
    xmax = MIN( 1.0, (1.0 + rsqr(rr) - rsqr(rmin))/(2.0*rr));
    if (rr>=0.5)
        xmax = 0.5/rr;
    Nx=10;
    xGL=dvector(1,Nx);
    wGL=dvector(1,Nx);
    gauleg(xmin,xmax,xGL,wGL,Nx);
    for (j=1; j<=Nx; j++) {
        xv = xGL[j];
        w = wGL[j];
        KQ13 = KQ13_function(ki, rr, xv);
        psl = psInterpolation_nr(ki*rsqrt(abskmq), kPS, pPS, nPSLT);
        Q13 += w*KQ13*psl;
    }
    
    free_dvector(wGL,1,Nx);
    free_dvector(xGL,1,Nx);
    
    return gd.p*PSLB*Q13;
}

local real Q13_function(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    result=2.0*rlog(10.0)*(rsqr(ki)/FOURPI2)
    *qromo(GaussLegendreQ13_func,ymin,ymax,midpnt,cmd.epsquad);
    
    return result;
}

// END DE Q13

// BEGIN QI

local real KQI_function(real ki, real rr, real xv)
{
    real k2, KA, KB, KQI;
    global_D2v2_ptr ptmp;
    
    k2 = ki * rsqrt(abskmq);
    ptmp = DsSecondOrder_func(ki, ki*rr, k2);
    if (model_int_flag==LCDM) {
        KQI = rsqr(rr) * (1.0 - rsqr(xv))/(1.0 + rsqr(rr - 2.0*rr*xv))
            * (
               KA_LCDM - KB_LCDM*(rsqr(xv) + rsqr(rr) - 2.0*rr*xv)/
                (1.0 + rsqr(rr) - 2.0*rr*xv)
               );
    } else {
        KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KQI = rsqr(rr) * (1.0 - rsqr(xv))/(1.0 + rsqr(rr - 2.0*rr*xv))
            * (
                KA - KB*(rsqr(xv) + rsqr(rr) - 2.0*rr*xv)/
                (1.0 + rsqr(rr) - 2.0*rr*xv)
            );
    }
    return KQI;
}

local real GaussLegendreQI_func(real y)
{
    int j;
    real *xGL, *wGL;
    real kmin, kmax, ki;
    real QI, KQI;
    real PSLA, PSLB;
    real rmin, rmax;
    real kk, rr;
    real xmin, xmax;
    real xv, w, psl;
    int Nx;
    
    gd.p = rpow(10.0,y);
    
    kk = gd.p;
    ki = gd.k;
    rr = kk/ki;
    PSLB = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
    QI = 0.0;
    kmin = kPos(PSLT+1);
    kmax = kPos(PSLT+nPSLT-1);
    rmax = kmax/ki;
    rmin = kmin/ki;
    
    xmin = MAX(-1.0, (1.0 + rsqr(rr) - rsqr(rmax))/(2.0*rr));
    xmax = MIN( 1.0, (1.0 + rsqr(rr) - rsqr(rmin))/(2.0*rr));
    if (rr>=0.5)
        xmax = 0.5/rr;
    Nx=10;
    xGL=dvector(1,Nx);
    wGL=dvector(1,Nx);
    gauleg(xmin,xmax,xGL,wGL,Nx);
    for (j=1; j<=Nx; j++) {
        xv = xGL[j];
        w = wGL[j];
        KQI = KQI_function(ki, rr, xv);
        psl = psInterpolation_nr(ki*rsqrt(abskmq), kPS, pPS, nPSLT);
        QI += w*KQI*psl;
    }
    
    free_dvector(wGL,1,Nx);
    free_dvector(xGL,1,Nx);
    
    return gd.p*PSLB*QI;
}

local real QI_function(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    result=2.0*rlog(10.0)*(rsqr(ki)/FOURPI2)
    *qromo(GaussLegendreQI_func,ymin,ymax,midpnt,cmd.epsquad);
    
    return result;
}

// END DE QI


// BEGIN Q5

local real KQ5_function(real ki, real rr, real xv)
{
    real k2, KA, KB, KQ5;
    global_D2v2_ptr ptmp;
    
    k2 = ki * rsqrt(abskmq);
    ptmp = DsSecondOrder_func(ki, ki*rr, k2);
    if (model_int_flag==LCDM) {
        KQ5 = rr*xv*(
                     KA_LCDM - KB_LCDM*rsqr(-rr+xv)/abskmq
                     );
    } else {
        KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KQ5 = rr*xv*(
                    KA - KB*rsqr(-rr+xv)/abskmq
                    );
    }
    return KQ5;
}

local real GaussLegendreQ5_func(real y)
{
    int j;
    real ki;
    real Q5, KQ5;
    real PSLB;
    real kk, rr;
    real xv, w, psl;
    
    gd.p = rpow(10.0,y);
    
    kk = gd.p;
    ki = gd.k;
    rr = kk/ki;
    PSLB = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
    Q5 = 0.0;
    
    for (j=1; j<=nGL(pGL); j++) {
        xv = xGL(pGL)[j];
        w = wGL(pGL)[j];
        KQ5 = KQ5_function(ki, rr, xv);
        psl = psInterpolation_nr(ki*rsqrt(abskmq), kPS, pPS, nPSLT);
        Q5 += w*KQ5*psl;
    }

    return gd.p*PSLB*Q5;
}

local real Q5_function(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);

    result=rlog(10.0)*(rsqr(ki)/FOURPI2)
    *qromo(GaussLegendreQ5_func,ymin,ymax,midpnt,cmd.epsquad);
    
    return result;
}

// END DE Q5

// BEGIN Q7

local real KQ7_function(real ki, real rr, real xv)
{
    real KQ7;
    
    KQ7 = rsqr(xv)*(1-rr*xv)/abskmq;
    return KQ7;
}

local real GaussLegendreQ7_func(real y)
{
    int j;
    real ki;
    real Q7, KQ7;
    real PSLB;
    real kk, rr;
    real xv, w, psl;
    
    gd.p = rpow(10.0,y);
    
    kk = gd.p;
    ki = gd.k;
    rr = kk/ki;
    PSLB = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
    Q7 = 0.0;
    
    for (j=1; j<=nGL(pGL); j++) {
        xv = xGL(pGL)[j];
        w = wGL(pGL)[j];
        KQ7 = KQ7_function(ki, rr, xv);
        psl = psInterpolation_nr(ki*rsqrt(abskmq), kPS, pPS, nPSLT);
        Q7 += w*KQ7*psl;
    }
    
    return gd.p*PSLB*Q7;
}

local real Q7_function(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    result=rlog(10.0)*(rsqr(ki)/FOURPI2)
    *qromo(GaussLegendreQ7_func,ymin,ymax,midpnt,cmd.epsquad);
    
    return result;
}

// END DE Q7

// BEGIN Q11

local real KQ11_function(real ki, real rr, real xv)
{
    real KQ11;
    
    KQ11 = rsqr(xv);
    return KQ11;
}

local real GaussLegendreQ11_func(real y)
{
    int j;
    real ki;
    real Q11, KQ11;
    real PSLB;
    real kk, rr;
    real xv, w, psl;
    
    gd.p = rpow(10.0,y);
    
    kk = gd.p;
    ki = gd.k;
    rr = kk/ki;
    PSLB = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
    Q11 = 0.0;
    
    for (j=1; j<=nGL(pGL); j++) {
        xv = xGL(pGL)[j];
        w = wGL(pGL)[j];
        KQ11 = KQ11_function(ki, rr, xv);
        psl = psInterpolation_nr(ki*rsqrt(abskmq), kPS, pPS, nPSLT);
        Q11 += w*KQ11*psl;
    }
    
    return gd.p*PSLB*Q11;
}

local real Q11_function(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    result=rlog(10.0)*(rsqr(ki)/FOURPI2)
    *qromo(GaussLegendreQ11_func,ymin,ymax,midpnt,cmd.epsquad);
    
    return result;
}

// END DE Q11

// BEGIN Q12

local real KQ12_function(real ki, real rr, real xv)
{
    real KQ12;
    
    KQ12 = rr*xv;
    return KQ12;
}

local real GaussLegendreQ12_func(real y)
{
    int j;
    real ki;
    real Q12, KQ12;
    real PSLB;
    real kk, rr;
    real xv, w, psl;
    
    gd.p = rpow(10.0,y);
    
    kk = gd.p;
    ki = gd.k;
    rr = kk/ki;
    PSLB = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
    Q12 = 0.0;
    
    for (j=1; j<=nGL(pGL); j++) {
        xv = xGL(pGL)[j];
        w = wGL(pGL)[j];
        KQ12 = KQ12_function(ki, rr, xv);
        psl = psInterpolation_nr(ki*rsqrt(abskmq), kPS, pPS, nPSLT);
        Q12 += w*KQ12*psl;
    }
    
    return gd.p*PSLB*Q12;
}

local real Q12_function(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    result=rlog(10.0)*(rsqr(ki)/FOURPI2)
    *qromo(GaussLegendreQ12_func,ymin,ymax,midpnt,cmd.epsquad);
    
    return result;
}

// END DE Q12


// BEGIN R1

global real funcR1int(real y)
{
    real ftmp,fac;
    global_D3v2_ptr ptmp;
    
    gd.p = rpow(10.0,y);
    fac = rpow(gd.p,3.0)*psInterpolation_nr(gd.p, kPS, pPS, nPSLT);

    ptmp = DsThirdOrder_func(gd.x, gd.k, gd.p);
    ftmp = (21.0/10.0)*D3symmD3v2(ptmp)/( DpkD3v2(ptmp)*DppD3v2(ptmp)*DppD3v2(ptmp) );

    return fac*ftmp;
}

local real R1_function(real eta, real ki)
{
    real result, fac, s;
    real pmin, pmax, ymin, ymax;
    int j;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);

    fac=(rlog(10.0)/FOURPI2)*psInterpolation_nr(gd.k, kPS, pPS, nPSLT);
    
    s=0;
    for (j=1;j<=nGL(pGL)/2;j++) {
        gd.x = xGL(pGL)[j];
        result = qromo(funcR1int,ymin,ymax,midpnt,cmd.epsquad);
        s += 2.0*wGL(pGL)[j]*result;
    }

    return fac*s;
}

// END DE R1


// BEGIN R2

local real GaussLegendreR2_func(real y)
{
    real ss, fac;
    real R2aB, rr, xv, w, k2, KA, KB, KR2;
    global_D2v2_ptr ptmp;
    int j;
    
    gd.p = rpow(10.0,y);
    rr = gd.p/gd.k;
    fac = gd.p*psInterpolation_nr(gd.p, kPS, pPS, nPSLT);

    R2aB = 0.0;
    for (j=1; j<=nGL(pGL); j++) {
        xv = xGL(pGL)[j];
        w = wGL(pGL)[j];
        k2 = gd.k * rsqrt(1.0 + rsqr(rr) - 2.0*rr*xv);
        ptmp = DsSecondOrder_func(k2, gd.k, gd.p);
        if (model_int_flag==LCDM) {
            KR2 = (rr*xv*(1.0-rr*xv)/abskmq)
                *(
                  KA_LCDM - KB_LCDM*rsqr(xv)
                  );
        } else {
            KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
            KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
            KR2 = (rr*xv*(1.0-rr*xv)/abskmq)
                *(
                  KA - KB*rsqr(xv)
                  );
        }
        R2aB += w*KR2;
    }

    return fac*R2aB;
}

local real R2_function(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);

    result=rlog(10.0)*(rsqr(gd.k)/FOURPI2)*psInterpolation_nr(gd.k, kPS, pPS, nPSLT)
            *qromo(GaussLegendreR2_func,ymin,ymax,midpnt,cmd.epsquad);

    return result;
}

// END DE R2

// BEGIN RI

local real KRI_function(real ki, real rr, real xv)
{
    real k2, KA, KB, KRI;
    global_D2v2_ptr ptmp;
    
    k2 = ki * rsqrt(abskmq);
    ptmp = DsSecondOrder_func(k2,ki,ki*rr);
    if (model_int_flag==LCDM) {
        KRI = (rsqr(rr)*(1.0-rsqr(xv))/abskmq)
            *(
              KA_LCDM - KB_LCDM*rsqr(xv)
              );
    } else {
        KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KRI = (rsqr(rr)*(1.0-rsqr(xv))/abskmq)
                *(
                    KA - KB*rsqr(xv)
                  );
    }
    return KRI;
}

local real GaussLegendreRI_func(real y)
{
    int j;
    real ki;
    real RI, KRI;
    real PSL;
    real kk, rr;
    real xv, w;
    
    gd.p = rpow(10.0,y);
    
    kk = gd.p;
    ki = gd.k;
    rr = kk/ki;
    PSL = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
    RI = 0.0;

    for (j=1; j<=nGL(pGL); j++) {
        xv = xGL(pGL)[j];
        w = wGL(pGL)[j];
        KRI = KRI_function(ki, rr, xv);
        RI += w*KRI;
    }
    
    return gd.p*RI*PSL;
}

local real RI_function(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    real PSL;
    
    gd.xstop = eta;
    gd.k = ki;
    PSL = psInterpolation_nr(ki, kPS, pPS, nPSLT);
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    result=rlog(10.0)*(rsqr(ki)/FOURPI2)*PSL
        *qromo(GaussLegendreRI_func,ymin,ymax,midpnt,cmd.epsquad);

    return result;
}

// END DE RI


// BEGIN R1p2

local real KR1p2_function(real ki, real rr, real xv)
{
    real k2, KA, KB, KR1p2;
    global_D2v2_ptr ptmp;
    
    k2 = ki * rsqrt(abskmq);
    ptmp = DsSecondOrder_func(k2,ki,ki*rr);
    if (model_int_flag==LCDM) {
        KR1p2 = (rsqr(rr)*(1.0-rr*xv)/abskmq)
                *(
                  KA_LCDM - KB_LCDM*rsqr(xv)
                  );
    } else {
        KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KR1p2 = (rsqr(rr)*(1.0-rr*xv)/abskmq)
                *(
                  KA - KB*rsqr(xv)
                 );
    }
    return KR1p2;
}

local real GaussLegendreR1p2_func(real y)
{
    int j;
    real ki;
    real R1p2, KR1p2;
    real PSL;
    real kk, rr;
    real xv, w;
    
    gd.p = rpow(10.0,y);
    
    kk = gd.p;
    ki = gd.k;
    rr = kk/ki;
    PSL = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
    R1p2 = 0.0;
    
    for (j=1; j<=nGL(pGL); j++) {
        xv = xGL(pGL)[j];
        w = wGL(pGL)[j];
        KR1p2 = KR1p2_function(ki, rr, xv);
        R1p2 += w*KR1p2;
    }
    
    return gd.p*R1p2*PSL;
}

local real R1p2_function(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    real PSL;

    gd.xstop = eta;
    gd.k = ki;
    PSL = psInterpolation_nr(ki, kPS, pPS, nPSLT);

    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);

    result=rlog(10.0)*(rsqr(ki)/FOURPI2)*PSL
        *qromo(GaussLegendreR1p2_func,ymin,ymax,midpnt,cmd.epsquad);

    return result;
}

// END DE R1p2



// BEGIN Qs and Rs

global_QRs QsRs_functions_trapezoid1(real eta, real ki)
{
    int i, j;
    real ftmpR1, fac;
    global_D2v2_ptr ptmp;
    global_D3v2_ptr ptmpR1;
    real Dpkmin, Dpk;
    real aTime;
    //
    real *xGL, *wGL;
    real kmin, kmax;
    real Q1p, Q2p, Q3p;
    real Q1aA, Q2aA, Q3aA;
    real Q1aB, Q2aB, Q3aB;
    real R1aA, R1aB, R1p;
    real R2p, R2aA, R2aB;
    real KR2;
    real PSLA, PSLB;
    real rmin, rmax;
    real kk, rr, deltar;
    real mumin, mumax;
    real xv, w, k2, psl;
    real KA, KB, KQ1, KQ2, KQ3;
    int Nx;
    //
    // NEW Qs AND Rs
    real Q8p, Q8aA, Q8aB, KQ8;
    real Q9p, Q9aA, Q9aB, KQ9;
    real Q13p, Q13aA, Q13aB, KQ13;
    real Q5p, Q5aA, Q5aB, KQ5;
    real Q7p, Q7aA, Q7aB, KQ7;
    real Q11p, Q11aA, Q11aB, KQ11;
    real Q12p, Q12aA, Q12aB, KQ12;
    real RIp, RIaA, RIaB, KRI;
    real R1p2p, R1p2aA, R1p2aB, KR1p2;
    
    //
    pointPSTableptr p;
    //
    global_QRs_ptr QRstmp;
    
    QRstmp = (global_QRs_ptr) allocate(1 * sizeof(global_QRs));
    
    aTime = cputime();
    kmin = kPos(PSLT+1);
    kmax = kPos(PSLT+nPSLT-1);
    fprintf(gd.outlog,"\ntrapezoid1 :: nPSLT, kmin and kmax :: %d %g %g",nPSLT,kmin,kmax);
    fac = psInterpolation_nr(ki, kPS, pPS, nPSLT);
    
    Q1p = 0.0; Q2p = 0.0; Q3p = 0.0;
    Q1aA = 0.0; Q2aA = 0.0; Q3aA = 0.0;
    Q1aB = 0.0; Q2aB = 0.0; Q3aB = 0.0;
    
    // NEW Qs AND Rs
    Q8p = 0.0;
    Q8aA = 0.0;
    Q8aB = 0.0;
    Q9p = 0.0;
    Q9aA = 0.0;
    Q9aB = 0.0;
    Q13p = 0.0;
    Q13aA = 0.0;
    Q13aB = 0.0;
    Q5p = 0.0;
    Q5aA = 0.0;
    Q5aB = 0.0;
    Q7p = 0.0;
    Q7aA = 0.0;
    Q7aB = 0.0;
    Q11p = 0.0;
    Q11aA = 0.0;
    Q11aB = 0.0;
    Q12p = 0.0;
    Q12aA = 0.0;
    Q12aB = 0.0;
    RIp = 0.0;
    RIaA = 0.0;
    RIaB = 0.0;
    R1p2p = 0.0;
    R1p2aA = 0.0;
    R1p2aB = 0.0;
    //
    R1p = 0.0;
    R1aA = 0.0;
    R1aB = 0.0;
    
    R2p = 0.0;
    R2aA = 0.0;
    R2aB = 0.0;
    PSLA = 0.0;
    rmax = kmax/ki;
    rmin = kmin/ki;
    p = PSLCDMtab;
    //
    for (i=1; i<nPSTable; i++) {
        kk = kPos(p+i);
        rr = kk/ki;
        PSLB = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
        mumin = MAX(-1.0, (1.0 + rsqr(rr) - rsqr(rmax))/(2.0*rr));
        mumax = MIN( 1.0, (1.0 + rsqr(rr) - rsqr(rmin))/(2.0*rr));
        if (rr>=0.5)
            mumax = 0.5/rr;
        Nx=10;
        xGL=dvector(1,Nx);
        wGL=dvector(1,Nx);
        gauleg(mumin,mumax,xGL,wGL,Nx);
        for (j=1; j<=Nx; j++) {
            xv = xGL[j];
            w = wGL[j];
            k2 = ki * rsqrt(1.0 + rsqr(rr) - 2.0*rr*xv);
            ptmp = DsSecondOrder_func(ki, ki*rr, k2);
            KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
            KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
            KQ1 = rsqr(rr)
            *rsqr(
                  KA + KB*(-1.0+(1.0-rsqr(xv))/abskmq)
                  );
            KQ2 = (rr*xv*(1.0-rr*xv)/abskmq)
            *(
              KA - KB*((rsqr(xv)+rsqr(rr)-2.0*rr*xv)/abskmq)
              );
            KQ3 = rsqr(xv)*rsqr(1.0-rr*xv)/rsqr(abskmq);
            
            // NEW Qs AND Rs
            KQ8 = rsqr(rr)*(
                            KA - KB*rsqr(-rr+xv)/abskmq
                            );
            KQ9 = rr*xv*(1-rr*xv)/abskmq;
            KQ13 = rsqr(rr);
            //
            psl = psInterpolation_nr(ki*rsqrt(abskmq), kPS, pPS, nPSLT);
            Q1aB += wGL[j]*KQ1*psl;
            Q2aB += wGL[j]*KQ2*psl;
            Q3aB += wGL[j]*KQ3*psl;
            
            // NEW Qs AND Rs
            Q8aB += w*KQ8*psl;
            Q9aB += w*KQ9*psl;
            Q13aB += w*KQ13*psl;
            //
        }
        //
        for (j=1; j<=nGL(pGL); j++) {
            xv = xGL(pGL)[j];
            w = wGL(pGL)[j];
            k2 = ki * rsqrt(1.0 + rsqr(rr) - 2.0*rr*xv);
            ptmp = DsSecondOrder_func(k2, ki, ki*rr);
            KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
            KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
            KR2 = (rr*xv*(1.0-rr*xv)/abskmq)
            *(
              KA - KB*rsqr(xv)
              );
            psl = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
            R2aB += w*KR2*psl;
            // NEW Qs AND Rs
            KRI = (rsqr(rr)*(1.0-rsqr(xv))/abskmq)
            *(
              KA - KB*rsqr(xv)
              );
            RIaB += w*KRI*psl;
            KR1p2 = (rsqr(rr)*(1.0-rr*xv)/abskmq)
            *(
              KA - KB*rsqr(xv)
              );
            R1p2aB += w*KR1p2*psl;
            
            psl = psInterpolation_nr(ki*rsqrt(abskmq), kPS, pPS, nPSLT);
            ptmp = DsSecondOrder_func(ki, ki*rr, k2);
            KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
            KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
            KQ5 = rr*xv*(
                         KA - KB*rsqr(-rr+xv)/abskmq
                         );
            Q5aB += w*KQ5*psl;
            
            KQ7 = rsqr(xv)*(1-rr*xv)/abskmq;
            Q7aB += w*KQ7*psl;
            
            KQ11 = rsqr(xv);
            Q11aB += w*KQ11*psl;
            
            KQ12 = rr*xv;
            Q12aB += w*KQ12*psl;
            
            //
            ptmpR1 = DsThirdOrder_func(xGL(pGL)[j], ki, kk);
            ftmpR1 = (21.0/10.0)*D3symmD3v2(ptmpR1)
            /( DpkD3v2(ptmpR1)*DppD3v2(ptmpR1)*DppD3v2(ptmpR1) );
            //
            R1aB += rsqr(rr)*psInterpolation_nr(kk, kPS, pPS, nPSLT)
            *wGL(pGL)[j]*ftmpR1;
        }
        //
        deltar = (kPos(p+i)-kPos(p+i-1))/ki;
        
        Q1p += deltar*(Q1aA*PSLA + Q1aB*PSLB)/2.0;
        Q2p += deltar*(Q2aA*PSLA + Q2aB*PSLB)/2.0;
        Q3p += deltar*(Q3aA*PSLA + Q3aB*PSLB)/2.0;
        
        // NEW Qs AND Rs
        Q8p += deltar*(Q8aA*PSLA + Q8aB*PSLB)/2.0;
        Q9p += deltar*(Q9aA*PSLA + Q9aB*PSLB)/2.0;
        Q13p += deltar*(Q13aA*PSLA + Q13aB*PSLB)/2.0;
        Q5p += deltar*(Q5aA*PSLA + Q5aB*PSLB)/2.0;
        Q7p += deltar*(Q7aA*PSLA + Q7aB*PSLB)/2.0;
        Q11p += deltar*(Q11aA*PSLA + Q11aB*PSLB)/2.0;
        Q12p += deltar*(Q12aA*PSLA + Q12aB*PSLB)/2.0;
        RIp += deltar*(RIaA + RIaB)/2.0;
        R1p2p += deltar*(R1p2aA + R1p2aB)/2.0;
        //
        Q1aA = Q1aB;
        Q2aA = Q2aB;
        Q3aA = Q3aB;
        
        // NEW Qs AND Rs
        Q8aA = Q8aB;
        Q9aA = Q9aB;
        Q13aA = Q13aB;
        Q5aA = Q5aB;
        Q7aA = Q7aB;
        Q11aA = Q11aB;
        Q12aA = Q12aB;
        RIaA = RIaB;
        R1p2aA = R1p2aB;
        //
        PSLA = PSLB;
        Q1aB = 0.0;
        Q2aB = 0.0;
        Q3aB = 0.0;
        
        // NEW Qs AND Rs
        Q8aB = 0.0;
        Q9aB = 0.0;
        Q13aB = 0.0;
        Q5aB = 0.0;
        Q7aB = 0.0;
        Q11aB = 0.0;
        Q12aB = 0.0;
        RIaB = 0.0;
        R1p2aB = 0.0;
        //
        R1p += deltar*( R1aA + R1aB )/2.0;
        R1aA = R1aB;
        R1aB = 0.0;
        
        R2p += deltar*(R2aA + R2aB)/2.0;
        R2aA = R2aB;
        R2aB = 0.0;
        
        free_dvector(wGL,1,Nx);
        free_dvector(xGL,1,Nx);
    }
    Q1p *= 2.0*(rpow(ki,3)/FOURPI2);
    Q2p *= 2.0*(rpow(ki,3)/FOURPI2);
    Q3p *= 2.0*(rpow(ki,3)/FOURPI2);
    
    // NEW Qs AND Rs
    Q8p *= 2.0*(rpow(ki,3)/FOURPI2);
    Q9p *= 2.0*(rpow(ki,3)/FOURPI2);
    Q13p *= 2.0*(rpow(ki,3)/FOURPI2);
    Q5p *= (rpow(ki,3)/FOURPI2);
    Q7p *= (rpow(ki,3)/FOURPI2);
    Q11p *= (rpow(ki,3)/FOURPI2);
    Q12p *= (rpow(ki,3)/FOURPI2);
    RIp *= (rpow(ki,3)/FOURPI2)*fac;
    R1p2p *= (rpow(ki,3)/FOURPI2)*fac;
    //
    
    R1p *=(rpow(gd.k,3)/FOURPI2)*fac;
    R2p *= (rpow(ki,3)/FOURPI2)*fac;
    
    etaQRs(QRstmp) = eta;
    kQRs(QRstmp)    = ki;
    Q1(QRstmp)      = Q1p;
    Q2(QRstmp)      = Q2p;
    Q3(QRstmp)      = Q3p;
    // NEW Qs AND Rs
    Q8(QRstmp)      = Q8p;
    Q9(QRstmp)      = Q9p;
    Q13(QRstmp)     = Q13p;
    Q5(QRstmp)      = Q5p;
    Q7(QRstmp)      = Q7p;
    Q11(QRstmp)      = Q11p;
    Q12(QRstmp)      = Q12p;
    RI(QRstmp)      = RIp;
    R1p2(QRstmp)      = R1p2p;
    //
    R1(QRstmp)      = R1p;
    R2(QRstmp)      = R2p;
    
    return *QRstmp;
}

global_QRs QsRs_functions_trapezoid2(real eta, real ki)
{
    int i, j;
    real ftmpR1, fac;
    global_D2v2_ptr ptmp;
    //
    real *xGL, *wGL;
    real kmin, kmax;
    real Q1p, Q2p, Q3p;
    real Q1aA, Q2aA, Q3aA;
    real Q1aB, Q2aB, Q3aB;
    real R1p, R2p;
    real PSLA, PSLB;
    real rmin, rmax;
    real kk, rr, deltar;
    real mumin, mumax;
    real xv, w, k2, psl;
    real KA, KB, KQ1, KQ2, KQ3;
    int Nx;
    //
    // NEW Qs AND Rs
    real Q8p, Q8aA, Q8aB, KQ8;
    real Q9p, Q9aA, Q9aB, KQ9;
    real Q13p, Q13aA, Q13aB, KQ13;
    real QIp, QIaA, QIaB, KQI;
    real Q5p;
    real Q7p;
    real Q11p;
    real Q12p;
    real RIp;
    real R1p2p;
    //
    pointPSTableptr p;
    //
    global_QRs_ptr QRstmp;
    
    QRstmp = (global_QRs_ptr) allocate(1 * sizeof(global_QRs));
    
    kmin = kPos(PSLT+1);
    kmax = kPos(PSLT+nPSLT-1);
    fprintf(gd.outlog,"\ntrapezoid2 :: nPSLT, kmin and kmax :: %d %g %g",nPSLT,kmin,kmax);
    fac = psInterpolation_nr(ki, kPS, pPS, nPSLT);
    
    Q1p = 0.0; Q2p = 0.0; Q3p = 0.0;
    Q1aA = 0.0; Q2aA = 0.0; Q3aA = 0.0;
    Q1aB = 0.0; Q2aB = 0.0; Q3aB = 0.0;
    // NEW Qs AND Rs
    Q8p = 0.0;
    Q8aA = 0.0;
    Q8aB = 0.0;
    Q9p = 0.0;
    Q9aA = 0.0;
    Q9aB = 0.0;
    Q13p = 0.0;
    Q13aA = 0.0;
    Q13aB = 0.0;
    //
    PSLA = 0.0;
    rmax = kmax/ki;
    rmin = kmin/ki;
    p = PSLCDMtab;
    //
    for (i=1; i<nPSTable; i++) {
        kk = kPos(p+i);
        rr = kk/ki;
        PSLB = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
        mumin = MAX(-1.0, (1.0 + rsqr(rr) - rsqr(rmax))/(2.0*rr));
        mumax = MIN( 1.0, (1.0 + rsqr(rr) - rsqr(rmin))/(2.0*rr));
        if (rr>=0.5)
            mumax = 0.5/rr;
        Nx=10;
        xGL=dvector(1,Nx);
        wGL=dvector(1,Nx);
        gauleg(mumin,mumax,xGL,wGL,Nx);
        for (j=1; j<=Nx; j++) {
            xv = xGL[j];
            w = wGL[j];
            k2 = ki * rsqrt(1.0 + rsqr(rr) - 2.0*rr*xv);
            ptmp = DsSecondOrder_func(ki, ki*rr, k2);
            KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
            KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
            KQ1 = rsqr(rr)
            *rsqr(
                  KA + KB*(-1.0+(1.0-rsqr(xv))/abskmq)
                  );
            KQ2 = (rr*xv*(1.0-rr*xv)/abskmq)
            *(
              KA - KB*((rsqr(xv)+rsqr(rr)-2.0*rr*xv)/abskmq)
              );
            KQ3 = rsqr(xv)*rsqr(1.0-rr*xv)/rsqr(abskmq);
            // NEW Qs AND Rs
            KQ8 = rsqr(rr)*(
                            KA - KB*rsqr(-rr+xv)/abskmq
                            );
            KQ9 = rr*xv*(1-rr*xv)/abskmq;
            KQ13 = rsqr(rr);
            //
            psl = psInterpolation_nr(ki*rsqrt(abskmq), kPS, pPS, nPSLT);
            Q1aB += wGL[j]*KQ1*psl;
            Q2aB += wGL[j]*KQ2*psl;
            Q3aB += wGL[j]*KQ3*psl;
            // NEW Qs AND Rs
            Q8aB += w*KQ8*psl;
            Q9aB += w*KQ9*psl;
            Q13aB += w*KQ13*psl;
            //
        }
        //
        deltar = (kPos(p+i)-kPos(p+i-1))/ki;
        
        Q1p += deltar*(Q1aA*PSLA + Q1aB*PSLB)/2.0;
        Q2p += deltar*(Q2aA*PSLA + Q2aB*PSLB)/2.0;
        Q3p += deltar*(Q3aA*PSLA + Q3aB*PSLB)/2.0;
        // NEW Qs AND Rs
        Q8p += deltar*(Q8aA*PSLA + Q8aB*PSLB)/2.0;
        Q9p += deltar*(Q9aA*PSLA + Q9aB*PSLB)/2.0;
        Q13p += deltar*(Q13aA*PSLA + Q13aB*PSLB)/2.0;
        //
        Q1aA = Q1aB;
        Q2aA = Q2aB;
        Q3aA = Q3aB;
        // NEW Qs AND Rs
        Q8aA = Q8aB;
        Q9aA = Q9aB;
        Q13aA = Q13aB;
        //
        PSLA = PSLB;
        Q1aB = 0.0;
        Q2aB = 0.0;
        Q3aB = 0.0;
        // NEW Qs AND Rs
        Q8aB = 0.0;
        Q9aB = 0.0;
        Q13aB = 0.0;
        QIaB = 0.0;
        //
        free_dvector(wGL,1,Nx);
        free_dvector(xGL,1,Nx);
    }
    Q1p *= 2.0*(rpow(ki,3)/FOURPI2);
    Q2p *= 2.0*(rpow(ki,3)/FOURPI2);
    Q3p *= 2.0*(rpow(ki,3)/FOURPI2);
    // NEW Qs AND Rs
    Q8p *= 2.0*(rpow(ki,3)/FOURPI2);
    Q9p *= 2.0*(rpow(ki,3)/FOURPI2);
    Q13p *= 2.0*(rpow(ki,3)/FOURPI2);
    QIp *= 2.0*(rpow(ki,3)/FOURPI2);
    //
    //
    Q5p = Q5_function(eta, ki);
    Q7p = Q7_function(eta, ki);
    Q11p = Q11_function(eta, ki);
    Q12p = Q12_function(eta, ki);
    RIp = RI_function(eta, ki);
    R1p2p = R1p2_function(eta, ki);
    //
    R1p = R1_function(eta, ki);        // CHECAR EL DESEMPEÑO DE ESTA VERSION
    R2p = R2_function(eta, ki);
    
    etaQRs(QRstmp) = eta;
    kQRs(QRstmp)    = ki;
    Q1(QRstmp)      = Q1p;
    Q2(QRstmp)      = Q2p;
    Q3(QRstmp)      = Q3p;
    // NEW Qs AND Rs
    Q8(QRstmp)      = Q8p;
    Q9(QRstmp)      = Q9p;
    Q13(QRstmp)     = Q13p;
    Q5(QRstmp)      = Q5p;
    Q7(QRstmp)      = Q7p;
    Q11(QRstmp)      = Q11p;
    Q12(QRstmp)      = Q12p;
    RI(QRstmp)      = RIp;
    R1p2(QRstmp)      = R1p2p;
    //
    R1(QRstmp)      = R1p;
    R2(QRstmp)      = R2p;
    
    return *QRstmp;
}



global_QRs QsRs_functions_romo(real eta, real ki)
{
    real Q1p, Q2p, Q3p, Q8p, Q9p, Q13p, QIp;
    real Q5p, Q7p, Q11p, Q12p;
    real R1p, R2p;
    real RI, R1p2;
    global_QRs_ptr QRstmp;

    QRstmp = (global_QRs_ptr) allocate(1 * sizeof(global_QRs));

    Q1p = Q1_function(eta, ki);
    Q2p = Q2_function(eta, ki);
    Q3p = Q3_function(eta, ki);
    Q8p = Q8_function(eta, ki);
    Q9p = Q9_function(eta, ki);
    Q13p = Q13_function(eta, ki);
    QIp = QI_function(eta, ki);
    Q5p = Q5_function(eta, ki);
    Q7p = Q7_function(eta, ki);
    Q11p = Q11_function(eta, ki);
    Q12p = Q12_function(eta, ki);
    RI = RI_function(eta, ki);
    R1p2 = R1p2_function(eta, ki);

    R1p = R1_function(eta, ki);
    R2p = R2_function(eta, ki);
    
    etaQRs(QRstmp) = eta;
    kQRs(QRstmp)    = ki;
    Q1(QRstmp)      = Q1p;
    Q2(QRstmp)      = Q2p;
    Q3(QRstmp)      = Q3p;
    Q8(QRstmp)      = Q8p;
    Q9(QRstmp)      = Q9p;
    Q13(QRstmp)      = Q13p;
    QI(QRstmp)      = QIp;
    Q5(QRstmp)      = Q5p;
    Q7(QRstmp)      = Q7p;
    Q11(QRstmp)      = Q11p;
    Q12(QRstmp)      = Q12p;
    RI(QRstmp)      = RI;
    R1p2(QRstmp)      = R1p2;
    R1(QRstmp)      = R1p;
    R2(QRstmp)      = R2p;

    return *QRstmp;
}

/*
global_QRs QsRs_functions_LCDM(real eta, real ki)
{
    real Q1p, Q2p, Q3p, Q8p, Q9p, Q13p, QIp;
    real Q5p, Q7p, Q11p, Q12p;
    real R1p, R2p;
    real RI, R1p2;
    global_QRs_ptr QRstmp;
    
    QRstmp = (global_QRs_ptr) allocate(1 * sizeof(global_QRs));
    
    Q1p = Q1_function(eta, ki);
    Q2p = Q2_function(eta, ki);
    Q3p = Q3_function(eta, ki);
    Q8p = Q8_function(eta, ki);
    Q9p = Q9_function(eta, ki);
    Q13p = Q13_function(eta, ki);
    QIp = QI_function(eta, ki);
    Q5p = Q5_function(eta, ki);
    Q7p = Q7_function(eta, ki);
    Q11p = Q11_function(eta, ki);
    Q12p = Q12_function(eta, ki);
    RI = RI_function(eta, ki);
    R1p2 = R1p2_function(eta, ki);
    
    R1p = RI;
    R2p = R2_function(eta, ki);
    
    etaQRs(QRstmp) = eta;
    kQRs(QRstmp)    = ki;
    Q1(QRstmp)      = Q1p;
    Q2(QRstmp)      = Q2p;
    Q3(QRstmp)      = Q3p;
    Q8(QRstmp)      = Q8p;
    Q9(QRstmp)      = Q9p;
    Q13(QRstmp)      = Q13p;
    QI(QRstmp)      = QIp;
    Q5(QRstmp)      = Q5p;
    Q7(QRstmp)      = Q7p;
    Q11(QRstmp)      = Q11p;
    Q12(QRstmp)      = Q12p;
    RI(QRstmp)      = RI;
    R1p2(QRstmp)      = R1p2;
    R1(QRstmp)      = R1p;
    R2(QRstmp)      = R2p;
    
    return *QRstmp;
}
*/

// END DE Qs and Rs

global_QRs qrs;

global global_QRs QsRs_functions_driver(real eta, real ki) // Remove eta
{
//    global_QRs qrs;

//    qrs = QsRs_functions(gd.xstop, ki);
    quadrature(ki);

    return qrs;
}

/*
global global_QRs QsRs_functions_driver_LCDM(real eta, real ki)
{
    global_QRs qrs;
    
    qrs = QsRs_functions_LCDM(gd.xstop, ki);
    
    return qrs;
}
*/

#define ROMO 1
#define NULLMETHOD 0
#define TRAPEZOID1 2
#define TRAPEZOID2 3
#define TRAPEZOIDAA 4

local void quadrature(real ki)
{
    //    global_QRs qrs;
    switch (gd.quadmethod_int) {
        case ROMO:
            qrs = QsRs_functions_romo(gd.xstop, ki);
            break;
            //
        case TRAPEZOID1:
            qrs = QsRs_functions_trapezoid1(gd.xstop, ki);
            break;
            //
        case TRAPEZOID2:
            qrs = QsRs_functions_trapezoid2(gd.xstop, ki);
            break;
            //
        case TRAPEZOIDAA:
            qrs = QsRs_functions_trapezoidAA(gd.xstop, ki);
            break;
            //
        case NULLMETHOD:
            qrs = QsRs_functions_romo(gd.xstop, ki);
            break;
            //
        default:
            qrs = QsRs_functions_romo(gd.xstop, ki);
            break;
    }
    //    return qrs;
}

void quadraturemethod_string_to_int(string method_str,int *method_int)
{
    *method_int=-1;
    if (strcmp(method_str,"romo") == 0) {
        *method_int = ROMO;
        strcpy(gd.quadraturemethod_comment, "romo quadrature method");
    }
    //
    if (strcmp(method_str,"trapezoid1") == 0) {
        *method_int = TRAPEZOID1;
        strcpy(gd.quadraturemethod_comment, "trapezoid1 quadrature method");
    }
    //
    if (strcmp(method_str,"trapezoid2") == 0) {
        *method_int = TRAPEZOID2;
        strcpy(gd.quadraturemethod_comment, "trapezoid2 quadrature method");
    }
    //
    if (strcmp(method_str,"trapezoidAA") == 0) {
        *method_int = TRAPEZOIDAA;
        strcpy(gd.quadraturemethod_comment, "trapezoidAA quadrature method");
    }
    //
    if (strnull(method_str)) {
        *method_int = NULLMETHOD;
        strcpy(gd.quadraturemethod_comment,
               "given null quadrature method ... running deafult (romo)");
        fprintf(stdout,"\n\tintegration: default integration method (romo)...\n");
    }
    //
    if (*method_int == -1) {
        *method_int = ROMO;
        strcpy(gd.quadraturemethod_comment,
               "Unknown quadrature method ... running deafult (romo)");
        fprintf(stdout,"\n\tquadrature: Unknown method... %s ",cmd.quadratureMethod);
        fprintf(stdout,
                "\n\trunning default quadrature method (romo)...\n");
    }
}

#undef ROMO
#undef TRAPEZOID1
#undef TRAPEZOID2
#undef TRAPEZOIDAA
#undef NULLMETHOD

#undef abskmq


// ==============================================================================
// AA TRAPEZOID METHOD:

//global_QRs QsRs_functions(real eta, real ki);
//global_QRs QsRs_functions_LCDM(real eta, real ki);

local real KQ8_function_AA(real k, real r, real x);
local real KQ9_function_AA(real k, real r, real x);
local real KQ13_function_AA(real k, real r, real x);
local real KQI_function_AA(real k, real r, real x);

local real KQ5_function_AA(real k, real r, real x);
local real KQ7_function_AA(real k, real r, real x);
local real KQ11_function_AA(real k, real r, real x);
local real KQ12_function_AA(real k, real r, real x);

local real KRI_function_AA(real k, real r, real x);
local real KR1p2_function_AA(real k, real r, real x);

local real GaussLegendreQ1_func_AA(real y);
local real GaussLegendreQ2_func_AA(real y);
local real GaussLegendreQ3_func_AA(real y);
local real GaussLegendreQ8_func_AA(real y);
local real GaussLegendreQ9_func_AA(real y);
local real GaussLegendreQ13_func_AA(real y);
local real GaussLegendreQI_func_AA(real y);
local real GaussLegendreQ5_func_AA(real y);
local real GaussLegendreQ7_func_AA(real y);
local real GaussLegendreQ11_func_AA(real y);
local real GaussLegendreQ12_func_AA(real y);
local real GaussLegendreRI_func_AA(real y);
local real GaussLegendreR1p2_func_AA(real y);
local real GaussLegendreR1_func_AA(real y);
local real GaussLegendreR2_func_AA(real y);
local real funcR1int_AA(real y);

local real Q1_function_AA(real eta, real ki);
local real Q2_function_AA(real eta, real ki);
local real Q3_function_AA(real eta, real ki);
local real Q8_function_AA(real eta, real ki);
local real Q9_function_AA(real eta, real ki);
local real Q13_function_AA(real eta, real ki);
local real QI_function_AA(real eta, real ki);
local real Q5_function_AA(real eta, real ki);
local real Q7_function_AA(real eta, real ki);
local real Q11_function_AA(real eta, real ki);
local real Q12_function_AA(real eta, real ki);
local real RI_function_AA(real eta, real ki);
local real R1p2_function_AA(real eta, real ki);
local real R1_function_AA(real eta, real ki);
local real R2_function_AA(real eta, real ki);

#define abskmq      (1.0+rsqr(rr)-2.0*rr*xv)
#define Nsteps      100   //Modificacion
//#define Nsteps      200   //Modificacion


// BEGIN Q1

local real GaussLegendreQ1_func_AA(real y)
{
    global_D2v2_ptr ptmp;
    int j;
    real *xGL, *wGL;
    real kmin, kmax, ki;
    real Q1p, Q1aA, Q1aB, KQ1;
    real PSLA, PSLB;
    real rmin, rmax;
    real kk, rr, deltar;
    real mumin, mumax;
    real xv, w, k2, psl;
    real KA, KB;
    int Nx;
    
    gd.p = rpow(10.0,y);
    
    kk = gd.p;
    ki = gd.k;
    rr = kk/ki;
    PSLB = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
    Q1aB = 0.0;
    kmin = kPos(PSLT+1);
    kmax = kPos(PSLT+nPSLT-1);
    rmax = kmax/ki;
    rmin = kmin/ki;
    
    mumin = MAX(-1.0, (1.0 + rsqr(rr) - rsqr(rmax))/(2.0*rr));
    mumax = MIN( 1.0, (1.0 + rsqr(rr) - rsqr(rmin))/(2.0*rr));
    if (rr>=0.5)
        mumax = 0.5/rr;
    Nx=10;
    xGL=dvector(1,Nx);
    wGL=dvector(1,Nx);
    gauleg(mumin,mumax,xGL,wGL,Nx);
    for (j=1; j<=Nx; j++) {
        xv = xGL[j];
        w = wGL[j];
        k2 = ki * rsqrt(1.0 + rsqr(rr) - 2.0*rr*xv);
        ptmp = DsSecondOrder_func(ki, ki*rr, k2);
        if (model_int_flag==LCDM) {
            KQ1 = rsqr(
                       KA_LCDM + KB_LCDM*(-1.0+(1.0-rsqr(xv))/abskmq)
                       );
        } else {
            KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
            KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
            KQ1 = rsqr(
                       KA + KB*(-1.0+(1.0-rsqr(xv))/abskmq)
                       );
        }
        psl = psInterpolation_nr(ki*rsqrt(abskmq), kPS, pPS, nPSLT);
        Q1aB += w*KQ1*psl;
    }
    
    free_dvector(wGL,1,Nx);
    free_dvector(xGL,1,Nx);
    
    return 2.0*rpow(gd.p,3.0)*PSLB*Q1aB;
}

local real Q1_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    real deltay;  //Modificacion
    //~ real Nsteps; //Modificacion
    int j; //Modificacion
    real Qa; //Modificacion
    real Qb;     //Modificacion
    real yp;    //Modificacion
    //~ Nsteps=400.0;  //Modificacion
    
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    
    deltay=(ymax-ymin)/Nsteps;  //Modificacion
    
    result=0; //Modificacion
    Qa=GaussLegendreQ1_func( ymin ); //Modificacion
    
    for (j=1; j<=Nsteps; j++) {     //Modificacion
        yp = ymin + deltay * j;
        Qb=GaussLegendreQ1_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }//Modificacion
    
    //~ result = (rlog(10.0)/FOURPI2)*result;  //Modificacion
    result = rlog(10.0)*result/FOURPI2;  //Modificacion
    
    
    //~ result=(rlog(10.0)/FOURPI2)
    //~ *qromo(GaussLegendreQ1_func,ymin,ymax,midpnt,cmd.epsquad);
    
    return result;
}

// END DE Q1


// BEGIN Q2

local real GaussLegendreQ2_func_AA(real y)
{
    global_D2v2_ptr ptmp;
    int j;
    real *xGL, *wGL;
    real kmin, kmax, ki;
    real Q2p, Q2aA, Q2aB, KQ2;
    real PSLA, PSLB;
    real rmin, rmax;
    real kk, rr, deltar;
    real mumin, mumax;
    real xv, w, k2, psl;
    real KA, KB;
    int Nx;
    
    gd.p = rpow(10.0,y);
    
    kk = gd.p;
    ki = gd.k;
    rr = kk/ki;
    PSLB = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
    Q2aB = 0.0;
    kmin = kPos(PSLT+1);
    kmax = kPos(PSLT+nPSLT-1);
    rmax = kmax/ki;
    rmin = kmin/ki;
    
    mumin = MAX(-1.0, (1.0 + rsqr(rr) - rsqr(rmax))/(2.0*rr));
    mumax = MIN( 1.0, (1.0 + rsqr(rr) - rsqr(rmin))/(2.0*rr));
    if (rr>=0.5)
        mumax = 0.5/rr;
    Nx=10;
    xGL=dvector(1,Nx);
    wGL=dvector(1,Nx);
    gauleg(mumin,mumax,xGL,wGL,Nx);
    for (j=1; j<=Nx; j++) {
        xv = xGL[j];
        w = wGL[j];
        k2 = ki * rsqrt(1.0 + rsqr(rr) - 2.0*rr*xv);
        ptmp = DsSecondOrder_func(ki, ki*rr, k2);
        if (model_int_flag==LCDM) {
            KQ2 = (rr*xv*(1.0-rr*xv)/abskmq)
            *(
              KA_LCDM - KB_LCDM*((rsqr(xv)+rsqr(rr)-2.0*rr*xv)/abskmq)
              );
        } else {
            KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
            KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
            KQ2 = (rr*xv*(1.0-rr*xv)/abskmq)
            *(
              KA - KB*((rsqr(xv)+rsqr(rr)-2.0*rr*xv)/abskmq)
              );
        }
        psl = psInterpolation_nr(ki*rsqrt(abskmq), kPS, pPS, nPSLT);
        Q2aB += w*KQ2*psl;
    }
    
    free_dvector(wGL,1,Nx);
    free_dvector(xGL,1,Nx);
    
    return 2.0*rpow(ki,2.0)*gd.p*PSLB*Q2aB;
}

local real Q2_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    real deltay;  //Modificacion
    //~ real Nsteps; //Modificacion
    int j; //Modificacion
    real Qa; //Modificacion
    real Qb;     //Modificacion
    real yp;    //Modificacion
    //~ Nsteps=400.0;  //Modificacion
    
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    
    deltay=(ymax-ymin)/Nsteps;  //Modificacion
    
    result=0; //Modificacion
    Qa=GaussLegendreQ2_func( ymin ); //Modificacion
    
    for (j=1; j<=Nsteps; j++) {     //Modificacion
        yp = ymin + deltay * j;
        Qb=GaussLegendreQ2_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }//Modificacion
    
    
    result = rlog(10.0)*result/FOURPI2;  //Modificacion
    
    
    //~ result=(rlog(10.0)/FOURPI2)
    //~ *qromo(GaussLegendreQ2_func,ymin,ymax,midpnt,cmd.epsquad);
    
    return result;
}

// END DE Q2


// BEGIN Q3

local real GaussLegendreQ3_func_AA(real y)
{
    int j;
    real *xGL, *wGL;
    real kmin, kmax, ki;
    real Q3p, Q3aA, Q3aB, KQ3;
    real PSLA, PSLB;
    real rmin, rmax;
    real kk, rr, deltar;
    real mumin, mumax;
    real xv, w, k2, psl;
    int Nx;
    
    gd.p = rpow(10.0,y);
    
    kk = gd.p;
    ki = gd.k;
    rr = kk/ki;
    PSLB = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
    Q3aB = 0.0;
    kmin = kPos(PSLT+1);
    kmax = kPos(PSLT+nPSLT-1);
    rmax = kmax/ki;
    rmin = kmin/ki;
    
    mumin = MAX(-1.0, (1.0 + rsqr(rr) - rsqr(rmax))/(2.0*rr));
    mumax = MIN( 1.0, (1.0 + rsqr(rr) - rsqr(rmin))/(2.0*rr));
    if (rr>=0.5)
        mumax = 0.5/rr;
    Nx=10;
    xGL=dvector(1,Nx);
    wGL=dvector(1,Nx);
    gauleg(mumin,mumax,xGL,wGL,Nx);
    for (j=1; j<=Nx; j++) {
        xv = xGL[j];
        w = wGL[j];
        k2 = ki * rsqrt(1.0 + rsqr(rr) - 2.0*rr*xv);
        KQ3 = rsqr(xv)*rsqr(1.0-rr*xv)/rsqr(abskmq);
        psl = psInterpolation_nr(ki*rsqrt(abskmq), kPS, pPS, nPSLT);
        Q3aB += w*KQ3*psl;
    }
    
    free_dvector(wGL,1,Nx);
    free_dvector(xGL,1,Nx);
    
    return 2.0*rpow(ki,2.0)*gd.p*PSLB*Q3aB;
}

local real Q3_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    real deltay;  //Modificacion
    //~ real Nsteps; //Modificacion
    int j; //Modificacion
    real Qa; //Modificacion
    real Qb;     //Modificacion
    real yp;    //Modificacion
    //~ Nsteps=400.0;  //Modificacion
    
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    
    deltay=(ymax-ymin)/Nsteps;  //Modificacion
    
    result=0; //Modificacion
    Qa=GaussLegendreQ3_func( ymin ); //Modificacion
    
    for (j=1; j<=Nsteps; j++) {     //Modificacion
        yp = ymin + deltay * j;
        Qb=GaussLegendreQ3_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }//Modificacion
    
    result = rlog(10.0)*result/FOURPI2;  //Modificacion
    
    
    //~ result=(rlog(10.0)/FOURPI2)
    //~ *qromo(GaussLegendreQ3_func,ymin,ymax,midpnt,cmd.epsquad);
    
    return result;
}

// END Q3


// BEGIN Q8

local real KQ8_function_AA(real ki, real rr, real xv)
{
    real k2, KA, KB, KQ8;
    global_D2v2_ptr ptmp;
    
    k2 = ki * rsqrt(abskmq);
    ptmp = DsSecondOrder_func(ki, ki*rr, k2);
    if (model_int_flag==LCDM) {
        KQ8 = rsqr(rr)*(
                        KA_LCDM - KB_LCDM*rsqr(-rr+xv)/abskmq
                        );
    } else {
        KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KQ8 = rsqr(rr)*(
                        KA - KB*rsqr(-rr+xv)/abskmq
                        );
    }
    return KQ8;
}

local real GaussLegendreQ8_func_AA(real y)
{
    int j;
    real *xGL, *wGL;
    real kmin, kmax, ki;
    real Q8, KQ8;
    real PSLA, PSLB;
    real rmin, rmax;
    real kk, rr;
    real xmin, xmax;
    real xv, w, psl;
    int Nx;
    
    gd.p = rpow(10.0,y);
    
    kk = gd.p;
    ki = gd.k;
    rr = kk/ki;
    PSLB = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
    Q8 = 0.0;
    kmin = kPos(PSLT+1);
    kmax = kPos(PSLT+nPSLT-1);
    rmax = kmax/ki;
    rmin = kmin/ki;
    
    xmin = MAX(-1.0, (1.0 + rsqr(rr) - rsqr(rmax))/(2.0*rr));
    xmax = MIN( 1.0, (1.0 + rsqr(rr) - rsqr(rmin))/(2.0*rr));
    if (rr>=0.5)
        xmax = 0.5/rr;
    Nx=10;
    xGL=dvector(1,Nx);
    wGL=dvector(1,Nx);
    gauleg(xmin,xmax,xGL,wGL,Nx);
    for (j=1; j<=Nx; j++) {
        xv = xGL[j];
        w = wGL[j];
        KQ8 = KQ8_function(ki, rr, xv);
        psl = psInterpolation_nr(ki*rsqrt(abskmq), kPS, pPS, nPSLT);
        Q8 += w*KQ8*psl;
    }
    
    free_dvector(wGL,1,Nx);
    free_dvector(xGL,1,Nx);
    
    return gd.p*PSLB*Q8;
}

local real Q8_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    real deltay;  //Modificacion
    //~ real Nsteps; //Modificacion
    int j; //Modificacion
    real Qa; //Modificacion
    real Qb;     //Modificacion
    real yp;    //Modificacion
    //~ Nsteps=400.0;  //Modificacion
    
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    
    deltay=(ymax-ymin)/Nsteps;  //Modificacion
    
    result=0; //Modificacion
    Qa=GaussLegendreQ8_func( ymin ); //Modificacion
    
    for (j=1; j<=Nsteps; j++) {     //Modificacion
        yp = ymin + deltay * j;
        Qb=GaussLegendreQ8_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }//Modificacion
    
    //~ result = (rlog(10.0)/FOURPI2)*result;  //Modificacion
    result = 2.0*rlog(10.0)*(rsqr(ki)/FOURPI2)*result;  //Modificacion
    
    
    //~ result=2.0*rlog(10.0)*(rsqr(ki)/FOURPI2)
    //~ *qromo(GaussLegendreQ8_func,ymin,ymax,midpnt,cmd.epsquad);
    
    return result;
}

// END DE Q8


// BEGIN Q9

local real KQ9_function_AA(real ki, real rr, real xv)
{
    real KQ9;
    
    KQ9 = rr*xv*(1-rr*xv)/abskmq;
    return KQ9;
}

local real GaussLegendreQ9_func_AA(real y)
{
    int j;
    real *xGL, *wGL;
    real kmin, kmax, ki;
    real Q9, KQ9;
    real PSLA, PSLB;
    real rmin, rmax;
    real kk, rr;
    real xmin, xmax;
    real xv, w, psl;
    int Nx;
    
    gd.p = rpow(10.0,y);
    
    kk = gd.p;
    ki = gd.k;
    rr = kk/ki;
    PSLB = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
    Q9 = 0.0;
    kmin = kPos(PSLT+1);
    kmax = kPos(PSLT+nPSLT-1);
    rmax = kmax/ki;
    rmin = kmin/ki;
    
    xmin = MAX(-1.0, (1.0 + rsqr(rr) - rsqr(rmax))/(2.0*rr));
    xmax = MIN( 1.0, (1.0 + rsqr(rr) - rsqr(rmin))/(2.0*rr));
    if (rr>=0.5)
        xmax = 0.5/rr;
    Nx=10;
    xGL=dvector(1,Nx);
    wGL=dvector(1,Nx);
    gauleg(xmin,xmax,xGL,wGL,Nx);
    for (j=1; j<=Nx; j++) {
        xv = xGL[j];
        w = wGL[j];
        KQ9 = KQ9_function(ki, rr, xv);
        psl = psInterpolation_nr(ki*rsqrt(abskmq), kPS, pPS, nPSLT);
        Q9 += w*KQ9*psl;
    }
    
    free_dvector(wGL,1,Nx);
    free_dvector(xGL,1,Nx);
    
    return gd.p*PSLB*Q9;
}

local real Q9_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    real deltay;  //Modificacion
    //~ real Nsteps; //Modificacion
    int j; //Modificacion
    real Qa; //Modificacion
    real Qb;     //Modificacion
    real yp;    //Modificacion
    //~ Nsteps=400.0;  //Modificacion
    
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    
    deltay=(ymax-ymin)/Nsteps;  //Modificacion
    
    result=0; //Modificacion
    Qa=GaussLegendreQ9_func( ymin ); //Modificacion
    
    for (j=1; j<=Nsteps; j++) {     //Modificacion
        yp = ymin + deltay * j;
        Qb=GaussLegendreQ9_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }//Modificacion
    
    
    result = 2.0*rlog(10.0)*(rsqr(ki)/FOURPI2)*result;  //Modificacion
    
    //~ result=2.0*rlog(10.0)*(rsqr(ki)/FOURPI2)
    //~ *qromo(GaussLegendreQ9_func,ymin,ymax,midpnt,cmd.epsquad);
    
    return result;
}

// END DE Q9

// BEGIN Q13

local real KQ13_function_AA(real ki, real rr, real xv)
{
    real KQ13;
    
    KQ13 = rsqr(rr);
    return KQ13;
}

local real GaussLegendreQ13_func_AA(real y)
{
    int j;
    real *xGL, *wGL;
    real kmin, kmax, ki;
    real Q13, KQ13;
    real PSLA, PSLB;
    real rmin, rmax;
    real kk, rr;
    real xmin, xmax;
    real xv, w, psl;
    int Nx;
    
    gd.p = rpow(10.0,y);
    
    kk = gd.p;
    ki = gd.k;
    rr = kk/ki;
    PSLB = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
    Q13 = 0.0;
    kmin = kPos(PSLT+1);
    kmax = kPos(PSLT+nPSLT-1);
    rmax = kmax/ki;
    rmin = kmin/ki;
    
    xmin = MAX(-1.0, (1.0 + rsqr(rr) - rsqr(rmax))/(2.0*rr));
    xmax = MIN( 1.0, (1.0 + rsqr(rr) - rsqr(rmin))/(2.0*rr));
    if (rr>=0.5)
        xmax = 0.5/rr;
    Nx=10;
    xGL=dvector(1,Nx);
    wGL=dvector(1,Nx);
    gauleg(xmin,xmax,xGL,wGL,Nx);
    for (j=1; j<=Nx; j++) {
        xv = xGL[j];
        w = wGL[j];
        KQ13 = KQ13_function(ki, rr, xv);
        psl = psInterpolation_nr(ki*rsqrt(abskmq), kPS, pPS, nPSLT);
        Q13 += w*KQ13*psl;
    }
    
    free_dvector(wGL,1,Nx);
    free_dvector(xGL,1,Nx);
    
    return gd.p*PSLB*Q13;
}

local real Q13_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    real deltay;  //Modificacion
    //~ real Nsteps; //Modificacion
    int j; //Modificacion
    real Qa; //Modificacion
    real Qb;     //Modificacion
    real yp;    //Modificacion
    //~ Nsteps=400.0;  //Modificacion
    
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    
    deltay=(ymax-ymin)/Nsteps;  //Modificacion
    
    result=0; //Modificacion
    Qa=GaussLegendreQ13_func( ymin ); //Modificacion
    
    for (j=1; j<=Nsteps; j++) {     //Modificacion
        yp = ymin + deltay * j;
        Qb=GaussLegendreQ13_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }//Modificacion
    
    result = 2.0*rlog(10.0)*(rsqr(ki)/FOURPI2)*result;  //Modificacion
    
    
    //~ result=2.0*rlog(10.0)*(rsqr(ki)/FOURPI2)
    //~ *qromo(GaussLegendreQ13_func,ymin,ymax,midpnt,cmd.epsquad);
    
    return result;
}

// END DE Q13

// BEGIN QI

local real KQI_function_AA(real ki, real rr, real xv)
{
    real k2, KA, KB, KQI;
    global_D2v2_ptr ptmp;
    
    k2 = ki * rsqrt(abskmq);
    ptmp = DsSecondOrder_func(ki, ki*rr, k2);
    if (model_int_flag==LCDM) {
        KQI = rsqr(rr) * (1.0 - rsqr(xv))/(1.0 + rsqr(rr - 2.0*rr*xv))
        * (
           KA_LCDM - KB_LCDM*(rsqr(xv) + rsqr(rr) - 2.0*rr*xv)/
           (1.0 + rsqr(rr) - 2.0*rr*xv)
           );
    } else {
        KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KQI = rsqr(rr) * (1.0 - rsqr(xv))/(1.0 + rsqr(rr - 2.0*rr*xv))
        * (
           KA - KB*(rsqr(xv) + rsqr(rr) - 2.0*rr*xv)/
           (1.0 + rsqr(rr) - 2.0*rr*xv)
           );
    }
    return KQI;
}

local real GaussLegendreQI_func_AA(real y)
{
    int j;
    real *xGL, *wGL;
    real kmin, kmax, ki;
    real QI, KQI;
    real PSLA, PSLB;
    real rmin, rmax;
    real kk, rr;
    real xmin, xmax;
    real xv, w, psl;
    int Nx;
    
    gd.p = rpow(10.0,y);
    
    kk = gd.p;
    ki = gd.k;
    rr = kk/ki;
    PSLB = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
    QI = 0.0;
    kmin = kPos(PSLT+1);
    kmax = kPos(PSLT+nPSLT-1);
    rmax = kmax/ki;
    rmin = kmin/ki;
    
    xmin = MAX(-1.0, (1.0 + rsqr(rr) - rsqr(rmax))/(2.0*rr));
    xmax = MIN( 1.0, (1.0 + rsqr(rr) - rsqr(rmin))/(2.0*rr));
    if (rr>=0.5)
        xmax = 0.5/rr;
    Nx=10;
    xGL=dvector(1,Nx);
    wGL=dvector(1,Nx);
    gauleg(xmin,xmax,xGL,wGL,Nx);
    for (j=1; j<=Nx; j++) {
        xv = xGL[j];
        w = wGL[j];
        KQI = KQI_function(ki, rr, xv);
        psl = psInterpolation_nr(ki*rsqrt(abskmq), kPS, pPS, nPSLT);
        QI += w*KQI*psl;
    }
    
    free_dvector(wGL,1,Nx);
    free_dvector(xGL,1,Nx);
    
    return gd.p*PSLB*QI;
}

local real QI_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    real deltay;  //Modificacion
    //~ real Nsteps; //Modificacion
    int j; //Modificacion
    real Qa; //Modificacion
    real Qb;     //Modificacion
    real yp;    //Modificacion
    //~ Nsteps=400.0;  //Modificacion
    
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    
    deltay=(ymax-ymin)/Nsteps;  //Modificacion
    
    result=0; //Modificacion
    Qa=GaussLegendreQI_func( ymin ); //Modificacion
    
    for (j=1; j<=Nsteps; j++) {     //Modificacion
        yp = ymin + deltay * j;
        Qb=GaussLegendreQI_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }//Modificacion
    
    
    result = 2.0*rlog(10.0)*(rsqr(ki)/FOURPI2)*result;  //Modificacion
    
    
    //~ result=2.0*rlog(10.0)*(rsqr(ki)/FOURPI2)
    //~ *qromo(GaussLegendreQI_func,ymin,ymax,midpnt,cmd.epsquad);
    
    return result;
}

// END DE QI


// BEGIN Q5

local real KQ5_function_AA(real ki, real rr, real xv)
{
    real k2, KA, KB, KQ5;
    global_D2v2_ptr ptmp;
    
    k2 = ki * rsqrt(abskmq);
    ptmp = DsSecondOrder_func(ki, ki*rr, k2);
    if (model_int_flag==LCDM) {
        KQ5 = rr*xv*(
                     KA_LCDM - KB_LCDM*rsqr(-rr+xv)/abskmq
                     );
    } else {
        KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KQ5 = rr*xv*(
                     KA - KB*rsqr(-rr+xv)/abskmq
                     );
    }
    return KQ5;
}

local real GaussLegendreQ5_func_AA(real y)
{
    int j;
    real ki;
    real Q5, KQ5;
    real PSLB;
    real kk, rr;
    real xv, w, psl;
    
    gd.p = rpow(10.0,y);
    
    kk = gd.p;
    ki = gd.k;
    rr = kk/ki;
    PSLB = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
    Q5 = 0.0;
    
    for (j=1; j<=nGL(pGL); j++) {
        xv = xGL(pGL)[j];
        w = wGL(pGL)[j];
        KQ5 = KQ5_function(ki, rr, xv);
        psl = psInterpolation_nr(ki*rsqrt(abskmq), kPS, pPS, nPSLT);
        Q5 += w*KQ5*psl;
    }
    
    return gd.p*PSLB*Q5;
}

local real Q5_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    real deltay;  //Modificacion
    //~ real Nsteps; //Modificacion
    int j; //Modificacion
    real Qa; //Modificacion
    real Qb;     //Modificacion
    real yp;    //Modificacion
    //~ Nsteps=400.0;  //Modificacion
    
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    
    deltay=(ymax-ymin)/Nsteps;  //Modificacion
    
    result=0; //Modificacion
    Qa=GaussLegendreQ5_func( ymin ); //Modificacion
    
    for (j=1; j<=Nsteps; j++) {     //Modificacion
        yp = ymin + deltay * j;
        Qb=GaussLegendreQ5_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }//Modificacion
    
    
    result = rlog(10.0)*(rsqr(ki)/FOURPI2)*result;  //Modificacion
    
    
    //~ result=rlog(10.0)*(rsqr(ki)/FOURPI2)
    //~ *qromo(GaussLegendreQ5_func,ymin,ymax,midpnt,cmd.epsquad);
    
    return result;
}

// END DE Q5

// BEGIN Q7

local real KQ7_function_AA(real ki, real rr, real xv)
{
    real KQ7;
    
    KQ7 = rsqr(xv)*(1-rr*xv)/abskmq;
    return KQ7;
}

local real GaussLegendreQ7_func_AA(real y)
{
    int j;
    real ki;
    real Q7, KQ7;
    real PSLB;
    real kk, rr;
    real xv, w, psl;
    
    gd.p = rpow(10.0,y);
    
    kk = gd.p;
    ki = gd.k;
    rr = kk/ki;
    PSLB = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
    Q7 = 0.0;
    
    for (j=1; j<=nGL(pGL); j++) {
        xv = xGL(pGL)[j];
        w = wGL(pGL)[j];
        KQ7 = KQ7_function(ki, rr, xv);
        psl = psInterpolation_nr(ki*rsqrt(abskmq), kPS, pPS, nPSLT);
        Q7 += w*KQ7*psl;
    }
    
    return gd.p*PSLB*Q7;
}

local real Q7_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    real deltay;  //Modificacion
    //~ real Nsteps; //Modificacion
    int j; //Modificacion
    real Qa; //Modificacion
    real Qb;     //Modificacion
    real yp;    //Modificacion
    //~ Nsteps=400.0;  //Modificacion
    
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    
    deltay=(ymax-ymin)/Nsteps;  //Modificacion
    
    result=0; //Modificacion
    Qa=GaussLegendreQ7_func( ymin ); //Modificacion
    
    for (j=1; j<=Nsteps; j++) {     //Modificacion
        yp = ymin + deltay * j;
        Qb=GaussLegendreQ7_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }//Modificacion
    
    //~ result = (rlog(10.0)/FOURPI2)*result;  //Modificacion
    result = rlog(10.0)*(rsqr(ki)/FOURPI2)*result;  //Modificacion
    
    
    //~ result=rlog(10.0)*(rsqr(ki)/FOURPI2)
    //~ *qromo(GaussLegendreQ7_func,ymin,ymax,midpnt,cmd.epsquad);
    
    return result;
}

// END DE Q7

// BEGIN Q11

local real KQ11_function_AA(real ki, real rr, real xv)
{
    real KQ11;
    
    KQ11 = rsqr(xv);
    return KQ11;
}

local real GaussLegendreQ11_func_AA(real y)
{
    int j;
    real ki;
    real Q11, KQ11;
    real PSLB;
    real kk, rr;
    real xv, w, psl;
    
    gd.p = rpow(10.0,y);
    
    kk = gd.p;
    ki = gd.k;
    rr = kk/ki;
    PSLB = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
    Q11 = 0.0;
    
    for (j=1; j<=nGL(pGL); j++) {
        xv = xGL(pGL)[j];
        w = wGL(pGL)[j];
        KQ11 = KQ11_function(ki, rr, xv);
        psl = psInterpolation_nr(ki*rsqrt(abskmq), kPS, pPS, nPSLT);
        Q11 += w*KQ11*psl;
    }
    
    return gd.p*PSLB*Q11;
}

local real Q11_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    real deltay;  //Modificacion
    //~ real Nsteps; //Modificacion
    int j; //Modificacion
    real Qa; //Modificacion
    real Qb;     //Modificacion
    real yp;    //Modificacion
    //~ Nsteps=400.0;  //Modificacion
    
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    
    deltay=(ymax-ymin)/Nsteps;  //Modificacion
    
    result=0; //Modificacion
    Qa=GaussLegendreQ11_func( ymin ); //Modificacion
    
    for (j=1; j<=Nsteps; j++) {     //Modificacion
        yp = ymin + deltay * j;
        Qb=GaussLegendreQ11_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }//Modificacion
    
    //~ result = (rlog(10.0)/FOURPI2)*result;  //Modificacion
    result = rlog(10.0)*(rsqr(ki)/FOURPI2)*result;  //Modificacion
    
    
    //~ result=rlog(10.0)*(rsqr(ki)/FOURPI2)
    //~ *qromo(GaussLegendreQ11_func,ymin,ymax,midpnt,cmd.epsquad);
    
    return result;
}

// END DE Q11

// BEGIN Q12

local real KQ12_function_AA(real ki, real rr, real xv)
{
    real KQ12;
    
    KQ12 = rr*xv;
    return KQ12;
}

local real GaussLegendreQ12_func_AA(real y)
{
    int j;
    real ki;
    real Q12, KQ12;
    real PSLB;
    real kk, rr;
    real xv, w, psl;
    
    gd.p = rpow(10.0,y);
    
    kk = gd.p;
    ki = gd.k;
    rr = kk/ki;
    PSLB = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
    Q12 = 0.0;
    
    for (j=1; j<=nGL(pGL); j++) {
        xv = xGL(pGL)[j];
        w = wGL(pGL)[j];
        KQ12 = KQ12_function(ki, rr, xv);
        psl = psInterpolation_nr(ki*rsqrt(abskmq), kPS, pPS, nPSLT);
        Q12 += w*KQ12*psl;
    }
    
    return gd.p*PSLB*Q12;
}

local real Q12_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    real deltay;  //Modificacion
    //~ real Nsteps; //Modificacion
    int j; //Modificacion
    real Qa; //Modificacion
    real Qb;     //Modificacion
    real yp;    //Modificacion
    //~ Nsteps=400.0;  //Modificacion
    
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    
    deltay=(ymax-ymin)/Nsteps;  //Modificacion
    
    result=0; //Modificacion
    Qa=GaussLegendreQ12_func( ymin ); //Modificacion
    
    for (j=1; j<=Nsteps; j++) {     //Modificacion
        yp = ymin + deltay * j;
        Qb=GaussLegendreQ12_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }//Modificacion
    
    //~ result = (rlog(10.0)/FOURPI2)*result;  //Modificacion
    result = rlog(10.0)*(rsqr(ki)/FOURPI2)*result;  //Modificacion
    
    
    //~ result=rlog(10.0)*(rsqr(ki)/FOURPI2)
    //~ *qromo(GaussLegendreQ12_func,ymin,ymax,midpnt,cmd.epsquad);
    
    return result;
}

// END DE Q12


// BEGIN R1

global real funcR1int_AA(real y)
{
    real ftmp,fac;
    global_D3v2_ptr ptmp;
    
    gd.p = rpow(10.0,y);
    fac = rpow(gd.p,3.0)*psInterpolation_nr(gd.p, kPS, pPS, nPSLT);
    
    ptmp = DsThirdOrder_func(gd.x, gd.k, gd.p);
    ftmp = (21.0/10.0)*D3symmD3v2(ptmp)/( DpkD3v2(ptmp)*DppD3v2(ptmp)*DppD3v2(ptmp) );
    
    return fac*ftmp;
}

local real R1_function_AA(real eta, real ki)
{
    real result, fac, s;
    real pmin, pmax, ymin, ymax;
    int j;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    //~ if (model_int_flag==LCDM) {
    //~ fac=0;
    //~ s=0;
    //~ }  else  {
    
    fac=(rlog(10.0)/FOURPI2)*psInterpolation_nr(gd.k, kPS, pPS, nPSLT);
    
    s=0;
    for (j=1;j<=nGL(pGL)/2;j++) {
        gd.x = xGL(pGL)[j];
        result = qromo(funcR1int,ymin,ymax,midpnt,cmd.epsquad);
        s += 2.0*wGL(pGL)[j]*result;
    }
    //~ }
    
    return fac*s;
}

// END DE R1


// BEGIN R2

local real GaussLegendreR2_func_AA(real y)
{
    real ss, fac;
    real R2aB, rr, xv, w, k2, KA, KB, KR2;
    global_D2v2_ptr ptmp;
    int j;
    
    gd.p = rpow(10.0,y);
    rr = gd.p/gd.k;
    fac = gd.p*psInterpolation_nr(gd.p, kPS, pPS, nPSLT);
    
    R2aB = 0.0;
    for (j=1; j<=nGL(pGL); j++) {
        xv = xGL(pGL)[j];
        w = wGL(pGL)[j];
        k2 = gd.k * rsqrt(1.0 + rsqr(rr) - 2.0*rr*xv);
        ptmp = DsSecondOrder_func(k2, gd.k, gd.p);
        if (model_int_flag==LCDM) {
            KR2 = (rr*xv*(1.0-rr*xv)/abskmq)
            *(
              KA_LCDM - KB_LCDM*rsqr(xv)
              );
        } else {
            KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
            KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
            KR2 = (rr*xv*(1.0-rr*xv)/abskmq)
            *(
              KA - KB*rsqr(xv)
              );
        }
        R2aB += w*KR2;
    }
    
    return fac*R2aB;
}

local real R2_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    
    
    real deltay;  //Modificacion
    //~ real Nsteps; //Modificacion
    int j; //Modificacion
    real Qa; //Modificacion
    real Qb;     //Modificacion
    real yp;    //Modificacion
    //~ Nsteps=400.0;  //Modificacion
    
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    
    deltay=(ymax-ymin)/Nsteps;  //Modificacion
    
    result=0; //Modificacion
    Qa=GaussLegendreR2_func( ymin ); //Modificacion
    
    for (j=1; j<=Nsteps; j++) {     //Modificacion
        yp = ymin + deltay * j;
        Qb=GaussLegendreR2_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }//Modificacion
    
    
    result = rlog(10.0)*(rsqr(gd.k)/FOURPI2)*psInterpolation_nr(gd.k, kPS, pPS, nPSLT)*result;  //Modificacion
    
    
    //~ result=rlog(10.0)*(rsqr(gd.k)/FOURPI2)*psInterpolation_nr(gd.k, kPS, pPS, nPSLT)
    //~ *qromo(GaussLegendreR2_func,ymin,ymax,midpnt,cmd.epsquad);
    
    return result;
}

// END DE R2

// BEGIN RI

local real KRI_function_AA(real ki, real rr, real xv)
{
    real k2, KA, KB, KRI;
    global_D2v2_ptr ptmp;
    
    k2 = ki * rsqrt(abskmq);
    ptmp = DsSecondOrder_func(k2,ki,ki*rr);
    if (model_int_flag==LCDM) {
        KRI = (rsqr(rr)*(1.0-rsqr(xv))/abskmq)
        *(
          KA_LCDM - KB_LCDM*rsqr(xv)
          );
    } else {
        KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KRI = (rsqr(rr)*(1.0-rsqr(xv))/abskmq)
        *(
          KA - KB*rsqr(xv)
          );
    }
    return KRI;
}

local real GaussLegendreRI_func_AA(real y)
{
    int j;
    real ki;
    real RI, KRI;
    real PSL;
    real kk, rr;
    real xv, w;
    
    gd.p = rpow(10.0,y);
    
    kk = gd.p;
    ki = gd.k;
    rr = kk/ki;
    PSL = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
    RI = 0.0;
    
    for (j=1; j<=nGL(pGL); j++) {
        xv = xGL(pGL)[j];
        w = wGL(pGL)[j];
        KRI = KRI_function(ki, rr, xv);
        RI += w*KRI;
    }
    
    return gd.p*RI*PSL;
}

local real RI_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    real PSL;
    
    
    real deltay;  //Modificacion
    int j; //Modificacion
    real Qa; //Modificacion
    real Qb;     //Modificacion
    real yp;    //Modificacion
    
    
    
    gd.xstop = eta;
    gd.k = ki;
    PSL = psInterpolation_nr(ki, kPS, pPS, nPSLT);
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    
    deltay=(ymax-ymin)/Nsteps;  //Modificacion
    
    result=0; //Modificacion
    Qa=GaussLegendreRI_func( ymin ); //Modificacion
    
    for (j=1; j<=Nsteps; j++) {     //Modificacion
        yp = ymin + deltay * j;
        Qb=GaussLegendreRI_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }//Modificacion
    
    
    result = rlog(10.0)*(rsqr(ki)/FOURPI2)*PSL*result;  //Modificacion
    
    
    
    
    //~ result=rlog(10.0)*(rsqr(ki)/FOURPI2)*PSL
    //~ *qromo(GaussLegendreRI_func,ymin,ymax,midpnt,cmd.epsquad);
    
    return result;
}

// END DE RI


// BEGIN R1p2

local real KR1p2_function_AA(real ki, real rr, real xv)
{
    real k2, KA, KB, KR1p2;
    global_D2v2_ptr ptmp;
    
    k2 = ki * rsqrt(abskmq);
    ptmp = DsSecondOrder_func(k2,ki,ki*rr);
    if (model_int_flag==LCDM) {
        KR1p2 = (rsqr(rr)*(1.0-rr*xv)/abskmq)
        *(
          KA_LCDM - KB_LCDM*rsqr(xv)
          );
    } else {
        KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KR1p2 = (rsqr(rr)*(1.0-rr*xv)/abskmq)
        *(
          KA - KB*rsqr(xv)
          );
    }
    return KR1p2;
}

local real GaussLegendreR1p2_func_AA(real y)
{
    int j;
    real ki;
    real R1p2, KR1p2;
    real PSL;
    real kk, rr;
    real xv, w;
    
    gd.p = rpow(10.0,y);
    
    kk = gd.p;
    ki = gd.k;
    rr = kk/ki;
    PSL = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
    R1p2 = 0.0;
    
    for (j=1; j<=nGL(pGL); j++) {
        xv = xGL(pGL)[j];
        w = wGL(pGL)[j];
        KR1p2 = KR1p2_function(ki, rr, xv);
        R1p2 += w*KR1p2;
    }
    
    return gd.p*R1p2*PSL;
}

local real R1p2_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    real PSL;
    
    
    real deltay;  //Modificacion
    int j; //Modificacion
    real Qa; //Modificacion
    real Qb;     //Modificacion
    real yp;    //Modificacion
    
    
    
    gd.xstop = eta;
    gd.k = ki;
    PSL = psInterpolation_nr(ki, kPS, pPS, nPSLT);
    
    pmin = kPos(PSLCDMtab+10);
    pmax = kPos(PSLCDMtab+nPSTable-10);
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    
    deltay=(ymax-ymin)/Nsteps;  //Modificacion
    
    result=0; //Modificacion
    Qa=GaussLegendreR1p2_func( ymin ); //Modificacion
    
    for (j=1; j<=Nsteps; j++) {     //Modificacion
        yp = ymin + deltay * j;
        Qb=GaussLegendreR1p2_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }//Modificacion
    
    
    result = rlog(10.0)*(rsqr(ki)/FOURPI2)*PSL*result;  //Modificacion
    
    //~ result=rlog(10.0)*(rsqr(ki)/FOURPI2)*PSL
    //~ *qromo(GaussLegendreR1p2_func,ymin,ymax,midpnt,cmd.epsquad);
    
    return result;
}

// END DE R1p2



// BEGIN Qs and Rs

global_QRs QsRs_functions_trapezoidAA(real eta, real ki)
{
    real Q1p, Q2p, Q3p, Q8p, Q9p, Q13p, QIp;
    real Q5p, Q7p, Q11p, Q12p;
    real R1p, R2p;
    real RI, R1p2;
    global_QRs_ptr QRstmp;
    
    fprintf(gd.outlog,"\nQuadrature method: %s",gd.quadraturemethod_comment);
    QRstmp = (global_QRs_ptr) allocate(1 * sizeof(global_QRs));
    
    Q1p = Q1_function_AA(eta, ki);
    Q2p = Q2_function_AA(eta, ki);
    Q3p = Q3_function_AA(eta, ki);
    Q8p = Q8_function_AA(eta, ki);
    Q9p = Q9_function_AA(eta, ki);
    Q13p = Q13_function_AA(eta, ki);
    QIp = QI_function_AA(eta, ki);
    Q5p = Q5_function_AA(eta, ki);
    Q7p = Q7_function_AA(eta, ki);
    Q11p = Q11_function_AA(eta, ki);
    Q12p = Q12_function_AA(eta, ki);
    RI = RI_function_AA(eta, ki);
    R1p2 = R1p2_function_AA(eta, ki);
    
    R1p = R1_function_AA(eta, ki);
    R2p = R2_function_AA(eta, ki);
    
    etaQRs(QRstmp) = eta;
    kQRs(QRstmp)    = ki;
    Q1(QRstmp)      = Q1p;
    Q2(QRstmp)      = Q2p;
    Q3(QRstmp)      = Q3p;
    Q8(QRstmp)      = Q8p;
    Q9(QRstmp)      = Q9p;
    Q13(QRstmp)      = Q13p;
    QI(QRstmp)      = QIp;
    Q5(QRstmp)      = Q5p;
    Q7(QRstmp)      = Q7p;
    Q11(QRstmp)      = Q11p;
    Q12(QRstmp)      = Q12p;
    RI(QRstmp)      = RI;
    R1p2(QRstmp)      = R1p2;
    R1(QRstmp)      = R1p;
    R2(QRstmp)      = R2p;
    
    return *QRstmp;
}

/*
 global_QRs QsRs_functions_trapezoidAA_LCDM(real eta, real ki)
 {
 real Q1p, Q2p, Q3p, Q8p, Q9p, Q13p, QIp;
 real Q5p, Q7p, Q11p, Q12p;
 real R1p, R2p;
 real RI, R1p2;
 global_QRs_ptr QRstmp;
 
 QRstmp = (global_QRs_ptr) allocate(1 * sizeof(global_QRs));
 
 Q1p = Q1_function_AA(eta, ki);
 Q2p = Q2_function_AA(eta, ki);
 Q3p = Q3_function_AA(eta, ki);
 Q8p = Q8_function_AA(eta, ki);
 Q9p = Q9_function_AA(eta, ki);
 Q13p = Q13_function_AA(eta, ki);
 QIp = QI_function_AA(eta, ki);
 Q5p = Q5_function_AA(eta, ki);
 Q7p = Q7_function_AA(eta, ki);
 Q11p = Q11_function_AA(eta, ki);
 Q12p = Q12_function_AA(eta, ki);
 RI = RI_function_AA(eta, ki);
 R1p2 = R1p2_function_AA(eta, ki);
 
 R1p = RI;
 R2p = R2_function_AA(eta, ki);
 
 etaQRs(QRstmp) = eta;
 kQRs(QRstmp)    = ki;
 Q1(QRstmp)      = Q1p;
 Q2(QRstmp)      = Q2p;
 Q3(QRstmp)      = Q3p;
 Q8(QRstmp)      = Q8p;
 Q9(QRstmp)      = Q9p;
 Q13(QRstmp)      = Q13p;
 QI(QRstmp)      = QIp;
 Q5(QRstmp)      = Q5p;
 Q7(QRstmp)      = Q7p;
 Q11(QRstmp)      = Q11p;
 Q12(QRstmp)      = Q12p;
 RI(QRstmp)      = RI;
 R1p2(QRstmp)      = R1p2;
 R1(QRstmp)      = R1p;
 R2(QRstmp)      = R2p;
 
 return *QRstmp;
 }
 */

// END DE Qs and Rs

/*
 global global_QRs QsRs_functions_driver(real eta, real ki)
 {
 global_QRs qrs;
 
 qrs = QsRs_functions(gd.xstop, ki);
 
 return qrs;
 }
 
 global global_QRs QsRs_functions_driver_LCDM(real eta, real ki)
 {
 global_QRs qrs;
 
 qrs = QsRs_functions_LCDM(gd.xstop, ki);
 
 return qrs;
 }
 */

#undef Nsteps
#undef abskmq


// ==============================================================================
