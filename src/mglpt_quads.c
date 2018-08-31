/*==============================================================================
 MODULE: mglpt_quads.c				[mglpt]
 ==============================================================================*/

#include "globaldefs.h"
#include "protodefs.h"

global_QRs QsRs_functions(real eta, real ki);

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
        KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KQ1 = rsqr(
                   KA + KB*(-1.0+(1.0-rsqr(xv))/abskmq)
                   );
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
            *qromo(GaussLegendreQ1_func,ymin,ymax,midpnt);

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
        KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KQ2 = (rr*xv*(1.0-rr*xv)/abskmq)
        *(
          KA - KB*((rsqr(xv)+rsqr(rr)-2.0*rr*xv)/abskmq)
          );
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
        *qromo(GaussLegendreQ2_func,ymin,ymax,midpnt);

    return result;
}

// END DE Q2


// BEGIN Q3

local real GaussLegendreQ3_func(real y)
{
    global_D2v2_ptr ptmp;
    int j;
    real *xGL, *wGL;
    real kmin, kmax, ki;
    real Q3p, Q3aA, Q3aB, KQ3;
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
        ptmp = DsSecondOrder_func(ki, ki*rr, k2);
        KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
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
    *qromo(GaussLegendreQ3_func,ymin,ymax,midpnt);
    
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
    KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
    KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
    KQ8 = rsqr(rr)*(
                    KA - KB*rsqr(-rr+xv)/abskmq
                    );
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
    *qromo(GaussLegendreQ8_func,ymin,ymax,midpnt);
    
    return result;
}

// END DE Q8


// BEGIN Q9

local real KQ9_function(real ki, real rr, real xv)
{
    real k2, KA, KB, KQ9;
    global_D2v2_ptr ptmp;
    
    k2 = ki * rsqrt(abskmq);
    ptmp = DsSecondOrder_func(ki, ki*rr, k2);
    KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
    KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
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
    *qromo(GaussLegendreQ9_func,ymin,ymax,midpnt);
    
    return result;
}

// END DE Q9

// BEGIN Q13

local real KQ13_function(real ki, real rr, real xv)
{
    real k2, KA, KB, KQ13;
    global_D2v2_ptr ptmp;
    
    k2 = ki * rsqrt(abskmq);
    ptmp = DsSecondOrder_func(ki, ki*rr, k2);
    KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
    KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
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
    *qromo(GaussLegendreQ13_func,ymin,ymax,midpnt);
    
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
    KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
    KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
    KQI = rsqr(rr) * (1.0 - rsqr(xv))/(1.0 + rsqr(rr - 2.0*rr*xv))
            * (
                KA - KB*(rsqr(xv) + rsqr(rr) - 2.0*rr*xv)/(1.0 + rsqr(rr) - 2.0*rr*xv)
            );
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
    *qromo(GaussLegendreQI_func,ymin,ymax,midpnt);
    
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
    KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
    KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
    KQ5 = rr*xv*(
                    KA - KB*rsqr(-rr+xv)/abskmq
                    );
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
    *qromo(GaussLegendreQ5_func,ymin,ymax,midpnt);
    
    return result;
}

// END DE Q5

// BEGIN Q7

local real KQ7_function(real ki, real rr, real xv)
{
    real k2, KA, KB, KQ7;
    global_D2v2_ptr ptmp;
    
    k2 = ki * rsqrt(abskmq);
    ptmp = DsSecondOrder_func(ki, ki*rr, k2);
    KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
    KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
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
    *qromo(GaussLegendreQ7_func,ymin,ymax,midpnt);
    
    return result;
}

// END DE Q7

// BEGIN Q11

local real KQ11_function(real ki, real rr, real xv)
{
    real k2, KA, KB, KQ11;
    global_D2v2_ptr ptmp;
    
    k2 = ki * rsqrt(abskmq);
    ptmp = DsSecondOrder_func(ki, ki*rr, k2);
    KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
    KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
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
    *qromo(GaussLegendreQ11_func,ymin,ymax,midpnt);
    
    return result;
}

// END DE Q11

// BEGIN Q12

local real KQ12_function(real ki, real rr, real xv)
{
    real k2, KA, KB, KQ12;
    global_D2v2_ptr ptmp;
    
    k2 = ki * rsqrt(abskmq);
    ptmp = DsSecondOrder_func(ki, ki*rr, k2);
    KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
    KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
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
    *qromo(GaussLegendreQ12_func,ymin,ymax,midpnt);
    
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
        result = qromo(funcR1int,ymin,ymax,midpnt);
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
        KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KR2 = (rr*xv*(1.0-rr*xv)/abskmq)
        *(
          KA - KB*rsqr(xv)
          );
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
            *qromo(GaussLegendreR2_func,ymin,ymax,midpnt);

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
    KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
    KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
    KRI = (rsqr(rr)*(1.0-rsqr(xv))/abskmq)
                *(
                    KA - KB*rsqr(xv)
                  );
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
        *qromo(GaussLegendreRI_func,ymin,ymax,midpnt);

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
    KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
    KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
    KR1p2 = (rsqr(rr)*(1.0-rr*xv)/abskmq)
                *(
                  KA - KB*rsqr(xv)
                 );
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
        *qromo(GaussLegendreR1p2_func,ymin,ymax,midpnt);

    return result;
}

// END DE R1p2



// BEGIN Qs and Rs

global_QRs QsRs_functions(real eta, real ki)
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

// END DE Qs and Rs

global global_QRs QsRs_functions_driver(real eta, real ki)
{
    global_QRs qrs;

    qrs = QsRs_functions(gd.xstop, ki);

    return qrs;
}

#undef abskmq
