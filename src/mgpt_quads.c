/*==============================================================================
 MODULE: mglpt_quads.c				[mgpt]
 ==============================================================================*/

#include "globaldefs.h"
#include "protodefs.h"
#include "models.h"

#define QROMBERG     qromo
#define KK  5

global_QRs QsRs_functions_trapezoid(real eta, real ki);
global_QRs QsRs_functions_trapezoid_LCDM(real eta, real ki);
global_QRs QsRs_functions_trapezoid3(real eta, real ki);
global_QRs QsRs_functions_trapezoid3_LCDM(real eta, real ki);
global_QRs QsRs_functions_romo(real eta, real ki);
global_QRs QsRs_functions_romo_LCDM(real eta, real ki);

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
local void quadrature_LCDM(real ki);

#define abskmq      (1.0+rsqr(rr)-2.0*rr*xv)
#define abskmqm     (1.0+rsqr(rr)+2.0*rr*xv)  //Modificacion

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
    kmin = kPS[2];
    kmax = kPS[nPSLT];
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
        if (model_int_flag==LCDM) {
            KQ1 = rsqr(
                       KA_LCDM + KB_LCDM*(-1.0+(1.0-rsqr(xv))/abskmq)
                       );
        } else {
            ptmp = DsSecondOrder_func(ki, ki*rr, k2);
            KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2(ptmp)*Dpk2D2(ptmp) );
            KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2(ptmp)*Dpk2D2(ptmp) );
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
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    result=(rlog(10.0)/FOURPI2)
    *QROMBERG(GaussLegendreQ1_func,ymin,ymax,midpnt,cmd.epsquad,KK);

    return result;
}

// END Q1


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
    kmin = kPS[2];
    kmax = kPS[nPSLT];
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
        if (model_int_flag==LCDM) {
            KQ2 = (rr*xv*(1.0-rr*xv)/abskmq)
            *(
              KA_LCDM - KB_LCDM*((rsqr(xv)+rsqr(rr)-2.0*rr*xv)/abskmq)
              );
        } else {
            ptmp = DsSecondOrder_func(ki, ki*rr, k2);
            KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2(ptmp)*Dpk2D2(ptmp) );
            KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2(ptmp)*Dpk2D2(ptmp) );
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
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    result=(rlog(10.0)/FOURPI2)
    *QROMBERG(GaussLegendreQ2_func,ymin,ymax,midpnt,cmd.epsquad,KK);

    return result;
}

// END Q2


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
    kmin = kPS[2];
    kmax = kPS[nPSLT];
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
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    result=(rlog(10.0)/FOURPI2)
    *QROMBERG(GaussLegendreQ3_func,ymin,ymax,midpnt,cmd.epsquad,KK);
    
    return result;
}

// END Q3


// BEGIN Q8

local real KQ8_function(real ki, real rr, real xv)
{
    real k2, KA, KB, KQ8;
    global_D2v2_ptr ptmp;

    k2 = ki * rsqrt(abskmq);
    if (model_int_flag==LCDM) {
        KQ8 = rsqr(rr)*(
                        KA_LCDM - KB_LCDM*rsqr(-rr+xv)/abskmq
                        );
    } else {
        ptmp = DsSecondOrder_func(ki, ki*rr, k2);
        KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2(ptmp)*Dpk2D2(ptmp) );
        KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2(ptmp)*Dpk2D2(ptmp) );
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
    kmin = kPS[2];
    kmax = kPS[nPSLT];
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
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    result=2.0*rlog(10.0)*(rsqr(ki)/FOURPI2)
    *QROMBERG(GaussLegendreQ8_func,ymin,ymax,midpnt,cmd.epsquad,KK);
    
    return result;
}

// END Q8


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
    kmin = kPS[2];
    kmax = kPS[nPSLT];
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
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    result=2.0*rlog(10.0)*(rsqr(ki)/FOURPI2)
    *QROMBERG(GaussLegendreQ9_func,ymin,ymax,midpnt,cmd.epsquad,KK);
    
    return result;
}

// END Q9

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
    kmin = kPS[2];
    kmax = kPS[nPSLT];
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
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    result=2.0*rlog(10.0)*(rsqr(ki)/FOURPI2)
    *QROMBERG(GaussLegendreQ13_func,ymin,ymax,midpnt,cmd.epsquad,KK);
    
    return result;
}

// END Q13

// BEGIN QI

local real KQI_function(real ki, real rr, real xv)
{
    real k2, KA, KB, KQI;
    global_D2v2_ptr ptmp;
    
    k2 = ki * rsqrt(abskmq);
    if (model_int_flag==LCDM) {
        KQI = rsqr(rr) * (1.0 - rsqr(xv))/(1.0 + rsqr(rr) - 2.0*rr*xv)
            * (
               KA_LCDM - KB_LCDM*(rsqr(xv) + rsqr(rr) - 2.0*rr*xv)/
                (1.0 + rsqr(rr) - 2.0*rr*xv)
               );
    } else {
        ptmp = DsSecondOrder_func(ki, ki*rr, k2);
        KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2(ptmp)*Dpk2D2(ptmp) );
        KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2(ptmp)*Dpk2D2(ptmp) );
        KQI = rsqr(rr) * (1.0 - rsqr(xv))/(1.0 + rsqr(rr) - 2.0*rr*xv)
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
    kmin = kPS[2];
    kmax = kPS[nPSLT];
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
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    result=2.0*rlog(10.0)*(rsqr(ki)/FOURPI2)
    *QROMBERG(GaussLegendreQI_func,ymin,ymax,midpnt,cmd.epsquad,KK);
    
    return result;
}

// END QI


// BEGIN Q5

local real KQ5_function(real ki, real rr, real xv)
{
    real k2, KA, KB, KQ5;
    global_D2v2_ptr ptmp;
    
    k2 = ki * rsqrt(abskmq);
    if (model_int_flag==LCDM) {
        KQ5 = rr*xv*(
                     KA_LCDM - KB_LCDM*rsqr(-rr+xv)/abskmq
                     );
    } else {
        ptmp = DsSecondOrder_func(ki, ki*rr, k2);
        KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2(ptmp)*Dpk2D2(ptmp) );
        KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2(ptmp)*Dpk2D2(ptmp) );
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
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);

    result=rlog(10.0)*(rsqr(ki)/FOURPI2)
    *QROMBERG(GaussLegendreQ5_func,ymin,ymax,midpnt,cmd.epsquad,KK);
    
    return result;
}

// END Q5

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
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    result=rlog(10.0)*(rsqr(ki)/FOURPI2)
    *QROMBERG(GaussLegendreQ7_func,ymin,ymax,midpnt,cmd.epsquad,KK);
    
    return result;
}

// END Q7

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
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    result=rlog(10.0)*(rsqr(ki)/FOURPI2)
    *QROMBERG(GaussLegendreQ11_func,ymin,ymax,midpnt,cmd.epsquad,KK);
    
    return result;
}

// END Q11

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
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    result=rlog(10.0)*(rsqr(ki)/FOURPI2)
    *QROMBERG(GaussLegendreQ12_func,ymin,ymax,midpnt,cmd.epsquad,KK);
    
    return result;
}

// END Q12


// BEGIN R1

global real funcR1int(real y)
{
    real ftmp,fac;
    global_D3v2_ptr ptmp;
    
    gd.p = rpow(10.0,y);
    fac = rpow(gd.p,3.0)*psInterpolation_nr(gd.p, kPS, pPS, nPSLT);

    if (model_int_flag==LCDM) {
        ftmp = KR1_LCDM;
    } else {
        ptmp = DsThirdOrder_func(gd.x, gd.k, gd.p);
        ftmp = (21.0/10.0)*D3symmD3(ptmp)/( DpkD3(ptmp)*DppD3(ptmp)*DppD3(ptmp) );
    }

    return fac*ftmp;
}

local real R1_function(real eta, real ki)
{
    real result, fac, s;
    real pmin, pmax, ymin, ymax;
    int j;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);

    fac=(rlog(10.0)/FOURPI2)*psInterpolation_nr(gd.k, kPS, pPS, nPSLT);
    
    s=0;
    for (j=1;j<=nGL(pGL)/2;j++) {
        gd.x = xGL(pGL)[j];
        result =
        QROMBERG(funcR1int,ymin,ymax,midpnt,cmd.epsquad,KK);
        s += 2.0*wGL(pGL)[j]*result;
    }

    return fac*s;
}

// END R1


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
        if (model_int_flag==LCDM) {
            KR2 = (rr*xv*(1.0-rr*xv)/abskmq)
                *(
                  KA_LCDM - KB_LCDM*rsqr(xv)
                  );
        } else {
            ptmp = DsSecondOrder_func(k2, gd.k, gd.p);
            KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2(ptmp)*Dpk2D2(ptmp) );
            KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2(ptmp)*Dpk2D2(ptmp) );
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
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);

    result=rlog(10.0)*(rsqr(gd.k)/FOURPI2)*psInterpolation_nr(gd.k, kPS, pPS, nPSLT)
        *QROMBERG(GaussLegendreR2_func,ymin,ymax,midpnt,cmd.epsquad,KK);

    return result;
}

// END R2

// BEGIN RI

local real KRI_function(real ki, real rr, real xv)
{
    real k2, KA, KB, KRI;
    global_D2v2_ptr ptmp;
    
    k2 = ki * rsqrt(abskmq);
    if (model_int_flag==LCDM) {
        KRI = (rsqr(rr)*(1.0-rsqr(xv))/abskmq)
            *(
              KA_LCDM - KB_LCDM*rsqr(xv)
              );
    } else {
        ptmp = DsSecondOrder_func(k2,ki,ki*rr);
        KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2(ptmp)*Dpk2D2(ptmp) );
        KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2(ptmp)*Dpk2D2(ptmp) );
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
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    result=rlog(10.0)*(rsqr(ki)/FOURPI2)*PSL
        *QROMBERG(GaussLegendreRI_func,ymin,ymax,midpnt,cmd.epsquad,KK);

    return result;
}

// END RI


// BEGIN R1p2

local real KR1p2_function(real ki, real rr, real xv)
{
    real k2, KA, KB, KR1p2;
    global_D2v2_ptr ptmp;
    
    k2 = ki * rsqrt(abskmq);
    if (model_int_flag==LCDM) {
        KR1p2 = (rsqr(rr)*(1.0-rr*xv)/abskmq)
                *(
                  KA_LCDM - KB_LCDM*rsqr(xv)
                  );
    } else {
        ptmp = DsSecondOrder_func(k2,ki,ki*rr);
        KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2(ptmp)*Dpk2D2(ptmp) );
        KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2(ptmp)*Dpk2D2(ptmp) );
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

    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);

    result=rlog(10.0)*(rsqr(ki)/FOURPI2)*PSL
        *QROMBERG(GaussLegendreR1p2_func,ymin,ymax,midpnt,cmd.epsquad,KK);

    return result;
}

// END R1p2



// BEGIN Qs and Rs

global_QRs QsRs_functions_trapezoid3(real eta, real ki)
{
    int i, j;
    real KR1, fac;
    global_D2v2_ptr ptmp;
    global_D3v2_ptr ptmpR1;
    real Dpkmin, Dpk;
    //
    real *xxGL, *wwGL;
    real kmin, kmax;
    
    real Q1p, Q2p, Q3p;
    real Q1aA, Q2aA, Q3aA;
    real Q1aB, Q2aB, Q3aB;
    
    real PSLA, PSLB;
    real rmin, rmax;
    real rr, deltar;
    real mumin, mumax;
    real xv, w, k2, psl;
    real psl1;
    real KA, KB, KQ1, KQ2, KQ3;
    real KAB;
    real KABm; //Modificacion
    int Nx;
    //
    real ypi, dk;
    //
    real Q8p, Q8aA, Q8aB, KQ8;
    real Q9p, Q9aA, Q9aB, KQ9;
    real Q13p, Q13aA, Q13aB, KQ13;
    
    real QIp, QIaA, QIaB, KQI;
    
    real Q5p, Q5aA, Q5aB, KQ5;
    real Q7p, Q7aA, Q7aB, KQ7;
    real Q11p, Q11aA, Q11aB, KQ11;
    real Q12p, Q12aA, Q12aB, KQ12;
    
    real R2p, R2aA, R2aB, KR2;
    real RIp, RIaA, RIaB, KRI;
    real R1p2p, R1p2aA, R1p2aB, KR1p2;
    
    real R1aA, R1aB, R1p;
    
    
    real *kk, *dkk;
    //
    pointPSTableptr p;
    //
    global_QRs_ptr QRstmp;
    
    QRstmp = (global_QRs_ptr) allocate(1 * sizeof(global_QRs));
    
    kmin = kPS[1];
    kmax = kPS[nPSLT];
    if (cmd.nquadSteps==1) {
        dk = 0.;
    } else {
        dk = (rlog10(kmax) - rlog10(kmin))/((real)(cmd.nquadSteps - 1));
    }
    
    kk=dvector(1,cmd.nquadSteps);
    dkk=dvector(1,cmd.nquadSteps);
    kk[1] = rpow(10.0,rlog10(kmin));
    for (i=2; i<cmd.nquadSteps; i++) {
        ypi = rlog10(kmin) + dk*((real)(i - 1));
        kk[i] = rpow(10.0,ypi);
        dkk[i] = (kk[i]-kk[i-1]);
    }
    
//
// LOOP FOR Q1, Q2, Q3, Q8, Q9, Q13 and QI and
// LOOP FOR Q5
    Q1p = 0.0; Q2p = 0.0; Q3p = 0.0;
    Q1aA = 0.0; Q2aA = 0.0; Q3aA = 0.0;
    Q1aB = 0.0; Q2aB = 0.0; Q3aB = 0.0;
    
    Q8p = 0.0; Q8aA = 0.0; Q8aB = 0.0;
    Q9p = 0.0; Q9aA = 0.0; Q9aB = 0.0;
    Q13p = 0.0; Q13aA = 0.0; Q13aB = 0.0;
    QIp = 0.0; QIaA = 0.0; QIaB = 0.0;
       
    //Modificacion
    Q5p = 0.0; Q5aA = 0.0; Q5aB = 0.0;
    //
    
    //
    PSLA = 0.0;
    rmax = kPS[nPSLT]/ki;
    rmin = kPS[1]/ki;
    Nx=10;
    xxGL=dvector(1,Nx);
    wwGL=dvector(1,Nx);
    for (i=2; i<cmd.nquadSteps; i++) {
        rr = kk[i]/ki;
        PSLB = psInterpolation_nr(kk[i], kPS, pPS, nPSLT);
        mumin = MAX(-1.0, (1.0 + rsqr(rr) - rsqr(rmax))/(2.0*rr));
        mumax = MIN( 1.0, (1.0 + rsqr(rr) - rsqr(rmin))/(2.0*rr));
        if (rr>=0.5)
            mumax = 0.5/rr;
        gauleg(mumin,mumax,xxGL,wwGL,Nx);
        for (j=1; j<=Nx; j++) {
            xv = xxGL[j];
            w = wwGL[j];
            k2 = ki * rsqrt(abskmq);
            ptmp = DsSecondOrder_func(ki, kk[i], k2);
            KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2(ptmp)*Dpk2D2(ptmp) );
            KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2(ptmp)*Dpk2D2(ptmp) );
            KQ1 = rsqr(rr)*rsqr(KA + KB*(-1.0+(1.0-rsqr(xv))/abskmq));
            KQ2 = (rr*xv*(1.0-rr*xv)/abskmq)
            *(KA - KB*((rsqr(xv)+rsqr(rr)-2.0*rr*xv)/abskmq));
            KQ3 = rsqr(xv)*rsqr(1.0-rr*xv)/rsqr(abskmq);
            KQI = (rsqr(rr)*(1.0 - rsqr(xv))/abskmq)
            *(KA - KB*( (rsqr(xv) + rsqr(rr) - 2.0*rr*xv)/abskmq));
            KQ8 = rsqr(rr)*(KA - KB*rsqr(-rr+xv)/abskmq);
            KQ9 = rr*xv*(1-rr*xv)/abskmq;
            KQ13 = rsqr(rr);
            
            
            //Modificacion
            KQ5 = (KA - KB*rsqr(-rr+xv)/abskmq)*0.5*(rr*xv + rr*rr*(1.-rr*xv)/abskmq );    
			//
            
//
            psl = psInterpolation_nr(k2, kPS, pPS, nPSLT);
            Q1aB += w*KQ1*psl;
            Q2aB += w*KQ2*psl;
            Q3aB += w*KQ3*psl;
            Q8aB += w*KQ8*psl;
            Q9aB += w*KQ9*psl;
            Q13aB += w*KQ13*psl;
            QIaB += w*KQI*psl;
            
            
            //Modificacion
            Q5aB += w*KQ5*psl;
			//
               
        }
        Q1p += dkk[i]*(Q1aA*PSLA + Q1aB*PSLB)/2.0;
        Q2p += dkk[i]*(Q2aA*PSLA + Q2aB*PSLB)/2.0;
        Q3p += dkk[i]*(Q3aA*PSLA + Q3aB*PSLB)/2.0;
        Q8p += dkk[i]*(Q8aA*PSLA + Q8aB*PSLB)/2.0;
        Q9p += dkk[i]*(Q9aA*PSLA + Q9aB*PSLB)/2.0;
        Q13p += dkk[i]*(Q13aA*PSLA + Q13aB*PSLB)/2.0;
        QIp += dkk[i]*(QIaA*PSLA + QIaB*PSLB)/2.0;
        Q1aA = Q1aB;
        Q2aA = Q2aB;
        Q3aA = Q3aB;
        Q8aA = Q8aB;
        Q9aA = Q9aB;
        Q13aA = Q13aB;
        QIaA = QIaB;
        Q1aB = 0.0;
        Q2aB = 0.0;
        Q3aB = 0.0;
        Q8aB = 0.0;
        Q9aB = 0.0;
        Q13aB = 0.0;
        QIaB = 0.0;
        
        
        //Modificacion
        Q5p     += dkk[i]*(Q5aA*PSLA + Q5aB*PSLB)/2.0;
        Q5aA = Q5aB; 
        Q5aB = 0.0; 
		//
		
        PSLA = PSLB;
    }
    Q1p *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
    Q2p *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
    Q3p *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
    Q8p *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
    Q9p *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
    Q13p *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
    QIp *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
    
    //Modificacion
    Q5p     *= 2.0*(rpow(ki,3.0)/FOURPI2)/ki;
    //~ Q7p     *= 2.0*(rpow(ki,3.0)/FOURPI2);
    //~ Q11p    *= 2.0*(rpow(ki,3.0)/FOURPI2);
    //~ Q12p    *= 2.0*(rpow(ki,3.0)/FOURPI2);
    //
        
    free_dvector(wwGL,1,Nx);
    free_dvector(xxGL,1,Nx);
//
//~ // LOOP FOR Q7, Q11 and Q12 and
    //~ Q5p = 0.0; Q5aA = 0.0; Q5aB = 0.0;
    Q7p = 0.0; Q7aA = 0.0; Q7aB = 0.0;
    Q11p = 0.0; Q11aA = 0.0; Q11aB = 0.0;
    Q12p = 0.0; Q12aA = 0.0; Q12aB = 0.0;
    PSLA = 0.0;
    for (i=2; i<cmd.nquadSteps; i++) {
        rr = kk[i]/ki;
        PSLB = psInterpolation_nr(kk[i], kPS, pPS, nPSLT);
        for (j=1; j<=nGL(pGL); j++) {
            xv = xGL(pGL)[j];
            w = wGL(pGL)[j];
            k2 = ki * rsqrt(abskmq);
            //~ ptmp = DsSecondOrder_func(ki, kk[i], k2);
            //~ KA = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2(ptmp)*Dpk2D2(ptmp) );
            //~ KB = DB2D2(ptmp)/( (3.0/7.0)*Dpk1D2(ptmp)*Dpk2D2(ptmp) );
            //~ KQ5 = rr*xv*(KA - KB*rsqr(-rr+xv)/abskmq);
            KQ7 = rsqr(xv)*(1-rr*xv)/abskmq;
            KQ11 = rsqr(xv);
            KQ12 = rr*xv;
            psl = psInterpolation_nr(k2, kPS, pPS, nPSLT);
            //~ Q5aB += w*KQ5*psl;
            Q7aB += w*KQ7*psl;
            Q11aB += w*KQ11*psl;
            Q12aB += w*KQ12*psl;
        }
        //~ Q5p     += dkk[i]*(Q5aA*PSLA + Q5aB*PSLB)/(2.0*ki);
        Q7p     += dkk[i]*(Q7aA*PSLA + Q7aB*PSLB)/(2.0*ki);
        Q11p    += dkk[i]*(Q11aA*PSLA + Q11aB*PSLB)/(2.0*ki);
        Q12p    += dkk[i]*(Q12aA*PSLA + Q12aB*PSLB)/(2.0*ki);
        //~ Q5aA = Q5aB; 
        Q7aA = Q7aB; Q11aA = Q11aB; Q12aA = Q12aB;
        PSLA = PSLB;
        //~ Q5aB = 0.0;
        Q7aB = 0.0; Q11aB = 0.0; Q12aB = 0.0;
    }
    //~ Q5p     *= (rpow(ki,3.0)/FOURPI2);
    Q7p     *= (rpow(ki,3.0)/FOURPI2);
    Q11p    *= (rpow(ki,3.0)/FOURPI2);
    Q12p    *= (rpow(ki,3.0)/FOURPI2);
//
// LOOP FOR R1
// and FOR R2, RI and R1p2
    R1p = 0.0; R1aA = 0.0; R1aB = 0.0;
    //Modificacion    
    R2p = 0.0; R2aA = 0.0; R2aB = 0.0;
    RIp = 0.0; RIaA = 0.0; RIaB = 0.0;
    R1p2p = 0.0; R1p2aA = 0.0; R1p2aB = 0.0;
    //
// Addition to local GL quad...
    Nx=10;
    xxGL=dvector(1,Nx);
    wwGL=dvector(1,Nx);
    gauleg(-1.,1.,xxGL,wwGL,Nx);
//
    for (i=2; i<cmd.nquadSteps; i++) {
        rr = kk[i]/ki;
//
//        for (j=1; j<=nGL(pGL)/2; j++) {
//            xv = xGL(pGL)[j];
//            w = wGL(pGL)[j];
// Instead these:
        for (j=1; j<=Nx/2; j++) {
        //~ for (j=1; j<=Nx; j++) {
            xv = xxGL[j];
            w = wwGL[j];
//
            //~ ptmpR1 = DsThirdOrder_func(xv, ki, kk[i]);
            ptmpR1 = DsThirdOrder_func(-xv, ki, kk[i]); //Modificacion
            KR1 = rsqr(rr)*(21.0/10.0)*D3symmD3(ptmpR1)
            /( DpkD3(ptmpR1)*DppD3(ptmpR1)*DppD3(ptmpR1) );
            
            //Modificacion
            KAB=  D2fD3(ptmpR1)/( (3.0/7.0)*DpkD3(ptmpR1)*DppD3(ptmpR1));
            KABm  = D2mfD3(ptmpR1)/( (3.0/7.0)*DpkD3(ptmpR1)*DppD3(ptmpR1));
            KR2   = ( rr*xv*(1.0-rr*xv)/abskmq )*KAB + ( -rr*xv*(1.0+rr*xv)/abskmqm )*KABm;
            KRI   = ( rsqr(rr)*(1.0-rsqr(xv))/abskmq )*KAB + ( rsqr(rr)*(1.0-rsqr(xv))/abskmqm )*KABm;
            KR1p2 = ( rsqr(rr)*(1.0-rr*xv)/abskmq )*KAB + ( rsqr(rr)*(1.0+rr*xv)/abskmqm )*KABm;
            //
            
            
            psl = psInterpolation_nr(kk[i], kPS, pPS, nPSLT);
            R1aB += w*KR1*psl;
            //Modificacion
            R2aB += w*KR2*psl;
            R1p2aB += w*KR1p2*psl;
            RIaB += w*KRI*psl;
            //
        }
        //~ R1p += 2.0*dkk[i]*( R1aA + R1aB )/(2.0*ki);
        R1p += 2.0*dkk[i]*( R1aA + R1aB )/(2.0*ki);
        R1aA = R1aB;
        R1aB = 0.0;
        //Modificacion
        R2p     += dkk[i]*(R2aA + R2aB)/(2.0*ki);
        RIp     += dkk[i]*(RIaA + RIaB)/(2.0*ki);
        R1p2p   += dkk[i]*(R1p2aA + R1p2aB)/(2.0*ki);
        R2aA = R2aB; RIaA = RIaB; R1p2aA = R1p2aB;
        R2aB = 0.0; RIaB = 0.0; R1p2aB = 0.0;   
        //     
    }
    psl1 = psInterpolation_nr(ki, kPS, pPS, nPSLT);
    R1p *= (rpow(ki,3.0)/FOURPI2)*psl1;
    //Modificacion
    R2p *= (rpow(ki,3.0)/FOURPI2)*psl1;
    RIp *= (rpow(ki,3.0)/FOURPI2)*psl1;
    R1p2p *= (rpow(ki,3.0)/FOURPI2)*psl1;
    //
// Addition to local GL quad...
    free_dvector(wwGL,1,Nx);
    free_dvector(xxGL,1,Nx);
//
//
    etaQRs(QRstmp) = eta;
    kQRs(QRstmp)    = ki;
    
    Q1(QRstmp)      = Q1p;
    Q2(QRstmp)      = Q2p;
    Q3(QRstmp)      = Q3p;
    Q8(QRstmp)      = Q8p;
    Q9(QRstmp)      = Q9p;
    Q13(QRstmp)     = Q13p;
    QI(QRstmp)      = QIp;
    
    Q5(QRstmp)      = Q5p;
    Q7(QRstmp)      = Q7p;
    Q11(QRstmp)     = Q11p;
    Q12(QRstmp)     = Q12p;
    
    R2(QRstmp)      = R2p;
    RI(QRstmp)      = RIp;
    R1p2(QRstmp)    = R1p2p;
    R1(QRstmp)      = R1p;
    
    free_dvector(dkk,1,cmd.nquadSteps);
    free_dvector(kk,1,cmd.nquadSteps);
    
    return *QRstmp;
}

global_QRs QsRs_functions_trapezoid3_LCDM(real eta, real ki)
{
    int i, j;
    real KR1, fac;
    global_D2v2_ptr ptmp;
//    global_D3v2_ptr ptmpR1;
    real Dpkmin, Dpk;
//
    real *xxGL, *wwGL;
    real kmin, kmax;
    
    real Q1p, Q2p, Q3p;
    real Q1aA, Q2aA, Q3aA;
    real Q1aB, Q2aB, Q3aB;
    
    real PSLA, PSLB;
    real rmin, rmax;
    real rr, deltar;
    real mumin, mumax;
    real xv, w, k2, psl;
    real psl1;
    real KA, KB, KQ1, KQ2, KQ3;
    real KAB;
    int Nx;
//
    real ypi, dk;
//
    real Q8p, Q8aA, Q8aB, KQ8;
    real Q9p, Q9aA, Q9aB, KQ9;
    real Q13p, Q13aA, Q13aB, KQ13;
    
    real QIp, QIaA, QIaB, KQI;
    
    real Q5p, Q5aA, Q5aB, KQ5;
    real Q7p, Q7aA, Q7aB, KQ7;
    real Q11p, Q11aA, Q11aB, KQ11;
    real Q12p, Q12aA, Q12aB, KQ12;
    
    real R2p, R2aA, R2aB, KR2;
    real RIp, RIaA, RIaB, KRI;
    real R1p2p, R1p2aA, R1p2aB, KR1p2;
    
    real R1aA, R1aB, R1p;
    
    real *kk, *dkk;
    //
    pointPSTableptr p;
    //
    global_QRs_ptr QRstmp;
    
    QRstmp = (global_QRs_ptr) allocate(1 * sizeof(global_QRs));
    
    kmin = kPS[1];
    kmax = kPS[nPSLT];
    if (cmd.nquadSteps==1) {
        dk = 0.;
    } else {
        dk = (rlog10(kmax) - rlog10(kmin))/((real)(cmd.nquadSteps - 1));
    }
    
    kk=dvector(1,cmd.nquadSteps);
    dkk=dvector(1,cmd.nquadSteps);
    kk[1] = rpow(10.0,rlog10(kmin));
    for (i=2; i<cmd.nquadSteps; i++) {
        ypi = rlog10(kmin) + dk*((real)(i - 1));
        kk[i] = rpow(10.0,ypi);
        dkk[i] = (kk[i]-kk[i-1]);
    }
//
// LOOP FOR Q1, Q2, Q3, Q8, Q9, Q13 and QI
    Q1p = 0.0; Q2p = 0.0; Q3p = 0.0;
    Q1aA = 0.0; Q2aA = 0.0; Q3aA = 0.0;
    Q1aB = 0.0; Q2aB = 0.0; Q3aB = 0.0;
    
    Q8p = 0.0; Q8aA = 0.0; Q8aB = 0.0;
    Q9p = 0.0; Q9aA = 0.0; Q9aB = 0.0;
    Q13p = 0.0; Q13aA = 0.0; Q13aB = 0.0;
    QIp = 0.0; QIaA = 0.0; QIaB = 0.0;
//
    PSLA = 0.0;
    rmax = kmax/ki;
    rmin = rmin/ki;
    Nx=10;
    xxGL=dvector(1,Nx);
    wwGL=dvector(1,Nx);
    for (i=2; i<cmd.nquadSteps; i++) {
        rr = kk[i]/ki;
        PSLB = psInterpolation_nr(kk[i], kPS, pPS, nPSLT);
        mumin = MAX(-1.0, (1.0 + rsqr(rr) - rsqr(rmax))/(2.0*rr));
        mumax = MIN( 1.0, (1.0 + rsqr(rr) - rsqr(rmin))/(2.0*rr));
        if (rr>=0.5)
            mumax = 0.5/rr;
        gauleg(mumin,mumax,xxGL,wwGL,Nx);
        for (j=1; j<=Nx; j++) {
            xv = xxGL[j];
            w = wwGL[j];
            k2 = ki * rsqrt(abskmq);
            KA = KA_LCDM;
            KB = KA;
            KQ1 = rsqr(rr)*rsqr(KA + KB*(-1.0+(1.0-rsqr(xv))/abskmq));
            KQ2 = (rr*xv*(1.0-rr*xv)/abskmq)
            *(KA - KB*((rsqr(xv)+rsqr(rr)-2.0*rr*xv)/abskmq));
            KQ3 = rsqr(xv)*rsqr(1.0-rr*xv)/rsqr(abskmq);
            KQI = (rsqr(rr)*(1.0 - rsqr(xv))/abskmq)
            *(KA - KB*( (rsqr(xv) + rsqr(rr) - 2.0*rr*xv)/abskmq));
            KQ8 = rsqr(rr)*(KA - KB*rsqr(-rr+xv)/abskmq);
            KQ9 = rr*xv*(1-rr*xv)/abskmq;
            KQ13 = rsqr(rr);
            //
            psl = psInterpolation_nr(k2, kPS, pPS, nPSLT);
            Q1aB += w*KQ1*psl;
            Q2aB += w*KQ2*psl;
            Q3aB += w*KQ3*psl;
            Q8aB += w*KQ8*psl;
            Q9aB += w*KQ9*psl;
            Q13aB += w*KQ13*psl;
            QIaB += w*KQI*psl;
        }
        Q1p += dkk[i]*(Q1aA*PSLA + Q1aB*PSLB)/2.0;
        Q2p += dkk[i]*(Q2aA*PSLA + Q2aB*PSLB)/2.0;
        Q3p += dkk[i]*(Q3aA*PSLA + Q3aB*PSLB)/2.0;
        Q8p += dkk[i]*(Q8aA*PSLA + Q8aB*PSLB)/2.0;
        Q9p += dkk[i]*(Q9aA*PSLA + Q9aB*PSLB)/2.0;
        Q13p += dkk[i]*(Q13aA*PSLA + Q13aB*PSLB)/2.0;
        QIp += dkk[i]*(QIaA*PSLA + QIaB*PSLB)/2.0;
        Q1aA = Q1aB;
        Q2aA = Q2aB;
        Q3aA = Q3aB;
        Q8aA = Q8aB;
        Q9aA = Q9aB;
        Q13aA = Q13aB;
        QIaA = QIaB;
        Q1aB = 0.0;
        Q2aB = 0.0;
        Q3aB = 0.0;
        Q8aB = 0.0;
        Q9aB = 0.0;
        Q13aB = 0.0;
        QIaB = 0.0;
        PSLA = PSLB;
    }
    Q1p *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
    Q2p *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
    Q3p *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
    Q8p *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
    Q9p *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
    Q13p *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
    QIp *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
    free_dvector(wwGL,1,Nx);
    free_dvector(xxGL,1,Nx);
//
// LOOP FOR Q5, Q7, Q11 and Q12 and
    Q5p = 0.0; Q5aA = 0.0; Q5aB = 0.0;
    Q7p = 0.0; Q7aA = 0.0; Q7aB = 0.0;
    Q11p = 0.0; Q11aA = 0.0; Q11aB = 0.0;
    Q12p = 0.0; Q12aA = 0.0; Q12aB = 0.0;
    PSLA = 0.0;
    for (i=2; i<cmd.nquadSteps; i++) {
        rr = kk[i]/ki;
        PSLB = psInterpolation_nr(kk[i], kPS, pPS, nPSLT);
        for (j=1; j<=nGL(pGL); j++) {
            xv = xGL(pGL)[j];
            w = wGL(pGL)[j];
            k2 = ki * rsqrt(abskmq);
            KA = KA_LCDM;
            KB = KA;
            KQ5 = rr*xv*(KA - KB*rsqr(-rr+xv)/abskmq);
            KQ7 = rsqr(xv)*(1-rr*xv)/abskmq;
            KQ11 = rsqr(xv);
            KQ12 = rr*xv;
            psl = psInterpolation_nr(k2, kPS, pPS, nPSLT);
            Q5aB += w*KQ5*psl;
            Q7aB += w*KQ7*psl;
            Q11aB += w*KQ11*psl;
            Q12aB += w*KQ12*psl;
        }
        Q5p     += dkk[i]*(Q5aA*PSLA + Q5aB*PSLB)/(2.0*ki);
        Q7p     += dkk[i]*(Q7aA*PSLA + Q7aB*PSLB)/(2.0*ki);
        Q11p    += dkk[i]*(Q11aA*PSLA + Q11aB*PSLB)/(2.0*ki);
        Q12p    += dkk[i]*(Q12aA*PSLA + Q12aB*PSLB)/(2.0*ki);
        Q5aA = Q5aB; Q7aA = Q7aB; Q11aA = Q11aB; Q12aA = Q12aB;
        PSLA = PSLB;
        Q5aB = 0.0; Q7aB = 0.0; Q11aB = 0.0; Q12aB = 0.0;
    }
    Q5p     *= (rpow(ki,3.0)/FOURPI2);
    Q7p     *= (rpow(ki,3.0)/FOURPI2);
    Q11p    *= (rpow(ki,3.0)/FOURPI2);
    Q12p    *= (rpow(ki,3.0)/FOURPI2);
//
// LOOP FOR R2, RI and R1p2
    R2p = 0.0; R2aA = 0.0; R2aB = 0.0;
    RIp = 0.0; RIaA = 0.0; RIaB = 0.0;
    R1p2p = 0.0; R1p2aA = 0.0; R1p2aB = 0.0;
    for (i=2; i<cmd.nquadSteps; i++) {
        rr = kk[i]/ki;
        for (j=1; j<=nGL(pGL); j++) {
            xv = xGL(pGL)[j];
            w = wGL(pGL)[j];
            k2 = ki * rsqrt(abskmq);
            KA = KA_LCDM;
            KB = KA;
            KAB = KA - KB*rsqr(xv);
            KR2 = ( rr*xv*(1.0-rr*xv)/abskmq )*KAB;
            KRI = ( rsqr(rr)*(1.0-rsqr(xv))/abskmq )*KAB;
            KR1p2 = ( rsqr(rr)*(1.0-rr*xv)/abskmq )*KAB;
            psl = psInterpolation_nr(kk[i], kPS, pPS, nPSLT);
            R2aB += w*KR2*psl;
            R1p2aB += w*KR1p2*psl;
            RIaB += w*KRI*psl;
        }
        R2p     += dkk[i]*(R2aA + R2aB)/(2.0*ki);
        RIp     += dkk[i]*(RIaA + RIaB)/(2.0*ki);
        R1p2p   += dkk[i]*(R1p2aA + R1p2aB)/(2.0*ki);
        R2aA = R2aB; RIaA = RIaB; R1p2aA = R1p2aB;
        R2aB = 0.0; RIaB = 0.0; R1p2aB = 0.0;
    }
    fac = psInterpolation_nr(ki, kPS, pPS, nPSLT);
    R2p *= (rpow(ki,3.0)/FOURPI2)*fac;
    RIp *= (rpow(ki,3.0)/FOURPI2)*fac;
    R1p2p *= (rpow(ki,3.0)/FOURPI2)*fac;
//
// LOOP FOR R1
    R1p = 0.0; R1aA = 0.0; R1aB = 0.0;
// Addition to local GL quad...
    Nx=10;
    xxGL=dvector(1,Nx);
    wwGL=dvector(1,Nx);
    gauleg(-1.,1.,xxGL,wwGL,Nx);
//
    for (i=2; i<cmd.nquadSteps; i++) {
        rr = kk[i]/ki;
//
//        for (j=1; j<=nGL(pGLRs)/2; j++) {
//            xv = xGL(pGLRs)[j];
//            w = wGL(pGLRs)[j];
// Instead these:
            for (j=1; j<=Nx/2; j++) {
                xv = xxGL[j];
                w = wwGL[j];
//
            KR1 = rsqr(rr)*KR1_LCDM;
            psl = psInterpolation_nr(kk[i], kPS, pPS, nPSLT);
            R1aB += w*KR1*psl;
        }
        R1p += 2.0*dkk[i]*( R1aA + R1aB )/(2.0*ki);
        R1aA = R1aB;
        R1aB = 0.0;
    }
    psl1 = psInterpolation_nr(ki, kPS, pPS, nPSLT);
    R1p *= (rpow(ki,3.0)/FOURPI2)*psl1;
//
//    R1p = RIp;
// Addition to local GL quad...
    free_dvector(wwGL,1,Nx);
    free_dvector(xxGL,1,Nx);
//
//
    etaQRs(QRstmp) = eta;
    kQRs(QRstmp)    = ki;
    
    Q1(QRstmp)      = Q1p;
    Q2(QRstmp)      = Q2p;
    Q3(QRstmp)      = Q3p;
    Q8(QRstmp)      = Q8p;
    Q9(QRstmp)      = Q9p;
    Q13(QRstmp)     = Q13p;
    QI(QRstmp)      = QIp;
    
    Q5(QRstmp)      = Q5p;
    Q7(QRstmp)      = Q7p;
    Q11(QRstmp)     = Q11p;
    Q12(QRstmp)     = Q12p;
    
    R2(QRstmp)      = R2p;
    RI(QRstmp)      = RIp;
    R1p2(QRstmp)    = R1p2p;
    R1(QRstmp)      = R1p;
    
    free_dvector(dkk,1,cmd.nquadSteps);
    free_dvector(kk,1,cmd.nquadSteps);
    
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


global_QRs QsRs_functions_romo_LCDM(real eta, real ki)
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


// END Qs and Rs

global_QRs qrs;

global global_QRs QsRs_functions_driver(real eta, real ki)
{
    quadrature(ki);
    return qrs;
}

global global_QRs QsRs_functions_driver_LCDM(real eta, real ki)
{
    quadrature_LCDM(ki);
    return qrs;
}

#define ROMO 1
#define NULLMETHOD 0
#define TRAPEZOID 2
#define TRAPEZOID3 5

local void quadrature(real ki)
{
    switch (gd.quadmethod_int) {
        case ROMO:
            qrs = QsRs_functions_romo(gd.xstop, ki);
            //~ qrs = QsRs_functions_trapezoid3(gd.xstop, ki);
            break;
//
        case TRAPEZOID:
            qrs = QsRs_functions_trapezoid(gd.xstop, ki);
            //~ qrs = QsRs_functions_trapezoid3(gd.xstop, ki);
            break;
//
        case TRAPEZOID3:
            qrs = QsRs_functions_trapezoid3(gd.xstop, ki);
            //~ qrs = QsRs_functions_trapezoid3(gd.xstop, ki);
            break;
//
        case NULLMETHOD:
            qrs = QsRs_functions_trapezoid(gd.xstop, ki);
            //~ qrs = QsRs_functions_trapezoid3(gd.xstop, ki);
            break;
//
        default:
            qrs = QsRs_functions_trapezoid(gd.xstop, ki);
            //~ qrs = QsRs_functions_trapezoid3(gd.xstop, ki);
            break;
    }
}

local void quadrature_LCDM(real ki)
{
    switch (gd.quadmethod_int) {
        case ROMO:
            qrs = QsRs_functions_romo_LCDM(gd.xstop, ki);
            //~ qrs = QsRs_functions_trapezoid3_LCDM(gd.xstop, ki);
            break;
//
        case TRAPEZOID:
            qrs = QsRs_functions_trapezoid_LCDM(gd.xstop, ki);
            //~ qrs = QsRs_functions_trapezoid3_LCDM(gd.xstop, ki);
            break;
//
        case TRAPEZOID3:
            qrs = QsRs_functions_trapezoid3_LCDM(gd.xstop, ki);
            break;
//
        case NULLMETHOD:
            qrs = QsRs_functions_trapezoid_LCDM(gd.xstop, ki);
            //~ qrs = QsRs_functions_trapezoid3_LCDM(gd.xstop, ki);
            break;
//
        default:
            qrs = QsRs_functions_trapezoid_LCDM(gd.xstop, ki);
            //~ qrs = QsRs_functions_trapezoid3_LCDM(gd.xstop, ki);
            break;
    }
}

void quadraturemethod_string_to_int(string method_str,int *method_int)
{
    *method_int=-1;
    if (strcmp(method_str,"romberg") == 0) {
        *method_int = ROMO;
        strcpy(gd.quadraturemethod_comment, "romberg open quadrature method");
    }
//
    if (strcmp(method_str,"trapezoid") == 0) {
        *method_int = TRAPEZOID;
        strcpy(gd.quadraturemethod_comment, "trapezoid quadrature method");
    }
//
    if (strcmp(method_str,"trapezoid3") == 0) {
        *method_int = TRAPEZOID3;
        strcpy(gd.quadraturemethod_comment, "trapezoid3 quadrature method");
    }
//
    if (strnull(method_str)) {
        *method_int = NULLMETHOD;
        strcpy(gd.quadraturemethod_comment,
               "given null quadrature method ... running deafult (trapezoid)");
        fprintf(stdout,"\n\tintegration: default integration method (trapezoid)...\n");
    }
//
    if (*method_int == -1) {
        *method_int = TRAPEZOID;
        strcpy(gd.quadraturemethod_comment,
               "Unknown quadrature method ... running deafult (trapezoid)");
        fprintf(stdout,"\n\tquadrature: Unknown method... %s ",cmd.quadratureMethod);
        fprintf(stdout,
                "\n\trunning default quadrature method (trapezoid)...\n");
    }
}

#undef ROMO
#undef TRAPEZOID
#undef TRAPEZOID3
#undef NULLMETHOD

// ==============================================================================
// AA TRAPEZOID METHOD:

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


// BEGIN Q1

local real Q1_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    real deltay;
    int j;
    real Qa;
    real Qb;
    real yp;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    deltay=(ymax-ymin)/cmd.nquadSteps;
    
    result=0;
    Qa=GaussLegendreQ1_func( ymin );
    
    for (j=1; j<=cmd.nquadSteps; j++) {
        yp = ymin + deltay * j;
        Qb=GaussLegendreQ1_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }
    
    result = rlog(10.0)*result/FOURPI2;

    return result;
}

// END Q1


// BEGIN Q2

local real Q2_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    real deltay;
    int j;
    real Qa;
    real Qb;
    real yp;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);

    deltay=(ymax-ymin)/cmd.nquadSteps;
    
    result=0;
    Qa=GaussLegendreQ2_func( ymin );
    
    for (j=1; j<=cmd.nquadSteps; j++) {
        yp = ymin + deltay * j;
        Qb=GaussLegendreQ2_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }
    
    
    result = rlog(10.0)*result/FOURPI2;

    return result;
}

// END Q2


// BEGIN Q3

local real Q3_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    real deltay;
    int j;
    real Qa;
    real Qb;
    real yp;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    deltay=(ymax-ymin)/cmd.nquadSteps;
    
    result=0;
    Qa=GaussLegendreQ3_func( ymin );
    
    for (j=1; j<=cmd.nquadSteps; j++) {
        yp = ymin + deltay * j;
        Qb=GaussLegendreQ3_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }
    
    result = rlog(10.0)*result/FOURPI2;

    return result;
}

// END Q3


// BEGIN Q8

local real Q8_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    real deltay;
    int j;
    real Qa;
    real Qb;
    real yp;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    
    deltay=(ymax-ymin)/cmd.nquadSteps;
    
    result=0;
    Qa=GaussLegendreQ8_func( ymin );
    
    for (j=1; j<=cmd.nquadSteps; j++) {
        yp = ymin + deltay * j;
        Qb=GaussLegendreQ8_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }
    
    result = 2.0*rlog(10.0)*(rsqr(ki)/FOURPI2)*result;

    return result;
}

// END Q8


// BEGIN Q9

local real Q9_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    real deltay;
    int j;
    real Qa;
    real Qb;
    real yp;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    
    deltay=(ymax-ymin)/cmd.nquadSteps;
    
    result=0;
    Qa=GaussLegendreQ9_func( ymin );
    
    for (j=1; j<=cmd.nquadSteps; j++) {
        yp = ymin + deltay * j;
        Qb=GaussLegendreQ9_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }
    
    result = 2.0*rlog(10.0)*(rsqr(ki)/FOURPI2)*result;
    
    return result;
}

// END Q9

// BEGIN Q13

local real Q13_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    real deltay;
    int j;
    real Qa;
    real Qb;
    real yp;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    deltay=(ymax-ymin)/cmd.nquadSteps;
    
    result=0;
    Qa=GaussLegendreQ13_func( ymin );
    
    for (j=1; j<=cmd.nquadSteps; j++) {
        yp = ymin + deltay * j;
        Qb=GaussLegendreQ13_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }
    
    result = 2.0*rlog(10.0)*(rsqr(ki)/FOURPI2)*result;

    return result;
}

// END Q13

// BEGIN QI

local real QI_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    real deltay;
    int j;
    real Qa;
    real Qb;
    real yp;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    deltay=(ymax-ymin)/cmd.nquadSteps;
    
    result=0;
    Qa=GaussLegendreQI_func( ymin );
    
    for (j=1; j<=cmd.nquadSteps; j++) {
        yp = ymin + deltay * j;
        Qb=GaussLegendreQI_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }
    
    result = 2.0*rlog(10.0)*(rsqr(ki)/FOURPI2)*result;

    return result;
}

// END QI


// BEGIN Q5

local real Q5_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    real deltay;
    int j;
    real Qa;
    real Qb;
    real yp;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    deltay=(ymax-ymin)/cmd.nquadSteps;
    
    result=0;
    Qa=GaussLegendreQ5_func( ymin );
    
    for (j=1; j<=cmd.nquadSteps; j++) {
        yp = ymin + deltay * j;
        Qb=GaussLegendreQ5_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }
    
    result = rlog(10.0)*(rsqr(ki)/FOURPI2)*result;

    return result;
}

// END Q5

// BEGIN Q7

local real Q7_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    real deltay;
    int j;
    real Qa;
    real Qb;
    real yp;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    deltay=(ymax-ymin)/cmd.nquadSteps;
    
    result=0;
    Qa=GaussLegendreQ7_func( ymin );
    
    for (j=1; j<=cmd.nquadSteps; j++) {
        yp = ymin + deltay * j;
        Qb=GaussLegendreQ7_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }
    
    result = rlog(10.0)*(rsqr(ki)/FOURPI2)*result;

    return result;
}

// END Q7

// BEGIN Q11

local real Q11_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    real deltay;
    int j;
    real Qa;
    real Qb;
    real yp;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    deltay=(ymax-ymin)/cmd.nquadSteps;
    
    result=0;
    Qa=GaussLegendreQ11_func( ymin );
    
    for (j=1; j<=cmd.nquadSteps; j++) {
        yp = ymin + deltay * j;
        Qb=GaussLegendreQ11_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }
    
    result = rlog(10.0)*(rsqr(ki)/FOURPI2)*result;

    return result;
}

// END Q11

// BEGIN Q12

local real Q12_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    real deltay;
    int j;
    real Qa;
    real Qb;
    real yp;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    deltay=(ymax-ymin)/cmd.nquadSteps;
    
    result=0;
    Qa=GaussLegendreQ12_func( ymin );
    
    for (j=1; j<=cmd.nquadSteps; j++) {
        yp = ymin + deltay * j;
        Qb=GaussLegendreQ12_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }
    
    result = rlog(10.0)*(rsqr(ki)/FOURPI2)*result;
    
    return result;
}

// END Q12


// BEGIN R1

local real R1_function_AA(real eta, real ki)
{
    real result, fac, s;
    real pmin, pmax, ymin, ymax;
    int j;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    fac=(rlog(10.0)/FOURPI2)*psInterpolation_nr(gd.k, kPS, pPS, nPSLT);
    
    s=0;
    for (j=1;j<=nGL(pGL)/2;j++) {
        gd.x = xGL(pGL)[j];
        result =QROMBERG(funcR1int,ymin,ymax,midpnt,cmd.epsquad,KK);
        s += 2.0*wGL(pGL)[j]*result;
    }
    
    return fac*s;
}

// END R1


// BEGIN R2

local real R2_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    
    real deltay;
    int j;
    real Qa;
    real Qb;
    real yp;
    
    gd.xstop = eta;
    gd.k = ki;
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    deltay=(ymax-ymin)/cmd.nquadSteps;
    
    result=0;
    Qa=GaussLegendreR2_func( ymin );
    
    for (j=1; j<=cmd.nquadSteps; j++) {
        yp = ymin + deltay * j;
        Qb=GaussLegendreR2_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }
    
    result = rlog(10.0)*(rsqr(gd.k)/FOURPI2)*psInterpolation_nr(gd.k, kPS, pPS, nPSLT)*result;
    
    return result;
}

// END R2

// BEGIN RI

local real RI_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    real PSL;
    
    
    real deltay;
    int j;
    real Qa;
    real Qb;
    real yp;
    
    
    
    gd.xstop = eta;
    gd.k = ki;
    PSL = psInterpolation_nr(ki, kPS, pPS, nPSLT);
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    
    deltay=(ymax-ymin)/cmd.nquadSteps;
    
    result=0;
    Qa=GaussLegendreRI_func( ymin );
    
    for (j=1; j<=cmd.nquadSteps; j++) {
        yp = ymin + deltay * j;
        Qb=GaussLegendreRI_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }
    
    
    result = rlog(10.0)*(rsqr(ki)/FOURPI2)*PSL*result;
    
    return result;
}

// END RI


// BEGIN R1p2

local real R1p2_function_AA(real eta, real ki)
{
    real result;
    real pmin, pmax, ymin, ymax;
    real PSL;
    
    real deltay;
    int j;
    real Qa;
    real Qb;
    real yp;
    
    gd.xstop = eta;
    gd.k = ki;
    PSL = psInterpolation_nr(ki, kPS, pPS, nPSLT);
    
    pmin = kPS[11];
    pmax = kPS[nPSLT-9];
    ymin = rlog10(pmin);
    ymax = rlog10(pmax);
    
    deltay=(ymax-ymin)/cmd.nquadSteps;
    
    result=0;
    Qa=GaussLegendreR1p2_func( ymin );
    
    for (j=1; j<=cmd.nquadSteps; j++) {
        yp = ymin + deltay * j;
        Qb=GaussLegendreR1p2_func( yp );
        result += (Qb + Qa)/2. * deltay;
        Qa = Qb;
    }

    result = rlog(10.0)*(rsqr(ki)/FOURPI2)*PSL*result;
    return result;
}

// END R1p2



// BEGIN Qs and Rs

global_QRs QsRs_functions_trapezoid(real eta, real ki)
{
    real Q1p, Q2p, Q3p, Q8p, Q9p, Q13p, QIp;
    real Q5p, Q7p, Q11p, Q12p;
    real R1p, R2p;
    real RI, R1p2;
    global_QRs_ptr QRstmp;
    
//    fprintf(gd.outlog,"\nQuadrature method: %s",gd.quadraturemethod_comment);
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


 global_QRs QsRs_functions_trapezoid_LCDM(real eta, real ki)
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

// END Qs and Rs

#undef abskmq
#undef KK
#undef QROMBERG

// ==============================================================================
