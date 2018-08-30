/*==============================================================================
 MODULE: mglpt_diffeqs.c			[mgpt]
 ==============================================================================*/

#include "globaldefs.h"
#include "protodefs.h"


local void derivsSecondOrder_ver2(double x,double y[],double dydx[]);
local void derivsThirdOrder(double x,double y[],double dydx[]);
local void derivsThirdOrder_ver2(double eta,double y[],double dydx[]);

local global_D3_ptr DsThirdOrder_func(real kf, real k1, real k2);

global void derivsFirstOrder(double x,double y[],double dydx[])
{
    nrhs++;

    dydx[1] = y[2];
    dydx[2]=f1(x)*mu(x,gd.kf)*y[1]-f2(x)*y[2];
}

global void derivsSecondOrder(double x,double y[],double dydx[])
{
    nrhs++;

    dydx[1] = y[2];
    dydx[2]=f1(x)*mu(x,gd.k1)*y[1]-f2(x)*y[2];
//
    dydx[3] = y[4];
    dydx[4]=f1(x)*mu(x,gd.k2)*y[3]-f2(x)*y[4];
//
    dydx[5] = y[6];
    dydx[6]=f1(x)*mu(x,gd.kf)*y[5]-f2(x)*y[6]
        + f1(x)*mu(x,gd.kf)*y[1]*y[3];
//
    dydx[7] = y[8];
    dydx[8]=f1(x)*mu(x,gd.kf)*y[7]-f2(x)*y[8]
    + f1(x)*(mu(x,gd.k1)+mu(x,gd.k2)-mu(x,gd.kf))*y[1]*y[3];
//
    dydx[9] = y[10];
    dydx[10]=f1(x)*mu(x,gd.kf)*y[9]-f2(x)*y[10]
    + f1(x)
    *(sqr(mass(x))/PiF(x,gd.kf))
    *KFL(x,gd.kf,gd.k1,gd.k2)
    *y[1]*y[3];
//
    dydx[11] = y[12];
    dydx[12]=f1(x)*mu(x,gd.kf)*y[11]-f2(x)*y[12]
            + (1.0/6.0)
                *sqr(
                        (OmM(x)*H(x))/(rexp(x)*H02)
                     )
                *(
                    (sqr(gd.kf)*M2(x))/( PiF(x,gd.kf)*PiF(x,gd.k1)*PiF(x,gd.k2) )
                )
                *y[1]*y[3];
}

local void derivsSecondOrder_ver2(double x,double y[],double dydx[])
{
    nrhs++;
    
    dydx[1] = y[2];
    dydx[2]=f1(x)*mu(x,gd.k1)*y[1]-f2(x)*y[2];
    //
    dydx[3] = y[4];
    dydx[4]=f1(x)*mu(x,gd.k2)*y[3]-f2(x)*y[4];
    //
    dydx[5] = y[6];
    dydx[6]=f1(x)*mu(x,gd.kf)*y[5]-f2(x)*y[6]
    + sourceA(x,gd.kf,gd.k1,gd.k2)*y[1]*y[3];
    //
    dydx[7] = y[8];
    dydx[8]=f1(x)*mu(x,gd.kf)*y[7]-f2(x)*y[8]
    + sourceb(x,gd.kf,gd.k1,gd.k2)*y[1]*y[3];
}

local void derivsThirdOrder(double x,double y[],double dydx[])
{
    nrhs++;
    
    dydx[1] = y[2];
    dydx[2]=f1(x)*mu(x,gd.k1)*y[1]-f2(x)*y[2];
//
    dydx[3] = y[4];
    dydx[4]=f1(x)*mu(x,gd.k2)*y[3]-f2(x)*y[4];
//
    dydx[5] = y[6];
    dydx[6]=f1(x)*mu(x,gd.kf)*y[5]-f2(x)*y[6]
            + f1(x)*mu(x,gd.kf)*y[1]*y[3];
//
    dydx[7] = y[8];
    dydx[8]=f1(x)*mu(x,gd.kf)*y[7]-f2(x)*y[8]
        + f1(x)*(mu(x,gd.k1)+mu(x,gd.k2)-mu(x,gd.kf))*y[1]*y[3];
//
    dydx[9] = y[10];
    dydx[10]=f1(x)*mu(x,gd.kf)*y[9]-f2(x)*y[10]
            + f1(x)
                *(sqr(mass(x))/PiF(x,gd.kf))
                *KFL(x,gd.kf,gd.k1,gd.k2)
                *y[1]*y[3];
//
    dydx[11] = y[12];
    dydx[12]=f1(x)*mu(x,gd.kf)*y[11]-f2(x)*y[12]
            + (1.0/6.0)
                *sqr(
                        (OmM(x)*H(x))/(rexp(x)*H02)
                     )
                *(
                    (sqr(gd.kf)*M2(x))/( PiF(x,gd.kf)*PiF(x,gd.k1)*PiF(x,gd.k2) )
                )
                *y[1]*y[3];
//
    dydx[13] = y[14];
    dydx[14]=f1(x)*mu(x,gd.k1)*y[13]-f2(x)*y[14]
            +f1(x)*(mu(x,gd.k2)-mu(x,gd.k1))*y[5]*y[3]
            +(dydx[6]+f2(x)*y[6])*y[3];
//
    dydx[15] = y[16];
    dydx[16]=f1(x)*mu(x,gd.k1)*y[15]-f2(x)*y[16]
            +f1(x)*(mu(x,gd.k2)-mu(x,gd.k1))*y[7]*y[3]
            +(dydx[8]+f2(x)*y[8])*y[3];
//
    dydx[17] = y[18];
    dydx[18]=f1(x)*mu(x,gd.k1)*y[17]-f2(x)*y[18]
            +f1(x)*(mu(x,gd.k2)-mu(x,gd.k1))*y[9]*y[3]
            +(dydx[10]+f2(x)*y[10])*y[3];
//
    dydx[19] = y[20];
    dydx[20]=f1(x)*mu(x,gd.k1)*y[19]-f2(x)*y[20]
            +f1(x)*(mu(x,gd.k2)-mu(x,gd.k1))*y[11]*y[3]
            +(dydx[12]+f2(x)*y[12])*y[3];
//
}

local void derivsThirdOrder_ver2(double eta,double y[],double dydx[])
{
    nrhs++;
    real Dpk, Dpp, D2f, D2mf, kplusp, kpluspm;

    kplusp = kpp(gd.x,gd.k,gd.p);
    kpluspm = kpp(-gd.x,gd.k,gd.p);
    Dpk  = y[1];
    Dpp  = y[3];
    D2f  = y[5];
    D2mf = y[7];

    dydx[1] = y[2];
    dydx[2] = f1(eta)*mu(eta,gd.k)*y[1]-f2(eta)*y[2];
//
    dydx[3] = y[4];
    dydx[4] = f1(eta)*mu(eta,gd.p)*y[3]-f2(eta)*y[4];
//
    dydx[5] = y[6];
    dydx[6] = f1(eta)*mu(eta,kplusp)*y[5]-f2(eta)*y[6] //;
    + SD2(eta,gd.x,gd.k,gd.p)*y[1]*y[3];
//
    dydx[7] = y[8];
    dydx[8] = f1(eta)*mu(eta,kpluspm)*y[7]-f2(eta)*y[8] //;
    + SD2(eta,-gd.x,gd.k,gd.p)*y[1]*y[3];
//
    dydx[9] = y[10];
    dydx[10] = f1(eta)*mu(eta,gd.k)*y[9]-f2(eta)*y[10]
    + S3I( eta,gd.x,gd.k,gd.p,Dpk,Dpp,D2f,D2mf)
    + S3II(eta,gd.x,gd.k,gd.p,Dpk,Dpp,D2f,D2mf)
    + S3FL(eta,gd.x,gd.k,gd.p,Dpk,Dpp,D2f,D2mf)
    + S3dI(eta,gd.x,gd.k,gd.p,Dpk,Dpp,D2f,D2mf);
}

global real DpFunction(real k)
{
    int nbad,nok;
    double *ystart;
    real Dptmp;

    gd.kf = k;
    
    ystart=dvector(1,NEQS1Order);
    xp=dvector(1,200);
    yp=dmatrix(1,NEQS1Order,1,200);

    ystart[1]=rexp(gd.xnow);
    ystart[2]=rexp(gd.xnow);

    nrhs=0;
    kmax=100;
    dxsav=(cmd.xstop-gd.xnow)/20.0;
//
    integration(ystart,NEQS1Order,gd.xnow,cmd.xstop,cmd.eps,gd.dx,cmd.dxmin,
                 &nok,&nbad,cmd.maxnsteps,derivsFirstOrder);

    Dptmp = yp[1][kount];
    
    free_dmatrix(yp,1,NEQS1Order,1,200);
    free_dvector(xp,1,200);
    free_dvector(ystart,1,NEQS1Order);

    return (Dptmp);
}

global global_D2_ptr DsSecondOrder_func(real kf, real k1, real k2)
{
    int nbad,nok;
    double *ystart;
    global_D2_ptr ptmp;
//
    gd.kf = kf;
    gd.k1 = k1;
    gd.k2 = k2;
//
    ystart=dvector(1,NEQS2Order);
    xp=dvector(1,200);
    yp=dmatrix(1,NEQS2Order,1,200);
    
    ystart[1]=rexp(gd.xnow);
    ystart[2]=rexp(gd.xnow);
    ystart[3]=rexp(gd.xnow);
    ystart[4]=rexp(gd.xnow);
    ystart[5]=3.0*rexp(2.0*gd.xnow)/7.0;
    ystart[6]=6.0*rexp(2.0*gd.xnow)/7.0;
    ystart[7]=3.0*rexp(2.0*gd.xnow)/7.0;
    ystart[8]=6.0*rexp(2.0*gd.xnow)/7.0;
    ystart[9]=0.0;
    ystart[10]=0.0;
    ystart[11]=0.0;
    ystart[12]=0.0;

    nrhs=0;
    kmax=100;
    dxsav=(cmd.xstop-gd.xnow)/100.0;

    integration(ystart,NEQS2Order,gd.xnow,cmd.xstop,cmd.eps,gd.dx,cmd.dxmin,
                &nok,&nbad,cmd.maxnsteps,derivsSecondOrder);

    ptmp = (global_D2_ptr) allocate(1 * sizeof(global_D2));

    etaD2(ptmp) = xp[kount];
    Dpk1D2(ptmp) = yp[1][kount];
    Dpk2D2(ptmp) = yp[3][kount];
    Dak1k2D2(ptmp) = yp[5][kount];
    Dbk1k2D2(ptmp) = yp[7][kount];
    DFLk1k2D2(ptmp) = yp[9][kount];
    DdIk1k2D2(ptmp) = yp[11][kount];
//
    free_dmatrix(yp,1,NEQS2Order,1,200);
    free_dvector(xp,1,200);
    free_dvector(ystart,1,NEQS2Order);
    
    return ptmp;
}

global global_D2v2_ptr DsSecondOrder_func_ver2(real kf, real k1, real k2)
{
    int nbad,nok;
    double *ystart;
    global_D2v2_ptr ptmp;
//
    gd.kf = kf;
    gd.k1 = k1;
    gd.k2 = k2;
//
    ystart=dvector(1,NEQS2Orderv2);
    xp=dvector(1,200);
    yp=dmatrix(1,NEQS2Orderv2,1,200);
    
    ystart[1]=rexp(gd.xnow);
    ystart[2]=rexp(gd.xnow);
    ystart[3]=rexp(gd.xnow);
    ystart[4]=rexp(gd.xnow);
    ystart[5]=3.0*rexp(2.0*gd.xnow)/7.0;
    ystart[6]=6.0*rexp(2.0*gd.xnow)/7.0;
    ystart[7]=3.0*rexp(2.0*gd.xnow)/7.0;
    ystart[8]=6.0*rexp(2.0*gd.xnow)/7.0;
    
    nrhs=0;
    kmax=100;
    dxsav=(cmd.xstop-gd.xnow)/20.0;

    integration(ystart,NEQS2Orderv2,gd.xnow,cmd.xstop,cmd.eps,gd.dx,cmd.dxmin,
                &nok,&nbad,cmd.maxnsteps,derivsSecondOrder_ver2);
    
    ptmp = (global_D2v2_ptr) allocate(1 * sizeof(global_D2v2));

    etaD2(ptmp) = xp[kount];
    Dpk1D2(ptmp) = yp[1][kount];
    Dpk2D2(ptmp) = yp[3][kount];
    DA2D2(ptmp) = yp[5][kount];
    DB2D2(ptmp) = yp[7][kount];
//
    free_dmatrix(yp,1,NEQS2Orderv2,1,200);
    free_dvector(xp,1,200);
    free_dvector(ystart,1,NEQS2Orderv2);
    
    return ptmp;
}


local global_D3_ptr DsThirdOrder_func(real kf, real k1, real k2)
{
    int nbad,nok;
    double *ystart;
    global_D3_ptr ptmp;
//
    gd.kf = kf;
    gd.k1 = k1;
    gd.k2 = k2;
//
    ystart=dvector(1,NEQS3Order);
    xp=dvector(1,200);
    yp=dmatrix(1,NEQS3Order,1,200);
    
    ystart[1]=rexp(gd.xnow);
    ystart[2]=rexp(gd.xnow);
    ystart[3]=rexp(gd.xnow);
    ystart[4]=rexp(gd.xnow);
//
    ystart[5]=3.0*rexp(2.0*gd.xnow)/7.0;
    ystart[6]=6.0*rexp(2.0*gd.xnow)/7.0;
    ystart[7]=3.0*rexp(2.0*gd.xnow)/7.0;
    ystart[8]=6.0*rexp(2.0*gd.xnow)/7.0;
//
    ystart[9]=0.0;
    ystart[10]=0.0;
    ystart[11]=0.0;
    ystart[12]=0.0;
//
    ystart[13]=5.0*rexp(2.0*gd.xnow)/21.0;;
    ystart[14]=15.0*rexp(2.0*gd.xnow)/21.0;;
    ystart[15]=5.0*rexp(2.0*gd.xnow)/21.0;;
    ystart[16]=15.0*rexp(2.0*gd.xnow)/21.0;;
//
    ystart[17]=0.0;
    ystart[18]=0.0;
    ystart[19]=0.0;
    ystart[20]=0.0;
//
    nrhs=0;
    kmax=100;
    dxsav=(cmd.xstop-gd.xnow)/100.0;
    integration(ystart,NEQS3Order,gd.xnow,cmd.xstop,cmd.eps,gd.dx,cmd.dxmin,
                &nok,&nbad,cmd.maxnsteps,derivsThirdOrder);
    
    ptmp = (global_D3_ptr) allocate(1 * sizeof(global_D3));

    etaD3(ptmp) = xp[kount];
    DpkD3(ptmp) = yp[1][kount];
    DppD3(ptmp) = yp[3][kount];
    Da2D3(ptmp) = yp[5][kount];
    Db2D3(ptmp) = yp[7][kount];
    DFL2D3(ptmp) = yp[9][kount];
    DdI2D3(ptmp) = yp[11][kount];
    DIA3D3(ptmp) = yp[13][kount];
    DIB3D3(ptmp) = yp[15][kount];
    DIFL3D3(ptmp) = yp[17][kount];
    DIdI3D3(ptmp) = yp[19][kount];
//
    free_dmatrix(yp,1,NEQS3Order,1,200);
    free_dvector(xp,1,200);
    free_dvector(ystart,1,NEQS3Order);
    
    return ptmp;
}


global global_D3v2_ptr DsThirdOrder_func_ver2(real x, real k, real p)
{
    int nbad,nok;
    double *ystart;
    global_D3v2_ptr ptmp;

    gd.x = x;
    gd.k = k;
    gd.p = p;
//
    ystart=dvector(1,NEQS3Orderv2);
    xp=dvector(1,200);
    yp=dmatrix(1,NEQS3Orderv2,1,200);
//
    ystart[1]=rexp(gd.xnow);
    ystart[2]=rexp(gd.xnow);
    ystart[3]=rexp(gd.xnow);
    ystart[4]=rexp(gd.xnow);
//
    ystart[5]=(3.0/7.0)*rexp(2.0*gd.xnow)*(1.0-rsqr(x));
    ystart[6]=(6.0/7.0)*rexp(2.0*gd.xnow)*(1.0-rsqr(x));
    ystart[7]=(3.0/7.0)*rexp(2.0*gd.xnow)*(1.0-rsqr(x));
    ystart[8]=(6.0/7.0)*rexp(2.0*gd.xnow)*(1.0-rsqr(x));
//
    ystart[9]=(5.0/(7.0*9.0))*rexp(3.0*gd.xnow)*rsqr(1.0-rsqr(x))
        *(
            1.0/(1.0 + rsqr(p/k) + 2.0 * (p/k) * x)
          + 1.0/(1.0 + rsqr(p/k) - 2.0 * (p/k) * x)
        );
//
    ystart[10]=(15.0/(7.0*9.0))*rexp(3.0*gd.xnow)*rsqr(1.0-rsqr(x))
        *(
            1.0/(1.0 + rsqr(p/k) + 2.0 * (p/k) * x)
          + 1.0/(1.0 + rsqr(p/k) - 2.0 * (p/k) * x)
          );
//
    nrhs=0;
    kmax=100;
    dxsav=(cmd.xstop-gd.xnow)/20.0;
    integration(ystart,NEQS3Orderv2,gd.xnow,cmd.xstop,cmd.eps,gd.dx,cmd.dxmin,
                &nok,&nbad,cmd.maxnsteps,derivsThirdOrder_ver2);

    ptmp = (global_D3v2_ptr) allocate(1 * sizeof(global_D3v2));

    etaD3v2(ptmp) = xp[kount];
    DpkD3v2(ptmp) = yp[1][kount];
    DppD3v2(ptmp) = yp[3][kount];
    D2fD3v2(ptmp) = yp[5][kount];
    D2mfD3v2(ptmp) = yp[7][kount];
    D3symmD3v2(ptmp) = yp[9][kount];
//
    free_dmatrix(yp,1,NEQS3Orderv2,1,200);
    free_dvector(xp,1,200);
    free_dvector(ystart,1,NEQS3Orderv2);
    
    return ptmp;
}

#define BSSTEP          0
#define NULLMETHOD      1
#define RKQS            2

global void integration(double ystart[], int nvar, double x1, double x2, double eps, double h1,
                       double hmin, int *nok, int *nbad, int maxnsteps,
                       void (*derivsin)(double, double [], double []))
{
    switch (gd.method_int) {
        case BSSTEP:
            odeint(ystart,nvar,gd.xnow,cmd.xstop,cmd.eps,gd.dx,cmd.dxmin,
                   nok,nbad,maxnsteps,derivsin,bsstep);
            break;
//
        case RKQS:
            odeint(ystart,nvar,gd.xnow,cmd.xstop,cmd.eps,gd.dx,cmd.dxmin,
                   nok,nbad,maxnsteps,derivsin,rkqs);
            break;
//
        case NULLMETHOD:
            odeint(ystart,nvar,gd.xnow,cmd.xstop,cmd.eps,gd.dx,cmd.dxmin,
                   nok,nbad,maxnsteps,derivsin,bsstep);
            break;
//
        default:
            odeint(ystart,nvar,gd.xnow,cmd.xstop,cmd.eps,gd.dx,cmd.dxmin,
                   nok,nbad,maxnsteps,derivsin,bsstep);
            break;
    }
}

void integration_method_string_to_int(string method_str,int *method_int)
{
    *method_int=-1;
    if (strcmp(method_str,"bsstep") == 0) {
        *method_int = BSSTEP;
        strcpy(gd.integration_method_comment, "bsstep integration method");
    }
//
    if (strcmp(method_str,"rkqs") == 0) {
        *method_int = RKQS;
        strcpy(gd.integration_method_comment, "rkqs integration method");
    }
//
    if (strnull(method_str)) {
        *method_int = NULLMETHOD;
        strcpy(gd.integration_method_comment,
               "null integration method ... running deafult (bsstep)");
        fprintf(stdout,"\n\tintegration: default integration method (bsstep)...\n");
    }
//
    if (*method_int == -1) {
        *method_int = BSSTEP;
        strcpy(gd.integration_method_comment,
               "Unknown integration method ... running deafult (bsstep)");
        fprintf(stdout,"\n\tintegration: Unknown method... %s ",cmd.integration_method);
        fprintf(stdout,
                "\n\trunning default integration method (bsstep)...\n");
    }
}

#undef BSSTEP
#undef RKQS
#undef NULLMETHOD






