/*==============================================================================
	MODULE: models.c			[mgpt]
	Written by: Mario A. Rodriguez-Meza
	Starting date: January 2018
	Purpose: Rutines to create several types of models
	Language: C
	Use: 'set_model();'
	Routines and functions:
	Modules, routines and external headers:
	Coments and notes:
	Info: Mario A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.inin.gob.mx

	Mayor revisions: January 22, 2018
	Copyright: (c) 2005-2018 Mario A. Rodriguez-Meza.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#define global // check the behaviour of this
#include "globaldefs.h"

#include "protodefs.h"
#include "models.h"

local void model_string_to_int(string, int *);

// Models
local void Model_HS(void);
local void Model_fR1(void);

local real sourcea(real eta, real kf);
local real sourceFL(real eta, real kf, real k1, real k2);
local real sourcedI(real eta, real kf, real k1, real k2);


local real S3IIplus(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f);
local real S3IIminus(real eta, real x, real k, real p, real Dpk, real Dpp, real D2mf);
local real S3FLplus(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f);
local real S3FLminus(real eta, real x, real k, real p, real Dpk, real Dpp, real D2mf);
local real D2phiplus(real eta, real x, real k, real p,
                     real Dpk, real Dpp, real D2f);
local real D2phiminus(real eta, real x, real k, real p,
                      real Dpk, real Dpp, real D2mf);
local real K3dI(real eta, real x, real k,  real p,
                 real Dpk, real Dpp, real D2f, real D2mf);


#define HS                          0
#define fR1                         1

global void set_model(void)
{
	int model_int;

	model_string_to_int(cmd.mgmodel, &model_int);
	switch (model_int){

	case HS: Model_HS(); break;
	case fR1: Model_fR1(); break;

	default: error("\nUnknown model type %s\n\n",cmd.mgmodel);
	}
}

local void model_string_to_int(string model_str,int *model_int)
{
	*model_int = -1;
	if (strcmp(model_str,"HS") == 0)				*model_int=HS;
	if (strcmp(model_str,"fR1") == 0)               *model_int=fR1;
}

#undef HS
#undef fR1


// MODELS -------------------

// ==========================================
// Begin: Hu-Sawicky Model


local void Model_HS(void)
{
    strcpy(gd.model_comment, "HS Model");
//
// Parameters set and default values:
//    nHS=1.0;
//    fR0=1.0e-5;
//    beta2=1.0/6.0;
//    omegaBD = 0.0;
//    screening = 1.0;
    gd.beta2 = 1.0/6.0;
    cmd.omegaBD = 0.0;
}

global real mass(real eta)
{
    real masstmp;
    
    masstmp = (1.0/H02)*rsqrt(1.0/(2.0*cmd.fR0))
                * rpow(cmd.om*rexp(-3.0*eta) + 4.0*(1.0-cmd.om),(2.0+cmd.nHS)/2.0)
                / rpow(cmd.om+4.0*(1.0-cmd.om),(1.0+cmd.nHS)/2.0);

    return (masstmp);
}

global real mu(real eta, real k)
{
    real mutmp;
    mutmp = 1.0 + (2.0*gd.beta2*k*k)/(k*k + rexp(2.0*eta)*rsqr(mass(eta)));

    return (mutmp);
}

global real PiF(real eta, real k)
{
    real PiFtmp;
    
    PiFtmp = k*k/rexp(2.0*eta) + rsqr(mass(eta));
    
    return (PiFtmp);
}


// ===========================================================================
// BEGIN :: SECOND ORDER (six second order differential equations)
//

global real M2(real eta)
{
    real M2tmp;
    
    M2tmp = (9.0/(4.0*rsqr(H02)))*rsqr(1.0/cmd.fR0)
    * rpow(cmd.om*rexp(-3.0*eta)+4.0*(1.0-cmd.om),5.0)
    / rpow(cmd.om+4*(1.0-cmd.om),4.0);
    
    return (M2tmp);
}


global real KFL(real eta, real k, real k1, real k2)
{
    real KFLtmp;

    KFLtmp = 0.5*(rsqr(sqr(k)-sqr(k1)-sqr(k2))/(sqr(k1)*sqr(k2)))*(mu(eta,k1)+mu(eta,k2)-2.0)
        + 0.5*((sqr(k)-sqr(k1)-sqr(k2))/sqr(k1))*(mu(eta,k1)-1.0)
        + 0.5*((sqr(k)-sqr(k1)-sqr(k2))/sqr(k2))*(mu(eta,k2)-1.0);

    return (KFLtmp);
}

global real sourceA(real eta, real kf, real k1, real k2)
{
    real Stmp;

    Stmp = sourcea(eta, kf)
            + sourceFL(eta, kf, k1, k2)
            - cmd.screening * sourcedI(eta, kf, k1, k2);

    return Stmp;
}

local real sourcea(real eta, real kf)
{
    real Stmp;
    
    Stmp = f1(eta)*mu(eta, kf);

    return Stmp;
}

global real sourceb(real eta, real kf, real k1, real k2)
{
    real Stmp;
    
    Stmp = f1(eta)*( mu(eta, k1) + mu(eta, k2) - mu(eta, kf));

    return Stmp;
}

local real sourceFL(real eta, real kf, real k1, real k2)
{
    real Stmp;
    
    Stmp = f1(eta)*(rsqr(mass(eta))/PiF(eta,kf))*KFL(eta, kf, k1, k2);

    return Stmp;
}

local real sourcedI(real eta, real kf, real k1, real k2)
{
    real Stmp;
    
    Stmp = (1.0/6.0)*rsqr(OmM(eta)*H(eta)/(rexp(eta)*H02))
            * rsqr(kf)*M2(eta)/ ( PiF(eta,kf)*PiF(eta,k1)*PiF(eta,k2) );
    
    return Stmp;
}

//
// END :: SECOND ORDER (six second order differential equations)
// ===========================================================================


// ===========================================================================
// BEGIN :: THIRD ORDER (Dsymmetric, five second order differential equations)
//

global real M1(real eta)
{
    real M1tmp;
    
    M1tmp = 3.0*rsqr(mass(eta));
    
    return (M1tmp);
}

global real M3(real eta)
{
    real M3tmp;
    
    M3tmp = (45.0/(8.0*rsqr(H02)))*rpow(1.0/cmd.fR0,3.0)
    * rpow(cmd.om*rexp(-3.0*eta)+4.0*(1.0-cmd.om),7.0)
    / rpow(cmd.om+4*(1.0-cmd.om),6.0);
    
    return (M3tmp);
}

global real kpp(real x, real k, real p)
{
    real kpptmp;
    
    kpptmp = rsqrt(rsqr(k) + rsqr(p) + 2.0*k*p*x);
    
    return kpptmp;
}

global real KFL2(real eta, real x, real k, real p)
{
    real KFLtmp;
    
    KFLtmp = 2.0*rsqr(x)*( mu(eta,k) + mu(eta,p)-2.0 )
            + (p*x/k)*(mu(eta,k) - 1.0)
            + (k*x/p)*(mu(eta,p) - 1.0);

    return (KFLtmp);
}

global real JFL(real eta, real x, real k, real p)
{
    real JFLtmp;
    
    JFLtmp = (9.0/(2.0*A0(eta)))
            * KFL2(eta, x, k, p) * PiF(eta, k) * PiF(eta, p);

    return (JFLtmp);
}

local real D2phiplus(real eta, real x, real k, real p,
                     real Dpk, real Dpp, real D2f)
{
    real D2tmp;
    
    D2tmp = (
             (1.0 + rsqr(x))
             -( 2.0*A0(eta)/3.0 )
             * (
                ( M2(eta) + JFL(eta,x,k,p)*(3.0+2.0*cmd.omegaBD) )
                / (3.0*PiF(eta,k)*PiF(eta,p))
                )
             ) * Dpk*Dpp + D2f;
    
    return (D2tmp);
}

local real D2phiminus(real eta, real x, real k, real p,
                     real Dpk, real Dpp, real D2mf)
{
    real D2tmp;
    
    D2tmp = (
             (1.0 + rsqr(x))
             -( 2.0*A0(eta)/3.0 )
             * (
                ( M2(eta) + JFL(eta,-x,k,p)*(3.0+2.0*cmd.omegaBD) )
                / (3.0*PiF(eta,k)*PiF(eta,p))
                )
             ) * Dpk*Dpp + D2mf;

    return (D2tmp);
}

local real K3dI(real eta, real x, real k,  real p,
                 real Dpk, real Dpp, real D2f, real D2mf)
{
    real K3tmp, kplusp, kpluspm;
    real t1, t2, t3, t4, t5, t6;

    kplusp = kpp(x,k,p);
    kpluspm = kpp(-x,k,p);

    t1 = 2.0*rsqr(OmM(eta)*H(eta)/H02)
            *(M2(eta)/(PiF(eta,k)*PiF(eta,0)));

    t2 = (1.0/3.0)*(rpow(OmM(eta),3.0)*rpow(H(eta),4.0)/rpow(H02,4) )
        *(
            M3(eta) - M2(eta)*(M2(eta) + JFL(eta,-1.0,p,p)*(3.0+2.0*cmd.omegaBD))
                        /(PiF(eta,0))
          ) / ( rsqr(PiF(eta,p)) * PiF(eta,k) );

    t3 = rsqr(OmM(eta)*H(eta)/H02)
        *(M2(eta)/(PiF(eta,p)*PiF(eta,kplusp)))
        *(
            1.0 + rsqr(x) + (D2f)/(Dpk*Dpp)
          );
    
    t4 = (1.0/3.0)*(rpow(OmM(eta),3.0)*rpow(H(eta),4.0)/rpow(H02,4) )
        *(
            M3(eta) - M2(eta)*(M2(eta) + JFL(eta,x,k,p)*(3.0+2.0*cmd.omegaBD))
                        /(PiF(eta,kplusp))
          ) / ( rsqr(PiF(eta,p)) * PiF(eta,k) );
    
    t5 = rsqr(OmM(eta)*H(eta)/H02)
        *(M2(eta)/(PiF(eta,p)*PiF(eta,kpluspm)))
        *(
            1.0 + rsqr(x) + (D2mf)/(Dpk*Dpp)
          );
    
    t6 = (1.0/3.0)*(rpow(OmM(eta),3.0)*rpow(H(eta),4.0)/rpow(H02,4) )
        *(
          M3(eta) - M2(eta)*(M2(eta) + JFL(eta,-x,k,p)*(3.0+2.0*cmd.omegaBD))
                /(PiF(eta,kpluspm))
          ) / ( rsqr(PiF(eta,p)) * PiF(eta,k) );

    K3tmp = t1 + t2 + t3 + t4 + t5 + t6;

    return (K3tmp);
}

global real S2a(real eta, real x, real k, real p)
{
    real Dtmp, kplusp;

    kplusp = kpp(x,k,p);
    Dtmp = f1(eta)*mu(eta,kplusp);

    return Dtmp;
}

global real S2b(real eta, real x, real k, real p)
{
    real Dtmp, kplusp;
    
    kplusp = kpp(x,k,p);
    Dtmp = f1(eta)*(mu(eta,k)+mu(eta,p)-mu(eta,kplusp));
    
    return Dtmp;
}

global real S2FL(real eta, real x, real k, real p)
{
    real Dtmp, kplusp;
    
    kplusp = kpp(x,k,p);
    Dtmp = f1(eta)*(
                    M1(eta)/(3.0*PiF(eta,kplusp))
                    *KFL2(eta,x,k,p)
                    );
    return Dtmp;
}

global real S2dI(real eta, real x, real k, real p)
{
    real Dtmp, kplusp;
    
    kplusp = kpp(x,k,p);
    Dtmp = (1.0/6.0)*
            rsqr(OmM(eta)*H(eta)/(rexp(eta)*H02))
        * ( (rsqr(kplusp)*M2(eta)) / (PiF(eta,kplusp)*PiF(eta,k)*PiF(eta,p)));

    return Dtmp;
}

global real SD2(real eta, real x, real k, real p)
{
    real Dtmp;

    Dtmp = S2a(eta, x, k, p) -  S2b(eta, x, k, p)*rsqr(x)
        + S2FL(eta, x, k, p) - S2dI(eta, x, k, p);

    return Dtmp;
}

global real S3I(real eta, real x, real k, real p, real Dpk, real Dpp,
                real D2f, real D2mf)
{
    real Stmp, kplusp, kpluspm;

    kplusp = kpp(x,k,p);
    kpluspm = kpp(-x,k,p);
    Stmp = (
            f1(eta)*(mu(eta,p)+mu(eta,kplusp)-mu(eta,k))*D2f*Dpp
                + SD2(eta,x,k,p)*Dpk*Dpp*Dpp
            )*(1.0 - rsqr(x))/(1.0 + rsqr(p/k) + 2.0*(p/k)*x)
        + (
           f1(eta)*(mu(eta,p)+mu(eta,kpluspm)-mu(eta,k))*D2mf*Dpp
            + SD2(eta,-x,k,p)*Dpk*Dpp*Dpp
           )*(1.0 - rsqr(x))/(1.0 + rsqr(p/k) - 2.0*(p/k)*x);

    return (Stmp);
}

local real S3IIplus(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f)
{
    real Stmp, kplusp;
    
    kplusp = kpp(x,k,p);
    
    Stmp =
    -f1(eta)*(mu(eta,p)+mu(eta,kplusp)-2.0*mu(eta,k))
    * Dpp*( D2f + Dpk*Dpp*rsqr(x) )
    
    -f1(eta)*(mu(eta,kplusp)-mu(eta,k))*Dpk*Dpp*Dpp
    
    -(
      (M1(eta)/(3.0*PiF(eta,kplusp))) * f1(eta)*KFL2(eta,x,k,p)
      -rsqr(OmM(eta)*H(eta)/H02)
      * (M2(eta)*kplusp*kplusp*rexp(-2.0*eta))
      / (6.0*PiF(eta,kplusp)*PiF(eta,k)*PiF(eta,p))
      )*Dpk*Dpp*Dpp;
    
    return (Stmp);
}

local real S3IIminus(real eta, real x, real k, real p, real Dpk, real Dpp, real D2mf)
{
    real Stmp, kpluspm;

    kpluspm = kpp(-x,k,p);

    Stmp =
    -f1(eta)*(mu(eta,p)+mu(eta,kpluspm)-2.0*mu(eta,k))
    * Dpp*( D2mf + Dpk*Dpp*rsqr(x) )
    
    -f1(eta)*(mu(eta,kpluspm)-mu(eta,k))*Dpk*Dpp*Dpp
    
    -(
      (M1(eta)/(3.0*PiF(eta,kpluspm))) * f1(eta)*KFL2(eta,-x,k,p)
      -rsqr(OmM(eta)*H(eta)/H02)
      * (M2(eta)*kpluspm*kpluspm*rexp(-2.0*eta))
      / (6.0*PiF(eta,kpluspm)*PiF(eta,k)*PiF(eta,p))
      )*Dpk*Dpp*Dpp;
    
    return (Stmp);
}

global real S3II(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f, real D2mf)
{
    real Stmp;
    
    Stmp =  S3IIplus(eta, x, k, p, Dpk, Dpp, D2f)
         + S3IIminus(eta, x, k, p, Dpk, Dpp, D2mf);

    return (Stmp);
}

local real S3FLplus(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f)
{
    real Stmp, kplusp;
    
    kplusp = kpp(x,k,p);
    
    Stmp = f1(eta)*(M1(eta)/(3.0*PiF(eta,k)))
    *(
        (2.0*rsqr(p+k*x)/rsqr(kplusp) - 1.0 - (k*x)/p )
        *( mu(eta,p)-1.0 )* D2f * Dpp
      
        + ( (rsqr(p) + 3.0*k*p*x + 2.0*k*k * x*x)/rsqr(kplusp) )
        *( mu(eta,kplusp) - 1.0) * D2phiplus(eta,x,k,p,Dpk,Dpp,D2f) * Dpp
      
      + 3.0*rsqr(x)*( mu(eta,k) + mu(eta,p) - 2.0 ) * Dpk * Dpp * Dpp
    );
    
    return (Stmp);
}

local real S3FLminus(real eta, real x, real k, real p, real Dpk, real Dpp, real D2mf)
{
    real Stmp, kpluspm;
    
    kpluspm = kpp(-x,k,p);
    
    Stmp = f1(eta)*(M1(eta)/(3.0*PiF(eta,k)))
    *(
      (2.0*rsqr(p-k*x)/rsqr(kpluspm) - 1.0 + (k*x)/p )
      *( mu(eta,p)-1.0 )* D2mf * Dpp
      
      + ( (rsqr(p) - 3.0*k*p*x + 2.0*k*k * x*x)/rsqr(kpluspm) )
      *( mu(eta,kpluspm) - 1.0) * D2phiminus(eta,x,k,p,Dpk,Dpp,D2mf) * Dpp
      
      + 3.0*rsqr(x)*( mu(eta,k) + mu(eta,p) - 2.0 ) * Dpk * Dpp * Dpp
      );
    
    return (Stmp);
}

global real S3FL(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f, real D2mf)
{
    real Stmp;

    Stmp = S3FLplus(eta, x, k, p, Dpk, Dpp, D2f)
        + S3FLminus(eta, x, k, p, Dpk, Dpp, D2mf);

    return (Stmp);
}

global real S3dI(real eta, real x, real k, real p, real Dpk, real Dpp,
                 real D2f, real D2mf)
{
    real Stmp;
    
    Stmp = -(rsqr(k)/rexp(2.0*eta))
        *(1.0/(6.0*PiF(eta,k)))
        *K3dI(eta,x,k,p,Dpk,Dpp,D2f,D2mf)*Dpk*Dpp*Dpp;

    return (Stmp);
}

//
// END :: THIRD ORDER (Dsymmetric, five second order differential equations)
// ===========================================================================


// End: Hu-Sawicky Model
// ==========================================



// ==========================================
// Begin: fR1 Model

local void Model_fR1(void)
{
    strcpy(gd.model_comment, "fR1 Model");
    
    //
    // Parameters set and default values:
    //    nHS=1.0;
    //    fR0=1.0e-5;
    //    beta2=1.0/6.0;
    //    omegaBD = 0.0;
    //    screening = 1.0;

    //
    // Array dictionary:
    //    param[0] =    nHS;
    //    param[1] =    fR0;
    //    param[2] =    beta2;
    //    param[3] =    omegaBD;
    //    param[4] =    screening;
}

// End: fR1 Model
// ==========================================



// ==========================================
// Begin: XXX Model

// End: XXX Model
// ==========================================



