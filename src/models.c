/*==============================================================================
	MODULE: models.c			[mgpt]
==============================================================================*/

#define global
#include "globaldefs.h"
#include "protodefs.h"
#include "models.h"

// ===========================================================================
// MODELS SECTION 1
// ==========================================
// Begin: HS Model global HEADERS -> local
local void set_Model_HS(void);
local real mu_HS(real eta, real k);
local real sourceA_HS(real eta, real kf, real k1, real k2);
local real sourceb_HS(real eta, real kf, real k1, real k2);
local real SD2_HS(real eta, real x, real k, real p);
local real S3I_HS(real eta, real x, real k, real p, real Dpk, real Dpp,
                real D2f, real D2mf);
local real S3II_HS(real eta, real x, real k, real p, real Dpk, real Dpp,
                 real D2f, real D2mf);
local real S3FL_HS(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f, real D2mf);
local real S3dI_HS(real eta, real x, real k, real p, real Dpk, real Dpp,
                 real D2f, real D2mf);
// End: HS Model global HEADERS
// ==========================================


// ==========================================
// Begin: fR1 Model global HEADERS -> local
local void set_Model_fR1(void);
// End: fR1 Model global HEADERS
// ==========================================


// ==========================================
// Begin: LCDM Model global HEADERS -> local
local void set_Model_LCDM(void);
local real mu_LCDM(real eta, real k);
local real sourceA_LCDM(real eta, real kf, real k1, real k2);
local real sourceb_LCDM(real eta, real kf, real k1, real k2);
local real SD2_LCDM(real eta, real x, real k, real p);
local real S3I_LCDM(real eta, real x, real k, real p, real Dpk, real Dpp,
                  real D2f, real D2mf);
local real S3II_LCDM(real eta, real x, real k, real p, real Dpk, real Dpp,
                   real D2f, real D2mf);
local real S3FL_LCDM(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f, real D2mf);
local real S3dI_LCDM(real eta, real x, real k, real p, real Dpk, real Dpp,
                   real D2f, real D2mf);
// End: LCDM Model global HEADERS
// ==========================================



local void model_string_to_int(string, int *);

#define HS                          0
#define fR1                         1

global void set_model(void)
{
    int model_int;
    
    model_string_to_int(cmd.mgmodel, &model_int);
    model_int_flag = model_int;
    switch (model_int){
        case HS: set_Model_HS(); break;
        case fR1: set_Model_fR1(); break;
        case LCDM: set_Model_LCDM(); break;
        default: error("\nUnknown model type %s\n\n",cmd.mgmodel);
    }
}

local void model_string_to_int(string model_str,int *model_int)
{
    *model_int = -1;
    if (strcmp(model_str,"HS") == 0)                *model_int=HS;
    if (strcmp(model_str,"fR1") == 0)               *model_int=fR1;
    if (strcmp(model_str,"LCDM") == 0)              *model_int=LCDM;
}

global real kpp(real x, real k, real p)
{
    real kpptmp;
    
    kpptmp = rsqrt(rsqr(k) + rsqr(p) + 2.0*k*p*x);
    
    return kpptmp;
}

// ===========================================================================
// MODELS SECTION 2
// ===========================================================================

// ===========================================================================
// BEGIN :: SWITCHING MODELS ROUTINES
// ===========================================================================

global real mu(real eta, real k)
{
    real tmp;
    switch (model_int_flag){
        case HS: tmp=mu_HS(eta, k); break;
        case LCDM: tmp=mu_LCDM(eta, k); break;
    }
    return tmp;
}

global real sourceA(real eta, real kf, real k1, real k2)
{
    real tmp;
    switch (model_int_flag){
        case HS: tmp=sourceA_HS(eta, kf, k1, k2); break;
        case LCDM: tmp=sourceA_LCDM(eta, kf, k1, k2); break;
    }
    return tmp;
}

global real sourceb(real eta, real kf, real k1, real k2)
{
    real tmp;
    switch (model_int_flag){
        case HS: tmp=sourceb_HS(eta, kf, k1, k2); break;
        case LCDM: tmp=sourceb_LCDM(eta, kf, k1, k2); break;
    }
    return tmp;
}

global real SD2(real eta, real x, real k, real p)
{
    real tmp;
    switch (model_int_flag){
        case HS: tmp=SD2_HS(eta, x, k, p); break;
        case LCDM: tmp=SD2_LCDM(eta, x, k, p); break;
    }
    return tmp;
}

global real S3I(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f, real D2mf)
{
    real tmp;
    switch (model_int_flag){
        case HS: tmp=S3I_HS(eta, x, k, p, Dpk, Dpp, D2f, D2mf); break;
        case LCDM: tmp=S3I_LCDM(eta, x, k, p, Dpk, Dpp, D2f, D2mf); break;
    }
    return tmp;
}

global real S3II(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f, real D2mf)
{
    real tmp;
    switch (model_int_flag){
        case HS: tmp=S3II_HS(eta, x, k, p, Dpk, Dpp, D2f, D2mf); break;
        case LCDM: tmp=S3II_LCDM(eta, x, k, p, Dpk, Dpp, D2f, D2mf); break;
    }
    return tmp;
}

global real S3FL(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f, real D2mf)
{
    real tmp;
    switch (model_int_flag){
        case HS: tmp=S3FL_HS(eta, x, k, p, Dpk, Dpp, D2f, D2mf); break;
        case LCDM: tmp=S3FL_LCDM(eta, x, k, p, Dpk, Dpp, D2f, D2mf); break;
    }
    return tmp;
}

global real S3dI(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f, real D2mf)
{
    real tmp;
    switch (model_int_flag){
        case HS: tmp=S3dI_HS(eta, x, k, p, Dpk, Dpp, D2f, D2mf); break;
        case LCDM: tmp=S3dI_LCDM(eta, x, k, p, Dpk, Dpp, D2f, D2mf); break;
    }
    return tmp;
}

#undef HS
#undef fR1

// ===========================================================================
// END :: SWITCHING MODELS ROUTINES
// ===========================================================================

// ===========================================================================
// MODELS SECTION 3
// ===========================================================================

// ===========================================================================
// Begin: Hu-Sawicky Model
// ===========================================================================

// Begin: Hu-Sawicky Model local HEADERS
local real mass_HS(real eta);
local real JFL_HS(real eta, real x, real k, real p);
local real KFL_HS(real eta, real k, real k1, real k2);
local real KFL2_HS(real eta, real x, real k, real p);
local real PiF_HS(real eta, real k);
local real M1_HS(real eta);
local real M2_HS(real eta);
local real M3_HS(real eta);
local real S2a_HS(real eta, real x, real k, real p);
local real S2b_HS(real eta, real x, real k, real p);
local real S2FL_HS(real eta, real x, real k, real p);
local real S2dI_HS(real eta, real x, real k, real p);

local real sourcea_HS(real eta, real kf);
local real sourceFL_HS(real eta, real kf, real k1, real k2);
local real sourcedI_HS(real eta, real kf, real k1, real k2);

local real S3IIplus_HS(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f);
local real S3IIminus_HS(real eta, real x, real k, real p, real Dpk, real Dpp, real D2mf);
local real S3FLplus_HS(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f);
local real S3FLminus_HS(real eta, real x, real k, real p, real Dpk, real Dpp, real D2mf);
local real D2phiplus_HS(real eta, real x, real k, real p,
                     real Dpk, real Dpp, real D2f);
local real D2phiminus_HS(real eta, real x, real k, real p,
                      real Dpk, real Dpp, real D2mf);
local real K3dI_HS(real eta, real x, real k,  real p,
                 real Dpk, real Dpp, real D2f, real D2mf);
// End: Hu-Sawicky Model local HEADERS

local void set_Model_HS(void)
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

local real mass_HS(real eta)
{
    real masstmp;
    
    masstmp = (1.0/H02)*rsqrt(1.0/(2.0*cmd.fR0))
                * rpow(cmd.om*rexp(-3.0*eta) + 4.0*(1.0-cmd.om),(2.0+cmd.nHS)/2.0)
                / rpow(cmd.om+4.0*(1.0-cmd.om),(1.0+cmd.nHS)/2.0);

    return (masstmp);
}

local real mu_HS(real eta, real k)
{
    real mutmp;
    mutmp = 1.0 + (2.0*gd.beta2*k*k)/(k*k + rexp(2.0*eta)*rsqr(mass_HS(eta)));

    return (mutmp);
}

local real PiF_HS(real eta, real k)
{
    real PiFtmp;
    
    PiFtmp = k*k/rexp(2.0*eta) + rsqr(mass_HS(eta));
    
    return (PiFtmp);
}


// ---------------------------------------------------------------------------
// BEGIN :: SECOND ORDER (six second order differential equations)
//

local real M2_HS(real eta)
{
    real M2tmp;

    M2tmp = cmd.screening;
    M2tmp *= (9.0/(4.0*rsqr(H02)))*rsqr(1.0/cmd.fR0)
    * rpow(cmd.om*rexp(-3.0*eta)+4.0*(1.0-cmd.om),5.0)
    / rpow(cmd.om+4*(1.0-cmd.om),4.0);
    
    return (M2tmp);
}


local real KFL_HS(real eta, real k, real k1, real k2)
{
    real KFLtmp;

    KFLtmp = 0.5*(rsqr(sqr(k)-sqr(k1)-sqr(k2))/(sqr(k1)*sqr(k2)))*(mu_HS(eta,k1)+mu_HS(eta,k2)-2.0)
        + 0.5*((sqr(k)-sqr(k1)-sqr(k2))/sqr(k1))*(mu_HS(eta,k1)-1.0)
        + 0.5*((sqr(k)-sqr(k1)-sqr(k2))/sqr(k2))*(mu_HS(eta,k2)-1.0);

    return (KFLtmp);
}

local real sourceA_HS(real eta, real kf, real k1, real k2)
{
    real Stmp;

    Stmp = sourcea_HS(eta, kf)
            + sourceFL_HS(eta, kf, k1, k2)
            - sourcedI_HS(eta, kf, k1, k2);
//    - cmd.screening * sourcedI_HS(eta, kf, k1, k2);

    return Stmp;
}

local real sourcea_HS(real eta, real kf)
{
    real Stmp;
    
    Stmp = f1(eta)*mu_HS(eta, kf);

    return Stmp;
}

local real sourceb_HS(real eta, real kf, real k1, real k2)
{
    real Stmp;
    
    Stmp = f1(eta)*( mu_HS(eta, k1) + mu_HS(eta, k2) - mu_HS(eta, kf));

    return Stmp;
}

local real sourceFL_HS(real eta, real kf, real k1, real k2)
{
    real Stmp;
    
    Stmp = f1(eta)*(rsqr(mass_HS(eta))/PiF_HS(eta,kf))*KFL_HS(eta, kf, k1, k2);

    return Stmp;
}

local real sourcedI_HS(real eta, real kf, real k1, real k2)
{
    real Stmp;
    
    Stmp = (1.0/6.0)*rsqr(OmM(eta)*H(eta)/(rexp(eta)*H02))
            * rsqr(kf)*M2_HS(eta)/ ( PiF_HS(eta,kf)*PiF_HS(eta,k1)*PiF_HS(eta,k2) );
    
    return Stmp;
}

//
// END :: SECOND ORDER (six second order differential equations)
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
// BEGIN :: THIRD ORDER (Dsymmetric, five second order differential equations)
//

local real M1_HS(real eta)
{
    real M1tmp;
    
    M1tmp = 3.0*rsqr(mass_HS(eta));
    
    return (M1tmp);
}

local real M3_HS(real eta)
{
    real M3tmp;

    M2tmp = cmd.screening;
    M3tmp *= (45.0/(8.0*rsqr(H02)))*rpow(1.0/cmd.fR0,3.0)
    * rpow(cmd.om*rexp(-3.0*eta)+4.0*(1.0-cmd.om),7.0)
    / rpow(cmd.om+4*(1.0-cmd.om),6.0);
    
    return (M3tmp);
}

local real KFL2_HS(real eta, real x, real k, real p)
{
    real KFLtmp;
    
    KFLtmp = 2.0*rsqr(x)*( mu_HS(eta,k) + mu_HS(eta,p)-2.0 )
            + (p*x/k)*(mu_HS(eta,k) - 1.0)
            + (k*x/p)*(mu_HS(eta,p) - 1.0);

    return (KFLtmp);
}

local real JFL_HS(real eta, real x, real k, real p)
{
    real JFLtmp;
    
    JFLtmp = (9.0/(2.0*A0(eta)))
            * KFL2_HS(eta, x, k, p) * PiF_HS(eta, k) * PiF_HS(eta, p);

    return (JFLtmp);
}

local real D2phiplus_HS(real eta, real x, real k, real p,
                     real Dpk, real Dpp, real D2f)
{
    real D2tmp;
    
    D2tmp = (
             (1.0 + rsqr(x))
             -( 2.0*A0(eta)/3.0 )
             * (
                ( M2_HS(eta) + JFL_HS(eta,x,k,p)*(3.0+2.0*cmd.omegaBD) )
                / (3.0*PiF_HS(eta,k)*PiF_HS(eta,p))
                )
             ) * Dpk*Dpp + D2f;
    
    return (D2tmp);
}

local real D2phiminus_HS(real eta, real x, real k, real p,
                     real Dpk, real Dpp, real D2mf)
{
    real D2tmp;
    
    D2tmp = (
             (1.0 + rsqr(x))
             -( 2.0*A0(eta)/3.0 )
             * (
                ( M2_HS(eta) + JFL_HS(eta,-x,k,p)*(3.0+2.0*cmd.omegaBD) )
                / (3.0*PiF_HS(eta,k)*PiF_HS(eta,p))
                )
             ) * Dpk*Dpp + D2mf;

    return (D2tmp);
}

local real K3dI_HS(real eta, real x, real k,  real p,
                 real Dpk, real Dpp, real D2f, real D2mf)
{
    real K3tmp, kplusp, kpluspm;
    real t1, t2, t3, t4, t5, t6;

    kplusp = kpp(x,k,p);
    kpluspm = kpp(-x,k,p);

    t1 = 2.0*rsqr(OmM(eta)*H(eta)/H02)
            *(M2_HS(eta)/(PiF_HS(eta,k)*PiF_HS(eta,0)));

    t2 = (1.0/3.0)*(rpow(OmM(eta),3.0)*rpow(H(eta),4.0)/rpow(H02,4) )
        *(
            M3_HS(eta) - M2_HS(eta)*(M2_HS(eta) + JFL_HS(eta,-1.0,p,p)*(3.0+2.0*cmd.omegaBD))
                        /(PiF_HS(eta,0))
          ) / ( rsqr(PiF_HS(eta,p)) * PiF_HS(eta,k) );

    t3 = rsqr(OmM(eta)*H(eta)/H02)
        *(M2_HS(eta)/(PiF_HS(eta,p)*PiF_HS(eta,kplusp)))
        *(
            1.0 + rsqr(x) + (D2f)/(Dpk*Dpp)
          );
    
    t4 = (1.0/3.0)*(rpow(OmM(eta),3.0)*rpow(H(eta),4.0)/rpow(H02,4) )
        *(
            M3_HS(eta) - M2_HS(eta)*(M2_HS(eta) + JFL_HS(eta,x,k,p)*(3.0+2.0*cmd.omegaBD))
                        /(PiF_HS(eta,kplusp))
          ) / ( rsqr(PiF_HS(eta,p)) * PiF_HS(eta,k) );
    
    t5 = rsqr(OmM(eta)*H(eta)/H02)
        *(M2_HS(eta)/(PiF_HS(eta,p)*PiF_HS(eta,kpluspm)))
        *(
            1.0 + rsqr(x) + (D2mf)/(Dpk*Dpp)
          );
    
    t6 = (1.0/3.0)*(rpow(OmM(eta),3.0)*rpow(H(eta),4.0)/rpow(H02,4) )
        *(
          M3_HS(eta) - M2_HS(eta)*(M2_HS(eta) + JFL_HS(eta,-x,k,p)*(3.0+2.0*cmd.omegaBD))
                /(PiF_HS(eta,kpluspm))
          ) / ( rsqr(PiF_HS(eta,p)) * PiF_HS(eta,k) );

    K3tmp = t1 + t2 + t3 + t4 + t5 + t6;

    return (K3tmp);
}

local real S2a_HS(real eta, real x, real k, real p)
{
    real Dtmp, kplusp;

    kplusp = kpp(x,k,p);
    Dtmp = f1(eta)*mu_HS(eta,kplusp);

    return Dtmp;
}

local real S2b_HS(real eta, real x, real k, real p)
{
    real Dtmp, kplusp;
    
    kplusp = kpp(x,k,p);
    Dtmp = f1(eta)*(mu_HS(eta,k)+mu_HS(eta,p)-mu_HS(eta,kplusp));
    
    return Dtmp;
}

local real S2FL_HS(real eta, real x, real k, real p)
{
    real Dtmp, kplusp;
    
    kplusp = kpp(x,k,p);
    Dtmp = f1(eta)*(
                    M1_HS(eta)/(3.0*PiF_HS(eta,kplusp))
                    *KFL2_HS(eta,x,k,p)
                    );
    return Dtmp;
}

local real S2dI_HS(real eta, real x, real k, real p)
{
    real Dtmp, kplusp;
    
    kplusp = kpp(x,k,p);
    Dtmp = (1.0/6.0)*
            rsqr(OmM(eta)*H(eta)/(rexp(eta)*H02))
        * ( (rsqr(kplusp)*M2_HS(eta)) / (PiF_HS(eta,kplusp)*PiF_HS(eta,k)*PiF_HS(eta,p)));

    return Dtmp;
}

local real SD2_HS(real eta, real x, real k, real p)
{
    real Dtmp;

    Dtmp = S2a_HS(eta, x, k, p) -  S2b_HS(eta, x, k, p)*rsqr(x)
        + S2FL_HS(eta, x, k, p) - S2dI_HS(eta, x, k, p);

    return Dtmp;
}

local real S3I_HS(real eta, real x, real k, real p, real Dpk, real Dpp,
                real D2f, real D2mf)
{
    real Stmp, kplusp, kpluspm;

    kplusp = kpp(x,k,p);
    kpluspm = kpp(-x,k,p);
    Stmp = (
            f1(eta)*(mu_HS(eta,p)+mu_HS(eta,kplusp)-mu_HS(eta,k))*D2f*Dpp
                + SD2_HS(eta,x,k,p)*Dpk*Dpp*Dpp
            )*(1.0 - rsqr(x))/(1.0 + rsqr(p/k) + 2.0*(p/k)*x)
        + (
           f1(eta)*(mu_HS(eta,p)+mu_HS(eta,kpluspm)-mu_HS(eta,k))*D2mf*Dpp
            + SD2_HS(eta,-x,k,p)*Dpk*Dpp*Dpp
           )*(1.0 - rsqr(x))/(1.0 + rsqr(p/k) - 2.0*(p/k)*x);

    return (Stmp);
}

local real S3IIplus_HS(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f)
{
    real Stmp, kplusp;
    
    kplusp = kpp(x,k,p);
    
    Stmp =
    -f1(eta)*(mu_HS(eta,p)+mu_HS(eta,kplusp)-2.0*mu_HS(eta,k))
    * Dpp*( D2f + Dpk*Dpp*rsqr(x) )
    
    -f1(eta)*(mu_HS(eta,kplusp)-mu_HS(eta,k))*Dpk*Dpp*Dpp
    
    -(
      (M1_HS(eta)/(3.0*PiF_HS(eta,kplusp))) * f1(eta)*KFL2_HS(eta,x,k,p)
      -rsqr(OmM(eta)*H(eta)/H02)
      * (M2_HS(eta)*kplusp*kplusp*rexp(-2.0*eta))
      / (6.0*PiF_HS(eta,kplusp)*PiF_HS(eta,k)*PiF_HS(eta,p))
      )*Dpk*Dpp*Dpp;
    
    return (Stmp);
}

local real S3IIminus_HS(real eta, real x, real k, real p, real Dpk, real Dpp, real D2mf)
{
    real Stmp, kpluspm;

    kpluspm = kpp(-x,k,p);

    Stmp =
    -f1(eta)*(mu_HS(eta,p)+mu_HS(eta,kpluspm)-2.0*mu_HS(eta,k))
    * Dpp*( D2mf + Dpk*Dpp*rsqr(x) )
    
    -f1(eta)*(mu_HS(eta,kpluspm)-mu_HS(eta,k))*Dpk*Dpp*Dpp
    
    -(
      (M1_HS(eta)/(3.0*PiF_HS(eta,kpluspm))) * f1(eta)*KFL2_HS(eta,-x,k,p)
      -rsqr(OmM(eta)*H(eta)/H02)
      * (M2_HS(eta)*kpluspm*kpluspm*rexp(-2.0*eta))
      / (6.0*PiF_HS(eta,kpluspm)*PiF_HS(eta,k)*PiF_HS(eta,p))
      )*Dpk*Dpp*Dpp;
    
    return (Stmp);
}

local real S3II_HS(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f, real D2mf)
{
    real Stmp;
    
    Stmp =  S3IIplus_HS(eta, x, k, p, Dpk, Dpp, D2f)
         + S3IIminus_HS(eta, x, k, p, Dpk, Dpp, D2mf);

    return (Stmp);
}

local real S3FLplus_HS(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f)
{
    real Stmp, kplusp;
    
    kplusp = kpp(x,k,p);
    
    Stmp = f1(eta)*(M1_HS(eta)/(3.0*PiF_HS(eta,k)))
    *(
        (2.0*rsqr(p+k*x)/rsqr(kplusp) - 1.0 - (k*x)/p )
        *( mu_HS(eta,p)-1.0 )* D2f * Dpp
      
        + ( (rsqr(p) + 3.0*k*p*x + 2.0*k*k * x*x)/rsqr(kplusp) )
        *( mu_HS(eta,kplusp) - 1.0) * D2phiplus_HS(eta,x,k,p,Dpk,Dpp,D2f) * Dpp
      
      + 3.0*rsqr(x)*( mu_HS(eta,k) + mu_HS(eta,p) - 2.0 ) * Dpk * Dpp * Dpp
    );
    
    return (Stmp);
}

local real S3FLminus_HS(real eta, real x, real k, real p, real Dpk, real Dpp, real D2mf)
{
    real Stmp, kpluspm;
    
    kpluspm = kpp(-x,k,p);
    
    Stmp = f1(eta)*(M1_HS(eta)/(3.0*PiF_HS(eta,k)))
    *(
      (2.0*rsqr(p-k*x)/rsqr(kpluspm) - 1.0 + (k*x)/p )
      *( mu_HS(eta,p)-1.0 )* D2mf * Dpp
      
      + ( (rsqr(p) - 3.0*k*p*x + 2.0*k*k * x*x)/rsqr(kpluspm) )
      *( mu_HS(eta,kpluspm) - 1.0) * D2phiminus_HS(eta,x,k,p,Dpk,Dpp,D2mf) * Dpp
      
      + 3.0*rsqr(x)*( mu_HS(eta,k) + mu_HS(eta,p) - 2.0 ) * Dpk * Dpp * Dpp
      );
    
    return (Stmp);
}

local real S3FL_HS(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f, real D2mf)
{
    real Stmp;

    Stmp = S3FLplus_HS(eta, x, k, p, Dpk, Dpp, D2f)
        + S3FLminus_HS(eta, x, k, p, Dpk, Dpp, D2mf);

    return (Stmp);
}

local real S3dI_HS(real eta, real x, real k, real p, real Dpk, real Dpp,
                 real D2f, real D2mf)
{
    real Stmp;
    
    Stmp = -(rsqr(k)/rexp(2.0*eta))
        *(1.0/(6.0*PiF_HS(eta,k)))
        *K3dI_HS(eta,x,k,p,Dpk,Dpp,D2f,D2mf)*Dpk*Dpp*Dpp;

    return (Stmp);
}

//
// END :: THIRD ORDER (Dsymmetric, five second order differential equations)
// ---------------------------------------------------------------------------


// ===========================================================================
// End: Hu-Sawicky Model
// ===========================================================================



// ===========================================================================
// Begin: fR1 Model
// ===========================================================================

local void set_Model_fR1(void)
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

// ===========================================================================
// End: fR1 Model
// ===========================================================================


// ===========================================================================
// Begin: LCDM Model
// ===========================================================================

// Begin: LCDM Model local HEADERS
local real mass_LCDM(real eta);
local real JFL_LCDM(real eta, real x, real k, real p);
local real KFL_LCDM(real eta, real k, real k1, real k2);
local real KFL2_LCDM(real eta, real x, real k, real p);
local real PiF_LCDM(real eta, real k);
local real M1_LCDM(real eta);
local real M2_LCDM(real eta);
local real M3_LCDM(real eta);
local real S2a_LCDM(real eta, real x, real k, real p);
local real S2b_LCDM(real eta, real x, real k, real p);
local real S2FL_LCDM(real eta, real x, real k, real p);
local real S2dI_LCDM(real eta, real x, real k, real p);

local real sourcea_LCDM(real eta, real kf);
local real sourceFL_LCDM(real eta, real kf, real k1, real k2);
local real sourcedI_LCDM(real eta, real kf, real k1, real k2);

local real S3IIplus_LCDM(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f);
local real S3IIminus_LCDM(real eta, real x, real k, real p, real Dpk, real Dpp, real D2mf);
local real S3FLplus_LCDM(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f);
local real S3FLminus_LCDM(real eta, real x, real k, real p, real Dpk, real Dpp, real D2mf);
local real D2phiplus_LCDM(real eta, real x, real k, real p,
                        real Dpk, real Dpp, real D2f);
local real D2phiminus_LCDM(real eta, real x, real k, real p,
                         real Dpk, real Dpp, real D2mf);
local real K3dI_LCDM(real eta, real x, real k,  real p,
                   real Dpk, real Dpp, real D2f, real D2mf);
// End: LCDM Model local HEADERS

local void set_Model_LCDM(void)
{
    strcpy(gd.model_comment, "LCDM Model");
    //
    // Parameters set and default values:
    //    nHS=1.0;
    //    fR0=1.0e-5;
    //    beta2=1.0/6.0;
    //    omegaBD = 0.0;
    //    screening = 1.0;
    gd.beta2 = 0.0;
    cmd.omegaBD = 0.0;
}

local real mass_LCDM(real eta)
{
    real masstmp;
    
    masstmp = (1.0/H02)*rsqrt(1.0/(2.0*cmd.fR0))
    * rpow(cmd.om*rexp(-3.0*eta) + 4.0*(1.0-cmd.om),(2.0+cmd.nHS)/2.0)
    / rpow(cmd.om+4.0*(1.0-cmd.om),(1.0+cmd.nHS)/2.0);
    
    return (masstmp);
}

local real mu_LCDM(real eta, real k)
{
    real mutmp;
    mutmp = 1.0;
//    + (2.0*gd.beta2*k*k)/(k*k + rexp(2.0*eta)*rsqr(mass_LCDM(eta)));
    
    return (mutmp);
}

local real PiF_LCDM(real eta, real k)
{
    real PiFtmp;
    
    PiFtmp = k*k/rexp(2.0*eta) + rsqr(mass_LCDM(eta));
    
    return (PiFtmp);
}


// ---------------------------------------------------------------------------
// BEGIN :: SECOND ORDER (six second order differential equations)
//

local real M2_LCDM(real eta)
{
    real M2tmp;
    
    M2tmp = (9.0/(4.0*rsqr(H02)))*rsqr(1.0/cmd.fR0)
    * rpow(cmd.om*rexp(-3.0*eta)+4.0*(1.0-cmd.om),5.0)
    / rpow(cmd.om+4*(1.0-cmd.om),4.0);
    
    return (M2tmp);
}


local real KFL_LCDM(real eta, real k, real k1, real k2)
{
    real KFLtmp;
    
    KFLtmp = 0.5*(rsqr(sqr(k)-sqr(k1)-sqr(k2))/(sqr(k1)*sqr(k2)))*(mu_LCDM(eta,k1)+mu_LCDM(eta,k2)-2.0)
    + 0.5*((sqr(k)-sqr(k1)-sqr(k2))/sqr(k1))*(mu_LCDM(eta,k1)-1.0)
    + 0.5*((sqr(k)-sqr(k1)-sqr(k2))/sqr(k2))*(mu_LCDM(eta,k2)-1.0);
    
    return (KFLtmp);
}

local real sourceA_LCDM(real eta, real kf, real k1, real k2)
{
    real Stmp;
    
    Stmp = sourcea_LCDM(eta, kf)
            + sourceFL_LCDM(eta, kf, k1, k2)
            - sourcedI_LCDM(eta, kf, k1, k2);
//    - cmd.screening * sourcedI_LCDM(eta, kf, k1, k2);

    return Stmp;
}

local real sourcea_LCDM(real eta, real kf)
{
    real Stmp;
    
    Stmp = f1(eta)*mu_LCDM(eta, kf);
    
    return Stmp;
}

local real sourceb_LCDM(real eta, real kf, real k1, real k2)
{
    real Stmp;
    
    Stmp = f1(eta)*( mu_LCDM(eta, k1) + mu_LCDM(eta, k2) - mu_LCDM(eta, kf));
    
    return Stmp;
}

local real sourceFL_LCDM(real eta, real kf, real k1, real k2)
{
    real Stmp;
    
    Stmp = f1(eta)*(rsqr(mass_LCDM(eta))/PiF_LCDM(eta,kf))*KFL_LCDM(eta, kf, k1, k2);
    
    return Stmp;
}

local real sourcedI_LCDM(real eta, real kf, real k1, real k2)
{
    real Stmp;
    
    Stmp = (1.0/6.0)*rsqr(OmM(eta)*H(eta)/(rexp(eta)*H02))
    * rsqr(kf)*M2_LCDM(eta)/ ( PiF_LCDM(eta,kf)*PiF_LCDM(eta,k1)*PiF_LCDM(eta,k2) );
    
    return Stmp;
}

//
// END :: SECOND ORDER (six second order differential equations)
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
// BEGIN :: THIRD ORDER (Dsymmetric, five second order differential equations)
//

local real M1_LCDM(real eta)
{
    real M1tmp;
    
    M1tmp = 3.0*rsqr(mass_LCDM(eta));
    
    return (M1tmp);
}

local real M3_LCDM(real eta)
{
    real M3tmp;
    
    M3tmp = (45.0/(8.0*rsqr(H02)))*rpow(1.0/cmd.fR0,3.0)
    * rpow(cmd.om*rexp(-3.0*eta)+4.0*(1.0-cmd.om),7.0)
    / rpow(cmd.om+4*(1.0-cmd.om),6.0);
    
    return (M3tmp);
}

local real KFL2_LCDM(real eta, real x, real k, real p)
{
    real KFLtmp;
    
    KFLtmp = 2.0*rsqr(x)*( mu_LCDM(eta,k) + mu_LCDM(eta,p)-2.0 )
    + (p*x/k)*(mu_LCDM(eta,k) - 1.0)
    + (k*x/p)*(mu_LCDM(eta,p) - 1.0);
    
    return (KFLtmp);
}

local real JFL_LCDM(real eta, real x, real k, real p)
{
    real JFLtmp;
    
    JFLtmp = (9.0/(2.0*A0(eta)))
    * KFL2_LCDM(eta, x, k, p) * PiF_LCDM(eta, k) * PiF_LCDM(eta, p);
    
    return (JFLtmp);
}

local real D2phiplus_LCDM(real eta, real x, real k, real p,
                        real Dpk, real Dpp, real D2f)
{
    real D2tmp;
    
    D2tmp = (
             (1.0 + rsqr(x))
             -( 2.0*A0(eta)/3.0 )
             * (
                ( M2_LCDM(eta) + JFL_LCDM(eta,x,k,p)*(3.0+2.0*cmd.omegaBD) )
                / (3.0*PiF_LCDM(eta,k)*PiF_LCDM(eta,p))
                )
             ) * Dpk*Dpp + D2f;
    
    return (D2tmp);
}

local real D2phiminus_LCDM(real eta, real x, real k, real p,
                         real Dpk, real Dpp, real D2mf)
{
    real D2tmp;
    
    D2tmp = (
             (1.0 + rsqr(x))
             -( 2.0*A0(eta)/3.0 )
             * (
                ( M2_LCDM(eta) + JFL_LCDM(eta,-x,k,p)*(3.0+2.0*cmd.omegaBD) )
                / (3.0*PiF_LCDM(eta,k)*PiF_LCDM(eta,p))
                )
             ) * Dpk*Dpp + D2mf;
    
    return (D2tmp);
}

local real K3dI_LCDM(real eta, real x, real k,  real p,
                   real Dpk, real Dpp, real D2f, real D2mf)
{
    real K3tmp, kplusp, kpluspm;
    real t1, t2, t3, t4, t5, t6;
    
    kplusp = kpp(x,k,p);
    kpluspm = kpp(-x,k,p);
    
    t1 = 2.0*rsqr(OmM(eta)*H(eta)/H02)
    *(M2_LCDM(eta)/(PiF_LCDM(eta,k)*PiF_LCDM(eta,0)));
    
    t2 = (1.0/3.0)*(rpow(OmM(eta),3.0)*rpow(H(eta),4.0)/rpow(H02,4) )
    *(
      M3_LCDM(eta) - M2_LCDM(eta)*(M2_LCDM(eta) + JFL_LCDM(eta,-1.0,p,p)*(3.0+2.0*cmd.omegaBD))
      /(PiF_LCDM(eta,0))
      ) / ( rsqr(PiF_LCDM(eta,p)) * PiF_LCDM(eta,k) );
    
    t3 = rsqr(OmM(eta)*H(eta)/H02)
    *(M2_LCDM(eta)/(PiF_LCDM(eta,p)*PiF_LCDM(eta,kplusp)))
    *(
      1.0 + rsqr(x) + (D2f)/(Dpk*Dpp)
      );
    
    t4 = (1.0/3.0)*(rpow(OmM(eta),3.0)*rpow(H(eta),4.0)/rpow(H02,4) )
    *(
      M3_LCDM(eta) - M2_LCDM(eta)*(M2_LCDM(eta) + JFL_LCDM(eta,x,k,p)*(3.0+2.0*cmd.omegaBD))
      /(PiF_LCDM(eta,kplusp))
      ) / ( rsqr(PiF_LCDM(eta,p)) * PiF_LCDM(eta,k) );
    
    t5 = rsqr(OmM(eta)*H(eta)/H02)
    *(M2_LCDM(eta)/(PiF_LCDM(eta,p)*PiF_LCDM(eta,kpluspm)))
    *(
      1.0 + rsqr(x) + (D2mf)/(Dpk*Dpp)
      );
    
    t6 = (1.0/3.0)*(rpow(OmM(eta),3.0)*rpow(H(eta),4.0)/rpow(H02,4) )
    *(
      M3_LCDM(eta) - M2_LCDM(eta)*(M2_LCDM(eta) + JFL_LCDM(eta,-x,k,p)*(3.0+2.0*cmd.omegaBD))
      /(PiF_LCDM(eta,kpluspm))
      ) / ( rsqr(PiF_LCDM(eta,p)) * PiF_LCDM(eta,k) );
    
    K3tmp = t1 + t2 + t3 + t4 + t5 + t6;
    
    return (K3tmp);
}

local real S2a_LCDM(real eta, real x, real k, real p)
{
    real Dtmp, kplusp;
    
    kplusp = kpp(x,k,p);
    Dtmp = f1(eta)*mu_LCDM(eta,kplusp);
    
    return Dtmp;
}

local real S2b_LCDM(real eta, real x, real k, real p)
{
    real Dtmp, kplusp;
    
    kplusp = kpp(x,k,p);
    Dtmp = f1(eta)*(mu_LCDM(eta,k)+mu_LCDM(eta,p)-mu_LCDM(eta,kplusp));
    
    return Dtmp;
}

local real S2FL_LCDM(real eta, real x, real k, real p)
{
    real Dtmp, kplusp;
    
    kplusp = kpp(x,k,p);
    Dtmp = f1(eta)*(
                    M1_LCDM(eta)/(3.0*PiF_LCDM(eta,kplusp))
                    *KFL2_LCDM(eta,x,k,p)
                    );
    return Dtmp;
}

local real S2dI_LCDM(real eta, real x, real k, real p)
{
    real Dtmp, kplusp;
    
    kplusp = kpp(x,k,p);
    Dtmp = (1.0/6.0)*
    rsqr(OmM(eta)*H(eta)/(rexp(eta)*H02))
    * ( (rsqr(kplusp)*M2_LCDM(eta)) / (PiF_LCDM(eta,kplusp)*PiF_LCDM(eta,k)*PiF_LCDM(eta,p)));
    
    return Dtmp;
}

local real SD2_LCDM(real eta, real x, real k, real p)
{
    real Dtmp;
    
    Dtmp = S2a_LCDM(eta, x, k, p) -  S2b_LCDM(eta, x, k, p)*rsqr(x)
    + S2FL_LCDM(eta, x, k, p) - S2dI_LCDM(eta, x, k, p);
    
    return Dtmp;
}

local real S3I_LCDM(real eta, real x, real k, real p, real Dpk, real Dpp,
                  real D2f, real D2mf)
{
    real Stmp, kplusp, kpluspm;
    
    kplusp = kpp(x,k,p);
    kpluspm = kpp(-x,k,p);
    Stmp = (
            f1(eta)*(mu_LCDM(eta,p)+mu_LCDM(eta,kplusp)-mu_LCDM(eta,k))*D2f*Dpp
            + SD2_LCDM(eta,x,k,p)*Dpk*Dpp*Dpp
            )*(1.0 - rsqr(x))/(1.0 + rsqr(p/k) + 2.0*(p/k)*x)
    + (
       f1(eta)*(mu_LCDM(eta,p)+mu_LCDM(eta,kpluspm)-mu_LCDM(eta,k))*D2mf*Dpp
       + SD2_LCDM(eta,-x,k,p)*Dpk*Dpp*Dpp
       )*(1.0 - rsqr(x))/(1.0 + rsqr(p/k) - 2.0*(p/k)*x);
    
    return (Stmp);
}

local real S3IIplus_LCDM(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f)
{
    real Stmp, kplusp;
    
    kplusp = kpp(x,k,p);
    
    Stmp =
    -f1(eta)*(mu_LCDM(eta,p)+mu_LCDM(eta,kplusp)-2.0*mu_LCDM(eta,k))
    * Dpp*( D2f + Dpk*Dpp*rsqr(x) )
    
    -f1(eta)*(mu_LCDM(eta,kplusp)-mu_LCDM(eta,k))*Dpk*Dpp*Dpp
    
    -(
      (M1_LCDM(eta)/(3.0*PiF_LCDM(eta,kplusp))) * f1(eta)*KFL2_LCDM(eta,x,k,p)
      -rsqr(OmM(eta)*H(eta)/H02)
      * (M2_LCDM(eta)*kplusp*kplusp*rexp(-2.0*eta))
      / (6.0*PiF_LCDM(eta,kplusp)*PiF_LCDM(eta,k)*PiF_LCDM(eta,p))
      )*Dpk*Dpp*Dpp;
    
    return (Stmp);
}

local real S3IIminus_LCDM(real eta, real x, real k, real p, real Dpk, real Dpp, real D2mf)
{
    real Stmp, kpluspm;
    
    kpluspm = kpp(-x,k,p);
    
    Stmp =
    -f1(eta)*(mu_LCDM(eta,p)+mu_LCDM(eta,kpluspm)-2.0*mu_LCDM(eta,k))
    * Dpp*( D2mf + Dpk*Dpp*rsqr(x) )
    
    -f1(eta)*(mu_LCDM(eta,kpluspm)-mu_LCDM(eta,k))*Dpk*Dpp*Dpp
    
    -(
      (M1_LCDM(eta)/(3.0*PiF_LCDM(eta,kpluspm))) * f1(eta)*KFL2_LCDM(eta,-x,k,p)
      -rsqr(OmM(eta)*H(eta)/H02)
      * (M2_LCDM(eta)*kpluspm*kpluspm*rexp(-2.0*eta))
      / (6.0*PiF_LCDM(eta,kpluspm)*PiF_LCDM(eta,k)*PiF_LCDM(eta,p))
      )*Dpk*Dpp*Dpp;
    
    return (Stmp);
}

local real S3II_LCDM(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f, real D2mf)
{
    real Stmp;
    
    Stmp =  S3IIplus_LCDM(eta, x, k, p, Dpk, Dpp, D2f)
    + S3IIminus_LCDM(eta, x, k, p, Dpk, Dpp, D2mf);
    
    return (Stmp);
}

local real S3FLplus_LCDM(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f)
{
    real Stmp, kplusp;
    
    kplusp = kpp(x,k,p);
    
    Stmp = f1(eta)*(M1_LCDM(eta)/(3.0*PiF_LCDM(eta,k)))
    *(
      (2.0*rsqr(p+k*x)/rsqr(kplusp) - 1.0 - (k*x)/p )
      *( mu_LCDM(eta,p)-1.0 )* D2f * Dpp
      
      + ( (rsqr(p) + 3.0*k*p*x + 2.0*k*k * x*x)/rsqr(kplusp) )
      *( mu_LCDM(eta,kplusp) - 1.0) * D2phiplus_LCDM(eta,x,k,p,Dpk,Dpp,D2f) * Dpp
      
      + 3.0*rsqr(x)*( mu_LCDM(eta,k) + mu_LCDM(eta,p) - 2.0 ) * Dpk * Dpp * Dpp
      );
    
    return (Stmp);
}

local real S3FLminus_LCDM(real eta, real x, real k, real p, real Dpk, real Dpp, real D2mf)
{
    real Stmp, kpluspm;
    
    kpluspm = kpp(-x,k,p);
    
    Stmp = f1(eta)*(M1_LCDM(eta)/(3.0*PiF_LCDM(eta,k)))
    *(
      (2.0*rsqr(p-k*x)/rsqr(kpluspm) - 1.0 + (k*x)/p )
      *( mu_LCDM(eta,p)-1.0 )* D2mf * Dpp
      
      + ( (rsqr(p) - 3.0*k*p*x + 2.0*k*k * x*x)/rsqr(kpluspm) )
      *( mu_LCDM(eta,kpluspm) - 1.0) * D2phiminus_LCDM(eta,x,k,p,Dpk,Dpp,D2mf) * Dpp
      
      + 3.0*rsqr(x)*( mu_LCDM(eta,k) + mu_LCDM(eta,p) - 2.0 ) * Dpk * Dpp * Dpp
      );
    
    return (Stmp);
}

local real S3FL_LCDM(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f, real D2mf)
{
    real Stmp;
    
    Stmp = S3FLplus_LCDM(eta, x, k, p, Dpk, Dpp, D2f)
    + S3FLminus_LCDM(eta, x, k, p, Dpk, Dpp, D2mf);
    
    return (Stmp);
}

local real S3dI_LCDM(real eta, real x, real k, real p, real Dpk, real Dpp,
                   real D2f, real D2mf)
{
    real Stmp;
    
    Stmp = -(rsqr(k)/rexp(2.0*eta))
    *(1.0/(6.0*PiF_LCDM(eta,k)))
    *K3dI_LCDM(eta,x,k,p,Dpk,Dpp,D2f,D2mf)*Dpk*Dpp*Dpp;
    
    return (Stmp);
}

//
// END :: THIRD ORDER (Dsymmetric, five second order differential equations)
// ---------------------------------------------------------------------------

// ===========================================================================
// End: LCDM Model
// ===========================================================================


// ===========================================================================
// Begin: XXX Model

// End: XXX Model
// ===========================================================================

