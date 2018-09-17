/*==============================================================================
	MODULE: models_user.h			[mgpt]
================================================================================*/

// AS IT IS HERE, IS EQUIVALENT TO HS MODEL... results must be the same.

#define USERMODEL 100

// User model parameters to appear in model parameters input file
typedef struct {
    int ipar1;       // nHS

    real par1;       // fR0
    real par2;       // screening

} usrparam_struct, *usrparam_struct_ptr;

local usrparam_struct usrp;

// Local parameters not to appear in model parameter input file
// (will be constant always for this model)

typedef struct {
    real par1;    // beta2
    real par2;    // omegaBD
} local_usrparam, *usrparam_ptr;

local local_usrparam usrpl;

// MACROS :: user parameters
#define nHS_u     usrp.ipar1
#define fR0_u     usrp.par1
#define screening_u   usrp.par2
#define beta2_u       usrpl.par1
#define omegaBD_u     usrpl.par2


local void set_Model_USER(void)
{
    strcpy(gd.model_comment, "USER Model");
//
// Parameters set and default values:
// Local in parameter file parameters:
    usrp.ipar1  =1.0;              //    nHS=1.0;
    usrp.par1   =1.0e-5;           //    fR0=1.0e-5;
    usrp.par2   =1.0;               //    screening = 1.0;
// Local parameters:
    usrpl.par1  =1.0/6.0;           //    beta2=1.0/6.0;
    usrpl.par2  =0.0;               //    omegaBD = 0.0;
}


// ===========================================================================
// MODELS SECTION 1
// ==========================================
// Begin: USER Model global HEADERS -> local
local void set_Model_USER(void);
local real OmM_USER(real eta);
local  real H_USER(real eta);
local  real f1_USER(real eta);
local  real f2_USER(real eta);
local  real A0_USER(real eta);
local real mu_USER(real eta, real k);
local real sourceA_USER(real eta, real kf, real k1, real k2);
local real sourceb_USER(real eta, real kf, real k1, real k2);
local real SD2_USER(real eta, real x, real k, real p);
local real S3I_USER(real eta, real x, real k, real p, real Dpk, real Dpp,
                  real D2f, real D2mf);
local real S3II_USER(real eta, real x, real k, real p, real Dpk, real Dpp,
                   real D2f, real D2mf);
local real S3FL_USER(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f, real D2mf);
local real S3dI_USER(real eta, real x, real k, real p, real Dpk, real Dpp,
                   real D2f, real D2mf);
// End: HS Model global HEADERS
// ==========================================


// ===========================================================================
// MODELS SECTION 2
// ===========================================================================

// NOTHING IN THIS SECTION...


// ===========================================================================
// MODELS SECTION 3
// ===========================================================================

// NOTHING IN THIS SECTION...


// ===========================================================================
// MODELS SECTION 4
// ===========================================================================

// ===========================================================================
// Begin: USER Model
// ===========================================================================

// Begin :: mglpt_fns
local real OmM_USER(real eta)
{
    real OmMtmp;

    OmMtmp=1.0/(1.0 + ( (gd.ol)/cmd.om)*rexp(3.0*eta) );
    return (OmMtmp);
}

local  real H_USER(real eta)
{
    real Htmp;
    
    Htmp=rsqrt(cmd.om*rexp(-3.0*eta)+(gd.ol));
    return (Htmp);
}

local  real f1_USER(real eta)
{
    real f1tmp;
    
    f1tmp=3.0/(2.0*(1.0
                    + ((gd.ol)/cmd.om)*rexp(3.0*eta)
                    )
               );
    return (f1tmp);
}

local  real f2_USER(real eta)
{
    real f2tmp;
    
    f2tmp=2.0 - 3.0/(2.0*(1.0
                          + ((gd.ol)/cmd.om)*rexp(3.0*eta)
                          )
                     );
    return (f2tmp);
}

local real A0_USER(real eta)
{
    real A0tmp;
    
    A0tmp = 1.5*OmM_USER(eta)*rsqr(H_USER(eta))/rsqr(invH0);
    
    return (A0tmp);
}
// End :: mglpt_fns

// Begin: Hu-Sawicky Model local HEADERS
local real mass_USER(real eta);
local real JFL_USER(real eta, real x, real k, real p);
local real KFL_USER(real eta, real k, real k1, real k2);
local real KFL2_USER(real eta, real x, real k, real p);
local real PiF_USER(real eta, real k);
local real M1_USER(real eta);
local real M2_USER(real eta);
local real M3_USER(real eta);
local real S2a_USER(real eta, real x, real k, real p);
local real S2b_USER(real eta, real x, real k, real p);
local real S2FL_USER(real eta, real x, real k, real p);
local real S2dI_USER(real eta, real x, real k, real p);

local real sourcea_USER(real eta, real kf);
local real sourceFL_USER(real eta, real kf, real k1, real k2);
local real sourcedI_USER(real eta, real kf, real k1, real k2);

local real S3IIplus_USER(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f);
local real S3IIminus_USER(real eta, real x, real k, real p, real Dpk, real Dpp, real D2mf);
local real S3FLplus_USER(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f);
local real S3FLminus_USER(real eta, real x, real k, real p, real Dpk, real Dpp, real D2mf);
local real D2phiplus_USER(real eta, real x, real k, real p,
                     real Dpk, real Dpp, real D2f);
local real D2phiminus_USER(real eta, real x, real k, real p,
                      real Dpk, real Dpp, real D2mf);
local real K3dI_USER(real eta, real x, real k,  real p,
                real Dpk, real Dpp, real D2f, real D2mf);
// End: Hu-Sawicky Model local HEADERS

/*
local void set_Model_USER(void)
{
    strcpy(gd.model_comment, "USER Model");
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
*/

local real mass_USER(real eta)
{
    real masstmp;
    
    masstmp = (1.0/invH0)*rsqrt(1.0/(2.0*fR0_u))
                * rpow(cmd.om*rexp(-3.0*eta) + 4.0*(gd.ol),(2.0+nHS_u)/2.0)
                / rpow(cmd.om+4.0*(gd.ol),(1.0+nHS_u)/2.0);

    return (masstmp);
}

local real mu_USER(real eta, real k)
{
    real mutmp;
    mutmp = 1.0 + (2.0*beta2_u*k*k)/(k*k + rexp(2.0*eta)*rsqr(mass_USER(eta)));

    return (mutmp);
}

local real PiF_USER(real eta, real k)
{
    real PiFtmp;
    
    PiFtmp = k*k/rexp(2.0*eta) + rsqr(mass_USER(eta));
    
    return (PiFtmp);
}


// ===========================================================================
// BEGIN :: SECOND ORDER (six second order differential equations)
//

local real M2_USER(real eta)
{
    real M2tmp;
    
    M2tmp = screening_u;
    M2tmp *= (9.0/(4.0*rsqr(invH0)))*rsqr(1.0/fR0_u)
            * rpow(cmd.om*rexp(-3.0*eta)+4.0*(gd.ol),5.0)
            / rpow(cmd.om+4*(gd.ol),4.0);
    
    return (M2tmp);
}

local real KFL_USER(real eta, real k, real k1, real k2)
{
    real KFLtmp;

    KFLtmp = 0.5*(rsqr(sqr(k)-sqr(k1)-sqr(k2))/(sqr(k1)*sqr(k2)))*(mu_USER(eta,k1)+mu_USER(eta,k2)-2.0)
        + 0.5*((sqr(k)-sqr(k1)-sqr(k2))/sqr(k1))*(mu_USER(eta,k1)-1.0)
        + 0.5*((sqr(k)-sqr(k1)-sqr(k2))/sqr(k2))*(mu_USER(eta,k2)-1.0);

    return (KFLtmp);
}

local real sourceA_USER(real eta, real kf, real k1, real k2)
{
    real Stmp;

    Stmp = sourcea_USER(eta, kf)
            + sourceFL_USER(eta, kf, k1, k2)
            - sourcedI_USER(eta, kf, k1, k2);

    return Stmp;
}

local real sourcea_USER(real eta, real kf)
{
    real Stmp;
    
    Stmp = f1(eta)*mu_USER(eta, kf);

    return Stmp;
}

local real sourceb_USER(real eta, real kf, real k1, real k2)
{
    real Stmp;
    
    Stmp = f1(eta)*( mu_USER(eta, k1) + mu_USER(eta, k2) - mu_USER(eta, kf));

    return Stmp;
}

local real sourceFL_USER(real eta, real kf, real k1, real k2)
{
    real Stmp;
    
    Stmp = f1(eta)*(rsqr(mass_USER(eta))/PiF_USER(eta,kf))*KFL_USER(eta, kf, k1, k2);

    return Stmp;
}

local real sourcedI_USER(real eta, real kf, real k1, real k2)
{
    real Stmp;
    
    Stmp = (1.0/6.0)*rsqr(OmM_USER(eta)*H_USER(eta)/(rexp(eta)*invH0))
            * rsqr(kf)*M2_USER(eta)/ ( PiF_USER(eta,kf)*PiF_USER(eta,k1)*PiF_USER(eta,k2) );
    
    return Stmp;
}

//
// END :: SECOND ORDER (six second order differential equations)
// ===========================================================================


// ===========================================================================
// BEGIN :: THIRD ORDER (Dsymmetric, five second order differential equations)
//

local real M1_USER(real eta)
{
    real M1tmp;
    
    M1tmp = 3.0*rsqr(mass_USER(eta));
    
    return (M1tmp);
}

local real M3_USER(real eta)
{
    real M3tmp;

    M3tmp = screening_u;
    M3tmp *= (45.0/(8.0*rsqr(invH0)))*rpow(1.0/fR0_u,3.0)
        * rpow(cmd.om*rexp(-3.0*eta)+4.0*(gd.ol),7.0)
        / rpow(cmd.om+4*(gd.ol),6.0);

    return (M3tmp);
}

local real KFL2_USER(real eta, real x, real k, real p)
{
    real KFLtmp;
    
    KFLtmp = 2.0*rsqr(x)*( mu_USER(eta,k) + mu_USER(eta,p)-2.0 )
            + (p*x/k)*(mu_USER(eta,k) - 1.0)
            + (k*x/p)*(mu_USER(eta,p) - 1.0);

    return (KFLtmp);
}

local real JFL_USER(real eta, real x, real k, real p)
{
    real JFLtmp;
    
    JFLtmp = (9.0/(2.0*A0(eta)))
            * KFL2_USER(eta, x, k, p) * PiF_USER(eta, k) * PiF_USER(eta, p);

    return (JFLtmp);
}

local real D2phiplus_USER(real eta, real x, real k, real p,
                     real Dpk, real Dpp, real D2f)
{
    real D2tmp;
    
    D2tmp = (
             (1.0 + rsqr(x))
             -( 2.0*A0(eta)/3.0 )
             * (
                ( M2_USER(eta) + JFL_USER(eta,x,k,p)*(3.0+2.0*omegaBD_u) )
                / (3.0*PiF_USER(eta,k)*PiF_USER(eta,p))
                )
             ) * Dpk*Dpp + D2f;
    
    return (D2tmp);
}

local real D2phiminus_USER(real eta, real x, real k, real p,
                     real Dpk, real Dpp, real D2mf)
{
    real D2tmp;
    
    D2tmp = (
             (1.0 + rsqr(x))
             -( 2.0*A0(eta)/3.0 )
             * (
                ( M2_USER(eta) + JFL_USER(eta,-x,k,p)*(3.0+2.0*omegaBD_u) )
                / (3.0*PiF_USER(eta,k)*PiF_USER(eta,p))
                )
             ) * Dpk*Dpp + D2mf;

    return (D2tmp);
}

local real K3dI_USER(real eta, real x, real k,  real p,
                 real Dpk, real Dpp, real D2f, real D2mf)
{
    real K3tmp, kplusp, kpluspm;
    real t1, t2, t3, t4, t5, t6;

    kplusp = kpp(x,k,p);
    kpluspm = kpp(-x,k,p);

    t1 = 2.0*rsqr(OmM(eta)*H(eta)/invH0)
            *(M2_USER(eta)/(PiF_USER(eta,k)*PiF_USER(eta,0)));

    t2 = (1.0/3.0)*(rpow(OmM(eta),3.0)*rpow(H(eta),4.0)/rpow(invH0,4) )
        *(
            M3_USER(eta) - M2_USER(eta)*(M2_USER(eta) + JFL_USER(eta,-1.0,p,p)*(3.0+2.0*omegaBD_u))
                        /(PiF_USER(eta,0))
          ) / ( rsqr(PiF_USER(eta,p)) * PiF_USER(eta,k) );

    t3 = rsqr(OmM(eta)*H(eta)/invH0)
        *(M2_USER(eta)/(PiF_USER(eta,p)*PiF_USER(eta,kplusp)))
        *(
            1.0 + rsqr(x) + (D2f)/(Dpk*Dpp)
          );
    
    t4 = (1.0/3.0)*(rpow(OmM(eta),3.0)*rpow(H(eta),4.0)/rpow(invH0,4) )
        *(
            M3_USER(eta) - M2_USER(eta)*(M2_USER(eta) + JFL_USER(eta,x,k,p)*(3.0+2.0*omegaBD_u))
                        /(PiF_USER(eta,kplusp))
          ) / ( rsqr(PiF_USER(eta,p)) * PiF_USER(eta,k) );
    
    t5 = rsqr(OmM(eta)*H(eta)/invH0)
        *(M2_USER(eta)/(PiF_USER(eta,p)*PiF_USER(eta,kpluspm)))
        *(
            1.0 + rsqr(x) + (D2mf)/(Dpk*Dpp)
          );
    
    t6 = (1.0/3.0)*(rpow(OmM(eta),3.0)*rpow(H(eta),4.0)/rpow(invH0,4) )
        *(
          M3_USER(eta) - M2_USER(eta)*(M2_USER(eta) + JFL_USER(eta,-x,k,p)*(3.0+2.0*omegaBD_u))
                /(PiF_USER(eta,kpluspm))
          ) / ( rsqr(PiF_USER(eta,p)) * PiF_USER(eta,k) );

    K3tmp = t1 + t2 + t3 + t4 + t5 + t6;

    return (K3tmp);
}

local real S2a_USER(real eta, real x, real k, real p)
{
    real Dtmp, kplusp;

    kplusp = kpp(x,k,p);
    Dtmp = f1(eta)*mu_USER(eta,kplusp);

    return Dtmp;
}

local real S2b_USER(real eta, real x, real k, real p)
{
    real Dtmp, kplusp;
    
    kplusp = kpp(x,k,p);
    Dtmp = f1(eta)*(mu_USER(eta,k)+mu_USER(eta,p)-mu_USER(eta,kplusp));
    
    return Dtmp;
}

local real S2FL_USER(real eta, real x, real k, real p)
{
    real Dtmp, kplusp;
    
    kplusp = kpp(x,k,p);
    Dtmp = f1(eta)*(
                    M1_USER(eta)/(3.0*PiF_USER(eta,kplusp))
                    *KFL2_USER(eta,x,k,p)
                    );
    return Dtmp;
}

local real S2dI_USER(real eta, real x, real k, real p)
{
    real Dtmp, kplusp;
    
    kplusp = kpp(x,k,p);
    Dtmp = (1.0/6.0)*
            rsqr(OmM(eta)*H(eta)/(rexp(eta)*invH0))
        * ( (rsqr(kplusp)*M2_USER(eta)) / (PiF_USER(eta,kplusp)*PiF_USER(eta,k)*PiF_USER(eta,p)));

    return Dtmp;
}

local real SD2_USER(real eta, real x, real k, real p)
{
    real Dtmp;

    Dtmp = S2a_USER(eta, x, k, p) -  S2b_USER(eta, x, k, p)*rsqr(x)
        + S2FL_USER(eta, x, k, p) - S2dI_USER(eta, x, k, p);

    return Dtmp;
}

local real S3I_USER(real eta, real x, real k, real p, real Dpk, real Dpp,
                real D2f, real D2mf)
{
    real Stmp, kplusp, kpluspm;

    kplusp = kpp(x,k,p);
    kpluspm = kpp(-x,k,p);
    Stmp = (
            f1(eta)*(mu_USER(eta,p)+mu_USER(eta,kplusp)-mu_USER(eta,k))*D2f*Dpp
                + SD2_USER(eta,x,k,p)*Dpk*Dpp*Dpp
            )*(1.0 - rsqr(x))/(1.0 + rsqr(p/k) + 2.0*(p/k)*x)
        + (
           f1(eta)*(mu_USER(eta,p)+mu_USER(eta,kpluspm)-mu_USER(eta,k))*D2mf*Dpp
            + SD2_USER(eta,-x,k,p)*Dpk*Dpp*Dpp
           )*(1.0 - rsqr(x))/(1.0 + rsqr(p/k) - 2.0*(p/k)*x);

    return (Stmp);
}

local real S3IIplus_USER(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f)
{
    real Stmp, kplusp;
    
    kplusp = kpp(x,k,p);
    
    Stmp =
    -f1(eta)*(mu_USER(eta,p)+mu_USER(eta,kplusp)-2.0*mu_USER(eta,k))
    * Dpp*( D2f + Dpk*Dpp*rsqr(x) )
    
    -f1(eta)*(mu_USER(eta,kplusp)-mu_USER(eta,k))*Dpk*Dpp*Dpp
    
    -(
      (M1_USER(eta)/(3.0*PiF_USER(eta,kplusp))) * f1(eta)*KFL2_USER(eta,x,k,p)
      -rsqr(OmM(eta)*H(eta)/invH0)
      * (M2_USER(eta)*kplusp*kplusp*rexp(-2.0*eta))
      / (6.0*PiF_USER(eta,kplusp)*PiF_USER(eta,k)*PiF_USER(eta,p))
      )*Dpk*Dpp*Dpp;
    
    return (Stmp);
}

local real S3IIminus_USER(real eta, real x, real k, real p, real Dpk, real Dpp, real D2mf)
{
    real Stmp, kpluspm;

    kpluspm = kpp(-x,k,p);

    Stmp =
    -f1(eta)*(mu_USER(eta,p)+mu_USER(eta,kpluspm)-2.0*mu_USER(eta,k))
    * Dpp*( D2mf + Dpk*Dpp*rsqr(x) )
    
    -f1(eta)*(mu_USER(eta,kpluspm)-mu_USER(eta,k))*Dpk*Dpp*Dpp
    
    -(
      (M1_USER(eta)/(3.0*PiF_USER(eta,kpluspm))) * f1(eta)*KFL2_USER(eta,-x,k,p)
      -rsqr(OmM(eta)*H(eta)/invH0)
      * (M2_USER(eta)*kpluspm*kpluspm*rexp(-2.0*eta))
      / (6.0*PiF_USER(eta,kpluspm)*PiF_USER(eta,k)*PiF_USER(eta,p))
      )*Dpk*Dpp*Dpp;
    
    return (Stmp);
}

local real S3II_USER(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f, real D2mf)
{
    real Stmp;
    
    Stmp =  S3IIplus_USER(eta, x, k, p, Dpk, Dpp, D2f)
         + S3IIminus_USER(eta, x, k, p, Dpk, Dpp, D2mf);

    return (Stmp);
}

local real S3FLplus_USER(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f)
{
    real Stmp, kplusp;
    
    kplusp = kpp(x,k,p);
    
    Stmp = f1(eta)*(M1_USER(eta)/(3.0*PiF_USER(eta,k)))
    *(
        (2.0*rsqr(p+k*x)/rsqr(kplusp) - 1.0 - (k*x)/p )
        *( mu_USER(eta,p)-1.0 )* D2f * Dpp
      
        + ( (rsqr(p) + 3.0*k*p*x + 2.0*k*k * x*x)/rsqr(kplusp) )
        *( mu_USER(eta,kplusp) - 1.0) * D2phiplus_USER(eta,x,k,p,Dpk,Dpp,D2f) * Dpp
      
      + 3.0*rsqr(x)*( mu_USER(eta,k) + mu_USER(eta,p) - 2.0 ) * Dpk * Dpp * Dpp
    );
    
    return (Stmp);
}

local real S3FLminus_USER(real eta, real x, real k, real p, real Dpk, real Dpp, real D2mf)
{
    real Stmp, kpluspm;
    
    kpluspm = kpp(-x,k,p);
    
    Stmp = f1(eta)*(M1_USER(eta)/(3.0*PiF_USER(eta,k)))
    *(
      (2.0*rsqr(p-k*x)/rsqr(kpluspm) - 1.0 + (k*x)/p )
      *( mu_USER(eta,p)-1.0 )* D2mf * Dpp
      
      + ( (rsqr(p) - 3.0*k*p*x + 2.0*k*k * x*x)/rsqr(kpluspm) )
      *( mu_USER(eta,kpluspm) - 1.0) * D2phiminus_USER(eta,x,k,p,Dpk,Dpp,D2mf) * Dpp
      
      + 3.0*rsqr(x)*( mu_USER(eta,k) + mu_USER(eta,p) - 2.0 ) * Dpk * Dpp * Dpp
      );
    
    return (Stmp);
}

local real S3FL_USER(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f, real D2mf)
{
    real Stmp;

    Stmp = S3FLplus_USER(eta, x, k, p, Dpk, Dpp, D2f)
        + S3FLminus_USER(eta, x, k, p, Dpk, Dpp, D2mf);

    return (Stmp);
}

local real S3dI_USER(real eta, real x, real k, real p, real Dpk, real Dpp,
                 real D2f, real D2mf)
{
    real Stmp;
    
    Stmp = -(rsqr(k)/rexp(2.0*eta))
        *(1.0/(6.0*PiF_USER(eta,k)))
        *K3dI_USER(eta,x,k,p,Dpk,Dpp,D2f,D2mf)*Dpk*Dpp*Dpp;

    return (Stmp);
}

//
// END :: THIRD ORDER (Dsymmetric, five second order differential equations)
// ===========================================================================

// ===========================================================================
// End: USER Model
// ===========================================================================



// ===========================================================================
// MODELS SECTION 5
// ==========================================

//=============================================================
// Begin: Modified gravity model reading and writing parameters

global void ReadMGModelParameterFile(char *fname)
{
#define DOUBLE 1
#define STRING 2
#define INT 3
#define BOOLEAN 4
#define MAXTAGS 300
    
    FILE *fd,*fdout;
    
    char buf[200],buf1[200],buf2[200],buf3[200];
    int  i,j,nt;
    int  id[MAXTAGS];
    void *addr[MAXTAGS];
    char tag[MAXTAGS][50];
    int  errorFlag=0;
    
    nt=0;
    
    // Modified gravity model parameters:
    IPName(nHS_u,"nHS");
    RPName(fR0_u,"fR0");
    RPName(screening_u,"screening");
    //
    if((fd=fopen(fname,"r"))) {
        while(!feof(fd)) {
            fgets(buf,200,fd);
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<1)
                continue;
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<2)
                *buf2='\0';
            if(buf1[0]=='%')
                continue;
            for(i=0,j=-1;i<nt;i++)
                if(strcmp(buf1,tag[i])==0) {
                    j=i;
                    tag[i][0]=0;
                    break;
                }
            if(j>=0) {
                switch(id[j]) {
                    case DOUBLE:
                        *((double*)addr[j])=atof(buf2);
                        break;
                    case STRING:
                        strcpy(addr[j],buf2);
                        break;
                    case INT:
                        *((int*)addr[j])=atoi(buf2);
                        break;
                    case BOOLEAN:
                        if (strchr("tTyY1", *buf2) != NULL) {
                            *((bool*)addr[j])=TRUE;
                        } else
                            if (strchr("fFnN0", *buf2) != NULL)  {
                                *((bool*)addr[j])=FALSE;
                            } else {
                                error("getbparam: %s=%s not bool\n",buf1,buf2);
                            }
                        break;
                }
            } else {
                fprintf(stdout, "Error in file %s: Tag '%s' %s.\n",
                        fname, buf1, "not allowed or multiple defined");
                errorFlag=1;
            }
        }
        fclose(fd);
    } else {
        fprintf(stdout,"Parameter file %s not found.\n", fname);
        errorFlag=1;
        exit(1);
    }
    
    for(i=0;i<nt;i++) {
        if(*tag[i]) {
            fprintf(stdout,
                    "Error. I miss a value for tag '%s' in parameter file '%s'.\n",
                    tag[i],fname);
            exit(0);
        }
    }
#undef DOUBLE
#undef STRING
#undef INT
#undef BOOLEAN
#undef MAXTAGS
}

#define FMTT    "%-35s%s\n"
#define FMTI    "%-35s%d\n"
#define FMTR    "%-35s%g\n"

global void PrintMGModelParameterFile(char *fname)
{
    FILE *fdout;
    char buf[200];
    
    sprintf(buf,"%s/%s%s%s",gd.tmpDir,fname,cmd.suffixModel,"-usedvalues");
    if(!(fdout=fopen(buf,"w"))) {
        fprintf(stdout,"error opening file '%s' \n",buf);
        exit(0);
    } else {
        fprintf(fdout,"%s\n",
                "%-------------------------------------------------------------------");
        fprintf(fdout,"%s %s\n","% Modified gravity parameter model input file for:",gd.headline0);
        fprintf(fdout,"%s\n","%");
        fprintf(fdout,"%s %s: %s\n%s\t    %s\n","%",gd.headline1,gd.headline2,"%",
                gd.headline3);
        fprintf(fdout,"%s\n%s\n",
                "%-------------------------------------------------------------------",
                "%");
        //
        // Modified gravity model parameters:
        fprintf(fdout,FMTI,"nHS",nHS_u);
        fprintf(fdout,FMTR,"fR0",fR0_u);
        fprintf(fdout,FMTR,"screening",screening_u);
        //
        fprintf(fdout,"\n\n");
    }
    fclose(fdout);
}

#undef FMTT
#undef FMTI
#undef FMTR

// End: Modified gravity model reading and writing parameters
//=============================================================
