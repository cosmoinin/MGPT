/*==============================================================================
 MODULE: mglpt_fns.c				[mglpt]
 ==============================================================================*/

#include "globaldefs.h"
#include "protodefs.h"

/*
global real OmM(real eta)
{
    real OmMtmp;
    
    OmMtmp=1.0/(1.0 + ( (1.0-cmd.om)/cmd.om)*rexp(3.0*eta) );
    return (OmMtmp);
}

global  real H(real eta)
{
    real Htmp;
    
    Htmp=rsqrt(cmd.om*rexp(-3.0*eta)+(1.0-cmd.om));
    return (Htmp);
}


global  real f1(real eta)  //This is 3/2 OmegaM
{
    real f1tmp;
    
    f1tmp=3.0/(2.0*(1.0
                          + ((1.0-cmd.om)/cmd.om)*rexp(3.0*eta)
                    )
               );
    return (f1tmp);
}

global  real f2(real eta) //This is 2 + H'(t)/H^2
{
    real f2tmp;
    
    f2tmp=2.0 - 3.0/(2.0*(1.0
                            + ((1.0-cmd.om)/cmd.om)*rexp(3.0*eta)
                          )
                     );
    return (f2tmp);
}

global real A0(real eta)
{
    real A0tmp;
    
    A0tmp = 1.5*OmM(eta)*rsqr(H(eta))/rsqr(H02);
    
    return (A0tmp);
}
*/

/*
global  real psInterpolation_nr(real k, double kPS[], double pPS[], int nPS)
{
    pointPSTableptr pf, pi;
    real psftmp;
    real dps;
    
    pi = PSLCDMtab;
    pf = PSLCDMtab+nPSTable-1;

    splint(kPS,pPS,pPS2,nPS,k,&psftmp);

    return (psftmp);
}
*/
