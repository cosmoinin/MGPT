/*==============================================================================
 MODULE: mglpt_fns.c				[mglpt]
 ==============================================================================*/

#include "globaldefs.h"
#include "protodefs.h"

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


global  real f1(real eta)
{
    real f1tmp;
    
    f1tmp=3.0/(2.0*(1.0
                          + ((1.0-cmd.om)/cmd.om)*rexp(3.0*eta)
                    )
               );
    return (f1tmp);
}

global  real f2(real eta)
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

global  real psLCDMf(real k)
{
    pointPSTableptr p, pf, pi;
    int jl, ju, jm;
    real dk, psftmp;
    bool ascnd;
    
    pi = PSLCDMtab;
    pf = PSLCDMtab+nPSTable-1;

    if ( k < kPos(pi) || k > kPos(pf) || nPSTable < 2 )
        error("\n\npsLCDMf: k is out of range or nPSTable is wrong... %g %g %g\n",
              k,kPos(pi),kPos(pf));

    ascnd = (kPos(pf) >= kPos(pi));
    
    jl=0;
    ju=nPSTable-1;
    while (ju-jl > 1) {
        jm = (ju+jl) >> 1;
        if (k >= kPos(pi+jm) == ascnd)
            jl=jm;
        else
            ju=jm;
    }
    
    p = PSLCDMtab + jl;
    dk = kPos(p+1)-kPos(p);
    psftmp = PS(p)+(PS(p+1)-PS(p))*(k-kPos(p))/dk;
    return (psftmp);
}

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

