/*==============================================================================
 MODULE: mglpt.c				[mgpt]
==============================================================================*/

#include "globaldefs.h"
#include "protodefs.h"
#include "models.h"

local void computingQRs(void);
local void loopQsRs(stream outstr, int imin, int imax, real dk);

void MainLoop(void)
{
    computingQRs();
    biasterms_processing();
    qfunctions_processing();
    CLPT_correlation_processing();
}


#define FMTQRDAT	\
"%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e \
%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n"

#define FMTQR       "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n"

local void computingQRs(void)
{
    char namebuf[256];
    stream outstrQsRs;
    real ki, kval, dk;
    real bTime;
    real kk, rr, xv, k2;
    global_D2v2_ptr ptmp;

    bTime = second();

    if (model_int_flag==LCDM) {
        gd.p = cmd.kmin;
        kk = gd.p;
        ki = cmd.kmin;
        rr = kk/ki;
        xv = xGL(pGL)[1];
        k2 = ki * rsqrt(1.0 + rsqr(rr) - 2.0*rr*xv);
        ptmp = DsSecondOrder_func(ki, ki*rr, k2);
        KA_LCDM = DA2D2(ptmp)/( (3.0/7.0)*Dpk1D2v2(ptmp)*Dpk2D2v2(ptmp) );
        KB_LCDM = KA_LCDM;
        fprintf(stdout,"\n\nKA_LCDM, KB_LCDM: %g %g\n",KA_LCDM, KB_LCDM);
    }

    outstrQsRs = stropen(gd.fpfnamekfun,"w!");

    fprintf(outstrQsRs,"%1s%5s%12s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%13s",
            "#","k","<Q1>","<Q2>","<Q3>",
            "<Q5>","<Q7>","<Q8>","<Q9>",
            "<Q11>","<Q12>","<Q13>","<QI>",
            "<R1>","<R2>","<R1p2>","<RI>","<Dpk>","<PSLMG>\n");
    
    fprintf(outstrQsRs,
            "%1s%6s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%12s",
            "#","<1>","<2>","<3>","<4>","<5>","<6>","<7>","<8>","<9>",
            "<10>","<11>","<12>","<13>","<14>","<15>","<16>","<17>","<18>\n");

    fprintf(stdout,"\nTesting Nk values from kmin to kmax of the power spectrum: %d %g %g\n\n",
            cmd.Nk, cmd.kmin, cmd.kmax);
    if (cmd.Nk==1) {
        dk = 0.;
    } else
        dk = (rlog10(cmd.kmax) - rlog10(cmd.kmin))/((real)(cmd.Nk - 1));
//
//    fprintf(stdout,"\nUsing %g tolerance quad...\n",cmd.epsquad);
    loopQsRs(outstrQsRs, 1, cmd.Nk, dk);

    fclose(outstrQsRs);

    fprintf(stdout,"\nTotal time to compute all k functions: %g sec.",second()-bTime);
}

local void loopQsRs(stream outstr, int imin, int imax, real dk)
{
    global_QRs qrs;
    real Q1, Q2, Q3, Q8, Q9, Q13, QI;
    real Q5, Q7, Q11, Q12;
    real RI, R1p2;
    real R1, R2;
    real aTime;
    real kval;
    real ki;
    real Dpk, PSLMG;
    int i;

    if (model_int_flag==LCDM) {
    for (i=imin; i<=imax; i++) {
        aTime = second();
        kval = rlog10(cmd.kmin) + dk*((real)(i - 1));
        ki = rpow(10.0,kval);
//        fprintf(stdout,"i: %d :: ki: %g :: ",i,ki);
        fprintf(stdout,"i: %d :: ki: %e :: ",i,ki);
        fflush(stdout);
        qrs = QsRs_functions_driver_LCDM(gd.xstop, ki);

        Q1 = qrs.Q1;
        Q2 = qrs.Q2;
        Q3 = qrs.Q3;
        Q8 = qrs.Q8;
        Q9 = qrs.Q9;
        Q13 = qrs.Q13;
        QI = qrs.QI;
        Q5 = qrs.Q5;
        Q7 = qrs.Q7;
        Q11 = qrs.Q11;
        Q12 = qrs.Q12;
        RI = qrs.RI;
        R1p2 = qrs.R1p2;
        R1 = qrs.RI;
        R2 = qrs.R2;
        
        Dpk = DpFunction(ki);
        PSLMG = psInterpolation_nr(ki, kPS, pPS, nPSLT);
/*
        fprintf(stdout,FMTQR,
                Q1, Q2, Q3,
                Q5, Q7, Q8, Q9,
                Q11, Q12, Q13, QI,
                R1, R2, R1p2, RI,
                second()-aTime);
*/
        fprintf(stdout,"time = %f\n",
                (second()-aTime));
        
        fprintf(outstr,FMTQRDAT,
                ki, Q1, Q2, Q3,
                Q5, Q7, Q8, Q9,
                Q11, Q12, Q13, QI,
                R1, R2, R1p2, RI,
                Dpk, PSLMG
                );

        fflush(outstr);
        fflush(stdout);
    }
    } else {
        for (i=imin; i<=imax; i++) {
            aTime = second();
            kval = rlog10(cmd.kmin) + dk*((real)(i - 1));
            ki = rpow(10.0,kval);
//            fprintf(stdout,"i: %d :: ki: %g :: ",i,ki);
            fprintf(stdout,"i: %d :: ki: %e :: ",i,ki);
            fflush(stdout);
            qrs = QsRs_functions_driver(gd.xstop, ki);
            
            Q1 = qrs.Q1;
            Q2 = qrs.Q2;
            Q3 = qrs.Q3;
            Q8 = qrs.Q8;
            Q9 = qrs.Q9;
            Q13 = qrs.Q13;
            QI = qrs.QI;
            Q5 = qrs.Q5;
            Q7 = qrs.Q7;
            Q11 = qrs.Q11;
            Q12 = qrs.Q12;
            RI = qrs.RI;
            R1p2 = qrs.R1p2;
            R1 = qrs.R1;
            R2 = qrs.R2;
            
            Dpk = DpFunction(ki);
            PSLMG = psInterpolation_nr(ki, kPS, pPS, nPSLT);
/*
            fprintf(stdout,FMTQR,
                    Q1, Q2, Q3,
                    Q5, Q7, Q8, Q9,
                    Q11, Q12, Q13, QI,
                    R1, R2, R1p2, RI,
                    second()-aTime);
*/
            fprintf(stdout,"time = %f\n",
                    (second()-aTime));
            
            fprintf(outstr,FMTQRDAT,
                    ki, Q1, Q2, Q3,
                    Q5, Q7, Q8, Q9,
                    Q11, Q12, Q13, QI,
                    R1, R2, R1p2, RI,
                    Dpk, PSLMG
                    );
            
            fflush(outstr);
            fflush(stdout);
        }
    }
}

#undef FMTQRDAT
#undef FMTQR

