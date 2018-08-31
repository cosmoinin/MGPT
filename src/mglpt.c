/*==============================================================================
 MODULE: mglpt.c				[mgpt]
==============================================================================*/


#include "globaldefs.h"
#include "protodefs.h"
#include "models.h"


local void computingQRs(void);

void MainLoop(void)
{
    computingQRs();
    biasterms_processing();

}

#define FMTQRDATHD    "%1s%4s%8s%8s%8s%8s%8s%8s%7s%8s%7s%7s%9s%7s%9s%7s%9s"

#define FMTQRDAT	\
"%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e \
%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n"

#define fpfnameQsRs   "CLPT/kfunctions.dat"
#define FMTQR       "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n"


local void computingQRs(void)
{
    char namebuf[256];
    stream outstrQsRs;
    global_QRs qrs;
    real Q1, Q2, Q3, Q8, Q9, Q13, QI;
    real Q5, Q7, Q11, Q12;
    real RI, R1p2;
    real R1, R2;
    real ki, kval, dk;
    int i;
    real aTime, bTime;
    real Dpk, PSLMG;

    bTime = cputime();

    outstrQsRs = stropen(fpfnameQsRs,"w!");

    fprintf(outstrQsRs,"%1s%5s%12s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%12s%12s%12s",
            "#","k","<Q1>","<Q2>","<Q3>",
            "<Q5>","<Q7>","<Q8>","<Q9>",
            "<Q11>","<Q12>","<Q13>","<QI>",
            "<R1>","<R2>","<R1p2>","<RI>","<Dpk>","<PSLMG>\n");

    fprintf(outstrQsRs,
            "%1s%6s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%12s%12s%12s",
            "#","<1>","<2>","<3>","<4>","<5>","<6>","<7>","<8>","<9>",
            "<10>","<11>","<12>","<13>","<14>","<15>","<16>","<17>","<18>\n");

    fprintf(stdout,"\nTesting Nk values from kmin to kmax of the power spectrum: %d %g %g\n\n",
            cmd.Nk, cmd.kmin, cmd.kmax);
    if (cmd.Nk==1) {
        dk = 0.;
    } else
        dk = (rlog10(cmd.kmax) - rlog10(cmd.kmin))/((real)(cmd.Nk - 1));
 
    for (i=1; i<=cmd.Nk; i++) {
        aTime = cputime();
        kval = rlog10(cmd.kmin) + dk*((real)(i - 1));
        ki = rpow(10.0,kval);
        fprintf(stdout,"i: %d :: ki: %g :: ",i,ki);
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

        fprintf(stdout,FMTQR,
                Q1, Q2, Q3,
                Q5, Q7, Q8, Q9,
                Q11, Q12, Q13, QI,
                R1, R2, R1p2, RI,
                cputime()-aTime);

        fprintf(outstrQsRs,FMTQRDAT,
                ki, Q1, Q2, Q3,
                Q5, Q7, Q8, Q9,
                Q11, Q12, Q13, QI,
                R1, R2, R1p2, RI,
                Dpk, PSLMG
                );

        fflush(outstrQsRs);
        fflush(stdout);
    }

    fclose(outstrQsRs);

    fprintf(stdout,"\nTotal time to compute all k's: %g",cputime()-bTime);
}

#undef FMTQRDATHD
#undef FMTQRDAT
#undef fpfnameQsRs

#undef FMTQR


