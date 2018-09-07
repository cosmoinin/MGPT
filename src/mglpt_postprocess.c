/*==============================================================================
 MODULE: mglpt_postprocess.c			[mgpt]
 ==============================================================================*/

#include "globaldefs.h"
#include "protodefs.h"

#define EPSQ 1.0e-6

local  real Interpolation_nr(real k, double kPS[], double pPS[], int nPS, double pPS2[]);

local real sigma2L_function_int(real y);
local real sigma2L_function(void);

local real PSLF(real k);
local real Q1F(real k);
local real Q2F(real k);
local real Q3F(real k);
local real Q5F(real k);
local real QIF(real k);
local real R1F(real k);
local real R2F(real k);
local real R1plus2F(real k);
local real RIF(real k);

local real tildeV(real k);
local real tildeT(real k);
local real xL(real k, real q);
local real xloop(real k, real q);
local real yL(real k, real q);
local real yloop(real k, real q);
local real x10(real k, real q);
local real y10(real k, real q);

local void InputQsRsTable(void);

global_qfunctions qfunctions(real qi);
global_corrfunctions correlation_functions(real qi);

local void postprocess_string_to_int(string, int *);

// Qs and Rs table structure
typedef struct _pointQsRsTable {
    real k;             //1
    real Q1;            //2
    real Q2;            //3
    real Q3;            //4
    real Q5;            //5
    real Q7;            //6
    real Q8;            //7
    real Q9;            //8
    real Q11;           //9
    real Q12;           //10
    real Q13;           //11
    real QI;            //12
    real R1;            //13
    real R2;            //14
    real R1plus2;       //15
    real RI;            //16
    real Dpk;           //17
    real PSMGL;         //18
} pointQsRsTable, *pointQsRsTableptr;

local int nQsRsTable;
local pointQsRsTableptr QsRstab;

#define kQsRs(x)    (((pointQsRsTableptr) (x))->k)
#define Q1QsRs(x)    (((pointQsRsTableptr) (x))->Q1)
#define Q2QsRs(x)    (((pointQsRsTableptr) (x))->Q2)
#define Q3QsRs(x)    (((pointQsRsTableptr) (x))->Q3)
#define Q5QsRs(x)    (((pointQsRsTableptr) (x))->Q5)
#define Q7QsRs(x)    (((pointQsRsTableptr) (x))->Q7)
#define Q8QsRs(x)    (((pointQsRsTableptr) (x))->Q8)
#define Q9QsRs(x)    (((pointQsRsTableptr) (x))->Q9)
#define Q11QsRs(x)    (((pointQsRsTableptr) (x))->Q11)
#define Q12QsRs(x)    (((pointQsRsTableptr) (x))->Q12)
#define Q13QsRs(x)    (((pointQsRsTableptr) (x))->Q13)
#define QIQsRs(x)    (((pointQsRsTableptr) (x))->QI)
#define R1QsRs(x)    (((pointQsRsTableptr) (x))->R1)
#define R2QsRs(x)    (((pointQsRsTableptr) (x))->R2)
#define R1plus2QsRs(x)    (((pointQsRsTableptr) (x))->R1plus2)
#define RIQsRs(x)    (((pointQsRsTableptr) (x))->RI)
#define DpkQsRs(x)    (((pointQsRsTableptr) (x))->Dpk)
#define PSMGLQsRs(x)    (((pointQsRsTableptr) (x))->PSMGL)

local pointQsRsTableptr PQsRstab;


#define BIASTERMS                   3

global void PostProcessing(void)
{
    int process_int;

    fprintf(stdout,"\n\nStarting postprocessing...\n");
    
    postprocess_string_to_int(cmd.options, &process_int);
    switch (process_int){
        case BIASTERMS: biasterms_processing(); break;
        default: error("\nUnknown postprocessing type %s\n\n",cmd.options);
    }
}

local void postprocess_string_to_int(string process_str,int *process_int)
{
    *process_int = -1;
    if (strcmp(process_str,"biasterms") == 0)               *process_int=BIASTERMS;
}

#undef BIASTERMS

local real *kTab;
local real *Q1T;
local real *Q1T2;
local real *Q2T;
local real *Q2T2;
local real *Q3T;
local real *Q3T2;
local real *Q5T;
local real *Q5T2;
local real *Q8T;
local real *Q8T2;
local real *QIT;
local real *QIT2;
local real *R1T;
local real *R1T2;
local real *R2T;
local real *R2T2;
local real *R1plus2T;
local real *R1plus2T2;
local real *RIT;
local real *RIT2;
local real *PSLMGT;
local real *PSLMGT2;


global void qfunctions_processing(void)
{
    stream outstr;
    pointQsRsTableptr p;
    real aTime;
    real qval;
    real qi;
    real dq;
    real qmin, qmax;
    int i, Nq;
    
    global_qfunctions qfun;
    global_corrfunctions corrfun;
    
    fprintf(stdout,"\n\nqfunctions processing...\n");
    InputQsRsTable();
    
    kTab = dvector(1,nQsRsTable);
    Q1T = dvector(1,nQsRsTable);
    Q1T2 = dvector(1,nQsRsTable);
    Q2T = dvector(1,nQsRsTable);
    Q2T2 = dvector(1,nQsRsTable);
    Q5T = dvector(1,nQsRsTable);
    Q5T2 = dvector(1,nQsRsTable);
    Q8T = dvector(1,nQsRsTable);
    Q8T2 = dvector(1,nQsRsTable);
    QIT = dvector(1,nQsRsTable);
    QIT2 = dvector(1,nQsRsTable);
    R1T = dvector(1,nQsRsTable);
    R1T2 = dvector(1,nQsRsTable);
    R2T = dvector(1,nQsRsTable);
    R2T2 = dvector(1,nQsRsTable);
    R1plus2T = dvector(1,nQsRsTable);
    R1plus2T2 = dvector(1,nQsRsTable);
    RIT = dvector(1,nQsRsTable);
    RIT2 = dvector(1,nQsRsTable);
    PSLMGT = dvector(1,nQsRsTable);
    PSLMGT2 = dvector(1,nQsRsTable);
    
    i=1;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        kTab[i] = kQsRs(p);
        Q1T[i] = Q1QsRs(p);
        Q2T[i] = Q2QsRs(p);
        Q5T[i] = Q5QsRs(p);
        Q8T[i] = Q8QsRs(p);
        QIT[i] = QIQsRs(p);
        R1T[i] = R1QsRs(p);
        R2T[i] = R2QsRs(p);
        R1plus2T[i] = R1plus2QsRs(p);
        RIT[i] = RIQsRs(p);
        PSLMGT[i] = PSMGLQsRs(p);
        i++;
    }
    
    spline(kTab,Q1T,nQsRsTable,1.0e30,1.0e30,Q1T2);
    spline(kTab,Q2T,nQsRsTable,1.0e30,1.0e30,Q2T2);
    spline(kTab,Q5T,nQsRsTable,1.0e30,1.0e30,Q5T2);
    spline(kTab,Q8T,nQsRsTable,1.0e30,1.0e30,Q8T2);
    spline(kTab,QIT,nQsRsTable,1.0e30,1.0e30,QIT2);
    spline(kTab,R1T,nQsRsTable,1.0e30,1.0e30,R1T2);
    spline(kTab,R2T,nQsRsTable,1.0e30,1.0e30,R2T2);
    spline(kTab,R1plus2T,nQsRsTable,1.0e30,1.0e30,R1plus2T2);
    spline(kTab,RIT,nQsRsTable,1.0e30,1.0e30,RIT2);
    spline(kTab,PSLMGT,nQsRsTable,1.0e30,1.0e30,PSLMGT2);
    //
    qmin = 0.0001;
    qmax = 300.0;
    Nq = 1200;
    fprintf(stdout,"\nTesting Nq values from qmin to qmax of the power spectrum: %d %g %g\n\n",
            Nq, qmin, qmax);
    if (Nq==1) {
        dq = 0.;
    } else
        dq = (rlog10(qmax) - rlog10(qmin))/((real)(Nq - 1));
    //
    fprintf(stdout,"\n\nWriting q functions to file %s...",gd.fpfnameqfunctions);
    outstr = stropen(gd.fpfnameqfunctions,"w!");
    for (i=1; i<=Nq; i++) {
        aTime = cputime();
        qval = rlog10(qmin) + dq*((real)(i - 1));
        qi = rpow(10.0,qval);
        fflush(stdout);
        qfun = qfunctions(qi);
        corrfun = correlation_functions(qi);
        fprintf(outstr,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
                qfun.q,
                qfun.XL,
                qfun.YL,
                qfun.Xloop,
                qfun.Yloop,
                qfun.VT,
                qfun.TT,
                qfun.X10,
                qfun.Y10,
                qfun.U10L,
                qfun.U10loop,
                qfun.U11,
                qfun.U20,
                corrfun.xi,
                corrfun.Lapxi,
                corrfun.nabla4xi
                );
    }
    fclose(outstr);
    fprintf(stdout," finished writing.\n");
    
    free_dvector(PSLMGT2,1,nQsRsTable);
    free_dvector(PSLMGT,1,nQsRsTable);
    free_dvector(RIT2,1,nQsRsTable);
    free_dvector(RIT,1,nQsRsTable);
    free_dvector(R1plus2T2,1,nQsRsTable);
    free_dvector(R1plus2T,1,nQsRsTable);
    free_dvector(R2T2,1,nQsRsTable);
    free_dvector(R2T,1,nQsRsTable);
    free_dvector(R1T2,1,nQsRsTable);
    free_dvector(R1T,1,nQsRsTable);
    free_dvector(QIT2,1,nQsRsTable);
    free_dvector(QIT,1,nQsRsTable);
    free_dvector(Q8T2,1,nQsRsTable);
    free_dvector(Q8T,1,nQsRsTable);
    free_dvector(Q5T2,1,nQsRsTable);
    free_dvector(Q5T,1,nQsRsTable);
    free_dvector(Q2T2,1,nQsRsTable);
    free_dvector(Q2T,1,nQsRsTable);
    free_dvector(Q1T2,1,nQsRsTable);
    free_dvector(Q1T,1,nQsRsTable);
    free_dvector(kTab,1,nQsRsTable);
    free(PQsRstab);
}


#define FMTBIASTERMDAT    \
"%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n"

global void biasterms_processing(void)
{
    stream outstr;
    pointQsRsTableptr p;
    real sigma2L;
    real P22, P13;
    real a10, a01, a20, a11, a02;
    real k, PSLT;
    real Ploop;
    real a02Off;
    int i;

    fprintf(gd.outlog,"\n\nBias terms processing...\n");
    InputQsRsTable();
    
    kTab = dvector(1,nQsRsTable);
    PSLMGT = dvector(1,nQsRsTable);
    PSLMGT2 = dvector(1,nQsRsTable);
    i=1;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        kTab[i] = kQsRs(p);
        PSLMGT[i] = PSMGLQsRs(p);
        i++;
    }
    spline(kTab,PSLMGT,nQsRsTable,1.0e30,1.0e30,PSLMGT2);

    sigma2L = sigma2L_function();
    fprintf(stdout,"\n\nsigma2L = %g\n",sigma2L);
    
// a02 offset
    a02Off = (1./2.)*Q13QsRs(PQsRstab + 2);

    fprintf(stdout,"\n\nWriting bias terms and the power spectrum to file %s...\n",gd.fpfnameSPTPowerSpectrum);
    outstr = stropen(gd.fpfnameSPTPowerSpectrum,"w!");

    fprintf(outstr,"%1s%5s%12s%11s%11s%11s%11s%11s%11s%11s",
            "#","k","<PSLT>","<P22>","<P13>",
            "<a10>","<a01>","<a20>","<a11>",
            "<a02>\n");
    
    fprintf(outstr,
            "%1s%6s%11s%11s%11s%11s%11s%11s%11s%11s",
            "#","<1>","<2>","<3>","<4>","<5>","<6>","<7>","<8>","<9>\n");


    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {

        k = kQsRs(p);
        PSLT = PSMGLQsRs(p);
        P22 = (9./98.)*Q1QsRs(p) + (3./7.)*Q2QsRs(p) + (1./2.)*Q3QsRs(p);
        
//        P13 = (10./21.)*R1QsRs(p) + (6./7.)*Q2QsRs(p) - sigma2L*k*k*PSLT;
        P13 = (10./21.)*R1QsRs(p) + (6./7.)*R2QsRs(p) - sigma2L*k*k*PSLT;
        
        a10 = (10./21.)*R1QsRs(p) + (6./7.)*R1plus2QsRs(p) + (6./7.)*R2QsRs(p)
                + (6./7.)*Q5QsRs(p) + 2.*Q7QsRs(p) + 2.*(1. - sigma2L*k*k)*PSLT;
        a01 = Q9QsRs(p) + (3./7.)*Q8QsRs(p);
        a20 = Q11QsRs(p) + (1. - sigma2L*k*k)*PSLT
              + Q9QsRs(p) + (6./7.)*R1plus2QsRs(p);
        a11 = 2.*Q12QsRs(p);
        a02 = (1./2.)*Q13QsRs(p) - a02Off;
        Ploop = PSLT + P22 + P13;
        fprintf(outstr,FMTBIASTERMDAT,
                k, PSLT, P22, P13, a10, a01, a20, a11, a02
                );
    }
    fclose(outstr);
//
    free_dvector(PSLMGT2,1,nQsRsTable);
    free_dvector(PSLMGT,1,nQsRsTable);
    free_dvector(kTab,1,nQsRsTable);
    free(PQsRstab);
}
#undef FMTBIASTERMDAT

local void InputQsRsTable(void)
{
    stream outstr;
    pointQsRsTableptr p;
    int i;
//
    fprintf(gd.outlog,"\n\nReading solution Q1 from file %s...\n",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 2, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputQsRsTable: nQsRsTable = %d is absurd\n\n", nQsRsTable);

    PQsRstab = (pointQsRsTableptr) allocate(nQsRsTable * sizeof(pointQsRsTable));
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        kQsRs(p) = inout_xval[i];
        Q1QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution Q2 from file %s...\n",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 3, &nQsRsTable);

    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q2QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution Q3 from file %s...\n",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 4, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q3QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution Q5 from file %s...\n",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 5, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q5QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution Q7 from file %s...\n",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 6, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q7QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution Q8 from file %s...\n",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 7, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q8QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution Q9 from file %s...\n",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 8, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q9QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution Q11 from file %s...\n",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 9, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q11QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution Q12 from file %s...\n",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 10, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q12QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution Q13 from file %s...\n",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 11, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q13QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution QI from file %s...\n",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 12, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        QIQsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution R1 from file %s...\n",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 13, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        R1QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution R2 from file %s...\n",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 14, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        R2QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution R1plus2 from file %s...\n",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 15, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        R1plus2QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution RI from file %s...\n",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 16, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        RIQsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution Dpk from file %s...\n",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 17, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        DpkQsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution PSMGL from file %s...\n",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 18, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        PSMGLQsRs(p) = inout_yval[i];
        ++i;
    }
//
}

//
// AUXILIARY FUNCTIONS FOR q FUNCTIONS COMPUTATION
//

local real tildeV(real k)
{
    // -(3./35.) (QI[k] - 3. Q2[k] + 2. RI[k] - 6. R2[k])
    real func;
    
    func= -(3./35.)*( QIF(k) - 3.*Q2F(k) + 2.*RIF(k) - 6.*R2F(k) );
    
    return (func);
}

local real tildeT(real k)
{
    // -(9./14.)  (QI[k] + 2. Q2[k] + 2. RI[k] + 4. R2[k])
    real func;
    
    func= -(9./14.)*(QIF(k) + 2.*Q2F(k) + 2.*RIF(k) + 4.*R2F(k) );
    
    return (func);
}

local real xL(real k, real q)
{
    // PSL[k] (1./3. - j1[k q]/(k q))
    real func;
    
    func= PSLF(k)*(1./3. - rj1Bessel(k*q)/(k*q) );
    
    return (func);
}

local real xloop(real k, real q)
{
    // (9./98. Q1[k] + 10./21. R1[k]) (1./3. - j1[k q]/(k q))
    real func;
    
    func= ( (9./98.)*Q1F(k) + (10./21.)*R1F(k) )*(1./3. - rj1Bessel(k*q)/(k*q) );
    
    return (func);
}

local real yL(real k, real q)
{
    // PSL[k] j2[k q]
    real func;
    
    func= PSLF(k)*rj2Bessel(k*q);
    
    return (func);
}

local real yloop(real k, real q)
{
    // (9./98. Q1[k] + 10./21. R1[k]) j2[k q]
    real func;
    
    func= ( (9./98.)*Q1F(k) + (10./21.)*R1F(k) )*rj2Bessel(k*q);
    
    return (func);
}

local real x10(real k, real q)
{
    // 1./14. (2. RI[k] - 2. R2[k] + 3. RI[k] j0[k q]
    //    -3 (RI[k] + 2. R2[k] + 2. R1plus2[k] + 2. Q5[k]) j1[k q]/(k q))
    real func;
    
    func= (1./14.)*(
                    2.*RIF(k) - 2.*R2F(k) + 3.*RIF(k)*rj0Bessel(k*q)
                    -3.*( RIF(k) + 2.*R2F(k) + 2.*R1plus2F(k) + 2.*Q5F(k) )
                    *rj1Bessel(k*q)/(k*q)
                    );
    
    return (func);
}

local real y10(real k, real q)
{
    // -3./14. (RI[k] + 2. R2[k] + 2. R1plus2[k] + 2. Q5[k]) (j0[k q]
    //    -3 j1[k q]/(k q))
    real func;
    
    func= (-3./14.)*(
                     RIF(k) + 2.*R2F(k) + 2.*R1plus2F(k) + 2.*Q5F(k)
                     )
    *(rj0Bessel(k*q) - 3.*rj1Bessel(k*q)/(k*q));
    
    return (func);
}

local real PSLF(real k)
{
    real func;
    //    func = psInterpolation_nr(k, kPS, pPS, nPSLT);
    func = Interpolation_nr(k, kTab, PSLMGT, nQsRsTable, PSLMGT2);
    return (func);
}

local real Q1F(real k)
{
    real func;
    func = Interpolation_nr(k, kTab, Q1T, nQsRsTable, Q1T2);
    return (func);
}

local real Q2F(real k)
{
    real func;
    func = Interpolation_nr(k, kTab, Q2T, nQsRsTable, Q2T2);
    return (func);
}

local real Q3F(real k)
{
    real func;
    func = Interpolation_nr(k, kTab, Q3T, nQsRsTable, Q3T2);
    return (func);
}

local real Q5F(real k)
{
    real func;
    func = Interpolation_nr(k, kTab, Q5T, nQsRsTable, Q5T2);
    return (func);
}

local real Q8F(real k)
{
    real func;
    func = Interpolation_nr(k, kTab, Q8T, nQsRsTable, Q8T2);
    return (func);
}

local real QIF(real k)
{
    real func;
    func = Interpolation_nr(k, kTab, QIT, nQsRsTable, QIT2);
    return (func);
}

local real R1F(real k)
{
    real func;
    func = Interpolation_nr(k, kTab, R1T, nQsRsTable, R1T2);
    return (func);
}

local real R2F(real k)
{
    real func;
    func = Interpolation_nr(k, kTab, R2T, nQsRsTable, R2T2);
    return (func);
}

local real R1plus2F(real k)
{
    real func;
    func = Interpolation_nr(k, kTab, R1plus2T, nQsRsTable, R1plus2T2);
    return (func);
}

local real RIF(real k)
{
    real func;
    func = Interpolation_nr(k, kTab, RIT, nQsRsTable, RIT2);
    return (func);
}

local  real Interpolation_nr(real k, double kPS[], double pPS[], int nPS, double pPS2[])
{
    real psftmp;
    
    if ( k < kPS[1] || k > kPS[nPS] )
        fprintf(gd.outlog,"\n\nInterpolation_nr: warning! :: k is out of range... %g\n",k);
    
    splint(kPS,pPS,pPS2,nPS,k,&psftmp);
    
    return (psftmp);
}


// q functions
global_qfunctions qfunctions(real qi)
{
    global_qfunctions_ptr qfunp;
    int i, Nk;
    real kmin, kmax, dk, kvali, kvalim1, ki, kim1;
    real kk;
    real deltak;
//
    real U10Lp, U10LA, U10LB;
    real U10loopp, U10loopA, U10loopB;
    real U11p, U11A, U11B;
    real U20p, U20A, U20B;
//
    real XLp, XLA, XLB;
    real Xloopp, XloopA, XloopB;
    real X10p, X10A, X10B;
    real YLp, YLA, YLB;
    real Yloopp, YloopA, YloopB;
    real Y10p, Y10A, Y10B;
//
    real Vp, preVp, preVA, preVB;
    real Tp, TA, TB;
    
    qfunp = (global_qfunctions_ptr) allocate(1 * sizeof(global_qfunctions));
    
    kmin = 0.0001;
    kmax = 100.0;
    Nk = 1200;
    if (Nk==1)
        dk = 0.;
    else
        dk = (rlog10(kmax) - rlog10(kmin))/((real)(Nk - 1));
    
    U10Lp = 0.;
    U10LA = -kmin*PSLF(kmin)*rj1Bessel(kmin*qi);
    U10loopp = 0.;
    U10loopA = -kmin*( (5./21.)*R1F(kmin) )*rj1Bessel(kmin*qi);
    U11p = 0.;
    U11A = -kmin*( (6./7.)*R1plus2F(kmin) )*rj1Bessel(kmin*qi);
    U20p = 0.;
    U20A = -kmin*( (3./7.)*Q8F(kmin) )*rj1Bessel(kmin*qi);
    //
    XLp = 0.;
    XLA = xL(kmin, qi);
    Xloopp = 0.;
    XloopA = xloop(kmin, qi);
    X10p = 0.;
    X10A = x10(kmin, qi);
    YLp = 0.;
    YLA = yL(kmin, qi);
    Yloopp = 0.;
    YloopA = yloop(kmin, qi);
    Y10p = 0.;
    Y10A = y10(kmin, qi);
//
    preVp = 0.;
    preVA = tildeV(kmin)*rj1Bessel(kmin*qi)/kmin;
    Tp = 0.;
    TA = tildeT(kmin)*rj3Bessel(kmin*qi)/kmin;
//
    for (i=2; i<=Nk; i++) {
        kvali = rlog10(kmin) + dk*((real)(i - 1));
        kvalim1 = rlog10(kmin) + dk*((real)(i - 2));
        ki = rpow(10.0,kvali);
        kim1 = rpow(10.0,kvalim1);
        deltak = (ki - kim1);
        kk = rpow(10.0,kvali);
//
        U10LB = -kk*PSLF(kk)*rj1Bessel(kk*qi);
        U10Lp = U10Lp + (U10LA + U10LB)*deltak/2.0;
        U10LA = U10LB;
//
        U10loopB = -kk*( (5./21.)*R1F(kk) )*rj1Bessel(kk*qi);
        U10loopp = U10loopp + (U10loopA + U10loopB)*deltak/2.0;
        U10loopA = U10loopB;
//
        U11B = -kk*( (6./7.) *R1plus2F(kk) )*rj1Bessel(kk*qi);
        U11p = U11p + (U11A + U11B)*deltak/2.0;
        U11A = U11B;
//
        U20B = -kk*( (3./7.)*Q8F(kk) )*rj1Bessel(kk*qi);
        U20p = U20p + (U20A + U20B)*deltak/2.0;
        U20A = U20B;
//
//
        XLB = xL(kk, qi);
        XLp = XLp + (XLA + XLB)*deltak/2.0;
        XLA = XLB;
//
        XloopB = xloop(kk, qi);
        Xloopp = Xloopp + (XloopA + XloopB)*deltak/2.0;
        XloopA = XloopB;
//
        X10B = x10(kk, qi);
        X10p = X10p + (X10A + X10B)*deltak/2.0;
        X10A = X10B;
//
        YLB = yL(kk, qi);
        YLp = YLp + (YLA + YLB)*deltak/2.0;
        YLA = YLB;
//
        YloopB = yloop(kk, qi);
        Yloopp = Yloopp + (YloopA + YloopB)*deltak/2.0;
        YloopA = YloopB;
//
        Y10B = y10(kk, qi);
        Y10p = Y10p + (Y10A + Y10B)*deltak/2.0;
        Y10A = Y10B;
//
        preVB = tildeV(kk)*rj1Bessel(kk*qi)/kk;
        preVp = preVp + (preVA + preVB)*deltak/2.0;
        preVA = preVB;
//
        TB = tildeT(kk)*rj3Bessel(kk*qi)/kk;
        Tp = Tp + (TA + TB)*deltak/2.0;
        TA = TB;
//
    }
//
    U10Lp /= TWOPI2;
    U10loopp /= TWOPI2;
    U11p /= TWOPI2;
    U20p /= TWOPI2;
    //
    XLp /= PI2;
    Xloopp /= PI2;
    X10p /= PI2;
    YLp /= PI2;
    Yloopp /= PI2;
    Y10p /= PI2;
//
    Vp = preVp/PI2 - (1./5.)*Tp/PI2;
    Tp = Tp/PI2;
    if (qi<0.05) {
        Vp = 0.;
        Tp = 0.;
    }
//
    qqfun(qfunp) = qi;
    U10Lqfun(qfunp) = U10Lp;
    U10loopqfun(qfunp) = U10loopp;
    U11qfun(qfunp) = U11p;
    U20qfun(qfunp) = U20p;
    XLqfun(qfunp) = XLp;
    Xloopqfun(qfunp) = Xloopp;
    X10qfun(qfunp) = X10p;
    YLqfun(qfunp) = YLp;
    Yloopqfun(qfunp) = Yloopp;
    Y10qfun(qfunp) = Y10p;
    VTqfun(qfunp) = Vp;
    TTqfun(qfunp) = Tp;
    
    return *qfunp;
}

// correlation functions
global_corrfunctions correlation_functions(real qi)
{
    global_corrfunctions_ptr corrfunp;
    int i, Nk;
    real kmin, kmax, dk, kvali, kvalim1, ki, kim1;
    real kk;
    real deltak;
    
    real xip, xiA, xiB;
    real Lapxip, LapxiA, LapxiB;
    real nabla4xip, nabla4xiA, nabla4xiB;
//
    corrfunp = (global_corrfunctions_ptr) allocate(1 * sizeof(global_corrfunctions));
    
    kmin = 0.0001;
    kmax = 20.0;
    Nk = 12000;
    if (Nk==1)
        dk = 0.;
    else
        dk = (rlog10(kmax) - rlog10(kmin))/((real)(Nk - 1));
    
    xip = 0.;
    xiA = rsqr(kmin)*PSLF(kmin)*rj0Bessel(kmin*qi);
    Lapxip = 0.;
    LapxiA = rpow(kmin,4.0)*PSLF(kmin)*rj0Bessel(kmin*qi)*rexp(-rsqr(kmin));
    nabla4xip = 0.;
    nabla4xiA = rpow(kmin,6.0)*PSLF(kmin)*rj0Bessel(kmin*qi)*rexp(-rsqr(2.0*kmin));
//
    for (i=2; i<=Nk; i++) {
        kvali = rlog10(kmin) + dk*((real)(i - 1));
        kvalim1 = rlog10(kmin) + dk*((real)(i - 2));
        ki = rpow(10.0,kvali);
        kim1 = rpow(10.0,kvalim1);
        deltak = (ki - kim1);
        kk = rpow(10.0,kvali);
//
        xiB = rsqr(kk)*PSLF(kk)*rj0Bessel(kk*qi);
        xip = xip + (xiA + xiB)*deltak/2.0;
        xiA = xiB;
//
        LapxiB = rsqr(kk)*xiB*rexp(-rsqr(kk));
        Lapxip = Lapxip + (LapxiA + LapxiB)*deltak/2.0;
        LapxiA = LapxiB;
//
        nabla4xiB = rpow(kk,4.0)*xiB*rexp(-rsqr(2*kk));
        nabla4xip = nabla4xip + (nabla4xiA + nabla4xiB)*deltak/2.0;
        nabla4xiA = nabla4xiB;
    }
//
    xip /= TWOPI2;
    Lapxip /= -TWOPI2;
    nabla4xip /= TWOPI2;
//
    qcorrfun(corrfunp) = qi;
    xicorrfun(corrfunp) = xip;
    Lapxicorrfun(corrfunp) = Lapxip;
    nabla4xicorrfun(corrfunp) = nabla4xip;
    
    return *corrfunp;
}

// sigma2L :: integration over the extended MG power spectrum (not using PSMG in QsRs table)
local real sigma2L_function_int(real y)
{
    real p;
    real PSL;
    
    p = rpow(10.0,y);

    PSL = psInterpolation_nr(p, kPS, pPS, nPSLT);

    return p*PSL;
}

local real sigma2L_function(void)
{
    real result;
    real kmin, kmax;
    real ymin, ymax;

    kmin = kPos(PSLT+1);
    kmax = kPos(PSLT+nPSLT-1);
    ymin = rlog10(kmin);
    ymax = rlog10(kmax);

    result= (1.0/SIXPI2)*rlog(10.0)
    *qromo(sigma2L_function_int,ymin,ymax,midpnt,EPSQ);

    return result;

}

#undef EPSQ
