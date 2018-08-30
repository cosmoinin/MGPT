/*==============================================================================
 MODULE: mglpt_postprocess.c			[mgpt]
 Written by: Mario A. Rodriguez-Meza
 Starting date:	January 2018
 Purpose:
 Language: C
 Use:
 Routines and functions:
 External modules, routines and headers:
 Comments and notes:
 Info: Mario A. Rodriguez-Meza
 Depto. de Fisica, ININ
 Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
 e-mail: marioalberto.rodriguez@inin.gob.mx
 http://www.inin.gob.mx/
 
 Major revisions:
 Copyright: (c) 2005-2018 Mario A. Rodriguez-Meza.  All Rights Reserved
 ================================================================================
 Legal matters:
 The author does not warrant that the program and routines it contains
 listed below are free from error or suitable for particular applications,
 and he disclaims all liability from any consequences arising from their	use.
 ==============================================================================*/

#include "globaldefs.h"
#include "protodefs.h"

#define SIXPI2  59.2176

local real Q12_function(real eta, real ki);

local  real Interpolation_nr(real k, double kPS[], double pPS[], int nPS, double pPS2[]);

local real sigma2L_function_int(real y);
local real sigma2L_function_ver2(void);
local real sigma2L_function(real ki);

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
local void InputQsRsTable_old(void);

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

#define fpfnameQsRs "CLPT/kfunctions.dat"
#define fpfnamebiasterms "SPTPowerSpectrum.dat"

local real *kTab;
local real *Q1T;
local real *Q1T2;
local real *Q2T;
local real *Q2T2;
local real *Q3T;
local real *Q3T2;
local real *Q5T;
local real *Q5T2;
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

    fprintf(stdout,"\n\nBias terms processing...\n");
    InputQsRsTable();
    
    sigma2L = sigma2L_function_ver2();
    fprintf(stdout,"\n\nsigma2L = %g\n",sigma2L);
    
// a02 offset
    a02Off = (1./2.)*Q13QsRs(PQsRstab + 2);

    fprintf(stdout,"\n\nWriting bias terms and the power spectrum to file %s...\n",fpfnamebiasterms);
    outstr = stropen(fpfnamebiasterms,"w!");

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
        P13 = (10./21.)*R1QsRs(p) + (6./7.)*Q2QsRs(p) - sigma2L*k*k*PSLT;
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
    free(PQsRstab);
}
#undef FMTBIASTERMDAT

local void InputQsRsTable(void)
{
    stream outstr;
    pointQsRsTableptr p;
    int i;
//
    fprintf(gd.outlog,"\n\nReading solution Q1 from file %s...\n",fpfnameQsRs);
    inout_InputData(fpfnameQsRs, 1, 2, &nQsRsTable);
    
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
    fprintf(gd.outlog,"\n\nReading solution Q2 from file %s...\n",fpfnameQsRs);
    inout_InputData(fpfnameQsRs, 1, 3, &nQsRsTable);

    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q2QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution Q3 from file %s...\n",fpfnameQsRs);
    inout_InputData(fpfnameQsRs, 1, 4, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q3QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution Q5 from file %s...\n",fpfnameQsRs);
    inout_InputData(fpfnameQsRs, 1, 5, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q5QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution Q7 from file %s...\n",fpfnameQsRs);
    inout_InputData(fpfnameQsRs, 1, 6, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q7QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution Q8 from file %s...\n",fpfnameQsRs);
    inout_InputData(fpfnameQsRs, 1, 7, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q8QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution Q9 from file %s...\n",fpfnameQsRs);
    inout_InputData(fpfnameQsRs, 1, 8, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q9QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution Q11 from file %s...\n",fpfnameQsRs);
    inout_InputData(fpfnameQsRs, 1, 9, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q11QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution Q12 from file %s...\n",fpfnameQsRs);
    inout_InputData(fpfnameQsRs, 1, 10, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q12QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution Q13 from file %s...\n",fpfnameQsRs);
    inout_InputData(fpfnameQsRs, 1, 11, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q13QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution QI from file %s...\n",fpfnameQsRs);
    inout_InputData(fpfnameQsRs, 1, 12, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        QIQsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution R1 from file %s...\n",fpfnameQsRs);
    inout_InputData(fpfnameQsRs, 1, 13, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        R1QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution R2 from file %s...\n",fpfnameQsRs);
    inout_InputData(fpfnameQsRs, 1, 14, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        R2QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution R1plus2 from file %s...\n",fpfnameQsRs);
    inout_InputData(fpfnameQsRs, 1, 15, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        R1plus2QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution RI from file %s...\n",fpfnameQsRs);
    inout_InputData(fpfnameQsRs, 1, 16, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        RIQsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution Dpk from file %s...\n",fpfnameQsRs);
    inout_InputData(fpfnameQsRs, 1, 17, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        DpkQsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\n\nReading solution PSMGL from file %s...\n",fpfnameQsRs);
    inout_InputData(fpfnameQsRs, 1, 18, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        PSMGLQsRs(p) = inout_yval[i];
        ++i;
    }
//
}

#undef fpfnameQsRs
#undef fpfnamebiasterms

// sigma2L
local real sigma2L_function(real ki)
{
    int i;
    
    real kmin, kmax;
    real PSLA, PSLB;
    real deltar, kk;

    real sigma2L;
    
    pointPSTableptr p;
    
    kmin = kPos(PSLT+1);
    kmax = kPos(PSLT+nPSLT-1);
    fprintf(gd.outlog,"nPSLT, kmin and kmax :: %d %g %g\n",nPSLT,kmin,kmax);
    
    sigma2L = 0.0;
    
    PSLA = 0.0;
    p = PSLCDMtab;
    //
    for (i=1; i<nPSTable; i++) {
        kk = kPos(p+i);
        PSLB = psInterpolation_nr(kk, kPS, pPS, nPSLT);
        deltar = (kPos(p+i)-kPos(p+i-1))/kk;
        sigma2L += deltar*(PSLA + PSLB)/2.0;
        PSLA = PSLB;
    }

    sigma2L *= (1.0/SIXPI2);
    
    return sigma2L;
}

local real sigma2L_function_int(real y)
{
    real p;
    real PSL;
    
    p = rpow(10.0,y);

    PSL = psInterpolation_nr(p, kPS, pPS, nPSLT);

    return p*PSL;
}

local real sigma2L_function_ver2(void)
{
    real result;
    real kmin, kmax;
    real ymin, ymax;

    kmin = kPos(PSLT+1);
    kmax = kPos(PSLT+nPSLT-1);
    ymin = rlog10(kmin);
    ymax = rlog10(kmax);

    result= (1.0/SIXPI2)*rlog(10.0)
    *qromo(sigma2L_function_int,ymin,ymax,midpnt);
    
    return result;

}
