/*==============================================================================
 MODULE: mglpt_postprocess.c			[mgpt]
 ==============================================================================*/

#include "globaldefs.h"
#include "protodefs.h"

#define KK  5
#define EPSQ 1.0e-6

local  real Interpolation_nr(real k, double kPS[], double pPS[], int nPS, double pPS2[]);

local real sigma2L_function_int(real y);
local real sigma2L_function(void);

// BEGIN :: CLPT correlation auxiliary functions and structures
local real XLF(real q);
local real YLF(real q);
local real XloopF(real q);
local real YloopF(real q);
local real VF(real q);
local real TF(real q);
local real X10F(real q);
local real Y10F(real q);
local real ULF(real q);
local real UloopF(real q);
local real U11F(real q);
local real U20F(real q);
local real xiLF(real q);
local real LapxiF(real q);
local real nabla4xiF(real q);

local real UF(real q);

local real fXF(real q);
local real hYF(real q);
local real BDetF(real ff, real hh);
local real MZAF(real q, real r, real mu);
local real preMZAF(real q, real r, real mu, real fx, real hy);

local real AijGijF(real q, real r, real mu, real ft, real ht, real Xloopt, real Yloopt);
local real WijkGammaijkF(real q, real r, real mu, real ft, real ht, real Vt, real Tt);
local real qigiF(real q, real r, real mu, real ft, real ht);
local real qiqjGijF(real q, real r, real mu, real ft, real ht);

local real MUF(real q, real r, real mu);
local real MA10F(real q, real r, real mu);
local real MxiRF(real q, real r, real mu);
local real MLapxiRF(real q, real r, real mu);
local real Mnabla4xiRF(real q, real r, real mu);
local real MUUF(real q, real r, real mu);
local real MU11F(real q, real r, real mu);
local real MU20F(real q, real r, real mu);
local real MxiR2F(real q, real r, real mu);
local real MUxiRF(real q, real r, real mu);

local real MAF(real q, real r, real mu);
local real MWF(real q, real r, real mu);

local real M10F(real q, real r, real mu);
local real M10LF(real q, real r, real mu);
local real M10loopF(real q, real r, real mu);
local real M01F(real q, real r, real mu);
local real M20F(real q, real r, real mu);
local real M20LF(real q, real r, real mu);
local real M20loopF(real q, real r, real mu);
local real M02F(real q, real r, real mu);
local real M11F(real q, real r, real mu);


local real *qTab;
local real *XLT;
local real *XLT2;
local real *YLT;
local real *YLT2;
local real *XloopT;
local real *XloopT2;
local real *YloopT;
local real *YloopT2;
local real *VT;
local real *VT2;
local real *TT;
local real *TT2;
local real *X10T;
local real *X10T2;
local real *Y10T;
local real *Y10T2;
local real *ULT;
local real *ULT2;
local real *UloopT;
local real *UloopT2;
local real *U11T;
local real *U11T2;
local real *U20T;
local real *U20T2;
local real *xiLT;
local real *xiLT2;
local real *LapxiT;
local real *LapxiT2;
local real *nabla4xiT;
local real *nabla4xiT2;

// qfunctions table structure
typedef struct _pointqfunctionsTable {
    real q;         //1
    real XL;        //2
    real YL;        //3
    real Xloop;     //4
    real Yloop;     //5
    real V;         //6
    real T;         //7
    real X10;       //8
    real Y10;       //9
    real UL;        //10
    real Uloop;     //11
    real U11;       //12
    real U20;       //13
    real xiL;       //14
    real Lapxi;     //15
    real nabla4xi;  //16
} pointqfunctionsTable, *pointqfunctionsTableptr;

local int nqfunctionsTable;
local pointqfunctionsTableptr qfunctionsstab;

#define qqfuncs(x)      (((pointqfunctionsTableptr) (x))->q)
#define XLqfuncs(x)     (((pointqfunctionsTableptr) (x))->XL)
#define YLqfuncs(x)     (((pointqfunctionsTableptr) (x))->YL)
#define Xloopqfuncs(x)     (((pointqfunctionsTableptr) (x))->Xloop)
#define Yloopqfuncs(x)     (((pointqfunctionsTableptr) (x))->Yloop)
#define Vqfuncs(x)     (((pointqfunctionsTableptr) (x))->V)
#define Tqfuncs(x)     (((pointqfunctionsTableptr) (x))->T)
#define X10qfuncs(x)     (((pointqfunctionsTableptr) (x))->X10)
#define Y10qfuncs(x)    (((pointqfunctionsTableptr) (x))->Y10)
#define ULqfuncs(x)    (((pointqfunctionsTableptr) (x))->UL)
#define Uloopqfuncs(x)    (((pointqfunctionsTableptr) (x))->Uloop)
#define U11qfuncs(x)     (((pointqfunctionsTableptr) (x))->U11)
#define U20qfuncs(x)     (((pointqfunctionsTableptr) (x))->U20)
#define xiLqfuncs(x)     (((pointqfunctionsTableptr) (x))->xiL)
#define Lapxiqfuncs(x)    (((pointqfunctionsTableptr) (x))->Lapxi)
#define nabla4xiqfuncs(x)     (((pointqfunctionsTableptr) (x))->nabla4xi)

local pointqfunctionsTableptr Pqfunctab;

global_zacorrfunctions zacorrelation_functions(real ri);
global_clptcorrfunctions clptcorrelation_functions(real ri);

local void InputqfunctionsTable(void);

// END :: CLPT correlation auxiliary functions and structures




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

#define FMTCORRELATIONFUNCTIONSDAT    \
"%e %e %e %e %e %e %e %e \
%e %e %e %e %e %e \
%e\n"

global void CLPT_correlation_processing(void)
{
    stream outstr;
    pointqfunctionsTableptr p;
    real aTime;
    real dr, ri;
    int i;
    
    global_zacorrfunctions zacorrfun;
    global_clptcorrfunctions clptcorrfun;
    
    fprintf(stdout,"\n\nCLPT correlation processing...\n");
    InputqfunctionsTable();
    
    if (nqfunctionsTable<=1) {
        error("\nNumber of q's in qfunctions table must be greater than 1...\n\n",nqfunctionsTable);
    }
    
    qTab = dvector(1,nqfunctionsTable);
    XLT = dvector(1,nqfunctionsTable);
    XLT2 = dvector(1,nqfunctionsTable);
    YLT = dvector(1,nqfunctionsTable);
    YLT2 = dvector(1,nqfunctionsTable);
    XloopT = dvector(1,nqfunctionsTable);
    XloopT2 = dvector(1,nqfunctionsTable);
    YloopT = dvector(1,nqfunctionsTable);
    YloopT2 = dvector(1,nqfunctionsTable);
    VT = dvector(1,nqfunctionsTable);
    VT2 = dvector(1,nqfunctionsTable);
    TT = dvector(1,nqfunctionsTable);
    TT2 = dvector(1,nqfunctionsTable);
    X10T = dvector(1,nqfunctionsTable);
    X10T2 = dvector(1,nqfunctionsTable);
    Y10T = dvector(1,nqfunctionsTable);
    Y10T2 = dvector(1,nqfunctionsTable);
    ULT = dvector(1,nqfunctionsTable);
    ULT2 = dvector(1,nqfunctionsTable);
    UloopT = dvector(1,nqfunctionsTable);
    UloopT2 = dvector(1,nqfunctionsTable);
    U11T = dvector(1,nqfunctionsTable);
    U11T2 = dvector(1,nqfunctionsTable);
    U20T = dvector(1,nqfunctionsTable);
    U20T2 = dvector(1,nqfunctionsTable);
    xiLT = dvector(1,nqfunctionsTable);
    xiLT2 = dvector(1,nqfunctionsTable);
    LapxiT = dvector(1,nqfunctionsTable);
    LapxiT2 = dvector(1,nqfunctionsTable);
    nabla4xiT = dvector(1,nqfunctionsTable);
    nabla4xiT2 = dvector(1,nqfunctionsTable);
    
    i=1;
    for (p=Pqfunctab; p<Pqfunctab+nqfunctionsTable; p++) {
        qTab[i] = qqfuncs(p);
        XLT[i] = XLqfuncs(p);
        YLT[i] = YLqfuncs(p);
        XloopT[i] = Xloopqfuncs(p);
        YloopT[i] = Yloopqfuncs(p);
        VT[i] = Vqfuncs(p);
        TT[i] = Tqfuncs(p);
        X10T[i] = X10qfuncs(p);
        Y10T[i] = Y10qfuncs(p);
        ULT[i] = ULqfuncs(p);
        UloopT[i] = Uloopqfuncs(p);
        U11T[i] = U11qfuncs(p);
        U20T[i] = U20qfuncs(p);
        xiLT[i] = xiLqfuncs(p);
        LapxiT[i] = Lapxiqfuncs(p);
        nabla4xiT[i] = nabla4xiqfuncs(p);
        i++;
    }
    
    spline(qTab,XLT,nqfunctionsTable,1.0e30,1.0e30,XLT2);
    spline(qTab,YLT,nqfunctionsTable,1.0e30,1.0e30,YLT2);
    spline(qTab,XloopT,nqfunctionsTable,1.0e30,1.0e30,XloopT2);
    spline(qTab,YloopT,nqfunctionsTable,1.0e30,1.0e30,YloopT2);
    spline(qTab,VT,nqfunctionsTable,1.0e30,1.0e30,VT2);
    spline(qTab,TT,nqfunctionsTable,1.0e30,1.0e30,TT2);
    spline(qTab,X10T,nqfunctionsTable,1.0e30,1.0e30,X10T2);
    spline(qTab,Y10T,nqfunctionsTable,1.0e30,1.0e30,Y10T2);
    spline(qTab,ULT,nqfunctionsTable,1.0e30,1.0e30,ULT2);
    spline(qTab,UloopT,nqfunctionsTable,1.0e30,1.0e30,UloopT2);
    spline(qTab,U11T,nqfunctionsTable,1.0e30,1.0e30,U11T2);
    spline(qTab,U20T,nqfunctionsTable,1.0e30,1.0e30,U20T2);
    spline(qTab,xiLT,nqfunctionsTable,1.0e30,1.0e30,xiLT2);
    spline(qTab,LapxiT,nqfunctionsTable,1.0e30,1.0e30,LapxiT2);
    spline(qTab,nabla4xiT,nqfunctionsTable,1.0e30,1.0e30,nabla4xiT2);
    
    fprintf(stdout,"\nTesting Nr values from rmin to rmax to compute CLPT: %d %g %g\n",
            cmd.Nr, cmd.rmin, cmd.rmax);
    if (cmd.Nr==1) {
        dr = 0.;
    } else
        dr = (cmd.rmax - cmd.rmin)/((real)(cmd.Nr - 1));
    
    fprintf(stdout,"\nWriting CLPT correlation functions to file %s...",
            gd.fpfnameclptfunctions);
    outstr = stropen(gd.fpfnameclptfunctions,"w!");

    fprintf(outstr,"%1s%5s%13s%12s%10s%11s%12s%12s%10s%12s%10s%11s%11s%11s%12s%11s",
            "#","r","<xiL>","<xiZA>","<xiA>",
            "<xiW>","<xi10L>","<xi10loop>","<xi20L>",
            "<xi20loop>","<xi01>","<xi02>","<xi11>",
            "<Lapxi>","<nabla4xi>",
            "<xiCLPT>\n");
    
    fprintf(outstr,
            "%1s%6s%11s%11s%11s%11s%11s%11s%11s%11s%12s%11s%11s%11s%11s%11s",
            "#","<1>","<2>","<3>","<4>",
            "<5>","<6>","<7>","<8>",
            "<9>","<10>","<11>","<12>",
            "<13>","<14>",
            "<15>\n");

    for (i=1; i<=cmd.Nr; i++) {
        aTime = second();
        ri = cmd.rmin + dr*((real)(i - 1));
        zacorrfun = zacorrelation_functions(ri);
        clptcorrfun = clptcorrelation_functions(ri);
        
        fprintf(outstr,FMTCORRELATIONFUNCTIONSDAT,
                zacorrfun.r,
                xiLF(ri),
                zacorrfun.xi,
                clptcorrfun.xiA,
                clptcorrfun.xiW,
                clptcorrfun.xi10L,
                clptcorrfun.xi10loop,
                clptcorrfun.xi20L,
                clptcorrfun.xi20loop,
                clptcorrfun.xi01,
                clptcorrfun.xi02,
                clptcorrfun.xi11,
                clptcorrfun.Lapxi,
                clptcorrfun.nabla4xi,
                zacorrfun.xi+clptcorrfun.xiA+clptcorrfun.xiW
                );
    }
    fclose(outstr);
    fprintf(stdout," finished writing.\n");
    
    free_dvector(nabla4xiT2,1,nqfunctionsTable);
    free_dvector(nabla4xiT,1,nqfunctionsTable);
    free_dvector(LapxiT2,1,nqfunctionsTable);
    free_dvector(LapxiT,1,nqfunctionsTable);
    free_dvector(xiLT2,1,nqfunctionsTable);
    free_dvector(xiLT,1,nqfunctionsTable);
    free_dvector(U20T2,1,nqfunctionsTable);
    free_dvector(U20T,1,nqfunctionsTable);
    free_dvector(U11T2,1,nqfunctionsTable);
    free_dvector(U11T,1,nqfunctionsTable);
    free_dvector(UloopT2,1,nqfunctionsTable);
    free_dvector(UloopT,1,nqfunctionsTable);
    free_dvector(ULT2,1,nqfunctionsTable);
    free_dvector(ULT,1,nqfunctionsTable);
    free_dvector(Y10T2,1,nqfunctionsTable);
    free_dvector(Y10T,1,nqfunctionsTable);
    free_dvector(X10T2,1,nqfunctionsTable);
    free_dvector(X10T,1,nqfunctionsTable);
    free_dvector(TT2,1,nqfunctionsTable);
    free_dvector(TT,1,nqfunctionsTable);
    free_dvector(VT2,1,nqfunctionsTable);
    free_dvector(VT,1,nqfunctionsTable);
    free_dvector(YloopT2,1,nqfunctionsTable);
    free_dvector(YloopT,1,nqfunctionsTable);
    free_dvector(XloopT2,1,nqfunctionsTable);
    free_dvector(XloopT,1,nqfunctionsTable);
    free_dvector(YLT2,1,nqfunctionsTable);
    free_dvector(YLT,1,nqfunctionsTable);
    free_dvector(XLT2,1,nqfunctionsTable);
    free_dvector(XLT,1,nqfunctionsTable);
    free_dvector(qTab,1,nqfunctionsTable);
    free(Pqfunctab);
}
#undef FMTCORRELATIONFUNCTIONSDAT

#define FMTQFUNCTIONSDAT    \
"%e %e %e %e %e %e %e %e \
%e %e %e %e %e %e %e \
%e\n"

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

    if (nQsRsTable<=1) {
        error("\nNumber of k's in kfunctions table must be greater than 1...\n\n",nQsRsTable);
    }

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

    fprintf(outstr,"%1s%5s%12s%11s%13s%11s%9s%11s%12s%11s%11s%13s%8s%11s%11s%13s%13s",
            "#","q","<XL>","<YL>","<Xloop>",
            "<Yloop>","<VT>","<TT>","<X10>",
            "<Y10>","<U10L>","<U10loop>","<U11>",
            "<U20>","<xi>","<Lapxi>",
            "<nabla4xi>\n");
    
    fprintf(outstr,
            "%1s%6s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%12s",
            "#","<1>","<2>","<3>","<4>",
            "<5>","<6>","<7>","<8>",
            "<9>","<10>","<11>","<12>",
            "<13>","<14>","<15>",
            "<16>\n");

    for (i=1; i<=Nq; i++) {
        aTime = second();
        qval = rlog10(qmin) + dq*((real)(i - 1));
        qi = rpow(10.0,qval);
        fflush(stdout);
        qfun = qfunctions(qi);
        corrfun = correlation_functions(qi);

        fprintf(outstr,FMTQFUNCTIONSDAT,
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
#undef FMTQFUNCTIONSDAT


//#define FMTBIASTERMDAT    \
//"%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n"

#define FMTBIASTERMDAT    \
"%e %e %e %e %e %e %e %e %e %e\n"

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

    if (nQsRsTable<=1) {
        error("\nNumber of k's in kfunctions table must be greater than 1...\n\n",nQsRsTable);
    }

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
    fprintf(stdout,"\nsigma2L = %g\n",sigma2L);
    
// a02 offset
    a02Off = (1./2.)*Q13QsRs(PQsRstab + 2);

    fprintf(stdout,"\nWriting bias terms and the power spectrum to file %s...",gd.fpfnameSPTPowerSpectrum);
    outstr = stropen(gd.fpfnameSPTPowerSpectrum,"w!");

    fprintf(outstr,"%1s%5s%13s%11s%11s%11s%11s%11s%11s%11s%13s",
            "#","k","<PSLT>","<P22>","<P13>",
            "<a10>","<a01>","<a20>","<a11>","<a02>",
            "<Ploop>\n");
    
    fprintf(outstr,
            "%1s%6s%11s%11s%11s%11s%11s%11s%11s%11s%12s",
            "#","<1>","<2>","<3>","<4>","<5>","<6>","<7>","<8>","<9>","<10>\n");


    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {

        k = kQsRs(p);
        PSLT = PSMGLQsRs(p);
        P22 = (9./98.)*Q1QsRs(p) + (3./7.)*Q2QsRs(p) + (1./2.)*Q3QsRs(p);
        
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
                k, PSLT, P22, P13, a10, a01, a20, a11, a02, Ploop
                );
    }
    fclose(outstr);
    fprintf(stdout," finished writing.\n");
//
    free_dvector(PSLMGT2,1,nQsRsTable);
    free_dvector(PSLMGT,1,nQsRsTable);
    free_dvector(kTab,1,nQsRsTable);
    free(PQsRstab);
}
#undef FMTBIASTERMDAT

local void InputqfunctionsTable(void)
{
    pointqfunctionsTableptr p;
    int i;
    
    fprintf(gd.outlog,"\n\nReading XL from file %s...",gd.fpfnameqfunctions);
    inout_InputData(gd.fpfnameqfunctions, 1, 2, &nqfunctionsTable);
    
    if (nqfunctionsTable < 1)
        error("\n\nInputqfunctionsTable: nqfunctionsTable = %d is absurd\n\n", nqfunctionsTable);
    
    Pqfunctab = (pointqfunctionsTableptr) allocate(nqfunctionsTable * sizeof(pointqfunctionsTable));
    
    i = 0;
    for (p=Pqfunctab; p<Pqfunctab+nqfunctionsTable; p++) {
        qqfuncs(p) = inout_xval[i];
        XLqfuncs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading YL from file %s...",gd.fpfnameqfunctions);
    inout_InputData(gd.fpfnameqfunctions, 1, 3, &nqfunctionsTable);
    
    if (nqfunctionsTable < 1)
        error("\n\nInputqfunctionsTable: nqfunctionsTable = %d is absurd\n\n", nqfunctionsTable);
    
    i = 0;
    for (p=Pqfunctab; p<Pqfunctab+nqfunctionsTable; p++) {
        YLqfuncs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading Xloop from file %s...",gd.fpfnameqfunctions);
    inout_InputData(gd.fpfnameqfunctions, 1, 4, &nqfunctionsTable);
    
    if (nqfunctionsTable < 1)
        error("\n\nInputqfunctionsTable: nqfunctionsTable = %d is absurd\n\n", nqfunctionsTable);
    
    i = 0;
    for (p=Pqfunctab; p<Pqfunctab+nqfunctionsTable; p++) {
        Xloopqfuncs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading Yloop from file %s...",gd.fpfnameqfunctions);
    inout_InputData(gd.fpfnameqfunctions, 1, 5, &nqfunctionsTable);
    
    if (nqfunctionsTable < 1)
        error("\n\nInputqfunctionsTable: nqfunctionsTable = %d is absurd\n\n", nqfunctionsTable);
    
    i = 0;
    for (p=Pqfunctab; p<Pqfunctab+nqfunctionsTable; p++) {
        Yloopqfuncs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading V from file %s...",gd.fpfnameqfunctions);
    inout_InputData(gd.fpfnameqfunctions, 1, 6, &nqfunctionsTable);
    
    if (nqfunctionsTable < 1)
        error("\n\nInputqfunctionsTable: nqfunctionsTable = %d is absurd\n\n", nqfunctionsTable);
    
    i = 0;
    for (p=Pqfunctab; p<Pqfunctab+nqfunctionsTable; p++) {
        Vqfuncs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading T from file %s...",gd.fpfnameqfunctions);
    inout_InputData(gd.fpfnameqfunctions, 1, 7, &nqfunctionsTable);
    
    if (nqfunctionsTable < 1)
        error("\n\nInputqfunctionsTable: nqfunctionsTable = %d is absurd\n\n", nqfunctionsTable);
    
    i = 0;
    for (p=Pqfunctab; p<Pqfunctab+nqfunctionsTable; p++) {
        Tqfuncs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading X10 from file %s...",gd.fpfnameqfunctions);
    inout_InputData(gd.fpfnameqfunctions, 1, 8, &nqfunctionsTable);
    
    if (nqfunctionsTable < 1)
        error("\n\nInputqfunctionsTable: nqfunctionsTable = %d is absurd\n\n", nqfunctionsTable);
    
    i = 0;
    for (p=Pqfunctab; p<Pqfunctab+nqfunctionsTable; p++) {
        X10qfuncs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading Y10 from file %s...",gd.fpfnameqfunctions);
    inout_InputData(gd.fpfnameqfunctions, 1, 9, &nqfunctionsTable);
    
    if (nqfunctionsTable < 1)
        error("\n\nInputqfunctionsTable: nqfunctionsTable = %d is absurd\n\n", nqfunctionsTable);
    
    i = 0;
    for (p=Pqfunctab; p<Pqfunctab+nqfunctionsTable; p++) {
        Y10qfuncs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading UL from file %s...",gd.fpfnameqfunctions);
    inout_InputData(gd.fpfnameqfunctions, 1, 10, &nqfunctionsTable);
    
    if (nqfunctionsTable < 1)
        error("\n\nInputqfunctionsTable: nqfunctionsTable = %d is absurd\n\n", nqfunctionsTable);
    
    i = 0;
    for (p=Pqfunctab; p<Pqfunctab+nqfunctionsTable; p++) {
        ULqfuncs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading Uloop from file %s...",gd.fpfnameqfunctions);
    inout_InputData(gd.fpfnameqfunctions, 1, 11, &nqfunctionsTable);
    
    if (nqfunctionsTable < 1)
        error("\n\nInputqfunctionsTable: nqfunctionsTable = %d is absurd\n\n", nqfunctionsTable);
    
    i = 0;
    for (p=Pqfunctab; p<Pqfunctab+nqfunctionsTable; p++) {
        Uloopqfuncs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading U11 from file %s...",gd.fpfnameqfunctions);
    inout_InputData(gd.fpfnameqfunctions, 1, 12, &nqfunctionsTable);
    
    if (nqfunctionsTable < 1)
        error("\n\nInputqfunctionsTable: nqfunctionsTable = %d is absurd\n\n", nqfunctionsTable);
    
    i = 0;
    for (p=Pqfunctab; p<Pqfunctab+nqfunctionsTable; p++) {
        U11qfuncs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading U20 from file %s...",gd.fpfnameqfunctions);
    inout_InputData(gd.fpfnameqfunctions, 1, 13, &nqfunctionsTable);
    
    if (nqfunctionsTable < 1)
        error("\n\nInputqfunctionsTable: nqfunctionsTable = %d is absurd\n\n", nqfunctionsTable);
    
    i = 0;
    for (p=Pqfunctab; p<Pqfunctab+nqfunctionsTable; p++) {
        U20qfuncs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading xiL from file %s...",gd.fpfnameqfunctions);
    inout_InputData(gd.fpfnameqfunctions, 1, 14, &nqfunctionsTable);
    
    if (nqfunctionsTable < 1)
        error("\n\nInputqfunctionsTable: nqfunctionsTable = %d is absurd\n\n", nqfunctionsTable);
    
    i = 0;
    for (p=Pqfunctab; p<Pqfunctab+nqfunctionsTable; p++) {
        xiLqfuncs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading Lapxi from file %s...",gd.fpfnameqfunctions);
    inout_InputData(gd.fpfnameqfunctions, 1, 15, &nqfunctionsTable);
    
    if (nqfunctionsTable < 1)
        error("\n\nInputqfunctionsTable: nqfunctionsTable = %d is absurd\n\n", nqfunctionsTable);
    
    i = 0;
    for (p=Pqfunctab; p<Pqfunctab+nqfunctionsTable; p++) {
        Lapxiqfuncs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading nabla4xi from file %s...",gd.fpfnameqfunctions);
    inout_InputData(gd.fpfnameqfunctions, 1, 16, &nqfunctionsTable);
    
    if (nqfunctionsTable < 1)
        error("\n\nInputqfunctionsTable: nqfunctionsTable = %d is absurd\n\n", nqfunctionsTable);
    
    i = 0;
    for (p=Pqfunctab; p<Pqfunctab+nqfunctionsTable; p++) {
        nabla4xiqfuncs(p) = inout_yval[i];
        ++i;
    }
}



local void InputQsRsTable(void)
{
    pointQsRsTableptr p;
    int i;
//
    fprintf(gd.outlog,"\n\nReading solution Q1 from file %s...",gd.fpfnamekfun);
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
    fprintf(gd.outlog,"\nReading solution Q2 from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 3, &nQsRsTable);

    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q2QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution Q3 from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 4, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q3QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution Q5 from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 5, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q5QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution Q7 from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 6, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q7QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution Q8 from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 7, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q8QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution Q9 from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 8, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q9QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution Q11 from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 9, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q11QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution Q12 from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 10, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q12QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution Q13 from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 11, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q13QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution QI from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 12, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        QIQsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution R1 from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 13, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        R1QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution R2 from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 14, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        R2QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution R1plus2 from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 15, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        R1plus2QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution RI from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 16, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        RIQsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution Dpk from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 17, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        DpkQsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution PSMGL from file %s...",gd.fpfnamekfun);
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


// BEGIN :: CLPT correlation auxiliary functions

local real XLF(real q)
{
    real func;
    func = Interpolation_nr(q, qTab, XLT, nqfunctionsTable, XLT2);
    return (func);
}

local real YLF(real q)
{
    real func;
    func = Interpolation_nr(q, qTab, YLT, nqfunctionsTable, YLT2);
    return (func);
}

local real XloopF(real q)
{
    real func;
    func = Interpolation_nr(q, qTab, XloopT, nqfunctionsTable, XloopT2);
    return (func);
}

local real YloopF(real q)
{
    real func;
    func = Interpolation_nr(q, qTab, YloopT, nqfunctionsTable, YloopT2);
    return (func);
}

local real VF(real q)
{
    real func;
    func = Interpolation_nr(q, qTab, VT, nqfunctionsTable, VT2);
    return (func);
}

local real TF(real q)
{
    real func;
    func = Interpolation_nr(q, qTab, TT, nqfunctionsTable, TT2);
    return (func);
}

local real X10F(real q)
{
    real func;
    func = Interpolation_nr(q, qTab, X10T, nqfunctionsTable, X10T2);
    return (func);
}

local real Y10F(real q)
{
    real func;
    func = Interpolation_nr(q, qTab, Y10T, nqfunctionsTable, Y10T2);
    return (func);
}

local real ULF(real q)
{
    real func;
    func = Interpolation_nr(q, qTab, ULT, nqfunctionsTable, ULT2);
    return (func);
}

local real UloopF(real q)
{
    real func;
    func = Interpolation_nr(q, qTab, UloopT, nqfunctionsTable, UloopT2);
    return (func);
}

local real U11F(real q)
{
    real func;
    func = Interpolation_nr(q, qTab, U11T, nqfunctionsTable, U11T2);
    return (func);
}

local real U20F(real q)
{
    real func;
    func = Interpolation_nr(q, qTab, U20T, nqfunctionsTable, U20T2);
    return (func);
}

local real xiLF(real q)
{
    real func;
    func = Interpolation_nr(q, qTab, xiLT, nqfunctionsTable, xiLT2);
    return (func);
}

local real LapxiF(real q)
{
    real func;
    func = Interpolation_nr(q, qTab, LapxiT, nqfunctionsTable, LapxiT2);
    return (func);
}

local real nabla4xiF(real q)
{
    real func;
    func = Interpolation_nr(q, qTab, nabla4xiT, nqfunctionsTable, nabla4xiT2);
    return (func);
}

// DERIVED FUNCTIONS:

local real UF(real q)
{
    real func;
    func = ULF(q) + UloopF(q);
    return (func);
}

local real fXF(real q)
{
    real func;
    func = 1.0/XLF(q);
    return (func);
}

local real hYF(real q)
{
    real func;
    func = -YLF(q)/(XLF(q)*XLF(q)+ XLF(q)*YLF(q));
    return (func);
}

local real BDetF(real ff, real hh)
{
    real func;
    func = ff*rsqrt(ff+hh);
    return (func);
}

local real MZAF(real q, real r, real mu)
{
    real func;
    func = preMZAF(q, r, mu, fXF(q), hYF(q));
    return (func);
}

local real preMZAF(real q, real r, real mu, real fx, real hy)
{
    real func;
    func = rsqr(q)*INVSQRTDTWOPI*BDetF(fx,hy)
    * rexp(
           -0.5*(
                 hy* rsqr(q - mu*r) + fx*(rsqr(q) - 2.0*mu*q*r + rsqr(r))
                 )
           );
    return (func);
}

local real AijGijF(real q, real r, real mu, real ft, real ht, real Xloopt, real Yloopt)
{
    real func;
    func = -(rsqr(ft)*(rsqr(q) - 2.0*mu*q*r + rsqr(r)) + ht*(-1.0 + ht*rsqr(q - mu*r))
             + ft*(-3.0 + 2.0*ht*rsqr(q - mu*r)))*Xloopt
    - (ft + ht)*(-1.0 + (ft + ht)*rsqr(q - mu*r) )*Yloopt;
    return (func);
}

local real WijkGammaijkF(real q, real r, real mu, real ft, real ht, real Vt, real Tt)
{
    real func;
    func = (ft + ht)*(-q + mu*r)
    *((ft + ht)*(-3.0 + (ft + ht)*rsqr(q - mu*r))*Tt
      + 3.0*(rsqr(ft)*(rsqr(q) - 2.0*mu*q*r + rsqr(r)) + ht*(-3.0 + ht*rsqr(q - mu*r))
             + ft*(-5.0 + 2.0*ht*rsqr(q - mu*r)))*Vt);
    return (func);
}

local real qigiF(real q, real r, real mu, real ft, real ht)
{
    real func;
    func = (ft + ht)*(q - mu*r);
    return (func);
}

local real qiqjGijF(real q, real r, real mu, real ft, real ht)
{
    real func;
    func = -(ft + ht) *(-1.0 + (ft + ht) *rsqr(q - mu*r));
    return (func);
}

local real MUF(real q, real r, real mu)
{
    real func;
    func = -2.0*MZAF(q, r, mu) * UF(q) * qigiF(q, r, mu, fXF(q), hYF(q));
    return (func);
}

local real MA10F(real q, real r, real mu)
{
    real func;
    func = -MZAF(q, r, mu) *AijGijF(q, r, mu, fXF(q), hYF(q), X10F(q), Y10F(q));
    return (func);
}

local real MxiRF(real q, real r, real mu)
{
    real func;
    func = MZAF(q, r, mu) *xiLF(q);
    return (func);
}

local real MLapxiRF(real q, real r, real mu)
{
    real func;
    func = MZAF(q, r, mu) *LapxiF(q);
    return (func);
}

local real Mnabla4xiRF(real q, real r, real mu)
{
    real func;
    func = MZAF(q, r, mu) *nabla4xiF(q);
    return (func);
}

local real MUUF(real q, real r, real mu)
{
    real func;
    func = -MZAF(q, r, mu) *ULF(q) *ULF(q) *qiqjGijF(q, r, mu, fXF(q), hYF(q));
    return (func);
}

local real MU11F(real q, real r, real mu)
{
    real func;
    func = -MZAF(q, r, mu) *U11F(q) *qigiF(q, r, mu, fXF(q), hYF(q));
    return (func);
}

local real MU20F(real q, real r, real mu)
{
    real func;
    func = -MZAF(q, r, mu) *U20F(q) *qigiF(q, r, mu, fXF(q), hYF(q));
    return (func);
}

local real MxiR2F(real q, real r, real mu)
{
    real func;
    func = 0.5 *MZAF(q, r, mu) *xiLF(q) *xiLF(q);
    return (func);
}

local real MUxiRF(real q, real r, real mu)
{
    real func;
    func = -2.0 *MZAF(q, r, mu) *ULF(q) *qigiF(q, r, mu, fXF(q), hYF(q)) *xiLF(q);
    return (func);
}

local real MAF(real q, real r, real mu)
{
    real func;
    func = -0.5*MZAF(q, r, mu)
    *AijGijF(q, r, mu, fXF(q), hYF(q), XloopF(q), YloopF(q));
    return (func);
}

local real MWF(real q, real r, real mu)
{
    real func;
    func = +(1./6.)*MZAF(q, r, mu)
    *WijkGammaijkF(q, r, mu, fXF(q), hYF(q), VF(q), TF(q));
    return (func);
}


local real M10F(real q, real r, real mu)
{
    real func;
    func = MUF(q, r, mu) + MA10F(q, r, mu);
    return (func);
}

local real M10LF(real q, real r, real mu)
{
    real func;
    func = MUF(q, r, mu);
    return (func);
}

local real M10loopF(real q, real r, real mu)
{
    real func;
    func = MA10F(q, r, mu);
    return (func);
}

local real M01F(real q, real r, real mu)
{
    real func;
    func = MUUF(q, r, mu) + MU20F(q, r, mu);
    return (func);
}

local real M20F(real q, real r, real mu)
{
    real func;
    func = MxiRF(q, r, mu) + MUUF(q, r, mu) + MU11F(q, r, mu);
    return (func);
}

local real M20LF(real q, real r, real mu)
{
    real func;
    func = MxiRF(q, r, mu);
    return (func);
}

local real M20loopF(real q, real r, real mu)
{
    real func;
    func = MUUF(q, r, mu) + MU11F(q, r, mu);
    return (func);
}

local real M02F(real q, real r, real mu)
{
    real func;
    func = MxiR2F(q, r, mu);
    return (func);
}

local real M11F(real q, real r, real mu)
{
    real func;
    func = MUxiRF(q, r, mu);  
    return (func);
}

// correlation functions

global_zacorrfunctions zacorrelation_functions(real ri)
{
    int i, j;
//
    real *muGL, *wGL;
    int Nx;
    real qmin, qmax, dq, q;
    int Nq;
    real xip, xiaA, xiaB;
    real mu, w;
    real mza;
//
    global_zacorrfunctions_ptr zacorrfuncp;
    
    zacorrfuncp = (global_zacorrfunctions_ptr) allocate(1 * sizeof(global_zacorrfunctions));
    
    qmin = 1;
    qmax = 250.0;
    Nq = 250;
    if (Nq==1)
        dq = 0.;
    else
        dq = (qmax - qmin)/((real)(Nq - 1));
    
    Nx=128;
    muGL=dvector(1,Nx);
    wGL=dvector(1,Nx);
    gauleg(-1.0,1.0,muGL,wGL,Nx);
    
    xip = 0.0;
    xiaA = 0.0;
    xiaB = 0.0;
    //
    for (i=1; i<Nq; i++) {
        q = qmin + dq*((real)(i - 1));
        for (j=1; j<=Nx; j++) {
            mu = muGL[j];
            w = wGL[j];
            mza = MZAF(q,ri,mu);
            xiaB += wGL[j]*mza;
        }
//
        xip += dq*(xiaA + xiaB)/2.0;
        xiaA = xiaB;
        xiaB = 0.0;
//
    }
    xip -= 1.0;
//
    rzacorrfun(zacorrfuncp)    = ri;
    xizacorrfun(zacorrfuncp)      = xip;
    
    free_dvector(wGL,1,Nx);
    free_dvector(muGL,1,Nx);
    
    return *zacorrfuncp;
}


global_clptcorrfunctions clptcorrelation_functions(real ri)
{
    int i, j;
//
    real *muGL, *wGL;
    int Nx;
    real qmin, qmax, dq, q;
    int Nq;
    real xiAp, xiAA, xiAB;
    real xiWp, xiWA, xiWB;
    real xi10Lp, xi10LA, xi10LB;
    real xi10loopp, xi10loopA, xi10loopB;
    real xi20Lp, xi20LA, xi20LB;
    real xi20loopp, xi20loopA, xi20loopB;
    real xi01p, xi01A, xi01B;
    real xi02p, xi02A, xi02B;
    real xi11p, xi11A, xi11B;
    real Lapxip, LapxiA, LapxiB;
    real nabla4xip, nabla4xiA, nabla4xiB;
    real mu, w;
    real mA, mW, m10L, m10loop, m20L, m20loop, m01, m02, m11, mLapxi, mnabla4xi;
    
    global_clptcorrfunctions_ptr clptcorrfunp;
    
    clptcorrfunp = (global_clptcorrfunctions_ptr) allocate(1 * sizeof(global_clptcorrfunctions));
    
    qmin = 3;
    qmax = 240.0;
    Nq = 240;
    if (Nq==1)
        dq = 0.;
    else
        dq = (qmax - qmin)/((real)(Nq - 1));
    
    Nx=64;
    muGL=dvector(1,Nx);
    wGL=dvector(1,Nx);
    gauleg(-1.0,1.0,muGL,wGL,Nx);
    
    xiAp = 0.0; xiAA = 0.0; xiAB = 0.0;
    xiWp = 0.0; xiWA = 0.0; xiWB = 0.0;
    xi10Lp = 0.0; xi10LA = 0.0; xi10LB = 0.0;
    xi10loopp = 0.0; xi10loopA = 0.0; xi10loopB = 0.0;
    xi20Lp = 0.0; xi20LA = 0.0; xi20LB = 0.0;
    xi20loopp = 0.0; xi20loopA = 0.0; xi20loopB = 0.0;
    xi01p = 0.0; xi01A = 0.0; xi01B = 0.0;
    xi02p = 0.0; xi02A = 0.0; xi02B = 0.0;
    xi11p = 0.0; xi11A = 0.0; xi11B = 0.0;
    Lapxip = 0.0; LapxiA = 0.0; LapxiB = 0.0;
    nabla4xip = 0.0; nabla4xiA = 0.0; nabla4xiB = 0.0;
//
    for (i=1; i<Nq; i++) {
        q = qmin + dq*((real)(i - 1));
        for (j=1; j<=Nx; j++) {
            mu = muGL[j];
            w = wGL[j];
            mA = MAF(q,ri,mu);
            mW = MWF(q,ri,mu);
            m10L = M10LF(q,ri,mu);
            m10loop = M10loopF(q,ri,mu);
            m20L = M20LF(q,ri,mu);
            m20loop = M20loopF(q,ri,mu);
            m01 = M01F(q,ri,mu);
            m02 = M02F(q,ri,mu);
            m11 = M11F(q,ri,mu);
            mLapxi = MLapxiRF(q,ri,mu);
            mnabla4xi = Mnabla4xiRF(q,ri,mu);
            xiAB += wGL[j]*mA;
            xiWB += wGL[j]*mW;
            xi10LB += wGL[j]*m10L;
            xi10loopB += wGL[j]*m10loop;
            xi20LB += wGL[j]*m20L;
            xi20loopB += wGL[j]*m20loop;
            xi01B += wGL[j]*m01;
            xi02B += wGL[j]*m02;
            xi11B += wGL[j]*m11;
            LapxiB += wGL[j]*mLapxi;
            nabla4xiB += wGL[j]*mnabla4xi;
        }
//
        xiAp += dq*(xiAA + xiAB)/2.0; xiAA = xiAB; xiAB = 0.0;
        xiWp += dq*(xiWA + xiWB)/2.0; xiWA = xiWB; xiWB = 0.0;
        xi10Lp += dq*(xi10LA + xi10LB)/2.0; xi10LA = xi10LB; xi10LB = 0.0;
        xi10loopp += dq*(xi10loopA + xi10loopB)/2.0; xi10loopA = xi10loopB; xi10loopB = 0.0;
        xi20Lp += dq*(xi20LA + xi20LB)/2.0; xi20LA = xi20LB; xi20LB = 0.0;
        xi20loopp += dq*(xi20loopA + xi20loopB)/2.0; xi20loopA = xi20loopB; xi20loopB = 0.0;
        xi01p += dq*(xi01A + xi01B)/2.0; xi01A = xi01B; xi01B = 0.0;
        xi02p += dq*(xi02A + xi02B)/2.0; xi02A = xi02B; xi02B = 0.0;
        xi11p += dq*(xi11A + xi11B)/2.0; xi11A = xi11B; xi11B = 0.0;
        Lapxip += dq*(LapxiA + LapxiB)/2.0; LapxiA = LapxiB; LapxiB = 0.0;
        nabla4xip += dq*(nabla4xiA + nabla4xiB)/2.0; nabla4xiA = nabla4xiB; nabla4xiB = 0.0;
//
    }
    
    rclptcorrfun(clptcorrfunp) = ri;
    xiAclptcorrfun(clptcorrfunp) = xiAp;
    xiWclptcorrfun(clptcorrfunp) = xiWp;
    xi10Lclptcorrfun(clptcorrfunp) = xi10Lp;
    xi10loopclptcorrfun(clptcorrfunp) = xi10loopp;
    xi20Lclptcorrfun(clptcorrfunp) = xi20Lp;
    xi20loopclptcorrfun(clptcorrfunp) = xi20loopp;
    xi01clptcorrfun(clptcorrfunp) = xi01p;
    xi02clptcorrfun(clptcorrfunp) = xi02p;
    xi11clptcorrfun(clptcorrfunp) = xi11p;
    Lapxiclptcorrfun(clptcorrfunp) = Lapxip;
    nabla4xiclptcorrfun(clptcorrfunp) = nabla4xip;
    
    free_dvector(wGL,1,Nx);
    free_dvector(muGL,1,Nx);
    
    return *clptcorrfunp;
}

// END :: CLPT correlation auxiliary functions


//
// AUXILIARY FUNCTIONS FOR q FUNCTIONS COMPUTATION
//

local real tildeV(real k)
{
    real func;
    
    func= -(3./35.)*( QIF(k) - 3.*Q2F(k) + 2.*RIF(k) - 6.*R2F(k) );
    
    return (func);
}

local real tildeT(real k)
{
    real func;
    
    func= -(9./14.)*(QIF(k) + 2.*Q2F(k) + 2.*RIF(k) + 4.*R2F(k) );
    
    return (func);
}

local real xL(real k, real q)
{
    real func;
    
    func= PSLF(k)*(1./3. - rj1Bessel(k*q)/(k*q) );
    
    return (func);
}

local real xloop(real k, real q)
{
    real func;
    
    func= ( (9./98.)*Q1F(k) + (10./21.)*R1F(k) )*(1./3. - rj1Bessel(k*q)/(k*q) );
    
    return (func);
}

local real yL(real k, real q)
{
    real func;
    
    func= PSLF(k)*rj2Bessel(k*q);
    
    return (func);
}

local real yloop(real k, real q)
{
    real func;
    
    func= ( (9./98.)*Q1F(k) + (10./21.)*R1F(k) )*rj2Bessel(k*q);
    
    return (func);
}

local real x10(real k, real q)
{
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
    
//    if ( k < kPS[1] || k > kPS[nPS] )
//        fprintf(gd.outlog,"\n\nInterpolation_nr: warning! :: k is out of range... %g\n",k);
    
    splint(kPS,pPS,pPS2,nPS,k,&psftmp);
    
    return (psftmp);
}


// q functions
global_qfunctions qfunctions(real qi)
{
    global_qfunctions_ptr qfunp;
    int i, Nk;
    real dk, kvali, kvalim1, ki, kim1;
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
    
    Nk = 1200;
    if (Nk==1)
        dk = 0.;
    else
        dk = (rlog10(cmd.kmax) - rlog10(cmd.kmin))/((real)(Nk - 1));
    
    U10Lp = 0.;
    U10LA = -cmd.kmin*PSLF(cmd.kmin)*rj1Bessel(cmd.kmin*qi);
    U10loopp = 0.;
    U10loopA = -cmd.kmin*( (5./21.)*R1F(cmd.kmin) )*rj1Bessel(cmd.kmin*qi);
    U11p = 0.;
    U11A = -cmd.kmin*( (6./7.)*R1plus2F(cmd.kmin) )*rj1Bessel(cmd.kmin*qi);
    U20p = 0.;
    U20A = -cmd.kmin*( (3./7.)*Q8F(cmd.kmin) )*rj1Bessel(cmd.kmin*qi);
//
    XLp = 0.;
    XLA = xL(cmd.kmin, qi);
    Xloopp = 0.;
    XloopA = xloop(cmd.kmin, qi);
    X10p = 0.;
    X10A = x10(cmd.kmin, qi);
    YLp = 0.;
    YLA = yL(cmd.kmin, qi);
    Yloopp = 0.;
    YloopA = yloop(cmd.kmin, qi);
    Y10p = 0.;
    Y10A = y10(cmd.kmin, qi);
//
    preVp = 0.;
    preVA = tildeV(cmd.kmin)*rj1Bessel(cmd.kmin*qi)/cmd.kmin;
    Tp = 0.;
    TA = tildeT(cmd.kmin)*rj3Bessel(cmd.kmin*qi)/cmd.kmin;
//
    for (i=2; i<=Nk; i++) {
        kvali = rlog10(cmd.kmin) + dk*((real)(i - 1));
        kvalim1 = rlog10(cmd.kmin) + dk*((real)(i - 2));
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

    kmin = kPS[1];
    kmax = kPS[nPSLT];
    ymin = rlog10(kmin);
    ymax = rlog10(kmax);

    result= (1.0/SIXPI2)*rlog(10.0)
    *qromo(sigma2L_function_int,ymin,ymax,midpnt,EPSQ,KK);

    return result;

}

#undef EPSQ
#undef KK
