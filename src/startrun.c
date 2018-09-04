/*==============================================================================
 MODULE: startrun.c				
 ==============================================================================*/

#include "globaldefs.h"
#include "protodefs.h"

local void ReadParameterFile(char *);
local void PrintParameterFile(char *);

local void startrun_parameterfile(void);
local void startrun_cmdline(void);
local void ReadParametersCmdline(void);
local void startrun_Common(void);
local void startrun_ParamStat(void);
local void CheckParameters(void);

local void setFilesDirs(void);

local void InputPSTable(void);
local void PSLTable(void);
local void GaussLegendrePoints(void);

void StartRun(string head0, string head1, string head2, string head3)
{
    real aTime;
    aTime = cputime();

    gd.headline0 = head0; gd.headline1 = head1;
    gd.headline2 = head2; gd.headline3 = head3;
    printf("\n%s\n%s: %s\n\t %s\n",
		gd.headline0, gd.headline1, gd.headline2, gd.headline3);

	gd.stopflag = 0;

    cmd.paramfile = GetParam("paramfile");
    if (!strnull(cmd.paramfile))
		startrun_parameterfile();
	else
		startrun_cmdline();

    StartOutput();

    fprintf(stdout,"\nTime to StartRun %g\n\n",cputime()-aTime);
    fflush(stdout);
}

local void startrun_parameterfile(void)
{
	ReadParameterFile(cmd.paramfile);
    startrun_ParamStat();
	startrun_Common();
	PrintParameterFile(cmd.paramfile);
}

#define parameter_null	"parameters_null-mgpt"

local void startrun_cmdline(void)
{
	ReadParametersCmdline();
	startrun_Common();
	PrintParameterFile(parameter_null);
}

local void ReadParametersCmdline(void)
{
// Modified gravity model parameters:
    cmd.mgmodel = GetParam("mgmodel");
    cmd.nHS = GetiParam("nHS");
    cmd.fR0 = GetdParam("fR0");
    cmd.screening = GetdParam("screening");
//
// Power spectrum table:
    cmd.fnamePS = GetParam("fnamePS");
    cmd.kmin = GetdParam("kmin");
    cmd.kmax = GetdParam("kmax");
    cmd.Nk = GetiParam("Nk");
//
    cmd.om = GetdParam("om");
    cmd.h = GetdParam("h");
//
// Differential equations evolution parameters:
    cmd.x = GetdParam("etaini");
    cmd.dxstr = GetParam("deta");
    cmd.dxmin = GetdParam("detamin");
    cmd.eps = GetdParam("eps");
    cmd.xstop = GetdParam("zout");
    cmd.maxnsteps = GetiParam("maxnsteps");
	cmd.integration_method = GetParam("integration_method");
//
// Integration parameters:
    cmd.ngausslegpoints = GetiParam("ngausslegpoints");
// Post processing parameters:
    cmd.postprocessing = GetbParam("postprocessing");
    cmd.options = GetParam("options");
}

#undef parameter_null

#define logfile			"mgpt.log"

local void startrun_Common(void)
{
	real dx1, dx2;

    setFilesDirs();
    strcpy(gd.mode,"w");
	if(!(gd.outlog=fopen(gd.logfilePath, gd.mode)))
		error("\nstart_Common: error opening file '%s' \n",gd.logfilePath);

//
    gd.dx = (sscanf(cmd.dxstr, "%lf/%lf", &dx1, &dx2) == 2 ?
				dx1/dx2 : atof(cmd.dxstr));
    if ( dx2 == 0. )
        error("\n\nstartrun_Common: dx : dx2 must be finite\n");

    CheckParameters();
    GaussLegendrePoints();
	integration_method_string_to_int(cmd.integration_method, &gd.method_int);
    gd.xnow = cmd.x;
    gd.xout = gd.xnow;
    gd.xoutinfo = gd.xnow;

    gd.xstop = rlog(1.0/(1.0+cmd.xstop));

    set_model();

    if (!strnull(cmd.fnamePS)) {
        sprintf(gd.fnamePSPath,"Input/%s",cmd.fnamePS);
        fprintf(gd.outlog,"PS file Path and file name: %s\n",gd.fnamePSPath);
        fflush(gd.outlog);
        InputPSTable();
        PSLTable();
    }

}

local void startrun_ParamStat(void)
{
	real dx1, dx2;

    if (GetParamStat("fnamePS") & ARGPARAM)
        cmd.fnamePS = GetParam("fnamePS");

	if (GetParamStat("etaini") & ARGPARAM)
		cmd.x = GetdParam("etaini");
	if (GetParamStat("deta") & ARGPARAM) {
		cmd.dxstr = GetParam("deta");
		gd.dx = (sscanf(cmd.dxstr, "%lf/%lf", &dx1, &dx2) == 2 ?
                    dx1/dx2 : atof(cmd.dxstr));
		if ( dx2 == 0. )
			error("\n\nstartrun_ParamStat: deta : deta2 must be finite\n");
	}
	if (GetParamStat("zout") & ARGPARAM)
		cmd.xstop = GetdParam("zout");

	if (GetParamStat("maxnsteps") & ARGPARAM) 
		cmd.maxnsteps = GetiParam("maxnsteps");

	if (GetParamStat("integration_method") & ARGPARAM) {
		cmd.integration_method = GetParam("integration_method");
		fprintf(gd.outlog,"\n\nrunning now %s integration method ...\n",
				cmd.integration_method);
	}

    if (GetParamStat("ngausslegpoints") & ARGPARAM)
        cmd.ngausslegpoints = GetiParam("ngausslegpoints");

    if (GetParamStat("postprocessing") & ARGPARAM)
        cmd.postprocessing = GetbParam("postprocessing");
    
    if (GetParamStat("options") & ARGPARAM)
        cmd.options = GetParam("options");
    
}

local void CheckParameters(void)
{
// Power spectrum table:
    if (strnull(cmd.fnamePS))
        error("CheckParameters: You should give a power spectrum filename\n");
    if (cmd.kmin < 0.0)
        error("CheckParameters: absurd value for kmin\n");
    if (cmd.kmax < 0.0)
        error("CheckParameters: absurd value for kmax\n");
    if (cmd.kmin > cmd.kmax)
        error("CheckParameters: kmin can not be greater than kmax\n");
    if (cmd.Nk < 0)
        error("CheckParameters: absurd value for Nk\n");
//
    if (cmd.om > 1.0 || cmd.om < 0.0)
        error("CheckParameters: absurd value for om\n");
    if (cmd.h < 0.0)
        error("CheckParameters: absurd value for h\n");

    if (gd.dx == 0)
        error("CheckParameters: absurd value for deta\n");
    if(cmd.x == rlog(1.0/(1.0+cmd.xstop)) )
        error("\n\nstartrun_Common: etaini and etaout=exp(-zout)-1 must be different\n");

    if (cmd.maxnsteps < 1)
        error("CheckParameters: absurd value for maxnsteps\n");

    if (cmd.ngausslegpoints <= 1)
        error("CheckParameters: absurd value for ngausslegpoints\n");
}

// I/O directories:
local void setFilesDirs(void)
{
    char buf[200];

    sprintf(buf,"if [ ! -d tmp ]; then mkdir tmp; fi");
//    fprintf(stdout,"system: %s\n",buf);
    system(buf);

    sprintf(gd.logfilePath,"tmp/%s",logfile);
//    fprintf(stdout,"Log file Path and file name: %s\n",gd.logfilePath);
}

#undef logfile

local void ReadParameterFile(char *fname)
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

// Power spectrum table:
    SPName(cmd.fnamePS,"fnamePS",100);
    RPName(cmd.kmin,"kmin");
    RPName(cmd.kmax,"kmax");
    IPName(cmd.Nk,"Nk");
//
// Modified gravity model parameters:
    SPName(cmd.mgmodel,"mgmodel",100);
    IPName(cmd.nHS,"nHS");
    RPName(cmd.fR0,"fR0");
    RPName(cmd.screening,"screening");
//
    RPName(cmd.om,"om");
    RPName(cmd.h,"h");
//
// Differential equations evolution parameters:
    RPName(cmd.x,"etaini");
	SPName(cmd.dxstr,"deta",100);
    RPName(cmd.dxmin,"detamin");
    RPName(cmd.eps,"eps");
	RPName(cmd.xstop,"zout");
    IPName(cmd.maxnsteps,"maxnsteps");
	SPName(cmd.integration_method,"integration_method",100);
//
// Integration parameters:
    IPName(cmd.ngausslegpoints,"ngausslegpoints");
//
// Post processing parameters:
    BPName(cmd.postprocessing,"postprocessing");
    SPName(cmd.options,"options",100);
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

#define FMTT	"%-35s%s\n"
#define FMTI	"%-35s%d\n"
#define FMTR	"%-35s%g\n"

local void PrintParameterFile(char *fname)
{
    FILE *fdout;
    char buf[200];
    
    sprintf(buf,"tmp/%s%s",fname,"-usedvalues");
    if(!(fdout=fopen(buf,"w"))) {
        fprintf(stdout,"error opening file '%s' \n",buf);
        exit(0);
    } else {
        fprintf(fdout,"%s\n",
		"%-------------------------------------------------------------------");
        fprintf(fdout,"%s %s\n","% Parameter input file for:",gd.headline0);
        fprintf(fdout,"%s\n","%");
        fprintf(fdout,"%s %s: %s\n%s\t    %s\n","%",gd.headline1,gd.headline2,"%",
            gd.headline3);
        fprintf(fdout,"%s\n%s\n",
		"%-------------------------------------------------------------------",
		"%");
// Power spectrum table:
        fprintf(fdout,FMTT,"fnamePS",cmd.fnamePS);
        fprintf(fdout,FMTR,"kmin",cmd.kmin);
        fprintf(fdout,FMTR,"kmax",cmd.kmax);
        fprintf(fdout,FMTI,"Nk",cmd.Nk);
//
// Modified gravity model parameters:
        fprintf(fdout,FMTT,"mgmodel",cmd.mgmodel);
        fprintf(fdout,FMTI,"nHS",cmd.nHS);
        fprintf(fdout,FMTR,"fR0",cmd.fR0);
        fprintf(fdout,FMTR,"screening",cmd.screening);
//
        fprintf(fdout,FMTR,"om",cmd.om);
        fprintf(fdout,FMTR,"h",cmd.h);
//
// Differential equations evolution parameters:
        fprintf(fdout,FMTR,"etaini",cmd.x);
        fprintf(fdout,FMTR,"detamin",cmd.dxmin);
        fprintf(fdout,FMTR,"eps",cmd.eps);
        fprintf(fdout,FMTT,"deta",cmd.dxstr);
        fprintf(fdout,FMTR,"zout",cmd.xstop);
        fprintf(fdout,FMTI,"maxnsteps",cmd.maxnsteps);
        fprintf(fdout,FMTT,"integration_method",cmd.integration_method);
//
        fprintf(fdout,FMTI,"ngausslegpoints",cmd.ngausslegpoints);
        fprintf(fdout,FMTT,"postprocessing",cmd.postprocessing ? "true" : "false");
        fprintf(fdout,FMTT,"options",cmd.options);
        fprintf(fdout,"\n\n");
    }
    fclose(fdout);
}

#undef FMTT
#undef FMTI
#undef FMTR


#define NPT 20
#define SPREAD 1.0
local void InputPSTable(void)
{
    stream outstr;
    pointPSTableptr p, plog, pn;
    pointPSTableptr PSLCDMtabtmp;
    int nPSTabletmp;
    int i;
    real dk, kval, PSval, kmin, kmax;
    int mwt;
    double al,bl,chi2,q,siga,sigb,*x,*y,*sig;
    double au, bu;
    real *kPStmp;
    real *pPStmp;
    real *pPS2tmp;
    char namebuf[256];
    real kminext, kmaxext, dktmp, kmn, kmx;
    int Nkext=600, NkL=50, NkU=50;
    real kminT=1.0e-6, kmaxT=350.0;


    fprintf(gd.outlog,"\n\nReading power spectrum from file %s...\n",gd.fnamePSPath);
    inout_InputData(gd.fnamePSPath, 1, 2, &nPSTabletmp);

    if (nPSTabletmp < 1)
        error("\n\nInputPSTable: nPSTable = %d is absurd\n\n", nPSTabletmp);

    PSLCDMtabtmp = (pointPSTableptr) allocate(nPSTabletmp * sizeof(pointPSTable));
    
    fprintf(gd.outlog,"nPSTable : %d\n", nPSTabletmp);

    i = 0;
    for (p=PSLCDMtabtmp; p<PSLCDMtabtmp+nPSTabletmp; p++) {
        kPos(p) = inout_xval[i];
        PS(p) = inout_yval[i];
        ++i;
    }

    fprintf(gd.outlog,"\n\nCreating log power spectrum...\n");

    PSLCDMLogtab = (pointPSTableptr) allocate(nPSTabletmp * sizeof(pointPSTable));
    plog = PSLCDMLogtab;
    nPSLogT=0;
    for (p=PSLCDMtabtmp; p<PSLCDMtabtmp+nPSTabletmp; p++) {
        kPos(plog) = rlog10(kPos(p));
        PS(plog) = rlog10(PS(p));
        plog++;
        nPSLogT++;
    }

    fprintf(gd.outlog,"\nTotal numbers in Log PS: %d %d\n",nPSLogT,plog-PSLCDMLogtab);
    fprintf(gd.outlog,"Total numbers in Normal PS: %d %d\n\n",nPSTabletmp,p-PSLCDMtabtmp);

    fprintf(gd.outlog,"\n\nLinear fit (a + b x) to log-log power spectrum at minset and maxset...\n");
    fprintf(gd.outlog,"\n\nLinear fit to log-log power spectrum at minset and maxset...\n");

    x=dvector(1,NPT);
    y=dvector(1,NPT);
    sig=dvector(1,NPT);
    
// Lower part of the PS
    fprintf(gd.outlog,"\nAt lower part of the spectrum...\n");

    plog = PSLCDMLogtab;
    for (i=1;i<=NPT;i++) {
        x[i]=kPos(plog);
        y[i]=PS(plog);
        sig[i]=SPREAD;
        plog++;
    }
    for (mwt=0;mwt<=1;mwt++) {
        fit(x,y,NPT,sig,mwt,&al,&bl,&siga,&sigb,&chi2,&q);
        if (mwt == 0)
            fprintf(gd.outlog,"\nIgnoring standard deviations\n");
        else
            fprintf(gd.outlog,"\nIncluding standard deviations\n");
        fprintf(gd.outlog,"%12s %9.6f %18s %9.6f \n",
               "a  =  ",al,"uncertainty:",siga);
        fprintf(gd.outlog,"%12s %9.6f %18s %9.6f \n",
               "b  =  ",bl,"uncertainty:",sigb);
        fprintf(gd.outlog,"%19s %14.6f \n","chi-squared: ",chi2);
        fprintf(gd.outlog,"%23s %10.6f \n","goodness-of-fit: ",q);
    }

// Upper part of the PS
    fprintf(gd.outlog,"\nAt upper part of the spectrum...\n");

    plog = PSLCDMLogtab+nPSLogT-1;
    for (i=1;i<=NPT;i++) {
        x[i]=kPos(plog);
        y[i]=PS(plog);
        sig[i]=SPREAD;
        plog--;
    }
    for (mwt=0;mwt<=1;mwt++) {
        fit(x,y,NPT,sig,mwt,&au,&bu,&siga,&sigb,&chi2,&q);
        if (mwt == 0)
            fprintf(gd.outlog,"\nIgnoring standard deviations\n");
        else
            fprintf(gd.outlog,"\nIncluding standard deviations\n");
        fprintf(gd.outlog,"%12s %9.6f %18s %9.6f \n",
               "a  =  ",au,"uncertainty:",siga);
        fprintf(gd.outlog,"%12s %9.6f %18s %9.6f \n",
               "b  =  ",bu,"uncertainty:",sigb);
        fprintf(gd.outlog,"%19s %14.6f \n","chi-squared: ",chi2);
        fprintf(gd.outlog,"%23s %10.6f \n","goodness-of-fit: ",q);
    }

// Extending power spectrum
    kPStmp = dvector(1,nPSLogT);
    pPStmp = dvector(1,nPSLogT);
    pPS2 = dvector(1,nPSLogT);
    pn = PSLCDMLogtab;
    i=1;
    for (pn = PSLCDMLogtab; pn<PSLCDMLogtab+nPSLogT; pn++) {
        kPStmp[i] = kPos(pn);
        pPStmp[i] = PS(pn);
        i++;
    }
//
    spline(kPStmp,pPStmp,nPSLogT,1.0e30,1.0e30,pPS2);

    kmin = kPos(PSLCDMtabtmp);
    kmax = kPos(PSLCDMtabtmp+nPSTabletmp-1);
    fprintf(gd.outlog,"\nkmin, kmax of the given power spectrum (with %d values): %g %g",kmin, kmax, nPSTabletmp);
    dktmp = (rlog10(kmax) - rlog10(kmin))/((real)(nPSTabletmp - 1));
    kminext = rpow(10.0, rlog10(kmin)-((real)NkL)*dktmp);
    kmaxext = rpow(10.0, rlog10(kmax)+((real)NkU)*dktmp);
//
    fprintf(gd.outlog,"\n\nNkL, NkU: %d %d\n",NkL, NkU);
    NkL = rlog10(kmin/kminT)/dktmp;
    NkU = rlog10(kmaxT/kmax)/dktmp;
    fprintf(gd.outlog,"\n\nNkL, NkU targets: %d %d\n",NkL,NkU);
    kminext = rpow(10.0, rlog10(kmin)-((real)NkL)*dktmp);
    kmaxext = rpow(10.0, rlog10(kmax)+((real)NkU)*dktmp);
//

    fprintf(gd.outlog,"\nkmin, kmax of the extended power spectrum (first try): %g %g",
            kminext, kmaxext);

    kmn = MIN(kminext,cmd.kmin);
    kmx = MAX(kmaxext,cmd.kmax);
    fprintf(gd.outlog,"\nkmin, kmax of the extended power spectrum (second try): %g %g\n",
            kmn, kmx);

    nPSTable = Nkext;
    fprintf(gd.outlog,"\n\nCreating new PSTable with %d values\n",nPSTable);
    PSLCDMtab = (pointPSTableptr) allocate(nPSTable * sizeof(pointPSTable));
    dk = (rlog10(kmx) - rlog10(kmn))/((real)(nPSTable - 1));
    p = PSLCDMtab;
    
    for (i=1; i<=nPSTable; i++) {
        kval = rlog10(kmn) + dk*((real)(i - 1));
        if (rpow(10.0,kval) > kmin && rpow(10.0,kval) < kmax)
            PSval = psInterpolation_nr(kval, kPStmp, pPStmp, nPSLogT);
        else
            if (rpow(10.0,kval) < kmin)
                PSval = al + bl*kval;
            else
                if (rpow(10.0,kval) > kmax)
                    PSval = au + bu*kval;
                else
                    error("\n\nError: InputPSTable :: kmin, kmax, kval: %g %g %g", kmin, kmax, kval);
//
        kPos(p) = rpow(10.0,kval);
        PS(p) = rpow(10.0,PSval);
        p++;
    }

    sprintf(namebuf,"tmp/%s_%s", cmd.fnamePS, "ext.dat");
    outstr = stropen(namebuf,"w!");
    for (p=PSLCDMtab; p<PSLCDMtab+nPSTable; p++) {
        fprintf(outstr,"%g %g\n",
                kPos(p),PS(p));
    }
    fclose(outstr);
//

    free_dvector(pPS2,1,nPSLogT);
    free_dvector(pPStmp,1,nPSLogT);
    free_dvector(kPStmp,1,nPSLogT);

    free_dvector(sig,1,NPT);
    free_dvector(y,1,NPT);
    free_dvector(x,1,NPT);
}
#undef NPT
#undef SPREAD


local void PSLTable(void)
{
    char namebuf[256];
    stream outstr;
    real kmin, Dpkmin, Dpk;
    pointPSTableptr p, pn;
    int i;
//
    real xstoptmp, Dp0, Dpzout, fac;

    xstoptmp = gd.xstop;
    gd.xstop = 0.;
    Dp0 = DpFunction(0.);
    fprintf(gd.outlog,"\n\n Dp(0) = %g\n",Dp0);
    gd.xstop = xstoptmp;
    Dpzout = DpFunction(0.);
    fprintf(gd.outlog,"\n\n Dp(zout) = %g\n",Dpzout);
    
    fac = rsqr(Dpzout/Dp0);
    for (p = PSLCDMtab; p<PSLCDMtab+nPSTable; p++) {
        PS(p) *= fac;
    }
    
    sprintf(namebuf,"tmp/%s_%s",cmd.fnamePS,"ext2.dat");
    outstr = stropen(namebuf,"w!");
    for (p=PSLCDMtab; p<PSLCDMtab+nPSTable; p++) {
        fprintf(outstr,"%g %g\n",
                kPos(p),PS(p));
    }
    fclose(outstr);

//

    kmin = kPos(PSLCDMtab);
    Dpkmin = DpFunction(kmin);
    fprintf(gd.outlog,"\n\n Dpkmin = %g\n",Dpkmin);

    PSLT = (pointPSTableptr) allocate(nPSTable * sizeof(pointPSTable));
    nPSLT = 0;
    pn = PSLT;
    for (p = PSLCDMtab; p<PSLCDMtab+nPSTable; p++) {
        kPos(pn) = kPos(p);
        Dpk = DpFunction(kPos(p));
        PS(pn) = rsqr(Dpk/Dpkmin)*PS(p);
        pn++;
        nPSLT++;
    }

    kPS = dvector(1,nPSLT);
    pPS = dvector(1,nPSLT);
    pPS2 = dvector(1,nPSLT);
    pn = PSLT;
    i=1;
    for (pn = PSLT; pn<PSLT+nPSLT; pn++) {
        kPS[i] = kPos(pn);
        pPS[i] = PS(pn);
        i++;
    }
//
    spline(kPS,pPS,nPSLT,1.0e30,1.0e30,pPS2);

}

local void GaussLegendrePoints(void)
{
    double dx,ss,xval;
    real x1=-1.0, x2=1.0;
    
    int i;
    double xx=0.0;
    
    pGL = (global_GL_ptr) allocate(sizeof(global_GL));

    nGL(pGL) = cmd.ngausslegpoints;
    xGL(pGL)=dvector(1,cmd.ngausslegpoints);
    wGL(pGL)=dvector(1,cmd.ngausslegpoints);
    x1GL(pGL) = -1.0;
    x2GL(pGL) = 1.0;

    fprintf(gd.outlog,"\nComputation of GLs...\n");
    gauleg(x1GL(pGL),x2GL(pGL),xGL(pGL),wGL(pGL),nGL(pGL));

    fprintf(gd.outlog,"\nEnd of Gauss-Legendre computing (StartRun)\n");

}

