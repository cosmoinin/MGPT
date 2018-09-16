/*==============================================================================
 MODULE: template_io.c		[mgpt]
 ==============================================================================*/

#include "globaldefs.h"
#include "protodefs.h"


void StartOutput(void)
{

    fprintf(stdout,"\n  \t -- %s --\n", gd.model_comment);
//
    fprintf(stdout,"  \t -- %s --\n", gd.integration_method_comment);
    fprintf(stdout,"  \t -- %s --\n\n", gd.quadraturemethod_comment);

    fprintf(gd.outlog,"\n%8s%8s%8s", "maxnsteps", "etaini", "deta");
    fprintf(gd.outlog,"%8s\n","etaout");
    fprintf(gd.outlog,"%8d%8.2f%8.4f%8.4f",cmd.maxnsteps,cmd.x,gd.dx,gd.xstop);
    if (! strnull(cmd.options))
        fprintf(stdout,"\n\toptions: %s\n", cmd.options);
    
}

/*
void output(void)
{
    real xeff;
    int j;

	xeff = gd.xnow + gd.dx/8.0;

}
*/

// I/O directories:
global void setFilesDirs_log(void)
{
    char buf[200];
    
    sprintf(gd.tmpDir,"tmp");
    
    sprintf(buf,"if [ ! -d %s ]; then mkdir %s; fi",gd.tmpDir,gd.tmpDir);
    system(buf);
    
    sprintf(gd.logfilePath,"%s/mgpt%s.log",gd.tmpDir,cmd.suffixModel);
}

global void setFilesDirs(void)
{
    char buf[200];
    
    sprintf(gd.clptDir,"CLPT");
    sprintf(buf,"if [ ! -d %s ]; then mkdir %s; fi",gd.clptDir,gd.clptDir);
    fprintf(gd.outlog,"system: %s\n",buf);
    system(buf);
    
    sprintf(gd.fpfnamekfun,"CLPT/kfunctions%s.dat",cmd.suffixModel);
    sprintf(gd.fpfnameSPTPowerSpectrum,"SPTPowerSpectrum%s.dat",cmd.suffixModel);
    sprintf(gd.fpfnameqfunctions,"CLPT/qfunctions%s.dat",cmd.suffixModel);
    sprintf(gd.fpfnameclptfunctions,"CorrelationFunction%s.dat",cmd.suffixModel);
}

void EndRun(void)
{
    char   buf[200];
	FILE *fd;

	fclose(gd.outlog);
    
    printf("\nFinal CPU time : %g min.\n\n", cputime() - gd.cpuinit);
}



