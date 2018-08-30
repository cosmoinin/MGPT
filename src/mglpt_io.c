/*==============================================================================
 MODULE: template_io.c		[mglpt]
 ==============================================================================*/

#include "globaldefs.h"
#include "protodefs.h"


void StartOutput(void)
{

    fprintf(stdout,"\n\n  \t -- %s --\n", gd.model_comment);
//
    fprintf(stdout,"  \t -- %s --\n\n", gd.integration_method_comment);

    fprintf(gd.outlog,"\n%8s%8s%8s", "maxnsteps", "etaini", "deta");
    fprintf(gd.outlog,"%8s\n","etaout");
    fprintf(gd.outlog,"%8d%8.2f%8.4f%8.4f",cmd.maxnsteps,cmd.x,gd.dx,gd.xstop);
    if (! strnull(cmd.options))
        fprintf(stdout,"\n\toptions: %s\n", cmd.options);
    
}

void output(void)
{
    real xeff;
    int j;

	xeff = gd.xnow + gd.dx/8.0;

}

void EndRun(void)
{
    char   buf[200];
	FILE *fd;

	fclose(gd.outlog);
    
    printf("\nFinal CPU time : %g\n\n", cputime() - gd.cpuinit);
}



