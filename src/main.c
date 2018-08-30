/*==============================================================================
 NAME: main.c				[mgpt]
 Mario A. Rodriguez-Meza
 Alejandro Aviles
 ================================================================================ 
 Use: mgpt -help
 Input: 	Command line parameters, Parameters file and/or icfile
 Output: ...
 Units:
 History:
 Comments and notes:
 References:
 ==============================================================================*/

#define global

#include "globaldefs.h"
#include "cmdline_defs.h"
#include "protodefs.h"

int main(int argc, string argv[])
{
    gd.cpuinit = cputime();
    InitParam(argv, defv);
    StartRun(argv[0], HEAD1, HEAD2, HEAD3);
	MainLoop();
	EndRun();
    return 0;
}

