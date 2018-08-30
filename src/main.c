/*==============================================================================
 NAME: main.c				[mgpt]
 Mario A. Rodriguez-Meza (marioalberto.rodriguez@inin.gob.mx)
 Alejandro Aviles (avilescervantes@gmail.com)
 ================================================================================ 
 Use: mgpt -help
 Input: 	Command line parameters or Parameters file
 Output: kfunctions and SPTPowerSpectrum.
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

