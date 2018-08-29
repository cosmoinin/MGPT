/*==============================================================================
 NAME: main.c				[mgpt]
 Written by: Mario A. Rodriguez-Meza
 Starting date: January 2018
 Purpose: Main routine
 Language: C
 Comments and notes:
 Info: Mario A. Rodriguez-Meza
 Depto. de Fisica, ININ
 Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
 e-mail: marioalberto.rodriguez@inin.gob.mx
 http://www.inin.mx/
 
 Major revision:
 Copyright: (c) 2005-2018 Mar.  All Rights Reserved
 ================================================================================
 
 Use: mgpt -help
 Input: 	Command line parameters, Parameters file and/or icfile
 Output: ...
 Units:
 History:
 Comments and notes:
 References:
 ================================================================================
 Legal matters:
 The author does not warrant that the program and routines it contains
 listed below are free from error or suitable for particular applications,
 and he disclaims all liability from any consequences arising from their use.
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

