/*==============================================================================
	HEADER: getparam.h			[General_libs]
	Written by: M.A. Rodriguez-Meza
	Starting date: May 2006
	Purpose: definitions and protodefinitions to initialize the parameters code
	Language: C
	Use: '#include "getparam.h"'
	Use in routines and functions:
	External headers:
	Coments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions: January 2007;
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their use.
==============================================================================*/


#ifndef _getparam_h
#define _getparam_h

void InitParam(string *, string *);

string GetParam(string);

int GetiParam(string);
bool GetbParam(string);
double GetdParam(string);

#define getargv0()    (GetParam("argv0"))
#define getversion()  (GetParam("Version"))

int GetParamStat(string);

#define DEFPARAM        001                     
#define REQPARAM        002                     
#define ARGPARAM        004                     
#define INPARAM         010                     
#define OUTPARAM        020                     
#define HIDPARAM        040                     

#endif  // ! _getparam_h

