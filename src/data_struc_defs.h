/*==============================================================================
 HEADER: data_struc_defs.h		[mglpt]
 Written by: Mario A. Rodriguez-Meza
 Starting date: January 2018
 Purpose: Definitions of global variables and parameters
 Language: C
 Use: '#include "global_defs.h"
 Use in routines and functions:
 External headers: stdinc.h, data_struc_defs.h
 Comments and notes:
 Info: M.A. Rodriguez-Meza
 Depto. de Fisica, ININ
 Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
 e-mail: marioalberto.rodriguez@inin.gob.mx
 http://www.inin.gob.mx/
 
 Major revisions:
 Copyright: (c) 2005-2018 Mar.  All Rights Reserved
 ================================================================================
 Legal matters:
 The author does not warrant that the program and routines it contains
 listed below are free from error or suitable for particular applications,
 and he disclaims all liability from any consequences arising from their	use.
 ==============================================================================*/

 
#ifndef _data_struc_defs_h
#define _data_struc_defs_h


#if !defined(global)					// global def question must be here
#  define global extern
#endif


#define IPName(param,paramtext)										\
  {strcpy(tag[nt],paramtext);										\
  addr[nt]=&(param);												\
  id[nt++]=INT;}

#define RPName(param,paramtext)										\
  {strcpy(tag[nt],paramtext);										\
  addr[nt]=&param;													\
  id[nt++]=DOUBLE;}

#define BPName(param,paramtext)										\
  {strcpy(tag[nt],paramtext);										\
  addr[nt]=&param;													\
  id[nt++]=BOOLEAN;}

#define SPName(param,paramtext,n)									\
  {strcpy(tag[nt],paramtext);										\
  param=(string) malloc(n);											\
  addr[nt]=param;													\
  id[nt++]=STRING;}

#endif // ! _data_struc_defs_h

