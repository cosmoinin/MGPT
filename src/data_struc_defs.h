/*==============================================================================
 HEADER: data_struc_defs.h		[mglpt]
 ==============================================================================*/

 
#ifndef _data_struc_defs_h
#define _data_struc_defs_h


#if !defined(global)	
#  define global extern
#endif


#define IPName(param,paramtext)										\
  {strcpy(tag[nt],paramtext);										  \
  addr[nt]=&(param);												         \
  id[nt++]=INT;}

#define RPName(param,paramtext)										\
  {strcpy(tag[nt],paramtext);										  \
  addr[nt]=&param;													          \
  id[nt++]=DOUBLE;}

#define BPName(param,paramtext)										\
  {strcpy(tag[nt],paramtext);										  \
  addr[nt]=&param;													          \
  id[nt++]=BOOLEAN;}

#define SPName(param,paramtext,n)								\
  {strcpy(tag[nt],paramtext);										  \
  param=(string) malloc(n);											   \
  addr[nt]=param;													           \
  id[nt++]=STRING;}

#endif // ! _data_struc_defs_h

