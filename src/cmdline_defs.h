/*==============================================================================
 HEADER: cmdline_defs.h		[mgpt]
*/

#ifndef _cmdline_defs_h
#define _cmdline_defs_h

#define HEAD1	"NagBody"
#define HEAD2	"mgpt Code for the Modified Gravity-Perturbation Theory."
#define HEAD3	"..."

string defv[] = {  ";"HEAD1": " HEAD2 "\n\t " HEAD3,
    "paramfile=",                   ";Parameter input file. Overwritten by what follows",
//
// Power spectrum table:
    "fnamePS=psLCDM.in",		    ";Filename with power spectrum table",
    "kmin=1e-4",                    ";kmin to analyse from the power spectrum table",
    "kmax=100",                      ";kmax to analyse from the power spectrum table",
    "Nk=300",                       ";Total number of k to analyse from the power spectrum table",
//
// Modified gravity model parameters:
    "mgmodel=HS",                   ";Modified gravity model to study, default f(R) Hu-Sawicki", ":mgm",
    "nHS=1",                        ";Hu-Sawicki index",
    "fR0=1.0e-5",                   ";Hu-Sawicki f_R0",
    "screening=1.0",                ";Hu-Sawicki screening", ":sc",
//
    "om=0.281",                     ";Omega matter value (z=0)",
    "h=0.697",                      ";Hubble parameter value (z=0)",
//
// Differential equations evolution parameters:
    "etaini=-4.0",                  ";Initial eta value :: Log[1/(1 + zini)]",
    "deta=2/5",                     ";deta integration step",
    "detamin=0.",                   ";Min eta integration step size",
    "eps=1.0e-4",                   ";Error parameter",
    "zout=0.0",                     ";redshift value to stop integration :: eta = exp[-eta] - 1",
    "maxnsteps=10000",              ";Maximum number of integration steps", ":maxn",
    "integration_method=bsstep",	";Integration method to use", ":im",
//
// Integration parameters:
    "ngausslegpoints=10",           ";Maximum number of steps", ":nglpts",
//
// Post processing parameters:
    "postprocessing=false",			";Post processing options", ":pp",
    "options=",                     ";Various control options", ":opt",
//
    "Version=1.0.0",                ";Mario A. Rodríguez-Meza/Alejandro Aviles 2018",
    NULL,
};

#endif // ! _cmdline_defs_h
