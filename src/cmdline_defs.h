/*==============================================================================
 HEADER: cmdline_defs.h		[mgpt]
*/

#ifndef _cmdline_defs_h
#define _cmdline_defs_h

#define HEAD1	""
#define HEAD2	"mgpt code for Modified gravity perturbation theory."
#define HEAD3	"..."

string defv[] = {  ";"HEAD1": " HEAD2 "\n\t " HEAD3,
    "paramfile=",                   ";Parameter input file. Overwritten by what follows",
//
// Power spectrum table:
    "fnamePS=psLCDM.in",		    ";Filename with power spectrum table",
    "kmin=1e-4",                    ";kmin to analyse from the power spectrum table",
    "kmax=100",                     ";kmax to analyse from the power spectrum table",
    "Nk=300",                       ";Total number of k to analyse from the power spectrum table",":nk",
//
// Modified gravity model parameters:
    "mgModel=LCDM",                 ";Modified gravity model to study, default f(R) Hu-Sawicki", ":mgm",
    "suffixModel=_LCDM",            ";Suffix model to add to output filenames", ":suffix",
    "nHS=1",                        ";Hu-Sawicky index",
    "fR0=1.0e-5",                   ";Hu-Sawicky f_R0",
    "screening=1.0",                ";Hu-Sawicky screening", ":sc",
// DGP:
    "epsDGP=-1.0",                  ";DGP parameter, use with mgModel=DGP",":epsdgp",
    "rcDGP=1.0",                    ";DGP parameter, use with mgModel=DGP",":rcdgp",
//
    "modelParamfile=",              ";If mgmodel=USER, you may give its parameter file name", ":mpf",
//
// Background cosmology:
    "Om=0.281",                     ";Omega matter value (z=0)",":om",
    "OL= 1 - Om",                   ";Omega Lambda value (z=0)",":ol",
    "h=0.697",                      ";Hubble parameter value (z=0)",
//
// Differential equations evolution parameters:
    "etaini=-4.0",                  ";Initial eta value :: Log[1/(1 + zini)]",
    "deta=2/5",                     ";deta integration step",
    "detamin=0.",                   ";Min eta integration step size",
    "eps=1.0e-4",                   ";Differential equations solver tolerance error parameter",
    "zout=0.0",                     ";Output redshift value",
    "maxnsteps=10000",              ";Maximum number of integration steps", ":maxn",
    "solverMethod=bsstep",	        ";Integration method to use", ":solver",
//
// Quadrature parameters:
    "quadratureMethod=trapezoid",   ";Quadrature method to use", ":quadm",
    "nquadSteps=100",               ";Number of k´s from the power spectrum table to integrate (trapezoid)",":nquad",
    "ngausslegpoints=10",           ";Number of Gauss-Legendre of integration points", ":nglpts",
    "epsquad=1.0e-5",               ";Quadrature tolerance error parameter (open Romberg method: romo)",
//
// Post processing parameters:
    "postprocessing=false",			";Post processing options", ":pp",
    "options=",                     ";Various control options", ":opt",
//
    "Version=1.0.0",                ";Mario A. Rodríguez-Meza/Alejandro Aviles 2018",
    NULL,
};

#endif // ! _cmdline_defs_h
