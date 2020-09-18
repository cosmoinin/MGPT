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
    "fnamePS=psLCDM.in",		    ";Input filename with the linear power spectrum table (k,P(k))",
    "kmin=1e-3",                    ";kmin to compute the power spectrum",
    "kmax=100",                     ";kmax to compute the power spectrum",
    "Nk=200",                       ";Total number of k in the power spectrum",":nk",
//
// CLPT correlation functions table:
    "rmin=50",                      ";rmin of the range for CLPT correlation functions",
    "rmax=130",                     ";rmax of the range for CLPT correlation functions",
    "Nr=100",                       ";Total number of r´s in the CLPT correlation function",":nr",
// Modified gravity model parameters:
    "mgModel=LCDM",                 ";Modified gravity model to study (HS or DGP), default is LCDM", ":mgm",
    "suffixModel=",                 ";Suffix model to add to output filenames", ":suffix",
    "fR0=1.0e-5",                   ";Hu-Sawicky f_R0",
    "screening=1.0",                ";set to =0 if you want no screenings", ":sc",
// DGP:
    "epsDGP=-1.0",                  ";epsilon DGP parameter, =-1 normal branch, =1 self-accelerating branch",":epsdgp",
    "rcDGP=1.0",                    ";crossover scale parameter in DGP, (in units of 1/H0)",":rcdgp",
//
    "modelParamfile=",              ";If mgmodel=USER, to use the model in models_user.h", ":mpf",
//
// Background cosmology:
    "Om=0.281",                     ";Omega matter value (z=0)",":om",
    "OL= 1 - Om",                   ";Omega Lambda value (z=0). Only works for DGP!",":ol",
    "h=0.697",                      ";Hubble parameter",
//
// Differential equations evolution parameters:
    "etaini=-4.0",                  ";Initial conformal time value :: Log[1/(1 + zini)]",
    "deta=2/5",                     ";Conformal time integration step",
    "detamin=0.",                   ";Min conformal time integration step size",
    "epsSolver=1.0e-4",             ";Differential equations solver tolerance error parameter",":epssolver",
    "zout=0.0",                     ";Output redshift value",
    "maxnsteps=10000",              ";Maximum number of integration steps", ":maxn",
    "solverMethod=rkqs",	        ";Integration method to use", ":solver",
//
// Quadrature parameters:
    "quadratureMethod=trapezoid3",   ";Quadrature method to use", ":quadm",
    "nquadSteps=200",               ";Number of k´s from the power spectrum table to integrate (trapezoid)",":nquad",
    "ngausslegpoints=16",           ";Number of Gauss-Legendre of integration points", ":nglpts",
    "epsquad=1.0e-6",               ";Quadrature tolerance error parameter (open Romberg method: romberg)",
//
// Post processing parameters:
    "postprocessing=false",			";Post processing options", ":pp",
    "options=",                     ";Various control options", ":opt",
//
    "Version=1.0.0",                ";Mario A. Rodríguez-Meza/Alejandro Aviles 2018",
    NULL,
};

#endif // ! _cmdline_defs_h
