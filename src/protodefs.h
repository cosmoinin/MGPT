/*==============================================================================
 HEADER: protodefs.h				[mglpt]
 ==============================================================================*/

#ifndef _protodefs_h
#define _protodefs_h

void integration_method_string_to_int(string,int *);

void output(void);

void MainLoop(void);
void StartRun(string, string, string, string);
void StartOutput(void);
void EndRun(void);

// Postprocessing
global void PostProcessing(void);
global void biasterms_processing(void);

// MGLPT DIFFEQS and QUAD
global void integration(double ystart[], int nvar, double x1, double x2, double eps, double h1,
                        double hmin, int *nok, int *nbad, int maxnsteps,
                        void (*derivsin)(double, double [], double []));
global real DpFunction(real k);
global global_D2v2_ptr DsSecondOrder_func(real kf, real k1, real k2);
global global_D3v2_ptr DsThirdOrder_func(real x, real k, real p);

// MGLPT QUADS
global_QRs QsRs_functions_driver(real eta, real ki);


// functions (mgpt_fns)
real OmM(real eta);
real H(real eta);
real f1(real eta);
real f2(real eta);
global real A0(real eta);
global  real psInterpolation_nr(real k, double kPS[], double pPS[], int nPS);
//

#endif // ! _protodefs_h
