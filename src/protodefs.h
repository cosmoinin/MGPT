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


real psLCDMf(real k);       // Interpolation of the power spectrum
global  real psInterpolation(real k, pointPSTableptr PSLtab, int nPSL);
global  real psInterpolation_nr(real k, double kPS[], double pPS[], int nPS);

// MGLPT DIFFEQS
global void integration(double ystart[], int nvar, double x1, double x2, double eps, double h1,
                        double hmin, int *nok, int *nbad, int maxnsteps,
                        void (*derivsin)(double, double [], double []));
global void derivsFirstOrder(double x,double y[],double dydx[]);
global real DpFunction(real k);
global void derivsSecondOrder(double x,double y[],double dydx[]);
global global_D2_ptr DsSecondOrder_func(real kf, real k1, real k2);
global global_D2v2_ptr DsSecondOrder_func_ver2(real kf, real k1, real k2);
global global_D3v2_ptr DsThirdOrder_func_ver2(real x, real k, real p);

// MGLPT QUADS
global_QRs QsR1R2_functions_driver(real eta, real ki);
global real GaussLegendreQ1_func_ver3(real y);
global real GaussLegendreQ2_func_ver3(real y);
global real GaussLegendreQ3_func_ver3(real y);
global real GaussLegendreQ8_func_ver3(real y);
global real GaussLegendreQ9_func_ver3(real y);
global real GaussLegendreQ13_func_ver3(real y);
global real GaussLegendreQI_func_ver3(real y);
global real GaussLegendreQ5_func_ver3(real y);
global real GaussLegendreQ7_func_ver3(real y);
global real GaussLegendreQ11_func_ver3(real y);
global real GaussLegendreQ12_func_ver3(real y);
global real GaussLegendreRI_func_ver3(real y);
global real GaussLegendreR1p2_func_ver3(real y);
global real GaussLegendreR1_func_ver3(real y);
global real GaussLegendreR2_func_ver3(real y);
global real funcR1int_ver4(real y);


// functions (mglpt_fns)
real OmM(real eta);
real H(real eta);
real f1(real eta);
real f2(real eta);
global real A0(real eta);
//

#endif // ! _protodefs_h
