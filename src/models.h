/*==============================================================================
	HEADER: models.h			[model]
	Written by: M.A. Rodriguez-Meza
	Starting date: January 2018
	Purpose: proto definitios of some model routines
	Language: C
	Use: '#include "...."
	Use in routines and functions: datanaly_md (main)
	External headers: None
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions:
	Copyright: (c) 2005-2018 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#ifndef _models_h
#define _models_h

global void set_model(void);


// ==========================================
// Begin: HS Model

global real mass(real eta);
global real mu(real eta, real k);

global real JFL(real eta, real x, real k, real p);

global real KFL(real eta, real k, real k1, real k2);
global real sourceA(real eta, real kf, real k1, real k2);
global real sourceb(real eta, real kf, real k1, real k2);


global real KFL2(real eta, real x, real k, real p);

global real PiF(real eta, real k);
global real M1(real eta);
global real M2(real eta);
global real M3(real eta);


global real kpp(real x, real k, real p);

global real S2a(real eta, real x, real k, real p);
global real S2b(real eta, real x, real k, real p);
global real S2FL(real eta, real x, real k, real p);
global real S2dI(real eta, real x, real k, real p);

global real SD2(real eta, real x, real k, real p);

global real S3I(real eta, real x, real k, real p, real Dpk, real Dpp,
                real D2f, real D2mf);
global real S3II(real eta, real x, real k, real p, real Dpk, real Dpp,
                 real D2f, real D2mf);
global real S3FL(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f, real D2mf);
global real S3dI(real eta, real x, real k, real p, real Dpk, real Dpp,
                 real D2f, real D2mf);

// End: HS Model
// ==========================================




// ==========================================
// Begin: HS Model

// End: HS Model
// ==========================================


// ==========================================
// Begin: HS Model

// End: HS Model
// ==========================================

#endif // !_models_h
