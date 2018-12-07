/*
 *  soltis_biophysJ2010_masterODEfile.h
 *  
 *
 *  Created by Mao-Tsuen Jeng on 12/22/11.
 *  Copyright 2011 __UCDavis__. All rights reserved.
 *
 */

// ************************************************************
//    Modified from the following matlab function:
//       DAEs are solved using iterative method.
//       ODEs are integrated using RK2 with adaptive step size.
// ************************************************************ 

#include "soltis_biophysJ2010_BARsignalling_odefile.h"
#include "soltis_biophysJ2010_CaMKII_signaling_ODEfile.h"
#include "soltis_biophysJ2010_camODEfile.h"
#include "soltis_biophysJ2010_eccODEfile.h"


/*
 function dydt = soltis_biophysJ2010_masterODEfile(t,y,p)
 % This function calls the ode files for EC coupling, CaM reactions, CaMKII
 % phosphorylation module, and PKA phosphorylation module
 
 % Dependencies: soltis_biophysJ2010_camODEfile.m,
 %               soltis_biophysJ2010_eccODEfile.m
 %
 % Author: Jeff Saucerman <jsaucerman@virginia.edu>
 % Copyright 2008, University of Virginia, All Rights Reserved
 %
 % References:
 % -----------
 % JJ Saucerman and BM Bers, Calmodulin mediates differential
 % sensitivity of CaMKII and calcineurin to local Ca2+ in cardiac myocytes. 
 % Biophys J. 2008 Aug 8. [Epub ahead of print] 
 % Please cite the above paper when using this model.
 */
void soltis_biophysJ2010_masterODEfile( double dt, double tt, double *y, double *p, double *dydt, pars_rec * pars1, double * allDrugs );

void soltis_biophysJ2010_masterODEfile( double dt, double tt, double *y, double *p, double *dydt, pars_rec * pars1, double * allDrugs ) {
	
	
	
	//// Collect params and ICs for each module
	// Allocate ICs for each moduel
	// y[59] -> y[64] are state transitions for mode 1 junctional LCCs
	// y[65] -> y[70] are state transitoins for mode 2 junctional LCCs
	// y[71] -> y[76] are state transitions for mode 1 sarcolemmal LCCs
	// y[77] -> y[82] are state transitions for mode 2 sarcolemmal LCCs
	
	double * y_ecc = &y[0]; // y_ecc = y[0:124];  // Ca_j is y[35], Ca_sl is y[36], Ca_cytosol is y[37]
	double * y_camDyad = &y[125]; // y_camDyad
	double * y_camSL = &y[125+15]; // y_camSL
	double * y_camCyt = &y[125+15+15] ; // y_camCyt
	double * y_CaMKII = &y[125+45]; // y_CaMKII     // 6 state vars in CaMKII module
	double * y_BAR = &y[125+45+6]; // y_BAR     // 30 state vars in BAR module
	
	double * dydt_ecc = &dydt[0]; 
	double * dydt_camDyad = &dydt[125]; 
	double * dydt_camSL = &dydt[125+15]; 
	double * dydt_camCyt = &dydt[125+15+15] ; 
	double * dydt_CaMKIIDyad = &dydt[125+45]; 
	double * dydt_BAR = &dydt[125+45+6]; 
	
	// break up parameters from master function into modules
	
	double cycleLength = p[0];
	double recoveryTime = p[1];
	double CaMtotDyad = p[2];
	double BtotDyad = p[3];
	double CaMKIItotDyad = p[4];
	double CaNtotDyad = p[5];
	double PP1totDyad = p[6];
	double CaMtotSL = p[7];
	double BtotSL = p[8];
	double CaMKIItotSL = p[9];
	double CaNtotSL = p[10];
	double PP1totSL = p[11];
	double CaMtotCyt = p[12];
	double BtotCyt = p[13];
	double CaMKIItotCyt = p[14];
	double CaNtotCyt = p[15];
	double PP1totCyt = p[16];
	double LCCtotDyad = p[17];
	double RyRtot = p[18];
	double PP1_dyad = p[19];
	double PP2A_dyad = p[20];
	double OA = p[21];
	double PLBtot = p[22];
	double LCCtotSL = p[23];
	double PP1_SL = p[24];
	double Ligtot = p[25];
	double LCCtotBA = p[26];
	double RyRtotBA = p[27];
	double PLBtotBA = p[28];
	double TnItotBA = p[29];
	double IKstotBA = p[30];
	double ICFTRtotBA = p[31];
	double PP1_PLBtot = p[32];
	double CKIIOE = p[33];
	
	double PLMtotBA = p[34];
	double factor_Ito = p[35];
	
	double K = 120; // [mM]
	double Mg = 1;  // [mM]
	
	//// Distribute parameters by module
	
	// CaM module
	double CaDyad = y[35]*1e3; // from ECC model, *** Converting from [mM] to [uM] ***
	double compart_dyad = 2;
	// ** NOTE: Btotdyad being sent to the dyad camODEfile is set to zero, but is used below for transfer between SL and dyad
	double pCaMDyad[] = { K, Mg, CaMtotDyad, 0, CaMKIItotDyad, CaNtotDyad, PP1totDyad, CaDyad, cycleLength, compart_dyad };
	double CaSL = y[36]*1e3; // from ECC model, *** Converting from [mM] to [uM] ***
	double compartSL = 1;
	double pCaMSL[] = { K, Mg, CaMtotSL, BtotSL, CaMKIItotSL, CaNtotSL, PP1totSL, CaSL, cycleLength, compartSL };
	double CaCyt = y[37]*1e3; // from ECC model, *** Converting from [mM] to [uM] ***
	double compartCyt = 0;
	double pCaMCyt[] = { K, Mg, CaMtotCyt, BtotCyt, CaMKIItotCyt, CaNtotCyt, PP1totCyt, CaCyt, cycleLength, compartCyt };
	
	// CaMKII phosphorylation module 
	double CaMKIIact_Dyad = CaMKIItotDyad * ( y[83+8+34-1 +8 ] + y[83+9+34-1 +8 ] + y[83+10+34-1 +8 ] + y[83+11+34-1 +8 ] ); // Multiply total by fraction
	double CaMKIIact_SL = CaMKIItotSL * ( y[83+15+8+34-1 +8 ] + y[83+15+9+34-1 +8 ] + y[83+15+10+34-1 +8 ] + y[83+15+11+34-1 +8 ] );
	double PP1_PLB_avail = y[83+45+6+22+34-1 +8 ] / PP1_PLBtot + .0091;  // Active PP1 near PLB / total PP1 conc + basal value
	double pCaMKII[] = { CaMKIIact_Dyad, LCCtotDyad, RyRtot, PP1_dyad, PP2A_dyad, OA, PLBtot, 
		CaMKIIact_SL, LCCtotSL, PP1_SL, 
		PP1_PLB_avail };
	
	double LCC_CKdyadp = y[83+45+2+34-1 +8 ] / LCCtotDyad;   // fractional CaMKII-dependent LCC dyad phosphorylation
	double RyR_CKp = y[83+45+4+34-1 +8 ] / RyRtot;           // fractional CaMKII-dependent RyR phosphorylation
	double PLB_CKp = y[83+45+5+34-1 +8 ] / PLBtot;           // fractional CaMKII-dependent PLB phosphorylation
	
	// BAR (PKA phosphorylation) module
	double pBAR[] = { Ligtot, LCCtotBA, RyRtotBA, PLBtotBA, TnItotBA, IKstotBA, ICFTRtotBA, PP1_PLBtot };
	double LCCa_PKAp = y[83+45+6+23+34-1 +8 ] / LCCtotBA;
	double LCCb_PKAp = y[83+45+6+24+34-1 +8 ] / LCCtotBA;
	double PLB_PKAn = ( PLBtotBA - y[83+45+6+19+34-1 +8 ] ) / PLBtotBA;
	double RyR_PKAp = y[83+45+6+25+34-1 +8 ] / RyRtotBA;
	double TnI_PKAp = y[83+45+6+26+34-1 +8 ] / TnItotBA;
	double IKs_PKAp = y[83+45+6+29+34-1 +8 ] / IKstotBA;
	double ICFTR_PKAp = y[83+45+6+30+34-1 +8 ] / ICFTRtotBA;
    double PLM_PKAn = ( PLMtotBA - y[83+45+6+31+34-1 +8]) / PLMtotBA;

		
	// ECC module
	double pECC[] = { cycleLength, LCC_CKdyadp, RyR_CKp, PLB_CKp, 
		LCCa_PKAp, LCCb_PKAp, PLB_PKAn, RyR_PKAp, TnI_PKAp, IKs_PKAp, ICFTR_PKAp, 
		CKIIOE, recoveryTime, PLM_PKAn, factor_Ito};
	
	//// Solve dydt in each module
	soltis_biophysJ2010_eccODEfile( dt, tt, y_ecc, pECC, dydt_ecc, pars1, allDrugs );
	soltis_biophysJ2010_camODEfile( y_camDyad, pCaMDyad, dydt_camDyad, pars1 );
	soltis_biophysJ2010_camODEfile( y_camSL, pCaMSL, dydt_camSL, pars1 );
	soltis_biophysJ2010_camODEfile( y_camCyt, pCaMCyt, dydt_camCyt, pars1 );
	soltis_biophysJ2010_CaMKII_signaling_ODEfile( y_CaMKII, pCaMKII, dydt_CaMKIIDyad );
	soltis_biophysJ2010_BARsignalling_odefile( tt, y_BAR, pBAR, dydt_BAR );
	
	// incorporate Ca buffering from CaM, convert JCaCyt from uM/msec to mM/msec
	// global JCaCyt JCaSL JCaDyad;
	dydt_ecc[35] = dydt_ecc[35] + 1e-3* pars1->JCaDyad;
	dydt_ecc[36] = dydt_ecc[36] + 1e-3* pars1->JCaSL;
	dydt_ecc[37] = dydt_ecc[37] + 1e-3* pars1->JCaCyt; 
	
	// incorporate CaM diffusion between compartments
	double Vmyo = 2.1454e-11;      // [L]
	double Vdyad = 1.7790e-014;    // [L]
	double VSL = 6.6013e-013;      // [L]
	double kDyadSL = 3.6363e-16;	// [L/msec]
	double kSLmyo = 8.587e-15;     // [L/msec]
	double k0Boff = 0.0014;        // [s^-1] 
	double k0Bon = k0Boff / 0.2;     // [uM^-1 s^-1] kon = koff/Kd
	double k2Boff = k0Boff / 100;    // [s^-1] 
	double k2Bon = k0Bon;          // [uM^-1 s^-1]
	double k4Boff = k2Boff;        // [s^-1]
	double k4Bon = k0Bon;          // [uM^-1 s^-1]
	
	//	CaMtotDyad = sum(y_camDyad(1:6)) + CaMKIItotDyad*sum(y_camDyad(7:10)) + sum(y_camDyad(13:15));
	double sum1, sum2, sum3;
	sum1 = 0;
	sum2 = 0; 
	sum3 = 0;
	for ( int it1 = 0; it1 < 6; it1++ ) {
		sum1 += y_camDyad[it1];
	}
	for ( int it1 = 6; it1 < 10; it1++ ) {
		sum2 += y_camDyad[it1];
	}
	sum2 = sum2 * CaMKIItotDyad; 
	for ( int it1 = 12; it1 < 15; it1++ ) {
		sum3 += y_camDyad[it1];
	}
	CaMtotDyad = sum1 + sum2 + sum3;
	
	double Bdyad = BtotDyad - CaMtotDyad; // [uM dyad]
	double J_cam_dyadSL = 1e-3 * ( k0Boff * y_camDyad[0] - k0Bon * Bdyad * y_camSL[0] ); // [uM/msec dyad]
	double J_ca2cam_dyadSL = 1e-3 * ( k2Boff * y_camDyad[1] - k2Bon * Bdyad * y_camSL[1] ); // [uM/msec dyad]
	double J_ca4cam_dyadSL = 1e-3 * ( k2Boff * y_camDyad[2] - k4Bon * Bdyad * y_camSL[2] ); // [uM/msec dyad]
	double J_cam_SLmyo = kSLmyo * ( y_camSL[0] - y_camCyt[0] ); // [umol/msec]
	double J_ca2cam_SLmyo = kSLmyo * ( y_camSL[1] - y_camCyt[1] ); // [umol/msec]
	double J_ca4cam_SLmyo = kSLmyo * ( y_camSL[2] - y_camCyt[2] ); // [umol/msec]
	dydt_camDyad[0] = dydt_camDyad[0] - J_cam_dyadSL;
	dydt_camDyad[1] = dydt_camDyad[1] - J_ca2cam_dyadSL;
	dydt_camDyad[2] = dydt_camDyad[2] - J_ca4cam_dyadSL;
	dydt_camSL[0] = dydt_camSL[0] + J_cam_dyadSL * Vdyad / VSL - J_cam_SLmyo / VSL;
	dydt_camSL[1] = dydt_camSL[1] + J_ca2cam_dyadSL * Vdyad / VSL - J_ca2cam_SLmyo / VSL;
	dydt_camSL[2] = dydt_camSL[2] + J_ca4cam_dyadSL * Vdyad / VSL - J_ca4cam_SLmyo / VSL;
	dydt_camCyt[0] = dydt_camCyt[0] + J_cam_SLmyo / Vmyo;
	dydt_camCyt[1] = dydt_camCyt[1] + J_ca2cam_SLmyo / Vmyo;
	dydt_camCyt[2] = dydt_camCyt[2] + J_ca4cam_SLmyo / Vmyo;
	
	// Collect all dydt terms
	// dydt = [dydt_ecc; dydt_camDyad; dydt_camSL; dydt_camCyt; dydt_CaMKIIDyad; dydt_BAR];
	
}
