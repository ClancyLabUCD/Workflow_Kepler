/*
 % This function calls the ode files for EC coupling, CaM reactions, CaMKII
 % phosphorylation module, and PKA phosphorylation module
 %
 % Re-implemented by Mao-Tsuen Jeng,
 % and modified by Pei-Chi Yang

 */

#include "morotti_et_al_mouse_BARsignalling_odefile.h"
#include "morotti_et_al_mouse_CaMKII_signaling_ODEfile.h"
#include "morotti_et_al_mouse_camODEfile.h"
#include "morotti_et_al_mouse_eccODEfile.h"

void morotti_et_al_mouse_masterODEfile( double dt, double tt, double *y, double *p, double *dydt , pars_rec * pars1, double * allDrugs );

void morotti_et_al_mouse_masterODEfile( double dt, double tt, double *y, double *p, double *dydt , pars_rec * pars1, double * allDrugs ) {
	
	
	
	//// Collect params and ICs for each module
	// Allocate ICs for each moduel
	// y(60) -> y(65) are state transitions for mode 1 junctional LCCs
	// y(66) -> y(71) are state transitoins for mode 2 junctional LCCs
	// y(72) -> y(77) are state transitions for mode 1 sarcolemmal LCCs
	// y(78) -> y(83) are state transitions for mode 2 sarcolemmal LCCs
	
	double * y_ecc = &y[0]; // y_ecc = y(1:131);  // Ca_j is y(36), Ca_sl is y(37), Ca_cytosol is y(38)
	double * y_camDyad = &y[130]; // y_camDyad = y(131+1:117+15);
	double * y_camSL = &y[130+15]; // y_camSL = y(131+15+1:117+15+15);
	double * y_camCyt = &y[130+15+15]; // y_camCyt = y(131+15+15+1:131+15+15+15);
	double * y_CaMKII = &y[130+45]; // y_CaMKII = y(131+45+1:117+45+6);      // 6 state vars in CaMKII module
	double * y_BAR = &y[130+45+6]; // y_BAR = y(131+45+6+1:131+45+6+36);    // 36 state vars in BAR module
    
	double * dydt_ecc = &dydt[0];
	double * dydt_camDyad = &dydt[130];
	double * dydt_camSL = &dydt[130+15];
	double * dydt_camCyt = &dydt[130+15+15];
	double * dydt_CaMKIIDyad = &dydt[130+45];
	double * dydt_BAR = &dydt[130+45+6];
	
	// break up parameters from master function into modules
	// paramsCell=mat2cell(p,ones(size(p,1),1),ones(size(p,2),1));
	double cycleLength = p[0];
	double recoveryTime = p[1];
    double variablePar = p[2];
    double CaMtotDyad = p[3];
    double BtotDyad = p[4];
    double CaMKIItotDyad = p[5];
    double CaNtotDyad = p[6];
    double PP1totDyad = p[7];
    double CaMtotSL = p[8];
    double BtotSL = p[9];
    double CaMKIItotSL = p[10];
    double CaNtotSL = p[11];
    double PP1totSL = p[12];
    double CaMtotCyt = p[13];
    double BtotCyt = p[14];
    double CaMKIItotCyt = p[15];
    double CaNtotCyt = p[16];
    double PP1totCyt = p[17];
    double LCCtotDyad = p[18];
    double RyRtot = p[19];
    double PP1_dyad = p[20];
    double PP2A_dyad = p[21];
    double OA = p[22];
    double PLBtot = p[23];
    double LCCtotSL = p[24];
    double PP1_SL = p[25];
    double Ligtot = p[26];
    double LCCtotBA = p[27];
    double RyRtotBA = p[28];
    double PLBtotBA = p[29];
    double TnItotBA = p[30];
    double IKstotBA = p[31];
    double ICFTRtotBA = p[32];
    double PP1_PLBtot = p[33];
    double IKurtotBA = p[34];
    double PLMtotBA = p[35];
    double CKIIOE = p[36];
    
    
	double K = 135; // [mM]
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
	double CaMKIIact_Dyad = CaMKIItotDyad * ( y[130+8 -1] + y[130+9 -1] + y[130+10 -1] + y[130+11 -1] ); // Multiply total by fraction
	double CaMKIIact_SL = CaMKIItotSL * ( y[130+15+8 -1] + y[130+15+9 -1] + y[130+15+10 -1] + y[130+15+11 -1] );
	double PP1_PLB_avail = 1 -  y[130+45+6+24 -1] / PP1_PLBtot + 0.081698;  // Active PP1 near PLB / total PP1 conc + basal value
	double pCaMKII[] = { CaMKIIact_Dyad, LCCtotDyad, RyRtot, PP1_dyad, PP2A_dyad, OA, PLBtot,
		CaMKIIact_SL, LCCtotSL, PP1_SL,
		PP1_PLB_avail };
	
	double LCC_CKdyadp = y[130+45+2 -1] / LCCtotDyad;   // fractional CaMKII-dependent LCC dyad phosphorylation
	double RyR_CKp = y[130+45+4 -1] / RyRtot;           // fractional CaMKII-dependent RyR phosphorylation
	double PLB_CKp = y[130+45+5 -1] / PLBtot;           // fractional CaMKII-dependent PLB phosphorylation
	
	// BAR (PKA phosphorylation) module
	double pBAR[] = { Ligtot,LCCtotBA,RyRtotBA,PLBtotBA,TnItotBA,IKstotBA,ICFTRtotBA,PP1_PLBtot,IKurtotBA,PLMtotBA };
    
	double LCCa_PKAp = y[130+45+6+28 -1]/LCCtotBA;
    double LCCb_PKAp = y[130+45+6+29 -1]/LCCtotBA;
    double PLB_PKAn = (PLBtotBA - y[130+45+6+26 -1])/PLBtotBA;
    double RyR_PKAp = y[130+45+6+30 -1]/RyRtotBA;
    double TnI_PKAp = y[130+45+6+31 -1]/TnItotBA;
    double IKs_PKAp = y[130+45+6+34 -1]/IKstotBA;
    double ICFTR_PKAp = y[130+45+6+35 -1]/ICFTRtotBA;
    double IKur_PKAp = y[130+45+6+36 -1]/IKurtotBA;
    
    double PLM_PKAp = y[130+45+6+27 -1]/PLMtotBA;
    
    // KO phospho-regulation
    // LCC_CKp = 0;
    // RyR_CKp = 0;
    // PLB_CKp = 0;
    
    // Temporary KO (resting values)
    // RyR_CKp = 57.3715/RyRtot;
    // PLB_CKp = 2.0751e-322/PLBtot;
    // LCC_CKp = 5.2773e-60/LCCtot;
    
    // Temporary Isolation OE (1 Hz) - fixed at normal phos levels - updated 05/10/10
    // RyR_CKp = 74.7073/RyRtot;
    // LCC_CKp = 15.4709/LCCtot;
    // PLB_CKp = .4568/PLBtot;
    
    // Temporary Isolation OE (2 Hz) - fixed at normal phos levels
    // USED FOR ANALYSIS OF CAMKII FEEDBACK
    // LCC_CKp = 28.0510/LCCtot;
    // RyR_CKp = 97.0553/RyRtot;
    // PLB_CKp = 1.1833/PLBtot;
    
    // ISO Test (2 Hz) - Fix phos levels at ss normal ISO values
    // USED FOR ARRHYTHMIA ANALYSIS
    // RyR_CKp = 101.9432/RyRtot;
    // LCC_CKp = 28.3810/LCCtot;
    // PLB_CKp = 1.3982/PLBtot;
    
    // ISO Test (1 Hz) - using KO value
    // RyR_CKp = .15;
    // LCC_CKp = 0;
    
    // ECC module
    double pECC[] = { cycleLength,LCC_CKdyadp,RyR_CKp,PLB_CKp, LCCa_PKAp,LCCb_PKAp,PLB_PKAn,RyR_PKAp,TnI_PKAp,IKs_PKAp,ICFTR_PKAp, CKIIOE,recoveryTime,variablePar,Ligtot,IKur_PKAp,PLM_PKAp };
    
    //// Solve dydt in each module
    morotti_et_al_mouse_eccODEfile( dt, tt, y_ecc, pECC, dydt_ecc , pars1, allDrugs  );
    morotti_et_al_mouse_camODEfile( y_camDyad, pCaMDyad, dydt_camDyad , pars1 );
    morotti_et_al_mouse_camODEfile( y_camSL, pCaMSL, dydt_camSL , pars1 );
    morotti_et_al_mouse_camODEfile( y_camCyt, pCaMCyt, dydt_camCyt , pars1 );
    morotti_et_al_mouse_CaMKII_signaling_ODEfile( y_CaMKII, pCaMKII, dydt_CaMKIIDyad );
    morotti_et_al_mouse_BARsignalling_odefile( y_BAR, pBAR, dydt_BAR );
    
    // incorporate Ca buffering from CaM, convert JCaCyt from uM/msec to mM/msec
    // global JCaCyt JCaSL JCaDyad;
    dydt_ecc[35] = dydt_ecc[35] + 1e-3 * pars1->JCaDyad;
    dydt_ecc[36] = dydt_ecc[36] + 1e-3 * pars1->JCaSL;
    dydt_ecc[37] = dydt_ecc[37] + 1e-3 * pars1->JCaCyt;
    
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
