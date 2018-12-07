/*
 *  soltis_biophysJ2010_CaMKII_signaling_ODEfile.h
 *  
 *
 *  Created by Mao-Tsuen Jeng on 12/23/11.
 *  Copyright 2011 __UCDavis__. All rights reserved.
 *
 */

/*
 function dydt = soltis_biophysJ2010_CaMKII_signaling_ODEfile(t,y,p)
 % This function computes the CaMKII-dependent phosphorylation profiles for
 % LCCs (dyadic and subsarcolemmal), RyRs, and PLB. 
 
 % Created by Anthony Soltis (ars7h@virginia.edu) on 01/23/09
 % Updated 12/08/09 as final signaling module
 % Final updates completed 07/21/10
 
 % NOTE: This module includes PKA-dependent phosphorylation of LCCs and
 % RyRs, but these outputs do not feed into any other moduel (PKA
 % phosphorylation is currently computed in the BAR module). The ODEs for
 % these interactions can be set to zero without any consequence to the
 % overall model. 
 % 
 % References:
 % -----------
 */
void soltis_biophysJ2010_CaMKII_signaling_ODEfile( double * y, double * p, double * dydt);
void soltis_biophysJ2010_CaMKII_signaling_ODEfile( double * y, double * p, double * dydt) {
	
	
	//// RATE CONSTANTS and KM VALUES
	// L-Type Ca Channel (LCC) parameters
	double k_ckLCC = 0.4;                  // [s^-1]
	double k_pp1LCC = 0.1103;              // [s^-1] 
	double k_pkaLCC = 13.5;                // [s^-1] 
	double k_pp2aLCC = 10.1;               // [s^-1] 
	
	double KmCK_LCC = 12.;                  // [uM] 
	double KmPKA_LCC = 21.;                 // [uM] 
	double KmPP2A_LCC = 47.;                // [uM] 
	double KmPP1_LCC = 9.;                  // [uM] 
	
	// Ryanodine Receptor (RyR) parameters
	double k_ckRyR = 0.4;                  // [s^-1] 
	double k_pkaRyR = 1.35;                // [s^-1] 
	double k_pp1RyR = 1.07;                // [s^-1] 
	double k_pp2aRyR = 0.481;              // [s^-1] 
	
	// Basal RyR phosphorylation (numbers based on param estimation)
	double kb_2809 = 0.51;                 // [uM/s] - PKA site
	double kb_2815 = 0.35;                 // [uM/s] - CaMKII site
	
	double KmCK_RyR = 12.;                  // [uM] 
	double KmPKA_RyR = 21.;                 // [uM] 
	double KmPP1_RyR = 9.;                  // [uM] 
	double KmPP2A_RyR = 47.;                // [uM] 
	
	// Phospholamban (PLB) parameters
	double k_ckPLB = 8e-3;                 // [s^-1]
	double k_pp1PLB = .0428;               // [s^-1]
	
	double KmCK_PLB = 12.;
	double KmPP1_PLB = 9.;
	
	// Okadaic Acid inhibition params (based on Huke/Bers [2008])
	// Want to treat OA as non-competitive inhibitor of PP1 and PP2A
	double Ki_OA_PP1 = 0.78;        // [uM] - Values from fit
	double Ki_OA_PP2A = 0.037;      // [uM] - Values from fit
	
	//// Assign input params and y0s 
	
	double CaMKIIactDyad = p[0];
	double LCCtotDyad = p[1];
	double RyRtot = p[2];
	double PP1_dyad = p[3];
	double PP2A_dyad = p[4];
	double OA = p[5];
	double PLBtot = p[6];
	double CaMKIIactSL = p[7];
	double LCCtotSL = p[8];
	double PP1_SL = p[9];
	double PP1_PLB_avail = p[10];
	
	//// Description of state variables
    double LCC_PKAp = y[0]; // [LCCp] by PKA (currently unused anywhere else)
	double LCC_CKdyadp = y[1];  // Dyadic [LCCp] by dyadic CaMKII
	double RyR2809p = y[2]; // [RyR-Ser2809p] by PKA (currently unused anywhere else)
	double RyR2815p = y[3]; // [RyR-Ser2815p] by CaMKII
	double PLBT17p = y[4];  // [PLB-Thr17p] by CaMKII
	double LCC_CKslp = y[5];    // Subsarcolemmal [LCCp] by subsarcolemmal CaMKII
	
	
	// Default PKA level
	double PKAc = 95.6 * .54;
	
	//// OA inhibition term (non-competitive) for PP1 and PP2A
	double OA_PP1 = 1. / ( 1 + pow( (OA/Ki_OA_PP1),3 ) );
	double OA_PP2A = 1. / ( 1 + pow( (OA/Ki_OA_PP2A),3 ) );
	
	//// ODE EQUATIONS
	//// LCC states (note: PP2A is acting on PKA site and PP1 on CKII site)
	// CaMKII phosphorylation of Dyadic LCCs
	double LCC_CKdyadn = LCCtotDyad - LCC_CKdyadp;
	double LCCDyad_PHOS = ( k_ckLCC * CaMKIIactDyad * LCC_CKdyadn ) / ( KmCK_LCC + LCC_CKdyadn );
	double LCCDyad_DEPHOS = ( k_pp1LCC * PP1_dyad * LCC_CKdyadp ) / ( KmPP1_LCC + LCC_CKdyadp ) * OA_PP1;
	double dLCC_CKdyadp = ( LCCDyad_PHOS - LCCDyad_DEPHOS );
	
	// CaMKII phosphorylation of Sub-sarcolemmal LCCs
	double LCC_CKsln = LCCtotSL - LCC_CKslp;
	double LCCSL_PHOS = ( k_ckLCC * CaMKIIactSL * LCC_CKsln ) / ( KmCK_LCC + LCC_CKsln ); 
	double LCCSL_DEPHOS = ( k_pp1LCC * PP1_SL * LCC_CKslp ) / ( KmPP1_LCC + LCC_CKslp ) * OA_PP1;
	double dLCC_CKslp = ( LCCSL_PHOS - LCCSL_DEPHOS ); 
	
	// PKA phosphorylation (currently unused elsewhere)
	double LCC_PKAn = LCCtotDyad - LCC_PKAp;
	double dLCC_PKAp = ( ( k_pkaLCC * PKAc * LCC_PKAn ) 
						/ ( KmPKA_LCC + LCC_PKAn ) 
						- ( k_pp2aLCC * PP2A_dyad * LCC_PKAp ) 
						/ ( KmPP2A_LCC + LCC_PKAp) 
						* OA_PP2A );
	//// RyR states
	double RyR2815n = RyRtot - RyR2815p;
	double RyR_BASAL = kb_2815 * RyR2815n;
	double RyR_PHOS = ( k_ckRyR * CaMKIIactDyad * RyR2815n ) / ( KmCK_RyR + RyR2815n );
	double RyR_PP1_DEPHOS = ( k_pp1RyR * PP1_dyad * RyR2815p ) / ( KmPP1_RyR + RyR2815p ) * OA_PP1;
	double RyR_PP2A_DEPHOS = ( k_pp2aRyR * PP2A_dyad * RyR2815p ) / ( KmPP2A_RyR + RyR2815p ) * OA_PP2A;
	double dRyR2815p = ( RyR_BASAL + RyR_PHOS - RyR_PP1_DEPHOS - RyR_PP2A_DEPHOS);
	
	// PKA phosphorylation of Ser 2809 on RyR (currently unused elsewhere)
	double RyR2809n = RyRtot - RyR2809p;
	double dRyR2809p = (kb_2809*RyR2809n + (k_pkaRyR*PKAc*RyR2809n)/(KmPKA_RyR+RyR2809n) - 
						(k_pp1RyR*PP1_dyad*RyR2809p)/(KmPP1_RyR+RyR2809p)*OA_PP1);        
	//// PLB states
	double PP1_PLB = PP1_dyad * PP1_PLB_avail;    // Inhibitor-1 regulation of PP1_dyad included here
	double PLBT17n = PLBtot - PLBT17p;
	double PLB_PHOS = ( k_ckPLB * PLBT17n * CaMKIIactDyad ) / ( KmCK_PLB + PLBT17n );
	double PLB_DEPHOS = ( k_pp1PLB * PP1_PLB * PLBT17p ) / ( KmPP1_PLB + PLBT17p ) * OA_PP1;
	double dPLBT17p = ( PLB_PHOS - PLB_DEPHOS ); 
	
	//// Collect ODEs and convert to uM/ms
	// dydt = [dLCC_PKAp; 
	//        dLCC_CKdyadp;
	//        dRyR2809p; 
	//        dRyR2815p;
	//        dPLBT17p; 
	//        dLCC_CKslp].*10^-3;  // Convert to uM/ms
	dydt[0] = dLCC_PKAp * 0.001; 
	dydt[1] = dLCC_CKdyadp * 0.001;
	dydt[2] = dRyR2809p * 0.001; 
	dydt[3] = dRyR2815p * 0.001;
	dydt[4] = dPLBT17p * 0.001; 
	dydt[5] = dLCC_CKslp * 0.001;  // Convert to uM/ms
	
}
