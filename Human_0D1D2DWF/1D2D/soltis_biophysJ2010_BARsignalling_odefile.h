/*
 *  soltis_biophsJ2010_BARsignalling_odefile.h
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

/*
 % function ydot = soltis_biophysJ2010_BARsignalling_odefile(t,y,pin)
 % Re-implemented by Anthony Soltis <ars7h@virginia.edu> for communication
 % with Shannon (2004) EC coupling model. Final version 07/21/10
 
 % saucerman_circres2004.m
 % coupled beta adrenergic signaling/EC for adult rabbit ventricular
 % myocytes, with extensions for yotiao interactions, KCNQ1-G589D mutation
 %
 % Copyright (2004) The Regents of the University of California
 % All Rights Reserved
 % Permission to use, copy, and modify, any part of this software for
 % academic purposes only, including non-profit  education and research,
 % without fee, and without a written agreement is hereby granted, provided
 % that the above copyright notice, this paragraph and the following three
 % paragraphs appear in all copies.
 % The receiving party agrees not to further distribute nor disclose the
 % source code to third parties without written permission and not to use
 % the software for research under commercial sponsorship or for any other
 % any commercial undertaking.  Those desiring to incorporate this software
 % into commercial products of to use it for commercial purposes should
 % contact the Technology Transfer Office, University of California, San
 % Diego, 9500 Gilman Drive, La Jolla, California, 92093-0910, Ph: (619)
 % 534 5815, referring to Case SDC98008.
 % IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
 % FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
 % INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 % THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
 % DAMAGE.
 % THE SOFTWARE PROVIDED HEREUNDER IS ON AN AS IS BASIS, AND THE
 % UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE,
 % SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.  THE UNIVERSITY OF
 % CALIFORNIA MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY
 % KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE
 % IMPLIED WARRANTIES OF MECHANT ABILITY OR FITNESS FOR A PARTICULAR
 % PURPOSE, OR THAT THE USE OF THE MATERIAL WILL NOT INFRINGE ANY PATENT,
 % TRADEMARK OR OTHER RIGHTS.
 %
 % Those using this code should acknowledge Dr. Andrew D. McCulloch and the
 % National Biomedical Computation Resource (NBCR), NIH grant #P41 RR08605.
 %
 % References:
 % -----------
 % Jeffrey J. Saucerman, Sarah N. Healy, Mary E. Belik, Jose L. Puglisi, and
 % Andrew D. McCulloch.  "Proarrhythmic Consequences of a KCNQ1 AKAP-Binding 
 % Domain Mutation: Computational Models of Whole Cells and Heterogeneous
 % Tissue", Circ. Res., Vol 95: 1216-1224 (2004).
 %
 % Jeffrey J. Saucerman and Andrew D. McCulloch
 % "Mechanistic systems models of cell signaling networks: a case study of
 % myocyte adrenergic regulation", Prog. Biophys. Mol. Bio., Vol 84: 261-278 (2004).
 %
 % Jeffrey J. Saucerman, Laurence L. Brunton, Anushka P. Michailova, and Andrew D. McCulloch
 % "Modeling beta-adrenergic control of cardiac myocyte contractility in
 % silico", J. Biol. Chem., Vol 278: 47977-48003 (2003).
 %
 % Last modified: 12/20/2004
 % Implemented by: Jeffrey Saucerman <jsaucer@ucsd.edu>
 %
 % Notes:
 % - This code was used for the single cell simulations in the Circ Res ms.
 % - In Matlab 7, you can enable cell mode to jump between modules.
 % - Please email me for any questions or comments.
 */
// ************************************************************************
// ************************************************************************
#ifndef soltis_biophysJ2010_BARsignalling_odefile_H
#define soltis_biophysJ2010_BARsignalling_odefile_H

void soltis_biophysJ2010_BARsignalling_odefile( double tt, double *y, double *pin, double *ydot );

void soltis_biophysJ2010_BARsignalling_odefile( double tt, double *y, double *pin, double *ydot ) {
	
	//// Assign passed in params
	double Ltot = pin[0];
	double LCCtot = pin[1];
	double RyRtot = pin[2];
	double PLBtot = pin[3];
	double TnItot = pin[4];
	double IKstot = pin[5];
	double ICFTRtot = pin[6];
	double PP1_PLBtot = pin[7];
	double PLMtot= pin[8];
	
	//// Parameters
	//// ----- Signaling model parameters -------
	// b-AR/Gs module
	double p[99];
	p[0] = Ltot;    // Ltotmax   [uM] ** apply agonist concentration here **
	p[1] = 0.028;   // sumb1AR   [uM]
	p[2] = 3.83;    // Gstot     [uM]
	p[3] = 0.285;   // Kl        [uM]
	p[4] = 0.062;   // Kr        [uM]
	p[5] = 33.0;    // Kc        [uM]
	p[6] = 1.1e-3;  // k_barkp   [1/sec]
	p[7] = 2.2e-3;  // k_barkm   [1/sec]
	p[8] = 3.6e-3;  // k_pkap    [1/sec/uM]
	p[9] = 2.2e-3; // k_pkam    [1/sec]
	p[10] = 16.0;   // k_gact    [1/sec]
	p[11] = 0.8;    // k_hyd     [1/sec]
	p[12] = 1.21e3; // k_reassoc [1/sec/uM]
					// cAMP module
	p[13] = 0.047;  // AC_tot    [uM]
	p[14] = 5.0e3;  // ATP       [uM]
	p[15] = 0.036;  // PDE3tot   [uM]   // Changed from .06 to .036
	p[16] = 0.036;  // PDE4tot   [uM]
	p[17] = 0.0;    // IBMXtot   [uM]
	p[18] = 0.0;    // Fsktot    [uM] (10 uM when used)
	p[19] = 0.2;    // k_ac_basal[1/sec]
	p[20] = 8.5;    // k_ac_gsa  [1/sec]
	p[21] = 7.3;    // k_ac_fsk  [1/sec]
	p[22] = 1.03e3; // Km_basal  [uM]
	p[23] = 315.0;  // Km_gsa    [uM]
	p[24] = 860.0;  // Km_fsk    [uM]
	p[25] = 0.4;    // Kgsa      [uM]
	p[26] = 44.0;   // Kfsk      [uM]
	p[27] = 3.5;    // k_pde3    [1/sec]
	p[28] = 0.15;   // Km_pde3   [uM]
	p[29] = 5.0;    // k_pde4    [1/sec]
	p[30] = 1.3;    // Km_pde4   [uM]
	p[31] = 30.0;   // Ki_ibmx   [uM]
					// PKA module
	p[32] = 0.46;   // PKAItot   [uM]
	p[33] = 0.084;  // PKAIItot  [uM]
	p[34] = 0.18;   // PKItot    [uM]
	p[35] = 9.14;   // Ka        [uM]
	p[36] = 1.64;   // Kb        [uM]
	p[37] = 4.375;  // Kd        [uM]
	p[38] = 0.2e-3; // Ki_pki    [uM]
					// PLB module
	p[39] = 10;     // epsilon   [none]
	p[40] = PLBtot; // PLBtot    [uM]
	p[41] = PP1_PLBtot;   // PP1tot    [uM]
	p[42] = 0.3;    // Inhib1tot [uM]
	p[43] = 54;     // k_pka_plb     [1/sec]
	p[44] = 21;     // Km_pka_plb    [uM]
	p[45] = 8.5;    // k_pp1_plb     [1/sec]
	p[46] = 7.0;    // Km_pp1_plb    [uM]
	p[47] = 60;     // k_pka_i1      [1/sec]
	p[48] = 1.0;    // Km_pka_i1     [uM]
	p[49] = 14.0;   // Vmax_pp2a_i1  [uM/sec]
	p[50] = 1.0;    // Km_pp2a_i1    [uM]
	p[51] = 1.0e-3; // Ki_inhib1     [uM]
					// LCC module
	p[52] = LCCtot; // LCCtot        [uM]
	p[53] = 0.025;  // PKAIIlcctot   [uM]
	p[54] = 0.025;  // PP1lcctot     [uM]
	p[55] = 0.025;  // PP2Alcctot    [uM]
	p[56] = 54;     // k_pka_lcc     [1/sec]
	p[57] = 21;     // Km_pka_lcc    [uM]
	p[58] = 8.52;   // k_pp1_lcc     [1/sec]
	p[59] = 3;      // Km_pp1_lcc    [uM]
	p[60] = 10.1;   // k_pp2a_lcc    [1/sec]
	p[61] = 3;      // Km_pp2a_lcc   [uM]
					// RyR module
	p[62] = RyRtot; // RyRtot        [uM]
	p[63] = 0.034;  // PKAIIryrtot   [uM]
	p[64] = 0.034;  // PP1ryr        [uM]
	p[65] = 0.034;  // PP2Aryr       [uM]
	p[66] = 54;     // kcat_pka_ryr  [1/sec]
	p[67] = 21;     // Km_pka_ryr    [uM]
	p[68] = 8.52;   // kcat_pp1_ryr  [1/sec]
	p[69] = 7;      // Km_pp1_ryr    [uM]
	p[70] = 10.1;   // kcat_pp2a_ryr [1/sec]
	p[71] = 4.1;    // Km_pp2a_ryr   [uM]
					// TnI module
	p[72] = TnItot; // TnItot        [uM]
	p[73] = 0.67;   // PP2Atni       [uM]
	p[74] = 54;     // kcat_pka_tni  [1/sec]
	p[75] = 21;     // Km_pka_tni    [uM]
	p[76] = 10.1;   // kcat_pp2a_tni [1/sec]
	p[77] = 4.1;    // Km_pp2a_tni   [uM]
					// Iks module
	p[78] = IKstot; // Iks_tot       [uM]
	p[79] = 0.025;  // Yotiao_tot    [uM]
	p[80] = 0.1e-3; // K_yotiao      [uM] ** apply G589D mutation here **
	p[81] = 0.025;  // PKAII_ikstot  [uM]
	p[82] = 0.025;  // PP1_ikstot    [uM]
	p[83] = 54;     // k_pka_iks     [1/sec]
	p[84] = 21;     // Km_pka_iks    [uM]
	p[85] = 8.52;   // k_pp1_iks     [1/sec]
	p[86] = 7;      // Km_pp1_iks    [uM]
					// Icftr Module - Added 04/30/10 by Anthony Soltis
	p[87] = ICFTRtot;  // CFTR_tot      [uM]
	p[88] = 0.025;  // PKAII_CFTRtot [uM]
	p[89] = 0.025;  // PP1_CFTRtot   [uM]
	p[90] = 54;     // k_pka_CFTR    [1/sec]
	p[91] = 8.5;    // Km_pka_CFTR   [uM]
	p[92] = 8.52;   // k_pp1_CFTR    [1/sec]
	p[93] = 7;      // Km_pp1_CFTR   [uM]
	
	// PLM Module - Added 09/13/12 by Ele Grandi
	p[94] = PLMtot; // PLMtot    [uM]
	p[95] = 54;     // k_pka_plm     [1/sec]
	p[96] = 21;     // Km_pka_plm    [uM]
	p[97] = 8.5;    // k_pp1_plm     [1/sec]
	p[98] = 7.0;    // Km_pp1_plm    [uM]
	
					//// -------- SIGNALING MODEL -----------
	
	// ydot = zeros(size(y));
	//// b-AR module
	double LR = y[0]*y[1]/p[3];
	double LRG = LR*y[2]/p[4];
	double RG = y[1]*y[2]/p[5];
//	double BARKDESENS = p[6]*(LR+LRG);
//	double BARKRESENS = p[7]*y[4];
//	double PKADESENS = p[8]*y[16]*y[3];  
//	double PKARESENS = p[9]*y[5];
//	double GACT = p[10]*(RG+LRG);
//	double HYD = p[11]*y[6];
//	double REASSOC = p[12]*y[7]*y[8];
	
	ydot[0] = p[0]-LR-LRG-y[0];
	ydot[1] = y[3]-LR-LRG-RG-y[1];
	ydot[2] = p[2]-LRG-RG-y[2];
	// Solving DAEs for y[0], y[1], y[2].
	double y0, y1, y2;
	double ip3 = 1. / p[3];
	double ip4 = 1. / p[4];
	double ip5 = 1. / p[5];

	double myeps = 1.e-14;

	int iter = 0;
	y0 = y[0];
	y1 = y[1];
	y2 = y[2];
	
	while ( iter < 1000 && ( fabs( ydot[0] ) > myeps || fabs( ydot[1] ) > myeps || fabs( ydot[2] ) > myeps ) )
	{
		// y[0] = p[0] / ( 1 + y[1] * ip3 + y[1] * y[2] * ip3 * ip4 );
		y0 = p[0] / ( 1. + y[1] * ip3 * ( 1 + y[2] * ip4 ) );
		y1 = y[3] / ( 1. + y0 * ip3 * ( 1 + y[2] * ip4 ) + y[2] * ip5 );
		y2 = p[2] / ( 1. + y0 * y1 * ip3 * ip4 +  y1 * ip5 );

		y[0] = p[0] / ( 1. + y1 * ip3 * ( 1 + y2 * ip4 ) );
		y[1] = y[3] / ( 1. + y[0] * ip3 * ( 1 + y2 * ip4 ) + y2 * ip5 );
		y[2] = p[2] / ( 1. + y[0] * y[1] * ip3 * ip4 +  y[1] * ip5 );
		
		ydot[0] = y[0] - y0;
		ydot[1] = y[1] - y1;
		ydot[2] = y[2] - y2;		
		
		iter++;
	}
//	cout << iter << "\t" << ydot[0] << "\t" << ydot[1] << "\t" << ydot[2] << "\t" << endl;
	// if ( iter > 6 ) {
	// cout << "iter = " << iter << "\n";
	// }
	ydot[0] = 0;
	ydot[1] = 0;
	ydot[2] = 0;	
	// End solving DAEs for y[0], y[1], y[2].
	
	LR = y[0]*y[1]/p[3];
	LRG = LR*y[2]/p[4];
	RG = y[1]*y[2]/p[5];
	double BARKDESENS = p[6]*(LR+LRG);
	double BARKRESENS = p[7]*y[4];
	double PKADESENS = p[8]*y[16]*y[3];  
	double PKARESENS = p[9]*y[5];
	double GACT = p[10]*(RG+LRG);
	double HYD = p[11]*y[6];
	double REASSOC = p[12]*y[7]*y[8];
	
	ydot[3] = (BARKRESENS-BARKDESENS)+(PKARESENS-PKADESENS);
	ydot[4] = BARKDESENS-BARKRESENS;
	ydot[5] = PKADESENS-PKARESENS;
	ydot[6] = GACT-HYD;
	ydot[7] = HYD-REASSOC;
	ydot[8] = GACT-REASSOC;
	// end b-AR module
	
	//// cAMP module
	double Gsa_gtp_AC = y[9]*y[11]/p[25];
	double Fsk_AC = y[10]*y[11]/p[26];
//	double AC_ACT_BASAL = p[19]*y[11]*p[14]/(p[22]+p[14]);	    
//	double AC_ACT_GSA = p[20]*Gsa_gtp_AC*p[14]/(p[23]+p[14]); 
//	double AC_ACT_FSK = p[21]*Fsk_AC*p[14]/(p[24]+p[14]);	   
//	// PDE3_ACT = p[27]*y[12]*y[15]/(p[28]+y[15]);	
//	// PDE4_ACT = p[29]*y[12]*y[15]/(p[30]+y[15]);	
//	double PDE3_ACT = p[27]*p[15]*y[15]/(p[28]*(1+p[17]/p[31])+y[15]);	// new PDE3 term w IBMX
//	double PDE4_ACT = p[29]*p[16]*y[15]/(p[30]*(1+p[17]/p[31])+y[15]);	// new PDE4 term w IBMX
//	double PDE_IBMX = y[12]*y[13]/p[31];
	ydot[9] = y[6]-Gsa_gtp_AC-y[9];
	ydot[10] = p[18]-Fsk_AC-y[10];
	ydot[11] = p[13]-Gsa_gtp_AC-y[11];  // note: assumes Fsk = 0.  Change Gsa_gtp_AC to Fsk_AC for Forskolin.
	
	// Solving DAEs for y[9], y[10], y[11].
	double ip25 = 1. / p[25];
	double ip26 = 1. / p[26];
	double y9, y10, y11;

	y9 = y[9];
	y10 = y[10];
	y11 = y[11];
	
	iter = 0;
	while ( iter < 1000 && ( fabs( ydot[9] ) > myeps || fabs( ydot[10] ) > myeps || fabs( ydot[11] ) > myeps ) )
	{
		y9 = y[6] / ( 1. + y[11] * ip25 );
		y10 = p[18] / ( 1. + y[11] * ip26 );
		y11 = p[13] / ( 1. + y9 * ip25 );		
		
		y[9] = y[6] / ( 1. + y11 * ip25 );
		y[10] = p[18] / ( 1. + y11 * ip26 );
		y[11] = p[13] / ( 1. + y[9] * ip25 );
		
		ydot[9] = y[9] - y9;
		ydot[10] = y[10] - y10;
		ydot[11] = y[11] - y11;		
		
		iter++;
	}
//	cout << iter << "\t" << ydot[9] << "\t" << ydot[10] << "\t" << ydot[11] << "\t" << endl;
	// if ( iter > 4 ) {
	// cout << "iter = " << iter << "\n";
	// }
	ydot[9] = 0;
	ydot[10] = 0;
	ydot[11] = 0;	
	// End solving DAEs for y[9], y[10], y[11].
	
	Gsa_gtp_AC = y[9]*y[11]/p[25];
	Fsk_AC = y[10]*y[11]/p[26];
	double AC_ACT_BASAL = p[19]*y[11]*p[14]/(p[22]+p[14]);	    
	double AC_ACT_GSA = p[20]*Gsa_gtp_AC*p[14]/(p[23]+p[14]); 
	double AC_ACT_FSK = p[21]*Fsk_AC*p[14]/(p[24]+p[14]);	   
	// PDE3_ACT = p[27]*y[12]*y[15]/(p[28]+y[15]);	
	// PDE4_ACT = p[29]*y[12]*y[15]/(p[30]+y[15]);	
	double PDE3_ACT = p[27]*p[15]*y[15]/(p[28]*(1+p[17]/p[31])+y[15]);	// new PDE3 term w IBMX
	double PDE4_ACT = p[29]*p[16]*y[15]/(p[30]*(1+p[17]/p[31])+y[15]);	// new PDE4 term w IBMX
	double PDE_IBMX = y[12]*y[13]/p[31];

	
	// ydot[12] = p[16]-PDE_IBMX-y[12];
	// ydot[13] = p[17]-PDE_IBMX-y[13];
	ydot[12] = 0;
	ydot[13] = 0;
	ydot[14] = AC_ACT_BASAL+AC_ACT_GSA+AC_ACT_FSK-PDE3_ACT-PDE4_ACT;
	// end cAMP module
	
	//// PKA module
	double PKI = p[34]*p[38]/(p[38]+y[16]+y[17]);
	double A2RC_I = (y[16]/p[37])*y[16]*(1+PKI/p[38]);
	double A2R_I = y[16]*(1+PKI/p[38]);
	double A2RC_II = (y[17]/p[37])*y[17]*(1+PKI/p[38]);
	double A2R_II = y[17]*(1+PKI/p[38]);
	double ARC_I = (p[35]/y[15])*A2RC_I;
	double ARC_II = (p[35]/y[15])*A2RC_II;
	ydot[15] = y[14]-(ARC_I+2*A2RC_I+2*A2R_I)-(ARC_II+2*A2RC_II+2*A2R_II)-y[15];
//	double PKAtemp = p[35]*p[36]/p[37]+p[35]*y[15]/p[37]+y[15]^2/p[37];
	double PKAtemp = p[35]*p[36]/p[37]+p[35]*y[15]/p[37]+y[15]*y[15]/p[37];
	ydot[16] = 2*p[32]*y[15]*y[15]-y[16]*(1+PKI/p[38])*(PKAtemp*y[16]+y[15]*y[15]);
	ydot[17] = 2*p[33]*y[15]*y[15]-y[17]*(1+PKI/p[38])*(PKAtemp*y[17]+y[15]*y[15]);
	double y15, y16, y17, y15sq;
	double ip37 = 1./p[37];
	double ip38 = 1./p[38];
	
	// Solving DAEs for y[15], y[16], y[17].
	iter = 0;
	double mk; 
	while ( iter < 1000 && ( fabs( ydot[15] ) > myeps || fabs( ydot[16] ) > myeps || fabs( ydot[17] ) > myeps ) )
	{
		mk = 0.35; // Increase convergent speed when close to fixed point.
		PKI = p[34]*p[38]/(p[38]+y[16]+y[17]);
		A2RC_I = (y[16]/p[37])*y[16]*(1+PKI/p[38]);
		A2R_I = y[16]*(1+PKI/p[38]);
		A2RC_II = (y[17]/p[37])*y[17]*(1+PKI/p[38]);
		A2R_II = y[17]*(1+PKI/p[38]);
		ARC_I = (p[35]/y[15])*A2RC_I;
		ARC_II = (p[35]/y[15])*A2RC_II;
		ydot[15] = y[14]-(ARC_I+2*A2RC_I+2*A2R_I)-(ARC_II+2*A2RC_II+2*A2R_II)-y[15];
//		while ( fabs( mk * ydot[15] / y[15] ) > 0.1 ) {
//			mk = 0.1 * mk;
//		}
		y[15] += mk * ydot[15];
		y15sq = y[15] * y[15];
		PKAtemp = ( p[35] * ( p[36] + y[15] ) + y15sq )* ip37;
		y16 = y[16];
		y[16] = 2*p[32]*y15sq / ( (1+PKI*ip38)*(PKAtemp*y[16]+y15sq) );
		ydot[16] = y[16] - y16;
		y17 = y[17];
		y[17] = 2*p[33]*y15sq / ( (1+PKI*ip38)*(PKAtemp*y[17]+y15sq) );
		ydot[17] = y[17] - y17;
		iter++;
	}
	// if ( iter > 60 ) {
	//  cout << "iter = " << iter << "\n";
	// }
	ydot[15] = 0;
	ydot[16] = 0;
	ydot[17] = 0;	
	// End solving DAEs for y[9], y[10], y[11].
	
	// end PKA module
	
	//// PLB module
	double PLB = p[40]-y[18];
	double PLB_PHOSPH = p[43]*y[16]*PLB/(p[44]+PLB);
	double PLB_DEPHOSPH = p[45]*y[21]*y[18]/(p[46]+y[18]);
	ydot[18] = PLB_PHOSPH-PLB_DEPHOSPH;
	
	double Inhib1 = p[42]-y[19];
	double Inhib1p_PP1 = y[20]*y[21]/p[51];
	double Inhib1_PHOSPH = p[47]*y[16]*Inhib1/(p[48]+Inhib1); 
	double Inhib1_DEPHOSPH = p[49]*y[19]/(p[50]+y[19]);
	ydot[19] = Inhib1_PHOSPH-Inhib1_DEPHOSPH;
	ydot[20] = y[19]-Inhib1p_PP1-y[20];
	ydot[21] = p[41]-Inhib1p_PP1-y[21];
	
	// Solving DAEs for y[20], y[21].
	double ip51 = 1. / p[51];
	double y20, y21;
	
	y20 = y[20];
	y21 = y[21];
		
	iter = 0;
	while ( iter < 1000 && ( fabs( ydot[20] ) > myeps || fabs( ydot[21] ) > myeps ) )
	{
		y20 = y[19] / ( 1 + y[21] * ip51 );
		y21 = p[41] / ( 1 + y20 * ip51 );
		
		y[20] = y[19] / ( 1 + y21 * ip51 );
		y[21] = p[41] / ( 1 + y[20] * ip51 );
		
		ydot[20] = y[20] - y20;
		ydot[21] = y[21] - y21;		
		
		iter++;
	}
	ydot[20] = 0;
	ydot[21] = 0;	
	// End solving DAEs for y[20], y[21].
   // if ( iter > 4 ) {
	// cout << "iter = " << iter << "\n";
	// }
	
	// end PLB module
	
	//// LCC module
	double PKAClcc = (p[53]/p[33])*y[17];
	double LCCa = p[52]-y[22];
	double LCCa_PHOSPH = p[39]*p[56]*PKAClcc*LCCa/(p[57] + p[39]*LCCa);
	double LCCa_DEPHOSPH = p[39]*p[60]*p[55]*y[22]/(p[61]+p[39]*y[22]);
	ydot[22] = LCCa_PHOSPH - LCCa_DEPHOSPH;
	
	double LCCb = p[52]-y[23];
	double LCCb_PHOSPH = p[39]*p[56]*PKAClcc*LCCb/(p[57]+p[39]*LCCb);   
	double LCCb_DEPHOSPH = p[39]*p[58]*p[54]*y[23]/(p[59]+p[39]*y[23]);
	ydot[23] = LCCb_PHOSPH-LCCb_DEPHOSPH;
	// end LCC module
	
	//// RyR module
	double PKACryr = (p[63]/p[33])*y[17];
	double RyR = p[62]-y[24];
	double RyRPHOSPH = p[39]*p[66]*PKACryr*RyR/(p[67]+p[39]*RyR);
	double RyRDEPHOSPH1 = p[39]*p[68]*p[64]*y[24]/(p[69]+p[39]*y[24]);
	double RyRDEPHOSPH2A = p[39]*p[70]*p[65]*y[24]/(p[71]+p[39]*y[24]);
	ydot[24] = RyRPHOSPH-RyRDEPHOSPH1-RyRDEPHOSPH2A;
	// end RyR module
	
	//// TnI module
	double TnI = p[72]-y[25];
	double TnIPHOSPH = p[74]*y[16]*TnI/(p[75]+TnI);
	double TnIDEPHOSPH = p[76]*p[73]*y[25]/(p[77]+y[25]);
	ydot[25] = TnIPHOSPH-TnIDEPHOSPH;
	// end TnI module
	
	//// Iks module
	double IksYot = y[26]*y[27]/p[80];           // [uM]
	ydot[26] = p[78] - IksYot - y[26];    // [uM]
	ydot[27] = p[79] - IksYot - y[27];    // [uM]
	
	// Solving DAEs for y[26], y[27].
	double ip80 = 1. / p[80];
	double y26, y27;
	
	y26 = y[26];
	y27 = y[27];
	
	iter = 0;
	while ( iter < 1000 && ( fabs( ydot[26] ) > myeps || fabs( ydot[27] ) > myeps ) )
	{
		y26 = p[78] / ( 1 + y[27] * ip80 );
		y27 = p[79] / ( 1 + y26 * ip80 );
				
		y[26] = p[78] / ( 1 + y27 * ip80 );
		y[27] = p[79] / ( 1 + y[26] * ip80 );
		
		ydot[26] = y[26] - y26;
		ydot[27] = y[27] - y27;		
		
		iter++;
	}
	ydot[26] = 0;
	ydot[27] = 0;	
	// End solving DAEs for y[26], y[27].	
	// if ( iter > 1 ) {
	// cout << "iter = " << iter << "\n";
	// }
	double PKACiks = (IksYot/p[78])*(p[81]/p[33])*y[17];
	double PP1iks = (IksYot/p[78])*p[82];
	double Iks = p[78]-y[28];
	double IKS_PHOSPH = p[39]*p[83]*PKACiks*Iks/(p[84]+p[39]*Iks);
	double IKS_DEPHOSPH = p[39]*p[85]*PP1iks*y[28]/(p[86]+p[39]*y[28]);
	ydot[28] = IKS_PHOSPH-IKS_DEPHOSPH;
	
	//// CFTR module (included 04/30/10)
	double CFTRn = p[87] - y[29];  // Non-phos = tot - phos
	double PKAC_CFTR = (p[88]/p[33])*y[17];    // (PKACFTRtot/PKAIItot)*PKAIIact
	double CFTRphos = p[39]*CFTRn*PKAC_CFTR*p[90]/(p[91]+p[39]*CFTRn);
	double CFTRdephos = p[89]*p[92]*p[39]*y[29]/(p[93] + p[39]*y[29]);
	ydot[29] = CFTRphos - CFTRdephos;
	
	// PLM module (included 09/18/12)
	double PLM = p[94]-y[30];
	double PLM_PHOSPH = p[95]*y[16]*PLM/(p[96]+PLM);
	double PLM_DEPHOSPH = p[97]*y[21]*y[30]/(p[98]+y[30]);
	ydot[30] = PLM_PHOSPH-PLM_DEPHOSPH;
	
	
	//// Gather odes
	// Need to convert all ydot terms that are ODEs (not DAEs) to miliseconds
	// odes = [4,5,6,7,8,9,13,14,15,19,20,23,24,25,26,29,30,31];
	// ydot(odes) = ydot(odes).*1e-3;
	for ( int i = 0; i < 31; i++ ) {
		ydot[i] = ydot[i] * 0.001;
	}
}
#endif