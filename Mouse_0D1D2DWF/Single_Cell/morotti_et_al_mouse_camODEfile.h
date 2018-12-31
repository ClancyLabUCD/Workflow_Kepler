/*
 function dydt = morotti_et_al_mouse_camODEfile(t,y,params)
 % This function describes the ODE's for CaM, CaMKII, and CaN
 %
 % Re-implemented by Mao-Tsuen Jeng for CPVT model
 % modified by Pei-Chi Yang on 06/22/2015
 //
*/

void morotti_et_al_mouse_camODEfile( double * y, double * params, double * dydt , pars_rec * pars1 );
void morotti_et_al_mouse_camODEfile( double * y, double * params, double * dydt , pars_rec * pars1 ) {
	
	//// Descriptions for state variables
	// Ca2CaM = y(1);  // 2 Ca bound to C terminal sites
	// Ca4CaM = y(2);  // 4 Ca bound
	// CaMB = y(3);   
	// Ca2CaMB = y(4);
	// Ca4CaMB = y(5);
	// Pb2 = y(6);     // probability of a Ca2CaM bound CaMKII subunit
	// Pb = y(7);      // probability of a Ca4CaM bound CaMKII subunit
	// Pt = y(8);      // probability of a Ca4CaM bound autophosphorylated CaMKII subunit
	// Pt2 = y(9);     // probability of a Ca2CaM bound autophosphorylated CaMKII subunit
	// Pa = y(10);     // probability of an autonomous autophosphorylated CaMKII subunit
	// CaMCaN = y(11);
	// Ca2CaMCaN = y(12);
	// Ca4CaMCaN = y(13);  // active calcineurin
	// description of intermediate variables
	// CaM- Ca free CaM
	////
	
	//paramsCell=mat2cell(params,ones(size(params,1),1),ones(size(params,2),1));
	//[K,Mg,CaMtot,Btot,CaMKIItot,CaNtot,PP1tot,Ca,cyclelength,compartment]=paramsCell{:};
	double K = params[0];
	double Mg = params[1];
	double CaMtot = params[2];
	double Btot = params[3];
	double CaMKIItot = params[4];
	double CaNtot = params[5];
	double PP1tot = params[6];
	double Ca = params[7];
	double cyclelength = params[8];
	double compartment = params[9];
	
	//yCell=mat2cell(y,ones(size(y,1),1),ones(size(y,2),1));
	//	[CaM,Ca2CaM,Ca4CaM,CaMB,Ca2CaMB,Ca4CaMB,...                 // variables 1-6
	//	 Pb2,Pb,Pt,Pt2,Pa,...                                    // variables 7-11
	//	 Ca4CaN,CaMCa4CaN,Ca2CaMCa4CaN,Ca4CaMCa4CaN]=yCell{:};   // variables 12-15
	double CaM = y[0];
	double Ca2CaM = y[1];
	double Ca4CaM = y[2];
	double CaMB = y[3];
	double Ca2CaMB = y[4];
	double Ca4CaMB = y[5];
	double Pb2 = y[6];
	double Pb = y[7];
	double Pt = y[8];
	double Pt2 = y[9];
	double Pa = y[10];             
	double Ca4CaN = y[11];
	double CaMCa4CaN = y[12];
	double Ca2CaMCa4CaN = y[13];
	double Ca4CaMCa4CaN = y[14];   
	
	// Parameters
	// Ca/CaM parameters
	double Kd02, Kd24;
	if ( Mg <= 1 ) {
		Kd02 = 0.0025 * ( 1. + K / 0.94 - Mg / 0.012 ) * ( 1. + K / 8.1 + Mg / 0.022 );  // [uM^2]
		Kd24 = 0.128 * ( 1. + K / 0.64 + Mg / 0.0014 ) * ( 1. + K / 13.0 - Mg / 0.153 ); // [uM^2]
	} else {
		Kd02 = 0.0025 * ( 1. + K / 0.94 - 1. / 0.012 + ( Mg - 1. ) / 0.060 ) * ( 1. + K / 8.1 + 1. / 0.022 + ( Mg - 1. ) / 0.068 );   // [uM^2]
		Kd24 = 0.128 * ( 1. + K / 0.64 + 1. / 0.0014 + ( Mg - 1. ) / 0.005 ) * ( 1. + K / 13.0 - 1. / 0.153 + ( Mg - 1. ) / 0.150 );  // [uM^2]
	}
	// end
	double k20 = 10.;               // [s^-1]      
	double k02 = k20 / Kd02;         // [uM^-2 s^-1]
	double k42 = 500.;              // [s^-1]      
	double k24 = k42 / Kd24;         // [uM^-2 s^-1]
	
	// CaM buffering (B) parameters
	double k0Boff = 0.0014;        // [s^-1] 
	double k0Bon = k0Boff/0.2;   // [uM^-1 s^-1] kon = koff/Kd
	double k2Boff = k0Boff/100;    // [s^-1] 
	double k2Bon = k0Bon;          // [uM^-1 s^-1]
	double k4Boff = k2Boff;        // [s^-1]
	double k4Bon = k0Bon;          // [uM^-1 s^-1]
								   // using thermodynamic constraints
	double k20B = k20/100; // [s^-1] thermo constraint on loop 1
	double k02B = k02;     // [uM^-2 s^-1] 
	double k42B = k42;     // [s^-1] thermo constraint on loop 2
	double k24B = k24;     // [uM^-2 s^-1]
	
	// CaMKII parameters
	// Wi Wa Wt Wp
	double kbi = 2.2;      // [s^-1] (Ca4CaM dissocation from Wb)
	double kib = kbi/33.5e-3; // [uM^-1 s^-1]
	double kib2 = kib;
	double kb2i = kib2*5;
	double kb24 = k24;
	double kb42 = k42*33.5e-3/5;
	double kpp1 = 1.72;    // [s^-1] (PP1-dep dephosphorylation rates)
	double Kmpp1 = 11.5;   // [uM]
	double kta = kbi/1000; // [s^-1] (Ca4CaM dissociation from Wt)
	double kat = kib;      // [uM^-1 s^-1] (Ca4CaM reassociation with Wa)
	double kt42 = k42*33.5e-6/5;
	double kt24 = k24;
	double kat2 = kib;
	double kt2a = kib*5;
	
	// CaN parameters
	double kcanCaoff = 1;              // [s^-1] 
	double kcanCaon = kcanCaoff/0.5;   // [uM^-1 s^-1] 
	double kcanCaM4on = 46;            // [uM^-1 s^-1]
	double kcanCaM4off = 1.3e-3;       // [s^-1]
	double kcanCaM2on = kcanCaM4on;
	double kcanCaM2off = 2508*kcanCaM4off;
	double kcanCaM0on = kcanCaM4on;
	double kcanCaM0off = 165*kcanCaM2off;
	double k02can = k02;
	double k20can = k20/165;
	double k24can = k24;
	double k42can = k20/2508;
	
	// CaM Reaction fluxes
	double rcn02 = k02*Ca*Ca*CaM - k20*Ca2CaM;
	double rcn24 = k24*Ca*Ca*Ca2CaM - k42*Ca4CaM;
	// CaM buffer fluxes
	double B = Btot - CaMB - Ca2CaMB - Ca4CaMB;
	double rcn02B = k02B*Ca*Ca*CaMB - k20B*Ca2CaMB;
	double rcn24B = k24B*Ca*Ca*Ca2CaMB - k42B*Ca4CaMB;
	double rcn0B = k0Bon*CaM*B - k0Boff*CaMB;
	double rcn2B = k2Bon*Ca2CaM*B - k2Boff*Ca2CaMB;
	double rcn4B = k4Bon*Ca4CaM*B - k4Boff*Ca4CaMB;
	// CaN reaction fluxes 
	double Ca2CaN = CaNtot - Ca4CaN - CaMCa4CaN - Ca2CaMCa4CaN - Ca4CaMCa4CaN;
	double rcnCa4CaN = kcanCaon*Ca*Ca*Ca2CaN - kcanCaoff*Ca4CaN;
	double rcn02CaN = k02can*Ca*Ca*CaMCa4CaN - k20can*Ca2CaMCa4CaN; 
	double rcn24CaN = k24can*Ca*Ca*Ca2CaMCa4CaN - k42can*Ca4CaMCa4CaN;
	double rcn0CaN = kcanCaM0on*CaM*Ca4CaN - kcanCaM0off*CaMCa4CaN;
	double rcn2CaN = kcanCaM2on*Ca2CaM*Ca4CaN - kcanCaM2off*Ca2CaMCa4CaN;
	double rcn4CaN = kcanCaM4on*Ca4CaM*Ca4CaN - kcanCaM4off*Ca4CaMCa4CaN;
	// CaMKII reaction fluxes
	double Pi = 1 - Pb2 - Pb - Pt - Pt2 - Pa;
	double rcnCKib2 = kib2*Ca2CaM*Pi - kb2i*Pb2;
	double rcnCKb2b = kb24*Ca*Ca*Pb2 - kb42*Pb;
	double rcnCKib = kib*Ca4CaM*Pi - kbi*Pb;
	double T = Pb + Pt + Pt2 + Pa;
	double kbt = 0.055*T + .0074*T*T + 0.015*T*T*T;
	double rcnCKbt = kbt*Pb - kpp1*PP1tot*Pt/(Kmpp1+CaMKIItot*Pt);
	double rcnCKtt2 = kt42*Pt - kt24*Ca*Ca*Pt2;
	double rcnCKta = kta*Pt - kat*Ca4CaM*Pa;
	double rcnCKt2a = kt2a*Pt2 - kat2*Ca2CaM*Pa;
	double rcnCKt2b2 = kpp1*PP1tot*Pt2/(Kmpp1+CaMKIItot*Pt2);
	double rcnCKai = kpp1*PP1tot*Pa/(Kmpp1+CaMKIItot*Pa);
	
	// CaM equations
	double dCaM = 1e-3*(-rcn02 - rcn0B - rcn0CaN);
	double dCa2CaM = 1e-3*(rcn02 - rcn24 - rcn2B - rcn2CaN + CaMKIItot * (-rcnCKib2 + rcnCKt2a) );
	double dCa4CaM = 1e-3*(rcn24 - rcn4B - rcn4CaN + CaMKIItot * (-rcnCKib+rcnCKta) );
	double dCaMB = 1e-3*(rcn0B-rcn02B);
	double dCa2CaMB = 1e-3*(rcn02B + rcn2B - rcn24B);
	double dCa4CaMB = 1e-3*(rcn24B + rcn4B);
	
	// CaMKII equations
	double dPb2 = 1e-3*(rcnCKib2 - rcnCKb2b + rcnCKt2b2); // Pb2
	double dPb = 1e-3*(rcnCKib + rcnCKb2b - rcnCKbt);    // Pb
	double dPt = 1e-3*(rcnCKbt-rcnCKta-rcnCKtt2);        // Pt
	double dPt2 = 1e-3*(rcnCKtt2-rcnCKt2a-rcnCKt2b2);     // Pt2
	double dPa = 1e-3*(rcnCKta+rcnCKt2a-rcnCKai);       // Pa
	
	// CaN equations
	double dCa4CaN = 1e-3*(rcnCa4CaN - rcn0CaN - rcn2CaN - rcn4CaN);                       // Ca4CaN
	double dCaMCa4CaN = 1e-3*(rcn0CaN - rcn02CaN);           // CaMCa4CaN
	double dCa2CaMCa4CaN = 1e-3*(rcn2CaN+rcn02CaN-rcn24CaN);    // Ca2CaMCa4CaN
	double dCa4CaMCa4CaN = 1e-3*(rcn4CaN+rcn24CaN);             // Ca4CaMCa4CaN
	
	//	dydt=[dCaM;dCa2CaM;dCa4CaM;dCaMB;dCa2CaMB;dCa4CaMB;...
	//		  dPb2;dPb;dPt;dPt2;dPa;dCa4CaN;dCaMCa4CaN;dCa2CaMCa4CaN;dCa4CaMCa4CaN];
	dydt[0] = dCaM;
	dydt[1] = dCa2CaM;
	dydt[2] = dCa4CaM;
	dydt[3] = dCaMB;
	dydt[4] = dCa2CaMB;
	dydt[5] = dCa4CaMB;
	dydt[6] = dPb2;
	dydt[7] = dPb;
	dydt[8] = dPt;
	dydt[9] = dPt2;
	dydt[10] = dPa;
	dydt[11] = dCa4CaN;
	dydt[12] = dCaMCa4CaN;
	dydt[13] = dCa2CaMCa4CaN;
	dydt[14] = dCa4CaMCa4CaN;
	
	// write to global variables for adjusting Ca buffering in EC coupling model
	double JCa = 1e-3 * ( 2. * CaMKIItot * ( rcnCKtt2 - rcnCKb2b ) 
						 - 2. * ( rcn02 + rcn24 + rcn02B + rcn24B 
								 + rcnCa4CaN + rcn02CaN + rcn24CaN ) ); // [uM/msec]
	
	// global JCaCyt JCaSL JCaDyad;
	if ( compartment == 0 ) {
		pars1->JCaCyt = JCa;
	} else if ( compartment == 1 ) {
		pars1->JCaSL = JCa;
	} else {
		pars1->JCaDyad = JCa;
	}
	// end
	
}
