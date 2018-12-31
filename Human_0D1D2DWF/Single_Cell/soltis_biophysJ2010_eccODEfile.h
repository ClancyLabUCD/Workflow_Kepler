/*
 *  soltis_biophysJ2010_eccODEfile.h
 *  Covert from Matlab code by Mao-Tsuen Jeng on 12/22/11.
 *
 *  Created by Pei-Chi Yang on 03/13/12 with Grandi Human model.
 *	& add cAMP effects on 01/03/12 (1)  - Iso && Ach
 *  (2) LCC model from SS model 
 *  Copyright 2011 __CCLab__. All rights reserved.
 *
 */



/*
 function ydot = soltis_biophysJ2010_eccODEfile(t,y,p)
 % This is the Shannon/Bers EC coupling model Biophys J. 2004, implemented
 % by Jeff Saucerman with help from Eleonora Grandi and Fei Wang.
 %
 % Author: Jeff Saucerman <jsaucerman@virginia.edu>
 % Copyright 2008, University of Virginia, All Rights Reserved
 %
 % Re-implemented by Anthony Soltis <ars7h@virginia.edu> to include 
 % regulation of Ca2+ fluxes and currents by CaMKII and PKA. Final version 
 % completed 07/21/10
 %
 % References:
 % -----------
 % JJ Saucerman and DM Bers, Calmodulin mediates differential
 % sensitivity of CaMKII and calcineurin to local Ca2+ in cardiac myocytes. 
 % Biophys J. 2008 Aug 8. [Epub ahead of print]
 */

void soltis_biophysJ2010_eccODEfile( double dt, double tt, double * y, double * p, double * ydot  , pars_rec * pars1 , double * allDrugs );

void soltis_biophysJ2010_eccODEfile( double dt, double tt, double * y, double * p, double * ydot  , pars_rec * pars1 , double * allDrugs ){
	
	dt = 0.5 * dt;
	
	//// Assign passed-in parameters
	double cycleLength = p[0];
	//double t = tt - cycleLength * floor( tt / cycleLength );
	
	// CaMKII phosphorylated targets (%)
	double LCC_CKp = p[1];
	double RyR_CKp = p[2];
	double PLB_CKp = p[3];
	
	// PKA phosphorylated targets (%)
	double LCCa_PKAp = p[4];
	double LCCb_PKAp = p[5];
	double PLB_PKAn = p[6];    // This is % non-phosphorylated PLB targets
	double RyR_PKAp = p[7];
	double TnI_PKAp = p[8];
	double IKs_PKAp = p[9];
	double ICFTR_PKAp = p[10];
	
	// Flag for CKII OE (tells code whether or not it should use Ito and INa
	// changes) - is either 0 (WT) or 1 (OE)
	double CKIIOE = p[11];
	
	double rest = p[12];
	
    double PLM_PKAn = p[13];
    
    double factor_Kr = p[14];
    
    // const double drug = 0. * (1e-6);
    const double drug = allDrugs[0];

    
    const double BlockGkr = allDrugs[1]; // ratio
    
    

	//// Model Parameters
	// Constants
	double R = 8314.;       // [J/kmol*K]  
	double Frdy = 96485.;   // [C/mol]  
	double Temp = 310.;     // [K] 310 K (37 C) for BT / 295 K (22 C) for RT
	double FoRT = Frdy / R / Temp;
	double Cmem = 1.3810e-10;   // [F] membrane capacitance
	double Qpow = ( Temp - 310 ) / 10 ;
	
	// Cell geometry
	double cellLength = 100.;     // cell length [um]
	double cellRadius = 10.25;   // cell radius [um]
	double junctionLength = 160e-3;  // junc length [um]
	double junctionRadius = 15e-3;   // junc radius [um]
	double distSLcyto = 0.45;    // dist. SL to cytosol [um]
	double distJuncSL = 0.5;  // dist. junc to SL [um]
	double DcaJuncSL = 1.64e-6;  // Dca junc to SL [cm^2/sec]
	double DcaSLcyto = 1.22e-6; // Dca SL to cyto [cm^2/sec]
	double DnaJuncSL = 1.09e-5;  // Dna junc to SL [cm^2/sec]
	double DnaSLcyto = 1.79e-5;  // Dna SL to cyto [cm^2/sec] 
	// double Vcell = pi * cellRadius^2 * cellLength * 1e-15;    // [L]
	double Vcell = pi * cellRadius * cellRadius * cellLength * 1e-15;    // [L]
	double Vmyo = 0.65 * Vcell; 
	double Vsr = 0.035 * Vcell; 
	double Vsl = 0.02 * Vcell; 
	double Vjunc = 0.0539 * .01 * Vcell; 
	double SAjunc = 20150. * pi * 2 * junctionLength * junctionRadius;  // [um^2]
	double SAsl = pi * 2 * cellRadius * cellLength;          // [um^2]
	double J_ca_juncsl = 1. / 1.2134e12; // [L/msec]
	double J_ca_slmyo = 1. / 2.68510e11; // [L/msec]
	double J_na_juncsl = 1. / (1.6382e12 / 3 * 100 ); // [L/msec] 
	double J_na_slmyo = 1. / (1.8308e10 / 3 * 100 );  // [L/msec] 
	
	// Fractional currents in compartments
	double Fjunc = 0.11;   
	double Fsl = 1. - Fjunc;
	double Fjunc_CaL = 0.9; 
	double Fsl_CaL = 1. - Fjunc_CaL;
	
	// Fixed ion concentrations     
	double Cli = 15.;   // Intracellular Cl  [mM]
	double Clo = 150.;  // Extracellular Cl  [mM]
	double Ko = 5.4;   // Extracellular K   [mM]
	double Nao = 140.;  // Extracellular Na  [mM]
	double Cao = 1.8;  // Extracellular Ca  [mM]
	double Mgi = 1.;    // Intracellular Mg  [mM]
	
	// Nernst Potentials
	double ena_junc = ( 1. / FoRT ) * log( Nao / y[31] );     // [mV]
	double ena_sl = ( 1. / FoRT ) * log( Nao / y[32] );       // [mV]
	double ek = ( 1. / FoRT ) * log( Ko / y[34] );	        // [mV]
	double eca_junc = ( 1. / FoRT / 2 ) * log( Cao / y[35] );   // [mV]
	double eca_sl = ( 1. / FoRT / 2 ) * log( Cao / y[36] );     // [mV]
	double ecl = ( 1. / FoRT ) * log( Cli / Clo );            // [mV]
	
	// Na transport parameters
	double GNa = 23;      // [mS/uF]
	double GNaB = 0.597e-3;    // [mS/uF] 
	
	
	
	double IbarNaK = 1.0 * 1.8;     // [uA/uF]
	
		
	double KmNaip = 11;         // [mM]
	double KmKo = 1.5;         // [mM]
	double Q10NaK = 1.63;  
	double Q10KmNai = 1.39;
	
	// K current parameters
	double pNaK = 0.01833;      
	double GtoSlow = 0.13*0.12;     // [mS/uF] 
	double GtoFast = 0.13*0.88;     // [mS/uF] 
	double gkp = 2 * 0.001;
	
	// Cl current parameters
	double GClCa = 0.5 * 0.109625;   // [mS/uF]
	double GClB = 0.6 * 9e-3;        // [mS/uF]
	double KdClCa = 100e-3;    // [mM]
	
	// I_Ca parameters 
	double pNa =  0.5 * 1.5e-8;       // [cm/sec]
	double pCa =   0.5 * 5.4e-4;       // [cm/sec]
	double pK =  0.5 * 2.7e-7;        // [cm/sec]
	double KmCa = 0.6e-3;      // [mM]
	double Q10CaL = 1.8;       
	
	// Ca transport parameters
	double IbarNCX = 1.0 * 4.5;      //  [uA/uF]5.5 before - 9 in rabbit
	

	const double KmCai = 3.59e-3;    //  [mM]
	const double KmCao = 1.3;        //  [mM]
	const double KmNai = 12.29;      //  [mM]
	const double KmNao = 87.5;       //  [mM]
	const double ksat = 0.32;        //  [none]  
	const double nu = 0.27;          //  [none]
	const double Kdact = 0.150e-3;   //  [mM] 
	const double Q10NCX = 1.57;      //  [none]
	double IbarSLCaP = 0.0673; //  IbarSLCaP FEI changed [uA/uF](2.2 umol/L cytosol/sec) jeff 0.093 [uA/uF]
	const double KmPCa = 0.5e-3;     //  [mM] 
	double GCaB = 5.513e-4;    //  [uA/uF] 3
	

	const double Q10SLCaP = 2.35;    //  [none]
	
	// SR flux parameters 
	double Q10SRCaP = 2.6;          // [none]
	double Vmax_SRCaP = 1.0 * 5.3114e-3; // [mM/msec] (mmol/L cytosol/msec)
	double Kmf = 0.246e-3;          // [mM]
	double Kmr = 1.7;               // [mM]L cytosol
	double hillSRCaP = 1.787;       // [mM]
	double ks = 25.;                 // [1/ms]      
	double koCa = 10.;               // [mM^-2 1/ms]    
	
	
	double kom = 0.06;              // [1/ms]     
	double kiCa = 0.5;              // [1/mM/ms]
	double kim = 0.005;             // [1/ms]
	double ec50SR = 0.45;           // [mM]
	
	// Buffering parameters
	const double Bmax_Naj = 7.561;       //  [mM] //  Bmax_Naj = 3.7; (c-code difference?)  //  Na buffering
	const double Bmax_Nasl = 1.65;       //  [mM]
	const double koff_na = 1e-3;         //  [1/ms]
	const double kon_na = 0.1e-3;        //  [1/mM/ms]
	const double Bmax_TnClow = 70e-3;    //  [mM]                      //  TnC low affinity
	double koff_tncl = 19.6e-3;    //  [1/ms] 
	const double kon_tncl = 32.7;        //  [1/mM/ms]
	const double Bmax_TnChigh = 140e-3;  //  [mM]                      //  TnC high affinity 
	const double koff_tnchca = 0.032e-3; //  [1/ms] 
	const double kon_tnchca = 2.37;      //  [1/mM/ms]
	const double koff_tnchmg = 3.33e-3;  //  [1/ms] 
	const double kon_tnchmg = 3e-3;      //  [1/mM/ms]
	const double Bmax_myosin = 140e-3;   //  [mM]                      //  Myosin buffering
	const double koff_myoca = 0.46e-3;   //  [1/ms]
	const double kon_myoca = 13.8;       //  [1/mM/ms]
	const double koff_myomg = 0.057e-3;  //  [1/ms]
	const double kon_myomg = 0.0157;     //  [1/mM/ms]
	const double Bmax_SR = 19. * .9e-3;     //  [mM] (Bers text says 47e-3) 19e-3
	const double koff_sr = 60e-3;        //  [1/ms]
	const double kon_sr = 100.;           //  [1/mM/ms]
	const double Bmax_SLlowsl = 37.4e-3 * Vmyo / Vsl;        //  [mM]    //  SL buffering
	const double Bmax_SLlowj = 4.6e-3 * Vmyo / Vjunc * 0.1;    //  [mM]    // Fei *0.1!!! junction reduction factor
	const double koff_sll = 1300e-3;     //  [1/ms]
	const double kon_sll = 100.;          //  [1/mM/ms]
	const double Bmax_SLhighsl = 13.35e-3 * Vmyo / Vsl; //13.4e-3 * Vmyo / Vsl;       //  [mM] 
	const double Bmax_SLhighj = 1.65e-3 * Vmyo / Vjunc * 0.1;  //  [mM] // Fei *0.1!!! junction reduction factor
	const double koff_slh = 30e-3;       //  [1/ms]
	const double kon_slh = 100.;          //  [1/mM/ms]
	const double Bmax_Csqn = 2.7; //140e-3 * Vmyo / Vsr;            //  [mM] //  Bmax_Csqn = 2.6;      //  Csqn buffering
	const double koff_csqn = 65.;         //  [1/ms] 
	const double kon_csqn = 100.;         //  [1/mM/ms] 
	
	
	// PKA-dependent phosphoregulation of TnI (increases Kd of TnC)
	double fracTnIpo = .0031;  // Derived quantity (TnI_PKAp(baseline)/TnItot)
	double fPKA_TnI = ( 1.45 - 0.45 * ( 1 - TnI_PKAp ) / ( 1 - fracTnIpo ) );
	koff_tncl = koff_tncl * fPKA_TnI;
 
	
	//// MEMBRANE CURRENTS
	//// Flags
	// Set CKIIflag to 1 for CKII OE, 0 otherwise
	int CKIIflag = CKIIOE;
	// Set ICa_MarkovFlag to 1 for Markov ICa, 0 otherwise
	int ICa_MarkovFlag = 1;
	// Set MarkovFlag to 1 for Markov INa, 0 otherwise
	int MarkovFlag = 0;
	// Set Ito to use either original params (=0) or Grandi Params (=1)
	int ItoFlag = 1;
	// I_Na: Fast Na Current
	// Max INa alterations with CKII hyperactivity as in Hund & Rudy 2008
	double inashift , alphaCKII , deltGbarNal_CKII ;  
	if ( CKIIflag == 1 ) {    
		inashift = -3.25;
		alphaCKII = -.18;
		deltGbarNal_CKII = 2.;
        
        
	} else {
		inashift = 0;
		alphaCKII = 0;
		deltGbarNal_CKII = 0;
		
	}
    //Unused in the Code
    //H-H Na model
    ydot[0] = 0;
    ydot[1] = 0;
    ydot[2] = 0;
    
    //---Replaced Markov-INa -----------------------------------------------------------------
    /////from "Jonathan D Moreno, Z Iris Zhu, Pei-Chi Yang, John R Bankston, Mao-Tsuen Jeng, Chaoyi Kang, Lianguo Wang, Jason D Bayer, David J Christini, Natalia A Trayanova, Crystal M Ripplinger, Robert S Kass, and Colleen E Clancy. A computational model to predict the effects of class I anti-arrhythmic drugs on ventricular rhythms. Sci Trans Med 3, 98ra83 (2011)"
	
	double Q10 = 3. ;
	
	double Tfactor = 1.0 / pow( Q10 , ( ( 37.0 - ( Temp - 273 ) ) / 10.0 ) );
	double pH = 7.4;
	double pKa = 9.3;
	
	double portion = 1.0 / ( 1 + pow( 10 , ( pH - pKa ) ) );
	double diffusion = 5500.;
	
	double drug_charged = drug * portion;
	double drug_neutral = drug * ( 1 - portion );
	double dd = -0.7;
	
	
	////Transition Rates for WT Channel
	double a11= Tfactor*8.5539/(7.4392e-2*exp(-y[38]/17.0)+ 2.0373e-1*exp(-y[38]/150));
	double a12= Tfactor*8.5539/(7.4392e-2*exp(-y[38]/15.0)+ 2.0373e-1*exp(-y[38]/150));
	double a13= Tfactor*8.5539/(7.4392e-2*exp(-y[38]/12.0)+ 2.0373e-1*exp(-y[38]/150));
	double b11= Tfactor*7.5215e-2*exp(-y[38]/20.3);
	double b12= Tfactor*2.7574*exp(-(y[38]-5)/20.3);
	double b13= Tfactor*4.7755e-1*exp(-(y[38]-10)/20.3);
	
	double a3 = Tfactor*5.1458e-6*exp(-y[38]/8.2471);
	double b3=Tfactor*6.1205*exp((y[38])/13.542);	
	
	// CaMKII shifting to hyperpolarized potentials ( 5mV)

	
		
	double a2= Tfactor*(13.370*exp(y[38]/43.749));
	double b2= ((a13*a2*a3)/(b13*b3));
	
	double a4 = 0*a2;
	double b4 = 0*a3;
	double a5= 0*a2;
	double b5 = 0*a3;
	
	double mu1;
	double mu2;
	
	
	
		mu1 = 2.0462e-07; mu2 = 8.9731e-04; 
	
	
	
	double ax = 3.4229e-2*a2;
	double bx = 1.7898e-2*a3;
	
    //WT Flec  bursting states, and just then optimized the kb0
    
    double ax1 =5.7839e-05 * ax;
    double bx1 =  1.6689e-08* bx;
    double a13c = 3.6324e-03 *a13;
    double a22 = 1.4847e+03 *a2;
    double b33 =  1.7352e-06* b3;
    double a33 = 6.7505e-05 * a3;
    double a44 =  2.4135e+00* a2;
    double b44 =  4.9001e-02* a3;
    double a55 = 0;
    double b55 = 0;
    
    
    
    double ax2 = 2.6126e-01 * ax;
    double a13n = 2.6452e+00 * a13;
    double a_22 =  4.2385e+01 * a2;
    double b_33 = 2.1181e+00 * b3;
    double a_44 =  1.0326e-03 * a2;
    double b_44 = 2.1378e-02 * a3;
    
    double a_55 = 0;
    double b_55 = 0;
	
	double kd0=11.2*(1e-6);
	double kd_open = kd0 * exp( dd * y[38] * Frdy / ( R * Temp ) );
	
	
	////charged drug
	double kon = drug_charged * diffusion;
	double koff = kd_open * diffusion;
	double kcoff = koff;
	double kcon = kon;
	
	//bursting drug	
	double kbon = kon;
	double kboff = 95.2165 *(1e-6)  *   exp( (dd*y[38]*Frdy) /(R*Temp))   *diffusion;	//This expression gives 6uM at -100mV (Canine VMs and Human HF myocytes)
	double kcbon = kbon;					
	double kcboff = kboff;
	
	double b13c, b22;
	if ( drug ==0 || drug_charged ==0 ) {
		b13c = 0;
	} else {
		b13c = (b13*kcon*koff*a13c)/(kon*kcoff*a13);
		// end
	}
	
	if ( b13c ==0 ) {
		b22 = 0;
	} else {
		b22 = ( a13c * a22 * a33 ) / ( b13c * b33 );
		// end
	}
	
	////neutral drug
	double k_on = drug_neutral * diffusion;
	double k_off = 400. * (1e-6) * diffusion;				
	double ki_on = k_on / 2;
	double ki_off = 5.4 * (1e-6) * diffusion;
	double kc_on = k_on / 2;
	double kc_off = 800 * (1e-6) * diffusion;				
	double a_33, b13n, b_22, bx2;
	
	if ( drug == 0 || drug_neutral == 0 ) { 
		a_33 = 0;
	} else {
		a_33 = ( ki_off * a3 * kc_on * b_33 ) / ( ki_on * kc_off * b3 );
		// end
	}
	
	if ( drug == 0 || drug_neutral == 0 ) {
		b13n = 0;
	} else {
		b13n = ( b13 * kc_on * a13n * k_off ) / ( kc_off * a13 * k_on );
		// end
	}
	
	if ( b13n == 0 ) {
		b_22 =0;
	} else {
		b_22 = ( a_33 * a13n * a_22 ) / ( b_33 * b13n );
		// end
	}
	
	if ( drug == 0 || drug_neutral == 0 ) {
		bx2 = 0;
	} else {
		bx2 = ( bx * k_on * ax2 * ki_off ) / ( ax * ki_on * k_off );
		// end
	}
	
	
	const double coef_O = (b13 + a2  + mu1 + kon + k_on + ax);	
	const double coef_C1 = (b12 + b3 + a13  + mu1 + kcon + kc_on);
	const double coef_C2 = (b11 + b3 + a12  + mu1 + kcon + kc_on);
	const double coef_C3 = (b3 + a11  + mu1 + kcon + kc_on);
	const double coef_IC3 = (a11 + a3 + ki_on);
	const double coef_IC2 = (b11 + a3 + a12 + ki_on);
	const double coef_IF = (b12 + b2 + a3 + a4 + ki_on);
	const double coef_IM1 = (b4 + a5);
	const double coef_IM2 = b5;
	const double coef_OS = (bx + ki_on);
	
	const double coef_BO = (mu2 + b13 + kbon + k_on);				//changed from kon to kbon
	const double coef_BC3 = (mu2 + a11 + kcbon + kc_on);			//changed from kcon to kcbon
	const double coef_BC2 = (mu2 + b11 + a12 + kcbon + kc_on);		//changed from kcon to kcbon
	const double coef_BC1 = (mu2 + b12 + a13 + kcbon + kc_on);		//changed from kcon to kcbon
	
	const double coef_DO = (koff + b13c + a22 + ax1 + mu1);
	const double coef_DC1 = (kcoff + b12 + b33 + a13c + mu1);
	const double coef_DC2 = (kcoff + b11 + b33 + a12 + mu1);
	const double coef_DC3 = (kcoff+ b33 + a11 + mu1);
	const double coef_DOS = bx1;
	const double coef_DIC3 = (a11 + a33);
	const double coef_DIC2 = (a33 + b11 + a12);
	const double coef_DIF = (a33 + b12 + a44 + b22);
	const double coef_DIM1 = ( b44 + a55 );
	const double coef_DIM2 = b55 ;
	
	const double coef_DBO = (kboff + b13 + mu2);
	const double coef_DBC3 = (kcboff + a11 + mu2);
	const double coef_DBC2 = (kcboff + b11 + a12 + mu2);
	const double coef_DBC1 = (kcboff + b12 + a13 + mu2);
	
	const double coef_D_O = (k_off + b13n + a_22 + ax2 + mu1);
	const double coef_D_C1 = (kc_off + b12 + b_33 + a13n + mu1);
	const double coef_D_C2 = (kc_off + b11 + b_33 + a12 + mu1);
	const double coef_D_C3 = (kc_off + b_33 + a11 + mu1);
	const double coef_D_OS = (bx2 + ki_off);
	const double coef_D_IC3 = (a_33 + a11 + ki_off);
	const double coef_D_IC2 = (a_33 + b11 + a12 + ki_off);
	const double coef_D_IF = (a_33 + a_44 + b_22 + b12 + ki_off);
	const double coef_D_IM1 = (b_44 + a_55);
	const double coef_D_IM2 = b_55;
	
	const double coef_D_BO = (k_off + b13 + mu2);
	const double coef_D_BC3 = (kc_off + a11 + mu2);
	const double coef_D_BC2 = (kc_off + b11 + a12+ mu2);
	const double coef_D_BC1 = (kc_off + b12 + a13 + mu2);	
	
	double yo83, yo84, yo85, yo86, yo87, yo88, yo89, yo90, yo91, yo92, yo93, yo94, yo95, yo96, yo97, yo98, yo99, yo100, yo101, yo102, yo103, yo104, yo105, yo106, yo107, yo108, yo109, yo110, yo111, yo112;
	double yo113, yo114, yo115, yo116, yo117, yo118, yo119, yo120, yo121, yo122, yo123, yo124;
	
	double O = y[83];
	double C1 = y[84];
	double C2 = y[85];
	double C3 = y[86];
	double IC3 = y[87];
	double IC2 = y[88];
	double IF = y[89];
	double IM1 = y[90];
	double IM2 = y[91];
	double OS = y[92];
	double BO = y[93];
	double BC3 = y[94];
	double BC2 = y[95];
	double BC1 = y[96];
	double DO = y[97];
	double DC1 = y[98];
	double DC2 = y[99];
	double DC3 = y[100];
	double DOS = y[101];
	double DIC3 = y[102];
	double DIC2 = y[103];
	double DIF = y[104];
	double DIM1 = y[105];
	double DIM2 = y[106];
	double DBO = y[107];
	double DBC3 = y[108];
	double DBC2 = y[109];
	double DBC1 = y[110];
	double D_O = y[111];
	double D_C1 = y[112];
	double D_C2 = y[113];
	double D_C3 = y[114];
	double D_OS = y[115];
	double D_IC3 = y[116];
	double D_IC2 = y[117];
	double D_IF = y[118];
	double D_IM1 = y[119];
	double D_IM2 = y[120];
	double D_BO = y[121];
	double D_BC3 = y[122];
	double D_BC2 = y[123];
	double D_BC1 = y[124];
	
	double dtinv = 1. / dt;
	double myeps = 1.e-100;
	double err1 = 1.;
	int iter1 = 0;
	
	while ( err1 > myeps && iter1 < 1000 ) {
		yo83 = O;
		yo84 = C1;
		yo85 = C2;
		yo86 = C3;
		yo87 = IC3;
		yo88 = IC2;
		yo89 = IF;
		yo90 = IM1;
		yo91 = IM2;
		yo92 = OS;
		yo93 = BO;
		yo94 = BC3;
		yo95 = BC2;
		yo96 = BC1;
		yo97 = DO;
		yo98 = DC1;
		yo99 = DC2;
		yo100 = DC3;
		yo101 = DOS;
		yo102 = DIC3;
		yo103 = DIC2;
		yo104 = DIF;
		yo105 = DIM1;
		yo106 = DIM2;
		yo107 = DBO;
		yo108 = DBC3;
		yo109 = DBC2;
		yo110 = DBC1;
		yo111 = D_O;
		yo112 = D_C1;
		yo113 = D_C2;
		yo114 = D_C3;
		yo115 = D_OS;
		yo116 = D_IC3;
		yo117 = D_IC2;
		yo118 = D_IF;
		yo119 = D_IM1;
		yo120 = D_IM2;
		yo121 = D_BO;
		yo122 = D_BC3;
		yo123 = D_BC2;
		yo124 = D_BC1;
		
		//Drug Free States
		O = ( y[83] + dt * ( a13 * C1 + b2 * IF + koff * DO + k_off * D_O + bx * OS + mu2 * BO) ) / ( 1. + dt * coef_O );		// O
		C1 = ( y[84] + dt * (a12 * C2 + a3 * IF + b13 * O  + kcoff * DC1 + kc_off * D_C1 + mu2 * BC1 ) ) / (1. + dt * coef_C1 );	// C1
		C2 = ( y[85] + dt * (a11 * C3 + a3 * IC2 + b12 * C1  + kcoff * DC2 + kc_off * D_C2 + mu2 * BC2) ) / ( 1. + dt * coef_C2 ); //C2	
		C3 = ( y[86] + dt * (a3 * IC3 + b11 * C2  + kcoff * DC3 + kc_off * D_C3 + mu2 * BC3 ) )/( 1.+ dt * coef_C3 ); //C3
		IC3 = ( y[87] + dt * (b3 * C3 + b11 * IC2 + ki_off * D_IC3 ) ) / ( 1. + dt * coef_IC3 );	
		IC2 = ( y[88] + dt * (a11 * IC3 + b3 * C2 + b12 * IF + ki_off * D_IC2 ) ) / ( 1. + dt * coef_IC2 );	
		IF = ( y[89] + dt * (a12 * IC2 + b3 * C1 + b4 * IM1 + a2 * O + ki_off * D_IF) ) / (1. + dt * coef_IF );	
		IM1 = ( y[90] + dt * (a4 * IF + b5 * IM2) ) / (1. + dt * coef_IM1 );	
		IM2 = ( y[91] + dt * (a5 * IM1 ) ) / ( 1. + dt * coef_IM2 );	
		OS = ( y[92] + dt * (ax * O + ki_off * D_OS) ) / (1. + dt * coef_OS );	
		BO = ( y[93] + dt * (mu1 * O + a13 * BC1 + kboff* DBO + k_off * D_BO) ) / (1. + dt * coef_BO );	
		BC3 = ( y[94] + dt * (mu1 * C3 + b11 * BC2 + kcboff * DBC3 + kc_off * D_BC3) ) / (1. + dt * coef_BC3 );	
		BC2 = ( y[95] + dt * (mu1 * C2 + a11 * BC3 + b12 * BC1 + kcboff * DBC2 + kc_off * D_BC2) ) / (1. + dt * coef_BC2 );	
		BC1 = ( y[96] + dt * (mu1 * C1 + a12 * BC2 + b13 * BO + kcboff * DBC1 + kc_off * D_BC1) ) / (1. + dt * coef_BC1 );
		
		
		//Charged Drug Bound States
		DO = ( y[97] + dt * (kon * O + a13c * DC1 + bx1 * DOS + b22 * DIF + mu2 * DBO ) ) / ( 1. + dt * coef_DO );	
		DC1 = ( y[98] + dt * (kcon * C1 + a12 * DC2 + a33 * DIF + b13c * DO  + mu2 * DBC1 ) ) /( 1.+ dt * coef_DC1 );	
		DC2 = ( y[99] + dt * (kcon * C2 + a11 * DC3 + a33 * DIC2 + b12 * DC1  + mu2 * DBC2 )) / ( 1. + dt * coef_DC2 );	
		DC3 = ( y[100] + dt * (kcon * C3 + b11 * DC2 + a33 * DIC3  + mu2 * DBC3 ) ) / ( 1. + dt * coef_DC3 );	
		DOS = ( y[101] + dt * (ax1 * DO )) / (1. + dt * coef_DOS );
		DIC3 = ( y[102] + dt * (b33 * DC3 + b11 * DIC2 ) ) / ( 1. + dt * coef_DIC3 );	
		DIC2 = ( y[103] + dt * (b33 * DC2 + a11 * DIC3 + b12 * DIF )) /(1. +dt * coef_DIC2 );	
		DIF = ( y[104] + dt * (b33 * DC1 + a12 * DIC2 + b44 * DIM1 + a22 * DO ))/( 1. + dt * coef_DIF );	
		DIM1 = ( y[105] + dt * (a44 * DIF + b55 * DIM2 )) / ( 1. + dt * coef_DIM1 );	
		DIM2 = ( y[106] + dt * (a55 * DIM1 )) /( 1. + dt * coef_DIM2 );	
		DBO = ( y[107] + dt * (kbon * BO + a13 * DBC1 + mu1 * DO )) / ( 1. + dt * coef_DBO );						//kbon stayed the same here	
		DBC3 = ( y[108] + dt * (kcbon * BC3 + b11 * DBC2 + mu1 * DC3 )) / (1. + dt * coef_DBC3 );
		DBC2 = ( y[109] + dt * (kcbon * BC2 + a11 * DBC3 + b12 * DBC1 + mu1 * DC2 )) / ( 1. + dt* coef_DBC2 );
		DBC1 = ( y[110] + dt * (kcbon * BC1 + a12 * DBC2 + b13 * DBO + mu1 * DC1 )) / (1. + dt* coef_DBC1 );
		
		//Need to check these	
		//Neutral Drug Bound States
		D_O = ( y[111] + dt * (k_on * O + a13n * D_C1 + b_22 * D_IF + bx2 * D_OS + mu2 * D_BO )) / ( 1. + dt * coef_D_O );
		D_C1 = ( y[112] + dt * (kc_on * C1 + a12 * D_C2 + a_33 * D_IF + b13n * D_O + mu2 * D_BC1  )) / ( 1. + dt * coef_D_C1 );	
		D_C2 = ( y[113] + dt * (kc_on * C2 + a11 * D_C3 + a_33 * D_IC2 + b12 * D_C1 + mu2 * D_BC2 )) / ( 1. + dt * coef_D_C2 );	
		D_C3 = ( y[114] + dt * (kc_on * C3 + a_33 * D_IC3 + b11 * D_C2 + mu2 * D_BC3 )) / ( 1. + dt * coef_D_C3 );
		D_OS = ( y[115] + dt * (ax2 * D_O + ki_on * OS )) / ( 1. + dt * coef_D_OS );	
		D_IC3 = ( y[116] + dt * (b_33 * D_C3 + b11 * D_IC2 + ki_on * IC3 )) / ( 1. + dt * coef_D_IC3 );	
		D_IC2 = ( y[117] + dt * (b_33 * D_C2 + a11 * D_IC3 + b12 * D_IF + ki_on * IC2 )) / ( 1. + dt * coef_D_IC2 );	
		D_IF = ( y[118] + dt * (b_33 * D_C1 + a12 * D_IC2 + b_44 * D_IM1 + a_22 * D_O + ki_on * IF )) / ( 1. + dt * coef_D_IF );
		D_IM1 = ( y[119] + dt * (a_44 * D_IF + b_55 * D_IM2 )) / ( 1. + dt * coef_D_IM1 );
		D_IM2 = ( y[120] + dt * (a_55 * D_IM1 )) / ( 1. + dt * coef_D_IM2 );
		D_BO = ( y[121] + dt * (k_on * BO + a13 * D_BC1 + mu1 * D_O )) / ( 1. + dt * coef_D_BO  );
		D_BC3 = ( y[122] + dt * (kc_on * BC3 + b11 * D_BC2 + mu1 * D_C3 )) / ( 1. + dt * coef_D_BC3  );
		D_BC2 = ( y[123] + dt * (kc_on * BC2 + a11 * D_BC3 + b12 * D_BC1 + mu1 * D_C2 )) / ( 1. + dt * coef_D_BC2 );
		D_BC1 = ( y[124] + dt * (kc_on * BC1 + a12 * D_BC2 + b13 * D_BO + mu1 * D_C1 )) / ( 1. + dt* coef_D_BC1 );
	
		err1 = ( fabs( yo83 - O ) 
				+ fabs( yo84 - C1 )
				+ fabs( yo85 - C2 )
				+ fabs( yo86 - C3 )
				+ fabs( yo87 - IC3 )
				+ fabs( yo88 - IC2 )
				+ fabs( yo89 - IF )
				+ fabs( yo90 - IM1 )
				+ fabs( yo91 - IM2 )
				+ fabs( yo92 - OS )
				+ fabs( yo93 - BO )
				+ fabs( yo94 - BC3 )
				+ fabs( yo95 - BC2 )
				+ fabs( yo96 - BC1 )
				+ fabs( yo97 - DO )
				+ fabs( yo98 - DC1 )
				+ fabs( yo99 - DC2 ) 
				+ fabs( yo100 - DC3 )
				+ fabs( yo101 - DOS )
				+ fabs( yo102 - DIC3 )
				+ fabs( yo103 - DIC2 )
				+ fabs( yo104 - DIF )
				+ fabs( yo105 - DIM1 )
				+ fabs( yo106 - DIM2 )
				+ fabs( yo107 - DBO )
				+ fabs( yo108 - DBC3 )
				+ fabs( yo109 - DBC2 )
				+ fabs( yo110 - DBC1 )
				+ fabs( yo111 - D_O )
				+ fabs( yo112 - D_C1 )
				+ fabs( yo113 - D_C2 )
				+ fabs( yo114 - D_C3 )
				+ fabs( yo115 - D_OS )
				+ fabs( yo116 - D_IC3 )
				+ fabs( yo117 - D_IC2 )
				+ fabs( yo118 - D_IF )
				+ fabs( yo119 - D_IM1 )
				+ fabs( yo120 - D_IM2 )
				+ fabs( yo121 - D_BO )
				+ fabs( yo122 - D_BC3 )
				+ fabs( yo123 - D_BC2 )
				+ fabs( yo124 - D_BC1 )
				
				);
		iter1++;
	}
	
	ydot[83] = ( O - y[83] ) * dtinv;//O
	ydot[84] = (C1 - y[84] ) * dtinv;//C1
	ydot[85] = ( C2 - y[85] ) * dtinv;//C2
	ydot[86] = ( C3 - y[86] ) * dtinv;//C3
	ydot[87] = ( IC3 - y[87] ) * dtinv;//IC3
	ydot[88] = ( IC2 - y[88] ) * dtinv;//IC2
	ydot[89] = ( IF - y[89] ) * dtinv;//IF
	ydot[90] = ( IM1 - y[90] ) * dtinv;//IM1
	ydot[91] = ( IM2 - y[91] ) * dtinv;//IM2
	ydot[92] = ( OS - y[92] ) * dtinv;//OS
	ydot[93] = ( BO - y[93] ) * dtinv;//BO
	ydot[94] = ( BC3 - y[94] ) * dtinv;//BC3
	ydot[95] = ( BC2 - y[95] ) * dtinv;//BC2
	ydot[96] = ( BC1 - y[96] ) * dtinv;//BC1
	ydot[97] = ( DO - y[97] ) * dtinv;//DO
	ydot[98] = ( DC1 - y[98] ) * dtinv;//DC1
	ydot[99] = ( DC2 - y[99] ) * dtinv;//DC2
	ydot[100] = ( DC3 - y[100] ) * dtinv;//DC3
	ydot[101] = ( DOS - y[101] ) * dtinv;//DOS
	ydot[102] = ( DIC3 - y[102] ) * dtinv;//DIC3
	ydot[103] = ( DIC2 - y[103] ) * dtinv;//DIC2
	ydot[104] = ( DIF - y[104] ) * dtinv;//DIF
	ydot[105] = ( DIM1 - y[105] ) * dtinv;//DIM1
	ydot[106] = ( DIM2 - y[106] ) * dtinv;//DIM2
	ydot[107] = ( DBO - y[107] ) * dtinv;//DBO
	ydot[108] = ( DBC3 - y[108] ) * dtinv;//DBC3
	ydot[109] = ( DBC2 - y[109] ) * dtinv;//DBC2
	ydot[110] = ( DBC1- y[110] ) * dtinv;//DBC1
	ydot[111] = ( D_O - y[111] ) * dtinv;//D_O
	ydot[112] = ( D_C1 - y[112] ) * dtinv;//D_C1
	ydot[113] = ( D_C2 - y[113] ) * dtinv;//D_C2
	ydot[114] = ( D_C3 - y[114] ) * dtinv;//D_C3
	ydot[115] = ( D_OS - y[115] ) * dtinv;//D_OS
	ydot[116] = ( D_IC3 - y[116] ) * dtinv;//D_IC3
	ydot[117] = ( D_IC2 - y[117] ) * dtinv;//D_IC2
	ydot[118] = ( D_IF - y[118] ) * dtinv;//D_IF
	ydot[119] = ( D_IM1 - y[119] ) * dtinv;//D_IM1
	ydot[120] = ( D_IM2 - y[120] ) * dtinv;//D_IM2
	ydot[121] = ( D_BO - y[121] ) * dtinv;//D_BO
	ydot[122] = ( D_BC3 - y[122] ) * dtinv;//D_BC3
	ydot[123] = ( D_BC2 - y[123] ) * dtinv;//D_BC2
	ydot[124] = ( D_BC1 - y[124] ) * dtinv;//D_BC1
	
	double I_Na_junc = Fjunc*GNa*(y[83]+y[93])*(y[38]-ena_junc);
	double I_Na_sl = Fsl*GNa*(y[83]+y[93])*(y[38]-ena_sl);
	double I_Na = I_Na_junc + I_Na_sl;
	

	ydot[46] = 0;
	
	ydot[57] = 0;
	ydot[58] = 0;
	
	
	// global I_Na_store
	pars1->I_Na_store = I_Na;
	
	
	double I_nabk_junc = Fjunc*GNaB*(y[38]-ena_junc) ;
	double I_nabk_sl = Fsl*GNaB*(y[38]-ena_sl) ;
	double I_nabk = I_nabk_junc+I_nabk_sl;
	
	
	//// I_nak: Na/K Pump Current
	const double sigma = (exp(Nao/67.3)-1)/7;
	const double fnak = 1./(1+0.1245*exp(-0.1*y[38]*FoRT)+0.0365*sigma*exp(-y[38]*FoRT));
	
	double fPKA_PLM = 1.0099 - 0.3551 * (PLM_PKAn/48);
	
	KmNaip=KmNaip*fPKA_PLM;
	
	const double I_nak_junc = ( Fjunc * IbarNaK * fnak * Ko 
							   / ( 1 + pow( (KmNaip / y[31]) , 4 ) ) 
							   / ( Ko + KmKo ) );
	const double I_nak_sl = ( Fsl * IbarNaK * fnak * Ko 
							 / ( 1 + pow( ( KmNaip / y[32] ), 4 ) ) 
							 / ( Ko + KmKo ) );
	const double I_nak = I_nak_junc+I_nak_sl;
	
	
	
	//// I_kr: Rapidly Activating K Current
	
	double IC50 = 1.5*(1e-6);
	double factor_flec = 1/(1+(drug/IC50));
	
	const double gkr =  factor_Kr * 0.035 * sqrt( Ko / 5.4 ) * factor_flec * 6 * BlockGkr ;
	
	
	const double xrss = 1. / ( 1 + exp( -(y[38]+10) / 5) );
	const double tauxr = 550./(1+exp((-22-y[38])/9))*6/(1+exp((y[38]-(-11))/9))+230/(1+exp((y[38]-(-40))/20));
	ydot[11] = (xrss-y[11])/tauxr;
	
	const double rkr = 1./(1+exp((y[38]+74)/24)); 
	
	const double I_kr = gkr * y[11] * rkr * (y[38]-ek);
	
    // Markov IKr from Lucia
//    const double GKr = 0.024 * 1.3 * 2.35;
//    
//    
//    double kO_D= 0.511e-6;
//    double rO_D= 1.6227e-5;
//    double kI_D= 0.511e-6;
//    double rI_D= 2.3181e-7;
//    
//    double ae=exp(24.335+0.0112*y[38]-25.914);
//    double be=exp(13.688-0.0603*y[38]-15.707);
//    double ai=exp(30.061-0.0312*y[38]-33.243);
//    double ain=exp(22.746-25.914);
//    double bin=exp(13.193-15.707);
//    double bi=exp(30.016+0.0223*y[38]-30.888)*pow((5.4/Ko),0.4);   //inactivation
//    double aa=exp(22.098+0.0365*y[38]-25.914);  //activation
//    double bb=exp(7.313-0.0399*y[38]-15.707);  //deactivation
//    
    
    //Calculation of k1 values using the initial differential equation (set as f(t,V)
    
    
    ydot[47] = 0;//(ai * y[51] + aa * y[48] + rO_D * y[52] - y[47] * (bi + bb + kO_D*Dofetilide));
    
    ydot[48] = 0;//(ain * y[49] + bb * y[47] - y[48] * (aa + bin));
    
    ydot[49] = 0;//(ae * y[50]+ bin * y[48] - y[49] * (be + ain));
    
    ydot[50] = 0;//(be * y[49] - ae* y[50]);
    
    ydot[51] = 0;//(bi * y[47] + rI_D * y[56] - y[51] * (ai + kI_D*Dofetilide));
    
    ydot[52] = 0;//( ai * y[56] + aa * y[53] + kO_D*Dofetilide * y[47] - y[52] * (bi + bb+rO_D));
    
    ydot[53] = 0;//(ain * y[54] + bb*y[52] - y[53] * (aa + bin));
    
    ydot[54] = 0;//(ae * y[55] + bin * y[53] - y[54] *(be + ain));
    
    ydot[55] = 0;//(be * y[54]- ae* y[55]);
    
    ydot[56] = 0;//(bi * y[52] + kI_D*Dofetilide*y[51] - y[56]*(ai +rI_D));
    
    
//    double I_kr = GKr * sqrt(Ko/5.4) * y[47] * (y[38] - ek); //!! NOT USING IN THIS CODE !!
    
    
    pars1->IKr_store = I_kr;

	//// I_ks: Slowly Activating K Current
	// Phosphoregulation of IKs by PKA parameters
	double fracIKspo = .0720;  // Derived quantity (IKs_PKAp(baseline)/IKstot)
	double fracIKsavail = ( 0.2 * (IKs_PKAp/fracIKspo)+0.8);
	
	double Xs05 = 1.5 * ( 2.0 - IKs_PKAp/fracIKspo); // Grandi version
	/// equations
	
	const double eks = (1./FoRT)*log((Ko+pNaK*Nao)/(y[34]+pNaK*y[33]));	
	
	const double gks_junc= fracIKsavail *0.0035 ;
	const double gks_sl= fracIKsavail *0.0035 ; // FRA
	const double xsss = 1. / (1+exp(-(y[38] + 3.8-Xs05)/14.25) ); //  fitting Fra
	const double tauxs=990.1/(1+exp(-(y[38]+2.436)/14.12) );
	ydot[12] = (xsss-y[12])/tauxs;
	const double I_ks_junc = Fjunc * gks_junc * y[12] * y[12] * (y[38]-eks);
	const double I_ks_sl = Fsl * gks_sl * y[12] * y[12] * (y[38]-eks);                                                                                                                                   
	
	double I_ks = (I_ks_junc+I_ks_sl);
	
	// global IKs_store
	pars1->IKs_store = I_ks;
	
	
	//// I_kp: Plateau K current
	const double kp_kp = 1./(1+exp(7.488-y[38]/5.98));
	const double I_kp_junc = Fjunc*gkp*kp_kp*(y[38]-ek);
	const double I_kp_sl = Fsl*gkp*kp_kp*(y[38]-ek);
	const double I_kp = I_kp_junc+I_kp_sl;
	
	//// I_to: Transient Outward K Current (slow and fast components)
	//  modified for human myocytes
	//	double GtoFast, GtoSlow;
		
	
	
	const double xtoss = 1./(1+exp(-(y[38]-19.0)/13));
	const double ytoss = 1./(1+exp((y[38]+19.5)/5));
	//double rtoss = 1. / ( 1 + exp( ( y[38] + 33.5 ) / 10 ) );
	const double tauxtos = 9./(1+exp((y[38]+3.0)/15))+0.5;
	double tauytos, taurtos, Py, Pr1, Pr2;
	
	if (ItoFlag == 0) {
		tauytos = 800./(1+exp((y[38]+60.0)/10))+30;
	} else if ( ItoFlag == 1 && CKIIflag == 0) {
		Py = 182; 
		Pr1 = 8085; 
		Pr2 = 313;                //Normal
		tauytos = Py/(1+exp((y[38]+33.5)/10))+1;
		taurtos = Pr1/(1+exp((y[38]+33.5)/10))+Pr2; 
		
	} else if (ItoFlag == 1 && CKIIflag == 1) {
		Py = 15; 
		Pr1 = 3600; 
		Pr2 = 500; 
		GtoSlow = GtoSlow*1.5;  // CKII OE
		tauytos = Py/(1+exp((y[38]+33.5)/10))+1;
		taurtos = Pr1/(1+exp((y[38]+33.5)/10))+Pr2;
	}
	
	ydot[7] = ( xtoss - y[7] ) / tauxtos;
	ydot[8] = ( ytoss - y[8] ) / tauytos;
	ydot[39]=0;
	const double I_tos = GtoSlow * y[7] * y[8] * (y[38] - ek);    //  [uA/uF]
	
	const double tauxtof = 8.5 * exp(- pow( ((y[38]+45)/50), 2) )+0.5;
	const double tauytof = 85 * exp( (- pow( (y[38]+40), 2) / 220 ) ) + 7;
	ydot[9] = ( xtoss - y[9] ) / tauxtof;
	ydot[10] = ( ytoss - y[10] ) / tauytof;
	const double I_tof = GtoFast * y[9] * y[10] * ( y[38] - ek );
	const double I_to = I_tos + I_tof;
	
	// global I_to_store
	pars1->I_to_store[0] = I_to;     // Total I_to
	pars1->I_to_store[1] = I_tof;    // Fast Component
	pars1->I_to_store[2] = I_tos;    // Slow component
	
	//// I_ki = I_k1: Time-Independent K Current
	const double aki = 1.02/(1+exp(0.2385*(y[38]-ek-59.215)));
	const double bki =(0.49124*exp(0.08032*(y[38]+5.476-ek)) + exp(0.06175*(y[38]-ek-594.31))) /(1 + exp(-0.5143*(y[38]-ek+4.753)));
	
	double kiss = aki/(aki+bki);
	
	double GK1 = 0.35 * sqrt(Ko/5.4);
	
		
	const double I_ki = GK1 *kiss*(y[38]-ek);	// global I_K1_store
	pars1->I_K1_store = I_ki;
	
	
	//// I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current
	const double I_ClCa_junc = Fjunc*GClCa/(1+KdClCa/y[35])*(y[38]-ecl);
	const double I_ClCa_sl = Fsl*GClCa/(1+KdClCa/y[36])*(y[38]-ecl);
	const double I_ClCa = I_ClCa_junc+I_ClCa_sl;
	const double I_Clbk = GClB*(y[38]-ecl);
	
	//// Original H-H formulation for LCC - unused if ICa_MarkovFlag = 1
	double dss = 1/(1+exp(-(y[38]+5)/6.0));
	double taud = dss*(1-exp(-(y[38]+5)/6.0))/(0.035*(y[38]+5));
	double fss = 1/(1+exp((y[38]+35)/9))+0.6/(1+exp((50-y[38])/20));
	double tauf = 1/(0.0197*exp( -pow( (0.0337*(y[38]+14.5)), 2) )+0.02);
	
	ydot[3] = 0;
	ydot[4] = 0;
	ydot[5] = 0;
	ydot[6] = 0;
	double ibarca_j = pCa*4*(y[38]*Frdy*FoRT) * (0.341*y[35]*exp(2*y[38]*FoRT)-0.341*Cao) /(exp(2*y[38]*FoRT)-1);
	double ibarca_sl = pCa*4*(y[38]*Frdy*FoRT) * (0.341*y[36]*exp(2*y[38]*FoRT)-0.341*Cao) /(exp(2*y[38]*FoRT)-1);
	double ibark = pK*(y[38]*Frdy*FoRT)*(0.75*y[34]*exp(y[38]*FoRT)-0.75*Ko) /(exp(y[38]*FoRT)-1);
	double ibarna_j = pNa*(y[38]*Frdy*FoRT) *(0.75*y[31]*exp(y[38]*FoRT)-0.75*Nao)  /(exp(y[38]*FoRT)-1);
	double ibarna_sl = pNa*(y[38]*Frdy*FoRT) *(0.75*y[32]*exp(y[38]*FoRT)-0.75*Nao)  /(exp(y[38]*FoRT)-1);
	//	
	double I_Ca_junc1 = (Fjunc_CaL*ibarca_j*y[3]*y[4]*(1-y[5])*pow(Q10CaL,Qpow))*0.45;
	double I_Ca_sl1 = (Fsl_CaL*ibarca_sl*y[3]*y[4]*(1-y[6])*pow(Q10CaL,Qpow))*0.45;
	double I_CaK1 = (ibark*y[3]*y[4]*(Fjunc_CaL*(1-y[5])+Fsl_CaL*(1-y[6]))*pow(Q10CaL,Qpow))*0.45;
	double I_CaNa_junc1 = (Fjunc_CaL*ibarna_j*y[3]*y[4]*(1-y[5])*pow(Q10CaL,Qpow))*0.45;
	double I_CaNa_sl1 = (Fsl_CaL*ibarna_sl*y[3]*y[4]*(1-y[6])*pow(Q10CaL,Qpow))*.45;
	
	
	//// LCC MARKOV MODEL - based on Mahajan et al. (2008) ////
	// This portion contains Markov state transitions for four channel types:
	// 'mode 1' channels in the junction and sl and 'mode 2' channels in the
	// same two compartments. Markov state transitions are computed for each
	// channel type independently - total currents are the sum of the two
	// channel types in each compartment (i.e. ICatot_junc = ICa_mode1_junc +
	// ICa_mode2_junc). Ca-dependent transition rates differ between juncitonal
	// and sl channels, whereas closing rate (r2) is adjusted to define mode1
	// vs. mode2 channels. Parameters determined through microscopic
	// reversibility are redefined to preserve constraint.
	
	// CaMKII shifts distribution of junctional and subsarcolemmal channels to 
	// either mode 2 at the expense of mode 1 channels (i.e. 10% mode 2 results 
	// in 90% mode 1).
	
	// PKA alters overall availability of channels (favail term that changes
	// overall scaling factor for currents) and also shifts distribution of
	// mode1/2 channels. PKA actions act on both junctional and sarcolemmal
	// channels.
	
	// To allow for CDI KO
	double cajLCC = y[35];
	double caslLCC = y[36];
	
	// LCC Current Fixed Parameters
	//pCa = 5.4e-4;       // [cm/sec] - Ca permeability
	double taupo = 1.;          // [ms] - Time constant of activation
	double TBa = 450.;          // [ms] - Time constant
	double s1o = .0221;
	double k1o = .03;
	double kop = 2.5e-3;       // [mM]
	double cpbar = 8e-3;       // [mM]
	double tca = 78.0312;
	double ICa_scale = 5.25; // to match human data
	double recoveryReduc = 3.;
	
	////// PKA PHOSPHOREGULATION OF LCC AVAILABLILITY (beta subunit phosph) ////////////
	double fracLCCbpo = .0328; // Derived quantity - (LCCbp(baseline)/LCCbtot)
	double favail = 1. * ( .017 * LCCb_PKAp / fracLCCbpo + 0.983);   // Test (max x1.5 pCa)
	ICa_scale =  ICa_scale * favail;
	
	// Voltage- and Ca-dependent Parameters
	double poss = 1. / ( 1. + exp( -y[38] / 8 ) );
	//	double fcaj = 1. / ( 1. + ( kop / cajLCC )^3 );            
	double fcaj = 1. / ( 1. + pow( (kop/cajLCC),3) );            
	double Rv = 10. + 4954. * exp( y[38] / 15.6 );
	double PrLCC = 1. - 1. / ( 1. + exp( -( y[38] + 40. ) / 4 ) );     
	double PsLCC = 1. / ( 1. + exp( -( y[38] + 40. ) / 11.32 ) );
	//	double TCaj = ( tca + 0.1 * ( 1. + ( cajLCC / cpbar )^2 ) ) / ( 1. + ( cajLCC / cpbar )^2 ); 
	double TCaj = ( tca + 0.1 * ( 1. + pow ( ( cajLCC / cpbar ),2) ) ) / ( 1. + pow ( ( cajLCC / cpbar ),2) ); 
	double tauCaj = ( Rv - TCaj ) * PrLCC + TCaj;     
	double tauBa = ( Rv - TBa ) * PrLCC + TBa;
	
	// Tranisition Rates (20 rates)
	double alphaLCC = poss / taupo;
	double betaLCC = ( 1. - poss ) / taupo;
	double r1 = 0.3;                               // [1/ms] - Opening rate
	double r2 = 3.;                                 // [1/ms] - closing rate
	double s1 = s1o * fcaj; 
	double s1p = .00195;                           // [ms] - Inactivation rate
	double k1 = k1o * fcaj;  
	double k1p = .00413;                           // [ms] - Inactivation rate
	double k2 = 1e-4;                              // [ms] - Inactivation rate
	double k2p = .00224;                           // [ms] - Inactivation rate
	double s2 = s1 * ( k2 / k1 ) * ( r1 / r2 );
	double s2p = s1p * ( k2p / k1p ) * ( r1 / r2 );
	double k3 = exp( -( y[38] + 40 ) / 3 ) / ( 3. * ( 1. + exp( -( y[38] + 40 ) / 3 ) ) );
	double k3p = k3;
	double k5 = ( 1. - PsLCC ) / tauCaj;
	double k6 = ( fcaj * PsLCC ) / tauCaj;
	double k5p = ( 1. - PsLCC ) / tauBa;
	
	// Recovery terms
	k5 = k5 / recoveryReduc;
	k5p = k5p / recoveryReduc;
	
	double k6p = PsLCC / tauBa;
	double k4 = k3 * ( alphaLCC / betaLCC ) * ( k1 / k2 ) * ( k5 / k6 );
	double k4p = k3p * ( alphaLCC / betaLCC ) * ( k1p / k2p ) * ( k5p / k6p );
	
	// global gates
	pars1->gates[1] = s1;
	pars1->gates[2] = k1;
	
	// State transitions for MODE 1 junctional LCCs //////
	// O = no differential; C2 = 60; C1 = 61; I1Ca = 62; I2Ca = 63;
	// I1Ba = 64; I2Ba = 65;
	double Po_LCCj_m1 = 1.0 - y[59] - y[60] - y[61] - y[62] - y[63] - y[64];                                           // O_m1j
	ydot[59] = betaLCC * y[60] + k5 * y[62] + k5p * y[64] - ( k6 + k6p + alphaLCC ) * y[59];                      // C2_m1j
	ydot[60] = alphaLCC * y[59] + k2 * y[61] + k2p * y[63] + r2 * Po_LCCj_m1 - ( r1 + betaLCC + k1 + k1p ) * y[60];   // C1_m1j
	ydot[61] = k1 * y[60] + k4 * y[62] + s1 * Po_LCCj_m1 - ( k2 + k3 + s2 ) * y[61];                              // I1Ca_m1j
	ydot[62] = k3 * y[61] + k6 * y[59] - ( k4 + k5 ) * y[62];                                                 // I2Ca_m1j
	ydot[63] = k1p * y[60] + k4p * y[64] + s1p * Po_LCCj_m1 - ( k2p + k3p + s2p ) * y[63];                        // I1Ba_m1j
	ydot[64] = k3p * y[63] + k6p * y[59] - ( k5p + k4p ) * y[64];                                             // I2Ba_m1j
	double ibarca_jm1 = ( 4. * pCa * y[38] * Frdy * FoRT ) * ( .001 * exp( 2 * y[38] * FoRT ) - 0.341 * Cao ) / ( exp( 2. * y[38] * FoRT ) - 1 );
	//	double I_Ca_junc_m1 = ( Fjunc_CaL * ibarca_jm1 * Po_LCCj_m1 * Q10CaL ^ Qpow ) * ICa_scale;
	double I_Ca_junc_m1 = ( Fjunc_CaL * ibarca_jm1 * Po_LCCj_m1 * pow( Q10CaL , Qpow ) ) * ICa_scale;
	
	////// Re-define all parameters as mode 2 specific parameters //////
	double s1om2 = .0221;
	double k1om2 = .03;
	double kopm2 = 2.5e-3;
	double cpbarm2 = 8e-3;
	double tcam2 = 78.0312;
	
	double possm2 = 1. / ( 1. + exp( -y[38] / 8 ) );
	//	double fcajm2 = 1. / ( 1. + ( kopm2 / cajLCC ) ^ 3 );    // Depends on junctional Ca
	double fcajm2 = 1. / ( 1. + pow( ( kopm2 / cajLCC ) , 3 ) );    // Depends on junctional Ca
	double Rvm2 = 10. + 4954 * exp( y[38] / 15.6 );
	double PrLCCm2 = 1. - 1. / ( 1. + exp( -( y[38] + 40 ) / 4 ) );     // Correct version I believe
	double PsLCCm2 = 1. / ( 1. + exp( -( y[38] + 40 ) / 11.32 ) );
	//	double TCajm2 = ( tcam2 + 0.1 * ( 1. + ( cajLCC / cpbarm2 ) ^ 2 ) ) / ( 1. + ( cajLCC / cpbarm2 ) ^ 2 ); // Caj dependent
	double TCajm2 = ( tcam2 + 0.1 * ( 1. +pow ( ( cajLCC / cpbarm2 ),2) ) ) / ( 1. + pow (( cajLCC / cpbarm2 ),2) ); // Caj dependent
	double tauCajm2 = ( Rvm2 - TCajm2 ) * PrLCCm2 + TCajm2;     // Caj dependence
	double tauBam2 = ( Rvm2 - TBa ) * PrLCCm2 + TBa;
	
	double alphaLCCm2 = possm2 / taupo;
	double betaLCCm2 = ( 1. - possm2 ) / taupo;
	double r1m2 = 0.3;                               // [1/ms] - Opening rate
	double r2m2 = 3. / 10;                                 // [1/ms] - closing rate
	double s1m2 = s1om2 * fcajm2; 
	double s1pm2 = .00195;                           // [ms] - Inactivation rate
	double k1m2 = k1om2 * fcajm2; 
	double k1pm2 = .00413;                           // [ms] - Inactivation rate
	double k2m2 = 1e-4;                              // [ms] - Inactivation rate
	double k2pm2 = .00224;                           // [ms] - Inactivation rate
	double s2m2 = s1m2*(k2m2/k1m2)*(r1m2/r2m2);
	double s2pm2 = s1pm2*(k2pm2/k1pm2)*(r1m2/r2m2);
	double k3m2 = exp( -( y[38] + 40. ) / 3. ) / ( 3. * ( 1. + exp( -( y[38] + 40. ) / 3. ) ) );
	double k3pm2 = k3m2;
	double k5m2 = ( 1. - PsLCCm2 ) / tauCajm2;
	double k6m2 = ( fcajm2 * PsLCCm2 ) / tauCajm2;
	double k5pm2 = ( 1. - PsLCCm2 ) / tauBam2;
	k5m2 = k5m2 / recoveryReduc;      // reduced for recovery
	k5pm2 = k5pm2 / recoveryReduc;    // reduced for recovery    
	double k6pm2 = PsLCCm2 / tauBam2;
	double k4m2 = k3m2 * ( alphaLCCm2 / betaLCCm2 ) * ( k1m2 / k2m2 ) * ( k5m2 / k6m2 );
	double k4pm2 = k3pm2 * ( alphaLCCm2 / betaLCCm2 ) * ( k1pm2 / k2pm2 ) * ( k5pm2 / k6pm2 );
	
	////// State transitions for MODE 2 junctional LCCs //////
	// O = no differential; C2 = 66; C1 = 67; I1Ca = 68; I2Ca = 69;
	// I1Ba = 70; I2Ba = 71;
	double Po_LCCj_m2 = 1.0-y[65]-y[66]-y[67]-y[68]-y[69]-y[70];                                                           // O_m2j
	ydot[65] = betaLCCm2 * y[66] + k5m2 * y[68] + k5pm2 * y[70] - ( k6m2 + k6pm2 + alphaLCCm2 ) * y[65];                          // C2_m2j
	ydot[66] = alphaLCCm2*y[65] + k2m2*y[67] + k2pm2*y[69] + r2m2*Po_LCCj_m2 - (r1m2+betaLCCm2+k1m2+k1pm2)*y[66];   // C1_m2j
	ydot[67] = k1m2*y[66] + k4m2*y[68] + s1m2*Po_LCCj_m2 - (k2m2+k3m2+s2m2)*y[67];                                  // I1Ca_m2j
	ydot[68] = k3m2*y[67] + k6m2*y[65] - (k4m2+k5m2)*y[68];                                                         // I2Ca_m2j
	ydot[69] = k1pm2*y[66] + k4pm2*y[70] + s1pm2*Po_LCCj_m2 - (k2pm2+k3pm2+s2pm2)*y[69];                            // I1Ba_m2j
	ydot[70] = k3pm2*y[69] + k6pm2*y[65] - (k5pm2+k4pm2)*y[70];                                                     // I2Ba_m2j
	double ibarca_jm2 = ( 4. * pCa * y[38] * Frdy * FoRT ) * ( .001 * exp( 2. * y[38] * FoRT ) - 0.341 * Cao ) / ( exp( 2. * y[38] * FoRT ) - 1 );
	//	double I_Ca_junc_m2 = ( Fjunc_CaL * ibarca_jm2 * ( Po_LCCj_m2 ) * Q10CaL ^ Qpow ) * ICa_scale;
	double I_Ca_junc_m2 = ( Fjunc_CaL * ibarca_jm2 * ( Po_LCCj_m2 ) * pow( Q10CaL , Qpow ) ) * ICa_scale;
	
	////// CaMKII AND PKA-DEPENDENT SHIFTING OF DYADIC LCCS TO MODE 2 ////////
	double fpkam2 = 0.1543 * LCCa_PKAp - .0043;  // Assumes max phosphorylation results in 15% mode 2 channels
	double fckiim2 = LCC_CKp * .1;               // Assumes max phosphorylation results in 10% mode 2 channels
	// Sum up total fraction of CKII and PKA-shifted mode 2 channels
	double junc_mode2 = fckiim2 + fpkam2;
	// Total junctional ICa
	double I_Ca_junc2 = ( 1. - junc_mode2 ) * I_Ca_junc_m1 + junc_mode2 * I_Ca_junc_m2;
	
	////// SUB-SARCOLEMMAL LCCs //////
	
	// Re-assign necessary params to be Casl sensitive
	//	double fcasl = 1. / ( 1. + ( kop / caslLCC ) ^ 3 );    // Depends on sl Ca
	double fcasl = 1. / ( 1. + pow( ( kop / caslLCC ) , 3 ) );    // Depends on sl Ca
	//	double TCasl = ( tca + 0.1 * ( 1. + ( caslLCC / cpbar ) ) ^ 2 ) / ( 1. + ( caslLCC / cpbar ) ^ 2 );
	double TCasl = ( tca + 0.1 * pow( ( 1. + ( caslLCC / cpbar ) ) , 2 ) ) / ( 1. + pow( ( caslLCC / cpbar ) , 2 ) );
	double tauCasl = ( Rv - TCasl ) * PrLCC + TCasl;
	
	// Re-assign necessary rates to be Casl sensitive
	double s1sl = s1o * fcasl;
	double k1sl = k1o * fcasl;
	double s2sl = s1sl * ( k2 / k1sl ) * ( r1 / r2 );
	double s2psl = s1p * ( k2p / k1p ) * ( r1 / r2 );
	double k5sl = ( 1. - PsLCC ) / tauCasl;
	k5sl = k5sl / recoveryReduc;  // Reduced for recovery
	double k6sl = ( fcasl * PsLCC ) / tauCasl;
	double k4sl = k3 * ( alphaLCC / betaLCC ) * ( k1sl / k2 ) * ( k5sl / k6sl );
	double k4psl = k3p * ( alphaLCC / betaLCC ) * ( k1p / k2p ) * ( k5p / k6p );
	
	// State transitions for 'mode 1' sarcolemmal LCCs
	// O = no differential; C2 = 72; C1 = 73; I1Ca = 74; I2Ca = 75;
	// I1Ba = 76; I2Ba = 77;
	double Po_LCCsl_m1 = 1. - y[71] - y[72] - y[73] - y[74] - y[75] - y[76];                                                // O_m1sl
	ydot[71] = betaLCC*y[72] + k5sl*y[74] + k5p*y[76] - (k6sl+k6p+alphaLCC)*y[71];                      // C2_m1sl
	ydot[72] = alphaLCC*y[71] + k2*y[73] + k2p*y[75] + r2*Po_LCCsl_m1 - (r1+betaLCC+k1sl+k1p)*y[72];    // C1_m1sl
	ydot[73] = k1sl*y[72] + k4sl*y[74] + s1sl*Po_LCCsl_m1 - (k2+k3+s2sl)*y[73];                         // I1Ca_m1sl
	ydot[74] = k3*y[73] + k6sl*y[71] - (k4sl+k5sl)*y[74];                                               // I2Ca_m1sl
	ydot[75] = k1p*y[72] + k4psl*y[76] + s1p*Po_LCCsl_m1 - (k2p+k3p+s2psl)*y[75];                       // I1Ba_m1sl
	ydot[76] = k3p*y[75] + k6p*y[71] - (k5p+k4psl)*y[76];                                               // I2Ba_m1sl
	double ibarca_slm1 = ( 4. * pCa * y[38] * Frdy * FoRT ) * ( .001 * exp( 2. * y[38] * FoRT ) - 0.341 * Cao ) / ( exp( 2. * y[38] * FoRT ) - 1 );
	//	double I_Casl_m1 = ( Fsl_CaL * ibarca_slm1 * Po_LCCsl_m1 * Q10CaL ^ Qpow ) * ICa_scale;
	double I_Casl_m1 = ( Fsl_CaL * ibarca_slm1 * Po_LCCsl_m1 * pow( Q10CaL , Qpow ) ) * ICa_scale;
	
	// Adjust closing rate for 'mode 2' sarcolemmal LCCs
	double r2slm2 = r2m2;
	double s2slm2 = s1sl*(k2/k1sl)*(r1/r2slm2);
	double s2pslm2 = s1p*(k2p/k1p)*(r1/r2slm2);
	
	////// State transitions for mode 2 sarcolemmal LCCs
	// O = no differential; C2 = 78; C1 = 79; I1Ca = 80; I2Ca = 81; I1Ba = 82; I2Ba = 83
	double Po_LCCsl_m2 = 1. - y[77]-y[78]-y[79]-y[80]-y[81]-y[82];                                                // O_m2sl
	ydot[77] = betaLCC*y[78] + k5sl*y[80] + k5p*y[82] - (k6sl+k6p+alphaLCC)*y[77];                      // C2_m2sl
	ydot[78] = alphaLCC*y[77] + k2*y[79] + k2p*y[81] + r2slm2*Po_LCCsl_m2 - (r1+betaLCC+k1sl+k1p)*y[78];// C1_m2sl
	ydot[79] = k1sl*y[78] + k4sl*y[80] + s1sl*Po_LCCsl_m2 - (k2+k3+s2slm2)*y[79];                       // I1Ca_m2sl
	ydot[80] = k3*y[79] + k6sl*y[77] - (k4sl+k5sl)*y[80];                                               // I2Ca_m2sl
	ydot[81] = k1p*y[78] + k4psl*y[82] + s1p*Po_LCCsl_m2 - (k2p+k3p+s2pslm2)*y[81];                     // I1Ba_m2sl
	ydot[82] = k3p*y[81] + k6p*y[77] - (k5p+k4psl)*y[82];                                               // I2Ba_m2sl
	double ibarca_slm2 = ( 4. * pCa * y[38] * Frdy * FoRT ) * ( .001 * exp( 2. * y[38] * FoRT ) - 0.341 * Cao ) / ( exp( 2. * y[38] * FoRT ) - 1 );
	//	cout << pCa << "\t" <<  y[38] << "\t" <<  Frdy << "\t" << FoRT << endl;
	//	cout << ( exp( 2. * y[38] * FoRT ) - 1 ) << endl;
	
	//	double I_Casl_m2 = ( Fsl_CaL * ibarca_slm2 * Po_LCCsl_m2 * Q10CaL ^ Qpow ) * ICa_scale;
	double I_Casl_m2 = ( Fsl_CaL * ibarca_slm2 * Po_LCCsl_m2 * pow( Q10CaL , Qpow ) ) * ICa_scale;
	//	cout << Fsl_CaL << "\t" <<  ibarca_slm2 << "\t" <<  Po_LCCsl_m2 << "\t" << ICa_scale << endl;
	
	// Sum mode 1 and mode 2 sl channels for total sl current
	double fckiim2_sl = 0; // Set to zero since SL LCCp by CaMKII is negligible
	double sl_mode2 = fckiim2_sl + fpkam2;
	double I_Ca_sl2 = ( 1. - sl_mode2 ) * I_Casl_m1 + sl_mode2 * I_Casl_m2; 
	//	cout << sl_mode2 << "\t" <<  I_Casl_m1 << "\t" <<  I_Casl_m2 << endl;
	
	// Na and K currents through LCC
	//	double I_CaKj2 = ibark * Fjunc_CaL * ( ( 1. - junc_mode2 ) * Po_LCCj_m1 + junc_mode2 * Po_LCCj_m2 ) * Q10CaL ^ Qpow * ICa_scale; 
	double I_CaKj2 = ibark * Fjunc_CaL * ( ( 1. - junc_mode2 ) * Po_LCCj_m1 + junc_mode2 * Po_LCCj_m2 ) * pow( Q10CaL , Qpow ) * ICa_scale; 
	//	double I_CaKsl2 = ibark * Fsl_CaL * ( ( 1. - sl_mode2 ) * Po_LCCsl_m1 + sl_mode2 * Po_LCCsl_m2) * Q10CaL ^ Qpow * ICa_scale;
	double I_CaKsl2 = ibark * Fsl_CaL * ( ( 1. - sl_mode2 ) * Po_LCCsl_m1 + sl_mode2 * Po_LCCsl_m2) * pow( Q10CaL , Qpow ) * ICa_scale;
	double I_CaK2 = I_CaKj2 + I_CaKsl2;
	//	double I_CaNa_junc2 = ( Fjunc_CaL * ibarna_j * ( ( 1. - junc_mode2 ) * Po_LCCj_m1 + junc_mode2 * Po_LCCj_m2 ) * Q10CaL ^ Qpow ) * ICa_scale;
	double I_CaNa_junc2 = ( Fjunc_CaL * ibarna_j * ( ( 1. - junc_mode2 ) * Po_LCCj_m1 + junc_mode2 * Po_LCCj_m2 ) * pow( Q10CaL , Qpow ) ) * ICa_scale;
	//	double I_CaNa_sl2 = Fsl_CaL * ibarna_sl * ( ( 1. - sl_mode2 ) * Po_LCCsl_m1 + sl_mode2 * Po_LCCsl_m2 ) * Q10CaL ^ Qpow * ICa_scale;
	double I_CaNa_sl2 = Fsl_CaL * ibarna_sl * ( ( 1. - sl_mode2 ) * Po_LCCsl_m1 + sl_mode2 * Po_LCCsl_m2 ) * pow( Q10CaL , Qpow ) * ICa_scale;
	
	// These are now able to switch depending on whether or not the flag to
	// switch to Markov model of ICa is ON
	double I_Ca_junc = ( 1. - ICa_MarkovFlag ) * I_Ca_junc1 + ICa_MarkovFlag*I_Ca_junc2;
	double I_Ca_sl = ( 1. - ICa_MarkovFlag ) * I_Ca_sl1 + ICa_MarkovFlag*I_Ca_sl2;
	double I_Ca = I_Ca_junc + I_Ca_sl;   // Total Ca curren throuhgh LCC
	double I_CaNa_junc = ( 1. - ICa_MarkovFlag ) * ( I_CaNa_junc1 ) + (ICa_MarkovFlag) * (I_CaNa_junc2);
	double I_CaNa_sl = ( 1. - ICa_MarkovFlag ) * ( I_CaNa_sl1 ) + (ICa_MarkovFlag) * (I_CaNa_sl2);
	double I_CaNa = I_CaNa_junc + I_CaNa_sl;   // Total Na current through LCC
	double I_CaK = ( 1. - ICa_MarkovFlag ) * ( I_CaK1 ) + ICa_MarkovFlag * (I_CaK2);  // Total K current through LCC
	
	// Collect all currents through LCC
	double I_Catot = I_Ca + I_CaK + I_CaNa;
	ydot[42]=-I_Ca*Cmem/(Vmyo*2*Frdy)*1e3;
	//	cout << ICa_MarkovFlag << "\t" << I_Ca_junc1 << "\t" << I_Ca_junc2 << "\t" << I_Ca_sl1 << "\t" << I_Ca_sl2 << endl;
	
	// global I_Ca_store ibar_store
	pars1->I_Ca_store = I_Catot;
	pars1->ibar_store = ibarca_j;
	// --- END LCC MARKOV MODEL --- //
	
	
	
	//// I_ncx: Na/Ca Exchanger flux
		const double Ka_junc = 1. / ( 1 + pow ( (Kdact / y[35]), 2) );
		const double Ka_sl = 1. / ( 1 + pow ( (Kdact /y[36]),2 ) );
		const double s1_junc = exp( nu * y[38] * FoRT ) * y[31] * y[31] * y[31] * Cao;
		const double s1_sl = exp( nu * y[38] * FoRT ) * y[32] * y[32] * y[32] * Cao;
		const double s2_junc = exp( (nu - 1) * y[38] * FoRT ) * Nao * Nao * Nao * y[35];
		
		const double s3_junc = KmCai * pow(Nao,3) * ( 1 + pow( ( y[31] / KmNai ), 3) ) + pow(KmNao,3) * y[35] * ( 1 + y[35] / KmCai ) + KmCao * y[31] * y[31] * y[31] + y[31] * y[31] * y[31] * Cao + Nao * Nao * Nao * y[35];
		const double s2_sl = exp( ( nu - 1 ) * y[38] * FoRT ) * Nao * Nao * Nao * y[36];
		const double s3_sl = KmCai * Nao * Nao * Nao * ( 1 + pow( (y[32] / KmNai), 3) ) + KmNao * KmNao * KmNao * y[36] * (1 + y[36]/KmCai) + KmCao * y[32] * y[32] * y[32] + y[32] * y[32] * y[32] * Cao + Nao * Nao * Nao * y[36];
		
		const double I_ncx_junc = Fjunc * IbarNCX * pow( Q10NCX , Qpow ) * Ka_junc * ( s1_junc - s2_junc ) / s3_junc / ( 1 + ksat * exp( (nu-1) * y[38] * FoRT ) );
		const double I_ncx_sl = Fsl * IbarNCX * pow( Q10NCX , Qpow ) * Ka_sl * ( s1_sl - s2_sl ) / s3_sl / ( 1 + ksat * exp( (nu-1) * y[38] * FoRT) );
		const double I_ncx = I_ncx_junc + I_ncx_sl;
	
		
		
	ydot[44] = 2. * I_ncx * Cmem / ( Vmyo * 2. * Frdy ) * 1e3;
	
	// global Incx
	pars1->Incx = I_ncx;
	//// I_pca: Sarcolemmal Ca Pump Current
	
	const double I_pca_junc = Fjunc * pow( Q10SLCaP , Qpow ) * IbarSLCaP * pow( y[35], 1.6 ) / ( pow( KmPCa, 1.6) + pow( y[35], 1.6 ) ) ;
	const double I_pca_sl = Fsl * pow( Q10SLCaP , Qpow ) * IbarSLCaP * pow( y[36], 1.6 ) / ( pow( KmPCa , 1.6 ) + pow( y[36] , 1.6 ) ) ;
	const double I_pca = I_pca_junc + I_pca_sl;
	
	ydot[43] = -I_pca * Cmem / ( Vmyo * 2 * Frdy ) * 1e3;
	
	//// I_cabk: Ca Background Current
	const double I_cabk_junc = Fjunc*GCaB*(y[38]-eca_junc);
	const double I_cabk_sl = Fsl*GCaB*(y[38]-eca_sl);
	const double I_cabk = I_cabk_junc+I_cabk_sl;
	
	ydot[45] = -I_cabk * Cmem / ( Vmyo * 2 * Frdy ) * 1e3;
	
	//// I_CFTR or I_cl_(cAMP) - Cystic Fibrosis Transmembrane Conductance Reg.
	// This is an Em- and time-independent current that is activated by PKA
	double fact_pka_cftr = 1.1933 * ICFTR_PKAp - 0.1933;
	double gCFTR = fact_pka_cftr * 4.9e-3;     // [A/F]  - Max value as in Shannon et al. (2005)
	double Icftr = gCFTR * ( y[38] - ecl );
	
	// global ICFTR
	pars1->ICFTR = Icftr;
	//// RyR model - SR release fluxes and leak
	const double MaxSR = 15.; 
	const double MinSR = 1.;
	const double kCaSR = MaxSR - ( MaxSR - MinSR ) / ( 1. + pow( ( ec50SR / y[30] ) , 2.5 ) ) ;
	double koSRCa = koCa / kCaSR;
	double kiSRCa = kiCa * kCaSR;
	double kleak = 5.348e-6;
	
	
	
	
	
	////// CaMKII and PKA-dependent phosphoregulation of RyR Po //////
	double fCKII_RyR = ( 20. * RyR_CKp / 3. - 1. / 3 );
	double fPKA_RyR = RyR_PKAp * 1.025 + 0.9750;
	koSRCa = ( fCKII_RyR + fPKA_RyR - 1 ) * koSRCa;
	
	// ODEs for RyR states and SR release through open RyRs
    
    kon = drug * 5500;
    
    double IC50_flecRyR = 2e-6; //test 5e-6, 11e-6, 17e-6
    double koff_flecRyR = IC50_flecRyR * 5500;
    
    double RI = 1. - y[13]-y[14]-y[15]-y[114];

    
	ydot[13] = ( kim * RI - kiSRCa * y[35] * y[13] ) - ( koSRCa * y[35] * y[35] * y[13] - kom * y[14]);   //  R
    ydot[14] = ((koSRCa*y[35]*y[35]*y[13]-kom*y[14])-(kiSRCa*y[35]*y[14]-kim*y[15])+koff_flecRyR*y[114]-kon*y[14]); //O

    ydot[15] = ((kiSRCa*y[35]*y[14]-kim*y[15])-(kom*y[15]-koSRCa*y[35]*y[35]*RI));   // I
    
    ydot[114] = ((kon*y[14]-koff_flecRyR*y[114]));

	
	
	const double J_SRCarel = ks* ( y[14]+0.2*y[114]  ) *(y[30]-y[35]);          //  [mM/ms]
	
	
	// Passive RyR leak - includes CaMKII regulation of leak flux
	kleak = ( 1. / 3 + 10. * RyR_CKp / 3. ) * kleak;
	double J_SRleak = kleak  * ( y[30] - y[35] );              //   [mmol/L cyt/ ms]
	
	// global Jleak 
	pars1->Jleak[0] = J_SRCarel * Vsr / Vmyo + J_SRleak;   // Total Jleak [mmol/L cyt/ms]
	pars1->Jleak[1] = J_SRleak;                        // Passive SR leak only [mmol/L cyt/ms]
	//// SERCA model - SR uptake fluxes
	// CaMKII and PKA-dependent phosphoregulation of PLB (changes to SERCA flux)
	double fCKII_PLB = ( 1. - .5 * PLB_CKp );
	const double fracPKA_PLBo = .9926;       // Derived quantity - ((PLBtot - PLBp(baseline))/PLBtot)
	double fPKA_PLB = ( PLB_PKAn / fracPKA_PLBo ) * 3. / 4 + 1. / 4;
	
	// Select smaller value (resulting in max reduction of Kmf)
	if ( fCKII_PLB < fPKA_PLB ) {
		Kmf = Kmf * fCKII_PLB;
	} else if ( fPKA_PLB < fCKII_PLB ) {
		Kmf = Kmf * fPKA_PLB;
		
	}
	

	
	double J_serca = ( pow( Q10SRCaP , Qpow ) 
					  * Vmax_SRCaP 
					  * ( pow( ( y[37] / Kmf ) , hillSRCaP ) 
						 - pow( ( y[30] / Kmr ) , hillSRCaP ) ) 
					  / ( 1. + pow( ( y[37] / Kmf ) , hillSRCaP ) + pow( ( y[30] / Kmr ) , hillSRCaP ) ) );
	
	// global Jserca  
	pars1->Jserca = J_serca;
	
	// //  Sodium and Calcium Buffering
	ydot[16] = kon_na*y[31]*(Bmax_Naj-y[16])-koff_na*y[16];        //  NaBj      [mM/ms]
	ydot[17] = kon_na*y[32]*(Bmax_Nasl-y[17])-koff_na*y[17];       //  NaBsl     [mM/ms]
	
	//  Cytosolic Ca Buffers
	
	ydot[18] = kon_tncl*y[37]*(Bmax_TnClow-y[18])-koff_tncl*y[18];            //  TnCL      [mM/ms]
	ydot[19] = kon_tnchca*y[37]*(Bmax_TnChigh-y[19]-y[20])-koff_tnchca*y[19]; //  TnCHc     [mM/ms]
	ydot[20] = kon_tnchmg*Mgi*(Bmax_TnChigh-y[19]-y[20])-koff_tnchmg*y[20];   //  TnCHm     [mM/ms]
	ydot[21] = 0;//*** commented b/c buffering done by CaM module kon_cam*y[37]*(Bmax_CaM-y[21])-koff_cam*y[21];                 //  CaM       [mM/ms]
	ydot[22] = kon_myoca*y[37]*(Bmax_myosin-y[22]-y[23])-koff_myoca*y[22];    //  Myosin_ca [mM/ms]
	ydot[23] = kon_myomg*Mgi*(Bmax_myosin-y[22]-y[23])-koff_myomg*y[23];      //  Myosin_mg [mM/ms]
	ydot[24] = kon_sr*y[37]*(Bmax_SR-y[24])-koff_sr*y[24];                    //  SRB       [mM/ms]
	
	
	// J_CaB_cytosol = sum(ydot(19:25));
	double J_CaB_cytosol = 0;
	for ( int it = 18; it <= 24; it ++ ) {
		J_CaB_cytosol += ydot[it];
	}
	//  Junctional and SL Ca Buffers
	ydot[25] = kon_sll*y[35]*(Bmax_SLlowj-y[25])-koff_sll*y[25];       //  SLLj      [mM/ms]
	ydot[26] = kon_sll*y[36]*(Bmax_SLlowsl-y[26])-koff_sll*y[26];      //  SLLsl     [mM/ms]
	ydot[27] = kon_slh*y[35]*(Bmax_SLhighj-y[27])-koff_slh*y[27];      //  SLHj      [mM/ms]
	ydot[28] = kon_slh*y[36]*(Bmax_SLhighsl-y[28])-koff_slh*y[28];     //  SLHsl     [mM/ms]
	const double J_CaB_junction = ydot[25] + ydot[27];
	const double J_CaB_sl = ydot[26] + ydot[28];
	
	
	// //  Ion concentrations
	//  SR Ca Concentrations
	ydot[29] = kon_csqn*y[30]*(Bmax_Csqn-y[29])-koff_csqn*y[29];       //  Csqn      [mM/ms]
	ydot[30] = J_serca-(J_SRleak*Vmyo/Vsr+J_SRCarel)-ydot[29];         //  Ca_sr     [mM/ms] // Ratio 3 leak current
	
	
	//  Sodium Concentrations
	const double I_Na_tot_junc = I_Na_junc+I_nabk_junc+3*I_ncx_junc+3*I_nak_junc+I_CaNa_junc;   /// [uA/uF]
	const double I_Na_tot_sl = I_Na_sl+I_nabk_sl+3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl;   //[uA/uF]
	
	ydot[31] = -I_Na_tot_junc*Cmem/(Vjunc*Frdy)+J_na_juncsl/Vjunc*(y[32]-y[31])-ydot[16];
	ydot[32] = -I_Na_tot_sl*Cmem/(Vsl*Frdy)+J_na_juncsl/Vsl*(y[31]-y[32]) + J_na_slmyo/Vsl*(y[33]-y[32])-ydot[17];
	
	ydot[33] = J_na_slmyo/Vmyo*(y[32]-y[33]);            // [mM/msec] 
	
	
	
	//  Potassium Concentration
	const double I_K_tot = I_to + I_kr + I_ks + I_ki - 2 * I_nak + I_CaK + I_kp;     //  [uA/uF]
	
	ydot[34] =0; //  -I_K_tot*Cmem/(Vmyo*Frdy);
	
	//  Calcium Concentrations
	const double I_Ca_tot_junc = I_Ca_junc+I_cabk_junc+I_pca_junc-2*I_ncx_junc;                   //  [uA/uF]
	const double I_Ca_tot_sl = I_Ca_sl+I_cabk_sl+I_pca_sl-2*I_ncx_sl;            //  [uA/uF]
	ydot[35] = ( -I_Ca_tot_junc*Cmem/(Vjunc*2*Frdy)+J_ca_juncsl/Vjunc*(y[36]-y[35]) 
				-J_CaB_junction+(J_SRCarel)*Vsr/Vjunc+J_SRleak*Vmyo/Vjunc );  //  Ca_j
	ydot[36] = ( -I_Ca_tot_sl*Cmem/(Vsl*2*Frdy)+J_ca_juncsl/Vsl*(y[35]-y[36]) 
				+ J_ca_slmyo/Vsl*(y[37]-y[36])-J_CaB_sl) ;   //  Ca_sl
	
	ydot[37] = -J_serca*Vsr/Vmyo-J_CaB_cytosol +J_ca_slmyo/Vmyo*(y[36]-y[37]);	
	double junc_sl = J_ca_juncsl / Vsl * ( y[35] - y[36] );
	double sl_junc = J_ca_juncsl / Vjunc * ( y[36] - y[35] );
	double sl_myo = J_ca_slmyo / Vsl * ( y[37] - y[36] );
	double myo_sl = J_ca_slmyo / Vmyo * ( y[36] - y[37] );
	
	
	
	//// Membrane Potential
	double I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;                 // [uA/uF]
	double I_Cl_tot = I_ClCa + I_Clbk + Icftr;                         // [uA/uF]
	double I_Ca_tot = I_Ca_tot_junc + I_Ca_tot_sl;
	double I_tot = I_Na_tot + I_Cl_tot + I_Ca_tot + I_K_tot;
	
    ydot[38] = -( I_tot );
	
}


