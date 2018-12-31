//
//  morotti_et_al_mouse_eccODEfile.h
//
//
//  Created by Mao-Tsuen on 6/9/15.
//  Modified by Pei-Chi Yang on 06/22/2015
//


void morotti_et_al_mouse_eccODEfile( double dt, double tt, double * y, double * p, double * ydot , pars_rec * pars1, double * allDrugs );

void morotti_et_al_mouse_eccODEfile( double dt, double tt, double * y, double * p, double * ydot , pars_rec * pars1, double * allDrugs ) {
    
    // function ydot = morotti_et_al_mouse_eccODEfile(t,y,p)
    // This file describes the EC coupling starting from the framework of the
    // Shannon-Bers model of rabbit ventricular EC coupling.
    // Reference: Shannon TR, Wang F, Puglisi J, Weber C & Bers DM. (2004).
    // A mathematical treatment of integrated Ca dynamics within the
    // ventricular myocyte. Biophysical Journal 87, 3351-3371.
    
    //// Assign passed-in parameters
    
    double cycleLength = p[0];
    // CaMKII phosphorylated targets (//)
    double LCC_CKp = p[1];
    double RyR_CKp = p[2];
    double PLB_CKp = p[3];
    // PKA phosphorylated targets (//)
    double LCCa_PKAp = p[4];
    double LCCb_PKAp = p[5];
    double PLB_PKAn = p[6]; // This is // non-phosphorylated PLB targets
    double RyR_PKAp = p[7];
    double TnI_PKAp = p[8];
    double IKs_PKAp = p[9]; // not present in MOUSE
    double ICFTR_PKAp = p[10]; // not present in MOUSE
    double IKur_PKAp = p[15]; // MOUSE
    double PLM_PKAp = p[16]; // MOUSE
    // Flag for CaMKII-OE
    double CKIIOE = p[11];
    // Protocols
    double rest = p[12];
    double variablePar = p[13];
    // PKA
    double Ligtot = p[14];
    
    const double drug = allDrugs[0];

    const double BlockGkr = allDrugs[1]; // ratio
    
    // ydot = zeros(size(y));
    //// Flags
    
    // Set CKIIflag to 1 for CaMKII-OE, 0 otherwise
    double  CKIIflag = CKIIOE;
    // Set ICa_MarkovFlag to 1 for Markov ICa, 0 otherwise
    double ICa_MarkovFlag = 1;
    // Set INa_MarkovFlag to 1 for Markov INa, 0 otherwise
    double  INa_MarkovFlag = 1;
    // Set Ito to use either origl params (=0) or Grandi Params (=1)
    double ItoFlag = 1;
    
    // Na clamp (set flag to 1, 0 otherwise)
    double NaClampFlag = 0;
    // PLM KO (set flag to 1, 0 otherwise)
    double PLMkoFlag = 0;
    // Strophanthidin (dyastolic Na influx) (set flag to 1, 0 otherwise)
    double StrophFlag = 0;
    // Caffeine-induced Ca transient (set flag to 1, 0 otherwise)
    double CaffeineFlag = 0;
    // Digitalis (set flag to 1, 0 otherwise)
    double DigitalisFlag = 0;
    //// Na loading and CaMKII-Na-Ca-CaMKII loop properties
    
    // Na loading parameters ON (set flag to 1, 0 otherwise)
    double NaGainFlag=0; // WT (default 0)
    if ( CKIIflag==1 ){
        NaGainFlag=1; // CaMKII-OE (default 1)
    } // end
    
    // CaMKII-Na-Ca-CaMKII loop closed (set flag to 1, 0 otherwise)
    double loop = 0; // (default 0)
    //// Model Parameters
    
    // Constants
    double R = 8314.;       // [J/kmol*K]
    double Frdy = 96485.;   // [C/mol]
    double Temp = 310.;     // [K] 310 K (37 C) for BT / 295 K (22 C) for RT
    double FoRT = Frdy / R / Temp;
    double Qpow = ( Temp - 310. ) / 10.;
    
    // Cell geometry
    double Acell = 20.e3; // [um^2] MOUSE
    //Cmem = 1.3810e-10; // [F] membrane capacitance RABBIT
    double Cmem = Acell * 1.e-14; // [F] 200 pF membrane capacitance MOUSE
    
    // Fractional currents in compartments
    //Fjunc = 0.11; Fsl = 1-Fjunc; // RABBIT
    //Fjunc = 0.3; Fsl = 1-Fjunc; // MOUSE
    double Fjunc = 17. / ( 17. + 31. ) * 7. / 17. + 31. / ( 17. + 31. ) * 2. / 31.;
    double Fsl = 1. - Fjunc; // MOUSE Fjunc = 0.1875;
    //Fjunc_nak = Fjunc; Fsl_nak = 1-Fjunc_nak; // RABBIT
    double Fjunc_nak = 1.6 * 17. / ( 1.6 * 17. + 31. ) * 7. / 17. + 31. / ( 1.6 * 17. + 31. ) * 2. / 31.;
    double Fsl_nak = 1. - Fjunc_nak; // MOUSE Fjunc = 0.2268;
    double Fjunc_ncx = Fjunc;
    double Fsl_ncx = 1. - Fjunc_ncx;
    //Fjunc_ncx = Fjunc_nak; Fsl_ncx = 1-Fjunc_ncx;
    //Fjunc_ncx = 0.5; Fsl_ncx = 1-Fjunc_ncx;
    //Fjunc_ncx = 0.9; Fsl_ncx = 1-Fjunc_ncx;
    double Fjunc_CaL = 0.9;
    double Fsl_CaL = 1. - Fjunc_CaL;
    
    double cellLength = 100.; // cell length [um]
    double cellRadius = 10.25; // cell radius [um]
    double junctionLength = 15.e-3; // junc length [um]
    double junctionRadius = 160.e-3; // junc radius [um]
    double distSLcyto = 0.45; // dist. SL to cytosol [um]
    //distJuncSL = 0.5; // dist. junc to SL [um] RABBIT
    double distJuncSL = 0.3; // dist. junc to SL [um] MOUSE
    double DcaJuncSL = 1.64e-6; // Dca junc to SL [cm^2/sec]
    double DcaSLcyto = 1.22e-6; // Dca SL to cyto [cm^2/sec]
    double DnaJuncSL = 1.09e-5; // Dna junc to SL [cm^2/sec]
    double DnaSLcyto = 1.79e-5; // Dna SL to cyto [cm^2/sec]
    double Vcell = pi * cellRadius * cellRadius * cellLength * 1.e-15; // [L]
    double Vmyo = 0.65 * Vcell;
    double Vsr = 0.035 * Vcell;
    double Vsl = 0.02 * Vcell;
    double Vjunc = 0.0539 * .01 * Vcell;
    //SAjunc = 20150*pi*2*junctionLength*junctionRadius; // [um^2] RABBIT
    //SAsl = pi*2*cellRadius*cellLength // [um^2] RABBIT
    double SAsl = Fsl * Acell; // [um^2]  MOUSE
    double Njunc = ( Fjunc * Acell ) / ( pi * junctionRadius * junctionRadius ); // [-]
    double SAjunc = Njunc * pi * 2. * junctionLength * junctionRadius; // [um^2] MOUSE
    //spacing=sqrt(2*Acell/Njunc); // [um] not used -> reduction in distJuncSL
    
    //J_ca_juncsl = 1/1.2134e12; // [L/msec] [m^2/sec] RABBIT
    //J_ca_slmyo = 1/2.68510e11; // RABBIT
    //J_na_juncsl = 1/(1.6382e12/3*100); // RABBIT
    //J_na_slmyo = 1/(1.8308e10/3*100);  // RABBIT
    double J_ca_juncsl = DcaJuncSL * SAjunc / distJuncSL * 1.e-10; // [L/msec] [m^2/sec] MOUSE
    double J_ca_slmyo = DcaSLcyto * SAsl / distSLcyto * 1.e-10; // MOUSE
    double J_na_juncsl = DnaJuncSL * SAjunc / distJuncSL * 1.e-10; // MOUSE
    double J_na_slmyo = DnaSLcyto * SAsl / distSLcyto * 1.e-10; // MOUSE
    
    // Fixed ion concentrations
    double Cli = 15.;   // Intracellular Cl  [mM]
    double Clo = 150.;  // Extracellular Cl  [mM]
    double Ko = 5.4;   // Extracellular K   [mM]
    double Nao = 140.;  // Extracellular Na  [mM]
    double Cao = 1.;    // Extracellular Ca  [mM] MOUSE // 1.8 mM in RABBIT
    double Mgi = 1.;    // Intracellular Mg  [mM]
    
    // Nernst Potentials
    double ena_junc = ( 1. / FoRT ) * log( Nao / y[31] );     // [mV]
    double ena_sl = ( 1. / FoRT ) * log( Nao / y[32] );       // [mV]
    double ek = ( 1. / FoRT ) * log( Ko / y[34] );	        // [mV]
    double eca_junc = ( 1. / FoRT / 2. ) * log( Cao / y[35] );   // [mV]
    double eca_sl = ( 1. / FoRT / 2. ) * log( Cao / y[36] );     // [mV]
    double ecl = ( 1. / FoRT ) * log( Cli / Clo );            // [mV]
    //// Na transport parameters
    
    double GNa = 16;// [mS/uF]
    
    double GNaB = 4.5 * 0.297e-3; // [mS/uF] changed from rabbit
    
    double IbarNaK = 5.; // [uA/uF] changed from rabbit (1.90719)
    if ( NaGainFlag == 1 ) {
        GNaB = GNaB * 4.;
        IbarNaK = IbarNaK * 0.9;
    } // end
    if ( DigitalisFlag == 1 ) {
        IbarNaK = IbarNaK * 0.5; // 50// block
    } // end
    double KmNaip = 19.; // [mM] changed from rabbit (11)
    double KmKo = 1.5; // [mM]
    double Q10NaK = 1.63;
    double Q10KmNai = 1.39;
    if ( PLMkoFlag==1 ) {
        PLM_PKAp = 1.;
        GNaB = GNaB * 48. / 20.;
        IbarNaK = IbarNaK * 0.8;
    } // end
    if ( StrophFlag==1 ) {
        IbarNaK = 0;
    } // end
    
    //    // INa Markov Model parameters
    //    double GNa2 = 10.64;//23;      // [mS/uF]
    //
    //    double P1a1 = 3.802;
    //    double P2a1 = 0.1027;
    //    double P3a1 = 2.5;
    //    double P4a1 = 17.;
    //    double P5a1 = 0.20;
    //    double P6a1 = 150.;
    //    double P4a2 = 15.;
    //    double P5a2 = 0.23;
    //    double P4a3 = 12.;
    //    double P5a3 = 0.25;
    //    double P1b1 = 0.1917;
    //    double P2b1 = 20.3;
    //    double P1b2 = 0.2;
    //    double P2b2 = 2.5;
    //    double P1b3 = 0.22;
    //    double P2b3 = 7.5;
    //    double P1a4 = 0.188495;
    //    double P2a4 = 16.6;
    //    double P3a4 = 0.393956;
    //    double P4a4 = 7.;
    //    double P1a5 = 7.e-7;
    //    double P2a5 = 7.2; // TG 7.8
    //    double P1b5 = 0.0044; // TG 0.00467
    //    double P2b5 = 2.e-5;
    //    double P1a6 = 100.;
    //    double P1b6 = 8.9554e-7; // TG 6.5449e-7
    //    double P2b6 = 11.3944;
    //    double P1a7 = 0.487e-4; // TG 0.3377e-3
    //    double P2a7 = 23.2696;
    //    double P1b7 = 0.2868e-3; // TG 1.868e-4
    //    double P2b7 = 35.9898;
    //    double P1a8 = 0.1e-7; // TG 6.5e-6
    //    double P1b8 = 9.8e-3; // TG 3.8e-3
    //    if ( CKIIflag == 1 ) { // MOUSE - CaMKII-OE
    //        P2a5 = 7.8;
    //        P1b5 = 0.00467;
    //        P1b6 = 6.5449e-7;
    //        P1a7 = 0.3377e-3;
    //        P1b7 = 1.868e-4;
    //        P1a8 = 6.5e-6;
    //        P1b8 = 3.8e-3;
    //    } // end
    //// K currents parameters
    
    double Kcoeff = 1.; // K current modulation
    
    double pNaK = 0.01833;
    double GtoSlow = 0.; // [mS/uF] changed from rabbit (0.06): NO ItoSlow in MOUSE
    double GtoFast = 0.44; // [mS/uF] changed from rabbit (0.02)
    if ( CKIIflag == 1 ) { // MOUSE
        GtoFast = GtoFast * 2. / 3.; // chronic CaMKII-OE effect
    } // end
    //Gkur = 0.3; // [mS/uF] only in MOUSE
    double Gkur1 = 1.1 * 0.16; // fast
    double Gkur2 = 0.14; // slow
    double Gss = 0.15; // [mS/uF] only in MOUSE
    
    //Flec block on Ikr
    double IC50 = 1.5 * (1e-6);
	double factor_flec = 1/(1+(drug/IC50));
    double gkr = 0.03 * pow( (Ko / 5.4) , 0.5 ) * factor_flec;
    
    double gkp = 0.001;
    
    // Cl current parameters
    double GClCa = 0.109625; // [mS/uF]
    double GClB = 9.e-3; // [mS/uF]
    double KdClCa = 100.e-3; // [mM]
    //// LTCC parameters
    
    double K_Ica = 1.65; // MOUSE
    double pNa = K_Ica * 1.5e-8; // [cm/sec]
    double pCa = K_Ica * 5.4e-4; // [cm/sec] - Ca permeability
    double pK = K_Ica * 2.7e-7; // [cm/sec]
    double KmCa = 0.6e-3; // [mM]
    double Q10CaL = 1.8;
    //// Ca transport parameters
    
    double IbarNCX = 1; // [uA/uF] changed from rabbit (9)
    if ( CKIIflag == 1 ) {
        IbarNCX = 1.5 * IbarNCX;
    } // end
    double KmCai = 3.59e-3;    // [mM]
    double KmCao = 1.3;        // [mM]
    double KmNai = 12.29;      // [mM]
    double KmNao = 87.5;       // [mM]
    double ksat = 0.27;        // [none]
    double nu = 0.35;          // [none]
    double Kdact = 1. / 2. * 0.256e-3; // [mM] changed from rabbit
    double Q10NCX = 1.57;      // [none]
    double IbarSLCaP = 0.0673; // [uA/uF]
    double KmPCa = 0.5e-3;     // [mM]
    double GCaB = 3. * 2.513e-4; // [uA/uF] changed from rabbit (2.513e-4)
    double Q10SLCaP = 2.35;    // [none]
    
    // SR flux parameters
    double Q10SRCaP = 2.6;          // [none]
    double Vmax_SRCaP = 1.15 * 1.15 * 2.86e-4; // [mM/msec] (mmol/L cytosol/msec) changed
    double Kmf = 0.3e-3; // [mM] changed from rabbit (0.246e-3) // from Yang-Saucerman
    double Kmr = 2.1; // [mM]L cytosol changed from rabbit (1.7) // from Yang-Saucerman
    double hillSRCaP = 1.787;       // [mM]
    double ks = 25.;                 // [1/ms]
    double koCa = 10.;               // [mM^-2 1/ms]
    double kom = 0.06;              // [1/ms]
    double kiCa = 0.5;              // [1/mM/ms]
    double kim = 0.005;             // [1/ms]
    double ec50SR = 0.45 + 0.05;      // [mM] changed from rabbit (0.45)
    
    if ( CaffeineFlag == 1 ) {
        koCa = koCa * 7.5;
        GCaB = 0;
        Vmax_SRCaP = 0;
    } // end
    //// Buffering parameters
    
    double Bmax_Naj = 7.561;       // [mM]
    double Bmax_Nasl = 1.65;       // [mM]
    double koff_na = 1.e-3;         // [1/ms]
    double kon_na = 0.1e-3;        // [1/mM/ms]
    double Bmax_TnClow = 70.e-3;    // [mM]                      // TnC low affinity
    double koff_tncl = 19.6e-3;    // [1/ms]
    double kon_tncl = 32.7;        // [1/mM/ms]
    double Bmax_TnChigh = 140.e-3;  // [mM]                      // TnC high affinity
    double koff_tnchca = 0.032e-3; // [1/ms]
    double kon_tnchca = 2.37;      // [1/mM/ms]
    double koff_tnchmg = 3.33e-3;  // [1/ms]
    double kon_tnchmg = 3.e-3;      // [1/mM/ms]
    double Bmax_myosin = 140.e-3;   // [mM]                      // Myosin buffering
    double koff_myoca = 0.46e-3;   // [1/ms]
    double kon_myoca = 13.8;       // [1/mM/ms]
    double koff_myomg = 0.057e-3;  // [1/ms]
    double kon_myomg = 0.0157;     // [1/mM/ms]
    double Bmax_SR = 19. * .9e-3;     // [mM]
    double koff_sr = 60.e-3;        // [1/ms]
    double kon_sr = 100.;           // [1/mM/ms]
    double Bmax_SLlowsl = 37.38e-3 * Vmyo / Vsl;        // [mM]    // SL buffering
    double Bmax_SLlowj = 4.62e-3 * Vmyo / Vjunc * 0.1;    // [mM]
    double koff_sll = 1300.e-3;     // [1/ms]
    double kon_sll = 100.;          // [1/mM/ms]
    double Bmax_SLhighsl = 13.35e-3 * Vmyo / Vsl;       // [mM]
    double Bmax_SLhighj = 1.65e-3 * Vmyo / Vjunc * 0.1;  // [mM]
    double koff_slh = 30.e-3;       // [1/ms]
    double kon_slh = 100.;          // [1/mM/ms]
    double Bmax_Csqn = 2.7;        //140e-3*Vmyo/Vsr; [mM]
    double koff_csqn = 65.;         // [1/ms]
    double kon_csqn = 100.;         // [1/mM/ms]
    
    // PKA-dependent phosphoregulation of TnI (increases Kd of TnC)
    double fracTnIpo = 0.062698;  // Derived quantity (TnI_PKAp(baseline)/TnItot)
    //fPKA_TnI = (1.45-0.45*(1-TnI_PKAp)/(1-fracTnIpo)); // Max effect +45//
    double fPKA_TnI = ( 1.61 - 0.61 * ( 1. - TnI_PKAp ) / ( 1. - fracTnIpo ) ); // Max effect +61//
    koff_tncl = koff_tncl * fPKA_TnI;
    //// Global Variable for Time
    
    // global tStep tArray
    // if t > tArray(tStep) // Roughly eliminates data from rejected time steps
    //    tStep = tStep + 1;
    // end
    // tArray(tStep) = t;
    
    //// I_Na: Fast Na Current
    
    // Max INa alterations with CaMKII hyperactivity as in Hund & Rudy 2008
    
    double inashift = 0;
    double alphaCKII = 0;
    double deltGbarNal_CKII = 0;
    
    if ( CKIIflag == 1 ) { // acute effects
        inashift = -3.25;
        alphaCKII = -.18;
        //deltGbarNal_CKII = 2;  // RABBIT
        if ( NaGainFlag == 1 ) {
            deltGbarNal_CKII = 3;  // MOUSE
        } else {
            deltGbarNal_CKII = 0;  // no Na Gain in OE
        } // end
        //    } else {
        //        inashift = 0;
        //        alphaCKII = 0;
        //        deltGbarNal_CKII = 0;
        //    } // end
    }
    double RyRp_WT_mean, RyRp_OE_mean, RyRp_OEloop_min, delta_loop, NaVsCaMKIIclamp;
    
    if ( loop == 1 ) {
        RyRp_WT_mean = 0.2101;
        RyRp_OE_mean = 0.7387; // Derived (1 Hz, no loop)
        RyRp_OEloop_min = 0.7033; // Derived (1 Hz, OE loop)
        delta_loop = ( ( 3. / ( RyRp_OE_mean - RyRp_WT_mean ) ) * RyR_CKp
                      - ( 3. / ( RyRp_OE_mean - RyRp_WT_mean ) ) * RyRp_WT_mean );
        NaVsCaMKIIclamp = 0; // if 1, CaMKII Clamp on NaV
        if ( NaVsCaMKIIclamp == 1 ) {
            delta_loop = ( ( 3. / ( RyRp_OE_mean - RyRp_WT_mean ) ) * RyRp_OEloop_min
                          - ( 3. / ( RyRp_OE_mean - RyRp_WT_mean ) ) * RyRp_WT_mean ) ;
        } // end
        GNaB = (4.5) * 0.297e-3 * ( 1. + delta_loop );
        if ( CKIIflag == 1 ) { // OE
            if ( NaGainFlag == 1 ) { // acute
                deltGbarNal_CKII = delta_loop;
            } else {
                deltGbarNal_CKII = 0;
            } // end
        } else { // WT
            deltGbarNal_CKII = 0;
        } // end
    } // end
    
    double am = 0.32 * ( y[38] + 47.13 ) / ( 1. - exp( -0.1 * ( y[38] + 47.13 ) ) );
    double bm = 0.08 * exp( -(y[38]) / 11. );
    
    double ah, aj, bh, bj;
    if ( ( y[38]-inashift ) >= -40. ) {
        ah = 0;
        aj = 0;
        //bh = 1/(0.13*(1+exp(-((y[38]-inashift)+10.66)/11.1))); // RABBIT
        bh = 0.66 * 1. / ( 0.13 * ( 1. + exp( -( ( y[38] - inashift ) + 10.66 ) / 11.1 ) ) ); // MOUSE
        bj = ( 0.3 * exp( -2.535e-7 * ( y[38] - inashift ) )
              / ( 1. + exp( -0.1 * ( ( y[38] - inashift ) + 32. ) ) ) );
    } else {
        ah = 0.135 * exp( ( 80. + ( y[38] - inashift ) ) / -6.8 );
        //bh = 3.56*exp(0.079*(y[38]-inashift))+3.1e5*exp(0.35*(y[38]-inashift)); // RABBIT
        bh = ( 1.1 * 3.56 * exp( 0.079 * ( y[38] - inashift - 2. ) )
              + 3.1e5 * exp( 0.35 * ( y[38] - inashift - 2. ) ) ); // MOUSE
        // Including alteration to aj as in Hund and Rudy 2008
        aj = ( ( 1. + alphaCKII )
              * ( ( -1.2714e5 * exp( 0.2444 * ( y[38] - inashift ) )
                   - 3.474e-5 * exp( -0.04391 * ( y[38] - inashift ) ) )
                 * ( ( y[38] - inashift ) + 37.78 )
                 / ( 1. + exp( 0.311 * ( ( y[38] - inashift ) + 79.23 ) ) ) ) );
        bj = ( 0.1212 * exp( -0.01052 * ( y[38] - inashift ) )
              / ( 1. + exp( -0.1378 * ( ( y[38] - inashift ) + 40.14 ) ) ) );
    } // end
    
    ydot[0] = 0;//am * ( 1. - y[0] ) - bm * y[0];
    ydot[1] = 0;//ah * ( 1. - y[1] ) - bh * y[1];
    ydot[2] = 0;//aj * ( 1. - y[2] ) - bj * y[2];
    
    double I_Na_junc1 = 0;//Fjunc * GNa * y[0]*y[0]*y[0] * y[1] * y[2] * ( y[38] - ena_junc );
    double I_Na_sl1 = 0;//Fsl * GNa * y[0]*y[0]*y[0] * y[1] * y[2] * ( y[38] - ena_sl );
    //// I_Na,L: Late INa current (as in Hund & Rudy 2008)
    
    //double GbarNal = .0065 * ( 1. + deltGbarNal_CKII ) * 2. ; // deltGbar assigned in 'Fast INa' section
    
    // h-gate (note: m-gate is same as INa m-gate -> using y[0] for this)
    //double hlss = 1. / ( 1. + exp( ( y[38] + 91. ) / 6.1 ) );
    //double tauhl = 600.; // ms
    ydot[46] = 0;//( hlss - y[46] ) / tauhl;
    double I_Nalj = 0;//Fjunc * GbarNal * y[0]*y[0]*y[0] * y[46] * ( y[38] - ena_junc );
    double I_Nalsl = 0;//Fsl * GbarNal * y[0]*y[0]*y[0] * y[46] * ( y[38] - ena_sl );
    //double I_Nal = I_Nalj+I_Nalsl;
    //// I_Na: alternative Markov Model - unused
    
    //    if ( INa_MarkovFlag == 1 ) {
    //        // State variables
    //        double CNa2 = y[47];
    //        double CNa1 = y[48];
    //        double ONa = y[49];
    //        double IFNa = y[50];
    //        double I1Na = y[51];
    //        double CNa3 = y[52];
    //        double ICNa2 = y[53];
    //        double ICNa3 = y[54];
    //        double LONa = y[55];
    //        double LCNa1 = y[56];
    //        double LCNa2 = y[57];
    //        double LCNa3 = y[58];
    //        double I2Na = ( 1. - ( ONa + CNa1 + CNa2 + CNa3 + IFNa + I1Na + ICNa2 + ICNa3 + LONa + LCNa1 + LCNa2 + LCNa3 ) );
    //        // Transition rates
    //        double alphaNa1 = P1a1 / ( P2a1 * exp( -( y[38] + P3a1 ) / P4a1 ) + P5a1 * exp( -( y[38] + P3a1 ) / P6a1 ) );
    //        double alphaNa2 = P1a1 / ( P2a1 * exp( -( y[38] + P3a1 ) / P4a2 ) + P5a2 * exp( -( y[38] + P3a1 ) / P6a1 ) );
    //        double alphaNa3 = P1a1 / ( P2a1 * exp( -( y[38] + P3a1 ) / P4a3 ) + P5a3 * exp( -( y[38] + P3a1 ) / P6a1 ) );
    //        double betaNa1 = P1b1 * exp( -( y[38] + P3a1 ) / P2b1 ); // shift
    //        double betaNa2 = P1b2 * exp( -( y[38] - P2b2 ) / P2b1 );
    //        double betaNa3 = P1b3 * exp( -( y[38] - P2b3 ) / P2b1 );
    //        double alphaNa4 = 1. / ( P1a4 * exp( -( y[38] + P4a4 ) / P2a4 ) + P3a4 );
    //        double alphaNa5 = P1a5 * exp( -( y[38] + P4a4 ) / P2a5 );
    //        double betaNa5 = ( P1b5 + P2b5 * ( y[38] + P4a4 ) );
    //        double betaNa6 = P1b6 * exp( -y[38] / P2b6 );
    //        double alphaNa7 = P1a7 * exp( y[38] / P2a7 );
    //        double betaNa7 = P1b7 * exp( -y[38] / P2b7 );
    //        double alphaNa8 = P1a8;
    //        double betaNa8 = P1b8;
    //        double betaNa4 = ( alphaNa3 * alphaNa4 * alphaNa5 ) / ( betaNa3 * betaNa5 );
    //        double alphaNa6 = alphaNa4 / P1a6;
    //        // ODEs
    //        double dCNa3  = betaNa8 * LCNa3 + betaNa1 * CNa2 + alphaNa5 * ICNa3 - ( alphaNa1 + betaNa5 + alphaNa8 ) * CNa3;
    //        double dCNa2  = betaNa8 * LCNa2 + alphaNa1 * CNa3 + betaNa2 * CNa1 + alphaNa5 * ICNa2 - ( betaNa1 + alphaNa2 + betaNa5 + alphaNa8 ) * CNa2;
    //        double dCNa1  = betaNa8 * LCNa1 + alphaNa2 * CNa2 + betaNa3 * ONa + alphaNa5 * IFNa - ( betaNa2 + alphaNa3 + betaNa5 + alphaNa8 ) * CNa1;
    //        double dONa   = betaNa8 * LONa + alphaNa3 * CNa1 + betaNa4 * IFNa - ( betaNa3 + alphaNa4 + alphaNa8 ) * ONa;
    //        double dIFNa  = alphaNa4 * ONa + betaNa5 * CNa1 + betaNa6 * I1Na + alphaNa2 * ICNa2 - ( betaNa4 + alphaNa5 + alphaNa6 + betaNa2 ) * IFNa;
    //        double dI1Na  = alphaNa6 * IFNa + betaNa7 * I2Na - ( betaNa6 + alphaNa7 ) * I1Na;
    //        double dICNa2 = alphaNa1 * ICNa3 + betaNa2 * IFNa + betaNa5 * CNa2 - ( betaNa1 + alphaNa2 + alphaNa5 ) * ICNa2;
    //        double dICNa3 = betaNa1 * ICNa2 + betaNa5 * CNa3 - ( alphaNa1 + alphaNa5 ) * ICNa3;
    //        double dLONa  = alphaNa3 * LCNa1 + alphaNa8 * ONa - ( betaNa8 + betaNa3 ) * LONa;
    //        double dLCNa1 = alphaNa8 * CNa1 + alphaNa2 * LCNa2 + betaNa3 * LONa - ( betaNa8 + betaNa2 + alphaNa3 ) * LCNa1;
    //        double dLCNa2 = betaNa2 * LCNa1 + alphaNa8 * CNa2 + alphaNa1 * LCNa3 - ( betaNa8 + betaNa1 + alphaNa2 ) * LCNa2;
    //        double dLCNa3 = alphaNa8 * CNa3 + betaNa1 * LCNa2 - ( betaNa8 + alphaNa1 ) * LCNa3;
    //        ydot[47] = dCNa2;
    //        ydot[48] = dCNa1;
    //        ydot[49] = dONa;
    //        ydot[50] = dIFNa;
    //        ydot[51] = dI1Na;
    //        ydot[52] = dCNa3;
    //        ydot[53] = dICNa2;
    //        ydot[54] = dICNa3;
    //        ydot[55] = dLONa;
    //        ydot[56] = dLCNa1;
    //        ydot[57] = dLCNa2;
    //        ydot[58] = dLCNa3;
    //    } else { // If not using INa Markov, set ODEs to zero to speed simulations
    ydot[47] = 0;
    ydot[48] = 0;
    ydot[49] = 0;
    ydot[50] = 0;
    ydot[51] = 0;
    ydot[52] = 0;
    ydot[53] = 0;
    ydot[54] = 0;
    ydot[55] = 0;
    ydot[56] = 0;
    ydot[57] = 0;
    ydot[58] = 0;
    //    } // end
    
    //    double I_Na_junc2 = Fjunc * GNa2 * ( y[49] + y[55] ) * ( y[38] - ena_junc ); // junc current
    //    double I_Na_sl2 = Fsl * GNa2 * ( y[49] + y[55] ) * ( y[38] - ena_sl ); // sl current
    //    //// I_Na: compute total current (fast and late components)
    //    double I_Na_junc;
    //    double I_Na_sl;
    //
    
    
    double sum;
	
	double O , OS , C1 , C2 , C3 , IC3 , IC2 , IF , IM1 , IM2 ;
	double DO , DOS , DC1 , DC2 , DC3 , DIC3 , DIC2 , DIF , DIM1 , DIM2 ;
	double D_O , D_OS , D_C1 , D_C2 , D_C3 , D_IC3 , D_IC2 , D_IF , D_IM1 , D_IM2 ;
	
	double BO, BC3, BC2, BC1;
	double DBO, DBC3, DBC2, DBC1;
	double D_BO, D_BC3, D_BC2, D_BC1;
	
	double O_n , OS_n , C1_n , C2_n , C3_n , IC3_n , IC2_n , IF_n , IM1_n , IM2_n ;
	double DO_n , DOS_n , DC1_n , DC2_n , DC3_n , DIC3_n , DIC2_n , DIF_n , DIM1_n , DIM2_n ;
	double D_O_n , D_OS_n , D_C1_n , D_C2_n , D_C3_n , D_IC3_n , D_IC2_n , D_IF_n , D_IM1_n , D_IM2_n ;
	
	double BO_n, BC3_n, BC2_n, BC1_n;
	double DBO_n, DBC3_n, DBC2_n, DBC1_n;
	double D_BO_n, D_BC3_n, D_BC2_n, D_BC1_n;
	
	double O_o , OS_o , C1_o , C2_o , C3_o , IC3_o , IC2_o , IF_o , IM1_o , IM2_o ;
	double DO_o , DOS_o , DC1_o , DC2_o , DC3_o , DIC3_o , DIC2_o , DIF_o , DIM1_o , DIM2_o ;
	double D_O_o , D_OS_o , D_C1_o , D_C2_o , D_C3_o , D_IC3_o , D_IC2_o , D_IF_o , D_IM1_o , D_IM2_o ;
	
	double BO_o, BC3_o, BC2_o, BC1_o;
	double DBO_o, DBC3_o, DBC2_o, DBC1_o;
	double D_BO_o, D_BC3_o, D_BC2_o, D_BC1_o;
	
	
	double a11, a12, a13, a2, a3, a4, a5;
	double b11, b12, b13, b2, b3, b4, b5;
	double ax, bx, ax1, bx1, ax2, bx2;
	double a13c, b13c, a13n, b13n;
	double a22, b22, a33, b33, a44, b44, a55, b55;
	double a_22, b_22, a_33, b_33, a_44, b_44, a_55, b_55;
	
	double mu1, mu2;
	
	double kd, kon, k_on, koff, k_off, kcon, kbon, kcbon, kc_on, kcoff, kc_off, ki_on, ki_off, kboff, kcboff;
	double kd_closed, kd_open, kd_openb;
	double hh;
	
	hh = dt;
	
	
	
	
	const double Q10=3;
	
	const double Tfactor = 1.0/(pow(Q10, (37.0-(Temp-273))/10.0));
	
	const double pH=7.4;
	const double pKa=7.2;
	const double portion = 1.0/(1+ pow(10, (pH-pKa)) );
	const double diffusion= 5500; // 500;
	
    
    
	const double drug_charged=drug*portion;
	const double drug_neutral=drug*(1-portion);
	const double dd= -0.7;
	
	//Rate Constants *************************************************************************************************************************************************************************
	//WT Fits Reduced Model (no IM1, IM2)
	a11= Tfactor*8.5539/(7.4392e-2*exp(-y[38]/17.0)+ 2.0373e-1*exp(-y[38]/150));
	a12= Tfactor*8.5539/(7.4392e-2*exp(-y[38]/15.0)+ 2.0373e-1*exp(-y[38]/150));
	a13= Tfactor*8.5539/(7.4392e-2*exp(-y[38]/12.0)+ 2.0373e-1*exp(-y[38]/150));
	b11= Tfactor*7.5215e-2*exp(-y[38]/20.3);
	b12= Tfactor*2.7574*exp(-(y[38]-5)/20.3);
	b13= Tfactor*4.7755e-1*exp(-(y[38]-10)/20.3) * 0.45; // MOUSE
	
	a3 = Tfactor*5.1458e-6*exp(-y[38]/8.2471);
	b3 = Tfactor*6.1205*exp((y[38])/13.542);
    
    
	a2 = Tfactor*(13.370*exp(y[38]/43.749)) * 0.45;  // MOUSE
    
    
	b2 = ((a13*a2*a3)/(b13*b3));
	
	a4 = 0*a2;
	b4 = 0*a3;
	a5 = 0*a2;
	b5 = 0*a3;
	
	
	
	mu1 = 2.0462e-7;
    mu2 = 8.9731e-4 * 1;
	
	
	ax = 3.4229e-2*a2 * 0.45;  // MOUSE
	bx = 1.7898e-2*a3;
    
	
	//************************************ Drug Rate Constants *******************************************************************
	//WT Flec  bursting states, and just then optimized the kb0
	
	ax1 = 5.7839e-05 * ax;
	bx1 =  1.6689e-08* bx;
	a13c = 3.6324e-03 *a13;
	a22 = 1.4847e+03 *a2;
	b33 =  1.7352e-06* b3;
	a33 = 6.7505e-05 * a3;
	a44 =  2.4135e+00* a2;
	b44 =  4.9001e-02* a3;
	a55 = b55 = 0;
	
	
	
	ax2 = 2.6126e-01 * ax;
	a13n = 2.6452e+00 * a13;
	a_22 =  4.2385e+01 * a2;
	b_33 = 2.1181e+00 * b3;
	a_44 =  1.0326e-03 * a2;
	b_44 = 2.1378e-02 * a3;
	
	a_55 = 0;
	b_55 = 0;
    
	
	
	
	const double kd0=11.2*(1e-6);
	kd_open=kd0*exp( (dd*y[38]*Frdy) /(R*Temp));
	
	// charged drug
	kon=drug_charged*diffusion;
	koff=kd_open*diffusion;
	kcoff = koff;
	kcon = kon;
	
	//bursting drug
	kbon = kon;
	kboff = 95.2165 *(1e-6)  *   exp( (dd*y[38]*Frdy) /(R*Temp))   *diffusion;	//This expression gives at 30uM, 70% block at -100mV Nagamoto, 2000
	kcbon = kbon;
	kcboff = kboff;
	
	
	if (drug ==0 || drug_charged ==0 ){b13c = 0;}
	else{b13c = (b13*kcon*koff*a13c)/(kon*kcoff*a13);}
	
	if (b13c ==0){b22 = 0;}
	else {b22=(a13c*a22*a33)/(b13c*b33);}
	
	// neutral drug
	k_on = drug_neutral*diffusion;
	k_off= 400*(1e-6)*diffusion;	//400
	ki_on=k_on/2;
	ki_off= 5.4*(1e-6)*diffusion;	//3.4
	kc_on=k_on/2;
	kc_off= 800*(1e-6)*diffusion;	//900
	
	if (drug ==0 || drug_neutral ==0 ){a_33 = 0;}
	else {a_33 = (ki_off*a3*kc_on*b_33)/(ki_on*kc_off*b3);}
	
	if (drug ==0 || drug_neutral ==0){b13n = 0;}
	else {b13n = (b13*kc_on*a13n*k_off)/(kc_off*a13*k_on);}
	
	if (b13n==0){b_22 =0;}
	else {b_22 = (a_33*a13n*a_22)/(b_33*b13n);}
	
	if (drug ==0 || drug_neutral ==0){bx2 = 0;}
	else{bx2 = (bx*k_on*ax2*ki_off)/(ax*ki_on*k_off);}
	
    
    IC3 = y[87];
    IC2 = y[88];
    IF = y[89];
    IM1 = y[90];
    IM2 = y[91];
    C3 = y[92];
    C2 = y[93];
    C1 = y[94];
    O = y[95];
    OS = y[96];
    BO = y[97];
    BC3 = y[98];
    BC2 = y[99];
    BC1 = y[100];
    
    DC3 = y[101];
    DC2 = y[102];
    DC1 = y[103];
    DO = y[104];
    DOS = y[105];
    DIC3 = y[106];
    DIC2 = y[107];
    DIF = y[108];
    DIM1 = y[109];
    DIM2 = y[110];
    
    DBO = y[111];
    DBC3 = y[112];
    DBC2 = y[113];
    DBC1 = y[114];
    
    D_C3 = y[115];
    D_C2 = y[116];
    D_C1 = y[117];
    D_O = y[118];
    D_OS = y[119];
    D_IC3 = y[120];
    D_IC2 = y[121];
    D_IF = y[122];
    D_IM1 = y[123];
    D_IM2 = y[124];
    
    D_BO = y[125];
    D_BC3 = y[126];
    D_BC2 = y[127];
    D_BC1 = y[128];
    
    
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
    
    double hh2 = hh/2;
    
    const double co_O = 1. / ( 1 + hh2 * ( b13 + a2  + mu1 + kon + k_on + ax ) );
    const double co_C1 = 1. / ( 1 + hh2 * ( b12 + b3 + a13  + mu1 + kcon + kc_on));
    const double co_C2 = 1. / ( 1 + hh2 * ( b11 + b3 + a12  + mu1 + kcon + kc_on));
    const double co_C3 = 1. / ( 1 + hh2 * ( b3 + a11  + mu1 +  kcon + kc_on));
    const double co_IC3 = 1. / ( 1 + hh2 * ( a11 + a3 + ki_on));
    const double co_IC2 = 1. / ( 1 + hh2 * ( b11 + a3 + a12 + ki_on));
    const double co_IF = 1. / ( 1 + hh2 * ( b12 + b2 + a3 + a4 + ki_on));
    const double co_IM1 = 1. / ( 1 + hh2 * ( b4 + a5));
    const double co_IM2 = 1. / ( 1 + hh2 * ( b5 ));
    const double co_OS = 1. / ( 1 + hh2 * ( bx + ki_on));
    
    const double co_BO = 1. / (1 + hh2 * (mu2 + b13 + kbon + k_on)  );				//changed from kon to kbon
    const double co_BC3 =1. / (1 + hh2 * (mu2 + a11 + kcbon + kc_on)  );			//changed from kcon to kcbon
    const double co_BC2 = 1. / (1 + hh2 * (mu2 + b11 + a12 + kcbon + kc_on)  );		//changed from kcon to kcbon
    const double co_BC1 = 1. / (1 + hh2 * (mu2 + b12 + a13 + kcbon + kc_on)  );		//changed from kcon to kcbon
    
    const double co_DO = 1. / ( 1 + hh2 * ( koff + b13c + a22 + ax1 + mu1 ));
    const double co_DC1 = 1. / ( 1 + hh2 * ( kcoff + b12 + b33 + a13c + mu1 ));
    const double co_DC2 = 1. / ( 1 + hh2 * ( kcoff + b11 + b33 + a12 + mu1));
    const double co_DC3 = 1. / ( 1 + hh2 * ( kcoff+ b33 + a11 + mu1 ));
    const double co_DOS = 1. / ( 1 + hh2 * ( bx1));
    const double co_DIC3 = 1. / ( 1 + hh2 * ( a11 + a33));
    const double co_DIC2 = 1. / ( 1 + hh2 * ( a33 + b11 + a12));
    const double co_DIF = 1. / ( 1 + hh2 * ( a33 + b12 + a44 + b22));
    const double co_DIM1 = 1. / ( 1 + hh2 * (  b44 + a55 ));
    const double co_DIM2 = 1. / ( 1 + hh2 * ( b55 ) );
    const double co_DBO = 1. / (1 + hh2 * (kboff + b13 + mu2)  );
    const double co_DBC3 =1. / (1 + hh2 * (kcboff + a11 + mu2)  );
    const double co_DBC2 = 1. / (1 + hh2 * (kcboff + b11 + a12 + mu2)  );
    const double co_DBC1 = 1. / (1 + hh2 * (kcboff + b12 + a13 + mu2)  );
    
    const double co_D_O = 1. / ( 1 + hh2 * ( k_off + b13n + a_22 + ax2 + mu1 ));
    const double co_D_C1 = 1. / ( 1 + hh2 * ( kc_off + b12 + b_33 + a13n + mu1 ));
    const double co_D_C2 = 1. / ( 1 + hh2 * ( kc_off + b11 + b_33 + a12 + mu1 ));
    const double co_D_C3 = 1. / ( 1 + hh2 * ( kc_off + b_33 + a11 + mu1 ));
    const double co_D_OS = 1. / ( 1 + hh2 * ( bx2 + ki_off));
    const double co_D_IC3 = 1. / ( 1 + hh2 * ( a_33 + a11 + ki_off));
    const double co_D_IC2 = 1. / ( 1 + hh2 * ( a_33 + b11 + a12 + ki_off));
    const double co_D_IF = 1. / ( 1 + hh2 * ( a_33 + a_44 + b_22 + b12 + ki_off));
    const double co_D_IM1 = 1. / ( 1 + hh2 * ( b_44 + a_55));
    const double co_D_IM2 = 1. / ( 1 + hh2 * ( b_55 ));
    
    const double co_D_BO = 1. / (1 + hh2 * (k_off + b13 + mu2)  );
    const double co_D_BC3 =1. / (1 + hh2 * (kc_off + a11 + mu2)  );
    const double co_D_BC2 = 1. / (1 + hh2 * (kc_off + b11 + a12+ mu2)  );
    const double co_D_BC1 = 1. / (1 + hh2 * (kc_off + b12 + a13 + mu2)  );
    
    
    //Drug Free States
    O_o = O + hh2 * ( a13 * C1 + b2 * IF + koff * DO + k_off * D_O + bx * OS + mu2 * BO  - O * coef_O );
    C1_o = C1 + hh2 * (a12 * C2 + a3 * IF + b13 * O  + kcoff * DC1 + kc_off * D_C1 + mu2 * BC1 - C1 * coef_C1 );
    C2_o = C2 + hh2 * (a11 * C3 + a3 * IC2 + b12 * C1  + kcoff * DC2 + kc_off * D_C2 + mu2 * BC2 - C2 * coef_C2 );
    C3_o = C3 + hh2 * (a3 * IC3 + b11 * C2  + kcoff * DC3 + kc_off * D_C3 + mu2 * BC3 - C3 * coef_C3 );
    IC3_o = IC3 + hh2 * (b3 * C3 + b11 * IC2 + ki_off * D_IC3 - IC3 * coef_IC3 );
    IC2_o = IC2 + hh2 * (a11 * IC3 + b3 * C2 + b12 * IF + ki_off * D_IC2 - IC2 * coef_IC2 );
    IF_o = IF + hh2 * (a12 * IC2 + b3 * C1 + b4 * IM1 + a2 * O + ki_off * D_IF - IF * coef_IF );
    IM1_o = IM1 + hh2 * (a4 * IF + b5 * IM2 - IM1 * coef_IM1 );
    IM2_o = IM2 + hh2 * (a5 * IM1 - IM2 * coef_IM2 );
    OS_o = OS + hh2 * (ax * O + ki_off * D_OS - OS * coef_OS );
    
    BO_o = BO + hh2 * (mu1 * O + a13 * BC1 + kboff* DBO + k_off * D_BO - BO * coef_BO );
    BC3_o = BC3 + hh2 * (mu1 * C3 + b11 * BC2 + kcboff * DBC3 + kc_off * D_BC3 - BC3 * coef_BC3 );
    BC2_o = BC2 + hh2 * (mu1 * C2 + a11 * BC3 + b12 * BC1 + kcboff * DBC2 + kc_off * D_BC2 - BC2 * coef_BC2 );
    BC1_o = BC1 + hh2 * (mu1 * C1 + a12 * BC2 + b13 * BO + kcboff * DBC1 + kc_off * D_BC1 - BC1 * coef_BC1 );
    
    //Charged Drug Bound States
    DO_o = DO + hh2 * (kon * O + a13c * DC1 + bx1 * DOS + b22 * DIF + mu2 * DBO - DO * coef_DO );
    DC1_o = DC1 + hh2 * (kcon * C1 + a12 * DC2 + a33 * DIF + b13c * DO  + mu2 * DBC1 - DC1 * coef_DC1 );
    DC2_o = DC2 + hh2 * (kcon * C2 + a11 * DC3 + a33 * DIC2 + b12 * DC1  + mu2 * DBC2 - DC2 * coef_DC2 );
    DC3_o = DC3 + hh2 * (kcon * C3 + b11 * DC2 + a33 * DIC3  + mu2 * DBC3 - DC3 * coef_DC3 );
    DOS_o = DOS + hh2 * (ax1 * DO -  DOS * coef_DOS );
    DIC3_o = DIC3 + hh2 * (b33 * DC3 + b11 * DIC2 - DIC3 * coef_DIC3 );
    DIC2_o = DIC2 + hh2 * (b33 * DC2 + a11 * DIC3 + b12 * DIF - DIC2 * coef_DIC2 );
    DIF_o = DIF + hh2 * (b33 * DC1 + a12 * DIC2 + b44 * DIM1 + a22 * DO - DIF * coef_DIF );
    DIM1_o = DIM1 + hh2 * (a44 * DIF + b55 * DIM2 - DIM1 * coef_DIM1 );
    DIM2_o = DIM2 + hh2 * (a55 * DIM1 - DIM2 * coef_DIM2 );
    
    DBO_o = DBO + hh2 * (kbon * BO + a13 * DBC1 + mu1 * DO - DBO * coef_DBO );						//kbon stayed the same here
    DBC3_o = DBC3 + hh2 * (kcbon * BC3 + b11 * DBC2 + mu1 * DC3 - DBC3 * coef_DBC3 );
    DBC2_o = DBC2 + hh2 * (kcbon * BC2 + a11 * DBC3 + b12 * DBC1 + mu1 * DC2 - DBC2* coef_DBC2 );
    DBC1_o = DBC1 + hh2 * (kcbon * BC1 + a12 * DBC2 + b13 * DBO + mu1 * DC1 - DBC1* coef_DBC1 );
    //Need to check these
    //Neutral Drug Bound States
    D_O_o = D_O + hh2 * (k_on * O + a13n * D_C1 + b_22 * D_IF + bx2 * D_OS + mu2 * D_BO - D_O * coef_D_O );
    D_C1_o = D_C1 + hh2 * (kc_on * C1 + a12 * D_C2 + a_33 * D_IF + b13n * D_O + mu2 * D_BC1  - D_C1 * coef_D_C1 );
    D_C2_o = D_C2 + hh2 * (kc_on * C2 + a11 * D_C3 + a_33 * D_IC2 + b12 * D_C1 + mu2 * D_BC2 - D_C2 * coef_D_C2 );
    D_C3_o = D_C3 + hh2 * (kc_on * C3 + a_33 * D_IC3 + b11 * D_C2 + mu2 * D_BC3 - D_C3 * coef_D_C3 );
    D_OS_o = D_OS + hh2 * (ax2 * D_O + ki_on * OS - D_OS * coef_D_OS );
    D_IC3_o = D_IC3 + hh2 * (b_33 * D_C3 + b11 * D_IC2 + ki_on * IC3 - D_IC3 * coef_D_IC3 );
    D_IC2_o = D_IC2 + hh2 * (b_33 * D_C2 + a11 * D_IC3 + b12 * D_IF + ki_on * IC2 - D_IC2 * coef_D_IC2 );
    D_IF_o = D_IF + hh2 * (b_33 * D_C1 + a12 * D_IC2 + b_44 * D_IM1 + a_22 * D_O + ki_on * IF - D_IF * coef_D_IF );
    D_IM1_o = D_IM1 + hh2 * (a_44 * D_IF + b_55 * D_IM2 - D_IM1 * coef_D_IM1 );
    D_IM2_o = D_IM2 + hh2 * (a_55 * D_IM1 - D_IM2 * coef_D_IM2 );
    D_BO_o = D_BO + hh2 * (k_on * BO + a13 * D_BC1 + mu1 * D_O - D_BO * coef_D_BO  );
    D_BC3_o = D_BC3 + hh2 * (kc_on * BC3 + b11 * D_BC2 + mu1 * D_C3 - D_BC3 * coef_D_BC3  );
    D_BC2_o = D_BC2 + hh2 * (kc_on * BC2 + a11 * D_BC3 + b12 * D_BC1 + mu1 * D_C2 - D_BC2 * coef_D_BC2 );
    D_BC1_o = D_BC1 + hh2 * (kc_on * BC1 + a12 * D_BC2 + b13 * D_BO + mu1 * D_C1 - D_BC1* coef_D_BC1 );
    
    int iter = 0;
    double err_sum = 1;
    while ( err_sum > 1E-100 && iter < 100 ) {
        
        
        //Drug Free States
        O_n = co_O * ( O_o + hh2 * ( a13 * C1 + b2 * IF + koff * DO + k_off * D_O + bx * OS + mu2 * BO) );
        C1_n = co_C1 * ( C1_o + hh2 * (a12 * C2 + a3 * IF + b13 * O  + kcoff * DC1 + kc_off * D_C1 + mu2 * BC1 ) );
        C2_n = co_C2 * ( C2_o + hh2 * (a11 * C3 + a3 * IC2 + b12 * C1  + kcoff * DC2 + kc_off * D_C2 + mu2 * BC2) );
        C3_n = co_C3 * ( C3_o + hh2 * (a3 * IC3 + b11 * C2  + kcoff * DC3 + kc_off * D_C3 + mu2 * BC3) );
        IC3_n = co_IC3 * ( IC3_o + hh2 * (b3 * C3 + b11 * IC2 + ki_off * D_IC3  ) );
        IC2_n = co_IC2 * ( IC2_o + hh2 * (a11 * IC3 + b3 * C2 + b12 * IF + ki_off * D_IC2  ) );
        IF_n = co_IF * ( IF_o + hh2 * (a12 * IC2 + b3 * C1 + b4 * IM1 + a2 * O + ki_off * D_IF  ) );
        IM1_n = co_IM1 * ( IM1_o + hh2 * (a4 * IF + b5 * IM2  ) );
        IM2_n = co_IM2 * ( IM2_o + hh2 * (a5 * IM1  ) );
        OS_n = co_OS * ( OS_o + hh2 * (ax * O + ki_off * D_OS  ) );
        
        BO_n = co_BO * (BO_o + hh2 * (mu1 * O + a13 * BC1 + kboff* DBO + k_off * D_BO)  );
        BC3_n = co_BC3 * (BC3_o + hh2 * (mu1 * C3 + b11 * BC2 + kcboff * DBC3 + kc_off * D_BC3) );
        BC2_n = co_BC2 * (BC2_o + hh2 * (mu1 * C2 + a11 * BC3 + b12 * BC1 + kcboff * DBC2 + kc_off * D_BC2) );
        BC1_n = co_BC1 * (BC1_o + hh2 * (mu1 * C1 + a12 * BC2 + b13 * BO + kcboff * DBC1 + kc_off * D_BC1) );
        
        //Charged Drug Bound States
        DO_n = co_DO * ( DO_o + hh2 * (kon * O + a13c * DC1 + bx1 * DOS + b22 * DIF  + mu2 * DBO ) );
        DC1_n = co_DC1 * ( DC1_o + hh2 * (kcon * C1 + a12 * DC2 + a33 * DIF + b13c * DO + mu2 * DBC1  ) );
        DC2_n = co_DC2 * ( DC2_o + hh2 * (kcon * C2 + a11 * DC3 + a33 * DIC2 + b12 * DC1 + mu2 * DBC2) );
        DC3_n = co_DC3 * ( DC3_o + hh2 * (kcon * C3 + b11 * DC2 + a33 * DIC3 +  mu2 * DBC3) );   //************TYPO
        DOS_n = co_DOS * ( DOS_o + hh2 * (ax1 * DO  ) );
        DIC3_n = co_DIC3 * ( DIC3_o + hh2 * (b33 * DC3 + b11 * DIC2 ) );
        DIC2_n = co_DIC2 * ( DIC2_o + hh2 * (b33 * DC2 + a11 * DIC3 + b12 * DIF  ) );
        DIF_n = co_DIF * ( DIF_o + hh2 * (b33 * DC1 + a12 * DIC2 + b44 * DIM1 + a22 * DO ) );
        DIM1_n = co_DIM1 * ( DIM1_o + hh2 * (a44 * DIF + b55 * DIM2  ) );
        DIM2_n = co_DIM2 * ( DIM2_o + hh2 * (a55 * DIM1 ) );
        
        DBO_n = co_DBO * (DBO_o + hh2 * (kbon * BO + a13 * DBC1 + mu1 * DO)  );								//kbon stayed the same here
        DBC3_n = co_DBC3 * (DBC3_o + hh2 * (kcbon * BC3 + b11 * DBC2 + mu1 * DC3) );
        DBC2_n = co_DBC2 * (DBC2_o + hh2 * (kcbon * BC2 + a11 * DBC3 + b12 * DBC1 + mu1 * DC2) );
        DBC1_n = co_DBC1 *(DBC1_o + hh2 * (kcbon * BC1 + a12 * DBC2 + b13 * DBO + mu1 * DC1) );
        //Need to check these
        //Neutral Drug Bound States
        D_O_n = co_D_O * ( D_O_o + hh2 * (k_on * O + a13n * D_C1 + b_22 * D_IF + bx2 * D_OS + mu2 * D_BO ) );
        D_C1_n = co_D_C1 * ( D_C1_o + hh2 * (kc_on * C1 + a12 * D_C2 + a_33 * D_IF + b13n * D_O + mu2 * D_BC1 ) );
        D_C2_n = co_D_C2 * ( D_C2_o + hh2 * (kc_on * C2 + a11 * D_C3 + a_33 * D_IC2 + b12 * D_C1 + mu2 * D_BC2 ) );
        D_C3_n = co_D_C3 * ( D_C3_o + hh2 * (kc_on * C3 + a_33 * D_IC3 + b11 * D_C2 + mu2 * D_BC3) );
        D_OS_n = co_D_OS * ( D_OS_o + hh2 * (ax2 * D_O + ki_on * OS ) );
        D_IC3_n = co_D_IC3 * ( D_IC3_o + hh2 * (b_33 * D_C3 + b11 * D_IC2 + ki_on * IC3 ) );
        D_IC2_n = co_D_IC2 * ( D_IC2_o + hh2 * (b_33 * D_C2 + a11 * D_IC3 + b12 * D_IF + ki_on * IC2  ) );
        D_IF_n = co_D_IF * ( D_IF_o + hh2 * (b_33 * D_C1 + a12 * D_IC2 + b_44 * D_IM1 + a_22 * D_O + ki_on * IF  ) );
        D_IM1_n = co_D_IM1 * ( D_IM1_o + hh2 * (a_44 * D_IF + b_55 * D_IM2  ) );
        D_IM2_n = co_D_IM2 * ( D_IM2_o + hh2 * (a_55 * D_IM1  ) );
        D_BO_n = co_D_BO * (D_BO_o + hh2 * (k_on * BO + a13 * D_BC1 + mu1 * D_O)  );
        D_BC3_n = co_D_BC3 * (D_BC3_o + hh2 * (kc_on * BC3 + b11 * D_BC2 + mu1 * D_C3)  );
        D_BC2_n = co_D_BC2 * (D_BC2_o + hh2 * (kc_on * BC2 + a11 * D_BC3 + b12 * D_BC1 + mu1 * D_C2) );
        D_BC1_n = co_D_BC1 * (D_BC1_o + hh2 * (kc_on * BC1 + a12 * D_BC2 + b13 * D_BO + mu1 * D_C1) );
        
        
        
        err_sum = fabs( IC3 - IC3_n ) + fabs( IC2 - IC2_n ) +  + fabs( IF - IF_n ) +  + fabs( IM1 - IM1_n ) +  + fabs( IM2 - IM2_n ) +  + fabs( C3 - C3_n ) +  + fabs( C2 - C2_n ) +  + fabs( C1 - C1_n ) +  + fabs( O - O_n ) +  + fabs( OS - OS_n ) +  + fabs( BO - BO_n) + + fabs(BC3 - BC3_n) + + fabs(BC2 - BC2_n) + + fabs(BC1 - BC1_n) + +
        fabs( DC3 - DC3_n ) +  + fabs( DC2 - DC2_n ) +  + fabs( DC1 - DC1_n ) +  + fabs( DO - DO_n ) +  + fabs( DOS - DOS_n ) +   + fabs( DIC3 - DIC3_n ) +  + fabs( DIC2 - DIC2_n ) +  + fabs( DIF - DIF_n ) +  + fabs( DIM1 - DIM1_n ) +  + fabs( DIM2 - DIM2_n ) +  + fabs( DBO - DBO_n) + + fabs(DBC3 - DBC3_n) + + fabs(DBC2 - DBC2_n) + + fabs(DBC1 - DBC1_n) + +
        fabs( D_C3 - D_C3_n ) +  + fabs( D_C2 - D_C2_n ) +  + fabs( D_C1 - D_C1_n ) +  + fabs( D_O - D_O_n ) +  + fabs( D_OS - D_OS_n ) +  + fabs( D_IC3 - D_IC3_n ) +  + fabs( D_IC2 - D_IC2_n ) +  + fabs( D_IF - D_IF_n ) +  + fabs( D_IM1 - D_IM1_n ) +  + fabs( D_IM2 - D_IM2_n ) + + fabs(D_BO - D_BO_n) + + fabs(D_BC3 - D_BC3_n) + + fabs(D_BC2 - D_BC2_n) + + fabs(D_BC1 - D_BC1_n);
        
        
        IC3 = IC3_n;
        IC2 = IC2_n;
        IF = IF_n;
        IM1 = IM1_n;
        IM2 = IM2_n;
        C3 = C3_n;
        C2 = C2_n;
        C1 = C1_n;
        O = O_n;
        OS = OS_n;
        BO = BO_n;
        BC3 = BC3_n;
        BC2 = BC2_n;
        BC1 = BC1_n;
        
        DC3 = DC3_n;
        DC2 = DC2_n;
        DC1 = DC1_n;
        DO = DO_n;
        DOS = DOS_n;
        DIC3 = DIC3_n;
        DIC2 = DIC2_n;
        DIF = DIF_n;
        DIM1 = DIM1_n;
        DIM2 = DIM2_n;
        DBO = DBO_n;
        DBC3 = DBC3_n;
        DBC2 = DBC2_n;
        DBC1 = DBC1_n;
        
        D_C3 = D_C3_n;
        D_C2 = D_C2_n;
        D_C1 = D_C1_n;
        D_O = D_O_n;
        D_OS = D_OS_n;
        D_IC3 = D_IC3_n;
        D_IC2 = D_IC2_n;
        D_IF = D_IF_n;
        D_IM1 = D_IM1_n;
        D_IM2 = D_IM2_n;
        D_BO = D_BO_n;
        D_BC3 = D_BC3_n;
        D_BC2 = D_BC2_n;
        D_BC1 = D_BC1_n;
        
        iter++;
    }
    
    
    
    y[87] = IC3_n;
    y[88] = IC2_n;
    y[89] = IF_n;
    y[90] = IM1_n;
    y[91] = IM2_n;
    y[92] = C3_n;
    y[93] = C2_n;
    y[94] = C1_n;
    y[95] = O_n;
    y[96] = OS_n;
    y[97] = BO_n;
    y[98] = BC3_n;
    y[99] = BC2_n;
    y[100] = BC1_n;
    
    y[101] = DC3_n;
    y[102] = DC2_n;
    y[103] = DC1_n;
    y[104] = DO_n;
    y[105] = DOS_n;
    y[106] = DIC3_n;
    y[107] = DIC2_n;
    y[108] = DIF_n;
    y[109] = DIM1_n;
    y[110] = DIM2_n;
    y[111] = DBO_n;
    y[112] = DBC3_n;
    y[113] = DBC2_n;
    y[114] = DBC1_n;
    
    y[115] = D_C3_n;
    y[116] = D_C2_n;
    y[117] = D_C1_n;
    y[118] = D_O_n;
    y[119] = D_OS_n;
    y[120] = D_IC3_n;
    y[121] = D_IC2_n;
    y[122] = D_IF_n;
    y[123] = D_IM1_n;
    y[124] = D_IM2_n;
    y[125] = D_BO_n;
    y[126] = D_BC3_n;
    y[127] = D_BC2_n;
    y[128] = D_BC1_n;
    
    
    
    ydot[87] = 0;
    ydot[88] = 0;
    ydot[89] = 0;
    ydot[90] = 0;
    ydot[91] = 0;
    ydot[92] = 0;
    ydot[93] = 0;
    ydot[94] = 0;
    ydot[95] = 0;
    ydot[96] = 0;
    ydot[97] = 0;
    ydot[98] = 0;
    ydot[99] = 0;
    ydot[100] = 0;
    
    ydot[101] = 0;
    ydot[102] = 0;
    ydot[103] = 0;
    ydot[104] = 0;
    ydot[105] = 0;
    ydot[106] = 0;
    ydot[107] = 0;
    ydot[108] = 0;
    ydot[109] = 0;
    ydot[110] = 0;
    ydot[111] = 0;
    ydot[112] = 0;
    ydot[113] = 0;
    ydot[114] = 0;
    
    ydot[115] = 0;
    ydot[116] = 0;
    ydot[117] = 0;
    ydot[118] = 0;
    ydot[119] = 0;
    ydot[120] = 0;
    ydot[121] = 0;
    ydot[122] = 0;
    ydot[123] = 0;
    ydot[124] = 0;
    ydot[125] = 0;
    ydot[126] = 0;
    ydot[127] = 0;
    ydot[128] = 0;
    
    //**************************************************************************************************************************************************************************************************/
    
    
    sum =  y[87] + y[88] + y[89] + y[90] + y[91] + y[92] + y[93] + y[94] + y[95] +
    y[96] + y[97]  + y[98] + y[99] + y[100] + y[101] + y[102] + y[103] + y[104] + y[105] +
    y[106] + y[107]  + y[108] + y[109] + y[110] + y[111]  + y[112] + y[113] +y[114] + y[115] +
    y[116] + y[117] + y[118] +y[119] + y[120] + y[121] + y[122] + y[123] + y[124] + y[125] +
    y[126] + y[127] + y[128];
    
    // cout << "sum = " << Cell_ptr->sum << endl;
    if (fabs (sum -1.0) > 0.0001 || sum != sum ) {cout << "Error in WT Sum " << sum << endl;}
    
	
	
	double I_Na_junc2 = Fjunc*GNa*(y[95]+y[97])*(y[38]-ena_junc);
	double I_Na_sl2 = Fsl*GNa*(y[95]+y[97])*(y[38]-ena_sl);
    
    double I_Na_junc, I_Na_sl;
    
    if ( INa_MarkovFlag == 1 ) {
        I_Na_junc = I_Na_junc2;
        I_Na_sl = I_Na_sl2;
    } else {
        I_Na_junc = I_Na_junc1 + I_Nalj;
        I_Na_sl = I_Na_sl1 + I_Nalsl;
    } // end
    
    double I_Na = I_Na_junc + I_Na_sl;
    // global I_Na_store
    pars1->I_Na_store = I_Na;
    // I_Na_store(tStep) = I_Na;
    //// I_nabk: Na Background Current
    
    double I_nabk_junc = Fjunc * GNaB * ( y[38] - ena_junc );
    double I_nabk_sl = Fsl * GNaB * ( y[38] - ena_sl );
    double I_nabk = I_nabk_junc + I_nabk_sl;
    
    // global I_Nabk_store
    pars1->I_Nabk_store = I_nabk;
    // I_Nabk_store(tStep) = I_nabk;
    //// I_nak: Na/K Pump Current
    
    double sigma = ( exp( Nao / 67.3 ) - 1. ) / 7.;
    double fnak = 1. / ( 1. + 0.1245 * exp( -0.1 * y[38] * FoRT ) + 0.0365 * sigma * exp( -y[38] * FoRT ) );
    //KmNaip = 14; // PKA effect - 1 mM (Despa 2005)
    double fracPKA_PLMo = 0.116738; // Derived quantity (PLM_PKAp(baseline)/PLMtot)
    double fracPKA_PLMiso = 0.859251; // Derived quantity (PLM_PKAp(ISO)/PLMtot)
    //kPKA_PLM=(KmNaip-14)/(fracPKA_PLMiso/fracPKA_PLMo-1); // PLM_PKAp ISO
    double kPKA_PLM = KmNaip * ( 1. - 0.7019 ) / ( fracPKA_PLMiso / fracPKA_PLMo - 1. ); // PLM_PKAp ISO
    double KmNaip_PKA = -kPKA_PLM + kPKA_PLM * ( PLM_PKAp / fracPKA_PLMo );
    KmNaip = KmNaip - KmNaip_PKA;
    
    double I_nak_junc = ( Fjunc_nak * IbarNaK * fnak * Ko
                         / ( 1. + ( KmNaip*KmNaip*KmNaip*KmNaip ) / ( y[31]*y[31]*y[31]*y[31] ) )
                         / ( Ko + KmKo ) );
    double I_nak_sl = ( Fsl_nak * IbarNaK * fnak * Ko
                       / ( 1. + ( KmNaip*KmNaip*KmNaip*KmNaip ) / ( y[32]*y[32]*y[32]*y[32] ) )
                       / ( Ko + KmKo ) );
    double I_nak = I_nak_junc + I_nak_sl;
    
    // global I_NaK_store
    pars1->I_NaK_store = I_nak;
    // I_NaK_store(tStep) = I_nak;
    //// I_kur - IK,slow
    
    double xurss = 1. / ( 1. + exp( -( y[38] + 15. ) / 14. ) );
    double yurss = 1. / ( 1. + exp( ( y[38] + 48 ) / 6.2 ) );
    
    double tauxur = 0.95 + 0.05 * exp( -0.08 * y[38] );
    
    ydot[83] = ( xurss - y[83] ) / tauxur;       // IKslow1
    double tauxur2 = ( 1. + 7. / ( 1. + exp( -( y[38] + 45. ) / 8. ) )
                      + 20. * exp( -( ( y[38] + 35. ) / 10. ) * ( ( y[38] + 35. ) / 10. ) ) );//8+20*exp(-((y[38]+35)/10)^2);
    ydot[12] = ( xurss - y[12] ) / tauxur2;      // IKslow2
    //tauyur = 400+900*exp(-((y[38]+55)/16)^2)+150/(1+exp(-(y[38]+60)/8));
    //ydot[84] = (yurss-y[84])/tauyur;
    double tauyur1 = ( 400. + 900. * exp( -( ( y[38] + 55. ) / 16. ) * ( ( y[38] + 55. ) / 16. ) )
                      - 250. / ( 1. + exp( -( y[38] + 60. ) / 8. ) ) ); // fast
    ydot[84] = ( yurss - y[84] ) / tauyur1;
    double tauyur2 = ( 400. + 900. * exp( -( ( y[38] + 55. ) / 16. ) * ( ( y[38] + 55. ) / 16. ) )
                      + 550. / ( 1. + exp( -( y[38] + 60 ) / 8. ) ) ); // slow
    ydot[86] = ( yurss - y[86] ) / tauyur2;
    
    // PKA-dependent phosphoregulation of Ik,slow1 (increases Gkur1)
    double fracIKurp0 = 0.437635;  // Derived quantity (IKur_PKAp(baseline)/IKurtot)
    double fracIKurpISO = 0.718207; // Derived quantity (IKur_PKAp(ISO)/IKurtot)
    double a_Kur = ( 1.20 - 1. ) / ( fracIKurpISO / fracIKurp0 - 1. );
    double fracIKuravail = ( 1. - a_Kur ) + a_Kur * ( IKur_PKAp / fracIKurp0 ); // +20// with 0.1 uM ISO
    
    //I_kur = Kcoeff*fracIKuravail*Gkur*y[83]*y[84]*(y[38]-ek);
    double I_kur1 = Kcoeff * fracIKuravail * Gkur1 * y[83] * y[84] * ( y[38] - ek ); // IKslow1
    //I_kur2 = Kcoeff*fracIKuravail*Gkur2*y[83]*y[86]*(y[38]-ek); // IKslow2
    double I_kur2 = Kcoeff * Gkur2 * y[83] * y[86] * ( y[38] - ek ); // IKslow2 // no PKA effect
    
    double I_kur = I_kur1 + I_kur2;
    
    // global I_kur1_store I_kur2_store
    pars1->I_kur1_store = I_kur1;
    pars1->I_kur2_store = I_kur2;
    // I_kur1_store(tStep) = I_kur1;
    // I_kur2_store(tStep) = I_kur2;
    //// I_ss
    
    double xssss = xurss; // = 1/(1+exp(-(y[38]+15)/14));
    //tauxss = 14+0.8*exp(-0.08*y[38]);
    double tauxss = 70. * exp( -( ( y[38] + 43. ) / 30. ) * ( ( y[38] + 43. ) / 30. ) ) + 14.; // Iss
    ydot[85] = ( xssss - y[85] ) / tauxss;
    double I_ss = Kcoeff * Gss * y[85] * ( y[38] - ek ); //store
    
    // global I_ss_store
    pars1->I_ss_store = I_ss;
    // I_ss_store(tStep) = I_ss;
    //// I_kr: Rapidly Activating K Current
    
    double xrss = 1. / ( 1. + exp( -( y[38] + 50. ) / 7.5 ) );
    double tauxr = ( 1. / ( 1.38e-3 * ( y[38] + 7. ) / ( 1. - exp( -0.123 * ( y[38] + 7. ) ) )
                           + 6.1e-4 * ( y[38] + 10. ) / ( exp( 0.145 * ( y[38] + 10. ) ) - 1. ) ) );
    ydot[11] = (xrss-y[11])/tauxr;
    double rkr = 1/(1+exp((y[38]+33)/22.4));
    double I_kr = Kcoeff*gkr*y[11]*rkr*(y[38]-ek) * BlockGkr;
    
    // global I_kr_store
    pars1->I_kr_store = I_kr;
    // I_kr_store(tStep) = I_kr;
    //// I_ks: Slowly Activating K Current
    
    // Phosphoregulation of IKs by PKA parameters
    // fracIKspo = 0.07344;  // Derived quantity (IKs_PKAp(baseline)/IKstot)
    // fracIKsavail = (0.2*(IKs_PKAp/fracIKspo)+0.8);
    // Xs05 = 1.5*(2.0 - IKs_PKAp/fracIKspo);
    //
    // pcaks_junc = -log10(y[35])+3.0;
    // pcaks_sl = -log10(y[36])+3.0;
    // gks_junc = fracIKsavail*0.07*(0.057 +0.19/(1+ exp((-7.2+pcaks_junc)/0.6))); // Now regulated by PKA
    // gks_sl = fracIKsavail*0.07*(0.057 +0.19/(1+ exp((-7.2+pcaks_sl)/0.6)));     // Now regulated by PKA
    // eks = (1/FoRT)*log((Ko+pNaK*Nao)/(y[34]+pNaK*y[33]));
    // // xsss = 1/(1+exp(-(y[38]-1.5)/16.7)); // Original version
    // xsss = 1/(1+exp(-(y[38]-Xs05)/16.7));   // Now regulated by PKA
    // tauxs = 1/(7.19e-5*(y[38]+30)/(1-exp(-0.148*(y[38]+30)))+1.31e-4*(y[38]+30)/(exp(0.0687*(y[38]+30))-1));
    //ydot[12] = 0;//(xsss-y[12])/tauxs;
    // state variable y[12] is now used for IKslow2 activation
    double I_ks_junc = 0.;//Fjunc*gks_junc*y[12]^2*(y[38]-eks); // No IKs in mouse
    double I_ks_sl = 0.;//Fsl*gks_sl*y[12]^2*(y[38]-eks); // No IKs in mouse
    double I_ks = I_ks_junc + I_ks_sl;
    
    // global IKs_store
    pars1->IKs_store = I_ks;
    // IKs_store(tStep) = I_ks;
    //// I_kp: Plateau K current
    
    double kp_kp = 1. / ( 1. + exp( 7.488 - y[38] / 5.98 ) );
    double I_kp_junc = Kcoeff * Fjunc * gkp * kp_kp * ( y[38] - ek );
    double I_kp_sl = Kcoeff * Fsl * gkp * kp_kp * ( y[38] - ek );
    double I_kp = I_kp_junc + I_kp_sl;
    //// I_to: Transient Outward K Current (slow and fast components)
    
    // Itos (ABSENT IN MOUSE)
    double xtoss = 1. / ( 1. + exp( -( y[38] + 3.0 ) / 13. ) );
    double ytoss = 1. / ( 1. + exp( ( y[38] + 48. ) / 5. ) );
    double rtoss = 1. / ( 1. + exp( ( y[38] + 33.5 ) /10. ) ); // Rto not used in MOUSE model
    
    double tauxtos = 0.08 + 0.7 * exp( -( ( y[38] + 25. ) / 30. ) * ( ( y[38] + 25. ) / 30. ) );
    
    double taurtos, tauytos, Py, Pr1, Pr2;
    if ( ItoFlag == 0 ) {// (not used)
        // Shannon Versions
        //tauytos = 3e3/(1+exp((y[38]+60.0)/10))+30;
        taurtos = 2.8e3 / ( 1. + exp( ( y[38] + 60.0 ) / 10. ) ) + 220.; // no Rto
        tauytos = 100. + 400. / ( 1. + exp( ( y[38] + 25. ) / 5. ) );
    } else if ( ItoFlag == 1 && CKIIflag == 0 ) { // WT
        // Grandi Versions
        Py = 182.;
        Pr1 = 8085.;
        Pr2 = 313.;            // Normal
        //tauytos = Py/(1+exp((y[38]+33.5)/10))+1;
        taurtos = Pr1 / ( 1. + exp( ( y[38] + 33.5 ) / 10. ) ) + Pr2; // no Rto
        tauytos = 100. + 400. / ( 1. + exp( ( y[38] + 25. ) / 5. ) );
    } else if ( ItoFlag == 1 && CKIIflag == 1 ) {            // CaMKII-OE acute effect
        Py = 15.;
        Pr1 = 3600.;
        Pr2 = 500.;
        //GtoSlow = GtoSlow*1.5; // Rabbit
        //tauytos = Py/(1+exp((y[38]+33.5)/10))+1;
        taurtos = Pr1 / ( 1. + exp( ( y[38] + 33.5 ) / 10. ) ) + Pr2; // no Rto
        //tauytos = 35+400/(1+exp((y[38]+25)/5));    // MOUSE tau rec -73//
        tauytos = 100. + 35. / ( 1. + exp( ( y[38] + 25. ) / 5. ) );     // MOUSE tau rec -73//
    } // end
    
    ydot[7] = ( xtoss - y[7] ) / tauxtos;
    ydot[8] = ( ytoss - y[8] ) / tauytos;
    ydot[39] = ( rtoss - y[39] ) / taurtos;                // no Rto
    //I_tos = GtoSlow*y[7]*(y[8]+0.5*y[39])*(y[38]-ek); // [uA/uF]
    double I_tos = 0. * Kcoeff * GtoSlow * y[7] * y[8] * ( y[38] - ek ); // N0 Itos in MOUSE // [uA/uF]
    
    // Itof
    double xtofs = 1. / ( 1. + exp( -( y[38] + 3.0 ) / 13. ) ); // = xtoss
    double ytofs = 1. / ( 1. + exp( ( y[38] + 48. ) / 5. ) ); // = ytoss
    
    //tauxtof = 3.5*exp(-y[38]*y[38]/30/30)+1.5; // Original
    //tauxtof = 3.5*exp(-((y[38]+3)/30)^2)+1.5; // Version in Grandi Code (does not change AP shape)
    double tauxtof = ( 0.08 + 0.7 * exp( -( ( y[38] + 25. ) / 30. )
                                        * ( ( y[38] + 25. ) / 30. ) ) ); // = tauxtos;
    
    double tauytof = ( 10. + 32. * exp( -( ( y[38] + 55. ) / 16. )
                                       * ( ( y[38] + 55. ) / 16. ) )
                      + 8. / ( 1. + exp( -( y[38] + 60. ) / 8. ) ) );
    
    if ( CKIIflag == 1 ) { // MOUSE (CaMKII-OE acute effect) // tau rec -38//
        tauytof = ( 5. + 32. * exp( -( ( y[38] + 55. ) / 16. )
                                   * ( ( y[38] + 55. ) / 16. ) )
                   + 12. / ( 1. + exp( -( y[38] + 60. ) / 8. ) ) );
    } // end
    
    ydot[9] = ( xtofs - y[9] ) / tauxtof;
    ydot[10] = ( ytofs - y[10] ) / tauytof;
    double I_tof = Kcoeff * GtoFast * y[9] * y[10] * ( y[38] - ek );
    
    double I_to = I_tos + I_tof;
    
    pars1->I_to_store = I_to;     // Total I_to
    pars1->I_tof_store = I_tof;    // Fast Component
    pars1->I_tos_store = I_tos;    // Slow component
    //        global I_to_store
    //        I_to_store(1,tStep) = I_to;     // Total I_to
    //        I_to_store(2,tStep) = I_tof;    // Fast Component
    //        I_to_store(3,tStep) = I_tos;    // Slow component
    //// I_k1: Time-Independent K Current (I_ki)
    
    double aki = 1.02 / ( 1. + exp( 0.2385 * ( y[38] - ek - 59.215 ) ) );
    double bki = ( ( 0.49124 * exp( 0.08032 * ( y[38] + 5.476 - ek ) )
                    + exp( 0.06175 * ( y[38] - ek - 594.31 ) ) )
                  / ( 1. + exp( -0.5143 * ( y[38] - ek + 4.753 ) ) ) );
    double kiss = aki / ( aki + bki );
    //I_ki = 0.9*sqrt(Ko/5.4)*kiss*(y[38]-ek); // RABBIT
    double I_ki;
    if ( CKIIflag == 1 ){ // MOUSE (chronic CaMKII-OE effect)
        I_ki = 1. / 2 * 0.3 * pow( ( Ko / 5.4 ), 0.5 ) * kiss * ( y[38] - ek ) * Kcoeff;
    } else {
        I_ki = 0.3 * pow( ( Ko / 5.4 ), 0.5 ) * kiss * ( y[38] - ek ) * Kcoeff;
    } // end
    
    // global I_K1_store
    pars1->I_K1_store = I_ki;
    // I_K1_store(tStep) = I_ki;
    //// I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current
    
    double I_ClCa_junc = Fjunc * GClCa / ( 1. + KdClCa / y[35] ) * ( y[38] - ecl );
    double I_ClCa_sl = Fsl * GClCa / ( 1. + KdClCa / y[36] ) * ( y[38] - ecl );
    double I_ClCa = I_ClCa_junc + I_ClCa_sl;
    double I_Clbk = GClB * ( y[38] - ecl );
    //// Original H-H formulation for LCC - unused if ICa_MarkovFlag = 1
    
    double dss = 1. / ( 1. + exp( -( y[38] + 14.5 ) / 6.0 ) );
    double taud = ( dss * ( 1. - exp( -( y[38] + 14.5 ) / 6.0 ) )
                   / ( 0.035 * ( y[38] + 14.5 ) ) );
    double fss = ( 1. / ( 1. + exp( ( y[38] + 35.06 ) / 3.6 ) )
                  + 0.6 / ( 1. + exp( ( 50. - y[38] ) / 20. ) ) );
    double tauf = ( 1. / ( 0.0197 * exp( -( 0.0337 * ( y[38] + 14.5 ) )
                                        * ( 0.0337 * ( y[38] + 14.5 ) ) ) + 0.02 ) );
    
    // ydot[3] = (dss-y[3])/taud;
    // ydot[4] = (fss-y[4])/(tauf);
    // ydot[5] = (1.7)*y[35]*(1-y[5])-11.9e-3*y[5]; // fCa_junc
    // ydot[6] = 1.7*y[36]*(1-y[6])-11.9e-3*y[6]; // fCa_sl
    ydot[3] = 0;
    ydot[4] = 0;
    ydot[5] = 0;
    ydot[6] = 0;
    
    double ibarca_j = ( pCa * 4. * ( y[38] * Frdy * FoRT )
                       * ( 0.341 * y[35] * exp( 2. * y[38] * FoRT ) - 0.341 * Cao )
                       / ( exp( 2. * y[38] * FoRT ) - 1. ) );
    double ibarca_sl = ( pCa * 4. * ( y[38] * Frdy * FoRT )
                        * ( 0.341 * y[36] * exp( 2. * y[38] * FoRT ) - 0.341 * Cao )
                        / ( exp( 2. * y[38] * FoRT ) - 1. ) );
    double ibark = ( pK * ( y[38] * Frdy * FoRT )
                    * ( 0.75 * y[34] * exp( y[38] * FoRT ) - 0.75 * Ko )
                    / ( exp( y[38] * FoRT ) - 1. ) );
    double ibarna_j = ( pNa * ( y[38] * Frdy * FoRT )
                       * ( 0.75 * y[31] * exp( y[38] * FoRT ) - 0.75 * Nao )
                       / ( exp( y[38] * FoRT ) - 1. ) );
    double ibarna_sl = ( pNa * ( y[38] * Frdy * FoRT )
                        * ( 0.75 * y[32] * exp( y[38] * FoRT ) - 0.75 * Nao)
                        / ( exp( y[38] * FoRT ) - 1. ) );
    
    double I_Ca_junc1 = ( Fjunc_CaL * ibarca_j * y[3] * y[4]
                         * ( 1. - y[5] ) * pow( Q10CaL, Qpow ) ) * 0.45 ;
    double I_Ca_sl1 = ( Fsl_CaL * ibarca_sl * y[3] * y[4]
                       * ( 1. - y[6] ) * pow( Q10CaL, Qpow ) ) * 0.45;
    double I_CaK1 = ( ibark * y[3] * y[4]
                     * ( Fjunc_CaL * ( 1. - y[5] ) + Fsl_CaL * ( 1 - y[6] ) )
                     * pow( Q10CaL, Qpow ) ) * 0.45;
    double I_CaNa_junc1 = ( Fjunc_CaL * ibarna_j * y[3] * y[4]
                           * ( 1. - y[5] ) * pow( Q10CaL, Qpow ) ) * 0.45;
    double I_CaNa_sl1 = ( Fsl_CaL * ibarna_sl * y[3] * y[4]
                         * ( 1. - y[6] ) * pow( Q10CaL, Qpow ) ) * 0.45;
    
    //// LCC MARKOV MODEL - based on Mahajan et al. (2008)
    
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
    // either mode 2 at the expense of mode 1 channels (i.e. 10// mode 2 results
    // in 90// mode 1).
    
    // PKA alters overall availability of channels (favail term that changes
    // overall scaling factor for currents) and also shifts distribution of
    // mode1/2 channels. PKA actions act on both junctional and sarcolemmal
    // channels.
    
    // To allow for CDI KO
    double cajLCC = y[35];
    double caslLCC = y[36];
    
    // LCC Current Fixed Parameters
    double taupo = 1.;          // [ms] - Time constant of activation
    double TBa = 450.;          // [ms] - Time constant
    double s1o = .0221;
    double k1o = .03;
    double kop = 2.5e-3;       // [mM]
    double cpbar = 8.e-3;       // [mM]
    double tca = 78.0312;
    double ICa_scale = 5.25;
    double recoveryReduc = 3.;
    
    // PKA PHOSPHOREGULATION OF LCC AVAILABLILITY (beta subunit phosph)
    double fracLCCbp0 = 0.250657; // Derived quantity - (LCCbp(baseline)/LCCbtot)
    double fracLCCbpISO = 0.525870; // Derived quantity - (LCCbp(ISO)/LCCbtot)
    //a_favail=(1.50-1)/(fracLCCbpISO/fracLCCbp0-1); // fracLCCbp ISO
    double a_favail = ( 1.56 - 1. ) / ( fracLCCbpISO / fracLCCbp0 - 1. ); // fracLCCbp ISO (x1.56 o.1 ISO)
    double favail = ( 1. - a_favail ) + a_favail * ( LCCb_PKAp / fracLCCbp0 ); // Test (max x2.52 100// phosph)
    //favail = 1; // no PKA effect on LTCCb
    ICa_scale =  ICa_scale * favail;
    
    double SSAshift = 0.;
    double SSIshift = 0.;
    // Voltage- and Ca-dependent Parameters
    double poss = 1. / ( 1. + exp( -( y[38] + SSAshift ) / 8. ) );
    double fcaj = 1. / ( 1. + pow( ( kop / cajLCC ), 3. ) );
    double Rv = 10. + 4954. * exp( y[38] / 15.6 );
    double PrLCC = 1. - 1. / ( 1. + exp( -( y[38] + 40. ) / 4. ) );
    double PsLCC = 1. / ( 1. + exp( -( y[38] + 40. + SSIshift ) / 11.32 ) );
    double TCaj = ( tca + 0.1 * ( 1. + pow( ( cajLCC / cpbar ), 2. ) ) ) / ( 1. + pow( ( cajLCC / cpbar ), 2. ) );
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
    double k2 = 1.e-4;                              // [ms] - Inactivation rate
    double k2p = .00224;                           // [ms] - Inactivation rate
    double s2 = s1 * ( k2 / k1 ) * ( r1 / r2 );
    double s2p = s1p * ( k2p / k1p ) * ( r1 / r2 );
    double k3 = exp( -( y[38] + 40. ) / 3. ) / ( 3. * ( 1. + exp( -( y[38] + 40. ) / 3. ) ) );
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
    pars1->gates[0] = s1;
    pars1->gates[1] = k1;
    // gates(1,tStep) = s1;
    // gates(2,tStep) = k1;
    
    // State transitions for MODE 1 junctional LCCs
    // O = no differential; C2 = 60; C1 = 61; I1Ca = 62; I2Ca = 63;
    // I1Ba = 64; I2Ba = 65;
    double Po_LCCj_m1 = 1.0 - y[59] - y[60] - y[61] - y[62] - y[63] - y[64];                                           // O_m1j
    ydot[59] = betaLCC * y[60] + k5 * y[62] + k5p * y[64] - ( k6 + k6p + alphaLCC ) * y[59];                      // C2_m1j
    ydot[60] = alphaLCC * y[59] + k2 * y[61] + k2p * y[63] + r2 * Po_LCCj_m1 - ( r1 + betaLCC + k1 + k1p ) * y[60];   // C1_m1j
    ydot[61] = k1 * y[60] + k4 * y[62] + s1 * Po_LCCj_m1 - ( k2 + k3 + s2 ) * y[61];                              // I1Ca_m1j
    ydot[62] = k3 * y[61] + k6 * y[59] - ( k4 + k5 ) * y[62];                                                 // I2Ca_m1j
    ydot[63] = k1p * y[60] + k4p * y[64] + s1p * Po_LCCj_m1 - ( k2p + k3p + s2p ) * y[63];                        // I1Ba_m1j
    ydot[64] = k3p * y[63] + k6p * y[59] - ( k5p + k4p ) * y[64];                                             // I2Ba_m1j
    double ibarca_jm1 = ( ( 4. * pCa * y[38] * Frdy * FoRT )
                         * ( .001 * exp( 2. * y[38] * FoRT ) - 0.341 * Cao )
                         / ( exp( 2. * y[38] * FoRT ) - 1. ) );
    double I_Ca_junc_m1 = ( Fjunc_CaL * ibarca_jm1 * Po_LCCj_m1 * pow( Q10CaL, Qpow ) ) * ICa_scale;
    
    // Re-define all parameters as mode 2 specific parameters
    double s1om2 = .0221;
    double k1om2 = .03;
    double kopm2 = 2.5e-3;
    double cpbarm2 = 8.e-3;
    double tcam2 = 78.0312;
    
    double possm2 = 1. / ( 1. + exp( -( y[38] + SSAshift ) / 8. ) );
    double fcajm2 = 1. / ( 1. + pow( ( kopm2 / cajLCC ), 3. ) ); // Depends on junctional Ca
    double Rvm2 = 10. + 4954. * exp( y[38] / 15.6 );
    double PrLCCm2 = 1. - 1. / ( 1. + exp( -( y[38] + 40. ) / 4. ) );
    double PsLCCm2 = 1. / ( 1. + exp( -( y[38] + 40. + SSIshift ) / 11.32 ) );
    double TCajm2 = ( ( tcam2 + 0.1 * ( 1. + pow( ( cajLCC / cpbarm2 ), 2. ) ) )
                     / ( 1. + pow( ( cajLCC / cpbarm2 ), 2. ) ) ); // Caj dependent
    double tauCajm2 = ( Rvm2 - TCajm2 ) * PrLCCm2 + TCajm2; // Caj dependence
    double tauBam2 = ( Rvm2 - TBa ) * PrLCCm2 + TBa;
    
    double alphaLCCm2 = possm2 / taupo;
    double betaLCCm2 = ( 1. - possm2 ) / taupo;
    double r1m2 = 0.3;                               // [1/ms] - Opening rate
    double r2m2 = 3. / 8.; // [1/ms] - closing rate,  changed from rabbit (3/10) - MOUSE
    double s1m2 = s1om2 * fcajm2;
    double s1pm2 = .00195;                           // [ms] - Inactivation rate
    double k1m2 = k1om2 * fcajm2;
    double k1pm2 = .00413;                           // [ms] - Inactivation rate
    double k2m2 = 1.e-4;                              // [ms] - Inactivation rate
    double k2pm2 = .00224;                           // [ms] - Inactivation rate
    double s2m2 = s1m2 * ( k2m2 / k1m2 ) * ( r1m2 / r2m2 );
    double s2pm2 = s1pm2 * ( k2pm2 / k1pm2 ) * ( r1m2 / r2m2 );
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
    
    // State transitions for MODE 2 junctional LCCs
    // O = no differential; C2 = 66; C1 = 67; I1Ca = 68; I2Ca = 69;
    // I1Ba = 70; I2Ba = 71;
    double Po_LCCj_m2 = 1.0 - y[65] - y[66] - y[67] - y[68] - y[69] - y[70];                                                           // O_m2j
    ydot[65] = betaLCCm2 * y[66] + k5m2 * y[68] + k5pm2 * y[70] - ( k6m2 + k6pm2 + alphaLCCm2 ) * y[65];                          // C2_m2j
    ydot[66] = alphaLCCm2 * y[65] + k2m2 * y[67] + k2pm2 * y[69] + r2m2 * Po_LCCj_m2 - ( r1m2 + betaLCCm2 + k1m2 + k1pm2 ) * y[66];   // C1_m2j
    ydot[67] = k1m2 * y[66] + k4m2 * y[68] + s1m2 * Po_LCCj_m2 - ( k2m2 + k3m2 + s2m2 ) * y[67];                                  // I1Ca_m2j
    ydot[68] = k3m2 * y[67] + k6m2 * y[65] - ( k4m2 + k5m2 ) * y[68];                                                         // I2Ca_m2j
    ydot[69] = k1pm2 * y[66] + k4pm2 * y[70] + s1pm2 * Po_LCCj_m2 - ( k2pm2 + k3pm2 + s2pm2 ) * y[69];                            // I1Ba_m2j
    ydot[70] = k3pm2 * y[69] + k6pm2 * y[65] - ( k5pm2 + k4pm2 ) * y[70];                                                     // I2Ba_m2j
    double ibarca_jm2 = ( ( 4. * pCa * y[38] * Frdy * FoRT )
                         * ( .001 * exp( 2. * y[38] * FoRT ) - 0.341 * Cao )
                         / ( exp( 2. * y[38] * FoRT ) - 1. ) );
    double I_Ca_junc_m2 = ( Fjunc_CaL * ibarca_jm2 * ( Po_LCCj_m2 ) * pow( Q10CaL, Qpow ) ) * ICa_scale;
    
    // CaMKII AND PKA-DEPENDENT SHIFTING OF DYADIC LCCS TO MODE 2
    //fpkam2 = 0.1543*LCCa_PKAp - .0043; // Assumes max phosphorylation results in 15// mode 2
    double fracLCCap0 = 0.219577; // Derived
    double frac_fpkam2 = ( 0.15 * fracLCCap0 ) / ( 1. - fracLCCap0 );
    double fpkam2 = ( 0.15 + frac_fpkam2 ) * LCCa_PKAp - frac_fpkam2; // Assumes max (100//) phosphorylation results in 15// mode 2 channels
    //(fpkam2 = 0 with NO ISO)
    //fpkam2 = 0; // no PKA effect on LTCCa
    //fpkam2 = 0.0765*0.8;
    double fckiim2 = LCC_CKp * .1; // Assumes max phosphorylation results in 10// mode 2 channels (max LCC_CKp = 1)
    // Sum up total fraction of CKII and PKA-shifted mode 2 channels
    double junc_mode2 = fckiim2 + fpkam2;
    // Total junctional ICa
    double I_Ca_junc2 = ( 1. - junc_mode2 ) * I_Ca_junc_m1 + junc_mode2 * I_Ca_junc_m2;
    
    // SUB-SARCOLEMMAL LCCs
    // Re-assign necessary params to be Casl sensitive
    double fcasl = 1. / ( 1. + pow( ( kop / caslLCC ), 3. ) );    // Depends on sl Ca
    double TCasl = ( ( tca + 0.1 * pow( ( 1. + ( caslLCC / cpbar ) ), 2. ) )
                    / ( 1. + pow( ( caslLCC / cpbar ), 2. ) ) );
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
    ydot[71] = betaLCC * y[72] + k5sl * y[74] + k5p * y[76] - ( k6sl + k6p + alphaLCC ) * y[71];                      // C2_m1sl
    ydot[72] = alphaLCC * y[71] + k2 * y[73] + k2p * y[75] + r2 * Po_LCCsl_m1 - ( r1 + betaLCC + k1sl + k1p ) * y[72];    // C1_m1sl
    ydot[73] = k1sl * y[72] + k4sl * y[74] + s1sl * Po_LCCsl_m1 - ( k2 + k3 + s2sl ) * y[73];                         // I1Ca_m1sl
    ydot[74] = k3 * y[73] + k6sl * y[71] - ( k4sl + k5sl ) * y[74];                                               // I2Ca_m1sl
    ydot[75] = k1p * y[72] + k4psl * y[76] + s1p * Po_LCCsl_m1 - ( k2p + k3p + s2psl ) * y[75];                       // I1Ba_m1sl
    ydot[76] = k3p * y[75] + k6p * y[71] - ( k5p + k4psl ) * y[76];                                               // I2Ba_m1sl
    double ibarca_slm1 = ( 4. * pCa * y[38] * Frdy * FoRT ) * ( .001 * exp( 2. * y[38] * FoRT ) - 0.341 * Cao ) / ( exp( 2. * y[38] * FoRT ) - 1. );
    double I_Casl_m1 = ( Fsl_CaL * ibarca_slm1 * Po_LCCsl_m1 * pow( Q10CaL, Qpow ) ) * ICa_scale;
    
    // Adjust closing rate for 'mode 2' sarcolemmal LCCs
    double r2slm2 = r2m2;
    double s2slm2 = s1sl * ( k2 / k1sl ) * ( r1 / r2slm2 );
    double s2pslm2 = s1p * ( k2p / k1p ) * ( r1 / r2slm2 );
    
    // State transitions for mode 2 sarcolemmal LCCs
    // O = no differential; C2 = 78; C1 = 79; I1Ca = 80; I2Ca = 81; I1Ba = 82; I2Ba = 83
    double Po_LCCsl_m2 = 1. - y[77] - y[78] - y[79] - y[80] - y[81] - y[82];                                                // O_m2sl
    ydot[77] = betaLCC * y[78] + k5sl * y[80] + k5p * y[82] - ( k6sl + k6p + alphaLCC ) * y[77];                      // C2_m2sl
    ydot[78] = alphaLCC * y[77] + k2 * y[79] + k2p * y[81] + r2slm2 * Po_LCCsl_m2 - ( r1 + betaLCC + k1sl + k1p ) * y[78];// C1_m2sl
    ydot[79] = k1sl * y[78] + k4sl * y[80] + s1sl * Po_LCCsl_m2 - ( k2 + k3 + s2slm2 ) * y[79];                       // I1Ca_m2sl
    ydot[80] = k3 * y[79] + k6sl * y[77] - ( k4sl + k5sl ) * y[80];                                               // I2Ca_m2sl
    ydot[81] = k1p * y[78] + k4psl * y[82] + s1p * Po_LCCsl_m2 - ( k2p + k3p + s2pslm2 ) * y[81];                     // I1Ba_m2sl
    ydot[82] = k3p * y[81] + k6p * y[77] - ( k5p + k4psl ) * y[82];                                               // I2Ba_m2sl
    double ibarca_slm2 = ( 4. * pCa * y[38] * Frdy * FoRT ) * ( .001 * exp( 2. * y[38] * FoRT ) - 0.341 * Cao ) / ( exp( 2. * y[38] * FoRT ) - 1. );
    double I_Casl_m2 = ( Fsl_CaL * ibarca_slm2 * Po_LCCsl_m2 * pow( Q10CaL, Qpow ) ) * ICa_scale;
    
    // Sum mode 1 and mode 2 sl channels for total sl current
    double fckiim2_sl = 0.; // Set to zero since SL LCCp by CaMKII is negligible
    double sl_mode2 = fckiim2_sl + fpkam2;
    double I_Ca_sl2 = ( 1. - sl_mode2 ) * I_Casl_m1 + sl_mode2 * I_Casl_m2;
    
    // Na and K currents through LCC
    double I_CaKj2 = ibark * Fjunc_CaL * ( ( 1. - junc_mode2 ) * Po_LCCj_m1 + junc_mode2 * Po_LCCj_m2 ) * pow( Q10CaL, Qpow ) * ICa_scale;
    double I_CaKsl2 = ibark * Fsl_CaL * ( ( 1. - sl_mode2 ) * Po_LCCsl_m1 + sl_mode2 * Po_LCCsl_m2 ) * pow( Q10CaL, Qpow ) * ICa_scale;
    double I_CaK2 = I_CaKj2 + I_CaKsl2;
    double I_CaNa_junc2 = ( Fjunc_CaL * ibarna_j * ( ( 1. - junc_mode2 ) * Po_LCCj_m1 + junc_mode2 * Po_LCCj_m2 ) * pow( Q10CaL, Qpow ) ) * ICa_scale;
    double I_CaNa_sl2 = Fsl_CaL * ibarna_sl * ( ( 1. - sl_mode2 ) * Po_LCCsl_m1 + sl_mode2 * Po_LCCsl_m2 ) * pow( Q10CaL, Qpow ) * ICa_scale;
    
    // These are now able to switch depending on whether or not the flag to
    // switch to Markov model of ICa is ON
    double I_Ca_junc = ( 1. - ICa_MarkovFlag ) * I_Ca_junc1 + ICa_MarkovFlag * I_Ca_junc2;
    double I_Ca_sl = ( 1. - ICa_MarkovFlag ) * I_Ca_sl1 + ICa_MarkovFlag * I_Ca_sl2;
    double I_Ca = I_Ca_junc + I_Ca_sl;   // Total Ca curren throuhgh LCC
    double I_CaNa_junc = ( 1. - ICa_MarkovFlag ) * ( I_CaNa_junc1 ) + ( ICa_MarkovFlag ) * ( I_CaNa_junc2 );
    double I_CaNa_sl = ( 1. - ICa_MarkovFlag ) * ( I_CaNa_sl1 ) + ( ICa_MarkovFlag ) * ( I_CaNa_sl2 );
    double I_CaNa = I_CaNa_junc + I_CaNa_sl;   // Total Na current through LCC
    double I_CaK = ( 1. - ICa_MarkovFlag ) * (I_CaK1) + ICa_MarkovFlag * (I_CaK2);  // Total K current through LCC
    
    // Collect all currents through LCC
    double I_Catot = I_Ca + I_CaK + I_CaNa;
    ydot[42] = -I_Ca * Cmem / ( Vmyo * 2. * Frdy ) * 1.e3;
    
    // global I_Ca_store ibar_store
    pars1->I_Ca_store = I_Catot;
    pars1->ibar_store = ibarca_j;
    // I_Ca_store(tStep) = I_Catot;
    // ibar_store(tStep) = ibarca_j;
    //// I_ncx: Na/Ca Exchanger flux
    
    double Ka_junc = 1. / ( 1. + pow( ( Kdact / y[35] ), 3. ) );
    double Ka_sl = 1. / ( 1. + pow( ( Kdact / y[36] ), 3. ) );
    double s1_junc = exp( nu * y[38] * FoRT ) * y[31]* y[31]* y[31] * Cao;
    double s1_sl = exp( nu * y[38] * FoRT ) * y[32]* y[32]* y[32] * Cao;
    double s2_junc = exp( ( nu - 1. ) * y[38] * FoRT ) * Nao* Nao* Nao * y[35];
    double s3_junc = ( ( KmCai * Nao* Nao* Nao * ( 1. + pow( ( y[31] / KmNai ), 3. ) )
                        + KmNao * KmNao * KmNao * y[35]
                        + KmNai * KmNai * KmNai * Cao * ( 1. + y[35] / KmCai )
                        + KmCao * y[31] * y[31] * y[31]
                        + y[31] * y[31] * y[31] * Cao
                        + Nao * Nao * Nao * y[35] )
                      * ( 1. + ksat * exp( ( nu - 1. ) * y[38] * FoRT ) ) );
    double s2_sl = exp( ( nu - 1. ) * y[38] * FoRT ) * Nao * Nao * Nao * y[36];
    double s3_sl = ( ( KmCai * Nao * Nao * Nao * ( 1. + pow( ( y[32] / KmNai ), 3 ) )
                      + KmNao * KmNao * KmNao * y[36]
                      + KmNai * KmNai * KmNai * Cao * ( 1. + y[36] / KmCai )
                      + KmCao * y[32] * y[32] * y[32]
                      + y[32] * y[32] * y[32] * Cao
                      + Nao *  Nao * Nao * y[36] )
                    * ( 1. + ksat * exp( ( nu - 1. ) * y[38] * FoRT ) ) );
    double I_ncx_junc = Fjunc_ncx * IbarNCX * pow( Q10NCX, Qpow ) * Ka_junc * ( s1_junc - s2_junc ) / s3_junc;
    double I_ncx_sl = Fsl_ncx * IbarNCX * pow( Q10NCX, Qpow ) * Ka_sl * ( s1_sl - s2_sl ) / s3_sl;
    double I_ncx = I_ncx_junc + I_ncx_sl;
    ydot[44] = 2. * I_ncx * Cmem / ( Vmyo * 2. * Frdy ) * 1.e3; //uM/ms
    
    // global Incx
    pars1->Incx = I_ncx;
    // Incx(tStep) = I_ncx;
    //// I_pca: Sarcolemmal Ca Pump Current
    
    double I_pca_junc = ( Fjunc * pow( Q10SLCaP, Qpow )
                         * IbarSLCaP * pow( y[35], 1.6 )
                         / ( pow( KmPCa, 1.6 ) + pow( y[35], 1.6 ) ) );
    double I_pca_sl = ( Fsl * pow( Q10SLCaP, Qpow )
                       * IbarSLCaP * pow( y[36], 1.6 )
                       / ( pow( KmPCa, 1.6 ) + pow( y[36], 1.6 ) ) );
    double I_pca = I_pca_junc + I_pca_sl;
    ydot[43] = -I_pca * Cmem / ( Vmyo * 2. * Frdy ) * 1.e3;
    
    // global Ipca_store
    pars1->Ipca_store = I_pca;
    // Ipca_store(tStep) = I_pca;
    //// I_cabk: Ca Background Current
    
    double I_cabk_junc = Fjunc * GCaB * ( y[38] - eca_junc );
    double I_cabk_sl = Fsl * GCaB * ( y[38] - eca_sl );
    double I_cabk = I_cabk_junc + I_cabk_sl;
    ydot[45] = -I_cabk * Cmem / ( Vmyo * 2. * Frdy ) * 1.e3;
    //// I_CFTR or I_cl_(cAMP) - Cystic Fibrosis Transmembrane Conductance Reg.
    // This is an Em- and time-independent current that is activated by PKA
    // fact_pka_cftr = 1.1933*ICFTR_PKAp - 0.1933; // Derived?
    // gCFTR = fact_pka_cftr*4.9e-3; // [A/F] - Max value as in Shannon et al. (2005)
    double Icftr = 0.; //gCFTR*(y[38] - ecl); // NO Icftr in MOUSE
    
    // global ICFTR
    pars1->ICFTR = Icftr;
    // ICFTR(tStep) = Icftr;
    //// RyR model - SR release fluxes and leak
    
    // CaMKII and PKA-dependent phosphoregulation of RyR Po
    double fCKII_ec50SR = 1.16 - 4. / 5. * RyR_CKp;
    ec50SR = fCKII_ec50SR * ec50SR; // MOUSE - 60//
    
    double MaxSR = 15.;
    double MinSR = 1.;
    double kCaSR = MaxSR - ( MaxSR - MinSR ) / ( 1. + pow( ( ec50SR / y[30] ), 2.5 ) );
    double koSRCa = koCa / kCaSR;
    double kiSRCa = kiCa * kCaSR;
    double kleak = 2. * 5.348e-6; // [1/ms] changed from rabbit (5.348e-6)
    
    //fCKII_RyR = (20*RyR_CKp/3 - 1/3); // 1 at basal condition - RABBIT
    double fCKII_RyR = ( 10. * RyR_CKp - 1. ); // 1 at basal condition - MOUSE
    
    //fPKA_RyR = RyR_PKAp*1.025 + 0.9750; // 1 with NO ISO
    double frac_RyRo = 0.204276; // Derived (RyR_PKAp(basal)/RyRtot)
    double a_RyR = ( 2. - 1. ) / ( 1. / frac_RyRo - 1. ); // Max effect: fPKA_RyR=2
    double fPKA_RyR = 1. - a_RyR + a_RyR * ( RyR_PKAp / frac_RyRo );
    koSRCa = ( fCKII_RyR + fPKA_RyR - 1. ) * koSRCa;
    
    // ODEs for RyR states and SR release through open RyRs
    //    double RI = 1. - y[13] - y[14] - y[15];
    //    ydot[13] = ( ( kim * RI - kiSRCa * y[35] * y[13] )
    //                - ( koSRCa * y[35] * y[35] * y[13] - kom * y[14] ) );   // R
    //    ydot[14] = ( ( koSRCa * y[35] * y[35] * y[13] - kom * y[14] )
    //                - ( kiSRCa * y[35] * y[14] - kim * y[15] ) );// O
    //    ydot[15] = ( ( kiSRCa * y[35] * y[14] - kim * y[15] )
    //                - ( kom * y[15] - koSRCa * y[35] * y[35] * RI ) );   // I
    
    //    drug = 0.4 *(1e-6);
    
    kon = drug * 5500;
	
	double IC50_flecRyR = 2e-6; //test 5e-6, 11e-6, 17e-6
    double koff_flecRyR = IC50_flecRyR * 5500;
    
	double RI = 1. - y[13]-y[14]-y[15]-y[129];
    
    ydot[13] = ( ( kim * RI - kiSRCa * y[35] * y[13] )
                - ( koSRCa * y[35] * y[35] * y[13] - kom * y[14] ) );   // R
    
    ydot[14] = ( ( koSRCa * y[35] * y[35] * y[13] - kom * y[14] )
                - ( kiSRCa * y[35] * y[14] - kim * y[15] ) + koff_flecRyR*y[129]-kon*y[14]);// O
    
    ydot[15] = ( ( kiSRCa * y[35] * y[14] - kim * y[15] )
                - ( kom * y[15] - koSRCa * y[35] * y[35] * RI ) );   // I
    
    ydot[129] = ((kon*y[14]- koff_flecRyR * y[129]));
    
    
    double J_SRCarel = ks * (y[14] + 0.2 * y[129] ) * ( y[30] - y[35] );          // [mmol/L SR/ ms]
    
    // Passive RyR leak - includes CaMKII regulation of leak flux
    //kleak = (1/3 + 10*RyR_CKp/3)*kleak; // RABBIT
    kleak = ( 1. / 2. + 5. * RyR_CKp / 2. ) * kleak; // MOUSE (reduced CaMKII effect on leak)
    double J_SRleak = kleak * ( y[30] - y[35] ); // [mmol/L cyt/ms]
    
    // global Jleak
    pars1->Jleak[0] = J_SRCarel*Vsr/Vmyo + J_SRleak; // Total leak [mmol/L cyt/ms]
    pars1->Jleak[1] = J_SRleak;                      // Passive leak only [mmol/L cyt/ms]
    // Jleak(tStep,1) = J_SRCarel*Vsr/Vmyo + J_SRleak; // Total leak [mmol/L cyt/ms]
    // Jleak(tStep,2) = J_SRleak;                      // Passive leak only [mmol/L cyt/ms]
    //// SERCA model - SR uptake fluxes
    
    // CaMKII and PKA-dependent phosphoregulation of PLB (changes to SERCA flux)
    double fCKII_PLB = ( 1. - .5 * PLB_CKp ); // Max effect: fCKII_PLB=0.5
    double fracPKA_PLBo = 1. - 0.079755; // Derived quantity - (1 - (PLBp(baseline)/PLBtot))
    //fPKA_PLB = (PLB_PKAn/fracPKA_PLBo)*3/4 + 1/4; // Max effect: fPKA_PLB=0.25
    double fPKA_PLB = ( PLB_PKAn / fracPKA_PLBo ) * ( 100. - 55.31 ) / 100. + 55.31 / 100.; // Max effect: fPKA_PLB=0.45
    
    // Select smaller value (resulting in max reduction of Kmf)
    
    if ( fCKII_PLB < fPKA_PLB ) {
        Kmf = Kmf * fCKII_PLB;//fCKII_PLB
    } else if ( fPKA_PLB < fCKII_PLB ) {
        Kmf = Kmf * fPKA_PLB;//fPKA_PLB
    } // end
    
    double J_serca = ( pow( Q10SRCaP, Qpow ) * Vmax_SRCaP
                      * ( pow( ( y[37] / Kmf ), hillSRCaP )
                         - pow( ( y[30] / Kmr ), hillSRCaP ) )
                      / ( 1. + pow( ( y[37] / Kmf ), hillSRCaP )
                         + pow( ( y[30] / Kmr ), hillSRCaP ) ) ); // [mM/msec]
    
    // global Jserca
    pars1->Jserca = J_serca; // [mM/msec] or [mmol/L cyt msec]
    //Jserca(tStep) = J_serca; // [mM/msec] or [mmol/L cyt msec]
    //// Na and Ca Buffering
    
    ydot[16] = kon_na * y[31] * ( Bmax_Naj - y[16] ) - koff_na * y[16];        // NaBj      [mM/ms]
    ydot[17] = kon_na * y[32] * ( Bmax_Nasl - y[17] ) - koff_na * y[17];       // NaBsl     [mM/ms]
    
    // Cytosolic Ca Buffers
    ydot[18] = kon_tncl * y[37] * ( Bmax_TnClow - y[18] ) - koff_tncl * y[18];            // TnCL      [mM/ms]
    ydot[19] = kon_tnchca * y[37] * ( Bmax_TnChigh - y[19] - y[20] ) - koff_tnchca * y[19]; // TnCHc     [mM/ms]
    ydot[20] = kon_tnchmg * Mgi * ( Bmax_TnChigh - y[19] - y[20] ) - koff_tnchmg * y[20];   // TnCHm     [mM/ms]
    ydot[21] = 0.; // commented b/c buffering done by CaM module
    // kon_cam*y[37]*(Bmax_CaM-y[21])-koff_cam*y[21]; // CaM       [mM/ms]
    ydot[22] = kon_myoca * y[37] * ( Bmax_myosin - y[22] - y[23] ) - koff_myoca * y[22];    // Myosin_ca [mM/ms]
    ydot[23] = kon_myomg * Mgi * ( Bmax_myosin - y[22] - y[23] ) - koff_myomg * y[23];      // Myosin_mg [mM/ms]
    ydot[24] = kon_sr * y[37] * ( Bmax_SR - y[24] ) - koff_sr * y[24];                    // SRB       [mM/ms]
    //J_CaB_cytosol = sum(ydot(19:25)); // wrong formulation
    double J_CaB_cytosol = ydot[18] + ydot[19] + ydot[21] + ydot[22] + ydot[24];
    
    // Junctional and SL Ca Buffers
    ydot[25] = kon_sll * y[35] * ( Bmax_SLlowj - y[25] ) - koff_sll * y[25];       // SLLj      [mM/ms]
    ydot[26] = kon_sll * y[36] * ( Bmax_SLlowsl - y[26] ) - koff_sll * y[26];      // SLLsl     [mM/ms]
    ydot[27] = kon_slh * y[35] * ( Bmax_SLhighj - y[27] ) - koff_slh * y[27];      // SLHj      [mM/ms]
    ydot[28] = kon_slh * y[36] * ( Bmax_SLhighsl - y[28] ) - koff_slh * y[28];     // SLHsl     [mM/ms]
    double J_CaB_junction = ydot[25] + ydot[27];
    double J_CaB_sl = ydot[26] + ydot[28];
    //// Ion concentrations
    
    // SR Ca Concentrations
    ydot[29] =  kon_csqn * y[30] * ( Bmax_Csqn - y[29] ) - koff_csqn * y[29]; // Csqn      [mM/ms]
    ydot[30] = J_serca * Vmyo / Vsr - ( J_SRleak * Vmyo / Vsr + J_SRCarel ) - ydot[29]; // Ca_sr     [mM/ms] // Ratio 3 leak current
    
    // Na Concentrations
    double I_Na_tot_junc = I_Na_junc + I_nabk_junc + 3. * I_ncx_junc + 3. * I_nak_junc + I_CaNa_junc;   // [uA/uF]
    double I_Na_tot_sl = I_Na_sl + I_nabk_sl + 3. * I_ncx_sl + 3. * I_nak_sl + I_CaNa_sl;   //[uA/uF]
    ydot[31] = -I_Na_tot_junc * Cmem / ( Vjunc * Frdy ) + J_na_juncsl / Vjunc * ( y[32] - y[31] ) - ydot[16];
    ydot[32] = ( -I_Na_tot_sl * Cmem / ( Vsl * Frdy )
                + J_na_juncsl / Vsl * ( y[31] - y[32] )
                + J_na_slmyo / Vsl * ( y[33] - y[32] ) - ydot[17] );
    
    if ( NaClampFlag == 1 ) {
        ydot[33] = 0; // Na clamp
    } else {
        ydot[33] = J_na_slmyo / Vmyo * ( y[32] - y[33] ); // [mM/msec]
    } // end
    
    // K Concentration
    double I_K_tot = I_to + I_kr + I_ks + I_ki - 2. * I_nak + I_CaK + I_kp + I_kur + I_ss;     // [uA/uF]
    ydot[34] = 0.; //-I_K_tot*Cmem/(Vmyo*Frdy);           // [mM/msec]
    
    // Ca Concentrations
    double I_Ca_tot_junc = I_Ca_junc + I_cabk_junc + I_pca_junc - 2. * I_ncx_junc; // [uA/uF]
    double I_Ca_tot_sl = I_Ca_sl + I_cabk_sl + I_pca_sl - 2. * I_ncx_sl;           // [uA/uF]
    
    ydot[35] = ( -I_Ca_tot_junc * Cmem / ( Vjunc * 2. * Frdy )
                + J_ca_juncsl / Vjunc * ( y[36] - y[35] )
                - J_CaB_junction
                + ( J_SRCarel ) * Vsr / Vjunc
                + J_SRleak * Vmyo / Vjunc );   // Ca_j
    
    ydot[36] = ( -I_Ca_tot_sl * Cmem / ( Vsl * 2. * Frdy )
                + J_ca_juncsl / Vsl * ( y[35] - y[36] )
                + J_ca_slmyo / Vsl * ( y[37] - y[36] )
                - J_CaB_sl );   // Ca_sl
    
    ydot[37] = ( -J_serca - J_CaB_cytosol + J_ca_slmyo / Vmyo * ( y[36] - y[37] ) ); // Cai
    
    double junc_sl = J_ca_juncsl / Vsl * ( y[35] - y[36] );
    double sl_junc = J_ca_juncsl / Vjunc * ( y[36] - y[35] );
    double sl_myo = J_ca_slmyo / Vsl * ( y[37] - y[36] );
    double myo_sl = J_ca_slmyo / Vmyo * ( y[36] - y[37] );
    //// Simulation type
    
    //        protocol = 'pace';
    //
    //        if CaffeineFlag==1,
    //            protocol = 'vcRest';
    //        end
    //        if StrophFlag==1,
    //            protocol = 'none';
    //        end
    //
    //        // AP Waveform for AP clamp
    //        global AP
    //
    //        switch protocol
    //    case {'none',''},
    //        I_app = 0;
    //    case 'pace',    // pace w/ current injection at cycleLength 'cycleLength'
    //        if mod(t,cycleLength) <= 5
    //            I_app = 9.5;
    //        else
    //            I_app = 0.0;
    //        end
    //    case 'vcRest'   // Resting vclamp to equilibrate cell at resting potential
    //        V_clamp = -83;//y[38];//-80;
    //        R_clamp = .01;
    //        I_app = (V_clamp-y[38])/R_clamp;
    //    case 'vclamp_ICa', // Based on Kohlhaas (2006), Maier (2003)
    //        V_hold = -90;
    //        V_test = 20; // Use 0.05 for 0 mV to avoid numerical problems
    //        V_step = -50;
    //        INa_inactT = 30; // ms
    //        if mod(t,cycleLength) <= INa_inactT
    //            V_clamp = V_step;
    //        elseif mod(t,cycleLength)>INa_inactT && mod(t,cycleLength)<= 200 + INa_inactT
    //        V_clamp = V_test;
    //        else
    //            V_clamp = V_hold;
    //        end
    //        R_clamp = 0.01;
    //        I_app = (V_clamp-y[38])/R_clamp;
    //    case 'vclamp_ICa_Gain'
    //        V_hold = -80;
    //        V_test = 20;//0.05;
    //                    //         if t>1e3
    //                    //             V_test = -10;
    //                    //         end
    //        if mod(t,cycleLength)>= 5 &&  mod(t,cycleLength) <= 35//205
    //            V_clamp = V_test;
    //        else
    //            V_clamp = V_hold;
    //        end
    //        R_clamp = 0.01;
    //        I_app = (V_clamp-y[38])/R_clamp;
    //    case 'I-V_ICa'
    //        // Start with 5 pre-pulses at 0.5 Hz (200 ms each) from -90 to 0 mV
    //        V_hold = -90;
    //        V_condition = 0;
    //        V_step = -50;
    //        V_test = -30;
    //        if t < 10e3
    //            if mod(t,cycleLength) <= 200
    //                V_clamp = V_condition;
    //            else
    //                V_clamp = V_hold;
    //        end
    //            else
    //                if mod(t,cycleLength) <= 50
    //                    V_clamp = V_step;
    //        elseif mod(t,cycleLength) > 50 && mod(t,cycleLength) <= 250
    //        V_clamp = V_test;
    //                else
    //                    V_clamp = V_hold;
    //        end
    //        end
    //        R_clamp = 0.01;
    //        I_app = (V_clamp-y[38])/R_clamp;
    //    case 'IV_SSI_ICa'
    //        V_hold1 = -90;
    //        V_condition = variablePar;
    //        V_hold2 = V_hold1;
    //        V_step = 0.05;
    //        if mod(t,cycleLength) <= 250
    //            V_clamp = V_hold1;
    //        elseif mod(t,cycleLength) > 250 && mod(t,cycleLength) <= 2250
    //        V_clamp = V_condition;
    //        elseif mod(t,cycleLength) > 2250 && mod(t,cycleLength) <= 2255
    //        V_clamp = V_hold2;
    //        elseif mod(t,cycleLength) > 2255 && mod(t,cycleLength) <= 2455
    //        V_clamp = V_step;
    //        else
    //            V_clamp = V_hold1;
    //        end
    //        R_clamp = 0.01;
    //        I_app = (V_clamp-y[38])/R_clamp;
    //    case 'recovery_ICa' // Based on protocol from Li 1997 in ferret
    //        V_hold1 = -70;
    //        V_hold2 = -90;
    //        V_hold3 = -90; // Vary this term to change holding potential of experiment
    //        V_test = 0;
    //        if t <= 10e3
    //            if mod(t,cycleLength) <= 200
    //                V_clamp = V_test;
    //            else
    //                V_clamp = V_hold1;
    //        end
    //        elseif t>10e3 && t<= 12e3 // 2 s at -90 mV
    //        V_clamp = V_hold2;
    //        elseif t>12e3 && t<= 12.5e3 // First 500 ms pulse
    //        V_clamp = V_test;
    //        elseif t > 12.5e3 && t <= 12.5e3+rest // Variable rest interval
    //        V_clamp = V_hold3;
    //        elseif t > 12.5e3+rest && t <= 12.5e3+rest+500 // Second, 500 ms test pulse
    //        V_clamp = V_test;
    //            else
    //                V_clamp = V_hold2;
    //        end
    //        R_clamp = 0.01;
    //        I_app = (V_clamp-y[38])/R_clamp;
    //    case 'recovery_Ito'
    //        V_hold = -80;
    //        V_test = +50;
    //        if t <= 500
    //            V_clamp = V_test;
    //        elseif t > 500 && t <= 500 + rest
    //        V_clamp = V_hold;
    //        elseif t > 500 + rest && t <= 1e3 + rest
    //        V_clamp = V_test;
    //        else
    //            V_clamp = V_hold;
    //        end
    //        R_clamp = .01;
    //        I_app = (V_clamp-y[38])/R_clamp;
    //    case 'vclamp_INa'
    //        V_hold = -140;
    //        V_test = -20;
    //        if mod(t,cycleLength) <= 500
    //            V_clamp = V_test;
    //        else
    //            V_clamp = V_hold;
    //        end
    //        R_clamp = .01;
    //        I_app = (V_clamp-y[38])/R_clamp;
    //    case 'recovery_INa'
    //        V_hold = -140;
    //        V_test = -20;
    //        if t <= 1e3
    //            V_clamp = V_test;
    //        elseif t > 1e3 && t <= 1e3 + rest
    //        V_clamp = V_hold;
    //        elseif t > 1e3 + rest && t <= 1e3 + rest + 10
    //        V_clamp = V_test;
    //        else
    //            V_clamp = V_hold;
    //        end
    //        R_clamp = .01;
    //        I_app = (V_clamp-y[38])/R_clamp;
    //    case 'step', // SM
    //        V_hold = -90;
    //        V_test = 0.05; // Use 0.05 for 0 mV to avoid numerical problems
    //        V_step = -90;
    //        INa_inactT = 20; // delay (ms)
    //        if mod(t,cycleLength) <= INa_inactT
    //            V_clamp = V_step;
    //        elseif mod(t,cycleLength)>INa_inactT && mod(t,cycleLength)<= 200 + INa_inactT
    //        V_clamp = V_test;
    //        else
    //            V_clamp = V_hold;
    //        end
    //        R_clamp = 0.01;
    //        I_app = (V_clamp-y[38])/R_clamp;
    //    case 'AP_clamp'
    //        // Determine appropriate voltage by interpolating between data points
    //        // Note - use mod(t,1e3) for 1 Hz data
    //        ind1 = find(AP(:,1) <= mod(t,1e3),1,'last');
    //        ind2 = find(AP(:,1) >= mod(t,1e3),1,'first');
    //        tint = [AP(ind1,1),mod(t,1e3),AP(ind2,1)];
    //        APint = interp1(AP(:,1),AP(:,2),tint);
    //        potential = APint[1];
    //
    //        R_clamp = .01;
    //        I_app = (potential-y[38])/R_clamp;
    //        end
    
    
    //// Membrane Potential
    
    double I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;                 // [uA/uF]
    pars1->Nai = y[33];
    
    double I_Cl_tot = I_ClCa+I_Clbk+Icftr;                         // [uA/uF]
    double I_Ca_tot = I_Ca_tot_junc+I_Ca_tot_sl;                   // [uA/uF]
    double I_tot = I_Na_tot+I_Cl_tot+I_Ca_tot+I_K_tot;             // [uA/uF]
    //ydot[38] = -(I_tot-I_app);
    ydot[38] = -(I_tot);
    // global dVm_store
    //dVm_store = ydot[38];
    // dVm_store(tStep) = ydot[38];
}

