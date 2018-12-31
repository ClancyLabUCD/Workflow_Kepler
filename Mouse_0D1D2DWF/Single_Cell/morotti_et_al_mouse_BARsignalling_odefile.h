/*
// function ydot = morotti_et_al_mouse_barODEfile(t,y,pin)
// % This module describes the beta-adrenergic signaling pathway in mouse
// % ventricular myocyte, and this file was built upon the code developeded
// % by Yang and Saucerman.
// % Reference: Yang JH & Saucerman JJ. (2012). Phospholemman is a negative
// % feed-forward regulator of Ca2+ in beta-adrenergic signaling,
// % accelerating beta-adrenergic inotropy. Journal of Molecular and Cellular
// % Cardiology 52, 1048-1055.
// Re-implemented by Pei-Chi Yang on 06/22/2015
 */
// ************************************************************************
// ************************************************************************
#ifndef morotti_et_al_mouse_BARsignalling_odefile_H
#define morotti_et_al_mouse_BARsignalling_odefile_H

void morotti_et_al_mouse_BARsignalling_odefile( double *y, double *pin, double *ydot );

void morotti_et_al_mouse_BARsignalling_odefile( double *y, double *pin, double *ydot ) {
	
    // Assign passed in params
    
    double LR=y[0];
    double LRG=y[1];
    double RG=y[2];
    double b1AR_S464=y[3];
    double b1AR_S301=y[4];
    double GsaGTPtot=y[5];
    double GsaGDP=y[6];
    double Gsby=y[7];
    double AC_GsaGTP=y[8];
    double PDEp=y[9];
    double cAMPtot=y[10];
    double RC_I=y[11];
    double RCcAMP_I=y[12];
    double RCcAMPcAMP_I=y[13];
    double RcAMPcAMP_I=y[14];
    double PKACI=y[15];
    double PKACI_PKI=y[16];
    double  RC_II=y[17];
    double RCcAMP_II=y[18];
    double  RCcAMPcAMP_II=y[19];
    double  RcAMPcAMP_II=y[20];
    double  PKACII=y[21];
    double  PKACII_PKI=y[22];
    double  I1p_PP1=y[23];
    double  I1ptot=y[24];
    double  PLBp=y[25]; // output
    double  PLMp=y[26]; // output
    double  LCCap=y[27]; // output
    double  LCCbp=y[28]; // output
    double  RyRp=y[29]; // output
    double   TnIp=y[30]; // output
    double  KS79=y[31];
    double  KS80=y[32];
    double   KSp=y[33]; // output -> 0
    double  CFTRp=y[34]; // output -> 0
    double   KURp=y[35]; // output
    // Drug Concentrations
    
    double   ISO = pin[0]; // (uM) isoproterenol concentration - Ltot
    double   FSK = 0; // (uM) forskolin concentration
    double  IBMX = 0; // (uM) IBMX concentration
    // -------- SIGNALING MODEL -----------
    
    
    // b-AR module
    
    double b1ARtot = 0.00528;        // (uM) total b1-AR protein % MOUSE
    //b1ARtot=0.028; % RABBIT
    
    double kf_LR   = 1;             // (1/[uM ms]) forward rate for ISO binding to b1AR
    double kr_LR   = 0.285;          // (1/ms) reverse rate for ISO binding to b1AR
    double kf_LRG  = 1;             // (1/[uM ms]) forward rate for ISO:b1AR association with Gs
    double kr_LRG  = 0.062;         // (1/ms) reverse rate for ISO:b1AR association with Gs
    double kf_RG   = 1;              // (1/[uM ms]) forward rate for b1AR association with Gs
    double kr_RG   = 33;             // (1/ms) reverse rate for b1AR association with Gs
    
    double Gstot          = 3.83;           // (uM) total Gs protein
    double k_G_act        = 16e-3;          // (1/ms) rate constant for Gs activation
    double k_G_hyd        = 0.8e-3;         // (1/ms) rate constant for G-protein hydrolysis
    double k_G_reassoc    = 1.21;           // (1/[uM ms]) rate constant for G-protein reassociation
    
    double kf_bARK        = 1.1e-6;         // (1/[uM ms]) forward rate for b1AR phosphorylation by b1ARK
    double kr_bARK        = 2.2e-6;         // (1/ms) reverse rate for b1AR phosphorylation by b1ARK
    double kf_PKA         = 3.6e-6;         // (1/[uM ms]) forward rate for b1AR phosphorylation by PKA
    double kr_PKA         = 2.2e-6;         // (1/ms) reverse rate for b1AR phosphorylation by PKA
    
    double b1ARact = b1ARtot - b1AR_S464 - b1AR_S301;
    double b1AR = b1ARact - LR - LRG - RG;
    double Gs = Gstot - LRG - RG - Gsby;
    
    double dLR = kf_LR*ISO*b1AR - kr_LR*LR + kr_LRG*LRG - kf_LRG*LR*Gs;
    double dLRG = kf_LRG*LR*Gs - kr_LRG*LRG - k_G_act*LRG;
    double dRG = kf_RG*b1AR*Gs - kr_RG*RG - k_G_act*RG;
    
    double bARK_desens = kf_bARK*(LR+LRG);
    double bARK_resens = kr_bARK*b1AR_S464;
    double PKA_desens = kf_PKA*PKACI*b1ARact;
    double PKA_resens = kr_PKA*b1AR_S301;
    double db1AR_S464 = bARK_desens - bARK_resens; // ydot(5)
    double db1AR_S301 = PKA_desens - PKA_resens; // ydot(6)
    
    double G_act = k_G_act*(RG+LRG);
    double G_hyd = k_G_hyd*GsaGTPtot;
    double G_reassoc = k_G_reassoc*GsaGDP*Gsby;
    double dGsaGTPtot = G_act - G_hyd; // ydot(7)
    double dGsaGDP = G_hyd - G_reassoc; // ydot(8)
    double dGsby = G_act - G_reassoc; // ydot(9)
    // end b-AR module
    
    // cAMP module
    
    double ACtot           = 70.57e-3;       // (uM) total adenylyl cyclase % MOUSE
    // ACtot=47e-3; // RABBIT
    double ATP             = 5e3;            // (uM) total ATP
    double k_AC_basal      = 0.2e-3;         // (1/ms) basal cAMP generation rate by AC
    double Km_AC_basal     = 1.03e3;         // (uM) basal AC affinity for ATP
    
    double Kd_AC_Gsa       = 0.4;            // (uM) Kd for AC association with Gsa
    double kf_AC_Gsa       = 1;              // (1/[uM ms]) forward rate for AC association with Gsa
    double kr_AC_Gsa       = Kd_AC_Gsa;      // (1/ms) reverse rate for AC association with Gsa
    
    double k_AC_Gsa        = 8.5e-3;         // (1/ms) basal cAMP generation rate by AC:Gsa
    double Km_AC_Gsa       = 315.0;          // (uM) AC:Gsa affinity for ATP
    
    double Kd_AC_FSK       = 44.0;           // (uM) Kd for FSK binding to AC
    double k_AC_FSK        = 7.3e-3;         // (1/ms) basal cAMP generation rate by AC:FSK
    double Km_AC_FSK       = 860.0;          // (uM) AC:FSK affinity for ATP
    
    double PDEtot          = 22.85e-3;      // (uM) total phosphodiesterase
    double k_cAMP_PDE      = 5e-3;           // (1/ms) cAMP hydrolysis rate by PDE
    double k_cAMP_PDEp     = 2*k_cAMP_PDE;   // (1/ms) cAMP hydrolysis rate by phosphorylated PDE
    double Km_PDE_cAMP     = 1.3;            // (uM) PDE affinity for cAMP
    
    double Kd_PDE_IBMX     = 30.0;          // (uM) Kd_R2cAMP_C for IBMX binding to PDE
    double k_PKA_PDE       = 7.5e-3;         // (1/ms) rate constant for PDE phosphorylation by type 1 PKA
    double k_PP_PDE        = 1.5e-3;         // (1/ms) rate constant for PDE dephosphorylation by phosphatases
    
    double cAMP = cAMPtot - (RCcAMP_I+2*RCcAMPcAMP_I+2*RcAMPcAMP_I) - (RCcAMP_II+2*RCcAMPcAMP_II+2*RcAMPcAMP_II);
    double AC = ACtot-AC_GsaGTP;
    double GsaGTP = GsaGTPtot - AC_GsaGTP;
    double dAC_GsaGTP = kf_AC_Gsa*GsaGTP*AC - kr_AC_Gsa*AC_GsaGTP;
    
    double AC_FSK = FSK*AC/Kd_AC_FSK;
    double AC_ACT_BASAL = k_AC_basal*AC*ATP/(Km_AC_basal+ATP);
    double AC_ACT_GSA = k_AC_Gsa*AC_GsaGTP*ATP/(Km_AC_Gsa+ATP);
    double AC_ACT_FSK = k_AC_FSK*AC_FSK*ATP/(Km_AC_FSK+ATP);
    
    double PDE_IBMX = PDEtot*IBMX/Kd_PDE_IBMX;
    double PDE = PDEtot - PDE_IBMX - PDEp;
    double dPDEp = k_PKA_PDE*PKACII*PDE - k_PP_PDE*PDEp;
    double PDE_ACT = k_cAMP_PDE*PDE*cAMP/(Km_PDE_cAMP+cAMP) + k_cAMP_PDEp*PDEp*cAMP/(Km_PDE_cAMP+cAMP);
    
    double dcAMPtot = AC_ACT_BASAL + AC_ACT_GSA + AC_ACT_FSK - PDE_ACT; // ydot(15)
    // end cAMP module
    
    // PKA module
    
    double PKItot          = 0.18;           // (uM) total PKI
    double kf_RC_cAMP      = 1;              // (1/[uM ms]) Kd for PKA RC binding to cAMP
    double kf_RCcAMP_cAMP  = 1;             //(1/[uM ms]) Kd for PKA RC:cAMP binding to cAMP
    double kf_RcAMPcAMP_C  = 4.375;         // (1/[uM ms]) Kd for PKA R:cAMPcAMP binding to C
    double kf_PKA_PKI      = 1;              // (1/[uM ms]) Ki for PKA inhibition by PKI
    double kr_RC_cAMP      = 1.64;           //(1/ms) Kd for PKA RC binding to cAMP
    double kr_RCcAMP_cAMP  = 9.14;           // (1/ms) Kd for PKA RC:cAMP binding to cAMP
    double kr_RcAMPcAMP_C  = 1;              // (1/ms) Kd for PKA R:cAMPcAMP binding to C
    double kr_PKA_PKI      = 2e-4;           // (1/ms) Ki for PKA inhibition by PKI
    double epsilon         = 10;             //(-) AKAP-mediated scaling factor
    
    double PKI = PKItot - PKACI_PKI - PKACII_PKI;
    
    double dRC_I = - kf_RC_cAMP*RC_I*cAMP + kr_RC_cAMP*RCcAMP_I;
    double dRCcAMP_I = - kr_RC_cAMP*RCcAMP_I + kf_RC_cAMP*RC_I*cAMP - kf_RCcAMP_cAMP*RCcAMP_I*cAMP + kr_RCcAMP_cAMP*RCcAMPcAMP_I;
    double dRCcAMPcAMP_I = - kr_RCcAMP_cAMP*RCcAMPcAMP_I + kf_RCcAMP_cAMP*RCcAMP_I*cAMP - kf_RcAMPcAMP_C*RCcAMPcAMP_I + kr_RcAMPcAMP_C*RcAMPcAMP_I*PKACI;
    double dRcAMPcAMP_I = - kr_RcAMPcAMP_C*RcAMPcAMP_I*PKACI + kf_RcAMPcAMP_C*RCcAMPcAMP_I;
    double dPKACI = - kr_RcAMPcAMP_C*RcAMPcAMP_I*PKACI + kf_RcAMPcAMP_C*RCcAMPcAMP_I - kf_PKA_PKI*PKACI*PKI + kr_PKA_PKI*PKACI_PKI; // ydot(17)
    double dPKACI_PKI = - kr_PKA_PKI*PKACI_PKI + kf_PKA_PKI*PKACI*PKI;
    
    double dRC_II = - kf_RC_cAMP*RC_II*cAMP + kr_RC_cAMP*RCcAMP_II;
    double dRCcAMP_II = - kr_RC_cAMP*RCcAMP_II + kf_RC_cAMP*RC_II*cAMP - kf_RCcAMP_cAMP*RCcAMP_II*cAMP + kr_RCcAMP_cAMP*RCcAMPcAMP_II;
    double dRCcAMPcAMP_II = - kr_RCcAMP_cAMP*RCcAMPcAMP_II + kf_RCcAMP_cAMP*RCcAMP_II*cAMP - kf_RcAMPcAMP_C*RCcAMPcAMP_II + kr_RcAMPcAMP_C*RcAMPcAMP_II*PKACII;
    double dRcAMPcAMP_II = - kr_RcAMPcAMP_C*RcAMPcAMP_II*PKACII + kf_RcAMPcAMP_C*RCcAMPcAMP_II;
    double dPKACII = - kr_RcAMPcAMP_C*RcAMPcAMP_II*PKACII + kf_RcAMPcAMP_C*RCcAMPcAMP_II - kf_PKA_PKI*PKACII*PKI + kr_PKA_PKI*PKACII_PKI; // ydot(18)
    double dPKACII_PKI = - kr_PKA_PKI*PKACII_PKI + kf_PKA_PKI*PKACII*PKI;
    // end PKA module
    
    // I-1/PP1 module
    
    double PP1tot = pin[7]; // PP1tot = 0.89; // (uM) total phosphatase 1
    double I1tot           = 0.3;            // (uM) total inhibitor 1
    double k_PKA_I1        = 60e-3;         // (1/ms) rate constant for I-1 phosphorylation by type 1 PKA
    double Km_PKA_I1       = 1.0;            // (uM) Km for I-1 phosphorylation by type 1 PKA
    double Vmax_PP2A_I1    = 14.0e-3;        // (uM/ms) Vmax for I-1 dephosphorylation by PP2A
    double Km_PP2A_I1      = 1.0;            // (uM) Km for I-1 dephosphorylation by PP2A
    
    double Ki_PP1_I1       = 1.0e-3;         // (uM) Ki for PP1 inhibition by I-1
    double kf_PP1_I1       = 1;              // (uM) Ki for PP1 inhibition by I-1
    double kr_PP1_I1       = Ki_PP1_I1;      // (uM) Ki for PP1 inhibition by I-1
    
    double I1 = I1tot - I1ptot;
    double PP1 = PP1tot - I1p_PP1;
    double I1p = I1ptot - I1p_PP1;
    double I1_phosph = k_PKA_I1*PKACI*I1/(Km_PKA_I1+I1);
    double I1_dephosph = Vmax_PP2A_I1*I1ptot/(Km_PP2A_I1+I1ptot);
    double
    dI1p_PP1 = kf_PP1_I1*PP1*I1p - kr_PP1_I1*I1p_PP1;
    double dI1ptot = I1_phosph - I1_dephosph; // ydot(20)
    // end I-1/PP1 module
    // PLB module
    
    double PLBtot = pin[3]; //p(41) = PLBtot; // PLBtot    [uM]
    double k_PKA_PLB = 54e-3;     //p(44) = 54;     // k_pka_plb     [1/ms]
    double  Km_PKA_PLB = 21;     //p(45) = 21;     //Km_pka_plb    [uM]
    double k_PP1_PLB = 8.5e-3;    //p(46) = 8.5;    // k_pp1_plb     [1/ms]
    double Km_PP1_PLB = 7.0;    //p(47) = 7.0;    % Km_pp1_plb    [uM]
    
    double PLB = PLBtot - PLBp;
    double PLB_phosph = k_PKA_PLB*PKACI*PLB/(Km_PKA_PLB+PLB);
    double PLB_dephosph = k_PP1_PLB*PP1*PLBp/(Km_PP1_PLB+PLBp);
    double dPLBp = PLB_phosph - PLB_dephosph; // ydot(19)
    // end PLB module
    
    // PLM module (included 09/18/12) MOUSE
    
    double PLMtot = pin[9]; // p(102) = PLMtot; % PLMtot    [uM]
    double k_PKA_PLM = 54e-3; // p(103) = 54;     % k_pka_plb     [1/ms]
    double Km_PKA_PLM = 21; // p(104) = 21;     % Km_pka_plb    [uM]
    double k_PP1_PLM = 8.5e-3; // p(105) = 8.5;    % k_pp1_plb     [1/ms]
    double Km_PP1_PLM = 7.0; // p(106) = 7.0;    % Km_pp1_plb    [uM]
    
    double PLM = PLMtot - PLMp;
    double PLM_phosph = k_PKA_PLM*PKACI*PLM/(Km_PKA_PLM+PLM);
    double PLM_dephosph = k_PP1_PLM*PP1*PLMp/(Km_PP1_PLM+PLMp);
    double dPLMp = PLM_phosph - PLM_dephosph; // ydot(32)
    // end PLM module
    // LCC module
    
    double PKAIItot        = 0.059;          // (uM) total type 2 PKA % MOUSE
    // PKAIItot=0.084; // RABBIT
    
    double LCCtot = pin[1]; //p(53) = LCCtot; % LCCtot        [uM]
    double PKACII_LCCtot = 0.025; //p(54) = 0.025;  % PKAIIlcctot   [uM]
    double PP1_LCC = 0.025;  //p(55) = 0.025;  % PP1lcctot     [uM]
    double PP2A_LCC = 0.025;  //p(56) = 0.025;  % PP2Alcctot    [uM]
    double k_PKA_LCC = 54e-3;     //p(57) = 54;     % k_pka_lcc     [1/ms]
    double Km_PKA_LCC = 21;     //p(58) = 21;%*1.6;     % Km_pka_lcc    [uM]
    double k_PP1_LCC = 8.52e-3;   //p(59) = 8.52;   % k_pp1_lcc     [1/ms] RABBIT, MOUSE
    //p(59) = 8.5;   // k_pp1_lcc     [1/sec] RAT
    double Km_PP1_LCC = 3;      //p(60) = 3;      % Km_pp1_lcc    [uM]
    double k_PP2A_LCC = 10.1e-3;   //p(61) = 10.1;   % k_pp2a_lcc    [1/ms]
    double Km_PP2A_LCC = 3;      //p(62) = 3;      % Km_pp2a_lcc   [uM]
    
    double  PKACII_LCC = (PKACII_LCCtot/PKAIItot)*PKACII;
    double  LCCa = LCCtot - LCCap;
    double  LCCa_phosph = epsilon*k_PKA_LCC*PKACII_LCC*LCCa/(Km_PKA_LCC+epsilon*LCCa);
    double  LCCa_dephosph = epsilon*k_PP2A_LCC*PP2A_LCC*LCCap/(Km_PP2A_LCC+epsilon*LCCap);
    double  dLCCap = LCCa_phosph - LCCa_dephosph; // ydot(23)
    double  LCCb = LCCtot - LCCbp;
    double  LCCb_phosph = epsilon*k_PKA_LCC*PKACII_LCC*LCCb/(Km_PKA_LCC+epsilon*LCCb);
    double  LCCb_dephosph = epsilon*k_PP1_LCC*PP1_LCC*LCCbp/(Km_PP1_LCC+epsilon*LCCbp);
    double  dLCCbp = LCCb_phosph - LCCb_dephosph; // ydot(24)
    // end LCC module
    // RyR module (not included in Yang-Saucerman)
    
    double  RyRtot = pin[2]; //p(63) = RyRtot; // RyRtot        [uM]
    double  PKAIIryrtot = 0.034;  //p(64) = 0.034;  // PKAIIryrtot   [uM]
    double  PP1ryr = 0.034;  //p(65) = 0.034;  // PP1ryr        [uM]
    double  PP2Aryr = 0.034;  //p(66) = 0.034;  // PP2Aryr       [uM]
    double  kcat_pka_ryr = 54e-3;     //p(67) = 54;     // kcat_pka_ryr  [1/ms]
    double  Km_pka_ryr = 21;     //p(68) = 21;     // Km_pka_ryr    [uM]
    double  kcat_pp1_ryr = 8.52e-3;   //p(69) = 8.52;   // kcat_pp1_ryr  [1/ms]
    double  Km_pp1_ryr = 7;      //p(70) = 7;      // Km_pp1_ryr    [uM]
    double  kcat_pp2a_ryr = 10.1e-3;   //p(71) = 10.1;   // kcat_pp2a_ryr [1/ms]
    double  Km_pp2a_ryr = 4.1;    //p(72) = 4.1;    // Km_pp2a_ryr   [uM]
    
    double  PKACryr = (PKAIIryrtot/PKAIItot)*PKACII;
    double  RyR = RyRtot-RyRp;
    double  RyRPHOSPH = epsilon*kcat_pka_ryr*PKACryr*RyR/(Km_pka_ryr+epsilon*RyR);
    double  RyRDEPHOSPH1 = epsilon*kcat_pp1_ryr*PP1ryr*RyRp/(Km_pp1_ryr+epsilon*RyRp);
    double  RyRDEPHOSPH2A = epsilon*kcat_pp2a_ryr*PP2Aryr*RyRp/(Km_pp2a_ryr+epsilon*RyRp);
    double  dRyRp = RyRPHOSPH-RyRDEPHOSPH1-RyRDEPHOSPH2A; // ydot(25)
    // end RyR module
    //// TnI module
    
    double  TnItot = pin[4]; //p(73) = TnItot; // TnItot        [uM]
    double  PP2A_TnI = 0.67;   // PP2Atni       [uM]
    double k_PKA_TnI = 54e-3;     // kcat_pka_tni  [1/ms]
    double  Km_PKA_TnI = 21;     // Km_pka_tni    [uM]
    double  k_PP2A_TnI = 10.1e-3;   // kcat_pp2a_tni [1/ms]
    double  Km_PP2A_TnI = 4.1;    // Km_pp2a_tni   [uM]
    
    double  TnI = TnItot - TnIp;
    double  TnI_phosph = k_PKA_TnI*PKACI*TnI/(Km_PKA_TnI+TnI);
    double  TnI_dephosph = k_PP2A_TnI*PP2A_TnI*TnIp/(Km_PP2A_TnI+TnIp);
    double  dTnIp = TnI_phosph - TnI_dephosph; // ydot(26)
    // end TnI module
    //// Iks module (not present in mouse)
    
    double  IKstot = pin[5];
    // p(79) = IKstot; // Iks_tot       [uM]
    // p(80) = 0.025;  // Yotiao_tot    [uM]
    // p(81) = 0.1e-3; // K_yotiao      [uM] ** apply G589D mutation here **
    // p(82) = 0.025;  // PKAII_ikstot  [uM]
    // p(83) = 0.025;  // PP1_ikstot    [uM]
    // p(84) = 54;     // k_pka_iks     [1/sec]
    // p(85) = 21;     // Km_pka_iks    [uM]
    // p(86) = 8.52;   // k_pp1_iks     [1/sec]
    // p(87) = 7;      // Km_pp1_iks    [uM]
    //
    // IksYot = y(27)*y(28)/p(81);           // [uM]
    // ydot(27) = p(79) - IksYot - y(27);    // [uM]
    // ydot(28) = p(80) - IksYot - y(28);    // [uM]
    // PKACiks = (IksYot/p(79))*(p(82)/p(34))*y(18);
    // PP1iks = (IksYot/p(79))*p(83);
    // Iks = p(79)-y(29);
    // IKS_PHOSPH = p(40)*p(84)*PKACiks*Iks/(p(85)+p(40)*Iks);
    // IKS_DEPHOSPH = p(40)*p(86)*PP1iks*y(29)/(p(87)+p(40)*y(29));
    // ydot(29) = IKS_PHOSPH-IKS_DEPHOSPH;
    // end Iks module
    double  dKS79 = 0; // ydot(27) not ODE
    double  dKS80 = 0; // ydot(28) not ODE
    double  dKSp = 0; // ydot(29)
    //// CFTR module (included 04/30/10)
    
    double  ICFTRtot = pin[6]; // p(88) = ICFTRtot;  // ICFTR_tot      [uM]
    // PKAII_CFTRtot = 0.025;  //p(89) = 0.025;  // PKAII_CFTRtot [uM]
    // PP1_CFTRtot = 0.025;  //p(90) = 0.025;  // PP1_CFTRtot   [uM]
    // k_pka_CFTR = 54e-3;     //p(91) = 54;     // k_pka_CFTR    [1/ms]
    // Km_pka_CFTR = 8.5;    //p(92) = 8.5;    // Km_pka_CFTR   [uM]
    // k_pp1_CFTR = 8.52e-3;   //p(93) = 8.52;   // k_pp1_CFTR    [1/ms]
    // Km_pp1_CFTR = 7;      //p(94) = 7;      // Km_pp1_CFTR   [uM]
    //
    // CFTRn = ICFTRtot - CFTRp;  // Non-phos = tot - phos
    // PKAC_CFTR = (PKAII_CFTRtot/PKAIItot)*PKACII;    // (PKACFTRtot/PKAIItot)*PKAIIact
    // CFTRphos = epsilon*CFTRn*PKAC_CFTR*k_pka_CFTR/(Km_pka_CFTR+epsilon*CFTRn);
    // CFTRdephos = PP1_CFTRtot*k_pp1_CFTR*epsilon*CFTRp/(Km_pp1_CFTR + epsilon*CFTRp);
    double dCFTRp = 0; //CFTRphos - CFTRdephos; // ydot(30)
    // end CFTR module
    //// Ikur module (included 04/10/12) MOUSE
    
    double IKurtot = pin[8];// p(95) = IKurtot;   // Ikur_tot      [uM]
    double  PKAII_KURtot = 0.025; // p (96) = 0.025;  // PKAII_KURtot [uM]
    double  PP1_KURtot = 0.025; // p(97) = 0.025;  // PP1_KURtot   [uM]
    double  k_pka_KUR = 54e-3; // p(98) = 54;    // k_pka_KUR    [1/ms]
    double  Km_pka_KUR = 21; // p(99) = 21;    // Km_pka_KUR   [uM]
    double  k_pp1_KUR = 8.52e-3; // p(100) = 8.52;  // k_pp1_KUR    [1/ms]
    double  Km_pp1_KUR = 7; // p(101) = 7;     // Km_pp1_KUR   [uM]
    
    double  KURn = IKurtot - KURp;  // Non-phos = tot - phos
    double  PKAC_KUR = (PKAII_KURtot/PKAIItot)*PKACII;    // (PKA_KURtot/PKAIItot)*PKAIIact
    double  KURphos = epsilon*KURn*PKAC_KUR*k_pka_KUR/(Km_pka_KUR+epsilon*KURn);
    double  KURdephos = PP1_KURtot*k_pp1_KUR*epsilon*KURp/(Km_pp1_KUR+epsilon*KURp);
    double  dKURp = KURphos - KURdephos; // ydot(31)
   
    // end Ikur module
    //// ydot
    
    ydot[0]=dLR;
    ydot[1]=dLRG;
    ydot[2]=dRG;
    ydot[3]=db1AR_S464;
    ydot[4]=db1AR_S301;
    ydot[5]=dGsaGTPtot;
    ydot[6]=dGsaGDP;
    ydot[7]=dGsby;
    ydot[8]=dAC_GsaGTP;
    ydot[9]=dPDEp;
    ydot[10]=dcAMPtot;
    ydot[11]=dRC_I;
    ydot[12]=dRCcAMP_I;
    ydot[13]=dRCcAMPcAMP_I;
    ydot[14]=dRcAMPcAMP_I;
    ydot[15]=dPKACI;
    ydot[16]=dPKACI_PKI;
    ydot[17]=dRC_II;
    ydot[18]=dRCcAMP_II;
    ydot[19]=dRCcAMPcAMP_II;
    ydot[20]=dRcAMPcAMP_II;
    ydot[21]=dPKACII;
    ydot[22]=dPKACII_PKI;
    ydot[23]=dI1p_PP1; // output CaMKII
    ydot[24]=dI1ptot;
    ydot[25]=dPLBp; // output
    ydot[26]=dPLMp; // output
    ydot[27]=dLCCap; // output
    ydot[28]=dLCCbp; // output
    ydot[29]=dRyRp; // output
    ydot[30]=dTnIp; // output
    ydot[31]=dKS79;
    ydot[32]=dKS80;
    ydot[33]=dKSp; // output -> 0
    ydot[34]=dCFTRp; // output -> 0
    ydot[35]=dKURp; // output
    
    
}

#endif