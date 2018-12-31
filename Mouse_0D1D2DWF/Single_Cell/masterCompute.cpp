/*
 % This function calls the ode solver and plots results
 %
 % Re-implemented by Mao-Tsuen Jeng and modified by Pei-Chi Yang
 % for CPVT simulation
 %
 %  %
 %
 */
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <math.h>
#include <time.h>
#include <string.h>

// for mkdir
#include <sys/types.h>
#include <sys/stat.h>

typedef struct {
    double I_Na_store, I_Nabk_store, I_NaK_store;
    double I_kur1_store, I_kur2_store, I_ss_store, I_kr_store, IKs_store;
    double I_to_store, I_tof_store, I_tos_store, I_K1_store;
    double I_Ca_store, ibar_store, Ipca_store, Incx, ICFTR, Jleak[2], Jserca, gates[2];
    double JCaCyt, JCaSL, JCaDyad;
    double Nai;
} pars_rec;

typedef struct {
    double y0n[217];
    pars_rec pars1;
    double p[37];
    double v, v_old, ryr_old, I_Total;
    double dvdt, dcaidt, cai_old, dcasrdt, casr_old;
    
    double ui1_v, uj1_v, ui2_v, uj2_v, dun_v, tu_v; // For computing change of direction.
    double ui1_cai, uj1_cai, ui2_cai, uj2_cai, dun_cai, tu_cai; // For computing change of direction.
    double ui1_casr, uj1_casr, ui2_casr, uj2_casr, dun_casr, tu_casr; // For computing change of direction.
    
} Cell;



typedef struct {
    double t, tt;
    int tStep, counter, beat;
    Cell cellData[1][1];
} SimState;

SimState theState;
SimState *S = &theState;

const double pi = 3.141592653589793; // pi = 3.1415926535897932384626433832795028841971...

#include "integrate_rk2.h"
#include "morotti_et_al_mouse_masterODEfile.h"

int main( int argc, char *argv[]) {
    
    cout.precision(16);
    
    
    char name[30];
    
    
    Cell * theCell;
    pars_rec * pars1; // = &(theCell.pars1);
    double * y0n; // = theCell.y0n;
    double * y0n2;
    S->tStep = 0;
    double * p;
    
    char *InputFile2 = argv[2];
    double input2[5];
    
    FILE *fp2;
    
    fp2 = fopen( InputFile2, "r");
    
    if ( fp2 == NULL ) {
        
        printf( "%s%s%s", "\nFile \'", InputFile2 , "\' NOT found. \n \n" );
        
        exit(1);
    }
    
    cout << "stimulus Parameters: " << endl;
    
    for ( int idy = 0; idy < 5; idy++ ) {
        fscanf ( fp2, "%lf", &input2[idy] );
        cout << idy << "\t" << input2[idy] << endl;
    }
    
    char expression[2];
    fscanf ( fp2, "%s", expression );
    // exit(0);
    fclose( fp2 );
    
    double allDrugs[2];
    allDrugs[0] = input2[0];
    allDrugs[1] = input2[1];
    
    
    double cycleLength = input2[2];     // [ms]
    double freq = 1000. / cycleLength ; // [Hz] CHANGE DEPENDING ON FREQUENCY
    
    double CaMtotDyad = 418.;           // [uM]
    double BtotDyad = 1.54/8.293e-4;    // [uM]
    double CaMKIItotDyad = 120.;        // [uM]
    double CaNtotDyad = 3.e-3/8.293e-4; // [uM]
    double PP1totDyad = 96.5;           // [uM]
    double CaMtotSL = 5.65;             // [uM]
    double BtotSL = 24.2;               // [uM]
    double CaMKIItotSL = 120.*8.293e-4; // [uM]
    double CaNtotSL = 3.e-3;            // [uM]
    double PP1totSL = 0.57;             // [uM]
    double CaMtotCyt = 5.65;            // [uM]
    double BtotCyt = 24.2;              // [uM]
    double CaMKIItotCyt = 120.*8.293e-4;// [uM]
    double CaNtotCyt = 3.e-3;           // [uM]
    double PP1totCyt = 0.57;            // [uM]
    
    // ADJUST CAMKII ACTIVITY LEVELS (expression = 'WT', 'OE', or 'KO')
    
    double CKIIOE = 0; // Should be zero during 'WT' and 'KO' runs
    
//    if ( expression == "OE" ) {
    if ( strcmp( expression , "OE" ) == 0 ) {
        int CKIIOE = 1; // Flag for CKII OE in ECC file (0=WT, 1=OE) - for Ito and INa
        CaMKIItotDyad = 120*6;        // [uM]
        CaMKIItotSL = 120*8.293e-4*6; // [uM]
        CaMKIItotCyt = 120*8.293e-4*6;// [uM]
//    }	else if ( expression == "KO" ) {
    } else if ( strcmp( expression , "KO" ) == 0 ) {
        CaMKIItotDyad = 0;          // [uM]
        CaMKIItotSL = 0;            // [uM]
        CaMKIItotCyt = 0;           // [uM]
    }
    // end
    
    // For Recovery from inactivation of LCC
    double recoveryTime = 10;  // initialize to smallest value
    double variablePar = 20; // initilization
    //plb_val=38; % RABBIT
    double plb_val=106; // MOUSE
    // Parameters for CaMKII module
    double LCCtotDyad = 31.4*.9;       // [uM] - Total Dyadic [LCC] - (umol/l dyad)
    double LCCtotSL = 0.0846;          // [uM] - Total Subsarcolemmal [LCC] (umol/l sl)
    double RyRtot = 382.6;             // [uM] - Total RyR (in Dyad)
    double PP1_dyad = 95.7;            // [uM] - Total dyadic [PP1]
    double PP1_SL = 0.57;              // [uM] - Total Subsarcolemmal [PP1]
    double PP2A_dyad = 95.76;          // [uM] - Total dyadic PP2A
    double OA = 0;                     // [uM] - PP1/PP2A inhibitor Okadaic Acid
    double PLBtot = plb_val;                // [uM] - Total [PLB] in cytosolic units
    
    // Parameters for BAR module
    double Ligtot = 0.;      // [uM] - SET LIGAND CONCENTRATION HERE
    double LCCtotBA = 0.025;           // [uM] - [umol/L cytosol]
    double RyRtotBA = 0.135;           // [uM] - [umol/L cytosol]
    double PLBtotBA = plb_val;              // [uM] - [umol/L cytosol]
    double TnItotBA = 70;              // [uM] - [umol/L cytosol]
    double IKstotBA = 0.025;           // [uM] - [umol/L cytosol]
    double ICFTRtotBA = 0.025;         // [uM] - [umol/L cytosol]
    double PP1_PLBtot = 0.89;          // [uM] - [umol/L cytosol]
    double IKurtotBA = 0.025;          // [uM] - [umol/L cytosol] MOUSE
    double PLMtotBA = 48;			   // [uM]
    
    
    
    int sy0 = 217;
    double dt = 0.01;               // [ms]  time step size
    int fold = cycleLength / dt;    // iterations in one beat.
    double t = 0;
    double tb = 0;
    double tt1;
    int iter_max = (int) ( input2[3]*cycleLength / dt );      // iterations = simulation time [ms] / dt //
    double endS1 = ( input2[3] / 2 ) * cycleLength ;           // End of S1 simulation

    int tcout = 0;
    time_t t_start, t_end;
    double dif, dv;
    int ww = 0;
    int ll = 0;
    double I_inj;
    double t0 = 0;
    
    
    theCell = &(S->cellData[0][0]);
    y0n = theCell->y0n;
    // char stateFileName[] = "initial_WTstates.txt";
    char *stateFileName = argv[1];
    
    char OutputFolder[60] ; // = "OutputFolder";
    if ( argc > 3 ) {
        sprintf( OutputFolder, argv[3] );
    } else {
        printf( "%s", "Please specify output folder.\n" );
        printf( "%s", "Example:\n" );
        printf( "%s", "./a.out inputFileName1 inputFileName2 OutputFolder\n\n");
        
        exit(1);
    }
    
    mkdir( OutputFolder, 0777 );
    
   
    
    theCell->tu_v = 0;
    theCell->tu_cai = 0;
    theCell->tu_casr = 0;
    
    
    // Initial values
    FILE *fp;
    
    fp = fopen( stateFileName, "r");
    
    if ( fp == NULL ) {
        
        printf( "%s%s%s", "\nFile \'", stateFileName , "\' NOT found. \n \n" );
        
        exit(1);
    }
    
    cout << "Initial States" << endl;
    for ( int idy = 0; idy < sy0; idy++ ) {
        fscanf ( fp, "%lf", &(theCell->y0n[idy]) );
        cout << idy << "\t" << y0n[idy] << endl;
    }
    
    fclose( fp );
    
    
    for ( int idy = 0; idy < sy0; idy++ ) {
        S->cellData[ww][ll].y0n[idy] = y0n[idy];
    }
    
            
            S->cellData[ww][ll].p[0] = cycleLength;
            S->cellData[ww][ll].p[1] = recoveryTime;
            S->cellData[ww][ll].p[2] = variablePar;
            S->cellData[ww][ll].p[3] = CaMtotDyad;
            S->cellData[ww][ll].p[4] = BtotDyad;
            S->cellData[ww][ll].p[5] = CaMKIItotDyad;
            S->cellData[ww][ll].p[6] = CaNtotDyad;
            S->cellData[ww][ll].p[7] = PP1totDyad;
            S->cellData[ww][ll].p[8] = CaMtotSL;
            S->cellData[ww][ll].p[9] = BtotSL;
            S->cellData[ww][ll].p[10] = CaMKIItotSL;
            S->cellData[ww][ll].p[11] = CaNtotSL;
            S->cellData[ww][ll].p[12] = PP1totSL;
            S->cellData[ww][ll].p[13] = CaMtotCyt;
            S->cellData[ww][ll].p[14] = BtotCyt;
            S->cellData[ww][ll].p[15] = CaMKIItotCyt;
            S->cellData[ww][ll].p[16] = CaNtotCyt;
            S->cellData[ww][ll].p[17] = PP1totCyt;
            S->cellData[ww][ll].p[18] = LCCtotDyad;
            S->cellData[ww][ll].p[19] = RyRtot;
            S->cellData[ww][ll].p[20] = PP1_dyad;
            S->cellData[ww][ll].p[21] = PP2A_dyad;
            S->cellData[ww][ll].p[22] = OA;
            S->cellData[ww][ll].p[23] = PLBtot;
            S->cellData[ww][ll].p[24] = LCCtotSL;
            S->cellData[ww][ll].p[25] = PP1_SL;
            S->cellData[ww][ll].p[26] = Ligtot;
            S->cellData[ww][ll].p[27] = LCCtotBA;
            S->cellData[ww][ll].p[28] = RyRtotBA;
            S->cellData[ww][ll].p[29] = PLBtotBA;
            S->cellData[ww][ll].p[30] = TnItotBA;
            S->cellData[ww][ll].p[31] = IKstotBA;
            S->cellData[ww][ll].p[32] = ICFTRtotBA;
            S->cellData[ww][ll].p[33] = PP1_PLBtot;
            S->cellData[ww][ll].p[34] = IKurtotBA;
            S->cellData[ww][ll].p[35] = PLMtotBA;
            S->cellData[ww][ll].p[36] = CKIIOE;
    
            
    
    char filename[90];
    
    sprintf( filename, "%s/%s", OutputFolder, "vm_1Hz.txt" );
    
    FILE *fpy = fopen( filename , "w" ); //
    
    sprintf( filename, "%s/%s", OutputFolder, "allresult_1Hz.txt" );
    
    FILE *fpall = fopen( filename, "w" );
    
    
    
   	
    
    int syid1 = 9; // size of yid1
    
    // index of variables in matlab that actually used in this program.
    int yid1[] = {  30, 31,
        32, 33, 34, 36, 37,
        38,	39};
    
    int yid;
    
    
    time_t startTime;
    time_t previousTime;
    
    const time_t timeSave = 1*60*60;
    startTime = time(NULL);
    previousTime = startTime;
    
    int beat = 0;
    
    
    S->t = 0;           // Initial time.
    
    
    for ( int iter = 1; iter <= iter_max; iter ++ ) {
        
        tcout++;
        
        t = iter * dt;
        
        t0 = 0;
        tt1 = t - t0;
        tb = tt1 - cycleLength * floor( tt1 / cycleLength );
        
        
        
        if ( tb <= 5  ) {
            I_inj = -9.5;
            
        } else {
            I_inj=0.0;
        }
        
        theCell = &(S->cellData[0][ll]);
        pars1 = &(theCell->pars1);
        y0n = theCell->y0n;
        p = theCell->p;
        
        theCell->v_old = y0n[38];
        theCell->ryr_old = y0n[14];
        theCell->cai_old = y0n[37];
        theCell->casr_old = y0n[30];// + y0n[29];
        
        
        
        integrate_rk2( morotti_et_al_mouse_masterODEfile, &t, sy0, y0n, dt, p, pars1, allDrugs ) ;
        theCell->I_Total = (theCell->v_old - y0n[38]) / dt + I_inj;
        theCell->v = y0n[38];
        
        
        
        theCell = &(S->cellData[ww][ll]);
        
        dv=dt*(-theCell->I_Total);
        
        theCell->y0n[38] = theCell->v_old + dv;
        
        
        
        S->t = t;
        theCell = &(S->cellData[0][0]);
        
        y0n = theCell->y0n;
        double voltage1=y0n[38];
        double cai1=y0n[37];
        double SRCa1=y0n[30];// + y0n[29];
        
        
        if ( (tcout % 100) == 0  ) {
            fprintf( fpy, "%8.6e\t", t );
            
            for( int idy = 0; idy < syid1; idy++ ) {
                yid = yid1[idy] - 1; // index in C starts from 0, index in Matlab starts from 1.
                fprintf( fpy, "%8.6e\t", y0n[yid] );
            }
            fprintf( fpy, "\n" );
        
            
            
            
            
            fprintf( fpall, "%16.14f\t%16.14f\t%16.14f\t%16.14f\t%16.14f\t%16.14f\t%16.14f\t%16.14f\t%16.14f\t%16.14f\t%16.14f\n",
                    S->t,
                    theCell->pars1.I_Na_store,
                    theCell->pars1.I_Ca_store,
                    theCell->pars1.I_kur1_store,
                    theCell->pars1.I_kur2_store,
                    theCell->pars1.I_ss_store,
                    theCell->pars1.Jserca,
                    theCell->pars1.Jleak[0],
                    theCell->pars1.Jleak[1],
                    theCell->pars1.Nai,
                    theCell->pars1.Incx );
            
           
        }
        
        S->t = t;
        
        if ( (tcout % fold ) == 0 ) {
            cout <<"Results: "<< endl;
            cout << "t: " << S->t<< ", beat: " << beat << ", runtime: " << (time(NULL) - startTime)/60 << " min " << endl;
            beat ++;
        }
        
        tb = 0;
    }
    
   
    fclose (fpy);
    fclose( fpall );
    
    FILE * frest =  fopen( "initial.txt", "w" );
    for ( int iter = 0; iter < 217; iter++ ) {
        fprintf( frest, "%16.14e\t", y0n[iter] );
    }
    fclose( frest );
    
    //}
}
// end

