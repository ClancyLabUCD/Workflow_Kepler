/*
 *  soltis_biophysJ2010_masterCompute.h
 *
 *  ***** Original code in Matlab ***
 *  ***** Coverted to C code by Mao-Tsuen Jeng
 *
 *  Created on 12/22/11.
 *  Copyright 2011 __UCDavis__. All rights reserved.
 %
 %
 %
 % function yfinal = soltis_biophysJ2010_masterCompute
 % This function calls the ode solver and plots results.
 % Dependency: soltis_biophysJ2010_masterODEfile.m
 %
 % Re-implemented by Anthony Soltis <ars7h@virginia.edu> for CaMKII/PKA
 % regulation of ECC model
 %
 % Author: Jeff Saucerman <jsaucerman@virginia.edu>
 % Copyright 2008, University of Virginia, All Rights Reserved
 %
 % Reference: JJ Saucerman and DM Bers, Calmodulin mediates differential
 % sensitivity of CaMKII and calcineurin to local Ca2+ in cardiac myocytes.
 % Biophys J. 2008 Aug 8. [Epub ahead of print]
 % Please cite the above paper when using this model.
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
	double I_Ca_store, I_to_store[3], I_Na_store, IKr_store, I_K1_store, ibar_store, gates[2];
	double Jserca, IKs_store, Jleak[2], ICFTR, Incx;
	
	double JCaCyt, JCaSL, JCaDyad;
} pars_rec;

typedef struct {
	double y0n[207];
	pars_rec pars1;
	double p[36];
	double v, v_old, I_Total;
    double dvdt;
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
#include "soltis_biophysJ2010_masterODEfile.h"


int main( int argc, char *argv[]) {
	
    cout.precision(16);
    
	char name[30];
	
	
	Cell * theCell;
	pars_rec * pars1; // = &(theCell.pars1);
	double * y0n; // = theCell.y0n;
	
	double * p;  // = theCell.p;
    
	
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
    
    //// Parameters for external modules
	// ECC and CaM modules
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
	// char expression[] = "WT";
	double CKIIOE = 0; // Should be zero during 'WT' and 'KO' runs
	
	if ( strcmp( expression, "OE" ) == 0 ) {
    //if ( expression == "OE" ) {
        int CKIIOE = 1; // Flag for CKII OE in ECC file (0=WT, 1=OE) - for Ito and INa
		CaMKIItotDyad = 120*6;        // [uM]
		CaMKIItotSL = 120*8.293e-4*6; // [uM]
		CaMKIItotCyt = 120*8.293e-4*6;// [uM]
	}	else if ( strcmp( expression, "KO" ) == 0 ) {
		CaMKIItotDyad = 0;          // [uM]
		CaMKIItotSL = 0;            // [uM]
		CaMKIItotCyt = 0;           // [uM]
	}
	// end
	
	// For Recovery from inactivation of LCC
	double recoveryTime = 10;  // initialize to smallest value
	
	// Parameters for CaMKII module
	double LCCtotDyad = 31.4*.9;       // [uM] - Total Dyadic [LCC] - (umol/l dyad)
	double LCCtotSL = 0.0846;          // [uM] - Total Subsarcolemmal [LCC] (umol/l sl)
	double RyRtot = 382.6;             // [uM] - Total RyR (in Dyad)
	double PP1_dyad = 95.7;            // [uM] - Total dyadic [PP1]
	double PP1_SL = 0.57;              // [uM] - Total Subsarcolemmal [PP1]
	double PP2A_dyad = 95.76;          // [uM] - Total dyadic PP2A
	double OA = 0;                     // [uM] - PP1/PP2A inhibitor Okadaic Acid
	double PLBtot = 38;                // [uM] - Total [PLB] in cytosolic units
	
	// Parameters for BAR module
	double Ligtot = input2[4];                 // [uM] - SET LIGAND CONCENTRATION HERE
	double LCCtotBA = 0.025;           // [uM] - [umol/L cytosol]
	double RyRtotBA = 0.135;           // [uM] - [umol/L cytosol]
	double PLBtotBA = 38;              // [uM] - [umol/L cytosol]
	double TnItotBA = 70;              // [uM] - [umol/L cytosol]
	double IKstotBA = 0.025;           // [uM] - [umol/L cytosol]
	double ICFTRtotBA = 0.025;         // [uM] - [umol/L cytosol]
	double PP1_PLBtot = 0.89;          // [uM] - [umol/L cytosol]
	double PLMtotBA = 48;
	
	
    
    
    
    
	int sy0 = 207;                  // size of variables for integration.
	double dt = 0.01;               // [ms]  time step size
	int fold = cycleLength / dt;    // iterations in one beat.
	double t = 0;                   // time from begining of simulation
	double tb = 0;                  // time from begining of beat.
    double tt1 = 0;                 // time in current stage ( S0, S1, ... ) of simulation
	int iter_max = (int) ( input2[3]*cycleLength / dt );      // iterations = simulation time [ms] / dt //
    double endS1 = ( input2[3] / 2 ) * cycleLength ;           // End of S1 simulation
	
	int tcout = 0;                  // output time
	time_t t_start, t_end;
	double dif, dv;
	int ww = 0;
	int ll = 0;
	double I_inj;   // inject current
	double t0 = 0;  // starting time of current simulation: S0, S1, ....
	double CL;      // Cycle Length
	
	double tt0,  tt2, highestV, lowestV, dvdt;
    
    theCell = &(S->cellData[0][0]);
	y0n = theCell->y0n;
    
    // char stateFileName[] = "initial_WTstates.txt";
    char *stateFileName = argv[1];
    
    char OutputFolder[60] ; // = "OutputFolder";
    if ( argc > 3 ) {
        sprintf( OutputFolder, "%s", argv[3] );
    } else {
        printf( "%s", "Please specify output folder.\n" );
        printf( "%s", "Example:\n" );
        printf( "%s", "./a.out inputFileName1 inputFileName2 OutputFolder\n\n");
        
        exit(1);
    }
    
    mkdir( OutputFolder, 0777 );

    
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
    S->cellData[ww][ll].p[2] = CaMtotDyad;
    S->cellData[ww][ll].p[3] = BtotDyad;
    S->cellData[ww][ll].p[4] = CaMKIItotDyad;
    S->cellData[ww][ll].p[5] = CaNtotDyad;
    S->cellData[ww][ll].p[6] = PP1totDyad;
    S->cellData[ww][ll].p[7] = CaMtotSL;
    S->cellData[ww][ll].p[8] = BtotSL;
    S->cellData[ww][ll].p[9] = CaMKIItotSL;
    S->cellData[ww][ll].p[10] = CaNtotSL;
    S->cellData[ww][ll].p[11] = PP1totSL;
    S->cellData[ww][ll].p[12] = CaMtotCyt;
    S->cellData[ww][ll].p[13] = BtotCyt;
    S->cellData[ww][ll].p[14] = CaMKIItotCyt;
    S->cellData[ww][ll].p[15] = CaNtotCyt;
    S->cellData[ww][ll].p[16] = PP1totCyt;
    S->cellData[ww][ll].p[17] = LCCtotDyad;
    S->cellData[ww][ll].p[18] = RyRtot;
    S->cellData[ww][ll].p[19] = PP1_dyad;
    S->cellData[ww][ll].p[20] = PP2A_dyad;
    S->cellData[ww][ll].p[21] = OA;
    S->cellData[ww][ll].p[22] = PLBtot;
    S->cellData[ww][ll].p[23] = LCCtotSL;
    S->cellData[ww][ll].p[24] = PP1_SL;
    S->cellData[ww][ll].p[25] = Ligtot;
    S->cellData[ww][ll].p[26] = LCCtotBA;
    S->cellData[ww][ll].p[27] = RyRtotBA;
    S->cellData[ww][ll].p[28] = PLBtotBA;
    S->cellData[ww][ll].p[29] = TnItotBA;
    S->cellData[ww][ll].p[30] = IKstotBA;
    S->cellData[ww][ll].p[31] = ICFTRtotBA;
    S->cellData[ww][ll].p[32] = PP1_PLBtot;
    S->cellData[ww][ll].p[33] = CKIIOE;
    S->cellData[ww][ll].p[34] = PLMtotBA;
    
   S->cellData[ww][ll].p[35] = 0; //Ito factor for cell type
    
    
    char filename[90];
    
    sprintf( filename, "%s/%s", OutputFolder, "vm_1Hz.txt" );
    
	FILE *fpy = fopen( filename , "w" ); //

    sprintf( filename, "%s/%s", OutputFolder, "allresult_1Hz.txt" );

	FILE *fpall = fopen( filename, "w" );
    
    sprintf( filename, "%s/%s", OutputFolder, "apds_1Hz.txt" );

	FILE *fpytwo = fopen( filename, "w" );
	
	
	int syid1 = 9; // size of yid1
    
    // index of variables in matlab that actually used in this program.
	int yid1[] = {  30, 31,
		32, 33, 34, 36, 37,
		38,	39};
    
	int yid;
    
    time_t startTime;
	time_t previousTime;
	
	const time_t timeSave = 1*60*60;  // Time to save the simuation
	startTime = time(NULL);
	previousTime = startTime;
	
    srand48(startTime);
    
    int beat = 0;
    tt0 = 6E20;  // Initial value for starting time of Apd90
    tt2 = 6E20;  // Initial value for end time of Apd90
    
    highestV = -86;     // Initial value for highest voltage
    lowestV = 0;        // Initial value for lowest volatge
    dvdt = 0;           // Initial value of dv/dt
    
	S->t = 0;           // Initial time.
    
	for ( int iter = 1; iter <= iter_max; iter ++ ) {
		
		
		if ( t < endS1) {
            
            CL = cycleLength;
            t0 = 0;
            
        } else {
            
			CL = 1e3;
			t0 = endS1;
		}
		
		tt1 = t - t0; // time in current stage
        
		tb = tt1 - CL * floor( tt1 / CL ); // time in current beat
		
		tcout++;
		
		t = iter * dt;  // total time
        
        theCell = &(S->cellData[0][ll]);
        pars1 = &(theCell->pars1);
        y0n = theCell->y0n;
        p = theCell->p;
        
        theCell->v_old = y0n[38];
        
        
        // Inject current
        if ( tb <= 0.5  ) { // stimilus duration [ms]
            I_inj = -80.0;
        } else {
            I_inj=0.0;
        }
        
        
        // integrate using Runge-Kutta order 2 method
        integrate_rk2( soltis_biophysJ2010_masterODEfile, &t, sy0, y0n, dt, p, pars1 , allDrugs ) ;
        theCell->I_Total = (theCell->v_old - y0n[38]) / dt + I_inj;
        theCell->v = y0n[38];
        
        
        dv = dt * (-theCell->I_Total);
        
        theCell->dvdt = (dv / dt);
        theCell->y0n[38] = theCell->v_old + dv;
        
		if (theCell->dvdt > dvdt){ // Find the maximum dv/dt
            dvdt = theCell->dvdt;
            tt0 = tb;  // Starting time of Apd90
        }
        
        if ( theCell->y0n[38] > highestV ){ // Find the highest voltage
            highestV = theCell->y0n[38];
            
        }
        
        if ( theCell->y0n[38] < lowestV) { // Find the lowest voltage.
            lowestV = theCell->y0n[38];
        }
        
        if ( theCell->y0n[38] > highestV - 0.9 * ( highestV - lowestV) ) {
            tt2 = tb; // End time of Apd90
        }
        
		
		if ( (tcout % 100) == 0 && t > endS1 ) {
			fprintf( fpy, "%8.6e\t", t );
			
			for( int idy = 0; idy < syid1; idy++ ) {
				yid = yid1[idy] - 1; // index in C starts from 0, index in Matlab starts from 1.
				fprintf( fpy, "%8.6e\t", y0n[yid] );
			}
			fprintf( fpy, "\n" );
			
			
			fprintf( fpall, "%16.14f\t%16.14f\t%16.14f\t%16.14f\t%16.14f\t%16.14f\t%16.14f\t%16.14f\t%16.14f\t%16.14f\t%16.14f\t%16.14f\t%16.14f\t%16.14f\t%16.14f\t%16.14f\n",
					S->t, theCell->pars1.I_Ca_store,
					theCell->pars1.I_to_store[0],
					theCell->pars1.I_to_store[1],
					theCell->pars1.I_to_store[2],
					theCell->pars1.I_Na_store,
					theCell->pars1.I_K1_store,
					theCell->pars1.gates[0],
					theCell->pars1.gates[1],
					theCell->pars1.Jserca,
					theCell->pars1.IKs_store,
                    theCell->pars1.IKr_store,
					theCell->pars1.Jleak[0],
					theCell->pars1.Jleak[1],
					theCell->pars1.ICFTR,
					theCell->pars1.Incx );
			
			
		}
		
        
		S->t = t;
		
		if ( (tcout % fold ) == 0 ) {
            cout <<"Results: "<< endl;
			cout << "t: " << S->t<< ", beat: " << beat << ", runtime: " << (time(NULL) - startTime)/60 << " min " << endl;
            beat ++;
            
            fprintf (fpytwo, "%d\t%8.6f\t%8.6f\t%8.6f\n", beat, tt2-tt0, tt0, tt2);
            
			highestV = -86;
            lowestV = 0;
            dvdt = 0;
            tt0 = 6E20;
            tt2 = 6E20;
		}
		tb = 0;
    }
	
	
	fclose (fpytwo);
	fclose( fpy );
	
	
	FILE * frest =  fopen( "finalStates.txt", "w" );
	for ( int iter = 0; iter < 207; iter++ ) {
		fprintf( frest, "%16.14e\t", y0n[iter] );
	}
	fclose( frest );
	
    fclose( fpall );
	
    
}
// end main


