/*
 *  soltis_biophysJ2010_masterCompute.h
 *  
 *
 *  Created by Mao-Tsuen Jeng on 12/22/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>



/*
 function yfinal = soltis_biophysJ2010_masterCompute
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

void masterCompute_1D( int argc, char *argv[] );

void masterCompute_1D( int argc, char *argv[] ) {
    
    char *InputFile2 = argv[2];
    double input2[5];
    
	FILE *fp2;
    char tissueType[20];
    int dimSize[2];
    char expression[2];
    char dim[2];

    
    fp2 = fopen( InputFile2, "r");
    
    if ( fp2 == NULL ) {
        
        printf( "%s%s%s", "\nFile \'", InputFile2 , "\' NOT found. \n \n" );
        
        exit(1);
    }
    
    try
    {
    
    for ( int idy = 0; idy < 5; idy++ ) {
        fscanf ( fp2, "%lf", &input2[idy] );
    }
    
    fscanf ( fp2, "%s", expression );
    
    fscanf( fp2, "%s", dim );
    
    fscanf( fp2, "%d", &dimSize[0] );   // tw
    fscanf( fp2, "%d", &dimSize[1] );   // tl
 
    fscanf( fp2, "%s", tissueType );
    }
    catch (exception& e)
    {
        cout << e.what() << '\n';
        fclose( fp2 );
        exit(1);
    }
    
    
    
    fclose( fp2 );
    
    if ( strcmp( tissueType , "homogeneous")!=0 && strcmp( tissueType , "heterogeneous")!=0 ) {
        printf( "%s", "Please enter 'homogeneous' or 'heterogeneous'.\n" );
        exit(1);
    }

    
    double allDrugs[2];
    allDrugs[0] = input2[0];
	allDrugs[1] = input2[1];
    
    
	
	char name[30];
    
    SimState theState;
    SimState *S = &theState;

    dimSize[0] = 1;
    const int tw = dimSize[0];
    const int tl = dimSize[1];
    Cell cellData[tw][tl];

    for ( int i = 0; i < tw ; i++ ) {
        S->cellData[i] = cellData[i];
    }
	
	Cell * theCell;
	pars_rec * pars1;
	double * y0n;
	S->tStep = 0;
	double * p;
	
    
	double Df = 0.00154;
	double dx = 0.01;
	double dy = 0.01;
	
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
	
	if ( strcmp( expression , "OE" ) ) {
		int CKIIOE = 1; // Flag for CKII OE in ECC file (0=WT, 1=OE) - for Ito and INa
		CaMKIItotDyad = 120*6;        // [uM] 
		CaMKIItotSL = 120*8.293e-4*6; // [uM]
		CaMKIItotCyt = 120*8.293e-4*6;// [uM]
	}	else if ( strcmp( expression , "KO" ) ) {
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
	

	// char stateFileName[] = "initial_WTstates.txt"; // from rest pacing
		
	int sy0 = 207; // size of variables for integration.
	
	double dt = 0.01; // [ms]  time step size
	int fold = 1 / dt;
	double t = 0;   // time from begining of simulation
	double tb = 0;  // time from begining of beat.
	int iter_max = (int) ( input2[3] * cycleLength / dt ); // simulation time [ms]/ dt;200e3/dt; //
	
	
	int tcout = 0;
	time_t t_start, t_end;
	double dif, dv;
	int ww = 0;
	int ll;
	double I_inj;
	double t0 = 0;
	double CL;
	double tt1;
	double Ev;
	
	theCell = &(S->cellData[0][0]);

    
	y0n = theCell->y0n;
    
    
    char *stateFileName = argv[1];
    
    char OutputFolder[600] ; // = "OutputFolder";
    if ( argc > 3 ) {
        sprintf( OutputFolder, argv[3] );
    } else {
        printf( "%s", "Please specify output folder.\n" );
        printf( "%s", "Example:\n" );
        printf( "%s", "./a.out inputFileName1 inputFileName2 OutputFolder\n\n");
        
        exit(1);
    }
    
    mkdir( OutputFolder, 0777 );
    
    
    
	FILE *fp = fopen( stateFileName, "r");
    
    if ( fp == NULL ) {
        
        printf( "%s%s%s", "\nFile \'", stateFileName , "\' NOT found. \n \n" );
        
        exit(1);
    }
    
    
	for ( int idy = 0; idy < sy0; idy++ ) {
		fscanf ( fp, "%lf", &y0n[idy] );
		
	}
	fclose( fp );
	
    
	
	ww = 0;
	for ( ll = 0; ll < tl; ll++ ) {
		if ( ll != 0 ) {
			for ( int idy = 0; idy < sy0; idy++ ) {
				S->cellData[ww][ll].y0n[idy] = y0n[idy];
			}
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
		
		S->cellData[ww][ll].p[35] = (ll-tl) * 3e-04;
		// cout << tissueType << endl;

		if ( strcmp( tissueType , "homogeneous" ) == 0 ) {
            S->cellData[ww][ll].p[35] = 0;
        }

	}
	//	}
	
    char filename[690];
    
    sprintf( filename, "%s/%s", OutputFolder, "y_1D.txt" );
    
	FILE *fpy = fopen( filename, "w" );
	
    
    sprintf( filename, "%s/%s", OutputFolder, "ECGs.txt" );
    
	FILE *ecg = fopen( filename, "w");
	

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
		
	
	
	for ( int iter = 1; iter <= iter_max; iter ++ ) {
		

            
            CL = cycleLength;
            t0 = 0;
		
		tt1 = t - t0;
		tb = tt1 - CL * floor( tt1 / CL );
		
		tcout++;
		
		t = iter * dt;  // total time
		
		
		
#pragma omp parallel for private( ll, theCell, I_inj, pars1, y0n, p)			
		
		
		for ( ll = 0; ll < tl; ll++ ) {
			
			if ( tb < 1.0 && ( ll== 0  )  ) {
				I_inj = -500.0;
			            
			} else {
				I_inj=0.0;
			}
			
			
			theCell = &(S->cellData[0][ll]);
			pars1 = &(theCell->pars1);
			y0n = theCell->y0n;
			p = theCell->p;
			
			theCell->v_old = y0n[38];
			integrate_rk2( soltis_biophysJ2010_masterODEfile, &t, sy0, y0n, dt, p, pars1, allDrugs ) ;
			theCell->I_Total = (theCell->v_old - y0n[38]) / dt + I_inj;
			theCell->v = y0n[38];
		}
		
		
		
#pragma omp parallel for private( ll, theCell, dv) 		
		
		
		for (ll=0; ll<=(tl-1); ll+=1){
			
			theCell = &(S->cellData[0][ll]);
			
			
			if(ll>0 && ll<(tl-1)) { dv=dt*(-theCell->I_Total+Df*(S->cellData[ww][ll-1].v - 2*theCell->v + S->cellData[ww][ll+1].v )/(dx*dx) ); }
			else if (ll==0) { dv=dt*(-theCell->I_Total+Df*(-theCell->v + S->cellData[ww][ll+1].v )/(dx*dx)); }
			else if (ll==(tl-1)) { dv=dt*(-theCell->I_Total+Df*(-theCell->v + S->cellData[ww][ll-1].v )/(dx*dx)); }
			
			theCell->y0n[38] = theCell->v_old + dv;
			
		}
		
		
		S->t = t;
		
		if ( (tcout % 100 ) == 0 ) {
			
			cout << "t: " << t << ", tStep: " << S->tStep << ", runtime: " << (time(NULL) - startTime)/60 << " min " << endl;
		}
		
		if ( (tcout % 100 ) == 0 ) {
			
                fprintf (fpy, "%8.6e\t", t);
            
			for (ll=0; ll<=(tl-1); ll+=1){
				
				
				
				fprintf(fpy, "%8.6f\t", S->cellData[ww][ll].y0n[38]);
			}
			fprintf(fpy, "\n");
			
			
			Ev=0;
			for (ll=0; ll<=(tl-2); ll++){
				Ev = Ev + ( S->cellData[ww][ll].y0n[38] - S->cellData[ww][ll+1].y0n[38] )*( 1 / ((tl - ll - 1)*dx+2)-1/((tl - ll)*dx+2));
			}
			fprintf(ecg, "%10.3f\t%10.6f\n", t, Ev);
		}
		tb = 0;
		
		
	}
	
	
	

	
	FILE * frest =  fopen( "fiber.cells", "w" );
    
    for (ll=0; ll<=(tl-1); ll+=1){
        
        theCell = &(S->cellData[0][ll]);
        fwrite( theCell, sizeof(Cell), 1, frest);
    }
    
	fclose( frest );
	
	fclose( fpy );
	fclose(ecg);
	
}
// end
