/*
 *  soltis_biophysJ2010_masterCompute.h
 *
 *
 *  Created by Mao-Tsuen Jeng on 12/22/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */


void masterCompute_2D( int argc, char *argv[] );

void masterCompute_2D( int argc, char *argv[] ) {
    
    char *InputFile2 = argv[2];
    double input2[5];
    
	FILE *fp2;
    char tissueType[20];
    int dimSize[2];
    char expression[2];
    char dim[2];
    double s2time;
    
    
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
    
    //    cout << tissueType << endl;
    //	 cout << strcmp( tissueType, "homogeneous" ) << endl;
    //	 cout << strcmp( tissueType , "heterogeneous" ) << endl;
    //	 cout << ( strcmp( tissueType , "homogeneous" ) != 0 && strcmp( tissueType , "heterogeneous" ) != 0 ) << endl;
    
//    if ( strcmp( tissueType , "homogeneous" ) != 0 && strcmp( tissueType , "heterogeneous" ) != 0 ) {
//        printf( "%s", "Please enter 'homogeneous' or 'heterogeneous'.\n" );
//        exit(1);
//    }
//    
    
    double allDrugs[2];
    allDrugs[0] = input2[0];
	allDrugs[1] = input2[1];
    
    
    
	char name[830];
    
    SimState theState;
    SimState *S = &theState;
    
    const int tw = dimSize[0];
    const int tl = dimSize[1];
    
    // cout << tw << "\t" << tl << "\t" << tw*tl*sizeof(Cell) << endl;
    
    // Cell cellData[tw][tl];
    Cell * cellData[tw];
    for ( int i = 0; i < tw; i++ ) {
        cellData[i] = (Cell*) calloc (tl,sizeof(Cell));
    }
    
    // cout << tw << "\t" << tl << endl;
    
    
    for ( int i = 0; i < tw ; i++ ) {
        S->cellData[i] = cellData[i];
    }
    
	Cell * theCell;
	pars_rec * pars1;
	double * y0n;
	S->tStep = 0;
	double * p;
    
    
    
	
	const double Dfx = 0.003/4; //tl direction
    const double Dfy = 0.00154; //tw direction
	const double dx = 0.01;
	const double dy = 0.01;
    
	// omp_set_num_threads( 20 );
    
	
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
    
    double CKIIOE = 0; // Should be zero during 'WT' and 'KO' runs
    
    // if ( expression == "OE" ) {
    if ( strcmp( expression, "OE" ) == 0 ) {
        int CKIIOE = 1; // Flag for CKII OE in ECC file (0=WT, 1=OE) - for Ito and INa
        CaMKIItotDyad = 120*6;        // [uM]
        CaMKIItotSL = 120*8.293e-4*6; // [uM]
        CaMKIItotCyt = 120*8.293e-4*6;// [uM]
    // }	else if ( expression == "KO" ) {
    } else if ( strcmp( expression, "KO" ) == 0 ) {
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
    
	
	
    char stateFileName[] = "fiber.cells"; // from fiber steady pacing
    // char *stateFileName = argv[1];
    
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
    
    
    
	int sy0 = 217;
	
	double dt = 0.01;
	int fold = 1 / dt;
	double t = 0;
	double tb = 0;
	int iter_max = (int) ( cycleLength / dt ) ; // simulation time [ms]/ dt;
    
	
	
	
	time_t t_start, t_end;
	double dif, dv;
	int ww = 0;
	int ll;
	double I_inj;
	double t0 = 0;
	double CL;
	double tt1;
	double Ev;
	
    
    char filename[690];
    
    
    FILE *fpytwo;
    
    
        
    
    const char *twoDstart = "twoDstart.SimState";
    
    
    // FILE *fp3 = fopen(twoDstart, "r");
    
    
    // if (fp3 == NULL) {
        cout << "From fiber initial file" << endl;
        
        FILE *fp = fopen( stateFileName, "r");
        
        if ( fp == NULL ){
            cout << "Fiber initial file not found." << endl;
            
            exit(1);
        } else {
            for ( ww = 0; ww < tw; ww++ ) {
                
                for ( ll = 0; ll < tl; ll++ ) {
                    
                    theCell = &(S->cellData[ww][ll]);
                    
                    fread( theCell, sizeof(Cell), 1, fp);
                }
                fclose( fp );
                fp = fopen( stateFileName, "r");
                
            }
            S->counter = 0;
            
        }
        
        fclose( fp );
        
        S->counter = 0;
    // }
//     else {
//        // fread(S, sizeof(SimState), 1, fp3);
//        fread(&(S->t), sizeof(double), 1, fp3);
//        fread(&(S->tt), sizeof(double), 1, fp3);
//        fread(&(S->tStep), sizeof(int), 1, fp3);
//        fread(&(S->counter), sizeof(int), 1, fp3);
//        fread(&(S->beat), sizeof(int), 1, fp3);
//        for ( ww = 0; ww < tw; ww++ ) {
//            for ( ll = 0; ll < tl; ll++ ) {
//                theCell = &(S->cellData[ww][ll]);
//                fread( theCell, sizeof( Cell ), 1, fp3 );
//            }
//        }
//        cout << "from 2D states " << endl;
//        fclose( fp3 );
//    }
    
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
        
        tt1 = t - t0;
        tb = tt1 - CL * floor( tt1 / CL );
        
        
        S->counter ++;
        t = iter * dt;  // total time
        
        
        
        
#pragma omp parallel for private( ww, ll, theCell, I_inj, pars1, y0n, p)
        for ( ww = 0; ww < tw; ww++ ) {
            
            for ( ll = 0; ll < tl; ll++ ) {
                
                if ( tb < 0.5 && ( ll==0 || ll==1 ) ) {
                    I_inj = -500.0;
                    
                } else {
                    I_inj=0.0;
                }
                
                theCell = &(S->cellData[ww][ll]);
                pars1 = &(theCell->pars1);
                y0n = theCell->y0n;
                p = theCell->p;
                
                theCell->v_old = y0n[38];
                integrate_rk2( morotti_et_al_mouse_masterODEfile, &t, sy0, y0n, dt, p, pars1, allDrugs ) ;
                theCell->I_Total = (theCell->v_old - y0n[38]) / dt + I_inj;
                theCell->v = theCell->v_old;
            }
        }
        
        
#pragma omp parallel for private( ww, ll, theCell, dv)
        
        
        for ( ww = 0; ww < tw; ww++ ) {
            for ( ll = 0; ll < tl; ll++ ) {
                
                theCell = &(S->cellData[ww][ll]);
                
                if(ww>0 && ww<(tw-1) && ll>0 && ll<(tl-1)) {
                    dv=dt*(-theCell->I_Total+Dfx*(S->cellData[ww][ll-1].v-2*theCell->v+S->cellData[ww][ll+1].v)/(dx*dx)+Dfy*(S->cellData[ww-1][ll].v-2*theCell->v+S->cellData[ww+1][ll].v)/(dy*dy) );
                }
                else if (ww==0 && ll==0) {
                    dv=dt*(-theCell->I_Total+Dfx*(-theCell->v+S->cellData[ww][ll+1].v)/(dx*dx)+Dfy*(-theCell->v+S->cellData[ww+1][ll].v)/(dy*dy) );
                }
                else if (ww==0 && ll==(tl-1)) {
                    dv=dt*(-theCell->I_Total+Dfx*(-theCell->v+S->cellData[ww][ll-1].v)/(dx*dx)+Dfy*(-theCell->v+S->cellData[ww+1][ll].v)/(dy*dy) );
                }
                else if (ww==(tw-1) && ll==0) {
                    dv=dt*(-theCell->I_Total+Dfx*(-theCell->v+S->cellData[ww][ll+1].v)/(dx*dx)+Dfy*(-theCell->v+S->cellData[ww-1][ll].v)/(dy*dy) );
                }
                else if (ww==(tw-1) && ll==(tl-1)) {
                    dv=dt*(-theCell->I_Total+Dfx*(-theCell->v+S->cellData[ww][ll-1].v)/(dx*dx)+Dfy*(-theCell->v+S->cellData[ww-1][ll].v)/(dy*dy) );
                }
                else if (ww==0 && ll>0 && ll<(tl-1)) {
                    dv=dt*(-theCell->I_Total+Dfx*(S->cellData[ww][ll-1].v-2*theCell->v+S->cellData[ww][ll+1].v)/(dx*dx)+Dfy*(-theCell->v+S->cellData[ww+1][ll].v)/(dy*dy) );
                }
                else if (ww==(tw-1) && ll>0 && ll<(tl-1)) {
                    dv=dt*(-theCell->I_Total+Dfx*(S->cellData[ww][ll-1].v-2*theCell->v+S->cellData[ww][ll+1].v)/(dx*dx)+Dfy*(-theCell->v+S->cellData[ww-1][ll].v)/(dy*dy) );
                }
                else if (ww>0 && ww<(tw-1) && ll==0) {
                    dv=dt*(-theCell->I_Total+Dfx*(-theCell->v+S->cellData[ww][ll+1].v)/(dx*dx)+Dfy*(S->cellData[ww-1][ll].v-2*theCell->v+S->cellData[ww+1][ll].v)/(dy*dy) );
                }
                else if (ww>0 && ww<(tw-1) && ll==(tl-1)) {
                    dv=dt*(-theCell->I_Total+Dfx*(-theCell->v+S->cellData[ww][ll-1].v)/(dx*dx)+Dfy*(S->cellData[ww-1][ll].v-2*theCell->v+S->cellData[ww+1][ll].v)/(dy*dy) );
                }
                
                theCell->y0n[38] = theCell->v_old + dv;
            }
        }
        
        
        
        
        S->t = t;
        
        if ( (S->counter % 100 ) == 0 ) {
            
            cout << "t: " << t << ", tStep: " << S->tStep << ", runtime: " << (time(NULL) - startTime)/60 << " min " << endl;
		}
		
		if ( (S->counter % 100 ) == 0 ) {
			
            
			sprintf(name, "%s/ap%d.dat", OutputFolder, S->counter/100);
            fpytwo = fopen(name, "w");
            cout << name << endl;
			
			for (ww=0; ww < tw; ww+=10) {
                
                for (ll=0; ll< tl; ll+=10){
                    
                    theCell = &(S->cellData[ww][ll]);
                    
                    fprintf(fpytwo, "%8.6f\t", S->cellData[ww][ll].y0n[38]);
                }
                fprintf(fpytwo, "\n");
            }
            fclose (fpytwo);
            
			
            
		}
		tb = 0;
		
		
    }
    
    fp2 = fopen(twoDstart, "w");
	// fwrite(S, sizeof(SimState), 1, fp2);
    fwrite(&(S->t), sizeof(double), 1, fp2);
    fwrite(&(S->tt), sizeof(double), 1, fp2);
    fwrite(&(S->tStep), sizeof(int), 1, fp2);
    fwrite(&(S->counter), sizeof(int), 1, fp2);
    fwrite(&(S->beat), sizeof(int), 1, fp2);
    
    for ( ww = 0; ww < tw; ww++ ) {
        for ( ll = 0; ll < tl; ll++ ) {
            theCell = &(S->cellData[ww][ll]);
            fwrite( theCell, sizeof( Cell ), 1, fp2 );
        }
    }
    
    
    // Close files
	fclose(fp2);
    
    
    
    // Free memory
    for ( int i = 0; i < tw; i++ ) {
        free( cellData[i] );
    }
    
    
}
// end
