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
#include <stdlib.h>
#include <cstring>

#include <sys/types.h>  // for mkdir using cigwin
#include <sys/stat.h>




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
} Cell;


typedef struct {
	double t, tt;
	int tStep, counter, beat;
	Cell * cellData[1000];
} SimState;


const double pi = 3.141592653589793; // pi = 3.1415926535897932384626433832795028841971...

#include "integrate_rk2.h"
#include "soltis_biophysJ2010_masterODEfile.h"

#include "masterCompute_1D.h"
#include "masterCompute_2D.h"




int main( int argc, char *argv[] ) {
    
    using namespace std;
    
    char *InputFile2 = argv[2];
    double input2[5];
    
	FILE *fp2;
    
    fp2 = fopen( InputFile2, "r");
    
    if ( fp2 == NULL ) {
        
        printf( "%s%s%s", "\nFile \'", InputFile2 , "\' NOT found. \n \n" );
        
        exit(1);
    }
    
    for ( int idy = 0; idy < 5; idy++ ) {
        fscanf ( fp2, "%lf", &input2[idy] );
    }
    
    char expression[2];
    fscanf ( fp2, "%s", expression );
    // exit(0);
    char dim[2];
    fscanf( fp2, "%s", dim );
    
    fclose( fp2 );
    
    if ( strcmp( dim , "1D" )==0 ) {
        cout << "Performing 1D simulation\n";
        masterCompute_1D( argc, argv );
    } else if ( strcmp( dim , "2D" )==0 ) {
        
        masterCompute_2D( argc, argv );
    } else {
        cout << "Please specify 1D or 2D.\n";
    }
    
    
}
// end
