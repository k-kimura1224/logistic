#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>
#include <math.h>

#include "svpsolver.h"
#include "cut_pool.h"
#include "probdata.h"
#include "vector.h"
#include "cut.h"
//#include "solution.h"
//#include "solution_pool.h"
//#include "node.h"

using namespace std;

#define debug	1
#define log		1

void	SVPsolver::gene_OAcuts( 
	double	*u,
	double 	*l,
	double 	*Q,
	double 	M
	)
{
	assert( u != NULL );
	assert( l != NULL );
	assert( Q != NULL );
	assert( M > 0.0 );
	assert( oa_cpool.get_size() <= 0 );
	assert( probdata.get_m() > 0 );

	int n = probdata.get_m();

	bool *z = new bool[n];

	for(int i=0; i<n; i++){
		if( u[i] == 0.0 && l[i] == 0.0 ){
			z[i] = false;
		}else{
			z[i] = true;
		}
	}

	int k = count( n, z);

#if debug && 0
	cout << "*gene_OAcuts: k = " << k << endl;
#endif

	double *subQ = new double[k*k];
	int ct = 0;

	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			if( z[i] == true && z[j] == true ){
				subQ[ct] = Q[(i*n)+j];
				ct++;
			}
		}
	}

	assert( ct == k*k );

	double *V = new double[k*k]; // eigenvectors
	double *lam = new double[k]; // eigenvalues
	Gen_ZeroVec( k, lam);

	int r;
	
	r = Com_eigen( subQ, k, lam, V);
		
	if( r != 0 ){
		cout << "error: cut_oa.cpp" << endl;
	}

	/* 
		cut_i (i=1,..,k) :
			- sqrt(M/lam[i]) <= V[i*k ~ (i*k)+k ]'x <= sqrt(M/lam[i])
	*/
	// generate the pool {
	oa_cpool.alloc( k );

	CUT_PLANE cut;
	double *coef = new double[n];
	double Mlam;

	for(int i=0; i<k; i++){
		ct = 0;
		for(int j=0; j<n; j++){
			if( z[j] == true ){
				coef[j] = V[(i*k)+ct];
				ct++;
			}else{
				coef[j] = 0.0;
			}
		}
		assert( ct == k );

		assert( lam[i] > 0 );
		Mlam = sqrt( M/lam[i] );
		
		cut.set_cut( n, coef, -Mlam, Mlam);
		oa_cpool.add_cut( &cut );
	}

	// } generate cut


	delete[] z;
	delete[] subQ;
	delete[] V;
	delete[] lam;
	delete[] coef;
}
