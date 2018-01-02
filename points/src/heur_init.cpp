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
#include "probdata.h"
#include "vector.h"
#include "solution.h"
#include "solution_pool.h"

using namespace std;

#define debug	0

void	SVPsolver::find_min_column()
{
	int		m;
	double 	*B_ = NULL;
	double	*nrm = NULL;
	double	min = 0.0;
	int		memo = -1;

	m = probdata.get_m();
	B_ = probdata.get_B_();
	nrm = new double[m];

	assert( m > 0 );
	assert( B_ != NULL );
	
	// compute norms
	for(int i=0; i<m; i++){
		nrm[i] = Com_nrm( &B_[i*m], m);
		//cout << i << ":" << nrm[i] << endl;
		if( i==0 ){
			min = nrm[i];
			memo = i;
		}else if( min > nrm[i] ){
			min = nrm[i];
			memo = i;
		}
	}

	assert( min > 0.0 );
	assert( memo >= 0 && memo < m );

	// copy
	assert( norm == NULL );
	norm = new double[m];
	Copy_vec( nrm, norm, m);

	// sort	
	if( BRANCHINGRULE_INT == 4 ){
		assert( order == NULL );
		order = new int[m];
		double *_null = NULL;
		Bubble_Sort( m, nrm, _null, order);
		printv( m, order);
	}

	// set
	for(int i=0; i<m; i++){
		if( i == memo )	nrm[i] = 1.0;
		else					nrm[i] = 0.0;
	}

	bestval = min*min;
	bestsol.set_sol( m, nrm, min*min);
	pool.add_solution( bestsol );

	assert( _Appfac != 0.0 );
	Appfac = sqrt( bestval ) / _Appfac;

	delete[] nrm;
}
