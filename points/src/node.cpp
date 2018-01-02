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

#include "node.h"
#include "vector.h"

using namespace std;

#define debug_class	0
#define debug	0

NODE::NODE(){	// default constructor
#if debug_class
	cout << "NODE: default constructor" << endl;
#endif
	m = -1;
	ub = NULL;
	lb = NULL;
	warm = NULL;
	relax_objval = - 1.0e+10;
	relax_solval = NULL;
	dpt = -1;
	zero = true;
	index = -1;
	solved = false;
	sumfixed = NULL;
}

NODE::NODE( const NODE &source )
{	// copy constructor
#if debug_class
	cout << "NODE: copy constructor" << endl;
#endif
	m = source.m;
	relax_objval = source.relax_objval;
	zero = source.zero;
	dpt = source.dpt;
	index = source.index;
	solved = source.solved;
	
	if( m > 0 ){
		assert( source.ub != NULL );
		assert( source.lb != NULL );
		assert( source.warm != NULL );
		
		ub = new double[m];
		lb = new double[m];
		warm = new double[m];
		
		Copy_vec( source.ub, ub, m);
		Copy_vec( source.lb, lb, m);
		Copy_vec( source.warm, warm, m);

		if( source.relax_solval != NULL ){
			relax_solval = new double[m];
			Copy_vec( source.relax_solval, relax_solval, m);
		}else{
			relax_solval = NULL;
		}
		if( source.sumfixed != NULL ){
			sumfixed = new double[m];
			Copy_vec( source.sumfixed, sumfixed, m);
		}else{
			sumfixed = NULL;
		}
	}else{
		ub = NULL;
		lb = NULL;
		warm = NULL;
		relax_solval = NULL;
		sumfixed = NULL;
	}

}

// assignment operator
NODE& NODE::operator=( const NODE& source )
{
#if debug_class
	cout << "NODE: assignment operator" << endl;
#endif

	if( this != &source ){ 	
		m = source.m;
		relax_objval = source.relax_objval;
		zero = source.zero;
		index = source.index;
		dpt = source.dpt;
		solved = source.solved;
		
		if( m > 0 ){
			assert( source.ub != NULL );
			assert( source.lb != NULL );
			assert( source.warm != NULL );
			
			delete[] ub;
			delete[] lb;
			delete[] warm;
			ub = new double[m];
			lb = new double[m];
			warm = new double[m];
			
			Copy_vec( source.ub, ub, m);
			Copy_vec( source.lb, lb, m);
			Copy_vec( source.warm, warm, m);

			if( source.relax_solval != NULL ){
				delete[] relax_solval;
				relax_solval = new double[m];
				Copy_vec( source.relax_solval, relax_solval, m);
			}else{
				relax_solval = NULL;
			}
			if( source.sumfixed != NULL ){
				delete[] sumfixed;
				sumfixed = new double[m];
				Copy_vec( source.sumfixed, sumfixed, m);
			}else{
				sumfixed = NULL;
			}
		}else{
			ub = NULL;
			lb = NULL;
			warm = NULL;
			relax_solval = NULL;
			sumfixed = NULL;
		}
	}

	return *this;
}

// destructor
NODE::~NODE()
{
#if debug_class
	cout << "NODE: destructor" << endl;
#endif
	delete[] ub;
	delete[] lb;
	delete[] warm;
	delete[] relax_solval;
	delete[] sumfixed;
	ub = NULL;
	lb = NULL;
	relax_solval = NULL;
	warm = NULL;
	sumfixed = NULL;
}

void NODE::set_vals(
	int		s_m,
	double	*s_ub,
	double	*s_lb,
	double	*s_warm,
	double	s_relax_objval,
	int		s_dpt,
	bool 		s_zero,
	int		s_index
	)
{
	assert( s_m > 0 );
	assert( ub == NULL );
	assert( lb == NULL );
	assert( warm == NULL );

	m = s_m;
	relax_objval = s_relax_objval;
	dpt = s_dpt;
	zero = s_zero;
	index = s_index;

	ub = new double[m];
	lb = new double[m];
	warm = new double[m];

	Copy_vec( s_ub, ub, m);
	Copy_vec( s_lb, lb, m);
	Copy_vec( s_warm, warm, m);
	
}

void NODE::set_relaxsolval(
	double	*solval
	)
{
	assert( solval != NULL );
	assert( m > 0 );

	delete[] relax_solval;
	relax_solval = new double[m];

	//Copy_vec( solval, relax_solval, m);
	double ep = 1e-10;
	for(int i=0; i<m ; i++){
		if( ceil(solval[i]) - solval[i] < ep ){
			relax_solval[i] = ceil(solval[i]);
		}else if(solval[i] - floor(solval[i]) <ep ){
			relax_solval[i] = floor(solval[i]);
		}else{
			relax_solval[i] = solval[i];
		}
	}
}

bool NODE::alloc_sumfixed()
{
	assert( m > 0 );

	bool result=true;

	if( sumfixed == NULL ){
		sumfixed = new double[m];
		result = true;
	}else{
		result = false;
	}

	return result; 
}

void NODE::set_sumfixed(
	double	c,
	double	*s_sumfixed
	)
{
	assert( s_sumfixed != NULL );
	assert( sumfixed != NULL );
	assert( m > 0 );
	assert( c != 0.0 );

	Com_scal( s_sumfixed, m, c, sumfixed);
}
	
void NODE::add_sumfixed(
	double	c,
	double	*s_sumfixed
	)
{
	assert( s_sumfixed != NULL );
	assert( sumfixed != NULL );
	assert( m > 0 );
	assert( c != 0.0 );

	for(int i=0; i<m; i++){
		sumfixed[i] += c*s_sumfixed[i];
	}
}
