#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>

#include "cut.h"
#include "vector.h"

using namespace std;

#define debug_class	0
#define debug	0

CUT_PLANE::CUT_PLANE(){	// default constructor
#if debug_class
	cout << "CUT_PLANE: default constructor" << endl;
#endif
	m = -1;
	coef = NULL;
	lb = NULL;
	ub = NULL;
}

CUT_PLANE::CUT_PLANE( const CUT_PLANE &source )
{	// copy constructor
#if debug_class
	cout << "CUT_PLANE: copy constructor" << endl;
#endif
	m = source.m;

	if( source.coef != NULL ){
		assert( m > 0 );
		coef = new double[m];
		Copy_vec( source.coef, coef, m);
	}else{
		coef = NULL;
	}

	if( source.lb != NULL ){
		assert( m > 0 );
		lb = new double;
		*lb = *(source.lb);
	}else{
		lb = NULL;
	}

	if( source.ub != NULL ){
		assert( m > 0 );
		ub = new double;
		*ub = *(source.ub);
	}else{
		ub = NULL;
	}

}

// assignment operator
CUT_PLANE& CUT_PLANE::operator=( const CUT_PLANE& source )
{
#if debug_class
	cout << "CUT_PLANE: assignment operator" << endl;
#endif

	if( this != &source ){ 	
		m = source.m;
	
		if( source.coef != NULL ){
			assert( m > 0 );
			delete[] coef;
			coef = new double[m];
			Copy_vec( source.coef, coef, m);
		}else{
			coef = NULL;
		}
	
		if( source.lb != NULL ){
			assert( m > 0 );
			delete[] lb;
			lb = new double;
			*lb = *(source.lb);
		}else{
			lb = NULL;
		}
	
		if( source.ub != NULL ){
			assert( m > 0 );
			delete[] ub;
			ub = new double;
			*ub = *(source.ub);
		}else{
			ub = NULL;
		}

	}

	return *this;
}

// destructor
CUT_PLANE::~CUT_PLANE()
{
#if debug_class
	cout << "CUT_PLANE: destructor" << endl;
#endif
	delete[] coef;
	delete[] lb;
	delete[] ub;
	coef = NULL;
	lb = NULL;
	ub = NULL;
}

void CUT_PLANE::set_cut(
	int		s_m,
	double	*s_coef,
	double	s_lb,
	double	s_ub
	)
{
	assert( s_m > 0 );
	assert( s_coef != NULL );
	assert( s_ub > s_lb );

	m = s_m;

	coef = new double[m];
	lb = new double;
	ub = new double;

	Copy_vec( s_coef, coef, m);
	*lb = s_lb;
	*ub = s_ub;

}
