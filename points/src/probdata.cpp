#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>

#include "probdata.h"
#include "vector.h"

using namespace std;

#define debug_class	0
#define debug	0

PROB_DATA::PROB_DATA(){	// default constructor
#if debug_class
	cout << "PROB_DATA: default constructor" << endl;
#endif
	m	=	-1;
	B = NULL;
	B_ = NULL;
	Q = NULL;
}

PROB_DATA::PROB_DATA( const PROB_DATA &source )
{	// copy constructor
#if debug_class
	cout << "PROB_DATA: copy constructor" << endl;
#endif
	m	=	source.m;

	if( m > 0 ){
		assert( source.B != NULL );
		assert( source.B_ != NULL );
		assert( source.Q != NULL );

		B = new double[m*m];
		B_ = new double[m*m];
		Q = new double[m*m];

		Copy_vec( source.B, B, m*m);
		Copy_vec( source.B_, B_, m*m);
		Copy_vec( source.Q, Q, m*m);
	}else{
		B = NULL;
		B_ = NULL;
		Q = NULL;
	}
}

// assignment operator
PROB_DATA& PROB_DATA::operator=( const PROB_DATA& source )
{
#if debug_class
	cout << "PROB_DATA: assignment operator" << endl;
#endif

	if( this != &source ){ 	
		m	=	source.m;

		if( m > 0 ){
			assert( source.B != NULL );
			assert( source.B_ != NULL );
			assert( source.Q != NULL );
			
			delete[] B;
			delete[] B_;
			delete[] Q;

			B = new double[m*m];
			B_ = new double[m*m];
			Q = new double[m*m];
	
			Copy_vec( source.B, B, m*m);
			Copy_vec( source.B_, B_, m*m);
			Copy_vec( source.Q, Q, m*m);
		}else{
			B = NULL;
			B_ = NULL;
			Q = NULL;
		}
	}
	return *this;
}

// destructor
PROB_DATA::~PROB_DATA()
{
#if debug_class
	cout << "PROB_DATA: destructor" << endl;
#endif
	delete[] B;
	delete[] B_;
	delete[] Q;
	B = NULL;
	B_ = NULL;
	Q = NULL;
}

void PROB_DATA::set_data(
	int		s_m,
	double	*s_B,
	double	*s_B_,
	double	*s_Q
	)
{
	assert( m == -1 );
	assert( s_m > 0 );
	assert( s_B != NULL );
	assert( s_B_ != NULL );
	assert( s_Q != NULL );

	m = s_m;

	delete[] B;
	delete[] B_;
	delete[] Q;

	B = new double[m*m];
	B_ = new double[m*m];
	Q = new double[m*m];
	
	Copy_vec( s_B, B, m*m);
	Copy_vec( s_B_, B_, m*m);
	Copy_vec( s_Q, Q, m*m);
}

