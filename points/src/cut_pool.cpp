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
#include "cut_pool.h"

using namespace std;

#define debug	0
#define debug_class	0

CUT_POOL::CUT_POOL(){	// default constructor
#if debug_class
	cout << "CUT_POOL: default constructor" << endl;
#endif
	size = -1;
	num = -1;
	cuts = NULL;
}

CUT_POOL::CUT_POOL( const CUT_POOL &source )
{	// copy constructor
#if debug_class
	cout << "CUT_POOL: copy constructor" << endl;
#endif
	size = source.size;
	num = source.num;
	
	if( size > 0 ){
		assert( source.cuts != NULL );
		assert( num >= 0 );
		cuts = new CUT_PLANE[size];

		for(int i=0; i<num; i++){
			cuts[i] = source.cuts[i];
		}
	}else{
		cuts = NULL;
	}

}

// assignment operator
CUT_POOL& CUT_POOL::operator=( const CUT_POOL& source )
{
#if debug_class
	cout << "CUT_POOL: assignment operator" << endl;
#endif

	if( this != &source ){ 	
		size = source.size;
		num = source.num;
		
		if( size > 0 ){
			assert( source.cuts != NULL );
			assert( num >= 0 );
			delete[] cuts;
			cuts = new CUT_PLANE[size];

			for(int i=0; i<num; i++){
				cuts[i] = source.cuts[i];
			}
		}else{
			cuts = NULL;
		}
	}
	return *this;
}

// destructor
CUT_POOL::~CUT_POOL()
{
#if debug_class
	cout << "CUT_POOL: destructor" << endl;
#endif
	delete[] cuts;
	cuts = NULL;
}

void	CUT_POOL::alloc(
	int	s
	)
{
	assert( s > 0 );
	assert( cuts == NULL );

	delete[] cuts;

	size = s;
	num = 0;
	cuts = new CUT_PLANE[size];
}

void CUT_POOL::add_cut(
	CUT_PLANE	*cut
	)
{
	assert( size > 0 );
	assert( num >= 0 );
	assert( num <= size );
	assert( cuts != NULL );

	if( num < size ){
		cuts[num] = *cut;
		num++;
	}else if( num >= size ){
		cout << "error: add_cut() " << endl;
		exit(-1);
	}
	
}



