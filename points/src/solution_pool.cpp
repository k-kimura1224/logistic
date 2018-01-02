#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>

#include "solution.h"
#include "solution_pool.h"

using namespace std;

#define debug	0
#define debug_class	0

SOLUTION_POOL::SOLUTION_POOL(){	// default constructor
#if debug_class
	cout << "SOLUTION_POOL: default constructor" << endl;
#endif
	size = -1;
	num = -1;
	sols = NULL;
}

SOLUTION_POOL::SOLUTION_POOL( const SOLUTION_POOL &source )
{	// copy constructor
#if debug_class
	cout << "SOLUTION_POOL: copy constructor" << endl;
#endif
	size = source.size;
	num = source.num;
	
	if( size > 0 ){
		assert( source.sols != NULL );
		assert( num >= 0 );
		sols = new SOLUTION[size];

		for(int i=0; i<num; i++){
			sols[i] = source.sols[i];
		}
	}else{
		sols = NULL;
	}

}

// assignment operator
SOLUTION_POOL& SOLUTION_POOL::operator=( const SOLUTION_POOL& source )
{
#if debug_class
	cout << "SOLUTION_POOL: assignment operator" << endl;
#endif

	if( this != &source ){ 	
		size = source.size;
		num = source.num;
		
		if( size > 0 ){
			assert( source.sols != NULL );
			assert( num >= 0 );

			delete[] sols;
			sols = new SOLUTION[size];
	
			for(int i=0; i<num; i++){
				sols[i] = source.sols[i];
			}
		}else{
			sols = NULL;
		}
	}
	return *this;
}

// destructor
SOLUTION_POOL::~SOLUTION_POOL()
{
#if debug_class
	cout << "SOLUTION_POOL: destructor" << endl;
#endif
	delete[] sols;
	sols = NULL;
}

void	SOLUTION_POOL::alloc(
	int	s
	)
{
	assert( s > 0 );
	assert( sols == NULL );

	delete[] sols;

	size = s;
	num = 0;
	sols = new SOLUTION[size];
}

void SOLUTION_POOL::add_solution(
	SOLUTION		sol
	)
{
	assert( size > 0 );
	assert( num >= 0 );
	assert( num <= size );
	assert( sols != NULL );

	if( num == 0 ){
		sols[0] = sol;
		num++;
	}else if( num < size ){
		for(int i=num-1; i>=0; i--){
			if( sols[i].get_solval() > sol.get_solval() ){
				sols[i+1] = sols[i];
			}else{
				sols[i+1] = sol;
				break;
			}
			if( i == 0 ){
				sols[i] = sol;
			}
		}
		num++;
	}else if( num == size ){
		if( sol.get_solval() < sols[size-1].get_solval() ){
			num--;
			for(int i=num-1; i>=0; i--){
				if( sols[i].get_solval() > sol.get_solval() ){
					sols[i+1] = sols[i];
				}else{
					sols[i+1] = sol;
					break;
				}
				if( i == 0 ){
					sols[i] = sol;
				}
			}
			num++;
		}
	}else{
		cout << "error: add_solution() " << endl;
		exit(-1);
	}
	
}



