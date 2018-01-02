#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>
#include <time.h>

#include "stopwatch.h"
#include "svpsolver.h"

using namespace std;

#define debug_class	0
#define debug	0

STOPWATCH::STOPWATCH(){	// default constructor
#if debug_class
	cout << "STOPWATCH: default constructor" << endl;
#endif
	run = false;
	limit = TIMELIMIT;
}

STOPWATCH::STOPWATCH( const STOPWATCH &source )
{	// copy constructor
#if debug_class
	cout << "STOPWATCH: copy constructor" << endl;
#endif
	run = source.run;
	limit = source.limit;

	if( run == true ){
		begin = source.begin;
		end = source.end;
	}

}

// assignment operator
STOPWATCH& STOPWATCH::operator=( const STOPWATCH& source )
{
#if debug_class
	cout << "STOPWATCH: assignment operator" << endl;
#endif

	if( this != &source ){ 	
		run = source.run;
		limit = source.limit;
	
		if( run == true ){
			begin = source.begin;
			end = source.end;
		}
	}

	return *this;
}

// destructor
STOPWATCH::~STOPWATCH()
{
#if debug_class
	cout << "STOPWATCH: destructor" << endl;
#endif
}

void STOPWATCH::start()
{
	begin = time(NULL);
	end = time(NULL);
	run = true;
}

void STOPWATCH::stop()
{
	end = time(NULL);
}

int STOPWATCH::get_time()
{
	time_t buf = time(NULL);
	return (int)(buf - begin);
}

bool STOPWATCH::check_time()
{
	time_t buf = time(NULL);
	bool result = true;
	if( (int)(buf - begin) > limit ){
		//cout << "----- Time over -----" << endl;
		result = false;
	}else{
		result = true;
	}
	return  result;
}
