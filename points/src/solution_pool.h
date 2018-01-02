#ifndef SOLUTION_POOL_H__
#define SOLUTION_POOL_H__

#include "solution.h"

class SOLUTION_POOL{
	int			size;
	int			num;
	SOLUTION 	*sols;
	public:

		SOLUTION_POOL();											// default constructor	
		SOLUTION_POOL( const SOLUTION_POOL &source );			// copy constructor
		SOLUTION_POOL& operator=( const SOLUTION_POOL& );		// assignment operator
		~SOLUTION_POOL();											// destructor

		void		alloc( int s );
		void		add_solution( SOLUTION sol );

};
#endif
