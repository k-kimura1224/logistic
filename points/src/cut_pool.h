#ifndef CUT_POOL_H__
#define CUT_POOL_H__

#include "cut.h"

class CUT_POOL{
	int			size;
	int			num;
	CUT_PLANE 	*cuts;
	public:

		CUT_POOL();											// default constructor	
		CUT_POOL( const CUT_POOL &source );			// copy constructor
		CUT_POOL& operator=( const CUT_POOL& );		// assignment operator
		~CUT_POOL();											// destructor

		void		alloc( int s );
		void		add_cut( CUT_PLANE *cut );

		int			get_size(){ return size; }
		CUT_PLANE*	get_cut(int i){
			assert( i>=0 && i<size && i<num );
			return &cuts[i];
		}

};
#endif
