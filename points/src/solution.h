#ifndef SOLUTION_H__
#define SOLUTION_H__

using namespace std;


class SOLUTION{
	int		m;
	double	*val;
	double	objval;
	
	public:

		SOLUTION();											// default constructor	
		SOLUTION( const SOLUTION &source );			// copy constructor
		SOLUTION& operator=( const SOLUTION& );		// assignment operator
		~SOLUTION();											// destructor

		void set_sol(int s_m, double *s_val, double s_objval);

		double* get_solval(){ return val; }
};

#endif
