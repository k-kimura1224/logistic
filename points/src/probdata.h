#ifndef PROBDATA_H__
#define PROBDATA_H__

using namespace std;

class PROB_DATA{
	int			m;			// the number of nodes in L
	double		*B;		// [m*m]
	double		*B_;		// [m*m], Colmajor
	double		*Q;		// [m*m], B'B
	public:

		PROB_DATA();											// default constructor	
		PROB_DATA( const PROB_DATA &source );			// copy constructor
		PROB_DATA& operator=( const PROB_DATA& );		// assignment operator
		~PROB_DATA();											// destructor

		void set_data(int s_m, double *s_B, double *s_B_, double *s_Q );

		int		get_m(){ return m; }
		double*	get_B(){ return B; }
		double*	get_B_(){ return B_; }
		double*	get_Q(){ return Q; }
		double*	get_bvec(int k){ return B_+(k*m); }
};

#endif
