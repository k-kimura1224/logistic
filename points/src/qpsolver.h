#ifndef QPSOLVER_H_
#define QPSOLVER_H_

#include <assert.h>
#include "stopwatch.h"

using namespace std;


#define	SOLVE_PROBLEM_ACTIVE	1
// 0:dgesv
// 1:dsysv
// 2:dspsv? unimplimented

/*********************************************
	min	x'Qx + p'x
	s.t.	Ax <= b
			Cx = d
			l <= x <= u
			x \in R^n
*********************************************/

class QPsolver{
	double	*Q;	// [n*n]
	double	*p;	// [n]
	double	*A;	// [m*n] 
	double	*b;	// [m]
	double	*C;	// [k*n]
	double	*d;	// [k]
	double	*l;	// [n]
	double	*u;	// [n]
	int		n;
	int		m;
	int		k;
	double	*warm;
	double	*bestsol;
	double	bestval;

	double	ep;
	public:
	//STOPWATCH stopwatch;

	QPsolver();											// default constructor	
	QPsolver( const QPsolver &source );			// copy constructor
	QPsolver& operator=( const QPsolver& );	// assignment operator
	~QPsolver();										// destructor

	// set 
	void	set_dim( int s_n ){ n = s_n; }
	void	set_obj( int s_n, double *s_Q, double *s_p ){
		assert( n == s_n );
		Q = s_Q;
		p = s_p;
	}
	void	set_obj_quad( int s_n, double *s_Q ){
		assert( n == s_n );
		Q = s_Q;
	}
	void	set_obj_line( int s_n, double *s_p ){
		assert( n == s_n );
		p = s_p;
	}
	void	set_ineq( int s_n, int s_m, double *s_A, double *s_b ){
		assert( n == s_n );
		m = s_m;
		A = s_A;
		b = s_b;
	}
	void	set_equ( int s_n, int s_k, double *s_C, double *s_d ){
		assert( n == s_n );
		k = s_k;
		C = s_C;
		d = s_d;
	}
	void	set_lb( int s_n, double *s_l ){
		assert( n == s_n );
		l = s_l;
	}
	void	set_ub( int s_n, double *s_u ){
		assert( n == s_n );
		u = s_u;
	}
	void	set_warm( int s_n, double *s_warm ){
		assert( n == s_n );
		warm = s_warm;
	}
	void	disp_prob();
	void	solve();
	int	solve_activeset( double *x_new, bool *W_ineq, bool *W_u, bool *W_l);
	double	compute_stepsize( double *x, double *dx, bool *W_ineq, bool *W_u, bool *W_l);
	bool	update_activeset( double *x_new, bool *W_ineq, bool *W_u, bool *W_l); 
	double	compute_objval( double *x );

	double	get_bestval(){ return bestval; }
	double*	get_bestsol(){ return bestsol; }

	void		set_ep( double s_ep ){ ep = s_ep; }
	double	get_ep(){ return ep; }
};

#endif
