#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>
#include <math.h>
#include <omp.h>

#include "qpsolver.h"
#include "vector.h"

using namespace std;


#define debug_class	0
#define debug	0

QPsolver::QPsolver(){	// default constructor
#if debug_class
	cout << "QPsolver: default constructor" << endl;
#endif
	Q = NULL;
	p = NULL;
	A = NULL;
	b = NULL;
	C = NULL;
	d = NULL;
	l = NULL;
	u = NULL;
	n = 0;
	m = 0;
	k = 0;
	warm = NULL;
	bestsol = NULL;
	bestval = 0.0;
}

QPsolver::QPsolver( const QPsolver &source )
{	// copy constructor
#if debug_class
	cout << "QPsolver: copy constructor" << endl;
#endif
	Q = source.Q;	
	p = source.p;
	A = source.A;
	b = source.b;
	C = source.C;
	d = source.d;
	l = source.l;
	u = source.u;
	n = source.n;
	m = source.m;
	k = source.k;
	warm = source.warm;
	bestval = source.bestval;
	if( source.bestsol != NULL ){
		assert( n > 0 );
		bestsol = new double[n];
		for(int i=0; i<n; i++){
			bestsol[i] = source.bestsol[i];
		}
	}else{
		bestsol = NULL;
	}
}

// assignment operator
QPsolver& QPsolver::operator=( const QPsolver& source )
{
#if debug_class
	cout << "QPsolver: assignment operator" << endl;
#endif

	if( this != &source ){ 	
		Q = source.Q;	
		p = source.p;
		A = source.A;
		b = source.b;
		C = source.C;
		d = source.d;
		l = source.l;
		u = source.u;
		n = source.n;
		m = source.m;
		k = source.k;
		warm = source.warm;
		bestval = source.bestval;
		if( source.bestsol != NULL ){
			assert( n > 0 );
			bestsol = new double[n];
			for(int i=0; i<n; i++){
				bestsol[i] = source.bestsol[i];
			}
		}else{
			bestsol = NULL;
		}
	}
	return *this;
}

// destructor
QPsolver::~QPsolver()
{
#if debug_class
	cout << "QPsolver: destructor" << endl;
#endif
	delete[] bestsol;
	bestsol = NULL;
}

void	QPsolver::disp_prob(){

	if( n > 0 ){
		cout << "dim: " << n << endl;
	
		if( Q ){
			cout << "Q:" << endl;
			printM( n, n, Q);
		}
		if( p ){
			cout << "p:" << endl;
			printv( n, p);
		}
		if( m ){
			cout << "m: " << m << endl;
			if( A ){
				cout << "A:" << endl;
				printM( m, n, A);
			}
			if( b ){
				cout << "b:" << endl;
				printv( m, b);
			}
		}
		if( k ){
			cout << "k: " << k << endl;
			if( C ){
				cout << "C:" << endl;
				printM( k, n, C);
			}
			if( d ){
				cout << "d:" << endl;
				printv( k, d);
			}
		}
		if( u ){
			cout << "u:" << endl;
			printv( n, u);
		}
		if( l ){
			cout << "l:" << endl;
			printv( n, l);
		}
		if( warm ){
			cout << "init:" << endl;
			printv( n, warm);
		}
	}

}

void	QPsolver::solve(){

	assert( n > 0 );
	assert( warm != NULL );
	assert( Q != NULL );
	
	double	ep = get_ep();
	double	*x;
	x = new double[n];

	for(int i=0; i<n; i++){
		x[i] = warm[i];
	}

#if debug
	cout << "initial point:";
	for(int i=0; i<n; i++){
		cout << x[i] << " ";
	}
	cout << endl;
#endif

	// find active constraint sets {{
	bool	*W_ineq = NULL;
	if( m > 0 ){
		assert( A != NULL );
		assert( b != NULL );
		W_ineq	=	new bool[m];

		for(int i=0; i<m; i++){
		//	double buf = b[i] - Com_dot( A+(i*n), x, n);

		//	if( fabs( buf ) < ep ){
		//		W_ineq[i] = true;
		//	}else if( buf > 0 ){
		//		W_ineq[i] = false;
		//	}else{ 
		//		cout << "the given initial point is not a feasible solution." << endl;
		//		cout << buf << endl;
		//		disp_prob();
		//		exit(-1);
		//	}
			double buf = Com_dot( A+(i*n), x, n);
			if( Equal( b[i], buf, ep ) == true ){
				W_ineq[i] = true;
			}else if( b[i] - buf > 0 ){
				W_ineq[i] = false;
			}else{
				cout << "the given initial point is not a feasible solution." << endl;
				cout << b[i] - buf << endl;
				disp_prob();
				exit(-1);
			}
		}

	}

#if debug
	cout << "active constraints:";
	if( W_ineq != NULL ){
		for(int i=0; i<m; i++){
			if( W_ineq[i] == true ){
				cout << i << " ";
			}
		}
		cout << endl;
	}else{
		cout << "There are not inequality constraints" << endl;
	}
#endif
	// }} find active constraint sets

	// check feasibility for equation constraints{
	if( k > 0 ){
		assert( C != NULL );
		assert( d != NULL );
		
		for(int i=0; i<k; i++){
			if( Equal( Com_dot( C+(i*n), x, n), d[i], ep) == false ){
				cout << "the given initial point is not a feasible solution." << endl;
				cout << "aa" << endl;
				exit(-1);
			}
		}
	}
	// }} check feasibility for equation constraints

	// find active bounds {{ 
	bool	*W_u = NULL;
	bool	*W_l = NULL;
	if( u != NULL ){
		W_u = new bool[n];
		for(int i=0; i<n; i++){
			//double buf = u[i] - x[i];

			//if( fabs( buf ) < ep ){
			//	W_u[i] = true;
			//}else if( buf > 0){
			//	W_u[i] = false;
			//}else{
			//	cout << "the given initial point is not a feasible solution." << endl;
			//	cout << "bb:" << buf << endl;
			//	cout << "point:";
			//	printv( n, warm);
			//	cout << "l:";
			//	printv( n, l);
			//	cout << "u:";
			//	printv( n, u);
			//	exit(-1);
			//}
			double buf = u[i] - x[i];
			if( Equal( u[i], x[i], ep) == true ){
				W_u[i] = true;
				if( buf != 0 ) x[i] = u[i];
			}else if( buf > 0){
				W_u[i] = false;
			}else{
				cout << "the given initial point is not a feasible solution." << endl;
				cout << "bb:" << buf << endl;
				cout << "point:";
				printv( n, warm);
				cout << "l:";
				printv( n, l);
				cout << "u:";
				printv( n, u);
				exit(-1);
			}
		
		}
	}
	if( l != NULL ){
		W_l = new bool[n];
		for(int i=0; i<n; i++){
	//		double buf = x[i] - l[i];
	//		if( fabs( buf ) < ep ){
	//			W_l[i] = true;
	//		}else if( buf > 0 ){
	//			W_l[i] = false;
	//		}else{
	//			cout << "the given initial point is not a feasible solution." << endl;
	//			cout << "cc:" << buf << endl;
	//			exit(-1);
	//		}
			double buf = x[i] - l[i];
			if( Equal( x[i], l[i], ep) ){
				W_l[i] = true;
				if( buf != 0 ) x[i] = l[i];
			}else if( buf > 0 ){
				W_l[i] = false;
			}else{
				cout << "the given initial point is not a feasible solution." << endl;
				cout << "cc:" << buf << endl;
				exit(-1);
			}
		
		}
	}

#if debug
	cout << "active upper bounds:";
	if( W_u != NULL ){
		for(int i=0; i<n; i++){
			if( W_u[i] == true ){
				cout << i << " ";
			}
		}
		cout << endl;
	}else{
		cout << "There are not upper bounds" << endl;
	}
	cout << "active lower bounds:";
	if( W_l != NULL ){
		for(int i=0; i<n; i++){
			if( W_l[i] == true ){
				cout << i << " ";
			}
		}
		cout << endl;
	}else{
		cout << "There are not lower bounds" << endl;
	}
#endif
	// }} find active bounds

	double	*x_new = NULL;
	double	*dx = NULL;
	bool		update;

	x_new = new double[n];
	dx = new double[n];

	int debug_ct = 0;
	while(1){
		if( debug_ct > 1000 ){
			cout << "error: qpsolver.cpp" << endl;
			cout << "th: " << omp_get_thread_num() << endl;	
			exit(1);
		}
		debug_ct++;

#if debug
	cout << "active constraints:";
	if( W_ineq != NULL ){
		for(int i=0; i<m; i++){
			if( W_ineq[i] == true ){
				cout << i << " ";
			}
		}
		cout << endl;
	}else{
		cout << "There are not inequality constraints" << endl;
	}
	cout << "active upper bounds:";
	if( W_u != NULL ){
		for(int i=0; i<n; i++){
			if( W_u[i] == true ){
				cout << i << " ";
			}
		}
		cout << endl;
	}else{
		cout << "There are not upper bounds" << endl;
	}
	cout << "active lower bounds:";
	if( W_l != NULL ){
		for(int i=0; i<n; i++){
			if( W_l[i] == true ){
				cout << i << " ";
			}
		}
		cout << endl;
	}else{
		cout << "There are not lower bounds" << endl;
	}
#endif

		// solve the problem with active set
		int argmin = solve_activeset( x_new, W_ineq, W_u, W_l);
		
#if debug
	cout << "argmin: " << argmin << endl;
#endif
			
		// dx = xnew - x 
		Com_linecomb( x_new, x, n, 1.0, -1.0, dx);

#if debug
	cout << "dx: ";
	printv( n, dx);
#endif
		if( Com_nrm( dx, n) > ep ){
			// compute stepsize
			double alpha = compute_stepsize( x, dx, W_ineq, W_u, W_l);	
#if debug
	cout << "alpha: " << alpha << endl;
#endif
			assert( alpha >= 0 );
			assert( alpha <= 1 );
			// x_new = x + alpha dx
			Com_linecomb( x, dx, n, 1.0, alpha, x_new);
#if debug
	cout << "x_new: ";
	printv( n, x_new);
#endif
		}
		
		/* x_new is feasible for the original problem now */
		
		// update W_ineq, W_u and W_l
		update = update_activeset( x_new, W_ineq, W_u, W_l);

		if( update == true ){
			Copy_vec( x_new, x, n);
			continue;
		}

		if( argmin == -1 ){
			break;
		}else{
			int ct = 0;
			update = false;
			if( W_ineq != NULL ){
				assert( m > 0 );
				if( count( m, W_ineq) > argmin ){
					for(int i=0; i<m; i++){
						if( W_ineq[i] == true ){
							if( ct == argmin ){
								W_ineq[i] = false;
								update = true;
								break;
							}else{
								ct++;
							}
						}
					}
				}else{
					ct += count( m, W_ineq);
				}
			}

			if( update == true ){
				Copy_vec( x_new, x, n);
				continue;
			}

			if( W_u != NULL ){
				if( ct + count( n, W_u) > argmin ){
					for(int i=0; i<n; i++){
						if( W_u[i] == true ){
							if( ct == argmin ){
								W_u[i] = false;
								update = true;
								break;
							}else{
								ct++;
							}
						}
					}
				}else{
					ct += count( n, W_u);
				}
			}

			if( update == true ){
				Copy_vec( x_new, x, n);
				continue;
			}

			if( W_l != NULL ){
				for(int i=0; i<n; i++){
					if( W_l[i] == true ){
						if( ct == argmin ){
							W_l[i] = false;
							update = true;
							break;
						}else{
							ct++;
						}
					}
				}
			}
			
			if( update == true ){
				Copy_vec( x_new, x, n);
				continue;
			}else{
				cout << "error: QPsolver" << endl;
				cout << "th: " << omp_get_thread_num() << endl;	
				exit(-1);
			}

		}
	}

	bestsol = new double[n];

	Copy_vec( x_new, bestsol, n);

	bestval = compute_objval( bestsol );

#if debug
	cout << "bestsol: ";
	printv( n, bestsol);
	cout << "bestval: " << bestval << endl;
#endif

	delete[] x;
	delete[] W_ineq;
	delete[] W_u;
	delete[] W_l;
	delete[] x_new;
	delete[] dx;
}

int QPsolver::solve_activeset(
	double	*x_new,
	bool		*W_ineq,
	bool		*W_u,
	bool		*W_l
	)
{
	/**********************************************
		min	x'Qx + p'x
		s.t.	F x = g

		by solving H z = r.

		This is

		| 2Q , F' || x | = | -p |
		| F  , 0  || y |   |  g |
		
	**********************************************/
	assert( x_new != NULL );

	int	w_equ = k;
	int	w_ineq = count( m, W_ineq);
	int	w_u = count( n, W_u);
	int	w_l = count( n, W_l);
	int	w = w_equ + w_ineq + w_u + w_l;
			
#if debug
	cout << "w: " << w << endl;
#endif

	double	*F = NULL;	// [w*n]
	double	*g = NULL;	// [w]

	if( w > 0 ){
		F = new double[w*n];
		g = new double[w];
	}

	/**********************************************
		1. C
		2. A
		3. u
		4. l
	**********************************************/

	// 1. C
	int		ct=0;
	if( w_equ > 0 ){
		assert( k > 0 );
		assert( C != NULL );
		assert( w >= k );
		assert( F != NULL );
		
		Copy_vec( C, F, k*n);
		Copy_vec( d, g, k);

		ct += k;

		assert( w >= ct );
	}

	// 2. A
	if( w_ineq > 0 ){
		assert( m > 0 );
		assert( A != NULL );
		assert( w >= w_ineq );
		assert( F != NULL );

		for(int i=0; i<m; i++){
			if( W_ineq[i] == true ){
				Copy_vec( &A[i*n], &F[ct*n], n);
				g[ct] = b[i];
				ct++;
			}
		}

		assert( w >= ct );
	}

	// 3. u
	if( w_u > 0 ){
		assert( u != NULL );
		assert( w >= w_u );
		assert( F != NULL );

		for(int i=0; i<n; i++){
			if( W_u[i] == true ){
				for(int j=0; j<n; j++){
					if( j == i ){
						F[(ct*n)+j] = 1.0;
					}else{
						F[(ct*n)+j] = 0.0;
					}
				}
				g[ct] = u[i];
				ct++;
			}
		}

		assert( w >= ct );
	}

	// 4. l
	if( w_l > 0 ){
		assert( l != NULL );
		assert( w >= w_l );
		assert( F != NULL );

		for(int i=0; i<n; i++){
			if( W_l[i] == true ){
				for(int j=0; j<n; j++){
					if( j == i ){
						F[(ct*n)+j] = -1.0;
					}else{
						F[(ct*n)+j] = 0.0;
					}
				}
				g[ct] = - l[i];
				ct++;
			}
		}

		assert( w >= ct );
	}

	assert( w == ct );	

#if debug
	cout << "F:" << endl;
	//printM( w, n, F);
	cout << "g: ";
	//printv( w, g);
#endif

	int		h = n+w;
	double	*H = NULL;	// [h*h]

	H = new double[h*h];

	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			H[(i*h)+j] = 2 * Q[(i*n)+j];
		}
	}

	if( w > 0 ){
		for(int i=n; i<h; i++){
			for(int j=0; j<n; j++){
				H[(i*h)+j] = F[((i-n)*n)+j];
				H[(j*h)+i] = F[((i-n)*n)+j];
			}
			for(int j=n; j<h; j++){
				H[(i*h)+j] = 0.0;
			}
		}
	}

	double	*z = NULL;	// [h]
	double	*r = NULL;	// [h]

	z = new double[h];
	r = new double[h];

	if( p == NULL ){
		for(int i=0; i<n; i++){
			r[i] = 0.0;
		}
	}else{
		for(int i=0; i<n; i++){
			r[i] = - p[i];
		}
	}

	if( w > 0 ){
		for(int i=n; i<h; i++){
			assert( w > (i-n) );
			r[i] = g[(i-n)];
		}
	}

#if debug
	cout << "H:" << endl;
	//printM( h, h, H);
	cout << "r: ";
	//printv( h, r);
#endif
	
	int	com;
	switch ( SOLVE_PROBLEM_ACTIVE ) {
		case 0:
		{
			com = Com_LS_dgesv( H, r, h, z);
			break;
		}
		case 1:
		{
			com = Com_LS_dsysv( H, r, h, z);
			break;
		}
		default:
		{
			cout << "error: SOLVE_PROBLEM_ACTIVE is " << SOLVE_PROBLEM_ACTIVE << endl;
			exit(-1);
			break;
		}
	}

#if debug
	cout << "com: " << com << endl;
	cout << "z: ";
	printv( h, z);
#endif
	
	if( com != 0 ){
		cout << "error: com = " << com << endl;
		exit(-1);
	}


	Copy_vec( z, x_new, n);

	double min = 0.0;
	int memo = -1;
	double ep = get_ep();

	for(int i=n+k; i<h; i++){
		if( z[i] < min && fabs( z[i] ) > ep ){
			min = z[i];
			memo = i - ( n + k );
		}
	}

	delete[] F;
	delete[] g;
	delete[] H;
	delete[] z;
	delete[] r;

	return memo;
}

double	QPsolver::compute_stepsize(
	double	*x,
	double	*dx,
	bool		*W_ineq, 
	bool		*W_u,
	bool		*W_l
	)
{
	assert( x != NULL );
	assert( dx != NULL );
	assert( Com_nrm( dx, n) >= get_ep() );

	double	alpha = 1.0;
	double ep = get_ep();
	
	if( W_ineq != NULL ){
		assert( m > 0 );
		assert( A != NULL );
		assert( b != NULL );
		for(int i=0; i<m; i++){
			if( W_ineq[i] == false ){
				double	dot = Com_dot( &A[i*n], dx, n);
				if( dot > ep ){
					double	a = ( b[i] - Com_dot( &A[i*n], x, n) ) / dot;
					if( alpha > a ){
						alpha = a;
					}
				}
			}
		}
	}


	if( W_u != NULL ){
		assert( n > 0 );
		assert( u != NULL );
		for(int i=0; i<n; i++){
			if( W_u[i] == false && dx[i] > ep ){
				double	a = ( u[i] - x[i] ) / dx[i];
				if( alpha > a ){
					alpha = a;
				}
			}
		}
	}

	if( W_l != NULL ){
		assert( n > 0 );
		assert( l != NULL );
		for(int i=0; i<n; i++){
			if( W_l[i] == false && dx[i] < -ep ){
				double	a = ( l[i] - x[i] ) / dx[i];
				if( alpha > a ){
					alpha = a;
				}
			}
		}
	}

	if( fabs( alpha ) < ep ){
		alpha = 0.0;
	}

	return alpha;
}

bool QPsolver::update_activeset(
	double	*x_new,
	bool		*W_ineq,
	bool		*W_u,
	bool		*W_l
	)
{
	double ep = get_ep();
	bool	update = false;

	if( W_ineq != NULL ){
		assert( m > 0 );
		assert( A != NULL );
		assert( b != NULL );
		for(int i=0; i<m; i++){
			if( W_ineq[i] == false ){
				double buf = Com_dot( &A[i*n], x_new, n);
				if( Equal( b[i], buf, ep) == true ){
					W_ineq[i] = true;
					update = true;
				}
				//if( fabs( b[i] - Com_dot( &A[i*n], x_new, n)) < ep ){
				//	W_ineq[i] = true;
				//	update = true;
				//}
			}
		}
	}

	if( W_u != NULL ){
		assert( n > 0 );
		assert( u != NULL );
		for(int i=0; i<n; i++){
			if( W_u[i] == false ){
				if( Equal( u[i], x_new[i], ep) == true ){
				//if( fabs( u[i] - x_new[i] ) < ep ){
					W_u[i] = true;
					update = true;
					if( u[i] - x_new[i] != 0 ) x_new[i] = u[i];
				}
			}
		}
	}

	if( W_l != NULL ){
		assert( n > 0 );
		assert( l != NULL );
		for(int i=0; i<n; i++){
			if( W_l[i] == false ){
				if( Equal( x_new[i], l[i], ep) == true ){
				//if( fabs( x_new[i] - l[i] ) < ep ){
					W_l[i] = true;
					update = true;
					if( x_new[i] - l[i] != 0 ) x_new[i] = l[i];
				}
			}
		}
	}

	return update;

}

double QPsolver::compute_objval(
	double	*x
	)
{
	assert( n > 0 );
	assert( Q != NULL );
	assert( x != NULL );

	double	obj = 0.0;
	
	double	*Qx; //[n]
	Qx = new double[n];

	Gen_ZeroVec( n, Qx);
	Com_mat_Ax( Q, n, n, x, Qx);

	obj += Com_dot( x, Qx, n);

	if( p != NULL ){
		obj += Com_dot( x, p, n);
	}

	delete[] Qx;

	return obj;
}
