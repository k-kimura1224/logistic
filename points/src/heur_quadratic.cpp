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

#include "svpsolver.h"
#include "probdata.h"
#include "vector.h"
#include "solution.h"
#include "solution_pool.h"
#include "node.h"

using namespace std;

#define debug	0
#define log		1

void	SVPsolver::heur_quadratic( int sel )
{
	assert( sel >= 0 );
	assert( sel < listsize );

	int		m = probdata.get_m();

	double	val;
	double	val_new;
	double	*solvals = NULL;
	double	*solvals_new = NULL;

	list<NODE>::iterator it = NodeList.begin();
	advance( it, sel);
	double	*relaxsolvals = it->get_relaxsolval();
	double	*u = it->get_ub();
	double	*l = it->get_lb();

	assert( m > 0 );
	assert( relaxsolvals != NULL );
	assert( u != NULL );
	assert( l != NULL );

	solvals = new double[m];
	solvals_new = new double[m];

	for(int i=0; i<m; i++){
		solvals[i] = round( relaxsolvals[i] );
	}
	
	val = compute_objval( solvals );

	double	lam, lam1, lam2;
	double	t, t1, t2;
	double s1, s2;
	double *normv = new double[m];
	double	*Bx = new double[m];
	double	*B_ = probdata.get_B_();
	
	assert( B_ != NULL );
	assert( norm != NULL );

#if debug
	cout << "best:" << bestval << endl;
	cout << "val:" << val << endl;
#endif
	for(int i=0; i<m; i++){
		normv[i] = norm[i] * norm[i];
	}
	while(1){

		Com_mat_Ax( B_, m, m, solvals, Bx);
		val_new = val;

		for(int i=0; i<m; i++){
			if( u[i] != l[i] ){
				lam = val;
				s1 = Com_dot( &B_[i*m], Bx, m);
				s2 = normv[i]; //norm[i] * norm[i];
				t = - s1 / s2;
	
				if( t < l[i] - solvals[i] ){
					t = l[i] - solvals[i];
					if( t != 0.0 ){
						lam += t * ( t*s2 + 2*s1 );
					}
				}else if( u[i] - solvals[i] < t ){
					t = u[i] - solvals[i];
					if( t != 0.0 ){
						lam += t * ( t*s2 + 2*s1 );
					}
				}else{
					t1 = ceil(t);
					lam1 = val;
					if( t1 != 0.0 ){
						lam1 += t1 * ( t1*s2 + 2*s1 );
					}
	
					t2 = floor(t);
					lam2 = val;
					if( t2 != 0.0 ){
						lam2 += t2 * ( t2*s2 + 2*s1 );
					}
	
					if( lam1 > lam2 ){
						lam = lam2;
						t = t2;
					}else{
						lam = lam1;
						t = t1;
					}
				}
		
	#if debug
		cout << "i:" << i << "-> " << lam << endl;
	#endif
	
				if( val_new > lam ){
					val_new = lam;
					Copy_vec( solvals, solvals_new, m);
					solvals_new[i] = solvals[i] + t;
				}
			}
		}

#if debug
		cout << "val_new:" << val_new << endl;
#endif
		// val and val_new are integer.
		if( val  - val_new > 0.1){
			//cout << "val:" << val << " val_new:" << val_new << endl;
			val = val_new;
			Copy_vec( solvals_new, solvals, m);
		}else{
#if debug
			cout << "break" << endl;
#endif
			break;
		}
	}
	
	SOLUTION	solution;
	solution.set_sol( m, solvals, val);
	pool.add_solution( solution );

	if( bestval > val ){
		bestval = val;
		bestsol = solution;
		Appfac = sqrt( bestval ) / _Appfac;
#if log
		cout << "*********************************" << endl;
		cout << "*   heur_quadratic  dpt: " << it->get_dpt() << endl;
		cout << "*********************************" << endl;
#endif
	}

	delete[] solvals;
	delete[] solvals_new;
	delete[] Bx;
	delete[] normv;
}
