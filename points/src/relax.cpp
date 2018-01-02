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
#include "node.h"
#include "probdata.h"
#include "qpsolver.h"
#include "solution.h"
#include "vector.h"
#include "Schmidt_manager.h"
#include "cut_pool.h"
#include "cut.h"

using namespace std;

#define debug	0


RelaxResult	SVPsolver::solve_relaxation(
	int	sel
	)
{
	RelaxResult result;
	list<NODE>::iterator it = NodeList.begin();
	advance( it, sel);
	
	//if( (bestval - it->get_lowerbound())/bestval < 1.0e-12 ){
	if( bestval - it->get_lowerbound() < 1.0 ){
		return INFEASIBLE;
	}

	if( it->get_zero() == true ){
		//result = solve_relaxation_BIN(sel);
		result = solve_relaxation_BIN_sch(sel);
	}else{
		result = solve_relaxation_INT(sel);
	}

	return result;
}

RelaxResult	SVPsolver::solve_relaxation_BIN(
	int	sel
	)
{
	int		m = probdata.get_m();
//	double	*Q = probdata.get_Q();
//	double	*u = NodeList[sel].get_ub();
//	double	*l = NodeList[sel].get_lb();
//
//	assert( m > 0 );
//
//	QPsolver	 qps;
//	int dim = m;
//
//	for(int i=0; i<m; i++){
//		for(int j=(int)lb[i]; j<=(int)ub[i]; j++){
//			dim++;
//		}
//	}
//
//	qps.set_dim( dim );
//
//	double	*P;
//	P = new double[dim*dim];
//
//	int	ct = 0;
//	for(int i=0; i<dim; i++){
//		for(int j=0; j<dim; j++){
//			if( i<m && j<m ){
//				P[ct] = Q[(i*m)+j];
//			}else{
//				P[ct] = 0.0;
//			}
//		}
//	}
//
//	qps.set_obj_quad( dim, P);
//
//	double	*A;
//	A = new double[dim];
//
//	for(int i=0; i<m; i++){
//		A[i] = 0.0;
//	}
//
//	ct = m;
//	for(int i=0; i<m; i++){
//		for(int j=(int)lb[i]; j<=(int)ub[i]; j++){
//			if( j==0 ){
//				A[ct] = 1.0;
//			}else{
//				A[ct] = 0.0;
//			}
//			ct++;
//		}
//	}
//
//	double	*b;
//	b = new double[1];
//	b[0] = (double)m - 1.0;
//
//	qps.set_ineq( dim, 1, A, b);
//
//	double	*C;
//	C = new double[2*m*dim];
//
//	ct = 0;
//	for(int k=0; k<m; k++){
//		for(int j=0; j<m; j++){
//			if( j==k ){
//				C[ct] = 1.0;
//			}else{
//				C[ct] = 0.0;
//			}
//			ct++;
//		}
//		for(int i=0; i<m; i++){
//			for(int j=(int)lb[i]; j<=(int)ub[i]; j++){
//				if( i==k ){
//					C[ct] = -(double)j;
//				}else{
//					C[ct] = 0.0;
//				}
//				ct++;
//			}
//		}
//	}
//
//	assert( ct == dim*m );
//
//	for(int k=0; k<m; k++){
//		for(int j=0; j<m; j++){
//			C[ct] = 0.0;
//			ct++;
//		}
//		for(int i=0; i<m; i++){
//			for(int j=(int)lb[i]; j<=(int)ub[i]; j++){
//				if( i==k ){
//					C[ct] = 1.0;
//				}else{
//					C[ct] = 0.0;
//				}
//				ct++;
//			}
//		}
//	}
//
//	assert( ct == m*2*dim );
//
//	double *d;
//	d = new double[2*m];
//
//	for(int i=0; i<m; i++){
//		d[i] = 0.0;
//	}
//	for(int i=m; i<2*m; i++){
//		d[i] = 1.0;
//	}
//
//	qps.set_equ( dim, 2*m, C, d);
//
//	double	*U;
//	double	*L;
//
//	U = new double[dim];
//	L = new double[dim];
//
//	for(int i=0; i<m; i++){
//		U[i] = u[i];
//		L[i] = l[i];
//	}
//
//	for(int i=m; i<dim; i++){
//			U[i] = 1.0;
//			L[i] = 0.0;
//	}
//
//	qps.set_lb( dim, L);
//	qps.set_ub( dim, U);
//
////	qps.disp_prob();
//
//	double	*base = NodeList[sel].get_warm();
//	double	*warm;
//	warm = new double[dim];
//
//	for(int i=0; i<m; 
////	qps.set_warm( m, warm);
////	qps.solve();
//
//	delete[] P;
//	delete[] A;
//	delete[] b;
//	delete[] C;
//	delete[] d;
//	delete[] U;
//	delete[] L;
//	delete[] warm;

	RelaxResult result;

	list<NODE>::iterator it = NodeList.begin();
	advance( it, sel);

	bool infe=false;
	double *u = it->get_ub();
	double *l = it->get_lb();

	for(int i=0; i<m; i++){
		if( u[i] != 0.0 || l[i] != 0.0 ){
			break;
		}
		if( i == m-1 ) infe = true;
	}

	if( infe == true ){
		result = INFEASIBLE;
	}else{
		it->set_lowerbound( 0.0 );
		result = FEASIBLE;
	}

	return result;
}

RelaxResult	SVPsolver::solve_relaxation_BIN_sch(
	int	sel
	)
{
	int		m = probdata.get_m();

	list<NODE>::iterator it = NodeList.begin();
	advance( it, sel);
	double *u = it->get_ub();
	double *l = it->get_lb();

	int ct = 0;

	for(int i=0; i<m; i++){
		if( u[i] == 0.0 && l[i] == 0.0 ){
			sch.set_z_i( i, false);
			ct++;
		}else{
			sch.set_z_i( i, true);
		}
	}

	if( ct == m ) return INFEASIBLE;

	assert( sch.get_n() > 0 );

	RelaxResult result;

	sch.compute_GS();

	if( it->get_lowerbound() < sch.get_min() ){
		it->set_lowerbound( sch.get_min() );
	}

	if( sch.get_min() >= bestval ){
		result = INFEASIBLE;
#if debug
		cout << "*********************************" << endl;
		cout << "*   Schmidt_manager  #zero: " << ct << endl;
		cout << "*********************************" << endl;
#endif
	}else{
		result = FEASIBLE;
	}

	return result;
}

RelaxResult	SVPsolver::solve_relaxation_INT(
	int	sel
	)
{
	list<NODE>::iterator it = NodeList.begin();
	advance( it, sel);

	int		m = probdata.get_m();
	double	*Q = probdata.get_Q();
	double	*u = it->get_ub();
	double	*l = it->get_lb();
	double	*warm = it->get_warm();
	
	int	n = 0;
	for(int i=0; i<m; i++){
		if( l[i] != u[i] ) n++;
	}

	if( n == 0 ){
		return INFEASIBLE;
	}
	
	double	*subQ = NULL;
	double	*subu = NULL;
	double	*subl = NULL;
	double	*subw = NULL;
	double	*p = NULL;
	double	c_term = 0;

	QPsolver	 qps;
	int	ct = 0;

	qps.set_dim( n );
	if( n == m ){
		qps.set_obj_quad( m, Q);
		qps.set_lb( m, l);
		qps.set_ub( m, u);
		qps.set_warm( m, warm);
	}else{
		subQ = new double[n*n];
		ct = 0;
		for(int i=0; i<m; i++){
			if( l[i] != u[i] ){
				for(int j=0; j<m; j++){
					if( l[j] != u[j] ){
						subQ[ct++] = Q[i+(j*m)];
					}
				}
			}
		}
		assert( ct == n*n );
	//	printM( n, n, subQ);
	//	cout << endl;
	//	cout << "[1,0]=" << subQ[0+(1*n)] <<endl;
	//	cout << "[0,1]=" << subQ[1+(0*n)] <<endl;
		assert( subQ[0+(1*n)] == subQ[1+(0*n)] );

		if( it->get_sumfixed() != NULL ){
			p = new double[n];
			ct = 0;
			for(int i=0; i<m; i++){
				if( l[i] != u[i] ){
					p[ct++] = 2*Com_dot( probdata.get_bvec(i), it->get_sumfixed(), m);
				}
			}
			assert( ct == n );

			qps.set_obj( n, subQ, p);
			c_term = Com_dot( it->get_sumfixed(), it->get_sumfixed(), m);
		}else{
			qps.set_obj_quad( n, subQ);
		}

		subu = new double[n];
		subl = new double[n];
		subw = new double[n];
		
		ct = 0;

		for(int i=0; i<m; i++){
			if( l[i] != u[i] ){
				subu[ct] = u[i];
				subl[ct] = l[i];
				subw[ct] = warm[i];
				ct++;
			}
		}
		assert( ct == n );

		qps.set_ub( n, subu);
		qps.set_lb( n, subl);
		qps.set_warm( n, subw );

		
	}




	// using oa_cut{{
	double	*A = NULL;
	double	*b = NULL;
	if( CUT_OA == true ){
		assert( oa_cpool.get_size() > 0 );
		cout << "stop!!!!!!!" << endl;
		exit(1);

		int	k 	= oa_cpool.get_size() * 2;
		A = new double[k*m];
		b = new double[k];

		CUT_PLANE *cut;
		ct = 0;
		for(int i=0; i<oa_cpool.get_size(); i++){
			cut = oa_cpool.get_cut( i );

			assert( cut->get_coef() != NULL );
			assert( cut->get_lb() != NULL );
			assert( cut->get_ub() != NULL );
			assert( cut->get_ub() > cut->get_lb() );

			Copy_vec( cut->get_coef(), &A[ct*m], m);
			for(int j=0; j<m; j++){
				A[(ct*m)+m+j] = - A[(ct*m)+j];
			}
			b[ct] = *(cut->get_ub());
			b[ct+1] = - *(cut->get_lb());
			
			ct += 2;
		}
		assert( ct == k );

		qps.set_ineq( m, k, A, b);
	}
	// }} using oa_cut

	//qps.disp_prob();

	clock_t start = clock();

	qps.set_ep( 1e-12 );

	qps.solve();
   clock_t end = clock();
	QP_time += (double)(end-start)/CLOCKS_PER_SEC;
	
	double *relaxvals = new double[m];
	double *qpbestvals =qps.get_bestsol();
	ct = 0;
	for(int i=0; i<m; i++){
		if( l[i] != u[i] ){
			relaxvals[i] = qpbestvals[ct++];
		}else{
			relaxvals[i] = l[i];
		}
	}
	assert( ct == n );
	
	it->set_relaxsolval( relaxvals );

	double qpbestval = qps.get_bestval() + c_term;

	if( it->get_lowerbound() < qpbestval ){
		it->set_lowerbound( qpbestval );
	}

	RelaxResult result;

	if( bestval - qpbestval  < 1.0 ){
		result = INFEASIBLE;
	}else{
		result = FEASIBLE;

		double *vals;

		vals = new double[m];

		for(int i=0; i<m; i++){
			vals[i] = round( relaxvals[i] );
		}

		double objval = compute_objval( vals );

		if( objval == qpbestval ){
			result = GETINTEGER;
		}
		
		SOLUTION	solution;
		solution.set_sol( m, vals, objval);
		pool.add_solution( solution );

		if( bestval > objval ){
			bestval = objval;
			bestsol = solution;
			result = UPDATE;
			Appfac = sqrt( bestval ) / _Appfac;
		}

		delete[] vals;
	}

	delete[] A;
	delete[] b;
	delete[] subQ;  
	delete[] subu;
	delete[] subl;
	delete[] subw;
	delete[] p;
	delete[] relaxvals;

	return result;
}

