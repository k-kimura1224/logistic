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

#include "svpsolver.h"
#include "solution.h"
#include "solution_pool.h"
#include "node.h"
#include "stopwatch.h"
#include "vector.h"
//#include "probdata.h"

#define debug 0

bool SVPsolver::p_solve()
{
	assert( nthreads > 0 );

	stopwatch.start();

	int m = probdata.get_m();
	assert( m > 0 );

	if( BRANCHINGRULE_INT == 5 ){
		double	*_val = new double[m];

		for(int i=0; i<m; i++){
			_val[i] = ub[i] - lb[i];
		}

		assert( order == NULL );
		order = new int[m];
		double *_null = NULL;
		Bubble_Sort( m, _val, _null, order);
		printv( m, order);

		delete[] _val;
		_val = NULL;
	}

	// output bounds
	cout << "Bounds: " << endl;
	for(int i=0; i<m; i++){
		cout << "x_" << i << ": [ " << lb[i] << ", " << ub[i] << "]" << endl; 	
	}

	// generate a root node
	NODE	root;
	int	index=0;
	root.set_vals( m, ub, lb, bestsol.get_solval(), GLB, 0, true, index);

	for(int i=0; i<m; i++){
		if( lb[i] - ub[i] == 0 && lb[i]!=0.0 ){
			if( root.alloc_sumfixed() == true ){
				root.set_sumfixed( lb[i], probdata.get_bvec(i) );
			}else{
				root.add_sumfixed( lb[i], probdata.get_bvec(i) );
			}
		}
	}

	NodeList.push_back( root );
	listsize++;
	index++;

	// generate oa_cpool
	if( CUT_OA == true ){
		gene_OAcuts( ub, lb, probdata.get_Q(), bestval);
	}

	int			sel;
	RelaxResult	r;
	QP_time  = 0.0;
	__time = 0.0;
	__start = clock();
	int	disp = index;
	int	cutoff=0;
	bool 	result;
	while(1){
		assert( (int)NodeList.size() > 0 );
		assert( listsize > 0 );
		assert( (int)NodeList.size() == listsize );
		
		if( listsize >= NUM_INITNODES ){
			result = false;
			break;
		}
		// select a node from the list
		sel = select_node(index, disp);

		assert( sel >= 0 );
		assert( sel < listsize );

		// solve a relaxation problem
		r = solve_relaxation(sel);

		if( r == INFEASIBLE || r == GETINTEGER ){
			cutoff++;
		}

		// run heuristics
		if( r == FEASIBLE && HEUR_APP < Appfac ){
			heur( sel, false);
		}

		// output
		if( (index-1) % 1000 == 0 && disp < index ){
			list<NODE>::iterator it = NodeList.begin();
			double min_lb = it->get_lowerbound();
			for(int i=1; i<listsize; i++){
				++it;
				if( min_lb > it->get_lowerbound() ){
					min_lb = it->get_lowerbound();
				}
			}
			GLB = min_lb;
			disp_log(sel, r, index, cutoff);
			disp = index;
			cutoff = 0;
		}

		// branch
		if( r == UPDATE || r == FEASIBLE ){
			branch( sel, index);
			index += 2;
		}

		// remove
		clock_t  start = clock();
		list<NODE>::iterator it = NodeList.begin();
		advance( it, sel);
		//NodeList.erase( NodeList.begin() + sel);
		NodeList.erase( it );
		listsize--;

		clock_t  end = clock();
		__time += (double)(end-start)/CLOCKS_PER_SEC;
		// break
		assert( (int)NodeList.size() == listsize );
		if( stopwatch.check_time() == false ){
			result = true;
			break;
		}
		if( listsize == 0 ){
			cout << "End" << endl;
			result = true;
			break;
		}

	}

	nnode = (unsigned long int)index -1;

	if( result == true ){
		stopwatch.stop();
		return true;
	}

	// parallel {{
	NodeList.sort();

// LLL n43 seed2
//	for(int i=0; i<listsize; i++){
//		list<NODE>::iterator it = NodeList.begin();
//		advance( it, i);
//
//		if( i==2541 ){
//			cout << i << "->" << it->get_lowerbound() << endl;
//			cout << "debug" << endl;
//			SVPsolver sub;
//			sub.create_probdata( m, probdata.get_B_());
//			if( it->get_zero() == true ) sub.create_sch( m, probdata.get_B_());
//			sub.set_bestval( bestval );
//			sub.set_bounds( it->get_ub(), it->get_lb());
//			if( BRANCHINGRULE_INT == 4 || BRANCHINGRULE_INT == 5 ) sub.set_order( order );
//			sub.set_norm( norm );
//			sub.set_Appfac( Appfac, _Appfac);
//			double testrun = sub.solve(true, it->get_zero(), it->get_lowerbound(), TIMELIMIT - stopwatch.get_time() );
//		}
//	}
//	exit(1);

	cout << "Parallel mode ---------------------------------------------";
	cout << "-------------------------------------------" << endl;
	assert( (int)NodeList.size() == listsize );

	bool			p_run;	
	double		p_min;

#pragma omp parallel num_threads(nthreads)
	{
	#pragma omp for schedule(dynamic, PARASIZE) private(p_run,p_min)
	for(int i=0; i<listsize; i++){
		//cout << "[debug|" << " th:" << omp_get_thread_num() << " i:" << i << "]" <<endl ;
		list<NODE>::iterator it = NodeList.begin();
		advance( it, i);
		SVPsolver sub;
		sub.create_probdata( m, probdata.get_B_());
		if( it->get_zero() == true ) sub.create_sch( m, probdata.get_B_());
		sub.set_bestval( bestval );
		sub.set_bounds( it->get_ub(), it->get_lb());
		if( BRANCHINGRULE_INT == 4 || BRANCHINGRULE_INT == 5 ) sub.set_order( order );
		sub.set_norm( norm );
		sub.set_Appfac( Appfac, _Appfac);
		p_run = sub.solve(true, it->get_zero(), it->get_lowerbound(), TIMELIMIT - stopwatch.get_time() );

		if( !p_run ){
			#pragma omp critical
			{
				nnode += sub.get_nnode();
			}
		}else{

			it->set_solved( true );
	
			int p_ct = 0;
			if( i%(int)PARASIZE == 0 ){
				p_min = bestval;
				list<NODE>::iterator itr = NodeList.begin();
				for(int j=0; j<listsize; ++j,++itr) {
					if( itr->get_solved() == true ){
						p_ct++;
					}
					if( itr->get_solved() == false && p_min > itr->get_lowerbound() ){
						p_min = itr->get_lowerbound();
					}
				}
				#pragma omp critical
				{
					GLB = p_min;
				}
			}
	
			#pragma omp critical
			{
				if( bestval > sub.get_bestval() ){
					bestval = sub.get_bestval();
					bestsol = sub.get_bestsol();
					pool.add_solution( bestsol );
					Appfac = sqrt( bestval ) / _Appfac;
				}
	
				nnode += sub.get_nnode();
	
				if( i%(int)PARASIZE == 0 ){
					cout << "*";
					cout << stopwatch.get_time() << "s: ";
					cout << "GLB=" << GLB << ", ";
					cout << "best=" << bestval << ", ";
					cout << "gap=" << 100*(bestval - GLB)/bestval << "%, ";
					cout << "AF=" << Appfac << ", ";
					cout << "#solved=" << p_ct;
					cout << endl;
				}
	
			}
		}
	}
	}
	cout << "End" << endl;
	// }} parallel
	
	double	min = bestval;
	list<NODE>::iterator itr = NodeList.begin();
	for(int j=0; j<listsize; ++j,++itr) {
		if( itr->get_solved() == false && min > itr->get_lowerbound() ){
			min = itr->get_lowerbound();
		}
	}
	GLB = min;

	if( GLB == bestval ){
		NodeList.clear();
		listsize = 0;
	}

	stopwatch.stop();
	return true;

}
