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

using namespace std;

#define debug	0

static int find_branchingvariable1(
		int		m,
		double 	*relax_solval
	)
{

	double	*eval;
	eval = new double[m];

	for(int i=0; i<m; i++){
		eval[i] = pow( relax_solval[i] - round(relax_solval[i]), 2);
	}

	double	max_eval = 0.0;	
	int		memo = -1;
	for(int i=0; i<m; i++){
		if( max_eval <= eval[i] ){
			max_eval = eval[i];
			memo = i;
		}
	}

	delete[] eval;
	return memo;
}

static int find_branchingvariable2(
		int		m,
		double 	*relax_solval,
		double	*ub,
		double	*lb
	)
{

	double	min = 1.0;	
	int		memo = -1;
	double	r;
	for(int i=0; i<m; i++){
		if( ub[i] != lb[i] 
			&& pow(relax_solval[i],2) < 1.0
			&& relax_solval[i] != 0.0 ){
			if( relax_solval[i] >  0.0 ){
				r = relax_solval[i];
			}else{
				r = -relax_solval[i];
			}
			if( min > r ){
				min = r;
				memo = i;
			}
		}
	}

	return memo;
}

static int find_branchingvariable3(
		int		m,
		double 	*relax_solval,
		double	*ub,
		double	*lb
	)
{

	for(int i=m-1; 0<=i; i--){
		if(  pow(relax_solval[i],2) < 1.0
			&& relax_solval[i] != 0.0 ){
			return i;
		}
	}

	return -1;
}

static int find_branchingvariable4(
		int		m,
		double 	*relax_solval
	)
{

	int memo;
	
	for(int i=m-1; 0<=i; i--){
		if( relax_solval[i] - round(relax_solval[i]) != 0){
			memo = i;
			break;
		}
	}

	return memo;
}

static int find_branchingvariable5(
		int		m,
		double 	*relax_solval,
		int		*order
	)
{

	assert( order != NULL );

	int memo = -1;
	for(int i=0; i<m; i++){
		assert( order[i] >= 0 && order[i] < m );
		if( relax_solval[order[i]] - round(relax_solval[order[i]]) != 0){
			memo = order[i];
			//cout << memo << endl;
			break;
		}
	}

	assert( memo >= 0 );
	assert( memo < m );

	return memo;
}

static int find_branchingvariable6(
		int		m,
		double 	*relax_solval,
		int		*order
	)
{

	assert( order != NULL );

	int memo = -1;
	for(int i=m-1; 0<=i; i--){
		assert( order[i] >= 0 && order[i] < m );
		if( relax_solval[order[i]] - round(relax_solval[order[i]]) != 0){
			memo = order[i];
			//cout << memo << endl;
			break;
		}
	}

	assert( memo >= 0 );
	assert( memo < m );

	return memo;
}

void	SVPsolver::branch(
	int	sel,
	int	index
	)
{
	list<NODE>::iterator it = NodeList.begin();
	advance( it, sel);
	if( it->get_zero() == true ){
		branch_BIN(sel, index);
	}else{
		branch_INT(sel, index);
	}
}

void	SVPsolver::branch_BIN(
	int	sel,
	int	index
	)
{
	//cout << "branch_start" << endl;
	NODE	C_BIN;
	NODE	C_INT;

	int		m = probdata.get_m();
	double	*set_ub;
	double	*set_lb;
	double	*set_warm;
	double	set_relax_objval;
	int		set_dpt;
	bool		set_zero;
	int		set_index;

	list<NODE>::iterator it = NodeList.begin();
	advance( it, sel);
	double	*warm = it->get_warm();

	set_ub = new double[m];
	set_lb = new double[m];

	// x_i = 0	
	int memo = -1;
	Copy_vec( it->get_ub(), set_ub, m);
	Copy_vec( it->get_lb(), set_lb, m);
	for(int i=0; i<m; i++){
		if( set_lb[i] != 0.0 || set_ub[i] != 0.0 ){
			set_lb[i] = 0.0;
			set_ub[i] = 0.0;
			memo = i;
			break;
		}
	}

	set_warm = it->get_warm();
	for(int i=0; i<m; i++){
		if( warm[i] >= set_lb[i] && warm[i] <= set_ub[i] ){
			set_warm[i] = warm[i];
		}else{
			set_warm[i] = (set_ub[i]+set_lb[i])/2.0;
		}
	}

	set_relax_objval = it->get_lowerbound();
	set_dpt = it->get_dpt() + 1;
	set_zero = true;
	set_index = index;

	C_BIN.set_vals( m, set_ub, set_lb,
						 set_warm, set_relax_objval,
						 set_dpt,  set_zero, set_index);

	if( it->get_sumfixed() != NULL ){
		bool r = C_BIN.alloc_sumfixed();
		assert( r == true );
		C_BIN.set_sumfixed( 1.0, it->get_sumfixed() );
	}

	if( set_lb[memo] == set_ub[memo] && set_lb[memo] != 0.0 ){
		if( C_BIN.alloc_sumfixed() == true ){
			C_BIN.set_sumfixed( set_lb[memo], probdata.get_bvec(memo) );
		}else{
			C_BIN.add_sumfixed( set_lb[memo], probdata.get_bvec(memo) );
		}
	}

	NodeList.push_back( C_BIN );
	listsize++;

	assert( memo >= 0 );
	// x_i >= 1
	Copy_vec( it->get_ub(), set_ub, m);
	Copy_vec( it->get_lb(), set_lb, m);
	for(int i=0; i<m; i++){
			set_lb[memo] = 1.0;
	}

	set_warm = it->get_warm();
	for(int i=0; i<m; i++){
		if( warm[i] >= set_lb[i] && warm[i] <= set_ub[i] ){
			set_warm[i] = warm[i];
		}else{
			set_warm[i] = set_lb[i];
		}
	}

	set_relax_objval = it->get_lowerbound();
	set_dpt = it->get_dpt() + 1;
	set_zero = false;
	set_index = index+1;

	C_INT.set_vals( m, set_ub, set_lb,
						set_warm, set_relax_objval,
						set_dpt,  set_zero, set_index);


	if( it->get_sumfixed() != NULL ){
		bool r = C_INT.alloc_sumfixed();
		assert( r == true );
		C_INT.set_sumfixed( 1.0, it->get_sumfixed() );
	}

	if( set_lb[memo] == set_ub[memo] && set_lb[memo] != 0.0 ){
		if( C_INT.alloc_sumfixed() == true ){
			C_INT.set_sumfixed( set_lb[memo], probdata.get_bvec(memo) );
		}else{
			C_INT.add_sumfixed( set_lb[memo], probdata.get_bvec(memo) );
		}
	}
	NodeList.push_back( C_INT );
	listsize++;

	delete[] set_ub;
	delete[] set_lb;
	//cout << "branch_end" << endl;
	//exit(1);

}
void	SVPsolver::branch_INT(
	int	sel,
	int	index
	)
{
	NODE	C_LEFT;
	NODE	C_RIGHT;

	int		m = probdata.get_m();
	double	*set_ub;
	double	*set_lb;
	double	*set_warm;
	double	set_relax_objval;
	int		set_dpt;
	bool		set_zero;
	int		set_index;

	list<NODE>::iterator it = NodeList.begin();
	advance( it, sel);
	double	*relax_solval = it->get_relaxsolval();

	set_ub = new double[m];
	set_lb = new double[m];
	set_warm = new double[m];

	assert( relax_solval != NULL );
	// find a branching variable
	int memo = -1;
	switch( BRANCHINGRULE_INT ){
		case 0:
		{
			memo = find_branchingvariable1(m, relax_solval);
			break;
		}
		case 1:
		{
			memo = find_branchingvariable2(m, relax_solval, it->get_ub(), it->get_lb());
			if( memo == -1 )
				memo = find_branchingvariable1(m, relax_solval);
			break;
		}
		case 2:
		{
			memo = find_branchingvariable3(m, relax_solval, it->get_ub(), it->get_lb());
			if( memo == -1 )
				memo = find_branchingvariable1(m, relax_solval);
			break;
		}
		case 3:
		{
			//printv( m, relax_solval);
			memo = find_branchingvariable4(m, relax_solval);
			break;
		}
		case 4:
		{
			memo = find_branchingvariable5( m, relax_solval, order);
			break;
		}
		case 5:
		{
			memo = find_branchingvariable6( m, relax_solval, order);
			break;
		}

		default:
		{
			cout << "error: branch.cpp" << endl;
			exit(-1);
			break;
		}
	}

	assert( memo >= 0 );
	assert( memo < m );

	// x_memo <= floor
	Copy_vec( it->get_ub(), set_ub, m);
	Copy_vec( it->get_lb(), set_lb, m);
	for(int i=0; i<m; i++){
		if( i == memo ){
			set_ub[i] = floor( relax_solval[i] );
		}
	}

	Copy_vec( relax_solval, set_warm, m);
	for(int i=0; i<m; i++){
		if( i== memo ){
			set_warm[i] = set_ub[i];
		}
	}


	set_relax_objval = it->get_lowerbound();
	set_dpt = it->get_dpt() + 1;
	set_zero = false;
	set_index = index;

	C_LEFT.set_vals( m, set_ub, set_lb,
						 set_warm, set_relax_objval,
						 set_dpt,  set_zero, set_index);

	if( it->get_sumfixed() != NULL ){
		bool r = C_LEFT.alloc_sumfixed();
		assert( r == true );
		C_LEFT.set_sumfixed( 1.0, it->get_sumfixed() );
	}

	if( set_lb[memo] == set_ub[memo] && set_lb[memo] != 0.0 ){
		if( C_LEFT.alloc_sumfixed() == true ){
			C_LEFT.set_sumfixed( set_lb[memo], probdata.get_bvec(memo) );
		}else{
			C_LEFT.add_sumfixed( set_lb[memo], probdata.get_bvec(memo) );
		}
	}
	
	NodeList.push_back( C_LEFT );
	listsize++;

	// ceil <= x_memo
	Copy_vec( it->get_ub(), set_ub, m);
	Copy_vec( it->get_lb(), set_lb, m);
	for(int i=0; i<m; i++){
		if( i == memo ){
			set_lb[i] = ceil( relax_solval[i] );
		}
	}

	Copy_vec( relax_solval, set_warm, m);
	for(int i=0; i<m; i++){
		if( i == memo ){
			set_warm[i] = set_lb[i];
		}
	}

	set_relax_objval = it->get_lowerbound();
	set_dpt = it->get_dpt() + 1;
	set_zero = false;
	set_index = index+1;

	C_RIGHT.set_vals( m, set_ub, set_lb,
						 set_warm, set_relax_objval,
						 set_dpt,  set_zero, set_index);

	if( it->get_sumfixed() != NULL ){
		bool r = C_RIGHT.alloc_sumfixed();
		assert( r == true );
		C_RIGHT.set_sumfixed( 1.0, it->get_sumfixed() );
	}

	if( set_lb[memo] == set_ub[memo] && set_lb[memo] != 0.0 ){
		if( C_RIGHT.alloc_sumfixed() == true ){
			C_RIGHT.set_sumfixed( set_lb[memo], probdata.get_bvec(memo) );
		}else{
			C_RIGHT.add_sumfixed( set_lb[memo], probdata.get_bvec(memo) );
		}
	}

	NodeList.push_back( C_RIGHT );
	//cout << "branching: " << memo << ", [" << floor(relax_solval[memo]) << ", " << ceil(relax_solval[memo]) << "]" << endl;
	listsize++;


	delete[] set_ub;
	delete[] set_lb;
	delete[] set_warm;
}


