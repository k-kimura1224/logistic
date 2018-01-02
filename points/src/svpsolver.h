#ifndef SVPSOLVER_H__
#define SVPSOLVER_H__

//#include <vector>
#include <list>
#include "probdata.h"
#include "solution.h"
#include "solution_pool.h"
#include "node.h"
#include "stopwatch.h"
#include "Schmidt_manager.h"
#include "cut_pool.h"

#define	NODESELECTION  2		// 0:depth fisrt search
										// 1: 0th
										// 2: 0th, sort

#define	BRANCHINGRULE_INT	3		// 0:
											// 1:
											// 2:
											// 3: ver 2.0
											// 4: ver 2.1
											// 5: ver 2.1

#define	CUT_OA	false

#define	HEUR_APP		0.95

#define	TIMELIMIT	86400

// for parallel
#define	NUM_INITNODES		30000
#define	PARASIZE				100	

#define	PARA_LOG				true


using namespace std;

enum RelaxResult {
	UPDATE,
	FEASIBLE,
	INFEASIBLE,
	GETINTEGER
};

class SVPsolver{
	PROB_DATA		probdata;

	double			*ub;
	double			*lb;

	double			GLB;
	double			bestval;
	SOLUTION			bestsol;
	SOLUTION_POOL	pool;
	double			_Appfac;
	double			Appfac;

	list<NODE>		NodeList;
	int				listsize;

	STOPWATCH		stopwatch;
	unsigned long int		nnode;

	int				*order;		// order of norm or bounds

	SCHMIDT_M		sch;

	double QP_time;
	clock_t __start;
	double __time;

	double	*norm;

	CUT_POOL		oa_cpool;	//a pool of cutting planes by using
									//outer approximation

	// for parallel
	int	nthreads;
	public:

		SVPsolver();											// default constructor	
		SVPsolver( const SVPsolver &source );			// copy constructor
		SVPsolver& operator=( const SVPsolver& );		// assignment operator
		~SVPsolver();											// destructor

		void			create_probdata( int m, double *B_);
		void			create_sch( int m, double *B_);
		bool			solve(bool para, bool s_zero, double s_GLB, int tlimit);
		bool			p_solve();
		void			find_min_column();
		void			compute_bounds();
		int			select_node(int index, int disp);
		RelaxResult	solve_relaxation(int sel);
		RelaxResult	solve_relaxation_BIN(int sel);
		RelaxResult	solve_relaxation_BIN_sch(int sel);
		RelaxResult	solve_relaxation_INT(int sel);
		void			heur(int sel, bool para);
		void			heur_unitsphere(int sel);
		void			heur_quadratic(int sel);
		void			disp_log(int	sel, RelaxResult r, int index, int cutoff);
		void			disp_bestsol();
		double		compute_objval( double *x );
		void			branch(int	sel, int index);
		void			branch_BIN(int	sel, int index);
		void			branch_INT(int	sel, int index);
		double		get_bestval(){ return bestval; }
		SOLUTION		get_bestsol(){ return bestsol; }
		PROB_DATA	get_probdata(){ return probdata; }
		unsigned long int	get_nnode(){ return nnode; }
		void			gene_OAcuts( double *u, double *l, double *Q, double M);

		// for parallel mode
		void			set_num_thread(int n){ nthreads = n; }
		void			set_bestval(double val){ bestval = val; }
		void			set_bounds(double *u, double *l){
			for(int i=0; i<probdata.get_m(); i++){
				ub[i] = u[i];
				lb[i] = l[i];
			}
		}
		void			set_order(int *s_order){
			order = new int[probdata.get_m()];
			for(int i=0; i<probdata.get_m(); i++) order[i] = s_order[i];
		}
		void			set_norm(double *s_norm){
			assert( probdata.get_m() > 0 );
			assert( s_norm != NULL );
			assert( norm == NULL );
			norm = new double[probdata.get_m()];
			Copy_vec( s_norm, norm, probdata.get_m());
		}
		void			set_Appfac(double a, double _a){
			Appfac = a;
			_Appfac = _a;
		}
};

#endif
