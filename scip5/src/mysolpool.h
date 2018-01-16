
/**@file     warmstart.h
 * @ingroup
 * @brief
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_MYSOLPOOL_H__
#define __SCIP_MYSOLPOOL_H__

#include <map>
#include <vector>

#include "scip/scip.h"
/**
 * When relaxator solves convex NLPs using the newton method,
 * relaxation solution of parent nodes is called
 * as initial point.
 *
**/

#define MP_PENAL 				2
#define MP_MAXPOOL		100000000	/* the maximal value of the number of solutions */
#define MP_SOLPOOLMODE	1			/* 1: normal mode
											 * 2:
											 * 3:
											 */

#if 1
#define MP_FB_PRI				200000	/* for branch_myfullstrong.cpp */
#define MP_MFB_PRI			0			/* for branch_mostfreq.cpp */
#elif 0
#define MP_FB_PRI				0
#define MP_MFB_PRI			200000
#else
#define MP_FB_PRI				-200000
#define MP_MFB_PRI			-200000
#endif

#ifdef __cplusplus
extern "C" {
#endif

using namespace std;

class MFB_Solution{
	int		N;

	public:
		int			*solval_01;
		SCIP_Real	objval;

		MFB_Solution();												// default constructor
		MFB_Solution( const MFB_Solution &source );			// copy constructor
		MFB_Solution& operator=( const MFB_Solution& );		// assignment operator
		~MFB_Solution();												// destructor

		int		MFBSol_N() const;
		void		MFBSol_free();

		void	MFBSol_initialize(int n);
		void	MFBSol_set( int*	val_01, SCIP_Real	val);
};

class MFB_SolPool{
	int		N;
	int		Nsols;

	public:
		MFB_Solution	sols[100];
		int				*freq;		// [N];

		MFB_SolPool();												// default constructor
		MFB_SolPool( const MFB_SolPool &source );			// copy constructor
		MFB_SolPool& operator=( const MFB_SolPool& );		// assignment operator
		~MFB_SolPool();												// destructor

		int			MFBSolPool_N() const;
		int			MFBSolPool_Nsols() const;
		void			MFBSolPool_free();

		void	MFBSolPool_initialize(int n);
		void	MFBSolPool_store(
						int	*solval_01,
						SCIP_Real	mL
					);
		void	MFBSolPool_update( int	*solval_01, SCIP_Real	val );
		int	MFBSolPool_maxfreq( int *list ) const; 	//return index

};

class Solution{
	int				N;
	int				n_key;
	int				ct;

	public:
		SCIP_Real		*val;		// [N]
		int				*key;		// [n_key]
		SCIP_Real		val_mL;

		Solution();											// default constructor
		Solution( int p1 );								// int constructor: alloc (N=p1)
		Solution( const Solution &source );			// copy constructor
		Solution& operator=( const Solution& );	// assignment operator
		~Solution();										// destructor

		void		MySol_set( SCIP_Real *solval, SCIP_Real set_mL, int *k, int c);
		void		MySol_info() const;
		void		MySol_free();

		int		MySol_N() const;
		int		MySol_key() const;
		int		MySol_ct() const;
		int		MySol_Nkey() const;
};

class SolPool{
	int					N;		// N = dimension of z
	int					max;	// the number of solutions
	int					mode;
	int					size;
	int					n_key;

	public:
		vector<Solution>	pool;
		MFB_SolPool			MFB_pool;

		SolPool();										// default constructor
		SolPool( int p1 );							// int constructor: (N=p1)
		SolPool( const SolPool &source );		// copy constructor
		SolPool& operator=( const SolPool& );	// assignment operator
		~SolPool();										// destructor

		int		MySolPool_N() const;
		int		MySolPool_maxsize() const;
		int		MySolPool_mode() const;
		int		MySolPool_poolsize() const;
		int		MySolPool_Nkey() const;
		void		MySolPool_genekey( int *solval_01, int *key );
		int		MySolPool_check( int *key ) const;
		void		MySolPool_updatesize();

		// store solution
		void		MySolPool_store(
						int	*solval_01,
						SCIP_Real	*solval,
						SCIP_Real	mL
					);
		void		MySolPool_info() const;

		Solution*	MySolPool_getsol( int n );		// get sol


};


#ifdef __cplusplus
}
#endif

#endif
