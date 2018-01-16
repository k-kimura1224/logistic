
/**
 *
**/
#include <map>
#include <vector>
#include <iostream>
#include <algorithm>		// copy
#include <math.h>			// pow
#include <iterator>		// back_inserter

#include "mysolpool.h"
//#include "set_myparameter.h"

#define 	debug 0

using namespace std;


// MFB_Solution::{{
MFB_Solution::MFB_Solution(){ // default constructor
#if debug
	cout << "MFB_Solution: default constructor" << endl;
#endif
	N = -1;
	solval_01 = NULL;
	objval = 1e+10;
}

MFB_Solution::MFB_Solution( const MFB_Solution &source )
{	// copy constructor
#if debug
	cout << "MFB_Solution: copy constructor" << endl;
#endif

	N	=	source.MFBSol_N();
	objval = source.objval;

	if( N > 0 ){

		int	i;
		solval_01 = new int[N];

		for(i=0; i<N; i++){
			solval_01[i] = source.solval_01[i];
		}

	}else{
		solval_01 = NULL;
	}
}

// assignment operator
MFB_Solution& MFB_Solution::operator=( const MFB_Solution& source )
{
#if debug
	cout << "MFB_Solution: assignment operator" << endl;
#endif

	if( this != &source ){
		N		=	source.MFBSol_N();
		objval = source.objval;

		if( N > 0 ){
			MFBSol_free();

			int	i;
			solval_01 = new int[N];

			for(i=0; i<N; i++){
				solval_01[i] = source.solval_01[i];
			}

		}else{
			solval_01 = NULL;
		}
	}
	return *this;
}

MFB_Solution::~MFB_Solution(){	// destructor
#if debug
	cout << "MFB_Solution: destructor" << endl;
#endif
	if( N > 0 ){
		MFBSol_free();
	}
}

int	MFB_Solution::MFBSol_N() const
{
	return N;
}

void	MFB_Solution::MFBSol_free()
{
#if debug
	cout << "MFB_Solution: free" << endl;
#endif

	assert( N>0 );
	delete[] solval_01;
}

void	MFB_Solution::MFBSol_initialize( int  n ){

	N = n;

	solval_01 = new int[N];

	int	i;
	for(i=0; i<N; i++){
		solval_01[i] = -1;
	}

}

void	MFB_Solution::MFBSol_set(
	int			*val_01,
	SCIP_Real	val
	)
{
	assert( N>0 );

	int	i;
	for(i=0; i<N; i++){

		assert( val_01[i]==1 || val_01[i]==0 );

		solval_01[i] = val_01[i];
	}

	objval = val;

}

// }} MFB_Solution

// MFB_SolPool::{{
MFB_SolPool::MFB_SolPool(){ // default constructor
#if debug
	cout << "MFB_SolPool: default constructor" << endl;
#endif
	N = -1;
	Nsols = 0;

	freq = NULL;
}

MFB_SolPool::MFB_SolPool( const MFB_SolPool &source )
{	// copy constructor
#if debug
	cout << "MFB_SolPool: copy constructor" << endl;
#endif

	N	=	source.MFBSolPool_N();
	Nsols = source.MFBSolPool_Nsols();

	if( N > 0 ){

		int	i;
		freq = new int[N];

		for(i=0; i<N; i++){
			freq[i] = source.freq[i];
		}

		for(i=0; i<100; i++){
			sols[i] = source.sols[i];
		}

	}else{
		freq = NULL;
	}
}

// assignment operator
MFB_SolPool& MFB_SolPool::operator=( const MFB_SolPool& source )
{
#if debug
	cout << "MFB_SolPool: assignment operator" << endl;
#endif

	if( this != &source ){
		N	=	source.MFBSolPool_N();
		Nsols = source.MFBSolPool_Nsols();

		if( N > 0 ){
			MFBSolPool_free();

			int	i;
			freq = new int[N];

			for(i=0; i<N; i++){
				freq[i] = source.freq[i];
			}

			for(i=0; i<100; i++){
				sols[i] = source.sols[i];
			}

		}else{
			freq = NULL;
		}
	}
	return *this;
}

MFB_SolPool::~MFB_SolPool(){	// destructor
#if debug
	cout << "MFB_SolPool: destructor" << endl;
#endif
	if( N > 0 ){
		MFBSolPool_free();
	}
}

int	MFB_SolPool::MFBSolPool_N() const
{
	return N;
}

int	MFB_SolPool::MFBSolPool_Nsols() const
{
	return Nsols;
}

void	MFB_SolPool::MFBSolPool_free()
{
#if debug
	cout << "MFB_SolPool: free" << endl;
#endif

	assert( N>0 );
	delete[] freq;
}

void	MFB_SolPool::MFBSolPool_initialize( int	n ){

	N = n;

	freq = new int[N];

	int i;
	for(i=0; i<N; i++){
		freq[i] = 0;
	}

	for(i=0; i<100; i++){
		sols[i].MFBSol_initialize( n );
	}
}

void	MFB_SolPool::MFBSolPool_store(
	int			*solval_01,
	SCIP_Real	mL
	)
{
	assert( N > 0 );

	int	store = -1;

	SCIP_Real	val;
	int			num,i;

	if( Nsols < 100 ){

		store = 1;

		num = 0;

		for(i=0; i<N; i++){
			assert( solval_01[i]==1 || solval_01[i]==0 );
			if( solval_01[i] == 1 ){
				num++;
			}
		}

		val = 2 * mL + (double)MP_PENAL * (double)num;
		//cout << "!!!!!!!!!!" <<endl;

	}else if( 2*mL >= sols[99].objval ){
		store = 0;
	}else{
		num = 0;

		for(i=0; i<N; i++){
			assert( solval_01[i]==1 || solval_01[i]==0 );
			if( solval_01[i] == 1 ){
				num++;
			}
		}

		val = 2 * mL + (double)MP_PENAL * (double)num;

		//cout << "val: " << val << ", 99: " << sols[99].objval << endl;
		if( val < sols[99].objval ){
			store = 1;
		}else{
			store = 0;
		}
	}

	assert( store==1 || store== 0 );

	if( store == 1 ){
		//cout << "store " << val << endl;
		MFBSolPool_update( solval_01, val );
		if( Nsols<100 ) Nsols++;
	}

}

void MFB_SolPool::MFBSolPool_update(
	int			*solval_01,
	SCIP_Real	val
	)
{
	assert( N>0 );

	int	i,memo;

	MFB_Solution	sol1,sol2;

	memo = -1;

	for(i=0; i<100; i++){
		if( val < sols[i].objval ){
			memo = i;
			break;
		}
	}

	assert( memo >= 0 );
	assert( memo < 100 );

	sol2 = sols[memo];

	// update
	sols[memo].MFBSol_set( solval_01, val);


	for(i=(memo+1); ( i<=Nsols ) && ( i<100 ); i++){
		sol1 = sol2;
		sol2 = sols[i];
		sols[i] = sol1;
	}

	for(i=0; i<N; i++){
		assert( solval_01[i]==0 || solval_01[i]==1 );
		if( solval_01[i] == 1 ){
			freq[i]++;
		}
	}

	if( Nsols==100 ){
		for(i=0; i<N; i++){
			assert( solval_01[i]==0 || solval_01[i]==1 );
			if( sol2.solval_01[i] == 1 ){
				freq[i]--;
			}
		}
	}
}
int	MFB_SolPool::MFBSolPool_maxfreq(
	int	*list
) const
{
	assert( N>0 );

	int	i;
	int	ind=-1;
	int	max=-1;

	for(i=0; i<N; i++){
		//cout << "freq_" << i << " = " << freq[i] << endl;
		assert( freq[i]>=0 );
		assert( freq[i]<=100 );
		assert( list[i]==0 || list[i]==1 );
		if( max < freq[i] && list[i]==1 ){
			max = freq[i];
			ind = i;
		}
	}

	assert( ind>=0 );
	assert( ind<N );

	return ind;

}
// }} MFB_SolPool



// Solution:: {{
Solution::Solution(){	// default constructor
#if debug
	cout << "Solution: default constructor" << endl;
#endif
	N		=	-1;
	ct		=	-1;
	val	=	NULL;
	key	=	NULL;
	val_mL = 0;
	n_key = -1;
}

Solution::Solution( int p1 ){	// int constructor
#if debug
	cout << "Solution: int constructor(allocation)" << endl;
#endif
	assert( p1 > 0 );
	N		=	p1;
	n_key = (int)ceil( (double)p1 / 20.0 );

	assert( n_key > 0 );

	val	=	new double[N];
	key	=	new int[n_key];
	ct		=	-1;
	val_mL = 0;
}

Solution::Solution( const Solution &source )
{	// copy constructor
#if debug
	cout << "Solution: copy constructor" << endl;
#endif

	N	=	source.MySol_N();
	ct	=	source.MySol_ct();
	val_mL = source.val_mL;
	n_key = source.MySol_Nkey();

	if( N > 0 ){

		assert( n_key > 0 );

		int	i;
		val =	new double[N];
		key = new int[n_key];

		for(i=0; i<N; i++){
			val[i] = source.val[i];
		}

		for(i=0; i<n_key; i++){
			key[i] = source.key[i];
		}
	}else{
		val = NULL;
		key = NULL;
	}
}

// assignment operator
Solution& Solution::operator=( const Solution& source )
{
#if debug
	cout << "Solution: assignment operator" << endl;
#endif

	if( this != &source ){
		N		=	source.MySol_N();
		ct		=	source.MySol_ct();
		val_mL = source.val_mL;
		n_key = source.MySol_Nkey();

		if( N > 0 ){
			MySol_free();

			assert( n_key > 0 );

			int	i;
			val = new double[N];
			key = new int[n_key];

			for(i=0; i<N; i++){
				val[i] = source.val[i];
			}

			for(i=0; i<n_key; i++){
				key[i] = source.key[i];
			}
		}else{
			val = NULL;
			key = NULL;
		}
	}
	return *this;
}

Solution::~Solution(){	// destructor
#if debug
	cout << "Solution: destructor" << endl;
#endif
	if( N > 0 ){
		MySol_free();
	}
}

void	Solution::MySol_set(
	SCIP_Real*	solval,
	SCIP_Real	set_mL,
	int*			k,
	int			c
	)
{
#if debug
	cout << "Solution: set" << endl;
#endif
	assert( N>0 && n_key>0 );

	int i;

	ct = c;

	for(i=0; i<N; i++){
		val[i] = solval[i];
	}

	val_mL = set_mL;

	for(i=0; i<n_key; i++){
		key[i] = k[i];
	}

}

void	Solution::MySol_info() const
{
	assert( N>0 );

	int i;

	cout << "dim: " << N << endl;
	cout << "ct: " << ct << endl;
	cout << "mL: " << val_mL << endl;
	cout << "n_key: " << n_key << endl;
	cout << "key: ";
	for(i=0; i<n_key; i++){
		cout << key[i] << " ";
	}

	cout << endl;

	for(i=0; i<N; i++){
		cout << "val[" << i << "]: " << val[i] << endl;
	}
	cout << endl;
}

void	Solution::MySol_free()
{
#if debug
	cout << "Solution: free" << endl;
#endif

	assert( N>0 );
	assert( n_key>0 );
	delete[] val;
	delete[] key;
}

int	Solution::MySol_N() const
{
	return N;
}

int	Solution::MySol_ct() const
{
	return ct;
}

int	Solution::MySol_Nkey() const{
	return n_key;
}

// }} Solution::

// SolPool:: {{
SolPool::SolPool(){	// default constructor
#if debug
	cout << "SolPool: default constructor" << endl;
#endif
	N = -1;
	max = MP_MAXPOOL;
	mode = MP_SOLPOOLMODE;
	size = 0;
	n_key = -1;
}

SolPool::SolPool(int p1){	// int constructor
#if debug
	cout << "SolPool: int constructor" << endl;
#endif
	assert( MP_MAXPOOL > 0 );
	N = p1;
	max = MP_MAXPOOL;
	mode = MP_SOLPOOLMODE;
	size = 0;
	n_key = (int)ceil( (double)p1 / 20.0 );

	if( MP_MFB_PRI > 0 ){
		MFB_pool.MFBSolPool_initialize( N );
	}
}

// copy constructor
SolPool::SolPool( const SolPool &source )
{
#if debug
	cout << "SolPool: copy constructor" << endl;
#endif

	N		=	source.MySolPool_N();
	max	=	source.MySolPool_maxsize();
	mode	=	source.MySolPool_mode();
	size	=	source.MySolPool_poolsize();
	n_key	=	source.MySolPool_Nkey();

	MFB_pool = source.MFB_pool;

	if( source.MySolPool_poolsize() > 0 ){
		copy( source.pool.begin(), source.pool.end(), back_inserter(pool));
	}
}

// assignment operator
SolPool& SolPool::operator=( const SolPool& source )
{
#if debug
	cout << "SolPool: assignment operator" << endl;
#endif

	if( this != &source ){
		N		=	source.MySolPool_N();
		max	=	source.MySolPool_maxsize();
		mode	=	source.MySolPool_mode();
		size	=	source.MySolPool_poolsize();
		n_key	=	source.MySolPool_Nkey();

		MFB_pool = source.MFB_pool;

		if( source.MySolPool_poolsize() > 0 ){
			copy( source.pool.begin(), source.pool.end(), back_inserter(pool));
		}
	}

	return *this;
}

// destructor
SolPool::~SolPool()
{
#if debug
	cout << "SolPool: destructor" << endl;
#endif

	if( size > 0 ){
		vector<Solution>().swap(pool);
	}
}

int	SolPool::MySolPool_N() const
{
	return	N;
}

int	SolPool::MySolPool_maxsize() const
{
	return	max;
}

int	SolPool::MySolPool_mode() const
{
	return	mode;
}

int	SolPool::MySolPool_poolsize() const
{
	return	size;
}

int	SolPool::MySolPool_Nkey() const
{
	return	n_key;
}

void	SolPool::MySolPool_updatesize(){
	size = (int)pool.size();
}

void	SolPool::MySolPool_genekey(
	int		*solval_01,		// [N]
	int		*key
	)
{
	assert( N > 0 );
	assert( n_key > 0 );

	int	i,j,ct1,ct2;

	ct2=0;

	for(i=0; i<n_key; i++){
		key [i] = 0;
		ct1=0;
		for(j=0; j<20; j++){

			key[i] += (int)pow( 2.0, (double)ct1 +1.0) * solval_01[ct2];

			//cout << "i,j=" << i<<"," << j << ", " << "key[0]=" << key[0];
			//cout << " key[1]=" << key[1] << endl;
			assert( key[i] >= 0);

			ct1++;
			ct2++;

			if( ct2 == N ){
				assert( i==(n_key -1) );
				break;
			}
		}
	}
}

int 	SolPool::MySolPool_check( int *key ) const
{
	int	result = -1;
	int	i,j;
	int	memo;

	for(i=size-1; i>=0; i--){

		assert( n_key == pool[i].MySol_Nkey() );

		memo = 0;

		for(j=0; j<n_key; j++){
			if( key[j] != pool[i].key[j] ){
				break;
			}else{
				memo++;
			}
		}

		if( memo == n_key ){
			result = i;
			break;
		}
	}

	return result;

}


// store solution
void	SolPool::MySolPool_store(
	int				*solval_01,
	SCIP_Real		*solval,
	SCIP_Real		mL
){
	assert( MP_MAXPOOL > 0 && max > 0 );
	assert( N > 0 );

	int	*key;

	key = new int[n_key];

	if( mode==1 ){
		MySolPool_genekey( solval_01, key);
		if( max > size && MySolPool_check( key ) == -1 ){
			Solution	sol(N);
			sol.MySol_set( solval, mL, key, 0);
			pool.push_back( sol );
			size++;

			if( MP_MFB_PRI > 0 ){
				MFB_pool.MFBSolPool_store( solval_01, mL);
			}
		}
	}

	delete[] key;
}

void	SolPool::MySolPool_info() const
{
	int	i;

	cout << "dim: " << N << endl;
	cout << "size: " << size << endl;
	cout << "max: " << max << endl;
	cout << "mode: " << mode << endl;
	cout << "n_key: " << n_key << endl;
	cout << endl;

	for(i=0; i<size; i++){
		pool[i].MySol_info();
		cout << endl;
	}
}


Solution*	SolPool::MySolPool_getsol(
	int				n
	)
{

	assert( n>=0 && n<size  );

	return &pool[n];

}

// }} SolPool::
