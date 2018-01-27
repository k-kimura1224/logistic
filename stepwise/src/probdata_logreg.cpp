
/**@file   probdata_logreg.c
 * @brief  Problem data for +***** problems
 * @author
 *
 *
 */

#include <iostream>
#include <vector>

#include <string.h>
#include <stdio.h>
#include <time.h>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include "probdata_logreg.h"
#include "convenient_tool.h"
#include "get_lineardependence.h"
#include "call_cblas.h"
#include "mysolpool.h"
#include "logistic_regression.h"
#include "stepsize.h"
#include "convenient_tool.h"
//#include "read_data.h"
//#include "normalization.h"
//#include "divide_data.h"
//#include "matrix.h"
//#include "vector.h"
//#include "cblapack.h"
//#include "mysolpool.h"
//#include "set_myparameter.h"

#define bigM				300.0
#define MP_WARM 1
#define debug 0
using namespace std;

struct SCIP_ProbData
{
	int				n;				/* the number of datas */
	int				p;				/*	the number of explain vars */
	SCIP_Real*		y;				/* [n] */
	SCIP_Real*		x;				/* [n*p1],   ColMajor,*/

	int				ndep;			/*	number of groups */
	int*				Mdep;			/* [ndep] max number in each groups */
	int*				groupX;		/* [ndep][p1] -> member of group[i]  */

	SCIP_Real*		coef_obj;	/* coefficients of obj.f */

	SCIP_VAR**		b;				/* [p+1] continuous variables */
	SCIP_VAR**		z;				/* [p+1] 01 variables */
	SCIP_VAR**		bx;			/* [n] continuous variables */
	SCIP_VAR**		EXP;			/* [n] continuous variables */
	SCIP_VAR**		LOG;			/* [n] continuous variables */
	int				nvars;		/* number of variables */

	SCIP_Real		penalcf;		/* penalty coefficient.
											if AIC, 2.
											if BIC, log(n). */
	SolPool*			pool;
};

/**@name Local methods
 *
 * @{
 */

/* check data */
static
SCIP_RETCODE checkData(
      int         n,
      int         p,
      SCIP_Real* data
      )
{
   assert( n > p );
   assert( p > 0 );

   if( !( n > p && p > 0 ) )
   {
      return SCIP_ERROR;
   }

   int i;
   int j;
   int p1 = p+1;
   SCIP_Real buf;

   for( i = 0; i < p1; i++ )
   {
      buf = data[i];

      for( j = 0; j < n; j++ )
      {
         if( fabs( buf - data[j*p1 + i] ) > 1.0e-6 )
            break;

         if( j == n - 1 )
         {
            cout << "error: colmun= " << i << endl;
            return SCIP_ERROR;
         }
      }
   }

   return SCIP_OKAY;
}


/** creates problem data */
static
SCIP_RETCODE probdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata           /**< pointer to problem data */
   )
{
   assert(scip != NULL);
   assert(probdata != NULL);

   /* allocate memory */
   SCIP_CALL( SCIPallocMemory(scip, probdata) );

   return SCIP_OKAY;
}

/** frees the memory of the given problem data */
static
SCIP_RETCODE probdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata            /**< pointer to problem data */
   )
{
	int	i;
	int	n = (*probdata)->n;
	int	p1 = (*probdata)->p+1;

	SCIPdebugMessage("probdataFree \n");
	assert(scip != NULL);
	assert(probdata != NULL);

	SCIPfreeMemoryArrayNull(scip, &(*probdata)->y);
	SCIPfreeMemoryArrayNull(scip, &(*probdata)->x);
	SCIPfreeMemoryArrayNull(scip, &(*probdata)->coef_obj);
	if( (*probdata)->ndep ){
		SCIPfreeMemoryArrayNull(scip, &(*probdata)->Mdep);
		SCIPfreeMemoryArrayNull(scip, &(*probdata)->groupX);
	}

	/* release variables */
	for(i=0; i<p1; i++){
		SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->b[i]) );
		SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->z[i]) );
	}
	for(i=0; i<n; i++){
		SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->bx[i]) );
		SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->EXP[i]) );
		SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->LOG[i]) );
	}

	SCIPfreeMemoryArrayNull(scip, &(*probdata)->b);
	SCIPfreeMemoryArrayNull(scip, &(*probdata)->z);
	SCIPfreeMemoryArrayNull(scip, &(*probdata)->bx);
	SCIPfreeMemoryArrayNull(scip, &(*probdata)->EXP);
	SCIPfreeMemoryArrayNull(scip, &(*probdata)->LOG);

	/* constraints */
	/*
	SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->defrss));
	SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->deflog));
	for(i=0; i<n; ++i) SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->defep[i]));
	for(i=0; i<p; ++i){
		SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->aM1[i]));
		SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->aM2[i]));
	}

	CIPfreeMemoryArrayNull(scip, &(*probdata)->defep);
	SCIPfreeMemoryArrayNull(scip, &(*probdata)->aM1);
	SCIPfreeMemoryArrayNull(scip, &(*probdata)->aM2);

	if( (*probdata)->ndep ){
		for(i=0; i<(*probdata)->ndep; i++){
			SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->ld[i]));
		}
		SCIPfreeMemoryArrayNull(scip, &(*probdata)->ld);
	}
	*/

	if( MP_MAXPOOL > 0 )
		delete ((*probdata)->pool);

	/* free probdata */
	SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}

/** forward selection */
static
SCIP_RETCODE forward(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,            /**< problem data */
	char*			probname	/*	output problemname	  */
   )
{

   int n        = SCIPprobdataGetNdatas( probdata );
   int p        = SCIPprobdataGetNexvars( probdata );
   int p1       = p + 1;
   int ndep     = SCIPprobdataGetNdep( probdata );
   SCIP_VAR** var_b    = SCIPprobdataGetVars_b( probdata );
   SCIP_VAR** var_z    = SCIPprobdataGetVars_z( probdata );
   SCIP_VAR** var_bx   = SCIPprobdataGetVars_bx( probdata );
   SCIP_VAR** var_EXP  = SCIPprobdataGetVars_EXP( probdata );
   SCIP_VAR** var_LOG  = SCIPprobdataGetVars_LOG( probdata );
   SCIP_Real penalcf  = SCIPprobdataGetPC( probdata );
   SCIP_Real* X        = SCIPprobdataGetX( probdata );
   SCIP_Real* coef_obj = SCIPprobdataGetCO( probdata );
   SolPool* pool     = SCIPprobdataGetPool( probdata );

   int*  Mdep;             /* [ndep] */
   int*  groupX;           /* [ndep*p] */

   int i,j,ct;

   /* for forward selection */
   int   dim;
   int*  list;             /* [p1] */

   int   t;
   int   memo;

   // for branching info
   int   *branchinfo;      /* [3*p1] */

   int   sum_branchinfo[3];

   if( ndep ){
      Mdep     =  SCIPprobdataGetMdep(probdata);
      groupX   =  SCIPprobdataGetgroupX(probdata);
   }else{
      Mdep     =  NULL;
      groupX   =  NULL;
   }


 // get branching info {{
   // alloc and initialize
   SCIP_CALL( SCIPallocBufferArray(scip, &branchinfo, 3*p1));
   SCIP_CALL( SCIPinitIntArrayZero( 3*p1, branchinfo ) );

   for(i=0; i<p1; ++i){
      branchinfo[p1 + i]  =  1;
   }
   // }} get branching info

   for( i = 0; i < 3; i++ )
      sum_branchinfo[i] = SCIPcalcIntSum( &branchinfo[p1*i], p1);

   clock_t start = clock();
   /* alloc */
   SCIP_CALL( SCIPallocBufferArray(scip, &list, p1));

   /* list */
   ct = 0;
   for( i = 0; i < p1; i++ )
   {
      if( branchinfo[i] == 1 )
      {
         list[i] = - 2;
      }
      else if( branchinfo[p1+i] == 1 )
      {
         list[i] = - 1;
      }
      else
      {
         list[i] = ct;
         ct++;
      }
   }

   dim = ct;

   // for forward selection
   SCIP_Real   AIC;
   SCIP_Real   AIC_new;
   SCIP_Real   mL;
   SCIP_Real   mL_new;
   SCIP_Real*  b;          /* [dim] */
   SCIP_Real*  b_old;      /* [dim] */

   int   *index_nonzero;
   int   *solval_01;
   SCIP_Real   *solval;
   int   *key;
   int   n_key = pool->MySolPool_Nkey();
   int   index;

   Solution*   mysol;

   // for the newton method
   SCIP_Real   *subX_;        // [n*dim], colMajor
   SCIP_Real   *subX_old_;    // [n*dim-1],  colMajor
   SCIP_Real   *subcoef;      // [dim],
   SCIP_Real   *subcoef_old;  // [dim-1],
   SCIP_Real   *A_;           // [dim*dim],colMajor
   SCIP_Real   *Y_;           // [n*dim],    colMajor
   SCIP_Real   *q;            // [dim],
   SCIP_Real   *pi;           // [n],
   SCIP_Real   *point;        // [dim], current point
   SCIP_Real   *d;            // [dim], direction vector
   SCIP_Real   *Xb;           // [n]

   SCIP_Real   objval;        // object value
   SCIP_Real   objval_new;    // new object value
   SCIP_Real   ss;            // step size
   SCIP_Real   newton_gap;
   SCIP_Real norm;

   int info;
   int buf_int;

   // alloc
   SCIP_CALL( SCIPallocBufferArray(scip, &solval_01, p1));
   SCIP_CALL( SCIPallocBufferArray(scip, &solval, p1));
   SCIP_CALL( SCIPallocBufferArray(scip, &key, n_key));
   SCIP_CALL( SCIPallocBufferArray(scip, &index_nonzero, p1));
   SCIP_CALL( SCIPallocBufferArray(scip, &pi, n));
   SCIP_CALL( SCIPallocBufferArray(scip, &Xb, n));

   for( i = 0; i < p1; i++ )
      index_nonzero[i] = -1;

   ct = 0;
   for( i = 0; i < p1; i++ )
   {
      if( list[i] >= 0 )
      {
         solval_01[i] = 1;
         index_nonzero[ct] = i;
         ct++;
      }
      else
      {
         solval_01[i] = 0;
      }
   }

   assert( dim >= 0 );
   assert( ct == dim );

#if debug
   cout << "index_nonzero:";
   SCIPprintIntVec( p1, index_nonzero);
   cout << "dim: " << dim << endl;
#endif

   if( dim >= 1 )
   {

      SCIP_CALL( SCIPallocBufferArray(scip, &b_old, dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &subX_old_, n*dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &subcoef_old, dim));

      if( MP_WARM == true )
      {
         pool->MySolPool_genekey( solval_01, key);
         index = pool->MySolPool_check( key );

         if( index > -1 )
         {
#if debug
   cout << "initial point is found from the pool" << endl;
#endif
            mysol = pool->MySolPool_getsol( index );

            for( i = 0; i < dim; i++ )
            {
               assert( index_nonzero[i] >= 0 );
               assert( index_nonzero[i] < p1 );

               b_old[i] = mysol->val[index_nonzero[i]];
            }

         }
         else
         {
            for( i = 0; i < dim; i++ )
               b_old[i] = 0.0;
         }

      }
      else
      {
         for( i = 0; i < dim; i++ )
            b_old[i] = 0.0;
      }

      ct = 0;
      for( i = 0; i < dim; i++ )
      {
         buf_int = index_nonzero[i] * n;
         for( j = 0; j < n; j++ )
            subX_old_[ct++]  = X[buf_int + j];
         subcoef_old[i] = coef_obj[index_nonzero[i]];
      }

      assert( ct == ( n * dim ) );

      for( j = 0; j < dim; j++ )
         subcoef_old[j] = coef_obj[index_nonzero[j]];
   }

   AIC = SCIPinfinity(scip);
   int ct_error;

   // forward selection
   while( 1 )
   {
      dim++;
      memo = -1;
      mL = SCIPinfinity(scip);

   cout << "[dim=" << dim << "] ";
#if debug
   cout << "[dim=" << dim << "] ";
   cout << "---------------------------" << endl;
#endif

      // alloc
      SCIP_CALL( SCIPallocBufferArray(scip, &d, dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &A_, dim*dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &Y_, n*dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &q, dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &subX_, n*dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &subcoef, dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &point, dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &b, dim));

      SCIP_CALL( SCIPinitArrayZero( dim*dim, A_) );

      // copy
      if( dim >= 2 )
      {
         assert( subX_old_ != NULL );
         assert( subcoef_old != NULL );

         SCIP_CALL( SCIPcblasCopy( subX_old_, subX_, n*(dim-1)) );
         SCIP_CALL( SCIPcblasCopy( subcoef_old, subcoef, dim-1) );

         SCIPfreeBufferArray(scip, &subX_old_);
         SCIPfreeBufferArray(scip, &subcoef_old);
      }

      SCIP_CALL( SCIPallocBufferArray(scip, &subX_old_, n*dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &subcoef_old, dim));

      for( i = 0; i < p1; i++ )
      {
#if debug
   cout << "b_" << i << " is added";
#endif
         if( list[i] == -1 )
         {
            if( ndep > 0 )
            {
               info = 0;
               for( j = 0; j < ndep; j++ )
               {
                  for( t = 0; t < p1; t++ )
                  {
                     buf_int = j * p1;
                     if( groupX[buf_int + t] == 1 )
                     {
                        if( list[t] < 0 && t != i )
                           break;
                     }

                     if( t == Mdep[j] )
                     {
                        info = 1;
                        break;
                     }
                  }

                  if( info == 1 )
                     break;

               }

               if( info == 1 )
                  continue;
            }

            for( j = 0; j < p1; j++ )
            {
               if( list[j] >= 0 || j == i )
                  solval_01[j] = 1;
               else
                  solval_01[j] = 0;
            }

            // newton method
            // step0: initialize
            // step0-1: generate sub matrix of X and sub vector of coef_obj
            SCIP_CALL( SCIPcblasCopy( &X[i*n], &subX_[n*(dim-1)], n) );
            subcoef[dim-1] = coef_obj[i];

            // step0-2: generate initial point
            info = 0;

            if( MP_WARM == 1 )
            {
               // 1. search the solution
               pool->MySolPool_genekey( solval_01, key);
               index = pool->MySolPool_check( key );

               if( index > -1 )
               {

                  mysol = pool->MySolPool_getsol( index );
                  mL_new = mysol->val_mL;
#if debug
   cout << " [the solution is found from the pool]";
   cout << ", mL: " << mL_new << endl;
#endif
                  if( mL > mL_new )
                  {
                     mL = mL_new;
                     memo = i;
                     // copy
                     for(j=0; j<dim-1; j++){
                        b[j] = mysol->val[index_nonzero[j]];
                     }
                     b[dim-1] = mysol->val[i];
                     SCIP_CALL( SCIPcblasCopy( subX_, subX_old_, n*dim) );
                     SCIP_CALL( SCIPcblasCopy( subcoef, subcoef_old, dim) );
                  }
                  continue;
               }

               if( dim == 1 )
               {
                  point[0] = 0.0;
                  info = 1;
               }
               else
               {
                  assert( dim >= 2 );
                  SCIP_CALL( SCIPcblasCopy( b_old, point, dim-1) );
                  point[dim-1] = 0.0;
                  info = 1;
               }
            }
            else
            {
               // then, warmstart is not used.
               for( j = 0; j < dim; j++ )
                  point[j] = 0.0;

               info = 1;
            }

            assert( info == 1 );

            // step0-3: set the objval
            SCIP_CALL( SCIPcblasDgemv2( subX_, n, dim, point, Xb) );
            objval = - Loglikelifood_( scip, n, dim, subcoef, Xb, point);
            assert( objval < SCIPinfinity(scip) );

            ct_error = 0;
            while( 1 )
            {
               // step1: find a descent direction by solving Ad = q
               // define pi
               for( j = 0; j < n; j++ )
               {
                  if( Xb[j] < 13.0 )
                     pi[j] = 1.0 - ( 1.0 / ( 1.0 + exp(Xb[j]) ) );
                  else
                     pi[j] = 1.0;
               }

               // define q:= c - X^t pi
               SCIP_CALL( SCIPcblasDgemv4( n, dim, subX_, pi, subcoef, -1.0, 1.0, q) );

               // define Y := P ( I_n - P ) X in ( n, dim)
               for( j = 0; j < n; j++ )
                  pi[j] *= 1 - pi[j];

               ct = 0;
               for( j = 0; j < dim; j++ )
               {
                  for( t = 0; t < n; t++ )
                  {
                     Y_[ct] = pi[t] * subX_[ct];
                     ct++;
                  }
               }

               assert( ct == n*dim );

               // define A:= X^t Y
               SCIP_CALL( SCIPcblasDgemm2( n, dim, dim, subX_, Y_, A_) );

               // solve Ad = q
               info = SCIPclapackDposv( scip, A_, q, dim, d);

               if( info != 0 )
               {
                  assert( ct_error == 1 || ct_error == 0 );
                  if( ct_error == 0 )
                  {
                     for( j = 0; j < dim; j++ )
                        point[j] = 0.0;

                     SCIP_CALL( SCIPcblasDgemv2( subX_, n, dim, point, Xb) );
                     objval = - Loglikelifood_( scip, n, dim, subcoef, Xb, point);
                     ct_error++;
                     continue;
                  }
                  else
                  {
                     exit(1);
                     info = SCIPclapackDgesv( scip, A_, q, dim, d);
                     if( info != 0 ) exit(1);
                  }
               }

               // step2: find the stepsize
               norm = SCIPcblasDnrm( d, dim );
               SCIP_CALL( SCIPcblasDscal( d, dim, 1.0/norm, d ) );
               ss = FindStepsize( scip, n, dim, subcoef, subX_, point, d);

               assert( ss >= 0 );

               if( ss <= 1e-05 )
                  ss = 0.0;

               // step3: update a new point and its objval
               SCIP_CALL( SCIPcblasDaxpy( d, point, dim, ss, 1.0, point) );

               SCIP_CALL( SCIPcblasDgemv2( subX_, n, dim, point, Xb) );
               objval_new = - Loglikelifood_( scip, n, dim, subcoef, Xb, point);

               assert( objval - objval_new >= - 1.0e-04 );

               newton_gap = pow( objval - objval_new, 2.0);

               if( newton_gap < 1e-08 ){
                  objval = objval_new;
                  break;
               }

               objval = objval_new;
            }// while

            // store to the pool
            if( MP_MAXPOOL > 0 )
            {
               for( j = 0; j < p1; j++ )
               {
                  if( i == j )
                     solval[j] = point[dim-1];
                  else if( list[j] < 0 )
                     solval[j] = 0.0;
                  else
                     solval[j] = point[list[j]];
               }
               pool->MySolPool_store( solval_01, solval, objval);
            }

            mL_new = objval;
#if debug
   cout << ", mL: " << mL_new;
#endif
            if( mL > mL_new )
            {
               mL = mL_new;
               memo = i;
               // copy
               SCIP_CALL( SCIPcblasCopy( point, b, dim) );
               SCIP_CALL( SCIPcblasCopy( subX_, subX_old_, n*dim) );
               SCIP_CALL( SCIPcblasCopy( subcoef, subcoef_old, dim) );
            }

         }
#if debug
   cout << endl;
#endif
      } // for

      assert( memo < p1 );

      if( memo == -1 )
      {
         // free
         SCIPfreeBufferArray(scip, &d);
         SCIPfreeBufferArray(scip, &A_);
         SCIPfreeBufferArray(scip, &Y_);
         SCIPfreeBufferArray(scip, &q);
         SCIPfreeBufferArray(scip, &subX_);
         SCIPfreeBufferArray(scip, &subcoef);
         SCIPfreeBufferArray(scip, &point);
         SCIPfreeBufferArray(scip, &b);

         SCIPfreeBufferArray(scip, &subX_old_);
         SCIPfreeBufferArray(scip, &subcoef_old);
#if debug
   cout << "break1" << endl;
#endif
         dim--;
         break;
      }

      AIC_new = 2 * mL + penalcf * (double)dim;

      if( AIC_new < AIC )
      {
         AIC = AIC_new;

         assert( list[memo] == -1 );
         list[memo] = dim-1;

         assert( index_nonzero[dim-1] == -1 );
         index_nonzero[dim-1] = memo;

         if( dim >= 2 )
            SCIPfreeBufferArray(scip, &b_old);

         SCIP_CALL(SCIPallocBufferArray(scip, &b_old, dim));
         SCIP_CALL( SCIPcblasCopy( b, b_old, dim) );

#if debug
   cout << "--->" << memo << "th variable, AIC: " << AIC << endl;
#endif
      }
      else
      {
         memo = -1;
#if debug
         cout << "---> no selection (AIC_new: " << AIC_new << ")" << endl;
#endif
      }

      // free
      SCIPfreeBufferArray(scip, &d);
      SCIPfreeBufferArray(scip, &A_);
      SCIPfreeBufferArray(scip, &Y_);
      SCIPfreeBufferArray(scip, &q);
      SCIPfreeBufferArray(scip, &subX_);
      SCIPfreeBufferArray(scip, &subcoef);
      SCIPfreeBufferArray(scip, &point);
      SCIPfreeBufferArray(scip, &b);

      if( memo == -1 )
      {
         dim--;
#if debug
   cout << "break2" << endl;
#endif
         SCIPfreeBufferArray(scip, &subX_old_);
         SCIPfreeBufferArray(scip, &subcoef_old);
         break;
      }
      else if( dim == sum_branchinfo[1] )
      {
         //dim--;
#if debug
   cout << "break3" << endl;
#endif
         SCIPfreeBufferArray(scip, &subX_old_);
         SCIPfreeBufferArray(scip, &subcoef_old);
         break;
      }
   }

   clock_t end = clock();
   int k = 0;
   {

      SCIP_SOL*   sol;
      SCIP_HEUR*  heur;
      int         nvars;
      SCIP_Real*  solvals;
      SCIP_VAR**  vars;
      //SCIP_Bool   success;

      nvars = SCIPprobdataGetNvars(probdata);
      assert( nvars == 2*p1 + 3*n );

      // alloc
      SCIP_CALL( SCIPallocBufferArray( scip, &solvals, nvars) );
      SCIP_CALL( SCIPallocBufferArray( scip, &vars, nvars) );

      ct = 0;

      // z
      for(i=0; i<p1; i++){
         vars[ct]    =  var_z[i];
         if( list[i] >= 0 ){
            solvals[ct] = 1.0;
            k++;
         }else{
            solvals[ct] = 0.0;
         }
         ct++;
      }

      // b
      for(i=0; i<p1; i++){
         vars[ct]    = var_b[i];
         if( list[i] < 0 ){
            solvals[ct] = 0.0;
            solval[i] = 0.0;
         }else{
            assert( list[i] >= 0 );
            assert( list[i] <= dim-1 );
            solvals[ct] = b_old[list[i]];
            solval[i] = b_old[list[i]];
         }
         ct++;
      }

      // bx
      for(i=0; i<n; i++){
         vars[ct]    =  var_bx[i];
         solvals[ct] =  0.0;

         for(j=0; j<p1; j++){
            solvals[ct] += solval[j] * X[i + j * n];
         }
         ct++;
      }

      // EXP
      for(i=0; i<n; i++){
         vars[ct]    =  var_EXP[i];
         solvals[ct] =  exp( solvals[ct-n] );
         ct++;
      }

      // LOG
      for(i=0; i<n; i++){
         vars[ct]    =  var_LOG[i];
         solvals[ct] =  log( 1.0 + solvals[ct-n] );
         ct++;
      }

      assert( ct == nvars );

      heur = SCIPfindHeur(scip, "trysol");
      SCIP_CALL( SCIPcreateSol(scip, &sol, heur));
      SCIP_CALL( SCIPsetSolVals(scip, sol, nvars, vars, solvals));

      FILE *file;
      char filename[SCIP_MAXSTRLEN];
      (void) SCIPsnprintf( filename, SCIP_MAXSTRLEN, "%s_f.sol", probname);
	   file = fopen(filename,"w");
      SCIP_CALL( SCIPprintSol(scip, sol, file, false));
	   fclose(file);
      //SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, FALSE, TRUE, TRUE, FALSE, &success) );

      SCIPfreeBufferArray(scip, &solvals);
      SCIPfreeBufferArray(scip, &vars);
   }

   cout << "Time of forward = " << (double)(end - start) / CLOCKS_PER_SEC << "sec.\n";
   cout << "AIC: " << AIC << endl;
   cout << "k: " << k << endl;

   SCIPfreeBufferArray(scip, &branchinfo);
   SCIPfreeBufferArray(scip, &list);
   SCIPfreeBufferArray(scip, &solval_01);
   SCIPfreeBufferArray(scip, &key);
   SCIPfreeBufferArray(scip, &index_nonzero);
   SCIPfreeBufferArray(scip, &pi);
   SCIPfreeBufferArray(scip, &Xb);
   SCIPfreeBufferArray(scip, &solval);
   SCIPfreeBufferArray(scip, &b_old);

   return SCIP_OKAY;
}

/** fix variable and check feasiblity by using linear dependence */
static
SCIP_Bool fixVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   p1,                 /**< the number of explanatory variables */
   int                   ndep,               /**< the number of linerly dependent sets */
   SCIP_PROBDATA*        probdata,           /**< user problem data */
   int*                  branchinfo          /**< branching information */
   )
{
   int* Mdep;             /* [ndep] */
   int* groupX;           /* [ndep*p1] */

   int i;
   int j;
   int buf;
   assert(scip != NULL);
   assert(p1 > 0);
   assert(ndep > 0);
   assert(probdata != NULL);
   assert(branchinfo != NULL);

   Mdep = SCIPprobdataGetMdep(probdata);
   groupX = SCIPprobdataGetgroupX(probdata);

   assert(Mdep != NULL);
   assert(groupX != NULL);

   // adjust fixed z from linear dependent {{
   bool ad = false;
   int   dummy;
   for( i = 0; i < ndep; ++i )
   {
      dummy = -1;
      buf = i * p1;
      for( j = 1; j < p1; ++j )
      {
         if( groupX[buf+j] == 1 )
         {
            if( branchinfo[j] == 1 )
               break;
            if( branchinfo[p1+j] == 1 )
               dummy = j;
            if( j == Mdep[i] )
            {
               if( dummy != -1 )
               {
                  branchinfo[p1+dummy] = 0;
                  branchinfo[dummy] = 1;
                  ad = true;
                  break;
               }
               else
               {
#if debug
                  cout << " cut-off!" << endl;
#endif
                  return FALSE;
               }
            }
         }
      }
   }

#if debug
   if( ad == true )
   {
      for(i=0; i<3; i++){
         for(j=0; j<p1; j++){
            printf("%d, ", *(branchinfo+(i*p1)+j));
         }
         cout << endl;
      }
   }
#endif
   // }} adjust fixed z from linear dependent

   return TRUE;
}
/** backward selection */
static
SCIP_RETCODE backward(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,            /**< problem data */
	char*			probname	/*	output problemname	  */
   )
{
   int   n;
   int   p;
   int   p1;
   int   ndep;
   SCIP_VAR**  var_b;      /* [p1] continuous variables */
   SCIP_VAR**  var_z;      /* [p1] 01 variables */
   SCIP_VAR**  var_bx;     /* [n] continuous variables */
   SCIP_VAR**  var_EXP;    /* [n] continuous variables */
   SCIP_VAR**  var_LOG;    /* [n] continuous variables */
   SCIP_Real   penalcf;
   SCIP_Real*  X;          /* [n*p1] */
   SCIP_Real*  coef_obj;   /* [n*p1] */
   SolPool*    pool;

   // for branching info
   int   *branchinfo;      /* [3*p1] */

   int   sum_branchinfo[3];

   int i,j,ct;

   /* for backward selection */
   int   dim;
   int*  list;             /* [p1] */
   SCIP_Real   AIC;
   SCIP_Real   AIC_new;
   SCIP_Real   mL;
   SCIP_Real   mL_new;
   SCIP_Real*  b;          /* [dim] */
   SCIP_Real*  b_old;      /* [dim] */

   int   t,ct_b;
   int   memo;

   assert(probdata != NULL);

   n        = SCIPprobdataGetNdatas( probdata );
   p        = SCIPprobdataGetNexvars( probdata );
   p1       = p + 1;
   ndep     = SCIPprobdataGetNdep( probdata );
   var_b    = SCIPprobdataGetVars_b( probdata );
   var_z    = SCIPprobdataGetVars_z( probdata );
   var_bx   = SCIPprobdataGetVars_bx( probdata );
   var_EXP  = SCIPprobdataGetVars_EXP( probdata );
   var_LOG  = SCIPprobdataGetVars_LOG( probdata );
   penalcf  = SCIPprobdataGetPC( probdata );
   X        = SCIPprobdataGetX( probdata );
   coef_obj = SCIPprobdataGetCO( probdata );
   pool     = SCIPprobdataGetPool( probdata );

   // get branching info {{
   // alloc and initialize
   SCIP_CALL( SCIPallocBufferArray(scip, &branchinfo, 3*p1));
   SCIP_CALL( SCIPinitIntArrayZero( 3*p1, branchinfo ) );

   for(i=0; i<p1; ++i){
      branchinfo[p1 +i ]  =  1;
   }
   // }} get branching info

   if( ndep )
   {
      /* fix non-fixed variable z and check feasiblity by using linear dependence */
      if( fixVariable(scip, p1, ndep, probdata, branchinfo) == FALSE )
      {
         /* cut-off */
         SCIPfreeBufferArray(scip, &branchinfo);
         return SCIP_OKAY;
      }
   }

   for(i=0; i<3; i++){
      sum_branchinfo[i] = SCIPcalcIntSum( &branchinfo[p1*i], p1);
   }

   clock_t start = clock();
   /* alloc */
   SCIP_CALL( SCIPallocBufferArray(scip, &list, p1));

   /* list */
   for(i=0; i<p1; i++){
      list[i] = 1 - branchinfo[i];
   }

   dim = sum_branchinfo[1] + sum_branchinfo[2];
   AIC = 1e+06;
   SCIP_CALL( SCIPallocBufferArray(scip, &b_old, dim));

   int   *index_nonzero;
   int   *solval_01;
   int   *val_01;
   SCIP_Real   *solval;
   SCIP_CALL( SCIPallocBufferArray(scip, &solval_01, p1));
   SCIP_CALL( SCIPallocBufferArray(scip, &solval, p1));
   SCIP_CALL( SCIPallocBufferArray(scip, &val_01, p1));

   // for the newton method
   SCIP_Real   *subX_;        // [n*dim], colMajor
   SCIP_Real   *subcoef;      // [dim],
   SCIP_Real   *A_;           // [dim*dim],colMajor
   SCIP_Real   *Y_;           // [n*dim],    colMajor
   SCIP_Real   *q;            // [dim],
   SCIP_Real   *pi;           // [n],
   SCIP_Real   *point;        // [dim], current point
   SCIP_Real   *d;            // [dim], direction vector
   SCIP_Real   *Xb;           // [n]

   SCIP_Real   objval;        // object value
   SCIP_Real   objval_new;    // new object value
   SCIP_Real   ss;            // step size
   SCIP_Real   newton_gap;
   SCIP_Real norm;

   int info;
   int buf_int;
   Solution*      mysol;

   SCIP_CALL( SCIPallocBufferArray(scip, &pi, n));
   SCIP_CALL( SCIPallocBufferArray(scip, &Xb, n));

   // backward selection
   // mark
   while(1)
   {
      dim--;
      memo = -1;
      mL = 1e+06;

   cout << "[dim=" << dim << "] ";
#if debug
   cout << "[dim=" << dim << "] ";
   cout << "-------------------------" << endl;
#endif

      // alloc
      SCIP_CALL( SCIPallocBufferArray(scip, &index_nonzero, dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &d, dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &A_, dim*dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &Y_, n*dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &q, dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &subX_, n*dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &subcoef, dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &point, dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &b, dim));

      SCIP_CALL( SCIPinitArrayZero( dim*dim, A_) );

      for(i=0; i<p1; i++){
#if debug
         cout << "b_" << i << "=0";
#endif

         if( (branchinfo[p1+i]==1) && (list[i]==1) )
         {
            ct = 0;
            for( j = 0; j < p1; j++ )
            {
               if( list[j] == 0 || j == i )
               {
                  solval_01[j] = 0;
               }
               else
               {
                  solval_01[j] = 1;
                  index_nonzero[ct++] = j;
               }
            }

            assert( ct==dim );

            // newton method
            // step0: initialize
            // step0-1: generate initial point
            info = 0;
            if( MP_WARM )
            {
               int *key;
               int n_key = pool->MySolPool_Nkey();
               int index;
               SCIP_CALL( SCIPallocBufferArray(scip, &key, n_key));
               // 1. search the solution
               pool->MySolPool_genekey( solval_01, key);
               index = pool->MySolPool_check( key );

               if( index > -1 ){
#if debug
                  cout << " [the solution is found from the pool]" << endl;
#endif
                  mysol = pool->MySolPool_getsol( index );
                  mL_new = mysol->val_mL;

#if debug
                  cout << " mL: " << mL_new << endl;
#endif
                  if( mL > mL_new ){
                     mL = mL_new;
                     memo = i;
                     // copy
                     for(j=0; j<dim; j++){
                        b[j] = mysol->val[index_nonzero[j]];
                     }
                  }
                  SCIPfreeBufferArray(scip, &key);
                  continue;
               }

               // 2. search
               for(j=0; j<p1; j++){
                  for(t=0; t<p1; t++){
                     val_01[t] = solval_01[t];
                  }

                  if( val_01[j] == 1 ){
                     val_01[j] = 0;
                  }else{
                     val_01[j] = 1;
                  }

                  pool->MySolPool_genekey( val_01, key);
                  index = pool->MySolPool_check( key );

                  if( index > -1 ){
#if debug
                     cout << " [generate the intial point]" << endl;
#endif
                     info = 1;
                     mysol = pool->MySolPool_getsol( index );

                     for(t=0; t<dim; t++){
                        point[t] = mysol->val[index_nonzero[t]];
                     }
                     break;
                  }
               }

               if( info == 0 ){
                  for(j=0; j<dim; j++){
                     point[j] = 0.0;
                  }
                  info = 1;
               }

               SCIPfreeBufferArray(scip, &key);
            }else{
               for(j=0; j<dim; j++){
                  point[j] = 0.0;
               }
               info = 1;

            }

            assert( info == 1 );

            // step0-2: generate sub matrix of X and sub vector of coef_obj
            ct = 0;
            for( j = 0; j < dim; j++ )
            {
               buf_int = index_nonzero[j] * n;
               for( t = 0; t < n; t++ )
                  subX_[ct++] = X[buf_int +t];

            }

            assert( ct == ( n * dim ));

            for( j = 0; j < dim; j++ )
               subcoef[j] = coef_obj[index_nonzero[j]];

            SCIP_CALL( SCIPcblasDgemv2( subX_, n, dim, point, Xb) );
            objval = - Loglikelifood_( scip, n, dim, subcoef, Xb, point);
            assert( objval < SCIPinfinity(scip) );

            while(1)
            {
               // step1: find a descent direction by solving Ad = q
               // define pi
               for( j = 0; j < n; j++ )
               {
                  if( Xb[j] < 13.0 )
                     pi[j] = 1.0 - ( 1.0 / ( 1.0 + exp(Xb[j]) ) );
                  else
                     pi[j] = 1.0;
               }

               // define q:= c - X^t pi
               SCIP_CALL( SCIPcblasDgemv4( n, dim, subX_, pi, subcoef, -1.0, 1.0, q) );

               // define Y:= P ( I_n - P ) X in ( n, dim)
               for( j = 0; j < n; j++ )
                  pi[j] *= 1 - pi[j];

               ct = 0;
               for( j = 0; j < dim; j++ )
               {
                  for( t = 0; t < n; t++ )
                  {
                     Y_[ct] = pi[t] * subX_[ct];
                     ct++;
                  }
               }
               assert( ct == n*dim );

               // define A:= X^t Y
               SCIP_CALL( SCIPcblasDgemm2( n, dim, dim, subX_, Y_, A_) );

               // solve Ad = q
               info = SCIPclapackDposv( scip, A_, q, dim, d);

               if( info != 0 )
               {
                  //mydcopy_( q, d, dim);
                  for( j = 0; j < dim; j++ )
                     point[j] = 0.0;

                  SCIP_CALL( SCIPcblasDgemv2( subX_, n, dim, point, Xb) );
                  objval = - Loglikelifood_( scip, n, dim, subcoef, Xb, point);
                  continue;

               }

               // step2: find the stepsize
               norm = SCIPcblasDnrm( d, dim );
               SCIP_CALL( SCIPcblasDscal( d, dim, 1.0/norm, d ) );
               ss = FindStepsize( scip, n, dim, subcoef, subX_, point, d);

               assert( ss >= 0 );

               if( ss <= 1e-05 ){
                  ss = 0.0;
               }

               // step3: update a new point and its objval
               SCIP_CALL( SCIPcblasDaxpy( d, point, dim, ss, 1.0, point) );

// test
               //for(j=0; j<dim; j++){
               //   if( fabs( point[j] ) > bigM ){
               //      //cout << "!!!" << endl;
               //      for(t=0; t<dim; t++){
               //         point[t] = 0.0;
               //      }
               //      break;
               //   }
               //}
// test

               SCIP_CALL( SCIPcblasDgemv2( subX_, n, dim, point, Xb) );
               objval_new = - Loglikelifood_( scip, n, dim, subcoef, Xb, point);

               assert( objval - objval_new >= - 1.0e-04 );

               newton_gap = pow( objval - objval_new, 2.0);

               if( newton_gap < 1e-08 ){
                  objval = objval_new;
                  break;
               }

               objval = objval_new;

            }//while

            if( MP_MAXPOOL > 0 ){
               ct=0;
               for(j=0; j<p1; j++){
                  if( solval_01[j] == 1 ){
                     solval[j] = point[ct];
                     ct++;
                  }else{
                     solval[j] = 0.0;
                  }
               }
               assert( ct==dim );
               pool->MySolPool_store( solval_01, solval, objval);
            }

// test
            //for(j=0; j<dim; j++){
            //   assert( fabs( point[j] ) <= bigM );
            //}

            mL_new = objval;
#if debug
            //printv( dim, point);
            cout << " mL: " << mL_new;
#endif
            if( mL > mL_new ){
               mL = mL_new;
               memo = i;
               // copy
               SCIP_CALL( SCIPcblasCopy( point, b, dim) );
            }

         }

#if debug
         cout << endl;
#endif
      }// for

      assert( memo < p1 );
      if( memo == -1 ){
         // free
         SCIPfreeBufferArray(scip, &index_nonzero);
         SCIPfreeBufferArray(scip, &d);
         SCIPfreeBufferArray(scip, &A_);
         SCIPfreeBufferArray(scip, &Y_);
         SCIPfreeBufferArray(scip, &q);
         SCIPfreeBufferArray(scip, &subX_);
         SCIPfreeBufferArray(scip, &subcoef);
         SCIPfreeBufferArray(scip, &point);
         SCIPfreeBufferArray(scip, &b);

#if debug
         cout << "break1" << endl;
#endif
         dim++;
         break;
      }

      AIC_new = 2 * mL + penalcf * (double)dim;

      if( AIC_new < AIC ){
         AIC = AIC_new;
         list[memo] = 0;
#if debug
         cout << "---> " << memo << "th variable, AIC:" << AIC << endl;
#endif

         SCIPfreeBufferArray(scip, &b_old);
         SCIP_CALL(SCIPallocBufferArray(scip, &b_old, dim));
         SCIP_CALL( SCIPcblasCopy( b, b_old, dim) );
      }else{
         memo = -1;
#if debug
         cout << "--> no selection (AIC_new: " << AIC_new << ")" << endl;
#endif
      }


      // free
      SCIPfreeBufferArray(scip, &index_nonzero);
      SCIPfreeBufferArray(scip, &d);
      SCIPfreeBufferArray(scip, &A_);
      SCIPfreeBufferArray(scip, &Y_);
      SCIPfreeBufferArray(scip, &q);
      SCIPfreeBufferArray(scip, &subX_);
      SCIPfreeBufferArray(scip, &subcoef);
      SCIPfreeBufferArray(scip, &point);
      SCIPfreeBufferArray(scip, &b);

      if( memo == -1 ){
         dim++;
#if debug
         cout << "break2" << endl;
#endif
         break;
      }else if( SCIPcalcIntSum( list, p1) == sum_branchinfo[2] ){
#if debug
         cout << "break3" << endl;
#endif
         dim++;
         break;
      }
   }

   clock_t end = clock();

   {

      SCIP_HEUR*  heur;
      SCIP_SOL*   sol;
      int         nvars;
      SCIP_Real*  solvals;
      SCIP_VAR**  vars;
      //SCIP_Bool   success;

      nvars = SCIPprobdataGetNvars(probdata);

      assert( nvars == 2*p1 + 3*n );

      // alloc
      SCIP_CALL( SCIPallocBufferArray( scip, &solvals, nvars) );
      SCIP_CALL( SCIPallocBufferArray( scip, &vars, nvars) );

      ct = 0;

      // z
      for(i=0; i<p1; i++){
         vars[ct]    =  var_z[i];
         solvals[ct] =  (double)list[i];
         ct++;
      }

      // b
      ct_b = 0;
      for(i=0; i<p1; i++){
         vars[ct]    = var_b[i];
         if( list[i] == 1 ){
            solvals[ct] = b_old[ct_b];
            solval[i] = b_old[ct_b];
            ct_b++;
         }else{
            solvals[ct] = 0.0;
            solval[i] = 0.0;
         }
         ct++;
      }

      // bx
      for(i=0; i<n; i++){
         vars[ct]    =  var_bx[i];
         solvals[ct] =  0.0;

         for(j=0; j<p1; j++){
            solvals[ct] += solval[j] * X[i + j * n];
         }
         ct++;
      }

      // EXP
      for(i=0; i<n; i++){
         vars[ct]    =  var_EXP[i];
         solvals[ct] =  exp( solvals[ct-n] );
         ct++;
      }

      // LOG
      for(i=0; i<n; i++){
         vars[ct]    =  var_LOG[i];
         solvals[ct] =  log( 1.0 + solvals[ct-n] );
         ct++;
      }

      assert( ct == nvars );

      heur = SCIPfindHeur(scip, "trysol");
      SCIP_CALL( SCIPcreateSol(scip, &sol, heur));
      SCIP_CALL( SCIPsetSolVals(scip, sol, nvars, vars, solvals));

      FILE *file;
      char filename[SCIP_MAXSTRLEN];
      (void) SCIPsnprintf( filename, SCIP_MAXSTRLEN, "%s_b.sol", probname);
	   file = fopen(filename,"w");
      SCIP_CALL( SCIPprintSol(scip, sol, file, false));
	   fclose(file);
      //SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, FALSE, TRUE, TRUE, FALSE, &success) );

      SCIPfreeBufferArray(scip, &solvals);
      SCIPfreeBufferArray(scip, &vars);

   }

   cout << "Time of backward = " << (double)(end - start) / CLOCKS_PER_SEC << "sec.\n";
   cout << "AIC: " << AIC << endl;
   cout << "k: " << ct_b << endl;

   // free
   SCIPfreeBufferArray(scip, &branchinfo);
   SCIPfreeBufferArray(scip, &solval_01);
   SCIPfreeBufferArray(scip, &solval);
   SCIPfreeBufferArray(scip, &val_01);
   SCIPfreeBufferArray(scip, &list);
   SCIPfreeBufferArray(scip, &b_old);
   SCIPfreeBufferArray(scip, &pi);
   SCIPfreeBufferArray(scip, &Xb);
   return SCIP_OKAY;
}

/** create constraints (in Flow or Price Mode) */
static
SCIP_RETCODE createConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,            /**< problem data */
	int*						 ng
   )
{
	int	i,j,ct;
	int	n = probdata->n;
	int	p = probdata->p;
	int	p1 = p+1;
	char	consname[SCIP_MAXSTRLEN];

	// create constraints
   SCIP_CONS** 	defbx;
   SCIP_CONS** 	defEXP;
   SCIP_CONS** 	defLOG;
   SCIP_CONS** 	bM1;
   SCIP_CONS** 	bM2;

	SCIP_Real	one = 1.0;
	SCIP_Real	minusone = -1.0;

	assert(scip != NULL);
	assert(probdata != NULL);
	assert( n >= 0 );
	assert( p >= 0 );

	SCIPdebugMessage("createConstraints \n");

	// alloc
	SCIP_CALL( SCIPallocMemoryArray( scip, &defbx, n));
	SCIP_CALL( SCIPallocMemoryArray( scip, &defEXP, n));
	SCIP_CALL( SCIPallocMemoryArray( scip, &defLOG, n));
	SCIP_CALL( SCIPallocMemoryArray( scip, &bM1, p1));
	SCIP_CALL( SCIPallocMemoryArray( scip, &bM2, p1));

 	/* linear constraint defbx
	 *		<b,x_i>  - bx[i] = 0
	 */
   {
      SCIP_Real* coef;

		SCIP_CALL( SCIPallocMemoryArray( scip, &coef, p1));

		for(i=0; i<n; ++i){
         (void) SCIPsnprintf( consname, SCIP_MAXSTRLEN, "def_bx_%d", i+1);
         for(j=0; j< p1; ++j){
            coef[j] = SCIPmatColMajor( probdata->x, n, i, j);
         }

         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &defbx[i], consname, p1, probdata->b, coef, 0.0, 0.0));
         SCIP_CALL( SCIPaddCoefLinear(scip, defbx[i], probdata->bx[i], minusone));
      }

		SCIPfreeMemoryArrayNull( scip, &coef);
   }

	/* nonlinear constraint defEXP
	 *		exp(bx[i]) - EXP[i] == 0.0
	 *		(? exp(bx[i]) <= EXP[i] )
	 */
   {

		for(i=0; i<n; ++i){
			SCIP_EXPR* bxexpr;
      	SCIP_EXPR* expr;
      	SCIP_VAR* var[1];
      	SCIP_VAR* var2[1];
      	SCIP_Real coef[1];
      	SCIP_EXPRTREE* exprtree;

      	// setup expression
      	SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &bxexpr, SCIP_EXPR_VARIDX, 0));

      	// expression for expr : exp(bx[i])
      	SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_EXP, bxexpr));

      	// expression tree from expr
      	SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, expr, 1, 0, NULL));
      	var[0] = probdata->bx[i];
      	SCIP_CALL( SCIPexprtreeSetVars(exprtree, 1, var));

      	// create nonlinear constraint for exprtree - EXP_btx_i  == 0.0
      	var2[0] = probdata->EXP[i];
      	coef[0] = -1.0;

			(void) SCIPsnprintf( consname, SCIP_MAXSTRLEN, "def_EXP_%d", i+1);
			//SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &defEXP[i], consname, 1, var2, coef, 1, &exprtree, &one, 0.0, 0.0));
			SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &defEXP[i], consname, 1, var2, coef, 1, &exprtree, &one, -SCIPinfinity(scip), 0.0));

			SCIP_CALL( SCIPexprtreeFree(&exprtree));
      }
   }

	/* nonlinear constraint def_LOG
	 * 	log(1+EXP[i]) - LOG[i] == 0.0
	 *		(? log(1+EXP[i]) <= LOG[i] )
	 */
   {

		for(i=0; i<n; ++i){
			SCIP_EXPR* expexpr;
      	SCIP_EXPR* expr;
      	SCIP_VAR* var[1];
      	SCIP_VAR* var2[1];
      	SCIP_Real coef[1];
      	SCIP_EXPRTREE* exprtree;

      	// setup expression
      	SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expexpr, SCIP_EXPR_VARIDX, 0));
      	SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), &expexpr, 1, &expexpr, &one, 1.0) );

      	// expression for expr : log(1+EXP[i])
      	SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_LOG, expexpr));

      	// expression tree from expr
      	SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, expr, 1, 0, NULL));
      	var[0] = probdata->EXP[i];
      	SCIP_CALL( SCIPexprtreeSetVars(exprtree, 1, var));

      	// create nonlinear constraint for exprtree - LOG[i]  == 0.0
      	var2[0] = probdata->LOG[i];
      	coef[0] = -1.0;

			(void) SCIPsnprintf( consname, SCIP_MAXSTRLEN, "def_LOG_%d", i+1);
			//SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &defLOG[i], consname, 1, var2, coef, 1, &exprtree, &one, 0.0, 0.0));
			SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &defLOG[i], consname, 1, var2, coef, 1, &exprtree, &one, -SCIPinfinity(scip), 0.0));

			SCIP_CALL( SCIPexprtreeFree(&exprtree));
      }
   }

   /* linear constraint bM1 : 0 <= b + Mz */
   {
      for(i=0; i<p1; ++i){
         (void) SCIPsnprintf( consname, SCIP_MAXSTRLEN, "bM1_%d", i);
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &bM1[i], consname, 0, NULL, NULL, 0.0, SCIPinfinity(scip)));
         SCIP_CALL( SCIPaddCoefLinear(scip, bM1[i], probdata->b[i], one));
         SCIP_CALL( SCIPaddCoefLinear(scip, bM1[i], probdata->z[i], bigM));
      }
   }

   /* linear constraint bM2 :  b - Mz <= 0 */
   {
      for(i=0; i<p1; ++i){
         (void) SCIPsnprintf( consname, SCIP_MAXSTRLEN, "bM2_%d", i);
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &bM2[i], consname, 0, NULL, NULL, -SCIPinfinity(scip), 0.0));
         SCIP_CALL( SCIPaddCoefLinear(scip, bM2[i], probdata->b[i], one));
         SCIP_CALL( SCIPaddCoefLinear(scip, bM2[i], probdata->z[i], -bigM));
      }
   }

 	// add constraints to problem
	for(i=0; i<p1; ++i){
      SCIP_CALL( SCIPaddCons(scip, bM1[i]) );
      SCIP_CALL( SCIPaddCons(scip, bM2[i]) );
   }

	for(i=0; i<n; ++i){
		SCIP_CALL( SCIPaddCons(scip, defbx[i]));
		SCIP_CALL( SCIPaddCons(scip, defEXP[i]));
		SCIP_CALL( SCIPaddCons(scip, defLOG[i]));
	}

	// release
   for (i=0; i<p1; ++i) {
	   SCIP_CALL( SCIPreleaseCons(scip, &bM1[i]) );
		SCIP_CALL( SCIPreleaseCons(scip, &bM2[i]) );
	}

	for(i=0; i<n; ++i){
		SCIP_CALL( SCIPreleaseCons(scip, &defbx[i]));
		SCIP_CALL( SCIPreleaseCons(scip, &defEXP[i]));
		SCIP_CALL( SCIPreleaseCons(scip, &defLOG[i]));
	}

	/**
	 * linearly dependent constraint
	 * ex.
	 *  dep = { 1, 2, 3, 4}
	 *  add constraint:
	 *  "ld"
	 *		0 <= z1 + z2 + z3 + z4 <= 3
	**/
	if( probdata->ndep ){
		/* constraint ld */

		SCIP_CONS**	ld;
		SCIP_CALL( SCIPallocMemoryArray( scip, &ld, probdata->ndep));

		ct=0;

		for(i=0; i<(probdata->ndep); ++i){
			(void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "ld_%d", i+1);
			SCIP_CALL( SCIPcreateConsBasicLinear( scip, &ld[i], consname,
							0, NULL, NULL, 0.0, ng[i]-1));
			for(j=0; j<p1; ++j){
				if( probdata->groupX[ct]==1 ){
					SCIP_CALL( SCIPaddCoefLinear( scip, ld[i], probdata->z[j], one));
				}
				ct++;
			}
			SCIP_CALL( SCIPaddCons(scip, ld[i]));
			SCIP_CALL( SCIPreleaseCons(scip, &ld[i]));
		}
		assert( ct == ( probdata->ndep * p1 ) );
		SCIPfreeMemoryArrayNull(scip, &ld);
	}

	SCIPfreeMemoryArrayNull( scip, &defbx);
	SCIPfreeMemoryArrayNull( scip, &defEXP);
	SCIPfreeMemoryArrayNull( scip, &defLOG);
	SCIPfreeMemoryArrayNull( scip, &bM1);
	SCIPfreeMemoryArrayNull( scip, &bM2);

   return SCIP_OKAY;
}

/** create initial columns */
static
SCIP_RETCODE createVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata           /**< problem data */
   )
{
	int	i;
	int	n = probdata->n;
	int	p = probdata->p;
	int	p1 = p+1;
	char	varname[SCIP_MAXSTRLEN];

	assert(scip != NULL);
	assert(probdata != NULL);
	assert( n >= 0 );
	assert( p >= 0 );

	/* create variables */
	SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->b, p1));
	SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->z, p1));
	SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->bx, n));
	SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->EXP, n));
	SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->LOG, n));


	for(i=0; i<p1; ++i){
		// b: parameters
		(void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "b%d", i);
		SCIP_CALL( SCIPcreateVarBasic( scip, &probdata->b[i], varname,
						-SCIPinfinity(scip), SCIPinfinity(scip),
						-2 * probdata->coef_obj[i], SCIP_VARTYPE_CONTINUOUS));
		// z: binary variables
		(void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "z%d", i);
		SCIP_CALL( SCIPcreateVarBasic( scip, &probdata->z[i], varname,
						0.0, 1.0, probdata->penalcf, SCIP_VARTYPE_BINARY));
	}

	for(i=0; i<n; ++i){
		// bx = b^t x
		(void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "bx(%d)", i+1);
		SCIP_CALL( SCIPcreateVarBasic( scip, &probdata->bx[i], varname,
						-SCIPinfinity(scip), SCIPinfinity(scip),
						0.0, SCIP_VARTYPE_CONTINUOUS));
		// EXP = exp( bx )
		(void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "EXP(%d)", i+1);
		SCIP_CALL( SCIPcreateVarBasic( scip, &probdata->EXP[i], varname,
						0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS));
		// LOG = log(1+EXP[i])
		(void) SCIPsnprintf( varname, SCIP_MAXSTRLEN, "LOG(%d)", i+1);
		SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->LOG[i], varname,
						0.0, SCIPinfinity(scip), 2.0, SCIP_VARTYPE_CONTINUOUS));

	}

	/* add variables to problem */
	for(i=0; i<p1; ++i){
		SCIP_CALL( SCIPaddVar(scip, probdata->b[i]));
		SCIP_CALL( SCIPaddVar(scip, probdata->z[i]));
	}

	for(i=0; i<n; i++){
		SCIP_CALL( SCIPaddVar( scip, probdata->bx[i]));
		SCIP_CALL( SCIPaddVar( scip, probdata->EXP[i]));
		SCIP_CALL( SCIPaddVar( scip, probdata->LOG[i]));
	}

	probdata->nvars =  p1 + p1 + n + n + n;

	/* branching variables */
	for(i=0; i<p1; ++i){
		SCIP_CALL( SCIPchgVarBranchPriority(scip, probdata->z[i], 1000));
	}

   return SCIP_OKAY;
}


/**@} */

/**@name Callback methods of problem data
 *
 * @{
 */

/** copies user data of source SCIP for the target SCIP */
static
SCIP_DECL_PROBCOPY(probcopyLogreg)
{

	int	n		=	sourcedata->n;
	int	p1		=	sourcedata->p + 1;
	int	ndep	=	sourcedata->ndep;
	int	i;
	SCIP_Bool success;

	SCIPdebugMessage("########################## probcopy ###########################\n");

	SCIP_CALL( probdataCreate(scip, targetdata) );

	(*targetdata)->n		=	sourcedata->n;
	(*targetdata)->p		=	sourcedata->p;
	(*targetdata)->ndep	=	sourcedata->ndep;
	(*targetdata)->nvars =	sourcedata->nvars;
	(*targetdata)->penalcf =	sourcedata->penalcf;

	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->y, n));
	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->x, n*p1));
	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->coef_obj, p1));

	if( ndep ){
		SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->Mdep, ndep));
		SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->groupX, ndep*p1));
	}

	for(i=0; i<n; i++)
		(*targetdata)->y[i] = sourcedata->y[i];

	for(i=0; i<n*p1; i++)
		(*targetdata)->x[i] = sourcedata->x[i];

	for(i=0; i<p1; i++)
		(*targetdata)->coef_obj[i] = sourcedata->coef_obj[i];

	if( ndep ){
		for(i=0; i<ndep; i++)
			(*targetdata)->Mdep[i] = sourcedata->Mdep[i];

		for(i=0; i<(ndep*p1); i++)
			(*targetdata)->groupX[i] = sourcedata->groupX[i];
	}else{
		(*targetdata)->Mdep = NULL;
		(*targetdata)->groupX = NULL;
	}

	/* variables */

	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->b, p1));
	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->z, p1));
	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->bx, n));
	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->EXP, n));
	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->LOG, n));

	for(i=0; i<p1; i++){
		SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->b[i],
			&((*targetdata)->b[i]), varmap, consmap, global, &success) );
		assert(success);
		SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->b[i]) );

		SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->z[i],
			&((*targetdata)->z[i]), varmap, consmap, global, &success) );
		assert(success);
		SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->z[i]) );
	}

	for(i=0; i<n; i++){
		SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->bx[i],
			&((*targetdata)->bx[i]), varmap, consmap, global, &success) );
		assert(success);
		SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->bx[i]) );

		SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->EXP[i],
			&((*targetdata)->EXP[i]), varmap, consmap, global, &success) );
		assert(success);
		SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->EXP[i]) );

		SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->LOG[i],
			&((*targetdata)->LOG[i]), varmap, consmap, global, &success) );
		assert(success);
		SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->LOG[i]) );
	}

	/* branching variables */
	for(i=0; i<p1; ++i){
		SCIP_CALL( SCIPchgVarBranchPriority(scip, (*targetdata)->z[i], 1000));
	}

	if( MP_MAXPOOL > 0 ){

		(*targetdata)->pool = new SolPool( p1 );
		if( sourcedata->pool->MySolPool_poolsize() > 0 ){

			copy( sourcedata->pool->pool.begin(), sourcedata->pool->pool.end(), back_inserter((*targetdata)->pool->pool));
			(*targetdata)->pool->MySolPool_updatesize();
		}
		(*targetdata)->pool->MFB_pool = sourcedata->pool->MFB_pool;

	}else{

		(*targetdata)->pool = NULL;

	}

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}


/** frees user data of original problem (called when the original problem is freed) */
static
SCIP_DECL_PROBDELORIG(probdelorigLogreg)
{
   SCIPdebugMessage("probdelorigLogreg \n");

   SCIPdebugMessage("free original problem data\n");

   /* free the (original) probdata */
   SCIP_CALL( probdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed) */
static
SCIP_DECL_PROBTRANS(probtransLogreg)
{
	SCIP_Real timelimit;
	SCIP_Bool update;

	int	n		=	sourcedata->n;
	int	p1		=	sourcedata->p + 1;
	int	ndep	=	sourcedata->ndep;
	int	i;

	SCIPdebugMessage("probtransLogreg \n");
   assert(scip != NULL);
   assert(sourcedata != NULL);
   assert(targetdata != NULL);
   assert(n > 0);
   assert(p1 > 1);
   assert(ndep >= 0);

	SCIP_CALL( SCIPgetBoolParam(scip, "logreg/countpresoltime", &update) );

	/* adjust time limit to take into account reading time */
	if( update )
	{
	   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
	   timelimit -= SCIPgetReadingTime(scip);
	   timelimit = MAX(0.0,timelimit);
	   SCIP_CALL( SCIPsetRealParam(scip, "limits/time", timelimit) );
	}

	/* create transform probdata */
	SCIP_CALL( probdataCreate(scip, targetdata) );

	(*targetdata)->n		=	sourcedata->n;
	(*targetdata)->p		=	sourcedata->p;
	(*targetdata)->ndep	=	sourcedata->ndep;
	(*targetdata)->nvars =	sourcedata->nvars;
	(*targetdata)->penalcf =	sourcedata->penalcf;

	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->y, n));
	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->x, n*p1));
	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->coef_obj, p1));

	if( ndep ){
		SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->Mdep, ndep));
		SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->groupX, ndep*p1));
	}

	for(i=0; i<n; i++)
		(*targetdata)->y[i] = sourcedata->y[i];

	for(i=0; i<n*p1; i++)
		(*targetdata)->x[i] = sourcedata->x[i];

	for(i=0; i<p1; i++)
		(*targetdata)->coef_obj[i] = sourcedata->coef_obj[i];

	if( ndep ){
		for(i=0; i<ndep; i++)
			(*targetdata)->Mdep[i] = sourcedata->Mdep[i];

		for(i=0; i<(ndep*p1); i++)
			(*targetdata)->groupX[i] = sourcedata->groupX[i];
	}else{
		(*targetdata)->Mdep = NULL;
		(*targetdata)->groupX = NULL;
	}

	/* variables */

	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->b, p1));
	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->z, p1));
	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->bx, n));
	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->EXP, n));
	SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->LOG, n));

	SCIP_CALL( SCIPtransformVars(scip, p1, sourcedata->b, (*targetdata)->b) );
	SCIP_CALL( SCIPtransformVars(scip, p1, sourcedata->z, (*targetdata)->z) );
	SCIP_CALL( SCIPtransformVars(scip, n, sourcedata->bx, (*targetdata)->bx) );
	SCIP_CALL( SCIPtransformVars(scip, n, sourcedata->EXP, (*targetdata)->EXP) );
	SCIP_CALL( SCIPtransformVars(scip, n, sourcedata->LOG, (*targetdata)->LOG) );

	/* branching variables */
	for(i=0; i<p1; ++i){
		SCIP_CALL( SCIPchgVarBranchPriority(scip, (*targetdata)->z[i], 1000));
	}

	if( MP_MAXPOOL > 0 ){

		(*targetdata)->pool = new SolPool( p1 );
		if( sourcedata->pool->MySolPool_poolsize() > 0 ){

			copy( sourcedata->pool->pool.begin(), sourcedata->pool->pool.end(), back_inserter((*targetdata)->pool->pool));
			(*targetdata)->pool->MySolPool_updatesize();
		}
		(*targetdata)->pool->MFB_pool = sourcedata->pool->MFB_pool;
	}else{

		(*targetdata)->pool = NULL;

	}

   return SCIP_OKAY;
}


static
SCIP_DECL_PROBEXITSOL(probexitsolLinreg)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(probdata != NULL);

   return SCIP_OKAY;
}

/** frees user data of transformed problem (called when the transformed problem is freed) */
static
SCIP_DECL_PROBDELTRANS(probdeltransLogreg)
{
   SCIPdebugMessage("free transformed problem data\n");

	SCIP_CALL( probdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** get problem name
 *
 *  Return NULL on error
 */
static
SCIP_RETCODE getProblemName(
	const char*	filename,	/*	input filename			  */
	char*			probname,	/*	output problemname	  */
	int			maxSize		/* maximum size of p.name */
	)
{
	int	i=0;
	int	j=0;
	int	l;

	/*	first find end of string */
	while( filename[i]!=0 )
		++i;
	l = i;

	/* go back until '.' or '/' or '\' appears */
	while( (i>0) && (filename[i]!='.') && (filename[i]!='/') && (filename[i]!='\\'))
		--i;

	/* if we found '.', search for '/' or '\\' */
	if( filename[i]=='.' ){
		l = i;
		while( (i>0) && (filename[i]!='/') && (filename[i]!='\\') )
			--i;
	}

	/* crrect counter */
	if( (filename[i]=='/') || (filename[i]=='\\') )
		++i;

	/* copy name */
	while( (i<l) && (filename[i]!=0) ){
		probname[j++] = filename[i++];
		if( j>maxSize-1)
			return SCIP_ERROR;
	}
	probname[j] = 0;

	return SCIP_OKAY;

}

/** exit a program if a input file has error */
static
SCIP_RETCODE readFail(
   int                   status
   )
{
   if( !status )
   {
      SCIPerrorMessage("Reading failed.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}


/** read dimension from the input file */
static
SCIP_RETCODE readDataDim(
   const char*           filename,           /**< name of the input file */
   int*                  n,                  /**< pointer to store the number of data points*/
   int*                  p,                  /**< pointer to store the  number of explanatory variables */
   int*                  i_ex                /**< pointer to store the index of explained varaible */
   )
{
   FILE *file;
   int status;
   char s[SCIP_MAXSTRLEN];

   assert(filename != NULL);
   assert(n != NULL);
   assert(p != NULL);
   assert(i_ex != NULL);

   /* open file */
   file = fopen(filename, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("Could not open file <%s>.\n", filename);
      return SCIP_ERROR;
   }

   /* skip one line */
   if( fgets(s, SCIP_MAXSTRLEN, file) == NULL )
   {
      SCIPerrorMessage("Error reading file <%s>.\n", filename);
      return SCIP_ERROR;
   }

   /* read n */
   status = fscanf(file, "%d", n);
   SCIP_CALL( readFail(status) );

   /* read p */
   status = fscanf(file, "%d", p);
   SCIP_CALL( readFail(status) );

   /* read i_ex */
   status = fscanf(file, "%d", i_ex);
   SCIP_CALL( readFail(status) );

   /* close file */
   fclose(file);

   return SCIP_OKAY;
}


/** read data points */
static
SCIP_RETCODE readData(
   const char*           filename,           /**< name of the input file */
   int                   n,                  /**< the number of data points */
   int                   p,                  /**< the number of explanatory varaibles */
   SCIP_Real*            data                /**< array to store data points */
   )
{
   int i;
   FILE *file;
   int buf;
   int status;
   char s[SCIP_MAXSTRLEN];

   /* open file */
   file = fopen(filename, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("Could not open file <%s>.\n", filename);
      return SCIP_ERROR;
   }

   /* skip one line */
   if( fgets(s, SCIP_MAXSTRLEN, file) == NULL )
   {
      SCIPerrorMessage("Error reading file <%s>.\n", filename);
      return SCIP_ERROR;
   }

   /* skip n */
   status = fscanf(file, "%d", &buf);
   SCIP_CALL( readFail(status) );

   /* skip p */
   status = fscanf(file, "%d", &buf);
   SCIP_CALL( readFail(status) );

   /* skip i_ex */
   status = fscanf(file, "%d", &buf);
   SCIP_CALL( readFail(status) );

   /* read data points */
   buf = n * (p + 1);
   for( i = 0; i < buf; i++ )
   {
      status = fscanf(file, "%lf", data + i);
      SCIP_CALL( readFail(status) );
   }

   /* close file */
   fclose(file);

   return SCIP_OKAY;
}

/** calculate the mean values for normalization */
static
SCIP_RETCODE calcMean(
   int                   n,                  /**< the number of data points */
   int                   p,                  /**< the number of explanatory varaibles */
   SCIP_Real*            data,               /**< data points */
   SCIP_Real*            mean                /**< array to store the mean values */
   )
{
   int i;
   int j;

   assert(n > 0);
   assert(p > 0);
   assert(data != NULL);
   assert(mean != NULL);

   for( i = 0; i < p + 1; i++ )
   {
      mean[i] = 0.0;
      for( j = 0; j < n; j++ )
      {
         mean[i] += data[j * (p + 1) + i];
      }
      mean[i] *= 1.0 / (SCIP_Real) n;
   }

   return SCIP_OKAY;
}


/** calculate the variance values for normalization */
static
SCIP_RETCODE calcVariance(
   int                   n,                  /**< the number of data points */
   int                   p,                  /**< the number of explanatory varaibles */
   SCIP_Real*            data,               /**< data points */
   SCIP_Real*            mean,               /**< the mean values */
   SCIP_Real*            variance            /**< array to store the variance values */
   )
{
   int i;
   int j;

   assert(n > 1);
   assert(p > 0);
   assert(data != NULL);
   assert(mean != NULL);
   assert(variance != NULL);

   for( i = 0; i < p + 1; i++ )
   {
      variance[i] = 0.0;
      for( j = 0; j < n; j++ )
      {
         variance[i] += pow(data[j * (p + 1) + i] - mean[i], 2.0);
      }
      variance[i] *= 1.0 / (SCIP_Real)(n - 1);
   }

   return SCIP_OKAY;
}


/** normalize data */
static
SCIP_RETCODE normalization(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   n,                  /**< the number of data points */
   int                   p,                  /**< the number of explanatory variables */
   int                   i_ex,
   SCIP_Real*            data                /**< data points */
   )
{
   int i;
   int j;
   SCIP_Real* mean;
   SCIP_Real* variance;

   assert(scip != NULL);
   assert(n > 0);
   assert(p > 0);
   assert(data != NULL);

   /* allocate memory for mean and variance */
   SCIP_CALL( SCIPallocMemoryArray(scip, &mean, p + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &variance, p + 1) );

   /* calculate the mean values and the variance values */
   SCIP_CALL( calcMean(n, p, data, mean) );
   SCIP_CALL( calcVariance(n, p, data, mean, variance) );

   for( i = 0; i < n; i++ )
   {
      for( j = 0; j < p + 1; j++ )
      {
         if( j != i_ex - 1 )
         {
            if( !SCIPisEQ(scip, variance[j], 0.0) )
               data[i * (p + 1) + j] = (data[i * (p + 1) + j] - mean[j]) / sqrt(variance[j]);
         }
      }
   }

   /* free */
   SCIPfreeMemoryArrayNull(scip, &mean);
   SCIPfreeMemoryArrayNull(scip, &variance);

   return SCIP_OKAY;
}

/** divide data into explained variable and explanatory variables */
static
SCIP_RETCODE divideData(
   int                   n,                  /**< the number of data points */
   int                   p,                  /**< the number of explanatory variables */
   int                   i_ex,               /**< index of the explained variable */
   SCIP_Real*            data,               /**< data points */
   SCIP_Real*            explained,          /**< array to store the explained variable */
   SCIP_Real*            explanatory         /**< array to store the explanatory variables */
   )
{
   int i;
   int j;
   int ct = 0;

   assert(n > 0);
   assert(p > 0);
   assert(i_ex > 0 && i_ex <= p+1);
   assert(data != NULL);
   assert(explained != NULL);
   assert(explanatory != NULL);

   for( i = 0; i < n; i++ )
   {
      explained[i] = data[i * (p + 1) + i_ex - 1];
      assert( explained[i] == 1 || explained[i] == 0 );

      if( !( explained[i] == 1 || explained[i] == 0 ) )
      {
         cout << "class[i] = " << explained[i] << endl;
         return SCIP_ERROR;
      }
   }

   for( i = 0; i < n; i++ )
   {
      explanatory[ct] = 1.0;
      ct++;
      for( j = 0; j < p + 1; j++ )
      {
         if( j != i_ex - 1 )
         {
            explanatory[ct] = data[i * (p + 1) + j];
            ct++;
         }
      }
   }

   assert( ct == n*(p+1) );

   return SCIP_OKAY;
}

/** sets up the problem data */
SCIP_RETCODE SCIPprobdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< file name */
   )
{
	SCIP_PROBDATA* probdata;
	int	n, p, i_ex, ndep;
	SCIP_Real*	data;	/* [n*(p+1)] */
	SCIP_Real*	X;		/* [n*(p+1)] */
	int*	D;				/* [p] */
	int*	ng;			/* ng: number of groupX[i] */


	int	i,j;
	int	ct=0;
	char	probname[SCIP_MAXSTRLEN];

	/* create problem data */
	SCIP_CALL( probdataCreate(scip, &probdata) );

	SCIP_CALL( getProblemName( filename, probname, SCIP_MAXSTRLEN));

   /* read dimension of data points */
   SCIP_CALL( readDataDim(filename, &n, &p, &i_ex) );

	SCIP_CALL( SCIPallocMemoryArray( scip, &data, n*(p+1)));

   /* read data points */
   SCIP_CALL( readData(filename, n, p, data) );

	probdata->n = n;
	probdata->p = p;

	assert( n > 0 );
	assert( p > 0 );
	assert( i_ex >= 0 );

	/* penalty coefficient */
	probdata->penalcf = 2.0;

   SCIP_CALL( checkData(n, p, data) );

   /* normalize data */
   SCIP_CALL( normalization(scip, n, p, i_ex, data) );

	SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->y, n));
	SCIP_CALL( SCIPallocMemoryArray(scip, &X, n*(p+1)));

   /* divide data into explained variable and explanatory variables */
   SCIP_CALL( divideData(n, p, i_ex, data, probdata->y, X) );

	SCIPfreeMemoryArrayNull( scip, &data);


	/* output information of problem */
	SCIPinfoMessage( scip, NULL, "File name\t:\t%s\n", filename);
	SCIPinfoMessage( scip, NULL, "Problem name \t:\t%s\n", probname);
	SCIPinfoMessage( scip, NULL, "Number of data\t:\t%d\n", n);
	SCIPinfoMessage( scip, NULL, "Number of var\t:\t%d\n", p);

	/**
	 * for using lapack and blas,
	 * define
	 *		x[n*(p+1)] := X[n][p+1]
	**/

	SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->x, n*(p+1)));
   SCIP_CALL( SCIPtransMat(n, p+1, X, probdata->x) );

	SCIPfreeMemoryArrayNull( scip, &X);

	/*	linear dependent */
	SCIP_CALL( SCIPallocMemoryArray(scip, &D, p+1));

   /* calculate the number of linearly dependent sets */
   SCIP_CALL( SCIPgetNLineDependSet(scip, probdata->x, n, p+1, D) );
   ndep = SCIPcalcIntSum(D, p+1);
   probdata->ndep = ndep;

	assert( ndep>=0 );

	if( ndep ){
		SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->Mdep, ndep));
		ct = 0;

		for(i=0; i<(p+1); ++i){
			if( D[i]==1 ){
				probdata->Mdep[ct++] = i;
			}
		}

		SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->groupX, ndep*(p+1)));
      SCIP_CALL( SCIPinitIntArrayZero(ndep*(p+1), probdata->groupX) );

		for(i=0; i<ndep; ++i){
			*(probdata->groupX+(i*(p+1))+probdata->Mdep[i]) = 1;
		}

      SCIP_Real *buf;
      SCIP_CALL( SCIPallocMemoryArray(scip, &buf, (p+1)*(p+1)));
      SCIP_CALL( SCIPcblasDgemm1(probdata->x, n, p+1, buf) );

      SCIP_CALL( SCIPgetLineDependSet(scip, buf, p+1, ndep, probdata->Mdep, D, probdata->groupX) );

      SCIPfreeMemoryArrayNull( scip, &buf);

		SCIP_CALL( SCIPallocBufferArray(scip, &ng, ndep));
		for(i=0; i<ndep; ++i){
			ng[i] = SCIPcalcIntSum( &(probdata->groupX[i*(p+1)]), p+1);
		}

	}else{
		probdata->Mdep		=	NULL;
		probdata->groupX	=	NULL;
		ng						=	NULL;
	}

	SCIPfreeMemoryArrayNull(scip, &D);

   /* print linearly dependent sets */
   SCIP_CALL( SCIPprintLineDependSet(scip, ndep, p+1, probdata->groupX) );

	// coef_obj: coefficients of obj.f
	SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->coef_obj, p+1));
   SCIP_CALL( SCIPinitArrayZero(p+1, probdata->coef_obj) );

	for(i=0; i<(p+1); i++){
		for(j=0; j<n; j++){
			probdata->coef_obj[i]	+=	probdata->y[j] * SCIPmatColMajor( probdata->x, n, j, i);
		}
	}

	/* create a problem in SCIP and add non-NULL callbacks via setter functions */
	SCIP_CALL( SCIPcreateProbBasic(scip, probname) );
	SCIP_CALL( SCIPsetProbDelorig(scip, probdelorigLogreg) );
	SCIP_CALL( SCIPsetProbTrans(scip, probtransLogreg) );
	SCIP_CALL( SCIPsetProbDeltrans(scip, probdeltransLogreg) );
	SCIP_CALL( SCIPsetProbExitsol(scip, probexitsolLinreg) );
	SCIP_CALL( SCIPsetProbCopy(scip, probcopyLogreg) );

	/* set objective sense */
	SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );

	/* set user problem data */
	SCIP_CALL( SCIPsetProbData(scip, probdata) );

	/* create and add variables */
	SCIP_CALL( createVariables(scip, probdata) );

	/* create and add constraints */
	SCIP_CALL( createConstraints(scip, probdata, ng) );


	/* free */
	if( ndep ){
		SCIPfreeBufferArray(scip, &ng);
	}

	// my solpool
	if( MP_MAXPOOL>0 ){

		probdata->pool = new SolPool( p+1 );

	}else{

		probdata->pool = NULL;

	}

   cout << endl;
   SCIP_CALL( forward(scip, probdata, probname) );
   cout << endl;
   SCIP_CALL( backward(scip, probdata, probname) );
   exit(0);

   return SCIP_OKAY;
}

/** returns the number of datas */
int SCIPprobdataGetNdatas(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->n;
}

/** returns the number of explain vars */
int SCIPprobdataGetNexvars(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->p;
}

/** returns the number of vars */
int SCIPprobdataGetNvars(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->nvars;
}

/** returns the number of linedep groups */
int SCIPprobdataGetNdep(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->ndep;
}

/** returns y */
extern
SCIP_Real*	SCIPprobdataGety(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->y;
}

/** returns X */
extern
SCIP_Real*	SCIPprobdataGetX(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->x;
}

/** returns coef_obj */
extern
SCIP_Real*	SCIPprobdataGetCO(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->coef_obj;
}

/** returns penalcf */
extern
SCIP_Real	SCIPprobdataGetPC(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->penalcf;
}

/** returns Mdep */
extern
int*	SCIPprobdataGetMdep(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->Mdep;
}

/** returns groupX */
extern
int*	SCIPprobdataGetgroupX(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->groupX;
}

/** returns var_b */
extern
SCIP_VAR**	SCIPprobdataGetVars_b(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->b;
}

/** returns var_z */
extern
SCIP_VAR**	SCIPprobdataGetVars_z(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->z;
}

/** returns var_bx */
extern
SCIP_VAR**	SCIPprobdataGetVars_bx(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->bx;
}

/** returns var_EXP */
extern
SCIP_VAR**	SCIPprobdataGetVars_EXP(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->EXP;
}

/** returns var_LOG */
extern
SCIP_VAR**	SCIPprobdataGetVars_LOG(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);

	return	probdata->LOG;
}

/* returns pool */
extern
SolPool*	SCIPprobdataGetPool(
   SCIP_PROBDATA*			probdata
	)
{
	assert(probdata != NULL);
	assert( probdata->pool != NULL );

	return probdata->pool;

}


/** writes end of log file */
/*
SCIP_RETCODE SCIPprobdataWriteLogfileEnd(
   SCIP*                 scip
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   if( probdata->logfile != NULL )
   {
      int success;
      SCIP_Real factor = 1.0;

      if( probdata->stp_type ==  STP_MAX_NODE_WEIGHT )
         factor = -1.0;

      SCIPprobdataWriteLogLine(scip, "End\n");
      SCIPprobdataWriteLogLine(scip, "\n");
      SCIPprobdataWriteLogLine(scip, "SECTION Run\n");
      if( probdata->ug )
      {
         SCIPprobdataWriteLogLine(scip, "Threads %d\n", probdata->nSolvers);
         SCIPprobdataWriteLogLine(scip, "Time %.1f\n", SCIPgetTotalTime(scip));
         SCIPprobdataWriteLogLine(scip, "Dual %16.9f\n", factor * probdata->ugDual);
      }
      else
      {
         SCIPprobdataWriteLogLine(scip, "Threads 1\n");
         SCIPprobdataWriteLogLine(scip, "Time %.1f\n", SCIPgetTotalTime(scip));
         SCIPprobdataWriteLogLine(scip, "Dual %16.9f\n", factor * SCIPgetDualbound(scip));
      }
      SCIPprobdataWriteLogLine(scip, "Primal %16.9f\n", factor * SCIPgetPrimalbound(scip));
      SCIPprobdataWriteLogLine(scip, "End\n");

      if( SCIPgetNSols(scip) > 0 )
      {
         SCIPprobdataWriteLogLine(scip, "\n");
         SCIPprobdataWriteLogLine(scip, "SECTION Finalsolution\n");

         SCIP_CALL( SCIPprobdataWriteSolution(scip, probdata->logfile) );
         SCIPprobdataWriteLogLine(scip, "End\n");
      }

      success = fclose(probdata->logfile);
      if( success != 0 )
      {
         SCIPerrorMessage("An error occurred while closing file <%s>\n", probdata->logfile);
         return SCIP_FILECREATEERROR;
      }

      probdata->logfile = NULL;
   }


   return SCIP_OKAY;
}
*/

