
/**@file   heur_forward.c
 * @brief  Forward selection
 *
 *
 *
 */

#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "heur_forward.h"
#include "probdata_logreg.h"
#include "mysolpool.h"
#include "logistic_regression.h"
#include "stepsize.h"
#include "call_cblas.h"
#include "convenient_tool.h"

#define HEUR_NAME             "forward"
#define HEUR_DESC             "primal heuristic using the forward selection"
#define HEUR_DISPCHAR         'f'
#define HEUR_PRIORITY         1000
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         10
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE

#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

#define debug                 0
#define MP_WARM               true

using namespace std;


/*
 * Local methods
 */

/** feasiblity by using linear dependence */
static
SCIP_Bool checkFeasiblity(
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

   return TRUE;
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyForward)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurForward(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeForward)
{   /*lint --e{715}*/

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecForward)
{  /*lint --e{715}*/

   // for probdata
   SCIP_PROBDATA* probdata;
   int   n;
   int   p;
   int   p1;
   int   ndep;
   int*  Mdep;             /* [ndep] */
   int*  groupX;           /* [ndep*p] */
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
   int   ublb;
   int   *branchinfo;      /* [3*p1] */

   int   sum_branchinfo[3];

   int i,j,ct;

   /* for forward selection */
   int   dim;
   int*  list;             /* [p1] */

   int   t;
   int   memo;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(result != NULL);

#if debug
   printf("forward selection!");
#endif

   /* get probdata */
   probdata = SCIPgetProbData(scip);
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
      ublb              =  SCIPround(scip, SCIPcomputeVarUbLocal(scip, var_z[i])
                        +  SCIPcomputeVarLbLocal(scip, var_z[i]));
      branchinfo[ublb*p1 + i]  =  1;
   }
   // }} get branching info

#if debug
   cout << "--------------------------------------" << endl;
   for(i=0; i<3; i++){
      for(j=0; j<p1; j++){
         printf("%d, ", *(branchinfo+(i*p1)+j));
      }
      cout << endl;
   }
#endif

   if( ndep )
   {
      /* check feasiblity by using linear dependence */
      if( checkFeasiblity(scip, p1, ndep, probdata, branchinfo) == FALSE )
      {
         /* cut-off */
         SCIPfreeBufferArray(scip, &branchinfo);
         *result = SCIP_DIDNOTFIND;
         return SCIP_OKAY;
      }
   }

   for( i = 0; i < 3; i++ )
      sum_branchinfo[i] = SCIPcalcIntSum( &branchinfo[p1*i], p1);

   if( sum_branchinfo[1] + sum_branchinfo[2] == 0 )
   {
      // if all z = 0
      SCIPfreeBufferArray(scip, &branchinfo);
      //*result = SCIP_DELAYED;
      //*result = SCIP_DIDNOTRUN;
      *result = SCIP_DIDNOTFIND;
      return SCIP_OKAY;
   }


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

#if debug
   cout << "list:";
   SCIPprintIntVec( p1, list);
#endif

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

   // forward selection
   while( 1 )
   {
      dim++;
      memo = -1;
      mL = SCIPinfinity(scip);

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

   SCIP_Bool store;
   int nsols = SCIPgetNSols(scip);
   int maxsol;

   SCIP_CALL( SCIPgetIntParam(scip, "limits/maxorigsol", &maxsol) );

   assert( nsols >= 0 );
   assert( maxsol > 0 );

   if( nsols < maxsol )
   {
      store = TRUE;
   }
   else
   {
      SCIP_Real bestval;
      bestval = SCIPgetSolOrigObj(scip, SCIPgetBestSol(scip));

      if( AIC <= bestval )
         store = TRUE;
      else
         store = FALSE;
   }

   if( store )
   {

      SCIP_SOL*   sol;
      int         nvars;
      SCIP_Real*  solvals;
      SCIP_VAR**  vars;
      SCIP_Bool   success;

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

      SCIP_CALL( SCIPcreateSol(scip, &sol, heur));
      SCIP_CALL( SCIPsetSolVals(scip, sol, nvars, vars, solvals));
      SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, FALSE, TRUE, TRUE, FALSE, &success) );

      SCIPfreeBufferArray(scip, &solvals);
      SCIPfreeBufferArray(scip, &vars);

   }

// // free
   SCIPfreeBufferArray(scip, &branchinfo);
   SCIPfreeBufferArray(scip, &list);
   SCIPfreeBufferArray(scip, &solval_01);
   SCIPfreeBufferArray(scip, &key);
   SCIPfreeBufferArray(scip, &index_nonzero);
   SCIPfreeBufferArray(scip, &pi);
   SCIPfreeBufferArray(scip, &Xb);
   SCIPfreeBufferArray(scip, &solval);
   SCIPfreeBufferArray(scip, &b_old);


   if( store == 1 ){
      *result = SCIP_FOUNDSOL;
   }else{
      //*result = SCIP_DELAYED;
      //*result = SCIP_DIDNOTRUN;
      *result = SCIP_DIDNOTFIND;
   }

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */


/** creates the local primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurForward(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEUR* heur;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecForward, NULL) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyForward) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeForward) );

   return SCIP_OKAY;
}
