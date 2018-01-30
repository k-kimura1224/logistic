
/**@file   heur_backward.c
 * @brief  Backward selection
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

#include "heur_backward.h"
#include "probdata_logreg.h"
#include "mysolpool.h"
#include "logistic_regression.h"
#include "stepsize.h"
#include "call_cblas.h"
#include "convenient_tool.h"

#define HEUR_NAME             "backward"
#define HEUR_DESC             "primal heuristic using the backward selection"
#define HEUR_DISPCHAR         'b'
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

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyBackward)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurBackward(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeBackward)
{   /*lint --e{715}*/

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecBackward)
{  /*lint --e{715}*/

   // for probdata
   SCIP_PROBDATA* probdata;
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
   int   ublb;
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

   assert(heur != NULL);
   assert(scip != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(result != NULL);

#if debug
   printf("backward selection!\n");
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

   // get branching info {{
   // alloc and initialize
   SCIP_CALL( SCIPallocBufferArray(scip, &branchinfo, 3*p1));
   SCIP_CALL( SCIPinitIntArrayZero( 3*p1, branchinfo ) );

   for(i=0; i<p1; ++i){
      ublb              =  SCIPround(scip, SCIPcomputeVarUbLocal(scip, var_z[i])
                        +  SCIPcomputeVarLbLocal(scip, var_z[i]));
      branchinfo[ublb * p1 +i ]  =  1;
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
      /* fix non-fixed variable z and check feasiblity by using linear dependence */
      if( fixVariable(scip, p1, ndep, probdata, branchinfo) == FALSE )
      {
         /* cut-off */
         SCIPfreeBufferArray(scip, &branchinfo);
         *result = SCIP_DIDNOTFIND;
         return SCIP_OKAY;
      }
   }

   for(i=0; i<3; i++){
      sum_branchinfo[i] = SCIPcalcIntSum( &branchinfo[p1*i], p1);
   }

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
   int ct_error;

   SCIP_CALL( SCIPallocBufferArray(scip, &pi, n));
   SCIP_CALL( SCIPallocBufferArray(scip, &Xb, n));

   // backward selection
   // mark
   while(1)
   {
      dim--;
      memo = -1;
      mL = 1e+06;

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

#if debug
            SCIPprintVec( dim, point);
#endif

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

            ct_error = 0;
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
                     objval = SCIPinfinity(scip);
                     ct_error++;
                     break;
                  }
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

            if( MP_MAXPOOL > 1 && ct_error <= 1 ){
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

   SCIP_Bool store;
   int nsols = SCIPgetNSols(scip);
   int maxsol;

   SCIP_CALL( SCIPgetIntParam(scip, "limits/maxorigsol", &maxsol) );

   assert(maxsol > 0);
   assert( nsols >= 0 );

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

   if( store ){

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

      SCIP_CALL( SCIPcreateSol(scip, &sol, heur));
      SCIP_CALL( SCIPsetSolVals(scip, sol, nvars, vars, solvals));
      SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, FALSE, TRUE, TRUE, FALSE, &success) );

      SCIPfreeBufferArray(scip, &solvals);
      SCIPfreeBufferArray(scip, &vars);

   }

   // free
   SCIPfreeBufferArray(scip, &branchinfo);
   SCIPfreeBufferArray(scip, &solval_01);
   SCIPfreeBufferArray(scip, &solval);
   SCIPfreeBufferArray(scip, &val_01);
   SCIPfreeBufferArray(scip, &list);
   SCIPfreeBufferArray(scip, &b_old);
   SCIPfreeBufferArray(scip, &pi);
   SCIPfreeBufferArray(scip, &Xb);

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
SCIP_RETCODE SCIPincludeHeurBackward(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEUR* heur;


   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecBackward, NULL) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyBackward) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeBackward) );

   return SCIP_OKAY;
}
