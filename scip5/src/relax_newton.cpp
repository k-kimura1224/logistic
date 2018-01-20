
/**@file   relax_dposv.c
 * @brief
 */

#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "relax_newton.h"
#include "probdata_logreg.h"
#include "call_cblas.h"
#include "convenient_tool.h"
#include "get_lineardependence.h"
#include "logistic_regression.h"
#include "stepsize.h"

#include "mysolpool.h"

#define RELAX_NAME             "myrelaxator_newton"
#define RELAX_DESC             "relaxator newton"
#define RELAX_PRIORITY         1000
#define RELAX_FREQ             1
#define RELAX_INCLUDESLP       TRUE

#define debug                 0
#define MP_WARM               true


using namespace std;

/*
 * Local methods
 */

/** return the last branching variable */
static
SCIP_RETCODE getLastBranchingVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< the current node */
   SCIP_VAR**            result              /**< pointer to store the last branching variable */
   )
{
   SCIP_VAR** branchvars;                    /* array of variables on which the branchings has been performed in all ancestors */
   SCIP_Real* branchbounds;                  /* array of bounds which the branchings in all ancestors set */
   SCIP_BOUNDTYPE* boundtypes;               /* array of boundtypes which the branchings in all ancestors set */
   int* nodeswitches;                        /* marks, where in the arrays the branching decisions of the next node on the path start
                                              * branchings performed at the parent of node always start at position 0. For single variable branching,
                                              * nodeswitches[i] = i holds */
   int nbranchvars;                          /* number of variables on which branchings have been performed in all ancestors
                                              *   if this is larger than the array size, arrays should be reallocated and method should be called again */
   int branchvarssize;                       /* available slots in arrays */
   int nnodes;                               /* number of nodes in the nodeswitch array */
   int nodeswitchsize;                       /* available slots in node switch array */

   SCIP_VAR* lastbranchvar;

   branchvarssize = SCIPnodeGetDepth(node);
   nodeswitchsize = branchvarssize;

   /* memory allocation */
   SCIP_CALL( SCIPallocBufferArray(scip, &branchvars, branchvarssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &branchbounds, branchvarssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &boundtypes, branchvarssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodeswitches, nodeswitchsize) );

   SCIPnodeGetAncestorBranchingPath(node, branchvars, branchbounds, boundtypes, &nbranchvars, branchvarssize, nodeswitches, &nnodes, nodeswitchsize);

   /* if the arrays were to small, we have to reallocate them and recall SCIPnodeGetAncestorBranchingPath */
   if( nbranchvars > branchvarssize || nnodes > nodeswitchsize )
   {
      branchvarssize = nbranchvars;
      nodeswitchsize = nnodes;

      /* memory reallocation */
      SCIP_CALL( SCIPreallocBufferArray(scip, &branchvars, branchvarssize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &branchbounds, branchvarssize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &boundtypes, branchvarssize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &nodeswitches, nodeswitchsize) );

      SCIPnodeGetAncestorBranchingPath(node, branchvars, branchbounds, boundtypes, &nbranchvars, branchvarssize, nodeswitches, &nnodes, nodeswitchsize);
      assert(nbranchvars == branchvarssize);
   }

   /* get the last branching */
   if( nbranchvars >= 1 )
      lastbranchvar = branchvars[nodeswitches[0]];
   else
      lastbranchvar = NULL;

   assert(lastbranchvar == NULL || SCIPvarGetType(lastbranchvar) == SCIP_VARTYPE_BINARY);

   /* free all local memory */
   SCIPfreeBufferArray(scip, &nodeswitches);
   SCIPfreeBufferArray(scip, &boundtypes);
   SCIPfreeBufferArray(scip, &branchbounds);
   SCIPfreeBufferArray(scip, &branchvars);

   *result = lastbranchvar;
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


/** calculate sum_branchinfo */
static
SCIP_RETCODE calcSumBranchinfo(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   p,                  /**< the number of explantory variables */
   int                   p2,                 /**< p * p */
   int*                  branchinfo,         /**< branching information */
   int*                  sum_branchinfo
   )
{
   int* ip0 = &branchinfo[0];
   int* ip1 = &branchinfo[p2];
   int i;

   assert(scip != NULL);
   assert(p > 0);
   assert(p2 == 2 * p);
   assert(branchinfo != NULL);
   assert(sum_branchinfo != NULL);

   sum_branchinfo[0] = 0;
   sum_branchinfo[2] = 0;

   for( i = 0; i < p; i++ )
      sum_branchinfo[0] += ip0[i];

   for( i = 0; i < p; i++ )
      sum_branchinfo[2] += ip1[i];

   sum_branchinfo[1] = p - sum_branchinfo[0] - sum_branchinfo[2];

   return SCIP_OKAY;
}

/*
 * Callback methods of relaxator
 */

/** copy method for relaxator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_RELAXCOPY(relaxCopyNewton)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(relax != NULL);
   assert(strcmp(SCIPrelaxGetName(relax), RELAX_NAME) == 0);

   /* call inclusion method of relaxator */
   SCIP_CALL( SCIPincludeRelaxNewton(scip));

   return SCIP_OKAY;
}

/** destructor of relaxator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_RELAXFREE(relaxFreeNewton)
{   /*lint --e{715}*/

   return SCIP_OKAY;
}

/** execution method of relaxator */
static
SCIP_DECL_RELAXEXEC(relaxExecNewton)
{

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
   int   *branchinfo;      /* [3*p1] */

   SCIP_NODE*  node;

   //int   root;

   int i,j,ct;

   assert(relax != NULL);
   assert(scip != NULL);
   assert(strcmp(SCIPrelaxGetName(relax), RELAX_NAME) == 0);
   assert(result != NULL);

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

   /* get current node */
   node     =  SCIPgetCurrentNode(scip);

   /* get the last branching variable */
   SCIP_VAR* lastbranchvar;

   SCIP_CALL( getLastBranchingVar(scip, node, &lastbranchvar) );
   assert(lastbranchvar == NULL || SCIPvarGetType(lastbranchvar) == SCIP_VARTYPE_BINARY);

   /* this means that lastbranchvar is fixed to 1 */
   if( lastbranchvar != NULL &&
         (int)SCIPround(scip, SCIPcomputeVarLbLocal(scip, lastbranchvar)) == 1 )
   {
      SCIP_NODE* parent;
      SCIP_Real parentdual;

      assert((int)SCIPround(scip, SCIPcomputeVarUbLocal(scip, lastbranchvar) ) == 1);

      /* get parent node */
      parent = SCIPnodeGetParent(node);

      /* get local lower bound */
      parentdual = SCIPgetNodeLowerbound(scip, parent);

      /* update */
      SCIP_CALL( SCIPupdateLocalLowerbound(scip, parentdual + penalcf) );

      *result = SCIP_CUTOFF;
      for( i = 0; i < p1; i++ )
      {
         if( SCIPround(scip, SCIPcomputeVarUbLocal(scip, var_z[i]) ) == 1
               && SCIPround(scip, SCIPcomputeVarLbLocal( scip, var_z[i]) ) == 0 )
         {
            *result = SCIP_SUCCESS;
            break;
         }
      }

      return SCIP_OKAY;
   }

   // get branching info {{
   // alloc and initialize
   SCIP_CALL( SCIPallocBufferArray(scip, &branchinfo, 3 * p1) );

   /* branchinfo[i] = 1 if z_i is fixed to 0 */
   int* ip0 = &branchinfo[0];
   for( i = 0; i < p1; i++ )
   {
      ip0[i] = 1 - SCIPround(scip, SCIPcomputeVarUbLocal(scip, var_z[i]) );
      assert(ip0[i] == 1 || ip0[i] == 0);
   }

   /* branchinfo[2p+i] = 1 if z_i is fixed to 1 */
   int* ip1 = &branchinfo[p1*2];
   for( i = 0; i < p1; i++ )
   {
      ip1[i] = SCIPround(scip, SCIPcomputeVarLbLocal(scip, var_z[i]) );
      assert(ip1[i] == 1 || ip1[i] == 0);
   }

   /* branchinfo[p+i] = 1 if z_i is not fixed */
   int* ip = &branchinfo[p1];
   for( i = 0; i < p1; i++ )
   {
      ip[i] = 1 - ip0[i] - ip1[i];
      assert(ip[i] == 1 || ip[i] == 0);
   }
   // }} get branching info

#if debug
   //cout << "(" << node->number << ")";
   cout << "---------------------------------------------";
   cout << "---------------------------------------------" << endl;
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
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
   }

   /* calculate sum_branchinfo */
   int   sum_branchinfo[3];
   SCIP_CALL( calcSumBranchinfo(scip, p1, p1 * 2, branchinfo, sum_branchinfo) );

   int dimb = sum_branchinfo[1] + sum_branchinfo[2];
   int dimz = sum_branchinfo[1];
   int dimz_ = dimz;

   if( dimb == 0 )
   {
      // if all z = 0
      SCIPfreeBufferArray(scip, &branchinfo);
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   int *index_nonzero;
   int *solval_01;
   SCIP_CALL( SCIPallocBufferArray(scip, &index_nonzero, dimb));
   SCIP_CALL( SCIPallocBufferArray(scip, &solval_01, p1));

   ct = 0;
   for( i = 0; i < p1; i++ )
   {
      if( branchinfo[i] == 0 )
      {
         index_nonzero[ct] = i;
         ct++;

         solval_01[i] = 1;
      }
      else
      {
         solval_01[i] = 0;
      }
   }

   assert( ct==dimb );

   // for the newton method
   SCIP_Real   *subX_;        // [n*dimb],   colMajor
   SCIP_Real   *subcoef;      // [dimb],
   SCIP_Real   *A_;           // [dimb*dimb],colMajor
   SCIP_Real   *Y_;           // [n*dim],    colMajor
   SCIP_Real   *q;            // [dimb],
   SCIP_Real   *pi;           // [n],
   SCIP_Real   *point;        // [dimb], current point
   SCIP_Real   *d;            // [dimb], direction vector
   SCIP_Real   *Xb;           // [n]

   SCIP_Real   objval;        // object value
   SCIP_Real   objval_new;    // new object value
   SCIP_Real   ss;            // step size
   SCIP_Real   newton_gap;
   SCIP_Real norm;

   int info;

   SCIP_Real   lb;
   Solution*      mysol;

   // alloc
   SCIP_CALL( SCIPallocBufferArray(scip, &point, dimb));

   // newton method
   // step0: initialize {{
   //    generate initial point {
   if( MP_WARM && lastbranchvar != NULL )
   {
      int *key;
      int n_key = pool->MySolPool_Nkey();
      int index;

      // 1. search the relaxsol of the current node
      SCIP_CALL( SCIPallocBufferArray(scip, &key, n_key));

      pool->MySolPool_genekey( solval_01, key);
      index = pool->MySolPool_check( key );

      if( index > -1 )
      {
#if debug
   cout << "the relaxsol is found from the pool" << endl;
#endif
         mysol = pool->MySolPool_getsol( index );
         lb = 2 * mysol->val_mL + penalcf * (double)sum_branchinfo[2];
         SCIP_CALL( SCIPupdateLocalLowerbound(scip, lb));

         SCIPfreeBufferArray(scip, &branchinfo);
         SCIPfreeBufferArray(scip, &index_nonzero);
         SCIPfreeBufferArray(scip, &solval_01);
         SCIPfreeBufferArray(scip, &point);
         SCIPfreeBufferArray(scip, &key);

         if( dimz == 0 || dimz_ == 0 )
         {
            *result = SCIP_CUTOFF;
         }else{
            *result=SCIP_SUCCESS;
         }

         return SCIP_OKAY;
      }

      // 2. search the relaxsol of the parent node
      int*  val_01;
      SCIP_CALL( SCIPallocBufferArray(scip, &val_01, p1));
      for(i=0; i<p1; i++){
         val_01[i] = solval_01[i];
      }

      int parentbr = -1;
      for( i = 0; i < p1; i++ )
      {
         if( var_z[i] == lastbranchvar )
         {
            parentbr = i;
            break;
         }
      }

      assert( parentbr >= 0 && parentbr < p1 );

      val_01[parentbr] = 1;

      pool->MySolPool_genekey( val_01, key);
      index = pool->MySolPool_check( key );

      info = 0;
      if( index > -1 )
      {
#if debug
   cout << "generate the an initial point from the relaxsol of the parent node" << endl;
#endif
         mysol = pool->MySolPool_getsol( index );

         for(i=0; i<dimb; i++){
            point[i] = mysol->val[index_nonzero[i]];
         }
         info = 1;
      }
      else
      {
      // 3. search
         for(i=0; i<p1; i++){
            for(j=0; j<p1; j++){
               val_01[j] = solval_01[j];
            }
            if( val_01[i] == 1 ){
               val_01[i] = 0;
            }else{
               val_01[i] = 1;
            }

            pool->MySolPool_genekey( val_01, key);
            index = pool->MySolPool_check( key );

            if( index > -1){
#if debug
   cout << "generate the intial point" << endl;
#endif
               info = 1;
               mysol = pool->MySolPool_getsol( index );

               for(j=0; j<dimb; j++){
                  point[j] = mysol->val[index_nonzero[j]];
               }
               break;
            }
         }
      }

      // 4.
      if( info == 0 )
      {
         for(i=0; i<dimb; i++)
            point[i] = 0.0;

         info = 1;
      }

      // free
      SCIPfreeBufferArray( scip, &val_01);
      SCIPfreeBufferArray(scip, &key);

   }else{

      for(i=0; i<dimb; i++){
         point[i] = 0.0;
      }

      info = 1;
   }

   assert( info == 1 );

#if debug
   cout << "dimb: " << dimb << endl;
   cout << "initial point: ";
   SCIPprintVec( dimb, point );
#endif

   // alloc
   SCIP_CALL( SCIPallocBufferArray(scip, &d, dimb));
   SCIP_CALL( SCIPallocBufferArray(scip, &A_, dimb*dimb));
   SCIP_CALL( SCIPallocBufferArray(scip, &Y_, n*dimb));
   SCIP_CALL( SCIPallocBufferArray(scip, &q, dimb));
   SCIP_CALL( SCIPallocBufferArray(scip, &subX_, n*dimb));
   SCIP_CALL( SCIPallocBufferArray(scip, &subcoef, dimb));
   SCIP_CALL( SCIPallocBufferArray(scip, &pi, n));
   SCIP_CALL( SCIPallocBufferArray(scip, &Xb, n));

   // } generate initial point
   // generate sub matrix of X and sub vector of coef_obj {
   int buf_int;
   ct = 0;
   for( i = 0; i < dimb; i++ )
   {
      buf_int = index_nonzero[i] * n;
      for( j = 0; j < n; j++ )
         subX_[ct++] = X[buf_int + j];
   }

   assert( ct == ( n * dimb ));

   for( i = 0; i < dimb; i++ )
      subcoef[i] = coef_obj[index_nonzero[i]];

   // } generate sub matrix of X and sub vector of coef_obj
   // }} step0

   SCIP_CALL( SCIPcblasDgemv2( subX_, n, dimb, point, Xb) );
   objval = - Loglikelifood_( scip, n, dimb, subcoef, Xb, point);
   assert( objval < SCIPinfinity(scip) );
//
//   GenerateZeroVec( dimb*dimb, A_);
//
#if debug
   cout << "objval: " << objval << endl;
#endif

   int ct_newton=0;

   while(1)
   {
      ct_newton++;

#if debug
      cout << "newton method(" << ct_newton << ")";
      cout << "--------------------------------" << endl;
#endif
      /* step1: find a descent direction by solving Ad = q
       *    where d is the descent direction,

       *    A:= X^t P ( I_n - P ) X   in ( dim, dim)
       *       X:= subX_  in ( n, dim)
       *       P(b):= diag(pi_1(b),.., pi_n(b))
       *          b is the current point
       *       I_n is (dim-n) identity matrix
       *
       *    q:= c - X^t pi  in ( dim, 1)
       *       c is the subcoef  in ( dim, 1)
       *       pi:= (pi_1(b),.., pi_n(b))^t
       *          b is the current point
       */

      // define pi
      for( i = 0; i < n; i++ )
      {
         if( Xb[i] < 13.0 )
            pi[i] = 1.0 - ( 1.0 / ( 1.0 + exp(Xb[i]) ) );
         else
            pi[i] = 1.0;
      }

      // define q:= c - X^t pi
      SCIP_CALL( SCIPcblasDgemv4( n, dimb, subX_, pi, subcoef, -1.0, 1.0, q) );

      // define Y:= P ( I_n - P ) X  in ( n, dim)
      for( i = 0; i < n; i++ )
         pi[i] *= 1 - pi[i];

      ct=0;
      for( j = 0; j < dimb; j++ )
      {
         for( i = 0; i < n; i++ )
         {
            Y_[ct] = pi[i] * subX_[ct];
            ct++;
         }
      }
      assert( ct == n*dimb );

      // define A:= X^t Y
      SCIP_CALL( SCIPcblasDgemm2( n, dimb, dimb, subX_, Y_, A_) );

#if 0
   cout << "A_:" << endl;
   printM_( A_, dimb, dimb);
#endif


#if 0
   cout << "q: " << endl;
   printv( dimb, q);
#endif

      // solve Ad = q
      info = SCIPclapackDposv( scip, A_, q, dimb, d);

      if( info != 0 )
      {
//#if debug
#if 0
         cout << "warning: info= " << info << " , objval= " << objval << endl;
         //cout << "A =";
         //printM_( A_, dimb, dimb);

         int *Dsub;
         Dsub = new int[dimb];
         LinearDependent( n, dimb, subX_, Dsub);

         if( sumint( Dsub, dimb) != 0 ){
            cout << "error: l.d." << endl;
            exit(1);
         }

         //printv( Dsub, dimb);
         //LinearDependent_( dimb, dimb, A_, Dsub);
         //printv( Dsub, dimb);
         delete[] Dsub;
#endif
//       ct_dposer++;
//
//       if( ct_dposer<3 ){
//          SCIP_Real max_point = -1.0;
//          int       max_num = -1;
//          for(i=0; i<dimb; i++){
//             if( max_point < fabs(point[i]) ){
//                max_point   =  fabs(point[i]);
//                max_num     =  i;
//             }
//          }
//
//          assert( max_num != -1 );
//
//          point[max_num] = point[max_num]/2.0;
//          objval   =  -Loglikelifood_( n, dimb, subcoef, subX_, point);
//#if debug
// cout << "point: ";
// printv( dimb, point);
// cout << "objval: " << objval << endl;
//#endif
//
//       }else{
//          cout << "error:---" << endl;
//          stop();
//       }

//    ]else[
         //SCIP_CALL( SCIPcblasCopy( q, d, dimb) );
         for( i = 0; i < dimb; i++ )
            point[i] = 0.0;

         SCIP_CALL( SCIPcblasDgemv2( subX_, n, dimb, point, Xb) );
         objval = - Loglikelifood_( scip, n, dimb, subcoef, Xb, point);
         continue;
      }

      // step2: find the stepsize
      norm = SCIPcblasDnrm( d, dimb );
      SCIP_CALL( SCIPcblasDscal( d, dimb, 1.0/norm, d ) );
      ss = FindStepsize( scip, n, dimb, subcoef, subX_, point , d);

      assert( ss >= 0  );

      if( ss <= 1e-05 ){
         ss = 0.0;
      }

      // step3: update a new point and its objval
      SCIP_CALL( SCIPcblasDaxpy( d, point, dimb, ss, 1.0, point) );
      SCIP_CALL( SCIPcblasDgemv2( subX_, n, dimb, point, Xb) );
      objval_new = - Loglikelifood_( scip, n, dimb, subcoef, Xb, point);
#if debug
   cout << "d: ";
   SCIPprintVec( dimb, d);
   cout << "stepsize: " << ss << endl;
   cout << "new point: ";
   SCIPprintVec( dimb, point);
   cout << "new objval: " << objval_new << endl;
#endif

      assert( objval - objval_new >= - 1.0e-04 );

      newton_gap = pow( objval - objval_new, 2.0);
#if debug
   cout << "gap: " << newton_gap << endl;
#endif

      if( newton_gap < 1e-08 ){
         objval = objval_new;
         break;
      }

      objval = objval_new;
   }

   // update the local bound
   lb = 2 * objval + penalcf * (double)sum_branchinfo[2];
   SCIP_CALL( SCIPupdateLocalLowerbound(scip, lb));

   // store the relaxation solution to the pool
   ct = 0;

   SCIP_Real*     solval;     // [p1]
   SCIP_CALL( SCIPallocBufferArray(scip, &solval, p1));

   for( i = 0; i < p1; i++ )
   {
      if( solval_01[i] == 1 )
      {
         solval[i] = point[ct];
         ct++;
      }
      else
      {
         solval[i] = 0.0;
      }
   }

   if( MP_WARM == true && MP_MAXPOOL > 0 ){
      pool->MySolPool_store( solval_01, solval, objval);
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

      SCIP_Real ub;
      ub = lb + penalcf * (double)sum_branchinfo[1];

      if( ub <= bestval )
         store = TRUE;
      else
         store = FALSE;

   }

   if( store )
   {
      // set primal solution {{

      SCIP_SOL*   sol;
      SCIP_HEUR*  heur;
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
         solvals[ct] =  (double)solval_01[i];
         ct++;
      }

      // b
      for(i=0; i<p1; i++){
         vars[ct]    = var_b[i];
         solvals[ct] = solval[i];
         ct++;
      }

      // bx
      for(i=0; i<n; i++){
         vars[ct]    =  var_bx[i];
         solvals[ct] =  Xb[i];
         ct++;
      }

      // EXP
      for(i=0; i<n; i++){
         vars[ct]    =  var_EXP[i];
         solvals[ct] =  exp( Xb[i] );
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
      SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, FALSE, TRUE, TRUE, FALSE, &success) );

      SCIPfreeBufferArray(scip, &solvals);
      SCIPfreeBufferArray(scip, &vars);

      // }} set primal solution
   }

   // free
   SCIPfreeBufferArray(scip, &branchinfo);
   SCIPfreeBufferArray(scip, &index_nonzero);
   SCIPfreeBufferArray(scip, &solval_01);
   SCIPfreeBufferArray(scip, &point);
   SCIPfreeBufferArray(scip, &d);
   SCIPfreeBufferArray(scip, &A_);
   SCIPfreeBufferArray(scip, &Y_);
   SCIPfreeBufferArray(scip, &q);
   SCIPfreeBufferArray(scip, &subX_);
   SCIPfreeBufferArray(scip, &subcoef);
   SCIPfreeBufferArray(scip, &pi);
   SCIPfreeBufferArray(scip, &Xb);
   SCIPfreeBufferArray(scip, &solval);

   if( dimz == 0 || dimz_ == 0 )
   {
      *result = SCIP_CUTOFF;
   }else{
      *result=SCIP_SUCCESS;
   }

   return SCIP_OKAY;
}

/*
 * relaxator  specific interface methods
 */


/** creates the myrelaxator and includes it in SCIP */
SCIP_RETCODE SCIPincludeRelaxNewton(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAX* relax;

   /* include relaxator */
   SCIP_CALL( SCIPincludeRelaxBasic(scip, &relax, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ,
         relaxExecNewton, NULL) );

   assert(relax != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetRelaxCopy(scip, relax, relaxCopyNewton) );
   SCIP_CALL( SCIPsetRelaxFree(scip, relax, relaxFreeNewton) );

   return SCIP_OKAY;
}
