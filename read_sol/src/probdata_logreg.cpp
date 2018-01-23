
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

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include "probdata_logreg.h"
#include "convenient_tool.h"
#include "get_lineardependence.h"
#include "call_cblas.h"
#include "logistic_regression.h"
#include "stepsize.h"
//#include "read_data.h"
//#include "normalization.h"
//#include "divide_data.h"
//#include "matrix.h"
//#include "vector.h"
//#include "cblapack.h"
//#include "mysolpool.h"
//#include "set_myparameter.h"

#define bigM				300.0
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

static
SCIP_RETCODE newton(
      SCIP*                 scip,
      SCIP_PROBDATA*        probdata,            /**< problem data */
      int*                  solvalint
      )
{
   int n = probdata->n;
   int p = probdata->p;
   int p1 = p + 1;
   SCIP_Real penalcf  = SCIPprobdataGetPC( probdata );
   SCIP_Real* X        = SCIPprobdataGetX( probdata );
   SCIP_Real* coef_obj = SCIPprobdataGetCO( probdata );

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
   int dimb = 0;
   int *index_nonzero;

   int i, j;
   int ct;

   int ndep = probdata->ndep;

   if( ndep > 0 )
   {
      int* Mdep = probdata->Mdep;
      int* groupX = probdata->groupX;
      int buf;
      int dummy;

      for( i = 0; i < ndep; ++i )
       {
          dummy = -1;
          buf = i * p1;
          for( j = 1; j < p1; ++j )
          {
             if( groupX[buf+j] == 1 )
             {
                if( solvalint[j] == 0 )
                   break;
                if( solvalint[j] == 1 )
                   dummy = j;
                if( j == Mdep[i] )
                {
                   if( dummy != -1 )
                   {
                      solvalint[dummy] = 0;
                      break;
                   }
                   else
                   {
                      cout << "error" << endl;
                      exit(1);
                   }
                }
             }
          }
       }
   }

   for( i = 0; i < p1; i++ )
      dimb += solvalint[i];

   cout << "k: " << dimb << endl;

   // alloc
   SCIP_CALL( SCIPallocBufferArray(scip, &index_nonzero, dimb));

   ct = 0;
   for( i = 0; i < p1; i++ )
   {
      if( solvalint[i] == 1 )
         index_nonzero[ct++] = i;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &point, dimb));

    for( i = 0; i < dimb; i++ )
      point[i] = 0.0;

   // alloc
   SCIP_CALL( SCIPallocBufferArray(scip, &d, dimb));
   SCIP_CALL( SCIPallocBufferArray(scip, &A_, dimb*dimb));
   SCIP_CALL( SCIPallocBufferArray(scip, &Y_, n*dimb));
   SCIP_CALL( SCIPallocBufferArray(scip, &q, dimb));
   SCIP_CALL( SCIPallocBufferArray(scip, &subX_, n*dimb));
   SCIP_CALL( SCIPallocBufferArray(scip, &subcoef, dimb));
   SCIP_CALL( SCIPallocBufferArray(scip, &pi, n));
   SCIP_CALL( SCIPallocBufferArray(scip, &Xb, n));

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

   SCIP_Real AIC = 2 * objval + penalcf * dimb;

   cout << "AIC: " << AIC << endl;


   SCIPfreeBufferArray(scip, &index_nonzero);
   SCIPfreeBufferArray(scip, &point);
   SCIPfreeBufferArray(scip, &d);
   SCIPfreeBufferArray(scip, &A_);
   SCIPfreeBufferArray(scip, &Y_);
   SCIPfreeBufferArray(scip, &q);
   SCIPfreeBufferArray(scip, &subX_);
   SCIPfreeBufferArray(scip, &subcoef);
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
		(void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "b_%d", i);
		SCIP_CALL( SCIPcreateVarBasic( scip, &probdata->b[i], varname,
						-SCIPinfinity(scip), SCIPinfinity(scip),
						-2 * probdata->coef_obj[i], SCIP_VARTYPE_CONTINUOUS));
		// z: binary variables
		(void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "z_%d", i);
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
      explained[i] = data[i * (p + 1) + i_ex - 1];

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

   {
      // read sol
      SCIP_SOL* readsol_05;
      SCIP_SOL* readsol_09;
      SCIP_SOL* readsol_17;
      SCIP_Bool error;
      SCIP_HEUR* heur;
      char varname[SCIP_MAXSTRLEN];
      char solutionfile[SCIP_MAXSTRLEN];
      SCIP_VAR* var;

      int* solvalint;

      heur = SCIPfindHeur(scip, "trysol");
      SCIP_CALL( SCIPcreateSol(scip, &readsol_05, heur));
      SCIP_CALL( SCIPcreateSol(scip, &readsol_09, heur));
      SCIP_CALL( SCIPcreateSol(scip, &readsol_17, heur));

      /* read solution */
      (void) SCIPsnprintf( solutionfile, SCIP_MAXSTRLEN, "%s_05.sol", probname);
      SCIP_CALL( SCIPreadSolFile(scip, solutionfile, readsol_05, TRUE, NULL, &error));

      (void) SCIPsnprintf( solutionfile, SCIP_MAXSTRLEN, "%s_09.sol", probname);
      SCIP_CALL( SCIPreadSolFile(scip, solutionfile, readsol_09, TRUE, NULL, &error));

      (void) SCIPsnprintf( solutionfile, SCIP_MAXSTRLEN, "%s_17.sol", probname);
      SCIP_CALL( SCIPreadSolFile(scip, solutionfile, readsol_17, TRUE, NULL, &error));

      assert( error == FALSE );

      SCIP_CALL( SCIPallocBufferArray(scip, &solvalint, 3*(p+1)));

      // readsol_05
      cout << "read_05" << endl;
      for( i = 0; i < p + 1; i++ )
      {
         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "z_%d", i);
         var = SCIPfindVar(scip, varname);
         solvalint[i] = SCIPgetSolVal(scip, readsol_05, var);
      }

      // readsol_09
      cout << "read_09" << endl;
      for( i = 0; i < p + 1; i++ )
      {
         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "z_%d", i);
         var = SCIPfindVar(scip, varname);
         solvalint[i+p+1] = SCIPgetSolVal(scip, readsol_09, var);
      }

      // readsol_17
      cout << "read_17" << endl;
      for( i = 0; i < p + 1; i++ )
      {
         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "z_%d", i);
         var = SCIPfindVar(scip, varname);
         solvalint[i+2*(p+1)] = SCIPgetSolVal(scip, readsol_17, var);
      }

      // print solvalint
      cout << endl;
      ct = 0;
      for( i = 0; i < 3; i++ )
      {
         for( j = 0; j < p + 1; j++ )
         {
            cout << solvalint[ct++] << " ";
         }
         cout << endl;
      }

      //newton
      for( i = 0; i < 3; i++ )
      {
         SCIP_CALL( newton( scip, probdata, &solvalint[i*(p+1)]) );
      }

      SCIPfreeBufferArray(scip, &solvalint);
      exit(1);
   }


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

