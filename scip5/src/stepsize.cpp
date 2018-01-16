/* stepsize.h */
/*********************************************/
/*********************************************/

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <assert.h>

#include "stepsize.h"
#include "call_cblas.h"
#include "logistic_regression.h"

#define debug  0

using namespace std;


static
SCIP_Real   StepsizeMyrule(
   SCIP*          scip,
   const int      n,
   const int      dim,
   SCIP_Real      *coef,      // coef[dim]
   SCIP_Real      *X_,        // X[n*dim] (colMajor)
   SCIP_Real      *x,         // x[dim] is current point
   SCIP_Real      *d       // d[dim] is descent direction
   )
{
   /* find step size
    *  for minimizing f(x+alpha*d)
    *
    */

   SCIP_Real      alpha;
   SCIP_Real      *x_new;
   SCIP_Real      objval = SCIPinfinity(scip);
   SCIP_Real      objval_new;
   SCIP_Real      old = 0.0;
   SCIP_Real      t = 0.5;
   SCIP_Real      *Xb;

   alpha =  10.0;
   SCIP_CALL( SCIPallocBufferArray(scip, &x_new, dim));
   SCIP_CALL( SCIPallocBufferArray(scip, &Xb, n));

   int ct = 0;

   while(1)
   {
      ct++;

      // x_new = x + alpha d
      SCIP_CALL( SCIPcblasDaxpy( x, d, dim, 1.0, alpha, x_new) );
      // Xb = X x_new
      SCIP_CALL( SCIPcblasDgemv2( X_, n, dim, x_new, Xb) );

      objval_new = - Loglikelifood_( scip, n, dim, coef, Xb, x_new);

#if debug
   cout << "stepsize: " << alpha << endl;
   cout << "new objval: " << objval_new << endl;
#endif

      // debug
      //cout << "alpha=" << alpha << endl;
      //cout << "objval=" << objval_new << endl;

      if( (objval_new > objval) && objval_new < 10.0e+08 )
      {
         break;
      }
      else
      {
         alpha    *= t;
         old      =  objval;
         objval   =  objval_new;
      }

      if( ct >= 1000 )
      {
         cout << "error" << ct << ": my rule" << endl;
         exit(1);
         //return (SCIP_Real)-ct;
      }

   }

   if( ct < 3 )
   {
      //cout << "error" << ct << ": my rule" << endl;
      //cout << "objval_new: " << objval_new << endl;
      //cout << "objval:" << objval << endl;
      //exit(1);
      SCIPfreeBufferArray(scip, &x_new);
      SCIPfreeBufferArray(scip, &Xb);
      return   alpha/t;
   }

   //cout << "ct : " << ct_total << endl;

   SCIP_Real   a,A,B,C,lag,L;
   a     =  (alpha/(t*t));
   A     =  old;
   B     =  objval;
   C     =  objval_new;
   lag   =  3*a*(A-5*B+4*C);
   lag   /= 8*(A-3*B+2*C);
   SCIP_CALL( SCIPcblasDaxpy( x, d, dim, 1.0, lag, x_new) );
   SCIP_CALL( SCIPcblasDgemv2( X_, n, dim, x_new, Xb) );
   L     =   - Loglikelifood_( scip, n, dim, coef, Xb, x_new);

   //cout << a << ": " << A << endl;
   //cout << a*t << ": " << B << endl;
   //cout << a*t*t << ": " << C  << endl;
   //cout << lag << ": " << L << endl;

   SCIPfreeBufferArray(scip, &x_new);
   SCIPfreeBufferArray(scip, &Xb);

   if( B > L ){
      //cout << "L" << endl;
      return lag;
   }else{
      //cout << "B" << endl;
      return a*t;
   }

}

// find the stepsize
SCIP_Real FindStepsize(
   SCIP*          scip,
   const int      n,
   const int      dim,
   SCIP_Real      *coef,      // coef[dim]
   SCIP_Real      *X_,        // X[n*dim] (colMajor)
   SCIP_Real      *x,         // x[dim] is the current point
   SCIP_Real      *d          // d[dim] is the descent direction
){
   SCIP_Real ss;

   if( MP_STEPSIZE_MODE == 0 )
   {
      // if 0, this mode is my rule
#if debug
   cout << "[StepsizeMyrule]start" << endl;
#endif

      ss =  StepsizeMyrule( scip, n, dim, coef, X_, x, d);

#if debug
   cout << "[StepsizeMyrule]end" << endl;
#endif

      return   ss;

   }else{
      cout << "error: stepsize.cpp" << endl;
      exit(1);
   }


   return   -1.0;
}

