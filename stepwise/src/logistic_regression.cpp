/* logistic_regression.h */
/*********************************************/
/*********************************************/

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include "call_cblas.h"

using namespace std;

// log likelifood function
double Loglikelifood_(
   SCIP*    scip,
   const int      n,
   const int      p,
   double         *c,      // [p]
   double         *Xb,     // [n]
   double         *b       // [p]
){
   /*
    * log likelifood function
    *
    *    b in R^p
    *       l(b) = <c,b> - sum_i=1^n( log(1+f_i(b)) )
    *
    *       where f_i(b) = exp( <b,x_i> )
    *             x_i is row vector of X (i=1,2,..,n)
    */

   int      i;
   double   f  = SCIPcblasDdot( c, b, p);

   for( i = 0; i < n; ++i )
   {
      if( Xb[i] < 13.0 )
         f -= log( 1.0 + exp( Xb[i] ));
      else
         f -= Xb[i];
   }

   return f;
}

