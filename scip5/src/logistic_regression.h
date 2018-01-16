/* logistic_regression.h */
/*********************************************/
/*********************************************/
#ifndef LOGISTIC_F_H
#define LOGISTIC_F_H

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>

//#ifdef __cplusplus
//extern "C" {
//#endif

// log likelifood function 
extern
double Loglikelifood_(	
   SCIP* scip,
	const int 		n,
	const int 		p,
	double		 	*c,		// [p]
	double 			*X_,		// [n*p] (colMajor)!!!!!!!
	double 			*b			// [p]
);


//#ifdef __cplusplus
//}
//#endif

#endif
