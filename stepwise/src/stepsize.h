/* stepsize.h */
/*********************************************/
/*********************************************/
#ifndef STEPSIZE_F_H
#define STEPSIZE_F_H

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>

#include "scip/scip.h"
#define MP_STEPSIZE_MODE 0
//#ifdef __cplusplus
//extern "C" {
//#endif

// log likelifood function 
extern
SCIP_Real FindStepsize(
   SCIP*          scip,
	const	int		n,
	const	int		dim,
	SCIP_Real			*coef,		// coef[dim]
	SCIP_Real			*X_,			// X[n*dim] (colMajor)
	SCIP_Real			*x,			// x[dim] is the current point
	SCIP_Real			*d				// d[dim] is the descent direction
);


//#ifdef __cplusplus
//}
//#endif

#endif
