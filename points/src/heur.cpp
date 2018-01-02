#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>

#include "svpsolver.h"
#include "probdata.h"
#include "vector.h"
#include "solution.h"
#include "solution_pool.h"

using namespace std;

// heur_unitsphere.cpp
#define	HEUR_UNITSPHERE				false
#define	HEUR_UNITSPHERE_FREQ			1
#define	HEUR_UNITSPHERE_FREQOFS		0
#define	HEUR_UNITSPHERE_DPT			40

#define	PARA_HEUR_UNITSPHERE					false
#define	PARA_HEUR_UNITSPHERE_FREQ			5
#define	PARA_HEUR_UNITSPHERE_FREQOFS		20
#define	PARA_HEUR_UNITSPHERE_DPT			20

// heur_quadratic.cpp
#define	HEUR_QUADRATIC				true
#define	HEUR_QUADRATIC_FREQ			1
#define	HEUR_QUADRATIC_FREQOFS		0
#define	HEUR_QUADRATIC_DPT			40

#define	PARA_HEUR_QUADRATIC					true
#define	PARA_HEUR_QUADRATIC_FREQ			5
#define	PARA_HEUR_QUADRATIC_FREQOFS		20
#define	PARA_HEUR_QUADRATIC_DPT			20

#define debug	0

void	SVPsolver::heur(int sel, bool para)
{

	assert( sel >= 0 );
	assert( sel < listsize );

	list<NODE>::iterator it = NodeList.begin();
	advance( it, sel);
	int	dpt = it->get_dpt();
	int	dpt2;

	assert( dpt >= 0 );

	if( para == false ){
		// single
		if( it->get_zero() == false ){
			// unit sphere algorithm
			dpt2 = dpt - (int)HEUR_UNITSPHERE_FREQOFS;
			if( dpt2 >= 0
				&& HEUR_UNITSPHERE == true 
				&& dpt2 % (int)HEUR_UNITSPHERE_FREQ == 0 
				&& dpt2 <= (int)HEUR_UNITSPHERE_DPT ){
				heur_unitsphere( sel );
			}
			// solve quadratic optimization with one variable
			dpt2 = dpt - (int)HEUR_QUADRATIC_FREQOFS;
			if( dpt2 >= 0
				&& HEUR_QUADRATIC == true
				&& dpt2 % (int)HEUR_QUADRATIC_FREQ == 0 
				&& dpt2 <= (int)HEUR_QUADRATIC_DPT ){
				heur_quadratic( sel );
			}
		}
	}else{
		// parallel
		if( it->get_zero() == false ){
			// unit sphere algorithm
			dpt2 = dpt - (int)PARA_HEUR_UNITSPHERE_FREQOFS;
			if( dpt2 >= 0 
				&& PARA_HEUR_UNITSPHERE == true 
				&& dpt2 % (int)PARA_HEUR_UNITSPHERE_FREQ == 0 
				&& dpt2 <= (int)PARA_HEUR_UNITSPHERE_DPT ){
				heur_unitsphere( sel );
			}
			// solve quadratic optimization with one variable
			dpt2 = dpt - (int)PARA_HEUR_QUADRATIC_FREQOFS;
			if( dpt2 >= 0
				&& PARA_HEUR_QUADRATIC == true
				&& dpt2 % (int)PARA_HEUR_QUADRATIC_FREQ == 0 
				&& dpt2 <= (int)PARA_HEUR_QUADRATIC_DPT ){
				heur_quadratic( sel );
			}
		}
	}
}

