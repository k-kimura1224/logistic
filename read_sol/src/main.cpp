
/**@file   main.cpp
 */
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "scip/scip.h"
#include "scip/scipshell.h"
#include "scip/scipdefplugins.h"

#include "reader_logreg.h"
#include "probdata_logreg.h"
#include "relax_newton.h"
#include "heur_backward.h"
#include "heur_forward.h"
#include "branch_frequent.h"


/* namespace usage */
using namespace std;
/** creates a SCIP instance with default plugins, evaluates command line parameters, runs SCIP appropriately,
 *  and frees the SCIP instance
 */

static
SCIP_RETCODE OutputSolution(
	SCIP*				scip
)
{

	SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
	//int	n	=	SCIPprobdataGetNdatas(probdata);
	int	p1	=	SCIPprobdataGetNexvars(probdata) + 1;
	SCIP_SOL*	sol;
	SCIP_VAR*	var;
	char			varname[SCIP_MAXSTRLEN];
	SCIP_Real 	solval;
	int			solvalint;
	SCIP_Real 	objval;
	int			ct=0;

	int i;

   assert(probdata != NULL);

   printf("\nSolution:\n");

	/* print a part of solution */

	objval = SCIPgetSolOrigObj( scip, SCIPgetBestSol(scip));
	sol = SCIPgetBestSol(scip);

	printf("\tz\t\t\tb\n\n");
	for(i=0; i<p1; i++){
		(void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "b%d", i);
		var			= SCIPfindVar(scip, varname);
		solval		= SCIPgetSolVal(scip, sol, var);
		(void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "z%d", i);
		var			= SCIPfindVar(scip, varname);
		solvalint 	= SCIPgetSolVal(scip, sol, var);
		printf("%2d:\t%d\t%20.15g\n", i, solvalint, solval);
		if( solvalint==1 ) ct++;

	}

	cout << endl;
	printf("k: %d\t", ct);
	printf("objval: %20.15g\n\n", objval);

	return SCIP_OKAY;
}

#if 0
static
SCIP_RETCODE OutputSolutionsN(
	SCIP*				scip
)
{

	SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
	int	p	=	SCIPprobdataGetNexvars(probdata);
	SCIP_VAR**	var_z	= SCIPprobdataGetVars_z(probdata);
	int	z_val;
	int*	n_z;			/* the number of times that z is 1 */
	SCIP_SOL**	sols;
	int			nsols = SCIPgetNSols(scip);

	int			N = 100;

	int i,j;

   assert(probdata != NULL);

   printf("\nSolution:\n");

	if( N > nsols ) N = nsols;

	if( N==0 ) return SCIP_OKAY;

	sols = SCIPgetSols(scip);

	SCIP_CALL( SCIPallocBufferArray(scip, &n_z, p));
	GenerateZeroVecInt( p, n_z);

	for(i=0; i<N; i++){
		for(j=0; j<p; j++){
				z_val = SCIPgetSolVal( scip, sols[i], var_z[j]);
				if( z_val == 1 ){
					n_z[j]++;
				}
		}
	}

	for(i=0; i<p; i++){
		printf("%d -- %d \n", i+1, n_z[i]);
	}

	stop();


	return SCIP_OKAY;
}
#endif

static
SCIP_RETCODE runShell(
   int                   argc,               /**< number of shell parameters */
   char**                argv,               /**< array containing shell parameters */
   const char*           defaultsetname      /**< name of default settings file */
   )
{
   SCIP* scip = NULL;

   /*********
    * Setup *
    *********/

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* include steiner tree reader */
   SCIP_CALL( SCIPincludeReaderLogreg(scip) );

   /* include relaxator */
	SCIP_CALL( SCIPincludeRelaxNewton(scip) );

	/* include primal heuristics */
	SCIP_CALL( SCIPincludeHeurBackward(scip));

	/* include primal heuristics */
	SCIP_CALL( SCIPincludeHeurForward(scip));

	///* include branching rule */
	//SCIP_CALL( SCIPincludeBranchruleMyfullstrong(scip));

	///* include branching rule */
	SCIP_CALL( SCIPincludeBranchruleFrequent(scip));

   /* set LOGREG-specific default parameters */
   /* propagating */
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/pseudoobj/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/vbounds/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/dualfix/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/linear/propfreq", -1) );
   /* presolving */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );
   /* if -1, disable the LP relaxation and only use my custom relaxation */
   SCIP_CALL( SCIPsetIntParam(scip, "lp/solvefreq", -1) );
    /**********************************
    * Process command line arguments *
    **********************************/

   SCIP_CALL( SCIPprocessShellArguments(scip, argc, argv, defaultsetname) );

   if( SCIPgetNSols(scip) > 0 )
   {
		SCIP_CALL( OutputSolution(scip));
		//if( 0 )
		//	OutputSolutionsN(scip);
	}

   /********************
    * Deinitialization *
    ********************/

   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

int
main(
   int                   argc,               /**< number of shell parameters */
   char**                argv                /**< array containing shell parameters */
   )
{
   SCIP_RETCODE retcode;

   retcode = runShell(argc, argv, "scip.set");
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
