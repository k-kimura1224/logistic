/**@file   branch_myrule.c
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "branch_frequent.h"
#include "probdata_logreg.h"
#include "convenient_tool.h"
#include "mysolpool.h"

#define BRANCHRULE_NAME            "frequent"
#define BRANCHRULE_DESC            "most frequent branching"
#define BRANCHRULE_PRIORITY        200000
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

#define debug  0

/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of branching rule
 */

/* TODO: Implement all necessary branching rule methods. The methods with an #if 0 ... #else #define ... are optional */


/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyFrequent)
{  /*lint --e{715}*/

   /* call inclusion method of branchrule */
   SCIP_CALL( SCIPincludeBranchruleFrequent(scip) ) ;

   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeFrequent)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitFrequent)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** deinitialization method of branching rule (called before transformed problem is freed) */
static
SCIP_DECL_BRANCHEXIT(branchExitFrequent)
{  /*lint --e{715}*/


   return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsFrequent)
{  /*lint --e{715}*/

   // for probdata
   SCIP_PROBDATA* probdata;
   int   p;
   int   p1;
   SCIP_VAR**  var_z;      /* [p] 01 variables */
   SolPool*    pool;

   // for branching
   SCIP_VAR**  cands;
   int   ncands;
   SCIP_NODE*  childnode_0;      /* z_j = 0 */
   SCIP_NODE*  childnode_1;      /* z_j = 1 */

   int*  list;       /* list of candidate variables */

   int   i,j;

#if debug
   printf("[myrule brnaching]");
   Longline();
#endif

   /* get problem data*/
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   p  =  SCIPprobdataGetNexvars(probdata);
   p1 =  p+1;
   var_z =  SCIPprobdataGetVars_z(probdata);
   pool  = SCIPprobdataGetPool( probdata );

   /* alloc */
   SCIP_CALL( SCIPallocBufferArray(scip, &list, p1));
   SCIP_CALL( SCIPinitIntArrayZero( p1, list ) );

   /* get pseudo candidates (non-fixed integer variables) */
   SCIP_CALL( SCIPgetPseudoBranchCands(scip, &cands, NULL, &ncands) );

   for( i = 0; i < ncands; i++ )
   {
      for( j = 0; j < p1; j++ )
      {
         if( cands[i] == var_z[j] )
         {
            list[j] = 1;
            break;
         }
         assert( j<=p );
      }
   }

   int ind = pool->MFB_pool.MFBSolPool_maxfreq( list );

#if debug
   printf("list:");
   printintv( p1, list);
   cout << "MFB: " << ind << "th var. " << endl;
#endif


   assert( ind >= 0 );
   assert( ind < p1 );
   assert( list[ind]==1 );

   SCIP_CALL( SCIPbranchVar( scip, var_z[ind], &childnode_0, NULL, &childnode_1));

   /* free */
   SCIPfreeBufferArray(scip, &list);

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/*
 * branching rule specific interface methods
 */

/** creates the myrule branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleFrequent(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULE* branchrule;



   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, NULL) );

   assert(branchrule != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyFrequent) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeFrequent) );
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitFrequent) );
   SCIP_CALL( SCIPsetBranchruleExit(scip, branchrule, branchExitFrequent) );
   SCIP_CALL( SCIPsetBranchruleExecPs(scip, branchrule, branchExecpsFrequent) );

   return SCIP_OKAY;
}
