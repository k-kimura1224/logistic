
/**@file   probdata_logreg.h
 * @brief  Problem data for stp problem
 * @author 
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROBDATA_LOGREG__
#define __SCIP_PROBDATA_LOGREG__

#include "scip/scip.h"
#include "mysolpool.h"

#ifdef __cplusplus
extern "C" {
#endif


/** sets up the problem data */
extern
SCIP_RETCODE SCIPprobdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< file name */
   );

/** returns the number of datas */
extern
int SCIPprobdataGetNdatas(
   SCIP_PROBDATA*			probdata
	);

/** returns the number of explain vars */
extern
int SCIPprobdataGetNexvars(
   SCIP_PROBDATA*			probdata
	);


/** returns the number of vars */
int SCIPprobdataGetNvars(
   SCIP_PROBDATA*			probdata
	);


/** returns the number of linedep groups */
extern
int SCIPprobdataGetNdep(
   SCIP_PROBDATA*			probdata
	);

/** returns y */
extern
SCIP_Real*	SCIPprobdataGety(
   SCIP_PROBDATA*			probdata
	);

/** returns X */
extern
SCIP_Real*	SCIPprobdataGetX(
   SCIP_PROBDATA*			probdata
	);

/** returns coef_obj */
extern
SCIP_Real*	SCIPprobdataGetCO(
   SCIP_PROBDATA*			probdata
	);

/** returns penalcf */
extern
SCIP_Real	SCIPprobdataGetPC(
   SCIP_PROBDATA*			probdata
	);

/** returns Mdep */
extern
int*	SCIPprobdataGetMdep(
   SCIP_PROBDATA*			probdata
	);

/** returns groupX */
extern
int*	SCIPprobdataGetgroupX(
   SCIP_PROBDATA*			probdata
	);

/** returns var_b */
extern
SCIP_VAR**	SCIPprobdataGetVars_b(
   SCIP_PROBDATA*			probdata
	);

/** returns var_z */
extern
SCIP_VAR**	SCIPprobdataGetVars_z(
   SCIP_PROBDATA*			probdata
	);

/** returns var_bx */
extern
SCIP_VAR**	SCIPprobdataGetVars_bx(
   SCIP_PROBDATA*			probdata
	);

/** returns var_EXP */
extern
SCIP_VAR**	SCIPprobdataGetVars_EXP(
   SCIP_PROBDATA*			probdata
	);

/** returns var_LOG */
extern
SCIP_VAR**	SCIPprobdataGetVars_LOG(
   SCIP_PROBDATA*			probdata
	);

/** returns pool */
extern
SolPool*	SCIPprobdataGetPool(
	SCIP_PROBDATA*			probdata
	);

/** set dual bound by ug */
extern
void SCIPprobdataSetDualBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             dual
   );

/** set the number of solvers */
extern
void SCIPprobdataSetNSolvers(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nSolvers            /**< the number of solvers */
   );

#ifdef __cplusplus
}
#endif

#endif
