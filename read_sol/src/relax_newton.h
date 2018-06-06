
/**@file   relax_dposv.h
 * @ingroup RELAXATORS
 * @brief  myrelaxator relaxator by dposv
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_RELAX_MYRELAXATOR_NEWTON_H__
#define __SCIP_RELAX_MYRELAXATOR_NEWTON_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the myrelaxator relaxator and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludeRelaxNewton(
	SCIP*						scip
   );

#ifdef __cplusplus
}
#endif

#endif
