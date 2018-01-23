
/**@file   heur_forward.h
 * @brief  Forward selection
 *
 * Forward selection involves starting with all variables in the model.
 * 
 *
 */

#ifndef __SCIP_HEUR_FORWARD_H__
#define __SCIP_HEUR_FORWARD_H__


#include "scip/scip.h"
#ifdef __cplusplus
extern "C" {
#endif

/** creates the local primal heuristic and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeHeurForward(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
