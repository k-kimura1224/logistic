
/**@file   reader_logreg.h
 * @author 
 *
 * This file implements the reader used to read and write ********** problems.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_LOGREG_H__
#define __SCIP_READER_LOGREG_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the stp file reader in SCIP */
extern
SCIP_RETCODE SCIPincludeReaderLogreg(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
