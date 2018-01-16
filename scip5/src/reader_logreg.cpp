
/**@file   reader_logreg.c
 * @brief  **** problem reader file reader
 * @author 
 *
 * This file implements the reader used to read and write ******* problems.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "probdata_logreg.h"
#include "reader_logreg.h"

#define READER_NAME             "logregreader"
#define READER_DESC             "file reader for linear regression data format"
#define READER_EXTENSION        "logreg"

#define   DEFAULT_COUNTPRESOLTIME  TRUE      /**< count presolving time as part of overall solution time? */

/**@name Callback methods
 *
 * @{
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyLogreg)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderLogreg(scip) );

   return SCIP_OKAY;
}

/** problem reading method of the reader */
static
SCIP_DECL_READERREAD(readerReadLogreg)
{  /*lint --e{715}*/
   SCIP_RETCODE          retcode;
   SCIP_PROBDATA*        probdata;

   *result = SCIP_DIDNOTRUN;

   retcode = SCIPprobdataCreate(scip, filename);

   if( retcode == SCIP_READERROR )
      return SCIP_READERROR;

   SCIP_CALL( retcode );

   probdata = SCIPgetProbData(scip);
   if( SCIPgetStage(scip) == SCIP_STAGE_INIT ||  probdata == NULL )
      return SCIP_READERROR;

#if 0
		SCIPinfoMessage(scip, NULL, "Original problem:\n");
		SCIP_CALL( SCIPprintOrigProblem(scip, NULL, NULL, FALSE) );
   	SCIPinfoMessage(scip, NULL, "\n");
#endif

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

/** problem writing method of the reader */
static
SCIP_DECL_READERWRITE(readerWriteLogreg)
{  /*lint --e{715}*/

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** includes the stp file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderLogreg(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create reader data */
   readerdata = NULL;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );
   assert(reader != NULL);

   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyLogreg) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadLogreg) );
   SCIP_CALL( SCIPsetReaderWrite(scip, reader, readerWriteLogreg) );

   /* include user parameters */

	SCIP_CALL( SCIPaddBoolParam(scip,
         "logreg/countpresoltime",
        "count presolving time to solving time?",
         NULL, FALSE, DEFAULT_COUNTPRESOLTIME, NULL, NULL) );

   return SCIP_OKAY;
}

/**@} */
