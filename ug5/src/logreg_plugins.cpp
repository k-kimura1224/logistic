
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "ug_scip/scipUserPlugins.h"
#include "ug_scip/scipParaSolver.h"
#include "ug_scip/scipParaInitiator.h"

#include "scip/scipdefplugins.h"
#include "reader_logreg.h"
#include "relax_newton.h"
#include "probdata_logreg.h"
#include "heur_backward.h"
//#include "heur_forward.h"
//#include "branch_frequent.h"
//#include "branch_myfullstrong.h"
//#include "set_myparameter.h"
using namespace UG;
using namespace ParaSCIP;

/* to provide rank ans size to user  */
static ParaComm *paraComm = 0;

class LineregUserPlugins : public ScipUserPlugins {
   void operator()(SCIP*	scip)
   {
      /* include steiner tree reader */
      SCIP_CALL_ABORT( SCIPincludeReaderLogreg(scip) );

		/* include relaxator */
		SCIP_CALL_ABORT( SCIPincludeRelaxNewton(scip) );

		/* include primal heuristics */
		//SCIP_CALL_ABORT( SCIPincludeHeurForward(scip));

		/* include primal heuristics */
		SCIP_CALL_ABORT( SCIPincludeHeurBackward(scip));

		/* include branching rule */
		//SCIP_CALL_ABORT( SCIPincludeBranchruleFrequent(scip));

		/* include branching rule */
		//SCIP_CALL_ABORT( SCIPincludeBranchruleMyfullstrong(scip));

		//SCIP_CALL_ABORT( SCIPsetMyParameter( scip ));
   }

   void writeUserSolution(SCIP *scip, int nSolvers, double dual)
   {
      //SCIPprobdataSetDualBound(scip, dual);
      //SCIPprobdataSetNSolvers(scip, nSolvers);
      //SCIPprobdataWriteLogfileEnd(scip);
   }

};


void
setUserPlugins(ParaInitiator *inInitiator)
{
   ScipParaInitiator *initiator = dynamic_cast<ScipParaInitiator *>(inInitiator);
   initiator->setUserPlugins(new LineregUserPlugins());
   paraComm = initiator->getParaComm();
}

void
setUserPlugins(ParaInstance *inInstance)
{
   ScipParaInstance *instance = dynamic_cast<ScipParaInstance *>(inInstance);
   instance->setUserPlugins(new LineregUserPlugins());
}

void
setUserPlugins(ParaSolver *inSolver)
{
   ScipParaSolver *solver = dynamic_cast<ScipParaSolver *>(inSolver);
   solver->setUserPlugins(new LineregUserPlugins());
   if( !paraComm )
   {
      paraComm = solver->getParaComm();
   }
}

extern "C"
int getUgRank()
{
   return paraComm->getRank();
}

extern "C"
int getUgSize()
{
   return paraComm->getSize();
}
