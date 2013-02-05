///////////////////////////////////////////////////////////////////////////////
// make demoC-start
///////////////////////////////////////////////////////////////////////////////
// This module handles main Monte Carlo loop and writing histogram onto disk.
// It is UNIVERSAL, works for any MC generator for event analysis.
// MC event generator and analysis is encapsulated inside Robol module.
// Provisions are implemented for running parallel jobs on a farm.
// In particular: random number seed is read from disk and Semaphore file
// remembers r.n. seed for restarting job at the point it was finished last time.
// NGroup variable decides how often events are dumped on the disk.
///////////////////////////////////////////////////////////////////////////////
// C++ headers
#include <stdlib.h>
#include <fstream.h>
#include <iomanip.h>
#include <iostream.h>
#include <math.h>
// ROOT headers
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
// OUR headers
#include "Semaph.h"
#include "ROBOL.h"

//=============================================================================
//  ROOT  ROOT ROOT   ROOT  ROOT  ROOT  ROOT  ROOT  ROOT  ROOT   ROOT   ROOT 
//=============================================================================
TROOT tuples("rmain","Main production program for B.E. exercise");
TFile DiskFile("rmain.root","RECREATE","histograms");
//=============================================================================
main()
{
  long NevTot = 1000000000; // total (limit) of events
  long iEvent = 0;          // serial number of event in the MC generation
  long NGroup = 50000;      // number of events in the group
  long iLoop;               // auxiliary event counter
  char chcyc[100];          // Cycle text variable
  Text_t *Tcycle;           // Cycle text variable

  cout << "   |--------------------| <<endl<<flush" ;
  cout << "   |  DemoC    Entering | <<endl<<flush" ;
  cout << "   |--------------------| <<endl<<flush" ;

  ROBOL Robol; // define main Maker (SuperFactory)
  ////////////

  ///////////////////////////////////////////////////////////////////////
  // Intializations of the random number generator, seed from the disk //
  ///////////////////////////////////////////////////////////////////////
  Semaph  MainSemaphore; // define our main semaphore utility
  TString Semaph1;
  MainSemaphore.Initialize(Semaph1);
  if( Semaph1 == "STOP") goto Finis;  // STOP for wrong semaphore
  

  Robol.Initialize(NevTot);
  /////////////////////////

  cout<< "||||||||||||||||||||||||||||||||||||||||||||||"<<endl;
  cout<< "||   Going to generate "<<NevTot<< "  events  "<<endl;
  cout<< "||||||||||||||||||||||||||||||||||||||||||||||"<<endl;

  for ( iLoop = 0; ; iLoop++)
    {
      for ( long iGroup = 0; iGroup< NGroup; iGroup++) 
	{
	  iEvent++;
	  
	  Robol.Production(iEvent);
	  /////////////////////////

	  if( iEvent  >= NevTot )  break;  // production stoped
	}; // for iGroup
      //
      Robol.KKMC_NORMA();                  // transfering normalization histog.
      ///////////////////
      cout <<"iEvent = "<<iEvent<<endl<<flush;
      ///////////////////////////////////////////////////////////////
      // Checks for EACH GROUP of NGroup events of semaphore etc.  //
      // Dump of histos on disk erasing trailing cycles            //
      ///////////////////////////////////////////////////////////////
      DiskFile.Write();
      sprintf(chcyc,"*;%i",iLoop);
      Tcycle = chcyc;
      DiskFile.Delete(Tcycle);
      //
      MainSemaphore.ReadStatus(Semaph1);
      if( iEvent  >= NevTot )      break;  // production stoped
      if( Semaph1 == "STOP" )      break;  // production stoped
      ///////////////////////////////////////////////////////////////
    } // for iLoop
  //
  Robol.KKMC_NORMA();                  // transfering normalization histog.
  Robol.Finalize();
  //---  last write after Robol
  DiskFile.Write();
  sprintf(chcyc,"*;%i",iLoop+1);
  Tcycle = chcyc;
  DiskFile.Delete(Tcycle);
  //
  DiskFile.ls();
  DiskFile.Close();
  
Finis:
  cout << "   |--------------------| \n" ;
  cout << "   |  DemoC    Exiting  | \n" ;
  cout << "   |--------------------| \n" ;
}
