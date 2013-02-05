///////////////////////////////////////////////////////////////////////////////
//             Weelcome to Bose-Einstein Exercise !!!                        //
///////////////////////////////////////////////////////////////////////////////
// To initiate enviroment for ROOT  source ./RootLogin
// To compile&link:                 make rmain.exe
// To execute:
//	make 172GeV.4J-start
//	make 172GeV.4J-stop
//	make 172GeV.4J-cont
//	make 172GeV.4J-debug
// To debug:                        xldb rmain.exe
//  set -g option in compilation everywhere before debug!!!
//  xldb -I "../korww ../interfaces ../ampli4f ../kwlund ../model" rmain.exe
//
//  NOTES:
//  FORTRAN is hidden in KoralwMaker.C and Semaph.C and koralw_fortface.f
//  this has to be reorganized!
///////////////////////////////////////////////////////////////////////////////
// C++ headers
#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include <math.h>
// ROOT headers
#include "TROOT.h"
#include "TFile.h"
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
  long NGroup = 500;        // number of events in the group
  long iLoop;               // auxiliary event counter

  cout << "   |--------------------| <<endl<<flush" ;
  cout << "   |  Rmain.C  Entering | <<endl<<flush" ;
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
      cout <<"iEvent = "<<iEvent<<endl<<flush;
      ///////////////////////////////////////////////////////////////
      // Checks for EACH GROUP of NGroup events of semaphore etc.  //
      // Dump of histos on disk erasing trailing cycles            //
      ///////////////////////////////////////////////////////////////
      DiskFile.Write();
      char chcyc[100]; sprintf(chcyc,"*;%i",iLoop);
      Text_t *Tcycle = chcyc;
      DiskFile.Delete(Tcycle);
      //
      MainSemaphore.ReadStatus(Semaph1);
      if( iEvent  >= NevTot )      break;  // production stoped
      if( Semaph1 == "STOP" )      break;  // production stoped
      ///////////////////////////////////////////////////////////////
    } // for iLoop

  DiskFile.ls();
  DiskFile.Close();
  
  Robol.Finalize();
  /////////////////
Finis:
  cout << "   |--------------------| \n" ;
  cout << "   |  Rmain.C  Exiting  | \n" ;
  cout << "   |--------------------| \n" ;
}
