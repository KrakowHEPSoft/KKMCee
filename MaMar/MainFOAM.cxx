///////////////////////////////////////////////////////////////////////////////
// This module handles main Monte Carlo loop and writing histograms onto disk.
// It is UNIVERSAL, works for any MC generator for event analysis.
// MC event generator and analysis is encapsulated inside RobolA module.
// Three different disk files are used
///////////////////////////////////////////////////////////////////////////////
// C++ headers
using namespace std;
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

// ROOT headers
#include "TROOT.h"
#include "TFile.h"
//#include "TH1.h"

// OUR headers
#include "TSemaf.h"
#include "TRobol.h"
#include "TSystem.h"


int main()
{
  gSystem->Load("../.libs/libKKfm.so");
  TFile *HistoFile;    // ROOT objects, histograms etc.
  TFile *MCgenFile;    // Monte Carlo generator objects NEVER CLOSE!
  TFile *SemafFile;    // Semaphore and loop params
  HistoFile =TFile::Open("histo.root","UPDATE","Histograms");
  MCgenFile =TFile::Open("mcgen.root","UPDATE","Generators");
  SemafFile =TFile::Open("semaf.root","UPDATE","Semaphor");
  ofstream   OutFile("pro.output",ios::out);  // Logfile output
  cout    <<endl<<flush;
  cout    << "  |--------------------| "<<endl<<flush;
  cout    << "  |  MainPr   Entering | "<<endl<<flush;
  cout    << "  |--------------------| "<<endl<<flush;
  OutFile <<endl<<flush;
  OutFile << "  |--------------------| "<<endl<<flush;
  OutFile << "  |  MainPr   Entering | "<<endl<<flush;
  OutFile << "  |--------------------| "<<endl<<flush;
  double NevTot =2e14; // Total (limit) of events, redefined below
  double NGroup =2e5;  // no. of events in the group, TO BE TUNED
  double iLoop;        // external loop event counter
  double iGroup;       // internal loop event counter
  double iEvent = 0;   // serial number of event in the MC generation
  //////////////////////////////////////////////////////////////////////////
  TString Status = "START";   // to be overwitten
  //------------------------------------------------------------------------
  SemafFile->cd();
  SemafFile->ls();
  TSemaf *Semafor = (TSemaf*)SemafFile->Get("Semafor"); /// read semafor from the disk
  Status   = Semafor->m_flag;
  NevTot   = Semafor->m_nevtot;
  NGroup   = Semafor->m_nevgrp;
  cout<<"  MainPr:  Status  = "<<  Status << endl;
  cout<<"  MainPr:  NevTot  = "<<  NevTot << endl;
  Semafor->m_flag = "CONTINUE";
  Semafor->Write("Semafor",TObject::kOverwrite);      /// write modified semaphore
  SemafFile->Write("",TObject::kOverwrite,0);         /// write modified semaphore
  SemafFile->ls();
  SemafFile->Close();
  //--------------------------------------------------------------------------
  if( Status == "STOP") cout<<"MainPr: +++Warning+++ Status=STOP in semaphore"<<endl;
  //////////////////////////////////////////////////////////////////////////
  TRobol *RoboT =(TRobol*)MCgenFile->Get("RoboT");   /// read analysis object
  RoboT->Initialize(
        &OutFile,   /// Central log-file for messages
        MCgenFile,  /// ROOT disk file for CRNG and MC gen.
        HistoFile); /// ROOT disk file for histograms
  /// immediate dump of results of initialization into the disk
  MCgenFile->Write("",TObject::kOverwrite,0);
  MCgenFile->Flush();
  cout<< "||||||||||||||||||||||||||||||||||||||||||||||||||"<<endl;
  cout<< "|| MainPr:  to be generated "<<NevTot<< " events  "<<endl;
  cout<< "||||||||||||||||||||||||||||||||||||||||||||||||||"<<endl;
  ///***********************************************************************
  ///          MAIN LOOP OVER MONTE CARLO EVENTS
  for ( iLoop = 0; ; iLoop++)
    {
      for ( iGroup = 0; iGroup< NGroup; iGroup++) 
	{
	  iEvent++;
	  ///***************************
	  RoboT->Production(iEvent); //
	  ///***************************
	  if( iEvent  >= NevTot )  break;  // production stoped
	}; // for iGroup
      //
      cout    <<"iEvent = "<<iEvent<<endl<<flush;
      OutFile <<"iEvent = "<<iEvent<<endl<<flush;
      ///////////////////////////////////////////////////////////////
      // Disk Write for all histos and r.n. generator
      RoboT->FileDump();
      ///////////////////////////////////////////////////////////////
      // Checks semaphore for each group of NGroup events
      SemafFile =TFile::Open("semaf.root","UPDATE","Semaphor");
      //SemafFile->ReOpen("UPDATE"); // Bug in ROOT?
      Semafor = (TSemaf*)SemafFile->Get("Semafor");  // read semafor from the disk
      Status   = Semafor->m_flag;                   // get flag
      Semafor->m_nevent = iEvent;                   // record no. of events
      Semafor->Write("Semafor",TObject::kOverwrite);// write modified semaphore
      SemafFile->Write("",TObject::kOverwrite);
      SemafFile->Close();
      //-------------------------------------------------------------
      if( iEvent  >= NevTot )      break;  // production stoped
      if( Status == "STOP" )       break;  // production stoped
      ///////////////////////////////////////////////////////////////
    } // for iLoop
  ///            End of MC LOOP
  ///************************************************************************
  SemafFile =TFile::Open("semaf.root","UPDATE","Semaphor");
  Semafor = (TSemaf*)SemafFile->Get("Semafor");
  Semafor->m_flag = "CONTINUE";  
  Semafor->Write("Semafor",TObject::kOverwrite);// write modified semaphore
  SemafFile->Write("",TObject::kOverwrite);
  SemafFile->Close();
  cout<< "||||||||||||||||||||||||||||||||||||||||||||||"<<endl;
  cout<< "||  MainPr: Generated "<<iEvent<< " events    "<<endl;
  cout<< "||||||||||||||||||||||||||||||||||||||||||||||"<<endl;
  //
  RoboT->Finalize();
  //////////////////////////////////////////////////////////////////////
  //  LAST WRITEs
  MCgenFile->cd();
  RoboT->Write("RoboT",TObject::kOverwrite);  // only for tests
  RoboT->FileDump();
  //////////////////////////////////////////////////////////////////////
  //      Some TESTS and control printouts
  cout<<"------------------------------ls----------------------------------"<<endl;
  HistoFile->ls();
  MCgenFile->Write();
  cout<<"------------------------------Map---------------------------------"<<endl;
  HistoFile->Map();
  cout<<"-------------------------ShowStreamerInfo-------------------------"<<endl;
  MCgenFile->ShowStreamerInfo();
  cout<<"------------------------GetListOfKeys-----------------------------"<<endl;
  HistoFile->GetListOfKeys()->Print();
  cout<<"------------------------------end---------------------------------"<<endl;
  //////////////////////////////////////////////////////////////////////
  //                 CLOSE ALL FILES
  HistoFile->Close();
  MCgenFile->Close();
  OutFile.close();
  cout    <<endl<<flush;
  cout    << "  |--------------------| "<<endl<<flush;
  cout    << "  |  MainPr   Ended    | "<<endl<<flush;
  cout    << "  |--------------------| "<<endl<<flush;
  OutFile <<endl<<flush;
  OutFile << "  |--------------------| "<<endl<<flush;
  OutFile << "  |  MainPr   Ended    | "<<endl<<flush;
  OutFile << "  |--------------------| "<<endl<<flush;
  return 0;
  }