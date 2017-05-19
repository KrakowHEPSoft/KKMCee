{
///============================================================================
///
/// This is configuration/initialization script for ../MainPr main program
/// To start MC run in the interactive mode just type "make start" 
///
/// This series of calculations integrates a given distribution 
/// by using MC methods
///============================================================================
gROOT->Reset();
cout<<"%%% ================== Start.C ================== %%%%"<<endl;
gSystem->Load("../.libs/libKKfm.so");
TFile HistoFile("histo.root","RECREATE","Histograms");
TFile GenFile(  "mcgen.root","RECREATE","r.n.generator, MCgens");
TFile SemFile(  "semaf.root","RECREATE","Semaphore");
///*****************************************************************
///   Create new instance of Semaphre object
///   and fill it with the MC run general parameters
TString semaf   = "START";
double nevtot   = 1e8;
double nevgrp   = 2e7; // 200k
///------------------------------------------------------------------
SemFile.cd();
TSemaf *Semafor = new TSemaf(semaf, nevtot, nevgrp);
Semafor->Write("Semafor",TObject::kOverwrite);
SemFile.Write();
SemFile.Close();
///*****************************************************************
///       Create new instance of RN generator and initialize it
GenFile.cd();
TRandom *RN_gen = new TRandom3();       // Central r.n.gen.
long    iniseed = 54217137;
RN_gen->SetSeed(iniseed);
RN_gen->Write("RN_gen",TObject::kOverwrite);
///*****************************************************************
///      Create new instance of MC generator
TMCgenFOAM *MCgen = new TMCgenFOAM("MCgen");
//########### Change some input parameters ###########
MCgen->m_kDim    =  3;       // energy spread  ON   ##
//MCgen->m_kDim    =  2;       // energy spread  OFF  ##
MCgen->m_KeyISR  =  2;       // Type of QED ISR corr##
MCgen->m_nCells  =  10000;    // No. of cells, (2k)  ##
MCgen->m_nSampl  =  10000;    // MC evts/cell (200)  ##
MCgen->m_sigE    =  0.008;   // energy spread [GeV] ##
MCgen->m_MH      =  125.7;   // Higgs energy  [GeV] ##
MCgen->m_GamH    =  0.0042;  // Higgs width   [GeV] ##
//####################################################
MCgen->ls();
MCgen->Write("MCgen",TObject::kOverwrite);
///*****************************************************************
///       Create new instance of the MC analysis object
TRobol *RoboT = new TRobolFOAM("RoboT");  /// base clase only
RoboT->Write("RoboT",TObject::kOverwrite);
///*****************************************************************
GenFile.Write();
cout<<"---------------------------------------------------------"<<endl;
GenFile.Close();
cout << "===========Output written in histo.root===========" << endl;
HistoFile.Write();
//HistoFile.ls();
HistoFile.Close();
cout<<"%%% ===============End Start.C ================== %%%%"<<endl;
return 0;
}
