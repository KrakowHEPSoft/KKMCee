{
///============================================================================
///
/// This is configuration/initialization script for ../MainFOAM main program
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
double nevtot   = 1e10;
//nevtot =  5e6;
//nevtot = 20e6;
//nevtot = 100e6;
double nevgrp   = 5e5; // 500k
nevgrp = 1e6;
///------------------------------------------------------------------
SemFile.cd();
TSemaf *Semafor = new TSemaf(semaf, nevtot, nevgrp);
Semafor->Write("Semafor",TObject::kOverwrite);
SemFile.Write();
SemFile.Close();
///*****************************************************************
///       Create new instance of RN generator and initialize it
GenFile.cd();
TRandom *RN_gen = new TRandom3();       // Central r.n.gen. Mersene
//TRandom *RN_gen = new TRandom1();       // Central r.n.gen. Ranlux
long    iniseed = 54217137;
RN_gen->SetSeed(iniseed);
RN_gen->Write("RN_gen",TObject::kOverwrite);
///*****************************************************************
///      Create new instance of MC generator
TMCgenFOAM *MCgen = new TMCgenFOAM("MCgen");
//########### Overwrite default input parameters ###########
MCgen->m_nSampl  = 100000;   // MC evts/cell (200)  ##
MCgen->m_nCells  = 10000;
MCgen->m_IsFoam5 = 0;   // Foam5 OFF
MCgen->m_IsFoam3 = 0;   // Foam3 OFF
MCgen->m_IsFoam3i= 1;   // Foam3i ON
MCgen->m_IsFoam1 = 0;   // Foam1 OFF
MCgen->m_IsFoam2 = 0;   // Foam0 OFF
//##########################################################
MCgen->ls();
MCgen->Write("MCgen",TObject::kOverwrite);
///*****************************************************************
///       Create new instance of the MC analysis object
TRobol *RoboT = new TRobolFOAM("RoboT");  /// base clase only
RoboT.f_HistNormName = "HST_FOAM_NORMA5";
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
