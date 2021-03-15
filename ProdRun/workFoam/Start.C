{
///============================================================================
///
/// This is configuration/initialization script for the main program
/// To start MC run in the interactive mode just type "make start" 
///============================================================================
gROOT->Reset();
cout<<"%%% ================== Start.C ================== %%%%"<<endl;
gSystem->Load("../.libs/libProdFoam.so");
gSystem->Load("../../SRChh/.libs/libKK2f.so");
gSystem->Load("../../SRChh/.libs/libKKfm.so");
//gSystem->Load("../../MCdev/.libs/libMCdev.so");
TFile HistoFile("histo.root","RECREATE","Histograms");
TFile GenFile(  "mcgen.root","RECREATE","r.n.generator, MCgens");
TFile SemFile(  "semaf.root","RECREATE","Semaphore");
///*****************************************************************
///  MC run general parameters
double nevtot   = 1e12; // 1000G
//nevtot = 2e5;
//nevtot =10;     // 10
//nevtot =100;     // 100
//nevtot =2000;      //2k
//nevtot =10000;     // 10k
//nevtot =100000;    // 100k
//nevtot =200000;    // 200k
//nevtot =1000000;   // 1M
//nevtot =4000000;   // 4M
//nevtot =100e6;     // 100M
nevtot = 2e9;        // 2G
double nevgrp   = 2e5; // 200k
//nevgrp = 1e6;          // 1M
//nevgrp = 5e4;
///------------------------------------------------------------------
cout<<"***   Create new instance of Semaphore object"<<endl;
TString semaf   = "START";
SemFile.cd();
TSemaf *Semafor = new TSemaf(semaf, nevtot, nevgrp);
Semafor->Write("Semafor",TObject::kOverwrite);
SemFile.Write();
SemFile.Close();
///*****************************************************************
cout<<"***   Create new instance of RN generator and initialize it"<<endl;
GenFile.cd();
TRandom *RN_gen = new TRandom3();  // RanMar Marsaglia
//TRandom *RN_gen = new TRandom1();  // Ranlux
long    iniseed = 54217137;
RN_gen->SetSeed(iniseed);
RN_gen->Write("RN_gen",TObject::kOverwrite);
///*****************************************************************
cout<<"***   Create new instance of MC generator"<<endl;
TMCgen *MCgen = new KKhhFoam("MCgen");
//########### Change some input parameters ###########
//MCgen->m_nSampl  = 100000;   // MC evts/cell (200)  ##
//####################################################
MCgen->ls();
MCgen->Write("MCgen",TObject::kOverwrite);
///*****************************************************************
cout<<"***   Create new instance of the MC analysis object"<<endl;
TRobol *RoboT = new TRobolFoam("RoboT");  /// base clase only
//RoboT.f_HistNormName = TString("HST_FOAM_NORMA7"); // instead of h_TMCgen_NORMA
RoboT->Write("RoboT",TObject::kOverwrite);
///*****************************************************************
GenFile.Write();
cout<<"---------------------------------------------------------"<<endl;
GenFile.Close();
cout << "===========Output written in histo.root===========" << endl;
HistoFile.Write();
HistoFile.ls();
HistoFile.Close();
cout<<"%%% ===============End Start.C ================== %%%%"<<endl;
return 0;
}
