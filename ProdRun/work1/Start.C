{
///============================================================================
///
/// This is configuration/initialization script for ../MainKKMC main program
/// To start MC run in the interactive mode just type "make start" 
///
/// This series of calculations integrates a given distribution 
/// by using MC methods
///============================================================================
gROOT->Reset();
cout<<"%%% ================== Start.C ================== %%%%"<<endl;
//gSystem->Load("../.libs/libKKmc.so"); // transferred to rootlogon.C
TFile HistoFile("histo.root","RECREATE","Histograms");
TFile GenFile(  "mcgen.root","RECREATE","r.n.generator, MCgens");
TFile SemFile(  "semaf.root","RECREATE","Semaphore");
///*****************************************************************
///   Create new instance of Semaphre object
///   and fill it with the MC run general parameters
TString semaf   = "START";
double nevtot   = 1e12; // 1000G
//nevtot = 2e5;
//nevtot =1000;
nevtot = 10e6;
double nevgrp   = 1e5; // 100k
nevgrp = 1e6;          // 1M
//nevgrp = 1e4;
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
KKee2f *MCgen = new KKee2f("MCgen");
//########### Change some input parameters ###########
//MCgen->m_nSampl  = 100000;   // MC evts/cell (200)  ##
//####################################################
MCgen->ls();
MCgen->Write("MCgen",TObject::kOverwrite);
///*****************************************************************
///       Create new instance of the MC analysis object
TRobol *RoboT = new TRobolKKMC("RoboT");  /// base clase only
RoboT->f_HistNormName = "HST_KKMC_NORMA";
RoboT->Write("RoboT",TObject::kOverwrite);
///*****************************************************************
GenFile.Write();
cout<<"--------------------------GenFile.ls-------------------------------"<<endl;
GenFile.ls();
cout<<"-------------------------GenFile.ShowStreamerInfo-------------------------"<<endl;
GenFile.ShowStreamerInfo();
GenFile.Close();
cout << "===========Output written in histo.root===========" << endl;
HistoFile.Write();
//HistoFile.ls();
HistoFile.Close();
cout<<"%%% ===============End Start.C ================== %%%%"<<endl;
return 0;
}
