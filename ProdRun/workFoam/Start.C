{
///============================================================================
///
/// This is configuration/initialization script for the main program
/// To start MC run in the interactive mode just type "make start" 
///============================================================================
gROOT->Reset();
cout<<"%%% ================== Start.C ================== %%%%"<<endl;
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
//nevtot =1000000;   // 1M
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
KKeeFoam *MCgen = new KKeeFoam("MCgen");
//####################################################
MCgen->ls();
MCgen->Write("MCgen",TObject::kOverwrite);
///*****************************************************************
cout<<"***   Create new instance of the MC analysis object"<<endl;
TRobol *RoboT = new TRobolFoam("RoboT");  /// base clase only
RoboT->f_HistNormName = TString("HST_FOAM_NORMA4"); // instead of h_TMCgen_NORMA
RoboT->Write("RoboT",TObject::kOverwrite);
///*****************************************************************
cout << "========================GenFile.Write (MCgen.root)=============" << endl;
GenFile.Write();
cout<<"--------------------------GenFile.ls-------------------------------"<<endl;
GenFile.ls();
//cout<<"--------------------------GenFile.ShowStreamerInfo-------------------------"<<endl;
//GenFile.ShowStreamerInfo();
GenFile.Close();
cout << "========================HistoFile.Write (histo.root)===========" << endl;
HistoFile.Write();
cout<<"--------------------------HistoFile.ls-------------------------------"<<endl;
HistoFile.ls();
cout<<"--------------------------HistoFile.ShowStreamerInfo-------------------------"<<endl;
//HistoFile.GetListOfKeys()->Print();
//HistoFile.Close();
cout<<"%%% ===============End Start.C ================== %%%%"<<endl;
return 0;
}
