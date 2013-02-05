//=============================================================================
TFile DiskFile("rmain.root","RECREATE","histograms");
//=============================================================================
//  histos  histos  histos  histos  histos  histos  histos  histos  histos
TH1F *hst_unQ2  = new TH1F("hst_unQ2" ,"Q2 2-pi"     ,100,0,1);
TH1F *hst_wtQ2  = new TH1F("hst_wtQ2" ,"Q2 2-pi wted",100,0,1);
TH1F *hst_wt2Q2 = new TH1F("hst_wt2Q2","Q2 2-pi wted",100,0,1);
TH1F *hst_unQ3  = new TH1F("hst_unQ3" ,"Q2 3-pi"     ,100,0,1);
TH1F *hst_wtQ3  = new TH1F("hst_wtQ3" ,"Q2 3-pi wted",100,0,1);
TH1F *hst_wt2Q3 = new TH1F("hst_wt2Q3","Q2 3-pi wted",100,0,1);
TH1F *hst_nclu  = new TH1F("hst_nclu","cluster dist",20,0,20);
//
TH1F *hst_pene = new TH1F("hst_pene" ,"energy parton",100,0,100);
TH1F *hst_jene = new TH1F("hst_jene" ,"energy jet   ",100,0,100);
//
//-----------------------------------------------------------------------------
// ntuples ntuples ntuples ntuples ntuples ntuples ntuples ntuples
//-----------------------------------------------------------------------------
TNtuple *jtuple = new TNtuple("jtuple","all",
	   "m0:m1:m2:m3:m4:m5:x0:x1:x2:x3:x4:x5:np:nj:wt:wt2:pan:jan:pen:jen");
//
TNtuple *ctuple2 = new TNtuple("ctuple2","Q2 ","Q:np:wt");
int   ctuple2_counter=0, ctuple2_max=10000;  // maximum entries in ntuple
TNtuple *ctuple3 = new TNtuple("ctuple3","Q3 ","Q:np:wt");
int   ctuple3_counter=0, ctuple3_max=10000;  // maximum entries in ntuple
//=============================================================================
//=============================================================================
