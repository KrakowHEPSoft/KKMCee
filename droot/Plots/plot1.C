{
// .x canvas.cxx
//
gROOT->Reset();
//------------------------------------------------------------------------
// Connect ROOT histogram/ntuple demonstration file
TFile f1("../demoC/rmain.root");

//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Float_t  WidPix, HeiPix;
WidPix = 500;  HeiPix = 500;
TCanvas *cPlot1 = new TCanvas("cPlot1","Plot1",440,90, WidPix,HeiPix);
cPlot1->SetFillColor(10);
//cPlot1->Divide(2, 1, 0.0, 0);
///////////////////////////////////////////////
cPlot1->Draw();
//cPlot1->cd(1);
hst_Q2kloe->DrawCopy("e");
cPlot1->Update();

//------------------------------------------------------------------------
}
