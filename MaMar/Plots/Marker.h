/////////////////////////////////////////////////////////////////////
void BlueBox(TH1D *Hst, Float_t msize){
  // define plotting style (ISR)
  Hst->SetMarkerStyle(21);//box
  Hst->SetMarkerColor(4); //blue
  Hst->SetMarkerSize(msize);
  Hst->SetLineColor(4);   // blue
}
/////////////////////////////////////////////////////////////////////
void RedTriangle(TH1D *Hst, Float_t msize){
  // define plotting style (ISR + Beamstrahl)
  Hst->SetMarkerStyle(22); //triangle
  Hst->SetMarkerColor(2); //blue
  Hst->SetMarkerSize(msize);
  Hst->SetLineColor(2);   // red
}
/////////////////////////////////////////////////////////////////////
void BlackBullet(TH1D *Hst, Float_t msize){
  // define plotting style (Born)
  Hst->SetMarkerStyle(8); //bullet
  Hst->SetMarkerColor(1); //black
  Hst->SetMarkerSize(msize);
  Hst->SetLineColor(1);    //black
}

/////////////////////////////////////////////////////////////////////
void MarkPad3(TPad *Lpad){
//  description of markers as a separate pad within the pad
  Lpad->SetFillStyle(4000);
  Lpad->Draw();
  Lpad->cd();
  TLatex *tex2 = new TLatex();
  tex2->SetTextAlign(21);
  tex2->SetTextSize(0.15);
  tex2->DrawLatex(0.50,0.85,"(#Delta d#sigma)/ d#sigma");
  //
  tex2->SetTextSize(0.12);
  tex2->SetTextAlign(12);
  tex2->DrawLatex(0.20,0.65,"NLL2, LL3");
  tex2->DrawLatex(0.20,0.40,"LL2,  NLL2");
  tex2->DrawLatex(0.20,0.15,"LL3");
  //
  TMarker *bullet = new TMarker(0.10,0.65, 8);
  bullet->SetMarkerSize(1.0);
  bullet->SetMarkerColor(1);
  bullet->Draw();
  //
  TMarker *box    = new TMarker(0.10,0.14, 21);
  box->SetMarkerSize(1.0);
  box->SetMarkerColor(4);
  box->Draw();
  //
  TMarker *triang = new TMarker(0.10,0.40, 22);
  triang->SetMarkerSize(1.3);
  triang->SetMarkerColor(2);
  triang->Draw();
}

/////////////////////////////////////////////////////////////////////
void MarkPad2(TPad *Lpad){
//  description of markers as a separate pad within the pad
  Lpad->SetFillStyle(4000);
  Lpad->Draw();
  Lpad->cd();
  TLatex *tex2 = new TLatex();
  tex2->SetTextAlign(21);
  tex2->SetTextSize(0.15);
  tex2->DrawLatex(0.50,0.85,"Q^{2}#times d#sigma/dQ^{2} [nb]");
  //
  tex2->SetTextSize(0.15);
  tex2->SetTextAlign(12);
  tex2->DrawLatex(0.20,0.65,"E_{#gamma}>10MeV");
  tex2->DrawLatex(0.20,0.40,"and #theta_{#gamma}");
  tex2->DrawLatex(0.20,0.15,"and #theta_{#pi}, P_{#pi}^{T}");
  //
  TMarker *bullet = new TMarker(0.10,0.65, 8);
  bullet->SetMarkerSize(1.0);
  bullet->SetMarkerColor(1);
  bullet->Draw();
  //
  TMarker *box    = new TMarker(0.10,0.14, 21);
  box->SetMarkerSize(1.0);
  box->SetMarkerColor(4);
  box->Draw();
  //
  TMarker *triang = new TMarker(0.10,0.40, 22);
  triang->SetMarkerSize(1.3);
  triang->SetMarkerColor(2);
  triang->Draw();
}

/////////////////////////////////////////////////////////////////////
void MarkPad1(TPad *Lpad){
//  description of markers as a separate pad within the pad
  Lpad->SetFillStyle(4000);
  Lpad->Draw();
  Lpad->cd();
  TLatex *tex2 = new TLatex();
  tex2->SetTextAlign(21);
  tex2->SetTextSize(0.15);
  tex2->DrawLatex(0.50,0.85,"Event selection");
  //
  tex2->SetTextSize(0.15);
  tex2->SetTextAlign(12);
  tex2->DrawLatex(0.20,0.65,"E_{#gamma}>10MeV");
  tex2->DrawLatex(0.20,0.40,"and #theta_{#gamma}");
  tex2->DrawLatex(0.20,0.15,"and #theta_{#pi}, P_{#pi}^{T}");
  //
  TMarker *bullet = new TMarker(0.10,0.65, 8);
  bullet->SetMarkerSize(1.0);
  bullet->SetMarkerColor(1);
  bullet->Draw();
  //
  TMarker *box    = new TMarker(0.10,0.14, 21);
  box->SetMarkerSize(1.0);
  box->SetMarkerColor(4);
  box->Draw();
  //
  TMarker *triang = new TMarker(0.10,0.40, 22);
  triang->SetMarkerSize(1.3);
  triang->SetMarkerColor(2);
  triang->Draw();
}
