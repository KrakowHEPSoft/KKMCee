#include "TRandom3.h"
#include "TH1.h"
#include "KKee2f.h"
/////////////////////////////////////////////////////////////
int main(){
ofstream   OutFile("pro.output",ios::out);  // Logfile output
TFile *MCgenFile = new TFile("mcgen.root","RECREATE","Generators");
MCgenFile->cd();
TRandom3 *RN_gen = new TRandom3();// Central random numb. gen.
long iniseed = 54217137; RN_gen->SetSeed(iniseed);
/////////////////////////////////////////////////////////////
KKee2f *KKMCgen = new KKee2f("MCgen"); // MC generator object
/////////////////////////////////////////////////////////////
int nb = 10000;
TH1D *h_NORMA= new TH1D("KKMCgen_NORMA","Normaliz. histo",nb,0,nb);
KKMCgen->Initialize( RN_gen, &OutFile, h_NORMA);
/////////////////////////////////////////////////////////////
int NevGen =100;
cout<<"MainMini: ********************************** "<<endl;
cout<<"MainMini: type in no. of MC events: (100?) ";
cin>>NevGen; cout<<" requested "<< NevGen <<" events"<<endl;
/////////////////////////////////////////////////////////////
// Loop over MC events
for(int iev=1; iev<=NevGen; iev++) {
   KKMCgen->Generate();
   if( (iev/1000)*1000 == iev || iev<=10) {
       cout<<"MainMini: iev="<<iev<<endl;
       KKMCgen->m_Event->EventPrintAll();
       KKMCgen->m_Event->EventPrintAll(&OutFile);
   }//if iev
}// for iev
/////////////////////////////////////////////////////////////
KKMCgen->Finalize(); // final printout
cout << "  |--------------------| "<<endl<<flush;
cout << "  |  MainMini  Ended   | "<<endl<<flush;
cout << "  |--------------------| "<<endl<<flush;
//
KKMCgen->Write("MCgen",TObject::kOverwrite);
RN_gen->Write("RN_gen",TObject::kOverwrite);
MCgenFile->Write();
MCgenFile->Close();
return 0;
}// main
