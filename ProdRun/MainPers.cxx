#include "TRandom3.h"
#include "KKee2f.h"
TFile MCgenFile("mcgen.root","UPDATE","Generators");
/////////////////////////////////////////////////////////////
int main(){
ofstream   OutFile("pro2.output",ios::out);  // Logfile output
cout << "  |--------------------| "<<endl<<flush;
cout << "  |  MainPers Begin    | "<<endl<<flush;
cout << "  |--------------------| "<<endl<<flush;

cout<<"------------------------MCgenFile.GetListOfKeys-----------------------------"<<endl;
MCgenFile.GetListOfKeys()->Print();
cout<<"-------------------------MCgenFile.ShowStreamerInfo-------------------------"<<endl;
MCgenFile.ShowStreamerInfo();
cout<<"-------------------------MCgenFile.ls---------------------------------------"<<endl;
MCgenFile.ls();

// Normalization histogram kept outside MC generator (because of farming)
TH1D *h_NORMA = (TH1D*)MCgenFile.Get("KKMCgen_NORMA");
cout<<"  h_NORMA[1]="  <<h_NORMA->GetBinContent(1)<<endl;
// Central r.n. generator
TRandom3 *RN_gen= (TRandom3*)MCgenFile.Get("RN_gen");
cout<<"  RN_gen="  <<RN_gen->Rndm()<<endl;

//
KKee2f *KKMCgen = (KKee2f*)MCgenFile.Get("MCgen");
cout<<"-------------------------KKMCgen.ls---------------------------------------"<<endl;
KKMCgen->ls();

KKMCgen->Redress(  RN_gen, &OutFile, h_NORMA);

/////////////////////////////////////////////////////////////
int NevGen =100;
cout<<"MainPers: ********************************** "<<endl;
cout<<"MainPers: type in no. of MC events: (100?) ";
cin>>NevGen; cout<<" requested "<< NevGen <<" events"<<endl;
/////////////////////////////////////////////////////////////
// Loop over MC events
for(int iev=1; iev<=NevGen; iev++) {
  KKMCgen->Generate();
  if( (iev/1000)*1000 == iev || iev<=10) {
      cout<<"MainPers: iev="<<iev<<endl;
      KKMCgen->m_Event->EventPrintAll();
      KKMCgen->m_Event->EventPrintAll(&OutFile);
  }//if iev
}//for iev
KKMCgen->Finalize();
cout << "  |--------------------| "<<endl<<flush;
cout << "  |  MainPers Ended    | "<<endl<<flush;
cout << "  |--------------------| "<<endl<<flush;
return 0;

}
