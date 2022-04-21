#include "TRandom1.h"
#include "TH1.h"
#include "KKee2f.h"
/////////////////////////////////////////////////////////////
int main(){
ofstream   OutFile("pro.output",ios::out);  // Logfile output
TRandom1 *RN_gen = new TRandom1();// Central random numb. gen.
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
   if( (iev/20000)*20000 == iev) cout<<" iev="<<iev<<endl;
   KKMCgen->Generate();
}// for iev
/////////////////////////////////////////////////////////////
KKMCgen->Finalize(); // final printout
cout << "  |--------------------| "<<endl<<flush;
cout << "  |  TestMini1 Ended   | "<<endl<<flush;
cout << "  |--------------------| "<<endl<<flush;
return 0;
}// main
