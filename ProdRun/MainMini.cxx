using namespace std;
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

// ROOT headers
#include "TROOT.h"
#include "TFile.h"
#include "TRandom1.h"
#include "TH1.h"

// OUR headers
#include "KKee2f.h"
#include "KKborn.h"

#define SW20 setw(20)<<setprecision(14)
#define SP15 setw(15)<<setprecision(9)

int main()
{
ofstream   OutFile("pro.output",ios::out);  // Logfile output

/////////////////////////////////////////////////////////////
// Random number object
TRandom1 *RN_gen = new TRandom1();       // Central r.n.gen.
long    iniseed = 54217137;
RN_gen->SetSeed(iniseed);

/////////////////////////////////////////////////////////////
// MC generator object
KKee2f *KKMCgen = new KKee2f("MCgen");
//KKMCgen->ls();

/////////////////////////////////////////////////////////////
TH1D *h_NORMA = new TH1D("KKMCgen_NORMA","Normalization histo",10000,0,10000);
KKMCgen->Initialize(  RN_gen, &OutFile, h_NORMA);

cout<<"MainMini: Testing Born distributions"<<endl;
// double KKborn::BornSimple(int KFi, int KFf, double svar, double costhe)
//double KKborn::Born_DizetS(int KFi, int KFf, double svar, double CosThe){
KKborn *Born = KKMCgen->m_BornDist;
double CMSene =189.0;
double svar = CMSene*CMSene;
double CosThe = 0.5;
int KFi=11, KFf=13;
cout<<"MainMini: KFi, KFf, svar, CosThe="<<KFi<<" "<<KFf<<" "<<svar<<" "<< CosThe<<endl;
double Dist0= Born->BornSimple(KFi, KFf, svar, CosThe);
cout<<"MainMini: BornSimple=  "<<SW20<<Dist0<<endl;
double Dist1= Born->Born_DizetS(KFi, KFf, svar, CosThe);
cout<<"MainMini: Born_DizetS= "<<SW20<<Dist1<<endl;

int NevLimPrt =100;
int NevGen =100;
cout<<"MainMini: ********************************** "<<endl;
cout<<"MainMini: type in no. of events: (100?) ";
cin>>NevGen;
cout<<"MainMini: requested "<< NevGen <<" MC events "<<endl;

/////////////////////////////////////////////////////////////
// Small loop over MC events
for(int iev=1; iev<=NevGen; iev++) {
   if( (iev/20000)*20000 == iev) cout<<" iev="<<iev<<endl;
   KKMCgen->Generate();
   if(iev <= NevLimPrt) KKMCgen->m_Event->PrintISR_FSR();
   if(iev <= NevLimPrt) KKMCgen->m_Event->PrintISR_FSR(&OutFile);
}

/////////////////////////////////////////////////////////////
// final printout
KKMCgen->Finalize();

cout << "  |--------------------| "<<endl<<flush;
cout << "  |  TestMini1 Ended   | "<<endl<<flush;
cout << "  |--------------------| "<<endl<<flush;
return 0;
}
