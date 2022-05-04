// make test1-start

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

/////////////////////////////////////////////////////////////
// Small loop over MC events
for(int iev=1; iev<=10000; iev++) {
   if(iev <= 10) cout<<" iev ="<< iev<<endl;
   KKMCgen->Generate();
   if(iev <= 10) KKMCgen->m_Event->PrintISR_FSR();
}// for iev
/////////////////////////////////////////////////////////////
// final printout from KKMCee goes to pro.output
KKMCgen->Finalize();


cout << "  |--------------------| "<<endl<<flush;
cout << "  |  TestMini1 Ended   | "<<endl<<flush;
cout << "  |--------------------| "<<endl<<flush;
return 0;
}
