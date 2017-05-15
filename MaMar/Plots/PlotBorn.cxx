//////////////////////////////////////////////////////////
//line-shape i Afb, with and without box corrections
//    w funkcji sqrt(s), zakres 70 - 150 GeV
//    fixed costheta = 0.0, 0.3, 0.6, 0.9
//    procesy: ee-> mumu, uubar, ddbar
//////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
using namespace std;

#include <math.h>
#include "TROOT.h"
#include "TCanvas.h"
#include "TF2.h"
#include "TH2.h"
#include "TApplication.h"
#include "TObjString.h"
#include "TFile.h"

#include "HisNorm.h"
#include "KKplot.h"

//
TFile DiskFileB("RhoSemi.root","RECREATE","histograms");
FILE *DFile;


void ReaData(char DiskFile[], int imax, double xpar[])
//////////////////////////////////////////////////////////////
//    subprogram reading input data file and packing        //
//    entries into matrix xpar                              //
//    WARNING: input file cannot include empty lines        //
//    it cannot handle entries like 1d-30, has to be 1e-30! //
//////////////////////////////////////////////////////////////
{
  char trail[200];
  char ch1;
  int  foundB=0, foundE=0, line, indx;
  int  line_max =2000;
  double value;
  cout<<"============================ReaData=============================="<<endl;
  cout<<"===                     "<< DiskFile <<"               =========="<<endl;
  cout<<"================================================================="<<endl;
  ifstream InputFile;
  InputFile.open(DiskFile);
  for(indx=0;indx<imax; indx++) xpar[indx]=0.0;
  for(line=0;line<line_max; line++){
    InputFile.get(ch1);
    if( ch1 == 'B') foundB=1;
    InputFile.getline(trail,200);
    if(foundB) break;
  }
  for(line=0;line<line_max; line++){
    InputFile.get(ch1);
    if( ch1 == 'E'){
      foundE=1;
      InputFile.getline(trail,200);
      cout<<ch1<<trail<<"["<<line<<"]"<<endl;
      break;
    }
    if( ch1 == '*'){
      InputFile.getline(trail,200);
      cout<<ch1<<trail<<endl;
    }else{
      InputFile>>indx>>value;
      if(indx<0 || indx>abs(imax) ){
	cout<<" ++++++++ReaData: wrong indx = "<<indx<<endl;
	exit(0);
      }
      xpar[indx] = value;
      //xpar[indx-1] = value; // correction for fortran indexing in input file
      InputFile.getline(trail,200);
      cout<<ch1;
      cout<<setw(4)<<indx<<setw(15)<<value<<" ";
      cout<<trail<<endl;
    }
  }
  cout<<"================================================================="<<endl;
  InputFile.close();
}

void TabBorn(){
cout<<"==========================================================="<<endl;
cout<<"================ TabBorn BEGIN ============================"<<endl;
//////////////////////////////////////////////////////////////
//   Initialize MC generator and analysis programs          //
//////////////////////////////////////////////////////////////
double   m_xpar[10001];      // complete input of KKMC run
const int jmax =10000;
//ReaData("./KK2f_defaults_ERW",    jmax, m_xpar);  // numbering as in input!!!
ReaData("../../.KK2f_defaults",     jmax, m_xpar);  // numbering as in input!!!
// User input data
//ReaData("../workAFB/workAFB_189GeV.input", -jmax, m_xpar);  // jmax<0 means no-zeroing
double ypar[jmax];
for(int j=0;j<jmax;j++) ypar[j]=m_xpar[j+1];    // ypar has c++ numbering
// KKMC and KKsem initialization
char *output_file = "./kksem.output";
long stl2 = strlen(output_file);
long mout =16;
kk2f_fort_open_(mout,output_file,stl2);
kk2f_initialize_(ypar);
kksem_initialize_(ypar);

//  ************* user histograms  *************
TH1D *hst_weight3 = new TH1D("hst_weight3" ,  "MC weight",      100, -1.0, 2.0);
hst_weight3->Sumw2();
//************************************
  DFile = fopen("TableBorn.txt","w");
//************************************
double svar2, sigma,  CosTheta;
double dSig0, dSig3, dSig6, dSig9;
int m_KFini = 11;
int m_KFf   = 13;
fprintf(DFile," e+e- --> mu+ mu- \n");
fprintf(DFile," d(sigma)/d(cos_theta) [nb], cos_theta = 0.0, 0.3, 0.6, 0.9 \n");
for( double Ene=70; Ene<=150; Ene += 2 ){
   cout<<"Ene="<<Ene<<endl;
   svar2 = Ene*Ene;
   //
   CosTheta=0; bornv_interpogsw_(m_KFf,svar2, CosTheta);
   dSig0 = bornv_dizet_( 1, m_KFini, m_KFf, svar2, CosTheta, 0.0, 0.0, 0.0, 0.0);
   //
   CosTheta=0.3; bornv_interpogsw_(m_KFf,svar2, CosTheta);
   dSig3 = bornv_dizet_( 1, m_KFini, m_KFf, svar2, CosTheta, 0.0, 0.0, 0.0, 0.0);
   //
   CosTheta=0.6; bornv_interpogsw_(m_KFf,svar2, CosTheta);
   dSig6 = bornv_dizet_( 1, m_KFini, m_KFf, svar2, CosTheta, 0.0, 0.0, 0.0, 0.0);
   //
   CosTheta=0.9; bornv_interpogsw_(m_KFf,svar2, CosTheta);
   dSig9 = bornv_dizet_( 1, m_KFini, m_KFf, svar2, CosTheta, 0.0, 0.0, 0.0, 0.0);
   fprintf(DFile,"Ene= %10.5f  dSig_dCos= %12.7f  %12.7f  %12.7f  %12.7f  \n", Ene, dSig0, dSig3, dSig6, dSig9);
 }//j

//************************************
  fclose(DFile);
//************************************

cout<<"================ TabBorn END   ============================"<<endl;

}//TabBorn

///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //++++++++++++++++++++++++++++++++++++++++
  //TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++
  //
  DiskFileB.cd();
  //
  TabBorn();
  //
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();

  //++++++++++++++++++++++++++++++++++++++++
  //theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}

