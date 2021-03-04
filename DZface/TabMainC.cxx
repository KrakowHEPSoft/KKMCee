///////////////////////////////////////////////////////////////////////////////
// Auxiliary program for testing writing EW formfactor tables into disk files
// Alternative for TabMain.f
///////////////////////////////////////////////////////////////////////////////
// C++ headers
using namespace std;
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

// ROOT headers
#include "TROOT.h"
#include "TFile.h"

extern "C" {
  void fort_open_( const int&, const char*, int);
  void fort_close_(const int&);
  void bornv_initialize_(double[]);
  void hhdizet_initialize_(double[]);
  void readataz_(const int&, const int&, double[]);
}

void ReaData(const char *DiskFile, int imax, double xpar[]);

int main()
{
//
  cout    << "  |--------------------| "<<endl<<flush;
  cout    << "  | MainProg  Entering | "<<endl<<flush;
  cout    << "  |--------------------| "<<endl<<flush;
//////////////////////////////////////////////////////////////
  int jmax = 10000;
  double m_xpar[10001];
  double m_ypar[10001];

  ReaData("./KKMCee_defaults", jmax, m_xpar);       // f77 indexing in xpar
  ReaData("./pro.input",      -jmax, m_xpar);       // jmax<0 means no-zeroing
  for(int j=0;j<jmax;j++) m_ypar[j]=m_xpar[j+1];    // c++ indexing in ypar

//=============================================================
//   opening disk fime for fortran part of code
  int m_out = 16;
  const char *output_file = "./pro77.output";
  int sl2 = strlen(output_file);
  fort_open_(m_out,output_file,sl2);
//=============================================================
  bornv_initialize_(m_ypar);
//=========================
  hhdizet_initialize_(m_ypar);
//==================================
  cout    <<endl<<flush;
  cout    << "  |--------------------| "<<endl<<flush;
  cout    << "  |  MainProg  Ended   | "<<endl<<flush;
  cout    << "  |--------------------| "<<endl<<flush;
  return 0;
  }

void ReaData(const char *DiskFile, int imax, double xpar[])
//////////////////////////////////////////////////////////////
//    subprogram reading input data file and packing        //
//    entries into matrix xpar                              //
//    WARNING: input file cannot include empty lines        //
//    it cannot handle entries like 1d-30, has to be 1e-30! //
//////////////////////////////////////////////////////////////
{
  char trail[200];
  char ch1;
  int  foundB=0, line, indx;
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
}// ReaData


