/////////////////////////////////////////////////////////////////////////////
//!               Class WtMonit                                            //
/////////////////////////////////////////////////////////////////////////////
#ifndef THwtMon_H
#define THwtMon_H

// C++ headers
using namespace std;
#include<stdlib.h>
#include <iostream>
#include<fstream>
#include<iomanip>
#include<math.h>

// ROOT headers
#include "TROOT.h"
#include "TH1.h"
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                Class  WtMonit                                             //
//                                                                           //
//      Small auxiliary class for detailed monitoring single MC weight.      //
//      It creates and uses 1 histograms.                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


class THwtMon: public TH1D {
 private:
  double  m_Ntot;
  double  m_Nacc;
  double  m_Nzer;
  double  m_Wtmax;
  double  m_Wtmin;
  double  m_Wtsum;
  double  m_Wtsu2;
  double  m_WtsuOv;
  double  m_Xup;
 public:
  THwtMon();
  THwtMon(const char*);
  THwtMon(const char*, const char*, int, const double);
  ~THwtMon(){};                         // Destructor
  void Reset();
  int Fill(const double);
  int Fill(const double, const double);
  void GetAver(double&,  double& );
  void GetAll( double&,  double&,  double&,  double&,  double&, double&,
               double&,  double&,  double&,  double&,  double&);
  void PrintAll();
  void PrintLine();
////////////////////////////////////////////////////////////////////////////
  ClassDef(THwtMon,2); // Info class for html documentation
};
/////////////////////////////////////////////////////////////////////////////
#endif
