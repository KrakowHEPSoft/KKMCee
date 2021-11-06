/////////////////////////////////////////////////////////////////////////////
//!               Class WtMonit                                            //
/////////////////////////////////////////////////////////////////////////////
#ifndef TWtMon_H
#define TWtMon_H

// C++ headers
using namespace std;
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

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


class TWtMon : public TObject{
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
  TH1D   *m_HST;
 public:
  TWtMon();
  TWtMon(const char*);
  TWtMon(const char*, const char*, int, const double);
  virtual ~TWtMon(){};                         // Destructor
  void Reset();
  void Fill(const double);
  void Fill(const double, const double);
  void GetAver(double&,  double& );
  void GetAll( double&,  double&,  double&,  double&,  double&, double&,
               double&,  double&,  double&,  double&,  double&);
  void PrintAll();
  void PrintLine();
////////////////////////////////////////////////////////////////////////////
  ClassDef(TWtMon,1); // Info class for html documentation
};
/////////////////////////////////////////////////////////////////////////////
#endif
