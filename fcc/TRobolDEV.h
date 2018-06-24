#ifndef TRobolDEV_H
#define TRobolDEV_H
///////////////////////////////////////////////////////////////////////////////
///                CLASS TRobol2pB
///     Three-stroke engine for analysis: Initialize-Production-Finalize
///     Daugter class specific to simple multidimensional integration

/// OUR headers
#include "TRobol.h"
#include "TMCgenDEV.h"

//[[[[[[[[[[[[[[[[[[[[
extern "C" void kk2f_fort_open_( const int&, const char*, int);
//]]]]]]]]]]]]]]]]]]]]


class TRobolDEV : public TRobol {
/// member functions
  public:
  TRobolDEV();                // explicit default constructor for streamer
  TRobolDEV(const char*);     // user constructor
  virtual ~TRobolDEV();               // explicit destructor
///---
  virtual void Initialize(ofstream*, TFile*, TFile*);
  virtual void Hbooker();
  virtual void Production(double &);
  virtual void Finalize();
/// data members
  double    m_xmin;          // dummy
  double    m_xmax;          // dummy
/// ============== Histograms follow =================================
  TH1D   *hst_weight1;           //! no streamer
////////////////////////////////////////////////////////////////////////////
  ClassDef(TRobolDEV,2); // Monte Carlo generator
};
/////////////////////////////////////////////////////////////////////////////
//                End of the class TRobolDEV                                 //
/////////////////////////////////////////////////////////////////////////////
#endif
