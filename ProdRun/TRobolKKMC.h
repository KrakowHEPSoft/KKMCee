#ifndef TRobolKKMC_H
#define TRobolKKMC_H
//                CLASS ROBOL                                                //
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"

/// OUR headers
#include "TRobol.h"
#include "KKee2f.h"

class TRobolKKMC : public TRobol
{
///--- data members
 public:
 long   m_NevGen;             // event serial number
 long   m_count1;             // auxiliary event counter (debug)
//
 TH1D   *hst_weight;          //!  No streamer!!!
///////////////////////////////////////////
/// mandatory constructors and destructors
 public:
 TRobolKKMC();                // explicit default constructor for streamer
 TRobolKKMC(const char*);     // user constructor
 virtual ~TRobolKKMC();       // explicit destructor
/// mandatory methods
 virtual void Initialize(ofstream*, TFile*, TFile*);
 virtual void Hbooker();
 virtual void Production(double &);
 void Finalize();
////////////////////////////////////////////////////////////////////////////
                        ClassDef(TRobolKKMC,1)
};
////////////////////////////////////////////////////////////////////////////
#endif
