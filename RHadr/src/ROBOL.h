#ifndef ROBOL_H
#define ROBOL_H
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                CLASS ROBOL                                                //
//                                                                           //
//    Factory of factories                                                   //
//    It contains Makers for MC generation and Analysis event per event      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include "TROOT.h"

#include "KorEvent.h"
#include "KoralwMaker.h"
#include "BEwtMaker.h"
#include "JetAnalyzer.h"

class ROBOL : public TNamed {
public:
  ///////////////////////////////
  //   Three basic Makers      //
  KoralwMaker   KoralW;    // define MC generator maker (manager)
  BEwtMaker     BEfactory; // define BEwt maker (calculator)
  JetAnalyzer   J4factory; // define 4jet maker (analyzer)
  /////////////////////////////////
  //  Current event     (static) //
  KorEvent      current_event;
  /////////////////////////////////
  //  Control Histograms         //
  TH1F       *hst_BEwt;    // BE weight before rejection
  TH1F       *hst_BEwtAR;  // BE weight After  Rejection
  /////////////////////////////////
  // rejection control variables //
  long       KeyRej;       // rejection key (=0 for rejection ON)
  double     WtMax;        // maximum rejection weight
public:
  ROBOL() {
    hst_BEwt  = new TH1F("hst_BEwt" ,  "BE weight before rejection",50,0,10);
    hst_BEwtAR= new TH1F("hst_BEwtAR" ,"BE weight after  rejection",50,0,10);
  }
  ~ROBOL() {;}
public:
// Methods
void ROBOL::Initialize(long &NevTot);
void ROBOL::Production(long &iEvent);
void ROBOL::Finalize();
////////////////////////////////////////////////////////////////////////////
                      ClassDef(ROBOL,1)
};
////////////////////////////////////////////////////////////////////////////
#endif
