#ifndef JetAnalyzer_H
#define JetAnalyzer_H
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                CLASS JetAnalyzer                                          //
//                                                                           //
//    Concstructs 4 Jets in double W decay final state                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TROOT.h"
#include "TFile.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TProfile.h"

#include <iostream.h>
#include <fstream.h>

#include "KorEvent.h"
#include "JetAnalyzer.h"

class JetAnalyzer : public TNamed {
public:
//
   TH1F      *hst_pene; // parton minimum energy
   TH1F      *hst_jene; // jet minimum energy
   TNtuple   *jtuple;   // big ntuple with jet masses etc
//
public:
  JetAnalyzer();
  ~JetAnalyzer() {;}
public:
//////////////////////////////////////////////////////////////////////////////
void JetAnalyzer::Book(KorEvent &event,
			long ntot, double wt, double wt2);
////////////////////////////////////////////////////////////////////////////
                      ClassDef(JetAnalyzer,1)
};
////////////////////////////////////////////////////////////////////////////
#endif
