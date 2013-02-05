#ifndef BEwtMaker_H
#define BEwtMaker_H
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                CLASS BEwtMaker                                            //
//                                                                           //
//    calculation of Bose-Einstein weight for a given event from Koralw      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TROOT.h"
#include "TObject.h"
#include "TFile.h"
#include "TProfile.h"
#include "TNtuple.h"

#include <iostream.h>
#include <fstream.h>

#include "KorEvent.h"
#include "partit1.h"


class BEwtMaker  : public TNamed {
public:
  ofstream BEoutfile;    //! output file for information printouts
  long     BEievent;     // Serial number of event (for printouts)
  long     BEprintlevel; // level of information printout
  long     BEprintlast;  // number of events printed
  double   BErange;    // Range, maximum sqrt(Q2) range for a single pair [GeV]
  long     BEFuncType; // Function type =1,2 for Gausian,exponential
  double   BEpp;       // p parameter
  double   BEradius;   // R radius in fermi units
  double   BElambda;   // lambda parameter for wt  (for wt renormalization)
  double   BEavewt;    // average weight   for wt  (for wt renormalization)
  double   BElambda2;  // lambda parameter for wt2 (for wt2 renormalization)
  double   BEavewt2;   // average weight   for wt2 (for wt2 renormalization)
	               // Histograms //
  TH1F    *hst_nclu;   // histo cluster multiplicity
  TH1F    *hst_unQ2;   // histo Q2 of like-sign-pion-pair
  TH1F    *hst_wtQ2;   // histo Q2 of like-sign-pion-pair, wt
  TH1F    *hst_wt2Q2;  // histo Q2 of like-sign-pion-pair, wt2
  TH1F    *hst_unQ3;   // histo Q2 of like-sign-pion-triplet
  TH1F    *hst_wtQ3;   // histo Q2 of like-sign-pion-triplet, wt
  TH1F    *hst_wt2Q3;  // histo Q2 of like-sign-pion-triplet, wt2
		       // Ntuples //
  double   ctuple2_counter; // counter of entries in ntuple
  double   ctuple2_max;     // maximum entries in ntuple
  TNtuple *ctuple2;         // ntuple Q2,wt,pn
  double   ctuple3_counter; // counter of entries in ntuple
  double   ctuple3_max;     // maximum entries in ntuple
  TNtuple *ctuple3;         // ntuple Q3,wt,pn
public:
  BEwtMaker();
  ~BEwtMaker() {;}
private:
  BEwtMaker(const BEwtMaker &org) { }
public:
SetModel(double range, long FuncType, double pp, double radius);
SetRenorm(double lambda,  double avewt, 
	  double lambda2, double avewt2);
private:
  double CorFun1(double Q2, int mult, double Rf, double pp);
  double CorFun2(double Q2, int mult, double Rf, double pp);
  void   WTcluster(int mult, PartLund *first, double &wt);
  void   BEchain2(KorEvent &event, int kfP, int &nP, double &wt);
/////////////  methods
public:
void   MakeWeight(KorEvent &event, 
			   long &ntot, double &wtot, double &wtot2);
void   BookLSP(KorEvent &event, 
			int kfP, int ntot, double wt, double wt2);
double Q2pair(PartLund *first, PartLund *second);
////////////////////////////////////////////////////////////////////////////
  ClassDef(BEwtMaker,1)
};
////////////////////////////////////////////////////////////////////////////
#endif
