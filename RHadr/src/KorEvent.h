#ifndef KorEvent_H
#define KorEvent_H
//////////////////////////////////////////////////////////////////////////////
//                     CLASS   KorEvent                                     //
//                                                                          //
// This is class of KORALW events                                           //
// after HADRONIZATION encoded in Lund/PDG coding style                     //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#include "TObject.h"

#include "PartLund.h"

class KorEvent : public TObject {
//class KorEvent {
// partonic level variables
  public:
  int          m_nphot;           // number of photons
  VLorenz      m_photmom[100];    // list of photon momenta
  VLorenz      m_wminus;          // 4-momentum of W-
  VLorenz      m_wplus;           // 4-momentum of W+
  VLorenz      m_ferm1;           // 4-momentum of fermion 1
  VLorenz      m_ferm2;           // 4-momentum of fermion 2
  VLorenz      m_ferm3;           // 4-momentum of fermion 3
  VLorenz      m_ferm4;           // 4-momentum of fermion 4
public:
// Lund/PDG style record, 
// to be moved to mother class from which KorEvent will inherit
   int         m_npart;           // total number of particles in part[]
PartLund       m_part[20000];     // all particles in Lund-like format
public:
   int         m_njet;            // total number of jets
   PartLund    m_jet[200];        // jets in Lund-like format
public:
//------ header of constructor/destructor -------
  KorEvent();
  ~KorEvent();
//------ printing -------
  void  Print( long level);
//------------------------------
void   GetPartonMass( double     m[6] );
double GetPartonAngle(double angle[6] );
void   GetJetMass(    double     m[6] );
double GetJetAngle(   double angle[6] );
  ClassDef(KorEvent,1)   // KorEvent  class, <--for dictionary
};
//////////////////////////////////////////////////////////////////////////////
//                   END OF CLASS KorEvent                                  //
//////////////////////////////////////////////////////////////////////////////
#endif
