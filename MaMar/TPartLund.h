#ifndef TPartLund_H
#define TPartLund_H
//////////////////////////////////////////////////////////////////////////////
//            Header   CLASS   TPartLund                                    //
//                                                                          //
// This is class for a single particle/object in Lund (HepEvt) coding style //
// It can be also string, jet, initial beam etc.                            //
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "TObject.h"

#include "TLorentzVector.h"

//class TPartLund  : public TObject {
class TPartLund  {
public:
  TPartLund  *fPrevious;   // pointer for constructing lists
  TPartLund  *fNext;       // pointer for constructing lists
  int     fSerial;       // Lund serial number
  int     fStatus;       // status
  int     fFlafor;       // flavour
  int     fParent;       // parent
  int     fFirstChild;   // first child
  int     fLastChild;    // last  child
  TLorentzVector fMom;   // 4-momentum, pmom[3] is energy [GeV]
  double  fMass;         // mass [GeV]
  TLorentzVector fVertex; // vertex position [mm], V[0] is time in Lab
  double  fLifeTime;      // lifetime [mm/sec]
public:
  TPartLund() {;}
  TPartLund( int lserial,    int kstatus,      int kflavor,
            int kparent,    int kFirstChild,  int kLastChild,
            double px,      double py,        double pz, 
            double en,      double pmass,
            double Vx,      double Vy,        double Vz, 
            double Vt,      double LifeTime):
            fSerial(lserial),
            fStatus(kstatus),
            fFlafor(kflavor),
            fParent(kparent),
            fFirstChild(kFirstChild),
            fLastChild(kLastChild),
            fMom(  px, py, pz, en),
            fMass(pmass),
            fVertex( Vx, Vy, Vz, Vt),
            fLifeTime(LifeTime)
{};
//------ Destructor------
  ~TPartLund() {;}
//------ header printing -------
  void Print(int mode = 1);
//------ header printing -------
  void ListPrint();
//------
  ClassDef(TPartLund,1)   // TPartLund  class, <--for dictionary
};
//============================================================================
//                   END OF CLASS TPartLund
//============================================================================
#endif
