//#ifndef PartLund_H
//#define PartLund_H
//////////////////////////////////////////////////////////////////////////////
//            Header   CLASS   PartLund                                     //
//                                                                          //
// This is class for a single particle/object in Lund (HepEvt) coding style //
// It can be also string, jet, initial beam etc.                            //
//////////////////////////////////////////////////////////////////////////////

#include "TObject.h"

#include "TLorentzVector.h"

//class PartLund  : public TObject {
class PartLund  {
public:
  PartLund  *fPrevious;   // pointer for constructing lists
  PartLund  *fNext;       // pointer for constructing lists
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
  PartLund() {;}
  PartLund( int lserial,    int kstatus,      int kflavor, 
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
  ~PartLund() {;}
//------ header printing -------
  void PartLund::Print(int mode = 1);
//------ header printing -------
  void PartLund::ListPrint();
//------
//  ClassDef(PartLund,1)   // PartLund  class, <--for dictionary
};
//============================================================================
//                   END OF CLASS PartLund
//============================================================================
//#endif
