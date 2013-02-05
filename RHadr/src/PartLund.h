#ifndef PartLund_H
#define PartLund_H
//////////////////////////////////////////////////////////////////////////////
//            Header   CLASS   PartLund                                     //
//                                                                          //
// This is class for single particle/object in Lund coding style            //
// It can be also string, jet, initial beam etc.                            //
//////////////////////////////////////////////////////////////////////////////

#include "TObject.h"

#include "VLorenz.h"

class PartLund  : public TObject {
//class PartLund  {
public:
  PartLund  *previous;   // pointer for constructing lists
  PartLund  *next;       // pointer for constructing lists
  int     m_lserial;       // Lund serial number
  int     m_kstatus;       // status
  int     m_kflavor;       // flavour
  int     m_kparent;       // parent
  int     m_kFirstChild;   // first child
  int     m_kLastChild;    // last  child
  VLorenz m_pmom;          // 4-momentum, pmom[0] is energy [GeV]
  double  m_pmass;         // mass [GeV]
  VLorenz m_Vertex;        // vertex position [mm], V[0] is time in Lab
  double  m_LifeTime;      // lifetime [mm/sec]
public:
  PartLund() {;}
  PartLund( int lserial,    int kstatus,      int kflavor, 
	    int kparent,    int kFirstChild,  int kLastChild,
	    double px,      double py,        double pz, 
	    double en,      double pmass,
	    double Vx,      double Vy,        double Vz, 
            double Vt,      double LifeTime):
                  m_lserial(lserial),
	          m_kstatus(kstatus),
		  m_kflavor(kflavor),
                  m_kparent(kparent),
                  m_kFirstChild(kFirstChild),
                  m_kLastChild(kLastChild),
                  m_pmom(  en, px, py, pz),
                  m_pmass(pmass),
                  m_Vertex(Vt, Vx, Vy, Vz),
                  m_LifeTime(LifeTime)
{};
//------ Destructor------
  ~PartLund() {;}
//------ header printing -------
  void PartLund::print(int mode = 1);
//------ header printing -------
  void PartLund::ListPrint();
  ClassDef(PartLund,1)   // PartLund  class, <--for dictionary
};
//============================================================================
//                   END OF CLASS PartLund
//============================================================================
#endif
