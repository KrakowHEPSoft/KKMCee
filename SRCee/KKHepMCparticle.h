#ifndef HEPMC_Particle
#define HEPMC_Particle


#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "KKpart.h"

class KKHepMCeven;


class KKHepMCparticle : public KKevent
{
  private:
   int id; // particle unique number
   HepMC::GenParticle *part
   KKHepMCeven * event;
  public:
  
  
  KKHepMCparticle(HepMC3::GenParticle& particle, KKevent * e, int Id);
  ~KKHepMCparticle() // destructor

  HEPEvent* GetEvent();
  /** returns the ID number of particle as used by MC-TESTER (not
      the same as GenParticle pdg_id or barcode).*/

  double const     GetE ()        ;
  double const     GetPx()        ;
  double const     GetPy()        ;
  double const     GetPz()        ; 
  double const     GetM ()        ;
  int    const     GetPDGId ()    ;
  int    const     GetStatus()    ;
  int    const     IsStable()     ;
  int    const     Decays();
  int    const     IsHistoryEntry();
  // vertex locations:
  double const     GetVx  ()      ;
  double const     GetVy  ()      ;
  double const     GetVz  ()      ;

  void   SetEvent(KKevent  *event); // seeting event that particle Belongs to
  void   SetId(int id);
  void   SetE(double E);
  void   SetPx(double px);
  void   SetPy(double py);
  void   SetPz(double pz);
  void   SetM(double m);
  void   SetPDGId(int pdg);
  void   SetStatus(int st);



};
#endif
