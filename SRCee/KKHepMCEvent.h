#ifndef KKHepMCEvent_H
#define KKHepMCEvent_H

#include "KKHepMCparticle.h"
#include "KKHepMCEvent.h"


class KKHepMCEvent: public HEPEvent
{
 private:
  KKHepMCparticle **particles;
  HepMC3::GenEvent *evt;
  int particle_count;
  
 public:
  KKHepMCparticle(HepMC3::GenEvent &e, int size);
  int GetNumOfParticles();
  int GetEventNumber();
  HEPParticle* GetParticle(int idx);
  KKHepMCparticle* GetParticleWithId(int id);
  std::vector<double> * Sum4Momentum();
  

  
  
};
#endif
