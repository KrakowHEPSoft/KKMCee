#include "KKHepMCEvent.h"

KKHepMCEvent::KKHepMCEvent(HepMC3::GenEvent &e, int size)
{
  evt=&e;
  particle_count=size;
  particles = new KKHepMCparticle*[particle_count];
  for(int i=0; i<m_particle_count; ++i) {
    particles[i] = new KKHepMCparticle(*e.particles()[i],this,i+1); // talk to Staszek about this, in KKMC event does not have particles inside?!
  }

int KKHepMCEvent::GetNumOfParticles() {
    return m_particle_count;
}

int KKHepMCEvent::GetEventNumber() {
    return evt->event_number();
}
HEPParticle* KKHepMCEvent::GetParticle(int idx) {
    if(idx < 1 || idx > GetNumOfParticles()) {
        std::cout << "Warning can not get particle "<< idx;
        std::cout <<", particle ID not valid" << std::endl;
        return 0;
    }
    return particles[idx-1]; //Particle ID starts at 1
}
KKHepMCparticle* KKHepMCEvent::GetParticleWithId( int id ) {
    for(int i=0; i <  GetNumOfParticles(); i++) {
        if(particles[i]->part->id()==id)
            return particles[i];
    }
    std::cout << "Could not find particle with id "<<id<<std::endl;
    return 0; //and have some error about not finding the
    //particle
}

std::vector<double> * KKHepMCEvent::Sum4Momentum() {
    std::vector<double> * sum = new std::vector<double>(4,0.0);
    for(int i=0; i < GetNumOfParticles(); i++) {
        if(particles[i]->IsStable()) {
            sum->at(0)+=particles[i]->GetPx();
            sum->at(1)+=particles[i]->GetPy();
            sum->at(2)+=particles[i]->GetPz();
            sum->at(3)+=particles[i]->GetE();
        }
    }
    return sum;
}

