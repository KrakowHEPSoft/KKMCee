#include "KKHepMCparticle.h"
#include "KKHepMCevent.h"

KKHepMCparticle::KKHepMCparticle()
{};
KKHepMCparticle::KKHepMCparticle(HepMC3::GenParticle& particle, KKevent * e, int Id)
{
   part = &particle;
   SetEvent(e);
   SetId(Id);
};
KKHepMCparticle::~KKHepMCparticle()
{};

KKHepMCparticle:::GetEvent() {
    return event;
}
int const KKHepMCparticle::GetId(){
    return id;
}

double const KKHepMCparticle::GetE(){
  return part->P[0]; // check with staszek the convention
}
double const KKHepMCparticle::GetPx(){
    return part->P[1]; // check with staszek the convention
}
double const KKHepMCparticle::GetPy(){
    return part->P[2]; // check with staszek the convention
}
double const KKHepMCparticle::GetPz(){
    return part->P[3]; // check with staszek the convention
}
double const KKHepMCparticle::GetM(){
    return part->M; 
}
int const KKHepMCparticle::GetPDGId(){ 
   return part->pdg_id(); 
}
int const KKHepMCparticle::GetStatus(){ 
    return part->status();
}
int const KKHepMCparticle::IsStable(){ 
    return (GetStatus() == 1 || !part->end_vertex());
}
int const KKHepMCparticle::Decays(){
    return (!IsHistoryEntry() && !IsStable());
}

int const KKHepMCparticle::IsHistoryEntry(){
    return (GetStatus() == 3);
}
double const KKHepMCparticle::GetVx(){
  return part->production_vertex[0];
}
double const KKHepMCparticle::GetVy(){
  return part->production_vertex[1];
}
double const KKHepMCparticle::GetVz(){
  return part->production_vertex[2];
}
void KKHepMCparticle::SetEvent(HEPEvent * event) {
    this->event=(HepMC3Event*)event;
}
void KKHepMCparticle::SetId(int id) {
    this->id = id;
}
void KKHepMCparticle::SetE(double E) {
    HepMC3::FourVector temp_mom(part->momentum());
    temp_mom.setE(E);
    part->set_momentum(temp_mom);
}
void KKHepMCparticle::SetPx(double px) {
    HepMC3::FourVector temp_mom(part->momentum());
    temp_mom.setPx(px);
    part->set_momentum(temp_mom);
}
void KKHepMCparticle::SetPy( double py ) {
    HepMC3::FourVector temp_mom(part->momentum());
    temp_mom.setPy(py);
    part->set_momentum(temp_mom);
}
void KKHepMCparticle::SetPz( double pz ) {
    HepMC3::FourVector temp_mom(part->momentum());
    temp_mom.setPz(pz);
    part->set_momentum(temp_mom);
}
void KKHepMCparticle::SetPDGId ( int pdg ) {
    part->set_pid( pdg );
}
void KKHepMCparticle::SetStatus( int st) {
    part->set_status( st );
}
