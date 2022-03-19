///////////////////////////////////////////////////////////////////////////////
#include "HepFace.h"

ClassImp(HepFace);


HepFace::HepFace()
{
  // This constructor is for ROOT streamers ONLY
  // all pointers has to be NULLed
  cout<< "----> HepFace Default Constructor (for ROOT only) "<<endl;
  m_Out= NULL;
}

///_____________________________________________________________
HepFace::HepFace(ofstream *OutFile)
{
  cout<< "----> HepFace USER Constructor "<<endl;
  m_Out = OutFile;
}//HepFace

///______________________________________________________________________________________
HepFace::~HepFace()
{
  //Explicit destructor
  cout<< "----> HepFace::HepFace !!!! DESTRUCTOR !!!! "<<endl;
}///destructor

double HepFace::sqr( const Double_t x ){ return x*x;};

///______________________________________________________________________________________
void HepFace::Initialize()
{
  cout  << "----> HepFace::Initialize, Entering "<<endl;
//=================================================================
// BX*** macros are in MCdev/BXFORMAT.h
  BXOPE(*m_Out);
  BXTXT(*m_Out,"========================================");
  BXTXT(*m_Out,"======    HepFace::Initialize     ======");
  BXTXT(*m_Out,"========================================");

  m_Hvent = new GenEvent(Units::GEV,Units::MM);

  ///////////////////////////////////////////////////
}// Initialize

///______________________________________________________________________________________
void HepFace::make1()
{
  m_Hvent->clear();

// status of outgoing fermions
  int status=3;                                 // stable mu, neutrina
  if(abs(m_Event->m_KFfin)==15) status=1;       // instable tau
  if(abs(m_Event->m_KFfin) >= 1 && abs(m_Event->m_KFfin)  < 5) status=1; // quarks
//-------------------------------
// in and out fermions
/*
  GenParticlePtr pe1 = std::make_shared<GenParticle>( FourVector(
     m_Event->m_Pf1.Px(),  m_Event->m_Pf1.Py(),m_Event->m_Pf1.Pz(),  m_Event->m_Pf1.E() ), m_Event->m_KFini,3);
*/
  GenParticlePtr pe1 = std::make_shared<GenParticle>( Vect4( m_Event->m_Pf1), m_Event->m_KFini,3);
  GenParticlePtr pe2 = std::make_shared<GenParticle>( Vect4( m_Event->m_Pf2), m_Event->m_KFini,3);
  GenParticlePtr pe3 = std::make_shared<GenParticle>( Vect4( m_Event->m_Qf1), m_Event->m_KFfin, status);
  GenParticlePtr pe4 = std::make_shared<GenParticle>( Vect4( m_Event->m_Qf2), m_Event->m_KFfin, status);
// ISR photons:
  auto ISR_photons=vector<GenParticlePtr>();
  for(int i=0; i< m_Event->m_nPhotISR; i++){
    GenParticlePtr tmp_photon=std::make_shared<GenParticle>( Vect4(m_Event->m_PhotISR[i]), 22, 3);
  }//for i
// FSR photons:
  auto FSR_photons=vector<GenParticlePtr>();
  for(int i=0; i< m_Event->m_nPhotFSR; i++){
    GenParticlePtr tmp_photon=std::make_shared<GenParticle>( Vect4(m_Event->m_PhotFSR[i]), 22, 3);
    FSR_photons.push_back(tmp_photon);
  }// for i
// intermediate Z bozon:
  GenParticlePtr pZ = std::make_shared<GenParticle>( Vect4(m_Event->m_PX), 23, 3);
//

/*
  GenVertexPtr v1 = std::make_shared<GenVertex>();
  v1->add_particle_in (pe1);
  v1->add_particle_out(pe2);
  m_Hvent->add_vertex(v1);
*/

// for(int i=0; i < m_Event->m_nPhotISR; i++){
//    if( m_Event->m_nPhotISR>0){
//    int i =0;
//    GenVertexPtr v_tmp = std::make_shared<GenVertex>();
//    v_tmp->add_particle_in (pe1);
//    v_tmp->add_particle_out(pe1); // here I am cheating, fix this later, talk to Staszek if Evt class has momentum after FS or before
//    for(int i=0; i < m_Event->m_nPhotISR; i++){
//    v_tmp->add_particle_out(ISR_photons[i]);
//    m_Hvent->add_vertex(v_tmp);
//    }

 // version without photons
  GenVertexPtr vZ =  std::make_shared<GenVertex>();
  vZ ->add_particle_in(pe1);
  vZ ->add_particle_in(pe2);
  vZ ->add_particle_out(pZ);
  m_Hvent->add_vertex(vZ);
// now decaying Z:
  GenVertexPtr vZ2 =  std::make_shared<GenVertex>();
  vZ2 ->add_particle_in(pZ);
  vZ2->add_particle_out(pe3);
  vZ2->add_particle_out(pe4);
  m_Hvent->add_vertex(vZ2);

  Print::listing(*m_Hvent);
  Print::content(*m_Hvent);

}//make1
