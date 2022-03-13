#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Print.h"



#include "HEPMC3.h"

using namespace HepMC3;



HEPMC3::HEPMC3(HEPEvent *evt)
{
  GenEvent evt(Units::GEV,Units::MM);
  // incoming particles:
  
  GenParticlePtr pe1 = std::make_shared<GenParticle>( FourVector( evt->m_Pf1->Px(),  evt->m_Pf1->Py(),
								  evt->m_Pf1->Pz(),  evt->m_Pf1->E() ), evt->m_KFini,
						      3);
  GenParticlePtr pe2 = std::make_shared<GenParticle>( FourVector( evt->m_Pf1->Px(),  evt->m_Pf1->Py(),
								  evt->m_Pf1->Pz(),  evt->m_Pf1->E() ), evt->m_KFini,
						      3);
  // discus with staszek which is positon wich is electon
  // ISR photons:
  auto ISR_photos=vector<GenParticlePtr>();
  for(int i=0; i < evt->m_nPhotISR, ++i)
    {
       GenParticlePtr tmp_photon=std::make_shared<GenParticle>( FourVector( evt->m_PhotISR[i]->Px(),  evt->m_PhotISR[i]->Py(),
									    evt->m_m_PhotISR[i]->Pz(),  evt->m_PhotISR[i]->E() ), 22,
								3);
       ISR_photos.push_back(tmp_photon);
    }
  // FSR photons:
  auto FSR_photos=vector<GenParticlePtr>();
  for(int i=0; i < evt->m_nPhotFSR, ++i)
    {
      GenParticlePtr tmp_photon=std::make_shared<GenParticle>( FourVector( evt->m_PhotFSR[i]->Px(),  evt->m_PhotFSR[i]->Py(),
                                                                          evt->m_m_PhotFSR[i]->Pz(),  evt->m_PhotFSR[i]->E() ), 22,
                                                              3);
      FSR_photos.push_back(tmp_photon);
    }
  // outgoing particles
  int status=1;
  if(abs(evt->m_KFfni) >= 11 && abs(evt->m_KFfni)  < 17)
    {
      if(abs(evt->m_KFfni)==15) continue;
      status=3;
			       
    }

      
  GenParticlePtr pe3 = std::make_shared<GenParticle>( FourVector( evt->m_Qf1->Px(),  evt->m_Qf1->Py(),
                                                                  evt->m_Qf1->Pz(),  evt->m_Qf1->E() ), evt->m_KFfni,
                                                      status);
  GenParticlePtr pe4 = std::make_shared<GenParticle>( FourVector( evt->m_Qf1->Px(),  evt->m_Qf1->Py(),
                                                                  evt->m_Qf1->Pz(),  evt->m_Qf1->E() ), evt->m_KFfni,
                                                      status);



  // here again we have a problem, we dont know which photon to wich lepton. Talk to Staszek
 for(int i=0; i < evt->m_nPhotISR, ++i)
    {
       GenVertexPtr v_tmp = std::make_shared<GenVertex>();
       v_tmp->add_particle_in (pe1);
       v_tmp->add_particle_out(pe1); // here I am cheating, fix this later, talk to Staszek if Evt class has momentum after FS or before
       v_tmp->add_particle_out(ISR_photos[i]);
       evt.add_vertex(v_tmp);
    }
 // now the Z bozon:
 GenParticlePtr pZ = std::make_shared<GenParticle>( FourVector( evt->m_PX->Px(),  evt->m_PX->Py(),
                                                                  evt->m_PX->Pz(),  evt->m_PX->E() ), 23,
                                                      3);
 GenVertexPtr vZ =  std::make_shared<GenVertex>();
 vZ ->add_particle_in(pe1);
 vZ ->add_particle_in(pe2);
 vZ ->add_particle_out(pZ);
 evt.add_vertex(vZ);
 // now decayign Z:
 GenVertexPtr vZ2 =  std::make_shared<GenVertex>();
 vZ2 ->add_particle_in(pZ);  
 vZ2->add_particle_out(pe3);
 vZ2->add_particle_out(pe4);
 evt.add_vertex(vZ2);

 // here again we have a problem, we dont know which photon to wich lepton. Talk to Staszek                                                    
 for(int i=0; i < evt->m_nPhotFSR, ++i)
    {
       GenVertexPtr v_tmp = std::make_shared<GenVertex>();
       v_tmp->add_particle_in (pe3);
       v_tmp->add_particle_out(pe3); // here I am cheating, fix this later, talk to Staszek if Evt class has momentum after FS or before
       v_tmp->add_particle_out(FSR_photos[i]);
       evt.add_vertex(v_tmp);
    }

};


 






 

 
}
