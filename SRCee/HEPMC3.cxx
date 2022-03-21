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

  FourVector pe1_v4(evt->m_Pf1->Px(),  evt->m_Pf1->Py(),  evt->m_Pf1->Pz(),  evt->m_Pf1->E() );
  FourVector pe2_v4(evt->m_Pf2->Px(),  evt->m_Pf2->Py(),  evt->m_Pf2->Pz(),  evt->m_Pf2->E() );

  
  GenParticlePtr pe1 = std::make_shared<GenParticle>( pe1_v4, evt->m_KFini, 4);
  GenParticlePtr pe2 = std::make_shared<GenParticle>( pe2_v4, evt->m_KFini, 4);
  

  // discus with staszek which is positon wich is electon
  // ISR photons:
  auto ISR_photos=vector<GenParticlePtr>();
  for(int i=0; i < evt->m_nPhotISR, ++i)
    {

      FourVector tmp_photon_v4=FourVector( evt->m_PhotISR[i]->Px(),  evt->m_PhotISR[i]->Py(),  evt->m_m_PhotISR[i]->Pz(),  evt->m_PhotISR[i]->E() ); 
      GenParticlePtr tmp_photon=std::make_shared<GenParticle>( tmp_photon_v4, 22,1);
       
      double scalar1=tmp_photon_v4*pe1_v4;
      double scalar2=tmp_photon_v4*pe2_v4;
      
      

      GenParticlePtr e1STAR = std::make_shared<GenParticle>( pe1_v4-tmp_photon_v4, evt->m_KFini, 2); 
      GenParticlePtr e2STAR = std::make_shared<GenParticle>( pe2_v4-tmp_photon_v4, evt->m_KFini, 2);
      
      GenVertexPtr v_tmp = std::make_shared<GenVertex>();
      v_tmp->add_particle_in (pe1);
      if( scalar1 < scalar2)
	{
	  v_tmp->add_particle_out(e1STAR); // here I am cheating, fix this later, talk to Staszek if Evt class has momentum after FS or before
	  pe1=e1STAR;
	  pe1_v4=e1STAR.momentum();
	}
      else
	{
	  v_tmp->add_particle_out(e2STAR); 
	  pe2=e2STAR;
	  pe2_v4=e2STAR.momentum();
	}
      
      v_tmp->add_particle_out(tmp_photon);
      evt.add_vertex(v_tmp);

    }
  GenParticlePtr pZ = std::make_shared<GenParticle>( pe2_v4+pe1_v4, 23,	 2);                                                                                                                                                                                    
  GenVertexPtr vZ =  std::make_shared<GenVertex>();                                                                                                                                                                                            vZ ->add_particle_in(pe1);                                                                                                                                                                                                                  vZ ->add_particle_in(pe2);                                                                                                                                                                                                                  vZ ->add_particle_out(pZ);      
  
  auto FSR_photos=vector<GenParticlePtr>();

  FourVector pe3_v4(evt->m_Qf1->Px(),  evt->m_Qf1->Py(),  evt->m_Qf1->Pz(),  evt->m_Qf1->E() );
  FourVector pe4_v4(evt->m_Qf2->Px(),  evt->m_Qf2->Py(),  evt->m_Qf2->Pz(),  evt->m_Qf2->E() );


  GenParticlePtr pe3 = std::make_shared<GenParticle>( pe3_v4, evt->m_KFfni, 1);
  GenParticlePtr pe4 = std::make_shared<GenParticle>( pe4_v4, evt->m_KFfni, 1);

  
  

  
  for(int i=0; i < evt->m_nPhotFSR, ++i)           
    {
      FourVector tmp_photon_v4=FourVector( evt->m_PhotFSR[i]->Px(),  evt->m_PhotFSR[i]->Py(),  evt->m_m_PhotFSR[i]->Pz(),  evt->m_PhotFSR[i]->E() );
      GenParticlePtr tmp_photon=std::make_shared<GenParticle>( tmp_photon_v4, 22,1);

      GenParticlePtr e3STAR = std::make_shared<GenParticle>( pe3_v4-tmp_photon_v4, evt->m_K, 2);
      GenParticlePtr e4STAR = std::make_shared<GenParticle>( pe4_v4-tmp_photon_v4, evt->m_KFini, 2);

      
      
      double scalar3=tmp_photon_v4*pe3_v4;
      double scalar4=tmp_photon_v4*pe4_v4;
      
      GenVertexPtr v_tmp = std::make_shared<GenVertex>();
      
      if( scalar3 < scalar4)
	{
	  v_tmp->add_particle_in(e3STAR); // here I am cheating, fix this later, talk to Staszek if Evt class has momentum after FS or before 
	  v_tmp->add_particle_out(pe3);
	 
	  pe3=e3STAR;
	  pe3_v4=e3STAR.momentum();
	}
      else
	{
	  v_tmp->add_particle_in(e4STAR);
	  v_tmp->add_particle_out(pe4);
	  
	  pe4=e4STAR;
	  pe4_v4=e4STAR.momentum();

	}

      v_tmp->add_particle_out(tmp_photon);
      evt.add_vertex(v_tmp);


    }
  /*
  
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

  */
 






 

 
}
