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

// beam particles  

  FourVector pe1_v4( Vect4( m_Event->m_Pf1) );
  FourVector pe2_v4( Vect4( m_Event->m_Pf2) );
  
  GenParticlePtr pe1 = std::make_shared<GenParticle>( pe1_v4, m_Event->m_KFini, 4);
  GenParticlePtr pe2 = std::make_shared<GenParticle>( pe2_v4, m_Event->m_KFini, 4);

// ISR photons:
  for(int i=0; i < m_Event->m_nPhotISR; ++i)
    {
      // create a HEPMC photon
      FourVector tmp_photon_v4=Vect4(m_Event->m_PhotISR[i]);
      GenParticlePtr tmp_photon=std::make_shared<GenParticle>( tmp_photon_v4, 22,1);

      // assign from which electron the photon was emmited 
      
      // we will check which electron is more parallel to the photon 
      // for now we check the R-rap-distance separation dR = sqrt(dphi^2 + drap^2) 
      double scalar1=tmp_photon_v4.delta_r_rap(pe1_v4);
      double scalar2=tmp_photon_v4.delta_r_rap(pe2_v4);

      // electrons after the emission of the photon - change of kinematics 	
      GenParticlePtr e1STAR = std::make_shared<GenParticle>( pe1_v4-tmp_photon_v4, m_Event->m_KFini, 2); 
      GenParticlePtr e2STAR = std::make_shared<GenParticle>( pe2_v4-tmp_photon_v4, m_Event->m_KFini, 2);

      // create a vertex with the selected electron	
      GenVertexPtr v_tmp = std::make_shared<GenVertex>();
      v_tmp->add_particle_in (pe1);
      if( scalar1 < scalar2)
        {
	  v_tmp->add_particle_out(e1STAR);
          pe1=e1STAR;
          pe1_v4=e1STAR->data().momentum;
        }
      else
	{
          v_tmp->add_particle_out(e2STAR); 
          pe2=e2STAR;
          pe2_v4=e2STAR->data().momentum;
        }

      v_tmp->add_particle_out(tmp_photon);
      m_Hvent->add_vertex(v_tmp);

    } // End ISR photons

// intermediate Z bozon:

  GenParticlePtr pZ = std::make_shared<GenParticle>( pe2_v4 + pe1_v4 , 23,  2);
  GenVertexPtr vZ =  std::make_shared<GenVertex>();  

  vZ ->add_particle_in(pe1);
  vZ ->add_particle_in(pe2);
  vZ ->add_particle_out(pZ);

  m_Hvent->add_vertex(vZ);

// FSR photons:

  // we will start from final state fermions
  FourVector pe3_v4( Vect4( m_Event->m_Qf1) );
  FourVector pe4_v4( Vect4( m_Event->m_Qf2) );

  GenParticlePtr pe3 = std::make_shared<GenParticle>( pe3_v4, m_Event->m_KFfin, status);
  GenParticlePtr pe4 = std::make_shared<GenParticle>( pe4_v4, m_Event->m_KFfin, status);





  for(int i=0; i< m_Event->m_nPhotFSR; i++){
      // create a HEPMC photon
      FourVector tmp_photon_v4=FourVector(Vect4(m_Event->m_PhotFSR[i]));
      GenParticlePtr tmp_photon=std::make_shared<GenParticle>( tmp_photon_v4, 22, 3);

      // assign from which electron the photon was emmited 

      // we will check which electron is more parallel to the photon 
      // for now we check the R-rap-distance separation dR = sqrt(dphi^2 + drap^2) 
      double scalar3=tmp_photon_v4.delta_r_rap(pe3_v4);
      double scalar4=tmp_photon_v4.delta_r_rap(pe4_v4);

      // electrons before the emission of the photon - change of kinematics
      GenParticlePtr e3STAR = std::make_shared<GenParticle>( pe3_v4+tmp_photon_v4, m_Event->m_KFini, 2);	
      GenParticlePtr e4STAR = std::make_shared<GenParticle>( pe4_v4+tmp_photon_v4, m_Event->m_KFini, 2);

      // create a vertex with the selected electron     
      GenVertexPtr v_tmp = std::make_shared<GenVertex>();
      if( scalar3 < scalar4)
        {
          v_tmp->add_particle_in(e3STAR); 
          v_tmp->add_particle_out(pe3);

          pe3=e3STAR;
          pe3_v4=e3STAR->data().momentum;
        }
      else
	{
          v_tmp->add_particle_in(e4STAR);
          v_tmp->add_particle_out(pe4);

          pe4=e4STAR;
          pe4_v4=e4STAR->data().momentum;

        }

      v_tmp->add_particle_out(tmp_photon);
      m_Hvent->add_vertex(v_tmp);

  }// for i


// now decaying Z:
  GenVertexPtr vZ2 =  std::make_shared<GenVertex>();
  vZ2 ->add_particle_in(pZ);
  vZ2->add_particle_out(pe3);
  vZ2->add_particle_out(pe4);
  m_Hvent->add_vertex(vZ2);

// Test print out
  Print::listing(*m_Hvent);
  Print::content(*m_Hvent);

}//make1
