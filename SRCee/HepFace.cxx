///////////////////////////////////////////////////////////////////////////////
#include "HepFace.h"

ClassImp(HepFace);

/*
double Scalar_product(FourVector a)
{
  double val=0.;
  val=a.e*a.e-(a.px*a.px+a.py*a.py+a.pz*a.pz);
  return val;

  
}
*/

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
  // clear tau container
  tauMdecay.clear();
  tauPdecay.clear();
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
  GenParticlePtr pe2 = std::make_shared<GenParticle>( pe2_v4, - m_Event->m_KFini, 4); // - PDG id antifermion 

// ISR photons:
  for(int i=0; i < m_Event->m_nPhotISR; ++i)
    {
      // create a HEPMC photon
      FourVector tmp_photon_v4=Vect4(m_Event->m_PhotISR[i]);
      GenParticlePtr tmp_photon=std::make_shared<GenParticle>( tmp_photon_v4, 22,1);

      // assign from which electron the photon was emmited 
      
      // we check which electron is closer to the photon 
      /////double scalar1=Scalar_product( tmp_photon_v4 + pe1_v4);
      /////double scalar2=Scalar_product( tmp_photon_v4 + pe2_v4);
      double scalar1= m_Event->m_Pf1 * m_Event->m_PhotISR[i]; // iloczyn skalarny TLorentzvector
      double scalar2= m_Event->m_Pf2 * m_Event->m_PhotISR[i]; // iloczyn skalarny TLorentzvector

      // electrons after the emission of the photon - change of kinematics 	
      GenParticlePtr e1STAR = std::make_shared<GenParticle>( pe1_v4-tmp_photon_v4, m_Event->m_KFini, 2); 
      GenParticlePtr e2STAR = std::make_shared<GenParticle>( pe2_v4-tmp_photon_v4, - m_Event->m_KFini, 2); // - PDG id antifermion 

      // create a vertex with the selected electron	
      GenVertexPtr v_tmp = std::make_shared<GenVertex>();
      if( scalar1 < scalar2)
        {
	  v_tmp->add_particle_out(e1STAR);
          v_tmp->add_particle_in(pe1);	  
          pe1=e1STAR;
          pe1_v4=e1STAR->data().momentum;
        }
      else
	{
          v_tmp->add_particle_out(e2STAR);
          v_tmp->add_particle_in(pe2); 
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
  GenParticlePtr pe4 = std::make_shared<GenParticle>( pe4_v4, - m_Event->m_KFfin, status);

  for(int i=0; i< m_Event->m_nPhotFSR; i++){
      // create a HEPMC photon
      FourVector tmp_photon_v4=FourVector(Vect4(m_Event->m_PhotFSR[i]));
      GenParticlePtr tmp_photon=std::make_shared<GenParticle>( tmp_photon_v4, 22, 3);

      // assign from which electron the photon was emmited 

      // we will check which electron is more parallel to the photon 
      // for now we check the R-rap-distance separation dR = sqrt(dphi^2 + drap^2) 

      /////double scalar3=Scalar_product( tmp_photon_v4 + pe3_v4);
      /////double scalar4=Scalar_product( tmp_photon_v4 + pe3_v4);
      double scalar3= m_Event->m_PhotFSR[i] * m_Event->m_Qf1;
      double scalar4= m_Event->m_PhotFSR[i] * m_Event->m_Qf2;
      
      // electrons before the emission of the photon - change of kinematics
      GenParticlePtr e3STAR = std::make_shared<GenParticle>( pe3_v4+tmp_photon_v4, m_Event->m_KFfin, 2);	
      GenParticlePtr e4STAR = std::make_shared<GenParticle>( pe4_v4+tmp_photon_v4, - m_Event->m_KFfin, 2);

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
//  Print::listing(*m_Hvent);
//  Print::content(*m_Hvent);

}//make1


////////////////////////////////////////////////////////////////////////////
//SUBROUTINE HepEvt_Fil1( n,ist,id,jmo1,jmo2,jda1,jda2,p4,pinv,phflag)
//INTEGER           n,ist,id,jmo1,jmo2,jda1,jda2
//DOUBLE PRECISION  p4(4),pinv
//LOGICAL           phflag
//SUBROUTINE FILHEP(N,IST,ID,JMO1,JMO2,JDA1,JDA2,P4,PINV,PHFLAG)
//LOGICAL PHFLAG
//DOUBLE PRECISION  PINVD,P4D(4)
//C anti electron neutrino (nu_e is 12)
//   CALL FILHEP(0,1,-12*ISGN,NPS,NPS,0,0,PNE,AM,.TRUE.)
//--------------------------------------------------------------------------

///______________________________________________________________________________________
void HepFace::FillHep3(int N, int IST, int ID, int JMO1, int JMO2, int JDA1, int JDA2, float P4[], float &PINV, bool PHFLAG)
{
  // Create a HEPMC3 particle from a tau decay
  FourVector p_v4( Vect4(P4) );
  GenParticlePtr ptemp = std::make_shared<GenParticle>( p_v4, ID, 2);
  
  //// store it in the vector which will be used later to fill the HEPMC record 
  // N = 1 is tau-
  if (N == 1) tauMdecay.push_back(ptemp);
  // N = -1 is tau+
  else if (N == -1) tauPdecay.push_back(ptemp);
  //
  else cout << "HepFace::FillHep3 something went wrong" << endl;
  /*
  cout<<"==================================FillHep3=============================================="<<endl;
  cout<<"   N= "<< N<<"   IST= "<< IST<<"   ID= "<< ID<<endl;
  cout<<"   JMO1= "<< JMO1<<"   JMO2= "<< JMO2<<"   JDA1= "<< JDA1<<"   JDA2= "<< JDA2<<endl;
  for(int i=0; i<4; i++) cout<<"  P4["<<i<<"]=  "<<P4[i]; cout<<endl;
  cout<<"  Mass PINV="<<PINV<< "   PHFLAG= "<<PHFLAG<<endl;
  cout<<"========================================================================================"<<endl;
  */
}//FillHep3


///////////////////////////////////////////////////////////////////////////////////
// tauolaToHEPMC3() after the tau decays are done this function saves it to HEPMC3 
void HepFace::tauolaToHEPMC3(){
  GenVertexPtr tauPlus =   std::make_shared<GenVertex>();   
  GenVertexPtr tauMinus =  std::make_shared<GenVertex>();

  for (auto p: m_Hvent->particles()){
     if (p->pid() == 15 && p->status()==1){ // is ustable tau-
        // add tau to the vertex
	tauMinus->add_particle_in(p);
        // add decay products of the tau to the vertex
        for (auto i : tauMdecay)
        {
	   tauMinus->add_particle_out(i);
        }
     } // end if tau -
     else if (p->pid() == -15 && p->status()==1 ){ // is unstable tau+
        // add tau to the vertex
        tauPlus->add_particle_in(p);
        // add decay products of the tau to the vertex
	for (auto i : tauPdecay)
        {
           tauPlus->add_particle_out(i);
	}

     }  // end if tau+
  } // end loop over particles
 
  // add vertex to the event
  m_Hvent->add_vertex(tauMinus);
  m_Hvent->add_vertex(tauPlus);

  // clear tau container
  tauMdecay.clear();
  tauPdecay.clear();

  // Test print out
  cout<<"==================================tauolaToHEPMC3==============================================";
  Print::listing(*m_Hvent);
  Print::content(*m_Hvent);

}


