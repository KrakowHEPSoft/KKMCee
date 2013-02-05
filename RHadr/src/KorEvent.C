//////////////////////////////////////////////////////////////////////////////
//                     CLASS   KorEvent                                     //
//                                                                          //
// This is class of KORALW events                                           //
// after HADRONIZATION encoded in Lund/PDG coding style                     //
//                                                                          //
//                 end of description   KorEvent                            //
//////////////////////////////////////////////////////////////////////////////
# define sw0 setw(3)
# define sw1 setw(7)
# define sw2 setprecision(10) << setw(18)

#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include <math.h>

#include "KorEvent.h"


ClassImp(KorEvent)
//////////////////////////////////////////////////////////////////////////////
KorEvent::KorEvent()
{
// My Constructor
}
//////////////////////////////////////////////////////////////////////////////
KorEvent::~KorEvent()
{
// My Destructor
}
//////////////////////////////////////////////////////////////////////////////
void KorEvent::Print( long Level)
{
// KorEvent printing
  int j;

  //print(); // inherited printout from class KorEvent 

  cout <<"|||||||||||||||||||||||KORALW partons |||||||||||||||||||||||"<<endl;
  cout << sw2 << "nphot= " << sw2 << m_nphot << endl;
  for ( j=0; j < m_nphot ; j++ )
    {
      cout << sw0<< j; m_photmom[j].print(); cout << endl;
    }
  cout << sw0<< j; j++;  m_wminus.print(); cout << endl;
  cout << sw0<< j; j++;  m_wplus.print();  cout << endl;
  cout << sw0<< j; j++;  m_ferm1.print();  cout << endl;
  cout << sw0<< j; j++;  m_ferm2.print();  cout << endl;
  cout << sw0<< j; j++;  m_ferm3.print();  cout << endl;
  cout << sw0<< j; j++;  m_ferm4.print();  cout << endl;


  cout << "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||";
  cout << "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||"<<endl;
  //  cout << "npart = " << m_npart << "\n";
# define sw0 setw(3)
  cout << sw0 <<  "j";
  cout << sw1 <<  "status" << sw1 <<  "flavor" << sw1 <<  "parent" ;
  cout << sw2 << "pmom[0]" << sw2 << "pmom[1]" << sw2 <<  "pmom[2]" 
       << sw2 << "pmom[3]" << sw2 << "pmass" << endl;
  for ( j=0; j < m_npart ; j++ )
    {
      //cout << sw0 << j;
      m_part[j].print();
    }
  for ( j=0; j < m_njet ; j++ )
    {
      //cout << sw0 << j;
      m_jet[j].print();
    }
}
//---------------------------------------------------------------------------//
//                    Functions for 4-jet analysis                           //
//---------------------------------------------------------------------------//
void KorEvent::GetPartonMass(double m[6] )
{
  VLorenz psum;
  for(int i=0; i<6; i++) m[i]=0;
  psum = m_ferm1 + m_ferm2; m[0] = sqrt(psum*psum);
  psum = m_ferm3 + m_ferm4; m[1] = sqrt(psum*psum);
  psum = m_ferm1 + m_ferm3; m[2] = sqrt(psum*psum);
  psum = m_ferm2 + m_ferm4; m[3] = sqrt(psum*psum);
  psum = m_ferm1 + m_ferm4; m[4] = sqrt(psum*psum);
  psum = m_ferm2 + m_ferm3; m[5] = sqrt(psum*psum);
};
//---------------------------------------------------------------------------
double KorEvent::GetPartonAngle(double angle[6] )
{
  for(int i=0; i<6; i++) angle[i]=0;
  angle[0] = m_ferm1 / m_ferm2;
  angle[1] = m_ferm3 / m_ferm4;
  angle[2] = m_ferm1 / m_ferm3;
  angle[3] = m_ferm2 / m_ferm4;
  angle[4] = m_ferm1 / m_ferm4;
  angle[5] = m_ferm2 / m_ferm3;
  double min_angle=10000;
  for(i=0; i<6; i++) 
    if(angle[i]< min_angle)  min_angle = angle[i];
  return min_angle;
};
//---------------------------------------------------------------------------
void KorEvent::GetJetMass(double m[6] )
{
  VLorenz psum;
  int nj = m_njet;
  for(int i=0; i<6; i++) m[i]=0;
  if(nj==4)
    {
      psum= m_jet[0].m_pmom +m_jet[1].m_pmom; m[0] =  sqrt(psum*psum);
      psum= m_jet[2].m_pmom +m_jet[3].m_pmom; m[1] =  sqrt(psum*psum);
      psum= m_jet[0].m_pmom +m_jet[2].m_pmom; m[2] =  sqrt(psum*psum);
      psum= m_jet[1].m_pmom +m_jet[3].m_pmom; m[3] =  sqrt(psum*psum);
      psum= m_jet[0].m_pmom +m_jet[3].m_pmom; m[4] =  sqrt(psum*psum);
      psum= m_jet[1].m_pmom +m_jet[2].m_pmom; m[5] =  sqrt(psum*psum);
    }
};
//---------------------------------------------------------------------------
double KorEvent::GetJetAngle(double angle[6] )
{
  int nj = m_njet;
  double min_angle=10000;
  for(int i=0; i<6; i++) angle[i]=0;
  if( nj == 4 )
    {
      angle[0]= m_jet[0].m_pmom /m_jet[1].m_pmom;
      angle[1]= m_jet[2].m_pmom /m_jet[3].m_pmom;
      angle[2]= m_jet[0].m_pmom /m_jet[2].m_pmom;
      angle[3]= m_jet[1].m_pmom /m_jet[3].m_pmom;
      angle[4]= m_jet[0].m_pmom /m_jet[3].m_pmom;
      angle[5]= m_jet[1].m_pmom /m_jet[2].m_pmom;
      for(i=0; i<6; i++) 
	if(angle[i]< min_angle)  min_angle = angle[i];
    }
  return min_angle;
};
///////////////////////////////////////////////////////////////////////////////
//                   END OF CLASS KorEvent                                   //
///////////////////////////////////////////////////////////////////////////////
