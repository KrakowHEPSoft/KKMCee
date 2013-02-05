///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                CLASS JetAnalyzer                                          //
//                                                                           //
//    Concstructs 4 Jets in double W decay final state                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "JetAnalyzer.h"



ClassImp(JetAnalyzer)

JetAnalyzer::JetAnalyzer()
{
  //
  // Constructor of JetAnalyzer
  //
  hst_pene  = new TH1F("hst_pene" ,"energy parton",100,0,100);
  hst_jene  = new TH1F("hst_jene" ,"energy jet   ",100,0,100);
  // ntuple
  jtuple    = new TNtuple("jtuple","all",
	   "m0:m1:m2:m3:m4:m5:x0:x1:x2:x3:x4:x5:np:nj:wt:wt2:pan:jan:pen:jen");
}
//////////////////////////////////////////////////////////////////////////////
void JetAnalyzer::Book(KorEvent &event,
			long ntot, double wtBE, double wtBE2)
{
  //
  // Introduce kinematical cuts and fill histograms/ntuples
  //
     int nj = event.m_njet;
     double mWW[6];
     event.GetPartonMass(mWW);
     double xWW[6];
     event.GetJetMass(xWW);
////////////////////////////////////////////////
//    energies and angles of partons/jets     //
////////////////////////////////////////////////
     float pene[100];
     pene[0]=event.m_ferm1[0];
     pene[1]=event.m_ferm2[0];
     pene[2]=event.m_ferm3[0];
     pene[3]=event.m_ferm4[0];
     float jene[100];
     for (int j=0; j<nj ;j++) jene[j]= event.m_jet[j].m_pmom[0];
// calculate minimum energy of jets/quarks
     double min_pene = 1e10;
     for (j=0; j<4 ;j++)  if(pene[j]<min_pene) min_pene = pene[j];
     double min_jene = 1e10;
     for (j=0; j<nj ;j++) if(jene[j]<min_jene) min_jene = jene[j];
//
     hst_pene->FillN( 4,pene,NULL);
     hst_jene->FillN(nj,jene,NULL);
// minimum angle among partons and jets (for jets only for nj==4)
     double angles[6];
     double jet_min_angle =    event.GetJetAngle(angles);
     double par_min_angle =    event.GetPartonAngle(angles);
/////////////////////////////////////////
//         filling big ntuple          //
/////////////////////////////////////////
     float xft[100];
     for (j=0; j<6 ;j++) xft[j  ] = mWW[j];
     for (j=0; j<6 ;j++) xft[j+6] = xWW[j];
     xft[12] = ntot;
     xft[13] = nj;
     xft[14] = wtBE;
     xft[15] = wtBE2;
     xft[16] = par_min_angle;
     xft[17] = jet_min_angle;
     xft[18] = min_pene;
     xft[19] = min_jene;
     jtuple->Fill(xft);
}
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//           End of CLASS JetAnalyzer                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
