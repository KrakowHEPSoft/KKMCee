///////////////////////////////////////////////////////////////////////////////
#include "KKevent.h"

ClassImp(KKevent);

#define SW20 setw(20)<<setprecision(14)
#define SP15 setw(15)<<setprecision(9)
#define SP10 setw(10)<<setprecision(7)
#define SW10 setw(10)

KKevent::KKevent()
{
  // This constructor is for ROOT streamers ONLY
  // all pointers has to be NULLed
  cout<< "----> KKevent Default Constructor (for ROOT only) "<<endl;
  m_Out= NULL;
}

///_____________________________________________________________
KKevent::KKevent(ofstream *OutFile)
{
  cout<< "----> KKevent USER Constructor "<<endl;
  m_Out = OutFile;
}//KKevent

///______________________________________________________________________________________
KKevent::~KKevent()
{
  //Explicit destructor
  cout<< "----> KKevent::KKevent !!!! DESTRUCTOR !!!! "<<endl;
}///destructor

double KKevent::sqr( const Double_t x ){ return x*x;};

///______________________________________________________________________________________
void KKevent::Initialize(double CMSene)
{
  cout  << "----> KKevent::Initialize, Entering "<<endl;
//=================================================================
// BX*** macros are in MCdev/BXFORMAT.h
  BXOPE(*m_Out);
  BXTXT(*m_Out,"========================================");
  BXTXT(*m_Out,"======    KKevent::Initialize     ======");
  BXTXT(*m_Out,"========================================");
  m_CMSene = CMSene;
  m_nPhotISR=0;
  m_nPhotFSR=0;
  for(int i=0; i<maxPhot;i++){
	  m_PhotAll[i].SetPxPyPzE(0,0,0,0);
	  m_PhotISR[i].SetPxPyPzE(0,0,0,0);
	  m_PhotFSR[i].SetPxPyPzE(0,0,0,0);
  }// for i
  ///////////////////////////////////////////////////
}// Initialize


void KKevent::Merge(){
//////////////////////////////////////////////////////
//   Merging photos, fermions unchanged
//////////////////////////////////////////////////////
m_nPhot  = m_nPhotISR +m_nPhotFSR;
double enex,eney;
int i1=0, i2=0;               // c++ indexing!!!
for(int i=0; i<m_nPhot;i++){  // i within c++ indexing!!!
   enex = 0;
   eney = 0;
// saveguard against Alex Read and Tiziano Camporesi effect (bug in old BHLUMI)
//   for(int k=0; k<4; k++) xphot[k]= (m_Event->m_PhotISR[i1+1])[k];
//   for(int k=0; k<4; k++) yphot[k]= (m_Event->m_PhotFSR[i2+1])[k];
   if(i1 < m_nPhotISR) enex = m_PhotISR[i1+1].E();
   if(i2 < m_nPhotFSR) eney = m_PhotFSR[i2+1].E();
   if(enex > eney) {
      m_PhotAll[i+1] = m_PhotISR[i1+1];// f77 indexing
//   	  m_PhotAll[i+1].SetPxPyPzE(xphot[0],xphot[1],xphot[2],xphot[3]); // f77 indexing
      m_isr[i+1] = 1;   // ISR origin, f77 indexing !!!
      i1=i1 +1;
   } else {
      m_PhotAll[i+1] = m_PhotFSR[i2+1]; // f77 indexing
//   	  m_PhotAll[i+1].SetPxPyPzE(yphot[0],yphot[1],yphot[2],yphot[3]);
      m_isr[i+1] = 0;   // FSR origin, f77 indexing!!!
      i2=i2 +1;
   }
}// for i
m_PhotAll[0] = m_PhotISR[0]+m_PhotFSR[0]; // sum of momenta
}//KKevent::Merge

void KKevent::GivePair(double CMSene, double qm1, double qm2,
		      TLorentzVector *Pf1, TLorentzVector *Pf2, double *beta, double *eta1, double *eta2 ){
// defining pair of momenta
//KinLib_givpair(qmsene,Mas1,Mas2,m_q1, m_q2 ,betn,eta1,eta2)
 double svar  =  sqr(CMSene);
 *beta  =  sqrt((svar-sqr(qm1-qm2))*(svar-sqr(qm1+qm2)))/svar;
 *eta1  =  (svar+qm1*qm1-qm2*qm2)/svar;
 *eta2  =  (svar-qm1*qm1+qm2*qm2)/svar;
 double EE    =  CMSene/2;
 (*Pf1).SetPxPyPzE(0,0, (*beta)*EE,(*eta1)*EE);  //
 (*Pf2).SetPxPyPzE(0,0,-(*beta)*EE,(*eta2)*EE);  //
}//GivePair

void KKevent::DefPair(double CMSene, double qm1, double qm2, TLorentzVector *Pf1, TLorentzVector *Pf2){
// defining pair of momenta
 double svar  =  sqr(CMSene);
 double beta  =  sqrt((svar-sqr(qm1-qm2))*(svar-sqr(qm1+qm2)))/svar;
 double eta1  =  (svar+qm1*qm1-qm2*qm2)/svar;
 double eta2  =  (svar-qm1*qm1+qm2*qm2)/svar;
 double EE    =  CMSene/2;
 (*Pf1).SetPxPyPzE(0,0, beta*EE,eta1*EE);
 (*Pf2).SetPxPyPzE(0,0,-beta*EE,eta2*EE);
}//DefPair


void KKevent::BoostEul(double the, double phi, TLorentzVector *QQk, TLorentzVector *PX, TLorentzVector *pvec){
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
// Three transformations:                                                                  //
// (1) Boost from final fermions rest frame to ferms+phots rest frame (Z frame).           //
// (2) Euler rotation erasing memory of fermion directions.                                //
// (3) Boost to laboratory system.                                                         //
// Note that this transformation generates 'flat 2-body Born' in PX frame.                 //
// This is very easy to implement.                                                         //
// Of course, more sophisticated angular distribution can be implemented.                  //
// In such a case BostEul will be replaced with some other transformation.                 //
// Photon removal procedure with Piatek will work for arbitrary transformation.            //
//                                                                                         //
// Note that the procedure is encapsulated functionaly by generating Euler                 //
// angles localy, during first call for mode=-1.                                           //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////

BoostQ( -1,QQk,pvec); // Boost along qqk
(*pvec).RotateX(the);
(*pvec).RotateZ(phi);
BoostQ(  1,PX,pvec); // Boost along qqk
//CALL KinLib_Rotor(2,3,the,qvec,qvec); ! Rotation y-z
//CALL KinLib_Rotor(1,2,phi,qvec,qvec); ! Rotation x-y
//CALL KinLib_BostQ( -1,PX,qvec,qvec);  ! Boost to CMS
}

void KKevent::RotEuler(double alf, double bet, double gam, TLorentzVector *Pf){
//   SUBROUTINE KinLib_RotEuler(alfa,beta,gamma,pvec,qvec)
//////////////////////////////////////////////////////////////////////////////////
//      Full Euler rotation R_3(alpha)*R_2(beta)*R_3(gamma)                     //
//////////////////////////////////////////////////////////////////////////////////
//    CALL KinLib_Rotor(1,2, alfa,  pvec,qvec) ! z-rotation
//    CALL KinLib_Rotor(3,1, beta,  qvec,qvec) ! y-rotation
//    CALL KinLib_Rotor(1,2, gamma, qvec,qvec) ! z-rotation
//////////////////////////////////////////////////////////////////////////////////
	(*Pf).RotateZ(alf);
	(*Pf).RotateY(bet);
	(*Pf).RotateZ(gam);
cout<<" KKevent::RotEuler NOT TESTED YET!!!"<<endl;
exit(50);
}//RotEuler

void KKevent::RotEul(double the, double phi, TLorentzVector *Pf){
// SUBROUTINE KinLib_RotEul(the,phi,pvec,qvec)
//////////////////////////////////////////////////////////////////////////////////
//      Theta-phi rotation, it turns vector r along z-axis into r(theta,phi)    //
//////////////////////////////////////////////////////////////////////////////////
// CALL KinLib_Rotor(3,1,the,pvec,qvec) ! y-rotation
// CALL KinLib_Rotor(1,2,phi,qvec,qvec) ! z-rotation
	(*Pf).RotateY(the);
	(*Pf).RotateZ(phi);
}//RotEul

void KKevent::BoostQ(int Mode, TLorentzVector *Q, TLorentzVector *Pf){
TVector3 b = (*Q).BoostVector();
if(Mode <0) b=-b;
(*Pf).Boost(b);
}//BoostQ

void KKevent::PhaSpac2(TLorentzVector *PX, double the, double phi, double amfin,
		               TLorentzVector *Qf1, TLorentzVector *Qf2){
// define fermion pair in its rest frame PX along z-axis
 double PXM= (*PX).M();
 DefPair( PXM,amfin,amfin, Qf1, Qf2);
// rotate and boost
 RotEul(the,phi,Qf1);
 BoostQ(1, PX,  Qf1);
 RotEul(the,phi,Qf2);
 BoostQ(1, PX,  Qf2);
}//


void KKevent::ThetaR( double *cth11, double *cth12, double *cth21, double *cth22){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//   Provides cos(theta) for four scattering angles among p1,p2 and q1,q2          //
//   Angles are calculated in the rest frame of Qtot                               //
//   Called in programs calculating Born(p1,p2,q1,q2) distribution.                //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////

TLorentzVector Qtot = m_Qf1 + m_Qf2 + m_PhotFSR[0];  // m_PhotFSR[0] is sum of all

TLorentzVector p1=m_Pf1;
TLorentzVector p2=m_Pf2;
TLorentzVector q1=m_Qf1;
TLorentzVector q2=m_Qf2;

BoostQ(-1,&Qtot,&p1);
BoostQ(-1,&Qtot,&p2);
BoostQ(-1,&Qtot,&q1);
BoostQ(-1,&Qtot,&q2);
// Calculate all four possible angles
double q1d=   (q1).P(); //     sqrt(q1(1)**2 +q1(2)**2 +q1(3)**2)
double q2d=   (q2).P(); //     sqrt(q2(1)**2 +q2(2)**2 +q2(3)**2)
double p1d=   (p1).P(); //     sqrt(p1(1)**2 +p1(2)**2 +p1(3)**2)
double p2d=   (p2).P(); //     sqrt(p2(1)**2 +p2(2)**2 +p2(3)**2)

*cth11 = (p1).Vect() * (q1).Vect() /p1d/q1d;
*cth12 =-(p1).Vect() * (q2).Vect() /p1d/q2d;
*cth21 =-(p2).Vect() * (q1).Vect() /p2d/q1d;
*cth22 = (p2).Vect() * (q2).Vect() /p2d/q2d;

// angle = v1.Angle(v2.Vect(v2));
}//ThetaR


void KKevent::ThetaD(TLorentzVector &PX, double &CosTheta){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//   Provides Svar and  ONE cos(theta) among p1,p2 and q1,q2                       //
//   Energy conservation NOT required!!                                            //
//   The angle is beteween 3-vectors (p1-p2) and (q1-q2) in PX rest frame.         //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////

TLorentzVector pd = m_Pf1 - m_Pf2;
TLorentzVector qd = m_Qf1 - m_Qf2;

double Svar = PX*PX;
if( Svar <=  0.0 ) {
   cout<< "++++++++++ KKevent::ThetaD : PX not timelike"<<endl;
   exit(90);
}//if
double a = PX*pd;
double b = PX*qd;
pd = pd - (a/Svar)*PX;
qd = qd - (b/Svar)*PX;
a = pd*pd;
b = qd*qd;
if( a*b <= 0.0 ) {
   cout<< "++++++++++ KKevent::ThetaD: a,b="<<a<<"  "<<b<<endl;
   exit(90);
}
CosTheta = -(qd*pd)/sqrt(abs(a*b));
//
if( abs(CosTheta) > 1.0 ) {
   cout<< "++++++++++ KKevent::ThetaD: CosTheta= "<<CosTheta<<endl;
}//if
}//ThetaD



void KKevent::ZBoostALL(){
double exe   = sqrt((1-m_r1)/(1-m_r2));
double b3= (exe -1/exe)/(exe+1/exe);
m_Pf1.Boost(0,0,b3);
m_Pf2.Boost(0,0,b3);
m_Qf1.Boost(0,0,b3);
m_Qf2.Boost(0,0,b3);
for(int i=1; i<=m_nPhotISR;i++) m_PhotISR[i].Boost(0,0,b3);
for(int i=1; i<=m_nPhotFSR;i++) m_PhotFSR[i].Boost(0,0,b3);
for(int i=1; i<=m_nPhot;i++)    m_PhotAll[i].Boost(0,0,b3);
//
}//ZBoostALL

void KKevent::MomPrint( TLorentzVector &Vect){
//////////////////////////////////////////////////////////////
// printing entire four-vector in one line (withot endline)
  for ( int k=0; k < 4 ; k++ )   cout << SW20 << Vect[k];
}//MomPrint


void KKevent::PrintISR(){
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Prints out four momenta of MC event ISR                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

cout<<"||||||||||||||||||Event::PrintISR|||||||||||||||>  m_EventCounter="<< m_EventCounter<<endl;
cout<<"m_nPhotISR= "<<m_nPhotISR<<endl;
cout<<"m_WT_ISR = "<< m_WT_ISR<<"      ";
cout<<"m_KFini  = "<< m_KFini<< "    m_KFfin = "<< m_KFfin <<endl;
cout<<"m_vv = "<< m_vv   << "   m_r1 = "<< m_r1 <<"   m_r2 = "<<m_r2<<endl;

cout<< "m_Pf1  "; MomPrint(m_Pf1); cout<<"   "<< m_Pf1.M() <<endl;
cout<< "m_Pf2  "; MomPrint(m_Pf2); cout<<"   "<< m_Pf2.M() <<endl;

TLorentzVector BeamSum = m_Pf1+m_Pf2;
cout<< "BeamSum"; MomPrint(BeamSum); cout<<"   "<< BeamSum.M() <<endl;

cout<< "m_Qf1  "; MomPrint(m_Qf1); cout<<"   "<< m_Qf1.M() <<endl;
cout<< "m_Qf2  "; MomPrint(m_Qf2); cout<<"   "<< m_Qf2.M() <<endl;

TLorentzVector MomtSum = m_Qf1+m_Qf2;

for(int i=1; i<= m_nPhotISR; i++){
	cout<< "PhotISR"; MomPrint( m_PhotISR[i] ); cout<<"   "<< m_PhotISR[i].M() <<endl;
	MomtSum += m_PhotISR[i];
}

cout<< "m_Sum  "; MomPrint(MomtSum); cout<<"   "<< MomtSum.M() <<endl;
cout<<"----------------------------------------------------------------------------------------------------------------"<<endl;
}//EventPrint

void KKevent::PrintISR_FSR(){
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Prints out four momenta of MC event ISR and FSR photons                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

cout<<"///////////////////////////////Event::PrintISR+FSR///////////////////////////>  m_EventCounter="<< m_EventCounter<<endl;
cout<<"m_nPhotISR= "<<m_nPhotISR<<"  m_nPhotFSR= "<<m_nPhotFSR<<endl;
cout<<"m_WT_ISR = "<< m_WT_ISR<<"      ";
cout<<"m_KFini  = "<< m_KFini<< "    m_KFfin = "<< m_KFfin <<endl;
cout<<"m_vv = "<< m_vv   << "   m_r1 = "<< m_r1 <<"   m_r2 = "<<m_r2<<endl;

cout<< "m_Pf1  "; MomPrint(m_Pf1); cout<<"   "<< m_Pf1.M() <<endl;
cout<< "m_Pf2  "; MomPrint(m_Pf2); cout<<"   "<< m_Pf2.M() <<endl;

TLorentzVector BeamSum = m_Pf1+m_Pf2;
cout<< "BeamSum"; MomPrint(BeamSum); cout<<"   "<< BeamSum.M() <<endl;

cout<<"-------------------------------  final state -------------------------------------------------------------------"<<endl;
cout<< "m_Qf1  "; MomPrint(m_Qf1); cout<<"   "<< m_Qf1.M() <<endl;
cout<< "m_Qf2  "; MomPrint(m_Qf2); cout<<"   "<< m_Qf2.M() <<endl;

TLorentzVector MomtSum = m_Qf1+m_Qf2;

for(int i=1; i<= m_nPhotISR; i++){
	cout<< "PhotISR"; MomPrint( m_PhotISR[i] ); cout<<"   "<< m_PhotISR[i].M() <<endl;
	MomtSum += m_PhotISR[i];
}

for(int i=1; i<= m_nPhotFSR; i++){
	cout<< "PhotFSR"; MomPrint( m_PhotFSR[i] ); cout<<"   "<< m_PhotFSR[i].M() <<endl;
	MomtSum += m_PhotFSR[i];
}

cout<< "m_Sum  "; MomPrint(MomtSum); cout<<"   "<< MomtSum.M() <<endl;
cout<<"----------------------------------------------------------------------------------------------------------------"<<endl;
}//EventPrint

void KKevent::MomPrint(ofstream *Out, TLorentzVector &Vect){
//////////////////////////////////////////////////////////////
// printing entire four-vector in one line (withot endline)
  for ( int k=0; k < 4 ; k++ )   *Out << SW20 << Vect[k];
}//MomPrint


void KKevent::PrintISR_FSR(ofstream *Out){
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Prints out four momenta of MC event ISR and FSR photons                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

*Out<<"//////////////////////////Event::PrintISR+FSR///////////////////////////////>  m_EventCounter="<< m_EventCounter<<endl;
*Out<<"m_nPhotISR= "<<m_nPhotISR<<"  m_nPhotFSR= "<<m_nPhotFSR<<endl;
*Out<<"m_WT_ISR = "<< m_WT_ISR<<"      ";
*Out<<"m_KFini  = "<< m_KFini<< "    m_KFfin = "<< m_KFfin <<endl;
*Out<<"m_vv = "<< m_vv   << "   m_r1 = "<< m_r1 <<"   m_r2 = "<<m_r2<<endl;

*Out<< "m_Pf1  "; MomPrint(Out,m_Pf1); *Out<<"   "<< m_Pf1.M() <<endl;
*Out<< "m_Pf2  "; MomPrint(Out,m_Pf2); *Out<<"   "<< m_Pf2.M() <<endl;

TLorentzVector BeamSum = m_Pf1+m_Pf2;
*Out<< "BeamSum"; MomPrint(Out,BeamSum); *Out<<"   "<< BeamSum.M() <<endl;

*Out<<"-------------------------------  final state -------------------------------------------------------------------"<<endl;

*Out<< "m_Qf1  "; MomPrint(Out,m_Qf1); *Out<<"   "<< m_Qf1.M() <<endl;
*Out<< "m_Qf2  "; MomPrint(Out,m_Qf2); *Out<<"   "<< m_Qf2.M() <<endl;
TLorentzVector MomtSum = m_Qf1+m_Qf2;

for(int i=1; i<= m_nPhotISR; i++){
	*Out<< "PhotISR"; MomPrint(Out, m_PhotISR[i] ); *Out<<"   "<< m_PhotISR[i].M() <<endl;
	MomtSum += m_PhotISR[i];
}

for(int i=1; i<= m_nPhotFSR; i++){
	*Out<< "PhotFSR"; MomPrint(Out, m_PhotFSR[i] ); *Out<<"   "<< m_PhotFSR[i].M() <<endl;
	MomtSum += m_PhotFSR[i];
}

*Out<< "m_Sum  "; MomPrint(Out,MomtSum); *Out<<"   "<< MomtSum.M() <<endl;
*Out<<"----------------------------------------------------------------------------------------------------------------"<<endl;
}//EventPrint

void KKevent::EventPrintAll(){
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Prints out four momenta of Beams and Final state particles,              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

cout<<"==================KKevent::EventPrintAll=====================>  m_EventCounter="<< m_EventCounter<<endl;
cout<<"m_nPhot   = "<<m_nPhot<<endl;
cout<<"m_WtMain = "<< m_WtMain<<"  m_WtCrude =  "<<m_WtCrude<<"      ";
cout<<"m_KFini  = "<< m_KFini<< "    m_KFfin = "<< m_KFfin <<endl;
cout<<"m_vv = "<< m_vv   << "   m_r1 = "<< m_r1 <<"   m_r2 = "<<m_r2<<endl;

cout<< "m_Pf1  "; MomPrint(m_Pf1); cout<<"   "<< m_Pf1.M() <<endl;
cout<< "m_Pf2  "; MomPrint(m_Pf2); cout<<"   "<< m_Pf2.M() <<endl;

TLorentzVector BeamSum = m_Pf1+m_Pf2;
cout<< "BeamSum"; MomPrint(BeamSum); cout<<"   "<< BeamSum.M() <<endl;

cout<< "m_Qf1  "; MomPrint(m_Qf1); cout<<"   "<< m_Qf1.M() <<endl;
cout<< "m_Qf2  "; MomPrint(m_Qf2); cout<<"   "<< m_Qf2.M() <<endl;

TLorentzVector MomtSum = m_Qf1+m_Qf2;
for(int i=1; i<= m_nPhot; i++){
	cout<< "m_Ph "<<m_isr[i]<<" ";
	    MomPrint( m_PhotAll[i] ); cout<<"   "<< m_PhotAll[i].M() <<endl;
	MomtSum += m_PhotAll[i];
}
cout<< "m_Sum  "; MomPrint(MomtSum); cout<<"   "<< MomtSum.M() <<endl;
cout<<"----------------------------------------------------------------------------------------------------------------"<<endl;
}//EventPrintAll



void KKevent::EventPrintAll(ofstream *Out){
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Prints out four momenta of Beams and Final state particles,              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

*Out<<"==================KKevent::EventPrintAll=====================>  m_EventCounter="<< m_EventCounter<<endl;
*Out<<"m_nPhot  =  "<<m_nPhot<<endl;
*Out<<"m_WtMain = "<< m_WtMain<<"  m_WtCrude = "<<m_WtCrude<<"      ";
*Out<<"m_KFini  = "<< m_KFini<< "    m_KFfin = "<< m_KFfin <<endl;
*Out<<"m_vv = "<< m_vv   << "   m_r1 = "<< m_r1 <<"   m_r2 = "<<m_r2<<endl;

*Out<< "m_Pf1  "; MomPrint(Out,m_Pf1); *Out<<"   "<< m_Pf1.M() <<endl;
*Out<< "m_Pf2  "; MomPrint(Out,m_Pf2); *Out<<"   "<< m_Pf2.M() <<endl;

*Out<< "m_Qf1  "; MomPrint(Out,m_Qf1); *Out<<"   "<< m_Qf1.M() <<endl;
*Out<< "m_Qf2  "; MomPrint(Out,m_Qf2); *Out<<"   "<< m_Qf2.M() <<endl;

TLorentzVector MomtSum = m_Qf1+m_Qf2;
for(int i=1; i<= m_nPhot; i++){
	*Out<< "m_Ph "<<m_isr[i]<<" ";
	    MomPrint(Out, m_PhotAll[i] ); *Out<<"   "<< m_PhotAll[i].M() <<endl;
	MomtSum += m_PhotAll[i];
}
*Out<< "m_Sum  "; MomPrint(Out,MomtSum); *Out<<"   "<< MomtSum.M() <<endl;
*Out<<"-------------------------------initial hadron remnants ---------------------------------------------------------"<<endl;
*Out<< "m_Rem1 "; MomPrint(Out,m_Rem1); *Out<<"   "<< m_Rem1.M() <<endl;
*Out<< "m_Rem2 "; MomPrint(Out,m_Rem2); *Out<<"   "<< m_Rem2.M() <<endl;
MomtSum += m_Rem1+m_Rem2;
*Out<< "SumAll "; MomPrint(Out,MomtSum); *Out<<"   "<< MomtSum.M() <<endl;
*Out<<"----------------------------------------------------------------------------------------------------------------"<<endl;
}//EventPrintAll



