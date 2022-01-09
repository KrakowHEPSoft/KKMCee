///////////////////////////////////////////////////////////////////////////////
#include "KKarLud.h"

ClassImp(KKarLud);

#define SW20 setw(20)<<setprecision(14)

KKarLud::KKarLud()
{
  // This constructor is for ROOT streamers ONLY
  // all pointers has to be NULLed
  cout<< "----> KKarLud Default Constructor (for ROOT only) "<<endl;
  m_Out   = NULL;
  DB      = NULL;
  m_RNgen = NULL;
  m_Event = NULL;
}

///_____________________________________________________________
KKarLud::KKarLud(ofstream *OutFile)
{
  cout<< "----> KKarLud USER Constructor "<<endl;
  m_Out   = OutFile;
  DB      = NULL;
  m_RNgen = NULL;
  m_Event = NULL;
}//KKarLud

///______________________________________________________________________________________
KKarLud::~KKarLud()
{
  //Explicit destructor
  cout<< "----> KKarLud::KKarLud !!!! DESTRUCTOR !!!! "<<endl;
}///destructor

double KKarLud::sqr( const Double_t x ){ return x*x;};

///______________________________________________________________________________________
void KKarLud::Initialize()
{
  cout  << "----> KKarLud::Initialize, Entering "<<endl;
  m_icont = 0;
  m_nfail = 0;
//=================================================================
// BX*** macros are in MCdev/BXFORMAT.h
  BXOPE(*m_Out);
  BXTXT(*m_Out,"========================================");
  BXTXT(*m_Out,"======    KKarLud::Initialize     ======");
  BXTXT(*m_Out,"========================================");
  m_vvmin = DB->vvmin;
  m_CMSene= DB->CMSene;  // to be moved

  ///////////////////////////////////////////////////
}// Initialize

void KKarLud::Make(TLorentzVector *PX, double *wt_ISR){
//////////////////////////////////////////////////////////////////////////
// generating ISR photons
//////////////////////////////////////////////////////////////////////////
  m_icont++;
// import variables generated in FOAM
  m_KFini  = m_Event->m_KFini;
  m_KFfin  = m_Event->m_KFfin;
  m_vv     = m_Event->m_vv;
  m_r1     = m_Event->m_r1;
  m_r2     = m_Event->m_r2;
  m_XXXene = m_Event->m_XXXene; // = sqrt(PX*PX)
  m_AvMult = m_Event->m_AvMult;
//     Redefine CMS energy and boost
//  m_XXXene = DB->CMSene*sqrt((1-m_r1)*(1-m_r2));
//  m_Event->m_XXXene = m_XXXene;

// usually beam is electron but muon is also allowed
  m_amel= DB->fmass[m_KFini];  // current masses

// defining incomming (beam) parton momenta
//  double  qm1 = m_amel, qm2 = m_amel;
//  m_Event->DefPair(m_XXXene,m_amel,m_amel, &(m_Event->m_Pf1), &(m_Event->m_Pf2));

  // Low-level multiphoton ISR generator (1986)
  if (DB->KeyISR == 1) {
       YFSgen(m_XXXene, m_vv, &(m_Event->m_PX), &(m_Event->m_WT_ISR) );
	   *PX     = m_Event->m_PX;      // temporary solution
	   *wt_ISR = m_Event->m_WT_ISR;  // temporary solution?
  } else if(DB->KeyISR == 0) {
	    m_Event->DefPair( m_XXXene,m_amel,m_amel, &(m_Event->m_Qf1), &(m_Event->m_Qf2));
	    double cth= 1 -2*m_RNgen->Rndm();
	    double the= acos(cth);
	    double phi= 2*M_PI *m_RNgen->Rndm();
	    m_Event->RotEul(the,phi,&(m_Event->m_Qf1));
	    m_Event->RotEul(the,phi,&(m_Event->m_Qf2));
	    m_nphot = 0;
	    *wt_ISR = 1;
	    *PX = m_Event->m_Qf1 + m_Event->m_Qf2;
   } else {
       cout<< " ++++ KKarLud: wrong KeyISR="<< DB->KeyISR <<endl;
       exit(95);
   }
//-------------------------------------------------------------
   if(*wt_ISR == 0 ) {
 // Set momenta to zero for WT=0 events
      m_Event->m_Qf1.SetPxPyPzE(0,0,0,0);
      m_Event->m_Qf2.SetPxPyPzE(0,0,0,0);
      m_nphot=0;
      for(int j=0; j<=maxPhot; j++)
    	  m_Event->m_PhotISR[j].SetPxPyPzE(0,0,0,0);
    } else {
  // Define final fermion momenta (NOT used in case of FSR!)
      double cth= 1 -2*m_RNgen->Rndm();
      double the= acos(cth);
      double phi= 2*M_PI*m_RNgen->Rndm();
      double amfi  =DB->fmass[m_KFfin];
      m_Event->PhaSpac2(PX,the,phi,amfi, &(m_Event->m_Qf1), &(m_Event->m_Qf2));
   }// wt_ISR == 0
}// Make



void KKarLud::YFSgen(double XXXene, double vv, TLorentzVector *PX, double *WtIni){
////////////////////////////////////////////////////////////////////////////////////
//                                                                                //
//  ======================================================================        //
//  ======================= Y F S G E N ==================================        //
//  ======================================================================        //
//  The algorithm in this subprogram was described in:                            //
//  ``Yennie-Frautschi-Suura soft photons in Monte Carlo event generators''       //
//             Unpublished report by S. Jadach,                                   //
//          MPI-Munchen, MPI-PAE/PTh 6/87, Jan. 1987.                             //
//                                                                                //
//  Later on used in YFS1,YFS2,YFS3, YFSWW, KORALZ, KORALW Monte Carlo programs   //
//                                                                                //
//  Purpose:  ISR photon emission, photon multiplicity and momenta                //
//                                                                                //
////////////////////////////////////////////////////////////////////////////////////
//   INPUT:    XXXene,vv                                                          //
//   OUTPUT:   PX,WtIni                                                           //
//                                                                                //
//   XXXene  = total cms energy                                                   //
//   amel    = beam mass                                                          //
//   MltISR  = flag normaly set to zero, for SPECIAL tests enforces               //
//             photon multiplicity to be exactly equal MltISR                     //
//   vv      = v=1-s'/s variable                                                  //
//   vmin    = minimum v variable (infrared cutoff)                               //
//   nmax    = maximum photon multiplicity                                        //
//   alfinv  = 1/apha_QED                                                         //
//   p1,2    = initial fermion momenta (along z-axix)                             //
//   nphot   = photon multiplicity                                                //
//   sphot   = photon four-momenta                                                //
//   sphum   = total photon four-momentum                                         //
//   ygr,zet = Sudakov variables                                                  //
//   PX      = 4-mmentum left after photon emission                               //
//   WtIni   = total weight from this sub-generator                               //
////////////////////////////////////////////////////////////////////////////////////

TLorentzVector   xphot[maxPhot];    // photon momenta before rescaling
double           xph[maxPhot], rr[maxPhot];
TLorentzVector   pp,pk;
//---------------------------------------
double Ene  = XXXene/2;

for(int i=0; i<maxPhot;i++){
   xph[i]=0;
   m_yini[i]=0;
   m_zini[i]=0;
   xphot[i].SetPxPyPzE(0,0,0,0);
   m_Event->m_PhotISR[i].SetPxPyPzE(0,0,0,0);
}// for i

double m_WtMass, m_WtDil, m_WtCut; // local!!!
if(vv <= m_vvmin) {
///////////////////////////////////////////////////
//    no photon above detectability threshold    //
///////////////////////////////////////////////////
   m_WtMass = 1;
   m_WtDil  = 1;
   m_WtCut  = 1;
   m_nphot  = 0;
   m_Event->m_nPhotISR = 0;
} else {
/////////////////////////////////////////////////////////
// one or more photons, generate photon multiplicity   //
// nphot = poisson(AvMult) + 1                         //
/////////////////////////////////////////////////////////
   double AvMult = m_AvMult;
//e100:
   PoissGen(AvMult,&m_nphot, rr);
  m_nphot = m_nphot+1;
   m_Event->m_nPhotISR = m_nphot;
// For special tests of program at fixed multiplicity (for advc. users)
//   if( (m_MltISR != 0) && (m_nphot != m_MltISR) ) goto e100;
   if(m_nphot == 1) {
      xph[1]=vv;
   } else {
      xph[1]=vv;
      for(int i=2; i<=m_nphot; i++) xph[i]=vv*exp(log(m_vvmin/vv)*rr[i-1]);
   }// if nphot
   m_WtMass=1;
   for(int i=1; i<=m_nphot; i++){
      double xk=xph[i];
      double am2  = sqr(m_amel/Ene);
      double del1,del2,cg,sg,dist0,dist1;
      AngBre(am2, &del1,&del2, &cg,&sg, &dist0,&dist1);
      dist0 = dist0 *m_Event->m_Xenph;
      m_WtMass    =m_WtMass *(dist1/dist0);
      double phi=2*M_PI*m_RNgen->Rndm();
      xphot[i].SetPxPyPzE(xk*sg*cos(phi), xk*sg*sin(phi), xk*cg, xk);
      m_yini[i]    =xk*del1/2;
      m_zini[i]    =xk*del2/2;
   }// for i
///////////////////////////////////////////////////////////////////////////
// Here we determine dilatation factor for rescaling 4-momenta           //
// of all photons such that total 4-momentum is conserved (for fixed v)  //
///////////////////////////////////////////////////////////////////////////
   double DilFac0,DilFac,DilJac, DilJac0;
//   double WtDil0;
   if(m_nphot == 1) {
      DilFac0 = 1;
      DilFac  = 1;
      DilJac  = 1;
   } else {
      pk.SetPxPyPzE(0,0,0,0);
      pp.SetPxPyPzE(0,0,0,2); //total energy in units of ene =2
      for(int i=1; i<=m_nphot; i++) pk=pk+xphot[i];
      double ppdpp = pp*pp;
      double pkdpk = pk*pk;
      double ppdpk = pp*pk;
      double AA    = ppdpp*pkdpk/sqr(ppdpk);
// Dilatation factor
      DilFac0 = 2*ppdpk/ppdpp/vv;
      DilFac  = DilFac0*.5*(1+sqrt(1-vv*AA));
// and the corresponding jacobian factor
      DilJac  = (1+1/sqrt(1-vv*AA))/2;
   }// if nphot
   DilJac0 = (1+1/sqrt(1-vv))/2;  // as in crude of vv-dist.
   m_WtDil  = DilJac/DilJac0;
   m_WtCut  = 1;
// scale down photon energies and momenta
   for(int i=1; i<=m_nphot;i++){
      m_yini[i] =m_yini[i]/DilFac;
      m_zini[i] =m_zini[i]/DilFac;
      m_Event->m_PhotISR[i]=xphot[i]*(1/DilFac);
   }// for i
//     Check on lower energy cut-off
   if(m_Event->m_PhotISR[m_nphot].E() < m_vvmin ) m_WtCut = 0;
} //if vv

//     Photon momenta rescaled into GEV units
for(int  i=1; i<=m_nphot; i++){
	  m_Event->m_PhotISR[i] *= Ene;
      m_Event->m_PhotISR[0] += m_Event->m_PhotISR[i];
}// for i
// 4-momentum left after photon emission
TLorentzVector PXX;
PXX.SetPxPyPzE(0,0,0,XXXene);
PXX= PXX -m_Event->m_PhotISR[0];
//
if(DB->KeyWtm == 1) m_WtMass=1; // old tests
// Total ISR weight
*WtIni = m_WtMass *m_WtDil *m_WtCut;
*PX =PXX;
m_Event->m_PX = PXX;
//----------------------------
}//KKarLud::YFSgen

void  KKarLud::PoissGen(double average, int *mult, double rr[]){
//////////////////////////////////////////////////////////////////////////////
// Last corr. nov. 91                                                       //
// This generates photon multipl. nphot according to poisson distr.         //
// Input:  average = average multiplicity                                   //
//         nmax  = maximum multiplicity                                     //
// Output: mult = generated multiplicity                                    //
//         rr(1:100) list of ordered uniform random numbers,                //
//         a byproduct result, to be eventually used for some further       //
//         purpose (i.e.  generation of photon energies).                   //
//////////////////////////////////////////////////////////////////////////////
double  rn,sum,y;
int     nn;
//------------------------------------------------------------------------------
e50:
nn=0;
sum=0;
for(int it=1; it< maxPhot; it++){
   rn = m_RNgen->Rndm();
   y= log(rn);
   sum=sum+y;
   nn=nn+1;
   rr[nn] = sum/(-average);
   if(sum < -average) goto e130;
}
m_nfail=m_nfail+1;
if(m_nfail > 100)
   {cout<< "KKarLud::poissg: EXIT, too small maxPhot"<<maxPhot<<endl; exit(20);}
goto e50;
e130:
*mult=nn-1;
}// KKarLud::PoissGen

void KKarLud::AngBre(double am2,
		double *ddel1, double *ddel2, double *costhg, double *sinthg, double *dist0, double *dist1){
//////////////////////////////////////////////////////////////////////////////
// This routine generates photon angular distribution                       //
// in the rest frame of the fermion pair.                                   //
// The distribution is the S-factor without mass term,                      //
// i.e. without terms 2p_1p_2/(kp_1)(kp_2)                                  //
// Fermion mass is treated exactly!                                         //
// INPUT:                                                                   //
//     am2 = 4*massf**2/s where massf is fermion mass                       //
//     and s is effective mass squared of the parent fermion-pair.          //
// OUTPUT:                                                                  //
//     del1= 1-beta*cos(theta)                                              //
//     del2= 1+beta*cos(theta)                                              //
//     costhg, sinthg, cos and sin of the photon                            //
//     angle with respect to fermions direction                             //
//     dist0 = distribution generated, without m**2/(kp)**2 terms           //
//     dist1 = distribution with m**2/(kp)**2 terms                         //
//////////////////////////////////////////////////////////////////////////////
double rn[2];
m_RNgen->RndmArray(2,rn);
double beta =sqrt(1-am2);
double eps  =am2/(1+beta);                      //= 1-beta
double del1 =(2-eps)* exp( log(eps/(2-eps))*rn[0] );  //= 1-beta*costhg
double del2 =2-del1;                                  //= 1+beta*costhg
// calculation of sin and cos theta from internal variables
*costhg=(del2-del1)/(2*beta);                   // exact
*sinthg=sqrt( del1*del2-am2*sqr(*costhg) );      // exact
// symmetrization
double a;
if(rn[1] < 0.5) {
  a=del1;
  del1=del2;
  del2=a;
  *costhg= -(*costhg);
}
*dist0=1/(del1*del2)*(1 -am2/2);
*dist1=1/(del1*del2)*(1 -am2/2 -am2/4*(del1/del2+del2/del1));
*ddel1=del1;
*ddel2=del2;
// totaly equivalent formula is the following
//     dist1=1d0/(del1*del2)   *beta*sinthg**2/(del1*del2)
//     WRITE(6,*) 'alt dist1 = ', 1d0/(del1*del2) *beta*sinthg**2/(del1*del2)
}//KKarLud::AngBre

void KKarLud::MomPrint1(TString text, TLorentzVector &Vect){
//void KKarLud::MomPrint1(TLorentzVector &Vect){
//////////////////////////////////////////////////////////////
// printing entire four-vector in one line (withot endline)
//  for ( int k=0; k < 4 ; k++ )   *m_Out<< SW20 << Vect[k];
  cout<<text;
  for ( int k=0; k < 4 ; k++ )   cout<< SW20 << Vect[k]; cout<<endl;
}//MomPrint
