///////////////////////////////////////////////////////////////////////////////
#include "KKarFin.h"

ClassImp(KKarFin);

extern "C" {
//
   void pseumar_initialize_(const int&, const int&, const int&);
   void pseumar_makevec_(float rvec[], const int&);
}//

#define SW20 setw(20)<<setprecision(14)

KKarFin::KKarFin()
{
  // This constructor is for ROOT streamers ONLY
  // all pointers has to be NULLed
  cout<< "----> KKarFin Default Constructor (for ROOT only) "<<endl;
  m_Out   = NULL;
  DB      = NULL;
  m_RNgen = NULL;
  m_Event = NULL;
  m_BVR   = NULL;
}

///_____________________________________________________________
KKarFin::KKarFin(ofstream *OutFile)
{
  cout<< "----> KKarFin USER Constructor "<<endl;
  m_Out   = OutFile;
  DB      = NULL;
  m_RNgen = NULL;
  m_Event = NULL;
  m_BVR   = NULL;
}//KKarFin

///______________________________________________________________________________________
KKarFin::~KKarFin()
{
  //Explicit destructor
  cout<< "----> KKarFin::KKarFin !!!! DESTRUCTOR !!!! "<<endl;
}///destructor

double KKarFin::sqr( const Double_t x ){ return x*x;};

///______________________________________________________________________________________
void KKarFin::Initialize()
{
  cout  << "----> KKarFin::Initialize, Entering "<<endl;
//=================================================================
// BX*** macros are in MCdev/BXFORMAT.h
  BXOPE(*m_Out);
  BXTXT(*m_Out,"========================================");
  BXTXT(*m_Out,"======    KKarFin::Initialize     ======");
  BXTXT(*m_Out,"========================================");

  m_NevGen = 0;
  m_MarTot = 0;

  m_vvmin     = DB->vvmin;
  m_Emin      = DB->CMSene/2 *m_vvmin;
  m_delta     = m_vvmin*DB->delfac;
  m_MasPhot   = 1e-60;  // dummy photon mass
  m_icont   =0;
  ///////////////////////////////////////////////////
}// Initialize


void  KKarFin::Make(TLorentzVector *PX, double *wt){
/////////////////////////////////////////////////////////////////////////////
// INPUT:
//     amfi1,2 = final masses, may change event per event (resonances)
//     CharSq  = charge squared
//     PX     = 4-momentum of the entire FSR system (fermions+photons)
//     KFfin     = final state fermion KF code
// OUTPUT:
//     m_qf1,2   = final state charged particle momenta
//     m_nphot   = FSR photon multiplicity
//     m_phot    = list of FSR photons
//     phsu    = total FSR photon 4-momentum
//     wt      = Crude mc weight
/////////////////////////////////////////////////////////////////////////////
m_icont++;
m_KFfin = m_Event->m_KFfin;
double amfi1  = DB->fmass[m_KFfin];
double amfi2  = amfi1;
double CharSq = sqr(DB->Qf[m_KFfin]);
//////////////////////
//     Generate photons and fermions in the rest frame of ferms Q=qf1+qf2
//     and transform to LAB system,
//     PX=qf1+qf2+phsu is total momentum of the entire FSR system,
//     as seen from the LAB system.
//////////////////////

m_IsFSR = DB->KeyFSR; // general FSR switch
// but exception for neutrinos
if( m_KFfin == 12 || m_KFfin == 12 || m_KFfin == 16) m_IsFSR = 0;
//-----------------------------------------------------------------------
double WtFin;
if(m_IsFSR == 1) {
   m_NevGen = m_NevGen+1;   // Only bremsstrahlung counts!!!
   YFSfin( PX,  amfi1, amfi2,  CharSq,  &WtFin);
   m_Event->m_Qf1 = m_q1;
   m_Event->m_Qf2 = m_q2;
   m_Event->m_nPhotFSR = m_nphot;
   (m_Event->m_PhotFSR)[0].SetPxPyPzE(0.,0.,0.,0.);
   for(int i =1; i<= m_nphot; i++) {
	   m_Event->m_PhotFSR[i] = m_phot[i];
	   (m_Event->m_PhotFSR)[0] += m_phot[i]; // sum in zero component
   }
} else {
//-----------------------------------------------------------------------
// Final state bremss, fermion momenta defined in Z frame if no ISR
// In case of ISR they are reassigned one more time.
    //[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
        float rvec[10];
        pseumar_makevec_(rvec,2);
        double cth= 1 -2*rvec[0];
        double the= acos(cth);
        double phi= 2*M_PI*rvec[1];
        //(*m_Out) <<"@@@ KKarFin::Make: the, phi= "<< the<<"  "<< phi <<endl;
    //]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
   //double cth= 1 -2*m_RNgen->Rndm();
   //double the= acos(cth);
   //double phi= 2*M_PI*m_RNgen->Rndm();
   m_Event->PhaSpac2(PX, the, phi, amfi1, &(m_Event->m_Qf1), &(m_Event->m_Qf2) );
   WtFin = 1;
   m_nphot = 0;
}
//-----------------------------------------------------------------------
//   FSR weight
*wt = WtFin;
//CALL GLK_Mfill(m_idyfs+69, wt,  1d0)
//CALL GLK_Mfill(m_idyfs+60, 1d0*m_nphot,  1d0)
// =============================================
}//!!!KKarFin_Make!!!



void  KKarFin::YFSfin( TLorentzVector *PX, double Mas1, double Mas2, double CharSq, double *WtFin){
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
//  Simulates final state bremsstrahlung out of pair of charged patricles.                 //
//  (algorithm of 11-th March 1989, by S. Jadach)                                          //
//  INPUT  :                                                                               //
//      MltFSR  = normaly set to zero, otherwise enforces photon                           //
//                multiplicity to be exactly equal MltFSR (special tests)                  //
//      KeyPia  = photon removal switch                                                    //
//      PX     = 4-momentum of FSR system as seen from LAB system                          //
//      Mas1,2 = masses of the final charged pair                                          //
//      delta   = lower energy bound (dimensionless)                                       //
//      emin    = lower energy bound (GeV) in LAB system                                   //
//      CharSq  = final state charge squared                                               //
//      alfinv  = 1/alpha_QED                                                              //
//  OUTPUT :                                                                               //
//      qf1,2    = final fermion four momentum (GeV)                                       //
//      nphot   = photon multiplicity                                                      //
//      phot    = photon four momenta (GeV) in cms                                         //
//      phsu    = sum of photon momenta                                                    //
//      ygr,zet = Sudakov variables                                                        //
//      WtFin   = MC weight at this stage                                                  //
//      martot  = control variable, no of marked photons                                   //
//  DEBUG:                                                                                 //
//      Mark    = marks on photons CLOSE to lower energy bound                             //
//      qf1c,2c = ficticious final fermion four momenta for crude S-factor                 //
//      wt1    = the weight - phase space limits for very hard phot.                       //
//      wt2    = the weight - translation jacobian.                                        //
//      wt3    = the weight - single final mass-weight                                     //
//      wtm    = the list/matrix of mass-weights.                                          //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////

double   WtMlist[maxPhot];
double   rr[maxPhot],xk[maxPhot],cgx[maxPhot],sgx[maxPhot];
double   dis0[maxPhot];
int      Mark[maxPhot];
double   rvec[maxPhot];
//-----------------------------------------------------------------------
m_NevGen = m_NevGen + 1;

double alf1 = 1/M_PI/DB->Alfinv0;
double wt1   =1;
double wt2   =1;
double wt3   =1;
//-----------------------------------------------------------------------
double svar = (*PX)*(*PX);
double amfin = min(Mas1,Mas2);
double amc2  = 4*sqr(amfin)/svar;  // overvalued svar for crude
double betc  = sqrt(1-amc2);      // overvalued svar for crude
for(int i=0; i<maxPhot; i++) {
   rr[i]=0;
   Mark[i]=0;
}// for i

// Initialize fermion and photon momenta to zero
m_q1.SetPxPyPzE(0,0,0,0);
m_q2.SetPxPyPzE(0,0,0,0);
m_nphot=0;
for(int i=0; i<maxPhot;i++){
   m_yfin[i]=0;
   m_zfin[i]=0;
   m_phot[i].SetPxPyPzE(0,0,0,0);
}//  for i
//-----------------------------------------------------------------------
//  generate photon multiplicity, average = average multiplicity (crude)
//-----------------------------------------------------------------------
double gamf2 = CharSq*alf1 *(1+sqrt(betc))/betc *log( sqr(1+betc)/amc2);
double average = gamf2*log(1/m_delta);
average = average *DB->Xenph;
//
e60:
PoissGen(average,&m_nphot, rr);
// This is for tests of program at fixed multiplicity (advanc. users)
// if((m_MltFSR != 0) && (m_nphot != m_MltFSR)) goto 5
double sprim;
if (m_nphot == 0) {
   sprim=svar;
} else {
//-----------------------------------------------------------------------
// begin with photon energy
   double xsum=0;
   for(int i=1; i<=m_nphot; i++){
      xk[i]=exp(log(m_delta)*rr[i]);  // f77 indexing in rr and xk
      if(xk[i] < sqrt(10.0)*m_delta) Mark[i]=1;
      xsum = xsum+xk[i];
   }//
// Event outside phase space (too hard photon)
//   if(xsum >= 1.0) goto e60; // BUG
   if(xsum >= 1.0) {
       m_nphot=0; *WtFin =0; return;
   }
   double xfact=1.0/(1.0-xsum);
   for(int i=1; i<=m_nphot; i++){
      xk[i]=xk[i]*xfact;
   }
   m_RNgen->RndmArray(m_nphot,rvec);
   for(int i=1; i<=m_nphot; i++){
//-----------------------------------------------------------------------
//     simplified photon angular distribution,
//     s'->s and m**2/(kp)**2 dropped
//     cg=cos(theta) and sg=sin(theta) memorized to avoid rounding err.
	  double dl1,dl2, cg, sg, dis1;
	  AngBre(amc2, &dl1,&dl2, &cg,&sg, &dis0[i],&dis1);
	 //-----------------------------------------------------------------------
//    define photon momenta (in units of sqrt(s')/2 )
      double phi=2.0*M_PI*rvec[i-1];
      m_phot[i].SetPxPyPzE(sg*cos(phi), sg*sin(phi), cg, 1);
      m_phot[i] *= xk[i];
      m_phot[0] = m_phot[0]+m_phot[i]; // accumulate sum in zero component
      cgx[i]=cg;
      sgx[i]=sg;
   }// for i
//-----------------------------------------------------------------------
//   determine rescaling factor and s', wt2 is dilatation jacobian
   double xmk2 = m_phot[0] *m_phot[0];
   double yy   = 1.0/(1.0  +m_phot[0].E() +xmk2/4.0 );
   wt2  = yy*(1.0+m_phot[0].E());
   sprim= svar*yy;
//-----------------------------------------------------------------------
// reject events with too hard photons
   double smini = sqr(Mas1+Mas2);
   if(sprim < smini)   {
	   m_nphot=0; *WtFin =0; return;
   }

//-----------------------------------------------------------------------
//  Recsale properly all photon momenta
//-----------------------------------------------------------------------
   double ener = sqrt(sprim)/2;
   for(int i=0; i<=m_nphot; i++){
         m_phot[i] *= ener;
   }
}// m_nphot=0
//[[[[[[[[[[[[[[[[[[[[[
//cout<<"##### final  wt2= "<<wt2<<endl;
//for(int i=1; i<=m_nphot; i++) MomPrint1("m_phot[*]=", m_phot[i]);
//]]]]]]]]]]]]]]]]]]]]]
//-----------------------------------------------------------------------
//     final fermion momenta
//-----------------------------------------------------------------------
double amcru  = amfin*sqrt(sprim/svar);
double qmsene = sqrt(sprim);
double ener   = qmsene/2.0;
double betn,eta1,eta2;

// define fermion pair in its rest frame along z-axix
m_Event->GivePair(qmsene, Mas1, Mas2,&m_q1, &m_q2 ,&betn,&eta1,&eta2);   //real
m_Event->GivePair(qmsene,amcru,amcru,&m_r1, &m_r2 ,&betc,&eta1,&eta2);   //real
//[[[[[[[[[[[[[
//[[[[[[[[[[[[[[[[[[[[[
//   MomPrint1("m_q1=", m_q1);
//   MomPrint1("m_q2=", m_q2);
//   MomPrint1("m_r1=", m_r1);
//   MomPrint1("m_r2=", m_r2);
//]]]]]]]]]]]]]]]]]]]]]
//-----------------------------------------------------------------------
//     Mass weight for theta distribution
//-----------------------------------------------------------------------
//     Mass weight compensates for s'->s and droping terms -m**2/(k.q)**2
//     Care is taken of machine rounding errors.
//     del1 and del2 RECALCULATED out of angles sgx(i),cgx(i)
//     with TRUE sprim, using EXACT formulas
double amd1 = sqr(Mas1/ener);
double amd2 = sqr(Mas2/ener);
double del1,del2;
for(int i=1; i<= m_nphot; i++){
   if( cgx[i] > 0.0 ) {
      del1 = amd1/(eta1+betn) +betn*sqr(sgx[i])/(1+cgx[i]);
      del2 = eta2 +betn*cgx[i];
   } else {
      del1 = eta1 -betn*cgx[i];
      del2 = amd2/(eta2+betn) +betn*sqr(sgx[i])/(1-cgx[i]);
   }
   double dist1=1.0/(del1*del2) *(1.0 -(amd1+amd2)/4.0
                    -amd2/4.0*del1/del2 -amd1/4.0*del2/del1);
   WtMlist[i]= dist1/dis0[i];
   if(WtMlist[i] < 1.0e-90) WtMlist[i]= 0.0;
//*****************************************
//** dist1x below is exactly the same as dist1 but for small masses is
//** prone to severe rounding errors (in contrast to dist1 used above)
//         IF((1-sprim/svar) .gt. 0.01d0) THEN
//            qf1qf2= m_q1(4)*m_q2(4) -m_q1(3)*m_q2(3) -m_q1(2)*m_q2(2) -m_q1(1)*m_q2(1)
//            qf1k = m_q1(4)*m_phot(i,4)-m_q1(3)*m_phot(i,3) -m_q1(2)*m_phot(i,2)-m_q1(1)*m_phot(i,1)
//            qf2k = m_q2(4)*m_phot(i,4)-m_q2(3)*m_phot(i,3) -m_q2(2)*m_phot(i,2)-m_q2(1)*m_phot(i,1)
//            dist1x = 2*qf1qf2/qf1k/qf2k -Mas1**2/qf1k**2 -Mas2**2/qf2k**2
//            dist1x = dist1x/4d0*m_phot(i,4)**2
//            WRITE(*,'(a,5f20.10)') '===>: ',dist1x/dist1
//         ENDIF
//*****************************************
//     finaly define Sudakov variables
   m_yfin[i]=del1*xk[i]/2;
   m_zfin[i]=del2*xk[i]/2;
//[[[[[[[[[[[[[[[[[[[[[
//   cout<<"##### final  WtMlist[i]= "<<WtMlist[i]<<endl;
//   cout<<"##### m_yfin[i]="<<m_yfin[i]<<"  m_zfin[i]="<<m_zfin[i]<<endl;
//]]]]]]]]]]]]]]]]]]]]]
}// for i
//-----------------------------------------------------------------------
// Transform from rest frame of Q=qf1+qf2 down to CMS=Lab,
// through intermediate rest frame of PX=qf1+qf2+phsu.
// KarFin_Kinf1(PX,m_q1,m_q2,m_r1,m_r2,m_nphot,m_phot,m_phsu)
//-----------------------------------------------------------------------
// Angles for Euler rotation of all FSR system
   double cth= 1.0 -2*m_RNgen->Rndm();
   double the= acos(cth);
   double phi= 2.0*M_PI*m_RNgen->Rndm();
   // QQk is combined momentum of q1+q2+photons in q1+q2 system
   TLorentzVector QQk= m_q1+m_q2+m_phot[0];
// Transform fermions back to LAB/CMS through intermediate Z frame
   m_Event->BoostEul(the,phi,&QQk,PX,&m_q1);
   m_Event->BoostEul(the,phi,&QQk,PX,&m_q2);
   m_Event->BoostEul(the,phi,&QQk,PX,&m_r1);
   m_Event->BoostEul(the,phi,&QQk,PX,&m_r2);
// Transform photons
   for( int i=0; i<=m_nphot; i++){
     m_Event->BoostEul(the,phi,&QQk,PX,&m_phot[i]);  // m_phot[0] included!
   }// for i
//-----------------------------------------------------------------------
// Calculate YFS formfactor (cut-off dependent part) and mass weights
// Optionally removing photons below emin from the list
Piatek( Mas1, Mas2, CharSq, WtMlist, &wt3);
//-----------------------------------------------------------------------
// Monitoring weights and other control variables,
// Non-essential for the MC generation itself.
//CALL GLK_Mfill(m_idyfs+64, 1d0*m_nphot  ,1d0)
//      uu = 1d0-sprim/svar
//      CALL GLK_Fil1(m_idyfs+31, uu  ,wctrl)
//      CALL GLK_Fil1(m_idyfs+32, uu  ,  1d0)
// marked photons
if(m_nphot >= 1) {  // spurious if ???
   for(int i=1; i<=m_nphot; i++) {
//            ul= log10(m_phot(i,4)/m_Emin)
//            IF(Mark(i) .EQ. 1)   CALL GLK_Fil1(m_idyfs+20,   ul,1.d0)
      if(Mark[i] == 1) m_MarTot =m_MarTot+1;
   }
}//if
// Final Monitoring weights
*WtFin = wt1*wt2*wt3;
//[[[[[[[[[[[[[[[[[[[[[[[[[
//if( m_icont <50 ) {
//	cout<<" WtFin ="<<*WtFin<<"  wt1="<<wt1<<"  wt2="<<wt2<<"  wt3="<<wt3<<endl;
//}//if( m_icont
//]]]]]]]]]]]]]]]]]]]]]]]]]
//-----------------------------------------------------------------------
//CALL GLK_Mfill(m_idyfs+61,     wt1  ,1d0)
//CALL GLK_Mfill(m_idyfs+62,     wt2  ,1d0)
//CALL GLK_Mfill(m_idyfs+63,     wt3  ,1d0)
//-----------------------------------------------------------------------
}// YFSfin


void  KKarFin::Piatek(double Mas1, double Mas2, double CharSq, double WtMlist[], double *Wt3){
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
// Written CERN, piatek, 22 sept. 1989  (S.J.)                                             //
// This routine calculates YFS form-factor and optionaly                                   //
// removes some of soft photons (below Emin in CMS), calculating compensating weight.      //
// Note the action of this routine is not Loretnz invariant !!!!                           //
// Input:                                                                                  //
//     KeyPia   = 0, NO removal of photons below Emin                                      //
//              = 1, with removal of photons below Emin                                    //
//     Mas1,2  = fermion masses                                                           //
//     delta    = infrared cut-off in generation (dimensionless)                           //
//     CharSq   = FSR charge squared                                                       //
//     alfinv   = 1/alpha_QED                                                              //
//     qf1,2    = fermion momenta                                                          //
//     qf1c,2c  = ghost fermion momenta in crude S-factor                                  //
//     phsu     = sum of photon momenta                                                    //
//     phot     = list of photon momenta                                                   //
//     WtMlist  = list of mass-weights for all photons                                     //
// OUTPUT:                                                                                 //
//     WtRem    = mass-weight of removed photons, for tests                                //
//              = 1,  for KeyPia=0                                                         //
//     WtMas    = mass-weight for all photons, see comments below,                         //
//                the MOST IMPORTANT OUTPUT of Piatek !!!                                  //
//     WtCtrl    = control-weight for remowed photons, for tests                           //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
// REMARK on r1, r2:                                                                       //
// They are ghost-fermion momenta in the definition of the crude                           //
// distribution, i.e. truncated S-factor is r1*r2/(k*r1)/(k*r2)                            //
// as generated in the crude MC.                                                           //
// In QMS frame 4-momenta r1, r2 have the same energy                                      //
// and directions as q1, q2                                                                //
// but different masses Mas1c, Mas2c and therefore longer 3-momenta.                       //
// They are usefull because we may use the same analytical                                 //
// expression for the exact Breal and  for 'crude Breal'                                   //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////
//-----------------------------
double WtCtrl =1.0;
double WtRem  =1.0;
double alfpi = 1/M_PI/DB->Alfinv0;
double alfch = alfpi*CharSq;
double alfCR = alfch* DB->Xenph;
double Mass = min(Mas1,Mas2);
TLorentzVector  PX = m_q1+m_q2 +m_phot[0];
TLorentzVector  QQ = m_q1+m_q2;
double svarX = PX*PX;
double svarQ = QQ*QQ;
double EQQ   = 0.5*sqrt(svarQ);
double q1q2  =  m_q1*m_q2;
double r1r2  =  m_r1*m_r2;
double QQk   =  QQ*m_phot[0];
double Eq1   = (svarQ +Mas1*Mas1 -Mas2*Mas2)/(2*sqrt(svarQ));
double Eq2   = (svarQ -Mas1*Mas1 +Mas2*Mas2)/(2*sqrt(svarQ));
double uu= 1-svarQ/svarX;
// Delta1 and Epsi1 are cutoffs located in YFS formfactor
// Note that xfact=(1+2*QQk/svarQ) is EXACTLY the same as in KarFin_YFSfin
// Delta1 = 2*Emin/sqrt(s'), where Emin is in QMS
double Delta1 = m_delta*(1+ 2*QQk/svarQ);
double Epsi1  = sqrt(m_Emin*m_Emin/m_q1.E()/m_q2.E());
double EminQ  = EQQ*Delta1;
// The total phase space integral for crude x-section and YFS formfactor cut-off dependend part
// Note that delta is a lower limit on y-z-variables in crude generation
double amc2  =  4*Mass*Mass/svarX;
double betc  =  sqrt(1-amc2);
//===========================================================================================
double VoluMC =  2.0*alfCR *(1.0+betc*betc)/(2*betc)* log(sqr(1+betc)/amc2) *log(1/m_delta);
double YFS_IR = -2.0*alfch *(     q1q2 *m_BVR->A( q1q2, Mas1,Mas2) -1.0  )    *log(1/Delta1);


// YFSkon is delegated/exported to QED3 (unused here).
double YFSkon =  1/4.0 *2*alfch*(log(svarQ/sqr(Mass))-1) + alfch*( -0.5  + sqr(M_PI)/3.0); // Mass<<sqrt(s)
//[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
//if(m_icont <50 && uu<0.01){
//cout<<" *******************************************************************************nphot="<<m_nphot<<" uu="<<uu<<endl;
//cout<<" KKarFin::Piatek VoluMC= "<<VoluMC<<" YFS_IR="<<YFS_IR<<" YFSkon="<<YFSkon<<endl;
//cout<<" alfch="<< alfch <<" q1q2="<<q1q2<<" Mas1="<<Mas1<<"  Delta1="<<Delta1<<endl;
//cout<<"  m_delta="<<m_delta<<"  QQk="<< QQk << " svarQ="<<svarQ<<" svarQ/svarX="<< svarQ/svarX <<endl;
//}//if(m_icont
//]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
//===========================================================================================
// Corrections necessary for photon remooval scenario (now default!)
double Mas1c  = Mass*sqrt(svarQ/svarX);
double Mas2c  = Mass*sqrt(svarQ/svarX);
double BtiXcru= m_BVR->Btildc(alfCR, r1r2, m_r1.E(),m_r2.E(), Mas1c, Mas2c, m_Emin, m_MasPhot); //crude
double BtiQcru= m_BVR->Btildc(alfCR, r1r2, EQQ,     EQQ,      Mas1c, Mas2c, EminQ,  m_MasPhot); //crude
double BtiXexa= m_BVR->Btilda(alfch, q1q2, m_q1.E(),m_q2.E(), Mas1,  Mas2,  m_Emin, m_MasPhot); //exact
double BtiQexa= m_BVR->Btilda(alfch, q1q2, Eq1,     Eq2,      Mas1 , Mas2,  EminQ,  m_MasPhot); //exact

double DelVol = BtiXcru -BtiQcru;   // positive
double DelYFS = BtiXexa -BtiQexa;   // positive
double DelB2  = -DelVol+DelYFS;
//------ oldies ------
// Total QMS-CMS, result should be negative ( delta<<epsilon )
double DelB  =   VoluMC +YFS_IR;
// Ultrarelativistic (small mass) old version is the following:
double DelB2u = -2*alfch*(log(svarX/svarQ)+1) *log(Epsi1/Delta1);
// The average mass-weight for removed photon = exp(DelB2)
// It can be calculated analyticaly as a  ratio of YFS formfactors
// On the other hand, it is checked by MC, see control weight WtCtrl.
//
//*********************************************************************************************
//*      IF(iCont.LE.10 ) THEN
//*         IF((1-svarQ/svarX) .GT. 0.1d0) THEN
//*            iCont = iCont+1
//*((( Approximate version of DelB2 without contr. with A4 terms for tests (surpisingly good!!!)
//*      DelB2w = -2*alfch*(  (1d0+betc**2)/(2d0*betc)* DLOG((1d0+betc)**2/amc2)
//*     $                    -( q1q2*BVR_A(q1q2,Mas1,Mas2) -1d0 )
//*     $                  ) *DLOG(Epsi1/Delta1)
//*         WRITE(*,'(a,5f20.10)') 'piatek: ',1-svarQ/svarX,DelB2,DelB2u/DelB2,DelB2w/DelB2
//*)))
//*            VoluMC2 =
//*     $            BVR_Btildc(alfch, r1r2, EQQ,EQQ,  Mas1c,Mas2c, EQQ,          m_MasPhot)
//*     $           -BVR_Btildc(alfch, r1r2, EQQ,EQQ,  Mas1c,Mas2c, EQQ*m_Delta,  m_MasPhot)
//*            WRITE(*,'(a,5f20.10)') '###Piatek: ',1-svarQ/svarX, VoluMC2,VoluMC,VoluMC2/VoluMC
//*         ENDIF
//*      ENDIF
//*********************************************************************************************
//--------------------------
double wtm1=1.0;
double wtm2=1.0;
// mass weight below and above Emin calculated separately
for(int i=1; i<=m_nphot; i++){
   if(m_phot[i].E() < m_Emin) {
      wtm1=wtm1*WtMlist[i] /DB->Xenph;
      if(wtm1 <= 1e-90) wtm1=0;
   } else {
      wtm2=wtm2*WtMlist[i] /DB->Xenph;
      if(wtm2 <= 1e-90) wtm2=0;
   }//if
}// for i
//------------------------------------------------------------------------------
// Control weight - its average should be =1 within statist. error!
WtCtrl =wtm1*exp(-DelB2);
if(DB->KeyPia == 0) {
   if( abs(DelB) > 100.0 ) cout<< "#### KarFin_Piatek: DelB= "<<DelB<<endl;
   WtRem    = 1.0;
   m_WtMass = wtm1*wtm2;      //!!! <--removal OFF
   m_VoluMC = exp( VoluMC);
   m_YFS_IR = exp( YFS_IR);
   m_YFSkon = exp( YFSkon);   //!!! <--finite part of YFS formfactor
// =====================================
} else {
// Optional removal of photons below Emin from the record
// in such a case WtMas includes exp(belb2)= <wt3> for removed ph.
//[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
//   cout<<"KKarFin::Piatek:  m_NevGen ="<<m_NevGen<<endl;
//   cout<<"//////KKarFin::Piatek: m_nphot="<<m_nphot<<"  m_Emin="<< m_Emin<<endl;
//   for(int i=1; i<=m_nphot;i++) { cout<<"i="<<i<<"  "; MomPrint1("m_phot[*] =", m_phot[i]);}
//]]]]]]]]
   int nph=m_nphot;
   for(int j=m_nphot; j>=1; j--){
//   cout<<"///// j="<<j<<" E[j]="<< m_phot[j].E()<<endl;
      if(m_phot[j].E() < m_Emin) {
    	 // WARNING: loop below is never executed for ordered energies??
//    	 if( j < nph)
            for(int i=j+1; i<=nph; i++){
//            cout<<"\\ YES! i="<<i<<endl;
               m_phot[i-1]=m_phot[i];
            }// for i
         nph=nph-1;
      }//if
   }//for j
// Correction of Alex Read, probably obsolete because KarFin_Merge is also corrected
   for(int j=nph+1; j<=m_nphot; j++) {
         m_phot[j].SetPxPyPzE(0,0,0,0);
   }
   m_nphot=nph;
//[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
//  for(int i=1; i<=m_nphot;i++) { cout<<"i="<<i<<"  "; MomPrint1("m_phot[*] =", m_phot[i]);}
//]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
// Important to remember: EXP(YFS_IR +YFSkon)         = full YFS formfactor in QMS
//                        EXP(YFS_IR +YFSkon +DelYFS) = full YFS formfactor in CMS
   WtRem    = wtm1;
   m_WtMass = wtm2;           //!!! <--removal ON
   m_VoluMC = exp(VoluMC -DelVol);
   m_YFS_IR = exp(YFS_IR +DelYFS);
// YFSkon is delegated/exported to QED3 (unused here).
   m_YFSkon = exp(YFSkon);    //!!! <--finite part of YFS formfactor in QMS frame!!!
//===============================================
}// if KeyPia
m_Event->m_YFS_IR_fin = m_YFS_IR;  // to be used in KKceex3
m_Event->m_YFSkon_fin = m_YFSkon;
//**********************************
*Wt3 = m_WtMass *m_VoluMC *m_YFS_IR;
//**********************************
//[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
//if(m_icont <50 && uu<0.01){
//cout<<" *******************************************************************************nphot="<<m_nphot<<" uu="<<uu<<endl;
//cout<<" KKarFin::Piatek VoluMC= "<<m_VoluMC<<" YFS_IR="<<m_YFS_IR<<" YFSkon="<<m_YFSkon<<" m_WtMass="<<m_WtMass<<endl;
//if(m_icont <100){
//cout<<" KKarFin::Piatek *************** m_YFS_IR= "<< m_YFS_IR <<endl;
//cout<<" KKarFin::Piatek m_WtMass= "<<m_WtMass<<"   Wt3="<< *Wt3  <<endl;
//cout<<" KKarFin::Piatek   WtCtrl= "<<WtCtrl<<  " WtRem="<< WtRem<<endl;
//}//if(m_icont <50
// Monitoring
//CALL GLK_Mfill(m_idyfs+65,       WtCtrl ,1d0)
//CALL GLK_Mfill(m_idyfs+66,       WtRem  ,1d0)
//
}//!!! KarFin_Piatek !!!

void  KKarFin::PoissGen(double average, int *mult, double rr[]){
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
   {cout<< "KKarFin::poissg: EXIT, too small maxPhot"<<maxPhot<<endl; exit(20);}
goto e50;
e130:
*mult=nn-1;
}// KKarFin::PoissGen

void KKarFin::AngBre(double am2,
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
double rn[10];
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
}//KKarFin::AngBre


void KKarFin::MomPrint1(TString text, TLorentzVector &Vect){
//void KKarLud::MomPrint1(TLorentzVector &Vect){
//////////////////////////////////////////////////////////////
// printing entire four-vector in one line (withot endline)
//  for ( int k=0; k < 4 ; k++ )   *m_Out<< SW20 << Vect[k];
  cout<<text;
  for ( int k=0; k < 4 ; k++ )   cout<< SW20 << Vect[k]; cout<<endl;
}//MomPrint
