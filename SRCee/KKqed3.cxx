///////////////////////////////////////////////////////////////////////////////
#include "KKqed3.h"

ClassImp(KKqed3);


KKqed3::KKqed3()
{
  // This constructor is for ROOT streamers ONLY
  // all pointers has to be NULLed
  cout<< "----> KKqed3 Default Constructor (for ROOT only) "<<endl;
  m_Out   = NULL;
  DB      = NULL;
  m_Event = NULL;
  m_BornV = NULL;
}

///_____________________________________________________________
KKqed3::KKqed3(ofstream *OutFile)
{
  cout<< "----> KKqed3 USER Constructor "<<endl;
  m_Out = OutFile;
  DB      = NULL;
  m_Event = NULL;
  m_BornV = NULL;
}//KKqed3

///______________________________________________________________________________________
KKqed3::~KKqed3()
{
  //Explicit destructor
  cout<< "----> KKqed3::KKqed3 !!!! DESTRUCTOR !!!! "<<endl;
}///destructor

//inline double KKqed3::sqr( double x ){ return x*x;};


///______________________________________________________________________________________
void KKqed3::Initialize()
{
  cout  << "----> KKqed3::Initialize, Entering "<<endl;
//=================================================================
// BX*** macros are in MCdev/BXFORMAT.h
  m_KeyOrd = 1;
  BXOPE(*m_Out);
  BXTXT(*m_Out,"========================================");
  BXTXT(*m_Out,"======    KKqed3::Initialize      ======");
  BXTXT(*m_Out,"========================================");
  // m_Event has to be allocated before Initialize
  m_vlim1= 1.0e-9;
  m_vlim2= 1.0e-9;
  m_icont = 0;
  ///////////////////////////////////////////////////
}// Initialize

void KKqed3::Make(){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//   main routine for calculation of long list of the weights                      //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
  m_icont++;

  for(int i=0; i< maxWT; i++) m_WtSet[i] = 0;
  int KFini = m_Event->m_KFini;
  int KFfin = m_Event->m_KFfin;
  double amini = DB->fmass[KFini];
  double amfin = DB->fmass[KFfin];
  double chain2 = sqr(DB->Qf[KFini] );
  double chafi2 = sqr(DB->Qf[KFfin] );
//
  int IsFSR = DB->KeyFSR; // general FSR switch
  // but exception for neutrinos
  if( KFfin == 12 || KFfin == 12 || KFfin == 16) IsFSR = 0;
//
  TLorentzVector pp = m_Event->m_Pf1 + m_Event->m_Pf2;
  TLorentzVector qq = m_Event->m_Qf1 + m_Event->m_Qf2;
  TLorentzVector xx = qq + m_Event->m_PhotFSR[0];  // m_PhotFSR[0] is sum of all
//
  double YFSkon_ini, YFSkon_fin;
  double svar  = pp*pp;
  double svar1 = xx*xx;
  double svar2 = qq*qq;
  double vv = 1 -svar1/svar;
  double uu = 1 -svar2/svar1;
  double gami = chain2* 2/DB->Alfinv0/M_PI*(log(svar/amini/amini)  -1); // Beware! Mass<<sqrt(s) version!
  double gamf = chafi2* 2/DB->Alfinv0/M_PI*(log(svar2/amfin/amfin) -1);
  // ISR
  double alfi   = chain2 /DB->Alfinv0 /M_PI;
  double alfch  = chafi2* 1/DB->Alfinv0/M_PI;
///////////////////////////////
//  YFSkon_ini = exp( 0.25*gami + alfi*(-0.5 + sqr(M_PI)/3.0) );
//  YFSkon_fin = exp( 1/4.0 *2*  alfch *(log(svar2/sqr(amfin))-1) + alfch*( -0.5  +sqr(M_PI)/3.0) ); // Mass<<sqrt(s)
///////////////////////////////
// Imported from MC generation KKarfin and Density, safer method
  YFSkon_ini = m_Event->m_YFSkon_ini;
  YFSkon_fin = m_Event->m_YFSkon_fin;
//  if(m_icont <100) cout<<" YFSkon_fin="<<YFSkon_fin<<endl;

  int nphox =  m_Event->m_nPhotISR;
  int nphoy =  m_Event->m_nPhotFSR;

  for(int j=1;j<=nphox;j++){
     m_yini[j] = 2* m_Event->m_PhotISR[j] *m_Event->m_Pf1 /svar;
     m_zini[j] = 2* m_Event->m_PhotISR[j] *m_Event->m_Pf2 /svar;
  }
  for(int j=1;j<=nphoy;j++){
     m_yfin[j] = 2* m_Event->m_PhotFSR[j] *m_Event->m_Qf1 /svar2;
     m_zfin[j] = 2* m_Event->m_PhotFSR[j] *m_Event->m_Qf2 /svar2;
  }
//////
  double DisCru;   // DisCru has to be the same as in FOAM!!!
  if( DB->KeyThe == 0) {
     DisCru = 4.0/3.0* m_BornV->BornSimple(KFini,KFfin,svar1,0.0);
  } else {
    double CosTheta = m_Event->m_CosTheta;
    DisCru = m_BornV->BornSimple(KFini,KFfin,svar1,CosTheta);
  }//KeyThe
//
  double delp=  sqr(amini)/svar;
  double delq=  sqr(amfin)/svar2;
//
  double cth11,cth12,cth21,cth22;
  m_Event->ThetaR( &cth11, &cth12, &cth21, &cth22);
  double andi11= m_BornV->Born_DizetS(  KFini,KFfin, svar1, cth11);
  double andi12= m_BornV->Born_DizetS(  KFini,KFfin, svar1, cth12);
  double andi21= m_BornV->Born_DizetS(  KFini,KFfin, svar1, cth21);
  double andi22= m_BornV->Born_DizetS(  KFini,KFfin, svar1, cth22);

  double deli1=0, deli2=0, deli3=0;
  if( DB->KeyISR == 1) bvirt0( DB->Alfinv0, chain2,  svar , amini, &deli1, &deli2, &deli3);
  double delf1=0, delf2=0, delf3=0;
  if(IsFSR == 1)       bvirt0( DB->Alfinv0, chafi2, svar2,  amfin, &delf1, &delf2, &delf3);

/////////////////////////////////////////////////////////////////////////////
  double gi1,gi2,ggi1,ggi2,gggi1,gggi2;
  double gf1,gf2,ggf1,ggf2;
  double hfac1,hfac2,hfac3;
  double y,z,sfacj,hfac,yy,zz;
  double y1,z1,y2,z2,y3,z3;
  double dist10, dist11, dist12, dist20, dist21, dist30;
  double bety12, betx12;
  double ntree;
//---------------------------------------------------------
 /////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------
// Beta0, initial+final, factorized form
  //------------------------------------------------------------
  double andis = (andi11 +andi12 +andi21 +andi22)/4;
//[[[
// Collins-Soper Angle calculated from final fermions (ATLAS definition)
//  double Mll  = (m_Event->m_Qf1 + m_Event->m_Qf2).Mag();
//  double PTll = (m_Event->m_Qf1 + m_Event->m_Qf2).Perp();
//  double CosThetaCS = (m_Event->m_Qf1).Plus()*(m_Event->m_Qf2).Minus() - (m_Event->m_Qf1).Minus()*(m_Event->m_Qf2).Plus();
//  CosThetaCS = CosThetaCS/Mll/sqrt(Mll*Mll + PTll*PTll);
//  //m_Event->ThetaR( &cth11, &cth12, &cth21, &cth22);
//  double cosThe = (cth11 + cth12 + cth21 + cth22)/4;  // for debugging
//  //andis =  m_BornV->BornSimple(  KFini,KFfin, svar1, cosThe);
//  andis =  m_BornV->BornSimple(  KFini,KFfin, svar1, CosThetaCS);
//]]]

// KKMC-hh SY note: Beta03 is missing delf3, presumably not fully tested yet.
  m_Beta03 = andis*(1+deli1+deli2+deli3)*(1+delf1+delf2); //O(alf3)
  m_Beta02 = andis*(1+deli1+deli2)      *(1+delf1+delf2); //O(alf2)
  m_Beta01 = andis*(1+deli1)            *(1+delf1);       //O(alf1)
  m_Beta00 = andis;                                       //O(alf0)
// Initial only
  m_beti03 = andis*(1+deli1+deli2+deli3);
  m_beti02 = andis*(1+deli1+deli2);
  m_beti01 = andis*(1+deli1);
  m_beti00 = andis;
// Auxiliary
  double betf01 = andis*(1+delf1);

  //---------------------------------------------------------
  //                beta1 initial
  //---------------------------------------------------------
  m_xBet10 = 0.0;
  m_xBet11 = 0.0;
  m_xBet12 = 0.0;
  m_sbti10 = 0.0;
  m_sbti11 = 0.0;
  m_sbti12 = 0.0;
  if(DB->KeyISR != 0  &&  vv > m_vlim1) {
     for(int jph=1; jph<=nphox; jph++){
        y = m_yini[jph];
        z = m_zini[jph];
        sfacj  =  2/(y*z) *wm0(delp,y,z);
        m_xSfac[jph] = sfacj;
        hfac = wmd(delp,y,z) *sfacj;
        gf1 = 0.5;
        gf2 = 0.5;
        if( m_KeyOrd == 0) {
            Disr1( gami,jph, &gi1,&gi2,&ggi1,&ggi2,&gggi1,&gggi2);
        } else {
            Disr1a(gami,jph, &gi1,&gi2,&ggi1,&ggi2,&gggi1,&gggi2);
        }// if KeyOrd
//-- O(alf1) ----,  tree_level --------
//    The unconventional (1+delf1) in betx10 helps ISR*FSR factorization
//    in the O(alf2) semi-analytical x-check
        dist10= (  gi1*gf1*andi11   +gi1*gf2*andi12
                  +gi2*gf1*andi21   +gi2*gf2*andi22)*hfac;
        m_betx10[jph]=(dist10 -m_Beta00*sfacj )*(1+delf1);
        m_xBet10 = m_xBet10 +m_betx10[jph] /sfacj;
//-- O(alf2) ----,  one_loop_level --------
        dist11= ( ggi1*gf1*andi11  +ggi1*gf2*andi12
           +ggi2*gf1*andi21  +ggi2*gf2*andi22)*hfac;
        m_betx11[jph]= dist11*(1+delf1)       -m_Beta01*sfacj;
        m_xBet11 = m_xBet11 +m_betx11[jph] /sfacj;
//-- O(alf3) ----,  two_loop_level -------- !!!NEW!!!
        dist12= (gggi1*gf1*andi11 +gggi1*gf2*andi12
                +gggi2*gf1*andi21 +gggi2*gf2*andi22)*hfac;
        betx12     = dist12*(1+delf1+delf2) -m_Beta02*sfacj;
        m_xBet12 = m_xBet12 +betx12 /sfacj;
// ***** pure ISR
        m_beti10[jph] =  dist10 -m_beti00*sfacj;   //O(alf1)
        m_sbti10      = m_sbti10 +m_beti10[jph] /sfacj;
        m_beti11[jph] =  dist11 -m_beti01*sfacj;   //O(alf2)
        m_sbti11      = m_sbti11 +m_beti11[jph] /sfacj;
        m_beti12      =  dist12 -m_beti02*sfacj;   //O(alf3) !!!NEW
        m_sbti12 = m_sbti12 +m_beti12      /sfacj;
     }// for jph
  } else {
     for(int jph=1; jph<=nphox; jph++){
        m_xSfac[jph]= -1.0;
        m_betx10[jph] =  0.0;
        m_betx11[jph] =  0.0;
        m_beti10[jph] =  0.0;
        m_beti11[jph] =  0.0;
     }
  }// if vlim
//---------------------------------------------------------
//               beta1 final
//---------------------------------------------------------
  m_yBet10=0.0;
  m_yBet11=0.0;
  m_yBet12=0.0;
  if(IsFSR != 0  &&  uu > m_vlim1) {
     for(int jph=1; jph<=nphoy; jph++){
        yy = m_yfin[jph];
        zz = m_zfin[jph];
        y  = yy/(1 +yy+zz);
        z  = zz/(1 +yy+zz);
        sfacj  =  2/(yy*zz)*wm0(delq,yy,zz);
        m_ySfac[jph] = sfacj;
        hfac = wmd(delq,y,z) *sfacj;
        gi1 = 0.5;
        gi2 = 0.5;
        Dfsr1(gamf, jph, &gf1, &gf2, &ggf1, &ggf2);
//-- O(alf1) ---,  tree level
//    unconventional (1+deli1) in bety10 helps ISR*FSR factorization
//    in the O(alf2) semi-analytical x-check
        dist10= (gi1*gf1*andi11   + gi1*gf2*andi12
                +gi2*gf1*andi21   + gi2*gf2*andi22)*hfac;
        m_bety10[jph] = (dist10 -m_Beta00*sfacj  )*(1+deli1);
        m_yBet10 = m_yBet10 +m_bety10[jph] /sfacj;
//-- O(alf2) ---, one loop level
        dist11= (gi1*ggf1*andi11 + gi1*ggf2*andi12
                +gi2*ggf1*andi21 + gi2*ggf2*andi22)*hfac;
        m_bety11[jph] =  dist11*(1+deli1) -m_Beta01*sfacj;
        m_yBet11 = m_yBet11 +m_bety11[jph] /sfacj;
//-- O(alf3) ---, two loop level for ISR !!!NEW
//   Additional O(alf2) ISR virtual correction deli2 only
        m_bety12 = (1+deli1+deli2)*(dist11 -m_Beta00*(1+delf1)*sfacj);
        m_yBet12 = m_yBet12 +m_bety12 /sfacj;
//*****  pure FSR *****, for construction of beta2, beta3
        m_betf10[jph] =  dist10 -m_Beta00*sfacj;
        m_betf11[jph] =  dist11 -betf01*sfacj;
     }//for jph
        } else {
   for(int jph=1; jph<=nphoy; jph++){
       m_ySfac[jph]  = -1.0;
       m_bety10[jph] =  0.0;
       m_bety11[jph] =  0.0;
   }// for jph
  }//if
//-----------------------------------------------------------
//-----------------------------------------------------------
//                beta2 initial-initial
//-----------------------------------------------------------
  m_xxBet20=0.0;
  m_xxBet21=0.0;
  m_sbti20=0.0;
  m_sbti21=0.0;
  if(DB->KeyISR != 0  &&  vv > m_vlim2) {
     for(int jph2=2; jph2<=nphox; jph2++) {
        for(int jph1=1;  jph1<=jph2-1; jph1++) {
           hfac1  =  m_xSfac[jph1]*wmd(delp,m_yini[jph1],m_zini[jph1]);
           hfac2  =  m_xSfac[jph2]*wmd(delp,m_yini[jph2],m_zini[jph2]);
        //     Summation over two LL fragmentation trees for 2 ISR ohotons,
        //     photon jph1 is always harder because of low level M.C. generation
           ntree = 2;  // for 2 ISR fragmentation trees
           gi1  = 0;
           gi2  = 0;
           ggi1 = 0;
           ggi2 = 0;
           if(m_KeyOrd == 0 ) {
              Disr2(gami,jph1, jph1,jph2, &gi1,&gi2,&ggi1,&ggi2); // 1-st tree
              Disr2(gami,jph1, jph2,jph1, &gi1,&gi2,&ggi1,&ggi2); // 2-nd tree
           } else {
              if( m_yini[jph1]*m_zini[jph1] > m_yini[jph2]*m_zini[jph2]) {
                 Disr2a(gami, jph1,jph2, &gi1,&gi2,&ggi1,&ggi2); // 1-st tree
              } else {
                 Disr2a(gami, jph2,jph1, &gi1,&gi2,&ggi1,&ggi2); // 2-nd tree
              }
           }
           gf1 = 0.5;      // 0.5d0 for averaging over 2 choices
           gf2 = 0.5;      // of Born angles in final state
        //---  O(alf2) ---, tree level---------
           dist20= (gi1*gf1*andi11 + gi1*gf2*andi12
                   +gi2*gf1*andi21 + gi2*gf2*andi22)*hfac1*hfac2/ntree;
        // In beta20 we use beti10 instead of betx10,
        // Reason: unusual definition of betx10, see relevant comment above
           m_beta20 = dist20
               -m_Beta00*m_xSfac[jph1]*m_xSfac[jph2]
               -m_beti10[jph1]*m_xSfac[jph2] -m_beti10[jph2]*m_xSfac[jph1];
           m_betxx20[jph1][jph2]=m_beta20;
           m_xxBet20=m_xxBet20 +m_beta20 /m_xSfac[jph1]/m_xSfac[jph2];
        //---  O(alf3) ---, one loop level ---------!!!!NEW!!!!
           dist21= (ggi1*gf1*andi11+ ggi1*gf2*andi12
                   +ggi2*gf1*andi21+ ggi2*gf2*andi22)
                   *hfac1*hfac2/ntree;
           m_beta21 = dist21*(1+delf1)
                -m_Beta01*m_xSfac[jph1]*m_xSfac[jph2]
                -m_betx11[jph1]*m_xSfac[jph2] -m_betx11[jph2]*m_xSfac[jph1];
           m_xxBet21=m_xxBet21 +m_beta21 /m_xSfac[jph1]/m_xSfac[jph2];
        //***** Pure ISR *****
           m_sbti20=m_sbti20 +m_beta20 /m_xSfac[jph1]/m_xSfac[jph2];
           m_beti21 = dist21
                -m_beti01*m_xSfac[jph1]*m_xSfac[jph2]
                -m_beti11[jph1]*m_xSfac[jph2] -m_beti11[jph2]*m_xSfac[jph1];
           m_sbti21=m_sbti21 +m_beti21 /m_xSfac[jph1]/m_xSfac[jph2];
        }//for jph1
     }//for jph2
  } else {
     for(int  jph2=2; jph2<= nphox; jph2++){
        for(int  jph1=1; jph1<= nphox; jph1++) {
           m_betxx20[jph1][jph2]= 0;
        }//for
     }//for
  }// if vlim
//-----------------------------------------------------------
//                beta2 final-final
//-----------------------------------------------------------
  m_yyBet20=0;
  m_yyBet21=0;
  if(IsFSR != 0  &&  uu > m_vlim2) {
     for(int  jph2=2; jph2<=nphoy; jph2++){
        for(int jph1=1; jph1<= jph2-1; jph1++) {
           y1  = m_yfin[jph1] /(1 +m_yfin[jph1] +m_zfin[jph1] );
           z1  = m_zfin[jph1] /(1 +m_yfin[jph1] +m_zfin[jph1] );
           y2  = m_yfin[jph2] /(1 +m_yfin[jph2] +m_zfin[jph2] );
           z2  = m_zfin[jph2] /(1 +m_yfin[jph2] +m_zfin[jph2] );
        //     Note y1,z1<1 (yfin,zfin cant be used directly in wmd(...))
           hfac1  =  m_ySfac[jph1]*wmd(delq,y1,z1);
           hfac2  =  m_ySfac[jph2]*wmd(delq,y2,z2);
        //     sum over two FSR fragmentation trees
           ntree = 2;  // for 2 FSR fragmentation trees
           gf1 = 0;
           gf2 = 0;
           Dfsr2(jph1,jph1,jph2,&gf1,&gf2); // 1-st tree
           Dfsr2(jph1,jph2,jph1,&gf1,&gf2); // 2-nd tree
           gi1 = 0.5;        // 0.5d0 for averaging over 2 choices
           gi2 = 0.5;        // of Born angles in initial state
        //---- O(alf2) ----, tree level
           dist20 =
                 (gi1*gf1*andi11+ gi1*gf2*andi12
                 +gi2*gf1*andi21+ gi2*gf2*andi22) *hfac1*hfac2/ntree;
           m_beta20 = dist20
                -m_Beta00*m_ySfac[jph1]*m_ySfac[jph2]
                -m_betf10[jph1]*m_ySfac[jph2] -m_betf10[jph2]*m_ySfac[jph1];
           m_betyy20[jph1][jph2]=m_beta20;
           m_yyBet20=m_yyBet20 +m_beta20 /m_ySfac[jph1]/m_ySfac[jph2];
        //---- O(alf3) ----, one loop level  !!!NEW
        // Primitive ISR virtual correction only
           m_beta21 = (1+deli1)*m_beta20;
           m_yyBet21=m_yyBet21 +m_beta21 /m_ySfac[jph1]/m_ySfac[jph2];
        }//for jph1
     }// for jph2
  } else {
     for(int  jph2=2; jph2<=nphoy; jph2++) {
        for(int jph1=1; jph1<= jph2-1;jph1++ ) {
           m_betyy20[jph1][jph2]= 0;
        }//for jph1
     }// for jph2
  }// if vlim
//-----------------------------------------------------------
//                beta2 initial-final
//-----------------------------------------------------------
// or in other terminology   beta1_init - beta1_final
  m_xyBet20=0;
  m_xyBet21=0;
  if(DB->KeyISR*IsFSR != 0 && vv > m_vlim1 && uu > m_vlim1) {
     for(int jph1=1; jph1 <= nphox; jph1++) {
        for(int jph2=1; jph2<= nphoy; jph2++) {
           y1  = m_yini[jph1];
           z1  = m_zini[jph1];
           y2  = m_yfin[jph2] /(1 +m_yfin[jph2] +m_zfin[jph2] );
           z2  = m_zfin[jph2] /(1 +m_yfin[jph2] +m_zfin[jph2] );
           hfac1 =    m_xSfac[jph1] *wmd(delp,y1,z1);
           hfac2 =    m_ySfac[jph2] *wmd(delq,y2,z2);
           if( m_KeyOrd == 0) {
              Disr1( gami,jph1, &gi1, &gi2, &ggi1, &ggi2, &gggi1, &gggi2);
           } else {
              Disr1a(gami,jph1, &gi1, &gi2, &ggi1, &ggi2, &gggi1, &gggi2);
           }
           Dfsr1(gamf,jph2, &gf1, &gf2, &ggf1, &ggf2);
        //---- O(alf2) -----, tree level
           dist20 =
                (gi1*gf1*andi11+ gi1*gf2*andi12
                +gi2*gf1*andi21+ gi2*gf2*andi22)*hfac1*hfac2;
           m_beta20 = dist20
                -m_Beta00*m_xSfac[jph1]*m_ySfac[jph2]
                -m_beti10[jph1]*m_ySfac[jph2] -m_betf10[jph2]*m_xSfac[jph1];
           m_betxy20[jph1][jph2]=m_beta20;
           m_xyBet20=m_xyBet20 +m_beta20 /m_xSfac[jph1]/m_ySfac[jph2];
        //---- O(alf3) -----, one loop level  !!!!!!!NEW
        // Note that virtual correction is factorized ISR*FSR, as usual
           dist21 =
                (ggi1*ggf1*andi11+ ggi1*ggf2*andi12
                +ggi2*ggf1*andi21+ ggi2*ggf2*andi22)*hfac1*hfac2;
           m_beta21 = dist21
                -m_Beta01*m_xSfac[jph1]*m_ySfac[jph2]
                -m_betx11[jph1]*m_ySfac[jph2] -m_bety11[jph2]*m_xSfac[jph1];
           m_xyBet21=m_xyBet21 +m_beta21 /m_xSfac[jph1]/m_ySfac[jph2];
        }//for jph2
     }// for jph1
  } else {
      for(int jph1=1; jph1 <= nphox; jph1++) {
         for(int jph2=1; jph2<= nphoy; jph2++) {
           m_betxy20[jph1][jph2]= 0;
         }//for
      }//for
  }// if vlim
//-----------------------------------------------------------
//-----------------------------------------------------------
//                beta3 initial-initial-initial
//-----------------------------------------------------------
      m_xxxBet30=0;
      m_sbti30=0;
      if( DB->KeyISR != 0  &&  vv > m_vlim2) {
         for(int jph3 = 3; jph3<= nphox; jph3++){
            for(int jph2 = 2; jph2<= jph3-1; jph2++){
               for(int jph1 = 1; jph1<= jph2-1; jph1++){
                  hfac1  =  m_xSfac[jph1]*wmd(delp,m_yini[jph1],m_zini[jph1]);
                  hfac2  =  m_xSfac[jph2]*wmd(delp,m_yini[jph2],m_zini[jph2]);
                  hfac3  =  m_xSfac[jph3]*wmd(delp,m_yini[jph3],m_zini[jph3]);
//      Summation over 6 LL fragmentation trees for 3 ISR photons,
//      photon jph1 is always harder because of low level M.C. generation
                  ntree = 6;      // for 2 ISR fragmentation trees
                  gi1 =  0;
                  gi2 =  0;
                  Disr3(jph1, jph1, jph2, jph3, &gi1,&gi2);
                  Disr3(jph1, jph2, jph1, jph3, &gi1,&gi2);
                  Disr3(jph1, jph1, jph3, jph2, &gi1,&gi2);
                  Disr3(jph1, jph2, jph3, jph1, &gi1,&gi2);
                  Disr3(jph1, jph3, jph1, jph2, &gi1,&gi2);
                  Disr3(jph1, jph3, jph2, jph1, &gi1,&gi2);
                  gf1 = 0.5;    // 0.5 for averaging over 2 choices
                  gf2 = 0.5;    // of Born angles in final state
//---- O(alf3) -----, tree level  !!!!!!!NEW
                  dist30 =
                       (gi1*gf1*andi11+ gi1*gf2*andi12
                       +gi2*gf1*andi21+ gi2*gf2*andi22) *hfac1*hfac2*hfac3/ntree;
                  m_beta30 = dist30
                       -m_Beta00 *m_xSfac[jph1] *m_xSfac[jph2] *m_xSfac[jph3]
                       -m_beti10[jph1] *m_xSfac[jph2] *m_xSfac[jph3]
                       -m_beti10[jph2] *m_xSfac[jph1] *m_xSfac[jph3]
                       -m_beti10[jph3] *m_xSfac[jph1] *m_xSfac[jph2]
                       -m_betxx20[jph1][jph2] *m_xSfac[jph3]
                       -m_betxx20[jph1][jph3] *m_xSfac[jph2]
                       -m_betxx20[jph2][jph3] *m_xSfac[jph1];
                  m_xxxBet30 = m_xxxBet30
                       +m_beta30/m_xSfac[jph1]/m_xSfac[jph2]/m_xSfac[jph3];
// Pure ISR, simply the same
                  m_sbti30 = m_sbti30
                       +m_beta30/m_xSfac[jph1]/m_xSfac[jph2]/m_xSfac[jph3];
               }// for
            }//for
         }//for
      }// if
//-----------------------------------------------------------
//                beta3 initial-initial-final
//-----------------------------------------------------------
      m_xxyBet30 = 0;
      if(  DB->KeyISR != 0  &&  vv > m_vlim2 && IsFSR != 0 &&  uu > m_vlim1) {
         for(int jph2=2; jph2<= nphox; jph2++){
            for(int jph1=1; jph1<= jph2-1; jph1++){
               for(int jph3=1; jph3 <= nphoy; jph3++){
                  hfac1  =  m_xSfac[jph1]*wmd(delp,m_yini[jph1],m_zini[jph1]);
                  hfac2  =  m_xSfac[jph2]*wmd(delp,m_yini[jph2],m_zini[jph2]);
                  y3  = m_yfin[jph3] /(1 +m_yfin[jph3] +m_zfin[jph3] );
                  z3  = m_zfin[jph3] /(1 +m_yfin[jph3] +m_zfin[jph3] );
                  hfac3  =  m_ySfac[jph3]*wmd(delq,y3,z3);
                  ntree = 2;      // for 2 ISR fragmentation trees
                  gi1  =  0;     // initialization
                  gi2  =  0;     // initialization
                  Disr2(gami,jph1, jph1,jph2,&gi1,&gi2,&ggi1,&ggi2);  // 1-st tree
                  Disr2(gami,jph1, jph2,jph1,&gi1,&gi2,&ggi1,&ggi2);  // 2-nd tree
                  Dfsr1(gamf,jph3,           &gf1,&gf2,&ggf1,&ggf2);
//---- O(alf3) -----, tree level  !!!!!!!NEW
                  dist30= (gi1*gf1*andi11 + gi1*gf2*andi12
                          +gi2*gf1*andi21 + gi2*gf2*andi22)*hfac1*hfac2*hfac3/ntree;
                  m_beta30 = dist30
                       -m_Beta00 *m_xSfac[jph1] *m_xSfac[jph2] *m_ySfac[jph3]
                       -m_beti10[jph1] *m_xSfac[jph2] *m_ySfac[jph3]
                       -m_beti10[jph2] *m_xSfac[jph1] *m_ySfac[jph3]
                       -m_betf10[jph3] *m_xSfac[jph1] *m_xSfac[jph2]
                       -m_betxx20[jph1][jph2] *m_ySfac[jph3]
                       -m_betxy20[jph1][jph3] *m_xSfac[jph2]
                       -m_betxy20[jph2][jph3] *m_xSfac[jph1];
                  m_xxyBet30 = m_xxyBet30
                       +m_beta30/m_xSfac[jph1]/m_xSfac[jph2]/m_ySfac[jph3];
               }//for
            }//for
         }//for
      }// if vlim
//-----------------------------------------------------------
//                beta3 initial-final-final
//-----------------------------------------------------------
      m_xyyBet30 =  0;
      if(  DB->KeyISR != 0 &&  vv > m_vlim2 && IsFSR != 0 &&  uu > m_vlim1) {
         for(int  jph2=2; jph2<= nphoy; jph2++){
            for(int jph1=1; jph1<= jph2-1; jph1++){
               for(int jph3=1;  jph3<= nphox; jph3++){
                  y1  = m_yfin[jph1] /(1 +m_yfin[jph1] +m_zfin[jph1] );
                  z1  = m_zfin[jph1] /(1 +m_yfin[jph1] +m_zfin[jph1] );
                  y2  = m_yfin[jph2] /(1 +m_yfin[jph2] +m_zfin[jph2] );
                  z2  = m_zfin[jph2] /(1 +m_yfin[jph2] +m_zfin[jph2] );
                  hfac1  =  m_ySfac[jph1]*wmd(delq,y1,z1);
                  hfac2  =  m_ySfac[jph2]*wmd(delq,y2,z2);
                  hfac3  =  m_xSfac[jph3]*wmd(delp,m_yini[jph3],m_zini[jph3]);
//     sum over two FSR fragmentation trees
                  ntree = 2;     // for 2 FSR fragmentation trees
                  gf1 =  0;     // initialization
                  gf2 =  0;     // initialization
                  Dfsr2(jph1, jph1,jph2, &gf1,&gf2);
                  Dfsr2(jph1, jph2,jph1, &gf1,&gf2);
                  if( m_KeyOrd == 0) {
                     Disr1( gami,jph3,&gi1,&gi2,&ggi1,&ggi2,&gggi1,&gggi2);
                  } else {
                     Disr1a(gami,jph3,&gi1,&gi2,&ggi1,&ggi2,&gggi1,&gggi2);
                  }//if
//---- O(alf3) -----, tree level  !!!!!!!NEW
                  dist30= (gi1*gf1*andi11 + gi1*gf2*andi12
                          +gi2*gf1*andi21 + gi2*gf2*andi22) *hfac1*hfac2*hfac3/ntree;
                  m_beta30 = dist30
                       -m_Beta00 *m_ySfac[jph1] *m_ySfac[jph2] *m_xSfac[jph3]
                       -m_betf10[jph1] *m_ySfac[jph2] *m_xSfac[jph3]
                       -m_betf10[jph2] *m_ySfac[jph1] *m_xSfac[jph3]
                       -m_beti10[jph3] *m_ySfac[jph1] *m_ySfac[jph2]
                       -m_betyy20[jph1][jph2] *m_xSfac[jph3]
                       -m_betxy20[jph3][jph1] *m_ySfac[jph2]
                       -m_betxy20[jph3][jph2] *m_ySfac[jph1];
                  m_xyyBet30 = m_xyyBet30
                       +m_beta30/m_ySfac[jph1]/m_ySfac[jph2]/m_xSfac[jph3];
               }//for
            }//for
         }//for
      }//if vlim
//-----------------------------------------------------------------
// Finite part of the YFS formfactor for the ISR/FSR
//-----------------------------------------------------------------
      double ForIni = YFSkon_ini;
      double ForFin = YFSkon_fin;
      if(DB->KeyISR == 0)   ForIni=1;
      if(     IsFSR == 0)   ForFin=1;
      double fYFS = ForIni*ForFin;
//-----------------------------------------------------------------
//     and the rejection weights = (new.distr/crude.distr)
//-----------------------------------------------------------------
//     ============================================
//     ========== INITIAL + FINAL =================
//     ============================================

// Note that m_xyBet20 (which is genuine O(alf2)) is added to O(alf1),
// because our semianalytical programs are only able to deal
// with factorized ini/fin, see also O(alf1) definitions of beta1's.
// Total's, all beta's ---------------------------------
   m_WtSet[71] =   fYFS*  m_Beta00/DisCru;
   m_WtSet[72] =   fYFS*( m_Beta01 +m_xBet10 +m_yBet10 +m_xyBet20)/DisCru;
   m_WtSet[73] =   fYFS*( m_Beta02 +m_xBet11 +m_yBet11
                         +m_xxBet20 +m_xyBet20 +m_yyBet20 )/DisCru;
// !!!NEW!!!
   m_WtSet[74] =   fYFS*( m_Beta03 +m_xBet12 +m_yBet12
                       +m_xxBet21 +m_xyBet21 +m_yyBet21
                       +m_xxxBet30 +m_xxyBet30 +m_xyyBet30 )/DisCru;
// First order, individual beta's -------------
   m_WtSet[80] =   fYFS*m_Beta01/DisCru;
   m_WtSet[81] =   fYFS*(m_xBet10+m_yBet10)/DisCru;
   m_WtSet[82] =   fYFS*(m_xBet10)/DisCru;
   m_WtSet[83] =   fYFS*(m_yBet10)/DisCru;
   m_WtSet[84] =   fYFS*(m_xyBet20)/DisCru;
// Second order, individual beta's ------------
   m_WtSet[90] =   fYFS*m_Beta02/DisCru;
   m_WtSet[91] =   fYFS*(m_xBet11+m_yBet11)/DisCru;
   m_WtSet[92] =   fYFS*(m_xxBet20+m_xyBet20+m_yyBet20)/DisCru;
   m_WtSet[93] =   fYFS*(m_xBet11)/DisCru;
   m_WtSet[94] =   fYFS*(m_yBet11)/DisCru;
   m_WtSet[95] =   fYFS*(m_xxBet20)/DisCru;
   m_WtSet[96] =   fYFS*(m_xyBet20)/DisCru;
   m_WtSet[97] =   fYFS*(m_yyBet20)/DisCru;
// Third order, individual beta's ------------!!!NEW!!!
   m_WtSet[100] =   fYFS*m_Beta03/DisCru;
   m_WtSet[101] =   fYFS*(m_xBet12 +m_yBet12)/DisCru;
   m_WtSet[102] =   fYFS*(m_xxBet21+m_xyBet21+m_yyBet21)/DisCru;
   m_WtSet[103] =   fYFS*m_xBet12/DisCru;
   m_WtSet[104] =   fYFS*m_yBet12/DisCru;
   m_WtSet[105] =   fYFS*m_xxBet21/DisCru;
   m_WtSet[106] =   fYFS*m_xyBet21/DisCru;
   m_WtSet[107] =   fYFS*m_yyBet21/DisCru;
   m_WtSet[108] =   fYFS*(m_xxxBet30+m_xxyBet30+m_xyyBet30)/DisCru;
   m_WtSet[109] =   fYFS*m_xxxBet30/DisCru;
   m_WtSet[110] =   fYFS*m_xxyBet30/DisCru;
   m_WtSet[111] =   fYFS*m_xyyBet30/DisCru;
//     ============================================
//     ========= INITIAL STATE ALONE ==============
//     ============================================
// Total's, all beta's ---------------------------------
   m_WtSet[ 1] =   ForIni* m_beti00/DisCru;
   m_WtSet[ 2] =   ForIni*(m_beti01+m_sbti10)/DisCru;
   m_WtSet[ 3] =   ForIni*(m_beti02+m_sbti11+m_sbti20)/DisCru;
//!!!NEW
   m_WtSet[ 4] =   ForIni*(m_beti03+m_sbti12+m_sbti21+m_sbti30)/DisCru;
// First order, individual beta's -------------
   m_WtSet[10] =   ForIni*m_beti01/DisCru;
   m_WtSet[11] =   ForIni*m_sbti10/DisCru;
// Second order, individual beta's ------------
   m_WtSet[20] =   ForIni*m_beti02/DisCru;
   m_WtSet[21] =   ForIni*m_sbti11/DisCru;
   m_WtSet[22] =   ForIni*m_sbti20/DisCru;
//!!!NEW
// Third order, individual beta's ------------
   m_WtSet[30] =   ForIni*m_beti03/DisCru;
   m_WtSet[31] =   ForIni*m_sbti12/DisCru;
   m_WtSet[32] =   ForIni*m_sbti21/DisCru;
   m_WtSet[33] =   ForIni*m_sbti30/DisCru;
//=================================================================//
//            ISR   Non-exponentiated version                      //
//            Not yet fully implemented......                      //
//=================================================================//
// Entire 0,1,2-photon distributions
   double m_dis0,m_dis1,m_dis2,m_dig1;
   double fYFSu = exp( -gami*log(1/DB->vvmin) );
   m_dis0   =0;
   m_dis1   =0;
   m_dis2   =0;
   m_dig1   =0;
   double Bor1;
   if( (nphox+nphoy) == 0) {
       m_dis0 = 1;
       m_dis1 = 1+gami*log(DB->vvmin);
       m_dis2 = 1+gami*log(DB->vvmin)+0.5* sqr(gami*log(DB->vvmin));
   } else if( nphox == 1) {
       y = m_yini[1];
       z = m_zini[1];
       gf1 = 0.5;
       gf2 = 0.5;
       gi1 = (sqr(1-y))/2;
       gi2 = (sqr(1-z))/2;
       Bor1=  gi1*gf1*andi11   +gi1*gf2*andi12
             +gi2*gf1*andi21   +gi2*gf2*andi22;
       m_dis1 = Bor1  *wmd(delp,y,z);            // S-factor divided off
// standard O(alf1) for comparisons with spinor methods
       m_dig1 = Bor1  *2/(y*z)*wm1(delp,y,z)*svar/svar1;
       m_dis2 = m_dis1*(1 +gami*log(DB->vvmin));
   }else  if( nphoy == 1) {
       yy = m_yfin[1];
       zz = m_zfin[1];
       y  = yy/(1 +yy+zz);
       z  = zz/(1 +yy+zz);
       gi1 = 0.5;
       gi2 = 0.5;
       gf2 = (sqr(1-y) ) /2;   // y,z are swapped! correct d_fsr !!!!
       gf1 = (sqr(1-z) ) /2;
       Bor1=    gi1*gf1*andi11   +gi1*gf2*andi12
               +gi2*gf1*andi21   +gi2*gf2*andi22;
// standard O(alf1) for comparisons with spinor methods
       m_dig1 = Bor1  *2/(y*z)*wm1(delq,y,z);
       m_dis1 = Bor1  *wmd(delq,y,z);            // S-factor divided off
   }else  if( nphox == 2) {
       m_dis2 = wm0(delp,m_yini[1],m_zini[1]) *wm0(delp,m_yini[2],m_zini[2]);
   }//if
// UNEXP Total O(alf0),O(alf1),O(alf2)
    m_WtSet[160] =    m_dis0 /fYFSu;
    m_WtSet[161] =    m_dis1 /fYFSu;
    m_WtSet[162] =    m_dis2 /fYFSu;
//|=================================================================|
//|        Model weight (the best)                                  |
//|=================================================================|
    //m_WtBest = m_WtSet(m_IdeWgt)
}// Make

void KKqed3::Disr1a(double gami, int j1, double *g1, double *g2, double *gg1, double *gg2, double *ggg1, double *ggg2){
/////////////////////////////////////////////////////////////////////////////////////
//     !!!!!! NEW VERSION with PT ordering !!!!!                                   //
//                                                                                 //
// Ingredients for O(alf3)NLL ISR matrix element.                                  //
// INPUT:                                                                          //
//     alfinv=  QED coupling                                                       //
//     charg2=  charge squared                                                     //
//     y,z=     Sudakov variables                                                  //
//     j1=      pointers to input-photon                                           //
// OUTPUT:                                                                         //
//     g's are set here: g=treelevel, gg=oneloop, ggg=twoloop                      //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
double a1 = m_yini[j1];
double b1 = m_zini[j1];
double zz = (1-a1)*(1-b1);
if(zz <= 0) cout<<"+++++ KKqed3::Disr1a zz="<<zz<<endl;
double dels1 = gami/2                   // LL constant part
          +1/DB->Alfinv0/M_PI*(
          +log(a1)*log(1-b1)  +log(b1)*log(1-a1)    //   LL part
          +Dilog(a1)          +Dilog(b1)          // NLL this and all the rest
          -0.5*sqr(log(1-a1))  -0.5*sqr(log(1-b1))
          +1.5*log(1-a1)       +1.5*log(1-b1)
          +0.5*a1*(1-a1)/(1+sqr(1-a1))
          +0.5*b1*(1-b1)/(1+sqr(1-b1))
       );
//****      dels1 =  gami/2 -gami/4*dlog(zz) !! averaged LL version
double dels2 =  sqr(gami)/8.0
        -sqr(gami)/8  *log(zz)
        +sqr(gami)/24 *sqr(log(zz));
// Exact O(alf1) matrix element for the hardest photon jhard
*g1   = (sqr(1-a1)            )/2;
*g2   = (            sqr(1-b1))/2;
*gg1  = (sqr(1-a1)            )/2 *(1+dels1);
*gg2  = (            sqr(1-b1))/2 *(1+dels1);
*ggg1 = (sqr(1-a1)            )/2 *(1+dels1+dels2);
*ggg2 = (            sqr(1-b1))/2 *(1+dels1+dels2);
}// KKqed3::Disr1a


void KKqed3::Disr1(double gami, int j1, double *g1, double *g2, double *gg1, double *gg2, double *ggg1, double *ggg2){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
// Ingredients for O(alf3)NLL ISR matrix element.                                  //
// INPUT:                                                                          //
//     alfinv=  QED coupling                                                       //
//     charg2=  charge squared                                                     //
//     y,z=     Sudakov variables                                                  //
//     j1=      pointers to input-photon                                           //
// OUTPUT:                                                                         //
//     g's are set here: g=treelevel, gg=oneloop, ggg=twoloop                      //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
double a1 = m_yini[j1];
double b1 = m_zini[j1];
double zz = (1-a1)*(1-b1);
if(zz <= 0) cout<<"+++++ KKqed3::Disr1 zz="<<zz<<endl;
double dels1 =  gami/2 -gami/4*log(zz);
double dels2 =  sqr(gami)/8
             -sqr(gami)/8  *log(zz)
             +sqr(gami)/24 *sqr(log(zz));
// Exact O(alf1) matrix element for the hardest photon jhard
*g1   = (sqr(1-a1)            )/2;
*g2   = (            sqr(1-b1))/2;
*gg1  = (sqr(1-a1)            )/2 *(1+dels1);
*gg2  = (            sqr(1-b1))/2 *(1+dels1);
*ggg1 = (sqr(1-a1)            )/2 *(1+dels1+dels2);
*ggg2 = (            sqr(1-b1))/2 *(1+dels1+dels2);
}//KKqed3::Disr1


void KKqed3::Disr2a(double gami, int j1,int j2, double *g1, double *g2, double *gg1, double *gg2){
/////////////////////////////////////////////////////////////////////////////////////
//     !!!!!! NEW VERSION with PT ordering !!!!!                                   //
// Ingredients for O(alf2)NLL ISR matrix element.                                  //
// INPUT:                                                                          //
//     gami  = 2*alfa/pi*(BigLog-1)                                                //
//     y,z   = Sudakov variables                                                   //
//     j1,j2 = pointers of two input-photons , j1 should have gigger pT            //
// OUTPUT:                                                                         //
//     g's, gg's are updated (have to be initialized in calling program)           //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
double  a1,z1,p2,delvir1,z1z2,a2,b1,p1,b2;
//-------------------------------------------------------------------------------------
a2 = m_yini[j2];
b2 = m_zini[j2];
a1 = m_yini[j1]/(1-m_yini[j2]);
b1 = m_zini[j1]/(1-m_zini[j2]);
if(a1  > 1) cout<<"+++++ KKqed3::Disr2a: a1="<<a1<<endl;
// Exact O(alf1) matrix element for the hardest photon jhard
p1= (sqr(1-a1)            ) *( sqr(1-a2) + sqr(1-b2) )/4;
p2= (            sqr(1-b1)) *( sqr(1-a2) + sqr(1-b2) )/4;
*g1 = *g1 +2*p1;
*g2 = *g2 +2*p2;
z1 =  (1-m_yini[j1])*(1-m_zini[j1]);
z1z2= (1-m_yini[j1]-m_yini[j2])*(1-m_zini[j1]-m_zini[j2]);
// soft limit to QED3_Disr1 OK, for 2 trees we get 3 terms gami/6*dlog(zz)
delvir1 = gami/2 -gami/6*log(z1) -gami/6*log(z1z2);
*gg1=*gg1 +2*p1*(1+delvir1);
*gg2=*gg2 +2*p2*(1+delvir1);
//
if(z1    <= 0) cout<<"+++++ KKqed3::Disr2a: z1="<<z1<<endl;
if(z1z2  <= 0) cout<<"+++++ KKqed3::Disr2a: z1z2="<<z1z2<<endl;
}//KKqed3::Disr2a



void KKqed3::Disr2(double gami, int jhard, int j1, int j2, double *g1, double *g2, double *gg1, double *gg2){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
// Ingredients for O(alf2)NLL ISR matrix element.                                  //
// INPUT:                                                                          //
//     gami  = 2*alfa/pi*(BigLog-1)                                                //
//     y,z   = Sudakov variables                                                   //
//     jhard = pointer of hardes photon                                            //
//     j1,j2 = pointers of two input-photons                                       //
// OUTPUT:                                                                         //
//     g's, gg's are updated (have to be initialized in calling program)           //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
double  a1,z1,p2,delvir1,z1z2,a2,b1,p1,b2;
//-------------------------------------------------------------------------------------
a1 = m_yini[j1];
b1 = m_zini[j1];
a2 = m_yini[j2]/(1-m_yini[j1]);
b2 = m_zini[j2]/(1-m_zini[j1]);
// Exact O(alf1) matrix element for the hardest photon jhard
if(jhard == j1) {
   p1= (sqr(1-a1)            ) *( sqr(1-a2) + sqr(1-b2) )/4;
   p2= (            sqr(1-b1)) *( sqr(1-a2) + sqr(1-b2) )/4;
} else {
   p1= (sqr(1-a1) +sqr(1-b1) ) *( sqr(1-a2)             )/4;
   p2= (sqr(1-a1) +sqr(1-b1) ) *(             sqr(1-b2) )/4;
}
*g1 = *g1 +p1;
*g2 = *g2 +p2;
z1 =  (1-m_yini[j1])*(1-m_zini[j1]);
z1z2= (1-m_yini[j1]-m_yini[j2])*(1-m_zini[j1]-m_zini[j2]);
// soft limit to QED3_Disr1 OK, for 2 trees we get 3 terms gami/6*dlog(zz)
delvir1 = gami/2 -gami/6*log(z1) -gami/6*log(z1z2);
*gg1=*gg1 +p1*(1+delvir1);
*gg2=*gg2 +p2*(1+delvir1);

if(z1   <= 0) cout<<"+++++ KKqed3::Disr2: z1="<<z1 << endl;
if(z1z2 <= 0) cout<<"+++++ KKqed3::Disr2: z1z2="<< z1z2<< endl;
} //KKqed3::Disr2


void KKqed3::Disr3(int jhard, int j1, int j2, int j3, double *g1, double *g2){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
// Ingredients for O(alf3)LL ISR matrix element.                                   //
// INPUT:                                                                          //
//     y,z Sudakov variables                                                       //
//     jhard pointer of hardes photon                                              //
//     j1,j2,j3 pointers of 3 input-photons                                        //
// OUTPUT:                                                                         //
//     g1,g2 are updated (have to be initialized in calling program)               //
/////////////////////////////////////////////////////////////////////////////////////
double  b3,a3,p2,p1,b1,a1,b2,a2,a,b;
//-------------------------------------------------------------------------------------
a1 = m_yini[j1];
b1 = m_zini[j1];
a2 = m_yini[j2]/(1-m_yini[j1]);
b2 = m_zini[j2]/(1-m_zini[j1]);
a3 = m_yini[j3]/(1-m_yini[j2]-m_yini[j1]);
b3 = m_zini[j3]/(1-m_zini[j2]-m_zini[j1]);
if( m_KeyOrd != 0 ) {
   a3 = m_yini[j3];
   b3 = m_zini[j3];
   a2 = m_yini[j2]/(1-m_yini[j3]);
   b2 = m_zini[j2]/(1-m_zini[j3]);
   a1 = m_yini[j1]/(1-m_yini[j2]-m_yini[j3]);
   b1 = m_zini[j1]/(1-m_zini[j2]-m_zini[j3]);
}
// Exact O(alf1) matrix element for the hardest photon jhard
if(jhard == j1) {
   p1= chi1(a1) *chi2(a2,b2) *chi2(a3,b3);
   p2= chi1(b1) *chi2(a2,b2) *chi2(a3,b3);
} else if(jhard == j2) {
   p1= chi2(a1,b1) *chi1(a2) *chi2(a3,b3);
   p2= chi2(a1,b1) *chi1(b2) *chi2(a3,b3);
} else {
   p1= chi2(a1,b1) *chi2(a2,b2) *chi1(a3);
   p2= chi2(a1,b1) *chi2(a2,b2) *chi1(b3);
}
*g1= *g1 +p1;
*g2= *g2 +p2;
}// KKqed3::Disr3


void KKqed3::Dfsr1(double gamf, int j1, double *g1, double *g2, double *gg1, double *gg2){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
// Ingredients for O(alf2)NLL FSR matrix element.                                  //
// INPUT:                                                                          //
//     y,z Sudakov variables                                                       //
//     j1 pointer of input-photons                                                 //
// OUTPUT:                                                                         //
//     g's are set here                                                            //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
double zz,dels1,a1,b1;
//-------------------------------------------------------------------------------------
// normal definition as in O(alf1) single-photon case
      a1 = m_yfin[j1]/( 1 +m_yfin[j1] +m_zfin[j1] );
      b1 = m_zfin[j1]/( 1 +m_yfin[j1] +m_zfin[j1] );
      zz = (1-a1)*(1-b1);
      if(zz  <= 0) cout<<"+++++ KKqed3::Dfsr1: zz="<<zz<<endl;
      dels1 =  gamf/2 +gamf/4*log(zz);
// Exact O(alf1) matrix element for the hardest photon jhard
      *g2 = (sqr(1-a1)            ) /2;              // corrected
      *g1 = (            sqr(1-b1)) /2;              // corrected
      *gg2= (sqr(1-a1)            ) /2 *(1+dels1);   // corrected
      *gg1= (            sqr(1-b1)) /2 *(1+dels1);   // corrected
}// KKqed3::Dfsr1

void KKqed3::Dfsr2(int jhard, int j1, int j2, double *g1, double *g2){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
// Ingredients for O(alf2)NLL FSR matrix element.                                  //
// INPUT:                                                                          //
//     y,z Sudakov variables                                                       //
//     jhard pointer of hardes photon                                              //
//     j1,j2 pointers of two input-photons                                         //
// OUTPUT:                                                                         //
//     g1,g2 are updated (have to be initialized in calling program)               //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
double  b2,a2,p1,p2,b1,a1,zp2,yp2;
//-------------------------------------------------------------------------------------
// normal definition as in O(alf1) single-photon case
      a1 = m_yfin[j1]/( 1 +m_yfin[j1] +m_zfin[j1] );
      b1 = m_zfin[j1]/( 1 +m_yfin[j1] +m_zfin[j1] );
// take into account primary photon emission
      yp2 = m_yfin[j2]/( 1 +m_yfin[j1] );
      zp2 = m_zfin[j2]/( 1 +m_zfin[j1] );
// as in O(alf1) single-photon case
      a2 = yp2/(1 + yp2 +zp2);
      b2 = zp2/(1 + yp2 +zp2);
// Exact O(alf1) matrix element for the hardest photon jhard
      if(jhard == j1) {
         p2= (sqr(1-a1)            ) *( sqr(1-a2) + sqr(1-b2) )/4; // corrected
         p1= (            sqr(1-b1)) *( sqr(1-a2) + sqr(1-b2) )/4; // corrected
      } else {
         p2= (sqr(1-a1) +sqr(1-b1) ) *( sqr(1-a2)             )/4; // corrected
         p1= (sqr(1-a1) +sqr(1-b1) ) *(             sqr(1-b2) )/4; // corrected
      }
      *g1= *g1 +p1;
      *g2= *g2 +p2;
}// KKqed3::Dfsr2


void KKqed3::bvirt0(double alfinv, double charg2, double svar, double am, double *dels1,  double *dels2,  double *dels3){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
// ISR/FSR virtual corrections to beta0                                            //
// beta0 is equal born*(1+dels1+dels2+dels3)                                       //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
double pi=3.1415926535897932e0;
double zet2= pi*pi/6e0;
double zet3= 1.2020569031595942854e0;
double gami,bilg,alfpi,Mlog;
//-------------------------------------------------------------------------------------
//DOUBLE COMPLEX     F1_1, F1_2, cL
//-------------------------------------------------------------------------------------
alfpi =  1/DB->Alfinv0/pi;
bilg  =  log(svar/am/am);
Mlog  =  log(svar/am/am) -1;
gami  =  2*charg2*alfpi*(bilg-1);
*dels1 =  gami/2;
// ISR with subleading terms from Berends, Burgers, Van Neerveen
// (the effect of including NLL is negligible, below 1d-4)
*dels2 =
      sqr(charg2)*sqr(alfpi)  *0.5*sqr(bilg);                    // LL
//           +charg2**2*alfpi**2*(
//             -(13d0/16d0 +1.5d0*zet2 -3d0*zet3)*bilg           ! NLL
//             -16d0/5d0*zet2*zet2 +51d0/8d0*zet2 +13d0/4d0      ! NNLL
//             -4.5d0*zet3 -6d0*zet2*log(2d0)                    ! NNLL
//            )
*dels3 = 1.0/6.0 *(gami/2)*sqr(gami/2);
/////////////////////////////////////////////////////////////////////////
// The assignements below will get together O(alf1)CEEX and O(alf1)EEX
// but it will spoil O(alf2)EEX because dels1 is also input for beta1
//      cL    = DCMPLX( DLOG(Svar/Am**2)-1d0, -1d0 )
//      F1_1  = Alfpi*charg2   *0.5d0*cL
//      dels1 = CDABS(1+ F1_1)**2 -1d0
//      F1_2 = F1_1
//           +(Alfpi*charg2)**2 *(
//                    +cL**2/8d0
//                    +cL*( 3d0/32 -3d0/4*zet2 +3d0/2*zet3 )
//           )
//      dels2 = CDABS(1+ F1_2)**2 -(1d0+dels1)
/////////////////////////////////////////////////////////////////////////
}//KKqed3::bvirt0


double KKqed3::Dilog(double x){
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
// dilogarithm FUNCTION: dilog(x)=int( -ln(1-z)/z ) , 0 < z < x .                          //
// this is the cernlib version.                                                            //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////
double  a,b,y,s,t,z, dilog;
z=-1.644934066848226e0;
if(x  <  -1.0) goto e1;
if(x  <=  0.5) goto e2;
if(x  ==  1.0) goto e3;
if(x  <=  2.0) goto e4;
z=3.289868133696453e0;
e1:
  t=1.e0/x;
  s=-0.5e0;
  z=z-0.5e0*sqr(log(abs(x)));
goto e5;
e2:
  t=x;
  s=0.5e0;
  z=0.e0;
goto e5;
e3:
 dilog = 1.644934066848226e0;
 return dilog;
e4:
  t=1.e0-x;
  s=-0.5e0;
  z=1.644934066848226e0-log(x)*log(abs(t));
e5:
y=2.666666666666667e0*t+0.666666666666667e0;
b=      0.000000000000001e0;
a=y*b  +0.000000000000004e0;
b=y*a-b+0.000000000000011e0;
a=y*b-a+0.000000000000037e0;
b=y*a-b+0.000000000000121e0;
a=y*b-a+0.000000000000398e0;
b=y*a-b+0.000000000001312e0;
a=y*b-a+0.000000000004342e0;
b=y*a-b+0.000000000014437e0;
a=y*b-a+0.000000000048274e0;
b=y*a-b+0.000000000162421e0;
a=y*b-a+0.000000000550291e0;
b=y*a-b+0.000000001879117e0;
a=y*b-a+0.000000006474338e0;
b=y*a-b+0.000000022536705e0;
a=y*b-a+0.000000079387055e0;
b=y*a-b+0.000000283575385e0;
a=y*b-a+0.000001029904264e0;
b=y*a-b+0.000003816329463e0;
a=y*b-a+0.000014496300557e0;
b=y*a-b+0.000056817822718e0;
a=y*b-a+0.000232002196094e0;
b=y*a-b+0.001001627496164e0;
a=y*b-a+0.004686361959447e0;
b=y*a-b+0.024879322924228e0;
a=y*b-a+0.166073032927855e0;
a=y*a-b+1.935064300869969e0;
dilog=s*t*(a-b)+z;
return dilog;
}//Dilog
