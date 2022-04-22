///////////////////////////////////////////////////////////////////////////////
#include "KKdizet.h"

ClassImp(KKdizet);

/*
// for testing disk tables
extern "C" {
//  Imports EW tables form kkDizet common block
//  DOUBLE PRECISION FUNCTION hhDizet_GetELW(it,md,KFi,KFf,i,j,k)
  double hhdizet_getelw_(const int&,const int&,const int&,const int&,const int&,const int&,const int&);
//  Imports EW tables form kkDizet common block
//  DOUBLE PRECISION FUNCTION hhDizet_GetQCD(it,KFi,KFf,i,k)
  double hhdizet_getqcd_(const int&,const int&,const int&,const int&,const int&);
// SUBROUTINE hhDizet_GetSMpar(swsq, gammz, MW, GammW)
  void   hhdizet_getsmpar_(const double&, const double&, const double&, const double&);
}
*/

// for testing only
//extern "C" {
//BornV_PrintGSW(nout,KFi,KFf,svar,CosThe)
//  void bornv_printgsw_(const int&,const int&,const int&,const double&, const double&);
//}

const double  KKdizet::m_WminLEP1 =   0.010;  // LEP1 basic range
const double  KKdizet::m_WmaxLEP1 =  95.000;  // LEP1 basic range
const double  KKdizet::m_WdelZ    =  60.000;  // Z range (amz +- m_WdelZ)
const double  KKdizet::m_WmaxLEP2 = 240.001;  // LEP2 interval (m_WmaxLEP1,m_WmaxLEP2)
const double  KKdizet::m_WmaxNLC  = 8000.00;  // LHC/NLC range (m_WmaxLEP2,m_WmaxNLC)


KKdizet::KKdizet()
{
  // This constructor is for ROOT streamers ONLY
  // all pointers has to be NULLed
  cout<< "----> KKdizet Default Constructor (for ROOT only) "<<endl;
  m_Out= NULL;
}

///_____________________________________________________________
KKdizet::KKdizet(ofstream *OutFile)
{
  cout<< "----> KKdizet USER Constructor "<<endl;
  m_Out = OutFile;
}//KKdizet

///______________________________________________________________________________________
KKdizet::~KKdizet()
{
  //Explicit destructor
  cout<< "----> KKdizet::KKdizet !!!! DESTRUCTOR !!!! "<<endl;
}///destructor

double KKdizet::sqr( const Double_t x ){ return x*x;};

///______________________________________________________________________________________
void KKdizet::Initialize()
{
  cout  << "----> KKdizet::Initialize, Entering "<<endl;
//=================================================================
// BX*** macros are in MCdev/BXFORMAT.h
  BXOPE(*m_Out);
  BXTXT(*m_Out,"========================================");
  BXTXT(*m_Out,"======    KKdizet::Initialize      =====");
  BXTXT(*m_Out,"========================================");

  m_KeyQCD = 0; //  QCD FSR for final leptonic states is ABSENT!

  // Import EW look-up tables from common blocks of DZface (temporary)
  // ImportEWtabs();

}// Initialize

void KKdizet::ReadEWtabs(){
///////////////////////////////////////////////////////////////////
// Reading EW and QCD form-factors of DIZET from disk files
///////////////////////////////////////////////////////////////////
  char trail[200];
  char ch1;
  cout<<"==========================KKdizet::ReadEWtabs=========================="<<endl;
  double MZ, amh, amtop, swsq, gammz, MW, GammW;
  ifstream InputFile;

  int KFi,KFf;
  InputFile.open("DIZET-table1");
  int nChanel;
  InputFile >>  nChanel;
  InputFile.getline(trail,200); //  cout<<trail<<endl;
  cout<<"KKdizet::ReadEWtabs:  nChanel="<<nChanel<<endl;
  for(int ich=0;ich<nChanel;ich++){
  char label;
  for(int i=0;i<33;i++) InputFile.get(label); // reading front text
  InputFile>>  KFi>> KFf;
  InputFile.getline(trail,200); //  cout<<trail<<endl;
  cout<<"KKdizet::ReadEWtabs:  KFi, KFf =="<<KFi<<"  "<<KFf<<endl;
  KFi =1; // this is always electron beam KF-10!!!
//------------------------------------------------------------------
//  Values of m_MZ, m_amh, m_amtop are the same as input in xpar.
//  DIZET calculates SM values of m_swsq, m_gammz, m_MW, M_GammW
//  Let us give the special names D_swsq, D_GamZ, D_MW, D_GammW
//------------------------------------------------------------------
// 1st line, SM input/output parameters
  InputFile>> m_MZ >> m_amh >> m_amtop >> D_swsq >> D_GamZ>>D_MW>> D_GammW;
  cout<<" m_MZ="<< m_MZ <<"  m_amh="<< m_amh <<"  m_amtop="<< m_amtop
      <<"  D_swsq="<< D_swsq <<"  D_GamZ="<< D_GamZ
      <<"  D_MW="<<D_MW<<"  D_GammW="<< D_GammW<< endl;
  InputFile.getline(trail,200);
/////////////////////////////////////////////////////////////////
//   basic s range LEP1 and below        m_cyy(m_poin1+1, 7)   //
/////////////////////////////////////////////////////////////////
//  DO i=0, m_poin1
//     READ(27,m_fmt1) chr,n,ww
//     READ(27,m_fmt2) (m_cyys(i+1,k,KFi,KFf),k=1,m_poinG) ! EW
//////////////////////////////////////////////////////////////////
  //
  for(int i=0;i<= m_poin1;i++){
	 InputFile.getline(trail,200); //cout<<trail<<endl;
     for(int k=1; k<= m_poinG;k++) InputFile >> m_cyys[i][k-1][KFi-1][KFf-1];
     InputFile.getline(trail,200);
     //
  }// for i
/////////////////////////////////////////////////////////////////
//             near Z0 resonance    m_czz(m_poin2+1, 7)        //
/////////////////////////////////////////////////////////////////
//  DO i=0,m_poin2
//     DO  j=0,m_poTh2
//       WRITE(m_ndisk,m_fmt1) 'b',i,ww,j,cosi
//       WRITE(m_ndisk,*)   (m_czzs(i+1,j+1,k,KFi,KFf),k=1,m_poinG)
/////////////////////////////////////////////////////////////////
  for(int i=0;i<= m_poin2;i++){
	  for(int j=0;j<= m_poTh2;j++){
	     InputFile.getline(trail,200); //cout<<trail<<endl;
	     for(int k=1; k<= m_poinG;k++) InputFile >> m_czzs[i][j][k-1][KFi-1][KFf-1];
	     InputFile.getline(trail,200);
	  }//j
  }//i
  /////////////////////////////////////////////////////////////////////
  //   the region of boxes, LEP2,   m_ctt(m_poin3+1, m_poTh3+1, 7)   //
  /////////////////////////////////////////////////////////////////////
  for(int i=0;i<= m_poin3;i++){
	  for(int j=0;j<= m_poTh3;j++){
	     InputFile.getline(trail,200); //cout<<trail<<endl;
	     for(int k=1; k<= m_poinG;k++) InputFile >> m_ctts[i][j][k-1][KFi-1][KFf-1];
	     InputFile.getline(trail,200);
	  }//j
  }//i
  /////////////////////////////////////////////////////////////////////
  //   the region of boxes, NLC,    m_clc(m_poin4+1, m_poTh4+1, 7)   //
  /////////////////////////////////////////////////////////////////////
  for(int i=0;i<= m_poin4;i++){
	  for(int j=0;j<= m_poTh4;j++){
	     InputFile.getline(trail,200); //cout<<trail<<endl;
	     for(int k=1; k<= m_poinG;k++) InputFile >> m_clcs[i][j][k-1][KFi-1][KFf-1];
	     InputFile.getline(trail,200);
	  }//j
  }//i
  //-----------
  }//ich
  InputFile.close();
  cout<<"==========================  ReadEWtabs END  ==========================="<<endl;
}//ReadEWtabs


void KKdizet::InterpoGSW(int KFi0, int KFf0, double svar, double CosThe){
//////////////////////////////////////////////////////////////////////////
//  Calculates GSW formfactors from tables using linear interpolation
//  of the matrix in KFi and KFf constructed in hhDizet.
//  For compatibility with KKMC, KFi is passed in the BornV common block.
//
//////////////////////////////////////////////////////////////////////////
//	Projecting into flavours recorded in the EW tables
    int KFi=1; // always electron
    int KFf=KFf0;
//  if(KFi0==1  || KFi0==3  || KFi0==5 )  KFi=1;  // d,s,b
//  if(KFi0==2  || KFi0==4 )              KFi=2;  // u,c
//    if(KFf0==11 || KFf0==13 || KFf0==15 ) KFf=13; // mu,tau
//    if(KFf0==12 || KFf0==14 || KFf0==16 ) KFf=14; // nue,numu,nutau
//
    double ww = sqrt(svar);
    double MZ = 91.1876;
    m_WminZ = MZ-m_WdelZ;
    m_WmaxZ = MZ+m_WdelZ;
    // Defaults, without FSR QCD. No need to interpolate if KeyQCD = 0.
    m_KeyQCD = 0; //  QCD FSR for final leptonic states is ABSENT!
    m_QCDcorR[0] = 1;
    m_QCDcorR[1] = 0;
    m_QCDcorR[2] = 0;
    m_QCDcorR[3] = 0;

    double x,h,y,hy;
    int i,j;
    if(  (ww >= m_WminZ) && (ww <= m_WmaxZ) ) {
    //  LEP1 near Z0 resonance
       x= (ww-m_WminZ)/(m_WmaxZ-m_WminZ);
       i= min( int(trunc(m_poin2*x))+1, m_poin2);
       h= x*m_poin2-double(i-1);
       y= (1+CosThe)/2;
       j= min( int(trunc(m_poTh2*y))+1, m_poTh2);
       hy= y*m_poTh2-double(j-1);
       for(int kk=1; kk<=m_poinG; kk++){
          m_GSW[kk-1]=m_czzs[i-1][j-1][kk-1][KFi-1][KFf-1]*(1-h)*(1-hy)
                     +m_czzs[i  ][j-1][kk-1][KFi-1][KFf-1]*h*(1-hy)
                     +m_czzs[i-1][j  ][kk-1][KFi-1][KFf-1]*(1-h)*hy
                     +m_czzs[i  ][j  ][kk-1][KFi-1][KFf-1]*h*hy;
       }
       if(m_KeyQCD != 0)
         for(int kk=1; kk<=m_poinQ; kk++){
           m_QCDcorR[kk-1]=m_szzs[i-1][kk-1][KFi-1][KFf-1]*(1-h)
                        +m_szzs[i  ][kk-1][KFi-1][KFf-1]*h;
       }
    } else if(  (ww >= m_WminLEP1) && (ww <= m_WmaxLEP1) ) {
    	// LEP1 outside Z0 and low energies
    	 x= log( ww/m_WminLEP1) / log(m_WmaxLEP1/m_WminLEP1);
    	 i= min( int(trunc(m_poin1*x))+1, m_poin1);
    	 h= x*m_poin1-double(i-1);
    	 for(int kk=1; kk<=m_poinG; kk++){
    	            m_GSW[kk-1]  =m_cyys[i-1][kk-1][KFi-1][KFf-1]*(1-h)
    	                         +m_cyys[i  ][kk-1][KFi-1][KFf-1]*h;
    	 }//for kk
    	 if(m_KeyQCD != 0) {
    	   for(int kk=1; kk<= m_poinQ; kk++)
    	          m_QCDcorR[kk-1]=m_syys[i-1][kk-1][KFi-1][KFf-1]*(1-h)
    	                         +m_syys[i  ][kk-1][KFi-1][KFf-1]*h;
    	 }
    } else if( (ww >= m_WmaxLEP1) && (ww <= m_WmaxLEP2) ) {
    // in the LEP2 region
        x= (ww-m_WmaxLEP1)/(m_WmaxLEP2-m_WmaxLEP1);
        i= min( int(trunc(m_poin3*x))+1, m_poin3);
        h= x*m_poin3-double(i-1);
        y= (1+CosThe)/2;
        j= min( int(trunc(m_poTh3*y))+1, m_poTh3);
        hy= y*m_poTh3-double(j-1);
        for(int kk=1; kk<=m_poinG; kk++){
           m_GSW[kk-1]=m_ctts[i-1][j-1][kk-1][KFi-1][KFf-1]*(1-h)*(1-hy)
                      +m_ctts[i  ][j-1][kk-1][KFi-1][KFf-1]*h*(1-hy)
                      +m_ctts[i-1][j  ][kk-1][KFi-1][KFf-1]*(1-h)*hy
                      +m_ctts[i  ][j  ][kk-1][KFi-1][KFf-1]*h*hy;
        }
        if(m_KeyQCD != 0)
          for(int kk=1; kk<=m_poinQ; kk++){
            m_QCDcorR[kk-1]= m_stts[i-1][kk-1][KFi-1][KFf-1]*(1-h)
                            +m_stts[i  ][kk-1][KFi-1][KFf-1]*h;
          }
    } else if( (ww >= m_WmaxLEP2) && (ww <= m_WmaxNLC) ) {
    //  above LEP2, ILC and LHC
        x= (ww-m_WmaxLEP2)/(m_WmaxNLC-m_WmaxLEP2);
        i= min( int(trunc(m_poin4*x))+1, m_poin4);
        h= x*m_poin4-double(i-1);
        y= (1+CosThe)/2;
        j= min( int(trunc(m_poTh4*y))+1, m_poTh4);
        hy= y*m_poTh4-double(j-1);
        for(int kk=1; kk<=m_poinG; kk++){
           m_GSW[kk-1]=m_clcs[i-1][j-1][kk-1][KFi-1][KFf-1]*(1-h)*(1-hy)
                      +m_clcs[i  ][j-1][kk-1][KFi-1][KFf-1]*h*(1-hy)
                      +m_clcs[i-1][j  ][kk-1][KFi-1][KFf-1]*(1-h)*hy
                      +m_clcs[i  ][j  ][kk-1][KFi-1][KFf-1]*h*hy;
        }
       if(m_KeyQCD != 0)
          for(int kk=1; kk<=m_poinQ; kk++){
            m_QCDcorR[kk-1]= m_slcs[i-1][kk-1][KFi-1][KFf-1]*(1-h)
                            +m_slcs[i  ][kk-1][KFi-1][KFf-1]*h;
          }
    } else {
    	cout<< "KKdizet::InterpoGSW: exit, out or range ww="<<ww<<endl;
    	exit(99);
    }// if ww

}//InterpoGSW

////////////////////////////////////////////////////////////////////////
// temporary interface to CEEX/GPS
//______________________________________________________________________
void KKdizet::GetGSWxy(double GSW_re[], double GSW_im[]){
	for(int i=0; i<10; i++){
		GSW_re[i] = m_GSW[i].real();
		GSW_im[i] = m_GSW[i].imag();
		//cout<<"KKdizet::GetGSWxy:  GSW_re[0], GSW_im[0]="<<GSW_re[0]<<"  "<< GSW_im[0]<<endl;
	}
}//GetGSWxy

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
// For testing EW tables
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
/*
//__________________________________________________________________
void KKdizet::ImportEWtabs(){
///////////////////////////////////////////////////////////////////
// for tests
// transferring EW and QCD form-factors from f77 dizet common block
///////////////////////////////////////////////////////////////////
hhdizet_getsmpar_(D_swsq, D_GamZ, D_MW, D_GammW);

int r;
for(int KFi = 1; KFi <= 2; KFi++){
	  for(int KFf = 13; KFf <= 14; KFf++){
		  r=1; // LEP1 range
		  for(int i=0; i<= m_poin1;i++){ // loop over energy points
			  // Electroweak
			  for(int j=0; j<= m_poTh1;j++){ // loop over theta points (1 point !!!)
				  for(int k=1; k<=m_poinG; k++){ // loop over form-factors
					  m_cyys[i][k-1][KFi-1][KFf-1].real() = hhdizet_getelw_(r,0,KFi,KFf,i+1,j+1,k);
					  m_cyys[i][k-1][KFi-1][KFf-1].imag() = hhdizet_getelw_(r,1,KFi,KFf,i+1,j+1,k);
				  }//k
			  }//j
			  // QCD
			  for(int k=1; k<=m_poinQ; k++){ // loop over form-factors
				  m_syys[i][k-1][KFi-1][KFf-1] = hhdizet_getqcd_(r,KFi,KFf,i+1,k);
			  }//k
		  }//i
		  r=2; // Near Z resonance
		  for(int i=0; i<= m_poin2;i++){ // loop over energy points
			  // Electroweak
			  for(int j=0; j<= m_poTh2;j++){ // loop over theta points
				  for(int k=1; k<=m_poinG; k++){ // loop over form-factors
					  m_czzs[i][j][k-1][KFi-1][KFf-1].real() = hhdizet_getelw_(r,0,KFi,KFf,i+1,j+1,k);
					  m_czzs[i][j][k-1][KFi-1][KFf-1].imag() = hhdizet_getelw_(r,1,KFi,KFf,i+1,j+1,k);
				  }//k
			  }//j
			  // QCD
			  for(int k=1; k<=m_poinQ; k++){ // loop over form-factors
				  m_szzs[i][k-1][KFi-1][KFf-1] = hhdizet_getqcd_(r,KFi,KFf,i+1,k);
			  }//k
		  }//i
		  r=3; // LEP2
		  for(int i=0; i<= m_poin3;i++){ // loop over energy points
			  // Electroweak
			  for(int j=0; j<= m_poTh3;j++){ // loop over theta points
				  for(int k=1; k<=m_poinG; k++){ // loop over form-factors
					  m_ctts[i][j][k-1][KFi-1][KFf-1].real() = hhdizet_getelw_(r,0,KFi,KFf,i+1,j+1,k);
					  m_ctts[i][j][k-1][KFi-1][KFf-1].imag() = hhdizet_getelw_(r,1,KFi,KFf,i+1,j+1,k);
				  }//k
			  }//j
			  // QCD
			  for(int k=1; k<=m_poinQ; k++){ // loop over form-factors
				  m_stts[i][k-1][KFi-1][KFf-1] = hhdizet_getqcd_(r,KFi,KFf,i+1,k);
			  }//k
		  }//i
		  r=4; // above LEP2
		  for(int i=0; i<= m_poin4;i++){ // loop over energy points
			  // Electroweak
			  for(int j=0; j<= m_poTh4;j++){ // loop over theta points
				  for(int k=1; k<=m_poinG; k++){ // loop over form-factors
					  m_clcs[i][j][k-1][KFi-1][KFf-1].real() = hhdizet_getelw_(r,0,KFi,KFf,i+1,j+1,k);
					  m_clcs[i][j][k-1][KFi-1][KFf-1].imag() = hhdizet_getelw_(r,1,KFi,KFf,i+1,j+1,k);
				  }//k
			  }//j
			  // QCD
			  for(int k=1; k<=m_poinQ; k++){ // loop over form-factors
				  m_slcs[i][k-1][KFi-1][KFf-1] = hhdizet_getqcd_(r,KFi,KFf,i+1,k);
			  }//k
		  }//i
	  }//KFf
}//KFi
}// ImportEWtabs
*/
