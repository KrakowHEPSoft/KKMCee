//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//               Class   KKabox                                             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
// KKabox is multipurpose toolbox for KKMC testing.
//  1. Interfaces (erappers) to KKMC and KKsem F77 subrograms
//  2. Integrand for Foam
//  3. A few routines for producing tatex table out of histograms
//////////////////////////////////////////////////////////////////////////////

#include "KKabox.h"

KKabox::KKabox(const char* Name)
{
		cout<< "----> KKabox USER Constructor "<<endl;
}

///////////////////////////////////////////////////////////////////////////////////
void KKabox::Initialize(TFile &DiskFileA){
  //------------------------------------------------------------------------
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  cout<<"==================================================================="<<endl;
  cout<<"================ KKsem initialization begin ==========================="<<endl;
  m_jmax =10000;
  for(int j=1; j<=m_jmax; j++)
    m_ypar[j-1]=HST_KKMC_NORMA->GetBinContent(j);    // xpar encoded
  //[[[[[

  m_CMSene  = m_ypar[ 1 -1];
  m_vvmax   = m_ypar[17 -1];

  m_alfinv   = m_ypar[30 -1];     // 1/alphaQED at Q^2=0
  m_gnanob   = m_ypar[31 -1];     // GeV^2 -> nanobarns

  m_pi      = 3.1415926535897932;
  m_ceuler  = 0.57721566;
//
  m_alfpi   = 1/m_alfinv/m_pi;
//
  m_beam    = 0.510999e-3;  // electron
  m_chini   = 1.0;          // electron

  m_fin     = 0.105;        // final ferm. muon

  m_KFini   = 11; // electron
  m_KFf     = 13; // muon

  m_kDim    =    3;         // No. of dim. for Foam, =2,3 Machine energy spread OFF/ON
  m_nCells  = 2000;         // No. of cells, optional, default=2000
  m_nSampl  =  200;         // No. of MC evts/cell in exploration, default=200

  m_KeyISR  = 2;            // Type of ISR/QED switch, 0,1,2

  m_Mode    = 3;            // Operation mode for Foam

  m_count   = 0;            // counter for debug

  cout<<"***********************KKsem::initialize****************************"<<endl;
  for(int j=0;j<30;j++)
    cout<<j+1<<"   "<<m_ypar[j]<<endl;
  //]]]]]
  char *output_file = "./kksem.output";
  long stl2 = strlen(output_file);
  long mout =16;
  kk2f_fort_open_(mout,output_file,stl2);
  kk2f_initialize_(m_ypar);
  kksem_initialize_(m_ypar);
  cout<<"================ KKsem initialization END   ==========================="<<endl;
  cout<<"==================================================================="<<endl;
  //long kdum; kk2f_getkeyfsr_(kdum); //<-- to avoid linker bug! (if no kk2f_initialize)
}

///////////////////////////////////////////////////////////////////////////////////
void KKabox::VVplot( TH1 *hstNew, long KF, char chak[5], long KeyDis, long KeyFob)
{
  long   nbin = hstNew->GetNbinsX();
  double xmin = hstNew->GetXaxis()->GetXmin();
  double xmax = hstNew->GetXaxis()->GetXmax();
  double Bin[nbin];
  //
  kksem_setkffin_(KF); // set m_KFfin in KKsem
  bornv_setkf_( KF );  // set SINGLE Final State
  //
  kksem_setkeyfob_(KeyFob);
  kksem_vvplot_vec_(KeyDis,chak,nbin,xmin,xmax,Bin);
  //
  hstNew->Reset();
  for(int ib=0;ib<nbin;ib++){
    cout<<"KKabox::VVplot: ib= "<<ib<<"  Bin(ib) =  "<<Bin[ib]<<endl;
    hstNew->SetBinContent(ib+1, Bin[ib]);
    hstNew->SetBinError(  ib+1, 0.0);
  }
}// KKabox::VVplot



///////////////////////////////////////////////////////////////////////////////////
void KKabox::Cplot( TH1 *hstNew,
		   long KF, char chak[5], long KeyDis, long KeyFob, double vmin, double vmax)
{
  long   nbc = hstNew->GetNbinsX();
  double cmin = hstNew->GetXaxis()->GetXmin();
  double cmax = hstNew->GetXaxis()->GetXmax();
  double cBin[nbc];
  long   nbv = 1;   // only one bin in v needed
  double vBin[100]; // only one bin in v needed
  //
  kksem_setkffin_(KF); // set m_KFfin in KKsem
  bornv_setkf_( KF );  // set SINGLE Final State
  //
  if( (KeyFob == 0) ||(KeyFob == -100)  ){
    kksem_setkeyfob_(KeyFob);
  }else{
    cout<<"++++ KKsem::Cplot wrong KeyFob= "<< KeyFob<<endl;
    exit(7);
  }
  //
  double c1,c2,dc,fnor,Dsig;
  dc = (cmax-cmin)/nbc;
  hstNew->Reset();
  for(int ib=0;ib<nbc;ib++){
    c1= cmin +ib*dc;
    c2= c1+dc;
    kksem_setcrange_(c1,c2);
    kksem_vvplot_vec_(KeyDis,chak,nbv,vmin,vmax,vBin); // vBin has only one bin!
    Dsig=vBin[0]/dc*(vmax-vmin); // integration over v!
    cout<<"KKabox::Cplot: ib= "<<ib<<" cBin(ib) =  "<<Dsig<<endl;
    hstNew->SetBinContent(ib+1, Dsig);
    hstNew->SetBinError(  ib+1, 0.0);
  }
  kksem_setcrange_(-1.0,1.0); // back to normal
}// KKabox::VVplot


void KKabox::PlInitialize(FILE *ltx, int lint)
{
//----------------------------------------------------------------------
// Lint =0     Normal mode, full LaTeX header
// Lint =1     For TeX file is used in \input, no  LaTeX header
// Lint =2     LaTeX header for one-page plot used as input for postscript
// Negative Lint only for debug, big frame around plot is added.
//----------------------------------------------------------------------
   m_lint=lint;
if( abs(lint) == 0){
// Normal mode, no colors!!!
   fprintf(ltx,"\\documentclass[12pt]{article}\n");
   fprintf(ltx,"\\textwidth  = 16cm\n");
   fprintf(ltx,"\\textheight = 18cm\n");
   fprintf(ltx,"\\begin{document}\n");
   fprintf(ltx,"  \n");
} else if( abs(lint) == 1) {
// For TeX file is used in \input
   fprintf(ltx,"  \n");
} else if( abs(lint) == 2){
// For one-page plot being input for postscript
   fprintf(ltx,"\\documentclass[12pt,dvips]{article}\n");
   fprintf(ltx,"\\usepackage{amsmath}\n");
   fprintf(ltx,"\\usepackage{amssymb}\n");
   fprintf(ltx,"\\usepackage{epsfig}\n");
   fprintf(ltx,"\\usepackage{epic}\n");
   fprintf(ltx,"\\usepackage{eepic}\n");
   fprintf(ltx,"\\usepackage{color}\n"); //<-for colors!!!
   fprintf(ltx,"\\begin{document}\n");
   fprintf(ltx,"\\pagestyle{empty}\n");
   fprintf(ltx,"  \n");
} else {
   cout<<"+++STOP in GLK_PlInt, wrong lint =" <<lint<< endl;
}// lint
}//GLK_PlCap

void KKabox::PlEnd(FILE *ltex)
{//---------------------------------------------------
// Note that TeX file is used in \input then you may not want
// to have header and \end{document}
if( m_lint |= 1){
   fprintf(ltex,"\\end{document} \nl");
}
}//GLK_PlEnd

void KKabox::PlTable2(int Ncol, TH1D *iHst[], FILE *ltex, Char_t *Capt[], Char_t Mcapt[] , const char *chr1, int k1,int k2,int dk)
//* Tables in TeX, up to 9 columns
//* Ncol          = numbers of columns/histograms
//* idl(1:Npl)    = list of histo id's
//* ccapt(1:Npl+1)= list of column-captions above each column
//* mcapt         = multicolumn header, none if mcapt=' ',
//* chr1          = ' ' normal default, = Header+Table+Ending
//*               = 'B' no page eject,  = Header+Table
//*               = 'E' no page eject,  =        Table+Ending
//*               = 'E' no page eject,  =        Table
//* k1,k2,dk      = range of bins is (k1,k2) with increment dk
{
  int Npl=abs(Ncol);
  if( chr1 == " " || chr1 == "B"){
	  //------------------------------!
	  //           Header
	  //------------------------------!
	  fprintf(ltex," \n");
	  fprintf(ltex,"% ========================================\n");
	  fprintf(ltex,"% ============ begin table ===============\n");
	  //
	  if (abs(m_lint) == 2 ){
	      fprintf(ltex,"\\noindent\n");
	  } else {
	      fprintf(ltex,"\\begin{table}[!ht] \n");
	      fprintf(ltex,"\\centering \n");
	  }

	  //------------------------------!
	  // Tabular header
	  //------------------------------!
	  //WRITE(m_ltx,'(20A)') m_BS,'begin{tabular} {|',  ('|r',j=1,Npl+1),  '||}'
	  //WRITE(m_ltx,'(4A)') m_BS,'hline',m_BS,'hline'
      fprintf(ltex,"\\begin{tabular}{|");
	  for(int i=0; i<=Npl; i++ ) fprintf(ltex,"|r");
	  fprintf(ltex,"||}\n");
	  fprintf(ltex,"\\hline\\hline\n");

	  //------------------------------!
	  // Captions in columns
	  //------------------------------!
	  //WRITE(m_ltx,'(2A)') ccapt(1),('&',ccapt(j+1),j=1,Npl)
	  fprintf(ltex,"%s  \n", Capt[0]);
	  for(int i=1; i<=Npl; i++ ) fprintf(ltex,"& %s \n", Capt[i]);

	  fprintf(ltex,"\\\\ \\hline\n");
  }

  //------------------------------------------!
  // Optional Multicolumn caption
  //------------------------------------------!
  if(Ncol>0){
     fprintf(ltex,"& \\multicolumn{ %i }{c||}{",Npl);
     fprintf(ltex,"  %s  } \\\\   \\hline\n", Mcapt);
  }//Mcapt

  // X range taken from 1-st histogram
  int      nbX  = (iHst[1])->GetNbinsX();
  Double_t Xmax = (iHst[1])->GetXaxis()->GetXmax();
  Double_t Xmin = (iHst[1])->GetXaxis()->GetXmin();
  double dx = (Xmax-Xmin)/nbX;
  cout<<"  nbX=  " <<nbX<<endl;
  cout<<"  Xmin= " <<Xmin<<endl;
  cout<<"  Xmax= " <<Xmax<<endl;
  // Raws

  double xk, yi, ei;
  for(int k=k1; k<=k2; k+=dk){    // loop over bin (raw) number
	cout<<" k="<<k<<endl;
	xk = Xmin +dx*k; // right edge
    fprintf(ltex,"$  %10.2f $", xk);
    for( int j=1; j<=Npl; j++ ){
	   yi = (iHst[j])->GetBinContent(k);
	   ei = (iHst[j])->GetBinError(k);
	   fprintf(ltex," & $ %10.4f \\pm %8.4f $ ", yi, ei);
     }//j
    fprintf(ltex,"\\\\ \n");
 }//k
 fprintf(ltex,"\\hline\n");

  //------------------------------!
  // Ending
  //------------------------------!
  if( chr1 == " " || chr1 == "E"){
	  fprintf(ltex,"\\end{tabular}\n");
	  fprintf(ltex,"% ============ end   table ===============\n");
	  fprintf(ltex,"% ========================================\n");
	  fprintf(ltex," \n");
  }//chr1

}//GLK_PlTable2

///------------------------------------------------------------------------
double KKabox::gamISR( double svar){
	  return  sqr(m_chini)*2*m_alfpi*( log(svar/sqr(m_beam)) -1);
}

///------------------------------------------------------------------------
double KKabox::gamFSR( double svar){
	  return              2*m_alfpi*( log(svar/sqr(m_fin)) -1);
}

///------------------------------------------------------------------------
double KKabox::Rho_isr(double svar, double vv){
/// ISR rho-function for ISR

  double alf1   = m_alfpi;
  double gami   = gamISR(svar);
  //gami = sqr(m_chini)*2*m_alfpi*( log(svar/sqr(m_beam)) -1);
///
  double gamfac = exp(-m_ceuler*gami)/TMath::Gamma(1+gami);
  double delb   = gami/4 +alf1*(-0.5  +sqr(m_pi)/3.0);
  double ffact  = gamfac*exp(delb);

  double rho,dels,delh;
  if(       m_KeyISR == 0){
/// zero   order exponentiated
	dels = 0;
	delh = 0;
	rho  = ffact*gami* exp( log(vv)*(gami-1) ) *(1 +dels +delh);
  }else if( m_KeyISR == 1){
/// first  order
	dels = gami/2;   /// NLO part =0 as for vector boson???
    delh = vv*(-1 +vv/2);
    rho = ffact*gami* exp( log(vv)*(gami-1) ) *(1 +dels +delh);
  }else if( m_KeyISR == 2){
/// second order without NLO part
    dels = gami/2 +sqr(gami)/8;
    delh = vv*(-1+vv/2.0)
          +gami*0.5*(-0.25*(4.0-6.0*vv+3.0*vv*vv)*log(1-vv)-vv);
    rho = ffact*gami* exp( log(vv)*(gami-1) ) *(1 +dels +delh);
  }else{
	  cout<<"+++++TMCgenH::KKdistr: Wrong KeyISR = " << m_KeyISR<<endl;
	  exit(5);
  }
///
  return rho;
}//Rho_isr



///------------------------------------------------------------------------
double KKabox::Rho_fsr(double svar, double uu){
/// ISR+FSR rho-function

  double alf1   = m_alfpi;
  double gamf   = gamFSR(svar*(1-uu));
///
  /////delb   = Chf2* alf1*(0.5d0*bilg -1d0  +pi**2/3d0)
  /////delb   = delb -betf/2 *dlog(1-uu)
  double gamfac = exp(-m_ceuler*gamf)/TMath::Gamma(1+gamf);
  double delb   = gamf/4 +alf1*(-0.5  +sqr(m_pi)/3.0)
		         -gamf/2 *log(1-uu);
  double ffact  = gamfac*exp(delb);

  double rho,dels,delh;
  if(       m_KeyISR == 0){
/// zero   order exponentiated
	dels = 0;
	delh = 0;
	rho  = ffact*gamf* exp( log(uu)*(gamf-1) ) *(1 +dels +delh);
  }else if( m_KeyISR == 1){
/// first  order
	dels = gamf/2;   /// NLO part =0 as for vector boson???
    delh = uu*(-1 +uu/2);
    rho = ffact*gamf* exp( log(uu)*(gamf-1) ) *(1 +dels +delh);
  }else if( m_KeyISR == 2){
/// second order without NLO part
    /////dels  = betf/2d0 +betf**2/8d0
    /////delh  = uu*(-1d0+uu/2d0)
    /////$        +betf*(-0.5d0*uu-0.25d0*uu*(-1d0+0.5d0*uu)*log(1d0-uu))
    dels = gamf/2 +sqr(gamf)/8;
    delh = uu*(-1+uu/2.0)
          +gamf*0.5*( -0.5*uu -0.25*uu*(-1.0 +0.5*uu)*log(1-uu));
    rho = ffact*gamf* exp( log(uu)*(gamf-1) ) *(1 +dels +delh);
  }else{
	  cout<<"+++++TMCgenH::KKdistr: Wrong KeyISR = " << m_KeyISR<<endl;
	  exit(5);
  }
///
  return rho;
}//Rho_fsr

///////////////////////////////////////////////////////////////
Double_t KKabox::Density(int nDim, Double_t *Xarg)
{ // density distribution for Foam
    if ( abs(m_Mode) == 3 )
    	Density3(nDim, Xarg);
    else {
    	cout<<"+++ KKabox::Density: Wrong m_Mode = "<< m_Mode<<endl;
    	exit(20);
    }
}// Density

///////////////////////////////////////////////////////////////
Double_t KKabox::Density3(int nDim, Double_t *Xarg)
{ // density distribution for Foam
	m_count++;  // counter for debug
	//
	Double_t Dist=1;

	double svar = sqr(m_CMSene);
	double svarCum = svar;

// ******** mapping for ISR *******
	double gamiCR,gami,alfi;
	double CMSene1= sqrt(svar);
	bornv_makegami_( CMSene1, gamiCR,gami,alfi);   // from KKMC
	//[[[[ debug
	//gami = gamISR(CMSene1);
	//]]]]
	// cout<<" CMSene1,gami= "<< CMSene1 <<"  "<< gami <<endl;
	double R= Xarg[0];
	m_vv = exp(1.0/gami *log(R)) *m_vvmax; // mapping
	Dist *= m_vv/R/gami ;                  // Jacobian
	if( gami < 0 )      return 0.0;    // temporary fix
    if( m_vv < 1e-200 ) return 0.0;    // temporary fix
    // ISR photonic distribution
	double Rho2 = Rho_isr(svar,m_vv);  // remember take care of m_mbeam!!!
	Dist *= Rho2;
	svarCum *= (1-m_vv);
	double svar2 = svar*(1-m_vv);

// ******** mapping for FSR *******
    double rr= Xarg[1];
    double gamf   = gamFSR(svar2);
    m_uu = exp(1.0/gamf *log(rr));     // mapping
    Dist *= m_uu/rr/gamf;              // Jacobian
// FSR photonic distribution
  	double Rho3 = Rho_fsr(svar2,m_uu);           // remember take care of m_mbeam!!!
  	Dist *= Rho3;
    svarCum *= (1-m_uu);

    // ******** mapping for polar angle *******
    m_CosTheta = -1.0 + 2.0* Xarg[2];
    Dist *= 2.0;

    double zz = (1-m_vv)*(1-m_uu);
    m_xx = 1-zz;
    m_xx = m_vv + m_uu - m_vv*m_uu;  // numerically more stable

// ******** finally Born factor *******
    long KeyFob;
    KeyFob=   10; // BornV_Dizet, with EW and without integration ???
    KeyFob=  -11; // BornV_Simple, for KeyLib=0, NO EW, NO integration OK
    KeyFob=  -10; // KKsem_BornV, NO EW, NO integration OK!
    KeyFob= -100; // KKsem_BornV, NO EW, WITH integration, OK
    KeyFob=    0; // With EW (BornV_Dizet) With integration OK!
//  -----------------
	kksem_setkeyfob_( KeyFob );

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//double xBorn;
	//Integrated Born from KKMC
	//kksem_makeborn_( svar2, xBorn);
	//Dist *= xBorn/2.0;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//    In BornV_Differential:
//    CALL BornV_InterpoGSW( ABS(KFf),  svar, CosThe)
//    Born= BornV_Dizet( 1,m_KFini,KFf, svar, CosThe, eps1,eps2,ta,tb)
//    Born = 4*pi*alfa**2/(3d0*svar )*BornY  *m_gnanob

//
	bornv_interpogsw_(m_KFf,svar2, m_CosTheta);
	double dSig_dCos = bornv_dizet_( 1, m_KFini, m_KFf, svar2, m_CosTheta, 0.0, 0.0, 0.0, 0.0);

// ******* effective masses *********
	m_Mll = sqrt(svarCum); // final after FSR
	m_Mka = sqrt(svar2);   // after ISR

// =============== Sigm/dOmega from spin amplitudes ===============
// Effective 4-momenta, KKMC convention: p={px,py,pz,E)
	double Ene = m_Mka/2;
	double Pmb  = sqrt( (Ene-m_beam)*(Ene+m_beam) ); // modulus
	m_p1 = { 0, 0 , Pmb, Ene};  // beam
	m_p2 = { 0, 0 ,-Pmb, Ene};  // beam
	double Pmf  =sqrt( (Ene-m_fin)*(Ene+m_fin) ); // modulus
	m_p3 = { Pmf*sqrt(1-sqr(m_CosTheta)), 0 , Pmf*m_CosTheta,  Ene}; // final
	m_p4 = {-Pmf*sqrt(1-sqr(m_CosTheta)), 0 ,-Pmf*m_CosTheta,  Ene}; // final
	double PX[4] = {0, 0, 0, 2*Ene};
	double dSigAngF;
	gps_bornfoam_( 0,m_KFini,m_KFf,m_Mka,m_CosTheta,dSigAngF);
//*****  Obsolete gps_bornf_()
	double dSigAng;
    gps_bornf_(m_KFini, m_KFf ,PX, m_CosTheta, m_p1,m_beam, m_p2, -m_beam,
                                               m_p3,m_fin,  m_p4, -m_fin,   dSigAng);
//*****GPS_Bornf_(KFi,KFf,PX,CosThe,p1,m1,p2,m2,p3,m3,p4,m4,Xborn)
//    double fleps = 1e-50;
//    gps_bornf_(m_KFini, m_KFf ,PX, m_CosTheta, m_p1,fleps,  m_p2, -fleps,
//    		                                     m_p3,fleps,  m_p4, -fleps,   dSigAng);
//************ Debug*** Debug*** Debug*** Debug*** Debug ***********
    if( m_count <1000 && fabs(svar/svar2-1)>0.20 ){  // debug
//    if( m_count <1000 ){  // debug
    	double Rat;
    	gps_bornfoam_( 1,m_KFini,m_KFf,m_Mka,m_CosTheta,dSigAngF);
        double dSigAngFF = gps_makerhofoam_(1.0);
//    	Rat = dSigAngF/( dSig_dCos );
    	Rat = dSigAngF/( dSigAngFF );
    	cout<<" Density debug m_count= "<< m_count<< endl;
    	cout<<" dSig_dCos  = "<< dSig_dCos;
    	cout<<" dSigAng    = "<< dSigAng;
    	cout<<" dSigAngF    = "<< dSigAngF;
    	cout<<" svar/svar2 = "<< svar/svar2;
    	cout<<" Rat = "<<Rat<<endl;
    } // end debug **********
//
	double sig0nb = 4*m_pi* sqr(1/m_alfinv)/(3.0*svar2 )*m_gnanob;
//	Dist *=  dSig_dCos *3.0/8.0 *sig0nb;
	Dist *=  dSigAngF *3.0/8.0 *sig0nb;

	return Dist;
}// Density


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                End of Class KKabox                                        //
///////////////////////////////////////////////////////////////////////////////

