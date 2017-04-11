///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                Class KKsem                                                //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include "KKsem.h"

///////////////////////////////////////////////////////////////////////////////////
void KKsem::Initialize(TFile &DiskFileA){
  //------------------------------------------------------------------------
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  cout<<"==================================================================="<<endl;
  cout<<"================ KKsem::initialize begin ==========================="<<endl;
  const int jmax =10000;
  double ypar[jmax];
  for(int j=1; j<=jmax; j++)
    ypar[j-1]=HST_KKMC_NORMA->GetBinContent(j);    // xpar encoded
  //[[[[[
  cout<<"***********************KKsem::initialize****************************"<<endl;
  for(int j=0;j<30;j++)
    cout<<j+1<<"   "<<ypar[j]<<endl;
  //]]]]]
  char *output_file = "./kksem.output";
  long stl2 = strlen(output_file);
  long mout =16;
  kk2f_fort_open_(mout,output_file,stl2);
  kk2f_initialize_(ypar);
  kksem_initialize_(ypar);
  cout<<"================ KKsem_initialize END   ==========================="<<endl;
  cout<<"==================================================================="<<endl;
  //long kdum; kk2f_getkeyfsr_(kdum); //<-- to avoid linker bug! (if no kk2f_initialize)
}

///////////////////////////////////////////////////////////////////////////////////
void KKsem::VVplot( TH1 *hstNew, long KF, char chak[5], long KeyDis, long KeyFob)
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
    cout<<"KKsem::VVplot: ib= "<<ib<<"  Bin(ib) =  "<<Bin[ib]<<endl;
    hstNew->SetBinContent(ib+1, Bin[ib]);
    hstNew->SetBinError(  ib+1, 0.0);
  }
}// KKsem::VVplot


///////////////////////////////////////////////////////////////////////////////////
void KKsem::Cplot( TH1 *hstNew, 
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
    cout<<"KKsem::Cplot: ib= "<<ib<<" cBin(ib) =  "<<Dsig<<endl;
    hstNew->SetBinContent(ib+1, Dsig);
    hstNew->SetBinError(  ib+1, 0.0);
  }
  kksem_setcrange_(-1.0,1.0); // back to normal
}// KKsem::VVplot


void KKsem::PlInitialize(FILE *ltx, int lint)
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

void KKsem::PlEnd(FILE *ltex)
{//---------------------------------------------------
// Note that TeX file is used in \input then you may not want
// to have header and \end{document}
if( m_lint |= 1){
   fprintf(ltex,"\\end{document} \nl");
}
}//GLK_PlEnd

//  SUBROUTINE GLK_PlTable2(Npl,idl,ccapt,mcapt,fmt,chr1,chr2,chr3)
void KKsem::PlTable2(int Npl, FILE *ltex, Char_t *Capt[], const char *chr1)
//* Tables in TeX, up to 9 columns
//* Npl           = numbers of columns/histograms
//* idl(1:Npl)    = list of histo id's
//* ccapt(1:Npl+1)= list of column-captions above each column
//* mcapt         = multicolumn header, none if mcapt=' ',
//* fmt(1:1)      = format to print x(i) in first columb,
//*                 h(i) and error he(i) in further columns
//* chr1          = ' ' normal default, ='S' the Same table continued
//* chr2          = ' ' midle of the bin for x(i) in the first column
//*               = 'R' right edge,     ='L' left edge of the bin
//* chr3          = ' ' no page eject,  ='E' with page eject at the end.
//* Furthermore:
//* Captions are defined by means of
//*    CALL GLK_PlCapt(capt) before CALL GLK_PlTable2
//*    where CHARACTER*80 capt(50) is content of
//*    caption, line by line, see also comments in GLK_PlCapt routine.
{

  //fprintf(ltex,"%s \n", Capt[1]);
  //fprintf(ltex,"%s \n", Capt[2]);

  //if( chr1 == "S" ) fprintf(ltex,"ch1= %s \n", chr1);

  if( chr1 == " "){
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
	  // Central Caption
	  //------------------------------!
	  //if(ABS(m_lint) |= 2 ) {
	  //  fprintf(ltex,"\\caption{\\footnotesize\\sf'
	  //  for(int i=1, i<=m_KeyTit; i++){
	  //  fprintf(ltex,"%s \n",    m_titch(i));
	  //}
	  //fprintf(ltex," \n");

	  //------------------------------!
	  // Tabular header
	  //------------------------------!
	  //WRITE(m_ltx,'(20A)') m_BS,'begin{tabular} {|',  ('|r',j=1,Npl+1),  '||}'
	  //WRITE(m_ltx,'(4A)') m_BS,'hline',m_BS,'hline'
      fprintf(ltex,"\\begin{tabular}{|");
	  for(int i=1; i<=Npl; i++ ) fprintf(ltex,"|r");
	  fprintf(ltex,"||}\n");
	  fprintf(ltex,"\\hline\\hline\n");
	  //------------------------------!
	  // Captions in columns
	  //------------------------------!
	  //WRITE(m_ltx,'(2A)') ccapt(1),('&',ccapt(j+1),j=1,Npl)
	  fprintf(ltex,"%s  \n", Capt[1]);
	  for(int i=2; i<=Npl; i++ ) fprintf(ltex,"& %s \n", Capt[i]);
  }//chr1
  fprintf(ltex,"\\\\ \\hline\n");




  //------------------------------!
  // Ending
  //------------------------------!
  fprintf(ltex,"\\hline\n");
  fprintf(ltex,"\\end{tabular}\n");
  if(abs(m_lint) == 2 ){
       fprintf(ltex,"% ========================================\n");
  } else {
       fprintf(ltex,"\\end{table}\n");
  }
  fprintf(ltex,"% ============= end   table ==============\n");
  fprintf(ltex,"% ========================================\n");

}//GLK_PlTable2


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                End of Class KKsem                                         //
//////////////////////////////////////////////////////////////////////////////

