///////////////////////////////////////////////////////////////////////////////
//     Class   KKdbase  of static input data including physics constants     //
///////////////////////////////////////////////////////////////////////////////
//                NAMING CONVENTION
// The class has no local variables hence prefix m_ for data members is
// not necessary. Moreover in all object refering to data members of THIS class
// prefix like DB->CMSene will be there anyway.
///////////////////////////////////////////////////////////////////////////////
#include "KKdbase.h"

ClassImp(KKdbase);


KKdbase::KKdbase()
{
  // This constructor is for ROOT streamers ONLY
  // all pointers has to be NULLed
  cout<< "----> KKdbase Default Constructor (for ROOT only) "<<endl;
  m_Out= NULL;
}

///_____________________________________________________________
KKdbase::KKdbase(ofstream *OutFile)
{
  cout<< "----> KKdbase USER Constructor "<<endl;
  m_Out = OutFile;
}//KKdbase

///______________________________________________________________________________________
KKdbase::~KKdbase()
{
  //Explicit destructor
  cout<< "----> KKdbase::KKdbase !!!! DESTRUCTOR !!!! "<<endl;
}///destructor


///______________________________________________________________________________________
void KKdbase::Initialize(double xpar[] )
// c++ indexing in xpar
{
  cout  << "----> KKdbase::Initialize, Entering "<<endl;
//=================================================================
// BX*** macros are in MCdev/BXFORMAT.h
  BXOPE(*m_Out);
  BXTXT(*m_Out,"========================================");
  BXTXT(*m_Out,"======    KKdbase::Initialize     ======");
  BXTXT(*m_Out,"========================================");
  ///////////////////////////////////////////////////
  m_xpar[0] = 0;
  for(int i=1; i< maxPar; i++)
  	  m_xpar[i] = xpar[i]; // f77 indexing in m_xpar

  BXTXT(*m_Out,"========== BASIC  INPUT ================");
  CMSene = m_xpar[ 1];
  BX1F(*m_Out,"   CMSene", CMSene, " Center of the Mass Energy    [ 1]     =");

  BXTXT(*m_Out,"========== PHYSICS SWITCHES ============");
  KeyISR = m_xpar[20];
  KeyFSR = m_xpar[21];
  KeyINT = m_xpar[27];
  KeyGPS = m_xpar[28];
  KeyWgt = m_xpar[10];
  KeyElw = m_xpar[12];
  KeyPia = m_xpar[22];
  BX1I(*m_Out,"   KeyISR", KeyISR,    " ISR initial state emission       [20]          =");
  BX1I(*m_Out,"   KeyFSR", KeyFSR,    " FSR final   state emission       [21]          =");
  BX1I(*m_Out,"   KeyINT", KeyINT,    " IFI interference on/off switch   [27]          =");
  BX1I(*m_Out,"   KeyGPS", KeyGPS,    " GPS/CEEX matrix element type     [28]          =");
  BX1I(*m_Out,"   KeyWgt", KeyWgt,    " WTed/unWTed events  switch       [10]          =");
  BX1I(*m_Out,"   KeyElw", KeyElw,    " Electroweak correction  switch   [12]          =");
  KeyFixAlf = m_xpar[3032];
  BX1I(*m_Out,"KeyFixAlf", KeyFixAlf, " Fixed/Running alpha QED          [3032]        =");
  KeyMasQ = m_xpar[3030];         // constituent vs current quarks in kinematics
  BX1I(*m_Out,"  KeyMasQ", KeyMasQ,   " Constituent/current quarks mass  [3030]        =");

  // default values
  KeyZet  = m_xpar[501];
  BX1I(*m_Out,"   KeyZet",  KeyZet,   " Z resonance features             [501]         =");
  MZ      = m_xpar[502];
  GamZ    = m_xpar[504];
  swsq    = m_xpar[503];
  MW      = m_xpar[505];
  GamW    = m_xpar[506];
  GFermi  = m_xpar[ 32];
  BXTXT(*m_Out,"=======ELECTROWEAK SPECIALS ============");
  BXTXT(*m_Out,"==========  DEFAULT ====================");
  BX1F(*m_Out,"       MZ",    MZ,  " Mass  of Z bozon GeV             [502] =");
  BX1F(*m_Out,"     GamZ",  GamZ,  " Width of Z bozon GeV             [504] =");
  BX1F(*m_Out,"     swsq",  swsq,  " Electroweak mixing angle         [503] =");
  BX1F(*m_Out,"       MW",    MW,  " Mass  of W bozon GeV             [502] =");
  BX1F(*m_Out,"     GamW",  GamW,  " Width of W bozon GeV             [504] =");
  // input of special benchmarks EW parameter overrides default
  if (KeyElw <= 0) {
      if (m_xpar[3502] != 0e0)  MZ     = m_xpar[3502];
      if (m_xpar[3503] != 0e0)  swsq   = m_xpar[3503];
      if (m_xpar[3504] != 0e0)  GamZ   = m_xpar[3504];
      if (m_xpar[3505] != 0e0)  MW     = m_xpar[3505];
      if (m_xpar[3506] != 0e0)  GamW   = m_xpar[3506];
      if (m_xpar[3532] != 0e0)  GFermi = m_xpar[3532];
  }// m_KeyElw
  // additional rescaling
  if (KeyZet == -2) {
      MZ   = MZ  /sqrt(1 + sqr(GamZ/MZ));
      GamZ = GamZ/sqrt(1 + sqr(GamZ/MZ));
      MZ   = MW  /sqrt(1 + sqr(GamW/MW));
      GamW = GamW/sqrt(1 + sqr(GamW/MW));
      swsq = 1 - sqr(GamW/GamZ); //   ! sin^2(Theta_W)
  }// m_KeyZet= -2
  if (KeyElw <= 0 || KeyZet <= 0) {
  BXTXT(*m_Out,"===== SPECIAL BENCHMARK ================");
  BX1F(*m_Out,"       MZ",    MZ,  " Mass  of Z bozon GeV             [502] =");
  BX1F(*m_Out,"     GamZ",  GamZ,  " Width of Z bozon GeV             [504] =");
  BX1F(*m_Out,"     swsq",  swsq,  " Electroweak mixing angle         [503] =");
  BX1F(*m_Out,"       MW",    MW,  " Mass  of W bozon GeV             [502] =");
  BX1F(*m_Out,"     GamW",  GamW,  " Width of W bozon GeV             [504] =");
  }

  BXTXT(*m_Out,"==== MONTE CARLO ALGORITHM PARAMS ======");
  WTmax  = m_xpar[ 9];
  vvmin  = m_xpar[16];                     // minimum v, infrared cut
  vvmax  = m_xpar[17];                     // default vvmax
  delfac = m_xpar[18];
  BX1F(*m_Out,"    WTmax", WTmax,  " Maximum weight in rejection      [ 9] =");
  BX1F(*m_Out,"    vvmin", vvmin,  " Infrared cutoff in photon energy [16] =");
  BX1F(*m_Out,"    vvmax", vvmax,  " Maximum photon energy            [17] =");
  BX1F(*m_Out,"   delfac", delfac, " FSR photon IR cut reduction      [18] =");
  Xenph  = xpar[40];
  if(KeyINT == 0)  Xenph  = 1;
  BX1F(*m_Out,"    Xenph", Xenph,  " Photon emission enhancement      [40] =");
  KeyWtm = xpar[26];
  BX1I(*m_Out,"   KeyWtm",  KeyWtm, " Photon emission without mass terms [26]        =");

  Vcut[0]= m_xpar[41] ;
  Vcut[1]= m_xpar[42] ;
  Vcut[2]= m_xpar[43] ;
  BX1F(*m_Out,"  Vcut[0]", Vcut[0]," Techn. parameter in CEEX/GPS   [41]   =");
  BX1F(*m_Out,"  Vcut[1]", Vcut[1]," Techn. parameter in CEEX/GPS   [42]   =");
  BX1F(*m_Out,"  Vcut[2]", Vcut[2]," Techn. parameter in CEEX/GPS   [43]   =");
  //
  VQcut = m_xpar[3015];
  BX1F(*m_Out,"    VQcut",   VQcut," Parton shower quark cutoff   [3015]   =");
  XXXmin = m_xpar[3003];
  XXXmax = m_xpar[3004];
  BX1F(*m_Out,"   XXXmin",  XXXmin," Techn. param. init. value    [3003]   =");
  BX1F(*m_Out,"   XXXmax",  XXXmax," Techn. param. init. value    [3004]   =");
  Nalgo  = xpar[3028];  // presently not used
  BX1I(*m_Out,"    Nalgo",  Nalgo,    " Kinematic algorithm type (unused)[3028]        =");


  BXTXT(*m_Out,"========= PHYSICS DATA/CONSTANTS =======");
  MasPhot= m_xpar[510];
  Alfinv0 = m_xpar[30];
  AlfinvZ = m_xpar[808];

  for(int KF=1; KF<20;KF++){ // f77 indexing!
    IsGenerated[KF] = xpar[400+KF];    // Generation flag
    Nc[KF]     = xpar[500+10*KF+2];    // color
    Qf[KF]     = xpar[500+10*KF+3]/3;  // electric charge
    T3f[KF]    = xpar[500+10*KF+4]/2;  // izospin, L-hand component
    fmcon[KF]  = xpar[500+10*KF+5];    // constituent mass
    fmass[KF]  = xpar[500+10*KF+6];    // current fermion mass
	if(IsGenerated[KF] ==1 )
	   BX1I(*m_Out,"       KF",     KF,    " Active fermion in the initial or final state   =");
  }// for j

  BX1F(*m_Out,"  MasPhot", MasPhot," Photon mass, IR regulator       [510] =");
  BX1F(*m_Out,"  Alfinv0", Alfinv0," 1/alpha QED at q^2 = 0          [ 30] =");
  BX1F(*m_Out,"  AlfinvZ", Alfinv0," 1/alpha QED at MZ               [808] =");
  gnanob  = m_xpar[31];
  BX1F(*m_Out,"   gnanob", gnanob, " Conversion GeV^(-2) to nanobarns [31] =");

  BXCLO(*m_Out);
  cout  << "----> KKdbase::Initialize, exiting "<<endl;
}//Initialize

