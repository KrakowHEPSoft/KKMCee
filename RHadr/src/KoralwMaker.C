//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//               CLASS   KoralwMaker                                        //
//                                                                          //
//              Interface  to MC event generator Koralw                     //
//                                                                          //
//      Notes:                                                              //
//      Present version of interface to fortran                             //
//      is based on programs of Piotr Golonka (1997)                        //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////



#include "KoralwMaker.h"
#include "string.h"      // to get strlen



//////////////////////////////////////////////////////////////////////////////
//       KORALW    external FORTRAN programs                                //
//////////////////////////////////////////////////////////////////////////////
extern "C" void koralw_(long &mode, double x[]);
extern "C" void glimit_(const long &n);
extern "C" void goutpu_(const long &n);
extern "C" void fort_open_( long &lunit, char *filename, long s1);
extern "C" void fort_close_(long &lunit);
//
// readatax obsolete
extern "C" void readatax_(long &lunit, double xpar[], long &imax);
//
// clone of KW_ReaDataX, in order to avoid diskfile
extern "C" void readatax2_(long &lunit, long &ireset, long &italk, long &imax, double xpar[]);
extern "C" void varran_(double rn[], long &n);
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//=============================================================
//      Common Block from KORALW
//      COMMON / momdec / q1(4),q2(4),p1(4),p2(4),p3(4),p4(4)
//                    q1 = W+        four-momentum
//                    q2 = W-        four-momentum
//                    p1 = fermion1  four-momentum
//                    p2 = fermion2  four-momentum
//                    p3 = fermion3  four-momentum
//                    p4 = fermion4  four-momentum
//=============================================================
typedef struct {
  double  q1[4];
  double  q2[4];
  double  p1[4];
  double  p2[4];
  double  p3[4];
  double  p4[4];
} CommonMOMDEC;
extern CommonMOMDEC &momdec ;
extern CommonMOMDEC momdec_ ;
#define cb_momdec momdec_
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
//===================================================================
//    COMMON / momset / qeff1(4),qeff2(4),sphum(4),sphot(100,4),nphot
//                   qeff1 = dummy 
//                   qeff2 = dummy 
//===================================================================
typedef struct {
  double  qeff1[4];
  double  qeff2[4];
  double  sphum[4];
  double  sphot[4][100];
  long    nphot;
} CommonMOMSET;
extern CommonMOMSET &momset ;
extern CommonMOMSET momset_ ;
#define cb_momset momset_
//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////
//========================================================
// COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
//--------------------------------------------------------
typedef struct {
    int    n; 
    int    k[5][4000];
    float  p[5][4000];
    float  v[5][4000];
                } CommonLUJETS;
extern CommonLUJETS &lujets ;
extern CommonLUJETS lujets_ ;
#define cb_lujets lujets_
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//============================================================
//      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
//------------------------------------------------------------
typedef struct {
    int    mstu[200]; 
    float  paru[200];
    int    mstj[200];
    float  parj[200];
                } CommonLUDAT1;
extern CommonLUDAT1 &ludat1 ;
extern CommonLUDAT1 ludat1_ ;
#define cb_ludat1 ludat1_
//////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////
//        LUND         External FORTRAN programs                             //
//                                                                           //
//  For strings we follow Piotr Golonka's example (see lugive)               //
//extern "C" void type_of_call pyinit(char *frame, char *beam, char *target,
//                                    float &win, Long_t l_frame, Long_t l_beam,
//                                    Long_t l_target);
///////////////////////////////////////////////////////////////////////////////
extern "C" void lulist_(const long&);
//
extern "C" void luclus_(long &njet);
//
extern "C" void lugive_(char *directive, long s1);
//
///////////////////////////////////////////////////////////////////////////////



ClassImp(KoralwMaker)
//////////////////////////////////////////////////////////////////////////////
KoralwMaker::KoralwMaker()
{
// Constructor of KoralwMaker
  NevTot = 0;
  for ( int j=0; j < 10001 ; j++ ) {
    xpar[j] =0.0;
  }
}
//////////////////////////////////////////////////////////////////////////////
KoralwMaker::~KoralwMaker()
{
// Destructor of KoralwMaker
}
//////////////////////////////////////////////////////////////////////////////
void KoralwMaker::ReadData( long &ntot )
{
  // NevTot is read in c++ from second line of "koralw.input"
  // Next two parsers written in fortran scan "koralw.input" and define
  // array xpar which contain input params of KoralW

  ///////////////////////////////////
  // this is part in c++
  int j;
  char txtdum [200];
#define SkipLine infile2.getline(txtdum,100,'\n') //<-- skip to the end of line
  ifstream infile2;
  infile2.open("koralw.input"); 
  SkipLine;
  infile2 >> NevTot; SkipLine;
  ntot = NevTot;
  cout<<"NevTot = "<<NevTot <<"\n";
  //[[[[[[[[[[[
  // this parser in c++ requires rewriting!!!
  //SkipLine;
  //for (j=0; j<30; j++)
  //  {
  //    infile2>> xpar[j+1];   SkipLine;
  //    //cout <<xpar[j+1]<<"\n";
  //  };
  //]]]]]]]]]]]
  infile2.close();

  ////////////////////////////////////
  // Now we use two fortran parsers
  long ninp = 12;
  long imax = 10000;
  long ireset,italk;

  FortOpen(ninp,"KW_defaults");
  ireset =1;
  italk  =0;
  readatax2_(ninp, ireset, italk, imax, xpar);
  FortClose(ninp);

  FortOpen(ninp,"koralw.input");
  ireset =0;
  italk  =1;
  readatax2_(ninp, ireset, italk, imax, xpar);
  FortClose(ninp);
}
///////////////////////////////////////////////////////////////////////////////
void KoralwMaker::Initialize()
{
  m_EvenCounter = 0;
// Initialization of LUND parameters in comon blocks
// Note shift -1 in the indexing

// Initialisation input data before generation
  cout << "*****************************" << endl;
  cout << "**   koralw initialize     **" << endl;
  cout << "*****************************" << endl;

  // Initialization of Glibk internal fortran library
  long lentot = 50000;                      // lenght in /cglib/ labeled block
  glimit_(lentot);                          // Initialize glibk
  long   nout = 16;                         // output unit number
  goutpu_(nout);                            // Initialize glibk
  FortOpen(nout,"koralw.output");           // open fortran output file!!!

// lenght 10000+1 of xpar is not accident! see actual koralw calls???????
// Note small trick on shift +1 of index in xpar???????
  long   mode =-1;
  //koralw_(mode,&xpar[1]);          //  Initialize koralw
  koralw_(mode,xpar);          //  Initialize koralw

  LuGive("MSTU(11)=16;");               // Output unit number
  LuGive("MSTU(51)=0;");                // LUBOEI, default=0
  LuGive("MDCY(C111,1)=0;");            // inhibit pi0 decay
  LuGive("PARU(44)=10.0;");             // jet thickness pT in GeV
}
///////////////////////////////////////////////////////////////////////////////
void KoralwMaker::Finalize()
{
// Initialization of LUND parameters in comon blocks
// Note shift -1 in the indexing

  cout << "*****************************" << endl;
  cout << "**   koralw finalize       **" << endl;
  cout << "*****************************" << endl;

  long   mode = 1;
  koralw_( mode,&xpar[1]);
  
  long nout = 16;    // output unit number
  FortClose(nout);
}
///////////////////////////////////////////////////////////////////////////////
void KoralwMaker::Generate(KorEvent &E)
{
// Generating Koralw event and storing it in E
  int i,j,k,l;
  long   mode =0;
  double xpar;
  
  m_EvenCounter++;

// in generation xpar is dummy, output transfered through commons

  koralw_(mode,&xpar);


  //////////////////////////////////////////////
  //      Partons from inner MC generator     //
  //////////////////////////////////////////////
  E.m_nphot= cb_momset.nphot;
  for ( j=0; j < E.m_nphot ; j++ ) {
    //E.m_photmom[j]  
    //  = VLorenz(ssphot[j][3],ssphot[j][0],ssphot[j][1],ssphot[j][2]);
    E.m_photmom[j] = 
      VLorenz(cb_momset.sphot[3][j],
	      cb_momset.sphot[0][j],
	      cb_momset.sphot[1][j],
	      cb_momset.sphot[2][j]);
  }
  E.m_wminus = 
    VLorenz(cb_momdec.q1[3],cb_momdec.q1[0],cb_momdec.q1[1],cb_momdec.q1[2]);
  E.m_wplus  =
    VLorenz(cb_momdec.q2[3],cb_momdec.q2[0],cb_momdec.q2[1],cb_momdec.q2[2]);
  E.m_ferm1 =
    VLorenz(cb_momdec.p1[3],cb_momdec.p1[0],cb_momdec.p1[1],cb_momdec.p1[2]);
  E.m_ferm2 =
    VLorenz(cb_momdec.p2[3],cb_momdec.p2[0],cb_momdec.p2[1],cb_momdec.p2[2]);
  E.m_ferm3 =
    VLorenz(cb_momdec.p3[3],cb_momdec.p3[0],cb_momdec.p3[1],cb_momdec.p3[2]);
  E.m_ferm4 =
    VLorenz(cb_momdec.p4[3],cb_momdec.p4[0],cb_momdec.p4[1],cb_momdec.p4[2]);

  //////////////////////////////////////////////
  //      Entire event as in Lund record      //
  //////////////////////////////////////////////
  E.m_npart = cb_lujets.n;
  for ( j=0; j < E.m_npart ; j++ ) {
    //////if( part[j] != NULL ) { cout << j <<" NULLLLLLL !!! \n";}
    E.m_part[j] = PartLund(j+1,
	 cb_lujets.k[0][j],    cb_lujets.k[1][j],
	 cb_lujets.k[2][j],    cb_lujets.k[3][j],
         cb_lujets.k[4][j],
	 cb_lujets.p[0][j],    cb_lujets.p[1][j],
	 cb_lujets.p[2][j],    cb_lujets.p[3][j],
         cb_lujets.p[4][j],
	 cb_lujets.v[0][j],    cb_lujets.v[1][j],
	 cb_lujets.v[2][j],    cb_lujets.v[3][j],
         cb_lujets.v[4][j]);
    ////////E.m_part[j].print();
  }
  //lulist_(2);
  //cout << "---------------- KorEvent done ----------------" << endl;
}
///////////////////////////////////////////////////////////////////////////////
void KoralwMaker::JetDefine(KorEvent &event)
{
// define jets in final state using luclus routine from Lund
  //cout << "|||||||||||||||| Jetology ||||||||||||||||" << endl;
  static int icont=0; icont++;
  int i,j,k,l;

  //==========================================================//
  //                     JETOLOGY                             //
  //==========================================================//

  long Njet;
  luclus_(Njet);

  if ( Njet < 0 ) {
    cout << "====== WARNING ====="<< endl;
    cout << "Njet= " << Njet << endl;
  };

  event.m_njet = cb_ludat1.mstu[3-1]; // jet multiplicity

  for ( int kj=0; kj < event.m_njet  ; kj++ )
    {
      int j = event.m_npart +kj;
      event.m_jet[kj] = PartLund(j+1,
	 cb_lujets.k[0][j],    cb_lujets.k[1][j],
	 cb_lujets.k[2][j],    cb_lujets.k[3][j],
         cb_lujets.k[4][j],
	 cb_lujets.p[0][j],    cb_lujets.p[1][j],
	 cb_lujets.p[2][j],    cb_lujets.p[3][j],
         cb_lujets.p[4][j],
	 cb_lujets.v[0][j],    cb_lujets.v[1][j],
	 cb_lujets.v[2][j],    cb_lujets.v[3][j],
         cb_lujets.v[4][j]);
      event.m_part[j]=event.m_jet[kj];
    };
  //==========================================================//
  //               end   JETOLOGY                             //
  //==========================================================//
}
///////////////////////////////////////////////////////////////////////////////
void KoralwMaker::LuGive(char *directive)
{
  long s1;
  s1 = strlen(directive);
  lugive_(directive, s1);
}
///////////////////////////////////////////////////////////////////////////////
void KoralwMaker::LuList(KorEvent &event,  long Level)
{
  lulist_( Level );
}
///////////////////////////////////////////////////////////////////////////////
void KoralwMaker::FortOpen(long &lunit, char *filename)
{
  long s1;
  s1 = strlen(filename);
  fort_open_(lunit,filename, s1);
}
///////////////////////////////////////////////////////////////////////////////
void KoralwMaker::FortClose(long &lunit)
{
  fort_close_(lunit);
}
///////////////////////////////////////////////////////////////////////////////
void KoralwMaker::VarRan(double rn[],long &n )
{
  varran_(rn,n);
}
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                 END OF   CLASS   KoralwMaker                             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
