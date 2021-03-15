#ifndef KKeeFoam_H
#define KKeeFoam_H
///     KKeeFoam class

#include <math.h>
#include "TMath.h"
#include "TH1D.h"

#include "TMCgen.h"

#include "KKdbase.h"
#include "KKborn.h"
#include "KKdizet.h"
#include "THwtMon.h"
#include "KKceex.h"
#include "KKbvir.h"
#include "KKevent.h"
#include "KKpart.h"

extern "C" {
// SRChh/ffff_aux.f
   void fort_open_( const int&, const char*, int);
   void fort_close_(const int&);
// LHAPDF
   void   hhpdf_initialize_(double[]);
   double hhpdf_strucfunc_(const int&, const int&, const double&, const double&);
}

class KKeeFoam :public TMCgen
{
/// member functions
  public:
  KKeeFoam();                // explicit default constructor for streamer
  KKeeFoam(const char*);     // user constructor
  ~KKeeFoam();               // explicit destructor
////////////////////////////////////////////////////////////
  public:
  KKdbase  *DB;                // Database
  KKdizet  *m_DZ;              // Dizet interface
  KKborn   *m_BornDist;        // Born differential distribution

  KKevent  *m_Event;           // MC event ISR+FSR in KKMC format
  KKbvir   *m_BVR;             // Library of virtual corrections
  KKceex   *m_GPS;             // CEEX matrix element

  TFOAM   *m_Foam9;        //  Additional Foam object
  double   m_Xsav9;        //  normalization

// Dimensionality
  static const int maxPar  = 10001;    // max. num. KKMC parameters +1
  double   m_ypar[maxPar];     // xpar input parameters of KKMC, c++ indexing
  double   m_xpar[maxPar];     // xpar input parameters of KKMC, f77 indexing
//
// data members
  int      m_count7;            // debug
  int      m_count9;            // debug
  int      m_out;              // f77 output unit no

  double   m_alfpi;            // alpha/pi
  double   m_ceuler;           // Euler const.
  //
  int      m_QuarkList[6];     // list of quarks KF indices to generate
  int      m_nQuarks;          // number of quark flavors selected to generate
  //
  int      m_LeptonList[6];    // list of final leptons KF indices to generate
  int      m_nLeptons;         // number of lepton flavors selected to generate

  /// Foam setup
  int      m_nCells;        // No. of cells, optional, default=2000
  int      m_nSampl;        // No. of MC evts/cell in exploration, default=200
  int      m_kDim;          // =2 for Bremss, =3 for energy spread
  int      m_FoamMode;      // operation mode for FOAM Density
  double   m_eps;           // for mapping
  double   m_del;           // for mapping
  double   m_Mffmin;        // minimum final V-bozon mass
  double   m_vvmax;         // for mapping
  double   m_Xnorm;         // Foam normalization
  int      m_nCallsFoam0;   // No of events in Foam initialization
  // additional Foam object
  double   m_Xnorm9;        // Foam normalization
  int      m_nCallsFoam9;   // No of events in Foam initialization
  TH1D    *h_TMCgen_NORMA9; //! Normalization histogram, no streamer!!!
//***********************************************
  double   m_CMSene;        //! no streamer!!!
  double   m_XXXmin;        //! no streamer!!!
  double   m_XXXmax;        //! no streamer!!!
//******** MC EVENT ********
  double   m_XXXene;        //! no streamer!!!
  int      m_KFini;         //! no streamer!!!
  int      m_KFfin;         //! no streamer!!!
  double   m_chini;         //! current init. fermion charge
  double   m_Mbeam;         //! current init. fermion mass
  double   m_chfin;         //! current fin.  fermion charge
  double   m_Mfin;          //! current fin.  fermion mass
  double   m_CosTheta;      //! no streamer!!
  int      m_AntiQ;         //! no streamer!!
  double   m_y1;            //! 1-x1, PDF1
  double   m_y2;            //! 1-x2, PDF2
  double   m_vv;            //! ISR
  double   m_uu;            //! FSR
  double   m_r1;            //! IFI
  double   m_r2;            //! IFI
  double   m_xx;            //! total = (1-m_y1)*(1-m_y2)*(1-m_vv)*(1-m_uu)
//
  double   m_WTfoam;        //! MC weight
  double   m_wt0;           //! MC weight  EEX0, Born
  double   m_WTset[1001];   //! MC correcting weights

  ///////////////////////////////////////////////////////////
  /// methods
  inline double sqr( double x ){ return x*x;};
  inline double qub( double x ){ return x*x*x;};
  void   Initialize(TRandom*, ofstream*, TH1D*);
  void   Finalize();
  void   Generate();
  double Density(int, double *);            // Method of the abstract class TFOAM_INTEGRAND
  double Density4(int nDim, double *Xarg);  // 7-dimensional FOAM integrand
  double Density6(int nDim, double *Xarg);  // 9-dimensional FOAM integrand
  void   MapPlus( double r, double gam, double &v, double &dJac);
  void   MapMinus( double r, double gam, double &v, double &dJac);
  void   MapIFI( double r, double gam, double &v, double &R);
  double RhoISR(int KeyISR, double svar, double vv, double eps);
  double RhoFSR(int KeyFSR, double svar, double uu, double eps);
  double RhoIFI(double costhe, double uu, double eps);
  double Fyfs(double gam);
  double gamISR( double svar);
  double gamFSR( double svar);
  double gamIFI( double costhe);
  void   SetEvent( double svarZ, double CosTheta);

  double Sig0nb(double &CMSene);
  void   ReaData(const char *DiskFile, int imax, double xpar[]);
  //
////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
    ClassDef(KKeeFoam,2); // Monte Carlo generator
  };
  /////////////////////////////////////////////////////////////////////////////
  //                End of the class KKeeFoam                                  //
  /////////////////////////////////////////////////////////////////////////////
 #endif
