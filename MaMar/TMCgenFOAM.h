#ifndef TMCgenFOAM_H
#define TMCgenFOAM_H
///////////////////////////////////////////////////////////////////////////////
///     TMCgenFOAM class
/// This is class for axiliary exercises, 
/// mainly for checking analytical integration with Monte Carlo

#include <math.h>
#include "TMath.h"
#include "TMCgen.h"
#include "TH1D.h"


//------------------------------------------------------------
//  wrappers to f77 routines in KKMC and KKsem
extern "C" void kk2f_fort_open_( const int&, const char*, int);
extern "C" void kk2f_fort_close_(const int&);
//      SUBROUTINE KK2f_Initialize(xpar)
extern "C" void kk2f_initialize_(double xpar[]);
//-----------------------
//      SUBROUTINE BornV_SetKF(KFferm)
extern "C" void bornv_setkf_( const long& ); // set SINGLE Final State
//---------------------------------
//      DOUBLE PRECISION  FUNCTION BornV_Sig0nb(CMSene)
extern "C" double bornv_sig0nb_(const double&);
//      SUBROUTINE BornV_MakeGami(CMSene,gamiCR,gami,alfi)
extern "C" void bornv_makegami_(const double&, double&, double&, double&);
//      DOUBLE PRECISION  FUNCTION BornV_Simple(KFi,KFf,svar,costhe)
extern "C" double bornv_simple_( const long&,  const long&, const double&, const double&);
//------------------------------------------------------------
//      SUBROUTINE BornV_InterpoGSW(KFf,svar,CosThe)
extern "C" double bornv_interpogsw_( const long&,  const double&, const double&);
//      DOUBLE PRECISION FUNCTION BornV_Dizet(Mode,KFi,KFf,svar,CosThe,eps1,eps2,ta,tb)
extern "C" double bornv_dizet_(const long&, const long&, const long&,
		const double&, const double&, const double&, const double&, const double&, const double& );
//------------------------------------------------------------
//      SUBROUTINE GPS_BornF(KFi,KFf,PX,CosThe,p1,m1,p2,m2,p3,m3,p4,m4,Xborn)
extern "C" void gps_bornf_(const long&, const long&, double[], const double&,
		    double[], const double&, double[], const double&, double[], const double&, double[], const double&,
		    const double&);
//------------------------------------------------------------
//      SUBROUTINE GPS_BornFoam(Mode,KFi,KFf,CMSene,CosThe,Xborn)
extern "C" void gps_bornfoam_(const long&,   const long&,   const long&,
		                      const double&, const double&, const double&);
//      DOUBLE PRECISION  FUNCTION GPS_MakeRhoFoam(XNorm)
extern "C" double gps_makerhofoam_(const double&);
//------------------------------------------------------------

//-----------------------
//      SUBROUTINE KKsem_Initialize(xpar_input)
extern "C" void kksem_initialize_(double xpar[]);
//      SUBROUTINE KKsem_SetKFfin(KFfin)
extern "C" void kksem_setkffin_( const int& );
//      SUBROUTINE KKsem_SetKeyFoB(KeyFoB)
extern "C" void kksem_setkeyfob_( const int& );
//      SUBROUTINE KKsem_VVplot_vec(key,chak,nbin,xmin,xmax,yy)
extern "C" void kksem_vvplot_vec_(const int&, char[5], const int&, const double&, const double&, double[]);
//      SUBROUTINE KKsem_SetCrange(Cmin,Cmax)
extern "C" void kksem_setcrange_(const double&, const double&);
//      SUBROUTINE KKsem_SetKeyZet(KeyZet)
extern "C" void kksem_setkeyzet_( const int& );
//      SUBROUTINE KKsem_MakeBorn(svar,Born)
extern "C" void kksem_makeborn_(const double&, double&);
//      SUBROUTINE KKsem_Ord1(KeyDist,KFi,KFf,CMSene,vv,Result)
extern "C" void kksem_ord1_(const int&, const int&, const int&, const double&, const double&, const double&);


class TMCgenFOAM :public TMCgen
{
/// member functions
  public:
  TMCgenFOAM();                // explicit default constructor for streamer
  TMCgenFOAM(const char*);     // user constructor
  ~TMCgenFOAM();               // explicit destructor
////////////////////////////////////////////////////////////
/// data members
/// obsolete part???
  double m_Xnorm;
  double m_WT;              //! MC weight
  double m_x;               //!
  double m_y;               //!
////////////////////////////////////////////////////////////
/// Physics
  double m_gnanob;     ///
  double m_pi;         ///
  double m_ceuler;     ///
  double m_alfinv;     ///
  double m_alfpi;      ///
  double m_amel;       ///
  ///
  double m_CMSene;
  double m_beam;
  double m_chini;
  double m_fin;
  double m_chfin;
  long   m_KFini;         // electron
  long   m_KFf;           // muon
  //
  int    m_KeyISR;        // Type of ISR/QED switch
  int    m_jmax;          // jmax=10000
  double m_xpar[10001];   // imput array for KKMC
  /// Foam setup
  int    m_nCells;        // No. of cells, optional, default=2000
  int    m_nSampl;        // No. of MC evts/cell in exploration, default=200
  int    m_kDim;          // =2 for Bremss, =3 for energy spread
  int    m_Mode;          // operation mode for Density
  double m_eps;           // for mapping
  double m_del;           // for mapping
  double m_vvmax;         // for mapping
////////////////////////////////////////////////////////////
// Additional Foam object for ISR+FSR without IFI
  int     m_IsFoam5;       //  ON/OFF switch
  int     m_IsFoam3;       //  ON/OFF switch
  int     m_IsFoam1;       //  ON/OFF switch
  double  m_Xsav3;        //  normalization
  double  m_Xsav1;        //  normalization
  TFOAM  *m_Foam3;        //  Additional Foam object
  TFOAM  *m_Foam1;        //  Additional Foam object
//******** MC EVENT ********
  double m_CosTheta;      //! no streamer!!!
  double m_vv;            //! ISR
  double m_uu;            //! FSR
  double m_r1;            //! IFI
  double m_r2;            //! IFI
  double m_xx;            //! total
  //
  double m_Mka;           //!
  double m_p1[4];         //!
  double m_p2[4];         //!
  double m_p3[4];         //!
  double m_p4[4];         //!
  // Model Weights
  double m_WTmodel[100];  //!
  long   m_count;
///////////////////////////////////////////////////////////
/// methods obligatory
  void Initialize(TRandom*, ofstream*, TH1D*);
  void Finalize();
  void Generate();
  double Density(int, double *);   /// Method of the abstract class TFOAM_INTEGRAND
////////////////////////////////////////////////////////////
// Foam integrand
  double Fyfs( double );
  double gamISR( double );
  double gamFSR( double );
  double gamIFI( double );
  double RhoFSR(double, double );
  double RhoISR(double, double );
  double RhoIFI(double, double );
  double RhoFSR1(double, double );
  double RhoISR1(double, double );

  double Soft_yfs(double );
  void MapPlus(  double, double, double &, double &);
  void MapMinus( double, double, double &, double &);
  void MapIFI(   double, double, double &, double &);
  double Density5(int, double *);   ///
  double Density3(int, double *);   ///
  double Density1(int, double *);   ///
  /// methods auxiliary
  void ReaData(const char*, int, double[]);
  void Vdef(double[4], const double, const double, const double, const double);
////////////////////////////////////////////////////////////////////////////
  ClassDef(TMCgenFOAM,2); // Monte Carlo generator
};
/////////////////////////////////////////////////////////////////////////////
//                End of the class TMCgenFOAM                                  //
/////////////////////////////////////////////////////////////////////////////
#endif
