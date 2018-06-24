#ifndef TMCgenDEV_H
#define TMCgenDEV_H
///////////////////////////////////////////////////////////////////////////////
///     TMCgenDEV class
/// This is class for axiliary exercises, 
/// mainly for checking analytical integration with Monte Carlo

#include <math.h>
#include "TMath.h"
#include "TMCgen.h"
#include "TH1D.h"


class TMCgenDEV :public TMCgen
{
/// member functions
  public:
  TMCgenDEV();                // explicit default constructor for streamer
  TMCgenDEV(const char*);     // user constructor
  ~TMCgenDEV();               // explicit destructor
////////////////////////////////////////////////////////////
/// data members
  double m_Xnorm;
  double m_WT;              //! MC weight
// Model Weights
double m_WTmodel[100];  //!
double m_count;
/// Foam setup
int    m_nCells;        // No. of cells, optional, default=2000
int    m_nSampl;        // No. of MC evts/cell in exploration, default=200
int    m_kDim;          // =2 for Bremss, =3 for energy spread

///////////////////////////////////////////////////////////
/// methods obligatory
void Initialize(TRandom*, ofstream*, TH1D*);
void Finalize();
void Generate();
double Density(int, double *);   /// Method of the abstract class TFOAM_INTEGRAND
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
  ClassDef(TMCgenDEV,2); // Monte Carlo generator
};
/////////////////////////////////////////////////////////////////////////////
//                End of the class TMCgenDEV                                  //
/////////////////////////////////////////////////////////////////////////////
#endif
