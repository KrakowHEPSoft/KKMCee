#ifndef TFOAM_H
#define TFOAM_H
#include "ROOT_DEF.h"

#include<stdlib.h>
#include<math.h>

#ifdef ROOT_DEF
#include "TROOT.h"
#include "TNamed.h"
#include "TObject.h"
#include "TH1.h"
#include "TObjArray.h"
#else
#include "TFHST.h"
#endif
#include "TRandom.h"

#include "TFMAXWT.h"
#include "TFOAM_INTEGRAND.h"

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;


///////////////////////////////////////////////////////////////////////////////
#ifdef ROOT_DEF
class TFVECT : public TObject {
#else
class TFVECT{
#endif
  // constructor
  private:
    int     m_Dim;                     // Dimension
    double *m_Coords;                  // [m_Dim] Coordinates
    TFVECT  *m_Next;                   // pointer for tree construction
    TFVECT  *m_Prev;                   // pointer for tree construction
    const static int m_chat =0;        // chatability (for debug)
  public:
    TFVECT();                           // Constructor
    TFVECT(const int );                 // USER Constructor
    TFVECT(const TFVECT &);              // Copy constructor
    ~TFVECT();                          // Destructor
//////////////////////////////////////////////////////////////////////////////
//                     Overloading operators                                //
//////////////////////////////////////////////////////////////////////////////
    TFVECT& operator =(const TFVECT&);   // = operator; Substitution (const ?)
    double &operator[](int);           // [] provides POINTER to coordinate
    TFVECT& operator =(double []);      // LOAD IN entire double vector
    TFVECT& operator =(double);         // LOAD IN double numer
//////////////////////////   OTHER METHODS  //////////////////////////////////
    TFVECT& operator+=(const  TFVECT&);       // +=; add vector u+=v      (FAST)
    TFVECT& operator-=(const  TFVECT&);       // +=; add vector u+=v      (FAST)
    TFVECT& operator*=(const double&);       // *=; mult. by scalar v*=x (FAST)
    TFVECT  operator+( const  TFVECT&);       // +;  u=v+s, NEVER USE IT, SLOW!!!
    TFVECT  operator-( const  TFVECT&);       // -;  u=v-s, NEVER USE IT, SLOW!!!
    void Print(void);                  // Prints vector
    void PrintList(void);              // Prints vector and the following linked list
    const int &GetDim(void);           // Returns dimension VALUE
    double GetCoord(int i){return m_Coords[i];};   // Returns coordinate
/////////////////////////////////////////////////////////////////////////////
#ifdef ROOT_DEF
    ClassDef(TFVECT,1) //n-dimensional vector with dynamical allocation
#endif
};


///////////////////////////////////////////////////////////////////////////////
//                         Class TFCELL  used in FOAM                        //
///////////////////////////////////////////////////////////////////////////////
#ifdef ROOT_DEF
class TFCELL : public TObject {
#else
class TFCELL{
#endif
  //   static, the same for all cells!
 private:
  short int     m_kDim;             // Dimension of h-cubical subspace
  short int     m_nDim;             // Dimension of simplical subspace
  short int     m_OptMCell;         // Option Mega-Cell
  short int     m_OptCu1st;         // m_OptCubFirst=1, Numbering of dims starts with h-cubic, =0 simpl.
  int           m_nVert;            // No. of vertices in the simplex = m_nDim+1
  //   MEMBERS
 private:
  //--- MC pool (not used) ---
  int    m_nPool;                   // No of events in the MC pool
  TFVECT *m_Pool;                   // [m_nPool] Linked list of MC wt=1 event pool
  //--- linked tree organization ---
  TFCELL **m_Cell0;                 //! Common root cell, for positioning of pointers in the tree
  TFVECT **m_Vert0;                 //! Common root cell, for positioning of pointers in the tree
  int    m_Serial;                  // Serial number (index in m_Cell0)
  int    m_Status;                  // Status (active, inactive)
  int    m_Parent;                  // Pointer to parent cell
  int    m_Daught0;                 // Pointer to daughter 1
  int    m_Daught1;                 // Pointer to daughter 2
  //--- M.C. sampling and choice of the best edge ---
 private:
  double m_Xdiv;                    // Factor for division
  int    m_Best;                    // Best Edge for division
  //--- Integrals of all kinds ---
  double m_Volume;                  // Cartesian Volume of cell
  double m_Integral;                // Integral over cell (estimate from exploration)
  double m_Drive;                   // Driver  integral, only for cell build-up
  double m_Primary;                 // Primary integral, only for MC generation
  //--- Geometry of the cell ---
  int   *m_Verts;                   // [m_nVert] Pointer to ARRAY of vertices in SIMPLEX subspace
  TFVECT *m_Posi;                   // Pointer to position vector,  in H-CUBIC subspace
  TFVECT *m_Size;                   // Pointer to size vector,      in H-CUBIC subspace
  //////////////////////////////////////////////////////////////////////////////////////
  //                           METHODS                                                //
  //////////////////////////////////////////////////////////////////////////////////////
 public:
  TFCELL();                          // Default Constructor for ROOT streamers
  TFCELL(int, int, int, int);        // User Constructor
  TFCELL(TFCELL &);                  // Copy Constructor
  ~TFCELL();                         // Destructor
  void  Fill(int, int, int, int, int*, TFVECT*,TFVECT*); // Assigns values of attributes
  TFCELL&  operator=(TFCELL&);       // Substitution operator (never used)
  TFVECT*& operator[](int);          // Pointer to one of vertices
  int  operator()(int);              // integer pointer to one of vertices
  //--------------- Geometry ----------------------------------
  double  GetXdiv(void){  return m_Xdiv;};          // Pointer to Xdiv
  int     GetBest(void){  return m_Best;};          // Pointer to Best
  void    SetBest(int    Best){ m_Best =Best;};     // Set Best edge candidate
  void    SetXdiv(double Xdiv){ m_Xdiv =Xdiv;};     // Set x-division for best edge cand.
  void    GetHcub(  TFVECT&, TFVECT&);   // Get position and size vectors (h-cubical subspace)
  void    GetHSize( TFVECT& );           // Get size only of cell vector  (h-cubical subspace)
  void    GetXSimp(TFVECT &, TFVECT & , int );
  //--------------- Integrals/Volumes -------------------------
  void    CalcVolume(void);                         // Calculates volume of cell
  double  GetVolume(void){ return m_Volume;};       // Volume of cell
  double  GetIntg(void){  return m_Integral;};      // Get Integral
  double  GetDriv(void){  return m_Drive;};         // Get Drive
  double  GetPrim(void){  return m_Primary;};       // Get Primary
  void    SetIntg(double Intg){ m_Integral=Intg;};  // Set true integral
  void    SetDriv(double Driv){ m_Drive   =Driv;};  // Set driver integral
  void    SetPrim(double Prim){ m_Primary =Prim;};  // Set primary integral
  //--------------- linked tree organization ------------------
  int     GetStat(void){ return m_Status;};         // Get Status
  void    SetStat(int Stat){ m_Status=Stat;};       // Set Status
  int     GetnVert(void){ return m_nVert;};         // Get No. of vrtices in the list
  TFCELL* GetPare(void){
    if(m_Parent<0)
      return NULL;
    else
      return m_Cell0[m_Parent];
  };  // Get Pointer to pointer of parent vertex
  TFCELL* GetDau0(void){
    if(m_Daught0<0)
      return NULL;
    else
      return m_Cell0[m_Daught0];
  }; // Get Pointer to 1-st daughter vertex
  TFCELL* GetDau1(void){
    if(m_Daught1<0)
      return NULL;
    else
      return m_Cell0[m_Daught1];
  }; // Get Pointer to 2-nd daughter vertex
  void   SetDau0(int Daug){ m_Daught0 = Daug;};    // Set pointer to 1-st daughter
  void   SetDau1(int Daug){ m_Daught1 = Daug;};    // Set pointer to 2-nd daughter
  void   SetSerial(int Serial){ m_Serial=Serial;}; // Set serial number
  int    GetSerial(void){ return m_Serial;};       // Get serial number
  //--- other ---
  void   SetCell0(TFCELL **Cell0){m_Cell0 =Cell0;}; // Set root cell pointer
  void   SetVert0(TFVECT **Vert0){m_Vert0 =Vert0;}; // Set root vertex pointer
  void   Print(void);               // Prints cell content
////////////////////////////////////////////////////////////////////////////
#ifdef ROOT_DEF
 ClassDef(TFCELL,1) //Single cell of FOAM
#endif
};


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//              Auxiliary class TFMATRIX of  square matrices                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#ifdef ROOT_DEF
#include "TROOT.h"
class TFMATRIX : public TObject {
#else
class TFMATRIX{
#endif
private:
  int Dim;                                // Dimension  
  int Dim2;                               // Dimension  
  double *Matrix;                         // [dim2] Table of elements
public:
  TFMATRIX(void );
  TFMATRIX(int dim );
  TFMATRIX(TFMATRIX &From);
  ~TFMATRIX();
  TFMATRIX& operator =(TFMATRIX &From);
  double& operator()(int i, int j);
  double Determinant(void);
  void Print(void);
/////////////////////////////////////////////////////////////////////////////
#ifdef ROOT_DEF
  ClassDef(TFMATRIX,1) //Square matrices and their determinant
#endif
};


///////////////////////////////////////////////////////////////////////////////
//                                                                            //
//           Auxiliary class for organizing loop over partitions              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#ifdef ROOT_DEF
class TFPARTITION : public TObject {
#else
class TFPARTITION{
#endif
private:
  // COMPONENTS //
  int m_len;                  // lenght=dimension
  int *m_digits;              // [m_len] Partition
  // METHODS //
public:
  TFPARTITION(int len);
  ~TFPARTITION(void);
  void Reset(void);
  int Next(void);
  const int& Digit(int i);
/////////////////////////////////////////////////////////////////////////////
#ifdef ROOT_DEF
  ClassDef(TFPARTITION,1) //sequential generation of partitions
#endif
};




///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                            Class  TFOAM                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#ifdef ROOT_DEF
class TFOAM : public TObject {
#else
class TFOAM{
#endif
 private:
  // COMPONENTS //
  //-------------- Input parameters
  char  m_Name[128];       // Name of a give instance of the FOAM class
  char m_Version[10];         // Actual VERSION of the Foam like (2.xx)
  char  m_Date[40];        // Release DATE   of the Foam
  int m_nDim;              // Dimension of the simplical subspace
  int m_kDim;              // Dimension of the hyper-cubical subspace
  int m_TotDim;            // Total dimension = m_nDim+m_kDim
  int m_nCells;            // Maximum number of cells (from input)
  int m_vMax;              // Maximum number of vertices (calculated)
  int m_LastVe;            // index of the last vertex
  int m_RNmax;             // maximum r.n. requested at once
  //-------------------
  int m_OptDrive;          // type of Drive =0,1,2 for TrueVol,Sigma,WtMax
  int m_OptEdge;           // decides whether vertices are included in the sampling
  int m_OptPeek;           // type of Peek =0,1,2 for maximum, random, random (new)
  int m_OptOrd;            // Single root simplex for OptOrd=1, otherwise root h-cubic
  int m_OptMCell;          // MegaCell option, slim memory for h-cubes
  int m_Chat;              // Chat=0,1,2 chat level in output, Chat=1 normal level
  int m_OptDebug;          // m_OptDebug=1, additional histogram, SetDirectory(1)      DIPSWITCH
  int m_OptCu1st;          // m_OptCubFirst=1,0, Numbering of dims starts with h-cubic DIPSWITCH
  int m_OptVert;          // m_OptVert=0,1 (Vertices of simplex are stored=1 or not=0)
  int m_OptRej;            // OptRej=0,1, rejection in MC off, on (wt=1 events).
  //-------------------
  int  m_nBin;             // binning of edge-histogram in cell exploration
  int  m_nProj;            // number of projection edges
  int  m_nSampl;           // No of sampling when dividing (exploring) cell
  int  m_EvPerBin;         // maximum number of EFFECTIVE event per bin
  //-------------------
  int  m_N0Si;             // Start of numbering dimensions  for simplex hyperspace (see m_OptCu1st)
  int  m_N0Cu;             // Start of numbering dimensions  for hypercubic  hyperspace
  int  m_P0Si;             // Start of numbering projections for simplex hyperspace
  int  m_P0Cu;             // Start of numbering projections for hypercubic  hyperspace
  //-------------------  MULTI-BRANCHING ---------------------
  int    *m_MaskDiv;       //! [m_nProj] Dynamic Mask for  cell division, h-cubic + simplical
  int    *m_InhiDiv;       //! [m_kDim]  Flags alowing to inhibit cell division, h-cubic projection/edge
  int     m_OptPRD;        // General Option switch for PRedefined Division, for quick check
  TFVECT **m_XdivPRD;      //! Lists of division values encoded in one vector per direction
  //-------------------  GEOMETRY ----------------------------
  int m_NoAct;             // number of active cells
  int m_LastCe;            // index of the last cell
  TFCELL **m_Cells;        // [m_nCells] array of ALL cells
  TFVECT **m_VerX;         // [m_vMax] array of pointers to vertex vectors
  //------------------ M.C. generation----------------------------
  TFMAXWT *m_MCMonit;      // Monitor of MC efficiency
  double m_MaxWtRej;       // Maximum weight in rejection for getting wt=1 events
  TFCELL **m_CellsAct;     //! Array of pointers to active cells, constructed at the end of build-up
  double *m_PrimAcu;       //! Array of Cumulative Primary of active cells (for fast MC generation!)
#ifdef ROOT_DEF
  TObjArray *m_HistEdg;    // histograms of wt, one for each cell edge
  TObjArray *m_HistDbg;    // histograms of wt, extra one for debug (m_OptDebug=1)
  TH1D      *m_HistWt;     // histograms of MC wt
#else
  TFHST    **m_HistEdg;    // Array of pointers to histograms
  TFHST     *m_HistWt;     // histograms of MC wt
#endif
  double *m_MCvect;        // [m_TotDim] Generated MC vector for the outside user
  double  m_MCwt;          // MC weight
  double *m_Rvec;          // [m_RNmax] random number vector from r.n. generator m_TotDim+1 maximum elements
  //----------- Procedures 
  TFOAM_INTEGRAND *m_Rho;   // pointer to abstract class providing function to integrate
  TRandom         *m_PseRan;// generator of pseudorandom numbers
  //----------- Statistics and MC results
  long m_nCalls;           // No of function calls
  long m_nEffev;           // Total no of effective events in build-up
  double m_SumWt, m_SumWt2;// Sum of wt and wt^2
  double m_SumOve;         // Sum of overveighted events
  double m_NevGen;         // No of MC events
  double m_WtMax, m_WtMin; // Max/Min wt
  double m_Prime;          // Primary integral Iprim (Itrue=Iprim<wt>)
  double m_MCresult;       // True Integral Itrue from MC series
  double m_MCerror;        // and its error
  //----------  working space for CELL exploration -------------
  double *m_Lambda;        // [m_nDim] Internal params of the simplex l_0+l_1+...+l_{nDim-1}<1
  double *m_Alpha;         // [m_kDim] Internal params of the hyp-cubic 0<a_i<1
  //////////////////////////////////////////////////////////////////////////////////////////////
  //                                     METHODS                                              //
  //////////////////////////////////////////////////////////////////////////////////////////////
 public:
  TFOAM();                          // Constructor (used by ROOT streamer)
  TFOAM(const char*);               // User constructor
  ~TFOAM();                         // Destructor
  TFOAM(const TFOAM&);              // Copy Constructor  NOT USED
  TFOAM& operator=(const TFOAM& );  // substitution      NOT USED 
  // Initialization
  void Initialize(TRandom*, TFOAM_INTEGRAND*); // Initialization of the FOAM (grid, cells, etc).
  void InitVertices(void);          // Initializes first vertices of the basic cube
  void InitCells(void);             // Initializes first n-factorial cells inside original cube
  int  CellFill(int, TFCELL*, int*,TFVECT*,TFVECT*); // fill next cell and return its index
  void LinkCells(void);             // lift up cells after re-read from disk
  void Explore(TFCELL *Cell);       // Exploration of new cell, determine <wt>, wtMax etc.
  void Carver(int&,double&,double&);// Determine the best edge, wtmax   reduction
  void Varedu(double [], int&, double&,double&); // Determine the best edge, variace reduction
  void MakeLambda(void);            // provides random point inside simplex
  void MakeAlpha(void);             // provides random point inside hypercubic
  void Grow(void);                  // Adds new cells to FOAM until buffer is full
  long PeekMax(void);               // choose one active cell, used by Grow and also in MC generation
  TFCELL* PeekRan(void);            // choose randomly one active cell, used only by Grow
  int  Divide(TFCELL *);            // Divide iCell into two daughters; iCell retained, taged as inactive
  void MakeActiveList(void);        // Creates table of active cells used by GenerCel2
  void GenerCell(TFCELL *&);        // Chose an active cell with probability ~ Primary integral
  void GenerCel2(TFCELL *&);        // Chose an active cell with probability ~ Primary integral
  // Generation
  void   MakeEvent(void);           // Make one MC event
  void   GetMCvect(double *);       // Provides random MC vector;
  void   GetMCwt(double &);         // Provides  obtained MC weight
  double GetMCwt(void);             // Provides  obtained MC weight
  double MCgenerate(double *MCvect);// All three above function in one
  // Finalization
  void GetIntegMC(double&, double&);// Provides Inegrand and abs. error from MC run
  void GetIntNorm(double&, double&);// Provides normalization Inegrand
  void GetWtParams(const double, double&, double&, double&);// Provides MC weight parameters
  void Finalize(  double&, double&);// Prints summary of MC integration
  TFOAM_INTEGRAND *GetRho(){return m_Rho;};    // Get pointer of the distribut. (after restoring from disk)
  TRandom      *GetPseRan(){return m_PseRan;}; // Get pointer of r.n. generator (after restoring from disk)
  void SetRho(TFOAM_INTEGRAND *Rho){m_Rho=Rho;};        // Set new integrand distr.
  void SetPseRan(TRandom  *PseRan){m_PseRan=PseRan;};  // Set new r.n. generator
  void ResetRho(TFOAM_INTEGRAND *Rho);                  // ReSet integrand distr.
  void ResetPseRan(TRandom  *PseRan);                  // ReSet r.n. generator
  // Getters and Setters
  void SetnDim(int nDim){m_nDim = nDim;};         // Sets dimension of simplical subspace
  void SetkDim(int kDim){m_kDim = kDim;};         // Sets dimension of hyper-cubical subspace
  void SetnCells(long nCells){m_nCells =nCells;}; // Sets maximum number of cells
  void SetnSampl(long nSampl){m_nSampl =nSampl;}; // Sets no of MC events in cell exploration
  void SetnBin(int nBin){m_nBin = nBin;};         // Sets no of bins in histogs in cell exploration
  void SetChat(int Chat){m_Chat = Chat;};         // Sets option Chat, chat level
  void SetOptRej(int OptRej){m_OptRej =OptRej;};  // Sets option for MC rejection 
  void SetOptDrive(int OptDrive){m_OptDrive =OptDrive;}; // Sets option Drive, type of driver integrand
  void SetOptPeek( int OptPeek){ m_OptPeek  =OptPeek;};  // Sets option Peek, type of choice of cell in build-up
  void SetOptEdge( int OptEdge){ m_OptEdge  =OptEdge;};  // Sets option Edge, include or not vertices
  void SetOptOrd(  int OptOrd){  m_OptOrd   =OptOrd;};   // Sets option Ord, single root simplex or h-cub.
  void SetOptMCell(int OptMCell){m_OptMCell =OptMCell;}; // Sets MegaCell option, slim memory for h-cubes
  void SetOptVert(int OptVert){m_OptVert =OptVert;}      //Sets OptVert option
  void SetEvPerBin(int EvPerBin){m_EvPerBin =EvPerBin;}; // Sets max. no. of eff. events per bin in exploration
  void SetMaxWtRej(const double MaxWtRej){m_MaxWtRej=MaxWtRej;}; // Sets max. weight for rejection (wt=1 evts)
  void SetInhiDiv(int, int );          // Set inhibition of cell division along certain edge
  void SetXdivPRD(int, int, double[]); // Set predefined division points
  // Getters and Setters
  char  *GetVersion(void ){return m_Version;};      // Get version of the foam
  int    GetTotDim(void ){ return m_TotDim;};       // Get total dimension
  double GetPrimary(void ){return m_Prime;};        // Get value of primary integral
  void GetPrimary(double &Prime){Prime = m_Prime;}; // Get value of primary integral
  long GetnCalls(void){return m_nCalls;};           // Get total no of function calls
  long GetnEffev(void){return m_nEffev;};           // Get total no of effective wt=1 events
  // Debug
  void CheckAll(const int);   // Checks corectness of FOAM structure
  void PrintCells(void);      // Prints all cells
  void PrintVertices(void);   // Prints all Vertices
  void LaTexPlot2dim(char*);  // Generating Latex document drawing FOAM
  void RootPlot2dim(char*);   // Generating c++ code for drawing FOAM
  // Inline
 private:
  double sqr( const double x ){ return x*x;};
  double dmax(double x, double y) {  if(x>y) return x; else return y; }
  double dmin(double x, double y) {  if(x>y) return y; else return x; }
  long silnia(int n){
    long s=n; for(int i=n-1; i>1; i--) s*=i; if(n==0) s=1; return s;} //factorial
  void StopM(char* message){
    cout <<"++++ Exit in TFOAM: "<< message << endl; exit(0);}    //Error message
  //////////////////////////////////////////////////////////////////////////////////////////////
#ifdef ROOT_DEF
  ClassDef(TFOAM,1)   // General purpos self-adapting Monte Carlo event generator
#endif
};
#endif
