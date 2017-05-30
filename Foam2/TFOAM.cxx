#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdio.h>

#include"TFOAM.h"


#ifdef ROOT_DEF
ClassImp(TFVECT)
ClassImp(TFCELL)
ClassImp(TFMATRIX)
ClassImp(TFPARTITION)
ClassImp(TFOAM)
#endif

  //FFFFFF  BoX-FORMATs for nice and flexible outputs
#define BXOPE cout<<\
"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"<<endl<<\
"F                                                                              F"<<endl
#define BXTXT(text) cout<<\
"F                   "<<setw(40)<<         text           <<"                   F"<<endl
#define BX1I(name,numb,text) cout<<\
"F "<<setw(10)<<name<<" = "<<setw(10)<<numb<<" = "          <<setw(50)<<text<<" F"<<endl
#define BX1F(name,numb,text)     cout<<"F "<<setw(10)<<name<<\
          " = "<<setw(15)<<setprecision(8)<<numb<<"   =    "<<setw(40)<<text<<" F"<<endl
#define BX2F(name,numb,err,text) cout<<"F "<<setw(10)<<name<<\
" = "<<setw(15)<<setprecision(8)<<numb<<" +- "<<setw(15)<<setprecision(8)<<err<<\
                                                      "  = "<<setw(25)<<text<<" F"<<endl
#define BXCLO cout<<\
"F                                                                              F"<<endl<<\
"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"<<endl
  //FFFFFF  BoX-FORMATs ends here

#define MIN -1.0e150
#define MAX  1.0e150
#define SW2 setprecision(7) << setw(12)


  ///////////////////////////////////////////////////////////////////////////////////////
  //                                                                                   //
  //                FOAM Version 2.06                                                  //
  //                                                                                   //
  //             November 2000 (2.01)                                                  //
  //             June     2001 (2.03)                                                  //
  //             December 2001 (2.04)                                                  //
  //             January  2002 (2.05)                                                  //
  //             October  2003 (2.06)                                                  //
  //                                                                                   //
  //  Author:                                                                          //
  //    S. Jadach                                                                      //
  //    Updgrade to version 2.06 by P. Sawicki                                         //
  //    Institute of Nucl. Physics, Cracow, Poland                                     //
  //    S. Jadach@ifj.edu.pl, S.Jadach@cern.ch                                         //
  //    pawel.sawicki@ifj.edu.pl                                                       //
  //                                                                                   //
  //  Multi-dimensional general purpose Monte Carlo event generator (integrator)       //
  //  with self-adapting simplical and/or hyper-cubical "foam of cells".               //
  //                                                                                   //
  //  The essential part of version 2.x was developed during the 2000 visit in         //
  //  IFH DESY Zeuthen and the 2001 visit in CERN TH Geneva.                           //
  //                                                                                   //
  //  Short history:                                                                   //
  //      - First version FOAM 1.0 by S. Jadach, May 1999 in Fortran 77                //
  //      - Older simplical algorithm described in  physics/9910004                    //
  //      - Early C++ version by M. Ciesla & M. Slusarczyk, May 2000 (translated 1.x)  //
  //      - The essential part of version 2.x was developed during the 2000 visit in   //
  //      - IFH DESY Zeuthen and the 2001 visit in CERN TH Geneva:                     //
  //      - new algorithm of the cell division (discussion with A. Para acknowledged)  //
  //      - In June 2001 many improvements using Insure++ and CodeWizard of PARASOFT   //
  //      - December 2001, correct variance reduction and predefined division points   //
  //      - January  2002, persistency using ROOT is realized and tested               //
  //      - April 2003, memory saving solution for simplicial grid, new option OptVert //
  //      - May   2003, new interface to random number generators                      //
  //      - May   2011, back to class naming of  CPC (S.J.)                             //
  //                                                                                   //
  ///////////////////////////////////////////////////////////////////////////////////////



  ///////////////////////////////////////////////////////////////////////////////
  //                                                                           //
  //      Auxiliary class TFVECT of n-dimensional vector (dynamic allocation)  //
  //      used in FOAM                                                         //
  //      (overloaded operators: +,-,*,/,[])                                   //
  //                                                                           //
  ///////////////////////////////////////////////////////////////////////////////
TFVECT::TFVECT(){
// Default constructor creates "empty vector of zero dimension"
  m_Dim    =0;
  m_Coords =NULL;
  m_Next   =NULL;
  m_Prev   =NULL;
  //exit(0);
}
///////////////////////////////////////////////////////////////////////////////
TFVECT::TFVECT(const int n){
// USER constructor creates n-densional vector
// and allocated dynamicaly array of components
  int i;
  m_Next=NULL;
  m_Prev=NULL;
  m_Dim=n;
  m_Coords = NULL;
  if (n>0){
    m_Coords = new double[m_Dim];
    if(m_chat) {
      if(m_Coords == NULL)
	{ cout<<"TFVECT constructor failed to allocate"<<endl; exit(0);}
    }
    for (i=0; i<n; i++) *(m_Coords+i)=0.0;
  }
  if(m_chat) cout<<" USER CONSTRUCTOR TFVECT(const int) "<<endl;
}
///////////////////////////////////////////////////////////////////////////////
TFVECT::TFVECT(const TFVECT &Vect){
// "copy constructor"
  m_Next=NULL;
  m_Prev=NULL;
  m_Dim=Vect.m_Dim;
  m_Coords = NULL;
  if(m_Dim>0)  m_Coords = new double[m_Dim];
  if(m_chat) {
    if(m_Coords == NULL){ cout<<"TFVECT constructor failed to allocate m_Coords"<<endl; exit(0); }
  }
  for(int i=0; i<m_Dim; i++)
    m_Coords[i] = Vect.m_Coords[i];
  //if(m_chat) cout<<"+++++ NEVER USE Copy constructor !!!!! "<<endl;
  cout<<"+++++ NEVER USE Copy constructor !!!!! "<<endl;
  //exit(0);
}
///////////////////////////////////////////////////////////////////////////////
TFVECT::~TFVECT(){
//  Explicit destructor
  if(m_chat) cout<<" DESTRUCTOR TFVECT~ "<<endl;
  delete [] m_Coords; //  free(m_Coords)
  m_Coords=NULL;
}
//////////////////////////////////////////////////////////////////////////////
//                     Overloading operators                                //
//////////////////////////////////////////////////////////////////////////////
TFVECT& TFVECT::operator =(const TFVECT& Vect){
// operator = ;  substitution
  int i;
  if (&Vect == this) return *this;
  if( m_Dim != Vect.m_Dim )
    cout<<"++++>> TFVECT::operator=    Dims. are different: "<<m_Dim<<" and "<<Vect.m_Dim<<endl;
  if( m_Dim != Vect.m_Dim ){  // cleanup
    delete [] m_Coords;
    m_Coords = new double[m_Dim];
  }
  m_Dim=Vect.m_Dim;
  for(i=0; i<m_Dim; i++)
    m_Coords[i] = Vect.m_Coords[i];
  m_Next=Vect.m_Next;
  m_Prev=Vect.m_Prev;
  if(m_chat)  cout<<" SUBSITUTE operator = "<<endl;
  return *this;
}
//============================================================================
double &TFVECT::operator[](int n){
// [] is for acces to elements as in ordinary matrix like a[j]=b[j]
// (Perhaps against some strict rules but very practical.)
// Range protection is built in, consequently for substitution
// one should use rather use a=b than explicit loop!
  if ((n<0) || (n>=m_Dim)){
    cout<<"++++>> TFVECT::operator[], out of range"<<endl;
    exit(1);
  }
  return m_Coords[n];
}
//============================================================================
TFVECT& TFVECT::operator*=(const double &x){
// operator *=; multiply by scalar; c*=x, tested
  for(int i=0;i<m_Dim;i++) 
    m_Coords[i] = m_Coords[i]*x;
  return *this;
}
//============================================================================
TFVECT& TFVECT::operator+=(const TFVECT& Shift){
// operator +=; add vector; c*=x, tested
  if( m_Dim != Shift.m_Dim){
    cout<<"++++>> TFVECT::operator+, different dimensions= "<<m_Dim<<" "<<Shift.m_Dim<<endl;
  }
  for(int i=0;i<m_Dim;i++) 
    m_Coords[i] = m_Coords[i]+Shift.m_Coords[i];
  return *this;
}
//============================================================================
TFVECT& TFVECT::operator-=(const TFVECT& Shift){
// operator -=; add vector; c*=x, tested
  if( m_Dim != Shift.m_Dim){
    cout<<"++++>> TFVECT::operator-, different dimensions= "<<m_Dim<<" "<<Shift.m_Dim<<endl;
  }
  for(int i=0;i<m_Dim;i++) 
    m_Coords[i] = m_Coords[i]-Shift.m_Coords[i];
  return *this;
}
//============================================================================
TFVECT TFVECT::operator+(const TFVECT &p2){
// operator + ; sum of 2 vectors; c=a+b, a=a+b,       NEVER USE IT, SLOW!!!
  TFVECT Temp(m_Dim);
  Temp  = (*this);
  Temp += p2;
  return Temp;
}
//============================================================================
TFVECT TFVECT::operator-(const TFVECT &p2){
// operator - ; difference of 2 vectors; c=a-b, a=a-b, NEVER USE IT, SLOW!!!
  TFVECT Temp(m_Dim);
  Temp  = (*this);
  Temp -= p2;
  return Temp;
}
//=============================================================================
TFVECT& TFVECT::operator =(double Vect[]){
// Load in double vector, sometimes can be usefull
  int i;
  for(i=0; i<m_Dim; i++)
    m_Coords[i] = Vect[i];
  return *this;
}
//=============================================================================
TFVECT& TFVECT::operator =(double x){
// Load in double number, sometimes can be usefull
  if(m_Coords != NULL){
    for(int i=0; i<m_Dim; i++)
      m_Coords[i] = x;
  }
  return *this;
}
//////////////////////////////////////////////////////////////////////////////
//                          OTHER METHODS                                   //
//////////////////////////////////////////////////////////////////////////////
void TFVECT::Print(void){
// Printout of all components on "cout"
  int i;
  cout << "(";
  for(i=0; i<m_Dim-1; i++) cout  << SW2 << *(m_Coords+i) << ",";
  cout  << SW2 << *(m_Coords+m_Dim-1);
  cout << ")";
}
///////////////////////////////////////////////////////////////////////////////
void TFVECT::PrintList(void){
// Printout of all member vectors in the list starting from "this"
  long i=0;
  if(this == NULL) return;
  TFVECT *current=this;
  while(current != NULL){
    cout<<"vec["<<i<<"]=";
    current->Print();
    cout<<endl;
    current = current->m_Next;
    i++;
  }
}
///////////////////////////////////////////////////////////////////////////////
const int &TFVECT::GetDim(void){
// Getter for vector dimension
  return m_Dim;
}
///////////////////////////////////////////////////////////////////////////////
//                End of Class TFVECT                                        //
///////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//        Class TFCELL  used in FOAM
// Implements single cell of foam in k dimensions hyper-cubical subspace
// and n dimensions in simplical subspace.
// 
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
TFCELL::TFCELL(){
// explicit default constructor
// cout<< "TFCELL constructor 1"<<endl;
  m_nVert   = 0;
  m_nPool   = 0;
  m_Pool    = NULL;
  m_Parent  =-1;
  m_Daught0 =-1;
  m_Daught1 =-1;
  m_Verts   = NULL;
  m_Posi    = NULL;
  m_Size    = NULL;
}
TFCELL::TFCELL(int nDim, int kDim, int OptMCell, int OptCu1st){
  // Constructor of a single Cell, verices assigned later on
  // cout<< "TFCELL constructor 2"<<endl;
  if (  nDim+kDim >0){
    //---------=========----------
    m_nDim     = nDim;
    m_kDim     = kDim;
    m_OptMCell = OptMCell;   // MegaCell option, slim memory for h-cubes
    m_OptCu1st = OptCu1st;   // =1 Numbering of dimens starts with h-cubic, DIPSWITCH
    //
    m_nPool    = 0;
    m_Pool     = NULL;
    m_Status   = 1;
    m_Parent   =-1;
    m_Daught0  =-1;
    m_Daught1  =-1;
    m_Xdiv     = 0.0;
    m_Best     = 0;
    m_Volume   = 0.0;
    m_Integral = 0.0;
    m_Drive    = 0.0;
    m_Primary  = 0.0;
    // m_Verts is allocated  and filled in latter on in Fill !!!!!
    m_nVert = 0;
    m_Verts = NULL;
    //--- hyp-cubical subspace ---
    m_Posi = NULL;
    m_Size = NULL;
    if( (m_kDim>0) && (m_OptMCell == 0)  ){
      m_Posi = new TFVECT(m_kDim);
      m_Size = new TFVECT(m_kDim);
    }
  }else{
    cout<<"TFCELL::TFCELL: total dimesion has to be >0 "<<endl;
    exit(0);
  }
}

////////////////////////////////////////////////////////////////////////////////
TFCELL::TFCELL(TFCELL &From){
//  Copy constructor (not tested!)                                            //
  int i;
  cout<<"+++++ NEVER USE Copy constructor for TFCELL"<<endl;
  m_nPool       = From.m_nPool;
  for(i=0; i<m_nPool+1; i++) m_Pool[i]=From.m_Pool[i];
  m_Status      = From.m_Status;
  m_Parent      = From.m_Parent;
  m_Daught0   = From.m_Daught0; 
  m_Daught1   = From.m_Daught1;
  m_Xdiv        = From.m_Xdiv;
  m_Best        = From.m_Best;
  m_Volume      = From.m_Volume;
  m_Integral    = From.m_Integral;
  m_Drive       = From.m_Drive;
  m_Primary     = From.m_Primary;
  m_nVert       = From.m_nVert;
  m_Verts       = NULL;
  if(m_nDim>0){
      m_Verts = new int[m_nVert];
      for(i=0; i<m_nVert; i++) m_Verts[i] = (From.m_Verts)[i]; // copy constructor of FVec
  }
  m_Posi        = From.m_Posi; // used copy constructor of FVec
  m_Size        = From.m_Size; // used copy constructor of FVec
}

///////////////////////////////////////////////////////////////////////////////
TFCELL::~TFCELL(){
// Explicit destructor
  delete [] m_Verts;
  if(m_Posi != NULL ) delete m_Posi;
  if(m_Size != NULL ) delete m_Size;
}

////////////////////////////////////////////////////////////////////////////////
TFVECT*& TFCELL::operator[](int i){
// Getter of vertex true pointer
  if ((i>=0) && (i<=m_nDim))
    return m_Vert0[ m_Verts[i] ];
  else{
    cout<<"TFCELL::operator[] out of range, m_nDim= "<<m_nDim<<",  index= "<<i<<endl;
    exit(1);
    return m_Vert0[0]; // this is to make insure happy
  }
}
////////////////////////////////////////////////////////////////////////////////
int TFCELL::operator()(int i){
// Getter of vertex integer pointer
  if ((i>=0) && (i<=m_nDim))
    return m_Verts[i];
  else{
    cout<<"TFCELL::operator() out of range, m_nDim= "<<m_nDim<<",  index= "<<i<<endl;
    exit(1);
    return -1; // this is to make insure happy
  }
}

////////////////////////////////////////////////////////////////////////////////
TFCELL& TFCELL::operator=(TFCELL &From){
//  operator =    is substitution (never used)
  int i;
  cout<<" TFCELL::operator= "<<endl;
  if (&From == this) return *this;
  m_nPool       = From.m_nPool;
  for(i=0; i<m_nPool+1; i++) m_Pool[i]=From.m_Pool[i];
  m_Status      = From.m_Status;
  m_Parent      = From.m_Parent;
  m_Daught0   = From.m_Daught0; 
  m_Daught1   = From.m_Daught1;
  m_Xdiv        = From.m_Xdiv;
  m_Best        = From.m_Best;
  m_Volume      = From.m_Volume;
  m_Integral    = From.m_Integral;
  m_Drive       = From.m_Drive;
  m_Primary     = From.m_Primary;
  for(i=0; i<m_nVert; i++) m_Verts[i] = From.m_Verts[i];
  m_Posi        = From.m_Posi;
  m_Size        = From.m_Size;
  return *this;
}


////////////////////////////////////////////////////////////////////////////////
void TFCELL::Fill(int Status, int Parent, int Daugh1, int Daugh2,
		    int *Vertices, TFVECT *Posi, TFVECT *Size){
// Fills Cell which has to exist already
  int i;
  m_Status  = Status;
  m_Parent  = Parent;
  m_Daught0 = Daugh1;
  m_Daught1 = Daugh2;
  if(Vertices!=NULL){
    if( (m_Verts != NULL) || (m_nVert != 0) ){
      cout<<"TFCELL::Fill ++++++  m_Verts already allocated m_nVert="<<m_nVert<<endl;
      exit(5);
    }
    m_nVert = m_nDim+1;
    m_Verts = new int[m_nVert];
    for(i=0; i<m_nVert; i++) m_Verts[i] = Vertices[i];  // Pure Pointers
  }
  if(m_OptMCell == 0){
    if(Posi!=NULL) (*m_Posi) = (*Posi);
    if(Size!=NULL) (*m_Size) = (*Size);
  }
}

////////////////////////////////////////////////////////////////////////////////
//              GETTERS/SETTERS
////////////////////////////////////////////////////////////////////////////////
void    TFCELL::GetHcub( TFVECT &Posi, TFVECT &Size){
// Sophisticated Getter for Hyper-Cubic dimensions and size
  if(m_kDim<1) return;
  if(m_OptMCell==0){
    Posi=*m_Posi;
    Size=*m_Size;
  }else{
    TFCELL *pCell,*dCell;
    Posi = 0.0; Size=1.0; // load all components

    dCell = this;
    while(dCell != NULL){
      pCell = dCell->GetPare();
      if( pCell== NULL) break;
      if((pCell->GetDau0()==NULL) && (pCell->GetDau1()==NULL)) break; //case nDim>0, kDim>0
      int    kBest = pCell->m_Best;
      double xDivi = pCell->m_Xdiv;
      int kShift; if(m_OptCu1st) kShift = 0; else  kShift=m_nDim*(m_nDim+1)/2 ;
      if( ( kShift <= kBest) && (kBest < kShift+m_kDim) ){    // only if division was in h-cubic subspace
	int kDiv = kBest-kShift;
	if(         dCell == pCell->GetDau0()  ){
	  Size[kDiv]=Size[kDiv]*xDivi;
	  Posi[kDiv]=Posi[kDiv]*xDivi;
	}else if(   dCell == pCell->GetDau1()  ){
	  Size[kDiv]=Size[kDiv]*(1.0-xDivi);
	  Posi[kDiv]=Posi[kDiv]*(1.0-xDivi)+xDivi;
	}else{
	  cout<<" +++++ STOP in TFCELL::GetHcub "<<endl;
	  exit(2);
	}
      }
      dCell=pCell;
    }//while
  }//OptMCell
}//GetHcub


////////////////////////////////////////////////////////////////////////////////
void    TFCELL::GetHSize( TFVECT &Size){
// Sophisticated Getter for size vector of h-cubic cell
  if(m_kDim<1) return;
  if(m_OptMCell==0){
    Size=*m_Size;
  }else{
    TFCELL *pCell,*dCell;
    Size=1.0; // load all components
    dCell = this;
    while(dCell != NULL){
      pCell = dCell->GetPare();
      if( pCell== NULL) break;
      if((pCell->GetDau0()==NULL) && (pCell->GetDau1()==NULL)) break; //case nDim>0, kDim>0
      int    kBest = pCell->m_Best;
      double xDivi = pCell->m_Xdiv;
      int kShift; if(m_OptCu1st) kShift = 0; else  kShift=m_nDim*(m_nDim+1)/2 ;
      if( ( kShift <= kBest) && (kBest < kShift+m_kDim) ){    // only if division was in h-cubic subspace
	int kDiv = kBest-kShift;
	if(        dCell == pCell->GetDau0() ){
	  Size[kDiv]=Size[kDiv]*xDivi;
	}else if(  dCell == pCell->GetDau1()  ){
	  Size[kDiv]=Size[kDiv]*(1.0-xDivi);
	}else{
	  cout<<" +++++ STOP in TFCELL::GetHSize "<<endl;
	  exit(2);
	}
      }
      dCell=pCell;
    }//while
  }//OptMCell
}//GetHSize

////////////////////////////////////////////////////////////////////////////////
void TFCELL::GetXSimp(TFVECT & X, TFVECT & lambda, int kVert){
////////////////////////////////////////////////////////////////////////////////
// Translates internal simplicial coordinates into absolute //////////////////// 
////////////////////////////////////////////////////////////////////////////////

  TFCELL *pCell,*dCell;
  TFVECT *p;
  int i,j,k,iDiv,jDiv,kVold,kBest;
  long *PtsDiv;
  double sum,xDivi;
  if(m_nDim<2){ 
    cout <<"No simplicial subspace"<<endl;  
    exit(2);
   }
  PtsDiv = new long [(m_nDim+1)*m_nDim/2];
  k=0;
  X=0.;
  for(i=0; i<m_nDim+1; i++)
    for(j=i+1; j<m_nDim+1; j++)
      PtsDiv[k++]=(m_nDim+1)*i+j;

  dCell = this;
  sum=0.;
  lambda[kVert]=0.;
  for(i=0; i< m_nDim+1; i++) sum+=lambda[i];
  //  lambda.Print();
  int kShift; 
  if(m_OptCu1st) kShift = m_kDim; else  kShift=0 ;

  while(dCell != NULL){
    pCell = dCell->GetPare();
    if( pCell== NULL) break;
    if((pCell->GetDau0()==NULL) && (pCell->GetDau1()==NULL))
     break; //case of starting cubic divided into simplices 
      kBest = pCell->m_Best;
      if(kBest>=kShift ){  //If division was in simplicial subspace
      xDivi = pCell->m_Xdiv;
      iDiv=PtsDiv[kBest-kShift]/(m_nDim+1);
      jDiv=PtsDiv[kBest-kShift] % (m_nDim+1);
      kVold=kVert;
     if(  dCell == pCell->GetDau0()  ){
       kVert=jDiv;
       if(kVert!=kVold)
	 {
          lambda[kVold]=1.-sum;
          sum=1.-lambda[kVert];
          lambda[kVert]=0.;
         }
          sum+=(xDivi-1)*lambda[iDiv];
          lambda[iDiv]*=xDivi;  
     }else if(   dCell == pCell->GetDau1()  ){
       kVert=iDiv;
       if(kVert!=kVold)
	 {
          lambda[kVold]=1.-sum;
          sum=1.-lambda[kVert];
          lambda[kVert]=0.;
         }
          sum+=-xDivi*lambda[jDiv];
          lambda[jDiv]*=(1.-xDivi); 
        }
      }//End of condition that division was in simplicial subspace
      if(dCell != pCell->GetDau0() && dCell != pCell->GetDau1())
        {
	  cout<<" +++++ STOP in TFCELL::GetXSimp "<<endl;
	  exit(2);
	}
	dCell=pCell;
     }// End while
  // cout <<endl << "Vertex = "<<kVert<<endl;
  //pCell->Print();

  for(i=0; i<m_nDim+1; i++){
     p=(*dCell)[i];
     if(i==kVert){
     for(j=0; j<m_nDim; j++)
       X[j]+=(1.-sum+lambda[kVert])*p->GetCoord(j);
     }else{
     for(j=0; j<m_nDim; j++)
       X[j]+=lambda[i]*p->GetCoord(j);
     }
  }
  delete PtsDiv;
  return;

} //End of GetXSimp


//////////////////////////////////////////////////////////////////////////////////
void TFCELL::CalcVolume(void){
/////////////////////////////////////////////////////////////////////////////////
// Another way to calculate volume
/////////////////////////////////////////////////////////////////////////////////
  TFCELL *pCell,*dCell;
  int i,k,kBest;
  double xDivi,volu;
  int kShift;
  volu=1.0;
  if(m_nDim>0){ 
  if(m_OptCu1st) kShift = m_kDim; else  kShift=0 ;
  dCell = this;
  while(dCell != NULL){
    pCell = dCell->GetPare();
    if( pCell== NULL) break;
    if((pCell->GetDau0()==NULL) && (pCell->GetDau1()==NULL))
     break; //case of starting cubic divided into simplices
      kBest = pCell->m_Best;
      if(kBest>=kShift ){  //If division was in simplicial subspace
         xDivi = pCell->m_Xdiv;
         if(  dCell == pCell->GetDau0()  ){
	   volu*=xDivi;
         }else if(   dCell == pCell->GetDau1()  ){
	   volu*=(1.-xDivi);
         }
        }//End of condition that division was in simplicial subspace
        if(dCell != pCell->GetDau0() && dCell != pCell->GetDau1())
        {
	  cout<<" +++++ STOP in TFCELL::CalcVolume "<<endl;
	  exit(2);
	}
   dCell=pCell;
  } //end while
  for(i=2; i<=m_nDim; i++) volu /= i;  // divide by factorial
  }                     //end of simplicial subspace
  if(m_kDim>0){         // h-cubical subspace
    if(m_OptMCell == 0){
      for(k=0; k<m_kDim; k++) volu *= (*m_Size)[k];
    }else{
      TFVECT Size(m_kDim);
      (this)->GetHSize(Size);
      for(k=0; k<m_kDim; k++) volu *= Size[k];
    }
  }
  m_Volume =volu;
}
////////////////////////////////////////////////////////////////////////////////
void TFCELL::Print(void){
// Printout for debug purpose
  int i;
  cout <<  " Status= "<<     m_Status   <<",";
  cout <<  " Volume= "<<     m_Volume   <<",";
  cout <<  " TrueInteg= " << m_Integral <<",";
  cout <<  " DriveInteg= "<< m_Drive    <<",";
  cout <<  " PrimInteg= " << m_Primary  <<",";
  cout<< endl;
  cout <<  " Xdiv= "<<m_Xdiv<<",";
  cout <<  " Best= "<<m_Best<<",";
  cout <<  " Parent=  {"<<m_Parent <<"} "; // extra DEBUG
  cout <<  " Daught0= {"<<m_Daught0<<"} "; // extra DEBUG
  cout <<  " Daught1= {"<<m_Daught1<<"} "; // extra DEBUG
  cout<< endl;
  //
  if( m_Verts == NULL ){
    cout<<"   Vertex list is NULL"<< endl;
  }else{
    for(i=0; i<m_nVert; i++){
      cout <<"   Vert["<<i<<"]= ";
      if( m_Verts==NULL ){
	cout<< " Vertices not assigned";
      }else{
	(*this)[i]->Print();
	//cout<<"  {"<<m_Verts[i]<<"} ";   // extra DEBUG
      }
      cout<<" "<< endl;
    }
  }
  //
  if(m_Posi!=NULL){ cout <<"   Posi= "; m_Posi->Print(); cout<<","<< endl;}
  if(m_Size!=NULL){ cout <<"   Size= "; m_Size->Print(); cout<<","<< endl;}
  //
  // Alternatively, Posi and Size has to be calculated for m_OptMCell=1 
  if(m_kDim>0 && m_OptMCell != 0 ){
    TFVECT Posi(m_kDim); TFVECT Size(m_kDim);
    (this)->GetHcub(Posi,Size);
    cout <<"   Posi= "; Posi.Print(); cout<<","<< endl;
    cout <<"   Size= "; Size.Print(); cout<<","<< endl;
  }
  // Also print momenta from MC pool
  //if( m_Pool != NULL ) m_Pool->PrintList();
}
///////////////////////////////////////////////////////////////////////////////////////
//                             End of  class  TFCELL                                 //
///////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//              Auxiliary class TFMATRIX of  square matrices                 //
//         Needed essentialy for calculating volume of the simplex           //
///////////////////////////////////////////////////////////////////////////////
TFMATRIX::TFMATRIX(void ){
// Default constructor                                                       //
  Dim=0;
  Dim2=0;
  Matrix=NULL;
  cout<<"TFMATRIX Default constructor called"<<endl;
}
///////////////////////////////////////////////////////////////////////////////
TFMATRIX::TFMATRIX(int dim ){
//  principal constructor !!!                                                //
  Dim = dim;
  Dim2 = dim*dim;
  Matrix = NULL;
  if (dim>0) Matrix = new double[dim*dim];
}
///////////////////////////////////////////////////////////////////////////////
TFMATRIX::TFMATRIX(TFMATRIX &From){
// "copy constructor"     (it is used!)                                      //
  int i, j;
  Dim= From.Dim;
  Matrix = NULL;
  if(Dim>0) Matrix = new double[Dim*Dim];
  for(i=0; i<Dim; i++)
    for(j=0; j<Dim; j++)
      *(Matrix + Dim*i + j) = *(From.Matrix + Dim*i + j);
}
///////////////////////////////////////////////////////////////////////////////
TFMATRIX::~TFMATRIX(){
//   destructor explicit
  if(Dim>0)
    delete [] Matrix;
  Matrix=NULL;
}
///////////////////////////////////////////////////////////////////////////////
TFMATRIX& TFMATRIX::operator =(TFMATRIX &From){
//  '='  substitution operator
  int i, j;
  if (&From == this) return *this;
  Dim=From.Dim;
  Matrix = new double[Dim*Dim];
  for(i=0; i<Dim; i++)
    for(j=0; j<Dim; j++)
      *(Matrix + Dim*i + j) = *(From.Matrix + Dim*i + j);
  return *this;
}
///////////////////////////////////////////////////////////////////////////////
double& TFMATRIX::operator()(int i, int j){
//  '+'  substitution operator
  return *(Matrix + Dim*i + j);
}
///////////////////////////////////////////////////////////////////////////////
double TFMATRIX::Determinant(void){
//  Determinan using moderately efficient method (Gauss)
  int i, j, k, l;
  double temp, det, elmax;
  TFMATRIX m = *this;
  //TFMATRIX m(Dim); m =*this;
  det=1;
  for(i=0; i<Dim; i++){
    for(j=i+1; j<Dim; j++){ // i-th raw will be subtracted from j-th such that M[j,i]=0
      l=i;
      elmax= 0.0;
      for(k=i; k<Dim;k++) // looking for a raw with maximum modulus of i-th element;
	if( fabs(m(k,i))>elmax){
	  elmax=fabs(m(k,i));
	  l=k;
	}
      if(elmax==0.0) return 0.0;
      if (l!=i){                      // swaping two raws
	for(k=0; k<Dim; k++){
	  temp = m(i,k);
	  m(i,k)=m(l,k); 
	  m(l,k)=temp;
	}
	det *= -1;
      }
      temp = m(j,i)/m(i,i); // subtract i-th from j-th such that M[j,i]=0
      for(k=0; k<Dim; k++) m(j,k) -= m(i,k) *temp;
    }
  } // matrix is now triangle type
  for(i=0; i<Dim; i++)
    det *= m(i,i);
  return det;
}
///////////////////////////////////////////////////////////////////////////////
void TFMATRIX::Print(void){
  // just printing
  int i, j;
  cout << "Matrix [" << endl;
  for(i=0; i<Dim; i++){
    cout << "(";
    for(j=0; j<Dim-1; j++)
      cout << SW2 << (*this)(i,j) << ",";
    cout << SW2 << (*this)(i,Dim-1) << ")" << endl;
  }
  cout << "] Matrix" << endl;
}
///////////////////////////////////////////////////////////////////////////////
//                End of Class TFMATRIX                                      //
///////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//           Auxiliary class for organizing loop over partitions              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
TFPARTITION::TFPARTITION(int len){
  // constructor
  int k;
  m_len = len;
  if(m_len>0){
    m_digits = new int[m_len];
    for(k=0;k<m_len;k++) m_digits[k]=0;
  }else{
    m_digits = NULL;
  }
}
////////////////////////////////////////////////////////////
TFPARTITION::~TFPARTITION(void){
  // destructor
  if(m_digits != NULL) delete [] m_digits;
}
////////////////////////////////////////////////////////////
void TFPARTITION::Reset(void){
  // reset
  for(int k=0; k<m_len; k++){
    m_digits[k]=0;
  }
}
////////////////////////////////////////////////////////////
int TFPARTITION::Next(void){
// next one
  int j,status;
  //cout<<"====>";for(j=0; j<m_len; j++) cout<<m_digits[j]<<" ";cout<<endl;
  //cout<<"[[[[";cout<<endl;
  if(m_len>0){
    status=0;
    m_digits[m_len-1]++;           // Basic increment
    for(j=m_len-1; j>=0; j--){
      if(m_digits[j]==2){          // Overflow goes to left-most digits
	m_digits[j]=0; 
	if(j>0) m_digits[j-1]++;
      }
      status +=m_digits[j];
    }
  }else{
    status=0;
  }
  //cout<<"]]]]";cout<<endl;
  return status;        // returns 0 AFTER last partition
}
////////////////////////////////////////////////////////////
const int& TFPARTITION::Digit(int i){
  // get digit
  if(m_len>0){
    return m_digits[i];
  }else{
    return m_len;
  }
}
///////////////////////////////////////////////////////////////////////////////////////
//                             End of  class  TFPARTITION                            //
///////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                           TFOAM class starts here!                                //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
TFOAM::TFOAM(){
  // Constructor for streamer
  m_nDim     = 0;
  m_kDim     = 0;
  m_TotDim   = 0;
  m_nProj    = 0;
  m_NoAct    = 0;
  m_nCells   = 0;
  m_vMax     = 0;
  m_RNmax    = 0;
  m_MaskDiv  = NULL;
  m_InhiDiv  = NULL;
  m_XdivPRD  = NULL;
  m_VerX     = NULL;
  m_Cells    = NULL;
  m_Lambda   = NULL;
  m_Alpha    = NULL;
  //
  m_CellsAct = NULL;
  m_PrimAcu  = NULL;
  m_HistEdg  = NULL;
  m_HistWt   = NULL;
#ifdef ROOT_DEF
  m_HistDbg  = NULL;
#endif
  m_MCMonit    = NULL;
  m_PseRan     = NULL;  // Random number generator
  m_Rho = NULL;  // Integrand function
  m_MCvect     = NULL;
  m_Rvec       = NULL;
  m_Rho = NULL;  // pointer to abstract class providing function to integrate
  m_PseRan     = NULL;  // generator of pseudorandom numbers
}
///////////////////////////////////////////////////////////////////////////////////////
TFOAM::TFOAM(const char* Name){
// Constructor
  //cout<< "====> TFOAM::TFOAM Constructor "<<endl;
  if(strlen(Name)  >129) {
    cout<< " ++++STOP in TFOAM::TFOAM, Name too long "<<strlen(Name)<<endl;
    exit(2);
  }
  sprintf(m_Name,"%s",Name);         // Class name
  sprintf(m_Date,"%s","  Release date:  2003.10.07    "); // Release date
  sprintf(m_Version,"%s", "2.06"); // Release version
  m_OptDebug = 1;                // =1 additional histogram, not SetDirectory(0), DIPSWITCH
  m_OptCu1st = 1;                // =1 Numbering of dimens starts with h-cubic, DIPSWITCH
  m_MaskDiv  = NULL;             // Dynamic Mask for  cell division, h-cubic + simplical
  m_InhiDiv  = NULL;             // Flag alowing to inhibit cell division in certain projection/edge
  m_XdivPRD  = NULL;             // Lists of division values encoded in one verctor per direction
  m_VerX     = NULL;
  m_Cells    = NULL;
  m_Lambda   = NULL;
  m_Alpha    = NULL;
  m_CellsAct = NULL;
  m_PrimAcu  = NULL;
  m_HistEdg  = NULL;
  m_HistWt   = NULL;
#ifdef ROOT_DEF
  m_HistDbg  = NULL;
#endif
  m_nDim     = 0;                // dimension of simplical sub-space
  m_kDim     = 0;                // dimension of hyp-cubical sub-space
  m_nCells   = 1000;             // Maximum number of Cells,    is usually re-set
  m_vMax     = 0;                // Maximum number of vertices, is re-calculated from m_nCells
  m_nSampl   = 200;              // No of sampling when dividing cell
  m_OptPRD   = 0;                // General Option switch for PRedefined Division, for quick check
  m_OptDrive = 2;                // type of Drive =1,2 for TrueVol,Sigma,WtMax
  m_OptEdge  = 0;                // decides whether vertices are included in the sampling
  m_OptPeek  = 0;                // type of Peek =0,1 for maximum, random
  m_OptOrd   = 0;                // root cell h-cub. or simplex, =0,1
  m_OptMCell = 1;                // MegaCell option, slim memory for h-cubes
  m_OptVert  = 0;                // Vertices of simplex are not stored
  m_Chat     = 1;                // Chat=0,1,2 chat level in output, Chat=1 normal level
  m_OptRej   = 0;                // OptRej=0, wted events; OptRej=0, wt-1 events
  //------------------------------------------------------
  m_nBin     = 8;                // binning of edge-histogram in cell exploration
  m_EvPerBin =25;                // maximum no. of EFFECTIVE event per bin, =0 option is inactive
  //------------------------------------------------------
  m_nCalls = 0;                  // No of function calls
  m_nEffev = 0;                  // Total no of eff. wt=1 events in build=up
  m_LastCe =-1;                  // Index of the last cell
  m_NoAct  = 0;                  // No of active cells (used in MC generation)
  m_LastVe =-1;                  // Index of the last vertex
  m_WtMin = MAX;                 // Minimal weight
  m_WtMax = MIN;                 // Maximal weight
  m_MaxWtRej =1.10;              // Maximum weight in rejection for getting wt=1 events
  m_PseRan   = NULL;             // Initialize private copy of random number generator
  m_MCMonit  = NULL;             // MC efficiency monitoring
}

///////////////////////////////////////////////////////////////////////////////
TFOAM::~TFOAM(){
//    Class destructor
//  cout<<" DESTRUCTOR entered "<<endl;
  if(m_Cells!= NULL) {
    for(int i=0; i<m_nCells; i++) delete m_Cells[i]; // TFCELL*[]
    delete [] m_Cells;
  }
  delete [] m_Rvec;    //double[]
  delete [] m_Alpha;   //double[]
  delete [] m_Lambda;  //double[]
  delete [] m_MCvect;  //double[]
  delete [] m_PrimAcu; //double[]
  delete [] m_MaskDiv; //int[]
  delete [] m_InhiDiv; //int[]
  //
  if( m_VerX!= NULL) {
    for(int i=0; i<m_vMax; i++) delete m_VerX[i]; // TFVECT*[]
    delete [] m_VerX;
  }
  //
  if( m_XdivPRD!= NULL) {
    for(int i=0; i<m_kDim; i++) delete m_XdivPRD[i]; // TFVECT*[]
    delete [] m_XdivPRD;
  }
#ifdef ROOT_DEF
  if( m_OptDebug == 0) delete m_HistEdg;   // ????
#else
  for(int i=0;i<m_nProj;i++) delete m_HistEdg[i];
  delete [] m_HistEdg;
#endif
  delete m_MCMonit;
  delete m_HistWt;
  //
}
////////////////////////////////////////////////////////////////////////////////
TFOAM::TFOAM(const TFOAM &From){
// Copy Constructor  NOT IMPLEMENTED (NEVER USED)
  cout<<"+++++ Stop in Copy Constructor  TFOAM::TFOAM  NOT IMPLEMENTED"<<endl;
  exit(1);
  m_nDim = From.m_nDim; // this is to apease insure
  
}
////////////////////////////////////////////////////////////////////////////////
TFOAM& TFOAM::operator=(const TFOAM &From){
//  operator =    is substitution NOT IMPLEMENTED (NEVER USED)
  cout<<"+++++ Stop in  TFOAM::operator= NOT IMPLEMENTED"<<endl;
  exit(1);
  m_nDim = From.m_nDim; // this is to apease insure
  return *this;         // this is to apease insure
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void TFOAM::Initialize(TRandom *PseRan, TFOAM_INTEGRAND *fun ){
//  Basic initialization, create "FOAM of cells"

  if(m_Chat>0){
    BXOPE;
    BXTXT("****************************************");
    BXTXT("******     TFOAM::Initialize      ******");
    BXTXT("****************************************");
    BXTXT(m_Name);
    BX1F("  Version",m_Version,  m_Date);
    BX1I("     nDim",m_nDim,     " Dimension of the simplical sub-space             ");
    BX1I("     kDim",m_kDim,     " Dimension of the hyper-cubical sub-space         ");
    BX1I("   nCells",m_nCells,   " Requested number of Cells (half of them active)  ");
    BX1I("   nSampl",m_nSampl,   " No of MC events in exploration of a cell cell    ");
    BX1I("     nBin",m_nBin,     " No of bins in histograms, MC exploration of cell ");
    BX1I(" EvPerBin",m_EvPerBin, " Maximum No effective_events/bin, MC exploration  ");
    BX1I(" OptDrive",m_OptDrive, " Type of Driver   =1,2 for Sigma,WtMax            ");
    BX1I("  OptEdge",m_OptEdge,  " Decides whether vertices are included in the MC  ");
    BX1I("  OptPeek",m_OptPeek,  " Type of the cell Peek =0,1 for maximum, random   ");
    BX1I("   OptOrd",m_OptOrd,   " Root cell hyp-cub. or simplex, =0,1              ");
    BX1I(" OptMCell",m_OptMCell, " MegaCell option, slim memory for hyp-cubes       ");
    BX1I(" OptVert" ,m_OptVert, " Decides whether to store vertices of simplices   ");
    BX1I(" OptDebug",m_OptDebug, " Additional debug histogram, SetDirectory(1)      ");
    BX1I(" OptCu1st",m_OptCu1st, " Numbering of dimensions  starts with h-cubic     ");
    BX1I("   OptRej",m_OptRej,   " MC rejection on/off for OptRej=0,1               ");
    BX1F(" MaxWtRej",m_MaxWtRej, " Maximum wt in rejection for wt=1 evts");
    BXCLO;
  }

  m_PseRan = PseRan; // use central random number generator
  //m_PseRan = new TPSEMAR(); // Initialize private copy of random number generator
  m_Rho = fun;  // Integrand function defined/stored

  m_TotDim=m_nDim+m_kDim;
  if(m_TotDim==0) StopM("TFOAM::Initialize: zero dimension not alowed");

  /////////////////////////////////////////////////////////////////////////
  //                   ALLOCATE SMALL LISTS                              //
  //   it is done globaly, not for each cell, to save on allocation time //
  /////////////////////////////////////////////////////////////////////////
  m_RNmax= m_TotDim+1;
  m_Rvec = new double[m_RNmax];   // Vector of random numbers
  if(m_Rvec==NULL)  StopM("TFOAM::Initialize: Cannot initialize buffer m_Rvec");

  if(m_nDim>0){
    m_Lambda = new double[m_nDim];   // sum<1 for internal parametrization of the simplex
    if(m_Lambda==NULL)  StopM("TFOAM::Initialize: Cannot initialize buffer m_Lambda");
  }
  if(m_kDim>0){
    m_Alpha = new double[m_kDim];    // sum<1 for internal parametrization of the simplex
    if(m_Alpha==NULL)  StopM("TFOAM::Initialize: Cannot initialize buffer m_Alpha");
  }
  m_MCvect = new double[m_TotDim]; // vector generated in the MC run
  if(m_MCvect==NULL)  StopM("TFOAM::Initialize: Cannot initialize buffer m_MCvect");

  //====== variables related to MC cell exploration (projections edges) =====
  if(m_OptCu1st){
    m_N0Cu=0;   m_N0Si=m_kDim;              // hcubic dimensions are first
    m_P0Cu=0;   m_P0Si=m_kDim;
  }else{
    m_N0Cu=m_nDim;                m_N0Si=0; // simplical dimensions are first
    m_P0Cu=m_nDim*(m_nDim+1)/2;   m_P0Si=0;
  }
  m_nProj = m_nDim*(m_nDim+1)/2 +m_kDim;        // number of the PROJECTION edges
  //====== List of directions inhibited for division
  if(m_InhiDiv == NULL){
    m_InhiDiv = new int[m_kDim];
    for(int i=0; i<m_kDim; i++) m_InhiDiv[i]=0;    // Hcubic space only!!!
  }
  //====== Dynamic mask used in Explore for edge determibation
  if(m_MaskDiv == NULL){
    m_MaskDiv = new int[m_nProj];
    for(int i=0; i<m_nProj; i++) m_MaskDiv[i]=1;    // Hcub+Simplical
  }
  //====== List of predefined division values in all directions (initialized as empty)
  if(m_XdivPRD == NULL){
    m_XdivPRD = new TFVECT*[m_kDim];
    for(int i=0; i<m_kDim; i++)  m_XdivPRD[i]=NULL; // Artificialy extended beyond m_kDim
  }
  //====== Initialize list of histograms
#ifdef ROOT_DEF
  char hname[100];
  char htitle[100];
  sprintf(hname,"%s_HistEdge_",m_Name);
  m_HistWt  = new TH1D(hname,"Histogram of MC weight",100,0.0, 1.5*m_MaxWtRej); // MC weight
  m_HistEdg = new TObjArray(m_nProj);           // Initialize list of histograms
  for(int i=0;i<m_nProj;i++){
    sprintf(hname,"%s_HistEdge_%1i",m_Name,i);
    sprintf(htitle,"Edge Histogram No. %1i",i);
    //cout<<"i= "<<i<<"  hname= "<<hname<<"  htitle= "<<htitle<<endl;
    (*m_HistEdg)[i] = new TH1D(hname,htitle,m_nBin,0.0, 1.0); // Initialize histogram for each edge
    ((TH1D*)(*m_HistEdg)[i])->Sumw2();
    if( m_OptDebug == 0) ((TH1D*)(*m_HistEdg)[i])->SetDirectory(0);// exclude from diskfile
  }
  //======  extra histograms for debug purposes
  if( m_OptDebug == 1) {
  m_HistDbg = new TObjArray(m_nProj);         // Initialize list of histograms
    for(int i=0;i<m_nProj;i++){
      sprintf(hname,"%s_HistDebug_%1i",m_Name,i);
      sprintf(htitle,"Debug Histogram %1i",i);
      (*m_HistDbg)[i] = new TH1D(hname,htitle,m_nBin,0.0, 1.0); // Initialize histogram for each edge
    }
  }
#else
  m_HistWt  = new TFHST(0.0, 1.5*m_MaxWtRej, 100);
  m_HistEdg = new TFHST*[m_nProj];      // Initialize list of histograms
  for(int i=0;i<m_nProj;i++) m_HistEdg[i]= new TFHST(0.0, 1.0, m_nBin); // Initialize histogram for each edge
  if(m_HistEdg==NULL) StopM("TFOAM::Initializes: Cannot allocate m_HistEdg");
#endif

  // ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| //
  //                     BUILD-UP of the FOAM                            //
  // ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| //
  //        Allocate and Initialize BIG list of vertices
  if(m_nDim>0) InitVertices();
  //        PrintVertices(); cout<<" ===== after InitVerices ====="<<endl;
  //
  //        Allocate BIG list of cells, within cMax limit
  //        Define and explore root cell(s)
  InitCells();
  //        PrintCells(); cout<<" ===== after InitCells ====="<<endl;
  Grow();
  //        PrintCells(); cout<<" ===== after Grow      ====="<<endl;

  MakeActiveList(); // Final Preperations for the M.C. generation

  // Preperations for the M.C. generation
  m_SumWt  = 0.0;               // M.C. generation sum of Wt
  m_SumWt2 = 0.0;               // M.C. generation sum of Wt**2
  m_SumOve = 0.0;               // M.C. generation sum of overweighed
  m_NevGen = 0.0;               // M.C. generation sum of 1d0
  m_WtMax  = MIN;               // M.C. generation maximum wt
  m_WtMin  = MAX;               // M.C. generation minimum wt
  m_MCresult=m_Cells[0]->GetIntg(); // M.C. Value of INTEGRAL,temporary asignment
  m_MCresult=m_Cells[0]->GetIntg(); // M.C. Value of INTEGRAL,temporary asignment
  m_MCerror =m_Cells[0]->GetIntg(); // M.C. Value of ERROR   ,temporary asignment
  m_MCMonit = new TFMAXWT(5.0,1000);  // monitoring M.C. efficiency
  //
  if(m_Chat>0){
    double Driver = m_Cells[0]->GetDriv();
    BXOPE;
    BXTXT("***  TFOAM::Initialize FINISHED!!!   ***");
    BX1I("    nCalls",m_nCalls,  "Total number of function calls         ");
    BX1F("    XPrime",m_Prime,   "Primary total integral                 ");
    BX1F("    XDiver",Driver,    "Driver  total integral                 ");
    BX1F("  MCresult",m_MCresult,"Estimate of the true MC Integral       ");
    BXCLO;
  }
  if(m_Chat==2){ PrintCells(); PrintVertices(); }
} // Initialize


///////////////////////////////////////////////////////////////////////////////
void TFOAM::InitVertices(void){
// Alocate and Define components of vertices in the simplical subspace
// OptOrd=0: 2^nDim  initial vertices of the root unit hypercube
//           like (0000),(0001),(0010),(0011),(0100)...(1111)
// OptOrd=1:  nDim+1 initial  vertices of the root unit simplex
//           like (0000),(1001),(1100),(1110),(1111)
//
  long i,j, NoVert, nDivi, iVe;
  int k;
  TFPARTITION partition(m_nDim);

  if( m_VerX != NULL ) delete [] m_VerX;

  if(m_nDim==0) StopM("TFOAM::InitVertices: m_nDim=0");
  switch(m_OptOrd){
  case 0:
    ///////////////////////////////////////////////////////////////////////
    // cMax = 1+ nDim! +2*nDivi, vMax = 2^nDim + nDivi, where  nDivi     //
    // is the no. of binary divisions after split of the unit h-cube     //
    ///////////////////////////////////////////////////////////////////////
    nDivi= (m_nCells -1 - silnia(m_nDim))/2;
    if(nDivi < 1 ) StopM("TFOAM::Initialize: too big m_nDim or too small m_nCells");
      m_vMax= (long)(pow(2, m_nDim))+nDivi;
    cout << "TFOAM::InitVertices: projected number of vertices = "<<m_vMax<<endl;
    if(m_OptVert==0)  m_vMax= (long)(pow(2, m_nDim));    
    /////////////////////////////////////////////////////////
    //          Allocate BIG list of vertices              //
    m_VerX = new TFVECT*[m_vMax];
    for(i=0;i<m_vMax;i++)
      m_VerX[i]= new  TFVECT(m_nDim);
    if(m_VerX==NULL) StopM("TFOAM::InitVertices: Cannot initialize buffer for vertices");
    /////////////////////////////////////////////////////////
    //   Define components of the first 2^nDim vertices    //
    NoVert = (long)pow(2,m_nDim);
    iVe=0;
    partition.Reset();
    do {
      for(k=0; k<m_nDim; k++) (*m_VerX[iVe])[k] = (partition.Digit(k));
      iVe++;
      //cout<<"****> "; for(k=0; k<m_nDim; k++) cout<< (partition.Digit(k)) <<" ";cout<<endl;
    } while ( partition.Next() != 0);
    if(iVe != NoVert ) StopM("TFOAM::InitVertices: Something wrong with sum ove partitions");
    break;
  case 1:
    ///////////////////////////////////////////////////////////////////////
    // cMax = 1+ 2*nDivi,  vMax = (nDim+1) + nDivi                       //
    ///////////////////////////////////////////////////////////////////////
    nDivi= (m_nCells -1)/2;
    m_vMax= (m_nDim+1)+nDivi;
    cout << "TFOAM::InitVertices: projected number of vertices = "<<m_vMax<<endl;
    if(m_OptVert==0)  m_vMax=  m_nDim+1;    
    /////////////////////////////////////////////////////////
    //          Allocate BIG list of vertices              //
    m_VerX = new TFVECT*[m_vMax];
    for(i=0;i<m_vMax;i++)
      m_VerX[i]= new  TFVECT(m_nDim);
    if(m_VerX==NULL) StopM("TFOAM::InitVertices: Cannot initialize buffer for vertices");
    /////////////////////////////////////////////////////////
    //   Define components of the first nDim+1 vertices    //
    NoVert = m_nDim+1;
    for(i=0; i<NoVert; i++){
      for(j=0; j<m_nDim; j++){
	if(j>=i)
	  (*m_VerX[i])[j] = 0;
	else
	  (*m_VerX[i])[j] = 1;
      }
      //for(j=0; j<m_nDim; j++) cout<< (*(m_VerX+i))[j]<<" ";cout<<endl;
    }
    break;
  default:
    StopM("STOP in TFOAM::InitVertices - wrong m_OptOrd = ");  
  }//switch
  m_LastVe = NoVert-1;
  //PrintVertices(); cout<<" ===== exiting InitVerices ====="<<endl;
}

///////////////////////////////////////////////////////////////////////////////
void TFOAM::InitCells(void){
// Initialize "root part" of the FOAM of cells
  int i,j,k,Mask;
  long iVe,jser,iter,npow;
  int  *Perm,*digit;
  int    *Vertices =NULL;
  TFVECT *Posi =NULL;
  TFVECT *Size =NULL;
  //
  m_LastCe =-1;                             // Index of the last cell
  if(m_Cells!= NULL) {
    for(i=0; i<m_nCells; i++) delete m_Cells[i];
    delete [] m_Cells;
  }
  //
  m_Cells = new TFCELL*[m_nCells];
  for(i=0;i<m_nCells;i++){
    m_Cells[i]= new TFCELL(m_nDim,m_kDim,m_OptMCell,m_OptCu1st); // Allocate BIG list of cells
    m_Cells[i]->SetCell0(m_Cells);
    m_Cells[i]->SetVert0(m_VerX);
    m_Cells[i]->SetSerial(i);
  }
  if(m_Cells==NULL) StopM("TFOAM::InitCells: Cannot initialize CELLS");

  // define initial values of Posi Size vectors, if h-cubic is present
  if(m_kDim>0){
     Posi = new TFVECT(m_kDim);
     Size = new TFVECT(m_kDim);
     if(Posi==NULL||Size==NULL) StopM("TFOAM::nitCells: Cannot allocate Posi or Size");
     (*Posi)=0.0;
     (*Size)=1.0;
  }

  if(m_nDim==0 && m_kDim>0){
    /////////////////////////////////////////////////////////////////////////////
    //              Single Root Hypercube                                      //
    /////////////////////////////////////////////////////////////////////////////
    // Status, Parent,     Verts,  Posi, Size
    CellFill(1,   NULL,  Vertices, Posi, Size);  //  0-th cell ACTIVE
  }else if(m_nDim>0){
    Vertices = new int[m_nDim+1];       // ALLOCATE array of pointers to vertices
    switch(m_OptOrd){
    case 0:
      // Status, Parent, Verts,  Posi, Size
      CellFill(0,   NULL,  NULL,  NULL, NULL); // 0-th cell INACTIVE
      /////////////////////////////////////////////////////////////////////////////
      //   Splits unit cube into nDim! SIMPLICES                                 //
      //   Method of constructing permutation is primitive/simple                //
      //   It is good enough for our n<10 problem                                //
      /////////////////////////////////////////////////////////////////////////////
      // LOOP OVER Permutations of (0,1,2,3,...,m_Ndim-1) starts here
      npow = (long)pow(m_nDim,m_nDim);
      Perm  = new int[m_nDim];
      digit = new int[m_nDim];
      for(i=0; i<m_nDim; i++) Perm[i]=0; // start
      for(iter=0; iter<npow; iter++){
	Mask=1;
	for(i=0; i<m_nDim; i++)
	  for(j=i+1; j<m_nDim; j++)
	    if(Perm[i]==Perm[j]) Mask=0;
	if(Mask==1){ //         NEW PERMUTATION is found, define NEW simplical CELL
	  for(iVe=0;iVe<m_nDim+1;iVe++){
	    // This digit represents single BASIC SIMPLEX, other obtained by Permuting dimensions
	    for(k=0;k<m_nDim;k++){
	      digit[k]=0;
	      if(k<iVe) digit[k]=1;
	    }
	    // Translation from "binary" digit to serial pointer of a given vertex
	    jser=0;
	    for( k=0;k<m_nDim;k++)
	      jser= 2*jser+  digit[Perm[k]];
	    Vertices[iVe] = jser;      // FILL in pointers to vertices
	  }//iVe
	  // Status,     Parent,     Verts,  Posi, Size
	  CellFill(1, m_Cells[0],  Vertices,  Posi, Size);
	} // end of permutation exploitation
	Perm[m_nDim-1]=Perm[m_nDim-1]+1; // increment
	for(k=m_nDim-1; k>0; k--)
	  if(Perm[k]==m_nDim){ Perm[k]=0; Perm[k-1]=Perm[k-1]+1;}
      } //END OF LOOP OVER Permutations
      delete [] Perm;  // insure delete_mismatch
      delete [] digit; // insure delete_mismatch
      break;
    case 1:
      /////////////////////////////////////////////////////////////////////////////
      //   ROOT 0-th Cell is just SINGLE SIMPLEX, it is ACTIVE here of course!   //
      /////////////////////////////////////////////////////////////////////////////
      for( k=0;k<m_nDim+1;k++)
	Vertices[k] = k;          // FILL in pointers to vertices
      // Status, Parent,     Verts, Posi, Size
      CellFill(1,   NULL,  Vertices, Posi, Size);  //  0-th cell ACTIVE
      break;
    default:
      StopM("STOP in TFOAM::InitCells - wrong m_OptOrd");
    }//switch(OptOrd)
  }//if(nDim>0)

  // Exploration of the root cell(s)
  int iStart=0;                           // normal case
  if( m_OptOrd==0 && m_nDim>0 ) iStart=1; // the case of nDim! root cells
  for(long iCell=iStart; iCell<=m_LastCe; iCell++){
    Explore( m_Cells[iCell] );               // Exploration of root cell(s)
  }
  //PrintCells(); cout<<" ===== After exploration of root cells  ====="<<endl;


  //========================= CLEANUP ==============================
  if(Vertices!=NULL)   delete [] Vertices;
  if(    Posi!=NULL)   delete Posi;
  if(    Size!=NULL)   delete Size;
  //PrintCells(); cout<<" ===== exiting InitCells ====="<<endl;
}//InitCells


///////////////////////////////////////////////////////////////////////////////
void TFOAM::LinkCells(void){
  // lift up certain paremeters of cells after re-read from disk
  for(int i=0;i<m_nCells;i++){
    m_Cells[i]->SetCell0(m_Cells);
    m_Cells[i]->SetVert0(m_VerX);
  }
  MakeActiveList(); // Final Preperations for the M.C. generation
}


///////////////////////////////////////////////////////////////////////////////
int TFOAM::CellFill(int Status, TFCELL *Parent, int *Vertices, TFVECT *Posi, TFVECT *Size){
// Fills in cell all content of the new active cell
  TFCELL *Cell;
  if (m_LastCe==m_nCells){
    cout << "TFOAM::CellFill: Too many cells" << endl;
    exit(0);
  }
  m_LastCe++;   // 0-th cell is the first
  if (Status==1) m_NoAct++;

  Cell = m_Cells[m_LastCe];

  int kParent = -1;
  if(Parent!=NULL) kParent =Parent->GetSerial();

  Cell->Fill(Status, kParent, -1, -1, Vertices, Posi, Size);

  Cell->SetBest( -1);         // pointer for planning division of the cell
  Cell->SetXdiv(0.5);         // factor for division
  double xInt2,xDri2;
  if(Parent!=NULL){
    xInt2  = 0.5*Parent->GetIntg();
    xDri2  = 0.5*Parent->GetDriv();
    Cell->SetIntg(xInt2);
    Cell->SetDriv(xDri2);
  }else{
    Cell->SetIntg(0.0);
    Cell->SetDriv(0.0);
  }
  return m_LastCe;
}


///////////////////////////////////////////////////////////////////////////////
void TFOAM::Explore(TFCELL *Cell){
//   Explore newly defined cell with help of special short MC sampling
//   As a result, estimetes of true and drive volume will be defined
//   Average and dispersion of the weight distribution will be found along
//   each edge and the best edge (minimum dispersion) is memorized
//   for future use.
//   Axerage x for eventual future cell division is also defined.
//   Recorded are aso minimum and maximu weight etc.
//   The volume estimate in all (inactive) parent cells is updated
//   Note that links to parents and initial volume = 1/2 parent has to be
//   already defined prior to calling this routine.
  double Factorial;
  double Wt, Dx, Dxx, Vsum, xBest, yBest;
  double IntOld, DriOld;

  long iev;
  double NevMC;
  int iv, jv, i, j, k;
  int nProj, kBest;
  double CeSum[5], Xproj;

  TFVECT  Size(m_kDim);
  TFVECT  Posi(m_kDim);
  TFVECT  VRand(m_nDim),Lambda(m_nDim+1), X(m_nDim);
  Cell->GetHcub(Posi,Size);

  TFCELL  *Parent;
  TFMATRIX Yrel(m_nDim), Xrel(m_nDim), XVert(m_nDim+1); // m_nDim=0 is handled internaly

  double *xRand = new double[m_TotDim];

  double *VolPart=NULL;
  if(m_nDim>0){
    VolPart= new double[m_nDim+1];   // partial volumes for handling simplex geometry
    if(VolPart==NULL)   StopM("TFOAM::Initialize: Cannot initialize buffer VolPart");
                         // Table of vertex position
    for(iv=0; iv<m_nDim+1; iv++)
      {
	if(m_OptVert){
            for(j=0; j< m_nDim; j++)
            XVert(iv,j)=(*(*Cell)[iv])[j];
	    }
          else
            {
             X=0.;
             Lambda=0.;
             Cell->GetXSimp(X,Lambda,iv);
             for(j=0; j< m_nDim; j++)
             XVert(iv,j)=X[j];
            }
      }
  }//endif m_nDim

  Cell->CalcVolume();
  Dx = Cell->GetVolume();   // Dx includes simplical and h-cubical parts

  Factorial=silnia(m_nDim);           // m_nDim=0 is handled internaly

  IntOld = Cell->GetIntg(); //memorize old values,
  DriOld = Cell->GetDriv(); //will be needed for correcting parent cells

  if(m_OptVert)
  {
  if(m_nDim>0){     // decode vertex vectors, relative last vertex
    for(iv=0; iv< m_nDim; iv++)
      for(j=0; j< m_nDim; j++)
	Xrel(iv,j) =  (*(*Cell)[iv])[j] - (*(*Cell)[m_nDim])[j];
  }
  }
  /////////////////////////////////////////////////////
  //    Special Short MC sampling to probe cell      //
  /////////////////////////////////////////////////////
  CeSum[0]=0;
  CeSum[1]=0;
  CeSum[2]=0;
  CeSum[3]=MAX;  //wtmin
  CeSum[4]=MIN;  //wtmax
  //
#ifdef ROOT_DEF
  if(m_nProj>0) for(i=0;i<m_nProj;i++) ((TH1D *)(*m_HistEdg)[i])->Reset(); // Reset histograms
  m_HistWt->Reset();
#else
  if(m_nProj>0) for(i=0;i<m_nProj;i++) m_HistEdg[i]->Reset();              // Reset histograms
  m_HistWt->Reset();
#endif
  //
  ///////////////////////////////////////////////////////////////////////
  // additional scan over vertices in order to improve max/min weights //
  int icont =0;
  TFPARTITION BiPart(m_kDim);
  if(m_OptEdge==1){
    for(iv=0; iv<m_nDim+1; iv++){      // loop over vertices in the simplex
      BiPart.Reset();
      do {                             // loop over vertices in hyp-cube
	for(j=0; j<m_nDim; j++)
	  xRand[m_N0Si+j]  = XVert(iv,j);
	for(k=0; k<m_kDim; k++) 
	  xRand[m_N0Cu+k] =  Posi[k] +(BiPart.Digit(k))*(Size[k]);
	Wt = m_Rho->Density(m_TotDim, xRand )*Dx; // <wt> normalised to integral over the cell
	m_nCalls++;
	icont++; if(icont>100) break; // protection against excesive scann
	if (CeSum[3]>Wt) CeSum[3]=Wt;
	if (CeSum[4]<Wt) CeSum[4]=Wt;
	//cout<<"****> "; for(k=0; k<m_kDim; k++) cout<< (BiPart.Digit(k)) <<" ";cout<<endl;
      } while ( BiPart.Next() != 0);
    }
  }

  // ||||||||||||||||||||||||||BEGIN MC LOOP|||||||||||||||||||||||||||||
  double NevEff;
  for(iev=0;iev<m_nSampl;iev++){
    MakeLambda();              // generate uniformly vector inside simplex
    MakeAlpha();               // generate uniformly vector inside hypercube

    if(!m_OptVert && m_nDim>0 ){
        for(j=0; j<m_nDim; j++) Lambda[j]=m_Lambda[j];
        Lambda[m_nDim]=0.;
        VRand=0.;
        Cell->GetXSimp(VRand,Lambda,m_nDim);
        for(j=0; j<m_nDim; j++) xRand[m_N0Si+j]=VRand[j];
    }else{
        for(j=0; j<m_nDim; j++){
        xRand[m_N0Si+j]=(*(*Cell)[m_nDim])[j];
        for(iv=0; iv<m_nDim; iv++)
	  xRand[m_N0Si+j] += m_Lambda[iv]*Xrel(iv,j);
       }
    }
    if(m_kDim>0){
      for(j=0; j<m_kDim; j++)
	xRand[m_N0Cu+j]= Posi[j] +m_Alpha[j]*(Size[j]);
    }
    
    //-------------------------------------------------------------------------
    Wt = m_Rho->Density(m_TotDim, xRand )*Dx;   // <wt> normalised to integral over the cell
    //-------------------------------------------------------------------------
    // calculate partial volumes, necessary for projecting in simplex edges
    if(m_nDim>0){
      Vsum = 0.0;
      for(jv=0; jv<m_nDim+1; jv++){ // all vertices relative to xRand, jv is omitted
	k=0;
	for(iv=0; iv<m_nDim+1; iv++){
	  if(iv!=jv){  // vertex jv will be replaced with the random vertex
	    for(j=0; j<m_nDim; j++)
	      Yrel(k,j) = XVert(iv,j) - xRand[m_N0Si+j];
	    k++;
	  }
	}
	Dxx = Yrel.Determinant();
	VolPart[jv] = fabs(Dxx/Factorial);
	Vsum += VolPart[jv];
      }
      // this x-check only makes sense for pure simplical case!
      if( (m_kDim*m_nDim ==0) && (fabs(Vsum-Dx) > 1.0e-7*(fabs(Vsum)+fabs(Dx))) ){
	cout<<"++++ WARNING: TFOAM::Explore: Vsum= "<< Vsum <<" Dx= "<< Dx <<endl;
	//StopM("TFOAM::Explore: something wrong with volume calculation ");
      }
    }
    //=============================================
    // $**$ IMPORTANT! Two Loops below determine indexing of edges (simplex and hypercube)
    nProj = m_P0Si;
    if(m_nDim>0){
      for(jv=0; jv<m_nDim+1; jv++){
	for(iv=jv+1; iv<m_nDim+1; iv++){
	  Xproj = VolPart[jv]/(VolPart[jv]+VolPart[iv]);
#ifdef ROOT_DEF
	  ((TH1D *)(*m_HistEdg)[nProj])->Fill(Xproj,Wt);
#else
	  m_HistEdg[nProj]->Fill(Xproj,Wt); // fill all histograms, search the best edge-candidate
#endif
	  nProj++;
	}
      }
    }
    nProj = m_P0Cu;
    if(m_kDim>0){
      for(k=0; k<m_kDim; k++){
	Xproj =m_Alpha[k];
#ifdef ROOT_DEF
	((TH1D *)(*m_HistEdg)[nProj])->Fill(Xproj,Wt);
#else
	m_HistEdg[nProj]->Fill(Xproj,Wt); // fill all histograms, search the best edge-candidate
#endif
	nProj++;
      }
    }
    //
    m_nCalls++;
    CeSum[0] += Wt;    // sum of weights
    CeSum[1] += Wt*Wt; // sum of weights squared
    CeSum[2]++;        // sum of 1
    if (CeSum[3]>Wt) CeSum[3]=Wt;  // minimum weight;
    if (CeSum[4]<Wt) CeSum[4]=Wt;  // maximum weight
    // test MC loop exit condition
    NevEff = CeSum[0]*CeSum[0]/CeSum[1];
    if( NevEff >= m_nBin*m_EvPerBin) break;
  }   // ||||||||||||||||||||||||||END MC LOOP|||||||||||||||||||||||||||||
  //------------------------------------------------------------------
  //---  predefine logics of searching for the best division edge ---
  for(k=0; k<m_kDim;k++){
    m_MaskDiv[m_P0Cu+k] =1;                       // default is all
    if( m_InhiDiv[k]==1) m_MaskDiv[m_P0Cu +k] =0; // inhibit some...
  }
  // Note that predefined division below overrule inhibition above
  kBest=-1;
  double rmin,rmax,rdiv;
  if(m_OptPRD){          // quick check
    for(k=0; k<m_kDim; k++){
      rmin= Posi[k];
      rmax= Posi[k] +Size[k];
      if( m_XdivPRD[k] != NULL){
	int N= (m_XdivPRD[k])->GetDim();
	for(j=0; j<N; j++){
	  rdiv=(*m_XdivPRD[k])[j];
	  // check predefined divisions is available in this cell
	  if( (rmin +1e-99 <rdiv) && (rdiv< rmax -1e-99)){
	    kBest=k;
	    xBest= (rdiv-Posi[k])/Size[k] ;
	    goto ee05;
	  }
	}
      }
    }//k
  }
  ee05:
  //------------------------------------------------------------------
  m_nEffev += (long)NevEff;
  NevMC          = CeSum[2];
  double IntTrue = CeSum[0]/(NevMC+0.000001);
  double IntDriv, IntPrim;
  switch(m_OptDrive){
  case 1:                       // VARIANCE REDUCTION
    if(kBest == -1) Varedu(CeSum,kBest,xBest,yBest); // determine the best edge,
    //IntDriv =sqrt( CeSum[1]/NevMC -IntTrue*IntTrue ); // Older ansatz, numericaly not bad
    IntDriv =sqrt(CeSum[1]/NevMC) -IntTrue; // Foam build-up, sqrt(<w**2>) -<w>
    IntPrim =sqrt(CeSum[1]/NevMC);          // MC gen. sqrt(<w**2>) =sqrt(<w>**2 +sigma**2)
    break;
  case 2:                       // WTMAX  REDUCTION
    if(kBest == -1) Carver(kBest,xBest,yBest);  // determine the best edge
    IntDriv =CeSum[4] -IntTrue; // Foam build-up, wtmax-<w>
    IntPrim =CeSum[4];          // MC generation, wtmax!
    break;
  default:
    StopM("STOP in TFOAM::Explore - wrong m_OptDrive = ");
  }//switch
  //=================================================================================
  //hist_Neff_distrib.Fill( m_LastCe/2.0+0.01, NevEff+0.01);  // 
  //hist_kBest_distrib.Fill( kBest+0.50, 1.0 ); //  debug 
  //hist_xBest_distrib.Fill( xBest+0.01, 1.0 ); //  debug 
  //=================================================================================
  Cell->SetBest(kBest);
  Cell->SetXdiv(xBest);
  Cell->SetIntg(IntTrue);
  Cell->SetDriv(IntDriv);
  Cell->SetPrim(IntPrim);
  // correct/update integrals in all parent cells to the top of the tree
  double  ParIntg, ParDriv;
  for(Parent = Cell->GetPare(); Parent!=NULL; Parent = Parent->GetPare()){
    ParIntg = Parent->GetIntg();  
    ParDriv = Parent->GetDriv();  
    Parent->SetIntg( ParIntg   +IntTrue -IntOld );
    Parent->SetDriv( ParDriv   +IntDriv -DriOld );
  }
  delete [] VolPart;
  delete [] xRand;
  //Cell->Print();
} // TFOAM::Explore

///////////////////////////////////////////////////////////////////////////////
void TFOAM::Varedu(double CeSum[5], int &kBest, double &xBest, double &yBest){
// The case of the optimization of the maximu weight.
// Determine the best edge-candidate for future cell division,
// using results of the MC exploration run within the cell stored in m_HistEdg
  double Nent   = CeSum[2];
  double swAll  = CeSum[0];
  double sswAll = CeSum[1];
  double SSw    = sqrt(sswAll)/sqrt(Nent);
  //
  double SwIn,SwOut,SSwIn,SSwOut,xLo,xUp;
  kBest =-1;
  xBest =0.5;
  yBest =1.0;
  double MaxGain=0.0;
  // Now go over all projections kProj
  for(int kProj=0; kProj<m_nProj; kProj++)
    if( m_MaskDiv[kProj]){
    // initialize search over bins
    double SigmIn =0.0; double SigmOut =0.0;
    double SSwtBest = MAX;
    double Gain =0.0;
    double xMin=0.0; double xMax=0.0; 
    // Double loop over all pairs jLo<jUp
    for(int jLo=1; jLo<=m_nBin; jLo++){
      double swIn=0;  double sswIn=0;
      for(int jUp=jLo; jUp<=m_nBin;jUp++){
#ifdef ROOT_DEF
	swIn  +=     ((TH1D *)(*m_HistEdg)[kProj])->GetBinContent(jUp);
	sswIn += sqr(((TH1D *)(*m_HistEdg)[kProj])->GetBinError(  jUp));
#else
	swIn  +=     (m_HistEdg[kProj])->GetBinContent(jUp);
	sswIn += sqr((m_HistEdg[kProj])->GetBinError(  jUp));
#endif
	xLo=(jLo-1.0)/m_nBin;
	xUp=(jUp*1.0)/m_nBin;
	SwIn  =        swIn/Nent;
	SwOut = (swAll-swIn)/Nent;
	SSwIn = sqrt(sswIn)       /sqrt(Nent*(xUp-xLo))     *(xUp-xLo);
	SSwOut= sqrt(sswAll-sswIn)/sqrt(Nent*(1.0-xUp+xLo)) *(1.0-xUp+xLo);
	if( (SSwIn+SSwOut) < SSwtBest){
	  SSwtBest = SSwIn+SSwOut;
	  Gain     = SSw-SSwtBest;
	  SigmIn   = SSwIn -SwIn;  // Debug
	  SigmOut  = SSwOut-SwOut; // Debug
	  xMin    = xLo;
	  xMax    = xUp;
	}
      }//jUp
    }//jLo
    int iLo = (int) (m_nBin*xMin);
    int iUp = (int) (m_nBin*xMax);
    //----------DEBUG printout
    //cout<<"@@@@@  xMin xMax = "<<xMin   <<" "<<xMax<<"  iLo= "<<iLo<<"  iUp= "<<iUp;
    //cout<<"  SSwtBest/SSw= "<<SSwtBest/SSw<<"  Gain/SSw= "<< Gain/SSw<<endl;
    //----------DEBUG auxilairy Plot
#ifdef ROOT_DEF
    if( m_OptDebug == 1) {
      for(int iBin=1;iBin<=m_nBin;iBin++)
	if( ((iBin-0.5)/m_nBin > xMin) && ((iBin-0.5)/m_nBin < xMax) ){
	  ((TH1D *)(*m_HistDbg)[kProj])->SetBinContent(iBin,SigmIn/(xMax-xMin));
	}else{
	  ((TH1D *)(*m_HistDbg)[kProj])->SetBinContent(iBin,SigmOut/(1-xMax+xMin));
	}
    }
#endif
    if(Gain>MaxGain){
      MaxGain=Gain;
      kBest=kProj; // <--- !!!!! The best edge
      xBest=xMin;
      yBest=xMax;
      if(iLo == 0   ) xBest=yBest; // The best division point
      if(iUp == m_nBin) yBest=xBest; // this is not realy used
    }
  }
  //----------DEBUG printout
  //cout<<"@@@@@@@>>>>> kBest= "<<kBest<<"  MaxGain/SSw= "<< MaxGain/SSw<<endl;
  if( (kBest >= m_nProj) || (kBest<0) ) StopM("TFOAM::Varedu: something wrong with kBest");
}          //TFOAM::Varedu

///////////////////////////////////////////////////////////////////////////////
void TFOAM::Carver(int &kBest, double &xBest, double &yBest){
// The case of the optimization of the maximu weight.
// Determine the best edge-candidate for future cell division,
// using results of the MC exploration run within the cell stored in m_HistEdg
  int    kProj,iBin;
  double Carve,CarvTot,CarvMax,CarvOne,BinMax,BinTot,PrimTot,PrimMax;
  int    jLow,jUp,iLow,iUp;
  double TheBin;
  int    jDivi; // TEST

  double *Bins  = new double[m_nBin];      // bins of histogram for single  PROJECTION
  if(Bins==NULL)    StopM("TFOAM::Carver: Cannot initialize buffer Bins");
  
  kBest =-1;
  xBest =0.5;
  yBest =1.0;
  CarvMax = MIN;
  PrimMax = MIN;
  for(kProj=0; kProj<m_nProj; kProj++)
    if( m_MaskDiv[kProj] ){
    //if( kProj==1 ){
    //cout<<"==================== Carver histogram: kProj ="<<kProj<<"==================="<<endl;
    //((TH1D *)(*m_HistEdg)[kProj])->Print("all");
    BinMax = MIN;
    for(iBin=0; iBin<m_nBin;iBin++){
#ifdef ROOT_DEF
      Bins[iBin]= ((TH1D *)(*m_HistEdg)[kProj])->GetBinContent(iBin+1);
#else
      Bins[iBin]= (m_HistEdg[kProj])->GetBinContent(iBin+1);      // Unload histogram
#endif
      BinMax = dmax( BinMax, Bins[iBin]);       // Maximum content/bin
    }
    if(BinMax < 0 ) {       //case of empty cell
      delete [] Bins;
      return;
    }
    CarvTot = 0.0;
    BinTot  = 0.0;
    for(iBin=0;iBin<m_nBin;iBin++){
      CarvTot = CarvTot + (BinMax-Bins[iBin]);     // Total Carve (more stable)
      BinTot  +=Bins[iBin];
    }
    PrimTot = BinMax*m_nBin;
     //cout <<"Carver:  CarvTot "<<CarvTot<< "    PrimTot "<<PrimTot<<endl;
    jLow =0;
    jUp  =m_nBin-1;
    CarvOne = MIN;
    float Ylevel = MIN;
    for(iBin=0; iBin<m_nBin;iBin++){
      TheBin = Bins[iBin];
      //-----  walk to the left and find first bin > TheBin
      iLow = iBin;
      for(int j=iBin; j>-1; j-- ) {
	if(TheBin< Bins[j]) break;
	iLow = j;
      }
      //iLow = iBin;
      //if(iLow>0)     while( (TheBin >= Bins[iLow-1])&&(iLow >0) ){iLow--;} // horror!!!
      //------ walk to the right and find first bin > TheBin
      iUp  = iBin;
      for(int j=iBin; j<m_nBin; j++){
	if(TheBin< Bins[j]) break;
	iUp = j;
      }
      //iUp  = iBin;
      //if(iUp<m_nBin-1) while( (TheBin >= Bins[iUp+1])&&( iUp<m_nBin-1 ) ){iUp++;} // horror!!!
      //
      Carve = (iUp-iLow+1)*(BinMax-TheBin);
      if( Carve > CarvOne){
	CarvOne = Carve;
	jLow = iLow;
	jUp  = iUp;
	Ylevel = TheBin;
      }
    }//iBin
    if( CarvTot > CarvMax){
      CarvMax   = CarvTot;
      PrimMax   = PrimTot;
      //cout <<"Carver:   PrimMax "<<PrimMax<<endl;
      kBest = kProj;    // Best edge
      xBest = ((double)(jLow))/m_nBin;
      yBest = ((double)(jUp+1))/m_nBin;
      if(jLow == 0 )       xBest = yBest;
      if(jUp  == m_nBin-1) yBest = xBest;
      // division ratio in units of 1/m_nBin, testing
      jDivi = jLow;
      if(jLow == 0 )     jDivi=jUp+1;
    }
    //======  extra histograms for debug purposes
    //cout<<"kProj= "<<kProj<<" jLow= "<<jLow<<" jUp= "<<jUp<<endl;
#ifdef ROOT_DEF
    if( m_OptDebug == 1) {
      for(iBin=0;    iBin<m_nBin;  iBin++)
	((TH1D *)(*m_HistDbg)[kProj])->SetBinContent(iBin+1,BinMax);
      for(iBin=jLow; iBin<jUp+1;   iBin++)
	((TH1D *)(*m_HistDbg)[kProj])->SetBinContent(iBin+1,Ylevel);
    }
#endif
  }//kProj
  if( (kBest >= m_nProj) || (kBest<0) ) StopM("TFOAM::Carver: something wrong with kBest");
  delete [] Bins;
}          //TFOAM::Carver


///////////////////////////////////////////////////////////////////////////////
void TFOAM::MakeAlpha(void){
// HYP-CUBICAL SUBSPACE
// Provides random vector Alpha  0< Alpha(i) < 1
  int k;
  if(m_kDim<1) return;

  // simply generate and load kDim uniform random numbers
  m_PseRan->RndmArray(m_kDim,m_Rvec);   // kDim random numbers needed
  for(k=0; k<m_kDim; k++) m_Alpha[k] = m_Rvec[k];
} //MakeAlpha


///////////////////////////////////////////////////////////////////////////////
void TFOAM::MakeLambda(void){
// SIMPLICAL SUBSPACE
// Provides random vector Lambda such that Sum Lamba(i) < 1, with uniform
// probability. This  vector is used to populate uniformly the interior
// of a simplex. The method is: generate point inside cube, order components
// (maping into simplex) and take differences of Lambda(i+1) - Lambda(i)
// For higher dimension random-walk method is used instead of ordering.
  int i,k;
  double x, sum;
  if(m_nDim<1) return;

  int nDimMax = 4; // maximum dimension for bubble-sort ordering method

  if (m_nDim>nDimMax){                   // faster random-walk algorithm for high dimensions
    m_PseRan->RndmArray(m_nDim+1,m_Rvec);  // nDim+1 random numbers needed
    sum = 0.0;
    for(i=0; i<=m_nDim; i++){
      sum += -log(m_Rvec[i]);
      m_Rvec[i]=sum;                     // ordering is here automatic
    }
    for(i=0; i<m_nDim; i++) m_Rvec[i] = m_Rvec[i] /m_Rvec[m_nDim]; // normalize to last one
  }else{                                 // bubble-sort ordering (maping cube->simplex)
    m_PseRan->RndmArray(m_nDim,m_Rvec);    // nDim random numbers needed
    for(i=m_nDim-1; i>=0; i--){
      for(k=1; k<=i; k++){
        if ( m_Rvec[k] < m_Rvec[k-1] ){
          x = m_Rvec[k]; m_Rvec[k] = m_Rvec[k-1]; m_Rvec[k-1] = x;} // get ordering
      }
    }
  }
  m_Lambda[0] = m_Rvec[0];
  for(k = 1 ; k < m_nDim ; k++)
    m_Lambda[k] = m_Rvec[k] - m_Rvec[k-1]; // Sum of lambda's should < 1 !!!
}   //MakeLambda



///////////////////////////////////////////////////////////////////////////////
void TFOAM::Grow(void){
//  Grow new cells by division
  long iCell;
  TFCELL* newCell;

  while ( (m_LastCe+2) < m_nCells ){  // this condition also checked inside Divide
    switch (m_OptPeek){
    case 0:
      iCell   = PeekMax();        // peek up cell with maximum driver integral
      if( (iCell<0) || (iCell>m_LastCe) ) StopM("TFOAM::Grow: wrong iCell");
      newCell = m_Cells[iCell];
      break;
    case 1:
      newCell = PeekRan();        // peek up randomly one cell, tree-wise algorithm
      break;
    default:
      StopM("TFOAM::Grow: wrong value of m_OptPeek");
    }
    if(m_LastCe !=0){
      int kEcho=10;
      if(m_LastCe>=10000) kEcho=100;
      if( (m_LastCe%kEcho)==0){
	if(m_nDim+m_kDim<10)
	  cout<<m_nDim+m_kDim<<flush;
	else
	  cout<<"."<<flush;
	if( (m_LastCe%(100*kEcho))==0)  cout<<"|"<<m_LastCe<<endl<<flush;
      }
    }
    //
    if( Divide( newCell )==0) break;  // and divide it into two
  }
  cout<<endl<<flush;
  CheckAll(0);   // set arg=1 for more info
}// Grow

///////////////////////////////////////////////////////////////////////////////
long  TFOAM::PeekMax(void){
// Initialization, Finds cell with maximal Driver integral
   long  i;
   long iCell = -1;
   double  DrivMax, Driv;

   DrivMax = -1E+150;
   for(i=0; i<=m_LastCe; i++) //without root
   {
      if( m_Cells[i]->GetStat() == 1 )
      {
         Driv =  fabs( m_Cells[i]->GetDriv());
	 //cout<<"PeekMax: Driv = "<<Driv<<endl;
         if(Driv > DrivMax)
         {
            DrivMax = Driv;
            iCell = i;
         }
      }
   }
//  cout << 'TFOAM_PeekMax: iCell=' << iCell << endl;
   if (iCell == -1)
      cout << "STOP in TFOAM::PeekMax: not found iCell=" <<  iCell << endl;
   return(iCell);
}                 // TFOAM_PeekMax


///////////////////////////////////////////////////////////////////////////////
TFCELL*  TFOAM::PeekRan(void){
// Initialization phase.
//       Tree-wise algorithm
//       Peek up randomly pointer iCell of an active cell
// We walk randomly from top of tree downwards until we find active
// cell m_CeStat=1
// At each step one of daugters cells is choosen randomly according
// to their volume estimates.
   long    iCell,iDau,nDaut;
   double  random,p1,volu1,volu2,sum;
   //
   TFCELL *newCell = m_Cells[0];
   //
   if( (m_nDim>1)&&(m_OptOrd==0)){
     // 0-th cell is special because has nDim! daughters
     double TotDriv = 0.0;
     double TotIntg = 0.0;
     nDaut = silnia(m_nDim);
     for(iCell= 1 ; iCell<=nDaut ; iCell++){ // without root
       //(m_Cells+iCell)->Print();
       TotDriv += m_Cells[iCell]->GetDriv();
       TotIntg += m_Cells[iCell]->GetIntg();
     }
     if( (fabs(TotDriv -m_Cells[0]->GetDriv() ) > 1.0e-5*TotDriv) )
       StopM("STOP in TFOAM::Peek -> something wrong with total driver integral");
     if( (fabs(TotIntg -m_Cells[0]->GetIntg() ) > 1.0e-5*TotDriv) )
       StopM("STOP in TFOAM::Peek -> something wrong with total true integral");
     random=m_PseRan->Rndm();
     iDau  = 0;
     sum   = 0.0;
     for(iCell = 1 ; iCell<=nDaut ; iCell++){     //without root
       iDau = iCell;
       sum += m_Cells[iCell]->GetDriv();
       if( random < sum/TotDriv ) break;
     }
     if (iDau<1 ) StopM("TFOAM::PeekRan: something went wrong, iDau<1 !!!!");
     newCell = m_Cells[iDau];
     if( newCell->GetStat() == 1 ) return newCell;
   }
   // now the other standard cells with 2 daughters
   while (   newCell->GetStat() != 1 ){
     volu1 = newCell->GetDau0()->GetDriv();
     volu2 = newCell->GetDau1()->GetDriv();
     p1 = volu1/(volu1+volu2);
     random=m_PseRan->Rndm();
     if( random < p1 )
       newCell = newCell->GetDau0();
     else
       newCell = newCell->GetDau1();
   }
   return newCell;
}        //Peek


///////////////////////////////////////////////////////////////////////////////
int TFOAM::Divide(TFCELL *Cell){
//  Divide cell iCell into two daughter cells.
//  The iCell is retained and taged as inactive, daughter cells are appended
//  at the end of the buffer.
//  New vertex is added to list of vertices.
//  List of active cells is updated, iCell remooved, two daughters added
//  and their properties set with help of MC sampling (TFOAM_Explore)
//  Return Code RC=-1 of buffer limit is reached,  m_LastCe=m_nBuf
  double Xdiv;
  int  *kVer1 =NULL;
  int  *kVer2 =NULL;
  TFVECT *Posi1  =NULL;
  TFVECT *Size1  =NULL;
  TFVECT *Posi2  =NULL;
  TFVECT *Size2  =NULL;
  TFVECT *NewVert, *OldVe1, *OldVe2;
  int   kBest, j, jv, Old1,Old2;

  if(m_LastCe+1 >= m_nCells) StopM("TFOAM::Divide: buffer limit is reached, m_LastCe=m_nBuf");

  if(m_nDim>0 && m_OptVert){
    kVer1 = new int[m_nDim+1];  // array of pointers to vertices for daughter1
    kVer2 = new int[m_nDim+1];  // array of pointers to vertices for daughter2
    if(kVer1==NULL||kVer2==NULL) StopM("TFOAM::Divide: Cannot initialize buffer kVer");
    for(jv=0; jv<m_nDim+1; jv++){
      kVer1[jv] = (*Cell)(jv);   //inherit vertex pointers from parent
      kVer2[jv] = (*Cell)(jv);   //inherit vertex pointers from parent
    }
  }
  if( (m_kDim>0)  && (m_OptMCell==0) ){
    Posi1 = new TFVECT(m_kDim);
    Posi2 = new TFVECT(m_kDim);
    Size1 = new TFVECT(m_kDim);
    Size2 = new TFVECT(m_kDim);
    Cell->GetHcub( *Posi1 , *Size1 );   //inherit from parent
    Cell->GetHcub( *Posi2 , *Size2 );   //inherit from parent
  }

  Cell->SetStat(0); // reset Cell as inactive
  m_NoAct--;

  Xdiv  = Cell->GetXdiv();
  kBest = Cell->GetBest();
  if( kBest<0 || kBest>=m_nProj ) StopM("TFOAM::Divide: wrong kBest");

  int nProj0 = (m_nDim+1)*m_nDim/2;
  if( ( m_P0Si <= kBest ) && ( kBest < m_P0Si+nProj0 ) ){
    if(m_OptVert){
    // decoding pair (Old1,Old2) from kBest (a bit tricky method, see $**$)
    Old1 =0; Old2 =0;
    for(j=m_P0Si; j<=kBest; j++){
      Old2++;
      if( Old2 > m_nDim ){ Old1++; Old2 =Old1+1; }
    }
    OldVe1 =(*Cell)[Old1];
    OldVe2 =(*Cell)[Old2];
    if (m_LastVe+1>m_vMax) StopM("STOP in TFOAM::Divide: too short list of vertices ");
    m_LastVe++;                  // add new vertex to the list
    NewVert  =m_VerX[m_LastVe];  // pointer of new vertex
    for(j=0; j<m_nDim; j++)      // define components of new vertex using Xdiv and (Old1,Old2)
      (*NewVert)[j] = Xdiv* (*OldVe1)[j] + (1.0-Xdiv)* (*OldVe2)[j];
    kVer1[Old1] = m_LastVe;     // old pointer is overwriten
    kVer2[Old2] = m_LastVe;     // old pointer is overwriten
    }
  }else{
    if( m_OptMCell==0 ) {
      ////////////////////////////////////////////////////////////////
      //   The best division is in hyp-cubical subspace             //
      ////////////////////////////////////////////////////////////////
      int kDiv = kBest-m_P0Cu;
      if(kDiv<0||kDiv>=m_kDim) StopM("Divide: wrong kDiv");  // kDiv is form zero to m_kDim-1 !!!
      (*Posi1)[kDiv]=                   (*Posi1)[kDiv];      // (1) Position unchanged
      (*Size1)[kDiv]=              Xdiv*(*Size1)[kDiv];      // (1) Size reduced by Xdiv
      (*Posi2)[kDiv]=   (*Posi1)[kDiv] +(*Size1)[kDiv];      // (2) Position shifted
      (*Size2)[kDiv]=        (1.0-Xdiv)*(*Size2)[kDiv];      // (2) Size reduced by (1-Xdiv)
    }
  }
  //////////////////////////////////////////////////////////////////
  //           define two daughter cells (active)                 //
  //////////////////////////////////////////////////////////////////
  //          Status, Parent, Verts,  Posi, Size
  int d1 = CellFill(1,   Cell, kVer1, Posi1, Size1);
  int d2 = CellFill(1,   Cell, kVer2, Posi2, Size2);
  Cell->SetDau0(d1);
  Cell->SetDau1(d2);
  Explore( (m_Cells[d1]) );
  Explore( (m_Cells[d2]) );
  // CLEANUP
  if(kVer1 != NULL)  delete [] kVer1;
  if(kVer2 != NULL)  delete [] kVer2;
  if(Posi1 != NULL)  delete Posi1;
  if(Posi2 != NULL)  delete Posi2;
  if(Size1 != NULL)  delete Size1;
  if(Size2 != NULL)  delete Size2;
  return 1;
} // TFOAM_Divide


///////////////////////////////////////////////////////////////////////////////
void TFOAM::MakeActiveList(void){
  //  Finds out number of active cells m_NoAct,
  //  Creates list of active cell m_CellsAct and primary cumulative m_PrimAcu
  //  They are used during MC to choose randomly an active cell
  long n, iCell;
  double sum;
  // fush previous result
  if(m_PrimAcu  != NULL) delete [] m_PrimAcu;
  if(m_CellsAct != NULL) delete [] m_CellsAct;
  // Count Active cells and find total Primary
  m_Prime = 0.0; n = 0;
  for(iCell=0; iCell<=m_LastCe; iCell++){
    if (m_Cells[iCell]->GetStat()==1){
      m_Prime += m_Cells[iCell]->GetPrim();
      n++;
    }
  }
  // Allocate tables of active cells
  m_NoAct = n;
  m_CellsAct = new TFCELL*[m_NoAct]; // list of pointers
  m_PrimAcu  = new  double[m_NoAct]; // acumulative primary for MC generation
  if( m_CellsAct==NULL || m_PrimAcu==NULL ) StopM("MakeActiveList - Cant allocate m_CellsAct or m_PrimAcu");

  // Fill-in tables of active cells
  n = 0;
  for(iCell=0; iCell<=m_LastCe; iCell++){
    if( m_Cells[iCell]->GetStat() == 1 ){
      m_CellsAct[n] = m_Cells[iCell]; // ?????
      n++;
    }
  }
  sum =0.0;
  for(iCell=0; iCell<m_NoAct; iCell++){
    sum = sum +(m_CellsAct[iCell])->GetPrim()/m_Prime;
    m_PrimAcu[iCell]=sum;
  }
  
} //MakeActiveList


///////////////////////////////////////////////////////////////////////////////
void TFOAM::GenerCell(TFCELL *&pCell){
//  Generation phase.
//  Return randomly chosen active cell with probability equal to its
//  contribution into total driver integral using binary search .
  long  lo=0, hi=m_NoAct-1, hit;
  double Random;

  Random=m_PseRan->Rndm();
  while(lo+1<hi){
    hit = (lo+hi)/2;
    if (m_PrimAcu[hit]>Random)
      hi = hit;
    else
      lo = hit;
  }
  if (m_PrimAcu[lo]>Random)
    pCell = m_CellsAct[lo];
  else
    pCell = m_CellsAct[hi];
}       // TFOAM::GenerCell

///////////////////////////////////////////////////////////////////////////////
void TFOAM::GenerCel2(TFCELL *&pCell){
// Generation phase
// Return randomly chosen active cell with probability equal to its
// contribution into total driver integral using interpolation search
  long  lo, hi, hit;
  double fhit, flo, fhi;
  double Random;

  Random=m_PseRan->Rndm();
  lo  = 0;              hi =m_NoAct-1;
  flo = m_PrimAcu[lo];  fhi=m_PrimAcu[hi];
  while(lo+1<hi){
    hit = lo + (int)( (hi-lo)*(Random-flo)/(fhi-flo)+0.5);
    if (hit<=lo)
      hit = lo+1;
    else if(hit>=hi)
      hit = hi-1;
    fhit=m_PrimAcu[hit];
    if (fhit>Random){
      hi = hit;
      fhi = fhit;
    }else{
      lo = hit;
      flo = fhit;
    }
  }
  if (m_PrimAcu[lo]>Random)
    pCell = m_CellsAct[lo];
  else
    pCell = m_CellsAct[hi];
}       // TFOAM::GenerCel2


///////////////////////////////////////////////////////////////////////////////
void TFOAM::MakeEvent(void){
// Generation phase.
// Generates point/vector Xrand with the weight MCwt.
  int      j,iv;
  double   Wt,Dx,MCwt;
  TFCELL *rCell;
  TFVECT  VRand(m_nDim),Lambda(m_nDim+1), X(m_nDim);

  int OptGener =2;      // local dip-switch for tests!!!
  //
  //********************** MC LOOP STARS HERE **********************
 ee0:
  if(OptGener==1){
    GenerCell(rCell);   // primitive - more effective for small number of cells (?)
  }else{
    GenerCel2(rCell);   // choose randomly one cell
  }
  MakeLambda();
  MakeAlpha();

    if(!m_OptVert && m_nDim>0 ){
        for(j=0; j<m_nDim; j++) Lambda[j]=m_Lambda[j];
        Lambda[m_nDim]=0.;
        VRand=0.;
        rCell->GetXSimp(VRand,Lambda,m_nDim);
        for(j=0; j<m_nDim; j++) 
           m_MCvect[m_N0Si+j]=VRand[j];
    }else{
      for(j=0 ;  j<m_nDim ; j++){
      m_MCvect[m_N0Si+j] = (*(*rCell)[m_nDim])[j];
      for(iv=0 ; iv <m_nDim ; iv++)
         m_MCvect[m_N0Si+j] += m_Lambda[iv] *( (*(*rCell)[iv])[j] - (*(*rCell)[m_nDim])[j] );
       }
    }

  TFVECT  Posi(m_kDim); TFVECT  Size(m_kDim);
  rCell->GetHcub(Posi,Size);
  for(j=0; j<m_kDim; j++)
    m_MCvect[m_N0Cu+j]= Posi[j] +m_Alpha[j]*Size[j];
  Dx = rCell->GetVolume();      // Cartesian volume of the Cell
  //  weight average normalised to PRIMARY integral over the cell
  Wt = m_Rho->Density(m_TotDim, m_MCvect )*Dx;
  MCwt = Wt / rCell->GetPrim();  // PRIMARY controls normalization
  m_nCalls++;
  m_MCwt   =  MCwt;
  // accumulation of statistics for the main MC weight
  m_SumWt  += MCwt;           // sum of Wt
  m_SumWt2 += MCwt*MCwt;      // sum of Wt**2
  m_NevGen++;                 // sum of 1d0
  m_WtMax  =  dmax(m_WtMax, MCwt);   // maximum wt
  m_WtMin  =  dmin(m_WtMin, MCwt);   // minimum wt
  m_MCMonit->Fill(MCwt);
  m_HistWt->Fill(MCwt,1.0);          // histogram
  //*******  Optional rejection ******
  if(m_OptRej == 1){
    double Random;
    Random=m_PseRan->Rndm();
    if( m_MaxWtRej*Random > m_MCwt) goto ee0;  // Wt=1 events, internal rejection 
    if( m_MCwt<m_MaxWtRej ){
      m_MCwt = 1.0;                  // normal Wt=1 event
    }else{
      m_MCwt = m_MCwt/m_MaxWtRej;    // weight for overveighted events! kept for debug
      m_SumOve += m_MCwt-m_MaxWtRej; // contribution of overveighted
    }
  }
  //********************** MC LOOP ENDS HERE **********************
} // MakeEvent

///////////////////////////////////////////////////////////////////////////////
void TFOAM::GetMCvect(double *MCvect){
//  Returns generated point/vector 
  for ( int k=0 ; k<m_TotDim ; k++) *(MCvect +k) = m_MCvect[k];
}//GetMCvect
///////////////////////////////////////////////////////////////////////////////
double TFOAM::GetMCwt(void){
//   Returns point/vector MCvect weight MCwt
  return(m_MCwt);
}
///////////////////////////////////////////////////////////////////////////////
void TFOAM::GetMCwt(double &MCwt){
//   Returns point/vector MCvect weight MCwt
  MCwt=m_MCwt;
}

///////////////////////////////////////////////////////////////////////////////
double TFOAM::MCgenerate(double *MCvect){
// Generates event and returns weight MCvect
   MakeEvent();
   GetMCvect(MCvect);
   return(m_MCwt);
}//MCgenerate

///////////////////////////////////////////////////////////////////////////////
void TFOAM::GetIntegMC(double &MCresult, double &MCerror){
// Provides integral resulting calculated from averages of the MC run
// To be called after (or during) MC run
  double MCerelat;
  MCresult = 0.0;
  MCerelat = 1.0;
  if (m_NevGen>0){
    MCresult = m_Prime*m_SumWt/m_NevGen;
    MCerelat = sqrt( m_SumWt2/(m_SumWt*m_SumWt) - 1/m_NevGen);
  }
  MCerror = MCresult *MCerelat;
}//GetIntegMC

///////////////////////////////////////////////////////////////////////////////
void  TFOAM::GetIntNorm(double& IntNorm, double& Errel ){
// Returns NORMALIZATION integral
// To be called after (or during) MC run
  if(m_OptRej == 1){    // Wt=1 events, internal rejection
    double IntMC,ErrMC;
    GetIntegMC(IntMC,ErrMC);
    IntNorm = IntMC;
    Errel   = ErrMC;
  }else{                // Wted events, NO internal rejection
    IntNorm = m_Prime;
    Errel   = 0;
  }
}//GetIntNorm

/////////////////////////////////////////////////////////////////////////////////////////
void  TFOAM::GetWtParams(const double eps, double &AveWt, double &WtMax, double &Sigma){
// Returns parameters of the MC weight for efficiency evaluation
  double MCeff, WtLim;
  m_MCMonit->GetMCeff(eps, MCeff, WtLim);
  WtMax = WtLim;
  AveWt = m_SumWt/m_NevGen;
  Sigma = sqrt( m_SumWt2/m_NevGen -AveWt*AveWt );
}//GetMCeff

/////////////////////////////////////////////////////////////////////////////////////////
void TFOAM::Finalize(double& IntNorm, double& Errel ){
// Finalization phase.
// After MC run is completed it provides normalization and
// also prints some information/statistics on the MC run.
  GetIntNorm(IntNorm,Errel);
  double MCresult,MCerror;
  GetIntegMC(MCresult,MCerror);
  double MCerelat= MCerror/MCresult;
  //
  if(m_Chat>0){
    double eps = 0.0005;
    double MCeff, MCef2, WtMax, AveWt, Sigma;
    GetWtParams(eps, AveWt, WtMax, Sigma);
    MCeff=0;
    if(WtMax>0.0) MCeff=AveWt/WtMax;
    MCef2 = Sigma/AveWt;
    double Driver = m_Cells[0]->GetDriv();
    //
    BXOPE;
    BXTXT("****************************************");
    BXTXT("******      TFOAM::Finalize       ******");
    BXTXT("****************************************");
    BX1I("    NevGen",m_NevGen, "Number of generated events in the MC generation   ");
    BX1I("    LastVe",m_LastVe, "Number of vertices (only for simplical option)    ");
    BX1I("    nCalls",m_nCalls, "Total number of function calls                    ");
    BXTXT("----------------------------------------");
    BX1F("     AveWt",AveWt,    "Average MC weight                      ");
    BX1F("     WtMin",m_WtMin,  "Minimum MC weight (absolute)           ");
    BX1F("     WtMax",m_WtMax,  "Maximum MC weight (absolute)           ");
    BXTXT("----------------------------------------");
    BX1F("    XPrime",m_Prime,  "Primary total integral, R_prime        ");
    BX1F("    XDiver",Driver,   "Driver  total integral, R_loss         ");
    BXTXT("----------------------------------------");
    BX2F("    IntMC", MCresult,  MCerror,      "Result of the MC Integral");
    BX1F(" MCerelat", MCerelat,  "Relative error of the MC intgral       ");
    BX1F(" <w>/WtMax",MCeff,     "MC efficiency, acceptance rate");
    BX1F(" Sigma/<w>",MCef2,     "MC efficiency, variance/ave_wt");
    BX1F("     WtMax",WtMax,     "WtMax(esp= 0.0005)            ");
    BX1F("     Sigma",Sigma,     "variance of MC weight         ");
    if(m_OptRej==1){
       double AvOve=m_SumOve/m_SumWt;
    BX1F("<OveW>/<W>",AvOve,     "Contrib. of events wt>MaxWtRej");
    }
    BXCLO;
  }
}  // Finalize

///////////////////////////////////////////////////////////////////////////////
void  TFOAM::SetInhiDiv(int iDim, int InhiDiv){
// This should be called before Initialize, after setting nDim and kDim
  if(m_kDim==0) StopM("TFOAM::SetInhiDiv: m_kDim=0");
  if(m_InhiDiv == NULL){
    m_InhiDiv = new int[ m_kDim ];
    for(int i=0; i<m_kDim; i++) m_InhiDiv[i]=0;
  }
  //
  if( ( 0<=iDim) && (iDim<m_kDim)){
    m_InhiDiv[iDim] = InhiDiv;
  }else
    StopM("TFOAM::SetInhiDiv: wrong iDim");
}//SetInhiDiv

///////////////////////////////////////////////////////////////////////////////
void  TFOAM::SetXdivPRD(int iDim, int len, double xDiv[]){
// This should be called before Initialize, after setting nDim and kDim
// 0<=iDim<_kDim is in H-cubic subspace!!!
  if(m_kDim<=0)  StopM("TFOAM::SetXdivPRD: m_kDim=0");
  if(   len<1 )  StopM("TFOAM::SetXdivPRD: len<1");
  // allocate list of pointers, if it was not done before
  if(m_XdivPRD == NULL){
    m_XdivPRD = new TFVECT*[m_kDim];
    for(int i=0; i<m_kDim; i++)  m_XdivPRD[i]=NULL; // only for hcubic subspace.
  }
  // set division list for direction iDim in H-cubic subspace!!!
  if( ( 0<=iDim) && (iDim<m_kDim)){
    m_OptPRD =1;      // !!!!
    if( m_XdivPRD[iDim] != NULL)
      StopM("TFOAM::SetXdivPRD: second alocation of XdivPRD not allowed");
    m_XdivPRD[iDim] = new TFVECT(len); // allocate list of division
    for(int i=0; i<len; i++){
      (*m_XdivPRD[iDim])[i]=xDiv[i]; // set list of division points
    }
  }else
    StopM("TFOAM::SetXdivPRD: wrong iDim");
  //[[[[[[[[[[[[[[[[[[[[[
  cout<<" SetXdivPRD, idim= "<<iDim<<"  len= "<<len<<"   "<<endl;
  for(int i=0; i<len; i++){
    cout<< (*m_XdivPRD[iDim])[i] <<"  ";
  }
  cout<<endl;
  for(int i=0; i<len; i++)  cout<< xDiv[i] <<"   ";
  cout<<endl;
  //]]]]]]]]]]]]]]]]]]]]]
}//SetXdivPRD

///////////////////////////////////////////////////////////////////////////////
void TFOAM::CheckAll(const int level){
//  Miscelaneous and debug.
//  Checks all pointers, this is useful autodiagnostic.
//  level=0, no printout, failures causes STOP
//  level=1, printout, failures lead to WARNINGS only
  int Errors, Warnings, Ref;
  TFCELL *Cell;
  long iCell,iVert,jCe;
  int k;

  Errors = 0; Warnings = 0;
  if (level==1) cout << "///////////////////////////// FOAM_Checks /////////////////////////////////" << endl;
  m_NoAct=0;
  for(iCell=1; iCell<=m_LastCe; iCell++){ // without root
    Cell = m_Cells[iCell];
    if (Cell->GetStat()==1)
      m_NoAct++;
//  checking general rules
    if( ((Cell->GetDau0()==NULL) && (Cell->GetDau1()!=NULL) ) ||
        ((Cell->GetDau1()==NULL) && (Cell->GetDau0()!=NULL) ) ){
      Errors++;
      if (level==1) cout <<  "ERROR: Cell's no " << iCell << " has only one doughter " << endl;
    }
    if( (Cell->GetDau0()==NULL) && (Cell->GetDau1()==NULL) && (Cell->GetStat()==0) ){
      Errors++;
      if (level==1) cout <<  "ERROR: Cell's no " << iCell << " has no dougter and is inactive " << endl;
    }
    if( (Cell->GetDau0()!=NULL) && (Cell->GetDau1()!=NULL) && (Cell->GetStat()==1) ){
      Errors++;
      if (level==1) cout <<  "ERROR: Cell's no " << iCell << " has two dougters and is active " << endl;
    }

// checking parents
    if( (Cell->GetPare())!=m_Cells[0] ){ // not child of the root
      if ( (Cell != Cell->GetPare()->GetDau0()) && (Cell != Cell->GetPare()->GetDau1()) ){
        Errors++;
        if (level==1) cout <<  "ERROR: Cell's no " << iCell << " parent not pointing to this cell " << endl;
      }
    }

// checking doughters
    if(Cell->GetDau0()!=NULL){
      if(Cell != (Cell->GetDau0())->GetPare()){
        Errors++;
        if (level==1) cout <<  "ERROR: Cell's no " << iCell << " doughter 0 not pointing to this cell " << endl;
      }
    }
    if(Cell->GetDau1()!=NULL){
      if(Cell != (Cell->GetDau1())->GetPare()){
        Errors++;
        if (level==1) cout <<  "ERROR: Cell's no " << iCell << " doughter 1 not pointing to this cell " << endl;
      }
    }
  }// loop after cells;

// check on vertices
  if(m_nDim>0){
    for(iVert=0; iVert<=m_LastVe; iVert++){
      Ref = 0; //False
      for(jCe=0; jCe<=m_LastCe; jCe++){
	Cell = m_Cells[jCe];
	int nVert = Cell->GetnVert();
	for(k=0; k<nVert; k++){
	  if( (*Cell)[k]==m_VerX[iVert] ){
	    Ref = 1; break; // True;
	  }
	}
	if (Ref==1) break;
      }
      if (Ref==0){
	Errors++;
	if (level==1) cout << "ERROR: Vertex no. " << iVert << " NOT referenced!" << endl;
      }
    }
  }

// Check for empty cells
  for(iCell=0; iCell<=m_LastCe; iCell++){
    Cell = m_Cells[iCell];
    if( (Cell->GetStat()==1) && (Cell->GetDriv()==0) ){
      Warnings++;
      if(level==1) cout << "Warning: Cell no. " << iCell << "is active but empty " << endl;
    }
  }
// summary
  if(level==1){
    cout << "TFOAM::Check has found " << Errors << " errors and " << Warnings << " warnings" << endl;
    cout << "/////////////////////////////////////////////////////////////////////////" << endl;
  }
  if(Errors>0){
    cout << "STOP in TFOAM::Check - found total " << Errors << " errors " << endl;
    exit(0);
  }
} // Check

///////////////////////////////////////////////////////////////////////////////

void TFOAM::PrintVertices(void){
// Prints ALL vertices of the FOAM.
  long iVert;

  if(m_nDim>0){
    for(iVert=0; iVert<=m_LastVe; iVert++){
      cout<<"Vert["<<iVert <<"]=  ";
      m_VerX[iVert]->Print();
      //cout<<"  {"<<m_VerX[iVert]<<"} ";   // extra DEBUG
      cout<<endl;
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
void TFOAM::PrintCells(void){
// Prints ALL cells of the FOAM
  long iCell;

  for(iCell=0; iCell<=m_LastCe; iCell++){
    cout<<"Cell["<<iCell<<"]={ ";
    //cout<<"  "<< m_Cells[iCell]<<"  ";  // extra DEBUG
    cout<<endl;
    m_Cells[iCell]->Print();
    cout<<"}"<<endl;
  }
}

///////////////////////////////////////////////////////////////////////////////
void TFOAM::RootPlot2dim(const char *filename){
// Debugging tool which plots cells as triangles or rectangles
// in C++ in format readible for root 
  ofstream outfile(filename, ios::out);
  double   x1,y1,x2,y2,x3,y3,x,y;
  long    iCell, iVe;
  long    k;
  double offs =0.1;
  double lpag   =1-2*offs;
  outfile<<"{" << endl;
  outfile<<"cMap = new TCanvas(\"Map1\",\" Cell Map \",600,600);"<<endl;
  //
  outfile<<"TBox*a=new TBox();"<<endl;
  outfile<<"a->SetFillStyle(0);"<<endl;  // big frame
  outfile<<"a->SetLineWidth(4);"<<endl;
  outfile<<"a->SetLineColor(2);"<<endl;
  outfile<<"a->DrawBox("<<offs<<","<<offs<<","<<(offs+lpag)<<","<<(offs+lpag)<<");"<<endl;
  //
  outfile<<"TText*t=new TText();"<<endl;  // text for numbering
  outfile<<"t->SetTextColor(4);"<<endl;
  if(m_LastCe<51)
    outfile<<"t->SetTextSize(0.025);"<<endl;  // text for numbering
  else if(m_LastCe<251)
    outfile<<"t->SetTextSize(0.015);"<<endl;
  else
    outfile<<"t->SetTextSize(0.008);"<<endl;
  //
  outfile<<"TBox*b=new TBox();"<<endl;  // single cell
  outfile <<"b->SetFillStyle(0);"<<endl;
  //
  if(m_kDim==2 && m_LastCe<=2000){
    TFVECT  Posi(m_kDim); TFVECT  Size(m_kDim);
    outfile << "// =========== Rectangular cells  ==========="<< endl;
    for(iCell=1; iCell<=m_LastCe; iCell++){
      if( m_Cells[iCell]->GetStat() == 1){
	m_Cells[iCell]->GetHcub(Posi,Size);
	x1 = offs+lpag*(        Posi[0]); y1 = offs+lpag*(        Posi[1]);
	x2 = offs+lpag*(Posi[0]+Size[0]); y2 = offs+lpag*(Posi[1]+Size[1]);
	//     cell rectangle
	if(m_LastCe<=2000)
	  outfile<<"b->DrawBox("<<x1<<","<<y1<<","<<x2<<","<<y2<<");"<<endl;
	//     cell number
	if(m_LastCe<=250){
	  x = offs+lpag*(Posi[0]+0.5*Size[0]); y = offs+lpag*(Posi[1]+0.5*Size[1]);
	  outfile<<"t->DrawText("<<x<<","<<y<<","<<"\""<<iCell<<"\""<<");"<<endl;
	}
      }
    }
    outfile<<"// ============== End Rectangles ==========="<< endl;
  }//kDim=2
  //
  if(m_nDim==2 && m_LastCe<=2000){
  if(m_OptVert==0) StopM("RootPlot2dim:: You can only use this method with OptVert=1");
    int  *NoRefsAc = new int[m_vMax];
    outfile<<"// =========== Vertices Vertices ==========="<< endl;
    outfile<<"TMarker*mstar=new TMarker(0.0,0.0,30);"<< endl;
    outfile<<"TMarker*mdisc=new TMarker(0.0,0.0,20);"<< endl;
    outfile<<"mdisc->SetMarkerColor(3);"<< endl;
    //Count references of vertices
    for(iVe = 0 ; iVe<=m_LastVe ; iVe++) NoRefsAc[iVe] = 0;
    for(iVe = 0 ; iVe<=m_LastVe ; iVe++)
      for(iCell = 0 ; iCell<m_LastCe ; iCell++){
	TFCELL *Cell = m_Cells[iCell];
	for(k = 0 ; k<m_nDim+1 ; k++){
	  if( Cell->GetStat() == 1)
	    if( (*Cell)[k] == m_VerX[iVe] ) NoRefsAc[iVe]++;
	}
      }
    //Loop over all vertices
    for(iVe = 0 ; iVe<=m_LastVe ; iVe++){
      x = offs+lpag*((*( m_VerX[iVe]))[0]);
      y = offs+lpag*((*( m_VerX[iVe]))[1]);
      if( NoRefsAc[iVe] <= 2 )
	outfile<<"mstar->DrawMarker("<<x<<","<<y<<");"<<endl;
      else
	outfile<<"mdisc->DrawMarker("<<x<<","<<y<<");"<<endl;
      if(m_LastVe<=250)
	outfile<<"t->DrawText("<<x+0.02<<","<<y<<","<<"\""<<iVe<<"\""<<");"<<endl;
    }
    outfile<<"// ============== End Vertices ============="<< endl;
    delete [] NoRefsAc;
    //
    outfile<<"// ==========Triangle  Cells  ==============="<<endl;
    int iCeStart=1;
    if(m_OptOrd==1) iCeStart=0;
    for(iCell=iCeStart;  iCell<=m_LastCe ; iCell++){
      TFVECT *v1 = (*m_Cells[iCell])[0];
      TFVECT *v2 = (*m_Cells[iCell])[1];
      TFVECT *v3 = (*m_Cells[iCell])[2];
      x1 = offs+lpag*(( *v1 )[0]);
      y1 = offs+lpag*(( *v1 )[1]);
      x2 = offs+lpag*(( *v2 )[0]);
      y2 = offs+lpag*(( *v2 )[1]);
      x3 = offs+lpag*(( *v3 )[0]);
      y3 = offs+lpag*(( *v3 )[1]);
      outfile<<"TGraph *triang = new TGraph(4);"<<endl;
      outfile<<"triang->SetPoint(0,"<< x1 <<","<< y1 <<");"<<endl;
      outfile<<"triang->SetPoint(1,"<< x2 <<","<< y2 <<");"<<endl;
      outfile<<"triang->SetPoint(2,"<< x3 <<","<< y3 <<");"<<endl;
      outfile<<"triang->SetPoint(3,"<< x1 <<","<< y1 <<");"<<endl;
      outfile<<"triang->Draw(\"L\");"<<endl;
    }
    outfile << "// =========== End Triangle  Cells ===========" << endl;
  }//nDim=2
  //
  outfile << "}" << endl;
  outfile.close();
}

///////////////////////////////////////////////////////////////////////////////
void TFOAM::LaTexPlot2dim(const char *filename){
// Debugging tool which plots cells as triangles or rectangles
// in the LeTeX picture environment. 
  int     icont = 0;
  long    iCell, iVe;
  long    k;
  int     kx1,ky1,kx2,ky2,kx3,ky3,kx,ky,lx,ly;
  const char star[]    = "\\makebox(0,0){\\Large\\color{red} $\\star$}";
  const char disc[]    = "\\circle*{20}";
  const char dot[]     = "\\circle*{4}";
  ofstream outfile(filename, ios::out);
  
  icont++;
  if (icont >= 1) outfile << "\\newpage" << endl;
//-------------------------------------------------------------------------
//           Header
//-------------------------------------------------------------------------
  outfile << "\\documentclass[12pt]{article}" << endl;
  outfile << "\\usepackage{color}" << endl;   // <-for colors!!!
  outfile << "\\usepackage{epic}"  << endl;   // <-for extended ploting
  outfile << "\\textwidth  = 16cm" << endl;
  outfile << "\\textheight = 18cm" << endl;
  outfile << "\\pagestyle{empty}"  << endl;
  outfile << "\\begin{document}"   << endl;
  outfile << endl;
//------------------------------!
  outfile << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  outfile << "\\begin{figure}[!ht]" << endl;
  outfile << "\\centering" << endl;
//------------------------------!
// Frames and labels
//------------------------------!
  outfile << "% =========== big frame, title etc. =======" << endl;
  outfile << "\\setlength{\\unitlength}{0.1mm}" << endl;
  outfile << "\\begin{picture}(1600,1600)" << endl;
  outfile << "\\put(0,0){\\framebox(1600,1600){ }}" << endl;

  if(m_kDim==2 && m_LastCe<=2000){
    TFVECT  Posi(m_kDim); TFVECT  Size(m_kDim);
    outfile << "% =========== Rectangulare cells  ===========" << endl;
    outfile << "\\newcommand{\\VD}[2]{\\put(#1,#2){" << dot << "}}" << endl;
    for(iCell=1; iCell<=m_LastCe; iCell++){
      if( m_Cells[iCell]->GetStat() == 1){
	m_Cells[iCell]->GetHcub(Posi,Size);
	kx = (int)(Posi[0]*1600); ky = (int)(Posi[1]*1600);
	lx = (int)(Size[0]*1600); ly = (int)(Size[1]*1600);
	kx1 = kx+lx/2; ky1 = ky+ly/2;
	//     cell rectangle
	if(m_LastCe<=5000) outfile<<"\\put("<<kx <<","<<ky <<
                "){\\framebox("<<lx<<","<<ly<<"){ }}"<< endl;
	//     cell number
	if(m_LastCe<=250) outfile<<"\\put("<<kx1<<","<<ky1<<
	        "){\\makebox(0,0)[b]{\\hbox{\\small\\color{magenta}\\scriptsize "<<iCell<<" }}}"<< endl;
      }
    }
    outfile << "% ============== End Rectangles ===========" << endl;
  }//kDim=2
//-----------------------------------------------------------------------------
  if(m_nDim==2 && m_LastCe<=2000){
    if(m_OptVert==0) StopM("LaTexPlot2dim:: You can only use this method with OptVert=1");
    outfile << "% =========== Vertices Vertices ===========" << endl;
    //Count references of vertices
    int  *NoRefsAc = new int[m_vMax];
    for(iVe = 0 ; iVe<=m_LastVe ; iVe++) NoRefsAc[iVe] = 0;
    for(iVe = 0 ; iVe<=m_LastVe ; iVe++)
      for(iCell = 0 ; iCell<m_LastCe ; iCell++){
	TFCELL *Cell = m_Cells[iCell];
	for(k = 0 ; k<m_nDim+1 ; k++){
	  if( Cell->GetStat() == 1)
	    if( (*Cell)[k] == m_VerX[iVe] ) NoRefsAc[iVe]++;
	}
      }
    //Plotting symbol
    outfile << "\\newcommand{\\VD}[2]{\\put(#1,#2){" << disc << "}}" << endl;
    outfile << "\\newcommand{\\VS}[2]{\\put(#1,#2){" << star << "}}" << endl;
    outfile << "\\newcommand{\\VN}[3]{\\put(#1,#2){\\makebox(0,0)[b]{\\hbox{\\small\\color{red} #3}}}}" << endl;
    //Loop[ over all vertices
    for(iVe = 0 ; iVe<=m_LastVe ; iVe++){
      kx = (int)((*( m_VerX[iVe]))[0]*1600);
      ky = (int)((*( m_VerX[iVe]))[1]*1600);
      if( NoRefsAc[iVe] <= 2 )
	outfile << "\\VD{" << kx << "}{" << ky << "}" << endl;
      else
	outfile << "\\VS{" << kx << "}{" << ky << "}" << endl;
      if(m_LastVe<=250) outfile << "\\VN{" << kx-8 << "}{" << ky+12 << "}{" << iVe << "}" << endl;
    }
    outfile << "% ============== End Vertices ===========" << endl;
    delete [] NoRefsAc;
  }
//-----------------------------------------------------------------------------
  if(m_nDim==2 && m_LastCe<=2000){
    outfile << "% ==========Triangle  Cells  ===============" << endl;
    int iCeStart=1;
    if(m_OptOrd==1) iCeStart=0;
    for(iCell=iCeStart;  iCell<=m_LastCe ; iCell++){
      TFVECT *v1 = (*m_Cells[iCell])[0];
      TFVECT *v2 = (*m_Cells[iCell])[1];
      TFVECT *v3 = (*m_Cells[iCell])[2];
      kx1 = (int)( ( *v1 )[0] *1600 );
      ky1 = (int)( ( *v1 )[1] *1600 );
      kx2 = (int)( ( *v2 )[0] *1600 );
      ky2 = (int)( ( *v2 )[1] *1600 );
      kx3 = (int)( ( *v3 )[0] *1600 );
      ky3 = (int)( ( *v3 )[1] *1600 );
      kx= (kx1+kx2+kx3)/3;
      ky= (ky1+ky2+ky3)/3;
      if( m_Cells[iCell]->GetStat() == 1)
	{
	  outfile << "\\drawline(" << kx1 << "," << ky1 << ")(" << kx2 << "," << ky2 << ")" << endl;
	  outfile << "\\drawline(" << kx2 << "," << ky2 << ")(" << kx3 << "," << ky3 << ")" << endl;
	  outfile << "\\drawline(" << kx3 << "," << ky3 << ")(" << kx1 << "," << ky1 << ")" << endl;
	  if(m_LastCe<=250)
	    outfile << "\\put("      << kx  << "," << ky  
		    << "){\\makebox(0,0)[b]{\\hbox{\\small\\color{magenta} " << iCell << " }}}" << endl;
	}
    }
  }
  // Close frame
  outfile << "% ============== End Cells ===========" << endl;
  
  outfile << "\\end{picture}" << endl;
  outfile << "\\end{figure}" << endl;
  outfile << "\\end{document}" << endl;
  
  outfile.close();
}

///////////////////////////////////////////////////////////////////////////////
void TFOAM::ResetRho(TFOAM_INTEGRAND *Rho){
//  Set new integrand distr.
  cout<<"TFOAM::ResetRho: Resetting distribution"<<endl;
  delete m_Rho;
  m_Rho=Rho;
}

void TFOAM::ResetPseRan(TRandom  *PseRan){
// Set new r.n. generator
  cout<<"TFOAM::ResetPseRan: Resetting random number generator"<<endl;
  delete m_PseRan;
  m_PseRan=PseRan;
}

///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//              End of Class TFOAM                                                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////
