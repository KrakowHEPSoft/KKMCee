//////////////////////////////////////////////////////////////////////////////
//                     CLASS   TPartLund                                     //
//                                                                          //
// This is class for single particle/object in Lund coding style            //
// It can be also string, jet, initial beam etc.                            //
//                                                                          //
//     COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)                //
//     COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)           //
//     COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)       // 
//                                                                          //
//     COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)                //
//     COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)           //
//     COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)       //
//                                                                          //
//struct {                                                                  //
//  int   mstu[200];                                                        //
//  float paru[200];                                                        //
//  int   mstj[200];                                                        //
//  float parj[200];  } ludat1_ ;                                           //
//                                                                          //
//                 end of description TPartLund                              //
//////////////////////////////////////////////////////////////////////////////
# define sw1 setw(7)
# define sw2 setprecision(10) << setw(18)

#include "TPartLund.h"

ClassImp(TPartLund);
//////////////////////////////////////////////////////////////////////////////
void TPartLund::Print(int mode)
{
// TPartLund printing
  //cout.setf(ios::scientific);

  cout << "[" << fSerial << "] ";
  if(mode == 1 )
    {
      cout << sw1 << fStatus     << sw1 << fFlafor    << sw1 <<  fParent;
      cout << sw1 << fFirstChild << sw1 << fLastChild;
      for ( int k=0; k < 4 ; k++ )   cout << sw2 << fMom[k];
      cout << sw2 << fMass << endl;
    };
}
//////////////////////////////////////////////////////////////////////////////
void TPartLund::ListPrint()
{
// printing list
  TPartLund *actual;
  for (actual = this; actual != NULL; actual = actual->fPrevious)
    (*actual).Print(0);
}
//////////////////////////////////////////////////////////////////////////////
//                   END OF CLASS TPartLund                                  //
//////////////////////////////////////////////////////////////////////////////
