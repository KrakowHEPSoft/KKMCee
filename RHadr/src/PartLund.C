//////////////////////////////////////////////////////////////////////////////
//                     CLASS   PartLund                                     //
//                                                                          //
// This is class for single particle/object in Lund coding style            //
// It can be also string, jet, initial beam etc.                            //
//                                                                          //
//     COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)                //
//     COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)           //
//     COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)       // 
//struct {                                                                  //
//  int   mstu[200];                                                        //
//  float paru[200];                                                        //
//  int   mstj[200];                                                        //
//  float parj[200];  } ludat1_ ;                                           //
//                                                                          //
//                 end of description PartLund                              //
//////////////////////////////////////////////////////////////////////////////
# define sw1 setw(7)
# define sw2 setprecision(10) << setw(18)

#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>

#include "PartLund.h"

ClassImp(PartLund)
//////////////////////////////////////////////////////////////////////////////
void PartLund::print(int mode)
{
// PartLund printing
  //cout.setf(ios::scientific);

  cout << "[" << m_lserial << "] ";
  if(mode == 1 )
    {
      cout << sw1 << m_kstatus << sw1 << m_kflavor << sw1 <<  m_kparent;
      m_pmom.print();
      cout << sw2 << m_pmass << endl;
    };
}
//////////////////////////////////////////////////////////////////////////////
void PartLund::ListPrint()
{
// printing list
  PartLund *actual;
  for (actual = this; actual != NULL; actual = actual->previous)
    (*actual).print(0);
}
//////////////////////////////////////////////////////////////////////////////
//                   END OF CLASS PartLund                                  //
//////////////////////////////////////////////////////////////////////////////
