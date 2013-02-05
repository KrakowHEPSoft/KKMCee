#ifndef KoralwMaker_H
#define KoralwMaker_H
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//               CLASS   KoralwMaker                                        //
//                                                                          //
//              Interface  to MC event generator Koralw                     //
//                                                                          //
//     Data members are input/output parameters common to ALL MC events     //
//     Member fuctions organize generation of MC events                     //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#include "TROOT.h"

#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>

#include "KorEvent.h"


class KoralwMaker: public TNamed {
//class KoralwMaker {
public:
  long      NevTot;         //   total numer of events to be generated
  long      m_EvenCounter;  //   event serial counter
  long      m_PrintLimit;   //   print limit for debug
  double    xpar[10001];    //   xpar input/output aprameters of koralw
public:
//------ header of event_lu constructor -------
  KoralwMaker();
  ~KoralwMaker();
public:
/////////////////////////////////////////////////////////////////////////////
//                    generation of basic MC event                         //
//------ read input data
void ReadData( long &ntot );
//------ generator initialization
void Initialize();
//------ generator finalization
void Finalize();
//------ generate single event E
void Generate(KorEvent &E);
////////////////////////////////////////////////////////////////////////////
//        Elements of data analysis to be moved out ??                    //
//------ define jets
void JetDefine(KorEvent &event);
////////////////////////////////////////////////////////////////////////////
//        interface to jetset
void LuGive(char *directive);
//------ lund printout
void LuList(KorEvent &event,  long Level);
////////////////////////////////////////////////////////////////////////////
void FortOpen(long &lunit, char *filename);
void FortClose(long &lunit);
void VarRan(double rn[],long &n );
////////////////////////////////////////////////////////////////////////////
ClassDef(KoralwMaker,1)   // KorEvent  class, <--for dictionary
};
////////////////////////////////////////////////////////////////////////////
#endif
