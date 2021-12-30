#ifndef TAUOLA_STRUCTS_H
#define TAUOLA_STRUCTS_H

// NOTE: these variables do not belong to any common in FORTRAN
//       (they cannot, as they are used as constants in FORTRAN code)
//       We have to explicitly repeat them here. This is a problem.
//
// WARNING: These values must match definitions in ../TAUDCDsize.inc
const int NLT=2, NMODE=196,NM1=40,NM2=71,NM3=19,NM4=32,NM5=21,NM6=13;

// COMMON / TAUDCD /IDFFIN(9,NMODE),MULPIK(NMODE)
// &                ,NAMES
extern "C" struct {
    int  IDFFIN[NMODE][9];
    int  MULPIK[NMODE];
    char NAMES[NMODE][31];
} taudcd_;
extern "C" struct {
    int  KEY0[2];
    int  KEY1[NM1];
    int  KEY2[NM2];
    int  KEY3[NM3];
    int  KEY4[NM4];
    int  KEY5[NM5];
    int  KEY6[NM6];
} metyp_;
extern "C" struct {
  float GAMPRT[500];
  int JLIST[500];
  int NCHAN;
  // for the C use only, this struct is  supplemented with
  // list of pointers defined in ChannelForTauolaInterface.h:
  // ChannelForTauola* taubra_userChannels[500];
} taubra_;

extern "C" struct {
  float PROB1[NM2];
  float PROB2[NM2];
  float AM2[NM2];
  float GAM2[NM2];
  float AM3[NM2];
  float GAM3[NM2];
} sampl2_;

extern "C" struct {
  float PROB1[NM3];
  float PROB2[NM3];
  float AMRX[NM3];
  float GAMRX[NM3];
  float AMRA[NM3];
  float GAMRA[NM3];
  float AMRB[NM3];
  float GAMRB[NM3];
} sampl3_;

extern "C" struct {
  float PROB1[NM4];
  float PROB2[NM4];
  float AMRX[NM4];
  float GAMRX[NM4];
  float AMRA[NM4];
  float GAMRA[NM4];
} sampl4_;

extern "C" struct {
  float  PROBa2[NM5];
  float  PROBOM[NM5];
  float  ama2[NM5];
  float  gama2[NM5];
  double AMOM[NM5];
  double GAMOM[NM5];
} sampl5_;

#endif
