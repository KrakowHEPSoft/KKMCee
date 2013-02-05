*-------------------------------------------
      SAVE
*------------------------------------------------------------------------------
      INTEGER    IdGenYFS3
      PARAMETER (IdGenYFS3 = 6)
*------------------------------------------------------------------------------
      INTEGER           m_iscale
      DOUBLE PRECISION  m_CMSene
      COMMON / c_afb  / m_iscale, m_CMSene
*------------------------------------------------------------------------------
      INTEGER           ninp,nout
      COMMON / inout  / ninp,nout
*------------------------------------------------------------------------------
      INTEGER    imax
      PARAMETER (imax=10000)
      REAL*8     xpar(imax)
      REAL*8     ymin,ymax
*------------------------------------------------------------------------------
* semianalytical
      INTEGER     iSemRef
      INTEGER     iO3best,iO3bestF,iO3bestB,iO3bestAFB
      PARAMETER(
     $  iSemRef    =  1000000+107777,
     $  iO3best    =  1000000+305302,
     $  iO3bestF   =  2000000+305302,
     $  iO3bestB   =  3000000+305302,
     $  iO3bestAFB =  4000000+305302)
*------------------------------------------------------------------------------
* basic histos 
      INTEGER idc,mdc,idb,idi,ids,mds,mdb,nbiv,idWork,idf
      INTEGER kdb
      PARAMETER(
     $  idWork =  999999,
     $  idi = 300000,
     $  idf =  20000,
     $  idc =  10000,
     $  mdc = 110000,
     $  ids  = 40000,
     $  mds  =140000,
     $  idb =  50000,
     $  mdb = 150000,
     $  kdb = 450000)
*------------------------------------------------------------------------------
* Angular distributions, with four cut-offs vmax
      INTEGER iangO2,iangO2x,iangO2xx,iangO2xxx,iangO2prop      !EEX
      INTEGER iangG2,iangG2x,iangG2xx,iangG2xxx,iangG2prop      !CEEX int-ON
      INTEGER iangN2,iangN2x,iangN2xx,iangN2xxx,iangN2prop      !CEEX int-OFF
      PARAMETER(
     $  iangO2prop   = mdc+5073,            !EEX, dsig/dtheta (vvPROP  <0.2775d0)
     $  iangO2       = mdc+4073,            !EEX, dsig/dtheta (vvBARE  <0.9000d0)
     $  iangO2x      = mdc+3073,            !EEX, dsig/dtheta (vvBARE  <0.2775d0)
     $  iangO2xx     = mdc+2073,            !EEX, dsig/dtheta (vvBARE  <0.1900d0)
     $  iangO2xxx    = mdc+1073)            !EEX, dsig/dtheta (vvAleph <0.2775d0)
      PARAMETER(
     $  iangG2prop   = mdc+5203,            !CEEX, dsig/dtheta (vvPROP  <0.2775d0)
     $  iangG2       = mdc+4203,            !CEEX, dsig/dtheta (vvBARE  <0.9000d0)
     $  iangG2x      = mdc+3203,            !CEEX, dsig/dtheta (vvBARE  <0.2775d0)
     $  iangG2xx     = mdc+2203,            !CEEX, dsig/dtheta (vvBARE  <0.1900d0)
     $  iangG2xxx    = mdc+1203)            !CEEX, dsig/dtheta (vvAleph <0.2775d0)
      PARAMETER(
     $  iangN2prop   = mdc+5253,            !CEEX, dsig/dtheta (vvPROP  <0.2775d0), inOFF
     $  iangN2       = mdc+4253,            !CEEX, dsig/dtheta (vvBARE  <0.9000d0), inOFF
     $  iangN2x      = mdc+3253,            !CEEX, dsig/dtheta (vvBARE  <0.2775d0), inOFF
     $  iangN2xx     = mdc+2253,            !CEEX, dsig/dtheta (vvBARE  <0.1900d0), inOFF
     $  iangN2xxx    = mdc+1253)            !CEEX, dsig/dtheta (vvAleph <0.2775d0), inOFF
      INTEGER igAfbG2,     igSigG2
      INTEGER igAfbG2x,    igSigG2x
      INTEGER igAfbG2xx,   igSigG2xx
      INTEGER igAfbG2xxx,  igSigG2xxx
      PARAMETER(
     $  igAfbG2      =  iangG2   +1000000,   igSigG2      =  iangG2   +3000000,
     $  igAfbG2x     =  iangG2x  +1000000,   igSigG2x     =  iangG2x  +3000000,
     $  igAfbG2xx    =  iangG2xx +1000000,   igSigG2xx    =  iangG2xx +3000000,
     $  igAfbG2xxx   =  iangG2xxx+1000000,   igSigG2xxx   =  iangG2xxx+3000000)
      INTEGER igAfb51,     igSig51
      INTEGER igAfb31,     igSig31
      INTEGER igAfb43,     igSig43
      PARAMETER(
     $  igAfb51      =  iangG2   +4000000,   igSig51      =  iangG2   +5000000,
     $  igAfb31      =  iangG2x  +4000000,   igSig31      =  iangG2x  +5000000,
     $  igAfb43      =  iangG2xx +4000000,   igSig43      =  iangG2xx +5000000)
*------------------------------------------------------------------------------
      INTEGER iDizSig,iDizSigNin,iDizSigInt
      INTEGER iDizAfb,iDizAfbNin,iDizAfbInt
      PARAMETER(
     $  iDizSig         = idi+1201,        !Sig Dizet int ON
     $  iDizSigNin      = idi+1251,        !Sig Dizet int OFF
     $  iDizSigInt      = idi+1281,        !Sig Dizet interf
     $  iDizAfb         = idi+2201,        !Afb Dizet int ON
     $  iDizAfbNin      = idi+2251,        !Afb Dizet int OFF
     $  iDizAfbInt      = idi+2281)        !Afb Dizet interf
      INTEGER iMustSig,     iMustAfb,     iMustSigF,    iMustSigB
      INTEGER iMustSigNin,  iMustAfbNin,  iMustSigNinF, iMustSigNinB
      INTEGER iMustSigInt,  iMustAfbInt
      INTEGER iMustDelSigF, iMustDelSigB, iMustRatSigF, iMustRatSigB
      INTEGER iKorzSig,iKorzAfb,iMustCor
      INTEGER iHyb1Sig,iHyb1Afb,iHyb1SigF,iHyb1SigB
      INTEGER iHyb2Sig,iHyb2Afb,iHyb2SigF,iHyb2SigB
      PARAMETER(
     $  iMustSig        = idi+ 200,        ! Sig KoralZ/Mustral with IFI on
     $  iMustAfb        = idi+ 230,        ! Afb
     $  iMustSigF       = idi+ 210,        ! Sig Forward
     $  iMustSigB       = idi+ 220,        ! Sig backward
     $  iMustSigNin     = idi+ 500,        ! IFI off, Sig
     $  iMustAfbNin     = idi+ 530,        ! IFI of,f Afb
     $  iMustSigNinF    = idi+ 510,        ! IFI off, Sig Forward
     $  iMustSigNinB    = idi+ 520,        ! IFI off, Sig backward
     $  iMustSigInt     = idi+ 300,        ! delta-IFI
     $  iMustAfbInt     = idi+ 330,        ! delta-IFI
     $  iMustCor        = idi+ 600,        ! Correction factor CEEX2/Mustral, obsolete
     $  iMustDelSigF    = idi+ 610,        ! IFI off, Sig Forward
     $  iMustDelSigB    = idi+ 620,        ! IFI off, Sig backward
     $  iMustRatSigF    = idi+ 630,        ! IFI off, Sig Forward
     $  iMustRatSigB    = idi+ 640,        ! IFI off, Sig backward
     $  iHyb1Sig        = idi+ 700,        ! Sig Hybrid1
     $  iHyb1Afb        = idi+ 730,        ! Afb
     $  iHyb1SigF       = idi+ 710,        ! Sig Forward
     $  iHyb1SigB       = idi+ 720,        ! Sig Forward
     $  iHyb2Sig        = idi+ 750,        ! Sig Hybrid2
     $  iHyb2Afb        = idi+ 780,        ! Afb
     $  iHyb2SigF       = idi+ 760,        ! Sig Forward
     $  iHyb2SigB       = idi+ 770,        ! Sig Forward
     $  iKorzSig        = idi+ 800,        ! Sig KORALZ/YFS3
     $  iKorzAfb        = idi+ 830)        ! Sig KORALZ/YFS3
      INTEGER iMustAng,iMustAngNint
      INTEGER jMustAfb,jMustSig
      PARAMETER(
     $  iMustAng        = idi+ 7000,        !dSig/dCoThe KoralZ/Mustral
     $  iMustAngNint    = idi+ 7001,        !dSig/dCoThe KoralZ/Mustral
     $  jMustAfb        = idi+ 8000,        !Afb(TheMax) KoralZ/Mustral
     $  jMustSig        = idi+ 8001)        !Sig(TheMax) KoralZ/Mustral
*------------------------------------------------------------------------------
*             Classical v-distributions cumulative
      INTEGER isigO3, isigG0, isigG0nin, isigG0int
      INTEGER         isigG1, isigG1nin, isigG1int
      INTEGER         isigG2, isigG2nin, isigG2int
      INTEGER         isigG2mG1, isigG2mG1nin, isigComb
      INTEGER         isigG2ninF,isigG2ninB
      PARAMETER(
     $  isigO3          = mdb+  74,        !sig(vmax) EEX3 
     $  isigG0          = mdb+ 201,        !sig(vmax) CEEX0 full
     $  isigG0nin       = mdb+ 251,        !sig(vmax) CEEX0 int OFF
     $  isigG0int       = mdb+ 281,        !sig(vmax) CEEX0 interf. only
     $  isigG1          = mdb+ 202,        !sig(vmax) CEEX1 full
     $  isigG1nin       = mdb+ 252,        !sig(vmax) CEEX1 int OFF
     $  isigG1int       = mdb+ 282,        !sig(vmax) CEEX1 interf. only
     $  isigG2          = mdb+ 203,        !sig(vmax) CEEX2 full
     $  isigG2nin       = mdb+ 253,        !sig(vmax) CEEX2 int OFF
     $  isigG2int       = mdb+ 283,        !sig(vmax) CEEX2 interf. only
     $  isigG2mG1       = mdb+1203,        !sig(vmax) CEEX2 CEEX O(alf2-alf1)
     $  isigG2mG1nin    = mdb+1253,        !sig(vmax) CEEX2 CEEX O(alf2-alf1) int OFF
     $  isigG2ninF      = mdb+5253,        !sig(vmax) CEEX2 int OFF Forward
     $  isigG2ninB      = mdb+6253,        !sig(vmax) CEEX2 int OFF Backward
     $  isigComb        = mdb+1903)        !sig(vmax) CEEX2 + IFI_Mustraal
*------------------------------------------------------------------------------
      INTEGER iafbO3, iafbG0, iafbG0nin, iafbG0int
      INTEGER         iafbG1, iafbG1nin, iafbG1int
      INTEGER         iafbG2, iafbG2nin, iafbG2int
      INTEGER         iafbG2mG1, iafbG2mG1nin, iafbComb
      PARAMETER(
     $  iafbO3          = mdb+2074,        !AFB(vmax) EEX O3
     $  iafbG0          = mdb+2201,        !AFB(vmax) CEEX O0
     $  iafbG0nin       = mdb+2251,        !AFB(vmax) CEEX O0
     $  iafbG0int       = mdb+2281,        !AFB(vmax) CEEX O0 interf. only
     $  iafbG1          = mdb+2202,        !AFB(vmax) CEEX O1
     $  iafbG1nin       = mdb+2252,        !AFB(vmax) CEEX O1 int OFF
     $  iafbG1int       = mdb+2282,        !AFB(vmax) CEEX O1 interf. only
     $  iafbG2          = mdb+2203,        !AFB(vmax) CEEX O2
     $  iafbG2nin       = mdb+2253,        !AFB(vmax) CEEX O2 int OFF
     $  iafbG2int       = mdb+2283,        !AFB(vmax) CEEX O2 interf. only
     $  iafbG2mG1       = mdb+3203,        !AFB(vmax) CEEX O2-O1
     $  iafbG2mG1nin    = mdb+3253,        !AFB(vmax) CEEX O2-O1 int OFF
     $  iafbComb        = mdb+3903)        !ABF(vmax) CEEX2 + IFI_Mustraal
*------------------------------------------------------------------------------
*  SEMIREALISTIC  v-like-distributios, all cumulative
      INTEGER         isigS2, isigS2nin, isigS2int
      INTEGER         iafbS2, iafbS2nin, iafbS2int
      PARAMETER(
     $  isigS2          = mds+ 202,        !sig(vmax) CEEX O1 full
     $  isigS2nin       = mds+ 252,        !sig(vmax) CEEX O1 int OFF
     $  isigS2int       = mds+ 282,        !sig(vmax) CEEX O1 interf. only
     $  iafbS2          = mds+2202,        !AFB(vmax) CEEX O1 full
     $  iafbS2nin       = mds+2252,        !AFB(vmax) CEEX O1 int OFF
     $  iafbS2int       = mds+2282)        !AFB(vmax) CEEX O1 interf. only
*---------------------------------------------------------------------------
      INTEGER         ksigG2, ksigG2int, ksigG2nin  ! bin-per-bin
      INTEGER         kafbG2, kafbG2nin, kafbG2int  ! bin-per-bin
      PARAMETER(
     $  ksigG2          = kdb+ 202,        !sig(v   ) CEEX O1 full
     $  ksigG2nin       = kdb+ 252,        !sig(v   ) CEEX O1 int OFF
     $  ksigG2int       = kdb+ 282,        !sig(v   ) CEEX O1 interf. only
     $  kafbG2          = kdb+2202,        !AFB(v   ) CEEX O1 full
     $  kafbG2nin       = kdb+2252,        !AFB(v   ) CEEX O1 int OFF
     $  kafbG2int       = kdb+2282)        !AFB(v   ) CEEX O1 interf. only
*---------------------------------------------------------------------------
* AFB as a function of the cut
      INTEGER iafbG0mO0, iafbG1mO1, iafbG1mG0, iafbG1mG0intoff
      PARAMETER(
     $  iafbG0mO0       = mdb+2271,        !AFB(vmax) CEEX-EEX
     $  iafbG1mO1       = mdb+2272,        !AFB(vmax) CEEX-EEX
     $  iafbG1mG0       = mdb+1202,        !AFB(vmax) CEEX O1-O0 Interf. on
     $  iafbG1mG0intoff = mdb+1252)        !AFB(vmax) CEEX O1-O0 Interf. off
*------------------------------------------------------------------------------
*     Mark plots for plots
      CHARACTER*32 star,diamond,circle,ring,times,disc,plus,box,dot
      PARAMETER (diamond ='\\makebox(0,0){\\Large $\\diamond$}')
      PARAMETER (star    ='\\makebox(0,0){\\Large $\\star$}')
      PARAMETER (circle  ='\\circle{30}')
      PARAMETER (ring    ='\\circle{20}')
      PARAMETER (times   ='\\makebox(0,0){\\Large $\\times$}')
      PARAMETER (disc    ='\\circle*{20}')
      PARAMETER (plus    ='\\makebox(0,0){\\Large $+$}')
      PARAMETER (box     ='\\makebox(0,0){\\Large $\\Box$}') !!! does not work???
      PARAMETER (dot     ='\\circle*{10}')
*---------------------------------------------------------------------------
      CHARACTER*60  TeXfile
      CHARACTER*16  fmtx,fmty
*---------------------------------------------------------------------------
      CHARACTER*80 LabelBasic(6)
      DATA LabelBasic /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\input{LabelMacros}',
     $'\\Pave{',
     $'\\PaveLt{600}{1170}{\\color{green}',
     $'        \\large $\\cal KK$ MC 1999, S.Jadach, Z. W\\c{a}s, B.F.L. Ward}',
     $'} % -- end Pave',
     $'% end-of-label'/
*-------------------------------------------
*---- robol4, table of KF-Cut
      INTEGER           nCut, iCut, nKF, iKF, iBin
      PARAMETER (       nCut = 11, nKF = 16 )
      DOUBLE PRECISION  VmaxTab(nCut)
      DATA VmaxTab / 0.01d0, 0.10d0, 0.20d0, 0.30d0, 0.40d0, 0.50d0,
     $                       0.60d0, 0.70d0, 0.80d0, 0.90d0, 0.99d0/
*-------------------------------------------
      INTEGER                 iKKsem,
     $                        iCeex2,     iCeex2dif,
     $                        iKor402sig, iKor402dif,
     $                        iKor404sig, iKor404dif,
     $                        iZft620sig, iZft620dif
      PARAMETER(   iKKsem      = idf+ 1000)
      PARAMETER(   iCeex2      = idf+ 2000,    iCeex2dif   = idf+ 3000 )
      PARAMETER(   iKor402sig  = idf+ 4000,    iKor402dif  = idf+ 5000 )
      PARAMETER(   iKor404sig  = idf+ 6000,    iKor404dif  = idf+ 7000 )
      PARAMETER(   iZft620sig  = idf+ 8000,    iZft620dif  = idf+ 9000 )
*-------------------------------------------

