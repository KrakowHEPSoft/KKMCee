
      SUBROUTINE DIZET(NPAR,AMW   !  NPAR : FLAGS; AMW : INPUT/OUTPUT
     &                     ,AMZ,AMT,AMH,DAL5H,V_TBA,ALSTR   !  INPUT
     &                     ,ALQED,ALSTRT,ZPAR,PARTZ,PARTW)  !  OUTPUT
*----------------------------------------------------------------------*
* DIZET is prepared to be called by ZFITTER.                           *
* For the program, description, updates etc. see (June 2008):          *
* http://www-zeuthen.desy.de/theory/research/zfitter/index.html 
*
* Update 15 Jan 2013:
* For the program, description, updates etc. see (January 2013):       *
* http://zfitter.com (responsible T. Riemann)
* and 
* http://zfitter.desy.de (S. Riemann, under responsibility of 
* DESY Board of Directors)
*
* For program licence and 'conditions of use' see: 
* http://cpc.cs.qub.ac.uk/licence/licence.html
* and the webpages; in case please contact the spokesperson
*----------------------------------------------------------------------*
* VERSION 6.45 (30 August 2019)
*
* # added JEGERLEHNER(2017) vacuum polarization code 
* # new corrections for SW_EFF implemented 
*----------------------------------------------------------------------*
* VERSION 6.44a (XX Yyyyyyy 2013)
*
* # two minor bugs in function DALPHL in QED vacuum polarization 
*   are fixed (thanks to Satoshi Mishima and A.A.)
* # a bug in O(alpha alpha_s) corrections to vacuum polarizations of  
*   gauge bosons is fixed according to B. Kniel, NPB 347 (1990) 86
* (thanks to Satoshi Mishima)
*----------------------------------------------------------------------*
* VERSION 6.44 (15 January 2013)
* 
* new flag IBAIKOV is used, is set fixed inside subroutine QCDCOF
* for the modifications see: "mod. 14 Jan 2013"
* The changes are made in order to implement the QDC corrections 
* as described in:
* P. Baikov, K. Chetyrkin, J. Kuehn, J. Rittinger
* "Complete QCD Corrections to Hadronic Z-Decays in Order alpha_s^4"
* following exactly: arXiv:1201.5804v3 (2 May 2012)
* published as Phys.Rev.Lett. 108 (2012) 222003
* contact: Tord Riemann tordriemann@gmail.com
*---------------------------------------------------------------------
*  VERSION 6.43 (12 June 2008) 
*  # correction to version 6.42 performed by T. Riemann, DESY, based on
*  an information by  Daisuke Nomura, KEK
*  Email  Sun, 04 May 2008 03:24:11 +0900 (JST)
*  "..I noticed that the length of the 1715-th line of the 
*  file "dizet6_42.f"
* "& ... +C6*DELTOP-C7*DELTOP**2-C8*DELHIG*DELTOP+C9*DELHSQ*DELTOP..."
*  exceeds the (notorious) 72-character limit of fortran 77..."
*  #
*  # Further, the new webpage of zfitter is indicated above
*----------------------------------------------------------------------*
* VERSION 6.42 (18 March 2005)                                         *
*  correction to version 6.41 by A. Freitas                            *
* In the treatment of universal (factorizable) higher-order Zff vertex *
*  corrections for the bb final state, a bug in the form factor rho in *
*  ROKANC was corrected.                                               *
*----------------------------------------------------------------------*
* VERSION 6.41 (15 October 2004)                                       *
*  preliminary changes by A. Freitas                                   *
* Changes made to ROKANC to include universal (factorizable) higher-   *
*  order Zff vertex corrections (for Zbb vertex still missing)         *
*----------------------------------------------------------------------*
* VERSION 6.40 (29 April 2004)                                         *
*  preliminary changes by A. Freitas                                   *
* NPAR(24)= IDMWW is now fully operative to simulate theoretical       *
*  uncertainties for new 2-loop result of MW                           *
*  -> vary IDMWW between -1 and 1 to simulate error range              *
*     (default: IDMWW=0)                                               *
* NPAR(25)= IDSWW is new flag to simulate theoretical                  *
*  uncertainties for new 2-loop result of kappa (sw_eff)               *
*  -> vary IDSWW between -1 and 1 to simulate error range              *
*     (default: IDSWW=0)                                               *
*----------------------------------------------------------------------*
* VERSION 6.36 (21 June 2001)                                          *
* NPAR(24)= IDMWW to incorporate Martin's wish to have MW TU-simulation*
* Attention!!!                                                         *
* Change of argument list compared to all previous versions, even 6_30 *
* VERSION 6.34 (21 January 2001)                                       *
* NPAR(22)= IAMW2 to incorporate Hollik, Weiglein at all results       *
* NPAR(23)= ISFSR to incorporate S.Jadach's wish to have a Switch'FSR  *
*----------------------------------------------------------------------*
* VERSION 6.30 (10  April   2000)                                      *
*--------------------------------                                      *
* A bug in \Gamma_WW is fixed: L. Kalinovskaya, March 2000             *
* V_TB is implemented into Z->bb and ee->bb channels:                  *
*                 L. Kalinovskaya, D. Bardin, A. Olshevsky, March 2000 *
* Warning!!!                                                           *
* DIZET-agument list is changed in order to treat V_TB on the same     *
*                               footing as other parameters to be fit! *
*----------------------------------------------------------------------*
* VERSION 6.05 (06 April 1999)                                         *
*----------------------------------------------------------------------*
* VERSION 5.20 (06 February 1999)                                      *
*----------------------------------------------------------------------*
* VERSION 5.10 (06 March 1998)                                         *
*----------------------------------------------------------------------*
* Recent versions:                                                     *
* version 5_0  - Reactivation of weak BOXES                            *
*       (\BETA)  March-September 1995                                  *
*                Used for "LEP2 Workshop"                              *
* version 4_9  - IMPLEMENTATION OR REACTIVATION OF `WORKING OPTIONS'   *
*                JULY-AUGUST 1994                                      *
*                USED FOR "Z-RESONANCE WORKSHOP"                       *
*----------------------------------------------------------------------*
* Insertions December 1996                                             *
*              - Implementation of more QCD corrections to R_V and R_A *
*                Carefully tested against TOPAZ0 at SQRT(s)=M_Z,15 GeV *
*              - Implementation of "running" Zbb vertex into ROKAPP    *
*                and ROKANC                                            *
* Insertions July 1997                                                 *    
*              - Implementation of Degrassi/Gambino/Sirlin/Vicini      *    
*                NLL corrections O(G^2_mu*m^2_t*M^2_Z)                 *    
*                for \Delta r and \sin^2\theta^{lept}_{eff}            *    
* Insertions February-March 1998                                       *
*              - Implementation of Degrassi/Gambino NLL corrections    *
*                for partial Z-widths and for all realistic observables*
*              - Reactivation of Kniehl's QCD library                  *
*              - Implementation of Czarnecki/Kuehn corrections         *
*======================================================================*
* THE PROGRAM DIZET IS A WEAK LIBRARY. IT CALCULATES ELECTROWEAK       *
* RADIATIVE CORRECTIONS IN THE STANDARD THEORY.                        *
* ---------------------------------------------------------------------*
* FOR THE CALCULATION OF                                               *
*       W-BOSON MASS AND WEAK MIXING ANGLE,                            *
*       Z-BOSON WIDTH,                                                 *
*       W-BOSON WIDTH,                                                 *
* --- CALL DIZET(...).                                                 *
* ---------------------------------------------------------------------*
* FOR THE CALCULATION OF                                               *
*       WEAK FORM FACTORS                                              *
*       FOR 2-FERMION INTO 2-FERMION NEUTRAL CURRENT PROCESSES         *
*       AND OF RUNNING ALPHA.QED,                                      *
* --- CALL ROKANC(...).                                                *
* FOR THE CALCULATION OF                                               *
*       WEAK FORM FACTORS                                              *
*       FOR 2-FERMION INTO 2-FERMION CHARGED CURRENT PROCESSES,        *
* --- CALL RHOCC(...).                                                 *
* IN BOTH CASES, AN EARLIER CALL OF DIZET(...) IS NECESSARY.           *
*======================================================================*
* IF THE CALL WAS 'CALL DIZET(...)':                                   *
*-------------------------------------                                 *
*         FLAGS TO BE SET BY THE USER, EXPLAINED BELOW                 *
*       ------------------------------------------------               *
* NPAR(1) = IHVP                                                       *
* NPAR(2) = IAMT4                                                      *
* NPAR(3) = IQCD                                                       *
* NPAR(4) = IMOMS                                                      *
* NPAR(5) = IMASS                                                      *
* NPAR(6) = ISCRE                                                      *
* NPAR(7) = IALEM                                                      *
* NPAR(8) = IMASK                                                      *
* NPAR(9) = ISCAL                                                      *
* NPAR(10)= IBARB                                                      *
* NPAR(11)= IFTJR                                                      *
* NPAR(12)= IFACR                                                      *
* NPAR(13)= IFACT                                                      *
* NPAR(14)= IHIGS                                                      *
* NPAR(15)= IAFMT                                                      *
* NPAR(16)= IEWLC                                                      *
* NPAR(17)= ICZAK                                                      *
* NPAR(18)= IHIG2                                                      *
* NPAR(19)= IALE2                                                      *
* NPAR(20)= IGFER                                                      *
* NPAR(21)= IDDZZ                                                      *
* NPAR(22)= IAMW2                                                      *
* NPAR(23)= ISFSR                                                      *
* NPAR(24)= IDMWW                                                      *
* NPAR(25)= IDSWW                                                      *
*----------------------------------------------------------------------*
*         INPUT PARAMETERS TO BE SET BY THE USER                       *
*       ------------------------------------------                     *
* AMW  -  W-BOSON MASS (BUT IS BEING CALCULATED FOR NPAR(4)=1)         *
* AMZ  -  Z-BOSON MASS (BUT IS BEING CALCULATED FOR NPAR(4)=2)         *
*    NOTE: DUE TO A POSSIBLE RECALCULATION, THE AMZ, AMW CANNOT BE     *
*         ASSIGNED BY A PARAMETER STATEMENT (INPUT/OUTPUT VARIABLES)   *
* AMT  -  T-QUARK MASS                                                 *
* AMH  -  HIGGS BOSON MASS                                             *
* DAL5H - \Delta_ALPHA^5_{had}(MZ)                                     *
* ALSTR - ALPHA_ST(MZ)                                                 *
*----------------------------------------------------------------------*
*         OUTPUT OF THE DIZET PACKAGE                                  *
*       --------------------------------                               *
* ALSTRT  = ALPHA_ST(MT)                                               *
* ZPAR(1) = DR                                                         *
* ZPAR(2) = DRREM                                                      *
* ZPAR(3) = SW2                                                        *
* ZPAR(4) = GMUC                                                       *
* ZPAR(5-14) = STORES EFFECTIVE SIN'S FOR ALL PARTIAL Z-DECAY CHANNELS *
* 5- NEUTRINO,  6-ELECTRON,  7-MUON, 8-TAU, 9-UP, 10-DOWN, 11-CHARM,   *
* 12-STRANGE , 13-TOP     , 14-BOTTOM.                                 *
* ZPAR(15)= ALPHST                                                     *
* ZPAR(16-30)= QCDCOR(0-14)                                            *
* 
*                                                                      *
* AMW  -  W-BOSON MASS (BUT IS INPUT IF NPAR(4)=2,3)                   *
* AMZ  -  Z-BOSON MASS (BUT IS INPUT IF NPAR(4)=1,3)                   *
* GMUC -  MUON DECAY CONSTANT (IT IS SET TO GMU IF NPAR(4)=1,2)        *
* GMUC -  MUON DECAY CONSTANT (IT IS CALCULATED IF NPAR(4)=3  )        *
*         IF GMU IS CALCULATED FROM AMZ, AMW, IT DEVIATES FROM THE     *
*         EXPERIMENTAL VALUE!                                          *
* DR   -  DELTA.R, THE LOOP CORRECTION TO THE MUON DECAY CONSTANT G.MU *
* DRREM - THE REMAINDER CONTRIBUTION OF THE ORDER ALPHA CALCULATION    *
*         OF DELTA.R AFTER SEPARATION OF THE RESUMMED TERMS            *
* SW2  -  WEAK MIXING ANGLE DEFINED BY WEAK BOSON MASSES               *
* ALPHST - THE QCD COUPLING CONSTANT AS USED IN THE HADRONIC (QUARK)   *
*         DECAY CHANNELS OF THE GAUGE BOSON WIDTHS.                    *
* QCDCOR(I) - QCD CORRECTION FACTOR FOR QUARK PRODUCTION PROCESSES AND *
*         Z - BOSON I-th PARTIAL WIDTHS INTO QUARKS.                   *
*                                                                      *
* PARTZ(I) - PARTIAL DECAY WIDTHS OF THE Z-BOSON FOR THE CHANNELS:     *
*        I=0: NEUTRINO              I= 7: STRANGE                      *
*        I=1: ELECTRON              I= 8: TOP (NOT PART OF WIDTH)      *
*        I=2: MUON                  I= 9: BOTTOM                       *
*        I=3: TAU                   I=10: ALL HADRONS                  *
*        I=4: UP                    I=11: TOTAL                        *
*        I=5: DOWN                                                     *
*        I=6: CHARM                                                    *
* PARTW(I) - PARTIAL DECAY WIDTHS OF THE W-BOSON FOR THE CHANNELS:     *
*        I=1: ONE OF LEPTONIC              I= 2: ONE OF QUARKONIC      *
*        I=3: TOTAL                                                    *
*======================================================================*
* THE OTHER TWO POSSIBLE CALLS HAVE FLAGS, INPUT AND OUTPUT WHICH ARE  *
* COMMENTED AT THE BEGINNING OF THE SUBROUTINES ROKANC AND RHOCC.      *
*======================================================================*
*                                                                      *
* THE FLAGS INTRODUCED ABOVE HAVE THE FOLLOWING MEANING:               *
*----------------------------------------------------------------------*
* CHOICE OF HADRONIC VACUUM POLARISATION:                              *
* IHVP :    IHVP=1: HADR. VACPOL. OF JEGERLEHNER - EIDELMAN REF. 11    *
*      I         2:       OF JEGERLEHNER(1988)                         *
*      I         3:       OF BURKHARDT ET AL., REF. 10                 *
*      I         4:       OF JEGERLEHNER(2016)                         *
*      I         5:       OF JEGERLEHNER(2017)                         *
*----------------------------------------------------------------------*
* HANDLING OF HIGHER ORDER (ALPHA*ALPHA.S) T-MASS TERMS:               *
* IQCD :    THE Z-WIDTH HAS FIXED QCD FACTOR (1+ALFAS/PI).             *
*      I    ADDITIONAL OPTIONS ARE FROM REF.8, WITH RUNNING ALPHA.S    *
*      I    IQCD=0: NO QCD CORRS. TO VECTOR BOSON SELF ENERGIES        *
*      I    IN DELTA-R, WIDTHS, CROSS SECTION                          *
*      I         1: APPROXIM. FAST QCD CORRS. (REALISED FOR  Z-WIDTH   *
*      I            AND LEP PROCESSES)                                 *
*      I            IMPORTANT NOTICE: THESE ARE RELIABLE ONLY FOR LEP-I*
*      I         2: EXACT FORMULAE (FROM BARDIN/CHIZHOV-LIBRARY)       *
*      I         3: EXACT FORMULAE (FROM B. KNIEHL-LIBRARY)            *
*----------------------------------------------------------------------*
* HANDLING OF ALPHA**2 T-MASS TERMS:                                   *
* IAMT4: IAMT4 = 0: NO  MT**4 CORRECTIONS,                             *
*      I       = 1: WITH MT**4 CORRECTIONS RESUMED, SEE REF. 7.        *
*                   IN ORDER TO HAVE COMPLETE BACKWARD COMPATIBILITY,  *
*                   NO COMMON RESUMMATION WITH IQCD.NE.0 TERMS IS DONE *
*      I       = 2: WITH RESUMMATION RECIPE DUE TO BARDIN/KNIEHL,      *
*                   SIMILAR TO THAT OF REF. 12                         *
*      I       = 3: WITH RESUMMATION RECIPE OF REFS. 13-15             *
*      I       = 4: WITH TWO LOOP SUBLEADING CORRECTIONS               *
*                   BY DEGRASSI, GAMBINO at al.                        *
*      I       = 5: WITH TWO LOOP FERMIONIC CORRECTIONS FOR M_W BY     *
*                   WEIGLEIN, HOLLIK ET AL. (JULY 2000)                *
*      I       = 6: WITH COMPLETE TWO LOOP CORRECTIONS FOR M_W AND     *
*                   FERMIONIC TWO-LOOP CORRECTION FOR KAPPA(SW_EFF) BY *
*                   AWRAMIK, CZAKON, FREITAS, WEIGLEIN (APRIL 2004)    *
*      I       = 7: WITH COMPLETE TWO LOOP CORRECTIONS FOR BOTTOM      *
*                   QUARK KAPPA(SW_EFF) BY DYBOVYK AT AL. (JULY 2016)  *
*                   AND OTHER FERMION KAPPA(SW_EFF) BY AWRAMIK AT AL.  *
*                   (NOVEMBER 2006)                                    *
*      I       = 8: WITH COMPLETE TWO LOOP CORRECTIONS FOR             *
*                   KAPPA(SW_EFF) BY DYBOVYK AT AL. (AUGUST 2019)      *
*----------------------------------------------------------------------*
* HANDLING OF HADRONIC VACUUM POLARISATION IN DELTA.R AND RUNNING ALPHA*
* IMASS: IMASS = 0: DEFAULT, USES A FIT TO DATA                        *
*      I       = 1: USES EFFECTIVE QUARK MASSES.OPTION EXISTS FOR TESTS*
*----------------------------------------------------------------------*
* IF IMASK=0: QUARK MASSES ARE USED EVERYWHERE                         *
*         =1: PHYSICAL THRESHOLD ARE USED IN THE PHASE SPACE           *
*----------------------------------------------------------------------*
* IF ISCRE=0: SCALE OF THE REMAINDER TERMS = 1                         *
*         =1: ------------------------------ RENORD                    *
*         =2: ------------------------------ RENORM (not recommended,  *
*                                                    = n.r.)           *
*----------------------------------------------------------------------*
* IF ISCAL=0: SCALET=1                                                 *
*         =1: SCALET=SCLAVR+SCTPLU     |                               *
*         =2: SCALET=SCLAVR            | Kniehl's                      *
*         =3: SCALET=SCLAVR-SCTMIN     |                               *
*         =4: SCALET=SIRLIN'S 0.204                                    *
*----------------------------------------------------------------------*
* IF IFACR=0: NON   -EXPANDED DR                                       *
*         =1: PARTLY-EXPANDED DR, RESPECTING DRL*DRR-ORDER             *
*         =2: PARTLY-EXPANDED DR, VIOLATING  DRL*DRR-ORDER       (n.r.)*
*         =3: FULLY -EXPANDED DR, DOUBLE VIOLATING DRL*DRR-ORDER (n.r.)*
*----------------------------------------------------------------------*
* IF IFACT=0: NON   -EXPANDED RHO AND KAPPA    analogs      ff level   *
*         =1: PARTLY-EXPANDED RHO AND KAPPA       of        ff level   *
*         =2: PARTLY-EXPANDED RHO AND KAPPA     those       ff level   *
*         =3: FULLY -EXPANDED RHO AND KAPPA    in IFACR     ff level   *
*         =4: FULLY -NEGLECTED REM^2 TERMS               width level   *
*         =5:       -FAMOUS EW/QCD NON FACTORIZATION     width level   *
*----------------------------------------------------------------------*
* IHIGS: IHIGS = 0: LEADING HIGGS CONTRIBUTION IS NOT RESUMMED         *
*      I       = 1: LEADING HIGGS CONTRIBUTION IS     RESUMMED         *
*----------------------------------------------------------------------*
* IALEM: IALEM = 0 OR 2: DALH5(AMZ) MUST BE SUPPLIED BY THE USER AS    *
*      I                            INPUT TO THE DIZET PACKAGE         *
*      I       = 1 OR 3: DALH5(AMZ) IS CALCULATED BY THE PROGRAM       *
*      I                            USING A PARAMETRIZATION (IHVP)     *
*----------------------------------------------------------------------*
* IFTJR  IFTJR = 0 WITHOUT FTJR CORRECTIONS                            * 
*        IFTJR > 1  WITH   FTJR CORRECTIONS                            * 
*----------------------------------------------------------------------*
* ICZAK: ICZAK = 0 WITHOUT CZARNECKI/KUEHN CORRECTIONS                 * 
*        ICZAK > 1  WITH   CZARNECKI/KUEHN CORRECTIONS                 * 
*----------------------------------------------------------------------*
* IHIG2: IHIG2 = 0 WITHOUT TWO-LOOP HIGGS  CORRECTIONS                 * 
*        IHIG2 = 1  WITH   --/--/--/--/--/--/--/--/--/                 * 
*----------------------------------------------------------------------*
* IALE2: IALE2 = 0 WITHOUT TWO-LOOP CONSTANT CORRECTIONS IN DELTA_ALPHA* 
*                  THIS IS FOR A BACK COMPATIBILITY ONLY               * 
*        IALE2 = 1 -- WITH   ONE-LOOP CORRECTIONS\                     * 
*        IALE2 = 2 -- WITH   TWO-LOOP CORRECTIONS | FOR LEPTONIC DALPHA* 
*        IALE2 = 3 -- WITH THREE-LOOP CORRECTIONS/                     * 
*----------------------------------------------------------------------*
* IGFER: IGFER = 0 FOR BACK COMPATIBILITY WITH 5.12                    *
*        IGFER = 1 ONE-LOOP QED CORRECTIONS FOR FERMI CONSTANT         *
*        IGFER = 2 TW0-LOOP QED CORRECTIONS FOR FERMI CONSTANT         *
*----------------------------------------------------------------------*
* CHOICE OF INPUT PARAMETERS BESIDES AMT, AMH, ALPHA.S:                *
* IMOMS: IMOMS = 1: (CONVENTIONAL)    ALPHA, GMU, AMZ (OUTPUT: AMW)    *
*      I       = 2: INPUT INSTEAD IS: ALPHA, GMU, AMW (OUTPUT: AMZ)    *
*      I       = 3:                   ALPHA, AMZ, AMW (OUTPUT: GMU)    *
* WHERE                                                                *
*      I GMU...... MUON DECAY CONSTANT                                 *
*      I AMT...... T-QUARK MASS                                        *
*      I AMH...... HIGGS BOSON MASS                                    *
*      I ALST..... STRONG COUPLING CONSTANT ALPHA.S (.11)              *
*======================================================================*
*                                                                      *
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 XFOTF3
*
      COMMON/CDZFLG/IHVP,IAMT4,IQCD,IMOMS,IMASS,IALEM,IMASK,IBARB,IFTJR
      COMMON/CDZSCT/ISCRE,ISCAL,IAFMT,IFACR,IFACT,IHIGS,IEWLC,ICZAK
     &             ,IHIG2,IALE2,IGFER      
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
      COMMON/CDZZWG/CAMZ,CAMH,GMU,A0,GAMZ,GAMW,CALSZ,CALST,CALXI,CALQED
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
*
      COMMON /CDZRKZ/ARROFZ(0:10),ARKAFZ(0:10),ARVEFZ(0:10),ARSEFZ(0:10)
     &              ,AROTFZ(0:10),AIROFZ(0:10),AIKAFZ(0:10),AIVEFZ(0:10)
      COMMON/CDZ513/CDAL5H
      COMMON/CDZIMW/IAMW2
      COMMON/CDZFSR/ISFSR
      COMMON/CDZDMW/IDMWW,IDSWW
*
* Internal COMMON to study \mt-additional terms
*
      COMMON/CDZ_LK/IMTADD
      COMMON/CDZDDZ/IDDZZ
      COMMON/CDZVTB/V_TB
*
      DIMENSION NPAR(25),PARTZ(0:11),PARTW(3),ZPAR(30),QCDCOR(0:14)
*
* Never change this flag
*
      IMTADD=0
      V_TB=V_TBA
*
* FLAGS SETTING
       IHVP=NPAR(1)
      IAMT4=NPAR(2)
       IQCD=NPAR(3)
      IMOMS=NPAR(4)
      IMASS=NPAR(5)
      ISCRE=NPAR(6)
      IALEM=NPAR(7)
      IMASK=NPAR(8)
      ISCAL=NPAR(9)
      IBARB=NPAR(10)
      IFTJR=NPAR(11)
      IFACR=NPAR(12)
      IFACT=NPAR(13)
      IHIGS=NPAR(14)
      IAFMT=NPAR(15)
      IEWLC=NPAR(16)
      ICZAK=MIN(1,NPAR(17))
      IHIG2=NPAR(18)
      IALE2=NPAR(19)
      IGFER=NPAR(20)
      IDDZZ=NPAR(21)
      IAMW2=NPAR(22)
      ISFSR=NPAR(23)
      IDMWW=NPAR(24)
      IDSWW=NPAR(25)
*
      CALL CONST1(IHVP,AMT,AMH)
*
      ALFAS = ALSTR
      AMZ2  = AMZ**2
*
*----------------------------------------------------------------------
* CALCULATION 0F ALPHST=ALPHAS(AMZ**2) AND ALPHTT=ALPHAS(AMT**2)
* FOR THE HIGHER ORDER ALPHA*ALPHA.S CORRECTIONS WITH T-MASS AS USED
* IN THE CALCULATIONS OF DELTA.R AND THE GAUGE BOSON WIDTHS
*
*-----------------------------------------------------------------------
      CAMZ=AMZ
      CAMH=AMH
*-----------------------------------------------------------------------
      IF(MOD(IALEM,2).EQ.0) THEN
        DAL5H=DAL5H
      ELSE
        DAL5H=DALH5(AMZ2,AMZ)     
      ENDIF
      CDAL5H=DAL5H
      DALFA=AL4PI*DREAL(XFOTF3(IALEM,IALE2,IHVP,1,1,DAL5H,-AMZ**2))
      ALQED=1D0/ALFAI/(1D0-DALFA)
*     print *,'DALFA=',DALFA
*     stop
      CALQED=ALQED
*-----------------------------------------------------------------------
      SW2=.2325D0
      ALPHST=ALFAS
      CALL QCDCOF(AMZ,AMT,SW2,ALQED,ALPHST,ALPHTT,ALPHXI,QCDCOR)
      CALSZ=ALPHST
      CALST=ALPHTT
      ALSTRT=ALPHTT
      CALXI=ALPHXI
       ALSZ=ALPHST
       ALST=ALPHTT
       ALSX=ALPHXI
*
* Second Call of XFOTF3, when ALSTR is SET
*
      DALFA=AL4PI*DREAL(XFOTF3(IALEM,IALE2,IHVP,1,1,DAL5H,-AMZ**2))
      ALQED=1D0/ALFAI/(1D0-DALFA)      
      CALQED=ALQED
*      print *,'ALQEDI=',1d0/ALQED
*      stop
*
* ITERATIVE PROCEDURE FOR THE CALCULATION OF IVB- MASSES
*
      CALL SETCON(AMW,AMZ,AMT,AMH,DR,DRBIG,DRREM)
*
      AMW=DSQRT(AMW2)
      SW2=R1
*
      IF(IMOMS.EQ.2) AMZ=AMW/SQRT(1D0-SW2)
*-----------------------------------------------------------------------
* CALCULATION OF FINAL STATE QCD-FACTORS FOR Z,W - DECAYS INTO QUARKS
*-----------------------------------------------------------------------
      IF(IAMT4.EQ.-1) THEN
        IF(ALFAS.LE.1D-10) THEN
          DO IQCDC=0,14
            QCDCOR(IQCDC)=1.000D0
          ENDDO
        ELSE
          QCDCOR(0)=1.000D0
          DO IQCDC=1,10
            QCDCOR(IQCDC)=1.040D0
          ENDDO
          QCDCOR(11)=1.045D0
          QCDCOR(12)=1.045D0
          QCDCOR(13)=0D0
          QCDCOR(14)=0D0
        ENDIF
      ELSE
        CALL QCDCOF(AMZ,AMT,SW2,ALQED,ALPHST,ALPHTT,ALFAXI,QCDCOR)
      ENDIF
*
* CALCULATION OF Z- AND W- WIDTHS
      CALL ZWRATE(DR,DRBIG,DRREM,QCDCOR,V_TB,PARTZ,PARTW)
* FILLING OF OUTPUT VECTOR
      ZPAR(1)=DR
      ZPAR(2)=DRREM
      ZPAR(3)=SW2
      GMUC   =PI/ALFAI/AMW2/SW2/(1D0-DR)*1D5/SQRT(2D0)
      ZPAR(4)=GMUC
      DO 3 IZ=5,14
      ZPAR(IZ)=ARSEFZ(IZ-5)
  3   CONTINUE
      ZPAR(15)=ALPHST
      DO IQCDC=0,14
      ZPAR(16+IQCDC)=QCDCOR(IQCDC)
      ENDDO
*-----------------------------------------------------------------------
      END
 
      FUNCTION DZEWBX(RACOS)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
      COMMON/CDZZWG/AMZ,AMH,GMU,A0,GAMZ,GAMW,CALSZ,CALST,CALXI,CALQED
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      COMMON/CDZFLG/IHVP,IAMT4,IQCD,IMOMS,IMASS,IALEM,IMASK,IBARB,IFTJR
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
      COMMON/CDZBOX/QE,QF,ALAM1,ALAM2,HELI1,HELI2,SC,VE,VF,CRFAC
      COMMON/FORCHI/XKAPP,XKAPPC,XMZ2,XMZ2C
      COMMON/XFORMZ/XFZ(4),XFZT,XFZTA
      COMMON/CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
      COMMON/CDZIBF/IBFLA
*
      S=SC
      GMU2=GMU*GMU
      AMZ4=AMZ2*AMZ2
      AMW4=AMW2*AMW2
      PI3=PI*PI2
      CNST=PI/(2.D0*S*ALFAI*ALFAI)
      CBXWW=CNST*GMU2*AMW4*ALFAI/(4.D0*PI3)
*IAMT4=-1
C       BACKWARD COMPATIBILITY TO YR 1989, LIN.ZWRATE,SOME BUGS
      IF(IAMT4.EQ.-1) THEN
      CBXZZ=CNST*GMU2*AMZ4*ALFAI/(64.D0*PI3)
      ELSE
      CBXZZ=CNST*GMU2*AMZ4*ALFAI/(256.D0*PI3)
      ENDIF
      COMPL=ALAM1+ALAM2
      COMMI=ALAM1-ALAM2
      HOMPL=HELI1+HELI2
      HOMMI=HELI1-HELI2
      QEM=DABS(QE)
      QFM=DABS(QF)
C AVOID THAT AI11,AI12 IS ZERO: MULTIPLY RACOS BY (1 - 2*ME**2/S)
      ACOSM=RACOS*(1.D0-0.0000005D0/S)
      IF (ABS(ACOSM) .EQ. 1D0) ACOSM = DSIGN(0.9999999999999D0,RACOS)
      AI11=S/2.D0*(1.D0+ACOSM)
      AI12=S/2.D0*(1.D0-ACOSM)
      VEP1=VE+1.D0
      VEM1=VE-1.D0
      VFP1=VF+1.D0
      VFM1=VF-1.D0
      VEP2=VEP1*VEP1
      VEM2=VEM1*VEM1
      VFP2=VFP1*VFP1
      VFM2=VFM1*VFM1
      VEP3=VEP1*VEP2
      VEM3=VEM1*VEM2
      VFP3=VFP1*VFP2
      VFM3=VFM1*VFM2
      VZP1=COMPL*HOMPL*VEP2*VFP2+COMMI*HOMMI*VEM2*VFM2
      VZP2=COMPL*HOMPL*VEP3*VFP3+COMMI*HOMMI*VEM3*VFM3
      VZM1=COMPL*HOMMI*VEP2*VFM2+COMMI*HOMPL*VEM2*VFP2
      VZM2=COMPL*HOMMI*VEP3*VFM3+COMMI*HOMPL*VEM3*VFP3
      XCHI=XKAPP*S/(S-XMZ2)
      XCHIG=DCONJG(XCHI)
cbard XFZTG=DCONJG(1.D0/(2.D0-XFZT))
      XFZTG=(1.D0,0.D0)
      SE=1D0
      SF=1D0
      IF(QE.NE.0D0) SE=QE/QEM
      IF(QF.NE.0D0) SF=QF/QFM
*
      XCHIGF=XCHIG*SF
      BOXWW=-COMPL*HOMPL*DREAL(
     & +(QF*XFZTG+XCHIGF*VEP1*VFP1)*
     &  (+(1D0+SE*SF)/2D0*(XBOX(IBFLA,AI11,AI12,AMW2)
     &   +4.D0*(AI11/S)**2*XJ3(S,AMW2))
     &   -(1D0-SE*SF)/2.D0*2.D0*(AI11/S)**2
     &   *(AI11/S*XJ4(S,AI11,AMW2)+2.D0*XJ3(S,AMW2)) ))
      BOXZZ=-DREAL(
     & +(QF*XFZTG*VZP1+XCHIGF*VZP2)
     & *(XBOX(0,AI11,AI12,AMZ2)-(2.D0*(AI11/S)**3)*XJ4(S,AI11,AMZ2))
     & -(QF*XFZTG*VZM1+XCHIGF*VZM2)
     & *(XBOX(0,AI12,AI11,AMZ2)-(2.D0*(AI12/S)**3)*XJ4(S,AI12,AMZ2)))
*
* BOX CROSS SECTION TAKES INTO ACCOUNT
* ARBITRARY BEAM POLARIZATIONS AND DEFINITE FINAL HELICITIES
* COLORF IS TAKEN INTO ACCOUNT IN ZFITTER VIA CORINT
*
      DZEWBX=CRFAC*(CBXWW*BOXWW+CBXZZ*BOXZZ)
*
      END
 
      FUNCTION DZEWBI(RACOS)
* 
* DZEWBI=DZEWBX (Linux problems)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
      COMMON/CDZZWG/AMZ,AMH,GMU,A0,GAMZ,GAMW,CALSZ,CALST,CALXI,CALQED
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      COMMON/CDZFLG/IHVP,IAMT4,IQCD,IMOMS,IMASS,IALEM,IMASK,IBARB,IFTJR
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
      COMMON/CDZBOX/QE,QF,ALAM1,ALAM2,HELI1,HELI2,SC,VE,VF,CRFAC
      COMMON/FORCHI/XKAPP,XKAPPC,XMZ2,XMZ2C
      COMMON/XFORMZ/XFZ(4),XFZT,XFZTA
      COMMON/CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
      COMMON/CDZIBF/IBFLA
*
      S=SC
      GMU2=GMU*GMU
      AMZ4=AMZ2*AMZ2
      AMW4=AMW2*AMW2
      PI3=PI*PI2
      CNST=PI/(2.D0*S*ALFAI*ALFAI)
      CBXWW=CNST*GMU2*AMW4*ALFAI/(4.D0*PI3)
*IAMT4=-1
C       BACKWARD COMPATIBILITY TO YR 1989, LIN.ZWRATE,SOME BUGS
      IF(IAMT4.EQ.-1) THEN
      CBXZZ=CNST*GMU2*AMZ4*ALFAI/(64.D0*PI3)
      ELSE
      CBXZZ=CNST*GMU2*AMZ4*ALFAI/(256.D0*PI3)
      ENDIF
      COMPL=ALAM1+ALAM2
      COMMI=ALAM1-ALAM2
      HOMPL=HELI1+HELI2
      HOMMI=HELI1-HELI2
      QEM=DABS(QE)
      QFM=DABS(QF)
C AVOID THAT AI11,AI12 IS ZERO: MULTIPLY RACOS BY (1 - 2*ME**2/S)
      ACOSM=RACOS*(1.D0-0.0000005D0/S)
      IF (ABS(ACOSM) .EQ. 1D0) ACOSM = DSIGN(0.9999999999999D0,RACOS)
      AI11=S/2.D0*(1.D0+ACOSM)
      AI12=S/2.D0*(1.D0-ACOSM)
      VEP1=VE+1.D0
      VEM1=VE-1.D0
      VFP1=VF+1.D0
      VFM1=VF-1.D0
      VEP2=VEP1*VEP1
      VEM2=VEM1*VEM1
      VFP2=VFP1*VFP1
      VFM2=VFM1*VFM1
      VEP3=VEP1*VEP2
      VEM3=VEM1*VEM2
      VFP3=VFP1*VFP2
      VFM3=VFM1*VFM2
      VZP1=COMPL*HOMPL*VEP2*VFP2+COMMI*HOMMI*VEM2*VFM2
      VZP2=COMPL*HOMPL*VEP3*VFP3+COMMI*HOMMI*VEM3*VFM3
      VZM1=COMPL*HOMMI*VEP2*VFM2+COMMI*HOMPL*VEM2*VFP2
      VZM2=COMPL*HOMMI*VEP3*VFM3+COMMI*HOMPL*VEM3*VFP3
      XCHI=XKAPP*S/(S-XMZ2)
      XCHIG=DCONJG(XCHI)
cbard XFZTG=DCONJG(1.D0/(2.D0-XFZT))
      XFZTG=(1.D0,0.D0)
      SE=1D0
      SF=1D0
      IF(QE.NE.0D0) SE=QE/QEM
      IF(QF.NE.0D0) SF=QF/QFM
*
      XCHIGF=XCHIG*SF
      BOXWW=-COMPL*HOMPL*DREAL(
     & +(QF*XFZTG+XCHIGF*VEP1*VFP1)*
     &  (+(1D0+SE*SF)/2D0*(XBOX(IBFLA,AI11,AI12,AMW2)
     &   +4.D0*(AI11/S)**2*XJ3(S,AMW2))
     &   -(1D0-SE*SF)/2.D0*2.D0*(AI11/S)**2
     &   *(AI11/S*XJ4(S,AI11,AMW2)+2.D0*XJ3(S,AMW2)) ))
      BOXZZ=-DREAL(
     & +(QF*XFZTG*VZP1+XCHIGF*VZP2)
     & *(XBOX(0,AI11,AI12,AMZ2)-(2.D0*(AI11/S)**3)*XJ4(S,AI11,AMZ2))
     & -(QF*XFZTG*VZM1+XCHIGF*VZM2)
     & *(XBOX(0,AI12,AI11,AMZ2)-(2.D0*(AI12/S)**3)*XJ4(S,AI12,AMZ2)))
*
* BOX CROSS SECTION TAKES INTO ACCOUNT
* ARBITRARY BEAM POLARIZATIONS AND DEFINITE FINAL HELICITIES
* COLORF IS TAKEN INTO ACCOUNT IN ZFITTER VIA CORINT
*
      DZEWBI=CRFAC*(CBXWW*BOXWW+CBXZZ*BOXZZ)
*
      END
 
      FUNCTION XJ2(S,AMV2)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
*
      RMV=AMV2/S
      REJ2=S*FJJ(-S,AMV2,AMV2)
      AIJ2=0.D0
      IF(1.D0.LE.4.D0*RMV.OR.S.LT.0D0) GO TO 1
      SLAMS=SQRT(1.D0-4.D0*RMV)
      AIJ2=2.D0*PI/SLAMS
1     XJ2=DCMPLX(REJ2,AIJ2)
      END
 
      FUNCTION XJ3(S,AMV2)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
      DATA EPS/1.D-20/
*
      XI=DCMPLX(0.D0,1.D0)
      XMV2E=AMV2-XI*EPS
      XALF=(-S)/XMV2E
      XS12=SQRT(1.D0+4.D0/XALF)
      XY1=(1.D0-XS12)/2.D0
      XY2=(1.D0+XS12)/2.D0
      XJ3=(LOG(-XY1/XY2))**2
      END
 
      FUNCTION XJ4(S,AI,AMV2)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
      DATA EPS/1.D-20/
*
      XI  =DCMPLX(0.D0,1.D0)
      XMV2=AMV2-XI*EPS
      XALF=(-S)/XMV2
      XBET=  AI/XMV2
      XS12=SQRT(1.D0+4.D0/XALF)
      XY1 =(1.D0-XS12)/2.D0
      XY2 =(1.D0+XS12)/2.D0
      XS34=SQRT(1.D0-4.D0*(1.D0-XBET)/(XALF*XBET))
      XY3 =(1.D0-XS34)/2.D0
      XY4 =(1.D0+XS34)/2.D0
      XNOR=(-XALF)*XBET*XS34
      XJ4A=2.D0*(S/AMV2)**2/XNOR*
     &(XSPENZ(-XY4/(XY2-XY4))-XSPENZ(-XY3/(XY1-XY3))
     &+XSPENZ(-XY4/(XY1-XY4))-XSPENZ(-XY3/(XY2-XY3)))
      REJ4=DREAL(XJ4A)
      AIJ4=0D0
      IF(S.GT.0D0) AIJ4=DIMAG(XJ4A)
      XJ4=DCMPLX(REJ4,AIJ4)
      END
 
      FUNCTION XBOX(IBFLA,AI11,AI12,AMV2)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
      COMMON/CDZFLG/IHVP,IAMT4,IQCD,IMOMS,IMASS,IALEM,IMASK,IBARB,IFTJR
      COMMON/CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
      COMMON/CDZ_LK/IMTADD
      COMMON/CDZVTB/V_TB
*
      S=AI11+AI12
      RI1S=AI11/S
      RI2S=AI12/S
      RMV=AMV2/S
C THE ABS SHOULD BE FURTHER INVESTIGATED
*     XBOX=2.D0*RI1S*LOG(RI2S/RMV)+RI1S*(1.D0-4.D0*RMV)*XJ2(S,AMV2)
*IAMT4=-1
      IF(IAMT4.EQ.-1) THEN
      XBOX=2.D0*RI1S*LOG(ABS(RI2S/RMV))+RI1S*(1.D0-4.D0*RMV)*XJ2(S,AMV2)
     *+2.D0*(RI2S-RI1S-2.D0*RMV)*(F1-SPENCE(1.D0-RI2S)+XJ3(S,AMV2))
     *+(2.D0*RMV**2*RI1S+RI1S**2*RI2S+RI2S*(RI2S-2.D0*RMV)**2)
     **XJ4(S,AI12,AMV2)
      ELSE
      XBOX=2.D0*RI1S*LOG(ABS(RI2S/RMV))+RI1S*(1.D0-4.D0*RMV)*XJ2(S,AMV2)
     * +2.D0*(RI2S-RI1S-2.D0*RMV)*(F1-SPENCE(1.D0-RI2S/RMV)+XJ3(S,AMV2))
     *+(2.D0*RMV**2*RI1S+RI1S**2*RI2S+RI2S*(RI2S-2.D0*RMV)**2)
     **XJ4(S,AI12,AMV2)
*
      IF(IBFLA.EQ.0) RETURN
*
* Genuine WW-box m_t corrections, p.12 of 1997 notes
*
      AMT2=(AMQ(5))**2
      AMW2=AMV2
      RTW =AMT2/AMW2
      ALRT=LOG(RTW)
      RWS =AMW2/S
*
      CALL J3WANA(AMT2,AMW2,AI12,XJ3WT)
      CALL J3WANA( 0D0,AMW2,AI12,XJ3W0)
      CALL S3WANA(AMT2,AMW2,-S,XJ0W,XS3WT,XS3W0)
      CALL S4WANA(AMT2,AMW2,-S,AI12,XS4WT)
      CALL S4WANA( 0D0,AMW2,-S,AI12,XS4W0)
*
      RJ2=AMT2/AI12
*
* WWbg_desc=1/M2W*rt*(RW*(2-1/2*RW-(1/2-RW)*rt-1/2*RW*rt**2)*s*Sw3
*           +1/4*([1/epsilon]-Lnw)
*           -1/2*(1/2-RW+RW*rt)*Jw
*           +1/2*RW*rt*[lnrt]+1/2*(3/2-RW+RW*rt)         
*                     );
*     WWbg =+RTW*(RWS*(2D0-.5D0*RWS-(.5D0-RWS)*RTW
*    &      -.5D0*RWS*RTW**2)*S*S3W     
*    &      -.5D0*(.5D0-RWS+RWS*RTW)*AJ0W
*    &      +.5D0*RWS*RTW*ALRT+.5D0*(1.5D0-RWS+RWS*RTW))
*
      XWWbg =+.5D0*RTW*RWS*(
     &        +(4D0-RWS-(1D0-2D0*RWS)*RTW-RWS*RTW**2)*S*XS3WT     
     &        +(1D0-RTW)*(XJ0W-1D0)+RTW*ALRT
     &                    )
*
      XBBX=AI11**2/S/AMW2*XWWbg
     &     +2D0/S*(
     &     +AI11*(LOG(1D0+RJ2)+RJ2*LOG(1D0+1D0/RJ2))
     &     +(2D0*S*AMW2-AI11**2-AI12**2)/2D0*(XS3WT-XS3W0)
     &     -S*AMT2/2D0*(XS3WT+XS3W0)
     &     +AI12*(AI12-AI11-2D0*AMW2)*(XJ3WT-XJ3W0)+AI12*AMT2*XJ3WT
     &     +(+AMW2*((2D0*S-AI11)*AMW2-2D0*AI12**2)       
     &       +AI12/2D0*(AI11**2+AI12**2))*(XS4WT-XS4W0)
     &     +(-AMW2*(2D0*S-AI11)+(AI11**2+AI12**2+S*AI12)/2D0
     &       +S*AMT2/2D0)*AMT2*XS4WT)
*
       IF(IMTADD.EQ.0) THEN
        XBOX=XBOX+XBBX*V_TB**2
       ELSE
        XBOX=XBOX
       ENDIF 
*
      ENDIF
*
      END
 
      SUBROUTINE ROKAPN(INL,T,QFER,DROV,ROFAC,AKFAC)
C
C INL=1 FOR ELECTRON NEUTRINO
C INL=2 FOR MUON     NEUTRINO
C T=Q^2
C QFER=-1   FOR NEUTRINO-ELECTRON SCATTERING
C QFER=-1/3 FOR NEUTRINO-DOWN     SCATTERING
C QFER= 2/3 FOR NEUTRINO-UP       SCATTERING
C ROFAC, AKFAC ARE CALCULATED EWFF RHO AND KAPPA
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      COMMON/CDZFLG/IHVP,IAMT4,IQCD,IMOMS,IMASS,IALEM,IMASK,IBARB,IFTJR
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
      COMMON/CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
      COMMON/CDZZWG/AMZ,AMH,GMU,A0,GAMZ,GAMW,CALSZ,CALST,CALXI,CALQED
      COMMON/CDZWSC/SL2,SQ2,W0,W0F,Z0,Z0F,DWZ0R1,DWZ0F,XWM1,XWM1F,XZM1,
     &      XZM1F,XWZ1R1,XDWZ1F,XZFM1,XZFM1F,XAMM1,XAMM1F,XWFM1,XWFM1F
      CALL FORMFN(INL,T,QFER,DDROV,RO1,AK1)
      DROV=AL4PI*DDROV
      AMT =AMQ(5)
      AMT2=AMT*AMT
      RMTW=AMT2/AMW2
      RTT=T/AMT2
      Q2M=ABS(T)
      ALSZ=CALSZ
      ALST=CALXI
      SW2=R1
      CW2=R
      AMT2=AMQ(5)**2
*
* Mixed QCD-corrections
*
      XZERO=DCMPLX(0.D0,0.D0)
      XRQCD=XZERO
      XKQCD=XZERO
*
      IF    (IQCD.EQ.1) THEN
        XRQCD=AL4PI*XRQCDS(ALSZ,ALST,AMZ2,AMW2,AMT2,Q2M)
        XKQCD=AL4PI*XKQCDS(ALST,AMZ2,AMW2,AMT2,Q2M)
      ELSEIF(IQCD.EQ.2) THEN
        XRQCD=AL4PI*XROQCD(ALSZ,ALST,AMZ2,AMW2,AMT2,Q2M)
        XKQCD=AL4PI*XKAQCD(ALST,AMZ2,AMW2,AMT2,Q2M)
      ELSEIF(IQCD.EQ.3) THEN
* not yet tested (25/02/1998)
        XRQCD=AL1PI*ALST/PI*XRMQCD(AMZ2,AMW2,AMT2,Q2M)
     &       -AL1PI*ALSZ/PI/8D0/SW2/CW2*(VT2+VB2+2D0)
     &                   *Q2M/(AMZ2-Q2M)*LOG(Q2M/AMZ2)
* light quarks not added (25/02/1998), 
* there is a problem for low Q2M, see mixed_QCD.frm
* also IM part is not yet added
        XKQCD=AL1PI*ALST/PI*XKMQCD(AMZ2,AMW2,AMT2,Q2M)
      ENDIF
*
      IF(RTT-.1D0)4,4,5
4     ROFACT=AL4PI/R1*3.D0*RMTW*(1.D0/4.D0-5.D0/12.D0*RTT
     *     +19.D0/120.D0*RTT*RTT)
      GO TO 10
5     ROFACT=AL4PI/R1*3.D0*RMTW*(.5D0*DREAL(XI0(AMW2,T,AMT2,AMT2))
     +     -DREAL(XI1(AMW2,T,AMT2,0.D0)))
10    CONTINUE
      ROFAC=1.D0+RO1*AL4PI/R1+ROFACT+DREAL(XRQCD)
      AKFAC=1.D0+AK1*AL4PI/R1       +DREAL(XKQCD)
C--------------------------------------------------------------------
      IF(IBARB.EQ.0.OR.IBARB.EQ.-1) THEN
       AMT4C=19-2D0*PI2
        ELSEIF(IBARB.EQ.1) THEN
       RBTH=AMT2/AMH2
       ALRB=LOG(RBTH)
       AMT4C=49D0/4D0+PI2+27D0/2D0*ALRB+3D0/2D0*ALRB**2
     &      +RBTH/3D0*(2D0-12D0*PI2+12D0*ALRB-27D0*ALRB**2)
     &  +RBTH**2/48D0*(1613-240*PI2-1500*ALRB-720 *ALRB**2)
        ELSEIF(IBARB.EQ.2) THEN
       RBARB=SQRT(AMH2/AMT2)
       AMT4C=FBARB(RBARB)
      ENDIF
C--------------------------------------------------------------------
      IF (IAMT4 .EQ. 1 ) THEN
       SW2 = R1
       AMT2  = AMQ(5)**2
       DRHOT = .75D0*AL4PI/SW2/R*AMT2/AMZ2
       TOPX2 = GMU*AMT2/DSQRT(2.D0)/8.D0/PI2
       DRHOT4= 3.D0*TOPX2*(1.D0+TOPX2*AMT4C)
       ROFAC =(ROFAC-DRHOT)/(1.D0-DRHOT4)
       AKFAC =(AKFAC-R/SW2*DRHOT)*(1.D0+R/SW2*DRHOT4)
      ELSEIF(IAMT4 .EQ. 2 ) THEN
       SW2 = R1
       AMT2  = AMQ(5)**2
       DRHOT = .75D0*AL4PI/SW2/R*AMT2/AMZ2
       TOPX2 = GMU*AMT2/DSQRT(2.D0)/8.D0/PI2
       DRHOT4= 3.D0*TOPX2*(1.D0+TOPX2*AMT4C)
       IF(IQCD.EQ.0) THEN
        TBQCD0= 0D0
        TBQCDL= 0D0
         ELSE
        TBQCD0= TOPX2*ALST/PI*2D0*(1D0+PI2/3D0)
        TBQCDL= AL4PI*ALST/PI*AMT2/AMW2/R1*(.5D0+PI2/6D0)
       ENDIF
       ROFAC =(ROFAC-DRHOT-TBQCDL)/(1.D0-DRHOT4-TBQCD0)
       AKFAC =(AKFAC-R/SW2*(DRHOT+TBQCDL))*(1D0+R/SW2*(DRHOT4+TBQCD0))
      ELSEIF(IAMT4 .GE. 3 ) THEN
       DWZ1AL=R/R1*DREAL(XWZ1R1+XDWZ1F)
       RENORM=SQRT(2D0)*GMU*AMZ2*R1*R/PI*ALFAI
       SCALE = AL4PI/R1*(41D0/6D0-11D0/3D0*R)*ALR
       CORKAP=(AL4PI*DWZ1AL+SCALE)+.75D0*AL4PI/SW2**2*AMT2/AMZ2
*
       SW2 = R1
       AMT2  = AMQ(5)**2
       DRHOT = .75D0*AL4PI/SW2/R*AMT2/AMZ2
       TOPX2 = GMU*AMT2/DSQRT(2.D0)/8.D0/PI2
       DRHOT4= 3.D0*TOPX2*(1.D0+TOPX2*AMT4C)
       IF(IQCD.EQ.0) THEN
        TBQCD0=0D0
        TBQCDL=0D0
         ELSE
        TBQCD0=TOPX2*ALST/PI*2D0*(1D0+PI2/3D0)
        TBQCDL=AL4PI*ALST/PI*AMT2/AMW2/R1*(.5D0+PI2/6D0)
       ENDIF
       ROFAC=(ROFAC-DRHOT-TBQCDL)/(1.D0-DRHOT4-TBQCD0)
       AKFAC=(AKFAC-R/SW2*(DRHOT+TBQCDL)+CORKAP)
     &      *(1D0+R/SW2*(DRHOT4+TBQCD0)-CORKAP*RENORM)
      ENDIF
C-------------------------------------------------------------------
      END
 
      SUBROUTINE FORMFN(INL,Q2,QFER,DDROV,RO1,AK1)
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      COMMON/CDZFLG/IHVP,IAMT4,IQCD,IMOMS,IMASS,IALEM,IMASK,IBARB,IFTJR
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
      COMMON/CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
      COMMON/CDZZWG/AMZ,AMH,GMU,A0,GAMZ,GAMW,CALSZ,CALST,CALXI,CALQED
      COMMON/CDZWSC/SL2,SQ2,W0,W0F,Z0,Z0F,DWZ0R1,DWZ0F,XWM1,XWM1F,XZM1,
     &      XZM1F,XWZ1R1,XDWZ1F,XZFM1,XZFM1F,XAMM1,XAMM1F,XWFM1,XWFM1F
      COMMON/CDZFNE/RNU,ZGM
      AQFER=ABS(QFER)
      QFER2=QFER*QFER
      SI=QFER/AQFER
      RO1=3.D0/4.D0*ALR/R1+9.D0/4.D0
     *   +3.D0/4.D0*RZ*(ALRZ/R/RZ1-ALRW/(R-RZ))
     *   -3.D0/2.D0*(1.D0+SI)
     *   -3.D0/2.D0/R*SI*(1.D0/2.D0-2.D0*R1*AQFER+4.D0*R12*QFER2)
      IND=2*INL
      V1B1W=4.D0*DREAL(XI3(AMW2,Q2,AML(IND)**2,AML(IND)**2))
      RNU=AL4PI/R1*V1B1W
      CHQ21=DREAL(XI3(AMW2,Q2,AML(2)**2,AML(2)**2)
     *     +XI3(AMW2,Q2,AML(4)**2,AML(4)**2)
     *     +XI3(AMW2,Q2,AML(6)**2,AML(6)**2))
      CHMQ1=CHQ21
      DO 201 I=1,6
      AMQ2=AMQ(I)**2
      CHMQ1=CHMQ1+3.D0*CQM(I)   *DREAL(XI3(AMW2,Q2,AMQ2,AMQ2))
      CHQ21=CHQ21+3.D0*CQM(I)**2*DREAL(XI3(AMW2,Q2,AMQ2,AMQ2))
201   CONTINUE
      XWZ1AL=R*XWZ1R1+R*XDWZ1F
      DWZ1AL=DREAL(XWZ1AL)
      ZGM=AL4PI/R1*(+8.D0*R1*CHQ21-2.D0*CHMQ1)
      AK1=-DWZ1AL+3.D0/2.D0*(1.D0+SI)-5.D0-2.D0/3.D0*R
     *    +3.D0/2.D0/R*SI*(1.D0/2.D0-3.D0*R1*AQFER+4.D0*R12*QFER2)+V1B1W
     *    +8.D0*R1*CHQ21-2.D0*CHMQ1
      DDROV=-DWZ1AL/R
      END
 
      FUNCTION XV1B(Q2,AMV2)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
*
      AL=Q2/AMV2
      ALA=LOG(ABS(AL))
      REV1B=-3.5D0+2.D0/AL+(3.D0-2.D0/AL)*ALA
     *     +2.D0*(1.D0-1.D0/AL)**2*(SPENCE(1.D0-AL)-F1)
      AIV1B=0.D0
      IF(Q2.LT.0.D0)
     &AIV1B=PI*(-3.D0+2.D0/AL+2.D0*(1.D0-1.D0/AL)**2*LOG(1.D0-AL))
      XV1B=DCMPLX(REV1B,AIV1B)
      END
 
      FUNCTION XA1B(Q2,AMV2)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      RMV=AMV2/Q2
      XA1B=-8.D0*RMV+32.D0/3.D0
     *    +(2.D0*RMV-17.D0/6.D0)/Q2*XL(Q2,AMV2,AMV2)
      END
 
      FUNCTION XV2B(Q2,AMV2)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      RMV=AMV2/Q2
      XV2B=5.D0/2.D0-8.D0/3.D0+2.D0*RMV
     *    +(5.D0/12.D0+3.D0/4.D0-RMV)/Q2*XL(Q2,AMV2,AMV2)
     *    +2.D0*(-2.D0*RMV+RMV**2)*XJ3(-Q2,AMV2)
*
      END
 
      FUNCTION XDZB(Q2)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
*
      XDRZH=XDL(Q2,-AMZ2,AMZ2,AMH2)
      XDRWW=XDL(Q2,-AMZ2,AMW2,AMW2)
      AQ=AMZ2/Q2
      XZH=-1.D0/12.D0*RZ12*AQ+1.D0/24.D0
     *   *(1.D0+RZ1*(-10.D0+5.D0*RZ-RZ2)*AQ+RZ1*RZ12*AQ**2)*ALRZ
     *   +(-11.D0+4.D0*RZ-RZ2+RZ12*AQ)/24.D0/Q2*XL(Q2,AMZ2,AMH2)
     *   +(1.D0/2.D0-RZ/6.D0+RZ2/24.D0)*XDRZH
      XZH1=-1.D0/12.D0*RZ12*AQ+1.D0/24.D0
     *    *(1.D0+RZ1*(-10.D0+5.D0*RZ-RZ2)*AQ+RZ1*RZ12*AQ**2)*ALRZ
     *    +(-11.D0+4.D0*RZ-RZ2+RZ12*AQ)/24.D0/Q2*XL(Q2,AMZ2,AMH2)
      XZH2=(1.D0/2.D0-RZ/6.D0+RZ2/24.D0)*XDRZH
      XZL=2.D0*R2/Q2*XL(Q2,AMW2,AMW2)
     *   +(-2.D0*R2-17.D0/6.D0*R+2.D0/3.D0+1.D0/24.D0/R)*XDRWW
      XDZB=34.D0/3.D0*R-35.D0/18.D0-4.D0/9.D0/R-ALR/R/12.D0+XZL+XZH/R
*
      END
 
      FUNCTION XDL(Q2,Q2SBT,AM12,AM22)
*
* XDL(Q2,AMQ2,AM12,AM22)=(L(Q2,AM12,AM22)-L(Q2SBT,AM12,AM22))/(Q2-Q2SBT)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      DATA EPS/1.D-3/
*
      Q2S=Q2SBT+AM12+AM22
      ALAM=Q2S**2-4.D0*AM12*AM22
      DSLAM=DSQRT(DABS(ALAM))
      QD=Q2-Q2SBT
      RQD=DABS(QD/DSLAM)
      IF(RQD.LE.EPS) GO TO 1
      XDL=(XL(Q2,AM12,AM22)-XL(Q2SBT,AM12,AM22))/QD
      RETURN
1     R=4.D0*AM12/AM22
      IF(R-1.D0)2,3,2
2     XJS=XJ(Q2SBT,AM12,AM22)
      XDL=2.D0+Q2S*XJS+QD/ALAM*(Q2S-2.D0*AM12*AM22*XJS)
     *   +(QD/ALAM)**2*(-Q2S**2/3.D0-8.D0/3.D0*AM12*AM22*Q2S*XJS)
      RETURN
3     CONTINUE
      RAT=QD/AM22
      XDL=4.D0+2.D0/3.D0*RAT-2.D0/15.D0*RAT*RAT
*
      END
 
      FUNCTION XJ(Q2,AM12,AM22)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
*
* XJ(Q2,M12,M22)=J(Q2,M12,M22)
*
      ALAM=(Q2+AM12+AM22)**2-4.D0*AM12*AM22
      REJ=FJJ(Q2,AM12,AM22)
      AIJ=0.D0
      TRES=(SQRT(AM12)+SQRT(AM22))**2
      IF(-Q2.LE.TRES) GO TO 1
      SLAM=SQRT(ALAM)
      AIJ=2.D0*PI/SLAM
1     XJ=DCMPLX(REJ,AIJ)
*
      END
 
      FUNCTION XL(Q2,AM12,AM22)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
*
* XL(Q2,M12,M22)=L(Q2,M12,M22)=ALAM(Q2,-M12,-M22)*J(Q2,M12,M22)
*
      ALAM=(Q2+AM12+AM22)**2-4.D0*AM12*AM22
      REL=ALAM*FJJ(Q2,AM12,AM22)
      AIL=0.D0
      TRES=(SQRT(AM12)+SQRT(AM22))**2
      IF(-Q2.LE.TRES.OR.ALAM.LT.0D0) GO TO 1
      SLAM=SQRT(ALAM)
      AIL=2.D0*PI*SLAM
1     XL=DCMPLX(REL,AIL)
*
      END
 
      FUNCTION FJJ(Q2,AM12,AM22)
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
*
      Q2M=Q2+AM12+AM22
      Q2L=4.D0*AM12*AM22
      ALAM=Q2M*Q2M-Q2L
      IF(ALAM)1,5,6
1     SLAM=SQRT(-ALAM)
      IF(Q2M)2,3,4
2     CONTINUE
      R1=2.D0/SLAM*ATAN(SLAM/Q2M)
      FJJ=R1+2.D0*PI/SLAM
      RETURN
3     FJJ=0.5D0*PI/SQRT(AM12*AM22)
      RETURN
4     CONTINUE
      R1=2.D0/SLAM*ATAN(SLAM/Q2M)
      FJJ=R1
      RETURN
5     FJJ=2.D0/Q2M
      RETURN
6     SLAM=SQRT(ALAM)
      IF(Q2M)7,8,8
7     FJJ=LOG(Q2L/(Q2M-SLAM)**2)/SLAM
      RETURN
8     FJJ=LOG((Q2M+SLAM)**2/Q2L)/SLAM
*
      END
 
      FUNCTION XI0(AMW2,Q2,AM12,AM22)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
*
      DATA EPS/1.D-4/
*
***** XI0(MW2,Q2,M12,M22)=I0(Q2,M12,M22)
*
      AL1W=LOG(AM12/AMW2)
      IF(AM22/AM12.LT.EPS) GO TO 1
      AL12=LOG(AM12/AM22)
      DM12=(AM12-AM22)/Q2
      XI0=AL1W-2.D0-(1.D0+DM12)/2.D0*AL12+XL(Q2,AM12,AM22)/2.D0/Q2
      RETURN
1     AQ=AM12/Q2
      RELQ=LOG(ABS(1.D0+1.D0/AQ))
      AILQ=0.D0
      TRES=(SQRT(AM12)+SQRT(AM22))**2
      IF(-Q2.GT.TRES) AILQ=-PI
      XLQ=DCMPLX(RELQ,AILQ)
      XI0=AL1W-2.D0+(1.D0+AQ)*XLQ
*
      END
 
      FUNCTION XI1(AMW2,Q2,AM12,AM22)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
*
      DATA EPS/1.D-4/
*
***** XI1(MW2,Q2,M12,M22)=I1(Q2,M12,M22)
*
      AL1W=LOG(AM12/AMW2)
      IF(AM22/AM12.LT.EPS) GO TO 1
      AL12=LOG(AM12/AM22)
      DM12=(AM12-AM22)/Q2
      XI1=AL1W/2D0-1.D0-DM12/2.D0
     *   -(1.D0+2.D0*AM12/Q2+DM12**2)/4.D0*AL12
     *   +(1.D0+DM12)/4.D0*XL(Q2,AM12,AM22)/Q2
      RETURN
1     AQ=AM12/Q2
      RELQ=LOG(ABS(1.D0+1.D0/AQ))
      AILQ=0.D0
      TRES=(SQRT(AM12)+SQRT(AM22))**2
      IF(-Q2.GT.TRES) AILQ=-PI
      XLQ=DCMPLX(RELQ,AILQ)
      XI1=AL1W/2D0-1.D0-AQ/2.D0+(1.D0+AQ)**2/2.D0*XLQ
*
      END
 
      FUNCTION XI3(AMW2,Q2,AM12,AM22)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
*
      DATA EPS/1.D-4/
*
***** XI3(MW2,Q2,M12,M22)=I3(Q2,M12,M22)=INT(Y*(1-Y)*LN...)
*
      AL1W=LOG(AM12/AMW2)
      IF(AM22/AM12.LT.EPS) GO TO 1
      AL12=LOG(AM12/AM22)
      DM12=(AM12-AM22)/Q2
      SM12=(AM12+AM22)/Q2
      XI3=AL1W/6D0-5.D0/18.D0+SM12/3.D0+DM12**2/3.D0
     *   +(-0.5D0+3.D0/2.D0*SM12*DM12+DM12**3)/6.D0*AL12
     *   +(0.5D0-SM12/2.D0-DM12**2)/6.D0*XL(Q2,AM12,AM22)/Q2
      RETURN
1     AQ=AM12/Q2
      RELQ=LOG(ABS(1.D0+1.D0/AQ))
      AILQ=0.D0
      TRES=(SQRT(AM12)+SQRT(AM22))**2
      IF(-Q2.GT.TRES) AILQ=-PI
      XLQ=DCMPLX(RELQ,AILQ)
      XI3=AL1W/6D0-5.D0/18.D0+AQ/3.D0+AQ**2/3.D0
     *   +(1.D0-2.D0*AQ)*(1.D0+AQ)**2/6.D0*XLQ
*
      END
 
      FUNCTION XDI0(Q2,AMV2,AM12,AM22)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      RAT=DABS(AM12/AMV2)
      AL12=LOG(AM12/AM22)
      XDI0=(XL(Q2,AM12,AM22)-(AM12-AM22)*AL12)/2.D0/Q2
     *    -XDL(Q2,-AMV2,AM12,AM22)/2.D0
*
      END
 
      FUNCTION XDI1(Q2,AMV2,AM12,AM22)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      DATA EPS/1.D-4/
*
      IF(AM22.LT.1.D-10) GO TO 1
      AL12=LOG(AM12/AM22)
      DM12=(AM12-AM22)/Q2
      AQ=Q2/AMV2
      SMV1=1.D0-AQ
      XDI1=-DM12/2.D0-(2.D0*AM12/Q2+DM12**2*SMV1)*AL12/4.D0
     *    +(1.D0+DM12*SMV1)/4.D0/Q2*XL(Q2,AM12,AM22)
     *    -(1.D0-(AM12-AM22)/AMV2)/4.D0*XDL(Q2,-AMV2,AM12,AM22)
      RETURN
* CHAIN1 WILL BE USED ONLY FOR W-WIDTH
1     QV=(Q2+AMV2)/AMV2
      IF(ABS(QV).LT.EPS) GO TO 2
      XDI1=(XI1(AMW2,Q2,AM12,AM22)-XI1(AMW2,-AMV2,AM12,AM22))/QV
      RETURN
2     R1V=AM12/AMV2
      VQ=AMV2/Q2
      AQ=AM12/Q2
      RDI1=0.5D0*(-1.D0+R1V*(1.D0-VQ)-0.5D0*QV
     *    +(2.D0*AQ+VQ*(VQ-1.D0)*R1V**2)*LOG(ABS(1.D0+1.D0/AQ)))
      XDI1=DCMPLX(RDI1,0.D0)
*
      END
 
      FUNCTION XDI3(Q2,AMV2,AM12,AM22)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      DATA EPS/1.D-4/
*
      IF(AM22.LT.1.D-10) GO TO 1
      AL12=LOG(AM12/AM22)
      DM12=(AM12-AM22)/Q2
      SM12=(AM12+AM22)/Q2
      AQ=Q2/AMV2
      SMV1=1.D0-AQ
      SMV2=1.D0-AQ+AQ*AQ
      XDI3=SM12/3.D0+DM12**2/3.D0*SMV1
     *    +(SM12*DM12/4.D0*SMV1+DM12**3/6.D0*SMV2)*AL12
     *    +(0.5D0-SM12/2.D0*SMV1-DM12**2*SMV2)/6.D0/Q2*XL(Q2,AM12,AM22)
     *    -(0.5D0+0.5D0*(AM12+AM22)/AMV2
     *    -((AM12-AM22)/AMV2)**2)/6.D0*XDL(Q2,-AMV2,AM12,AM22)
      RETURN
* CHAIN1 WILL BE USED ONLY FOR W-WIDTH
1     QV=(Q2+AMV2)/AMV2
      IF(ABS(QV).LT.EPS) GO TO 2
      XDI3=(XI3(AMW2,Q2,AM12,AM22)-XI3(AMW2,-AMV2,AM12,AM22))/QV
      RETURN
2     R1V=AM12/AMV2
      VQ=AMV2/Q2
      VQ2=1.D0-VQ+VQ**2
      AQ=AM12/Q2
      RDI3=-1.D0/6.D0+(VQ-0.5D0)/3.D0*R1V+VQ2/3.D0*R1V**2
     *    -QV*(0.5D0+R1V)/6.D0-VQ*R1V**2*((VQ-1.D0)/2.D0
     *    +VQ2/3.D0*R1V)*LOG(ABS(1.D0+1.D0/AQ))
      XDI3=DCMPLX(RDI3,0.D0)
*
      END
 
      FUNCTION XDWF(Q2)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      COMMON/CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
      COMMON/CDZWSC/SL2,SQ2,W0,W0F,Z0,Z0F,DWZ0R1,DWZ0F,XWM1,XWM1F,XZM1,
     &      XZM1F,XWZ1R1,XDWZ1F,XZFM1,XZFM1F,XAMM1,XAMM1F,XWFM1,XWFM1F
*
      NG=3
      XSI3  =DCMPLX(0D0,0D0)
      XSDI1U=DCMPLX(0D0,0D0)
      XSDI1D=DCMPLX(0D0,0D0)
      DO 1 I=1,NG
      AML2=AML(2*I)**2
      XI31=XI3 (AMW2,Q2,AML2,0D0)
      XI32=XDI3(Q2,AMW2,AML2,0D0)
      XSI3=XSI3+XI3(AMW2,Q2,AML2,0D0)-XDI3(Q2,AMW2,AML2,0D0)
      XSDI1D=XSDI1D+AML2/AMW2*XDI1(Q2,AMW2,AML2,0D0)
1     CONTINUE
      DO 2 I=1,NG
      AMU2=AMQ(2*I-1)**2
      AMD2=AMQ(2*I  )**2
      XSI3=XSI3+3D0*(XI3(AMW2,Q2,AMU2,AMD2)-XDI3(Q2,AMW2,AMU2,AMD2))
      XSDI1U=XSDI1U+3D0*AMU2/AMW2*XDI1(Q2,AMW2,AMU2,AMD2)
      XSDI1D=XSDI1D+3D0*AMD2/AMW2*XDI1(Q2,AMW2,AMD2,AMU2)
2     CONTINUE
      XDWF=2D0*XSI3+XSDI1U+XSDI1D
*
      END
 
      FUNCTION XDZF(Q2)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      COMMON/CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
      DATA EPS/1.D-3/
*
      NG=3
      ALQ=LOG(ABS(Q2/AMZ2))
      AIL1=0.D0
      AIL2=0.D0
      IF(Q2.LT.0.D0) AIL1=-PI
      IF(Q2.GT.0.D0) AIL2=-PI
      XL1=DCMPLX(ALQ,AIL1)
      XL2=DCMPLX(ALQ,AIL2)
      QD=AMZ2+Q2
      RQD=QD/AMZ2
      XLQ=XL1+1.D0+RQD/2.D0+RQD*RQD/3.D0
      IF(DABS(RQD).GT.EPS) XLQ=XL1-AMZ2/(AMZ2+Q2)*XL2
      XSNU=NG/R/6.D0*(-5.D0/3.D0-ALR+XLQ)
      XSI3=DCMPLX(0.D0,0.D0)
      XSQMI3=DCMPLX(0.D0,0.D0)
      XSQ2I3=DCMPLX(0.D0,0.D0)
      XSDI0=DCMPLX(0.D0,0.D0)
      DO 1 I=1,NG
      AML2=AML(2*I)**2
      XSI3=XSI3+XI3(AMW2,Q2,AML2,AML2)-XDI3(Q2,AMZ2,AML2,AML2)
      XSDI0=XSDI0+AML2/AMW2*XDI0(Q2,AMZ2,AML2,AML2)
1     CONTINUE
      XSQMI3=XSI3
      XSQ2I3=XSI3
      INQ=2*NG
      DO 2 I=1,INQ
      AMQ2=AMQ(I)**2
      XSI3=XSI3+3.D0*(XI3(AMW2,Q2,AMQ2,AMQ2)-XDI3(Q2,AMZ2,AMQ2,AMQ2))
      XSQMI3=XSQMI3+3.D0*CQM(I)*(XI3(AMW2,Q2,AMQ2,AMQ2)
     *      -XDI3(Q2,AMZ2,AMQ2,AMQ2))
      XSQ2I3=XSQ2I3+3.D0*CQM(I)**2*(XI3(AMW2,Q2,AMQ2,AMQ2)
     *      -XDI3(Q2,AMZ2,AMQ2,AMQ2))
      XSDI0=XSDI0+3.D0*AMQ2/AMW2*XDI0(Q2,AMZ2,AMQ2,AMQ2)
2     CONTINUE
      XDZF=(8.D0*R12*XSQ2I3-4.D0*R1*XSQMI3+XSI3)/R+XSDI0/2.D0+XSNU
*
      END
 
      FUNCTION XAMF(Q2)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      COMMON/CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
*
      NG=3
      XCHMQ1=DCMPLX(0.D0,0.D0)
      DO 1 I=1,NG
      AML2=AML(2*I)**2
      XCHMQ1=XCHMQ1+XI3(AMW2,Q2,AML2,AML2)
1     CONTINUE
      XCHQ21=XCHMQ1
* equal for leptons
      INQ=2*NG
      DO 2 I=1,INQ
      AMQ2=AMQ(I)**2
      XCHMQ1=XCHMQ1+3.D0*CQM(I)*XI3(AMW2,Q2,AMQ2,AMQ2)
      XCHQ21=XCHQ21+3.D0*CQM(I)**2*XI3(AMW2,Q2,AMQ2,AMQ2)
2     CONTINUE
      XAMF=8.D0*R1*XCHQ21-2.D0*XCHMQ1
*
      END
 
      FUNCTION XAMFEF(Q2,SEFF2)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      COMMON/CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
*
      NG=3
      XCHMQ1=DCMPLX(0.D0,0.D0)
      DO 1 I=1,NG
      AML2=AML(2*I)**2
      XCHMQ1=XCHMQ1+XI3(AMW2,Q2,AML2,AML2)
1     CONTINUE
      XCHQ21=XCHMQ1
      INQ=2*NG
      DO 2 I=1,INQ
      AMQ2=AMQ(I)**2
      XCHMQ1=XCHMQ1+3.D0*CQM(I)*XI3(AMW2,Q2,AMQ2,AMQ2)
      XCHQ21=XCHQ21+3.D0*CQM(I)**2*XI3(AMW2,Q2,AMQ2,AMQ2)
2     CONTINUE
      XAMFEF=8.D0*SEFF2*XCHQ21-2.D0*XCHMQ1
      END
 
      FUNCTION XFOTF(Q2)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      COMMON/CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
*
      NG=3
      XCHQ21=DCMPLX(0.D0,0.D0)
      DO 1 I=1,NG
      AML2=AML(2*I)**2
      XCHQ21=XCHQ21+XI3(AML2,Q2,AML2,AML2)
1     CONTINUE
      INQ=2*NG
      DO 2 I=1,INQ
      AMQ2=AMQ(I)**2
      XCHQ21=XCHQ21+3.D0*CQM(I)**2*XI3(AMQ2,Q2,AMQ2,AMQ2)
2     CONTINUE
      XFOTF=8.D0*XCHQ21
      END
 
      FUNCTION XFOTF1(IALEM,IHVP,IQCD,ITOP,Q2)
*
* NEG.Q2 IS S-CHANNEL, POS. T-CHANNEL
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      COMMON/CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
      COMMON/CDZZWG/AMZ,AMH,GMU,A0,GAMZ,GAMW,CALSZ,CALST,CALXI,CALQED
      COMMON/CDZ513/DAL5H
*
      NG=3
      XCHQ21=DCMPLX(0D0,0D0)
C********* LEPTONIC PART OF VACUUM POLARISATION ********************
      DO 1 I=1,NG
      AML2=AML(2*I)**2
      XCHQLL=XI3(AML2,Q2,AML2,AML2)
      XCHQ21=XCHQ21+XCHQLL
 1    CONTINUE
*
* Leptonic one-loop contribution via DALPHL
*
      DALEP=DALPHL(1,ABS(Q2))
      XCHQ21=DCMPLX(DALEP/AL4PI/8,DIMAG(XCHQ21))
*      
      AMT2=AMQ(5)**2
      ALSZ=CALSZ
      ALST=CALXI
      SW2=R1
* AL4PI is removed from ALFQCD as Tord suggested
      IF(IQCD-1) 2,3,4
 2    ALFQCD=0D0
      GO TO 5
 3    ALFQCD=ALQCDS(ALST,AMZ2,AMW2,AMT2)
      GO TO 5
 4    ALFQCD=ALQCD (ALST,AMZ2,AMW2,AMT2)
 5    CONTINUE
C********* HADRONIC PART OF VACUUM POLARISATION ********************
      IF (IHVP - 2) 6,7,9
 6    CONTINUE
      XCHQ25=3D0*CQM(5)**2*XI3(AMT2,Q2,AMT2,AMT2)
      XCHQ21=XCHQ21+XCHQ25*ITOP
*
      AMQ2=ABS(Q2)
C***** JEGERLEHNER/EIDELMAN
      IF(MOD(IALEM,2).EQ.0) THEN
        UDCSB=DALH5(-Q2,AMZ)-DALH5(AMZ2,AMZ)+DAL5H
      ELSE
        UDCSB=DALH5(-Q2,AMZ)     
      ENDIF
C*****
      INQ=2*NG
      XCHQQQ=DCMPLX(0.D0,0.D0)
      DO 61 I=1,INQ
      IF(ITOP.EQ.0.AND.I.EQ.5) GOTO 61
      AMQ2=AMQ(I)**2
      XCHQQQ=XCHQQQ+3D0*CQM(I)**2*XI3(AMQ2,Q2,AMQ2,AMQ2)
 61   CONTINUE
      XFOTF1=8*XCHQ21+DCMPLX(UDCSB/AL4PI,8*DIMAG(XCHQQQ))+ALFQCD*ITOP
*
* this is a correction to Im part of \alpha_em, which improved
* the agreement with TOPAZ0
*     XFOTF1=8*XCHQ21
*    &     +DCMPLX(UDCSB/AL4PI,8*DIMAG(XCHQQQ*(1+CALSZ/PI)))+ALFQCD*ITOP
      RETURN
C****** THIS HADRONIC PART WAS USED IN THE 1989 Z PHYSICS WORKSHOP
 7    INQ=2*NG
      XCHQQQ=DCMPLX(0.D0,0.D0)
      DO 8 I=1,INQ
      IF(ITOP.EQ.0.AND.I.EQ.5) GOTO8
      AMQ2=AMQ(I)**2
      XCHQQQ=XCHQQQ+3D0*CQM(I)**2*XI3(AMQ2,Q2,AMQ2,AMQ2)
 8    CONTINUE
      XCHQ21=XCHQ21+XCHQQQ
      XFOTF1=8D0*XCHQ21+ALFQCD*ITOP
      RETURN
C****** THIS IS BURKHARDT ET AL. HADRONIC VACUUM POLARIZATION
 9    CONTINUE
      XCHQ25=+3.D0*CQM(5)**2*XI3(AMT2,Q2,AMT2,AMT2)
      XCHQ21=XCHQ21+XCHQ25*ITOP
      XUDSCB = XADRQQ(-Q2)/AL4PI
      XFOTF1 = 8D0*XCHQ21 + XUDSCB +ALFQCD*ITOP
*
      END
 
      FUNCTION XFOTF3(IALEM,IALE2,IHVP,IQCD,ITOP,DAL5H,Q2)
*
* Contrary to XFOTF1, this function may return the leptonic part up to
* three-loop corrections, which is governed by the flag IALE2 as desribed. 
*
* NEG.Q2 IS S-CHANNEL, POS. T-CHANNEL
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      COMMON/CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
      COMMON/CDZZWG/AMZ,AMH,GMU,A0,GAMZ,GAMW,CALSZ,CALST,CALXI,CALQED
* For print controll
*      DATA IPRC/0/
*      IPRC=IPRC+1
*
      NG=3
      XCHQ21=DCMPLX(0D0,0D0)
C********* LEPTONIC PART OF VACUUM POLARISATION ********************
      DO 1 I=1,NG
      AML2=AML(2*I)**2
      XCHQLL=XI3(AML2,Q2,AML2,AML2)
      XCHQ21=XCHQ21+XCHQLL
 1    CONTINUE
*
* Leptonic one-loop contribution via DALPHL
*
      DALEP=DALPHL(IALE2,ABS(Q2))
*
* Control print of DALEP
*     DALEP1=DALPHL(1,ABS(Q2))
*     DALEP2=DALPHL(2,ABS(Q2))-DALPHL(1,ABS(Q2))
*     DALEP3=DALPHL(3,ABS(Q2))-DALPHL(2,ABS(Q2))
*     DALEP =DALPHL(3,ABS(Q2))
*     print *,'MZ,SQ =',AMZ,SQRT(ABS(Q2)),DALEP
*     stop 
      XCHQ21=DCMPLX(DALEP/AL4PI/8,DIMAG(XCHQ21))
*
      AMT2=AMQ(5)**2
      ALSZ=CALSZ
      ALST=CALXI
      SW2=R1
* AL4PI is removed from ALFQCD as Tord suggested
      IF(IQCD-1) 2,3,4
 2    ALFQCD=0D0
      GO TO 5
 3    ALFQCD=ALQCDS(ALST,AMZ2,AMW2,AMT2)
      GO TO 5
 4    CONTINUE
      IF(IQCD.EQ.2) ALFQCD=ALQCD (ALST,AMZ2,AMW2,AMT2)
      IF(IQCD.EQ.3) ALFQCD=AQCDBK(ALST,AMZ2,AMW2,AMT2)
 5    CONTINUE
*
* Control print of the QCDcontribution for the book
*     if(iprc.lt.10) print *,'ALST,QCD=',ALST,AL4PI*ALFQCD
C********* HADRONIC PART OF VACUUM POLARISATION ********************
      SELECT CASE (IHVP)
      CASE (1, 4:)
         XCHQ25=3D0*CQM(5)**2*XI3(AMT2,Q2,AMT2,AMT2)
         XCHQ21=XCHQ21+XCHQ25*ITOP
*
* Control print of the top contribution for the book
*      if(iprc.lt.10) print *,'top=',AL4PI*8*DREAL(XCHQ25)
*
         AMQ2=ABS(Q2)
C***** JEGERLEHNER/EIDELMAN
         IF(MOD(IALEM,2).EQ.0) THEN
            UDCSB=DALH5(-Q2,AMZ)-DALH5(AMZ2,AMZ)+DAL5H
         ELSE
            UDCSB=DALH5(-Q2,AMZ)
* Control print of UDCSB
*     print *,'UDCSB=',UDCSB,Q2,AMZ
         ENDIF
C*****
         INQ=2*NG
         XCHQQQ=DCMPLX(0.D0,0.D0)
         DO 61 I=1,INQ
         IF(ITOP.EQ.0.AND.I.EQ.5) GOTO 61
         AMQ2=AMQ(I)**2
         XCHQQQ=XCHQQQ+3D0*CQM(I)**2*XI3(AMQ2,Q2,AMQ2,AMQ2)
 61      CONTINUE
         XFOTF3=8*XCHQ21+DCMPLX(UDCSB/AL4PI,8*DIMAG(XCHQQQ))+ALFQCD*ITOP
* Control print
*     print *,'Before RETURN, XFOTF3=',XFOTF3,ALFQCD
*     stop
*
* this is a correction to Im part of \alpha_em, which improved
* the agreement with TOPAZ0
*     XFOTF3=8*XCHQ21
*    &     +DCMPLX(UDCSB/AL4PI,8*DIMAG(XCHQQQ*(1+CALSZ/PI)))+ALFQCD*ITOP
      CASE (2)
C****** THIS HADRONIC PART WAS USED IN THE 1989 Z PHYSICS WORKSHOP
         INQ=2*NG
         XCHQQQ=DCMPLX(0.D0,0.D0)
         DO 8 I=1,INQ
         IF(ITOP.EQ.0.AND.I.EQ.5) GOTO8
         AMQ2=AMQ(I)**2
         XCHQQQ=XCHQQQ+3D0*CQM(I)**2*XI3(AMQ2,Q2,AMQ2,AMQ2)
 8       CONTINUE
         XCHQ21=XCHQ21+XCHQQQ
         XFOTF3=8D0*XCHQ21+ALFQCD*ITOP
      CASE (3)
C****** THIS IS BURKHARDT ET AL. HADRONIC VACUUM POLARIZATION
         XCHQ25=+3.D0*CQM(5)**2*XI3(AMT2,Q2,AMT2,AMT2)
         XCHQ21=XCHQ21+XCHQ25*ITOP
         XUDSCB = XADRQQ(-Q2)/AL4PI
         XFOTF3 = 8D0*XCHQ21 + XUDSCB +ALFQCD*ITOP
      END SELECT
      RETURN
*
      END
 
      SUBROUTINE SETCON(AMW,AMZ,AMT,AMH,DR,DRBIG,DRREM)
*
************************************************************************
*   THIS SUBROUTINE                                                    *
*   CALCULATES GAUGE BOSON MASSES, DELTA.R, AND THE WEAK MIXING ANGLE  *
*   OF THE ELECTROWEAK THEORY FOR GIVEN INPUT VALUES, DEPENDING ON IMOMS
*   ( ALL QUANTITIES SHOULD BE GIVEN IN GEV )                          *
************************************************************************
*
      IMPLICIT REAL*8(A-H,O-X)
      COMPLEX*16 XFOTF3
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
      COMMON/CDZZWG/FMZ,FMH,GMU,A0,GAMZ,GAMW,CALSZ,CALST,CALXI,CALQED
      COMMON/CDZFLG/IHVP,IAMT4,IQCD,IMOMS,IMASS,IALEM,IMASK,IBARB,IFTJR
*
      COMMON/CDZSCT/ISCRE,ISCAL,IAFMT,IFACR,IFACT,IHIGS,IEWLC,ICZAK
     &             ,IHIG2,IALE2,IGFER      
      COMMON/CDZ513/CDAL5H      
      COMMON/CDZDMW/IDMWW,IDSWW
*
      IF(IMOMS.EQ.1) THEN
       AMZ2=AMZ**2
       AMW2=AMZ2*(1D0+SQRT(1D0-4D0*A0**2/AMZ2))/2D0
       DRP=0D0
       DO 1 ITER=1,20
       CALL CONST2(3,AMW2,AMZ2)
       CALL
     &    SEARCH(IHVP,IAMT4,IQCD,IMASS,IALEM,IBARB,AAFAC,DR,DRBIG,DRREM)
       DRF=DR
       IF(ABS(DRF-DRP).LE.1D-11) GO TO 11
       DRP=DRF
       AMW2=AMZ2*(1D0+SQRT(1D0-4D0*AAFAC**2/AMZ2))/2D0
 1     CONTINUE
 11    CONTINUE
       CALL CONST2(3,AMW2,AMZ2)
       SW2=1D0-AMW2/AMZ2
       IF(IAMT4.EQ.5) THEN
*       PRINT *,'DEGRASSIs M_W IS REPLACED BY COMPLETE FERMIONIC 2-LOOP FROM '
*       PRINT *,'FREITAS ET AL.'
*       print *,'AMW-OLD=',SQRT(AMW2)
*--
*- M_w block:
*--
        DAL5HH=CDAL5H
        DALPHA=AL4PI*DREAL(XFOTF3(IALEM,IALE2,IHVP,1,0,DAL5HH,-AMZ2))
        DELALP=DALPHA/0.05924D0-1D0
        DELALS=CALSZ/0.119D0-1D0
        DELHIG=LOG(AMH/100D0)
        DELTOP=(AMT/174.3D0)**2-1D0
        DELAMZ=AMZ-91.1875D0
*
* Improved expansion including MZ and MT**4 dependence
* result for fermionic 2-loop corrections as of June 2001:
*
        AMW0=80.3768D0
        C1  = 0.05619D0
        C2  = 1.078D0
        C3  = 0.5236D0
        C4  = 0.0765D0
        C5  = 0.009305D0
        C6  = 0.0005365D0
        C7  = 0.00544D0
        C8  = 1.261D0
        C9  = 0.0727D0
*
        AMWNEW=AMW0-C1*DELHIG-C5*DELHIG**2+C6*DELHIG**4
     &             +C3*DELTOP-C9*DELTOP**2-C7*DELHIG*DELTOP
     &             -C2*DELALP-C4*DELALS+C8*DELAMZ
*
* Theoretical uncertainty from the authors
*
*        IDMWW=0 ! +-1 should be new external flag with values 0,+-1 [sigma]
        E0  = 0.00566D0
        E1  = 0.000315D0
        E2  = 0.000328D0
        AMWNEW=AMWNEW+IDMWW*(E0+E1*DELHIG+E2*DELHIG**2)
*       print *,'AMW-NEW=',AMWNEW
        AMW2=AMWNEW**2
        CALL CONST2(3,AMW2,AMZ2)
        SW2=1D0-AMW2/AMZ2
       ELSEIF(IAMT4.GE.6) THEN
*       PRINT *,'DEGRASSIs M_W IS REPLACED BY COMPLETE 2-LOOP FROM '
*       PRINT *,'AWRAMIK ET AL.'
*--
*- M_w block:
*--
        DAL5HH=CDAL5H
        DALPHA=AL4PI*DREAL(XFOTF3(IALEM,IALE2,IHVP,1,0,DAL5HH,-AMZ2))
        DELALP=DALPHA/0.05907D0-1D0
        DELALS=CALSZ/0.119D0-1D0
        DELHIG=LOG(AMH/100D0)
	DELHSQ=(AMH/100D0)**2
        DELTOP=(AMT/174.3D0)**2-1D0
        DELAMZ=(AMZ/91.1875D0)-1D0
*
* Improved expansion including MZ, MT**4 and MH**2 dependence
* result for complete 2-loop corrections as of Nov. 2003:
*
        AMW0=80.3799D0
        C1  = 0.05429D0
        C2  = 0.008939D0
	C3  = 0.0000890D0
	C4  = 0.000161D0
	C5  = 1.070D0
	C6  = 0.5256D0
	C7  = 0.0678D0
	C8  = 0.00179D0
	C9  = 0.0000659D0
	C10 = 0.0737D0
	C11 = 114.9D0
*
        AMWNEW=AMW0-C1*DELHIG-C2*DELHIG**2+C3*DELHIG**4+C4*(DELHSQ-1D0)
* 06 12 2008: next line was too long in v.6.42
     &        +C6*DELTOP-C7*DELTOP**2-C8*DELHIG*DELTOP+C9*DELHSQ*DELTOP
     &             -C5*DELALP-C10*DELALS+C11*DELAMZ
*
* Theoretical uncertainty from the authors
*
*        IDMWW=0 ! +-1 should be new external flag with values 0,+-1 [sigma]
        E0  = 0.004D0
        AMWNEW=AMWNEW+IDMWW*E0
*       print *,'AMW-NEW=',AMWNEW
        AMW2=AMWNEW**2
        CALL CONST2(3,AMW2,AMZ2)
        SW2=1D0-AMW2/AMZ2
       ENDIF
      ELSEIF(IMOMS.EQ.2) THEN
       AMW2=AMW**2
       AMZ2=AMW2/(1D0-A0**2/AMW2)
       DO 2 ITER=1,20
       CALL CONST2(3,AMW2,AMZ2)
       CALL
     &    SEARCH(IHVP,IAMT4,IQCD,IMASS,IALEM,IBARB,AAFAC,DR,DRBIG,DRREM)
       DRF=DR
       IF(ABS(DRF-DRP).LE.1D-11) GO TO 21
       DRP=DRF
       AMZ2=AMW2/(1D0-AAFAC**2/AMW2)
 2     CONTINUE
 21    CONTINUE
       CALL CONST2(3,AMW2,AMZ2)
       SW2=1D0-AMW2/AMZ2
      ELSE
       AMW2=AMW**2
       AMZ2=AMZ**2
       CALL CONST2(3,AMW2,AMZ2)
       CALL
     &    SEARCH(IHVP,IAMT4,IQCD,IMASS,IALEM,IBARB,AAFAC,DR,DRBIG,DRREM)
       SW2=1D0-AMW2/AMZ2
      ENDIF
      END
 
      SUBROUTINE CONST1(IHVP,FPMT,FPMH)
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      COMMON/CDZSCT/ISCRE,ISCAL,IAFMT,IFACR,IFACT,IHIGS,IEWLC,ICZAK
     &             ,IHIG2,IALE2,IGFER      
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
      COMMON/CDZZWG/AMZ,AMH,GMU,A0,GAMZ,GAMW,CALSZ,CALST,CALXI,CALQED
      COMMON/CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
      COMMON/CDZTHR/AMTH(6)
      DIMENSION AMHF(4),CLMI(8),AMLI(8),AMTI(6)
      DATA CLMI/0.D0,1.D0,0.D0,1.D0,0.D0,1.D0,0.D0,1.D0/
      DATA AMLI/0D0,.51099907D-3,0D0,.105658389D0,0D0,1.77705D0,2*0D0/
      DATA AMHF/  1.D0,5.D1,  2.D2,  7.D1/
*                  HNU HLE    TP     BP
      DATA AMTI/2*.134974D0,1.548465D0,.493646D0,0D0,4.73016D0/
*
* NUMERICAL CONSTANTS
      PI=ATAN(1D0)*4D0
      ALFAI=137.035 989 5 D0
      D3=1.2020569031596D0
      PI2=PI**2
      F1=PI2/6D0
      AL1PI=1D0/PI/ALFAI
      AL2PI=AL1PI/2D0
      AL4PI=AL1PI/4D0
* WS-PARAMETERS
      AMH=FPMH
      IF(IGFER.LE.1) THEN
        GMU=1.16639D-5
      ELSE
        GMU=1.16637D-5
      ENDIF
      A0=SQRT(PI/ALFAI/SQRT(2D0)/GMU)
*
C     FERMION PARAMETERS (SEE ALSO DATA)
*
      DO 2 I2=1,6
      CLM(I2) =CLMI(I2)
      AML(I2) =AMLI(I2)
      AMTH(I2)=AMTI(I2)
2     CONTINUE
      AML(7)=AMHF(1)
      AML(8)=AMHF(2)
      DO 1 I=1,4
      CQM(2*I-1)=2.D0/3.D0
1     CQM(2*I)=1.D0/3.D0
      SELECT CASE (IHVP)
         CASE (1)
C IHVP=1 USES THIS SET TOGETHER WITH THE JEGERLEHNER FIT WITH KNIEHL
            AMQ(1)=.062D0
            AMQ(2)=.083D0
            AMQ(3)=1.50D0
            AMQ(4)=.215D0
            AMQ(5)=FPMT
            AMQ(6)=4.70D0
            AMQ(7)=AMHF(3)
            AMQ(8)=AMHF(4)
         CASE (2)
C FOR YB. MASSES ARE INFLUENTIAL
            AMQ(1)=.04145D0
            AMQ(2)=.04146D0
            AMQ(3)=1.50D0
            AMQ(4)=0.15D0
            AMQ(5)=FPMT
            AMQ(6)=4.70D0
            AMQ(7)=AMHF(3)
            AMQ(8)=AMHF(4)
         CASE (3)
C USED WITH BURKHARDT'S ROUTINE HADRQQ, XADRQQ
            AMQ(1)=.04145D0
            AMQ(2)=.04146D0
            AMQ(3)=1.50D0
            AMQ(4)=0.15D0
            AMQ(5)=FPMT
            AMQ(6)=4.70D0
            AMQ(7)=AMHF(3)
            AMQ(8)=AMHF(4)
         CASE (4 : )
C USED WITH 2016,2017 FITS OF JEGERLEHNER
            AMQ(1)=.062D0
            AMQ(2)=.083D0
            AMQ(3)=1.666D0
            AMQ(4)=0.15D0
            AMQ(5)=FPMT
            AMQ(6)=4.800D0
            AMQ(7)=AMHF(3)
            AMQ(8)=AMHF(4)
      END SELECT
*
      END
 
      SUBROUTINE CONST2 (NG,FPMW2,FPMZ2)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      COMMON/CDZFLG/IHVP,IAMT4,IQCD,IMOMS,IMASS,IALEM,IMASK,IBARB,IFTJR
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
      COMMON/CDZZWG/AMZ,AMH,GMU,A0,GAMZ,GAMW,CALSZ,CALST,CALXI,CALQED
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      COMMON/CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
      COMMON/CDZWSC/SL2,SQ2,W0,W0F,Z0,Z0F,DWZ0R1,DWZ0F,XWM1,XWM1F,XZM1,
     &      XZM1F,XWZ1R1,XDWZ1F,XZFM1,XZFM1F,XAMM1,XAMM1F,XWFM1,XWFM1F
*
* FILL CDZWSM
*
      AMW2=FPMW2
      AMZ2=FPMZ2
      AMH2=AMH**2
      AKSX=AMH/AMZ
      CW2M=AMW2/AMZ2
      SW2M=1.D0-CW2M
      R=AMW2/AMZ2
      R1=1.D0-R
* VECTOR COUPLINGS FOR QCDCOR
      VB=1D0-4D0*CQM(6)*R1
      VT=1D0-4D0*CQM(5)*R1
      VB2=VB**2
      VB2T=1D0-8D0*CQM(6)*R1
      VT2=VT**2
      VT2T=1D0-8D0*CQM(5)*R1
* VECTOR COUPLINGS FOR QCDCOR
      R12=R1**2
      R2=R**2
      R1W=1.D0+R
      R1W2=R1W*R1W
      RW=AMH2/AMW2
      RW1=1.D0-RW
      RW12=RW1**2
      RW2=RW**2
      RZ=AMH2/AMZ2
      RZ1=1.D0-RZ
      RZ12=RZ1**2
      RZ2=RZ**2
      ALR=LOG(R)
      ALRW=LOG(RW)
      ALRZ=LOG(RZ)
*
* FILL FERMIONIC PARTS
*
      SL2=0D0
      SQ2=0D0
      NQ=2*NG
      W0F=0D0
      Z0F=0D0
      XWM1F=(0D0,0D0)
      XZM1F=(0D0,0D0)
*     T'HOOFT SCALE NOW APPLIED ONLY FOR FERMIONS, CANCELS IN DELTA_R
cbardin  Why 'ONLY FOR FERMIONS' but not for BOSONS??? This I forgot.
      AMW2MU=AMW2
*
      DO 1 I=1,NG
*     NEUTRINO PART
      ALRNOR=LOG(AMW2MU/AMZ2)
      XALR=DCMPLX(ALRNOR,PI)
      XZM1F=XZM1F+1D0/R/6D0*(XALR+5D0/3D0)

      AL=LOG(AMW2MU/AML(2*I)**2)
      SL2=SL2+CLM(2*I)**2*AL
      AL=LOG(AMW2MU/AMQ(2*I)**2)
      SQ2=SQ2+3D0*CQM(2*I)**2*AL
      AL=LOG(AMW2MU/AMQ(2*I-1)**2)
      SQ2=SQ2+3D0*CQM(2*I-1)**2*AL
*
      AML2=AML(2*I  )**2
      AMT2=AMQ(2*I-1)**2
      AMB2=AMQ(2*I  )**2
*
      RMLW=AML2/AMW2
      RMTW=AMT2/AMW2
      RMBW=AMB2/AMW2
      ALLW=LOG(AML2/AMW2MU)
      ALTW=LOG(AMT2/AMW2MU)
      ALBW=LOG(AMB2/AMW2MU)
*
      W0F=W0F+1D0/2D0*(RMLW*ALLW-RMLW/2D0)
      Z0F=Z0F+1D0/2D0*RMLW*ALLW
      IF(RMTW.NE.RMBW) THEN
       W0F=W0F+3D0/2D0*((RMTW**2*ALTW-RMBW**2*ALBW)/(RMTW-RMBW)
     *    -(RMTW+RMBW)/2D0)
      ELSE
       W0F=W0F+3D0*(RMTW*ALTW-RMTW/2D0)
      ENDIF
      Z0F=Z0F+3D0/2D0*(RMTW*ALTW+RMBW*ALBW)
*
      XWM1F=XWM1F-2D0*XI3(AMW2MU,-AMW2,AML2,0D0 )
     *     +     RMLW*XI1(AMW2MU,-AMW2,AML2,0D0 )
      XWM1F=XWM1F-6D0*XI3(AMW2MU,-AMW2,AMT2,AMB2)
     *     +3D0*(RMTW*XI1(AMW2MU,-AMW2,AMT2,AMB2)
     *     +     RMBW*XI1(AMW2MU,-AMW2,AMB2,AMT2))

*
      V2PA2L=1D0+(1D0-4D0*R1*CLM(2*I  ))**2
      XZM1F=XZM1F-1D0/2D0*V2PA2L/R*XI3(AMW2MU,-AMZ2,AML2,AML2)
     *     +          1D0/2D0*RMLW*XI0(AMW2MU,-AMZ2,AML2,AML2)
      V2PA2T=1D0+(1D0-4D0*R1*CQM(2*I-1))**2
      V2PA2B=1D0+(1D0-4D0*R1*CQM(2*I  ))**2
      XZM1F=XZM1F-3D0/2D0*V2PA2T/R*XI3(AMW2MU,-AMZ2,AMT2,AMT2)
     *     +          3D0/2D0*RMTW*XI0(AMW2MU,-AMZ2,AMT2,AMT2)
      XZM1F=XZM1F-3D0/2D0*V2PA2B/R*XI3(AMW2MU,-AMZ2,AMB2,AMB2)
     *     +          3D0/2D0*RMBW*XI0(AMW2MU,-AMZ2,AMB2,AMB2)
 1    CONTINUE
*
*     DERIVATIVES, USED ONLY IN FORMFACTORS AND PARTIAL WIDTHS
*
      XWFM1F=XDWF(-AMW2)
      XZFM1F=XDZF(-AMZ2)
      XAMM1F=XAMF(-AMZ2)
*
* Gambino's modification of XZM1F, include mixing ZA squared in the mass 
*                                  counterterm (Dyson resummation)
      IF(IAMT4.GE.4) XZM1F=XZM1F-AL4PI*XAMM1F**2
*
      DWZ0F =( W0F - Z0F )/R1
      XDWZ1F=(XWM1F-XZM1F)/R1
*
*     FILL BOSONIC PARTS
*
      XL1=XL(-AMW2,AMH2,AMW2)/AMW2
      XJ1=XJ(-AMW2,AMH2,AMW2)*AMH2
      XL2=XL(-AMW2,AMW2,AMZ2)/AMW2
      XL3=XL(-AMZ2,AMH2,AMZ2)/AMW2
      XJ3=XJ(-AMZ2,AMH2,AMZ2)*AMH2/R
      XL4=XL(-AMZ2,AMW2,AMW2)/AMW2
      R3=R2*R
      W0=5.D0/8.D0/R-17.D0/4.D0+5.D0/8.D0*R*(1.D0+R)-RW/8.D0
     *  +3.D0/4.D0*RW/RW1*ALRW+(3.D0/4.D0/R+9.D0/4.D0-3.D0/R1)*ALR
      Z0=5.D0/8.D0/R-RW/8.D0+3.D0/4.D0/R*ALR+3.D0/4.D0*RW/RZ1*ALRZ
      XWM1=1.D0/12.D0/R2+23.D0/12.D0/R-157.D0/9.D0-RW/2.D0+RW2/12.D0
     *    -RW*(3.D0/4.D0-RW/4.D0+RW2/24.D0)*ALRW
     *    +(1.D0/24.D0/R3+7.D0/12.D0/R2-7.D0/2.D0/R)*ALR
     *    +(0.5D0-RW/6.D0+RW2/24.D0)*XL1
     *    +(1.D0/24.D0/R2+2.D0/3.D0/R-17.D0/6.D0-2.D0*R)*XL2
      XZM1=35.D0/18.D0/R+35.D0/18.D0-34.D0/3.D0*R-8.D0*R2-RW/2.D0
     *    +RW2*R/12.D0+RW*(-3.D0/4.D0+RZ/4.D0-RZ2/24.D0)*ALRZ
     *    +5.D0/6.D0/R*ALR+(0.5D0-RZ/6.D0+RZ2/24.D0)*XL3
     *    +(1.D0/24.D0+2.D0/3.D0*R-17.D0/6.D0*R2-2.D0*R3)*XL4
      DWZ0R1=(W0-Z0)/R1
      XWZ1R1=(XWM1-XZM1)/R1
      XZFM1=-4.D0*R2+17.D0/3.D0*R-23.D0/9.D0+5.D0/18.D0/R-RW/2.D0
     *     +RW*RZ/6.D0-ALR/12.D0/R
     *     +RW*(-3.D0/4.D0+3.D0/8.D0*RZ-RZ2/12.D0)*ALRZ+0.5D0/R*ALRZ
     *     +(-R*R2+7.D0/6.D0*R2-17.D0/12.D0*R-1.D0/8.D0)*XL4
     *     +(0.5D0-5.D0/24.D0*RZ+1.D0/12.D0*RZ2)*XL3+0.5D0*XJ3
      XAMM1=2.D0/9.D0/R+35.D0/18.D0-34.D0/3.D0*R-8.D0*R2
     *     +(1.D0/24.D0+2.D0/3.D0*R-17.D0/6.D0*R2-2.D0*R*R2)*XL4
      XWFM1=R-34.D0/9.D0+2.D0/R+1.D0/6.D0/R2-RW/2.D0+RW**2/6.D0
     *     +(3.D0*R+5.D0/2.D0-17.D0/4.D0/R+7.D0/8.D0/R2+1.D0/12.D0/R3)
     *     *ALR+(0.5D0-3.D0*RW/4.D0+3.D0*RW2/8.D0-RW**3/12.D0)*ALRW
     *     +(-R/2.D0-2.D0+25.D0/24.D0/R+1.D0/12.D0/R2)*XL2
     *     +(0.5D0-5.D0*RW/24.D0+RW2/12.D0)*XL1+0.5D0*XJ1
*
      END
 
      SUBROUTINE
     &    SEARCH(IHVP,IAMT4,IQCD,IMASS,IALEM,IBARB,AAFAC,DR,DRBIG,DRREM)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
      COMMON/CDZSCT/ISCRE,ISCAL,IAFMT,IFACR,IFACT,IHIGS,IEWLC,ICZAK
     &             ,IHIG2,IALE2,IGFER      
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
      COMMON/CDZZWG/AMZ,AMH,GMU,A0,GAMZ,GAMW,CALSZ,CALST,CALXI,CALQED
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      COMMON/CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
      COMMON/CDZWSC/SL2,SQ2,W0,W0F,Z0,Z0F,DWZ0R1,DWZ0F,XWM1,XWM1F,XZM1,
     &      XZM1F,XWZ1R1,XDWZ1F,XZFM1,XZFM1F,XAMM1,XAMM1F,XWFM1,XWFM1F
      COMMON/CDZDEG/DROBAR,DROBLO,DRREMD
      COMMON/CDZ513/DAL5H
      COMMON/CDZAPV/DRHOD,DKDREM,DROREM
*
      W0AL=W0+W0F
      WM1AL=DREAL(XWM1+XWM1F)
      DWZ1F =R/R1*DREAL(XDWZ1F)
      DWZ1B =R/R1*DREAL(XWZ1R1)
      DWZ1AL=R/R1*DREAL(XWZ1R1+XDWZ1F)
      RXX=-2D0/3D0+4D0/3D0*(SL2+SQ2)+DWZ1AL
     *   +(W0AL-WM1AL-5D0/8D0*R2-5D0/8D0*R+11D0/2D0+9D0/4D0*R/R1*ALR)/R1
      RXXFER=4D0/3D0*(SL2+SQ2)+R/R1*DREAL(XDWZ1F)+(W0F-DREAL(XWM1F))/R1
      RXXBOS=-2D0/3D0+R/R1*DREAL(XWZ1R1)+(W0
     &   -DREAL(XWM1)-5D0/8D0*R2-5D0/8D0*R+11D0/2D0+9D0/4D0*R/R1*ALR)/R1
      AMH=SQRT(AMH2)
      AMW=SQRT(AMW2)
*
      CLQQCD=0D0
      XTBQCD=DCMPLX(0D0,0D0)
      ALFQCD=0D0
      AMT2=AMQ(5)**2
      SW2=R1
*
* Mixed QCD-corrections to \dr
*
      IF(IQCD.EQ.0) GO TO 4
      ALSZ=CALSZ
      ALST=CALST
*
* Mixed QCD-corrections from two light doublets
*
      CLQQCD=-AL4PI*(R/R1-1D0)/R1*ALR*ALSZ/PI*(1D0+1.409*ALSZ/PI
     &       -12.805*(ALSZ/PI)**2)
*
* Mixed QCD-corrections from t-b doublet
*
      IF    (IQCD.EQ.0) THEN
        XTBQCD=DCMPLX(0D0,0D0)
        ALFQCD=0D0
      ELSEIF(IQCD.EQ.1) THEN
        XTBQCD=AL4PI*DCMPLX(RXQCDS(ALST,AMZ2,AMW2,AMT2),0D0)
        XTBQCD=AL4PI*DCMPLX(RXQCDS(ALST,AMZ2,AMW2,AMT2),0D0)
        ALFQCD=AL4PI*ALQCDS(ALST,AMZ2,AMW2,AMT2)
      ELSEIF(IQCD.EQ.2) THEN
        XTBQCD=AL4PI*DCMPLX( RXQCD(ALST,AMZ2,AMW2,AMT2),0D0)
        ALFQCD=AL4PI*ALQCD (ALST,AMZ2,AMW2,AMT2)
      ELSEIF(IQCD.EQ.3) THEN
        XTBQCD=AL1PI*ALST/PI*DRMQCD(AMZ2,AMW2,AMT2)
        ALFQCD=0D0
* for IQCD=3, ALFQCD is a part of \delta_r - reminder  
      ENDIF
*
* Hadronic vacuum polarization
*
 4    SELECT CASE(IHVP)
      CASE (1,4:)
         XQQ15=(0D0,0D0)
         DO 6 IQ=1,6
         AMQ2 = AMQ(IQ)*AMQ(IQ)
         IF(IQ.EQ.5) GO TO 6
         XQQ15=XQQ15 + 6D0*XI3(AMQ2,-AMZ2,AMQ2,AMQ2) * 3D0 *CQM(IQ)**2
 6       CONTINUE
C***** JEGERLEHNER/EIDELMAN
         IF(MOD(IALEM,2).EQ.0) THEN
           UDCSB=DAL5H
         ELSE
           UDCSB=DALH5(AMZ2,AMZ)
         ENDIF
C*****
         IF(IMASS.EQ.0) THEN
          DR1FER=AL4PI*RXXFER-AL4PI*4D0/3D0*DREAL(XQQ15)+UDCSB
          DR1BOS=AL4PI*RXXBOS
         ELSE
          DR1FER=AL4PI*RXXFER
          DR1BOS=AL4PI*RXXBOS
         ENDIF
      CASE (2)
         DR1FER=AL4PI*RXXFER
         DR1BOS=AL4PI*RXXBOS
      CASE (3)
         XQQ15=(0D0,0D0)
         DO 9 IQ=1,6
         AMQ2 = AMQ(IQ)*AMQ(IQ)
         IF(IQ.EQ.5) GO TO 9
         XQQ15=XQQ15 + 6D0*XI3(AMQ2,-AMZ2,AMQ2,AMQ2) * 3D0 *CQM(IQ)**2
 9       CONTINUE
         XUDSCB = XADRQQ(AMZ2)
*IN XADRQQ FROM H.BURKHARDT: NEG. ARGUMENT=T-CHANL, POSIT. ARG.=S-CHANNE
         IF(IMASS.EQ.0) THEN
          DR1FER=DREAL(AL4PI*RXXFER - AL4PI*4D0/3D0*XQQ15 + XUDSCB)
          DR1BOS=      AL4PI*RXXBOS
         ELSE
          DR1FER=AL4PI*RXXFER
          DR1BOS=AL4PI*RXXBOS
         ENDIF
      END SELECT
*
      TBQCDL=0D0
      TBQCD0=0D0
      TBQCD3=0D0
      CORRDR=0D0
      IF(IQCD.NE.0) THEN
       IF(IAFMT.EQ.0) THEN
* Below is the leading term as it is given in \zf (3.73)
        TBQCD0=-CALXI/PI*2D0/3D0*(1D0+PI2/3D0)
        TBQCD3=0D0
        CORRDR=AL4PI/R1**2*AMT2/AMZ2*
     &  (-ALST/PI*(.5D0+PI2/6D0)-.75D0*TBQCD0)
         ELSE
* This coincides with TBQCD0 above if O(\alpha\alpha^2_s)=0
        TBQCD0=AFMT3(CALST,AMT2,AMZ2,SW2)
* Below is pure AFMT-correction to \Delta r (see p.18 PCWG-notebook)
        TBQCD3=-.75D0*AL4PI/R1**2*AMT2/AMZ2
     &          *(TBQCD0-(-ALST/PI*2D0/3D0*(1D0+PI2/3D0)))
        CORRDR=0D0
       ENDIF
       TBQCD = DREAL(XTBQCD)
* Below is -c^2_W/s^2_W \Delta \rho
       TBQCDL= AL4PI*ALST/PI*AMT2/AMW2*R/R1**2*(.5D0+PI2/6D0)
       DR1FER=DR1FER+TBQCD+2D0*CLQQCD+ALFQCD
      ENDIF
*
      DR1=DR1FER+DR1BOS
      DR =DR1
      DRREM = 0D0
C--------------------------------------------------------------------
      IF(IBARB.EQ.0.OR.IBARB.EQ.-1) THEN
       AMT4C=19-2D0*PI2
        ELSEIF(IBARB.EQ.1) THEN
       RBTH=AMT2/AMH2
       ALRB=LOG(RBTH)
       AMT4C=49D0/4D0+PI2+27D0/2D0*ALRB+3D0/2D0*ALRB**2
     &      +RBTH/3D0*(2D0-12D0*PI2+12D0*ALRB-27D0*ALRB**2)
     &  +RBTH**2/48D0*(1613-240*PI2-1500*ALRB-720 *ALRB**2)
        ELSEIF(IBARB.EQ.2) THEN
       RBARB=SQRT(AMH2/AMT2)
       AMT4C=FBARB(RBARB)
      ENDIF
*
* XFOTF3-CALL to compute DALPH's
*
      DALFA1=DREAL(XFOTF3(IALEM,    1,IHVP,IQCD,0,DAL5H,-AMZ2))*AL4PI
      DALFA =DREAL(XFOTF3(IALEM,IALE2,IHVP,IQCD,0,DAL5H,-AMZ2))*AL4PI
*
*-----------------------------------------------------------------------
      IF (IAMT4 .EQ. 1) THEN
*-----------------------------------------------------------------------
       DRHO1 = .75D0*AL4PI/SW2/R*AMT2/AMZ2
       AXF   = AMT2*GMU/(8D0*PI2*SQRT(2D0))
       DRIRR = 3D0*AXF*(1D0+AMT4C*AXF)
* the leading HIGGS-contribution
       DRHIG1= 0d0
       DRHIGS= 0d0
       fachig= log(amh2/amw2)-5d0/6
       if(ihigs.eq.1.and.fachig.gt.0d0) then
        DRHIG1= al4pi/4/sw2                *11d0/3*fachig
        DRHIGS= sqrt(2d0)*gmu*amw2/16/pi**2*11d0/3*fachig
       endif
* See p.18 of PCWG-notebook
       DRREMF= DR1FER+R/SW2*DRHO1-DALFA1+TBQCD3+CORRDR
       DRREMB= DR1BOS-DRHIG1
       DRREM = DRREMF+DRREMB
       DRBIG=1-(1+R/SW2*DRIRR-DRHIGS)*(1D0-DALFA)
* end of HIGGS-modification
       DR=DRBIG+DRREM
       RENORD=SQRT(2D0)*GMU*AMZ2*R1*R*(1D0-DRBIG)/PI*ALFAI
       RENORM=SQRT(2D0)*GMU*AMZ2*R1*R*(1D0-DALFA)/PI*ALFAI
       IF(ISCRE.EQ.0) SCALER=1.00D0
       IF(ISCRE.EQ.1) SCALER=RENORD
       IF(ISCRE.EQ.2) SCALER=RENORM
       DR=DRBIG+SCALER*DRREM
*
*-----------------------------------------------------------------------
      ELSEIF(IAMT4 .EQ. 2) THEN
*-----------------------------------------------------------------------
       DRHO1 = .75D0*AL4PI/SW2/R*AMT2/AMZ2
       AXF   = AMT2*GMU/(8D0*PI2*SQRT(2D0))
       DRIRR = 3D0*AXF*(1D0+AMT4C*AXF+TBQCD0)
* the leading HIGGS-contribution
       DRHIG1= 0d0
       DRHIGS= 0d0
       fachig= log(amh2/amw2)-5d0/6
       if(ihigs.eq.1.and.fachig.gt.0d0) then
        DRHIG1= al4pi/4/sw2                *11d0/3*fachig
        DRHIGS= sqrt(2d0)*gmu*amw2/16/pi**2*11d0/3*fachig
       endif
       DRREMF= DR1FER+R/SW2*DRHO1-DALFA1-TBQCDL
       DRREMB= DR1BOS-DRHIG1
       DRREM = DRREMF+DRREMB
       DRHHS = -.005832*(AL1PI)**2/SW2**2*AMH2/AMW2
       DRREM = DRREM+DRHHS
       DRBIG=1D0-(1D0+R/SW2*DRIRR-DRHIGS)*(1D0-DALFA)
* end of HIGGS-modification
       DR=DRBIG+DRREM
       RENORM=SQRT(2D0)*GMU*AMZ2*R1*R*(1D0-DALFA)/PI*ALFAI
       RENORD=SQRT(2D0)*GMU*AMZ2*R1*R*(1D0-DRBIG)/PI*ALFAI
       IF(ISCRE.EQ.0) SCALER=1.00D0
       IF(ISCRE.EQ.1) SCALER=RENORD
       IF(ISCRE.EQ.2) SCALER=RENORM
       DR=DRBIG+SCALER*DRREM
*
*-----------------------------------------------------------------------
      ELSEIF(IAMT4 .EQ. 3) THEN
*-----------------------------------------------------------------------
       DRHO1 = .75D0*AL4PI/SW2/R*AMT2/AMZ2
       AXF   = AMT2*GMU/(8D0*PI2*SQRT(2D0))
       DRIRR = 3D0*AXF*(1D0+AMT4C*AXF+TBQCD0)
       SCALE = AL4PI/R1*(41D0/6D0-11D0/3D0*R)*ALR
* the leading HIGGS-contribution  (FOR IAMT4=3 11/12 --> 1/12)
       DRHIG1= 0d0
       DRHIGS= 0d0
       fachig= log(amh2/amw2)-5d0/6
       if(ihigs.eq.1.and.fachig.gt.0d0) then
        DRHIG1= al4pi/4/sw2                * 1d0/3*fachig
        DRHIGS= sqrt(2d0)*gmu*amw2/16/pi**2* 1d0/3*fachig
       endif
       DRREM = DR1FER+DR1BOS-DALFA1-TBQCDL-(DWZ1AL*AL4PI+SCALE)-DRHIG1
       DRHHS = -.005832D0*(AL1PI)**2/SW2**2*AMH2/AMW2
       DRREM = DRREM+DRHHS
       RENORM=SQRT(2D0)*GMU*AMZ2*R1*R*(1D0-DALFA)/PI*ALFAI
       DRBIG=1D0-(1D0+R/SW2*DRIRR-DRHIGS)*(1D0-DALFA)
     &   +((AL4PI*DWZ1AL+SCALE)+R/SW2*DRHO1)*RENORM
* end of HIGGS-modification
       DR=DRBIG+SCALER*DRREM
       RENORD=SQRT(2D0)*GMU*AMZ2*R1*R*(1D0-DRBIG)/PI*ALFAI
       IF(ISCRE.EQ.0) SCALER=1.00D0
       IF(ISCRE.EQ.1) SCALER=RENORD
       IF(ISCRE.EQ.2) SCALER=RENORM
       DR=DRBIG+SCALER*DRREM
*
*-----------------------------------------------------------------------
      ELSEIF(IAMT4 .GE. 4) THEN
*-----------------------------------------------------------------------
       DRHO1 = .75D0*AL4PI/SW2/R*AMT2/AMZ2
       AXF   = AMT2*GMU/(8D0*PI2*SQRT(2D0))
       SCALEB = AL4PI/R1*(1D0/6D0+7D0*R)*ALR
       ANUMF  = 24D0
       TRQ2F  =  8D0
       SCALEF = -AL4PI/R1*(ANUMF/6-4D0/3*R1*TRQ2F)*ALR
       SCALE=SCALEB+SCALEF
* the leading HIGGS-contribution  (FOR IAMT4=3 11/12 --> 1/12)
       DRHIG1= 0d0
       DRHIGS= 0d0
       fachig= log(amh2/amw2)-5d0/6
       if(ihigs.eq.1.and.fachig.gt.0d0) then
        DRHIG1= al4pi/4/sw2                *1d0/3*fachig
        DRHIGS= sqrt(2d0)*gmu*amw2/16/pi**2*1d0/3*fachig
       endif
*
       IF(IHIG2.EQ.1) THEN
         DRHHS=-.005832*(AL1PI)**2/SW2**2*AMH2/AMW2
       ELSE
         DRHHS = 0D0
       ENDIF
*
* "Old" reminder DRREM, commented
*
*      DRREM=DR1FER+DR1BOS-DALFA1-TBQCDL-(DWZ1AL*AL4PI+SCALE)
*    &      -DRHIG1+DRHHS
*       print *,'DALFA1=',DALFA1
*       print *,'DALFA =',DALFA 
*       print *,'DRREM =',DRREM
*       print *,'DRLEA =',(DWZ1AL*AL4PI+SCALE)
*
* "New" reminder DRREM from NEWDR \equiv "Old" reminder DRREM
*
       CALL NEWDR(DALFAN,DRLEAN,DRREMN,DRREMK)
       DRREMD=DRREMK
       DRREM =DRREMN+TBQCD+2D0*CLQQCD+ALFQCD-TBQCDL+DRHHS
*
*       print *,'DALFAN=',DALFAN
*       print *,'DRREMN=',DRREM
*       print *,'DRLEAN=',DRLEAN
*       print *,'DRREM =',DRREM
*
       AMT=AMQ(5)
       PI3QF=1D0
*
* PI3QF=|QF|, new parameter of GDEGNL
*
       CALL GDEGNL
     & (GMU,AMZ,AMT,AMH,AMW,PI3QF,AMZ2,DRDREM,DRHOD,DKDREM,DROREM)
*
       DROBAR=3D0*AXF*TBQCD0-(AL4PI*DWZ1AL+SCALE)
     &       *SQRT(2D0)*GMU*AMZ2*R1**2/PI*ALFAI+DRHOD
       DROBLO=              -(AL4PI*DWZ1AL+SCALE)
     &       *SQRT(2D0)*GMU*AMZ2*R1**2/PI*ALFAI
*      
* Activation of old options (SCRE, EXPR) for DR:
*
* New game with SCALER, April 99
*
       RENORM=SQRT(2D0)*GMU*AMZ2*R1*R/PI*ALFAI
       IF(ISCRE.EQ.0) SCALER2=1.00D0
       IF(ISCRE.EQ.1) SCALER2=1D0/RENORM**2
       IF(ISCRE.EQ.2) SCALER2=1D0*RENORM**2
       IF(IFACR.EQ.0) THEN
* DRHHS added to main option IFACR=0
        DR=1D0-(1D0+R/SW2*DROBAR-DRHIGS)*(1D0-DALFA-DRREM 
     &        -(DRHHS+DRDREM)*SCALER2) 
        AAFAC=A0/SQRT(1D0-DR)
       ELSEIF(IFACR.EQ.1) THEN
        DR=1D0-(1D0+R/SW2*DROBAR-DRHIGS)*(1D0-DALFA-DRREM)
     &        +(DRHHS+DRDREM)*SCALER2
        AAFAC=A0/SQRT(1D0-DR)
       ELSEIF(IFACR.EQ.2) THEN
        DR=1D0-(1D0+R/SW2*DROBAR-DRHIGS)*(1D0-DALFA)
     &        +(1D0+R/SW2*DROBLO)*DRREM 
     &        +(DRHHS+DRDREM)*SCALER2
        AAFAC=A0/SQRT(1D0-DR)
       ENDIF
      DRBIG=DR
      ENDIF
*
* DR-options of YR(1994)
*
      IF(IAMT4.LE.3) THEN
       IF(IFACR.EQ.0) THEN
        AAFAC=A0/SQRT(1D0-DR)
       ELSEIF(IFACR.EQ.1) THEN
        AAFAC=A0*SQRT(1/(1-DRBIG)*(1+SCALER*DRREM/(1-DRBIG)))
       ELSEIF(IFACR.EQ.2) THEN
        AAFAC=A0*SQRT((1+SCALER*DRREM)/(1D0-DRBIG))
       ELSEIF(IFACR.EQ.3) THEN
        AAFAC=A0*SQRT(1/(1D0-DRBIG)+SCALER*DRREM)
       ELSE
        STOP
       ENDIF
      ENDIF
*
      END

      SUBROUTINE NEWDR(DALFA,DRLEA,DRREM,DRREMK)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      COMMON/CDZFLG/IHVP,IAMT4,IQCD,IMOMS,IMASS,IALEM,IMASK,IBARB,IFTJR
      COMMON/CDZSCT/ISCRE,ISCAL,IAFMT,IFACR,IFACT,IHIGS,IEWLC,ICZAK
     &             ,IHIG2,IALE2,IGFER      
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
      COMMON/CDZZWG/AMZ,AMH,GMU,A0,GAMZ,GAMW,CALSZ,CALST,CALXI,CALQED
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      COMMON/CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
      COMMON/CDZWSC/SL2,SQ2,W0,W0F,Z0,Z0F,DWZ0R1,DWZ0F,XWM1,XWM1F,XZM1,
     &      XZM1F,XWZ1R1,XDWZ1F,XZFM1,XZFM1F,XAMM1,XAMM1F,XWFM1,XWFM1F
      COMMON/CDZ513/DAL5H
*
* This is done for a cross-check using the results from DB/GP
* It calculates leading, DRLEA, and reminder, DRREM, contributions to \dr
* using exactly final equation (7.405) from DB/GP
*
      ANUMF = 24D0
      TRQ2F =  8D0
*
      XLLA=(0D0,0D0)
      XLL =(0D0,0D0)
      DO 1 IL=1,3
      AML2=AML(2*IL)*AML(2*IL)
      XLLA=XLLA+6D0*XI3(AML2,-AMZ2,AML2,AML2)
      XLL =XLL +6D0*XI3(AMZ2,-AMZ2,AML2,AML2)
 1    CONTINUE
*
      XQQ5=(0D0,0D0)
      DO 2 IQ=1,6
      AMQ2=AMQ(IQ)*AMQ(IQ)
      IF(IQ.EQ.5) GO TO 2
      XQQ5=XQQ5+6D0*XI3(AMZ2,-AMZ2,AMQ2,AMQ2)*3D0*CQM(IQ)**2
 2    CONTINUE
C***** JEGERLEHNER/EIDELMAN
      IF(MOD(IALEM,2).EQ.0) THEN
        UDCSB=DAL5H
      ELSE
        UDCSB=DALH5(AMZ2,AMZ)     
      ENDIF
C*****
*
* NEWDR is unused again, DALFA1 is not modified
*
      DALFA1=UDCSB+AL4PI*4D0/3*DREAL(XLLA)
      DALFA =DALFA1
      DR1FER=DALFA1
     &      +AL1PI/3*(-4D0/3*LOG(AMQ(5)**2/AMZ2)-DREAL(XQQ5+XLL))
     &      +AL4PI/R1*(W0F-DREAL(XWM1F))
     &      +AL1PI/R1*ANUMF/24*LOG(R)
      DR1BOS=AL4PI*(-2D0/3D0
     &      +1D0/R1*(W0-DREAL(XWM1)
     &      -5D0/8D0*R*(1D0+R)+11D0/2D0+9D0/4D0*R/R1*ALR))
     &      -AL4PI/R1*(1D0/6+7D0*R)*LOG(R)
*
      DRREM=DR1FER+DR1BOS-DALFA1
      DRLEA=
     &      +AL4PI*R/R1*DREAL(XDWZ1F)
     &      +AL4PI/R1*(ANUMF/6-4D0/3*R1*TRQ2F)*(-LOG(R))
     &      +AL4PI*R/R1*DREAL(XWZ1R1)
     &      +AL4PI/R1*(1D0/6+7D0*R)*LOG(R)
*
      DRREMK=DR1BOS
* excluded following Guiseppe Degrassi
*     &      +AL1PI/3*(-4D0/3*LOG(AMQ(5)**2/AMZ2))
     &      +AL4PI/R1*(W0F-DREAL(XWM1F))
     &      +AL1PI/R1*ANUMF/24*LOG(R)
*    
      END
 
      FUNCTION XADRQQ(S)
C  HADRONIC IRREDUCIBLE QQ SELF-ENERGY: TRANSVERSE
C     PARAMETRIZE THE REAL PART OF THE PHOTON SELF ENERGY FUNCTION
C     BY  A + B LN(1+C*:S:) , AS IN MY 1981 TASSO NOTE BUT USING
C     UPDATED VALUES, EXTENDED USING RQCD UP TO 100 TEV
C     FOR DETAILS SEE:
C     H.BURKHARDT, F.JEGERLEHNER, G.PENSO AND C.VERZEGNASSI
C     IN CERN YELLOW REPORT ON "POLARIZATION AT LEP" 1988
C     H.BURKHARDT, CERN/ALEPH, AUGUST 1988
C     NEGATIVE VALUES MEAN T - CHANNEL (SPACELIKE)
C     POSITIVE VALUES MEAN S - CHANNEL (TIMELIKE )
C     IN THE SPACE LIKE VALUES AROUND 1 GEV ARE TYPICAL FOR LUMINOSITY
C     THE VALUES AT 92 GEV ( Z MASS ) GIVE THE LIGHT QUARK CONTRIBUTION
C     TO DELTA R
C     TAKE CARE OF THE SIGN OF REPI WHEN USING THIS IN DIFFERENT
C     PROGRAMS
C     HERE REPI WAS CHOSEN TO
C     BE POSITIVE (SO THAT IT CORRESPONDS DIRECTLY TO DELTA ALPHA)
C     OFTEN ITS ASSUMED TO BE NEGATIVE.
C
C     THE IMAGINARY PART IS PROPORTIONAL TO R (HAD / MU CROSS SECTION)
C     AND IS THEREFORE 0 BELOW THRESHOLD ( IN ALL THE SPACELIKE REGION)
C     NOTE ALSO THAT ALPHA_S USUALLY HAS BEEN DERIVED FROM THE MEASURED
C     VALUES OF R.
C     CHANGING ALPHA_S TO VALUES INCOMPATIBLE WITH CURRENT DATA
C     WOULD IMPLY TO BE ALSO INCONSISTENT WITH RE,IM PI
C     DEFINED HERE
C
C     H.BURKHARDT
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 XADRQQ
C
      DATA A1,B1,C1/   0.0   ,  0.00835,  1.0  /
      DATA A2,B2,C2/   0.0   ,  0.00238,  3.927 /
      DATA A3,B3,C3/ 0.00165 ,  0.00300,  1.0  /
      DATA A4,B4,C4/ 0.00221 ,  0.00293,  1.0  /
C
      DATA PI/3.141592653589793D0/,ALFAIN/137.0359895D0/,INIT/0/
C
      IF(INIT.EQ.0) THEN
        INIT=1
        ALFA=1.D0/ALFAIN
        ALFAPI=1.D0/PI/ALFAIN
      ENDIF
      T=ABS(S)
      IF(T.LT.0.3D0**2) THEN
        REPIAA=A1+B1*LOG(1.D0+C1*T)
      ELSEIF(T.LT.3.D0**2) THEN
        REPIAA=A2+B2*LOG(1.D0+C2*T)
      ELSEIF(T.LT.100.D0**2) THEN
        REPIAA=A3+B3*LOG(1.D0+C3*T)
      ELSE
        REPIAA=A4+B4*LOG(1.D0+C4*T)
      ENDIF
C     AS IMAGINARY PART TAKE -I ALFA/3 REXP
      XADRQQ=REPIAA-(0.D0,1.D0)*ALFA/3.*REXP(S)
CEXPO HADRQQ=HADRQQ/(4.D0*PI*ALFA)  ! EXPOSTAR DIVIDES BY 4 PI ALFA
      END
 
      FUNCTION REXP(S)
C  HADRONIC IRREDUCIBLE QQ SELF-ENERGY: IMAGINARY
      IMPLICIT REAL*8(A-H,O-Z)
C     CONTINUUM R = AI+BI W ,  THIS + RESONANCES WAS USED TO CALCULATE
C     THE DISPERSION INTEGRAL. USED IN THE IMAG PART OF HADRQQ
      PARAMETER (NDIM=18)
      DIMENSION WW(NDIM),RR(NDIM),AA(NDIM),BB(NDIM)
      DATA WW/1.,1.5,2.0,2.3,3.73,4.0,4.5,5.0,7.0,8.0,9.,10.55,
     . 12.,50.,100.,1000.,10 000.,100 000./
      DATA RR/0.,2.3,1.5,2.7,2.7,3.6,3.6,4.0,4.0,3.66,3.66,3.66,
     .  4.,3.87,3.84, 3.79, 3.76,    3.75/
      DATA INIT/0/
      IF(INIT.EQ.0) THEN
        INIT=1
C CALCULATE A,B FROM STRAIGHT LINES BETWEEN R MEASUREMENTS
        BB(NDIM)=0.
        DO 4 I=1,NDIM
        IF(I.LT.NDIM) BB(I)=(RR(I)-RR(I+1))/(WW(I)-WW(I+1))
        AA(I)=RR(I)-BB(I)*WW(I)
    4   CONTINUE
       ENDIF
       REXP=0.D0
       IF(S.GT.0.D0) THEN
        W=SQRT(S)
       IF(W.GT.WW(1)) THEN
       DO 2 I=1,NDIM
C      FIND OUT BETWEEN WHICH POINTS OF THE RR ARRAY W IS
       K=I
       IF(I.LT.NDIM) THEN
       IF(W.LT.WW(I+1)) GOTO 3
       ENDIF
    2  CONTINUE
    3  CONTINUE
       REXP=AA(K)+BB(K)*W
C   WRITE(6,'('' K='',I2,'' AA='',F10.2,'' BB='',F10.3)')
C    .   K,AA(K),BB(K)
       ENDIF
      ENDIF
      END
 
      SUBROUTINE ZWRATE(DR,DRBIG,DRREM,QCDCOR,V_TB,PARTZ,PARTW)
*
************************************************************************
*  ZWRATE - Z- AND W- BOSONS DECAY RATES, CALCULATES PARTIAL AND TOTAL *
*  WIDTHS OF Z- AND W- BOSONS WITH ACCOUNT OF ALL 1-LOOP ELECTROWEAK   *
*  AND QED CORRECTIONS (QCD CORRECTIONS ARE ALSO INCLUDED).            *
************************************************************************
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
      COMMON/CDZFLG/IHVP,IAMT4,IQCD,IMOMS,IMASS,IALEM,IMASK,IBARB,IFTJR
      COMMON/CDZSCT/ISCRE,ISCAL,IAFMT,IFACR,IFACT,IHIGS,IEWLC,ICZAK
     &             ,IHIG2,IALE2,IGFER      
      COMMON/CDZZWG/AMZ,AMH,GMU,A0,GAMZ,GAMW,CALSZ,CALST,CALXI,CALQED
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     &      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      COMMON/CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
      COMMON/CDZTHR/AMTH(6)
      COMMON/CDZ513/DAL5H
      COMMON/CDZDMW/IDMWW,IDSWW
*
      COMMON /CDZRKZ/ARROFZ(0:10),ARKAFZ(0:10),ARVEFZ(0:10),ARSEFZ(0:10)
     &              ,AROTFZ(0:10),AIROFZ(0:10),AIKAFZ(0:10),AIVEFZ(0:10)
      COMMON /CDZAUX/PARTZA(0:10),PARTZI(0:10),RENFAC(0:10),SRENFC(0:10)
*
      COMMON/CDZAQF/AQFI(10)
*
      COMMON/CDZDDZ/IDDZZ
*
      DIMENSION MCOLFZ(10),PARTZ(0:11),PARTW(3),QCDCOR(0:14)
      DIMENSION INDF(10),INDL(10),INDQ(10)
      DIMENSION MWFAC(2),AQFW(2)
      DIMENSION ARCZAK(0:6)
*
      DATA MWFAC/3,6/
      DATA INDF /2,1,1,1,4,3,4,3,4,5/
      DATA INDL /1,2,4,6,0,0,0,0,0,0/
      DATA INDQ /0,0,0,0,1,2,3,4,5,6/
      DATA MCOLFZ/1,1,1,1,3,3,3,3,3,3/
*
* MWFAC AND MZFAC - FLAVOUR*COLOR FACTORS FOR W- AND Z- DECAYS
*
      AQFI(1)=0.D0
      AQFI(2)=1.D0
      AQFI(3)=1.D0
      AQFI(4)=1.D0
      AQFI(5)=2.D0/3.D0
      AQFI(6)=1.D0/3.D0
      AQFI(7)=2.D0/3.D0
      AQFI(8)=1.D0/3.D0
      AQFI(9)=2.D0/3.D0
      AQFI(10)=1.D0/3.D0
*
* AQFI - ARRAY OF FINAL PARTICLE CHARGES FOR PARTIAL Z- WIDTHS%
* T,TBAR DECAY CHANNEL IS ASSUMED TO BE ABOVE Z- THRESHOLD AND IS
* NOT ADDED TO THE TOTAL Z- WIDTH
*
      AQFW(1)=1.D0
      AQFW(2)=2.D0/3.D0
*
* Numerical implementation of Czarnecki-Kuehn's corrections
*
      IF(CALSZ.LE.1D-4) THEN
        DO ICZ=0,6
           ARCZAK(ICZ)=0D0
        ENDDO
      ELSE
        ARCZAK(0)= 0.D0
        ARCZAK(1)=-0.113D-3/3
        ARCZAK(2)=-0.160D-3/3
        ARCZAK(3)=-0.113D-3/3
        ARCZAK(4)=-0.160D-3/3
        ARCZAK(5)= 0.D0
        ARCZAK(6)=-0.040D-3/3
      ENDIF
*
* THE SAME FOR PARTIAL W- WIDTHS, HERE ONLY TWO CHANNELS EXIST IF ONE
* NEGLECTS FERMION MASSES. AGAIN T,BBBAR DECAY CHANNEL IS ASSUMED TO BE
*
      GAM0T=0.D0
      GAM0H=0.D0
      GAMWT=0.D0
      GAM1H=0.D0
      GAM1HA=0.D0
      GAM1HI=0.D0
      GAM1T =0.D0
      CONSTZ=GMU*AMZ**3/12.D0/PI/SQRT(2.D0)
*
*     NEW RESULT FOR FERMIONIC 2-LOOP CORRECTIONS TO KAPPA BY AWRAMIK ET AL.
*     Fit formula only computes LEPTONIC eff. weak mixing angle;
*     to compute other flavours we need to take the differences with
*     respect to the leptonic case from the old computation.
*     For this reason the kappa_lept is stored as reference value
*     in AKFACREF
*
      IF((IAMT4.EQ.6).OR.(IAMT4.EQ.8)) THEN
       CALL VERTZW(1,1,1D0)
       CALL
     &  FOKAPP(AQFI(2),DR,DRBIG,DRREM,SCALER,ROFACI,AKFACI,AR1IM,AK1IM)
       AKFACREF = AKFACI
      ENDIF
*
* LOOP ON FERMIONS IN Z DECAY
*
      DO 3 INF=1,10
*
      IBFLA=0
      IF(IAMT4.GE.4.AND.INF.NE.10) THEN
       CALL VERTZW(1,INDF(INF),1D0)
       CALL 
     & FOKAPP(AQFI(INF),DR,DRBIG,DRREM,SCALER,ROFACI,AKFACI,AR1IM,AK1IM)
*
       IF(IAMT4.EQ.6) THEN
*       PRINT *,'DEGRASSIs KAPPA IS REPLACED BY FERMIONIC 2-LOOP FROM '
*       PRINT *,'AWRAMIK ET AL.'
*--
*- Sw_eff block:
*--
        DALPHA=AL4PI*DREAL(XFOTF3(IALEM,IALE2,IHVP,1,0,DAL5H,-AMZ2))
        DELALP=DALPHA/0.05907D0-1D0
        DELALS=CALSZ/0.117D0-1D0
        DELHIG=LOG(AMH/100D0)
	    DELHSP=(AMH/100D0)
        DELTOP=(AMQ(5)/174.3D0)**2-1D0
        DELAMZ=(AMZ/91.1875D0)-1D0
*
* Expansion including MZ, MT**4 and MH**2 dependence
*
        ASW0=0.2313714D0
        C1  = 0.0004734D0
        C2  = 0.0000205D0
        C3  = 3.89D-6
        C4  = -1.83D-6
        C5  = 0.0207D0
        C6  = -0.002749D0
        C7  = 0.000180D0
        C9  = -9.25D-6
        C10 = 0.000387D0
        C11 = -0.656D0
*
        ASWNEW=ASW0+C1*DELHIG+C2*DELHIG**2+C3*DELHIG**4+C4*
     &			(DELHSP**2-1.D0) 
     &             +C6*DELTOP+C7*DELTOP**2+C9*(DELHSP-1D0)*DELTOP 
     &             +C5*DELALP+C10*DELALS+C11*DELAMZ
*
* Theoretical uncertainty from the authors
*
*        IDSWW=0 ! +-1 should be new external flag with values 0,+-1 [sigma]
        E0  = 4.9D-5
        ASWNEW=ASWNEW+IDSWW*E0
	AKFACI=AKFACI-AKFACREF+ASWNEW/SW2M     
       ENDIF
*     NEW RESULT FOR 2-LOOP CORRECTIONS TO KAPPA BY DUBOVYK ET AL.
       IF(IAMT4.EQ.7) THEN
*--
*- Sw_eff block:
*--
        SELECT CASE (INF)
           CASE (1)
              ASW0=  0.2308772D0
              C1  =  4.713D-4
              C2  =  2.05D-5
              C3  =  3.85D-6
              C4  = -1.85D-6
              C5  =  2.06D-2
              C6  = -2.850D-3
              C7  =  1.82D-4
              C8  = -9.71D-6
              C9  =  3.96D-4
              C10 = -6.54D-1
           CASE (2 : 4)
              ASW0=  0.2312527D0
              C1  =  4.739D-4
              C2  =  2.07D-5
              C3  =  3.85D-6
              C4  = -1.85D-6
              C5  =  2.07D-2
              C6  = -2.851D-3
              C7  =  1.82D-4
              C8  = -9.74D-6
              C9  =  3.98D-4
              C10 = -6.55D-1
           CASE (5, 7)
              ASW0=  0.2311395D0
              C1  =  4.726D-4
              C2  =  2.07D-5
              C3  =  3.85D-6
              C4  = -1.85D-6
              C5  =  2.07D-2
              C6  = -2.853D-3
              C7  =  1.83D-4
              C8  = -9.73D-6
              C9  =  3.98D-4
              C10 = -6.55D-1
           CASE (6, 8)
              ASW0=  0.2310286D0
              C1  =  4.720D-4
              C2  =  2.06D-5
              C3  =  3.85D-6
              C4  = -1.85D-6
              C5  =  2.07D-2
              C6  = -2.848D-3
              C7  =  1.81D-4
              C8  = -9.73D-6
              C9  =  3.97D-4
              C10 = -6.55D-1
        END SELECT
        DALPHA=AL4PI*DREAL(XFOTF3(IALEM,IALE2,IHVP,1,0,DAL5H,-AMZ2))
        
        DELALP=DALPHA/0.05907D0-1D0
        DELALS=CALSZ/0.117D0-1D0
        DELHIG=LOG(AMH/100D0)
        DELHSP=(AMH/100D0)
        DELTOP=(AMQ(5)/178.0D0)**2-1D0
        DELAMZ=(AMZ/91.1876D0)-1D0
*
        ASWNEW=ASW0+C1*DELHIG+C2*DELHIG**2+C3*DELHIG**4+C4*
     &   	   (DELHSP**2-1.D0)+C5*DELALP
     &          +C6*DELTOP+C7*DELTOP**2+C8*(DELHSP-1D0)*DELTOP
     &          +C9*DELALS+C10*DELAMZ
        E0  = 4.5D-6

*
* Theoretical uncertainty from the authors
*
*        IDSWW=0 ! +-1 should be new external flag with values 0,+-1 [sigma]
        ASWNEW=ASWNEW+IDSWW*E0
        AKFACI=ASWNEW/SW2M

       ENDIF
       IF(IAMT4.EQ.8) THEN
*--
*- Sw_eff block:
*--
         ASW0=  0.231464D0
         C1  =  4.616D-4
         C2  =  0.539D-4
         C3  =  -0.0737D-4
         C4  =  206D-4
         C5  =  -25.71D-4
         C6  =  4.00D-4
         C7  =  0.288D-4
         C8  =  3.88D-4
         C9  =  -6.49D-4
         C10 = -6560D-4

         DALPHA=AL4PI*DREAL(XFOTF3(IALEM,IALE2,IHVP,1,0,DAL5H,-AMZ2))

         DELALP=DALPHA/0.059D0-1D0
         DELALS=CALSZ/0.1184D0-1D0
         DELHIG=LOG(AMH/125.7D0)
         DELTOP=(AMQ(5)/173.2D0)**2-1D0
         DELAMZ=(AMZ/91.1876D0)-1D0
*
         ASWNEW=ASW0+C1*DELHIG+C2*DELHIG**2+C3*DELHIG**4+C4*
     &           DELALP+C5*DELTOP+C6*DELTOP**2+C7*DELHIG*DELTOP
     &           +C8*DELALS+C9*DELALS*DELTOP+C10*DELAMZ
         E0  = 5.6D-6

*
* Theoretical uncertainty from the authors
*
*        IDSWW=0 ! +-1 should be new external flag with values 0,+-1 [sigma]
         ASWNEW=ASWNEW+IDSWW*E0
         AKFACI=AKFACI-AKFACREF+ASWNEW/SW2M

        ENDIF
       ELSE
       IF(INF.EQ.10) IBFLA=1
       CALL VERTZW(1,INDF(INF),1D0)
       CALL ROKAPP(AQFI(INF),IBFLA,1D0,DR,DRBIG,DRREM,
     &      ROFACI,ROFACL,ROFACR,SCALER,AKFACI,AKFACL,AKFACR)
       IF(INF.EQ.10) THEN
        IBFLA=1
        CALL VERTZW(1,INDF(INF),V_TB)
        CALL ROKAPP(AQFI(INF),IBFLA,V_TB,DR,DRBIG,DRREM,
     &       ROFABI,ROFABL,ROFABR,SCALER,AKFABI,AKFABL,AKFABR)
       ENDIF
       IF((IAMT4.EQ.7).AND.(INF.EQ.10)) THEN
        ASW0=  0.232704D0
        C1  =  4.723D-4
        C2  =  1.97D-4
        C3  =  2.07D-2
        C4  = -9.733D-4
        C5  =  3.93D-4
        C6  = -1.38D-4
        C7  =  2.42D-4
        C8  = -8.10D-4
        C9  = -6.64D-1
!       Attention! Wrong value in publication, arXiv:1607.08375 formula (21)
        DELALP=DALPHA/0.059D0-1D0
        DELALS=CALSZ/0.1184D0-1D0
        DELHIG=LOG(AMH/125.7D0)
        DELTOP=(AMQ(5)/173.2D0)**2-1D0
        DELAMZ=(AMZ/91.1876D0)-1D0
*
        ASWNEW=ASW0+C1*DELHIG+C2*DELHIG**2+C3*DELALP+C4*DELTOP
     &          +C5*DELTOP**2+C6*DELTOP*DELHIG+C7*DELALS
     &          +C8*DELTOP*DELALS+C9*DELAMZ
        E0  = 1.3D-6
        ASWNEW=ASWNEW+IDSWW*E0
        AKFACI=ASWNEW/SW2M
       ENDIF
       IF((IAMT4.EQ.8).AND.(INF.EQ.10)) THEN
        ASW0=  0.232704D0
        C1  =  4.638D-4
        C2  =  0.558D-4
        C3  =  -0.0700D-4
        C4  =  207D-4
        C5  =  -9.554D-4
        C6  =  3.83D-4
        C7  =  0.179D-4
        C8  =  2.41D-4
        C9  =  -8.24D-4
        C10 = -6630D-4
        DELALP=DALPHA/0.059D0-1D0
        DELALS=CALSZ/0.1184D0-1D0
        DELHIG=LOG(AMH/125.7D0)
        DELTOP=(AMQ(5)/173.2D0)**2-1D0
        DELAMZ=(AMZ/91.1876D0)-1D0
*
        ASWNEW=ASW0+C1*DELHIG+C2*DELHIG**2+C3*DELHIG**4+C4*
     &          DELALP+C5*DELTOP+C6*DELTOP**2+C7*DELHIG*DELTOP
     &          +C8*DELALS+C9*DELALS*DELTOP+C10*DELAMZ
        E0  = 2.5D-6
        ASWNEW=ASWNEW+IDSWW*E0
        AKFACI=ASWNEW/SW2M
       ENDIF      
       AK1IM=0D0
       AR1IM=0D0
      ENDIF
*
* SUBROUTINE RETURNS EW-FORMFACTORS RO AND KAPPA 
*
      RAT=0D0
      SQR=1D0
      IF(INF.LE.4) THEN
       RAT=AML(INDL(INF))**2/AMZ2
      ELSEIF(INF.NE.9.AND.
     &  INF.GT.4.AND.ABS(QCDCOR((MAX(0,2*INDQ(INF)-1)))-1).LT.1D-8) THEN
       RAT=AMQ(INDQ(INF))**2/AMZ2
      ENDIF
      SQR=SQRT(1-4*RAT)
      IF(INF.EQ.9) THEN
        SQR=0D0
        ROFACI=0D0
        AKFACI=0D0
      ENDIF
*
      SINEFF=AKFACI*SW2M
*
      IF(INF.LE.4.OR.ABS(QCDCOR((MAX(0,2*INDQ(INF)-1)))-1).LT.1D-8) THEN
        RQCDV=1D0+0.75D0*CALQED/PI*AQFI(INF)**2
        RQCDA=RQCDV
      ELSE
        RQCDV=QCDCOR(MAX(0,2*INDQ(INF)-1))
        RQCDA=QCDCOR(MAX(0,2*INDQ(INF)  ))
      ENDIF
*
* DD-ZZ game, internal flag, not for 
*
*      print *,'IDDZZ=',IDDZZ
      IF(IDDZZ.EQ.0) THEN
        RQCDV=1D0
        RQCDA=1D0
      ENDIF
*
* GAM0I - PARTIAL WIDTH FOR I-TH CHANNEL IN THE BORN APPROXIMATION
*
       VF0L=1-4*SW2M*AQFI(INF)
       GAM0I=CONSTZ*SQR*((1+2*RAT)*(VF0L**2*RQCDV+RQCDA)/2-3*RAT*RQCDA)
*
* GAMWI - THE SAME BUT INCLUDING NON QED ELECTROWEAK 1-LOOP CORRECTIONS
* GAM1I - THE SAME BUT INCLUDING QED AND QCD CORRECTIONS TOO
*
      IF(IFACT.LE.3) THEN
       VF1L=1-4*SW2M*AQFI(INF)*AKFACI
       GAMWI=CONSTZ*ROFACI*SQR*((1+2*RAT)*(VF1L**2+1)/2-3*RAT)
       GAM1I=CONSTZ*ROFACI*SQR*((1+2*RAT)
     &     *((VF1L**2+16*SW2M**2*(AQFI(INF))**2*AK1IM**2)*RQCDV+RQCDA)/2
     &                            -3*RAT*RQCDA)
     &      +ICZAK*ARCZAK(INDQ(INF))
       GAM1IA=CONSTZ*ROFACI*SQR*((1+2*RAT)
     &      *(VF1L**2*RQCDV+RQCDA)/2-3*RAT*RQCDA)
     &      +ICZAK*ARCZAK(INDQ(INF))
       GAM1II=CONSTZ*ROFACI*SQR*(1+2*RAT)
     &      *(16*SW2M**2*(AQFI(INF))**2*AK1IM**2*RQCDV)/2
* FOR Z->BB APPLIED ONLY FOR IFACT.LE.3
       BF1L=1-4*SW2M*AQFI(INF)*AKFABI
       BAMWI=CONSTZ*ROFABI*SQR*((1+2*RAT)*(BF1L**2+1)/2-3*RAT)
       BAM1I=CONSTZ*ROFABI*SQR*((1+2*RAT)
     &     *((BF1L**2+16*SW2M**2*(AQFI(INF))**2*AK1IM**2)*RQCDV+RQCDA)/2
     &                            -3*RAT*RQCDA)
     &      +ICZAK*ARCZAK(INDQ(INF))
       BAM1IA=CONSTZ*ROFABI*SQR*((1+2*RAT)
     &      *(BF1L**2*RQCDV+RQCDA)/2-3*RAT*RQCDA)
     &      +ICZAK*ARCZAK(INDQ(INF))
       BAM1II=CONSTZ*ROFABI*SQR*(1+2*RAT)
     &      *(16*SW2M**2*(AQFI(INF))**2*AK1IM**2*RQCDV)/2
      ELSEIF(IFACT.EQ.4) THEN
       ROFACL=ROFACI
       AKFACL=AKFACI
       ROFACR=SCALER*ROFACR
       AKFACR=SCALER*AKFACR
       AKFACI=AKFACL+AKFACR
       VF1LL=1-4*SW2M*AQFI(INF)*AKFACL
       VF1LR= -4*SW2M*AQFI(INF)*AKFACR
       GAMWI=CONSTZ*SQR*(ROFACL*((1+2*RAT)*( VF1LL**2+2*VF1LL*VF1LR+1)/2
     &                                                      -3*RAT)
     &                  +ROFACR*((1+2*RAT)*( VF1LL**2+1 )/2 -3*RAT)    )
       GAM1I=CONSTZ*SQR*(ROFACL*((1+2*RAT)
     & *(((VF1LL**2+2*VF1LL*VF1LR)+16*SW2M**2*(AQFI(INF))**2*AK1IM**2)/2
     &                                    *RQCDV+RQCDA/2)-3*RAT*RQCDA)
     &                  +ROFACR*((1+2*RAT)*( VF1LL**2/2
     &                                    *RQCDV+RQCDA/2)-3*RAT*RQCDA) )
     &      +ICZAK*ARCZAK(INDQ(INF))
      ELSEIF(IFACT.EQ.5) THEN
       ROFACL=ROFACI
       AKFACL=AKFACI
       ROFACR=SCALER*ROFACR
       AKFACR=SCALER*AKFACR
       AKFACI=AKFACL+AKFACR
       VF1LL=1-4*SW2M*AQFI(INF)*AKFACL
       VF1LR= -4*SW2M*AQFI(INF)*AKFACR
       GAMWI=CONSTZ*SQR*(ROFACL*((1+2*RAT)*( VF1LL**2+2*VF1LL*VF1LR+1)/2
     &                                                      -3*RAT)
     &                  +ROFACR*((1+2*RAT)*( VF1LL**2+1 )/2 -3*RAT)    )
       GAM1I=CONSTZ*SQR*(ROFACL*((1+2*RAT)
     &      *((VF1LL**2+16*SW2M**2*(AQFI(INF))**2*AK1IM**2)/2
     &                                      *RQCDV+RQCDA/2)-3*RAT*RQCDA)
     &                  +ROFACL*((1+2*RAT)*VF1LL*VF1LR)
     &                  +ROFACR*((1+2*RAT)*( VF1LL**2+1)/2 -3*RAT)     )
     &      +ICZAK*ARCZAK(INDQ(INF))
      ENDIF
*
      NCF=1
      IF(INF.EQ.1) NCF=3
      GAM0T=GAM0T+GAM0I*MCOLFZ(INF)*NCF
      GAMWT=GAMWT+GAMWI*MCOLFZ(INF)*NCF
      IF(INF.GT.4) GAM0H =GAM0H +GAM0I *MCOLFZ(INF)*NCF
      IF(INF.GT.4) GAM1H =GAM1H +GAM1I *MCOLFZ(INF)*NCF
      IF(INF.GT.4) GAM1HA=GAM1HA+GAM1IA*MCOLFZ(INF)*NCF
      IF(INF.GT.4) GAM1HI=GAM1HI+GAM1II*MCOLFZ(INF)*NCF
      GAM1T=GAM1T+GAM1I*MCOLFZ(INF)*NCF
*
      PARTZ (INF-1)=GAM1I *1D3*MCOLFZ(INF)
      PARTZA(INF-1)=GAM1I *1D3*MCOLFZ(INF)
cbard    PARTZI(INF-1)=GAM1II*1D3*MCOLFZ(INF)
      PARTZI(INF-1)=0D0
*
      IF(INF.NE.9) THEN
        RENFAC(INF-1)=GAM1I/GAM1IA
      ELSE
        RENFAC(9)=1D0
      ENDIF
      SRENFC(INF-1)=SQRT(RENFAC(INF-1))
      ARROFZ(INF-1)=ROFACI*RENFAC(INF-1)
      AROTFZ(INF-1)=ROFACI
      ARKAFZ(INF-1)=AKFACI
      ARVEFZ(INF-1)=1D0-4D0*SW2M*AKFACI*AQFI(INF)
      ARSEFZ(INF-1)=SINEFF
      AIROFZ(INF-1)=AR1IM
      AIKAFZ(INF-1)=AK1IM
      AIVEFZ(INF-1)=-4D0*SW2M*AK1IM*AQFI(INF)
*
* GAM.T - CORRESPONDING TOTAL Z- WIDTHS WITHIN DIFFERENT APPROXIMATIONS
*
3     CONTINUE
*
* END LOOP ON FERMION FLAVORS IN Z DECAY
*
      GAMZ = GAM1T
      GAM1TZ=GAM1T
      PARTZ (9) =BAM1I*1D3*MCOLFZ(10)
      PARTZ (10)=GAM1H*1D3
      PARTZA(10)=GAM1H*1D3
cbard      PARTZI(10)=GAM1HI*1D3
      PARTZI(10)=0D0
      PARTZ (11)=GAM1T *1D3
*
* END OF Z- WIDTHS CALCULATION
*
************************************************************************
*
* W- CHAIN STARTS HERE, IT IS QUITE SIMILAR TO Z- CHAIN, FOR THIS
* REASON ONLY BRIEF ADDITIONAL COMMENTS ARE ADDED BELOW
*
      CALL VERTZW(0,0,1D0)
      AMW=SQRT(AMW2)
      CONSTW=GMU*AMW**3/6.D0/PI/SQRT(2.D0)
      GAM0T=0.D0
      GAM1T=0.D0
      DO 7 IND=1,2
      CALL PROW (AQFW(IND),ROW)
* THIS SUBROUTINE RETURNS THE ONLY ONE ELECTROWEAK FORMFACTOR ROW
* EXISTING IN THE W- DECAY CASE.
* THE OTHER IMPORTANT DIFFERENCE FROM Z- CASE IS THAT IT IS IMPOSSIBLE T
* DEFINE HERE QED- GAUGE INVARIANT SUBSET OF DIAGRAMS, FOR THIS REASON
* ONLY TOTAL 1-LOOP GAMMAS AND WIDTHS ARE CALCULATED FOR W- DECAY
* 
      GAM0I=CONSTW
      GAM1I=CONSTW*ROW*QCDCOR(IND-1)
      DELT1I=(GAM1I/GAM0I-1.D0)*100.D0
      GAM0T=GAM0T+GAM0I*MWFAC(IND)
      GAM1T=GAM1T+GAM1I*MWFAC(IND)
      PARTW(IND)= GAM1I*MWFAC(IND)*1D3
7     CONTINUE
      DELT1T=(GAM1T/GAM0T-1.D0)*100.D0
      GAMW = GAM1T
      PARTW(3)=GAM1T*1D3
*
      END
 
      SUBROUTINE ROKAPP(CH,IBFLA,V_TB,DR,DRBIG,DRREM,
     &           ROFACI,ROFACL,ROFACR,SCALER,AKFACI,AKFACL,AKFACR)
*
* Calculates RHO's and KAPPA's for partial Z-widths
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
      COMMON/CDZFLG/IHVP,IAMT4,IQCD,IMOMS,IMASS,IALEM,IMASK,IBARB,IFTJR
      COMMON/CDZSCT/ISCRE,ISCAL,IAFMT,IFACR,IFACT,IHIGS,IEWLC,ICZAK
     &             ,IHIG2,IALE2,IGFER      
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      COMMON/CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
      COMMON/CDZZWG/AMZ,AMH,GMU,A0,GAMZ,GAMW,CALSZ,CALST,CALXI,CALQED
      COMMON/CDZWSC/SL2,SQ2,W0,W0F,Z0,Z0F,DWZ0R1,DWZ0F,XWM1,XWM1F,XZM1,
     &      XZM1F,XWZ1R1,XDWZ1F,XZFM1,XZFM1F,XAMM1,XAMM1F,XWFM1,XWFM1F
      COMMON/CDZVZW/V1ZZ,V1ZW,V2ZWW,V1WZ,V2WWZ,VTB
      COMMON/CDZ513/DAL5H
*
* Basic EW RHO and KAPPA
*
      CH2=CH*CH
      W0A  =W0+W0F
      ZM1A =DREAL(XZM1 +XZM1F )
      ZFM1A=DREAL(XZFM1+XZFM1F)
      WM1A =DREAL(XWM1 +XWM1F )
      AMM1A=DREAL(XAMM1+XAMM1F)
*
      UFF=(0.5D0/R-3.D0/R*CH*R1+6.D0*R12/R*CH2)*V1ZZ
     &   +(1.D0-2.D0*R-2.D0*CH*R1)*V1ZW
     &   +2.D0*R*V2ZWW+2.D0*VTB
*
      RO1=AL4PI/R1*(ZM1A+ZFM1A-W0A+5.D0/8.D0*R2+5.D0/8.D0*R
     &   -11.D0/2.D0-9.D0/4.D0*R/R1*ALR+UFF)
*
      AK1=AL4PI/R1*(R/R1*(ZM1A-WM1A)+AMM1A+R12/R*CH2*V1ZZ-.5D0*UFF)
*
* MIXED QCD-CORRECTIONS TO Z-WIDTH
*
      ALSZ=CALSZ
      ALST=CALST
      SW2=R1
      CW2=R
      SSZ=AMZ2
      ROQCD = 0D0
      AKQCD = 0D0
      AMT2=AMQ(5)**2
*
* Mixed QCD-corrections
*
      ROQCD=0D0
      AKQCD=0D0
*
      IF    (IQCD.EQ.1) THEN
        ROQCD=AL4PI*DREAL(XRQCDS(ALSZ,ALST,AMZ2,AMW2,AMT2,SSZ))
        AKQCD=AL4PI*DREAL(XKQCDS(ALST,AMZ2,AMW2,AMT2,SSZ))
      ELSEIF(IQCD.EQ.2) THEN
        ROQCD=AL4PI*DREAL(XROQCD(ALSZ,ALST,AMZ2,AMW2,AMT2,SSZ))
        AKQCD=AL4PI*DREAL(XKAQCD(ALST,AMZ2,AMW2,AMT2,SSZ))
      ELSEIF(IQCD.EQ.3) THEN
        ROQCD=AL1PI*ALST/PI*DREAL(XRMQCD(AMZ2,AMW2,AMT2,SSZ))
     &       +AL1PI*ALSZ/PI/8D0/SW2/CW2*(VT2+VB2+2D0)
* light quarks are added (25/02/1998), TWO doublets --> 2
        AKQCD=AL1PI*ALST/PI*DREAL(XKMQCD(AMZ2,AMW2,AMT2,SSZ))
     &       +AL1PI*ALSZ/PI*CW2/2D0/SW2**2*LOG(CW2)
      ENDIF
*
C-----------------------------------------------------------------------
C 13/10/1992 - Barbieri's m_t^4 are implemented
      IF(IBARB.EQ.-1) THEN
       AMT4C=0D0
       AMT4B=0D0
        ELSEIF(IBARB.EQ.0) THEN
       AMT4C=19-2D0*PI2
       AMT4B=(27-PI2)/3
        ELSEIF(IBARB.EQ.1) THEN
       RBTH=AMT2/AMH2
       ALRB=LOG(RBTH)
       AMT4C=49D0/4D0+PI2+27D0/2D0*ALRB+3D0/2D0*ALRB**2
     &      +RBTH/3D0*(2D0-12D0*PI2+12D0*ALRB-27D0*ALRB**2)
     &  +RBTH**2/48D0*(1613-240*PI2-1500*ALRB-720 *ALRB**2)
       AMT4B=1D0/144*(311D0+24*PI2+282*ALRB+90*ALRB**2
     &      -4D0*RBTH*(40D0+ 6*PI2+ 15*ALRB+18*ALRB**2)
     &      +3D0*RBTH**2*(242.09D0-60*PI2-454.2D0*ALRB-180*ALRB**2))
        ELSEIF(IBARB.EQ.2) THEN
       RBARB=SQRT(AMH2/AMT2)
       AMT4C=FBARB (RBARB)
       AMT4B=FBARBB(RBARB)
      ENDIF
*
      TOPX2 = GMU*AMT2/DSQRT(2.D0)/8.D0/PI2
      IF(IFTJR.EQ.0) THEN
        TAUBB1=-2*TOPX2
      ELSE
        TAUBB1=-2*TOPX2*(1-PI/3*CALST)
      ENDIF
      TAUBB2=-2*TOPX2*(TOPX2*AMT4B)*V_TB**2
      CORBB =+AL1PI/8/SW2*AMT2/AMW2*V_TB**2
*
      ROFACI=1D0+RO1+ROQCD
      AKFACI=1D0+AK1+AKQCD
      ROFACM=ROFACI
      AKFACM=AKFACI
*
      DALFA=AL4PI*DREAL(XFOTF3(IALEM,IALE2,IHVP,IQCD,0,DAL5H,-AMZ2))
*
      RENORD=SQRT(2D0)*GMU*AMZ2*R1*R/PI*ALFAI
      RENORM=SQRT(2D0)*GMU*AMZ2*R1*R/PI*ALFAI
      IF(ISCRE.EQ.0) SCALER=1.00D0
      IF(ISCRE.EQ.1) SCALER=RENORD
      IF(ISCRE.EQ.2) SCALER=RENORM
*
*---- without Zbb-vertex -----------------------------------------------
*
      ROFACI=1D0+RO1+ROQCD+2*CORBB*IBFLA
      AKFACI=1D0+AK1+AKQCD-  CORBB*IBFLA
*
* IF(AMT4)
*
      IF (IAMT4 .EQ. 1) THEN
       DRHOT = .75D0*AL4PI/SW2/R*AMT2/AMZ2
       TOPX2 = GMU*AMT2/DSQRT(2.D0)/8.D0/PI2
       DRHOT4=3D0*TOPX2*(1D0+TOPX2*AMT4C)
       IF(IQCD.EQ.0) THEN
        TBQCD0=0D0
        TBQCD3=0D0
        CORRXI=0D0
         ELSE
          IF(IAFMT.EQ.0) THEN
* Below is the leading term as it is given in \zf (3.71)
           TBQCD0=.75D0*AL4PI/R/SW2*AMT2/AMZ2
     &           *(-CALXI/PI*2D0/3D0*(1D0+PI2/3D0))
           TBQCDL=.75D0*AL4PI/R/SW2*AMT2/AMZ2
     &           *(-CALST/PI*2D0/3D0*(1D0+PI2/3D0))
           CORRXI=-TBQCDL+TBQCD0
           TBQCD3=0D0
            ELSE
* This coincides with TBQCD0 above if O(\alpha\alpha^2_s)=0
           TBQCD0=AFMT3(CALST,AMT2,AMZ2,SW2)
* Here TBQCDR has to be called!!! ROKAPP - 1
* Below is pure AFMT-correction to \Delta \rho (see p.18 PCWG-notebook)
           TBQCD3=.75D0*AL4PI/R/SW2*AMT2/AMZ2
     &          *(TBQCD0-(-CALST/PI*2D0/3D0*(1D0+PI2/3D0)))
           CORRXI=0D0
          ENDIF
       ENDIF
       ROFACL=1/(1-DRHOT4)
       AKFACL=(1+R/SW2*DRHOT4)
       ROFACR=ROFACI-      DRHOT-1 +       TBQCD3+CORRXI
       AKFACR=AKFACI-R/SW2*DRHOT-1 +R/SW2*(TBQCD3+CORRXI)
*
      ELSEIF(IAMT4 .EQ. 2) THEN
*
       DRHOT = .75D0*AL4PI/SW2/R*AMT2/AMZ2
       TOPX2 = GMU*AMT2/DSQRT(2.D0)/8.D0/PI2
       DRHOT4=3D0*TOPX2*(1D0+TOPX2*AMT4C)
       IF(IQCD.EQ.0) THEN
        TBQCD0=0D0
        TBQCDL=0D0
         ELSE
          IF(IAFMT.EQ.0) THEN
           TBQCD0=-TOPX2*CALXI/PI*2D0*(1D0+PI2/3D0)
            ELSE
           TBQCD0=3*TOPX2*AFMT3(CALST,AMT2,AMZ2,SW2)
* Here TBQCDR has to be called!!! ROKAPP - 2
          ENDIF
        TBQCDL=-AL4PI*ALST/PI*AMT2/AMW2/R1*(.5D0+PI2/6D0)
       ENDIF
       ROFACL=1/(1-DRHOT4-TBQCD0)
       AKFACL=(1+R/SW2*(DRHOT4+TBQCD0))
       ROFACR=ROFACI       -DRHOT-TBQCDL -1
       AKFACR=AKFACI-R/SW2*(DRHOT+TBQCDL)-1
*
      ELSEIF(IAMT4 .GE. 3) THEN
*
       DWZ1AL=R/R1*DREAL(XWZ1R1+XDWZ1F)
       RENORM=SQRT(2D0)*GMU*AMZ2*R1*R/PI*ALFAI
       SCALE = AL4PI/R1*(41D0/6D0-11D0/3D0*R)*ALR
       CORKAP=(AL4PI*DWZ1AL+SCALE)+.75D0*AL4PI/SW2**2*AMT2/AMZ2
       DRHOT =.75D0*AL4PI/SW2/R*AMT2/AMZ2
       TOPX2 =GMU*AMT2/DSQRT(2.D0)/8.D0/PI2
       DRHOT4=3D0*TOPX2*(1D0+TOPX2*AMT4C)
       IF(IQCD.EQ.0) THEN
        TBQCD0=0D0
        TBQCDL=0D0
         ELSE
          IF(IAFMT.EQ.0) THEN
           TBQCD0=-TOPX2*CALXI/PI*2D0*(1D0+PI2/3D0)
            ELSE
           TBQCD0=3*TOPX2*AFMT3(CALST,AMT2,AMZ2,SW2)
* Here TBQCDR has to be called!!! ROKAPP - 3
          ENDIF
        TBQCDL=-AL4PI*ALST/PI*AMT2/AMW2/R1*(.5D0+PI2/6D0)
       ENDIF
       ROFACL=1/(1-DRHOT4-TBQCD0)
       AKFACL=(1+R/SW2*(DRHOT4+TBQCD0)-CORKAP*RENORM)
       ROFACR=ROFACI       -DRHOT-TBQCDL        -1
       AKFACR=AKFACI-R/SW2*(DRHOT+TBQCDL)+CORKAP-1
      ENDIF
C-----------------------------------------------------------------------
      IF(IBFLA.EQ.1) THEN
       ROFACL=ROFACL*(1+(TAUBB1+TAUBB2)*IBFLA)**2
       AKFACL=AKFACL/(1+(TAUBB1+TAUBB2)*IBFLA)
      ENDIF
C-----------------------------------------------------------------------
      IF(IAMT4.LE.3.OR.IBFLA.EQ.1) THEN
      IF(IEWLC.EQ.1) THEN
       IF(IFACT.EQ.0)THEN
        ROFACI=1/(1/ROFACL-SCALER*ROFACR)
        AKFACI=AKFACL*(1+SCALER*AKFACR)
       ELSEIF(IFACT.EQ.1) THEN
        ROFACI=ROFACL*(1+SCALER*ROFACR*ROFACL)
        AKFACI=AKFACL+SCALER*AKFACR
       ELSEIF(IFACT.EQ.2) THEN
        ROFACI=ROFACL*(1+SCALER*ROFACR)
        AKFACI=AKFACL+SCALER*AKFACR
       ELSEIF(IFACT.EQ.3) THEN
        ROFACI=ROFACL+SCALER*ROFACR
        AKFACI=AKFACL+SCALER*AKFACR
         ELSE
        ROFACI=ROFACL
        AKFACI=AKFACL
       ENDIF
      ELSEIF(IEWLC.EQ.0) THEN
       ROFACI=ROFACL
       AKFACI=AKFACL
      ELSE
       STOP
      ENDIF
      ENDIF
*
      END
 
      SUBROUTINE 
     &        FOKAPP(CH,DR,DRBIG,DRREM,SCALER,ROFACI,AKFACI,AR1IM,AK1IM)
*
* Calculates KAPPA's for Degrassi's sin^2\theta^{lept}_{eff}
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
      COMMON/CDZFLG/IHVP,IAMT4,IQCD,IMOMS,IMASS,IALEM,IMASK,IBARB,IFTJR
      COMMON/CDZSCT/ISCRE,ISCAL,IAFMT,IFACR,IFACT,IHIGS,IEWLC,ICZAK
     &             ,IHIG2,IALE2,IGFER      
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      COMMON/CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
      COMMON/CDZZWG/AMZ,AMH,GMU,A0,GAMZ,GAMW,CALSZ,CALST,CALXI,CALQED
      COMMON/CDZWSC/SL2,SQ2,W0,W0F,Z0,Z0F,DWZ0R1,DWZ0F,XWM1,XWM1F,XZM1,
     &      XZM1F,XWZ1R1,XDWZ1F,XZFM1,XZFM1F,XAMM1,XAMM1F,XWFM1,XWFM1F
      COMMON/CDZVZW/V1ZZ,V1ZW,V2ZWW,V1WZ,V2WWZ,VTB
      COMMON/CDZDEG/DROBAR,DROBLO,DRREMD
      COMMON/CDZ513/DAL5H
*
* Basic EW RHO and KAPPA
*
      CH2=CH*CH
      W0A  =W0+W0F
      XZM1A =XZM1 +XZM1F 
      XZFM1A=XZFM1+XZFM1F
      XWM1A =XWM1 +XWM1F 
      XAMM1A=XAMM1+XAMM1F
*
      V1ZIM=DIMAG(XV1B(-AMZ2,AMZ2))
      V1WIM=DIMAG(XV1B(-AMZ2,AMW2))
      XV1ZZ=DCMPLX(V1ZZ,V1ZIM)
      XV1ZW=DCMPLX(V1ZW,V1WIM)
*
      XUFF=0.25D0/R*(1D0-6D0*CH*R1+12D0*R12*CH2)*XV1ZZ
     &    +(0.5D0-R-CH*R1)*XV1ZW+R*V2ZWW+VTB
*
      XRO1=AL4PI/R1*(DREAL(XZM1A+XZFM1A)-W0A+5.D0/8.D0*R*(1D0+R)
     &   -11.D0/2.D0-9.D0/4.D0*R/R1*ALR+2D0*XUFF)
*
      XAK1=AL4PI/R1*(R/R1*DREAL(XZM1A-XWM1A)
     &    +XAMM1A+R12/R*CH2*XV1ZZ-XUFF)
*
      RO1=SQRT((1D0+DREAL(XRO1))**2+(DIMAG(XRO1))**2)
      AK1=1D0+DREAL(XAK1)
*
* Mixed QCD-corrections to Z-widths
*
      ALSZ=CALSZ
      ALST=CALST
      SW2=R1
      CW2=R
      SSZ=AMZ2
      ROQCD = 0D0
      AKQCD = 0D0
      AMT2=AMQ(5)**2
*
      AKMIX=AL1PI*ALSZ/4D0/SW2*(7D0/3D0-44D0/9D0*SW2) 
      AR1IM=DIMAG(XRO1)
      AK1IM=DIMAG(XAK1)+AKMIX        
*
      ROQCD=0D0
      AKQCD=0D0
*
      IF    (IQCD.EQ.1) THEN
        ROQCD=AL4PI*DREAL(XRQCDS(ALSZ,ALST,AMZ2,AMW2,AMT2,SSZ))
        AKQCD=AL4PI*DREAL(XKQCDS(ALST,AMZ2,AMW2,AMT2,SSZ))
      ELSEIF(IQCD.EQ.2) THEN
        ROQCD=AL4PI*DREAL(XROQCD(ALSZ,ALST,AMZ2,AMW2,AMT2,SSZ))
        AKQCD=AL4PI*DREAL(XKAQCD(ALST,AMZ2,AMW2,AMT2,SSZ))
      ELSEIF(IQCD.EQ.3) THEN
        ROQCD=AL1PI*ALST/PI*DREAL(XRMQCD(AMZ2,AMW2,AMT2,SSZ))
     &       +AL1PI*ALSZ/PI/8D0/SW2/CW2*(VT2+VB2+2D0)
* light quarks are added (25/02/1998), TWO doublets --> 2
        AKQCD=AL1PI*ALST/PI*DREAL(XKMQCD(AMZ2,AMW2,AMT2,SSZ))
     &       +AL1PI*ALSZ/PI*CW2/2D0/SW2**2*LOG(CW2)
      ENDIF
*
*-----------------------------------------------------------------------
*
      ROFACI=RO1+ROQCD
      AKFACI=AK1+AKQCD
      ROFACM=ROFACI
      AKFACM=AKFACI
*
      DALFA=DREAL(XFOTF3(IALEM,IALE2,IHVP,IQCD,0,DAL5H,-AMZ2))*AL4PI
*
* DALFA is used only in the scale of ADDIM, elsewhere one uses RENORM
*
*---- without Zbb-vertex -----------------------------------------------
*
       RENORM=SQRT(2D0)*GMU*AMZ2*R1*R/PI*ALFAI
* the same scale as in SEARCH
       SCALEB = AL4PI/R1*(1D0/6D0+7D0*R)*ALR
       ANUMF  = 24D0
       TRQ2F  =  8D0
       SCALEF = -AL4PI/R1*(ANUMF/6-4D0/3*R1*TRQ2F)*ALR
       SCALE=SCALEB+SCALEF
*
* to CORrect reminders:
*
       CORKAP=AL4PI*R/R1*(XDWZ1F+XWZ1R1)+SCALE
       CORRHO=R1/R*CORKAP
       DRHOT =.75D0*AL4PI/SW2/R*AMT2/AMZ2
       TOPX2 =GMU*AMT2/DSQRT(2.D0)/8.D0/PI2
*
       IF(IQCD.EQ.0) THEN
        TBQCD0=0D0
        TBQCDL=0D0
        DKTB3R=0D0
         ELSE
          IF(IAFMT.EQ.0) THEN
           TBQCD0=-TOPX2*CALXI/PI*2D0*(1D0+PI2/3D0)
            ELSE
           TBQCD0=3*TOPX2*AFMT3(CALST,AMT2,AMZ2,SW2)
* Here TBQCDR is called for FOKAPP
           DKTB3R=3*TOPX2*TBQCDR(CALST,AMT2,AMZ2,SW2)
          ENDIF
        TBQCDL=-AL4PI*ALST/PI*AMT2/AMW2/R1*(.5D0+PI2/6D0)
       ENDIF
*
       ROFACR=ROFACI      -TBQCDL+CORRHO-1
       AKFACR=AKFACI-R/SW2*TBQCDL+CORKAP-1+DKTB3R
*
       AMW=SQRT(AMW2)
       AMT=AMQ(5)
       PI3QF=ABS(CH)
*    
       CALL GDEGNL
     & (GMU,AMZ,AMT,AMH,AMW,PI3QF,AMZ2,DRDREM,DRHOD,DKDREM,DROREM)
*
       DF1BAR=RENORM*DRREMD
*
* New game with SCALER, April 99  
*
      IF(ISCRE.EQ.0) SCALER2=1D0
      IF(ISCRE.EQ.1) SCALER2=1D0/RENORM**2
      IF(ISCRE.EQ.2) SCALER2=1D0*RENORM**2
*
      IF(IFACT.EQ.0)     THEN
       ROFAC=(1D0+RENORM*ROFACR+DROREM*SCALER2)
     &                          /(1D0-DROBAR*(1D0-DF1BAR))       
      ELSEIF(IFACT.EQ.1) THEN
       ROFAC=(1D0+RENORM*ROFACR)/(1D0-DROBAR*(1D0-DF1BAR))
     &                         +DROREM*SCALER2       
      ELSEIF(IFACT.EQ.2) THEN
       ROFAC=1D0+DROBAR-DROBLO*DF1BAR+DROBLO**2
     &           +RENORM*ROFACR*(1D0+DROBLO)
     &                         +DROREM*SCALER2       
      ENDIF
      ROFACI=ROFAC
      IF(IFACT.EQ.0) THEN
       AKFAC=(1D0+RENORM*AKFACR+DKDREM*SCALER2)
     &                          *(1D0+R/SW2*DROBAR*(1D0-DF1BAR))
      ELSEIF(IFACT.EQ.1) THEN
       AKFAC=(1D0+RENORM*AKFACR)*(1D0+R/SW2*DROBAR*(1D0-DF1BAR))
     &                         +DKDREM*SCALER2
      ELSEIF(IFACT.EQ.2) THEN
       AKFAC=1D0+R/SW2*DROBAR-R/SW2*DROBLO*DF1BAR
     &           +RENORM*AKFACR *(1D0+R/SW2*DROBLO)
     &                         +DKDREM*SCALER2
      ENDIF
*
      ADDIM=1D0/ALFAI**2/(1D0-DALFA)**2*35D0/18*(1D0-8D0/3*AKFAC*SW2)
      AKFACI=AKFAC+ADDIM/SW2 
      SINEFF=AKFACI*SW2
*cbardprint *,'CH,ADDIM,SINEFF=',CH,ADDIM,SINEFF 
*
* iterations are abandoned
*
      NOITER=0
      IF(NOITER.EQ.0) RETURN
      ADDIM=1D0/ALFAI**2/(1D0-DALFA)**2*35D0/18*(1D0-8D0/3*SINEFF)
      AKK=AL4PI*(DREAL(XAMFEF(-AMZ2,SINEFF))-DREAL(XAMM1F))+ADDIM
      AKFACI=AKFAC+AKK/SW2
      SINEFF=AKFACI*SW2
*
      ADDIM=1D0/ALFAI**2/(1D0-DALFA)**2*35D0/18*(1D0-8D0/3*SINEFF)
      AKK=AL4PI*(DREAL(XAMFEF(-AMZ2,SINEFF))-DREAL(XAMM1F))+ADDIM
      AKFACI=AKFAC+AKK/SW2
      SINEFF=AKFACI*SW2
*
      END
 
      SUBROUTINE VERTZW(MZ,INDF,V_TB)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      COMMON/CDZVZW/V1ZZ,V1ZW,V2ZWW,V1WZ,V2WWZ,VTB
*
      IF(MZ) 5,8,5
*
* Z-BOSON CHAIN *********************
* FILLS CDZVZW (Z PART)
*
5     SR=SQRT(4.D0*R-1.D0)
      AT=ATAN(SR/(2.D0*R-1.D0))
      V1ZZ=-5.5D0-8.D0*(F1-SPENCE(2.D0))
      SPERR=SPENCE(1.D0+1.D0/R)
      V1ZW=-3.5D0-2.D0*R-(3.D0+2.D0*R)*ALR-2.D0*(1.D0+R)**2*(F1-SPERR)
      V2ZWW=2.D0/9.D0/R2+43.D0/18.D0/R-1.D0/6.D0-2.D0*R
     *     +(-1.D0/12.D0/R2-1.5D0/R+7.D0/3.D0+2.D0*R)*SR*AT
     *     -2.D0*R*(2.D0+R)*AT**2
*
      IF(INDF.EQ.5) THEN
*
       CALL VTBANA(1,AMZ2,WWv2,WWv11,WWv12)
       QBM=1D0/3D0
      VTB=V_TB**2*(R*WWv2-.5D0*(1D0-2D0*R1*(1D0-QBM))*WWv11-.5D0*WWv12)      
*
      ELSE
       VTB=0.D0
      ENDIF
      GO TO 9
*
* W-BOSON CHAIN *********************
* FILLS CDZVZW (W PART)
*
8     ALAM=AMZ2*AMZ2-4.D0*AMW2*AMZ2
      V1WZ=-5.D0-2.D0/R+(3.D0+2.D0/R)*ALR
     *    -2.D0*R1W2/R2*(SPENCE(1.D0)-SPENCE(R1W))
      V2WWZ=-9.D0/4.D0/R-1.D0/12.D0/R2+23.D0/18.D0
     *     +(1.D0/2.D0/R-3.D0/4.D0/R2
     *     -1.D0/24.D0/R/R2+1.D0)*ALR
     *     -DREAL(XL(-AMW2,AMW2,AMZ2))
     *     *(5.D0/6.D0/R+1.D0/24.D0/R2+1.D0/2.D0)/AMW2
     *     +(1.D0/2.D0+1.D0/R)*ALAM*DREAL(XJ(-AMW2,AMW2,AMZ2))
     *     *DREAL(XJ(-AMW2,AMW2,AMZ2))-(1.D0/2.D0+1.D0/R)*ALR*ALR
9     CONTINUE
*
      END

      SUBROUTINE VTBANA(NUNI,S,WWv2,WWv11,WWv12)   
*
* Supplies all ingredients for off-resonance EW finite m_t corrections
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      COMMON/CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
*
      RWS =AMW2/S
      AMT2=AMQ(5)**2
      RTW =AMT2/AMW2
      RTW1=RTW-1D0
      ALRT=LOG(RTW)
*
      CALL S3WANA(AMT2,AMW2,-S,XJ0W,XS3W,XS3W0)
      CALL S3WANA(AMW2,AMT2,-S,XJ0T,XS3T,XS3T0)
*
      AJ0W=DREAL(XJ0W)
      AJ0T=DREAL(XJ0T)
      S3W =DREAL(XS3W)
      S3W0=DREAL(XS3W0)
      S3T =DREAL(XS3T)
      S3T0=DREAL(XS3T0)
*
      WWv2 =-2D0*RWS*(2D0+RWS)*S*(S3w-S3w0)
     &      +RTW*((3D0*RWS**2+2.5D0*RWS-2D0-(2D0*RWS-.5D0)*RTW
     &           +RWS*(.5D0-RWS)*RTW**2)*S*S3w
     &           -(RWS+1D0-(.5D0-RWS)*RTW)*(AJ0w-2D0)
     &           +(2D0*RWS+3D0/2/RTW1**2-2D0/RTW1+1D0/2
     &           -(.5D0-RWS)*RTW)*ALRT
     &           -(RWS+3D0/2/RTW1+3D0/4-(.5D0-RWS)*RTW)
     &           +.25D0/RWS*(AJ0w-3D0)*NUNI
     &           )
      WWv11=+2D0*(1D0+RWS)**2*S*(S3t-S3t0)
     &      +(2D0*RWS+3D0)*(AJ0t+ALRT+LOG(RWS))
     &      -RTW*(RWS*(3D0*RWS+2D0-RTW-RWS*RTW**2)*S*S3t
     &           +(RWS+.5D0+RWS*RTW)*(AJ0t+ALRT-2D0)
     &           -(2D0*RWS+3D0/2/RTW1**2-2D0/RTW1+1D0/2+RWS*RTW)*ALRT
     &           +RWS+3D0/2/RTW1+5D0/4+RWS*RTW
     &           ) 
      WWv12=-RTW*(RWS*(2D0+RWS-2D0*RWS*RTW+RWS*RTW**2)*S*S3t
     &           -(.5D0-RWS+RWS*RTW)*(AJ0t+ALRT-1D0)+RWS*RTW*ALRT
     &           )
*
      END
 
      SUBROUTINE PROW (QI,ROW)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      COMMON/CDZWSC/SL2,SQ2,W0,W0F,Z0,Z0F,DWZ0R1,DWZ0F,XWM1,XWM1F,XZM1,
     &      XZM1F,XWZ1R1,XDWZ1F,XZFM1,XZFM1F,XAMM1,XAMM1F,XWFM1,XWFM1F
      COMMON/CDZVZW/V1ZZ,V1ZW,V2ZWW,V1WZ,V2WWZ,VTB
*
      QIQJ=QI*(1.D0-QI)
      WM1A=DREAL(XWM1+XWM1F)
      W0A=W0+W0F
      WFM1A=DREAL(XWFM1+XWFM1F)
      ROW=1.D0+AL4PI/R1*(WM1A-W0A+WFM1A-7.D0/1.D0+5.D0/8.D0*R*R1W
     *   -9.D0/4.D0*R/R1*ALR+3.D0/4.D0/R+3.D0*R-3.D0/R*R12*QIQJ
     *   +(1.D0/2.D0/R-1.D0-2.D0*R12/R*QIQJ)*V1WZ
     *   +2.D0*R*V2WWZ+2.D0*R1*(77.D0/12.D0-2.D0/3.D0*PI2+109.D0/36.D0
     *   -3.D0/2.D0*QIQJ))
      PROW1=100.D0*(ROW-1.D0)
*
      END

      SUBROUTINE ROKANC(IBOXF,IBFLA,S,Q2,U,QI,QJ,XROK,XFOT,XFOT5)
*
* BEFORE USE OF ROKANC, AT LEAST ONE CALL OF DIZET MUST BE DONE.
* SEE ALSO THE COMMENTS THERE.
*---------------------------------------------------------------------
* THIS ROUTINE CALCULATES THE WEAK NEUTRAL CURRENT FORM FACTORS FOR
* THE 4-FERMION SCATTERING CROSS SECTION. ALSO: THE RUNNING ALPHA.QED.
*----------------------------------------------------------------------
* INPUT FROM USER:
*            S,Q2,U - THE KINEMATIC INVARIANTS FOR THE QUARK PROCESS
*                     (S+T-U=0)
*             QI,QJ - THE CHARGES OF THE FERMION PAIRS IN THE PROCESS
*             IBOXF -    FLAG FOR THE WW,ZZ-BOX CONTRIBUTIONS
*             IBOXF = 0: THEY ARE SET EQUAL ZERO. NOT BAD FOR LEP100.
*                     1: THEY ARE CALCULATED.
* SPECIAL HANDLING OF WEAK CROSS SECTION FORM FACTORS IN CASE OF B-QUARK
*             IBFLA = 0: ALL OTHER CHANNELS
*                   = 1: THE B-QUARK PRODUCTION CHANNEL IN ANNIHILATION
*                      WITH THIS FLAG, THE ADDITIONAL VERTEX CORRECTIONS
*                      DUE TO THE T-QUARK MASS ARE TAKEN INTO ACCOUNT
*                      FOR LEP PHYSICS. SEE REF. 2.
*-----------------------------------------------------------------------
* OUTPUT OF THE ROUTINE:
*               XROK - THE FOUR COMPLEX NEUTRAL CURRENT FORM FACTORS
*                      RHO, KAPPA.I, KAPPA.J, KAPPA.IJ
*               XFOT - THE COMPLEX RUNNING QED COUPLING CONSTANT AT
*                      SCALE Q2
*              XFOT5 - THE XFOT, BUT TAKING INTO ACCOUNT ONLY 5 QUARKS
*-----------------------------------------------------------------------
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      COMMON/CDZFLG/IHVP,IAMT4,IQCD,IMOMS,IMASS,IALEM,IMASK,IBARB,IFTJR
      COMMON/CDZSCT/ISCRE,ISCAL,IAFMT,IFACR,IFACT,IHIGS,IEWLC,ICZAK
     &             ,IHIG2,IALE2,IGFER      
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
      COMMON/CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
      COMMON/CDZZWG/AMZ,AMH,GMU,A0,GAMZ,GAMW,CALSZ,CALST,CALXI,CALQED
      COMMON/CDZWSC/SL2,SQ2,W0,W0F,Z0,Z0F,DWZ0R1,DWZ0F,XWM1,XWM1F,XZM1,
     &      XZM1F,XWZ1R1,XDWZ1F,XZFM1,XZFM1F,XAMM1,XAMM1F,XWFM1,XWFM1F
      COMMON/CDZVZW/V1ZZ,V1ZW,V2ZWW,V1WZ,V2WWZ,VTB
      COMMON/CDZXKF/XROKF
      COMMON/CDZRLR/ROKL(4),ROKR(4),AROK(4)
      COMMON/CDZDEG/DROBAR,DROBLO,DRREMD
      COMMON/CDZ513/DAL5H
      COMMON/CDZ_LK/IMTADD
      COMMON/CDZVTB/V_TB
      COMMON /CDZRKZ/ARROFZ(0:10),ARKAFZ(0:10),ARVEFZ(0:10),ARSEFZ(0:10)
     &              ,AROTFZ(0:10),AIROFZ(0:10),AIKAFZ(0:10),AIVEFZ(0:10)
*
      DIMENSION XROK(4),XROKIN(4)
*
      AMT2=AMQ(5)**2
      QIM=ABS(QI)
      QJM=ABS(QJ)
      SI=1.D0
      SJ=1.D0
      IF(QIM.NE.0.D0)  SI=QI/QIM
      IF(QJM.NE.0.D0)  SJ=QJ/QJM
      VI=1.D0-4.D0*R1*QIM
      VJ=1.D0-4.D0*R1*QJM
      XDF=DREAL(XDZF(Q2))
      XDB=XDZB(Q2)
      XDBAL=XDF+XDB
      XZM1AL=XZM1+XZM1F
      XWZ1AL=R*XWZ1R1+R*XDWZ1F
      XFMF=XAMF(Q2)
      XV1BW=XV1B(Q2,AMW2)
      XV1BZ=XV1B(Q2,AMZ2)
      XA1BW=XA1B(Q2,AMW2)
      XV2BW=XV2B(Q2,AMW2)
      XRFL3=XL(Q2,AMW2,AMW2)/Q2
      XROBZ=17D0/6D0-5D0/6D0*XRFL3
      IF(IBFLA.EQ.1) THEN
        XROBT=AMT2/AMW2/4*(3D0-1D0/2*XRFL3)
      ELSE
        XROBT=(0D0,0D0)
      ENDIF
      Q2M=ABS(Q2)
      ALSZ=CALSZ
      ALST=CALST
      SW2=R1
      CW2=R
*
* Mixed QCD-corrections
*
      XZERO=DCMPLX(0.D0,0.D0)
      XRQCD=XZERO
      XKQCD=XZERO
*
      IF    (IQCD.EQ.1) THEN
        XRQCD=AL4PI*XRQCDS(ALSZ,ALST,AMZ2,AMW2,AMT2,Q2M)
        XKQCD=AL4PI*XKQCDS(ALST,AMZ2,AMW2,AMT2,Q2M)
      ELSEIF(IQCD.EQ.2) THEN
        XRQCD=AL4PI*XROQCD(ALSZ,ALST,AMZ2,AMW2,AMT2,Q2M)
        XKQCD=AL4PI*XKAQCD(ALST,AMZ2,AMW2,AMT2,Q2M)
      ELSEIF(IQCD.EQ.3) THEN
* are added (25/02/1998) and tested (05/12/1998)
       IF(ABS(1D0-Q2M/AMZ2).LT.1D-3) THEN
        XRQCD=AL1PI*ALST/PI*XRMQCD(AMZ2,AMW2,AMT2,Q2M)
     &       +AL1PI*ALSZ/PI/8D0/SW2/CW2*(VT2+VB2+2D0)
       ELSE
        XRQCD=AL1PI*ALST/PI*XRMQCD(AMZ2,AMW2,AMT2,Q2M)
     &       -AL1PI*ALSZ/PI/8D0/SW2/CW2*(VT2+VB2+2D0)
     &                   *Q2M/(AMZ2-Q2M)*LOG(Q2M/AMZ2)
       ENDIF
* light quarks are added (25/02/1998) and tested (05/12/1998)  
* Qum^2+Qdm^2=5/9 and two doublets
        XKQCD=AL1PI*ALST/PI*XKMQCD(AMZ2,AMW2,AMT2,Q2M)
     &       +DCMPLX(AL1PI*ALSZ/PI/2D0/SW2**2*(CW2*LOG(CW2)
     &                -SW2*(1D0-20D0/9D0*SW2)*LOG(Q2M/AMZ2))
     &              ,AL1PI*ALSZ/8D0/3D0/SW2**2*(2D0*SW2-1D0))      
      ENDIF
*
*  XROK(1)=RO WITH INDEXES I AND J
*
*  FROW BOXZZ
      AI11=-S
      AI12=-U
      SB=-Q2
      IF (IBOXF.EQ.0) THEN
       XWWRO=0.D0
       XZZRO=0.D0
       XZZI =0.D0
       XZZJ =0.D0
       XZZIJ=0.D0
      ELSE
      XWWP= (1D0+SI*SJ)/2D0*(XBOX(IBFLA,AI11,AI12,AMW2)
     &     +4.D0*(AI11/SB)**2*XJ3(SB,AMW2) )
* This is with March'95 corrections
* Still could be wrong! Attention of sign!!!
      XWWM=-(1D0-SI*SJ)/2D0*2.D0*(AI11/SB)**2
     &     *(AI11/SB*XJ4(SB,AI11,AMW2)+2.D0*XJ3(SB,AMW2))
      XZZP=  XBOX(0,AI11,AI12,AMZ2)
     &     -2.D0*(AI11/SB)**3*XJ4(SB,AI11,AMZ2)
      XZZM=-(XBOX(0,AI12,AI11,AMZ2)
     &     -2.D0*(AI12/SB)**3*XJ4(SB,AI12,AMZ2))
      XWWPL=(-Q2)/AI11**2*XWWP
      XWWMI=(-Q2)/AI11**2*XWWM
      XZZPL=(-Q2)/AI11**2*XZZP
      XZZMI=(-Q2)/AI12**2*XZZM
      XWWRO=-SI*SJ*R*(Q2+AMZ2)*(XWWMI+XWWPL)
      XZZRO=-SI*SJ/32D0/R*(Q2+AMZ2)*
     &                 (((1D0+VI**2)*(1D0+VJ**2)+4D0*VI*VJ)*XZZPL
     &                 -((1D0+VI**2)*(1D0+VJ**2)-4D0*VI*VJ)*XZZMI)
      XZZI =-SI*SJ/32D0/R*(Q2+AMZ2)*(VI-1D0)*
     &         ((1D0+VJ)**2*XZZMI-(1D0-VJ)**2*XZZPL)       -XZZRO
      XZZJ =-SI*SJ/32D0/R*(Q2+AMZ2)*(VJ-1D0)*
     &         ((1D0+VI)**2*XZZMI-(1D0-VI)**2*XZZPL)       -XZZRO
      XZZIJ=-SI*SJ/16D0/R*(Q2+AMZ2)*(VI-1D0)*(VJ-1D0)*XZZPL-XZZRO
      ENDIF
*
      CONST=DREAL(XZM1AL)
     *     -W0F-5.D0/4.D0-5.D0/8.D0/R+RW/8.D0
     *     +3.D0/4.D0*((1.D0/R1-1.D0/R)*ALR-RW*ALRW/RW1)
      XROK(1)=1.D0+AL4PI/R1*(CONST+XDBAL+2.D0*R*XV2BW+XROBZ-XROBT
     *       +(-2.D0*R+0.5D0+(VI+VJ)/4.D0)*XV1BW
     *       +(1.D0+3.D0/2.D0*(VI**2+VJ**2))/8.D0/R*XV1BZ)
     *       +XRQCD
      XROK(1)=XROK(1) + IBOXF*AL4PI/R1*(XZZRO+XWWRO)
*-------------------
      GAUGE=AL4PI/R1*(-ALR*(41D0/6D0-11D0/3D0*R)+2D0/3D0*R1)
      XROKF=1D0+AL4PI/R1*(-XWZ1AL)+GAUGE+XKQCD
*-------------------
*  XROK(2)=KAPPA WITH INDEX I
      WZ1AL=DREAL(XWZ1AL)
      XROK(2)=1.D0+AL4PI/R1*(-WZ1AL+XFMF-R*XA1BW-2.D0/3.D0*R
     *       +43.D0/18.D0-3.D0/4.D0*XRFL3-(2.D0*R+AMW2/Q2)*XV2BW
     *       -XROBZ+(2.D0*R-QJM-(VI+VJ)/4.D0+(1.D0-QJM)*AMW2/Q2)*XV1BW
     *       +(-VI*(1.D0+VI)/8.D0/R-QJM*VJ/2.D0*(1.D0+AMZ2/Q2))*XV1BZ)
     *       +XKQCD
      XROK(2)=XROK(2) + IBOXF*AL4PI/R1*(XZZI-XWWRO)
*
*  XROK(3)=KAPPA WITH INDEX J
      XROK(3)=1.D0+AL4PI/R1*(-WZ1AL+XFMF-R*XA1BW-2.D0/3.D0*R
     *       +43.D0/18.D0-3.D0/4.D0*XRFL3-(2.D0*R+AMW2/Q2)*XV2BW
     *       -XROBZ+XROBT
     *       +(2.D0*R-QIM-(VI+VJ)/4.D0+(1.D0-QIM)*AMW2/Q2)*XV1BW
     *       +(-VJ*(1.D0+VJ)/8.D0/R-QIM*VI/2.D0*(1.D0+AMZ2/Q2))*XV1BZ)
     *       +XKQCD
      XROK(3)=XROK(3) + IBOXF*AL4PI/R1*(XZZJ-XWWRO)
*
*  XROK(4)=KAPPA WITH INDEXES I AND J
      XROK(4)=1.D0+AL4PI/R1*(-2.D0*WZ1AL+2.D0*XFMF+(-R+AMW2/Q2)*XA1BW
     *       -4.D0/3.D0*R+35.D0/18.D0-2.D0/3.D0*XRFL3-2.D0*R*XV2BW
     *       -XROBZ+XROBT+(2.D0*R-0.5D0-(VI+VJ)/4.D0)*XV1BW
     *       +(-1.D0/8.D0/R-3.D0*(VI**2+VJ**2)/16.D0/R
     *       +(QI**2+QJ**2)*R1/R*(1.D0+AMW2/Q2))*XV1BZ)
     *       +2.D0*XKQCD
      XROK(4)=XROK(4) + IBOXF*AL4PI/R1*(XZZIJ-XWWRO)
*
C-----------------------------------------------------------------------
      IF(IBARB.EQ.0.OR.IBARB.EQ.-1) THEN
       AMT4C=19-2D0*PI2
        ELSEIF(IBARB.EQ.1) THEN
       RBTH=AMT2/AMH2
       ALRB=LOG(RBTH)
       AMT4C=49D0/4D0+PI2+27D0/2D0*ALRB+3D0/2D0*ALRB**2
     &      +RBTH/3D0*(2D0-12D0*PI2+12D0*ALRB-27D0*ALRB**2)
     &  +RBTH**2/48D0*(1613-240*PI2-1500*ALRB-720 *ALRB**2)
        ELSEIF(IBARB.EQ.2) THEN
       RBARB=SQRT(AMH2/AMT2)
       AMT4C=FBARB(RBARB)
      ENDIF
C--------------------------------------------------------------------
      RENORM=SQRT(2D0)*GMU*AMZ2*R1*R/PI*ALFAI
      IF(ISCRE.EQ.0) SCALER=1.00D0
      IF(ISCRE.GE.1) SCALER=RENORM
C--------------------------------------------------------------------
      IF (IAMT4 .EQ. 1 ) THEN
       DRHOT = .75D0*AL4PI/SW2/R*AMT2/AMZ2
       TOPX2 = GMU*AMT2/DSQRT(2.D0)/8.D0/PI2
       DRHOT4=3D0*TOPX2*(1D0+TOPX2*AMT4C)
       IF(IQCD.EQ.0) THEN
        TBQCD0=0D0
        TBQCD3=0D0
        CORRXI=0D0
         ELSE
          IF(IAFMT.EQ.0) THEN
* Below is the leading term as it is given in \zf (3.71)
           TBQCD0=.75D0*AL4PI/R/SW2*AMT2/AMZ2
     &           *(-CALXI/PI*2D0/3D0*(1D0+PI2/3D0))
           TBQCDL=.75D0*AL4PI/R/SW2*AMT2/AMZ2
     &           *(-CALST/PI*2D0/3D0*(1D0+PI2/3D0))
           CORRXI=-TBQCDL+TBQCD0
           TBQCD3=0D0
            ELSE
* This coincides with TBQCD0 above if O(\alpha\alpha^2_s)=0
           TBQCD0=AFMT3(CALST,AMT2,AMZ2,SW2)
* Here TBQCDR has to be called!!!
* Below is pure AFMT-correction to \Delta \rho (see p.18 PCWG-notebook)
           TBQCD3=.75D0*AL4PI/R/SW2*AMT2/AMZ2
     &          *(TBQCD0-(-CALST/PI*2D0/3D0*(1D0+PI2/3D0)))
           CORRXI=0D0
          ENDIF
       ENDIF
*
       ROKR(1)=DREAL(XROK(1))-1-        DRHOT+         TBQCD3+CORRXI
       ROKR(2)=DREAL(XROK(2))-1-  R/SW2*DRHOT+  R/SW2*(TBQCD3+CORRXI)
       ROKR(3)=DREAL(XROK(3))-1-  R/SW2*DRHOT+  R/SW2*(TBQCD3+CORRXI)
       ROKR(4)=DREAL(XROK(4))-1-2*R/SW2*DRHOT+2*R/SW2*(TBQCD3+CORRXI)
       AROK(1)=DIMAG(XROK(1))
       AROK(2)=DIMAG(XROK(2))
       AROK(3)=DIMAG(XROK(3))
       AROK(4)=DIMAG(XROK(4))
       ROKL(1)= 1/(1-DRHOT4)
       ROKL(2)=(1+R/SW2*DRHOT4)
       ROKL(3)=(1+R/SW2*DRHOT4)
       ROKL(4)=(1+R/SW2*DRHOT4)**2
C      XROK(1)= DCMPLX(DREAL(XROK(1)-DRHOT)/
C    &         (1D0-DRHOT4),DIMAG(XROK(1)))
C      XROKF  = DCMPLX(DREAL(XROKF  -R/SW2*DRHOT)*
C    &         (1D0+R/SW2*DRHOT4),DIMAG(XROK(2)))
C      XROK(2)= DCMPLX(DREAL(XROK(2)-R/SW2*DRHOT)*
C    &         (1D0+R/SW2*DRHOT4),DIMAG(XROK(2)))
C      XROK(3)= DCMPLX(DREAL(XROK(3)-R/SW2*DRHOT)*
C    &         (1D0+R/SW2*DRHOT4),DIMAG(XROK(3)))
C      XROK(4)= DCMPLX(DREAL(XROK(4)-2D0*R/SW2*DRHOT)*
C    &         (1D0+R/SW2*DRHOT4)**2,DIMAG(XROK(4)))
      ELSEIF(IAMT4 .EQ. 2 ) THEN
       DRHOT = .75D0*AL4PI/SW2/R*AMT2/AMZ2
       TOPX2 = GMU*AMT2/DSQRT(2.D0)/8.D0/PI2
       DRHOT4=3D0*TOPX2*(1D0+TOPX2*AMT4C)
       IF(IQCD.EQ.0) THEN
        TBQCD0= 0D0
        TBQCDL= 0D0
         ELSE
          IF(IAFMT.EQ.0) THEN
           TBQCD0=-TOPX2*CALXI/PI*2D0*(1D0+PI2/3D0)
            ELSE
           TBQCD0=3*TOPX2*AFMT3(CALST,AMT2,AMZ2,SW2)
* Here TBQCDR has to be called!!!
          ENDIF
        TBQCDL=-AL4PI*ALST/PI*AMT2/AMW2/R1*(.5D0+PI2/6D0)
       ENDIF
       ROKR(1)=DREAL(XROK(1))-1-         DRHOT-TBQCDL
       ROKR(2)=DREAL(XROK(2))-1-  R/SW2*(DRHOT+TBQCDL)
       ROKR(3)=DREAL(XROK(3))-1-  R/SW2*(DRHOT+TBQCDL)
       ROKR(4)=DREAL(XROK(4))-1-2*R/SW2*(DRHOT+TBQCDL)
       AROK(1)=DIMAG(XROK(1))
       AROK(2)=DIMAG(XROK(2))
       AROK(3)=DIMAG(XROK(3))
       AROK(4)=DIMAG(XROK(4))
       ROKL(1)=1/(1-DRHOT4-TBQCD0)
       ROKL(2)=(1+R/SW2*(DRHOT4+TBQCD0))
       ROKL(3)=(1+R/SW2*(DRHOT4+TBQCD0))
       ROKL(4)=(1+R/SW2*(DRHOT4+TBQCD0))**2
C      XROK(1) = DCMPLX(DREAL(XROK(1)-DRHOT-TBQCDL)/
C    &          (1.D0-DRHOT4-TBQCD0),DIMAG(XROK(1)))
C      XROK(2) = DCMPLX(DREAL(XROK(2)-R/SW2*(DRHOT+TBQCDL))*
C    &          (1.D0+R/SW2*(DRHOT4+TBQCD0)),DIMAG(XROK(2)))
C      XROKF   = DCMPLX(DREAL(XROKF  -R/SW2*(DRHOT+TBQCDL))*
C    &          (1.D0+R/SW2*(DRHOT4+TBQCD0)),DIMAG(XROK(2)))
C      XROK(3) = DCMPLX(DREAL(XROK(3)-R/SW2*(DRHOT+TBQCDL))*
C    &          (1.D0+R/SW2*(DRHOT4+TBQCD0)),DIMAG(XROK(3)))
C      XROK(4) = DCMPLX(DREAL(XROK(4)-2D0*R/SW2*(DRHOT+TBQCDL))*
C    &          (1.D0+R/SW2*(DRHOT4+TBQCD0))**2,DIMAG(XROK(4)))
*
* AMT4=3 option should work for bb-channel at AMT4=4, since 
*        Degrassi's corrections are not known for bb-channel       
*
      ELSEIF(IAMT4.EQ.3.OR.(IAMT4.GE.4.AND.IBFLA.EQ.1)) THEN
*
       DWZ1AL=R/R1*DREAL(XWZ1R1+XDWZ1F)
       RENORM=SQRT(2D0)*GMU*AMZ2*R1*R/PI*ALFAI
       SCALE = AL4PI/R1*(41D0/6D0-11D0/3D0*R)*ALR
       CORKAP=(AL4PI*DWZ1AL+SCALE)+.75D0*AL4PI/SW2**2*AMT2/AMZ2
*
       DRHOT = .75D0*AL4PI/SW2/R*AMT2/AMZ2
       TOPX2 = GMU*AMT2/DSQRT(2.D0)/8.D0/PI2
       DRHOT4=3D0*TOPX2*(1D0+TOPX2*AMT4C)
       IF(IQCD.EQ.0) THEN
        TBQCD0= 0D0
        TBQCDL= 0D0
         ELSE
          IF(IAFMT.EQ.0) THEN
           TBQCD0=-TOPX2*CALXI/PI*2D0*(1D0+PI2/3D0)
            ELSE
           TBQCD0=3*TOPX2*AFMT3(CALST,AMT2,AMZ2,SW2)
* Here TBQCDR has to be called!!!
          ENDIF
        TBQCDL=-AL4PI*ALST/PI*AMT2/AMW2/R1*(.5D0+PI2/6D0)
       ENDIF
       ROKR(1)=DREAL(XROK(1))-1-         DRHOT-TBQCDL
       ROKR(2)=DREAL(XROK(2))-1-  R/SW2*(DRHOT+TBQCDL)+  CORKAP
       ROKR(3)=DREAL(XROK(3))-1-  R/SW2*(DRHOT+TBQCDL)+  CORKAP
       ROKR(4)=DREAL(XROK(4))-1-2*R/SW2*(DRHOT+TBQCDL)+2*CORKAP
       AROK(1)=DIMAG(XROK(1))
       AROK(2)=DIMAG(XROK(2))
       AROK(3)=DIMAG(XROK(3))
       AROK(4)=DIMAG(XROK(4))
       ROKL(1)=1/(1-DRHOT4-TBQCD0)
       ROKL(2)=(1+R/SW2*(DRHOT4+TBQCD0)-CORKAP*RENORM)
       ROKL(3)=(1+R/SW2*(DRHOT4+TBQCD0)-CORKAP*RENORM)
       ROKL(4)=(1+R/SW2*(DRHOT4+TBQCD0)-CORKAP*RENORM)**2
C      XROK(1)=DCMPLX(DREAL(XROK(1)-DRHOT-TBQCDL)/
C    &        (1D0-DRHOT4-TBQCD0),DIMAG(XROK(1)))
C      XROK(2)=DCMPLX(DREAL(XROK(2)-R/SW2*(DRHOT+TBQCDL)+CORKAP)*
C    &        (1D0+R/SW2*(DRHOT4+TBQCD0)-CORKAP*RENORM),DIMAG(XROK(2)))
C      XROKF  =DCMPLX(DREAL(XROKF-R/SW2*(DRHOT+TBQCDL)+CORKAP)*
C    &     (1D0+R/SW2*(DRHOT4+TBQCD0)-CORKAP*RENORM),DIMAG(XROK(2)))
C      XROK(3)=DCMPLX(DREAL(XROK(3)-R/SW2*(DRHOT+TBQCDL)+CORKAP)*
C    &        (1D0+R/SW2*(DRHOT4+TBQCD0)-CORKAP*RENORM),DIMAG(XROK(3)))
C      XROK(4)=DCMPLX(DREAL(XROK(4)-2D0*R/SW2*(DRHOT+TBQCDL)+2D0*CORKAP)
C    &     *(1D0+R/SW2*(DRHOT4+TBQCD0)-CORKAP*RENORM)**2,DIMAG(XROK(4)))
      ENDIF
*-----------------------------------------------------------------------
*
* bb-channel modification
*
      IF(IBFLA.EQ.1) THEN
************************************************************************
* Obsolete comment:
*     APPROX. CORRECTION FOR FINITE T-MASS IN B CHANNEL
*     I.E. : NO T-QUARK MASS IN BOXES AND IN PHOTON VERTICES
*          THE T- QUARK MASS EFFECT IN THE Z VERTICES HAS BEEN
*          CALCULATED PRELIMINARY FOR S = MZ**2 ONLY.THIS IS
*          ACCURATE AT LEP I.
*     WE ASSUME THAT MT IS LARGER THAN SQRT(S)/2
*     COMMENT : COMMON CDZVZW IS FILLED IN Z- DECAY CHAIN
*          THE VALUE VTB IS GIVEN BY SR F1ZBT
************************************************************************
C 13/10/1992 - Barbieri's m_t^4 are implemented
      IF(IBARB.EQ.0.OR.IBARB.EQ.-1) THEN
       AMT4B=(27-PI2)/3
        ELSEIF(IBARB.EQ.1) THEN
       RBTH=AMT2/AMH2
       ALRB=LOG(RBTH)
       AMT4B=1D0/144*(311D0+24*PI2+282*ALRB+90*ALRB**2
     &      -4D0*RBTH*(40D0+ 6*PI2+ 15*ALRB+18*ALRB**2)
     &      +3D0*RBTH**2*(242.09D0-60*PI2-454.2D0*ALRB-180*ALRB**2))
        ELSEIF(IBARB.EQ.2) THEN
       RBARB=SQRT(AMH2/AMT2)
       AMT4B=FBARBB(RBARB)
      ENDIF
      TOPX2 = GMU*AMT2/DSQRT(2.D0)/8.D0/PI2
      IF(IFTJR.EQ.1) THEN
        TAUBB1=-2*TOPX2*(1-PI/3*CALST)
      ELSE
        TAUBB1=-2*TOPX2
      ENDIF
      TAUBB2=-2*TOPX2*TOPX2*AMT4B*V_TB**2
*                       
* New, s-dependent bb-corrections
*
      IF(IMTADD.EQ.1) THEN
        CALL VTBANA(0,AMZ2,WWv2,WWv11,WWv12)
      ELSE
        CALL VTBANA(0,SB,WWv2,WWv11,WWv12)
      ENDIF
*
* WW-box via XBOX
*
      QBM=1D0/3D0
      VTBX1=R*WWv2-.5D0*(1D0-2D0*R1*(1-QBM))*WWv11-.5D0*WWv12
      VTBX2=R*(AMZ2/SB-1D0)*(-(1D0-QBM)*WWv11+WWv2)
*
* VTBX1 is asymptotically -AMT2/AMW2/2
*
      DVTBB1=(AL4PI*VTBX1/R1-AL1PI/4/SW2*(-AMT2/AMW2/2))*V_TB**2
      DVTBB2= AL4PI*VTBX2/R1*V_TB**2
*
      ROKR(1)=ROKR(1)+DVTBB1
      ROKR(2)=ROKR(2)+DVTBB2
      ROKR(3)=ROKR(3)-DVTBB1
      ROKR(4)=ROKR(4)-DVTBB1
*
      ROKL(1)=ROKL(1)*(1+(TAUBB1+TAUBB2)*IBFLA)
* 25/11/1998, confirm correctness of commenting, db
*     ROKL(2)=ROKL(2)/(1+(TAUBB1+TAUBB2)*IBFLA)
      ROKL(3)=ROKL(3)/(1+(TAUBB1+TAUBB2)*IBFLA)
      ROKL(4)=ROKL(4)/(1+(TAUBB1+TAUBB2)*IBFLA)
*
* End of bb-channel modification
*
      ENDIF
*
* Game with options for AMT4 < 4 only
*
      IF(IAMT4.LT.4.OR.IBFLA.EQ.1) THEN
      IF(IFACT.EQ.0) THEN
       XROKIN(1)=DCMPLX(1/(1/ROKL(1)-SCALER*ROKR(1)),AROK(1))
       XROKIN(2)=DCMPLX(ROKL(2)*(1+SCALER*ROKR(2)),AROK(2))
       XROKIN(3)=DCMPLX(ROKL(3)*(1+SCALER*ROKR(3)),AROK(3))
       XROKIN(4)=DCMPLX(ROKL(4)*(1+SCALER*ROKR(4)),AROK(4))
      ELSEIF(IFACT.EQ.1) THEN
       XROKIN(1)=DCMPLX(ROKL(1)*(1+SCALER*ROKR(1)*ROKL(1)),AROK(1))
       XROKIN(2)=DCMPLX(ROKL(2)+SCALER*ROKR(2),AROK(2))
       XROKIN(3)=DCMPLX(ROKL(3)+SCALER*ROKR(3),AROK(3))
       XROKIN(4)=DCMPLX(ROKL(4)+SCALER*ROKR(4),AROK(4))
      ELSEIF(IFACT.EQ.2) THEN
       XROKIN(1)=DCMPLX(ROKL(1)*(1+SCALER*ROKR(1)),AROK(1))
       XROKIN(2)=DCMPLX(ROKL(2)+SCALER*ROKR(2),AROK(2))
       XROKIN(3)=DCMPLX(ROKL(3)+SCALER*ROKR(3),AROK(3))
       XROKIN(4)=DCMPLX(ROKL(4)+SCALER*ROKR(4),AROK(4))
      ELSEIF(IFACT.EQ.3) THEN
       XROKIN(1)=DCMPLX(ROKL(1)+SCALER*ROKR(1),AROK(1))
       XROKIN(2)=DCMPLX(ROKL(2)+SCALER*ROKR(2),AROK(2))
       XROKIN(3)=DCMPLX(ROKL(3)+SCALER*ROKR(3),AROK(3))
       XROKIN(4)=DCMPLX(ROKL(4)+SCALER*ROKR(4),AROK(4))
      ELSE
       XROKIN(1)=DCMPLX(ROKL(1),AROK(1))
       XROKIN(2)=DCMPLX(ROKL(2),AROK(2))
       XROKIN(3)=DCMPLX(ROKL(3),AROK(3))
       XROKIN(4)=DCMPLX(ROKL(4),AROK(4))
         ROKR(1)=SCALER*ROKR(1)
         ROKR(2)=SCALER*ROKR(2)
         ROKR(3)=SCALER*ROKR(3)
         ROKR(4)=SCALER*ROKR(4)
      ENDIF
      ENDIF
*
* End of options
*
*-----------------------------------------------------------------------
*
* New option AMT4=4,5
*
*	work-around to include factorizable higher-order inital state
*	 corrections also for IBFLA>0 (bb final state):
*	 [A. Freitas, Oct 1, 2004]
*
      IF(IAMT4.GE.4) THEN
*
       RENORM=SQRT(2D0)*GMU*AMZ2*R1*R/PI*ALFAI
* the same scale as in SEARCH
       SCALEB = AL4PI/R1*(1D0/6D0+7D0*R)*ALR
       ANUMF  = 24D0
       TRQ2F  =  8D0
       SCALEF = -AL4PI/R1*(ANUMF/6-4D0/3*R1*TRQ2F)*ALR
       SCALE=SCALEB+SCALEF       
       CORKAP=AL4PI*R/R1*(XDWZ1F+XWZ1R1)+SCALE
       CORRHO=R1/R*CORKAP
       DRHOT = .75D0*AL4PI/SW2/R*AMT2/AMZ2
       TOPX2 = GMU*AMT2/DSQRT(2.D0)/8.D0/PI2
*
       IF(IQCD.EQ.0) THEN
        TBQCD0= 0D0
        TBQCDL= 0D0
         ELSE
          IF(IAFMT.EQ.0) THEN
           TBQCD0=-TOPX2*CALXI/PI*2D0*(1D0+PI2/3D0)
            ELSE
           TBQCD0=3*TOPX2*AFMT3(CALST,AMT2,AMZ2,SW2)
* Here TBQCDR has to be called!!!
          ENDIF
        TBQCDL=-AL4PI*ALST/PI*AMT2/AMW2/R1*(.5D0+PI2/6D0)
       ENDIF
*
       AMW=SQRT(AMW2)
       AMT=AMQ(5)
       PI3QE=ABS(QI)
       PI3QF=ABS(QJ)
*    
       CALL GDEGNL
     & (GMU,AMZ,AMT,AMH,AMW,PI3QE,Q2M,DRDREM,DRHOD,DKDREE,DROREE)
       CALL GDEGNL
     & (GMU,AMZ,AMT,AMH,AMW,PI3QF,Q2M,DRDREM,DRHOD,DKDREF,DROREF)
*    
       DF1BAR=RENORM*DRREMD
       DROREM=.5D0*(DROREE+DROREF)
*
cb     ROFACR=ROFACI      -TBQCDL+CORRHO-1
cb     AKFACR=AKFACI-R/SW2*TBQCDL+CORKAP-1
cb     ROFAC=(1D0+RENORM*ROFACR+DROREM)/(1D0      -DROBAR*(1D0-DF1BAR))
cb     AKFAC=(1D0+RENORM*AKFACR+DKDREM)*(1D0+R/SW2*DROBAR*(1D0-DF1BAR))
*
*	work-around to include factorizable higher-order inital state
*	 corrections also for IBFLA>0 (bb final state):
*	 [A. Freitas, Oct 1, 2004]
*
       ROKR(1)=1D0+DREAL(RENORM*(XROK(1)-TBQCDL+CORRHO-1D0)+DROREM)
       IF(IBFLA.EQ.1) THEN
        ROKL(1)=ROKL(1)
     &		*SQRT((1-DRHOT4-TBQCD0)/(1D0-DROBAR*(1D0-DF1BAR)))
*	correction to work-around of Oct 2004
*	 this part for a factorized treatment of the rho form factor
*	 was still missing for IBFLA>0 (bb final state):
*	 [A. Freitas, March 2005]
        ROKR(1)=1D0+DREAL(SQRT(RENORM)*(XROK(1)-TBQCDL-DRHOT/2D0
     &   +CORRHO/2D0-1D0) +DROREM/2D0)
	ROKR(1)=ROKR(1)+DVTBB1
       ELSE
        ROKL(1)=1D0/(1D0-DROBAR*(1D0-DF1BAR))
       ENDIF
       XROKIN(1)=DCMPLX(ROKL(1)*ROKR(1),DIMAG(XROK(1)))
*
       ROKL(2)=1D0+R/SW2*DROBAR*(1D0-DF1BAR)
       ROKR(2)=1D0
     &        +DREAL(RENORM*(XROK(2)-R/SW2*TBQCDL+CORKAP-1D0)+DKDREE)
       IF(IBFLA.EQ.1) THEN
        ROKR(2)=ROKR(2)+DVTBB2
       ENDIF
       XROKIN(2)=DCMPLX(ROKL(2)*ROKR(2),DIMAG(XROK(2)))
*
       IF(IBFLA.EQ.0) THEN
        ROKL(3)=1D0+R/SW2*DROBAR*(1D0-DF1BAR)
        ROKR(3)=1D0
     &        +DREAL(RENORM*(XROK(3)-R/SW2*TBQCDL+CORKAP-1D0)+DKDREF)
        XROKIN(3)=DCMPLX(ROKL(3)*ROKR(3),DIMAG(XROK(3)))
       ENDIF
*
*       ROKL(4)=(1D0+R/SW2*DROBAR*(1D0-DF1BAR))**2
       ROKL(4)=ROKL(2)*ROKL(3)
       ROKR(4)=1D0
     &        +DREAL(RENORM*(XROK(4)-2D0*R/SW2*TBQCDL+2D0*CORKAP-1D0)
     &        +DKDREE+DKDREF)
       IF(IBFLA.EQ.1) THEN
	ROKR(4)=ROKR(4)
     &        -DREAL(RENORM*(XROK(3)-R/SW2*TBQCDL+CORKAP-1D0)+DKDREF)
     &        +SCALER*ROKR(3)
       ENDIF	
       XROKIN(4)=DCMPLX(ROKL(4)*ROKR(4),DIMAG(XROK(4)))
        
*
      ENDIF
*
*	work-around to include new fermionic 2-loop corrections to kappa
*	into ROKANC: [A. Freitas, Oct 1, 2004]
*	The lept. eff. weak mixing angle has been calculated before by ZWRATE
*	and stored into ARSEFZ(2). The flavour dependence is taken from
*	the old calculation for AMT4=4
*
      IF(IAMT4.GE.6) THEN
       CALL VERTZW(1,1,1D0)
       CALL
     &  FOKAPP(1D0,DKDREF,DRHOD,DRDREM,SCALER,DROREE,DKDREE,AR1IM,AK1IM)
       DKADD=ARSEFZ(2)/SW2M-DKDREE
       IF(IBFLA.EQ.0) THEN
        XROKIN(2)=XROKIN(2)+DKADD
        XROKIN(3)=XROKIN(3)+DKADD
        XROKIN(4)=XROKIN(4)+2*DKADD
       ELSE
        XROKIN(4)=XROKIN(4)+DKADD
        XROKIN(2)=XROKIN(2)+DKADD
       ENDIF
      ENDIF
*
      DO IRK=1,4
        XROK(IRK)=XROKIN(IRK)
      ENDDO
*
*-----------------------------------------------------------------------
* 
* PHOTON FORMFACTOR
*
      IF(IALE2.EQ.0) THEN
        XFOT =1D0+AL4PI*XFOTF3(IALEM,    1,IHVP,IQCD,1,DAL5H,Q2)
        XFOT5=1D0+AL4PI*XFOTF3(IALEM,    1,IHVP,   0,0,DAL5H,Q2)
      ELSE
        XFOT =1D0+AL4PI*XFOTF3(IALEM,IALE2,IHVP,IQCD,1,DAL5H,Q2)
        XFOT5=1D0+AL4PI*XFOTF3(IALEM,IALE2,IHVP,   0,0,DAL5H,Q2)
      ENDIF
*
* Release for running quantities inside XFOTF3
*
      END

      SUBROUTINE J3WANA(MT2,MW2,AMI2,J3W)
*
      IMPLICIT NONE
      REAL*8 MT2,MW2,AMI2
      COMPLEX*16 AMT2,AMW2,J3W,DCMPLX,XSPENZ
*
      AMT2=DCMPLX(MT2,-1D-10)
      AMW2=DCMPLX(MW2,-1D-10)
*
      J3W=1D0/AMI2*(XSPENZ(1D0-AMT2/AMW2)
     &             -XSPENZ(1D0-AMT2/AMW2-AMI2/AMW2))
*
      RETURN
      END

      SUBROUTINE S3WANA(MT2,MW2,AMQ2,J0,S3,S30)
*
* Supplies real parts of J0,S3,S30 (w,t) in analytic presentation
* CALL S3ANA(MT2,MW2,-S,...) supplies `w' indices
* CALL S3ANA(MW2,MT2,-S,...) supplies `t' indices
*
      IMPLICIT NONE
      REAL*8 MT2,MW2,AMQ2
      REAL*8 PI,PI2,D2,D3,ALFAI,AL4PI,AL2PI,AL1PI
      COMPLEX*16 AMT2,AMW2
      COMPLEX*16 SQR,LQR,X0,X1,X2,X3
      COMPLEX*16 Y1,Y2,Y3,Y4,Y5,Y6,J0,S3,S30
      COMPLEX*16 XSPENZ,LOG,SQRT
c      COMPLEX*16 A1,A2,A3,A4,A5,A6
*
      COMMON/CDZCON/PI,PI2,D2,D3,ALFAI,AL4PI,AL2PI,AL1PI
*
      AMT2=DCMPLX(MT2,-1D-10)
      AMW2=DCMPLX(MW2,-1D-10)
*
      SQR=SQRT(1D0+4D0*AMW2/AMQ2)
*
      LQR=LOG((SQR+1D0)/(SQR-1D0))
      J0 =SQR*LQR
      IF(MT2.GT.MW2) THEN
        S30=1D0/AMQ2*LQR**2
      ELSE
        S30=1D0/AMQ2*(D2-XSPENZ(1D0-AMQ2/AMT2))
      ENDIF
*
      X1=(1D0-SQR)/2D0
      X2=(1D0+SQR)/2D0
      X0=(AMT2-AMW2)/AMQ2
      X3=AMT2/(AMT2-AMW2)
*
      Y1=X1/X0
      Y2=(1D0-X1)/(1D0-X0)
      Y3=X2/X0
      Y4=(1D0-X2)/(1D0-X0)
      Y5=X3/X0
      Y6=(1D0-X3)/(1D0-X0)
*
c      A1=1D0/(1D0-Y1)
c      A2=1D0/(1D0-Y2)
c      A3=1D0/(1D0-Y3)
c      A4=1D0/(1D0-Y4)
c      A5=1D0/(1D0-Y5)
c      A6=1D0/(1D0-Y6)
*
      S3=1D0/AMQ2*(XSPENZ(1D0/(1D0-Y1))-XSPENZ(1D0/(1D0-Y2))
     &            +XSPENZ(1D0/(1D0-Y3))-XSPENZ(1D0/(1D0-Y4))
     &            -XSPENZ(1D0/(1D0-Y5))+XSPENZ(1D0/(1D0-Y6)))
c      S3=1D0/AMQ2*(XSPENZ(A1)-XSPENZ(A2)
c     &            +XSPENZ(A3)-XSPENZ(A4)
c     &            -XSPENZ(A5)+XSPENZ(A6))
*
      RETURN
      END
 
      SUBROUTINE S4WANA(MT2,MW2,AMQ2,AMI2,S4W)
*
      IMPLICIT NONE
      REAL*8 MT2,MW2,AMQ2,AMI2
      COMPLEX*16 AMT2,AMW2,S4W
      COMPLEX*16 SQW,SQB,X1,X2,X3,X4,X1B,X2B
      COMPLEX*16 Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,Y11,Y12,Y13,Y14,Y15,Y16
      COMPLEX*16 DCMPLX,SQRT,XSPENZ
*
      AMT2=DCMPLX(MT2,-1D-10)
      AMW2=DCMPLX(MW2,-1D-10)
*
      SQW=SQRT(1D0+4D0*AMW2/AMQ2)
      SQB=SQRT(1D0+4D0*AMW2/AMQ2*AMI2*(AMT2+AMI2-AMW2)/(AMT2+AMI2)**2)
      X1 =1D0/2D0*(1D0-SQW)
      X2 =1D0/2D0*(1D0+SQW)
      X3 =1D0/(1D0-AMW2/AMT2)
      X4 =1D0/(1D0-AMW2/(AMT2+AMI2))
      X1B=X4/2D0*(1D0-SQB)
      X2B=X4/2D0*(1D0+SQB)
*
      Y1 =X1/X1B
      Y3 =X2/X1B
      Y5 =X4/X1B
      Y7 =X3/X1B
      Y9 =X1/X2B
      Y11=X2/X2B
      Y13=X4/X2B
      Y15=X3/X2B
*
      Y2 =(X1-1D0)/(X1B-1D0)
      Y4 =(X2-1D0)/(X1B-1D0)
      Y6 =(X4-1D0)/(X1B-1D0)
      Y8 =(X3-1D0)/(X1B-1D0)
      Y10=(X1-1D0)/(X2B-1D0)
      Y12=(X2-1D0)/(X2B-1D0)
      Y14=(X4-1D0)/(X2B-1D0)
      Y16=(X3-1D0)/(X2B-1D0)
*
      S4W=1D0/(AMQ2*(AMT2+AMI2)*SQB)*(
     &     +XSPENZ(1D0/(1D0-Y1 ))-XSPENZ(1D0/(1D0-Y2 ))
     &     +XSPENZ(1D0/(1D0-Y3 ))-XSPENZ(1D0/(1D0-Y4 ))
     &     +XSPENZ(1D0/(1D0-Y5 ))-XSPENZ(1D0/(1D0-Y6 ))
     &     -XSPENZ(1D0/(1D0-Y7 ))+XSPENZ(1D0/(1D0-Y8 ))
     &     -XSPENZ(1D0/(1D0-Y9 ))+XSPENZ(1D0/(1D0-Y10))
     &     -XSPENZ(1D0/(1D0-Y11))+XSPENZ(1D0/(1D0-Y12))
     &     -XSPENZ(1D0/(1D0-Y13))+XSPENZ(1D0/(1D0-Y14))
     &     +XSPENZ(1D0/(1D0-Y15))-XSPENZ(1D0/(1D0-Y16))
     &     )
*
      RETURN
      END
 
      SUBROUTINE RHOCC(S,Q2,U,QI,QJ,QK,QL,ROW)
* ------ FORMER NAME OF THE ROUTINE: ROWAL
* BEFORE USE OF RHOCC AT LEAST ONE CALL OF DIZET MUST BE DONE.
* SEE ALSO THE COMMENTS THERE.
*---------------------------------------------------------------------
* THIS ROUTINE CALCULATES THE WEAK CHARGED CURRENT FORM FACTOR FOR
* THE 4-FERMION SCATTERING CROSS SECTION. 
*----------------------------------------------------------------------
* EXAMPLES OF THE USE OF THIS ROUTINE MAY BE FOUND IN THE PACKAGE
* HECTOR, Comp. Phys. Commun. 94 (1996) 128 [hep-ph/9511434] 
*----------------------------------------------------------------------
* INPUT FROM USER:
*            S,Q2,U - THE KINEMATIC INVARIANTS FOR THE QUARK PROCESS
*                     (S+T-U=0)
*       QI,QJ,QK,QL - THE CHARGES OF THE FOUR FERMIONS IN THE PROCESS
* OUTPUT OF THE ROUTINE:
*               ROW - THE FORM FACTOR
*---------------------------------------------------------------------
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      COMMON/CDZFLG/IHVP,IAMT4,IQCD,IMOMS,IMASS,IALEM,IMASK,IBARB,IFTJR
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      COMMON/CDZWSC/SL2,SQ2,W0,W0F,Z0,Z0F,DWZ0R1,DWZ0F,XWM1,XWM1F,XZM1,
     &      XZM1F,XWZ1R1,XDWZ1F,XZFM1,XZFM1F,XAMM1,XAMM1F,XWFM1,XWFM1F
      COMMON/CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
      COMMON/CDZZWG/AMZ,AMH,GMU,A0,GAMZ,GAMW,CALSZ,CALST,CALXI,CALQED
*
      SI=1D0
      SJ=1D0
      SK=1D0
      SL=1D0
      QIM=ABS(QI)
       IF(QIM.NE.0)  SI=QI/QIM
      QJM=ABS(QJ)
       IF(QJM.NE.0)  SJ=QJ/QJM
      QKM=ABS(QK)
       IF(QKM.NE.0)  SK=QK/QKM
      QLM=ABS(QL)
       IF(QLM.NE.0)  SL=QL/QLM
      RDWB=DREAL(XDWB(Q2))
      RDWF=DREAL(XDWF(Q2))
      DBAL=RDWB+RDWF
      W0AL =W0+W0F
      WM1AL=DREAL(XWM1+XWM1F)
      QMIJ =QIM*QJM+QKM*QLM
      QMIK =QIM*QKM+QJM*QLM
      QMIL =QIM*QLM+QJM*QKM
      V1BZ =DREAL(XV1B(Q2,AMZ2))
      UBW  =UB(Q2,AMW2)
      V2BWZ=V2B(Q2,AMW2,AMZ2)
      ROBW =DREAL(XROBW(Q2,AMW2,AMZ2))
      BQSWZ=BF(Q2,S,AMW2,AMZ2)
      BQUWZ=BF(Q2,U,AMW2,AMZ2)
      AQUSWZ=AF(Q2,U,S,AMW2,AMZ2)
      AL4PI=1D0/ALFAI/PI/4D0
      ROW=1D0+AL4PI/R1*(DBAL-W0AL+WM1AL+5D0/8D0*R*(1D0+R)-11D0/2D0
     &-9D0/4D0*R*ALR/R1+(-1D0+1D0/2D0/R-R12/R*QMIJ)*V1BZ+2D0*R*V2BWZ
     &-R1*UBW+ROBW+(2D0-1D0/R+2D0*R12/R*QMIK)*S*(Q2+AMW2)*BQSWZ
     &+(2D0-1D0/R+2D0*R12/R*QMIL)*(Q2+AMW2)*(U*BQUWZ-AQUSWZ))
C****************************************************************
      SM=ABS(S)
      UM=ABS(U)
      Q2M =ABS(Q2)
      ALQW=LOG(Q2M/AMW2)
      ALSW=LOG(SM/AMW2)
      ALUW=LOG(UM/AMW2)
C****************************************************************
      PI2=PI**2
      FF1=PI2/6D0
      SW=S/AMW2
      UW=U/AMW2
      QW=Q2/AMW2
      ROWADD=AL4PI*(QMIJ*(4D0 -2D0*FF1-PI2*TET(-Q2))
     & -QMIK*((LOG(SM/Q2M))**2-2D0*FF1-PI2*TET(S))
     & -QMIL*((LOG(UM/Q2M))**2-2D0*FF1-PI2*TET(U))
     & -.5D0+2D0*(FF1-SPENCE(1D0+SW)-4D0*SPENCE(-QW)
     & -4D0*LOG(ABS(1D0+QW))*LOG(ABS(SW)))
     & -2D0*QMIL*(SPENCE(1D0+UW)-SPENCE(1D0+SW)+2D0*
     & LOG(ABS(1D0+QW))*LOG(ABS(U/S))+(Q2+AMW2)*AA00(Q2,U,S,AMW2)))
      ROW=ROW+ROWADD
*
      Q2M=-Q2
      ALST=CALXI
      SW2=R1
      AMT2=AMQ(5)**2
      XZERO=DCMPLX(0D0,0D0)
      XRQCD=XZERO
      IF(IQCD-1) 1,2,3
2     XRQCD=AL4PI*XCQCDS(ALST,SW2,AMT2,Q2M)
      GOTO 1
3     XRQCD=AL4PI*XRCQCD(ALST,SW2,AMT2,Q2M)
1     CONTINUE
      ROW=ROW+DREAL(XRQCD)
*
      END
 
      FUNCTION XCQCDS(ALST,SW2,AMT2,S)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      COMMON/CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
*
* CALPI/4 IS OMITTED
*
      DATA EPS/1.D-3/
      ALTW=-AMT2/AMW2
      ALTS=-AMT2/S
      SMW2=S/AMW2
      DMW2=1D0-SMW2
      XPWFTS=XPWFI(ALTS)
      XPWFTW=XPWFI(ALTW)
      IF(ABS(DMW2).LT.EPS) GO TO 1
      XCQCDS=ALST/(3D0*PI*SW2)*(
     *      +AMT2/4D0/AMW2*(1D0/DMW2*(XPWFTW-XPWFTS)-XPWFTW)
     *      -AMT2/AMW2*(PI2/2D0+105D0/8D0))
      RETURN
1     XDWFTW=XDPWFI(ALTW)
      XCQCDS=ALST/(3D0*PI*SW2)*(AMT2/4D0/AMW2*(
     *      +XDWFTW/ALTW-XPWFTW)-AMT2/AMW2*(PI2/2D0+105D0/8.D0))
*
      END
 
      FUNCTION XRCQCD(ALST,SW2,AMT2,S)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      COMMON/CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
*
* CALPI/4 IS OMITTED
*
      DATA EPS/1.D-3/
      ALTW=-AMT2/AMW2
      ALTS=-AMT2/S
      SMW2=S/AMW2
      DMW2=1D0-SMW2
      XPWFTS=XPWF(ALTS)
      XPWFTW=XPWF(ALTW)
      IF(ABS(DMW2).LT.EPS) GO TO 1
      XRCQCD=ALST/(3D0*PI*SW2)*
     &      (AMT2/4D0/AMW2*(1D0/DMW2*(XPWFTW-XPWFTS)-XPWFTW)
     *      -AMT2/AMW2*(PI2/2D0+105D0/8D0))
      RETURN
1     XDWFTW=XDPWF(ALTW)
      XRCQCD=ALST/(3D0*PI*SW2)*(AMT2/4D0/AMW2*(
     *      +XDWFTW/ALTW-XPWFTW)-AMT2/AMW2*(PI2/2D0+105D0/8D0))
*
      END
 
      FUNCTION XDPWF(ALTW)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
* DERIVATIVE OF XPDF  IS STILL MISSING 
*
      XDPWF=(0D0,0D0)
*
      END
 
      FUNCTION XDPWFI(ALTW)
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
* DERIVATIVE OF XPDFI IS STILL MISSING
*
      XDPWFI=(0D0,0D0)
*
      END
 
      FUNCTION AA00(Q2,U,S,AMW2)
      IMPLICIT REAL*8(A-H,O-Z)
*
      SW=S/AMW2
      UW=U/AMW2
      QW=Q2/AMW2
      AA00=1D0/S*(-LOG(ABS(UW))+(1D0+1D0/QW)*LOG(ABS(1D0+QW))
     & +(2D0-Q2/S-1D0/SW)*(SPENCE(1D0+QW)-SPENCE(1D0+UW)
     & +LOG(ABS(1D0+QW))*LOG(ABS(Q2/U))))
*
      END
 
      FUNCTION XK3(Q2,AM12,AM22)
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      PI=4D0*ATAN(1D0)
      PI2=PI**2
      Q2M=Q2+AM12+AM22
      Q2TR=Q2M+2D0*SQRT(AM12*AM22)
      ALAM=Q2M*Q2M-4D0*AM12*AM22
      XDL2=ALAM*(XJ(Q2,AM12,AM22))**2-(LOG(AM12/AM22))**2
      IF(Q2TR)1,1,2
1     XK3=.25D0*XDL2-PI2
      RETURN
2     XK3=.25D0*XDL2
*
      END
 
      FUNCTION UB(Q2,AMV2)
      IMPLICIT REAL*8(A-H,O-Z)
      DATA EPS/1D-3/
*
C  FUNCTION UB(Q2,AMW2) IS EQUAL TO 2.*UBAR+(1.+AIN)*LOG(ABS(1.+A))
      A=Q2/AMV2
      AIN=1D0/A
      IF(ABS(A).LT.EPS) GO TO 3
      UB=-43D0/6D0-13D0/6D0*AIN
     &+(-2D0/3D0+13D0/6D0*AIN)*(1D0+AIN)*LOG(ABS(1D0+A))
      RETURN
3     UB=-27D0/4D0-25D0/36D0*A
*
      END
 
      FUNCTION BF(Q2,X,AM12,AM22)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CDZFBF/CQ2,CX,CM12,CM22
      EXTERNAL BFIND
      DATA EPS/.001D0/
*
      RQ1=Q2/AM12
      RX1=X/AM12
      IF(ABS(RQ1).LT.EPS.AND.ABS(RX1).LT.EPS) GO TO 11
      CQ2=Q2
      CX=X
      CM12=AM12
      CM22=AM22
      CALL SIMPS(0D0,1D0,.1D0,EPS,1D-30,BFIND,Y,R1,R2,R3)
      BF=R1
      RETURN
11    AR=AM12/AM22
      AR1=1D0-AR
      ALR=LOG(AR)
      ALX=LOG(ABS(RX1))
      IF(AR1)12,13,12
12    BF=
     &(1D0-ALX+AR/AR1*ALR+RQ1/(AR1**2)*(-.5D0*AR*(1D0+AR)-AR**2/AR1*ALR)
     *+RX1*((1D0+AR)*(-.5D0+ALX)/2D0-.5D0*AR**2/AR1*ALR))/AM12/AM22
      RETURN
13    BF=(-ALX+RX1*ALX-RQ1/6D0)/(AM12**2)
*
      END
 
      FUNCTION BFIND(Y)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CDZFBF/Q2,X,AM12,AM22
      DATA EPS/1D-8/
*
      Y1=1D0-Y
      AMY2=Y*AM12+Y1*AM22
      AMY4=AMY2**2
      AKY2=Y*Y1*Q2+AMY2
      D=X*AKY2+AMY4
      R=X*AKY2/AMY4
      IF(ABS(R+1D0).GT.EPS) GO TO 1
      BFIND=1D0/AMY4
      RETURN
1     BFIND=-LOG(ABS(R))/D
*
      END
 
      FUNCTION AF(Q2,AX,AY,AM12,AM22)
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
*
      DATA EPS/1D-3/,EPSVAX/5D-39/
      F1=PI**2/6D0
      RQ1=Q2/AM12
      RX1=AX/AM12
      AR=AM12/AM22
      AR1=1D0-AR
      ALR=LOG(AR)
      IF(ABS(RQ1).LT.EPS.AND.ABS(RX1).LT.EPSVAX) GO TO 1
C--
      ALX=LOG(ABS(RX1))
      BFX=BF(Q2,AX,AM12,AM22)
      CFX=2D0*F1-SPENCE(1D0+AX/AM12)-SPENCE(1D0+AX/AM22)
     &+2D0*DREAL(XK3(Q2,AM12,AM22))+(Q2+AM12+AM22)*AX*BFX
      AF=(-ALX+(-1D0+AR1/AR/RQ1)/2D0*ALR+DREAL(XL(Q2,AM12,AM22))/Q2/2D0
     &-AM12*AM22*BFX+(1D0-(Q2+AM12+AM22)/AY/2D0)*CFX)/AY
      RETURN
1     IF(AR1)2,3,2
2     AF=(-3D0/2D0*ALR/AR1+RX1*(7D0/18D0 + 2D0/3D0*AR/AR1*ALR)+
     *5D0/6D0*RQ1*AR/(AR1**2)*(2D0+(1D0+AR)/AR1*ALR))/AM22
      RETURN
3     AF=(3D0/2D0-5D0/18D0*RX1-5D0/36D0*RQ1)/AM12
*
      END
 
      FUNCTION XROBW(Q2,AM12,AM22)
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
C  ROB(Q2,AM12,AM22)=-R*(Q2+AMV2)*OMEGA(Q2,AM12,AM22)-BAR
      COMMON/CDZWSM/AMW2,AMZ2,RC,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      DATA EPS/1D-3/
*
      R=AM12/AM22
      R1=1D0-R
      R12=R1**2
      R13=R1**3
      A1=Q2/AM12
      IF(ABS(A1).LT.EPS) GO TO 1
      XROBW=4D0/3D0*R**2+19D0/12D0*R-1D0/12D0-1D0/12D0*R12*AM12/Q2
     &+(3D0/8D0*R**2+R/3D0+1D0/24D0
     &+R1*(-R/3D0-11D0/24D0+1D0/24D0/R)*AM12/Q2+1D0/24D0/R*R13*
     &(AM12/Q2)**2)*LOG(R)+(-3D0/8D0*R**2-R/2D0+1D0/24D0+1D0/24D0*R12*
     &AM12/Q2)/Q2*XL(Q2,AM12,AM22)+4D0*(RC+R*AM12/Q2)*XK3(Q2,AM12,AM22)
      RETURN
1     IF(R1)2,3,2
2     XROBW=5D0/8D0*R*(1D0+R)
     &+A1*(-7D0/18D0*R**2+R/24D0+13D0/6D0*R**3/R12)
     &+(3D0*R-9D0/4D0*R/R1+A1*(-4D0*RC*R/R1+R**2/R1/12D0
     &+13D0/12D0*R**3*(1D0+R)/R13))*LOG(R)
      RETURN
3     XROBW=3.5D0+A1*(4D0*RC-11D0/18D0)
*
      END
 
      FUNCTION V2B(Q2,AM12,AM22)
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
      DATA EPS/1D-3/
*
      R=AM12/AM22
      R1=1D0-R
      R12=R1*R1
      R13=R12*R1
      RA=(1D0+10D0/R+1D0/R**2)*(1D0+R)/24D0
      AQ=AM12/Q2
      A1=Q2/AM12
      IF(ABS(A1).LT.EPS) GO TO 1
      V2B=-4D0/3D0*R+5D0/2D0-4D0/3D0/R+2*RA*AQ
     &+R1*(-3D0/8D0*R1/R+1D0/3D0*(1D0-1D0/2D0/R+1D0/R**2)*AQ-RA/R*AQ**2)
     &*LOG(R)+(3D0/8D0*R+5D0/12D0+3D0/8D0/R-RA*AQ)/Q2
     &*DREAL(XL(Q2,AM12,AM22))
     &+2D0/R*(-(1D0+R)*AQ+AQ**2)*DREAL(XK3(Q2,AM12,AM22))
      RETURN
1     IF(R1)2,3,2
2     V2B=-5D0/8D0*R+11D0/4D0-5D0/8D0/R+A1*(1D0+R)/R12/18D0*(7D0*R**2
     &-23D0*R+7D0)+(-3D0/4D0/R+3D0/2D0/R1-R*(R**2+4D0*R+1D0)/R13/6D0*A1)
     &*LOG(R)
      RETURN
3     V2B=7D0/9D0*A1
*
      END
 
      FUNCTION XDWB(Q2)
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
      COMMON/CDZWSM/AMW2,AMZ2,R,R1,R12,R2,AMH2,RW,RW1,RW12,RW2,RZ,RZ1,
     *      RZ12,RZ2,ALR,ALRW,ALRZ,SW2M,CW2M,AKSX,R1W,R1W2
      COMMON/CDZWSC/SL2,SQ2,W0,W0F,Z0,Z0F,DWZ0R1,DWZ0F,XWM1,XWM1F,XZM1,
     &      XZM1F,XWZ1R1,XDWZ1F,XZFM1,XZFM1F,XAMM1,XAMM1F,XWFM1,XWFM1F
      DATA EPS/1D-3/
*
      QW=(Q2+AMW2)/AMW2
      IF(ABS(QW-1D0).LT.EPS) GO TO 3
      IF(ABS(QW).LT.EPS) GO TO 1
      XDRWH=(XL(Q2,AMW2,AMH2)-XL(-AMW2,AMW2,AMH2))/QW/AMW2
      XDRWZ=(XL(Q2,AMW2,AMZ2)-XL(-AMW2,AMW2,AMZ2))/QW/AMW2
      GO TO 2
1     XDRWH=2D0+QW/(RW-4D0)
     &         +(RW-2D0*QW/(RW-4D0))*AMW2*XJ(-AMW2,AMW2,AMH2)
      XDRWZ=2D0+QW*R/(1D0-4D0*R)
     &         +(1D0/R-2D0*QW*R/(1D0-4D0*R))*AMW2*XJ(-AMW2,AMW2,AMZ2)
2     CONTINUE
      AQ=AMW2/Q2
      XWH=-1D0/12D0*RW12*AQ+1D0/24D0*(1D0+RW1*(-10D0+5D0*RW-RW2)*AQ
     &+RW1*RW12*AQ**2)*ALRW+(-11D0+4D0*RW-RW2+RW12*AQ)/24D0/Q2
     &*XL(Q2,AMW2,AMH2)+(1D0/2D0-RW/6D0+RW2/24D0)*XDRWH
      XWL=(-3D0/8D0*R2+R+25D0/12D0-2D0/3D0/R-1D0/24D0/R2
     &+R12*(1D0+10D0/R+1D0/R2)*AQ/24D0)/Q2*XL(Q2,AMW2,AMZ2)
     &+(-2D0*R-17D0/6D0+2D0/3D0/R+1D0/24D0/R2)*XDRWZ
      XDWB=4D0/3D0*R2+7D0/12D0*R+253D0/36D0
     &+(-R2+2D0*R+8D0-8D0/R-1D0/R2)/12D0*AQ
     &+(3D0/8D0*R2+11D0/6D0*R+2D0/3D0
     &+R1*(-8D0*R+35D0+61D0/R-15D0/R2-1D0/R**3)/24D0*AQ
     &+R1*R12*(1D0/R+10D0/R2+1D0/R**3)/24D0*AQ**2)*ALR
     &+R1*(-2D0+17D0/6D0*AQ+5D0/6D0*AQ**2)*LOG(ABS(1D0+1D0/AQ))+XWL+XWH
      RETURN
3     XDWB=W0-XWM1+Q2/AMW2*(XWM1-7D0/18D0*R2-31D0/24D0*R+319D0/36D0
     &-5D0/8D0/R+RW/8D0+2D0*R2/R12-(RW/3D0+.5D0)/RW12
     &+(5D0/6D0*R-43D0/12D0-3D0/4D0/R+23D0/4D0/R1
     &+(1D0+R)*R2/R12/R1)*ALR-(7D0/4D0+5D0*RW2/6D0/RW1)*RW/RW12*ALRW)
*
      END
 
      FUNCTION FBARB(X)
      IMPLICIT REAL*8(A-Z)
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
      DATA P1/-0.74141D0/,P2/ -11.483D0  /,P3/  9.6577D0/,
     &     P4/ -6.7270D0/,P5/  3.0659D0  /,P6/-0.82053D0/,
     &     P7/ 0.11659D0/,P8/-0.67712D-02/
      IF(X.LE.4D0) THEN
        FBARB=P1+P2*X+P3*X**2+P4*X**3+P5*X**4+P6*X**5+P7*X**6+P8*X**7
         ELSE
        RBTH=1/X**2
        ALRB=LOG(RBTH)
        FBARB=49D0/4D0+PI2+27D0/2D0*ALRB+3D0/2D0*ALRB**2
     &       +RBTH/3D0*(2D0-12D0*PI2+12D0*ALRB-27D0*ALRB**2)
     &       +RBTH**2/48D0*(1613-240*PI2-1500*ALRB-720 *ALRB**2)
      ENDIF
*
      END
 
      FUNCTION FBARBB(X)
C 13/10/1992 - Barbieri's m_t^4 are implemented
      IMPLICIT REAL*8(A-Z)
      COMMON/CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
* Approximation from 0 to 4 (Mhiggs/mtop)
      DATA P1/ 5.6807D0/,P2/ -11.015D0  /,P3/ 12.814D0/,
     &     P4/-9.2954D0/,P5/  4.3305D0  /,P6/-1.2125D0/,
     &     P7/0.18402D0/,P8/-0.11582D-01/
      IF(X.LE.4D0) THEN
        FBARBB=P1+P2*X+P3*X**2+P4*X**3+P5*X**4+P6*X**5+P7*X**6+P8*X**7
         ELSE
        RBTH=1/X**2
        ALRB=LOG(RBTH)
        FBARBB=1D0/144*(311D0+24*PI2+282*ALRB+90*ALRB**2
     &        -4D0*RBTH*(40D0+ 6*PI2+ 15*ALRB+18*ALRB**2)
     &        +3D0*RBTH**2*(242.09D0-60*PI2-454.2D0*ALRB-180*ALRB**2))
      ENDIF
*
      END
 
      SUBROUTINE QCDCOF(SQS,AMT,SW2,ALQED,ALFAS,ALFAT,ALFAXI,QCDCOR)
*
* mod. 14 Jan 2013 new flag: IBAIKOV
* if IBAIKOV=2012, then arXiv:1201.5804v3 (02 May 2012) is implemented
* if IBAIKOV=2008, then arXiv:0801.1821v2 (04 Jul 2008) is implemented
* if IBAIKOV=2005, then status of QCD radiators is that of dizet6_42 and
* dizet6_43 as described in CPC133 and CPC174
*
      IMPLICIT REAL*8(A-H,J-Z)
* 
* mod. 14 Jan 2013
      INTEGER IBAIKOV
*
      COMMON/CDZSCT/ISCRE,ISCAL,IAFMT,IFACR,IFACT,IHIGS,IEWLC,ICZAK
     &             ,IHIG2,IALE2,IGFER      
      COMMON/CDZZWG/CAMZ,CAMH,GMU,A0,GAMZ,GAMW,CALSZ,CALST,CALXI,CALQED
      COMMON /CDZZER/ PI,ALPHAQ,ANF,AKC,AKB
      COMMON /CDZPHM/ CMASS,BMASS,TMASS,SMRUN1
      COMMON /CDZINP/ BMASSI
      COMMON /CDZRUN/ CMQRUN(8)
      COMMON /CDZBGC/ BETA03,BETA13,BETA23,GAMA03,GAMA13,GAMA23,
     &                BETA04,BETA14,BETA24,GAMA04,GAMA14,GAMA24,
     &                BETA05,BETA15,BETA25,GAMA05,GAMA15,GAMA25,
     &                COEF13,COEF23,COEF14,COEF24,COEF15,COEF25,
     &                ALMSB3,ALMSB4,ALMSB5,ALMSB6,ALF1PI,ALFCPI,ALFBPI,
     &                AMSAMC,AMSAMB,AMCAMC,AMCAMB,AMBAMB
      COMMON /CDZBG3/ BETA33,BETA34,BETA35,GAMA33,GAMA34,GAMA35
      COMMON /CDZFSR/ ISFSR
*
      DIMENSION QCDCOR(0:14),QCDCON(0:14)
      DIMENSION AMQRUN(6),CHARQU(6)
      DIMENSION QCDC2V(6),QCDC2A(6),QCDC4V(6),QCDC4A(6)
      DIMENSION SCTAVR(6),SCTPLU(6),SCTMIN(6)
*
      EXTERNAL ZALPHA,ZALSFL,ZRMCIN,ZRMCMC,ZRMBMB
*
      DATA AMQRUN/0D0,0D0,.750D0,0D0,0D0,2.40D0/
      DATA SCTAVR/.20D0,.19D0,.18D0,.17D0,.15D0,.13D0/
      DATA SCTPLU/.09D0,.10D0,.10D0,.10D0,.09D0,.08D0/
      DATA SCTMIN/.05D0,.04D0,.04D0,.04D0,.04D0,.04D0/
*
* HEAVY QUARK POLE MASSES (DB/GP-CONVENTION)
*
      DATA SMASS/0.3D0/,CMASS/1.5D0/,BMASS/4.7D0/
*
      ZQED=1D0
      ZMIX=1D0
      IF(ISFSR.EQ.0) ZQED=0D0
      IF(ISFSR.EQ.-1) THEN
        ZQED=0D0
        ZMIX=0D0
      ENDIF
*
* ISCAL=1,2,3 ARE NOT UP TO DATE, SEE KNIEHL's CONTRIBUTION:
* GAMES WITH THE SCALE FOR THE ALFA*ALFAS(SCALET*AMT) CORRECTIONS
*
      AMZ =CAMZ
      AMZ2=AMZ**2
      TMASS=AMT
*
* S-QUARK RUNNING MASS AT 1 GEV
*
      SMRUN1=.189D0
*
      SCALET=1D0
      IF(ISCAL.NE.0) THEN
       IF(AMT.LE.126D0                 ) INDSCT=1
       IF(AMT.GT.126D0.AND.AMT.LE.151D0) INDSCT=2
       IF(AMT.GT.151D0.AND.AMT.LE.176D0) INDSCT=3
       IF(AMT.GT.176D0.AND.AMT.LE.201D0) INDSCT=4
       IF(AMT.GT.201D0.AND.AMT.LE.226D0) INDSCT=5
       IF(AMT.GT.226D0                 ) INDSCT=6
      ENDIF
      IF(ISCAL.EQ.1    ) THEN
       SCALET=SCTAVR(INDSCT)+SCTPLU(INDSCT)
      ELSEIF(ISCAL.EQ.2) THEN
       SCALET=SCTAVR(INDSCT)
      ELSEIF(ISCAL.EQ.3) THEN
       SCALET=SCTAVR(INDSCT)-SCTMIN(INDSCT)
      ELSEIF(ISCAL.EQ.4) THEN
       SCALET=.204D0
      ENDIF
*
      ALFAT=0D0
      ALFAXI=0D0
*
      DO IQ=0,14
       QCDCOR(IQ)=1D0
      ENDDO
*
      DO IQ=1,3
       CHARQU(2*IQ-1)=2D0/3
       CHARQU(2*IQ  )=1D0/3
      ENDDO
*
      DO IMQ=1,6
       CMQRUN(IMQ)=AMQRUN(IMQ)
      ENDDO
      CMQRUN(3)=CMASS
      CMQRUN(6)=BMASS
*
      IF(ALFAS.LT.1D-10) RETURN
*
* NUMERICAL CONSTANTS
*
      PI=4D0*ATAN(1D0)
      PI2=PI**2
      D2=PI2/6D0
      D3=1.2020569031596D0
      D4=PI2**2/90
      D5=1.0369277551434D0
*
* Check of SQRT(S), its redefinition or STOP
*
      IF(SQS.LT.13D0) THEN
        PRINT *,'Warning: you have requested SQRT(S).LT.13 GeV.'
        PRINT *,'Program resets SQRT(S) to 13 GeV and continues.'
        SQS=13D0
      ENDIF
*
      IF(SQS.GE.350D0) THEN
        PRINT *,'You have requested SQRT(S).GE.350 GeV.'
        PRINT *,'Program STORs, since it is not foreseen'
        PRINT *,'to run above the t-tbar threshold.'
      ENDIF
*
* COEFFICIENTS OF BETA AND GAMA FUNCTIONS
*
      IF (SQS.LE.CMASS) THEN
       ANF=3D0
      ELSEIF(SQS.LE.BMASS) THEN
       ANF=4D0
      ELSEIF(SQS.LE.TMASS) THEN
       ANF=5D0
      ELSE
       ANF=6D0
      END IF
      ANF3=3D0
      ANF4=4D0
      ANF5=5D0
*
* First MG's insertion
*
      AMIN=LOG(AMZ2/1D-3)
      AMAX=LOG(AMZ2/1D+2)
      ALPHAQ=ALFAS            
      DO I=1,10
      FMIN=ZALPHA(AMIN)
      FMAX=ZALPHA(AMAX)
C     PRINT*,1,I,AMIN,AMAX,FMIN,FMAX
      IF ((FMIN.LE.0D0.AND.FMAX.LE.0D0).OR.
     &    (FMIN.GE.0D0.AND.FMAX.GE.0D0)) THEN
        AMIN=AMIN+2D0
        AMAX=AMAX-2D0
      ELSE
        GOTO 97
      ENDIF
      ENDDO
 97   CONTINUE
      CALL DZERO(AMIN,AMAX,ROOT,DUM,1D-8,100,ZALPHA)
C     PRINT*,1,I,AMIN,AMAX,FMIN,FMAX
*     CALL DZERO(AMIN,AMAX,ROOT,DUM,1D-8,100,ZALSFL)
      ALMSB5=AMZ*EXP(-ROOT/2)
*
*      ALMSB4=ALMSB5*(BMASS/ALMSB5)**(2D0/25)
*     &         *(LOG(BMASS**2/ALMSB5**2))**(963D0/14375)
cbardin
*      ALMSB4=.32D0
*      ALMSB3=ALMSB4*(CMASS/ALMSB4)**(2D0/27)
*     &         *(LOG(CMASS**2/ALMSB4**2))**(107D0/2025)
*      ALMSB6=ALMSB5*(ALMSB5/TMASS)**(2D0/21)
*     &         *(LOG(TMASS**2/ALMSB5**2))**(-107D0/1127)
*
* BETA/GAMA - COEFFICIENTS FOR ARBITRARY ANF
*
      B0=11D0-2D0/3D0*ANF
      B1=102D0-38D0/3D0*ANF
      B2=.5D0*(2857D0-5033D0/9D0*ANF+325D0/27D0*ANF**2)
      BETA0=B0/4
      BETA1=B1/16
      BETA2=B2/64
      GAMA0=1D0
      GAMA1=(202D0/3-20D0/9*ANF)/16
      GAMA2=(1249-(2216D0/27+160D0/3*D3)*ANF-140D0/81*ANF**2)/64
*
* BETA/GAMMA - COEFFICIENTS FOR ANF=3
*
      BETA03=(11D0-2D0/3D0*ANF3)/4
      BETA13=(102D0-38D0/3D0*ANF3)/16
      BETA23=(.5D0*(2857D0-5033D0/9D0*ANF3+325D0/27D0*ANF3**2))/64
      GAMA03=1D0
      GAMA13=(202D0/3-20D0/9*ANF3)/16
      GAMA23=(1249-(2216D0/27+160D0/3*D3)*ANF3-140D0/81*ANF3**2)/64
      COEF13=GAMA13/BETA03-BETA13*GAMA03/BETA03**2
      COEF23=GAMA23/BETA03-BETA13*GAMA13/BETA03**2
     &      -BETA23*GAMA03/BETA03**2+BETA13**2*GAMA03/BETA03**3
*
* BETA/GAMA - COEFFICIENTS FOR ANF=4
*
      BETA04=(11D0-2D0/3D0*ANF4)/4
      BETA14=(102D0-38D0/3D0*ANF4)/16
      BETA24=(.5D0*(2857D0-5033D0/9D0*ANF4+325D0/27D0*ANF4**2))/64
      GAMA04=1D0
      GAMA14=(202D0/3-20D0/9*ANF4)/16
      GAMA24=(1249-(2216D0/27+160D0/3*D3)*ANF4-140D0/81*ANF4**2)/64
      COEF14=GAMA14/BETA04-BETA14*GAMA04/BETA04**2
      COEF24=GAMA24/BETA04-BETA14*GAMA14/BETA04**2
     &      -BETA24*GAMA04/BETA04**2+BETA14**2*GAMA04/BETA04**3
*
* BETA/GAMA - COEFFICIENTS FOR ANF=5
*
      BETA05=(11D0-2D0/3D0*ANF5)/4
      BETA15=(102D0-38D0/3D0*ANF5)/16
      BETA25=(.5D0*(2857D0-5033D0/9D0*ANF5+325D0/27D0*ANF5**2))/64
      GAMA05=1D0
      GAMA15=(202D0/3-20D0/9*ANF5)/16
      GAMA25=(1249-(2216D0/27+160D0/3*D3)*ANF5-140D0/81*ANF5**2)/64
      COEF15=GAMA15/BETA05-BETA15*GAMA05/BETA05**2
      COEF25=GAMA25/BETA05-BETA15*GAMA15/BETA05**2
     &      -BETA25*GAMA05/BETA05**2+BETA15**2*GAMA05/BETA05**3
*
* FOUR LOOP BETA/GAMA
*
      BETA33=1D0/256D0*(149753D0/6D0+3564D0*D3
     &      -(1078361D0/162D0+6508D0/27D0*D3)*ANF3
     &  +(50065D0/162D0+6472D0/81D0*D3)*ANF3**2+1093D0/729D0*ANF3**3)
* test
*      print *,'0=',+1D0/256D0*(149753D0/6D0+3564D0*D3)
*      print *,'1=',+1D0/256D0*(-(1078361D0/162D0+6508D0/27D0*D3))
*      print *,'2=',+1D0/256D0*(+(50065D0/162D0+6472D0/81D0*D3))
*      print *,'3=',+1D0/256D0*(+1093D0/729D0)
*      stop
      BETA34=1D0/256D0*(149753D0/6D0+3564D0*D3
     &      -(1078361D0/162D0+6508D0/27D0*D3)*ANF4
     &  +(50065D0/162D0+6472D0/81D0*D3)*ANF4**2+1093D0/729D0*ANF4**3)
      BETA35=1D0/256D0*(149753D0/6D0+3564D0*D3
     &      -(1078361D0/162D0+6508D0/27D0*D3)*ANF5
     &  +(50065D0/162D0+6472D0/81D0*D3)*ANF5**2+1093D0/729D0*ANF5**3)
*
      GAMA33=1D0/256D0*(4603055D0/162D0+135680D0/27D0*D3-8800D0*D5
     &  +(-91723D0/27D0-34192D0/9D0*D3+880D0*D4+18400D0/9D0*D5)*ANF3
     &      +(5242D0/243D0+800D0/9D0*D3-160D0/3D0*D4)*ANF3**2
     &      +(-332D0/243D0+64D0/27D0*D3)*ANF3**3)
      GAMA34=1D0/256D0*(4603055D0/162D0+135680D0/27D0*D3-8800D0*D5
     &  +(-91723D0/27D0-34192D0/9D0*D3+880D0*D4+18400D0/9D0*D5)*ANF4
     &      +(5242D0/243D0+800D0/9D0*D3-160D0/3D0*D4)*ANF4**2
     &      +(-332D0/243D0+64D0/27D0*D3)*ANF4**3)
      GAMA35=1D0/256D0*(4603055D0/162D0+135680D0/27D0*D3-8800D0*D5
     &  +(-91723D0/27D0-34192D0/9D0*D3+880D0*D4+18400D0/9D0*D5)*ANF5
     &      +(5242D0/243D0+800D0/9D0*D3-160D0/3D0*D4)*ANF5**2
     &      +(-332D0/243D0+64D0/27D0*D3)*ANF5**3)
*
* Calculation of other \Lambda's
*
      ALMSB4=ALNFM1(BMASS,ALMSB5,5)
      ALMSBN=ANFM1N(BMASS,ALMSB5,5)
      ALMSB3=ALNFM1(CMASS,ALMSB4,4)
*
* Test of matching
*
*      print *,'ALMSB4=',ALMSB4
*      print *,'ALMSBN=',ALMSBN
*      ANF=5D0
*      ALMS=LOG(BMASS**2/ALMSB5**2)
*      ALP5=ALPHAS(ALMS)
*      print *,'ALP5=',ALP5
*      ANF=4D0
*      ALMS=LOG(BMASS**2/ALMSB4**2)
*      ALP4=ALPHAS(ALMS)
*      print *,'ALP4=',ALP4
*      ALMS=LOG(BMASS**2/ALMSBN**2)
*      ALP4=ALPHAS(ALMS)
*      print *,'ALP4=',ALP4
*      stop
*
* ALPHAS at 1 GeV, Mc-pole and Mb-pole
*
      ANF=3D0
      ALMS1=LOG((1D0)**2/ALMSB3**2)
      ALF1PI=ALPHAS(ALMS1)/PI
*
      ANF=4D0
      ALMC21=LOG(CMASS**2/ALMSB4**2)
      ALFCPI=ALPHAS(ALMC21)/PI
      ALMC22=LOG(4D0*CMASS**2/ALMSB4**2)
      CLFCPI=ALPHAS(ALMC22)/PI
      ANF=3D0
      ALMC21=LOG(CMASS**2/ALMSB3**2)
      BLFCPI=ALPHAS(ALMC21)/PI
      ALMC22=LOG(4D0*CMASS**2/ALMSB3**2)
      DLFCPI=ALPHAS(ALMC22)/PI
*
      ANF=5D0
      ALMB21=LOG(BMASS**2/ALMSB5**2)
      ALFBPI=ALPHAS(ALMB21)/PI
      ALMB22=LOG(4D0*BMASS**2/ALMSB5**2)
      CLFBPI=ALPHAS(ALMB22)/PI
      ANF=4D0
      ALMB21=LOG(BMASS**2/ALMSB4**2)
      BLFBPI=ALPHAS(ALMB21)/PI
      ALMB22=LOG(4D0*BMASS**2/ALMSB4**2)
      DLFBPI=ALPHAS(ALMB22)/PI
*
      ANF=5D0
      ALMS2=LOG( SQS**2          /ALMSB5**2)
      ALMT2=LOG((       TMASS)**2/ALMSB5**2)
      ALMX2=LOG((SCALET*TMASS)**2/ALMSB5**2)
      ALFSPI=ALPHAS(ALMS2)/PI
      ALFTPI=ALPHAS(ALMT2)/PI
      ALFXPI=ALPHAS(ALMX2)/PI
      ALFAT =ALFTPI*PI
      ALFAXI=ALFXPI*PI
*
      AMSAMC=SMRUN1*(ALFCPI/ALF1PI)**(GAMA03/BETA03)
     &      *(1+COEF13*(ALFCPI-ALF1PI)+.5D0*COEF13**2*(ALFCPI-ALF1PI)**2
     &         +.5D0*COEF23*(ALFCPI**2-ALF1PI**2))
*
*      BMSAMC=SMRUN1*(BLFCPI/ALF1PI)**(GAMA03/BETA03)
*     &      *(1+COEF13*(BLFCPI-ALF1PI)+.5D0*COEF13**2*(BLFCPI-ALF1PI)**2
*     &         +.5D0*COEF23*(BLFCPI**2-ALF1PI**2))
*
      AMSAMB=AMSAMC*(ALFBPI/ALFCPI)**(GAMA04/BETA04)
     &      *(1+COEF14*(ALFBPI-ALFCPI)+.5D0*COEF14**2*(ALFBPI-ALFCPI)**2
     &         +.5D0*COEF24*(ALFBPI**2-ALFCPI**2))
*
* RUNNING C-MASS, NF=4 EFFECTIVE THEORY
*
* Calculation of K_c
*
      AKC=2905D0/288+(7D0/3+2D0/3*LOG(2D0))*D2-1D0/6*D3
     &   -1D0/3*(71D0/48+D2)*ANF4
*
* Three ways of running m_c calculation:
*
* 1) Running m_c(m_c)
*
*      CALL DZERO(.70D0,1.30D0,ROOT,DUM,1D-8,100,ZRMCIN)
*      AMCSIN=ROOT
*      print *,'AMCSIN=',AMCSIN
*
* RGE-running from m_c to M_c
*
*      ANF=4D0
*      ALMSIN=LOG(AMCSIN**2/ALMSB4**2)
*      ALFSIN=ALPHAS(ALMSIN)/PI
*      AMCAMC=AMCSIN*(ALFCPI/ALFSIN)**(GAMA03/BETA03)
*     &      *(1+COEF13*(ALFCPI-ALFSIN)+.5D0*COEF13**2*(ALFCPI-ALFSIN)**2
*     &         +.5D0*COEF23*(ALFCPI**2-ALFSIN**2))
*      print *,'AMCAMC_1(3)=',AMCAMC
*      AMCAMC=AMCSIN*(ALFCPI/ALFSIN)**(GAMA04/BETA04)
*     &      *(1+COEF14*(ALFCPI-ALFSIN)+.5D0*COEF14**2*(ALFCPI-ALFSIN)**2
*     &         +.5D0*COEF24*(ALFCPI**2-ALFSIN**2))
*      print *,'AMCAMC_1(4)=',AMCAMC
*
* 2) Running m_c(M_c)
*
* Second MG's insertion
*
      AMIN=.70D0
      AMAX=1.30D0
      DO I=1,10
      FMIN=ZRMCMC(AMIN)
      FMAX=ZRMCMC(AMAX)
C     PRINT*,2,I,AMIN,AMAX,FMIN,FMAX
      IF ((FMIN.LE.0D0.AND.FMAX.LE.0D0).OR.
     &    (FMIN.GE.0D0.AND.FMAX.GE.0D0)) THEN
        AMIN=AMIN/2D0
        AMAX=AMAX*2D0
      ELSE
        GOTO 98
      ENDIF
      ENDDO
 98   CONTINUE
      CALL DZERO(AMIN,AMAX,ROOT,DUM,1D-8,100,ZRMCMC)
C     PRINT*,2,I,AMIN,AMAX,FMIN,FMAX
      AMCAMC=ROOT
*      print *,'AMCAMC_2=',AMCAM2
*
*      AMCAMC=CMASS/(1+4D0/3*ALFCPI+13.3D0*ALFCPI**2)
*      print *,'AMCAMC_O=',AMCAMC
*
      AMCAMB=AMCAMC*(ALFBPI/ALFCPI)**(GAMA04/BETA04)
     &      *(1+COEF14*(ALFBPI-ALFCPI)+.5D0*COEF14**2*(ALFBPI-ALFCPI)**2
     &         +.5D0*COEF24*(ALFBPI**2-ALFCPI**2))
*      print *,'AMCAMB_2(4)=',AMCAMB
*      AMCAMB=AMCAMC*(ALFBPI/ALFCPI)**(GAMA05/BETA05)
*     &      *(1+COEF15*(ALFBPI-ALFCPI)+.5D0*COEF15**2*(ALFBPI-ALFCPI)**2
*     &         +.5D0*COEF25*(ALFBPI**2-ALFCPI**2))
*      print *,'AMCAMB_2(5)=',AMCAMB
*
      AMCRUN=AMCAMB*(ALFSPI/ALFBPI)**(GAMA05/BETA05)
     &      *(1+COEF15*(ALFSPI-ALFBPI)+.5D0*COEF15**2*(ALFSPI-ALFBPI)**2
     &         +.5D0*COEF25*(ALFSPI**2-ALFBPI**2))
      AMQRUN(3)=AMCRUN
*
* RUNNING B-MASS, NF=5 EFFECTIVE THEORY
*
* Calculation of K_b
*
      AKB=3817D0/288+2D0/3*(2D0+LOG(2D0))*D2-1D0/6*D3-8D0/3
     &   +D2-1D0/2-1D0/3*(71D0/48+D2)*ANF5
*
* Third MG's insertion
*
      AMIN=3.5D0
      AMAX=4.5D0
      DO I=1,10
      FMIN=ZRMBMB(AMIN)
      FMAX=ZRMBMB(AMAX)
C     PRINT*,3,I,AMIN,AMAX,FMIN,FMAX
      IF ((FMIN.LE.0D0.AND.FMAX.LE.0D0).OR.
     &    (FMIN.GE.0D0.AND.FMAX.GE.0D0)) THEN
        AMIN=AMIN/2D0
        AMAX=AMAX*2D0
      ELSE
        GOTO 99
      ENDIF
      ENDDO
 99   CONTINUE
      CALL DZERO(AMIN,AMAX,ROOT,DUM,1D-8,100,ZRMBMB)
C     PRINT*,3,I,AMIN,AMAX,FMIN,FMAX
      AMBAMB=ROOT
*
*      AMBAMB=BMASS*(1-4D0/3*ALFBPI-(12.4D0-16D0/9)*ALFBPI**2)
*
      AMBRUN=AMBAMB*(ALFSPI/ALFBPI)**(GAMA05/BETA05)
     &      *(1+COEF15*(ALFSPI-ALFBPI)+.5D0*COEF15**2*(ALFSPI-ALFBPI)**2
     &         +.5D0*COEF25*(ALFSPI**2-ALFBPI**2))
      AMQRUN(6)=AMBRUN
*
      DO IMQ=1,6
       CMQRUN(IMQ)=AMQRUN(IMQ)
      ENDDO
*
      RATU=-(1D0+4D0/3*SW2)/(1-8D0/3*SW2)
      RATD=+(1D0+4D0/3*SW2)/(1-4D0/3*SW2)
      RATV=+(1D0+4D0/3*SW2)**2
*
* Re-implementation of QCD-corrections
* References: CKK, CERN 95-03
*              CK, MPI/PhT/96-84, hep-ph/9609202
*                  
cbardin
cb        print *,'SQRTs =',SQS
cb        print *,'SMASS =',SMASS
cb        print *,'CMASS =',CMASS
cb        print *,'BMASS =',BMASS
cb        print *,'TMASS =',AMT
cb        print *,'ALMSB3=',ALMSB3
cb        print *,'ALMSB4=',ALMSB4
cb        print *,'ALMSB5=',ALMSB5
cb        print *,'ALFAS =',ALFAS
cb        print *,'ALFAS1=',ALF1PI*PI
cb        print *,'ALFAMC=',ALFCPI*PI
cb        print *,'BLFAMC=',BLFCPI*PI
cb        print *,'CLFAMC=',CLFCPI*PI
cb        print *,'DLFAMC=',DLFCPI*PI
cb        print *,'ALFAMB=',ALFBPI*PI
cb        print *,'BLFAMB=',BLFBPI*PI
cb        print *,'CLFAMB=',CLFBPI*PI
cb        print *,'DLFAMB=',DLFBPI*PI
cb        print *,'ALFASS=',ALFSPI*PI
cb        print *,'AMSAMC=',AMSAMC
cb        print *,'AMSAMB=',AMSAMB
cb        print *,'AMCAMC=',AMCAMC
cb        print *,'AMCAMB=',AMCAMB
cb        print *,'AMCRUN=',AMCRUN
cb        print *,'AMBAMB=',AMBAMB
cb        print *,'AMBRUN=',AMBRUN
*
* Massless corrections
***
* mod. 14 Jan 2013
* decide what you want to evaluate:
* ORIGINAL MASSLESS QCD CORRECTIONS:
**********************************     GOTO 2005
* BAIKOV ET AL. 0801.1821v2, 2008:     GOTO 2008
* BAIKOV ET AL. 1201.5804v3, 2012:     GOTO 2012
*
      IBAIKOV=2012
      print *, '   IBAIKOV =  ',IBAIKOV
      IF(IBAIKOV.EQ.2005) GOTO 2005
      IF(IBAIKOV.EQ.2008) GOTO 2008
      IF(IBAIKOV.EQ.2012) GOTO 2012
*
** the original massless corrections of dizet6_42.f and older:
*
 2005 CONTINUE
      ANF=5D0
      COEF01=1D0
      COEF02=365D0/24-11D0*D3+(-11D0/12+2D0/3*D3)*ANF
      COEF03=87029D0/288 -121D0/8*D2-1103D0/4*D3+275D0/6*D5
     &      +(-7847D0/216+ 11D0/6*D2+ 262D0/9*D3- 25D0/9*D5)*ANF
     &      +(151D0/162  - 1D0/18*D2- 19D0/27*D3           )*ANF**2
***
* mod. 14 Jan 2013
      COEF04=0D0
      GOTO 2020
* The above COEF01,COEF02,COEF03 contain only massless non-singlet 
* corrections, see CKK in 95-03, appendix A.3, eq. (227) for R1
* In BCK 2008 eq.(5) and BCKR 2012v3 eq.(6), the additional COEF04 
* is added:
* 
 2008 CONTINUE
      ANF=5D0
      COEF01=1D0
* the -0.1152*n_f in COEF02 from 2008 should be more precisely 
* -0.115295 = -0.1153
* the latter value was used in COEF02 of 2005. So there is a small numerical 
* difference between the two values of COEF02 with no essential origin.
      COEF02=1.9857D0+(-0.1152D0)*ANF
*           = 1.4097D0
      COEF03=-6.63694D0+(-1.20013D0)*ANF+(0.00518D0)*ANF**2
*           = -12.7671D0
      COEF04=-156.61D0+18.77D0*ANF+(-0.7974D0)*ANF**2+0.0215D0*ANF**3
*           = -79.9806D0
      GOTO 2020
*
* In Baikov Chetyrkin Kuehn Rittinger 2012v3, there are additionally 
* vector and axial singlet corrections, so not only COEF04
* is added, but also CAI4, and the RVSING is modified. 
*
 2012 CONTINUE
      ANF=5D0
      COEF01=1D0
*     COEF02=1.9857D0+(-0.1152D0)*ANF
      COEF02=1.4092D0
*     COEF03=-6.63694D0+(-1.20013D0)*ANF+(0.00518D0)*ANF**2
      COEF03=-12.7671D0
*     COEF04=-156.61D0+18.77D0*ANF+(-0.7974D0)*ANF**2+0.0215D0*ANF**3
      COEF04=-79.9806D0 
*
 2020 CONTINUE
*
*****************
*
* Quadratic Corrections
*
* We consider the case when only AMCRUN and AMBRUN are retained, 
* i.e these corrections should be applied only for sqrt(S).GE.13-15 GeV
*
* Light quarks
*
      COEFL1=0D0
      COEFL2=0D0
      COEFL3=-80D0+60D0*D3+(32D0/9-8D0/3*D3)*ANF
*
* Heavy quarks
*
      COEFV1=12D0
      COEFV2=253D0/2-13D0/3*ANF
      COEFV3=  2522D0   - 855D0/2*D2+ 310D0/3*D3- 5225D0/6*D5
     &      +(-4942D0/27+    34D0*D2-394D0/27*D3+1045D0/27*D5)*ANF
     &      +(125D0/54  -   2D0/3*D2                      )*ANF**2
*
      COEFA0=-6D0
      COEFA1=-22D0
      COEFA2=-8221D0/24+57D0*D2+117D0*D3
     &      +(151D0/12 - 2D0*D2-  4D0*D3)*ANF   
      COEFA3=-4544045D0/864+   1340*D2+118915D0/36*D3  -1270D0*D5
     &      +(71621D0/162  -209D0/2*D2   -216D0*D3+5D0*D4+55D0*D5)*ANF
     &      +(-13171D0/1944+ 16D0/9*D2  +26D0/9*D3            )*ANF**2
*
*     PRINT *,'COEF01,2,3=',COEF01,COEF02,COEF03
*     PRINT *,'COEFL1,2,3=',COEFL1,COEFL2,COEFL3
*     PRINT *,'(COEFV3+COEFL3)=',(COEFV3+COEFL3)/12 
*
      RLMC=AMCRUN**2/SQS**2
      RLMB=AMBRUN**2/SQS**2
      R2LU=(RLMC+RLMB)*COEFL3*ALFSPI**3
      RV20=COEFV1*ALFSPI+COEFV2*ALFSPI**2+COEFV3*ALFSPI**3
      RA20=COEFA0+COEFA1*ALFSPI+COEFA2*ALFSPI**2+COEFA3*ALFSPI**3
*
      QCDC2V(1)=R2LU
      QCDC2V(2)=R2LU
      QCDC2V(3)=R2LU+RLMC*RV20
      QCDC2V(4)=R2LU
      QCDC2V(5)=0D0
      QCDC2V(6)=R2LU+RLMB*RV20
*
      QCDC2A(1)=R2LU
      QCDC2A(2)=R2LU
      QCDC2A(3)=R2LU+RLMC*RA20
      QCDC2A(4)=R2LU
      QCDC2A(5)=0D0
      QCDC2A(6)=R2LU+RLMB*RA20
*
* Quartic Corrections. They contain ln(M^2/S) and can't be coded 
*                      as above as consts
*
* For bb-channel, the vector m^6_b correction is added
*                                                                   
* Light quarks                                             
*
      ALMC=LOG(RLMC)
      R4LC=RLMC**2*(13D0/3-ALMC-4D0*D3)*ALFSPI**2
      ALMB=LOG(RLMB)
      R4LB=RLMB**2*(13D0/3-ALMB-4D0*D3)*ALFSPI**2
*
* Heavy quarks (actually b-quark)                 
*                                
      RV40=-6D0-22D0*ALFSPI+(+12D0-3173D0/12+27D0*PI2+112D0*D3
     &     +(143D0/18-2D0/3*PI2- 8D0/3*D3)*ANF)*ALFSPI**2                      
*                                
*     PRINT *,'RV40=',12D0-3173D0/12+27D0*PI2+112D0*D3
*    &     +(143D0/18-2D0/3*PI2- 8D0/3*D3)*ANF                      
*
      RA40=+6D0+10D0*ALFSPI+(-12D0+3533D0/12-27D0*PI2-220D0*D3
     &     +(-41D0/6 +2D0/3*PI2+16D0/3*D3)*ANF)*ALFSPI**2
      RV4L=(-11D0/2+1D0/3*ANF)*ALFSPI**2
      RA4L=(+77D0/2-7D0/3*ANF)*ALFSPI**2
*
      QCDC4V(1)=R4LC+R4LB
      QCDC4V(2)=R4LC+R4LB
      QCDC4V(3)=R4LC+R4LB+RLMC**2*(RV40+RV4L*ALMC)+12*RLMB**2*ALFSPI**2
      QCDC4V(4)=R4LC+R4LB
      QCDC4V(5)=0D0
      QCDC4V(6)=R4LC+R4LB+RLMB**2*(RV40+RV4L*ALMB)+12*RLMC**2*ALFSPI**2
     &                   -RLMB**3*(8D0+16D0/27*(155D0+6D0*ALMB)*ALFSPI)
*
      QCDC4A(1)=R4LC+R4LB
      QCDC4A(2)=R4LC+R4LB
      QCDC4A(3)=R4LC+R4LB+RLMC**2*(RA40+RA4L*ALMC)-12*RLMB**2*ALFSPI**2
      QCDC4A(4)=R4LC+R4LB
      QCDC4A(5)=0D0
      QCDC4A(6)=R4LC+R4LB+RLMB**2*(RA40+RA4L*ALMB)-12*RLMC**2*ALFSPI**2
*
      RLMT=SQS**2/AMT**2
      ALMT=LOG(RLMT)
      CAI2=-37D0/12+ALMT+7D0/81*RLMT+.0132D0*RLMT**2
* comment 14 Jan 2013, see Baikov et al. 2012 eq. (3):
*     R_{s:T,B}^a=(-3.08333 + lt)a_s^2 +...=-4.35248a_s^2 +..., 
* this is CAI2 with RLTM=0
* we leave CAI2 as it was in dizet6_42,43
* comment 14 Jan 2013, see Baikov et al. 2012 eqs. (3,4):
*     R_{s:T,B}^a=..+(..)a_s^3+..=..-17.6245a_s^3.., this is CAI3 
* for M_Z=91.1875 GeV, m_t=172 GeV
      CAI3=-5651D0/216+8D0/3+23D0/6*D2+D3+67D0/18*ALMT+23D0/12*ALMT**2
* mod. 14 Jan 2013
      CAI4=0D0
* see Baikov et al. 2012v3 eqs. (3,4):
* R_{s:T,B}^a= ... +()a_s^$, this is CAI4
      IF(IBAIKOV.EQ.2012)  CAI4=+87.5520D0
*
      DO 10 IQ=1,6
      SINS=(-1)**IQ
*           
* R_V
*
      QCDCON(2*IQ-1)=1+ALFSPI
     &              +1D0/4*ALQED/PI*CHARQU(IQ)**2*(3*ZQED-ALFSPI*ZMIX)
     &              +(1.40923D0+(44D0/675-2D0/135*ALMT)*RLMT)*ALFSPI**2
     &              +(-12.76706D0)*ALFSPI**3
     &+12*AMQRUN(IQ)**2/SQS**2*ALFSPI*(1+8.7D0*ALFSPI+45.15D0*ALFSPI**2)
*
* mod. 14 Jan 2013: COEF04 ADDED
      QCDCOR(2*IQ-1)=1+ALFSPI
     &              +1D0/4*ALQED/PI*CHARQU(IQ)**2*(3*ZQED-ALFSPI*ZMIX)
     &              +(COEF02+(44D0/675-2D0/135*ALMT)*RLMT)*ALFSPI**2
     &              + COEF03*ALFSPI**3
     &              + COEF04*ALFSPI**4
     &              +QCDC2V(IQ)+QCDC4V(IQ)
*
* R_A
*
      QCDCON(2*IQ  )=1+ALFSPI
     &   +1D0/4*ALQED/PI*CHARQU(IQ)**2*(3*ZQED-ALFSPI*ZMIX)
     &   +(1.40923D0+(44D0/675-2D0/135*ALMT)*RLMT+SINS*CAI2)*ALFSPI**2
     &   +(-12.76706D0+SINS*CAI3)*ALFSPI**3
     &   -6*AMQRUN(IQ)**2/SQS**2*
     &                      (1+11D0/3*ALFSPI+(11.286D0+ALMT)*ALFSPI**2)
     &         -10*AMQRUN(IQ)**2/AMT**2*(8D0/81-1D0/54*ALMT)*ALFSPI**2
*
* Axial singlet contributions
*
      IF    (IQ.EQ.3) THEN
       RASING=-6D0*AMCRUN**2/SQS**2*(-3D0+ALMT)*ALFSPI**2
     &       -10D0*AMCRUN**2/AMT**2*(8D0/81-1D0/54*ALMT)*ALFSPI**2
      ELSEIF(IQ.EQ.6) THEN
       RASING=-6D0*AMBRUN**2/SQS**2*(-3D0+ALMT)*ALFSPI**2
     &       -10D0*AMBRUN**2/AMT**2*(8D0/81-1D0/54*ALMT)*ALFSPI**2
      ELSE
       RASING=0D0
      ENDIF
*
* mod. 14 Jan 2013: COEF04 and CAI4 ADDED
      QCDCOR(2*IQ  )=1+ALFSPI
     &        +1D0/4*ALQED/PI*CHARQU(IQ)**2*(3*ZQED-ALFSPI*ZMIX)
     &        +(COEF02+(44D0/675-2D0/135*ALMT)*RLMT+SINS*CAI2)*ALFSPI**2
     &        +(COEF03+SINS*CAI3)*ALFSPI**3                            
     &        +(COEF04+SINS*CAI4)*ALFSPI**4                            
     &        +QCDC2A(IQ)+QCDC4A(IQ)                                   
     &        +RASING
  10  CONTINUE   
*
* Vector singlet contribution
*
*     RVSING=1D0/192*(176D0/3-128D0*D3)
*
* comment 14 Jan 2013: this is = -0.495816
*     -.41317D0 = 5/6*1D0/192*(176D0/3-128D0*D3)
* so may assume that the analytical formula is wrong by lacking 5/6
* In CKK in YR 95-03 eq. (180) with dabc^2 (?=40/3) one gets these 
* numbers
      RVSING=-.41317D0
* mod. 14 Jan 2013:
      IF(IBAIKOV.EQ.2012) RVSING=-.41317D0 - 4.9841*ALFSPI 
*
      QCDCOR(13)=RVSING*ALFSPI**3
* QCDCOR(14) - Arbuzov, Bardin, Leike correction
      QCDCOR(14)=16D0/3*ALFSPI
*
* Prints for Shirkov
*
*     PRINT *,'NF  =',ANF
*     PRINT *,'FV_2=',COEF02
*     PRINT *,'FV_3=',COEF03
*     PRINT *,'FA_2=',COEF02+(44D0/675-2D0/135*ALMT)*RLMT
*     PRINT *,'FA_3=',COEF03
*     PRINT *,'C2,3=',CAI2,CAI3
*     PRINT *,'FA_2=',COEF02+(44D0/675-2D0/135*ALMT)*RLMT+SINS*CAI2
*     PRINT *,'FA_3=',COEF03+SINS*CAI3
*     D2EF03=(-121D0/8+11D0/6*ANF- 1D0/18*ANF**2)*D2
*     PRINT *,'D2EF03 =',D2EF03,-1D0/3*(PI*(11-2D0/3*ANF)/4)**2
*     D2EF04=COEF02*(-(PI*(11-2D0/3*ANF)/4)**2)
*     D2EF05=+1D0/5*(PI*(11-2D0/3*ANF)/4)**4
*     PRINT *,'D2EF04 =',D2EF04,D2EF04*ALFSPI**4
*     PRINT *,'D2EF05 =',D2EF05,D2EF05*ALFSPI**5
*     PRINT *,'CV1,2,3=',COEFV1*RLMB,COEFV2*RLMB,COEFV3*RLMB
*     PRINT *,'CA1,2,3=',COEFA1*RLMB,COEFA2*RLMB,COEFA3*RLMB
*     STOP
*
      END
 
      FUNCTION ZRMCIN(RMASS)
      IMPLICIT REAL*8(A-H,O-Z)
*
      COMMON /CDZZER/ PI,ALPHAQ,ANF,AKC,AKB
      COMMON /CDZPHM/ CMASS,BMASS,TMASS,SMRUN1
      COMMON /CDZBGC/ BETA03,BETA13,BETA23,GAMA03,GAMA13,GAMA23,
     &                BETA04,BETA14,BETA24,GAMA04,GAMA14,GAMA24,
     &                BETA05,BETA15,BETA25,GAMA05,GAMA15,GAMA25,
     &                COEF13,COEF23,COEF14,COEF24,COEF15,COEF25,
     &                ALMSB3,ALMSB4,ALMSB5,ALMSB6,ALF1PI,ALFCPI,ALFBPI,
     &                AMSAMC,AMSAMB,AMCAMC,AMCAMB,AMBAMB
*
      IF(RMASS.LE.1D0.OR.RMASS.LE.ALMSB3) THEN
       AMSRMC=SMRUN1
      ELSE
       ANF=3D0
       ALMC2=LOG(RMASS**2/ALMSB3**2)
       ALFRPI=ALPHAS(ALMC2)/PI
       AMSRMC=SMRUN1*(ALFRPI/ALF1PI)**(GAMA03/BETA03)
     &      *(1+COEF13*(ALFRPI-ALF1PI)+.5D0*COEF13**2*(ALFRPI-ALF1PI)**2
     &         +.5D0*COEF23*(ALFRPI**2-ALF1PI**2))
      ENDIF
*
      ANF=4D0
      ALMC2 =LOG(RMASS**2/ALMSB4**2)
      ALFRPI=ALPHAS(ALMC2)/PI
      RAT=AMSRMC/RMASS
      DEL=(PI**2/8-.597D0*RAT+.23D0*RAT**2)*RAT
*
      ZRMCIN=CMASS-RMASS*(1D0+4D0/3*ALFRPI+(AKC+4D0/3*DEL)*ALFRPI**2)
*
      END
 
      FUNCTION ZRMCMC(RMASS)
      IMPLICIT REAL*8(A-H,O-Z)
*
      COMMON /CDZZER/ PI,ALPHAQ,ANF,AKC,AKB
      COMMON /CDZPHM/ CMASS,BMASS,TMASS,SMRUN1
      COMMON /CDZBGC/ BETA03,BETA13,BETA23,GAMA03,GAMA13,GAMA23,
     &                BETA04,BETA14,BETA24,GAMA04,GAMA14,GAMA24,
     &                BETA05,BETA15,BETA25,GAMA05,GAMA15,GAMA25,
     &                COEF13,COEF23,COEF14,COEF24,COEF15,COEF25,
     &                ALMSB3,ALMSB4,ALMSB5,ALMSB6,ALF1PI,ALFCPI,ALFBPI,
     &                AMSAMC,AMSAMB,AMCAMC,AMCAMB,AMBAMB
*
      ANFE=4D0
      RAT=AMSAMC/RMASS
      DEL=(PI**2/8-.597D0*RAT+.23D0*RAT**2)*RAT
      ALRMC=LOG(CMASS**2/RMASS**2)
*
      ZRMCMC=CMASS-RMASS*(1D0+(4D0/3+ALRMC)*ALFCPI
     &+(AKC+4D0/3*DEL+(173D0/24-13D0/36*ANFE)*ALRMC
     &+(15D0/8-1D0/12*ANFE)*ALRMC**2)*ALFCPI**2)
cbardin
cb      print *,'First =',(4D0/3+ALRMC)*ALFCPI
cb      print *,'Second=',
cb     &+(AKC+4D0/3*DEL+(173D0/24-13D0/36*ANFE)*ALRMC
cb     &+(15D0/8-1D0/12*ANFE)*ALRMC**2)*ALFCPI**2
*
      END
 
      FUNCTION ZRMBMB(RMASS)
      IMPLICIT REAL*8(A-H,O-Z)
*
      COMMON /CDZZER/ PI,ALPHAQ,ANF,AKC,AKB
      COMMON /CDZPHM/ CMASS,BMASS,TMASS,SMRUN1
      COMMON /CDZBGC/ BETA03,BETA13,BETA23,GAMA03,GAMA13,GAMA23,
     &                BETA04,BETA14,BETA24,GAMA04,GAMA14,GAMA24,
     &                BETA05,BETA15,BETA25,GAMA05,GAMA15,GAMA25,
     &                COEF13,COEF23,COEF14,COEF24,COEF15,COEF25,
     &                ALMSB3,ALMSB4,ALMSB5,ALMSB6,ALF1PI,ALFCPI,ALFBPI,
     &                AMSAMC,AMSAMB,AMCAMC,AMCAMB,AMBAMB
*
      ANFE=5D0
      RATS=AMSAMB/RMASS
      RATC=AMCAMB/RMASS
      DELS=(PI**2/8-.597D0*RATS+.23D0*RATS**2)*RATS
      DELC=(PI**2/8-.597D0*RATC+.23D0*RATC**2)*RATC
      DEL =DELS+DELC
      ALRMB=LOG(BMASS**2/RMASS**2)
*
      ZRMBMB=BMASS-RMASS*(1D0+(4D0/3+ALRMB)*ALFBPI
     &+(AKB+4D0/3*DEL+(173D0/24-13D0/36*ANFE)*ALRMB
     &+(15D0/8-1D0/12*ANFE)*ALRMB**2)*ALFBPI**2)
cbardin
cb      print *,'First =',(4D0/3+ALRMB)*ALFBPI
cb      print *,'Second=',
cb     &+(AKB+4D0/3*DEL+(173D0/24-13D0/36*ANFE)*ALRMB
cb     &+(15D0/8-1D0/12*ANFE)*ALRMB**2)*ALFBPI**2
*
      END
 
      FUNCTION ZALPHA(A)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /CDZZER/ PI,ALPHAQ,ANF,AKC,AKB
      B=LOG(A)
      B0=11D0-2D0/3D0*ANF
      B1=102D0-38D0/3D0*ANF
      B2=.5D0*(2857D0-5033D0/9D0*ANF+325D0/27D0*ANF**2)
      C=B1/(B0**2*A)
      ZALPHA=ALPHAQ-
     &     4D0*PI/(B0*A)*(1D0-C*B+C**2*((B-.5D0)**2+B2*B0/B1**2-1.25D0))
      END
 
      FUNCTION ZALSFL(A)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /CDZZER/UPI,ALPHAQ,ANF,AKC,AKB
      COMMON /CDZCON/PI,PI2,F1,D3,ALFAI,AL4PI,AL2PI,AL1PI
      B=LOG(A)
      B0=1D0/4D0 *(11D0-2D0/3D0*ANF)
      B1=1D0/16D0*(102D0-38D0/3D0*ANF)/B0
      B2=1D0/64D0*(2857D0/2D0-5033D0/18D0*ANF+325D0/54D0*ANF**2)/B0
      B3=1D0/256D0*(149753D0/6D0+3564D0*D3
     &  -(1078361D0/162D0+6508D0/27D0*D3)*ANF
     &  +(50065D0/162D0+6472D0/81D0*D3)*ANF**2+1093D0/729D0*ANF**3)/B0
*
      ZALSFL=ALPHAQ-PI*(1D0/(B0*A)
     &                 -1D0/(B0*A)**2* B1*B
     &                 +1D0/(B0*A)**3*(B1**2*(B**2-B-1D0)+B2))
     &      -1D0/(B0*A)**4*(B1**3*(B**3-2.5D0*B**2-2D0*B+.5D0)
     &                      +3D0*B1*B2*B-.5D0*B3)    
*
      END
 
      FUNCTION ALPHAS(A)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /CDZZER/ PI,ALPHAQ,ANF,AKC,AKB
      B=LOG(A)
      B0=11D0-2D0/3D0*ANF
      B1=102D0-38D0/3D0*ANF
      B2=.5D0*(2857D0-5033D0/9D0*ANF+325D0/27D0*ANF**2)
      C=B1/(B0**2*A)
      ALPHAS=
     &     4D0*PI/(B0*A)*(1D0-C*B+C**2*((B-.5D0)**2+B2*B0/B1**2-1.25D0))
      END

      FUNCTION ALNFM1(AMQ,ALNFN,NF)
      IMPLICIT REAL*8(A-H,O-Z)
*
      COMMON /CDZBGC/ BETA03,BETA13,BETA23,GAMA03,GAMA13,GAMA23,
     &                BETA04,BETA14,BETA24,GAMA04,GAMA14,GAMA24,
     &                BETA05,BETA15,BETA25,GAMA05,GAMA15,GAMA25,
     &                COEF13,COEF23,COEF14,COEF24,COEF15,COEF25,
     &                ALMSB3,ALMSB4,ALMSB5,ALMSB6,ALF1PI,ALFCPI,ALFBPI,
     &                AMSAMC,AMSAMB,AMCAMC,AMCAMB,AMBAMB
*
      IF(NF.EQ.5) THEN
        BETA0N=BETA05
        BETA1N=BETA15 
        BETA2N=BETA25
        BETA0M=BETA04
        BETA1M=BETA14
        BETA2M=BETA24 
      ELSEIF(NF.EQ.4) THEN
        BETA0N=BETA04
        BETA1N=BETA14
        BETA2N=BETA24
        BETA0M=BETA03
        BETA1M=BETA13
        BETA2M=BETA23
      ELSE
        PRINT *,'ERROR IN ALNFM1 CALL'
        STOP
      ENDIF
*
      ALG=2D0*LOG(AMQ/ALNFN)
      ALGALG=LOG(ALG)
*
      CO1=BETA0N-BETA0M
      CO2=BETA1N/BETA0N-BETA1M/BETA0M
      CO3=-BETA1M/BETA0M
      CO4=BETA1N/BETA0N**2*(BETA1N/BETA0N-BETA1M/BETA0M)
      CO5=1D0/BETA0N*((BETA1N/BETA0N)**2-(BETA1M/BETA0M)**2
     &                -BETA2N/BETA0N+BETA2M/BETA0M+7D0/24D0)
*     &                -BETA2N/BETA0N+BETA2M/BETA0M-7D0/72D0)
*
      RHS=(CO1*ALG+CO2*ALGALG+CO3*LOG(BETA0N/BETA0M)
     &    +CO4*ALGALG/ALG+CO5/ALG)/BETA0M
      ALNFM1=ALNFN/(EXP(RHS/2D0))
*
      END

      FUNCTION ANFM1N(AMQ,ALNFN,NF)
      IMPLICIT REAL*8(A-H,O-Z)
*
* This is 1997 version of ALNFM1 with possible inclusion of four loops
*
      COMMON/CDZCON/PI,PI2,D2,D3,ALFAI,AL4PI,AL2PI,AL1PI
      COMMON /CDZBGC/ BETA03,BETA13,BETA23,GAMA03,GAMA13,GAMA23,
     &                BETA04,BETA14,BETA24,GAMA04,GAMA14,GAMA24,
     &                BETA05,BETA15,BETA25,GAMA05,GAMA15,GAMA25,
     &                COEF13,COEF23,COEF14,COEF24,COEF15,COEF25,
     &                ALMSB3,ALMSB4,ALMSB5,ALMSB6,ALF1PI,ALFCPI,ALFBPI,
     &                AMSAMC,AMSAMB,AMCAMC,AMCAMB,AMBAMB
      COMMON /CDZBG3/ BETA33,BETA34,BETA35,GAMA33,GAMA34,GAMA35
*
* Returns \Lambda(n_f-1)=ANFM1N, given \Lambda(n_f)=ALNFN
*
      IF(NF.EQ.5) THEN
        BETA0N=BETA05
        BETA1N=BETA15 
        BETA2N=BETA25
        BETA3N=BETA35
        BETA0M=BETA04
        BETA1M=BETA14
        BETA2M=BETA24 
        BETA3M=BETA34 
      ELSEIF(NF.EQ.4) THEN
        BETA0N=BETA04
        BETA1N=BETA14
        BETA2N=BETA24
        BETA3N=BETA34
        BETA0M=BETA03
        BETA1M=BETA13
        BETA2M=BETA23
        BETA3M=BETA33
      ELSE
        PRINT *,'ERROR IN ALNFM1 CALL'
        STOP
      ENDIF
*
* calculation of b_i and b'_i
*
      b1 =BETA1N/BETA0N    
      b2 =BETA2N/BETA0N    
      b3 =BETA3N/BETA0N    
      b1p=BETA1M/BETA0M    
      b2p=BETA2M/BETA0M    
      b3p=BETA3M/BETA0M    
*
      ALG=2D0*LOG(AMQ/ALNFN)
      ALGALG=LOG(ALG)
*
      Bc2=-7D0/24D0
      Bc3=-80507D0/27648D0*D3-2D0/3D0*(1D0/3D0*LOG(2D0)+1D0)*D2
     &   -58933D0/124416D0+1D0/9D0*(D2+2479D0/3456D0)*4D0
*
      CO1=BETA0N-BETA0M
      CO2=b1-b1p
      CO3=-b1p
      CO4=b1*(b1-b1p)
      CO5=b1**2-b1p**2-b2+b2p-Bc2
      CO6=b1**2/2D0*(b1-b1p)
      CO7=b1*(b1p*(b1-b1p)-b2+b2p-Bc2)
      CO8=.5D0*(-b1**3-b1p**3+b3-b3p)
     &   +b1p*(b1**2+b2p-b2-Bc2)+Bc3
*
      RHS=(CO1*ALG+CO2*ALGALG+CO3*LOG(BETA0N/BETA0M)
     &    +1D0/(BETA0N*ALG)*(CO4*ALGALG+CO5)
*     &    -1D0/(BETA0N*ALG)**2*(CO6*ALGALG**2+CO7*ALGALG+CO8)
     &    )/BETA0M
      ANFM1N=ALNFN/(EXP(RHS/2D0))
*
      END
 
      FUNCTION AFMT3(ALST,AMT2,AMZ2,SW2)
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      COMMON/CDZSCT/ISCRE,ISCAL,IAFMT,IFACR,IFACT,IHIGS,IEWLC,ICZAK
     &             ,IHIG2,IALE2,IGFER      
*
* NUMERICAL CONSTANTS
*
      PI=ATAN(1D0)*4D0
      PI2=PI**2
      D2=PI2/6D0
      D3=1.2020569031596D0
      D4=PI2**2/90
      AL2=LOG(2D0)
      TS2=+0.2604341376322D0
      TD3=-3.0270094939877D0
      TB4=-1.7628000870738D0
      AMU2=AMT2
      ALMU=LOG(AMU2/AMT2)
      NF=6
      CA1=-2D0/3*(1+2*D2)
      CA2L=157D0/648-3313D0/162*D2-308D0/27*D3+143D0/18*D4-4D0/3*D2*AL2
     & +441D0/8*TS2-1D0/9*TB4-1D0/18*TD3-(1D0/18-13D0/9*D2+4D0/9*D3)*NF
     & -(11D0/6-1D0/9*NF)*(1+2*D2)*ALMU
      ALZT=LOG(AMZ2/AMT2)
      CA2C=AMZ2/AMT2*(-17.2240D0+0.08829D0*ALZT+0.4722D0*ALZT**2
     &               +(22.6367D0+1.25270D0*ALZT-0.8519D0*ALZT**2)*SW2)
      CA2I=(AMZ2/AMT2)**2*(-7.7781D0-0.072263D0*ALZT+0.004938D0*ALZT**2+
     &(21.497D0+0.05794D0*ALZT-0.006584D0*ALZT**2)*SW2-21.0799D0*SW2**2)
      IF(IAFMT.EQ.1) CA2=CA2L
      IF(IAFMT.EQ.2) CA2=CA2L+CA2C
      IF(IAFMT.EQ.3) CA2=CA2L+CA2C+CA2I
      AFMT3=CA1*ALST/PI+CA2*(ALST/PI)**2
*
      END 

      FUNCTION TBQCDR(ALST,AMT2,AMZ2,SW2)
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      COMMON/CDZSCT/ISCRE,ISCAL,IAFMT,IFACR,IFACT,IHIGS,IEWLC,ICZAK
     &             ,IHIG2,IALE2,IGFER      
*
* NUMERICAL CONSTANTS
*
      PI=ATAN(1D0)*4D0
      CW2=1D0-SW2
*
      ALZT=LOG(AMZ2/AMT2)
      CA2C=AMZ2/AMT2*(
     &             +(-11.3184D0-0.62630D0*ALZT+0.4259D0*ALZT**2)*SW2
     &             +(+22.6367D0+1.25270D0*ALZT-0.8519D0*ALZT**2)*SW2)
      CA2I=(AMZ2/AMT2)**2*(
     &(-16.01860-0.02897D0*ALZT+0.003292D0*ALZT**2)*SW2+10.54D0*SW2**2+
     &(21.497D0+0.05794D0*ALZT-0.006584D0*ALZT**2)*SW2-21.0799D0*SW2**2)
      IF(IAFMT.EQ.1) CA2=0D0
      IF(IAFMT.EQ.2) CA2=CA2C
      IF(IAFMT.EQ.3) CA2=CA2C+CA2I
      TBQCDR=-CW2/SW2*CA2*(ALST/PI)**2
*
      END
 
      FUNCTION DALPHL(IALE2,S)
*
* LEPTONIC \Delta\alpha up to three loops
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      COMMON/CDZCON/PI,PI2,D2,D3,ALFAI,AL4PI,AL2PI,AL1PI
      COMMON/CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
      DIMENSION PI0(3),PI1(3),PI2A(3),PI2L(3,3),PI2F(3),PI2H(3)
     &                       ,PI1L(3),PI2B(3,3),PI2S(3),PISU(3)
*
      D5=1.0369277551434D0
      ALPHA0=ALFAI
      CEILER=.577216D0
C...  CONST=0D0 ! that was a bug. Correction fixed in V.6.44.
C... A.A.: for s-channel:
      CONST=D2
*
      DO I=1,3
      AML2=(AML(2*I))**2
      RML2=AML2/S
      ALG1=LOG(S/AML2)
      PI0(I)=1D0/ALFAI/PI*(1D0/9+1D0/3*(1D0+2D0*RML2)
     &      *DREAL(XI0(AML2,-S,AML2,AML2)))
      PI1(I) =(1D0/ALFAI/PI)**2*((1D0+12D0*RML2)/4*ALG1+D3-5D0/24)
      PI1L(I)=(1D0/ALFAI/PI)**2*( 1D0/4*ALG1)
      PI2A(I)=(1D0/ALFAI/PI)**3*(-1D0/4*(-121D0/48
     &       +(-5D0+8D0*LOG(2D0))*D2-99D0/16*D3+10D0*D5+1D0/8*ALG1))
      ENDDO
*
      DO I=1,3
      SUM=0D0
      DO J=1,3
      ALG1=LOG(S/(AML(2*I))**2)
      ALG2=LOG(S/(AML(2*J))**2)
      PI2F(I)=(1D0/ALFAI/PI)**3*(-1D0/4*(-307D0/216-8D0/3*D2
     &       +545D0/144*D3+(11D0/6-4D0/3*D3)*ALG1-1D0/6*ALG1**2+CONST))
      PI2L(I,J)=(1D0/ALFAI/PI)**3*(-1D0/4*(-116D0/27+4D0/3*D2
     &         +38D0/9*D3+14D0/9*ALG1+(5D0/18-4D0/3*D3)*ALG2
     &         +1D0/6*ALG1**2-1D0/3*ALG1*ALG2+CONST))
      PI2H(J)=(1D0/ALFAI/PI)**3*(-1D0/4*(-37D0/6+38D0/9*D3
     &       +(11D0/6-4D0/3*D3)*ALG2-1D0/6*ALG2**2+CONST))
      IF(I.LT.J) THEN
        PI2B(I,J)=PI2H(J)  
      ELSEIF(I.EQ.J) THEN 
        PI2B(I,J)=PI2F(I)  
      ELSEIF(I.GT.J) THEN
        PI2B(I,J)=PI2L(I,J)
      ENDIF
      SUM=SUM+PI2B(I,J)
      ENDDO
      PI2S(I)=PI2A(I)+SUM
*
      IF(IALE2.EQ.0) THEN
        PISU(I)=PI0(I)+PI1L(I)
      ELSEIF(IALE2.EQ.1) THEN
        PISU(I)=PI0(I)
      ELSEIF(IALE2.EQ.2) THEN
        PISU(I)=PI0(I)+PI1(I)
      ELSEIF(IALE2.EQ.3) THEN
        PISU(I)=PI0(I)+PI1(I)+PI2S(I)
      ENDIF
*
      ENDDO
*
      DALPHL=PISU(1)+PISU(2)+PISU(3)
*
      END

      double precision function dalh5(s,argmz)
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      COMMON/CDZFLG/IHVP,IAMT4,IQCD,IMOMS,IMASS,IALEM,IMASK,IBARB,IFTJR
      COMMON/CDZSCT/ISCRE,ISCAL,IAFMT,IFACR,IFACT,IHIGS,IEWLC,ICZAK
     &             ,IHIG2,IALE2,IGFER      
*
* for back-compatibility only
*
       IF(IALE2.EQ.0) THEN
        e=-s/dabs(s)*dsqrt(dabs(s))
       ELSE
        e=+s/dabs(s)*dsqrt(dabs(s))
       ENDIF
* st2=0.2322 is the reference value
       st2= 0.2322d0
       IF(IHVP.EQ.1) THEN
         CALL hadr5(e,argmz,st2,der,errder,deg,errdeg)
       ELSEIF(IHVP.EQ.4) THEN
         CALL dhadr5n(e,st2,der,errder,deg,errdeg)
       ELSEIF(IHVP.EQ.5) THEN
         CALL dhadr5x(e,st2,der,errdersta,errdersys,deg,errdegsta,
     &              errdegsys)
       ENDIF
       IF ((e.LT.4.e1).AND.(e.GT.0.25)) THEN
        dalh5=0d0
       ELSE
          dalh5=der
*          dalh5=der-0.0001
*          dalh5=der+0.0001
       ENDIF
*
      END
