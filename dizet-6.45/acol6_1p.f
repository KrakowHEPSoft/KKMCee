*
* package acol.f 
* v.61  (26 June  1999)
* v.603 (01 April 1999)
* v.601 (08 March 1999) developed to v.603 without numerical changes in any 
*                       relevant digits 
* v.600 contained a bug yet
*
* Authors of the package:
* P. Christova, M. Jack, T. Riemann
* emails:   penka@ifh.de
*            jack@ifh.de
*    tord.riemann@ifh.de
* see:
* DESY 99-015 [hep-ph/9902408] and material to be published
* DESY 99-037: contains some numerical examples and comparisons of old and new
* treatments of acollinearity cut
*-------------------------------------------------------------------------
* zfitter with acol6_1.f is enabled:
*- to also calculate angular distributions in \cos\theta
*  with or without acollinearity/fermion energy cuts 
*  (ICUT=2,3), replacing old code from v.5_20 and earlier (ICUT=0)
*
*- also available now: the final, optimized analytical expressions    
*
*-------------------------------------------------------------------------
* zfitter with acol6_03.f is enabled:
*- to calculate with short expressions for acollinearity and fermion
*  energy cuts *
*  this corresponds to SFAST and is a new feature (ICUT=2)
*- to calculate with acollinearity and fermion energy and acceptance
*  cuts; this replaces the old acol coding (ICUT=3 now, instead of ICUT=0)
*
* have in mind:
*- not yet available: the correct angular distribution with
*  acollinearity cut (ZCUTCOS part)
*- not yet available: the final, optimized analytical expressions 
*- the code calculates final state O(alpha) corrections; they are not
*  made accessible   
*
* now follow 
* NEW EXTRA SUBROUTINES FOR A CALCULATION WITH ACOLLINEARITY/ENERGY CUTS :
*-------------------------------------------------------------------------
*
* for backward compatibility we introduced a flag
* ifunfin (with its own common block) with 2 settings:
* ifunfin = 0: bugs in FUNFIN and in FAL2 are NOT corrected
* ifunfin = 0: these bugs are corrected
* we FIXED the flag to value 1 inside function FAL2 and subroutine FUNFIN.
*
*==========================================================================*
*                                                                          *
* Changes for v. 6_1 (07-05-99):                                           *
*                                                                          *
* Main changes in comparison to older version 6_05:                        *
*                                                                          *
* - Inclusion of corrected hard contributions to initial state and         *
*   interference contributions to angular distributions in \cos\theta      *
*                                                                          *
* - Modified function HARD                                                 *
*                                                                          *   
* - Replacement of subroutine HCUT by 2 new subroutines:                   *
*                                                                          *
*   * HCUTACOL: complete results of hard flux functions                    * 
*               to angular distributions with all cuts                     *
*   * HCUT:     results of hard flux functions                             *
*               to angular distributions without                           *
*               acollinearity/minimal energies' cut                        *
*                                                                          *
* - New subroutines called by HCUTACOL/HCUT:                               *
*   SAVAR,DLOGCP,DLOGCM (logarithms,polynomials)                           * 
*   DRADTINI  <--- DRTFUN0,DRTFUN1      (sym. flux functions in c)         *
*   DRADFBINI <--- DRFBFUN0,DRFBFUN1    (anti-sym. flux functions in c)    *
*                                                                          *
*==========================================================================*
*                                                                          *
* Changes for v. 6_05 (07-05-99):                                          *
*                                                                          *
* Main changes in comparison to older versions 6_04 and 6_03:              *
*                                                                          *
* - Simplification of expressions in acol6_05.f:                           *
*                                                                          *
*   much shorter expressions for O(alpha) hard interference                *
*   corrections for branch with acollinearity and acceptance cut (ICUT=3); *
*   former subr. HARDINTACOL in acol6_05.f removed and hard intf.          *
*   expressions with cuts now treated parallel to initial-state terms      *
*   in subr. RADTIN and RADFBIN and functions called there;                *
*   hard final-state corrections with cuts (not used yet) treated          * 
*   separately in subr. RADTFIN and RADFBFIN                               *
*                                                                          *
* - Removing of some Fortran bugs in acoll. branch (2 undefined logarithms *
*   [1 bug for ICUT=2 had already been present in acol6_04.f],             *
*   one missing exponent in acol6_05.f)                                    *
*                                                                          *
* - Just some optical 'improvements' (change of the order of the           *
*   different subroutines in acol6_05.f; subr. which are not used          *
*   at the end of acol6_05.f; shorter names for some subr. in acol6_05.f)  *
*==========================================================================*
 
C ZUCUTS 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C BORN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C FUNFIN,FAL2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C SETCUT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C ZCUT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C SFAST
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C SCUT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C FCROS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C FASYM
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C SHARD, AHARD
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C ZANCUT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C COSCUT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C HARD
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C!!!!!!!!!!!!!!!!!!!!!!!!!!C
C                          C
C       PACKAGE I:         C
C                          C
C!!!!!!!!!!!!!!!!!!!!!!!!!!C

C====================================================================
C 
C     Package containing hard flux functions (radiators) 
C   to total cross section and forward backward cross sections
C
C          $\sigma_{T,FB}(R,A(R),c)$ 
C 
C   with cuts on one of the final state fermions'
C   minimal scattering angles $\cos\vartheta$, 
C   on their maximal acollinearity angle $\xi^{max}$, 
C   on their minimal energies $E_{min}$, 
C   and/or on their invariant mass squared $s'$ 
C   
C   2 possible options in the code:
C
C   a. general case with all cuts (see above);  
C   
C   b. special case without acceptance cut ($c=1$);
C      only cuts on acollinearity/energies and on $s'$
C      --> short formulae, no or respectively only two cases for 
C          phase space splitting 
C 
C====================================================================
C      ** Introduced by M. JACK on 03/03/99 05:00pm **
C====================================================================



C-------------------------------------------------
C Subroutines preparing phase space boundaries 
C and defining logarithms and polynomials 
C-------------------------------------------------

      SUBROUTINE PHASEREGC
*=======================================================================
* Introduced by M. JACK on 03/03/99 05:00pm
C-----------------------------------------------
C Distinction of phase space regions depending
C on \cos\theta or cut-off C respeCtively
C-----------------------------------------------

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON /SVAR  / S
      COMMON /MASSZ / AME,AMF,AME2,AMF2
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /DELVAR/ DEL,DELM,DIF,DIF3,ALD,ALDM,ALDPL,ALDMI
     +               ,RD,RDI,RDI2,RDI3,RDM,RDMI,RDMI2,RDMI3
     +               ,ALRD,ALRDM
      COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM
C
      COMMON /VAR/ A,V,RP,V2,RP2,AP,AM,R1PI,R1PI2,R1PI3
      COMMON/INDEX/ MF,IIF,JF,KF
C
      XC1 = 0D0
      XC2 = C
C
      V = 1D0-R
      RP = 1D0+R
      DIFA = DIF
      BT2 = 1D0-4D0*AME2/S
      BT  = SQRT(BT2)

C------------------------------------------------------------
C Polynomials and "zero-conditions" :

      G1A = RP-V*DIFA
      G2A = RP+V*DIFA
      Y1A = V-RP*DIFA
      Y2A = V+RP*DIFA
      
      XCM1 = Y1A/G1A
      XCP1 = Y2A/G2A
C     
C      XCM1 = MIN(XCM1,1D0)
C      XCM1 = MAX(XCM1,-1D0)
C      XCP1 = MIN(XCP1,1D0)
C      XCP1 = MAX(XCP1,-1D0)
C     
      XCM2 = -XCM1
      XCP2 = -XCP1   

C----------------------------------------------------------
C Abbreviation of distinction between different cases in
C phase space depending on value of s' by introduction 
C of following indiCes I,J,K : 
C----------------------------------------------------------

C---------------------------------------------------
C 1. XC2=C >= 0:
C---------------------------------------------------
      IF (XC2.GT.0D0) THEN
         
         IF (XCM2.LT.XCM1) THEN
            
            KF = 1
            IF (XC2.GT.XCP1) THEN
               IIF = 1 
               JF = -1
            ELSE
               IF ((XC2.LE.XCP1)
     &              .AND.(XC2.GT.XCM1)) THEN               
                  IIF = 1 
                  JF = 0
               ELSE
                  IIF = 1 
                  JF = 1
               ENDIF
            ENDIF
            
         ELSE
            
            KF = 0
            IF (XC2.GT.XCP1) THEN
               IIF = 1 
               JF = -1
            ELSE
               IF ((XC2.LE.XCP1)
     &              .AND.(XC2.GT.XCM2)) THEN               
                  IIF = 1 
                  JF = 0
               ELSE
                  IIF = 0 
                  JF = 0
               ENDIF
            ENDIF
            
         ENDIF
         
      ELSE
         
C---------------------------------------------------
C 2. XC2=C < 0:
C---------------------------------------------------
     
         IF (XCM2.LT.XCM1) THEN
            
            KF = 1
            IF (XC2.LT.XCP2) THEN
               IIF = -1 
               JF = 1
            ELSE
               IF ((XC2.GE.XCP2)
     &              .AND.(XC2.LT.XCM2)) THEN               
                  IIF = 0 
                  JF = 1
               ELSE
                  IIF = 1 
                  JF = 1
               ENDIF
            ENDIF

         ELSE
            
            KF = 0
            IF (XC2.LT.XCP2) THEN
               IIF = -1 
               JF = 1
            ELSE
               IF ((XC2.GE.XCP2)
     &              .AND.(XC2.LT.XCM1)) THEN               
                  IIF = 0 
                  JF = 1
               ELSE
                  IIF = 0 
                  JF = 0
               ENDIF
            ENDIF
            
         ENDIF
         
      ENDIF

C------------------------------------------------------------
C       WRITE(*,*) '...........................',
C     &  '...........................'
C       WRITE(*,*)
C       WRITE(*,*) 'I1,J1,K1  = ',IIF,JF,KF
C       WRITE(*,*) 'XCP2,XCM2 = ',XCP2,XCM2               
C       WRITE(*,*) 'XCM1,XCP1 = ',XCM1,XCP1
C       WRITE(*,*) 'C, XC2, XC1 = ',C,XC2,XC1
C       WRITE(*,*) 'A=', DIFA
C       WRITE(*,*) 'V,R,RP=', V,R,RP
C       WRITE(*,*) '...........................',
C     &  '...........................' 

      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  
      SUBROUTINE SAVAR
*====================================================================
** Introduced by M. JACK on 22/06/99 07:00pm
*
C Correction by M. Jack (07/05/99): 
C Logarithm XLNG12Q had not been defined.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

*---------------------------------------------------------------
      COMMON /MASSZ / AME,AMF,AME2,AMF2
      COMMON /SVAR  / S
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /DELVAR/ DEL,DELM,DIF,DIF3,ALD,ALDM,ALDPL,ALDMI
     +               ,RD,RDI,RDI2,RDI3,RDM,RDMI,RDMI2,RDMI3
     +               ,ALRD,ALRDM
      COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM
*---------------------------------------------------------------

       COMMON /VAR/ A,V,RP,V2,RP2,AP,AM,R1PI,R1PI2,R1PI3
       COMMON /XMASSES/ BETA,BETA2
       COMMON /XLIMITS/ XC1,XC2,epsc,epsccut,epsr,
     & XCM1,XCM2,XCP1,XCP2,epsl,epsm,ep,epsa,epsa1,epsa2 

C--------------------------------------------------

       COMMON /XLOGSRM/ XLNRP,XLNR,XLNM,XLLPR2,XLN2

       COMMON /XLOGS2/ XLNG1A,XLNG2A,XLNY1A,XLNY2A,
     & XLNY12P,XLNY12Q,XLNG12P,XLNG12Q 
       COMMON /XLOGSA/ XLNAP,XLNAM,XLNA2,XLNAQ
       COMMON /ALOGS/  ALNA2,ALNAM,ALNAP
       COMMON /YLOGSA/ YLNAS,YLNAD

C--------------------------------------------------

       COMMON /POLYNOMIALS/ G1C,G2C,Y1C,Y2C,
     & G1A,G2A,Y1A,Y2A,YG1CP,YG1CM,YG2CP,YG2CM

       COMMON /XLOGS1/ XLNG1C,XLNG2C,XLNG1C2,XLNG2C2,
     & XLNYG1CP,XLNYG1CM,XLNYG2CP,XLNYG2CM

       COMMON /XLOGS3/ CLNC1M,CLNC1P,CLNC2M,CLNC2P,
     & YLNC1S,YLNC1D,YLNC2S,YLNC2D

       COMMON /CLOGS/  XLNCP,XLNCM,XLNC2,XLNCQ,CLNC2

C--------------------------------------------------

       COMMON /DLOGC1/ DLNYQ1,DLNYP1,DLNCP1,DLNCM1,DLNCQ1,
     &                 DLNGAMC1
       COMMON /DLOGC2/ DLNYQ2,DLNYP2,DLNCP2,DLNCM2,DLNCQ2,
     &                 DLNGAMC2 
       COMMON /DPOLY/  GAMC,GAMC2,GAMC3,GAMC4
       COMMON /DLOGC/  DLNYQ,DLNYP,DLNCP,DLNCM,DLNCQ,
     &                 DLNGAMC,DLNC2
C--------------------------------------------------

       COMMON /INDEXA/ INDA,INDM
       COMMON /INDEXACOL/ INDACOL

C--------------------------------------------------

       V = 1D0-R
       RP = 1D0+R
C
       XAMU2 = 4D0*AMF2/S
       XAME2 = 4D0*AME2/S
       XAM   = SQRT(XAME2)
       BETA2 = 1D0-XAME2
       BETA  = SQRT(BETA2)

C------------------------------------------------------
C       IF (XAM.LT.1D-6) THEN
C          BETA = 1D0-XAME2/2D0-XAME2*XAME2/8D0 
C     &    -XAME2*XAME2*XAME2/16D0 
C       ENDIF
C------------------------------------------------------

       A  = DIF
       AP = 1D0+A
       AM = 1D0-A

       IF ((INDA.EQ.1).AND.(XAMU2/R.LT.1D-5)) THEN
         A  = 1D0-XAMU2/R/2D0-(XAMU2/R)**2/8D0
         AP = 2D0-XAMU2/R/2D0-(XAMU2/R)**2/8D0
         AM = XAMU2/R/2D0+(XAMU2/R)**2/8D0   
       ENDIF
C     
       CB = BETA*C

       R2 = R*R
       R3 = R*R2
       R4 = R*R3
C     
       V2 = V*V
       RP2 = RP*RP
       RF = 1D0+R2
       R1PI = 1D0/RP
       R1PI2 = R1PI*R1PI
       R1PI3 = R1PI*R1PI2
C
       AD = AP*AM
       A2 = A*A
       A3 = A*A2
       A4 = A*A3
                
C------------------------------------------------------
C Polynomials and "zero-conditions" :

       G1A = RP-V*A
       G2A = RP+V*A
       Y1A = V-RP*A
       Y2A = V+RP*A
C
       XCM1 = Y1A/G1A/BETA
       XCP1 = Y2A/G2A/BETA

       XCM1 = MIN(XCM1,1D0)
       XCM1 = MAX(XCM1,-1D0)
       XCP1 = MIN(XCP1,1D0)
       XCP1 = MAX(XCP1,-1D0)
       XCM2 = -XCM1
       XCP2 = -XCP1   
C
       G1C = RP-V*CB
       G2C = RP+V*CB
       Y1C = V-RP*CB
       Y2C = V+RP*CB
C
       YG1CM = Y1C-A*G1C
       YG1CP = Y1C+A*G1C
       YG2CM = Y2C-A*G2C
       YG2CP = Y2C+A*G2C

C------------------------------------------------------------
C Basic logarithms :

       XLN2   = LOG(2D0)
       XLNM   = LOG(S/AME2)
       XLNR   = LOG(MAX(1.D-32,ABS(R)))
       XLLPR2 = XLNM-2D0*XLNR
       XLNRP  = LOG(RP)
C
       XLNAP  = LOG(MAX(1.D-32,ABS(AP)))
       XLNAM  = LOG(MAX(1.D-32,ABS(AM)))
C
       XLNG1A = LOG(MAX(1.D-32,ABS(G1A)))
       XLNG2A = LOG(MAX(1.D-32,ABS(G2A)))
       XLNG12P = LOG(MAX(1.D-32,ABS(G1A*G2A)))
       XLNG12Q = LOG(MAX(1.D-32,ABS(G2A/G1A)))
C
       XLNY1A = LOG(MAX(1.D-32,ABS(Y1A)))
       XLNY2A = LOG(MAX(1.D-32,ABS(Y2A))) 
       XLNY12P = LOG(MAX(1.D-32,ABS(Y1A*Y2A)))
       XLNY12Q = LOG(MAX(1.D-32,ABS(Y2A/Y1A)))
C
       XLNCM = LOG(MAX(1.D-32,ABS(1D0-CB)))
       XLNCP = LOG(MAX(1.D-32,ABS(1D0+CB)))
       XLNCM = MAX(-XLNM+XLN2,XLNCM)
       XLNCP = MAX(-XLNM+XLN2,XLNCP)
C
       XLNG1C = LOG(MAX(1.D-32,ABS(G1C)))
       XLNG2C = LOG(MAX(1.D-32,ABS(G2C)))
       XLNG1C2 = LOG(MAX(1.D-32,ABS(G1C/2D0)))
       XLNG2C2 = LOG(MAX(1.D-32,ABS(G2C/2D0)))
C
       XLNYG1CM = LOG(MAX(1.D-32,ABS(YG1CM)))
       XLNYG1CP = LOG(MAX(1.D-32,ABS(YG1CP)))
       XLNYG2CM = LOG(MAX(1.D-32,ABS(YG2CM)))
       XLNYG2CP = LOG(MAX(1.D-32,ABS(YG2CP)))

C------------------------------------------------------------
C Combined logarithms:

       XLNCQ  = LOG(MAX(1.D-32,ABS((1D0+CB)/(1D0-CB))))
       XLNC2  = LOG(MAX(1.D-32,ABS((1D0+CB)*(1D0-CB))))
       CLNC1M = (1D0+CB)*XLNCP 
       CLNC1P = (1D0-CB)*XLNCM 
       CLNC2M = (1D0-CB)*XLNCM
       CLNC2P = (1D0+CB)*XLNCP
       CLNC2  = (1D0-CB)*(1D0+CB)*XLNC2
C
       YLNC1S = YG1CP*XLNYG1CP+YG1CM*XLNYG1CM
       YLNC1D = YG1CP*XLNYG1CP-YG1CM*XLNYG1CM
       YLNC2S = YG2CP*XLNYG2CP+YG2CM*XLNYG2CM
       YLNC2D = YG2CP*XLNYG2CP-YG2CM*XLNYG2CM
C     
       XLNAQ = LOG(MAX(1.D-32,ABS(AP/AM)))
       XLNA2 = LOG(MAX(1.D-32,ABS(AP*AM)))
       ALNA2 = AP*AM*XLNA2
       ALNAP = AP*XLNAP
       ALNAM = AM*XLNAM              
       YLNAS = Y2A*XLNY2A+Y1A*XLNY1A
       YLNAD = Y2A*XLNY2A-Y1A*XLNY1A


C------------------------------------------------------------
C New logarithms for angular distributions:

       CCP = ABS(1D0+CB)
       CCM = ABS(1D0-CB)

       CCP = MAX(CCP,2D0*AME2/S)       
       CCM = MAX(CCM,2D0*AME2/S)  

       DLNCCP = LOG(CCP)
       DLNCCM = LOG(CCM)

       DLNCCM = MAX(-XLNM+XLN2,DLNCCM)
       DLNCCP = MAX(-XLNM+XLN2,DLNCCP)
        
       DLNCCQ = DLNCCP-DLNCCM
       DLNC2  = DLNCCP+DLNCCM
C
       DLNCP1 = DLNCCM
       DLNCM1 = DLNCCP
       DLNCQ1 =-DLNCCQ

       DLNCP2 = DLNCCP
       DLNCM2 = DLNCCM
       DLNCQ2 = DLNCCQ
            
       DLNGAMC1 = LOG(MAX(1.D-32,ABS(G1C/2D0)))
       DLNGAMC2 = LOG(MAX(1.D-32,ABS(G2C/2D0))) 

C       CCP1 = ABS((1D0-R)/2D0/R*(1D0+CB))
C       CCM1 = ABS((1D0-R)/2D0/R*(1D0-CB))
C       CCP2 = ABS((1D0-R)/2D0*(1D0+CB))
C       CCM2 = ABS((1D0-R)/2D0*(1D0-CB))
C
C       IF (CCP1.LT.1D-4) THEN 
C         DLNGAMC2 = ALR + CCP1 - CCP1**2/2D0
C         DLNGAMC1 = - CCP2 - CCP2**2/2D0
C       ELSE
C       IF (CCM1.LT.1D-4) THEN
C         DLNGAMC2 = - CCM2 - CCM2**2/2D0
C         DLNGAMC1 = ALR + CCM1 - CCM1**2/2D0
C       ENDIF
C       ENDIF

       ALMIN1=1D0/2D0*(DLNC2-XLNM)-DLNGAMC1+ALR+XLN2
       ALMIN2=1D0/2D0*(DLNC2-XLNM)-DLNGAMC2+ALR+XLN2
C
       DLNYG1CP = LOG(MAX(1.D-32,ABS(YG1CP/2D0))) 
       DLNYG1CM = LOG(MAX(1.D-32,ABS(YG1CM/2D0))) 
       DLNYG2CP = LOG(MAX(1.D-32,ABS(YG2CP/2D0))) 
       DLNYG2CM = LOG(MAX(1.D-32,ABS(YG2CM/2D0))) 

       DLNYG1CP = MAX(ALMIN1, DLNYG1CP)  
       DLNYG1CM = MAX(ALMIN1, DLNYG1CM)  
       DLNYG2CP = MAX(ALMIN2, DLNYG2CP)  
       DLNYG2CM = MAX(ALMIN2, DLNYG2CM)  

       IF (CCM.LT.1D-4) THEN
          DLNYG1CP = XLNAM + ALR
          DLNYG1CM = XLNAP + ALR
          DLNYG2CP = XLNAP
          DLNYG2CM = XLNAM      
       ELSE
       IF (CCP.LT.1D-4) THEN
          DLNYG2CP = XLNAM + ALR
          DLNYG2CM = XLNAP + ALR
          DLNYG1CP = XLNAP
          DLNYG1CM = XLNAM   
       ENDIF
       ENDIF

       DLNYP1 = DLNYG1CP+DLNYG1CM
       DLNYP2 = DLNYG2CP+DLNYG2CM
       DLNYP1 = DLNYP1-DLNC2+2D0*DLNGAMC1-2D0*ALR+XLNM
       DLNYP2 = DLNYP2-DLNC2+2D0*DLNGAMC2-2D0*ALR+XLNM

       DLNYQ1 = DLNYG1CP-DLNYG1CM
       DLNYQ2 = DLNYG2CP-DLNYG2CM

C------------------------------------------------------------
C 
C Special case: A=\sqrt{1-4 m_f^2/s'} --> 1
C
       IF ((INDA.EQ.1).AND.(XAMU2/R.LT.1D-5)) THEN
         XLNAP = XLN2 - XAMU2/R/4D0
         XLNAM = XLN2 - XLNR + LOG(XAMU2/4D0)
         XLNAQ = XLNAP-XLNAM
         XLNA2 = XLNAP+XLNAM
         ALNAP = AP*XLNAP
         ALNAM = AM*XLNAM
         ALNA2 = XAMU2/R*XLNA2
C
         XLNG12Q = -ALR
C
C         DLNYP1 = 2D0*DLNGAMC1-ALR+XLNM
C         DLNYP2 = 2D0*DLNGAMC2-ALR+XLNM
C         DLNYQ1 = -DLNCCQ-ALR 
C         DLNYQ2 =  DLNCCQ-ALR
C
       ENDIF
C------------------------------------------------------------


C------------------------------------------------------------
C       WRITE(*,*) '...........................',
C     &  '...........................'
C       WRITE(*,*)
C       WRITE(*,*) 'I1,J1,K1  = ',IIF,JF,KF
C       WRITE(*,*) 'XCP2,XCM2 = ',XCP2,XCM2               
C       WRITE(*,*) 'XCM1,XCP1 = ',XCM1,XCP1
C       WRITE(*,*) 'A,AP,AM=', A,AP,AM
C       WRITE(*,*) 'V,R,RP=', V,R,RP
C       WRITE(*,*) 'C,BETA*C=', C,CB
C       WRITE(*,*) 'S,AME,AMF,XLNM=',S,AME,AMF,XLNM
C       WRITE(*,*) 'AME2,AMF2,TE,TF=',AME2,AMF2,TE,TF
C       WRITE(*,*) 'XLNR=',XLNR
C       WRITE(*,*) 'XLNRP=',XLNRP
C       WRITE(*,*) 'G1C=', G1C
C       WRITE(*,*) 'G2C=', G2C
C       WRITE(*,*) 'YG1CP=', YG1CP
C       WRITE(*,*) 'YG1CM=', YG1CM
C       WRITE(*,*) 'YG2CP=', YG2CP
C       WRITE(*,*) 'YG2CM=', YG2CM
C       WRITE(*,*) 'LNG1C  =', XLNG1C
C       WRITE(*,*) 'LNG2C  =', XLNG2C
C       WRITE(*,*) 'LnYG1CP=', XLNYG1CP
C       WRITE(*,*) 'LnYG1CM=', XLNYG1CM
C       WRITE(*,*) 'LnYG2CP=', XLNYG2CP
C       WRITE(*,*) 'LnYG2CM=', XLNYG2CM
C       WRITE(*,*) 'YLNC1S=',YLNC1S 
C       WRITE(*,*) 'YLNC1D=',YLNC1D 
C       WRITE(*,*) 'YLNC2S=',YLNC2S 
C       WRITE(*,*) 'YLNC2D=',YLNC2D
C
C       WRITE(*,*) 'LNCP  =', XLNCP
C       WRITE(*,*) 'LNCM  =', XLNCM
C       WRITE(*,*) 'XLNCQ =',XLNCQ  
C       WRITE(*,*) 'XLNC2 =',XLNC2  
C       WRITE(*,*) 'CLNC2 =',CLNC2  
C       WRITE(*,*) 'CLNC1M=',CLNC1M 
C       WRITE(*,*) 'CLNC1P=',CLNC1P 
C       WRITE(*,*) 'CLNC2M=',CLNC2M 
C       WRITE(*,*) 'CLNC2P=',CLNC2P
C       WRITE(*,*) 'Y1A=', Y1A
C       WRITE(*,*) 'G1A=', G1A
C       WRITE(*,*) 'Y2A=', Y2A
C       WRITE(*,*) 'G2A=', G2A
C       WRITE(*,*) 'XLNA2=',XLNA2 
C       WRITE(*,*) 'XLNAP=',XLNAP 
C       WRITE(*,*) 'XLNAM=',XLNAM 
C       WRITE(*,*) 'XLNAQ=',XLNAQ 
C       WRITE(*,*) 'LnY1A=', XLNY1A
C       WRITE(*,*) 'LnY2A=', XLNY2A
C       WRITE(*,*) 'LNG1A=', XLNG1A
C       WRITE(*,*) 'LNG2A=', XLNG2A 
C       WRITE(*,*) 'YLNAS=',YLNAS 
C       WRITE(*,*) 'YLNAD=',YLNAD 
C       WRITE(*,*) 'XLNAQ=',XLNAQ
C       WRITE(*,*) 'ALNA2=',ALNA2
C       WRITE(*,*) 'ALNAP=',ALNAP
C       WRITE(*,*) 'ALNAM=',ALNAM
C
C       WRITE(*,*) 'DLNC2   =', DLNC2
C       WRITE(*,*) 'DLNCCQ  =', DLNCCQ
C       WRITE(*,*) 'DLNYP1  =', DLNYP1
C       WRITE(*,*) 'DLNYQ1  =', DLNYQ1
C       WRITE(*,*) 'DLNYP2  =', DLNYP2
C       WRITE(*,*) 'DLNYQ2  =', DLNYQ2
C       WRITE(*,*) 'DLNGAMC1 =',DLNGAMC1 
C       WRITE(*,*) 'DLNGAMC2 =',DLNGAMC2
C       WRITE(*,*) 'ALR,XLN2,XLNM =',ALR,XLN2,XLNM
C       WRITE(*,*) '...........................',
C     &  '...........................' 

      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


 
      SUBROUTINE SLOGSP
C=======================================================================
** Introduced by M. JACK on 03/03/99 05:00pm

       IMPLICIT DOUBLE PRECISION (a-h,o-z)

       COMMON /POLYNOMIALS/ G1C,G2C,Y1C,Y2C,
     & G1A,G2A,Y1A,Y2A,YG1CP,YG1CM,YG2CP,YG2CM
       COMMON /XLOGS1/ XLNG1C,XLNG2C,XLNG1C2,XLNG2C2,
     & XLNYG1CP,XLNYG1CM,XLNYG2CP,XLNYG2CM
       COMMON /XLOGS3/ CLNC1M,CLNC1P,CLNC2M,CLNC2P,
     & YLNC1S,YLNC1D,YLNC2S,YLNC2D

       COMMON /POLYC/  GC,YC,GC2,YC2,GC3,YC3
       COMMON /CLOGS/  XLNCP,XLNCM,XLNC2,XLNCQ,CLNC2
       COMMON /CLOGSC/ CLNCM,CLNCP
       COMMON /XLOGSC/ XLNCCP,XLNCCM,XLNCCQ,XLNGC
       COMMON /YLOGSC/ YLNCS,YLNCD
C
       GC = G2C
       YC = Y2C
       GC2 = GC*GC
       YC2 = YC*YC
       GC3 = GC*GC2
       YC3 = YC*YC2
       XLNGC = XLNG2C
       XLNCCP = XLNCP
       XLNCCM = XLNCM
       XLNCCQ = XLNCP-XLNCM
       CLNCM = CLNC2M 
       CLNCP = CLNC2P        
       YLNCS = YLNC2S
       YLNCD = YLNC2D

       END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


 
      SUBROUTINE SLOGSM
C=======================================================================
** Introduced by M. JACK on 03/03/99 05:00pm

       IMPLICIT DOUBLE PRECISION (a-h,o-z)

       COMMON /POLYNOMIALS/ G1C,G2C,Y1C,Y2C,
     & G1A,G2A,Y1A,Y2A,YG1CP,YG1CM,YG2CP,YG2CM
       COMMON /XLOGS1/ XLNG1C,XLNG2C,XLNG1C2,XLNG2C2,
     & XLNYG1CP,XLNYG1CM,XLNYG2CP,XLNYG2CM
       COMMON /XLOGS3/ CLNC1M,CLNC1P,CLNC2M,CLNC2P,
     & YLNC1S,YLNC1D,YLNC2S,YLNC2D

       COMMON /POLYC/  GC,YC,GC2,YC2,GC3,YC3
       COMMON /CLOGS/  XLNCP,XLNCM,XLNC2,XLNCQ,CLNC2
       COMMON /CLOGSC/ CLNCM,CLNCP
       COMMON /XLOGSC/ XLNCCP,XLNCCM,XLNCCQ,XLNGC
       COMMON /YLOGSC/ YLNCS,YLNCD
C
       GC = G1C
       YC = Y1C
       GC2 = GC*GC
       YC2 = YC*YC
       GC3 = GC*GC2
       YC3 = YC*YC2
       XLNGC = XLNG1C
       XLNCCP = XLNCM
       XLNCCM = XLNCP
       XLNCCQ = XLNCM-XLNCP
       CLNCM = CLNC1M 
       CLNCP = CLNC1P        
       YLNCS = YLNC1S
       YLNCD = YLNC1D

       END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C-----------------------C
C Total cross-section:  C
C-----------------------C
 
      SUBROUTINE SHACOL(H0,H4,H5)
*     ========== =============
* Introduced by M. JACK on 03/03/99 05:00pm

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /MASSZ / AME,AMF,AME2,AMF2
      COMMON / SVAR/ S
      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /DELVAR/ DEL,DELM,DIF,DIF3,ALD,ALDM,ALDPL,ALDMI
     +               ,RD,RDI,RDI2,RDI3,RDM,RDMI,RDMI2,RDMI3
     +               ,ALRD,ALRDM
      COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
      COMMON /CORINT/ CORINT
C
      COMMON /INDEXA/ INDA,INDM
C------------------------------------------------------------------------
      
      RR=R
C      
      V  = 1D0-R
      RP = 1D0+R
      XMU2 = 4D0*AMF2/S
C            
      IF (INDA.EQ.1) THEN

*******************************************************
* acollinearity cut region I (''s'-cut case'') (A=1)
*******************************************************

         H5 = (1D0+R2)*(TF+ALR)*R1MI

*******************************************************        
* Results with or without final-state mASses:
*******************************************************

         IF (INDM.EQ.0) THEN

            H0 =  4D0/3D0*(1D0+R2)*TE
            H4 = -4D0*R*(1D0+R)*R1MI

         ELSE          

            DIFA2 = 1D0-XMU2/R
            DIFA  = SQRT(DIFA2)

            H0 =  DIFA*((1D0+DIFA2/3D0)*(1D0+R2)*TE
     &           +2D0*XMU2)
            H4 = -4D0*SQRT(ABS(R*(R-XMU2)))*RP*R1MI

         ENDIF

*****************************************************       
* Regularize initial-state term by divergent 
* initial-state soft+vertex Contribution: 
*****************************************************

         IF (R.GE.RCUT) THEN
            H0=(H0-8D0/3D0*TE)*R1MI
            IF (IPHOT2.GE.0) H0=H0+4D0/3D0*SH2(RR,1)
         ELSE
            H0=H0*R1MI
         END IF

*****************************************************       
* Regularize final-state term by divergent 
* final-state soft+vertex Contribution: 
*****************************************************
C         IF (R.GE.RCUT) THEN
C            H5 = H5- ...
C            IF (IPHOT2.GE.0) H5=H5 + ....
C         ELSE
C            H5=H5*R1MI
C         END IF
C
      ELSE

***************************************
* acollinearity cut region II or III 
***************************************
          
         DELA  = DEL
         DELMA = DELM
         DIFA  = DIF

         DIFA2 = DIFA*DIFA
         AP = 2D0*DELMA
         AM = 2D0*DELA
         XLAP = LOG(MAX(1.D-32,ABS(AP)))
         XLAM = LOG(MAX(1.D-32,ABS(AM)))
         XLOGA = LOG(MAX(1.D-32,ABS(AP/AM)))
C          
         H0 = DIFA*((1D0+DIFA2/3D0)*(1D0+R2)*TE
     &        +2D0*(1D0-DIFA2)*R)*R1MI
         H4 = -4D0*DIFA*R*(1D0+R)*R1MI         
         H5 = ((1D0+R2)*XLOGA-2D0*DIFA*XMU2/AP/AM)*R1MI
     &        -DIFA*V

      ENDIF
C
C      IF (INDA.EQ.1) THEN
C      WRITE(*,*) 'INDA,INDM =',INDA,INDM
C      WRITE(*,*) 'DIF =',DIFA
C      WRITE(*,*) 'R,XMU2,AMF,S =',R,XMU2,AMF,S
C      WRITE(*,*) 'H0,H4 =', H0,H4
C      WRITE(*,*) '---------------------------------'
C      ENDIF
C
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


 
      SUBROUTINE SHCUTACOL(H0,H4,H5)
*     ========== ============
* Introduced by M. JACK on 03/03/99 05:00pm,
* modified by M. JACK on 22/04/99 13:00pm 

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /DELVAR/ DEL,DELM,DIF,DIF3,ALD,ALDM,ALDPL,ALDMI
     +               ,RD,RDI,RDI2,RDI3,RDM,RDMI,RDMI2,RDMI3
     +               ,ALRD,ALRDM
      COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM
      COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM

      COMMON /INDEXA/ INDA,INDM
      COMMON /INDEX/ MF,IIF,JF,KF

C------------------------------------------------
C With acceptance cut:
    
      CALL RADTIN(H0M,H4M)
      CALL RADTFIN(H5M)

c      WRITE(*,*) '---------------------------'
c      WRITE(*,*) 'I,J,K =',IIF,JF,KF
c      WRITE(*,*) 'R,C,A =',R,C,DIF
c      WRITE(*,*) 'H5,H0,H4 =',H5M,H0M,H4M      

C------------------------------------------------
C No acceptance cut (C=1) (just for checks):     
C     CALL RTACOL(H0M)
C------------------------------------------------

      H0=H0M
      H4=H4M
      H5=H5M      

C------------------------------------------------
C regularize in s' cut region:
      IF (R.GE.RCUT) THEN
        H0 = H0-2*TE*COPL3*R1MI
        IF (IPHOT2.GE.0) H0=H0+SH2(R,1)*COPL3
      ENDIF
C------------------------------------------------

C------------------------------------------------
C Re-Create old results by M. Bilenky:
C
C      IF (IIF.EQ.JF) RETURN
C      I = IIF
C      J = JF
C      K = KF
C
C------------------------------------------------
C With acceptance cut:
C      CALL SCONST (H0M,H4M,I,J,K)   
C
C No acceptance cut (C=1):
C      CALL SCACOL (H0M,H4M,K)
C------------------------------------------------
C
C      H0=H0-H0M
C      H4=H4-H4M
C
C*------------------------------------------------
C* Insertion in order to subtraCt non-logarithmiC 
C* terms prop. to 1, C, C^2 for Comparisons with 
C* Bilenky's original Code:
C*------------------------------------------------
C
C      CALL SETRC ( 1)
C      CALL SHCUTM (H0MB,H4MB,I)
C      H0B=H0MB
C      H4B=H4MB
C      CALL SETRC (-1)
C      CALL SHCUTM (H0MB,H4MB,J)
C      H0B=H0B-H0MB
C      H4B=H4B-H4MB
C
C      IF (I.EQ.J) RETURN
C
C      CALL SJUMP (H0MB,H4MB,I,J)
C      H0B=H0B+H0MB
C      H4B=H4B+H4MB
C      IB = I
C      JB = J
C      I = IIF
C      J = JF
C      K = KF
C
C       IF (I.EQ.J) RETURN
C       IF (DELC.GT.epsccut) THEN 
C         CALL SCONST (H0M,H4M,I,J)
C       ELSE
C         CALL SCACOL (H0M,H4M)
C       ENDIF 
C       H0=H0-H0M
C       H4=H4-H4M 
C       
C       DELH0 = ABS((H0-H0B)/H0B)
C       IF (DELH0.GE.EPSDEL) THEN
C          WRITE(*,*)
C          WRITE(*,*) '=========== R_T =============='
C          WRITE(*,*) 'Comparison of Codes: Bilenky -- JACK'
C          WRITE(*,*) 'IB,JB = ',IB,JB
C          WRITE(*,*) 'I,J,K = ',IIF,JF,KF
C          WRITE(*,*) 'M.B.: R , RTINI = ',R,H0B
C          WRITE(*,*) 'M.J.: R , RTINI = ',R,H0
C          WRITE(*,*) 'M.J.-M.B.: ',DELH0
C          WRITE(*,*) 'M.J.:  RTNONLOG = ',H0M
C          WRITE(*,*) '====================================',
C     &         '===================================='
C          WRITE(*,*)
C       ENDIF

      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C-----------------------------------------------
C Initial-state and interference contributions :
C-----------------------------------------------
 
      SUBROUTINE RADTIN(RADT,RADTI)
*=======================================================================
* Introduced by M. JACK on 03/03/99 05:00pm,
* modified by M. JACK on 22/04/99 13:00pm 

       IMPLICIT DOUBLE PRECISION (a-h,o-z)

       COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
       COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM

       COMMON /VAR/ A,V,RP,V2,RP2,AP,AM,R1PI,R1PI2,R1PI3
       COMMON /XMASSES/ BETA,BETA2

       COMMON /POLYNOMIALS/ G1C,G2C,Y1C,Y2C,
     & G1A,G2A,Y1A,Y2A,YG1CP,YG1CM,YG2CP,YG2CM

       COMMON/INDEX/ MF,IIF,JF,KF

C-----------------------------------------------
     
       I = IIF
       J = JF
       K = KF
 
       CT = BETA*C
       CT2 = CT*CT

C-----------------------------------------------
C Initial-state terms:
       F10 = -2D0/3D0*(1D0+R2)*R1MI * ALR
       C0  = -2D0*A*CT*V

C Interference terms:       
       F10I = 2D0*(1D0-CT2)*R1MI * ALR 
       C0I  = 2D0*A*CT*V*RP

C-----------------------------------------------

       RADT = 0D0
       RADTI = 0D0

       IF ((I.EQ.0).AND.(J.EQ.0)) THEN

          CALL SLOGSP 
          CALL RTFUN0(CT,F002,F002I)
          CALL SLOGSM 
          CALL RTFUN0(-CT,F001,F001I)

          RT00 =  F002 - F001 + C0
          RADT =  RT00

          RT00I =  F002I - F001I + C0I
          RADTI =  RT00I

       ELSE
       IF ((I.EQ.1).AND.(J.EQ.1)) THEN

          CALL SLOGSP 
          CALL RTFUN1(CT,F112,F112I)
          CALL SLOGSM 
          CALL RTFUN1(-CT,F111,F111I)
          
          RT11 =  F112 - F111 + C0
          RADT =  RT11

          RT11I =  F112I - F111I + C0I
          RADTI =  RT11I

       ELSE
       IF ((I.EQ.1).AND.(J.EQ.0)) THEN

          CALL SLOGSP 
          CALL RTFUN1(CT,F112,F112I)
          CALL SLOGSM 
          CALL RTFUN0(-CT,F001,F001I)

          RT10 =  F112 - F001 + F10 + C0
          RADT =  RT10

          RT10I =  F112I - F001I + F10I + C0I
          RADTI =  RT10I

       ELSE
       IF ((I.EQ.1).AND.(J.EQ.-1)) THEN 

          CALL SLOGSP 
          CALL RTFUN1(CT,F112,F112I)
          CALL SLOGSM 
          CALL RTFUN1(-CT,F111,F111I)
 
          RT11M =  F112 + F111 + C0
          RADT  =  RT11M 
     
          RT11MI =  F112I + F111I + C0I
          RADTI  =  RT11MI 

       ELSE
       IF ((I.EQ.0).AND.(J.EQ.1)) THEN

          CALL SLOGSP 
          CALL RTFUN0(CT,F002,F002I)
          CALL SLOGSM 
          CALL RTFUN1(-CT,F111,F111I)

          RT01 =- ( F111 - F002 + F10 ) + C0
          RADT =  RT01

          RT01I =- ( F111I - F002I + F10I ) + C0I
          RADTI =  RT01I

       ELSE
       IF ((I.EQ.-1).AND.(J.EQ.1)) THEN

          CALL SLOGSP 
          CALL RTFUN1(CT,F112,F112I)
          CALL SLOGSM 
          CALL RTFUN1(-CT,F111,F111I)

          RT1M1 =- ( F111 + F112 ) + C0
          RADT  = RT1M1

          RT1M1I =- ( F111I + F112I ) + C0I
          RADTI  = RT1M1I
       
       ENDIF
       ENDIF
       ENDIF
       ENDIF
       ENDIF
       ENDIF
C
C       RADT = RADT/BETA
C       RADTI = RADTI/BETA
C
       END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


 
      SUBROUTINE RTFUN0(CTH,RT,RTI)
*=======================================================================
* Introduced by M. JACK on 03/03/99 05:00pm,
* modified by M. JACK on 22/04/99 13:00pm 

       IMPLICIT DOUBLE PRECISION (a-h,o-z)

       COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR

       COMMON /VAR/ A,V,RP,V2,RP2,AP,AM,R1PI,R1PI2,R1PI3

       COMMON /XMASSES/ BETA,BETA2

       COMMON /XLOGSRM/ XLNRP,XLNR,XLNM,XLLPR2,XLN2

       COMMON /POLYNOMIALS/ G1C,G2C,Y1C,Y2C,
     & G1A,G2A,Y1A,Y2A,YG1CP,YG1CM,YG2CP,YG2CM

       COMMON /XLOGS2/ XLNG1A,XLNG2A,XLNY1A,XLNY2A,
     & XLNY12P,XLNY12Q,XLNG12P,XLNG12Q 
       COMMON /XLOGSA/ XLNAP,XLNAM,XLNA2,XLNAQ
       COMMON /ALOGS/  ALNA2,ALNAM,ALNAP
       COMMON /YLOGSA/ YLNAS,YLNAD

       COMMON /POLYC/  GC,YC,GC2,YC2,GC3,YC3
       COMMON /CLOGS/  XLNCP,XLNCM,XLNC2,XLNCQ,CLNC2
       COMMON /CLOGSC/ CLNCM,CLNCP
       COMMON /XLOGSC/ XLNCCP,XLNCCM,XLNCCQ,XLNGC
       COMMON /YLOGSC/ YLNCS,YLNCD
C
       A2 = A*A
       A3 = A*A2
c
       CT = CTH
       CT2 = CT*CT       
       CT3 = CT*CT2       
       CT4 = CT*CT3
       CM0 = 1D0-CT 
       CP0 = 1D0+CT
C
       XLNGCR = 2D0*XLNGC-ALR-2D0*XLN2
       XLNMR  = XLNM-1D0-ALR-2D0*XLN2
       XLNGCA = 2D0*XLNGC-2D0*XLN2+XLNA2+1D0
C
       GLG = CM0*CP0/GC*(CT+(YC+2D0*R*(1D0+R))/GC)
       GY1 = 1D0/2D0/GC*(CM0*CP0-(CT+R)* YC/GC)
       GY2 = 1D0/2D0/GC*(CT+R)
       FY1 = 1D0/6D0/GC*(2D0*GLG-(CM0*CM0*CP0
     & +4D0*(CT+R))-A2*(CT+R))
       FY2 = -1D0/3D0*GY1 

C-----------------------------------------------
C Initial-state term:

       RT = 
     & + 1D0/V*(CT*(1D0+CT2/3D0)*(XLNGCR+XLNMR)
     & + CT/3D0*CLNC2 - 4D0/3D0*(CLNCP-CLNCM))
     & + 1D0/6D0/V*((3D0+CT2+A2)*YLNCS - A*CT*YLNCD)
     & + FY1 * YLNCS + A * FY2 * YLNCD
     & + 2D0/3D0* (1D0+R) * XLNCCQ
     & + ( XLNGCR+XLNMR-XLNC2 )
     & *( 2D0/3D0*(4D0+R)
     &  - 4D0/3D0/GC3 
     &  * ( 7D0+12D0*CT+8D0*CT2+4D0*CT3+CT4+8D0*R)
     &  + 2D0/3D0/GC2 
     &  * (24D0+25D0*CT+11D0*CT2+CT3-CT4+12D0*R)
     &  - 1D0/3D0/GC 
     &  * (26D0+19D0*CT+3D0*CT2-CT3+CT4+12D0*R) )
     & + 4D0/3D0/GC3 
     &  * (15D0+28D0*CT+22D0*CT2+12D0*CT3+3D0*CT4+16D0*R)
     & - 2D0/3D0/GC2 
     &  * (52D0+63D0*CT+37D0*CT2+9D0*CT3-CT4+24D0*R)
     & + 1D0/3D0/GC 
     &  * (52D0+41D0*CT+10D0*CT2-3D0*CT3+2D0*R)       
     & - 1D0/12D0
     &  * (23D0-9D0*R+2D0*(11D0-R)*CT-(5D0+R)*CT2)
     & + A2/12D0
     &  * (3D0-R+12D0*CT+3D0*(1D0+R)*CT2-2D0*V*CT3)
     & - A2/3D0/GC
     &  * (CP0*CP0+2D0*R)
     & + 2D0/3D0*CT*(1D0-A)*(1D0+A)/V
     & + CT/2D0*(1D0+1D0/3D0*A2*CT2)*V

C-----------------------------------------------
C Interference term:

       RTI = 
     & + 1D0/V*( 2D0 * (1D0-CT2) * XLNCCQ - 4D0*CT )
     & - XLNCCP * ( 1D0 + 2D0*R - CT2 + 2D0*R*CT )
     & +(XLNCCM + ALR)* ( 2D0 + R2 - CT2*(R2+2D0*R+2D0)+2D0*R*CT )
     & - 1D0/2D0*(V - CT*RP)*(V + CT*RP)*XLNGCA
     & + CT*(3D0 + 2D0*R + R2)
     & - 1D0 - R2


       END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


 
      SUBROUTINE RTFUN1(CTH,RT,RTI)
*=======================================================================
* Introduced by M. JACK on 03/03/99 05:00pm,
* modified by M. JACK on 22/04/99 13:00pm 

       IMPLICIT DOUBLE PRECISION (a-h,o-z)

       COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR

       COMMON /VAR/ A,V,RP,V2,RP2,AP,AM,R1PI,R1PI2,R1PI3
       COMMON /XMASSES/ BETA,BETA2

       COMMON /XLOGSRM/ XLNRP,XLNR,XLNM,XLLPR2,XLN2

       COMMON /XLOGS2/ XLNG1A,XLNG2A,XLNY1A,XLNY2A,
     & XLNY12P,XLNY12Q,XLNG12P,XLNG12Q 
       COMMON /XLOGSA/ XLNAP,XLNAM,XLNA2,XLNAQ
       COMMON /ALOGS/  ALNA2,ALNAM,ALNAP
       COMMON /YLOGSA/ YLNAS,YLNAD

       COMMON /POLYC/  GC,YC,GC2,YC2,GC3,YC3
       COMMON /CLOGS/  XLNCP,XLNCM,XLNC2,XLNCQ,CLNC2
       COMMON /CLOGSC/ CLNCM,CLNCP
       COMMON /XLOGSC/ XLNCCP,XLNCCM,XLNCCQ,XLNGC
       COMMON /YLOGSC/ YLNCS,YLNCD
C 
       A2 = A*A

       CT = CTH       
       CT2 = CT*CT       
       CT3 = CT*CT2       

       CM0 = 1D0-CT 
       CP0 = 1D0+CT
C
       XLNGCR = 2D0*XLNGC-ALR-2D0*XLN2
       XLNMR = XLNM-1D0-ALR-2D0*XLN2
C
       GLG = CM0*CP0/GC*(CT+(YC+2D0*R*(1D0+R))/GC)
       GY1 = 1D0/2D0/GC*(CM0*CP0-(CT+R)* YC/GC)
       GY2 = 1D0/2D0/GC*(CT+R)
       FY1 = 1D0/6D0/GC*(2D0*GLG-(CM0*CM0*CP0
     & +4D0*(CT+R))-A2*(CT+R))
       FY2 = -1D0/3D0*GY1 

C-----------------------------------------------
C Initial-state term:

       RT = 
     & + (1D0+R2)/2D0/V*( A*(1D0+A2/3D0)*XLNMR 
     & + A/3D0*ALNA2 - 4D0/3D0* ( ALNAP - ALNAM) )
     & + 1D0/6D0/V *( (3D0+CT2+A2)*YLNCD - A*CT*YLNCS )       
     & + FY1 * YLNCD + A * FY2 * YLNCS      
     & - 2D0/3D0* A * GLG       
     & + A/2D0*(1D0+A2*CT2/3D0)*V 
     & + A/4D0/V*(YC2-A2/3D0*(GC2+8D0*R))
     & + 2D0*A/3D0*CM0*CP0/V

C-----------------------------------------------
C Interference term:

       RTI =
     & + 1D0/V*(2D0*(1D0-CT2)*XLNAQ - 4D0*A) 
     & + XLNAQ*(-1D0/2D0*(3D0 + 2D0*R + R2)*(1D0-CT2)-2D0*R*CT)
     & + A*(2D0*(1D0+R+R2)+1D0/2D0*(1D0-R2)*(1D0+CT2))
     & - A*CT*V2


       END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C------------------------------C
C Forward-backward asymmetry:  C
C------------------------------C
 
      SUBROUTINE AHACOL(H3,H1,H6)
*     ========== =============
* Introduced by M. JACK on 08/03/99 05:00pm
*
C Correction by M. Jack (07/05/99): 
C Logarithm XLGA had not been defined.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /MASSZ / AME,AMF,AME2,AMF2
      COMMON / SVAR/ S
      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /DELVAR/ DEL,DELM,DIF,DIF3,ALD,ALDM,ALDPL,ALDMI
     +               ,RD,RDI,RDI2,RDI3,RDM,RDMI,RDMI2,RDMI3
     +               ,ALRD,ALRDM
      COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
      COMMON /CORINT/ CORINT

      COMMON /INDEXA/ INDA,INDM      
C------------------------------------------------------------------------

      RR=R
C
      RP = 1D0+R
      RP2 = RP*RP
      V  = 1D0-R
      R1MI  = 1D0/V
      R1PI  = 1D0/RP
      R1PI2 = R1PI/RP
      DIFA0 = V/RP
      XMU2 = 4D0*AMF2/S
      BT2 = 1D0-4D0*AME2/S
      BT  = SQRT(BT2)
C      
      IF (INDA.EQ.1) THEN

*******************************************************
* acollinearity cut region I (''s'-cut case'') (A=1)
*******************************************************

         H6 = (1D0+R2)*TF + 2D0*ALR + V*V 

*******************************************************        
* Results without or with final-state mASses:
*******************************************************

         IF (INDM.EQ.0) THEN
 
            H3 = 4D0*R*(1D0+R2)/RP2
     &         * (TE-ALR-2D0*LOG(2D0)+2D0*LOG(RP))
            H1 = 2D0/3D0*R1MI*(2D0*R*(1D0-4D0*R+R2)/RP
     &         + R*(3D0+5D0*R2)*ALR
     &         - (1D0+R)*(5D0-2D0*R+5D0*R2)*LOG(RP))
            
         ELSE

            DELA  = 1D0/2D0*(1D0-SQRT(1D0-XMU2/R))
            DELMA = 1D0/2D0*(1D0+SQRT(1D0-XMU2/R))
            DIFA  = SQRT(1D0-XMU2/R)
            DIFA2 = 1D0-XMU2/R
C
            AP = 2D0*DELMA
            AM = 2D0*DELA
C                    
            YAPR = V+DIFA*RP
            YAMR = V-DIFA*RP
            GAPR = RP+DIFA*V
            GAMR = RP-DIFA*V
            YLAP = YAPR*LOG(ABS(YAPR))+YAMR*LOG(ABS(YAMR))
            YLAM = YAPR*LOG(ABS(YAPR))-YAMR*LOG(ABS(YAMR))
            XLAP = LOG(ABS(AP))
            XLAM = LOG(ABS(AM))

            XLOGA  = LOG(MAX(1.D-32,ABS(AP/AM)))
            XLGA   = LOG(MAX(1.D-32,ABS(GAPR/GAMR)))
            XLGAS  = LOG(MAX(1.D-32,ABS(GAPR*GAMR)))
            XYLGAM = LOG(MAX(1.D-32,ABS(YAPR/YAMR)))
            XYLGAP = LOG(MAX(1.D-32,ABS(YAPR*YAMR)))             
C
            IF (DIFA.LT.DIFA0) THEN 

C            H3 = (1D0+R2)*(R1PI*(DIFA*YLAP-DIFA0*YLAM)
C     &      + 4D0*DELA*DELMA*(XLAP-XLAM)) 
C     &      + 8D0*DIFA*R2/RP

C-------------------------------------------------------------
C alternative:
            H3 = (1D0+R2)*(-YAPR*YAMR*XYLGAM/RP2
     &      + 4D0*DELA*DELMA*XLOGA) + 8D0*DIFA*R2/RP
C-------------------------------------------------------------

            H1 = 2D0*R*(XLGA+DIFA*(1D0-R)
     &      - (2D0-R+5D0/3D0*R2)*R1MI*(XLAP-XLAM))

            ELSE

C            H3 = (1D0+R2)*(R1PI*(DIFA*YLAM-DIFA0*YLAP)
C     &      - 4D0*DELA*DELMA*(TE-ALR-2D0*LOG(2D0)-2D0*XLAP)) 
C     &      + 8D0*DIFA*DELA*R
C     &      + 4D0*R*(1D0+R2)/RP2
C     &      * (TE-2D0*ALR-4D0*LOG(2D0)+2D0*LOG(1D0+R))

C-------------------------------------------------------------
C alternative:
            H3 = (1D0+R2)*((4D0*R*(TE-ALR-2D0*LOG(2D0)
     &      + LOG(RP2/R/4D0)) - YAPR*YAMR*XYLGAP)/RP2
     &      - 4D0*DELA*DELMA*(TE-ALR-2D0*LOG(2D0)-2D0*XLAP)) 
     &      + 8D0*DIFA*DELA*R
C-------------------------------------------------------------

            H1 = 2D0*R*XLGA+2D0*R*(2D0-R+5D0/3D0*R2)*R1MI*ALR
     &      - 2D0/3D0*(1D0+R)*(5D0-2D0*R+5D0*R2)*R1MI
     &      * (LOG(1D0+R)-LOG(2D0)+XLAP)
     &      - 1D0/3D0*(1D0-R)/RP*(1D0-4D0*R+R2)
     &      + DIFA*(1D0-R)*(1D0+R)-2D0/3D0*DIFA2*(1D0+R3)*R1MI  

            ENDIF

         ENDIF  

*****************************************************       
* Regularize initial-state term by initial-state 
* soft+vertex contribution: 
*****************************************************
         IF (R.GE.RCUT) THEN
            H3=(H3-2D0*TE)*R1MI
            IF (IPHOT2.GE.0) H3=H3+AH2(RR)
         ELSE
            H3=H3*R1MI
         END IF

*****************************************************       
* Regularize final-state term by divergent 
* final-state soft+vertex contribution: 
*****************************************************
C
C         IF (R.GE.RCUT) THEN
C            H6 = H6- ...
C            IF (IPHOT2.GE.0) H6=H6 + ....
C         ELSE
            H6=H6*R1MI
C         END IF
C  
      ELSE

***************************************
* acollinearity cut region II or III 
***************************************          

         DELA  = DEL
         DELMA = DELM
         DIFA  = DIF

         DIFA2 = DIFA*DIFA
         AP = 2D0*DELMA
         AM = 2D0*DELA
C                
         YAPR = V+RP*DIFA
         YAMR = V-RP*DIFA
         GAPR = RP+V*DIFA
         GAMR = RP-V*DIFA
c
         YAPR0 = YAPR/RP
         YAMR0 = YAMR/RP
         GAPR0 = GAPR/RP
         GAMR0 = GAMR/RP
c     
         YLAP = YAPR*LOG(MAX(1.D-32,ABS(YAPR)))
     &        + YAMR*LOG(MAX(1.D-32,ABS(YAMR)))
         YLAM = YAPR*LOG(MAX(1.D-32,ABS(YAPR)))
     &        - YAMR*LOG(MAX(1.D-32,ABS(YAMR)))
c
         XLAP = LOG(MAX(1.D-32,ABS(AP)))
         XLAM = LOG(MAX(1.D-32,ABS(AM)))
c
         XLOGA  = LOG(MAX(1.D-32,ABS(AP/AM)))
         XLOGAS = LOG(MAX(1.D-32,ABS(AP*AM)))
         XLGA   = LOG(MAX(1.D-32,ABS(GAPR/GAMR)))
         XLGAS  = LOG(MAX(1.D-32,ABS(GAPR*GAMR)))       
         XYLGAM = LOG(MAX(1.D-32,ABS(YAPR/YAMR)))
         XYLGAP = LOG(MAX(1.D-32,ABS(YAPR*YAMR))) 

C-------------------------------------------------------------
         H6 = ((1D0+R2)*XLOGA - 2D0*DIFA*XMU2/AP/AM)*R1MI 
     &  - RP*XLGA

         IF (DIFA.LT.DIFA0) THEN 

C          H3 = (1D0+R2)*(R1PI*(R1MI*DIFA*YLAP-R1PI*YLAM)
C     &  + 4D0*DELA*DELMA*R1MI*(XLAP-XLAM)) 
C     &  + 8D0*DIFA*R2/RP*R1MI
C-------------------------------------------------------------
C alternative:
C          H3 = (1D0+R2)*(-YAPR*YAMR*XYLGAM/RP2
C     &  + 4D0*DELA*DELMA*XLOGA) + 8D0*DIFA*R2/RP
C          H3 = H3*R1MI
C-------------------------------------------------------------
C-------------------------------------------------------------
C alternative:
          H3 = (1D0+R2)*(-YAPR0*YAMR0*XYLGAM
     &  + 4D0*DELA*DELMA*XLOGA) + 8D0*DIFA*R2/RP
          H3 = H3*R1MI
C-------------------------------------------------------------

          H1 = 2D0*R*(XLGA+DIFA*V
     &       - (2D0-R+5D0/3D0*R2)*XLOGA/V)

         ELSE

C          H3 = (1D0+R2)*(R1PI*(R1MI*DIFA*YLAM-R1PI*YLAP)
C     &  - 4D0*R1MI*DELA*DELMA*(TE-ALR-2D0*LOG(2D0)-2D0*XLAP)) 
C     &  + 8D0*DIFA*DELA*R*R1MI
C     &  + 4D0*R*(1D0+R2)/RP2*R1MI
C     &  * (TE-2D0*ALR-4D0*LOG(2D0)+2D0*LOG(1D0+R))
C-------------------------------------------------------------
C alternative:
C          H3 = (1D0+R2)*((4D0*R*(TE-ALR-2D0*LOG(2D0)
C     &  + LOG(RP2/R/4D0)) - YAPR*YAMR*XYLGAP)/RP2
C     &  - 4D0*DELA*DELMA*(TE-ALR-2D0*LOG(2D0)-2D0*XLAP)) 
C     &  + 8D0*DIFA*DELA*R
C          H3 = H3*R1MI
C-------------------------------------------------------------
C-------------------------------------------------------------
C alternative:
          H3 = (1D0+R2)*((4D0*R*(TE-ALR-2D0*LOG(2D0)
     &  + LOG(RP2/R/4D0))/RP2 - YAPR0*YAMR0*XYLGAP)
     &  - 4D0*DELA*DELMA*(TE-ALR-2D0*LOG(2D0)-2D0*XLAP)) 
     &  + 8D0*DIFA*DELA*R
          H3 = H3*R1MI
C-------------------------------------------------------------

C-------------------------------------------------------------
C          H1 = 2D0*R*XLGA+2D0*R*(2D0-R+5D0/3D0*R2)*R1MI*ALR
C     &  - 2D0/3D0*(1D0+R)*(5D0-2D0*R+5D0*R2)*R1MI
C     &  * (LOG(1D0+R)-LOG(2D0)+XLAP)
C     &  - 1D0/3D0*(1D0-R)/RP*(1D0-4D0*R+R2)
C     &  + DIFA*(1D0-R)*(1D0+R)-2D0/3D0*DIFA2*(1D0+R3)*R1MI  
C-------------------------------------------------------------

          H1 = 2D0*R*XLGA*V
     &  + 2D0*R*(2D0-R+5D0/3D0*R2)*ALR
     &  - 2D0/3D0*(5D0-2D0*R+5D0*R2)*RP
     &  * (LOG(RP)-LOG(2D0)+XLAP)
     &  - 1D0/3D0*(1D0-4D0*R+R2)/RP*V*V
     &  + DIFA*RP*V*V
     &  - 2D0/3D0*DIFA2*(1D0+R3)  

          H1 = H1*R1MI

         ENDIF

      ENDIF
C
C       WRITE(*,*)
C       WRITE(*,*) 'INDA,INDM =',INDA,INDM
C       WRITE(*,*) 'DIF,DEL,DELM =',DIFA,DELA,DELMA
C       WRITE(*,*) 'XMU2,S,BT =',XMU2,S,BT
C       WRITE(*,*) 'RP2,RP,V =',RP2,RP,V
C       WRITE(*,*) 'TF,ALR =',TF,ALR
C       WRITE(*,*) 'YAPR,YAMR,GAPR,GAMR= ',YAPR,YAMR,GAPR,GAMR
C       WRITE(*,*) 'YLAP,YLAM,XLAP,XLAM= ',YLAP,YLAM,XLAP,XLAM
C       WRITE(*,*) 'XYLGAM,XYLGAP= ',XYLGAM,XYLGAP
C       WRITE(*,*) 'XLOGA,XLGA,XLGAS= ',XLOGA,XLGA,XLGAS
C       WRITE(*,*)
C
      END




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


 
      SUBROUTINE AHCUTACOL(H3,H1,H6)
*     ========== ============
* Introduced by M. JACK on 08/03/99 05:00pm,
* modified by M. JACK on 22/04/99 13:00pm 
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
       COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
       COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM
       COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT
       COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM

       COMMON /DELVAR/ DEL,DELM,DIF,DIF3,ALD,ALDM,ALDPL,ALDMI
     +               ,RD,RDI,RDI2,RDI3,RDM,RDMI,RDMI2,RDMI3
     +               ,ALRD,ALRDM

       COMMON /INDEXA/ INDA,INDM
       COMMON /INDEX/ MF,IIF,JF,KF

C------------------------------------------------
C With acceptance cut:    
      CALL RADFBIN(H3M,H1M)
      CALL RADFBFIN(H6M)

c      WRITE(*,*) '---------------------------'
c      WRITE(*,*) 'I,J,K =',IIF,JF,KF
c      WRITE(*,*) 'R,C,A =',R,C,DIF
c      WRITE(*,*) 'H6,H3,H1 =',H6M,H3M,H1M 

C------------------------------------------------
C No acceptance cut (C=1) (just for checks):
C     CALL RFBACOL(H3M)
C------------------------------------------------

      H3=H3M
      H1=H1M
      H6=H6M

C------------------------------------------------
C regularize in s' cut region:
      IF (R.GE.RCUT) THEN
        H3=H3-TE*2D0*C2*R1MI
        IF (IPHOT2.GE.0) H3=H3+AH2(R)*C2
      ENDIF
C------------------------------------------------

C------------------------------------------------
C Re-create old results by M. Bilenky:
C      IF (IIF.EQ.JF) RETURN
C      I = IIF
C      J = JF
C      K = KF
C
C------------------------------------------------
C With acceptance Cut:
C      CALL ACONST (H3M,H1M,I,J,K)   
C
C No acceptance cut (C=1):
C      CALL ACACOL (H3M,H1M,K)
C------------------------------------------------
C 
C      H3=H3-H3M
C      H1=H1-H1M
C
C*------------------------------------------------
C* Insertion in order to subtract non-logarithmiC 
C* terms prop. to 1, C, C^2 for Comparisons with 
C* Bilenky's original Code:
C*------------------------------------------------
C
C      CALL SETRC ( 1)
C      CALL AHCUTM (H3MB,H1MB,I)
C      H3B=H3MB
C      H1B=H1MB
C      CALL SETRC (-1)
C      CALL AHCUTM (H3MB,H1MB,J)
C      H3B=H3B+H3MB
C      H1B=H1B+H1MB
C      CALL SETRC ( 0)
C      CALL AHCUTM (H3MB,H1MB,K)
C      H3B=H3B-2*H3MB
C      H1B=H1B-2*H1MB
C
C      IF (I.NE.K) THEN
C       CALL AJUMP (H3MB,H1MB,I,K)
C       H3B=H3B+H3MB
C       H1B=H1B+H1MB
C      END IF
C
C      IF (J.NE.K) THEN
C       CALL AJUMP (H3MB,H1MB,J,K)
C       H3B=H3B+H3MB
C       H1B=H1B+H1MB
C      END IF 
C
C       IB = I
C       JB = J
C       KB = K
C       I = IIF
C       J = JF
C       K = KF
C 
C       IF (I.EQ.J) RETURN
C       IF (DELC.GT.epsccut) THEN 
C         CALL ACONST (H3M,H1M,I,J,K)
C       ELSE
C         CALL ACACOL (H3M,H1M,K)
C       ENDIF 
C       H3=H3-H3M
C       H1=H1-H1M
C
C       DELH3 = ABS((H3-H3B)/H3B)
C       IF (DELH3.GE.EPSDEL) THEN
C          WRITE(*,*)
C          WRITE(*,*) '----------- R_FB -------------'  
C          WRITE(*,*) 'Comparison of Codes: Bilenky -- JACK'
C          WRITE(*,*) 'IB,JB,KB = ',IB,JB,KB
C          WRITE(*,*) 'I,J,K    = ',IIF,JF,KF 
C          WRITE(*,*) 'M.B.: R , RFBINI = ',R,H3B
C          WRITE(*,*) 'M.J.: R , RFBINI = ',R,H3
C          WRITE(*,*) 'M.J.-M.B.: ',DELH3
C          WRITE(*,*) 'M.J.:  RFBNONLOG = ',H3M
C          WRITE(*,*) '::::::::::::::::::::::::::::::::::::'
C          WRITE(*,*)
C       ENDIF

      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C-----------------------------------------------
C Initial-state and interference contributions :
C-----------------------------------------------
 
      SUBROUTINE RADFBIN(RADFB,RADFBI)
*=======================================================================
* Introduced by M. JACK on 03/03/99 05:00pm,
* modified by M. JACK on 22/04/99 13:00pm 
*
C Correction by M. Jack (07/05/99): missing exponent (**2) 
C in term R2/6D0 of expression G0I inserted  

       IMPLICIT DOUBLE PRECISION (a-h,o-z)

       COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
       COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM

       COMMON /VAR/ A,V,RP,V2,RP2,AP,AM,R1PI,R1PI2,R1PI3
       COMMON /XMASSES/ BETA,BETA2
       COMMON /XLOGSRM/ XLNRP,XLNR,XLNM,XLLPR2,XLN2

       COMMON /POLYNOMIALS/ G1C,G2C,Y1C,Y2C,
     & G1A,G2A,Y1A,Y2A,YG1CP,YG1CM,YG2CP,YG2CM

       COMMON /XLOGS2/ XLNG1A,XLNG2A,XLNY1A,XLNY2A,
     & XLNY12P,XLNY12Q,XLNG12P,XLNG12Q 
       COMMON /XLOGSA/ XLNAP,XLNAM,XLNA2,XLNAQ
       COMMON /ALOGS/  ALNA2,ALNAM,ALNAP
       COMMON /YLOGSA/ YLNAS,YLNAD
C
       COMMON/INDEX/ MF,IIF,JF,KF
C----------------------------------------------------------

       I = IIF
       J = JF
       K = KF

       CT = BETA*C
       CT2 = CT*CT

       A2 = A*A

C----------------------------------------------------------
C Initial-state terms:

       XLNMR = XLNM-1D0-2D0*XLN2-ALR
       XLNRPR = 2D0*XLNRP-2D0*XLN2-ALR

C----------------------------------------------------------
C       G0 = (1D0+R2)/RP*(1D0/V*A*YLNAD-1D0/RP*YLNAS)
C     &    + 4D0*R/RP*(1D0+R2)/RP/V*( XLNRPR+XLNMR )
C
C       G1 = (1D0+R2)/RP*(1D0/V*A*YLNAS-1D0/RP*YLNAD)
C     &    - 4D0*R/RP * A
C
C       G10 = 1D0/2D0*(1D0+R2)/V*((1D0-A2)*XLNMR-ALNA2) 
C     &     + 2D0*A2*R/V
C----------------------------------------------------------

       G0 = (1D0+R2)/RP2*R1MI
     &    * (-Y1A*Y2A*XLNY12P+4D0*R*(XLNRPR+XLNMR))     

       G1 = (-(1D0+R2)/RP*R1MI*Y1A*Y2A*XLNY12Q-4D0*R*A)/RP

       G10 = (1D0+R2)*R1MI/2D0*((1D0-A2)*XLNMR-ALNA2)
     &     + 2D0*A2*R*R1MI

C----------------------------------------------------------
C Interference terms:
 
       G0I = - 2D0/3D0*R1MI 
     & * ( 4D0*( XLNA2 + XLNRPR ) + A2 )
     & + ( 11D0 + 8D0*R + 5D0*R2 )/6D0 
     & * ( XLNA2 + XLNRPR - ALR ) 
     & + ( 1D0 + 2D0*R ) * ALR
     & + A2/3D0* (1D0 + R + R2)
     & + 11D0/6D0 - R + R2/6D0 - 2D0/RP
     & - A/2D0*CT2*V*RP

       G1I = ( 5D0 - 4D0*R + 5D0*R2 )/3D0*XLNAQ
     & - A*V*(V + CT2*RP)

       G10I = CT2*(1D0/2D0*V*RP*XLNAQ + R*XLNG12Q) 

C----------------------------------------------------------

       RADFB = 0D0
       RADFBI = 0D0

       IF ((I.EQ.0).AND.(J.EQ.0)) THEN

          CALL SLOGSP 
          CALL RFBFUN0(CT,G002,G002I)
          CALL SLOGSM 
          CALL RFBFUN0(-CT,G001,G001I)

          RFB00 =  G002 + G001 + G0       
          RADFB =  RFB00

          RFB00I =  G002I + G001I        
          RADFBI =  RFB00I

       ELSE
       IF ((I.EQ.1).AND.(J.EQ.1)) THEN

          CALL SLOGSP 
          CALL RFBFUN1(CT,G112,G112I)
          CALL SLOGSM 
          CALL RFBFUN1(-CT,G111,G111I)

          RFB11 =  G112 + G111 + G1 
          RADFB =  RFB11

          RFB11I =  G112I + G111I + G1I + 2D0*G10I
          RADFBI =  RFB11I

       ELSE
       IF ((I.EQ.1).AND.(J.EQ.0).AND.(K.EQ.1)) THEN

          CALL SLOGSP 
          CALL RFBFUN1(C,G112,G112I)
          CALL SLOGSM 
          CALL RFBFUN0(-C,G001,G001I)

          RFB10K1= G112 + G001 + G1 + G10
          RADFB  = RFB10K1

          RFB10K1I= G112I + G001I + G1I - G0I + G10I
          RADFBI  = RFB10K1I

       ELSE
       IF ((I.EQ.1).AND.(J.EQ.0).AND.(K.EQ.0)) THEN

          CALL SLOGSP 
          CALL RFBFUN1(CT,G112,G112I)
          CALL SLOGSM 
          CALL RFBFUN0(-CT,G001,G001I)

          RFB10K0= G112 + G001 + G0 - G10
          RADFB  = RFB10K0

          RFB10K0I= G112I + G001I + G0I + G10I 
          RADFBI  = RFB10K0I

       ELSE
       IF ((I.EQ.1).AND.(J.EQ.-1).AND.(K.EQ.1)) THEN 

          CALL SLOGSP 
          CALL RFBFUN1(CT,G112,G112I)
          CALL SLOGSM 
          CALL RFBFUN1(-CT,G111,G111I)

          RFB11MK1= G112 - G111 + G1 
          RADFB   = RFB11MK1

          RFB11MK1I= G112I - G111I + G1I + 2D0*G10I 
          RADFBI   = RFB11MK1I

       ELSE
       IF ((I.EQ.1).AND.(J.EQ.-1).AND.(K.EQ.0)) THEN

          CALL SLOGSP 
          CALL RFBFUN1(CT,G112,G112I)
          CALL SLOGSM 
          CALL RFBFUN1(-CT,G111,G111I)

          RFB11MK0= G112 - G111 + G0 - 2*G10
          RADFB   = RFB11MK0

          RFB11MK0I= G112I - G111I + 2D0*(G0I+G10I) 
          RADFBI   = RFB11MK0I

       ELSE
       IF ((I.EQ.0).AND.(J.EQ.1).AND.(K.EQ.1)) THEN

          CALL SLOGSP 
          CALL RFBFUN0(CT,G002,G002I)
          CALL SLOGSM 
          CALL RFBFUN1(-CT,G111,G111I)

          RFB01K1= G002 + G111 + G1 + G10
          RADFB  = RFB01K1

          RFB01K1I= G002I + G111I + G1I - G0I + G10I
          RADFBI  = RFB01K1I

       ELSE
       IF ((I.EQ.0).AND.(J.EQ.1).AND.(K.EQ.0)) THEN

          CALL SLOGSP 
          CALL RFBFUN0(CT,G002,G002I)
          CALL SLOGSM 
          CALL RFBFUN1(-CT,G111,G111I)

          RFB01K0= G002 + G111 + G0 - G10
          RADFB  = RFB01K0

          RFB01K0I= G002I + G111I + G0I + G10I
          RADFBI  = RFB01K0I

       ELSE
       IF ((I.EQ.-1).AND.(J.EQ.1).AND.(K.EQ.1)) THEN

          CALL SLOGSP 
          CALL RFBFUN1(CT,G112,G112I)
          CALL SLOGSM 
          CALL RFBFUN1(-CT,G111,G111I)

          RFB1M1K1= G111 - G112 + G1 
          RADFB   = RFB1M1K1

          RFB1M1K1I= G111I - G112I + G1I + 2D0*G10I
          RADFBI   = RFB1M1K1I

       ELSE
       IF ((I.EQ.-1).AND.(J.EQ.1).AND.(K.EQ.0)) THEN

          CALL SLOGSP 
          CALL RFBFUN1(CT,G112,G112I)
          CALL SLOGSM
          CALL RFBFUN1(-CT,G111,G111I)

          RFB1M1K0= G111 - G112 + G0 - 2*G10
          RADFB   = RFB1M1K0

          RFB1M1K0I= G111I - G112I + 2D0*(G0I+G10I)
          RADFBI   = RFB1M1K0I

       ENDIF
       ENDIF
       ENDIF
       ENDIF
       ENDIF
       ENDIF
       ENDIF
       ENDIF
       ENDIF
       ENDIF
C
C       RADFB = RADFB/BETA
C       RADFBI = RADFBI/BETA
C
       END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


 
      SUBROUTINE RFBFUN0(CTH,RFB,RFBI)
*=======================================================================
* Introduced by M. JACK on 03/03/99 05:00pm,
* modified by M. JACK on 22/04/99 13:00pm 

       IMPLICIT DOUBLE PRECISION (a-h,o-z)

       COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      
       COMMON /VAR/ A,V,RP,V2,RP2,AP,AM,R1PI,R1PI2,R1PI3
       COMMON /XMASSES/ BETA,BETA2
       COMMON /XLOGSRM/ XLNRP,XLNR,XLNM,XLLPR2,XLN2

       COMMON /POLYNOMIALS/ G1C,G2C,Y1C,Y2C,
     & G1A,G2A,Y1A,Y2A,YG1CP,YG1CM,YG2CP,YG2CM

       COMMON /XLOGS2/ XLNG1A,XLNG2A,XLNY1A,XLNY2A,
     & XLNY12P,XLNY12Q,XLNG12P,XLNG12Q 
       COMMON /XLOGSA/ XLNAP,XLNAM,XLNA2,XLNAQ
       COMMON /ALOGS/  ALNA2,ALNAM,ALNAP
       COMMON /YLOGSA/ YLNAS,YLNAD

       COMMON /POLYC/  GC,YC,GC2,YC2,GC3,YC3
       COMMON /CLOGS/  XLNCP,XLNCM,XLNC2,XLNCQ,CLNC2
       COMMON /CLOGSC/ CLNCM,CLNCP
       COMMON /XLOGSC/ XLNCCP,XLNCCM,XLNCCQ,XLNGC
       COMMON /YLOGSC/ YLNCS,YLNCD
C
       A2 = A*A
       A3 = A*A2

       CT = CTH 
       CT2 = CT*CT       
       CT3 = CT*CT2       

       CM0 = 1D0-CT 
       CP0 = 1D0+CT
       CM02 = CM0*CM0
       CP02 = CP0*CP0
C
       XLNGCR  = 2D0*XLNGC-2D0*XLN2-ALR
       XLNRPR  = 2D0*XLNRP-2D0*XLN2-ALR
       XLNMR   = XLNM-1D0-ALR-2D0*XLN2

       XLNGCRM = 2D0*XLNGC-2D0*XLN2-2D0*XLNCCM+XLNA2
       XLNGCM  = XLNGC - XLNRP - XLNCCM
C
       GLG = CM0*CP0/GC*(CT+(YC+2D0*R*(1D0+R))/GC)
       GY1 = 1D0/2D0/GC*(CM0*CP0-(CT+R)* YC/GC)
       GY2 = 1D0/2D0/GC*(CT+R)

C----------------------------------------------------------
C Initial-state terms:

       RFB =
     & + 1D0/2D0/V * ( CT * YLNCS - A * YLNCD
     & - 2D0*CM0*CP0 * ( XLNGCR+XLNMR-XLNC2 ) )
     & + GY1 * YLNCS + A * GY2 * YLNCD
     & + GLG * ( XLNGCR+XLNMR-XLNC2 )
     & + 1D0/2D0*CT*(1+R)*( A2 + YC2/GC2 )

C----------------------------------------------------------
C Interference terms:

       RFBI = 
     & + 1D0/V*( - 2D0*(CT+CT3/3D0)*(XLNCCP-XLNCCM-ALR)
     & - 8D0/3D0*XLNC2 + 16D0/3D0*(XLNGC-XLNRP)
     & - 2D0/3D0*CT2 )
     & + ( 1D0 + 2D0*R - CT2 + 2D0*R*CT )*(XLNCCP-XLNCCM)
     & + (1D0/3D0*CT3*(1D0+R+R2)+1D0/2D0*CT2*(1D0-R2)+CT*(1D0-R+R2))
     & * XLNGCRM
     & - ( 11D0/3D0 + 8D0/3D0*R + 5D0/3D0*R2 )
     & * XLNGCM
     & - 4D0/GC*(2D0*RP+4D0*CT+CT/RP+5D0*CT2/RP+2D0*CT3/RP2)
     & + CT2/2D0*(1D0-R2)*XLNAQ
     & + R*CT2*XLNG12Q
     & + ( - 2D0/3D0*CT3*(1D0+R+R2) + CT2*R2 - 2D0*CT*(1D0+R2) )*ALR
     & + CT3/12D0* ( A2*V2 + 5D0*RP2 )
     & + CT2*( - 1D0/2D0*A*(1D0-R2) + 4D0/3D0*R - 1D0/2D0*R2 
     & + 19D0/6D0 + 2/RP + 4D0/RP2 )
     & + CT*( A2/4D0*V2 - 1D0/6D0*R + 17D0/4D0 
     & + 1D0/4D0*R2 + 8D0/RP ) + 8D0 


       END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


 
      SUBROUTINE RFBFUN1(CTH,RFB,RFBI)
*=======================================================================
* Introduced by M. JACK on 03/03/99 05:00pm,
* modified by M. JACK on 22/04/99 13:00pm 

       IMPLICIT DOUBLE PRECISION (a-h,o-z)

       COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR

       COMMON /VAR/ A,V,RP,V2,RP2,AP,AM,R1PI,R1PI2,R1PI3
       COMMON /XMASSES/ BETA,BETA2
       COMMON /XLOGSRM/ XLNRP,XLNR,XLNM,XLLPR2,XLN2

       COMMON /POLYNOMIALS/ G1C,G2C,Y1C,Y2C,
     & G1A,G2A,Y1A,Y2A,YG1CP,YG1CM,YG2CP,YG2CM

       COMMON /XLOGS2/ XLNG1A,XLNG2A,XLNY1A,XLNY2A,
     & XLNY12P,XLNY12Q,XLNG12P,XLNG12Q 
       COMMON /XLOGSA/ XLNAP,XLNAM,XLNA2,XLNAQ
       COMMON /ALOGS/  ALNA2,ALNAM,ALNAP
       COMMON /YLOGSA/ YLNAS,YLNAD

       COMMON /POLYC/  GC,YC,GC2,YC2,GC3,YC3
       COMMON /CLOGS/  XLNCP,XLNCM,XLNC2,XLNCQ,CLNC2
       COMMON /CLOGSC/ CLNCM,CLNCP
       COMMON /XLOGSC/ XLNCCP,XLNCCM,XLNCCQ,XLNGC
       COMMON /YLOGSC/ YLNCS,YLNCD
C
       A2 = A*A
       A3 = A*A2

       CT = CTH 
       CT2 = CT*CT       
       CT3 = CT*CT2       

       CM0 = 1D0-CT 
       CP0 = 1D0+CT
       CM02 = CM0*CM0
       CP02 = CP0*CP0
C
       XLNGCR  = 2D0*XLNGC-2D0*XLN2-ALR
       XLNRPR  = 2D0*XLNRP-2D0*XLN2-ALR
       XLNMR   = XLNM-1D0-ALR-2D0*XLN2
C
       GLG = CM0*CP0/GC*(CT+(YC+2D0*R*(1D0+R))/GC)
       GY1 = 1D0/2D0/GC*(CM0*CP0-(CT+R)* YC/GC)
       GY2 = 1D0/2D0/GC*(CT+R)

C----------------------------------------------------------
C Initial-state terms:

       RFB = 
     & + 1D0/2D0/V * (CT * YLNCD - A * YLNCS)
     & + GY1 * YLNCD + A * GY2 * YLNCS
     & + A * (2D0*R/V+(1D0+R)*CT) * YC/GC

C----------------------------------------------------------
C Interference terms:

       RFBI =
     & - 2D0/V*CT*(1D0 + CT2/3D0)*XLNAQ
     & + XLNAQ*( CT*(1D0 + CT2/3D0)*(1D0 + R + R2) 
     & - CT2/2D0*(1D0+R2) - 5D0/6D0 + 2D0/3D0*R - 5D0/6D0*R2 )
     & + A/2D0*(1D0 + CT2)*V*( V + CT*RP )


       END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C----------------------------------------
C Final-state contributions :
C----------------------------------------

C-----------------------C
C Total cross-section:  C
C-----------------------C
 
      SUBROUTINE RADTFIN(RTFIN)
*=======================================================================
* Introduced by M. JACK on 22/04/99 13:00pm

C-----------------------------------------------------
C Hard interference contributions for acollinearity
C + acceptance cut (not regularized)
C-----------------------------------------------------

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       COMMON /SVAR  / S
       COMMON /MASSZ / AME,AMF,AME2,AMF2
       COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
       COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     &                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM
       COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR

       COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT

       COMMON /VAR/ A,V,RP,V2,RP2,AP,AM,R1PI,R1PI2,R1PI3
       COMMON /XMASSES/ BETA,BETA2

       COMMON /POLYNOMIALS/ G1C,G2C,Y1C,Y2C,
     & G1A,G2A,Y1A,Y2A,YG1CP,YG1CM,YG2CP,YG2CM

       COMMON /XLOGSRM/ XLNRP,XLNR,XLNM,XLLPR2,XLN2
       COMMON /XLOGS2/ XLNG1A,XLNG2A,XLNY1A,XLNY2A,
     & XLNY12P,XLNY12Q,XLNG12P,XLNG12Q 
       COMMON /XLOGSA/ XLNAP,XLNAM,XLNA2,XLNAQ

       COMMON /INDEXA/ INDA,INDM
C------------------------------------------------------------

       PI=4D0*ATAN(1D0)  
C
       RTFIN = 0D0
C
       CT = BETA*C
       CT2 = CT*CT
       CT3 = CT*CT2
       COPL4 = CT*(1D0+CT2/3D0)
       COPL5 = CT*(1D0-CT)*(1D0+CT)  
C
       XMU2 = 4D0*AMF2/S
       XME2 = 4D0*AME2/S

C-----------------------------------------------------------------
C Hard final-state contribution:

       IF (INDA.EQ.1) THEN

        RTFIN= 
     &  + 3D0/4D0*(COPL4*(1D0+R2)*(TF+ALR)/V+COPL5*V)
       ELSE

        RTFIN= 
     &  + 3D0/4D0*( COPL4*(((1D0+R2)*XLNAQ
     &  - 2D0*A*XMU2/AP/AM)/V-A*V)
     &  + 4D0*COPL5*A*V*R/G1A/G2A ) 
       ENDIF

      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C------------------------------C
C Forward-backward asymmetry:  C
C------------------------------C

      SUBROUTINE RADFBFIN(RFBFIN)
*=======================================================================
* Introduced by M. JACK on 22/04/99 13:00pm

C-----------------------------------------------------
C Hard interference contributions for acollinearity
C + acceptance cut (not regularized)
C-----------------------------------------------------

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       COMMON /SVAR  / S
       COMMON /MASSZ / AME,AMF,AME2,AMF2
       COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
       COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     &                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM
       COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR

       COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT

       COMMON /VAR/ A,V,RP,V2,RP2,AP,AM,R1PI,R1PI2,R1PI3
       
       COMMON /XMASSES/ BETA,BETA2

       COMMON /POLYNOMIALS/ G1C,G2C,Y1C,Y2C,
     & G1A,G2A,Y1A,Y2A,YG1CP,YG1CM,YG2CP,YG2CM

       COMMON /XLOGSRM/ XLNRP,XLNR,XLNM,XLLPR2,XLN2
       COMMON /XLOGS2/ XLNG1A,XLNG2A,XLNY1A,XLNY2A,
     & XLNY12P,XLNY12Q,XLNG12P,XLNG12Q 
       COMMON /XLOGSA/ XLNAP,XLNAM,XLNA2,XLNAQ

       COMMON /INDEXA/ INDA,INDM
C------------------------------------------------------------

       PI=4D0*ATAN(1D0)  
C
       RFBFIN = 0D0
C
       CT = BETA*C
       CT2 = CT*CT
       CT3 = CT*CT2
C
       XMU2 = 4D0*AMF2/S
       XME2 = 4D0*AME2/S

C-----------------------------------------------------------------
C Hard final-state contribution:

       IF (INDA.EQ.1) THEN

        RFBFIN =
     +  + CT2*(((1D0+R2)*TF+2D0*ALR)*R1MI+V)
       ELSE

        RFBFIN =
     +  + CT2*(((1D0+R2)*XLNAQ-2D0*A*XMU2/AP/AM)*R1MI-RP*XLNG12Q)
       ENDIF

      END


C===========================================================================

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Extra Subroutines just used for checks (compatibility with former code):C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
      SUBROUTINE SCONST(H0,H4,I,J)
*     ========== =================
* Introduced by M. JACK on 03/03/99 05:00pm
* (not used, just for checks !)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM

      COMMON /VAR/ A,V,RP,V2,RP2,AP,AM,R1PI,R1PI2,R1PI3
* 
      A2 = A*A
      A3 = A*A2
*
      C2 = C*C
*
      RP  = 1D0+R
      RP2 = RP*RP
      RP3 = RP*RP2
*
      V  = 1D0-R
      V2 = V*V
      V3 = V*V2
*
      H0 = 0D0
      H4 = 0D0
*
      IF (((I.EQ.0).AND.(J.EQ.1)).OR.((I.EQ.1).AND.(J.EQ.0))) THEN
*
       H0 = - 1D0/12D0*C2*RP3  + 1D0/12D0*RP*V2
     +      + A/4D0*C2*RP2*V  - A/12D0*V*(3D0+2D0*R+3D0*R2)
     +      - A2/4D0*C2*RP*V2 + A2/12D0*RP*(3D0-2D0*R+3D0*R2)
     +      + A3/12D0*C2*V3   - A3/12D0*RP2*V    
       H4 = 0D0
       H0 = H0*R1MI2
       H4 = H4*R1MI3
      ELSE
       IF (((I.EQ.-1).AND.(J.EQ.1)).OR.((I.EQ.1).AND.(J.EQ.-1))) THEN
*
       H0 = A*( 1D0/2D0*C2*RP2 - 1D0/6D0*(3D0+2D0*R+3D0*R2)
     +        + A2/6D0*C2*V2  - A2/6D0*RP2 )
       H4 = 0D0
       H0 = H0*R1MI
       H4 = H4*R1MI3
       ELSE
         H0 = 0D0
         H4 = 0D0
       ENDIF
      ENDIF
*
      H0 = H0*ISIGN(1,I-J)
      H4 = H4*ISIGN(1,I-J)
*     
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


 
      SUBROUTINE RTACOL(RADT)
*=======================================================================
* Introduced by M. JACK on 03/03/99 05:00pm
* (not used, just for checks !)

       IMPLICIT DOUBLE PRECISION (a-h,o-z)
C
       COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
       COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM

       COMMON /VAR/ A,V,RP,V2,RP2,AP,AM,R1PI,R1PI2,R1PI3

       COMMON /XLOGSRM/ XLNRP,XLNR,XLNM,XLLPR2,XLN2
C
       A2 = A*A

       R2 = R*R
       V  = 1D0-R
       RP = 1D0+R
C
       RADT = A*((1D0+A2/3D0)*(1D0+R2)/V*(XLNM-1D0)
     &      + 2D0*(1D0-A2)*R/V)
C
       IF (C.LT.0D0) THEN
          RADT = -RADT
       ENDIF
C

       END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


 
      SUBROUTINE SCACOL(H0,H4)
*=======================================================================
* Introduced by M. JACK on 03/03/99 05:00pm
* (not used, just for checks !)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM

      COMMON /VAR/ A,V,RP,V2,RP2,AP,AM,R1PI,R1PI2,R1PI3
*
      A2 = A*A

      RP2 = RP*RP
      V  = 1D0-R
      RP = 1D0+R 
*
      H4 = 0D0
      H0 = 2D0/3D0*A*(1D0-A2)*R/V 
*
      IF (C.LT.0D0) THEN
         H0 = -H0
      ENDIF

      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


 
      SUBROUTINE ACONST(H3,H1,I,J,K)
*=======================================================================
* Introduced by M. JACK on 03/03/99 05:00pm
* ( not used, just for checks !)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM

      COMMON /VAR/ A,V,RP,V2,RP2,AP,AM,R1PI,R1PI2,R1PI3
*
      A2 = A*A
      A3 = A*A2
*
      RP  = 1D0+R
      RP2 = RP*RP
      RP3 = RP*RP2
      RP4 = RP*RP3
*
      V  = 1D0-R
      V2 = V*V
      V3 = V*V2
*
      H1 = 0D0
      H3 = 0D0
*
      IF ((I.EQ.1).AND.(J.EQ.0)) THEN
       
        IF (K.EQ.0) THEN 
         H3= - 1D0/2D0*C*RP2 - 1D0/2D0*RP*V
     + + A*C*RP*V + A*(1D0+R2)
     + - A2/2D0*C*V2 - A2/2D0*RP*V 
         H3 = RP*H3*R1MI2
         H1 = 0D0

        ELSE
        IF (K.EQ.1) THEN
        H3 = - 1D0/2D0*C*RP2 + 1D0/2D0*RP*V
     + +A*C*RP*V - A*(1D0+R2)
     + -A2/2D0*C*V2 + A2/2D0*RP*V 
        H3 = RP*H3*R1MI2
        H1 = 0D0
        ELSE        
         WRITE(*,*) '1. Something is wrong !'
         WRITE(*,*) 'I,J,K=',I,J,K        
        ENDIF
        ENDIF
      ELSE
      IF ((I.EQ.0).AND.(J.EQ.1)) THEN
 
        IF (K.EQ.0) THEN 
         H3 = - 1D0/2D0*C*RP2 + 1D0/2D0*RP*V
     + + A*C*RP*V - A*(1D0+R2)
     + - A2/2D0*C*V2 + A2/2D0*RP*V 
         H3 = -RP*H3*R1MI2
         H1 = 0D0

        ELSE
        IF (K.EQ.1) THEN 
         H3= - 1D0/2D0*C*RP2 - 1D0/2D0*RP*V
     + + A*C*RP*V + A*(1D0+R2)
     + - A2/2D0*C*V2 - A2/2D0*RP*V 
         H3 = -RP*H3*R1MI2
         H1 = 0D0
        ELSE        
         WRITE(*,*) '2. Something is wrong !'
         WRITE(*,*) 'I,J,K=',I,J,K 
        ENDIF
        ENDIF
      ELSE
      IF ((I.EQ.1).AND.(J.EQ.-1)) THEN
       
        IF (K.EQ.0) THEN 
         H3 = -RP2*(1D0-2D0*A*C+A2)
         H3 = H3*R1MI
         H1 = 0D0
        ELSE
        IF (K.EQ.1) THEN
         H3 = -2D0*A*RP*(1D0+R2-C*RP*V)
         H3 = H3*R1MI2
         H1 = 0D0
        ELSE        
         WRITE(*,*) '3. Something is wrong !'
         WRITE(*,*) 'I,J,K=',I,J,K 
        ENDIF
        ENDIF
      ELSE
      IF ((I.EQ.-1).AND.(J.EQ.1)) THEN
       
        IF (K.EQ.0) THEN 
         H3 = -RP2*(1D0+2D0*A*C+A2)
         H3 = H3*R1MI
         H1 = 0D0
        ELSE
        IF (K.EQ.1) THEN
         H3 = -2D0*A*RP*(1D0+R2+C*RP*V)
         H3 = H3*R1MI2
         H1 = 0D0
        ELSE        
         WRITE(*,*) '4. Something is wrong !'
         WRITE(*,*) 'I,J,K=',I,J,K 
        ENDIF
        ENDIF
      ELSE
         WRITE(*,*) '5. Something is wrong !'
         WRITE(*,*) 'I,J,K=',I,J,K 
      ENDIF
      ENDIF
      ENDIF
      ENDIF
C
C      WRITE(*,*) 'CC,BETA,R| ',CC,BETA,R
C      WRITE(*,*) 'A,A2,A3| ',A,A2,A3
C      WRITE(*,*) 'V,V2,V3| ',V,V2,V3
C      WRITE(*,*) 'RP,RP2,RP3,RP4| ', RP,RP2,RP3,RP4
C
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


 
      SUBROUTINE RFBACOL(RADFB)
*=======================================================================
* Introduced by M. JACK on 03/03/99 05:00pm
* ( not used, just for checks !)

       IMPLICIT DOUBLE PRECISION (a-h,o-z)

       COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR

       COMMON /VAR/ A,V,RP,V2,RP2,AP,AM,R1PI,R1PI2,R1PI3

       COMMON /POLYNOMIALS/ G1C,G2C,Y1C,Y2C,
     & G1A,G2A,Y1A,Y2A,YG1CP,YG1CM,YG2CP,YG2CM

       COMMON /XLOGSRM/ XLNRP,XLNR,XLNM,XLLPR2,XLN2

       COMMON /XLOGS2/ XLNG1A,XLNG2A,XLNY1A,XLNY2A,
     & XLNY12P,XLNY12Q,XLNG12P,XLNG12Q 
       COMMON /XLOGSA/ XLNAP,XLNAM,XLNA2,XLNAQ
       COMMON /ALOGS/  ALNA2,ALNAM,ALNAP
       COMMON /YLOGSA/ YLNAS,YLNAD

       COMMON/INDEX/ MF,IIF,JF,KF
C
       A2 = A*A
       AP = 1D0+A
       AM = 1D0-A

       R2 = R*R
       RP = 1D0+R
       V  = 1D0-R
       RP2 = RP*RP
C       
       I = IIF
       J = JF
       K = KF
C
       XLNMR = XLNM-1D0-ALR-2D0*XLN2
       XLNRPR = 2D0*XLNRP-ALR-2D0*XLN2

       RADFB = 0D0

       IF (K.EQ.1) THEN
C----------------------------------------------------------
C        RADFB = (1D0+R2)/RP*(1D0/V*A*YLNAS-1D0/RP*YLNAD)
C     &  + (1D0+R2)/V*AP*AM*XLNAQ + 8D0*A*R2/RP/V
C----------------------------------------------------------

        RADFB = (1D0+R2)*(-Y1A*Y2A*XLNY12Q/RP/RP+(1D0-A2)*XLNAQ) 
     &       + 8D0*A*R2/RP
        RADFB = RADFB/V
   
       ELSE          
       IF (K.EQ.0) THEN
C----------------------------------------------------------
C        RADFB = (1D0+R2)/RP*(1D0/V*A*YLNAD-1D0/RP*YLNAS)
C     &  - (1D0+R2)/V*((1D0-A2)*XLNMR-2D0*(1D0-A)*ALNAP) 
C     &  + 4D0*A*(1D0-A)*R/V
C     &  + 4D0*R*(1D0+R2)/RP/RP/V*(XLNMR+XLNRPR)
C----------------------------------------------------------

        RADFB = (1D0+R2)*((4D0*R*(XLNMR+XLNRPR)
     &  - Y1A*Y2A*XLNY12P)/RP/RP
     &  - AP*AM*(XLNMR-2D0*XLNAP)) 
     &  + 4D0*A*AM*R
        RADFB = RADFB/V

       ELSE
          WRITE(*,*) 'Something is wrong !'
       ENDIF
       ENDIF
C
       END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


 
      SUBROUTINE ACACOL(H3,H1,K)
*=======================================================================
* Introduced by M. JACK on 03/03/99 05:00pm
* ( not used, just for checks !)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR

      COMMON /VAR/ A,V,RP,V2,RP2,AP,AM,R1PI,R1PI2,R1PI3
*
      AM  = 1D0-A

      RP  = 1D0+R
      V  = 1D0-R
      RP2 = RP*RP 
      R2  = R*R
      V2  = V*V
*
      H1 = 0D0
      H3 = 0D0
*
      IF (K.EQ.1) THEN
         H3 = -4D0*A*R2*RP/V2
      ELSE                       
      IF (K.EQ.0) THEN
         H3 = -AM*AM*RP2/V
      ELSE
         WRITE(*,*) 'Something is wrong !'
      ENDIF
      ENDIF

      END




C!!!!!!!!!!!!!!!!!!!!!!!!!!C
C                          C
C       PACKAGE II:        C
C                          C
C!!!!!!!!!!!!!!!!!!!!!!!!!!C


C====================================================================
C 
C     Package containing hard flux functions (radiators) 
C              to the angular distributions 
C
C   $\frac{d{\sigma}}{d{\cos\vartheta}}(R,A(R),\cos\vartheta)$
C 
C   with cuts on one of the final state fermions'
C   minimal scattering angles $\cos\vartheta$, 
C   on their maximal acollinearity angle $\xi^{max}$, 
C   on their minimal energies $E_{min}$, 
C   and/or on their invariant mass squared $s'$ 
C   
C   2 possible options in the code:
C
C   a. general case with all cuts (see above);  
C   
C   b. special case without cuts on acollinearity 
C      and minimal energies ($A=1$);
C      only $s'$ cut (and cut on $\cos\vartheta$)
C 
C====================================================================
C      ** Introduced by M. JACK on 22/06/99 07:00pm **
C====================================================================


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE DLOGCP
C=======================================================================
** Introduced by M. JACK on 22/06/99 07:00pm

       IMPLICIT DOUBLE PRECISION (a-h,o-z)

       COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
       COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM

       COMMON /VAR/ A,V,RP,V2,RP2,AP,AM,R1PI,R1PI2,R1PI3

       COMMON /DLOGC2/ DLNYQ2,DLNYP2,DLNCP2,DLNCM2,DLNCQ2,
     &                 DLNGAMC2

       COMMON /DPOLY/  GAMC,GAMC2,GAMC3,GAMC4
       COMMON /DLOGC/  DLNYQ,DLNYP,DLNCP,DLNCM,DLNCQ,
     &                 DLNGAMC,DLNC2
C-----------------------------------------------------------

       GAMC  = (RP+C*V)/2D0
       GAMC2 = GAMC*GAMC     
       GAMC3 = GAMC*GAMC2
       GAMC4 = GAMC*GAMC3      
C
       DLNYQ = DLNYQ2
       DLNYP = DLNYP2
       DLNCP = DLNCP2
       DLNCM = DLNCM2
       DLNCQ = DLNCQ2
       DLNGAMC = DLNGAMC2

       END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
      SUBROUTINE DLOGCM
C=======================================================================
** Introduced by M. JACK on 22/06/99 07:00pm

       IMPLICIT DOUBLE PRECISION (a-h,o-z)

       COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
       COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM

       COMMON /VAR/ A,V,RP,V2,RP2,AP,AM,R1PI,R1PI2,R1PI3
       
       COMMON /DLOGC1/ DLNYQ1,DLNYP1,DLNCP1,DLNCM1,DLNCQ1,
     &                 DLNGAMC1

       COMMON /DPOLY/  GAMC,GAMC2,GAMC3,GAMC4
       COMMON /DLOGC/  DLNYQ,DLNYP,DLNCP,DLNCM,DLNCQ,
     &                 DLNGAMC,DLNC2
C-----------------------------------------------------------

       GAMC  = (RP-C*V)/2D0
       GAMC2 = GAMC*GAMC     
       GAMC3 = GAMC*GAMC2
       GAMC4 = GAMC*GAMC3      
C
       DLNYQ = DLNYQ1
       DLNYP = DLNYP1
       DLNCP = DLNCP1
       DLNCM = DLNCM1
       DLNCQ = DLNCQ1
       DLNGAMC = DLNGAMC1

       END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE HCUTACOL(H0,H4,H3,H1)
*     ========== ============
* Introduced by M. JACK on 22/06/99 07:00pm

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /DELVAR/ DEL,DELM,DIF,DIF3,ALD,ALDM,ALDPL,ALDMI
     +               ,RD,RDI,RDI2,RDI3,RDM,RDMI,RDMI2,RDMI3
     +               ,ALRD,ALRDM
      COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM
      COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM

      COMMON /VAR/ A,V,RP,V2,RP2,AP,AM,R1PI,R1PI2,R1PI3
      COMMON /XMASSES/ BETA,BETA2

C------------------------------------------------
C With accolinearity cut:

      CALL DRADTIN(H0M,H4M)
      CALL DRADFBIN(H3M,H1M)

C------------------------------------------------
C      WRITE(*,*) '---------------------------'
C      WRITE(*,*) 'I,J,K =',IIF,JF,KF
C      WRITE(*,*) 'R,C,A =',R,C,DIF
C      WRITE(*,*) 'H0,H4 =',H0M,H4M
C      WRITE(*,*) 'H3,H1 =',H3M,H1M        
C------------------------------------------------

      H0=H0M/BETA
      H4=H4M/BETA
      H3=H3M/BETA
      H1=H1M/BETA

C------------------------------------------------
C regularize in s' cut region:
C        IF (R.GE.RCUT) THEN
C          H0=H0-2*TE*COPL2*R1MI
C          H3=H3-TE*4D0*C*R1MI
C          IF (IPHOT2.GE.0) THEN
C           H0=H0+SH2(R,2)*COPL2
C           H3=H3+AH2(R)*2D0*C
C          END IF
C        END IF
C------------------------------------------------

      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C SPECIAL CASE:
C No cuts on final state fermions' maximal 
C acollinearity angle and on minimal energies 
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      SUBROUTINE HCUTC(H0,H4,H3,H1)
*     ========== ============
* Introduced by M. JACK on 22/06/99 07:00pm

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /DELVAR/ DEL,DELM,DIF,DIF3,ALD,ALDM,ALDPL,ALDMI
     +               ,RD,RDI,RDI2,RDI3,RDM,RDMI,RDMI2,RDMI3
     +               ,ALRD,ALRDM
      COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM
      COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM

      COMMON /VAR/ A,V,RP,V2,RP2,AP,AM,R1PI,R1PI2,R1PI3
      COMMON /XMASSES/ BETA,BETA2

      COMMON /INDEXA/ INDA,INDM
      COMMON /INDEX/ MF,IIF,JF,KF

C--------------------------------------------------------
C No accolinearity cut: A=\sqrt{1-4 m_f^2/s'} --> 1

      C0 = -2D0*A*V
      C0I = 2D0*A*RP*V
C
      CT = BETA*C
C    
      CALL DLOGCP 
      CALL DRTFUN0(CT,RTA,RTAI)
      CALL DRFBFUN0(CT,RFBA,RFBAI)
      CALL DLOGCM 
      CALL DRTFUN0(-CT,RTB,RTBI)
      CALL DRFBFUN0(-CT,RFBB,RFBBI)
C
      H0M = RTA+RTB+C0
      H4M = RTAI+RTBI+C0I
      H3M = RFBB-RFBA
      H1M = RFBBI-RFBAI

C------------------------------------------------
c      WRITE(*,*) '---------------------------'
c      WRITE(*,*) 'I,J,K =',IIF,JF,KF
c      WRITE(*,*) 'R,C,A =',R,C,DIF
c      WRITE(*,*) 'H0,H4 =',H0M,H4M
c      WRITE(*,*) 'H3,H1 =',H3M,H1M        
C------------------------------------------------

      H0=H0M/BETA
      H4=H4M/BETA
      H3=H3M/BETA
      H1=H1M/BETA

C------------------------------------------------
C regularize in s' cut region:
C        IF (R.GE.RCUT) THEN
C          H0=H0-2*TE*COPL2*R1MI
C          H3=H3-TE*4D0*C*R1MI
C          IF (IPHOT2.GE.0) THEN
C           H0=H0+SH2(R,2)*COPL2
C           H3=H3+AH2(R)*2D0*C
C          END IF
C        END IF
C------------------------------------------------

      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C-----------------------------------------------
C Initial-state and Interference contributions :
C-----------------------------------------------
 
      SUBROUTINE DRADTIN(DRADT,DRADTI)
*=======================================================================
* Introduced by M. JACK on 22/06/99 07:00pm

       IMPLICIT DOUBLE PRECISION (a-h,o-z)

       COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
       COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM

       COMMON /VAR/ A,V,RP,V2,RP2,AP,AM,R1PI,R1PI2,R1PI3
       COMMON /XMASSES/ BETA,BETA2

       COMMON /XLOGSRM/ XLNRP,XLNR,XLNM,XLLPR2,XLN2

       COMMON /XLOGS2/ XLNG1A,XLNG2A,XLNY1A,XLNY2A,
     & XLNY12P,XLNY12Q,XLNG12P,XLNG12Q 
       COMMON /XLOGSA/ XLNAP,XLNAM,XLNA2,XLNAQ

       COMMON /DPOLY/  GAMC,GAMC2,GAMC3,GAMC4
       COMMON /DLOGC/  DLNYQ,DLNYP,DLNCP,DLNCM,DLNCQ,
     &                 DLNGAMC,DLNC2

       COMMON/INDEX/ MF,IIF,JF,KF

C-----------------------------------------------
     
       I = IIF
       J = JF
       K = KF
 
       CT = BETA*C

C-----------------------------------------------
C Initial-state terms:
       C0 = -2D0*A*V

C Interference terms:       
       C0I = 2D0*A*V*RP

C-----------------------------------------------

       RADT = 0D0
       RADTI = 0D0

       IF ((I.EQ.0).AND.(J.EQ.0)) THEN

          CALL DLOGCP 
          CALL DRTFUN0(CT,F002,F002I)
          CALL DLOGCM 
          CALL DRTFUN0(-CT,F001,F001I)

          RT00 =  F001 + F002 + C0
          DRADT =  RT00

          RT00I =  F001I + F002I + C0I
          DRADTI =  RT00I

       ELSE
       IF ((I.EQ.1).AND.(J.EQ.1)) THEN

          CALL DLOGCP 
          CALL DRTFUN1(CT,F112,F112I)
          CALL DLOGCM 
          CALL DRTFUN1(-CT,F111,F111I)
          
          RT11 =  F111 + F112 + C0
          DRADT =  RT11

          RT11I =  F111I + F112I + C0I
          DRADTI =  RT11I

       ELSE
       IF ((I.EQ.1).AND.(J.EQ.0)) THEN

          CALL DLOGCP 
          CALL DRTFUN1(CT,F112,F112I)
          CALL DLOGCM 
          CALL DRTFUN0(-CT,F001,F001I)

          RT10 =  F001 + F112 + C0
          DRADT =  RT10

          RT10I =  F001I + F112I + C0I
          DRADTI =  RT10I

       ELSE
       IF ((I.EQ.1).AND.(J.EQ.-1)) THEN 

          CALL DLOGCP 
          CALL DRTFUN1(CT,F112,F112I)
          CALL DLOGCM 
          CALL DRTFUN1(-CT,F111,F111I)
 
          RT11M =  F112 - F111 + C0
          DRADT  =  RT11M 
     
          RT11MI =  F112I - F111I + C0I
          DRADTI  =  RT11MI 

       ELSE
       IF ((I.EQ.0).AND.(J.EQ.1)) THEN

          CALL DLOGCP 
          CALL DRTFUN0(CT,F002,F002I)
          CALL DLOGCM 
          CALL DRTFUN1(-CT,F111,F111I)

          RT01 = F111 + F002 + C0
          DRADT =  RT01

          RT01I = F111I + F002I + C0I
          DRADTI =  RT01I

       ELSE
       IF ((I.EQ.-1).AND.(J.EQ.1)) THEN

          CALL DLOGCP 
          CALL DRTFUN1(CT,F112,F112I)
          CALL DLOGCM 
          CALL DRTFUN1(-CT,F111,F111I)

          RT1M1 = F111 - F112  + C0
          DRADT  = RT1M1

          RT1M1I = F111I - F112I  + C0I
          DRADTI  = RT1M1I
       
       ENDIF
       ENDIF
       ENDIF
       ENDIF
       ENDIF
       ENDIF
C
       END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
      SUBROUTINE DRTFUN0(CTH,RT,RTI)
*=======================================================================
* Introduced by M. JACK on 22/06/99 07:00pm

       IMPLICIT DOUBLE PRECISION (a-h,o-z)

       COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR

       COMMON /VAR/ A,V,RP,V2,RP2,AP,AM,R1PI,R1PI2,R1PI3
       COMMON /XMASSES/ BETA,BETA2

       COMMON /XLOGSRM/ XLNRP,XLNR,XLNM,XLLPR2,XLN2

       COMMON /XLOGS2/ XLNG1A,XLNG2A,XLNY1A,XLNY2A,
     & XLNY12P,XLNY12Q,XLNG12P,XLNG12Q 
       COMMON /XLOGSA/ XLNAP,XLNAM,XLNA2,XLNAQ

       COMMON /DPOLY/  GAMC,GAMC2,GAMC3,GAMC4
       COMMON /DLOGC/  DLNYQ,DLNYP,DLNCP,DLNCM,DLNCQ,
     &                 DLNGAMC,DLNC2

C-----------------------------------------------

       A2 = A*A
       A3 = A*A2
C
       C = CTH
       C2 = C*C
       C3 = C*C2

C-----------------------------------------------
C Initial-state term:

c       RT = 
c     &  + 1D0/V * (1D0+C2)*(DLNYP-1D0)
c     &  + DLNYP* (
c     &  - 1D0/2D0/GAMC4* ( 19D0 + 28D0*R + 12D0*R2 + 4D0*R3 
c     &                   + 25D0*C + 7D0*C2 + C3 )
c     &  + 1D0/2D0/GAMC3* ( 53D0 + 40D0*R + 16D0*R2 + 4D0*R3 
c     &                   + 27D0*C + C2 - C3 )
c     &  + 1D0/2D0/GAMC2* (  - 2D0*( 16D0 + 8D0*R + 3D0*R2 + R3 ) 
c     &                   + 1D0 - C + C2 - C3 )
c     &  + 1D0/2D0/GAMC * ( 1D0 - C + C2 - C3 ) )
c     &  + 1D0/3D0/GAMC4* ( 130D0 + 190D0*R + 78D0*R2 + 22D0*R3 
c     &                   + 172D0*C + 49D0*C2 + 7D0*C3 )
c     &  + 1D0/3D0/GAMC3* (- 4D0*( 100D0 + 80D0*R + 29D0*R2 + 6D0*R3 ) 
c     &                   - 239D0*C - 25D0*C2 + 4D0*C3 )
c     &  + 1D0/6D0/GAMC2* ( 604D0 + 310D0*R + 105D0*R2 + 21D0*R3 
c     &                   + 92D0*C - 23D0*C2 + 5D0*C3 )
c     &  + 1D0/6D0/GAMC * ( - 4D0*(19D0 + 7D0*R + 3D0*R2 + R3) 
c     &                   + 30D0*C - 13D0*C2 + 3D0*C3 )
c     &  + A2*( 
c     &  - 1D0/2D0/GAMC2*( 1D0 + 2D0*R + R2 + R3 + C ) 
c     &  + 1D0/2D0/GAMC * (1D0-C) )
c     &  - 1D0/3D0 * (9D0 + 4D0*R + R2)
c     &  + A2/2D0*RP*( 2D0 + R*C )  
c     &  + C/6D0*(R2 + 5D0*R + 12D0)
c
c       RT = RT/R
       

       RT = 
     &   DLNYP*(
     & + 1D0/GAMC4*( 5D0 + 8D0*R + 4D0*R2 + 2D0*R3 + 6D0*C + C2 )
     & - 2D0/GAMC3*(7D0 + 6D0*R + 3D0*R2 + R3 + 3D0*C )
     & + 1D0/GAMC2*( 8D0 + 5D0*R + 2D0*R2 + R3 )  )
     & - 2D0/GAMC4/3D0*( 35D0 + 56D0*R + 28D0*R2 + 11D0*R3 
     &                 + 42D0*C + 7D0*C2)
     & + 2D0/GAMC3*( 110D0/3D0 + 34D0*R + 46D0/3D0*R2 + 4D0*R3 
     &                 + 19D0*C + C2 )
     & + 1D0/GAMC2*(- (169D0/3D0 + 205D0/6D0*R + 14D0*R2 + 7D0/2D0*R3) 
     &          + A2*R/2D0*(1D0+R2) - 6*C + C2) 
     & + 2D0/GAMC/3D0*( 9D0 + 4D0*R + 2D0*R2 + R3 - 3D0*C + C2 )
     & + (5D0 + 3D0*R + R2)/3D0  
     & - C/6D0*(7D0 + 4D0*R + R2) 
     & - A2*(R - RP*V*C/2D0)

       RT = RT/V

C-----------------------------------------------
C Interference term:

       RTI = 
     &  - 4D0*C/V*(DLNCQ-ALR)
     &  + 2D0 *( C - R )*(DLNCQ- ALR ) 
     &  + 2D0*C*RP2*(DLNGAMC - ALR - DLNCM) 
     &  - 4D0*( R2 + R + 2D0 + 2D0*C*R1PI )/GAMC 
     &  + C*RP2*XLNA2 
     &  + C*RP2 
     &  + R2 + 4D0*R + 7D0 + 8D0*R1PI


       END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
      SUBROUTINE DRTFUN1(CTH,RT,RTI)
*=======================================================================
* Introduced by M. JACK on 22/06/99 07:00pm

       IMPLICIT DOUBLE PRECISION (a-h,o-z)

       COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR

       COMMON /VAR/ A,V,RP,V2,RP2,AP,AM,R1PI,R1PI2,R1PI3
       COMMON /XMASSES/ BETA,BETA2

       COMMON /XLOGSRM/ XLNRP,XLNR,XLNM,XLLPR2,XLN2

       COMMON /XLOGS2/ XLNG1A,XLNG2A,XLNY1A,XLNY2A,
     & XLNY12P,XLNY12Q,XLNG12P,XLNG12Q 
       COMMON /XLOGSA/ XLNAP,XLNAM,XLNA2,XLNAQ

       COMMON /DPOLY/  GAMC,GAMC2,GAMC3,GAMC4
       COMMON /DLOGC/  DLNYQ,DLNYP,DLNCP,DLNCM,DLNCQ,
     &                 DLNGAMC,DLNC2

C-----------------------------------------------

       A2 = A*A
       A3 = A*A2

       C = CTH
       C2 = C*C
       C3 = C*C2

C-----------------------------------------------
C Initial-state term:

c       RT = 
c     &  + 1D0/V*DLNYQ * ( 1D0 + C2 )
c     &  + DLNYQ*(
c     &  + 1D0/2D0/GAMC4* (-( 19D0 + 28D0*R + 12D0*R2 + 4D0*R3 ) 
c     &                    - 25D0*C - 7D0*C2 - C3 )
c     &  + 1D0/2D0/GAMC3* ( 53D0 + 40D0*R + 16D0*R2 + 4D0*R3 
c     &                    + 27D0*C + C2 - C3 )     
c     &  + 1D0/2D0/GAMC2* (- ( 32D0 + 16D0*R + 6D0*R2 + 2D0*R3 ) 
c     &                    + 1D0 - C + C2 - C3 )
c     &  + 1D0/2D0/GAMC* ( 1D0 - C + C2 - C3 ) )
c     &  + A*( 
c     &  + 1D0/GAMC3 * ( 5D0 + 8D0*R + 4D0*R2 + 2D0*R3 
c     &                + 6D0*C + C2 )
c     &  + 1D0/GAMC2 * ( -10D0 - 6D0*R - 3D0*R2 - R3 - C + C2 )
c     &  + 1D0/2D0/GAMC * ( R*(1D0+R2) - 2D0*( C - C2 )) )
c     &  + A3/6D0/GAMC *R*(1D0+R2) 
c     &  - 1D0/6D0*A3*R*(1D0 - C + R*( 1D0+C ) )
c     &  + 1D0/2D0*A*R*(1D0 + C + R*( 1D0-C ))
c     &  - 2D0*A*C*RP
c
c       RT = RT/R


       RT = 
     & DLNYQ*( 
     & + 1D0/GAMC4*( 5D0 + 8D0*R + 4D0*R2 + 2D0*R3 
     &             + 6D0*C + C2 )
     & - 2D0/GAMC3*(7D0 + 6D0*R + 3D0*R2 + R3 + 3D0*C )
     & + 1D0/GAMC2*( 8D0 + 5D0*R + 2D0*R2 + R3 )  )
     & - 2D0*A/GAMC3*(1D0 + 2D0*R + R2 + R3 + C )
     & + A/GAMC2*( 4D0 + 3D0*R + 2D0*R2 + R3 )
     & + A/GAMC/2D0*(1D0+A2/3D0)*V*(1D0+R2)
     & + A/2D0*RP*(C*RP + V)
     & + A3/6D0*V*(C*V - RP) 

       RT = RT/V

C-----------------------------------------------
C Interference term:

       RTI =
     &  - 4D0*C/V*XLNAQ 
     &  + XLNAQ * ( C*RP2 + 2D0*(C - R) )
     &  - A*V*( 1D0 - R - C*RP )


       END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C-----------------------------------------------
C Initial-state and Interference Contributions :
C-----------------------------------------------
 
      SUBROUTINE DRADFBIN(DRADFB,DRADFBI)
*=======================================================================
* Introduced by M. JACK on 22/06/99 07:00pm

       IMPLICIT DOUBLE PRECISION (a-h,o-z)

       COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
       COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM

       COMMON /VAR/ A,V,RP,V2,RP2,AP,AM,R1PI,R1PI2,R1PI3
       COMMON /XMASSES/ BETA,BETA2

       COMMON /XLOGSRM/ XLNRP,XLNR,XLNM,XLLPR2,XLN2

       COMMON /XLOGS2/ XLNG1A,XLNG2A,XLNY1A,XLNY2A,
     & XLNY12P,XLNY12Q,XLNG12P,XLNG12Q 
       COMMON /XLOGSA/ XLNAP,XLNAM,XLNA2,XLNAQ

       COMMON /DPOLY/  GAMC,GAMC2,GAMC3,GAMC4
       COMMON /DLOGC/  DLNYQ,DLNYP,DLNCP,DLNCM,DLNCQ,
     &                 DLNGAMC,DLNC2

       COMMON/INDEX/ MF,IIF,JF,KF
C----------------------------------------------------------

       I = IIF
       J = JF
       K = KF

       CT = BETA*C

C----------------------------------------------------------
C Initial-state terms:

       G1 = 0D0

C----------------------------------------------------------
C Interference terms:
 
       G1I = 2D0*CT*(2D0*R*XLNG12Q + V*RP*XLNAQ 
     &      - A*V*RP)

C----------------------------------------------------------

       DRADFB  = 0D0
       DRADFBI = 0D0

       IF ((I.EQ.0).AND.(J.EQ.0)) THEN

          CALL DLOGCP 
          CALL DRFBFUN0(CT,G002,G002I)
          CALL DLOGCM 
          CALL DRFBFUN0(-CT,G001,G001I)

          RFB00  =  G001 - G002 
          DRADFB =  RFB00

          RFB00I  =  G001I - G002I        
          DRADFBI =  RFB00I

       ELSE
       IF ((I.EQ.1).AND.(J.EQ.1)) THEN

          CALL DLOGCP 
          CALL DRFBFUN1(CT,G112,G112I)
          CALL DLOGCM 
          CALL DRFBFUN1(-CT,G111,G111I)

          RFB11  =  G111 - G112
          DRADFB =  RFB11

          RFB11I  =  G111I - G112I 
          DRADFBI =  RFB11I

       ELSE
       IF ((I.EQ.1).AND.(J.EQ.0)) THEN

          CALL DLOGCP 
          CALL DRFBFUN1(C,G112,G112I)
          CALL DLOGCM 
          CALL DRFBFUN0(-C,G001,G001I)

          RFB10 = G001 - G112 
          DRADFB  = RFB10

          RFB10I = G001I - G112I 
          DRADFBI  = RFB10I

       ELSE
       IF ((I.EQ.1).AND.(J.EQ.-1)) THEN 

          CALL DLOGCP 
          CALL DRFBFUN1(CT,G112,G112I)
          CALL DLOGCM 
          CALL DRFBFUN1(-CT,G111,G111I)

          RFB11M = - G111 - G112 + G1 
          DRADFB   = RFB11M

          RFB11MI = - G111I - G112I + G1I 
          DRADFBI   = RFB11MI

       ELSE
       IF ((I.EQ.0).AND.(J.EQ.1)) THEN

          CALL DLOGCP 
          CALL DRFBFUN0(CT,G002,G002I)
          CALL DLOGCM 
          CALL DRFBFUN1(-CT,G111,G111I)

          RFB01 = G111 - G002  
          DRADFB  = RFB01

          RFB01I = G111I - G002I
          DRADFBI  = RFB01I

       ELSE
       IF ((I.EQ.-1).AND.(J.EQ.1)) THEN

          CALL DLOGCP 
          CALL DRFBFUN1(CT,G112,G112I)
          CALL DLOGCM
          CALL DRFBFUN1(-CT,G111,G111I)

          RFB1M1 = G111 + G112 + G1
          DRADFB   = RFB1M1

          RFB1M1I = G111I + G112I + G1I
          DRADFBI   = RFB1M1I

       ENDIF
       ENDIF
       ENDIF
       ENDIF
       ENDIF
       ENDIF
C
       END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
      SUBROUTINE DRFBFUN0(CTH,RFB,RFBI)
*=======================================================================
* Introduced by M. JACK on 22/06/99 07:00pm

       IMPLICIT DOUBLE PRECISION (a-h,o-z)

       COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR

       COMMON /VAR/ A,V,RP,V2,RP2,AP,AM,R1PI,R1PI2,R1PI3
       COMMON /XMASSES/ BETA,BETA2

       COMMON /XLOGSRM/ XLNRP,XLNR,XLNM,XLLPR2,XLN2

       COMMON /XLOGS2/ XLNG1A,XLNG2A,XLNY1A,XLNY2A,
     & XLNY12P,XLNY12Q,XLNG12P,XLNG12Q 
       COMMON /XLOGSA/ XLNAP,XLNAM,XLNA2,XLNAQ

       COMMON /DPOLY/  GAMC,GAMC2,GAMC3,GAMC4
       COMMON /DLOGC/  DLNYQ,DLNYP,DLNCP,DLNCM,DLNCQ,
     &                 DLNGAMC,DLNC2

C----------------------------------------------------------

       A2 = A*A
       A3 = A*A2

       C = CTH
       C2 = C*C
       C3 = C*C2

C----------------------------------------------------------
C Initial-state terms:

c       RFB =
c     &  - 1D0/V * ( 2D0*C*(DLNYP-1D0))
c     &  + DLNYP*(
c     &  + 1D0/GAMC3*( 5D0 + 8D0*R + 4D0*R2 + 2D0*R3 
c     &               + 6D0*C + C2 )
c     &  + 1D0/GAMC2*(-(10D0 + 6D0*R + 3D0*R2 + R3) 
c     &               - C + C2 )
c     &  - 1D0/GAMC*C*(1D0-C) )
c     &  + 2D0/GAMC3*(-(10D0 + 16D0*R + 8D0*R2 + 3D0*R3 ) 
c     &               - 2D0*(6D0*C + C2) )
c     &  + 2D0/GAMC2*( 24D0 + 18D0*R + 8D0*R2 + 2D0*R3 
c     &               + 7D0*C - C2 )
c     &  + 1D0/2D0/GAMC*(- (26D0 + 9D0*R + 4D0*R2 + R3) 
c     &               + 2D0*(6D0*C - C2))
c     &  + A2/2D0*R*(1D0+R2)/GAMC 
c     &  - 1D0/2D0*(12D0 + 5D0*R + R2)
c     &  - 1D0/2D0*A2*R*RP
c
c       RFB = RFB/R


       RFB = 
     &  DLNYP*(-2D0/GAMC3* ( 1D0 + 2D0*R + R2 + R3 + C )
     &  + 1D0/GAMC2*( 4D0 + 3D0*R + 2D0*R2 + R3 ))
     &  + 2D0/GAMC3*( 4D0 + 8D0*R + 5D0*R2 + 3D0*R3 + 4D0*C )
     &  - 4D0/GAMC2*( 5D0 + 5D0*R + 3D0*R2 + R3 + C )
     &  + 1D0/GAMC/2D0*( 11D0 + 5D0*R + 3D0*R2 + R3 
     &                 + A2*V*(1D0+R2) - 4D0*C )
     &  + (7D0 + 4D0*R + R2 - A2*V*RP)/2D0
 
       RFB = RFB/V

C----------------------------------------------------------
C Interference terms:

       RFBI = 
     &  + 2D0/V* ( 1D0 + C2 )*( DLNCQ - ALR) 
     &  - 2D0*(DLNGAMC - ALR - DLNCM) 
     &       * ( 1D0 + R2 - C*R2 + C2*(1D0+R+R2) )
     &  - 2D0*(DLNGAMC - DLNCP) *( C - R )
     &  - XLNA2* ( 1D0 - R + R2 + C2*(1D0+R+R2) )
     &  - 2D0*C* XLNAP * V * RP
     &  - 2D0*C*R *XLNG12Q 
     &  + 2D0/GAMC2 *( 3D0 + 6D0*C*R1PI + 2D0*C2*R1PI2 + RP2 )
     &  + 4D0/GAMC * ( 3D0 + R + C*(5D0*R1PI + 6D0*R1PI2) 
     &             + 2D0*C2*(R1PI2 + 2D0*R1PI3) )
     &  + A*C*V*RP
     &  - A2/4D0*(1D0+C2)*V2  
     &  - 5D0/4D0*C2*RP2 
     &  - 8D0*C*R1PI*( 1D0 + 2D0*R1PI + 2D0*R1PI2 ) 
     &  + C*(R2 - 4D0*R - 9D0)
     &  - 4D0*( 3D0 + 6D0*R1PI + 5D0*R1PI2 ) - 1D0/4D0*RP2


       END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
      SUBROUTINE DRFBFUN1(CTH,RFB,RFBI)
*=======================================================================
* Introduced by M. JACK on 22/06/99 07:00pm

       IMPLICIT DOUBLE PRECISION (a-h,o-z)

       COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR

       COMMON /VAR/ A,V,RP,V2,RP2,AP,AM,R1PI,R1PI2,R1PI3
       COMMON /XMASSES/ BETA,BETA2

       COMMON /XLOGSRM/ XLNRP,XLNR,XLNM,XLLPR2,XLN2

       COMMON /XLOGS2/ XLNG1A,XLNG2A,XLNY1A,XLNY2A,
     & XLNY12P,XLNY12Q,XLNG12P,XLNG12Q 
       COMMON /XLOGSA/ XLNAP,XLNAM,XLNA2,XLNAQ

       COMMON /DPOLY/  GAMC,GAMC2,GAMC3,GAMC4
       COMMON /DLOGC/  DLNYQ,DLNYP,DLNCP,DLNCM,DLNCQ,
     &                 DLNGAMC,DLNC2

C----------------------------------------------------------

       A2 = A*A
       A3 = A*A2

       C = CTH
       C2 = C*C
       C3 = C*C2
 
C----------------------------------------------------------
C Initial-state terms:

c       RFB = 
c     &  - 2D0*C/V*DLNYQ 
c     &  + DLNYQ*(
c     &  + 1D0/GAMC3*( 5D0 + 8D0*R + 4D0*R2 + 2D0*R3 
c     &              + 6D0*C + C2 )
c     &  + 1D0/GAMC2*(- (10D0 + 6D0*R + 3D0*R2 + R3) 
c     &               - C + C2 )
c     &  - 1D0/GAMC*C*(1D0 - C) )
c     &  + 2D0*A*(- 1D0/GAMC2*( 1D0 + 2D0*R + R2 + R3 + C )
c     &           + 1D0/GAMC*( 1D0 - C ) )
c     &  + A*(R2 + 3D0*R + 4D0)
c
c       RFB = RFB/R


       RFB = 
     &  DLNYQ*(- 2D0/GAMC3*(1D0 + 2D0*R + R2 + R3 + C )
     &         + 1D0/GAMC2*(4D0 + 3D0*R + 2D0*R2 + R3 ) )        
     & + 2D0*A*R*( 1D0 + R2 )/GAMC2 
     & - A*RP2

       RFB = RFB/V

C----------------------------------------------------------
C Interference terms:

       RFBI = 
     &  + 2D0 * (1D0 + C2 )/V*XLNAQ 
     &  - 2D0*C*R * XLNG12Q 
     &  + (2D0*C*R2 - (1D0+C2)*(1D0+R+R2) )*XLNAQ
     &  + A*V*( 2D0*C*R - (1D0+3D0*C2)*RP/2D0 )


       END


