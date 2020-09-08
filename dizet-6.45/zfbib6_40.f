      SUBROUTINE ZCUT(INTRF,IFAST,INDF,
     &      S,AMZ,GAMZ,WIDTHS,SW2,NPAR,ZPAR,SIGBRN,SIGQED,AFBBRN,AFBQED)
*     ==================================================================
************************************************************************
*                                                                      *
*  PROGRAM CALCULATES CROSS-SECTION                                    *
*  FOR ANNIHILATION OF E+E-  INTO  F+F-                                *
*  WITH CUTS ON ACOLLINEARITY AND FERMIONS' ENERGY                     *
*                          OR                                          *
*  WITH CUT ON S'(SQUARE OF FERMIONS' INVARIANT MASS)                  *
*  THE CROSS-SECTION IS INTEGRATED OVER ANGLE BETWEEN E+ AND F+        *
*  FROM THET1 TO THET2                                                 *
*                                                                      *
************************************************************************
      IMPLICIT COMPLEX*16(X)
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
*
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
      COMMON /CORINT/ CORINT
      COMMON /CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
      COMMON /CDZBOX/ QE,QF,ALAM1,ALAM2,HELI1,HELI2,SC,VE,VF,CRFAC
      COMMON /INTRFS/ INTRFC
      COMMON /POLHEL/ COMB1,COMB2,HOMB1,HOMB2
      COMMON /CDZIBF/ IBFLA
      COMMON /ZFCHMS/ ALLCH(0:11),ALLMS(0:11)
      COMMON /SMATRS/ AMZS,GAMZS,RESR,RESI,RES0,RES1,RES2,RESG,ISMA
      COMMON /CDZRUN/ CMQRUN(8)
*
      DIMENSION NPAR(30),ZPAR(31),WIDTHS(0:11)
*
      EXTERNAL DZEWBX
*
* Setting of IBFLA for DZEWBX (per-channel mode)
*
      IBFLA=0
      IF(INDF.EQ.9) IBFLA=1
*   
* OTHER FLAGS
*
      INTRFC=INTRF
      IWEAK =NPAR(1)
      INTERF=NPAR(8)
      IFINAL=NPAR(9)
      IPHOT2=NPAR(10)
      IRCUT =NPAR(11)
      IAFB  =NPAR(13)
      IBORN =NPAR(14)
      ANG1  =ZPAR(29)
      ANG2  =ZPAR(30)
      ZPRIPP=ZPAR(31)
*
***********************  FLAGS *****************************************
*  IAFB =1 CALCULATION OF CROSS-SECTION AND FORWARD-BACKWARD ASYMMETRY
*       =0 CALCULATION CROSS-SECTION ONLY
*-----------------------------------------------------------------------
*   IRCUT - CHOICE BETWEEN : S'-CUT AND MUONS' ACOL.+ENERGY CUT  CASES
***********************  SETS  *****************************************
*-----------------------------------------------------------------------
* S-INDEPENDENT INITIALIZATION
*-----------------------------------------------------------------------
      CALL EWINIT(INTRF,INDF,AMZ,GAMZ,SW2)
*-----------------------------------------------------------------------
* SET S-DEPENDENT FUNCTIONS
*-----------------------------------------------------------------------
      CALL SETFUN(S)
      CALL SETCUT(S,ZPAR,SPR,ECUT,ACUT)
*-----------------------------------------------------------------------
* SET RUNNING COUPLING CONSTANTS
*-----------------------------------------------------------------------
      CALL EWCOUP(INTRF,INDF,S)
*-----------------------------------------------------------------------
* CALCULATION OF CROSS-SECTION INTEGRATED FROM THET TO PI-THET
*-----------------------------------------------------------------------
      ISYM=0
      IF(ABS(ANG1+ANG2-180D0).LT.1D-4.OR.IFAST.EQ.1) ISYM=1
*
      IF(INDF.EQ.9.AND.S.LE.4*AMQ(6)**2) THEN
        SBRN=0D0
        ABRN=0D0
        SQED=0D0
        AQED=0D0
      ELSE
        IF (IFAST.EQ.1) THEN
         CALL SFAST(INDF,S,SPR,ZPRIPP,ECUT,ACUT,SBRN,ABRN,SQED,AQED)
        ELSE
         CALL SCUT(INDF,S,ANG2,SPR,ZPRIPP,ECUT,ACUT,SBRN,ABRN,SQED,AQED)
        ENDIF
      ENDIF
*
      IF (ISYM.EQ.1)THEN
        SIGQED=2D0*SQED
        SIGBRN=2D0*SBRN
        IF (IAFB.EQ.1) THEN
          AFBQED=2D0*AQED
          AFBBRN=2D0*ABRN
        ENDIF
      ELSE
        SIGQED=SQED+AQED
        SIGBRN=SBRN+ABRN
        IF (IAFB.EQ.1) THEN
          AFBQED=SIGQED
          AFBBRN=SIGBRN
        ENDIF
        IF(INDF.EQ.9.AND.S.LE.4*AMQ(6)**2) THEN
          SBRN=0D0
          ABRN=0D0
          SQED=0D0
          AQED=0D0
        ELSE
         CALL SCUT(INDF,S,ANG1,SPR,ZPRIPP,ECUT,ACUT,SBRN,ABRN,SQED,AQED)
        ENDIF
        SIGQED=SIGQED-SQED-AQED
        SIGBRN=SIGBRN-SBRN-ABRN
        IF (IAFB.EQ.1) THEN
         AFBQED=AFBQED+SQED+AQED
         AFBBRN=AFBBRN+SBRN+ABRN
        ENDIF
      ENDIF
*
* Adding of EW-boxes for IBOX=1
*
      IBOX = NPAR(4)
      SBOX = 0D0
      ABOX = 0D0
*
      IF(IBOX.NE.1) GOTO 90
      IF((INTRF.EQ.2.OR.ISMA.GT.0).AND.INDF.EQ.10) THEN 
* Special treatment of EW-boxes for ZUXSEC and S-Matrix for INDF=10
        QE =ZPAR(1)
        QEM=DABS(QE)
        VE=1D0-4D0*SW2*QEM
        SC=S
        CNANOB=.38937966D6
        CRFAC=CNANOB
        ALAM1=1D0
        ALAM2=0D0
        HELI1=1D0
        HELI2=0D0
        IF (IFAST.EQ.1) THEN
          C1=-1D0
          C2=+1D0
        ELSE
          C1=COS(ANG1*PI/180D0)
          C2=COS(ANG2*PI/180D0)
        ENDIF
        eps = 1d-4
        amin= c1+1d-8
        amax= c2-1d-8
        st  = (amax-amin)/2
*
* u-box
*
        IBFLA=0
        QF =+2D0/3
        QFM=DABS(QF)
        VF=1D0-4D0*SW2*QFM
        call simps(amin,0d0,st,eps,1d-20,DZEWBX,ac,resb,resb2,resb3)
        call simps(0d0,amax,st,eps,1d-20,DZEWBX,ac,resf,resf2,resf3)
        UNORM=DSQRT(1D0-4D0* ALLMS(4)**2/SC)
     &       +DSQRT(1D0-4D0*CMQRUN(3)**2/SC)
        SBOX = SBOX+UNORM*(RESF+RESB)
*
* d-box
*
        QF =-1D0/3
        QFM=DABS(QF)
        VF=1D0-4D0*SW2*QFM
        call simps(amin,0d0,st,eps,1d-20,DZEWBX,ac,resb,resb2,resb3)
        call simps(0d0,amax,st,eps,1d-20,DZEWBX,ac,resf,resf2,resf3)
        DNORM=DSQRT(1D0-4D0*ALLMS(5)**2/SC)
     &       +DSQRT(1D0-4D0*ALLMS(7)**2/SC)
        SBOX = SBOX+DNORM*(RESF+RESB)
*
* b-box
*
        IBFLA=1
        call simps(amin,0d0,st,eps,1d-20,DZEWBX,ac,resb,resb2,resb3)
        call simps(0d0,amax,st,eps,1d-20,DZEWBX,ac,resf,resf2,resf3)
        DNORM=DSQRT(1D0-4D0*CMQRUN(6)**2/SC)
        SBOX = SBOX+DNORM*(RESF+RESB)
      ELSE
        QE   = ZPAR(1)
        QF   = ALLCH(INDF)
        IF(INDF.EQ.9) THEN
          AMF= CMQRUN(6)
        ELSEIF(INDF.EQ.6) THEN
          AMF= CMQRUN(3)
        ELSE
          AMF= ALLMS(INDF)
        ENDIF
        ALAM1=COMB1
        ALAM2=COMB2
        HELI1=HOMB1
        HELI2=HOMB2
        QFM=DABS(QF)
        QEM=DABS(QE)
        VE=1D0-4D0*SW2*QEM
        VF=1D0-4D0*SW2*QFM
        SC=S
        FMUSQ= DSQRT(1D0-4D0*AMF*AMF/SC)
        CNANOB=.38937966D6
        CRFAC=CNANOB*FMUSQ
        IF (IFAST.EQ.1) THEN
          C1=-1D0
          C2=+1D0
        ELSE
          C1=COS(ANG1*PI/180D0)
          C2=COS(ANG2*PI/180D0)
        ENDIF
        eps = 1d-4
        amin = c1+1d-8
        amax = c2-1d-8
        st   = (amax-amin)/2
        call simps(amin,0d0,st,eps,1d-20,DZEWBX,ac,resb,resb2,resb3)
        call simps(0d0,amax,st,eps,1d-20,DZEWBX,ac,resf,resf2,resf3)
        sbox = resf+resb
        abox = resf-resb
      ENDIF
*
 90   continue
*
      sigbrn=sigbrn + sbox*corint
      afbbrn=afbbrn + abox*corint
      sigqed=sigqed + sbox*corint
      afbqed=afbqed + abox*corint
*
      IF(IAFB.EQ.1.AND.SIGBRN.NE.0D0) THEN
        AFBBRN=AFBBRN/SIGBRN
      ENDIF
*
      IF(IAFB.EQ.1.AND.SIGQED.NE.0D0) THEN
        AFBQED=AFBQED/SIGQED
      ENDIF
*                                                               END ZCUT
      END
 
      SUBROUTINE ZANCUT(INTRF,IFAST,INDF,COSTHR,
     &                    S,AMZ,GAMZ,WIDTHS,SW2,NPAR,ZPAR,SIGBRN,SIGQED)
*     ==================================================================
************************************************************************
*     IT IS A FORMAL ANALOG OF ZCUT                                    *
*                                                                      *
*  PROGRAM CALCULATES DIFFERENTIAL IN COSTHR CROSS-SECTION             *
*  FOR ANNIHILATION OF E+E-  INTO  F+F-                                *
*  WITH CUTS ON ACOLLINEARITY AND FERMIONS' ENERGY                     *
*                          OR                                          *
*  WITH CUT ON S'(SQUARE OF FERMIONS' INVARIANT MASS)                  *
*  THE CROSS-SECTION IS DIFFERENTIAL OVER ANGLE BETWEEN E+ AND F+      *
*                                                                      *
************************************************************************
      IMPLICIT COMPLEX*16(X)
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
*
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
      COMMON /CORINT/ CORINT
      COMMON /CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
      COMMON /CDZBOX/ QE,QF,ALAM1,ALAM2,HELI1,HELI2,SC,VE,VF,CRFAC
      COMMON /INTRFS/ INTRFC
      COMMON /POLHEL/ COMB1,COMB2,HOMB1,HOMB2
      COMMON /COCOST/ COCOST
      COMMON /CDZIBF/ IBFLA
      COMMON /ZFCHMS/ ALLCH(0:11),ALLMS(0:11)
      COMMON /SMATRS/ AMZS,GAMZS,RESR,RESI,RES0,RES1,RES2,RESG,ISMA
      COMMON /CDZRUN/ CMQRUN(8)
*
      DIMENSION NPAR(30),ZPAR(31),WIDTHS(0:11)
*
      COCOST=COSTHR
*
* Setting of IBFLA for DZEWBX (per-channel mode)
*
      IBFLA=0
      IF(INDF.EQ.9) IBFLA=1
*   
* OTHER FLAGS
*
      INTRFC=INTRF
      IWEAK =NPAR(1)
      INTERF=NPAR(8)
      IFINAL=NPAR(9)
      IPHOT2=NPAR(10)
      IRCUT =NPAR(11)
      IAFB  =NPAR(13)
      IBORN =NPAR(14)
      ANG1  =ZPAR(29)
      ANG2  =ZPAR(30)
      ZPRIPP=ZPAR(31)
*
************************ FLAGS *****************************************
*  IAFB =1 CALCULATION OF CROSS-SECTION AND FORWARD-BACKWARD ASYMMETRY
*       =0 CALCULATION CROSS-SECTION ONLY
*-----------------------------------------------------------------------
*   IRCUT - CHOICE BETWEEN : S'-CUT AND MUONS' ACOL.+ENERGY CUT  CASES
***********************  SETS  *****************************************
*-----------------------------------------------------------------------
* S-INDEPENDENT INITIALIZATION
*-----------------------------------------------------------------------
      CALL EWINIT(INTRF,INDF,AMZ,GAMZ,SW2)
*-----------------------------------------------------------------------
* SET S-DEPENDENT FUNCTIONS
*-----------------------------------------------------------------------
      CALL SETFUN(S)
      CALL SETCUT(S,ZPAR,SPR,ECUT,ACUT)
*-----------------------------------------------------------------------
* SET RUNNING COUPLING CONSTANTS
*-----------------------------------------------------------------------
      CALL EWCOUP(INTRF,INDF,S)
*-----------------------------------------------------------------------
* CALCULATION OF THE DIFFERENTIAL IN COSTHR CROSS-SECTION
*-----------------------------------------------------------------------
*
      IF(INDF.EQ.9.AND.S.LE.4*AMQ(6)**2) THEN
        SIGBRN=0D0
        SIGQED=0D0
      ELSE
        CALL COSCUT(INDF,S,COSTHR,SPR,ZPRIPP,ECUT,ACUT,SIGBRN,SIGQED)
      ENDIF
*
* Adding of EW-boxes for IBOX=1
*
      IBOX = NPAR(4)
      SBOX = 0D0
*
      IF(IBOX.NE.1) GOTO 90
      IF((INTRF.EQ.2.OR.ISMA.GT.0).AND.INDF.EQ.10) THEN 
* Special treatment of EW-boxes for ZUXSEC and S-Matrix for INDF=10
        QE =ZPAR(1)
        QEM=DABS(QE)
        VE=1D0-4D0*SW2*QEM
        SC=S
        CNANOB=.38937966D6
        CRFAC=CNANOB
        ALAM1=1D0
        ALAM2=0D0
        HELI1=1D0
        HELI2=0D0
*
* u-box
*
        IBFLA=0
        QF =+2D0/3
        QFM=DABS(QF)
        VF=1D0-4D0*SW2*QFM
        UNORM=DSQRT(1D0-4D0*ALLMS(4)**2/SC)
     &       +DSQRT(1D0-4D0*ALLMS(6)**2/SC)
        SBOX = SBOX+UNORM*DZEWBI(COSTHR)
*
* d-box
*
        QF =-1D0/3
        QFM=DABS(QF)
        VF=1D0-4D0*SW2*QFM
        DNORM=DSQRT(1D0-4D0*ALLMS(5)**2/SC)
     &       +DSQRT(1D0-4D0*ALLMS(7)**2/SC)
        SBOX = SBOX+DNORM*DZEWBI(COSTHR)
*
* b-box
*
        IBFLA=1
        DNORM=DSQRT(1D0-4D0*ALLMS(9)**2/SC)
        SBOX = SBOX+DNORM*DZEWBI(COSTHR)
      ELSE
        QE   = ZPAR(1)
        QF   = ALLCH(INDF)
        IF(INDF.EQ.9) THEN
          AMF= CMQRUN(6)
        ELSEIF(INDF.EQ.6) THEN
          AMF= CMQRUN(3)
        ELSE
          AMF= ALLMS(INDF)
        ENDIF
        ALAM1=COMB1
        ALAM2=COMB2
        HELI1=HOMB1
        HELI2=HOMB2
        QFM=DABS(QF)
        QEM=DABS(QE)
        VE=1D0-4D0*SW2*QEM
        VF=1D0-4D0*SW2*QFM
        SC=S
        FMUSQ =DSQRT(1D0-4D0*AMF*AMF/SC)
        CNANOB=.38937966D6
        CRFAC=CNANOB*FMUSQ
        SBOX = DZEWBI(COSTHR)
      ENDIF
 90   CONTINUE
*
      SIGBRN=SIGBRN + SBOX
      SIGQED=SIGQED + SBOX
*                                                             END ZANCUT
      END
 
      SUBROUTINE SETCUT(S,ZPAR,SPR,ECUT,ACUT)
*     ========== ============================
************************************************************************
* Modified by M. Jack on 05/03/99 03:30pm                              *
************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
*
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /MASSZ / AME,AMF,AME2,AMF2
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
      COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT
      COMMON /FLAGZP/ ISRPPR,IFSPPR,IFUNAN
      COMMON /SOFTPR/ SOFTPR,SFPRFB
*
      DIMENSION ZPAR(31)
*
      PCUT=ZPAR(31)/S
*
      IF(IRCUT.EQ.1) THEN
        SPR =S*(1D0-ZPAR(26))
        ECUT=0D0
        ACUT=0D0
          ELSE
        SPR =0D0
        ACUT=ZPAR(27)
        ECUT=ZPAR(28)
      ENDIF
*
      RMIN =MAX(1D-10,4D0*AMF2/S)
*
      IF(IRCUT.EQ.1) THEN
        RCUT=SPR/S
        IF(IBORN.EQ.0.AND.RCUT.GE.1D0) STOP
     & 'WRONG S-PRIME CUT!'
        RCUT=MAX(RMIN,RCUT)
        RECUT1=RCUT
        RECUT2=RCUT
        RACUT =RCUT
      ELSE
        IF(IBORN.EQ.0.AND.(ACUT.LT.0D0.OR.ACUT.GE.180D0)) STOP 
     & 'WRONG ACOLLINEARITY CUT!'
        SINAC=SIN(ACUT/2D0*PI/180)
        SINAC2=SINAC**2
        COSAC2=1D0-SINAC2
        RACUT=(1D0-SINAC)/(1D0+SINAC)
* KINEMATICAL LIMIT
        ECUT=MAX(ECUT,AMF)
        RECUT2=2D0*ECUT/SQRT(S)
        IF(IBORN.EQ.0.AND.RECUT2.GE.1D0) STOP
     & 'WRONG ENERGY CUTS!'
        IF(RACUT.LE.2D0*RECUT2-1) THEN
          RECUT1=2D0*RECUT2-1D0
        ELSE
          RECUT1=RECUT2*(1D0-SINAC2/(1D0-RECUT2*COSAC2))
        ENDIF
        RECUT1=MAX(RMIN,RECUT1)
        RECUT2=MAX(RMIN,RECUT2)
        RACUT =MAX(RMIN, RACUT)
        RCUT  =MAX(RACUT,RECUT2)
        RECUTA=MAX(RECUT2,(1D0-SIN(ACUT*PI/180))/(1D0+SIN(ACUT*PI/180)))
      ENDIF
      IF (IFINAL.EQ.1) CALL SETFIN
*                                                             END SETCUT
      END

      SUBROUTINE 
     &     SCUT(INDF,SARG,THET,SPR,ZPRIPP,ECUT,ACUT,SBRN,ABRN,SQED,AQED)
*======================================================================*
* Modified by M. Jack on 03/03/99 05:00pm: comment added               *
************************************************************************
* ROUTINE CALCULATES INTEGRALS OVER SCATERING ANGLES                   *
* FROM 0 TO THAT OF COS-EVEN AND COS-ODD PARTS OF ANGULAR DISTRIBUTION *
*   FOR ANNIHILATION OF E+E- INTO F^+ F^- AND N*GAM                    *
*   S'-CUT   OR   COMBINED ACOLINEARITY + MUON ENERGIES CUTS           *
*                                                                      *
*   I N P U T:                                                         *
* - INDF - THE STANDART ZFITTER 'CHANNEL INDEX'                        *
* - SARG - S-INVARIANT IN GEV**2                                       *
* - THET - LIMIT OF COS INTEGRATION (ACCEPTANCE) IN DEGREE             *
* - SPR  - MIN OF FINAL FERMIONS INVARIANT MASS IN GEV**2              *
* -ZPRIPP- Z_min CUT FOR INITIAL STATE PAIRS                           *
* - ECUT - MIN OF FINAL FERMIONS ENERGIES IN GEV                       *
* - ACUT - MAX OF FINAL FERMIONS ACOLLINEARITY ANGLE IN DEGREE         *
*                                                                      *
*   O U T P U T:                                                       *
* - SBRN, ABRN - IMPROVED BORN  CROSS-SECTION AND ASYMMETRY            *
* - SQED, AQED - QED-CONVOLUTED CROSS-SECTION AND ASYMMETRY            *
************************************************************************
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /MASSZ / AME,AMF,AME2,AMF2
      COMMON /PCONST/ ALFAI,AL1PI,ALQE2,ALQF2,ALQEF,GMU,CSIGNB
      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
      COMMON /SVAR  / S
      COMMON /BORN0I/ SBORN0,ABORN0,SBORNI,ABORNI
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
      COMMON /INTRFS/ INTRF
      COMMON /FLAGZP/ ISRPPR,IFSPPR,IFUNAN
      COMMON /SOFTPR/ SOFTPR,SFPRFB
      COMMON /CPAIRS/ FPI2,RINF,RRR0,RRR1,PAIRDL,CORFAC
      COMMON /PRECIS/ NPREC
      COMMON /CORINT/ CORINT
      COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT
      COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     &                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM
*
*------------------------------------------------------------------
      COMMON /INDCUT/ INDC
      COMMON /IMISDC/ IMISD
*------------------------------------------------------------------
*
* From ''SUBROUTINE ZUCUTS'' :
*
*     ICUT (int/read)  = -1 no cuts to be used                         
*                      =  0 ACOL+EMIN cuts                             
*                      =  1 S_PR cut  
*
*     IRCUT (int/read) =  0 ACOL+EMIN cuts                    (original code)
*                      =  1 S_PR cut                          (original code)
*                      =  2 ACOL+EMIN cuts, no ACCEPTANCE CUT (new code) 
*                      =  3 ACOL+EMIN+ACCEPTANCE cuts         (new code)
*     SCUT  uses IRCUT =  3
*
*----------------------------------------------------------------------------
*
      EXTERNAL SSOFT,ASOFT,FACT,FACINV,SHARD,AHARD,PH2,AJMSB
*
      DATA EPSH/1D-5/,EPSS/1D-14/
*
      INDC = IRCUT
*
      S=SARG
      C=COS(THET*PI/180)
      CALL ZETCOS(S,C)
*
* BORNN should be accessed only via COSCUT
*
      CALL BORN(IFINAL,1D0,1D0,SBORN0,ABORN0,SBORNS,ABORNS)
      CALL BORN(-1    ,1D0,1D0,SBORNI,ABORNI,SBORNS,ABORNS)
      IF(INDF.EQ.10.AND.INTRF.EQ.2.AND.IMISD.EQ.1) THEN
        SBORNI=SBORNS
        ABORNI=ABORNS
      ENDIF
*
      SBRN=3D0/8*SBORN0*COPL3*CSIGNB/S
      ABRN=3D0/8*ABORN0*C2   *CSIGNB/S
      SQED=0D0
      AQED=0D0
      IF(IBORN.EQ.1) RETURN
*
CAC...now it is possible to use pair rc without the photonic ones
      IF(IPHOT2.EQ.-1) THEN
        SQED=3D0/8*SBORN0*COPL3*CSIGNB/S
        AQED=3D0/8*ABORN0*C2   *CSIGNB/S
CAC     RETURN
        GOTO 120
      ENDIF
*----------------------------------------------------------------------*
*     SOFT+BOXES (INITIAL-FINAL INTERFERENCE O(ALPHA))                 *
*----------------------------------------------------------------------*
      IF (INTERF.GE.1) THEN
        CALL BOXINT(BOXIS,BOXIA)
        CALL SFTINT(SFTIS,SFTIA)
      ELSE
        BOXIS=0D0
        BOXIA=0D0
        SFTIS=0D0
        SFTIA=0D0
      ENDIF
      IF(IMISD.EQ.0.AND.INTRF.EQ.2.AND.INDF.EQ.10) THEN
        BOXIS=0D0
        BOXIA=0D0
        SFTIS=0D0
        SFTIA=0D0
      ENDIF        
*----------------------------------------------------------------------*
*     CONVOLUTION                                                      *
*----------------------------------------------------------------------*
*     TOTAL: SOFT+VERTEX (INITIAL + FINAL) WITH EXPONENTIATION         *
*----------------------------------------------------------------------*
      EPSHAR=1D-2/NPREC
      IF(INDF.EQ.1) THEN
        EPSHAR=EPSHAR/1D2
      ELSEIF(INDF.EQ.2) THEN
        EPSHAR=EPSHAR/1D1
      ELSEIF(INDF.EQ.4.OR.INDF.EQ.5) THEN
        EPSHAR=EPSHAR/1D1
      ENDIF
      REPS1=1D-5/NPREC
      AEPS1=1D-20            
*
      DCUT=1D0-RCUT
      H=.25D0*(DCUT-EPSS)
      CALL FDSIMP
     &     (EPSS,DCUT,H,REPS1,AEPS1,SSOFT,FACT,FACINV,RR,SSFT,AIH,AIA)
*
CAC
C.... IFI EXPONENTIATION: 
      SIFI = 0D0
      IF(INTERF.EQ.2.AND.ABS(SQRT(SCOM)-91.2D0).GE.10D0) 
     & CALL SXPIFI(C,SIFI)
CAC
      IF(IPHOT2.EQ.5) THEN
        SQED=(SIFI+SSFT+EPSS**BETTAE*SBORN0*FYFS(0D0))*COPL3
      ELSE
        SQED=(SIFI+SSFT+EPSS**BETTAE*SBORN0)*COPL3*(1D0+SOFTER)
      ENDIF
*
      IF(IPHOT2.NE.5.OR.INTERF.GE.1) THEN
*----------------------------------------------------------------------*
*     TOTAL: HARD (INITAL + FINAL + INITIAL-FINAL INTERFERENCE)        *
*----------------------------------------------------------------------*
        REPS1 = 1D-4*EPSHAR
        AEPS1 = 1D-3*ABS(SSFT)*EPSHAR
        BOUNDH=1D0-EPSH
        H=.25D0*(BOUNDH-RCUT)
        CALL SIMPS(BOUNDH,RECUT1,H,REPS1,AEPS1,SHARD,RR,SHRD,AIHH,AIHA)
        SQED=SQED-SHRD
      ENDIF
*
      SQED=3D0/8*SQED*CSIGNB/S
*
 120  CONTINUE
*
* Soft Initial State pairs and Final State Soft+Hard ones
*
      IF(IFSPPR.EQ.0) THEN
        SFTPR=SBORN0*SOFTPR
      ELSEIF(IFSPPR.EQ.2.AND.ISRPPR.NE.1.AND.COPL3.NE.0D0) THEN
        SIGQA=SQED*8D0/3D0*S/CSIGNB/COPL3
        SFTPR=SBORN0*SOFTPR + SIGQA*FSRPR(S,RCUT,PCUT)
      ELSE
        SFTPR=SBORN0*( SOFTPR + FSRPR(S,RCUT,PCUT) )
      ENDIF
*
* Hard Initial pairs
*
      IF (ISRPPR.GE.1) THEN 
        REPS1=1D-5/NPREC
        XMIN=PAIRDL*(2D0-PAIRDL)
        ZMIN=RCUT
        XMAX=1D0-MAX(ZMIN,RCUT)
        H=0.25D0*(XMAX-XMIN)
        PHRD=0D0
        IF(IFLAGS(IFIPTO).NE.-1.AND.H.GT.0D0) 
     &  CALL SIMPS(XMIN,XMAX,H,REPS1,1D-20,PH2,AX,PHRD,AIHH,AIHA)  
        HRDPR=PHRD
        IF(ISRPPR.EQ.1) THEN
          SQED=SQED*(1D0+3D0/8*(SFTPR+HRDPR)*COPL3*CSIGNB/S/SBRN)
        ELSEIF(ISRPPR.EQ.2) THEN
          SQED=SQED+3D0/8D0*(SFTPR+HRDPR)*COPL3*CSIGNB/S
        ELSEIF(ISRPPR.GE.3) THEN
          BETE= 2D0*ALQE2*(TEE-1D0)
          TMIN= 0D0
          TMAX= XMAX**BETE
          H=0.25D0*(TMAX-TMIN)
          CALL SIMPS(TMIN,TMAX,H,REPS1,1D-20,AJMSB,AX,AJMSI,AIHH,AIHA)  
          IF(ISRPPR.EQ.3) SQED= SQED 
     &        +3D0/8D0*(AJMSI+(SFTPR+HRDPR)*0.7D0)*COPL3*CSIGNB/S
          IF(ISRPPR.EQ.4) SQED= SQED 
     &        +3D0/8D0*(AJMSI+(SFTPR+HRDPR)*1.0D0)*COPL3*CSIGNB/S
        ENDIF
      ENDIF
*
CAC   SQED=SQED+3D0/8*(SFTIS+BOXIS)*CSIGNB/S
      IF(IPHOT2.NE.-1) SQED=SQED+3D0/8*(SFTIS+BOXIS)*CSIGNB/S
*
      IF ((IAFB.EQ.0).AND.(ISYM.EQ.1)) RETURN
C.... PAIR RC TO AQED CAN NOT BE APPLIED SEPARATELY FROM PHOTONIC ONES:
      IF(IPHOT2.EQ.-1) RETURN
*
* Continue with forward-backward asymmetry
*
*----------------------------------------------------------------------*
*     FORW.-BACKW.: SOFT+VERTEX (INITAL + FINAL WITH EXPONENTIATION)   *
*     LLA pairs are added 14.10.99 (ABA)
*----------------------------------------------------------------------*
      REPS1 = 1D-4*EPSHAR
      CALL FDSIMP
     &     (EPSS,DCUT,H,REPS1,AEPS1,ASOFT,FACT,FACINV,RR,ASFT,AIH,AIA)

CAC
C.... IFI EXPONENTIATION: 
      AIFI = 0D0
      IF(INTERF.EQ.2.AND.ABS(SQRT(SCOM)-91.2D0).GE.10D0) 
     & CALL AXPIFI(C,AIFI)
CAC
      AQED=(AIFI+ASFT+EPSS**BETTAE*ABORN0)*C2*(1D0+SOFTER+SFPRFB)
*----------------------------------------------------------------------*
*     YFS-FOR ASYMMETRY DOES NOT EXIST                                 *
*----------------------------------------------------------------------*
*     FORW.-BACKW.: HARD (INITAL + FINAL + INITIAL-FINAL INTERFERENCE) *
*----------------------------------------------------------------------*
      REPS1 = 1D-4*EPSHAR
      AEPS1 = 1D-3*MAX(ABS(AQED),1D-5)*EPSHAR
      BOUNDH=1D0-EPSH
      H=.25D0*(BOUNDH-RCUT)
      CALL SIMPS(RECUT1,BOUNDH,H,REPS1,AEPS1,AHARD,RR,AHRD1,AIHH,AIHA)
      AQED=AQED+AHRD1+(SFTIA+BOXIA)
*
      AQED=3D0/8*AQED*CSIGNB/S
*                                                               END SCUT
      END

      SUBROUTINE 
     &       COSCUT(INDF,SARG,COSTHR,SPR,ZPRIPP,ECUT,ACUT,SIGBRN,SIGQED)
*======================================================================*
*         Modified by M. JACK on 26/06/99 07:00pm                      *
*======================================================================*
* To be checked yet                                                    *
************************************************************************
*    ROUTINE CALCULATES ANGULAR DISTRIBUTION IN THE SCATTERING ANGLE   *
*    FOR ANNIHILATION OF E+E- INTO MU+MU- AND N*GAM                    *
*                               WITH                                   *
*    S'-CUT   OR   COMBINE ACOLINEARITY + MUON ENERGIES CUTS.          *
*    (IT IS A FORMAL ANALOG OF SCUT)                                   *
*                                                                      *
* FOR -0.9999 < COS(THET) < 0.9999  (THET IS ANGLE BETWEEN E+ AND MU+):*
*                                                                      *
*   I N P U T:                                                         *
* - INDF - THE STANDART ZFITTER 'CHANNEL INDEX'                        *
* - SARG - S-INVARIANT IN GEV**2                                       *
* - COSTHR - COS OF ANGLE BETWEEN e^+ and f                            *
* - SPR  - MIN OF FINAL FERMIONS INVARIANT MASS IN GEV**2              *
* -ZPRIPP- Z_min CUT FOR INITIAL STATE PAIRS                           *
* - ECUT - MIN OF FINAL FERMIONS ENERGIES IN GEV                       *
* - ACUT - MAX OF FINAL FERMIONS ACOLLINEARITY ANGLE IN DEGREE         *
*                                                                      *
*   O U T P U T:                                                       *
* - SIGBRN, SIGQED                                                     *
* - SIGBRN - IMPROVED BORN  DIFFERENTIAL CROSS-SECTION                 *
* - SIGQED - QED-CONVOLUTED DIFFERENTIAL CROSS-SECTION                 *
************************************************************************
*
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 XVEFI,XAEFI,XVPOL,XKAPP,XKAPPC,XMZ2,XMZ2C
      COMPLEX*16 XROZ,XVEZ
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /MASSZ / AME,AMF,AME2,AMF2
      COMMON /PCONST/ ALFAI,AL1PI,ALQE2,ALQF2,ALQEF,GMU,CSIGNB
      COMMON /PSCONS/ SW2,AMZ,GAMZ
      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
      COMMON /SVAR  / S
      COMMON /BORN0I/ SBORN0,ABORN0,SBORNI,ABORNI
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
      COMMON /INTRFS/ INTRF
      COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT
      COMMON /FLAGZP/ ISRPPR,IFSPPR,IFUNAN
      COMMON /SOFTPR/ SOFTPR,SFPRFB
      COMMON /CPAIRS/ FPI2,RINF,RRR0,RRR1,PAIRDL,CORFAC
      COMMON /PRECIS/ NPREC
      COMMON /CORINT/ CORINT
      COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     &                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM
      COMMON /COUPL / VEFA,XVEFI,VEFZ,AEFA,XAEFI,AEFZ,VEEZ,XVPOL,VPOL2
      COMMON /COUPL2/ XROZ,XVEZ
      COMMON /FORCHI/ XKAPP,XKAPPC,XMZ2,XMZ2C
*
*------------------------------------------------------------------
      COMMON /INDCUT/ INDC
      COMMON /IMISDC/ IMISD
*------------------------------------------------------------------
*
* From ''SUBROUTINE ZUCUTS'' :
*
*     ICUT (int/read)  = -1 no cuts to be used                         
*                      =  0 ACOL+EMIN cuts                             
*                      =  1 S_PR cut  
*
*     IRCUT (int/read) =  0 ACOL+EMIN cuts                    (original code)
*                      =  1 S_PR cut                          (original code)
*                      =  2 ACOL+EMIN cuts, no ACCEPTANCE CUT (new code) 
*                      =  3 ACOL+EMIN+ACCEPTANCE cuts         (new code)
*     COSCUT uses IRCUT=  3
*
*----------------------------------------------------------------------------
*
      COMMON /INDFOS/ INDFOS
*
      EXTERNAL SOFT,FACT,FACINV,HARD,PH2,AJMSB
*
      DATA EPSH/1D-5/,EPSS/1D-14/
*
      INDC = IRCUT
      INDFOS=INDF
*
      S=SARG
      C=COSTHR
      CALL ZETCOS(S,C)
*
      CALL BORN(IFINAL,1D0,1D0,SBORN0,ABORN0,SBORNS,ABORNS)
      CALL BORN(-1    ,1D0,1D0,SBORNI,ABORNI,SBORNS,ABORNS)
      IF(INDF.EQ.10.AND.INTRF.EQ.2.AND.IMISD.EQ.1) THEN
        SBORNI=SBORNS
        ABORNI=ABORNS
      ENDIF
*
* Access to BORNN without ISR convolution
*
      IF(INDF.NE.-1) THEN
        SIGBRN=3D0/8*(SBORN0*COPL2+ABORN0*2D0*C)*CSIGNB/S
      ELSEIF(INDF.EQ.-1) THEN
        CALL BORNN(1D0,1D0,C,SBORNN)
        SIGBRN=SBORNN
      ENDIF
      SIGQED=0D0
      IF(IBORN.EQ.1) RETURN
*
CAC...now it is possible to use pair rc without the photonic ones
CDB...but not for EENN
      IF(IPHOT2.EQ.-1) THEN
        SIGQED=3D0/8*(SBORN0*COPL2+ABORN0*2D0*C)*CSIGNB/S
CAC     RETURN
        GOTO 130
      ENDIF
*----------------------------------------------------------------------*
*     SOFT+BOXES (INITIAL-FINAL INTERFERENCE O(ALPHA))                 *
*----------------------------------------------------------------------*
      IF (INTERF.GE.1) THEN
        CALL BOXIN (BOXI)
        CALL SOFTIN(SFTI)
      ELSE 
        BOXI=0D0
        SFTI=0D0
      ENDIF
      IF(IMISD.EQ.0.AND.INTRF.EQ.2.AND.INDF.EQ.10) THEN
        BOXI=0D0
        SFTI=0D0
      ENDIF
*----------------------------------------------------------------------*
*     CONVOLUTION                                                      *
*----------------------------------------------------------------------*
*     DIFFERENTIAL: SOFT+VERTEX (INITIAL + FINAL) WITH EXPONENTIATION  *
*     EENN IS ACCESSED HERE VIA IPHOT2.EQ.5 AND IFI=0 FOR EENN         *
*----------------------------------------------------------------------*
      EPSHAR=1D-2/NPREC
      IF(INDF.EQ.1) THEN
        EPSHAR=EPSHAR/1D2
      ELSEIF(INDF.EQ.2) THEN
        EPSHAR=EPSHAR/1D1
      ELSEIF(INDF.EQ.4.OR.INDF.EQ.5) THEN
        EPSHAR=EPSHAR/1D1
      ENDIF
      REPS1=1D-5/NPREC
      AEPS1=1D-20 
*
      DCUT=1D0-RCUT
      H=.25D0*(DCUT-EPSS)
      CALL FDSIMP
     &     (EPSS,DCUT,H,REPS1,AEPS1,SOFT,FACT,FACINV,RR,SFT,AIH,AIA)
*
CAC
C.... IFI EXPONENTIATION: 
      CIFI = 0D0
      IF(INTERF.EQ.2.AND.ABS(SQRT(SCOM)-91.2D0).GE.10D0) 
     & CALL CXPIFI(C,CIFI)
CAC
      IF(IPHOT2.EQ.5) THEN
       IF(INDF.NE.-1) THEN
        SIGQED=(CIFI+SFT+(EPSS**BETTAE*SBORN0*FYFS(0D0))*COPL2
     &                  +(EPSS**BETTAE*ABORN0*FYFS(0D0))*2D0*C)
       ELSEIF(INDF.EQ.-1) THEN
        SIGQED=SFT+EPSS**BETTAE*SBORNN*FYFS(0D0)
       ELSE
        PRINT *,'Wrong flag setting'
        STOP
       ENDIF
      ELSE
        SIGQED=(CIFI+SFT+(EPSS**BETTAE*SBORN0)*COPL2
     &                  +(EPSS**BETTAE*ABORN0)*2D0*C)*(1D0+SOFTER)
      ENDIF
      IF(INDF.EQ.-1) RETURN
*
      IF(IPHOT2.NE.5.OR.INTERF.GE.1) THEN
*----------------------------------------------------------------------*
*     DIFFERENTIAL HARD (INITAL + FINAL + INITIAL-FINAL INTERFERENCE)  *
*----------------------------------------------------------------------*
        REPS1 = 1D-4*EPSHAR
        AEPS1 = 1D-5*ABS(SFT)
        BOUNDH=1D0-EPSH
        H=.25D0*(BOUNDH-RCUT)/20  !  to remember why?
        CALL SIMPS(BOUNDH,RECUT1,H,REPS1,AEPS1,HARD,RR,HRD,AIHH,AIHA)
        SIGQED=SIGQED-HRD
        SISQED=SIGQED-HRD
      ENDIF
*
      SIGQED=3D0/8*SIGQED*CSIGNB/S
*
 130  CONTINUE
*
* Soft Initial State and Final State pairs  
*
      IF(IFSPPR.EQ.0) THEN
        SFTPR=SBORN0*SOFTPR
      ELSEIF(IFSPPR.EQ.2.AND.ISRPPR.NE.1.AND.COPL2.NE.0D0) THEN
        SIGQAA=SIGQED*8D0/3D0*S/CSIGNB/COPL2
        SFTPR=SBORN0*SOFTPR + SIGQAA*FSRPR(S,RCUT,PCUT)
      ELSE
        SFTPR=SBORN0*( SOFTPR + FSRPR(S,RCUT,PCUT) )
      ENDIF
*
* Hard Initial pairs
*
      IF (ISRPPR.GE.1) THEN 
        REPS1=1D-5/NPREC
        XMIN=PAIRDL*(2D0-PAIRDL)
        ZMIN=RCUT
        XMAX=1D0-MAX(ZMIN,RCUT)
        H=0.25D0*(XMAX-XMIN)
        PHRD=0D0
        IF(IFLAGS(IFIPTO).NE.-1.AND.H.GT.0D0) 
     &  CALL SIMPS(XMIN,XMAX,H,REPS1,1D-20,PH2,AX,PHRD,AIHH,AIHA)  
        HRDPR=PHRD
        IF(ISRPPR.EQ.1) THEN
          SIGQED=SIGQED*(1D0+(SFTPR+HRDPR)/SBORN0)
        ELSEIF(ISRPPR.EQ.2) THEN
          SIGQED=SIGQED+3D0/8D0*(SFTPR+HRDPR)*COPL2*CSIGNB/S
        ELSEIF(ISRPPR.GE.3) THEN
          BETE= 2D0*ALQE2*(TEE-1D0)
          TMIN= 0D0
          TMAX= XMAX**BETE
          H=0.25D0*(TMAX-TMIN)
          CALL SIMPS(TMIN,TMAX,H,REPS1,1D-20,AJMSB,AX,AJMSI,AIHH,AIHA)  
          IF(ISRPPR.EQ.3) SIGQED=SIGQED
     &          +3D0/8D0*(AJMSI+(SFTPR+HRDPR)*0.7D0)*COPL2*CSIGNB/S
          IF(ISRPPR.EQ.4) SIGQED=SIGQED
     &          +3D0/8D0*(AJMSI+(SFTPR+HRDPR)*1.0D0)*COPL2*CSIGNB/S
        ENDIF
      ENDIF
*
CAC
      IF(IPHOT2.NE.-1) SIGQED=SIGQED+3D0/8*(SFTI+BOXI)*CSIGNB/S      
CAC
*                                                             END COSCUT
      END
 
      SUBROUTINE 
     &         SFAST(INDF,SARG,SPR,ZPRIPP,ECUT,ACUT,SBRN,ABRN,SQED,AQED)
*======================================================================*
* Modified by M. Jack on 03/03/99 05:00pm: comment added               *
************************************************************************
*  ROUTINE CALCULATES TOTAL X-SECTION AND ASYMMETRY                    *
*                                                                      *
*   I N P U T:                                                         *
* - INDF - THE STANDART ZFITTER 'CHANNEL INDEX'                        *
* - SARG - S-INVARIANT IN GEV**2                                       *
* - SPR  - MIN OF FINAL FERMIONS INVARIANT MASS IN GEV**2              *
*                                                                      *
*   O U T P U T:                                                       *
* - SBRN, ABRN - IMPROVED BORN  CROSS-SECTION AND ASYMMETRY            *
* - SQED, AQED - QED-CONVOLUTED CROSS-SECTION AND ASYMMETRY            *
************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /MASSZ / AME,AMF,AME2,AMF2
      COMMON /PCONST/ ALFAI,AL1PI,ALQE2,ALQF2,ALQEF,GMU,CSIGNB
      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
      COMMON /SVAR  / S
      COMMON /BORN0I/ SBORN0,ABORN0,SBORNI,ABORNI
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
      COMMON /INTRFS/ INTRF
      COMMON /INDFIT/ IND,INDFC
      COMMON /FLAGZP/ ISRPPR,IFSPPR,IFUNAN
      COMMON /SOFTPR/ SOFTPR,SFPRFB
      COMMON /CPAIRS/ FPI2,RINF,RRR0,RRR1,PAIRDL,CORFAC
      COMMON /PRECIS/ NPREC
      COMMON /CORINT/ CORINT
      COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT
      COMMON /RMINV / EPSHR,EPSSR
*
*------------------------------------------------------------------
      COMMON /INDCUT/ INDC
      COMMON /IMISDC/ IMISD
*------------------------------------------------------------------
*
* From ''SUBROUTINE ZUCUTS'' :
*
*     ICUT (int/read)  = -1 no cuts to be used                         
*                      =  0 ACOL+EMIN cuts                             
*                      =  1 S_PR cut  
*
*     IRCUT (int/read) =  0 ACOL+EMIN cuts                    (original code)
*                      =  1 S_PR cut                          (original code)
*                      =  2 ACOL+EMIN cuts, no ACCEPTANCE CUT (new code) 
*                      =  3 ACOL+EMIN+ACCEPTANCE cuts         (new code)
*     SFAST uses IRCUT =  2
*
*----------------------------------------------------------------------------
*
      EXTERNAL FCROS,FASYM,FACT,FACINV,PH2,AJMSB
*
      DATA EPSH/1D-5/,EPSS/1D-14/
*
      INDC = IRCUT
*
      EPSHR=EPSH/NPREC
      EPSSR=EPSS/NPREC
*
      S=SARG
      CALL ZETCOS(S,1D0)
*
      CALL BORN(IFINAL,1D0,1D0,SBORN0,ABORN0,SBORNS,ABORNS)
      CALL BORN(-1    ,1D0,1D0,SBORNI,ABORNI,SBORNS,ABORNS)
      IF(INDF.EQ.10.AND.INTRF.EQ.2.AND.IMISD.EQ.1) THEN
         SBORNI=SBORNS
         ABORNI=ABORNS
      ENDIF
*     
      SBRN=0.5D0*SBORN0*CSIGNB/S
      ABRN=3D0/8*ABORN0*CSIGNB/S
      SQED=0D0
      AQED=0D0
      IF(IBORN.EQ.1) RETURN
*
CAC...now it is possible to use pair rc without the photonic ones
      IF(IPHOT2.EQ.-1) THEN
        SQED=0.5D0*SBORN0*CSIGNB/S
        AQED=3D0/8*ABORN0*CSIGNB/S
CAC     RETURN
        GOTO 140
      ENDIF
*----------------------------------------------------------------------*
*     SOFT+BOXES (INITIAL-FINAL INTERFERENCE O(ALPHA))                 *
*----------------------------------------------------------------------*
      IF (INTERF.GE.1) THEN
        CALL BOXINT(BOXIS,BOXIA)
        CALL SFTINT(SFTIS,SFTIA)
      ELSE
        BOXIS=0D0
        BOXIA=0D0
        SFTIS=0D0
        SFTIA=0D0
      ENDIF
      IF(IMISD.EQ.0.AND.INTRF.EQ.2.AND.INDF.EQ.10) THEN
        BOXIS=0D0
        BOXIA=0D0
        SFTIS=0D0
        SFTIA=0D0
      ENDIF        
*----------------------------------------------------------------------*
*     CONVOLUTION                                                      *
*----------------------------------------------------------------------*
*     TOTAL: SOFT+VERTEX (INITIAL + FINAL) WITH EXPONENTIATION         *
*----------------------------------------------------------------------*
      REPS1=1D-5/NPREC
      AEPS1=1D-20
*
      DCUT1=1D0-RCUT
      IF(INDC.EQ.2) DCUT1=1D0-RECUT1
      H=.25D0*(DCUT1-EPSS)
      CALL FDSIMP
     &     (EPSS,DCUT1,H,REPS1,AEPS1,FCROS,FACT,FACINV,RR,CROS,AIH,AIA)
*
CAC
C.... IFI EXPONENTIATION: 
      SIFI = 0D0
      IF(INTERF.EQ.2.AND.ABS(SQRT(SCOM)-91.2D0).GE.10D0) 
     & CALL SXPIFI(1D0,SIFI)
      SIFI = SIFI*4D0/3D0
CAC

      IF(IPHOT2.EQ.5) THEN
        SQED=SIFI+CROS+EPSS**BETTAE*SBORN0*4D0/3*FYFS(0D0)
      ELSE
        SQED=SIFI+CROS+EPSS**BETTAE*SBORN0*4D0/3*(1D0+SOFTER)
      ENDIF
*----------------------------------------------------------------------*
*     HARD (INITIAL + FINAL) IN SFAST DOES NOT EXIST                   *
*----------------------------------------------------------------------*
*
      SQED=3D0/8*SQED*CSIGNB/S
*
 140  CONTINUE
*
* Soft Initial State and Final State pairs  
*
      IF(IFSPPR.EQ.0) THEN
        SFTPR=SBORN0*SOFTPR
      ELSEIF(IFSPPR.EQ.2.AND.ISRPPR.NE.1) THEN
        SIGQA=SQED*2D0*S/CSIGNB
        SFTPR=SBORN0*SOFTPR + SIGQA*FSRPR(S,RCUT,PCUT)
      ELSE
        SFTPR=SBORN0*( SOFTPR + FSRPR(S,RCUT,PCUT) )
      ENDIF

*
* Hard Initial pairs
*
      IF (ISRPPR.GE.1) THEN
        REPS1=1D-5/NPREC
        XMIN=PAIRDL*(2.-PAIRDL)
        ZMIN=RCUT
        XMAX=1D0-MAX(ZMIN,RCUT)
        H=0.25D0*(XMAX-XMIN)
        PHRD=0D0
        IF(IFLAGS(IFIPTO).NE.-1.AND.H.GT.0D0) 
     &  CALL SIMPS(XMIN,XMAX,H,REPS1,1D-20,PH2,AX,PHRD,AIHH,AIHA)  
        HRDPR=PHRD
        IF (ISRPPR.EQ.1) THEN
          SQED=SQED*(1D0+.5D0*(SFTPR+HRDPR)*CSIGNB/S/SBRN)
        ELSEIF(ISRPPR.EQ.2) THEN
          SQED=SQED+.5D0*(SFTPR+HRDPR)*CSIGNB/S
        ELSEIF(ISRPPR.GE.3) THEN
          BETE= 2D0*ALQE2*(TEE-1D0)
          TMIN= 0D0
          TMAX= XMAX**BETE
          H=0.25D0*(TMAX-TMIN)
          CALL SIMPS(TMIN,TMAX,H,REPS1,1D-20,AJMSB,AX,AJMSI,AIHH,AIHA)  
          IF(ISRPPR.EQ.3) SQED= SQED 
     &        + 0.5D0*(AJMSI+(SFTPR+HRDPR)*0.7D0)*CSIGNB/S
          IF(ISRPPR.EQ.4) SQED= SQED 
     &        + 0.5D0*(AJMSI+(SFTPR+HRDPR)*1.0D0)*CSIGNB/S
        ENDIF
      ENDIF
*
      IF(IPHOT2.NE.-1) SQED=SQED+3D0/8*(SFTIS+BOXIS)*CSIGNB/S
*
      IF ((IAFB.EQ.0).AND.(ISYM.EQ.1)) RETURN
CAC...PAIR RC TO AQED CAN NOT BE APPLIED SEPARATELY FROM PHOTONIC ONES
      IF(IPHOT2.EQ.-1) RETURN
*
* Continue with forward-backward asymmetry
*
*----------------------------------------------------------------------*
*     FORW.-BACKW.: SOFT+VERTEX (INITAL + FINAL WITH EXPONENTIATION)   *
*----------------------------------------------------------------------*
      CALL FDSIMP
     &     (EPSS,DCUT1,H,REPS1,AEPS1,FASYM,FACT,FACINV,RR,ASYM,AIH,AIA)
CAC
C.... IFI EXPONENTIATION:
      AIFI = 0D0
      IF(INTERF.EQ.2.AND.ABS(SQRT(SCOM)-91.2D0).GE.10D0) 
     & CALL AXPIFI(1D0,AIFI)
CAC
      AQED=AIFI+ASYM+EPSS**BETTAE*ABORN0*(1D0+SOFTER+SFPRFB)
*----------------------------------------------------------------------*
*     HARD (INITIAL + FINAL) IN SFAST DOES NOT EXIST                   *
*----------------------------------------------------------------------*
*     FORW.-BACKW.: HARD (INITIAL-FINAL INTERFERENCE O(ALPHA))         *
*----------------------------------------------------------------------*
      AQED=AQED+(SFTIA+BOXIA)
*
      AQED=3D0/8*AQED*CSIGNB/S
*                                                              END SFAST
      END

      SUBROUTINE EWINIT(INTRF,INDF,AMZC,GZC,SW2C)
*     ========== ================================
*----------------------------------------------------------------------*
*     S-INDEPENDENT INITIALIZATION ROUTINE FOR ZCUTCOS AND ZCUT        *
*----------------------------------------------------------------------*
      IMPLICIT COMPLEX*16(X)
      IMPLICIT REAL*8(A-H,O-W,Z)
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
      COMMON /FRINIT/ NPAR(30),ZPAR(31)
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /CHARGZ/ QE,QF,QEM,QFM,QEF,QEFM,QE2,QF2
      COMMON /MASSZ / AME,AMF,AME2,AMF2
      COMMON /PCONST/ ALFAI,AL1PI,ALQE2,ALQF2,ALQEF,GMU,CSIGNB
      COMMON /PSCONS/ SW2,AMZ,GAMZ
      COMMON /POLAR / ALAMP,ALAMM,HELP,HELM
      COMMON /FORCHI/ XKAPP,XKAPPC,XMZ2,XMZ2C
      COMMON /KAPPAC/ AKAPPA
*
* NCONST
      PI=ATAN(1.D0)*4.D0
      F1=PI**2/6
      AL2=LOG(2.D0)
      ZET3=1.202 056 903 130 49D0
* CHARGES
      QE= ZPAR(1)
      QF= ZPAR(2)
* CHARGZ
      QEM=ABS(QE)
      QFM=ABS(QF)
      QEF=QE*QF
      QEFM=ABS(QEF)
      QE2=QE**2
      QF2=QF**2
* MASSES
      AME=.510 999 06 D-3
      AMF=ZPAR(4)
      AME2=AME**2
      AMF2=AMF**2
      IF(AMF.LE.0D0) AMF=1D-10
* MASSZ
      AME2=AME**2
      AMF2=AMF**2
* PCONST
          IF(IFLAGS(IFGFER).EQ.0) THEN
       GMU=1.166388D-5
      ELSEIF(IFLAGS(IFGFER).EQ.1) THEN
       GMU=1.16639D-5
      ELSEIF(IFLAGS(IFGFER).EQ.2) THEN
       GMU=1.16637D-5
      ENDIF
*
      ALFAI=137.035 989 5 D0
      AL1PI=1.D0/(PI*ALFAI)
      ALQE2=QE2*AL1PI
      ALQF2=QF2*AL1PI
      ALQEF=QEF*AL1PI
* CONHC - NB CONVERSION FACTOR
      CONHC=.389 379 66 D+6
      CSIGNB=4D0/3D0*PI/ALFAI**2*CONHC
*
* ADDITIONAL INITIALIZATION FOR ZUXSEC
      COLOR= 1.D0
      IF(INDF.GE.4) COLOR=3D0
      ZPAR(3) = COLOR
* END OF ADDITIONAL INITIALIZATION
* INITIAL POLARISATIONS
      ALAMP= ZPAR(22)
      ALAMM= ZPAR(23)
* FINAL HELICITIES
      HELP= ZPAR(24)
      HELM= ZPAR(25)
* ELECTROWEAK PARAMETERS
      ZPAR(5)=AMZC
      AMZ  = ZPAR(5)
      GAMZ = GZC
      SW2  = SW2C
*
*    FOR CHI
*    TAKING INTO ACCOUNT S-DEPENDENCE IN Z-WIDTH
*
      AKAPPA=GMU*AMZ**2/(SQRT(2D0)*8D0*PI)*ALFAI
      IF(NPAR(5).EQ.0) THEN
        XMZ2 =DCMPLX(AMZ**2,-AMZ*GAMZ)
        XMZ2C=DCMPLX(AMZ**2, AMZ*GAMZ)
        GAMMA=GAMZ/AMZ
        XKAPP =AKAPPA
        XKAPPC=AKAPPA
          ELSE
        XMZ2 =DCMPLX(AMZ**2,-AMZ*GAMZ)
        XMZ2C=DCMPLX(AMZ**2, AMZ*GAMZ)
        GAMMA=GAMZ/AMZ
        GAM2=1.D0+GAMMA**2
        XMZ2 =XMZ2 /GAM2
        XMZ2C=XMZ2C/GAM2
        XKAPP =AKAPPA/(1.D0+(0.D0, 1.D0)*GAMMA)
        XKAPPC=AKAPPA/(1.D0+(0.D0,-1.D0)*GAMMA)
      ENDIF
*
      END
 
      SUBROUTINE ZETCOS(S,CIN)
*=============================
      IMPLICIT REAL*8(A-H,O-Z)
*
      COMMON /MASSZ / AME,AMF,AME2,AMF2
      COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     &                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM
* 
      C=CIN
      C2=C*C
      C3=C*C2
      COPL2=1D0+C2
      COPL3=1D0/3*C3+C
*
      CP=(1D0+C)/2D0
      CP=MAX(CP,AME2/S)
      CM=(1D0-C)/2D0
      CM=MAX(CM,AME2/S)
      CP2=CP*CP
      CM2=CM*CM
      CP3=CP2*CP
      CM3=CM2*CM
*
      CPM=CP*CM
      CPM2=CPM*CPM
      CPM3=CPM*CPM2
*
      ALCP=LOG(CP)
      ALCM=LOG(CM)
      ALPL=ALCP+ALCM
      ALMI=ALCP-ALCM
*
      DLCP=DDILOG(CP)
      DLCM=DDILOG(CM)
*
      END
 
      SUBROUTINE SETFUN(SARG)
*     ========== ============
      IMPLICIT COMPLEX*16(X)
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
*
* flags
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /PCONST/ ALFAI,AL1PI,ALQE2,ALQF2,ALQEF,GMU,CSIGNB
      COMMON /MASSZ / AME,AMF,AME2,AMF2
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
      COMMON /FLAGZP/ ISRPPR,IFSPPR,IFUNAN
      COMMON /SOFTPR/ SOFTPR,SFPRFB
      COMMON /CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
      COMMON /CPAIRS/ FPI2,RINF,RRR0,RRR1,PAIRDL,CORFAC
*
      S=SARG
      SCOM=S
      TE=LOG(S/AME2)-1
      TEE=LOG(S/(AML(2)**2))
      TMU=LOG(S/(AML(4)**2))
      TTU=LOG(S/(AML(6)**2))
      FPI2=4*(.1396D0)**2
      TPI=LOG(S/FPI2)
      TF=LOG(S/AMF2)-1
      BETTAE=TE*2*ALQE2
      BETTAF=TF*2*ALQF2
      SOFTER=3.D0/4*BETTAE+ALQE2*(F1*2-1D0/2)
      SOFTFR=3.D0/4*BETTAF+ALQF2*(F1*2-1D0/2)
      SOFTPR=0D0
       TL=TE+1
*SECOND ORDER INITIAL (VIRTUAL+SOFT) CORRECTIONS
      IF(IPHOT2.GE.0) THEN
       SOFTER=SOFTER+ALQE2**2*TL**2*(9.D0/8-2*F1)
      ENDIF
      IF(IPHOT2.GE.1) THEN
       SOFTER=SOFTER+ALQE2**2*(
     +       +TL*(-45D0/16+11D0/2*F1+3D0*ZET3) )
      ENDIF
      IF(IPHOT2.GE.2) THEN
       SOFTER=SOFTER+ALQE2**2*(
     +       -6D0/5*F1**2-9D0/2*ZET3-6*F1*AL2+3D0/8*F1+19D0/4 )
      ENDIF
*
* V+S are in 1\to 1 agreement with YR(1989), Eq.(3.17)
*
      IF(IPHOT2.EQ.3) THEN
       PSI21=-2D0*ZET3 ! to be confirmed
       SOFTER=SOFTER+ALQE2**3*TE**3*(9D0/16-3D0*F1-4D0/3*PSI21) 
      ENDIF
*
* YFS is treated separately 
*
      IF(IPHOT2.EQ.5) THEN
       SOFTER=SOFTER
      ENDIF
*
* PAIRS
*
      ISRPPR=IFLAGS(IFISPP)
      IFSPPR=IFLAGS(IFFSPP)
*
* FUNANG by Arbuzov
*
      IFUNAN=IFLAGS(IFFUNA) 
*
      IF(ISRPPR.EQ.-1) THEN
        PAIRDL=.0055D0
        CORFAC=2.6D0
      ELSEIF(ISRPPR.NE.-1) THEN
CAC     PAIRDL=1D-9
        PAIRDL=1D-6
        CORFAC=1D0
      ENDIF
*
      IF(ISRPPR.NE.0.AND.IFLAGS(IFIPTO).NE.-1) THEN
       RINF=4.0D0
       RRR0=-8.31D0
       RRR1= 13.1D0
       DELL=2D0*DLOG(2D0*PAIRDL)
* Thebook
CAC    SOFTPR=+ALQE2**2*(
C    & +(1D0/2*TEE**2-5D0/3*TEE+28D0/9-2*F1)*(1D0/3*DELL+1D0/2)
C    & +(TEE-5D0/3)*(1D0/6*DELL**2-7D0/12)+1D0/18*DELL**3+2D0/3*ZET3
C    & +5D0/8+2D0/3*ZET3+8D0/9*F1-1241D0/648)
C    &        +ALQE2**2*(
C    & +(1D0/2*TMU**2-5D0/3*TMU+28D0/9-2*F1)*(1D0/3*DELL+1D0/2)
C    & +(TMU-5D0/3)*(1D0/6*DELL**2-7D0/12)+1D0/18*DELL**3+2D0/3*ZET3
C    & +5D0/8)*CORFAC
C    &        +ALQE2**2*(
C    & +(1D0/2*TTU**2-5D0/3*TTU+28D0/9-2*F1)*(1D0/3*DELL+1D0/2)
C    & +(TTU-5D0/3)*(1D0/6*DELL**2-7D0/12)+1D0/18*DELL**3+2D0/3*ZET3
C    & +5D0/8)
C    &        +ALQE2**2*(
C    & +(RINF*(1D0/2*TPI**2-F1)+RRR0*TPI+RRR1)*(1D0/3*DELL+1D0/2)
C    & +(RINF*TPI+RRR0)*(1D0/6*DELL**2-7D0/12)
CAC  & +RINF*(1D0/18*DELL**3+2D0/3*ZET3+5D0/8) )
       ALEADE = (1D0/2.*TEE**2)*(1D0/3.*DELL+1D0/2.)
       ANLLAE =
     & +(-5D0/3.*TEE+28D0/9.-2.*F1)*(1D0/3.*DELL+1D0/2.)
     & +(TEE-5D0/3)*(1D0/6*DELL**2-7D0/12)+1D0/18*DELL**3+2D0/3*ZET3
     & +5D0/8. + 2D0/3.*ZET3 + 8D0/9.*F1 + 1241D0/648.
       ALEADM = (1D0/2.*TMU**2)*(1D0/3.*DELL+1D0/2.)*CORFAC
       ANLLAM =
     & ((-5D0/3.*TMU+28D0/9.-2.*F1)*(1D0/3.*DELL+1D0/2.)
     & +(TMU-5D0/3)*(1D0/6*DELL**2-7D0/12)+1D0/18*DELL**3+2D0/3*ZET3
     & +5D0/8. )*CORFAC
       ALEADT = (1D0/2.*TTU**2)*(1D0/3.*DELL+1D0/2.)
       ANLLAT =
     & +(-5D0/3*TTU+28D0/9-2*F1)*(1D0/3*DELL+1D0/2)
     & +(TTU-5D0/3)*(1D0/6*DELL**2-7D0/12)+1D0/18*DELL**3+2D0/3*ZET3
     & +5D0/8
       ALEADH = (RINF*(1D0/2.*TPI**2))*(1D0/3.*DELL+1D0/2.)
       ANLLAH = (RINF*(-F1)+RRR0*TPI+RRR1)*(1D0/3.*DELL+1D0/2.)
     & +(RINF*TPI+RRR0)*(1D0/6.*DELL**2-7D0/12.)
     & +RINF*(1D0/18.*DELL**3+2D0/3.*ZET3+5D0/8.) 

       IF(IFLAGS(IFISPP).LE.2) THEN
        IF(IFLAGS(IFIPFC).EQ.1) ALEAD = ALEADE
        IF(IFLAGS(IFIPFC).EQ.2) ALEAD = ALEADM
        IF(IFLAGS(IFIPFC).EQ.3) ALEAD = ALEADT
        IF(IFLAGS(IFIPFC).EQ.4) ALEAD = ALEADH
        IF(IFLAGS(IFIPFC).EQ.5) ALEAD = ALEADE+ALEADM+ALEADT+ALEADH
        IF(IFLAGS(IFIPFC).EQ.6) ALEAD = ALEADE+ALEADM+ALEADT
        IF(IFLAGS(IFIPFC).EQ.1) ANLLA = ANLLAE
        IF(IFLAGS(IFIPFC).EQ.2) ANLLA = ANLLAM
        IF(IFLAGS(IFIPFC).EQ.3) ANLLA = ANLLAT
        IF(IFLAGS(IFIPFC).EQ.4) ANLLA = ANLLAH
        IF(IFLAGS(IFIPFC).EQ.5) ANLLA = ANLLAE+ANLLAM+ANLLAT+ANLLAH
        IF(IFLAGS(IFIPFC).EQ.6) ANLLA = ANLLAE+ANLLAM+ANLLAT
       ELSEIF(IFLAGS(IFISPP).GE.3) THEN
        ALEAD = 0D0
        IF(IFLAGS(IFIPFC).EQ.4.OR.IFLAGS(IFIPFC).EQ.5) ALEAD = ALEADH
        ANLLA = 0D0
        IF(IFLAGS(IFIPFC).EQ.4.OR.IFLAGS(IFIPFC).EQ.5) ANLLA = ANLLAH
       ENDIF
       SOFTPR= ALQE2**2*( ALEAD + ANLLA )
      ELSEIF(ISRPPR.NE.0.AND.IFLAGS(IFIPTO).EQ.-1) THEN
       RINF= 4.00D0
       RRR0=-8.31D0
       RRR1= 13.1D0
       RRR2=-15.6D0
       AVIRTE= 2D0*( - TEE**3/36D0 + 19D0/72D0*TEE**2
     & + (F1/3D0-265D0/216D0)*TEE - 11D0*F1/6D0 + 383D0/108D0 )
       AVIRTM= 2D0*( - TMU**3/36D0 + 19D0/72D0*TMU**2
     & + (F1/3D0-265D0/216D0)*TMU - ZET3/3D0 - 19D0*F1/18D0
     & + 3355D0/1296D0 )
       AVIRTT= 2D0*( - TTU**3/36D0 + 19D0/72D0*TTU**2
     & + (F1/3D0-265D0/216D0)*TTU - ZET3/3D0 - 19D0*F1/18D0
     & + 3355D0/1296D0 )
       AVIRTH= 2D0*( 
     &   RINF*( - TPI**3/36D0 + 1D0/8D0*TPI**2 + (F1/6D0-7D0/24D0)*TPI
     & - F1/4D0 + 15D0/16D0 )
     & + RRR0*( - TPI**2/12D0 + TPI/4D0 + F1/6D0 - 7D0/24D0 )
     & + RRR1*( - TPI/6D0 + 0.25D0 ) - RRR2/6D0 )
       IF(IFLAGS(IFIPFC).EQ.1) AVIRT = AVIRTE
       IF(IFLAGS(IFIPFC).EQ.2) AVIRT = AVIRTM
       IF(IFLAGS(IFIPFC).EQ.3) AVIRT = AVIRTT
       IF(IFLAGS(IFIPFC).EQ.4) AVIRT = AVIRTH
       IF(IFLAGS(IFIPFC).EQ.5) AVIRT = AVIRTE+AVIRTM+AVIRTT+AVIRTH
       IF(IFLAGS(IFIPFC).EQ.6) AVIRT = AVIRTE+AVIRTM+AVIRTT
       SOFTPR= ALQE2**2*AVIRT
       GOTO 1
      ENDIF
C.... NEW THIRD ORDER CONTRIBUTION
      IF(ISRPPR.EQ.2.AND.IFLAGS(IFIPTO).GE.1) THEN
       IF(IFLAGS(IFIPFC).EQ.1) FCON3 = PPH0(TEE,DELL)
       IF(IFLAGS(IFIPFC).EQ.2) FCON3 = PPH0(TMU,DELL)
       IF(IFLAGS(IFIPFC).EQ.3) FCON3 = PPH0(TTU,DELL)
       IF(IFLAGS(IFIPFC).EQ.4) FCON3 = PHA0(TPI,DELL)
       IF(IFLAGS(IFIPFC).EQ.5) FCON3 = PPH0(TEE,DELL) 
     &                               + PPH0(TMU,DELL)
     &                               + PHA0(TPI,DELL)
     &                               + PPH0(TTU,DELL)
       IF(IFLAGS(IFIPFC).EQ.6) FCON3 = PPH0(TEE,DELL) 
     &                               + PPH0(TMU,DELL)
     &                               + PPH0(TTU,DELL)

       SOFTPR= SOFTPR + ALQE2**3*(TEE-1D0)*FCON3
       IF(IFLAGS(IFIPTO).GE.2) THEN
        IF(IFLAGS(IFIPFC).EQ.1.OR.IFLAGS(IFIPFC).EQ.5) SOFTPR= SOFTPR
     &       + ALQE2**3*TEE**3*(DELL+1.5D0)/27D0
       ENDIF
      ENDIF
C.... NEW! THE FOURTH ORDER
      IF(ISRPPR.EQ.2.AND.IFLAGS(IFIPTO).EQ.3) THEN 
       IF(IFLAGS(IFIPFC).EQ.1.OR.IFLAGS(IFIPFC).EQ.5) THEN
       P30   = 
     &      48D0*( ZET3/3D0 - F1/2D0*(DELL/2D0+3D0/4D0)
     &    + 1D0/6D0*(DELL/2D0+3D0/4D0)**3 )
       P20   = (DELL+1.5D0)**2 - 4D0*F1  
       P10   = DELL + 1.5D0
       FCOE4 = P30/12D0 + P20*22D0/(16D0*27D0) + P10/(4D0*27D0)
       SOFTPR= SOFTPR + ALQE2**4*(TEE-1D0)**4*FCOE4
       ENDIF
      ENDIF

C.... THE OPTION BELOW IS A MODIFICATION OF JMS IN ORDER
C     TO INCLUDE THIRD ORDER HADRONIC PAIRS
      IF(IFLAGS(IFISPP).EQ.4) THEN
       FCON3 = 0D0
       IF(IFLAGS(IFIPFC).EQ.4.OR.IFLAGS(IFIPFC).EQ.5)
     & FCON3 = PHA0(TPI,DELL)
       SOFTPR= SOFTPR + ALQE2**3*(TEE-1D0)*FCON3
      ENDIF

C.... SOFT PAIRS FOR A_FB AND ANGULAR DISTRIBUTIONS
      SFPRFB = 0D0
      IF(ISRPPR.GE.2.AND.IFLAGS(IFFBHO).NE.0) THEN
C....  ADD SECOND ORDER LLA PAIRS FOR CORRECTIONS TO A_FB
       IF(IFLAGS(IFIPFC).EQ.1) SFPRFB = (ALQE2/2D0*TEE)**2
       IF(IFLAGS(IFIPFC).EQ.2) SFPRFB = (ALQE2/2D0*TMU)**2
       IF(IFLAGS(IFIPFC).EQ.3) SFPRFB = (ALQE2/2D0*TTU)**2
       IF(IFLAGS(IFIPFC).EQ.4) SFPRFB = (ALQE2/2D0*TPI)**2*RINF
       IF(IFLAGS(IFIPFC).EQ.5) SFPRFB = (ALQE2/2D0*TEE)**2
     &                                + (ALQE2/2D0*TMU)**2
     &                                + (ALQE2/2D0*TTU)**2
     &                                + (ALQE2/2D0*TPI)**2*RINF
       IF(IFLAGS(IFIPFC).EQ.6) SFPRFB = (ALQE2/2D0*TEE)**2
     &                                + (ALQE2/2D0*TMU)**2
     &                                + (ALQE2/2D0*TTU)**2
      ENDIF

 1    CONTINUE
*
      IF(ISRPPR.EQ.-1) SOFTER=SOFTER+SOFTPR
* 
      END
  
      SUBROUTINE EWCOUP(INTRF,INDF,S)
*     ========== ====================
*----------------------------------------------------------------------------
*   S-DEPENDENT INITIALIZATION OF ELECTROWEAK FORMFACTORS VIA CALL TO ROKANC
*----------------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
* flags
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
      COMMON /CDZRKZ/ARROFZ(0:10),ARKAFZ(0:10),ARVEFZ(0:10),ARSEFZ(0:10)
     &              ,AROTFZ(0:10),AIROFZ(0:10),AIKAFZ(0:10),AIVEFZ(0:10)
      COMMON /CDZAUX/PARTZA(0:10),PARTZI(0:10),RENFAC(0:10),SRENFC(0:10)
      COMMON /COUPL / VEFA,XVEFI,VEFZ,AEFA,XAEFI,AEFZ,VEEZ,XVPOL,VPOL2
      COMMON /COUPL2/ XROZ,XVEZ
      COMMON /COUPL0/ VEFZ0,AVEFZ0(6)
      COMMON /FRINIT/ NPAR(30),ZPAR(31)
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /PCONST/ ALFAI,AL1PI,ALQE2,ALQF2,ALQEF,GMU,CSIGNB
      COMMON /MASSZ / AME,AMF,AME2,AMF2
      COMMON /PSCONS/ SW2,AMZ,GAMZ
      COMMON /CHARGZ/ QE,QF,QEM,QFM,QEF,QEFM,QE2,QF2
      COMMON /XFORMZ/ XFZ(4),XFZT,XFZTA
      COMMON /INDFIT/ IND,INDFC
      COMMON /EWFORM/ XALLCH(5,4),XFOTF
      COMMON /COCOST/ COCOST
      COMMON /CDAL5H/ DAL5H
      COMMON /POLAR / ALAMP,ALAMM,HELP,HELM
      COMMON /POLHEL/ COMB1,COMB2,HOMB1,HOMB2
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
      COMMON /FORCHI/ XKAPP,XKAPPC,XMZ2,XMZ2C
      COMMON /GAMFIT/ GAMLE,GAMFI
      COMMON /ASSFIT/ ROEC,GVEC,GAEC,ROFC,GVFC,GAFC
      COMMON /AS2FIT/ ROL2,GVL2,GAL2
      COMMON /AFBFIT/ FOUR,VAE2,VAF2
      COMMON /ZFCHMS/ ALLCH(0:11),ALLMS(0:11)
      COMMON /HADRON/ XXVEFI(6),XXAEFI(6),AVEFA(6),AAEFA(6),
     &                AVEEZ(6),AVEFZ(6),AAEFZ(6)
      COMMON /CORINT/ CORINT
      COMMON /CALQED/ ALQEDZ,ALQEDS
      COMMON /CDZRLR/ ROKL(4),ROKR(4),AROK(4)
      COMMON /KAPPAC/ AKAPPA
      COMMON /CZAKCO/ CZAKFF
      COMMON /CDZRUN/ CMQRUN(8)
      COMMON /IMISDC/ IMISD
*
* parameters
*
      COMMON /ZUPARS/QDF,QCDCOR(0:14),ALPHST,SIN2TW,S2TEFF(0:11),
     & WIDTHS(0:11)
*
      DIMENSION XFZu(4),XFZd(4)
      DIMENSION INDARR(0:11)
      DIMENSION INDQUA(4:9),INDQ(10)
      DIMENSION ARCZAK(0:10)
      DIMENSION QCDCCR(0:14)
*
* flavour-dependent CZAKUE-corrections
*
      DATA ARCZAK/
     &4*0D0,-.376D-3,-.418D-3,-.376D-3,-.418D-3,.0D-3,-.106D-3,-.336D-3/
*
* end of CZAKUE-corrections
*                                                                            
      DATA INDQ /0,0,0,0,1,2,3,4,5,6/
      DATA INDQUA/3,4,3,4,3,5/
      DATA FAA/1D0/,FZA/1D0/,FZZ/1D0/
*
      DATA INDARR/1,2,2,2,3,4,3,4,3,5,1,2/
*
      IMISD=IFLAGS(IFMISD)
      INDFC=INDF
      IHVP =NPAR(2)
      IQCD =NPAR(3)
      IALEM=NPAR(20)
      IALE2=NPAR(21)
*
* Z, Z+G, Z,G modification
*
      IF(IFLAGS(IFDIAG).EQ.-1)    THEN
        FAA=0D0
        FZA=0D0
      ELSEIF(IFLAGS(IFDIAG).EQ.0) THEN
        FZA=0D0
      ENDIF
*
*     INTRF = 1 STANDARD MODEL INTERF.
*
* New option, MISD=1. Here simple bypass of CALL ROKANC 
*
      IF(INTRF.GT.1.AND.IFLAGS(IFMISD).EQ.0) GOTO 11
*
* to have common logic in the SUBROUTINE BORN
*
      IND=0
*
      XFZ(1)=DCMPLX(1D0,0D0)
      XFZ(2)=DCMPLX(1D0,0D0)
      XFZ(3)=DCMPLX(1D0,0D0)
      XFZ(4)=DCMPLX(1D0,0D0)
*
* EW FORM-FACTORS INITIALIZATION, SM BRANCH
* IF(NPAR(1).EQ.-1) the same as 1 and 2 here
      IF(NPAR(1).GE.1.OR.NPAR(1).EQ.-1) THEN
        IBOXF=0
        IBFLA=0
* Attention IBFLA setting !!!
        IF(INDF.EQ.9) IBFLA=1
        Q2=S/2D0
        IF(NPAR(4).EQ.2) THEN
          IBOXF=1
          Q2=S/2D0*(1D0-COCOST)
        ENDIF
        U2=Q2-S
        CALL ROKANC(IBOXF,IBFLA,U2,-S,-Q2,QE,QF,XFZ,XFZT,XFZTA)
      ENDIF
*
      IND=INDARR(INDF)
*
      GOTO 20
*
 11   CONTINUE
*
*     INTRF > 1 CROSS SECTION AND ASYMMETRY, MI BRANCHES
*
      IND=INDARR(INDF)
*
      XFZ(1)=XALLCH(IND,1)
      XFZ(2)=XALLCH(IND,2)
      XFZ(3)=XALLCH(IND,3)
      XFZ(4)=XALLCH(IND,4)
*
 20   CONTINUE
*
* CORRECTIONS FOR RUNNING ALPHA_QED THE SAME FOR ALL INTERFACES.
*
      AL4PI=AL1PI/4
      AMZ2=AMZ*AMZ
      IF(IALEM.GE.2) THEN
        IF(IALE2.EQ.0)THEN
          XFZT =1D0+AL4PI*XFOTF3(IALEM,1    ,IHVP,IQCD,1,DAL5H,-S)
        ELSE
          XFZT =1D0+AL4PI*XFOTF3(IALEM,IALE2,IHVP,IQCD,1,DAL5H,-S)
        ENDIF  
* It returns XFZT(S) and makes rescaling if \alpha(M^2_Z) is released
      ELSE
        XFZT=XFOTF
* It returns XFZT(M^2_Z) and the release is inside XFOTF3 
      ENDIF
*
* Here XVPOL, VPOL2 and ALQEDS are filled: IF ALEM=0,1 at M^2_Z; 
*                                          IF ALEM=2,3 at S.
*
      XVPOL =1D0/(2D0-XFZT)
      VPOL2 =DREAL(XVPOL)**2+DIMAG(XVPOL)**2
      ALQEDS=1D0/ALFAI/(2D0-DREAL(XFZT))
*
* call QCDCOF
*
      IF(IFLAGS(IFMISD).EQ.1) THEN
        SQRS=SQRT(S)
        TMASS=ALLMS(8)
        SIN2TW=SIN2TW
        ALQED=ALQEDS
        ALFAST=ALPHST
      CALL QCDCOF(SQRS,TMASS,SIN2TW,ALQED,ALFAST,ALFATT,ALPHXI,QCDCCR)
      ENDIF
*
* CZAKUE corrections
*
      CZAKFF=ARCZAK(INDF) 
*
*--------------------------------------------------------
*  VECTOR AND AXIAL COUPLINGS INITIALIZATION
*--------------------------------------------------------
*
      CORQEE=1D0+.75D0*ALQEDZ/PI*QE2
* IF(NPAR(1).EQ.-1) the same as 2 here
      IF(NPAR(1).EQ.-1.OR.NPAR(1).EQ.0.OR.NPAR(1).EQ.2) THEN
        ADDIME=0D0
        ADDIMF=0D0
      ELSE
       IF(INDF.NE.9) THEN
        ADDIME=ALQEDS**2*35D0/18*(1D0-8D0/3*DREAL(XFZ(2))*SW2)
        ADDIMF=ALQEDS**2*35D0/18*(1D0-8D0/3*DREAL(XFZ(3))*SW2)
       ELSE
        ADDIME=0D0
        ADDIMF=0D0
       ENDIF
      ENDIF
*
      IF(INTRF.EQ.2.AND.INDF.EQ.10) GOTO 150
*
      XRO =XFZ(1)
      XROZ=XRO
      XVEZ=1D0-4D0*(SW2*XFZ(2)+ADDIME)*QEM
      XVFZ=1D0-4D0*(SW2*XFZ(3)+ADDIMF)*QFM
      XVEFZ=-1D0+XVEZ+XVFZ
     & +16D0*(SW2*SW2*XFZ(4)+SW2*(XFZ(2)*ADDIMF+XFZ(3)*ADDIME))*QEM*QFM
*cbardin: Very useful prints to study PO/RO-at-resonance mismatch
*cb      print *,'XRO  =',XRO
*cb      print *,'Re,Im=',SQRT(DCMPLX(AROTFZ(2)   ,AIROFZ(2)   )
*cb     &                     *DCMPLX(AROTFZ(INDF),AIROFZ(INDF)))
*cb      print *,'XVE  =',XVEZ
*cb      print *,'Re,Im=   ',ARVEFZ(2),   ',',AIVEFZ(2)
*cb      print *,'XVF  =',XVFZ
*cb      print *,'Re,Im=   ',ARVEFZ(INDF),',',AIVEFZ(INDF)
*cb      print *,'XVE*F=',XVEZ*XVFZ
*cb      print *,'XVEF =',XVEFZ
*cb      print *,'RVEF =   ',DREAL(XVEZ)*DREAL(XVFZ)
      RO2 =DREAL(XRO)**2+DIMAG(XRO)**2
*
      COLOR=ZPAR(3)
      AMF=ZPAR(4)
      IFACT=NPAR(19)
*
* because of new QCD-corrections the code below is heavily rewritten
*
      IF(INDF.GE.4.AND.INDF.LE.9) THEN
        RATG=1D0/9
        RATI=7D0/3-44D0/9*SW2
        RATZ=(1+4D0/3*SW2)**2
        IF(INDF.EQ.6)     THEN
          RMASSQ=CMQRUN(3)
        ELSEIF(INDF.EQ.9) THEN
          RMASSQ=CMQRUN(6)
        ELSE
          RMASSQ=0D0
        ENDIF
        IF(IFLAGS(IFMISD).EQ.0) THEN
          RQCDV=ZPAR(7+MAX(0,2*INDQ(INDF+1)-1))
          RQCDA=ZPAR(7+MAX(0,2*INDQ(INDF+1)  ))
          VSNGAA=RATG*ZPAR(20)
          VSNGZA=RATI*ZPAR(20)
          VSNGZZ=RATZ*ZPAR(20)
          RQCDAS=ZPAR(21)*RMASSQ/SQRT(S)
        ELSE
          RQCDV=QCDCCR(MAX(0,2*INDQ(INDF+1)-1))
          RQCDA=QCDCCR(MAX(0,2*INDQ(INDF+1)  ))
          VSNGAA=RATG*QCDCCR(13)
          VSNGZA=RATI*QCDCCR(13)
          VSNGZZ=RATZ*QCDCCR(13)
          RQCDAS=QCDCCR(14)*RMASSQ/SQRT(S)
        ENDIF
        IF(IFINAL.GE.1) THEN
* QUARKS, STILL, IF(IFINAL.EQ.0) THE ALQEDZ IS USED
          RQCDV=RQCDV-.75D0*ALQEDZ/PI*QF2
          RQCDA=RQCDA-.75D0*ALQEDZ/PI*QF2
        ENDIF
      ELSE
        RQCDV =1D0
        RQCDA =1D0
        VSNGAA=0D0
        VSNGZA=0D0
        VSNGZZ=0D0
        RQCDAS=0D0
      ENDIF
*
      IF(IFINAL.EQ.-1) THEN
        RQCDV =1D0
        RQCDA =1D0
        VSNGAA=0D0
        VSNGZA=0D0
        VSNGZZ=0D0
        RQCDAS=0D0
      ENDIF
*
      CORINT=COLOR*RQCDV
*
* POLAR
*
      COMB1=1D0-ALAMP*ALAMM
      COMB2=ALAMP-ALAMM
*
      HELPL = HELP
      HELMI = HELM
* PREPARATION FOR DIFFERENT HELICITY STATES
      IHP = INT(HELPL)
      IHM = INT(HELMI)
      IF(IHP.NE.0.AND.IHM.NE.0) THEN
       HELI1 = (1.D0-HELPL*HELMI)/4.D0
       HELI2 = (HELPL-HELMI)/4.D0
      ENDIF
      IF(IHP.EQ.0.AND.IHM.EQ.0) THEN
       HELI1 = 1.D0
       HELI2 = 0.D0
      ENDIF
      IF(IHP.EQ.0.AND.IHM.NE.0) THEN
       HELI1 = .5D0
       HELI2 = -.5D0*HELMI
      ENDIF
      IF(IHP.NE.0.AND.IHM.EQ.0) THEN
       HELI1 = .5D0
       HELI2 = .5D0*HELPL
      ENDIF
       HOMB1 = HELI1
       HOMB2 = HELI2
*
      GOTO (100,200,300,400,500,600) INTRF
*
*  INTRF = 1  STANDARD MODEL INTERF.
*
100   CONTINUE
*
*COUPLINGS (VACUUM POLARIZATION IS INCLUDED IN THE SUBROUTINE BORN)
*Controll helicities and polarisations for VEEZ!
*
* Vector singlet QCD-corrections, assigned democratically 1/5 per channel
*
      DXVEZ =XVEZ*DCONJG(XVEZ)
      DXVFZ =XVFZ*DCONJG(XVFZ)
      DXVEFZ=XVEFZ*DCONJG(XVEFZ)
      VZ1 =  (1+DXVEZ)*(RQCDA+1D0/5*VSNGZZ)
     &    +  (DXVFZ+DXVEFZ)*RQCDV
      VZ2 = 2*DREAL(XVEZ*RQCDA+XVFZ*DCONJG(XVEFZ)*RQCDV)
      VZ10= (1+XVEZ*DCONJG(XVEZ))
     &    + (XVFZ*DCONJG(XVFZ)+XVEFZ*DCONJG(XVEFZ))
      VZ20= 2*DREAL(XVEZ      +XVFZ*DCONJG(XVEFZ)      )
      IF(IFLAGS(IFASCR).EQ.0) THEN
       AZ1= 2*(DREAL(XVEZ)*DREAL(XVFZ)+DREAL(XVEFZ))*(1D0+RQCDAS)
      ELSE
       AZ1= 2*DREAL(XVEZ*DCONJG(XVFZ)+XVEFZ)*(1D0+RQCDAS) 
      ENDIF
      AZ2 = 2*DREAL(XVFZ+XVEZ*DCONJG(XVEFZ))
*
      VEEZ=COMB1*HOMB1*(1D0+XVEZ*DCONJG(XVEZ))*RO2*RQCDA*COLOR      *FZZ
      VEFZ=(COMB1*HOMB1*VZ1+COMB2*HOMB1*VZ2+COMB1*HOMB2*AZ2
     &     +COMB2*HOMB2*AZ1)*RO2*COLOR                              *FZZ
      VEFZ0=(COMB1*HOMB1*VZ10+COMB2*HOMB1*VZ20+COMB1*HOMB2*AZ2
     &     +COMB2*HOMB2*AZ1)*RO2*COLOR                              *FZZ
*
      AEFZ=(COMB1*HOMB1*AZ1+COMB2*HOMB1*AZ2+COMB1*HOMB2*VZ2+
     & COMB2*HOMB2*VZ1)*RO2*COLOR                                   *FZZ
* 
      XVEFI=(COMB1*HOMB1*(QEFM*XVEFZ*RQCDV+QEM*XVEZ*1D0/5*VSNGZA)
     &     +QEFM*(COMB2*HOMB1*XVFZ+COMB1*HOMB2*XVEZ+COMB2*HOMB2)*RQCDV)
     &     *                     XRO              *COLOR            *FZA
*
      XAEFI=QEFM*(COMB1*HOMB1 +COMB2*HOMB1*XVEZ+COMB1*HOMB2*XVFZ
     &     +COMB2*HOMB2*XVEFZ)*  XRO*(1D0+RQCDAS) *COLOR            *FZA
*
      VEFA=COMB1*HOMB1*(QEF**2*RQCDV+1D0/5*VSNGAA)*COLOR            *FAA
      AEFA=COMB2*HOMB2* QEF**2                    *COLOR            *FAA
*
      RETURN
*
*  INTRF = 2  CROSS SECTION ONLY (no initial and final polarizations)
*
200   CONTINUE
*
      IF(IFLAGS(IFMISC).EQ.0) THEN
        ROFZI=AROTFZ(INDF)
      ELSE
        ROFZI=ARROFZ(INDF)/RENFAC(INDF)
      ENDIF
*
      XVEFI=(QEM*QFM*XVEFZ*RQCDV+QEM*XVEZ*1D0/5*VSNGZA)*XRO*COLOR   *FZA
       VEFA=(QEF**2       *RQCDV+1D0/5*VSNGAA)             *COLOR   *FAA
      IF(IFLAGS(IFASCR).EQ.0) THEN 
       AEFZ=2*(DREAL(XVEZ)*DREAL(XVFZ)+DREAL(XVEFZ))*RO2
     &                                      *(1D0+RQCDAS)  *COLOR   *FZZ
      ELSE
       AEFZ=2*DREAL(XVEZ*DCONJG(XVFZ)+XVEFZ)*RO2
     &                                      *(1D0+RQCDAS)  *COLOR   *FZZ
      ENDIF
      XAEFI=QEM*QFM*(1D0+RQCDAS)*                       XRO*COLOR   *FZA
       AEFA=0D0                                                     *FAA
*
      VZ1 = (1+XVEZ*DCONJG(XVEZ))*(RQCDA+1D0/5*VSNGZZ)
     &    + (XVFZ*DCONJG(XVFZ)+XVEFZ*DCONJG(XVEFZ))*RQCDV
      VZ10= (1+XVEZ*DCONJG(XVEZ))       
     &    + (XVFZ*DCONJG(XVFZ)+XVEFZ*DCONJG(XVEFZ))       
*
      CORF=(AMF/AMZ)**2
      IF(INDF.GE.4) CORF=0D0
      CORF2=1D0+2D0*CORF
      THRES=SQRT(1D0-4D0*CORF)
      ANORM=GMU*AMZ**3/(24D0*PI*SQRT(2D0))
*
      VEEZ=((GAMLE+PARTZI(1)/1D3)/ANORM/CORQEE-WIDTHS(1)/1D3/ANORM
     &/CORQEE+(1D0+XVEZ*DCONJG(XVEZ))*RO2)*RQCDA*COLOR              *FZZ
      VEFZ=((GAMLE+PARTZI(1)/1D3)/ANORM/CORQEE*
     &((GAMFI+PARTZI(INDF)/1D3)/ANORM/THRES+6D0*ROFZI*CORF*RQCDA*COLOR)
     &/CORF2 -WIDTHS(1)/1D3/ANORM/CORQEE*
     &(WIDTHS(INDF)/1D3/ANORM/THRES+6D0*ROFZI*CORF*RQCDA*COLOR)/CORF2
     &     + VZ1*RO2*COLOR)                                         *FZZ
      VEFZ0=((GAMLE+PARTZI(1)/1D3)/ANORM/CORQEE*
     &((GAMFI+PARTZI(INDF)/1D3)/ANORM/THRES+6D0*ROFZI*CORF      *COLOR)
     &/CORF2 -WIDTHS(1)/1D3/ANORM/CORQEE*
     &(WIDTHS(INDF)/1D3/ANORM/THRES+6D0*ROFZI*CORF      *COLOR)/CORF2
     &     + VZ10*RO2*COLOR)                                        *FZZ
      RETURN
*
150   CONTINUE
*
* SPECIAL CHAIN TO CALCULATE THE TOTAL HADRONIC CROSS SECTION
*
      COMB1=1D0
      COMB2=0D0
      HOMB1=1D0
      HOMB2=0D0
*
      GAMHAD=0D0
      GAMHAI=0D0
      DO 9  I=4,9
      GAMHAD=GAMHAD+WIDTHS(I)/1D3
 9    GAMHAI=GAMHAI+PARTZI(I)/1D3
*
* For construction of "summed couplings"
*
      VEFA=0D0
      AEFA=0D0
      XVEFI=(0D0,0D0)
      XAEFI=(0D0,0D0)
      VEFZ=0D0
      VEFZ0=0D0
      AEFZ=0D0
*
      DO 10 I=4,9
      IF(I.EQ.8) GOTO 10
      QFH=ALLCH(I)
      QFHM=ABS(QFH)
       II=INDQUA(I)
*
* New option, MISD=1. Here important new CALL ROKANC
*
      IF(IFLAGS(IFMISD).EQ.0) GOTO 12
*
* EW FORM-FACTORS INITIALIZATION, SM BRANCH
* IF(NPAR(1).EQ.-1) the same as 1 and 2 here
      IF(NPAR(1).GE.1.OR.NPAR(1).EQ.-1) THEN
        IBOXF=0
        IBFLA=0
* Attention IBFLA setting !!!
        IF(I.EQ.9) IBFLA=1
        Q2=S/2D0
        IF(NPAR(4).EQ.2) THEN
          IBOXF=1
          Q2=S/2D0*(1D0-COCOST)
        ENDIF
        U2=Q2-S
*
* CPU saving arrangement: 2*u+2*d+b
*
        IF(I.EQ.4) THEN
          CALL ROKANC(IBOXF,IBFLA,U2,-S,-Q2,QE,QFH,XFZu,XFZT,XFZTA)
          DO IXZ=1,4
             XFZ(IXZ)=XFZu(IXZ)
          ENDDO
        ELSEIF(I.EQ.5) THEN
          CALL ROKANC(IBOXF,IBFLA,U2,-S,-Q2,QE,QFH,XFZd,XFZT,XFZTA)
          DO IXZ=1,4
             XFZ(IXZ)=XFZd(IXZ)
          ENDDO
        ELSEIF(I.EQ.6) THEN
          DO IXZ=1,4
             XFZ(IXZ)=XFZu(IXZ)
          ENDDO
        ELSEIF(I.EQ.7) THEN
          DO IXZ=1,4
             XFZ(IXZ)=XFZd(IXZ)
          ENDDO
        ELSEIF(I.EQ.9) THEN
          CALL ROKANC(IBOXF,IBFLA,U2,-S,-Q2,QE,QFH,XFZ,XFZT,XFZTA)
        ENDIF
      ENDIF
*
 12   CONTINUE
*
      IF(IFLAGS(IFMISD).EQ.0) THEN
        XRO=XALLCH(II,1)
      ELSE
        XRO=XFZ(1)
      ENDIF
      IF(I.NE.9) THEN
        ADDIME=ALQEDS**2*35D0/18*(1D0-8D0/3*DREAL(XFZ(2))*SW2)
        ADDIMF=ALQEDS**2*35D0/18*(1D0-8D0/3*DREAL(XFZ(3))*SW2)
      ELSE
        ADDIME=0D0
        ADDIMF=0D0
      ENDIF
      RO2=DREAL(XRO)**2 + DIMAG(XRO)**2
      IF(IFLAGS(IFMISD).EQ.0) THEN
        XVE =1D0-4D0*(SW2*XALLCH(II,2)+ADDIME)*QEM
        XVF =1D0-4D0*(SW2*XALLCH(II,3)+ADDIMF)*QFHM
        XVEF=-1D0+XVE+XVF
     &       +16D0*(SW2*SW2*XALLCH(II,4)
     &       +SW2*(XALLCH(II,2)*ADDIMF+XALLCH(II,3)*ADDIME))*QEM*QFHM
      ELSE
        XVE =1D0-4D0*(SW2*XFZ(2)+ADDIME)*QEM
        XVF =1D0-4D0*(SW2*XFZ(3)+ADDIMF)*QFHM
        XVEF=-1D0+XVE+XVF
     &       +16D0*(SW2*SW2*XFZ(4)
     &       +SW2*(XFZ(2)*ADDIMF+XFZ(3)*ADDIME))*QEM*QFHM
      ENDIF
*
      COLOR=3D0
      QFH2=QFH**2
      AMQ=ALLMS(I)
*
      RATG=1D0/9
      RATI=7D0/3-44D0/9*SW2
      RATZ=(1+4D0/3*SW2)**2
      IF(IFLAGS(IFMISD).EQ.1) THEN ! bug emulation
         RMASSQ=0D0
      ENDIF
      IF(I.EQ.6)     THEN
        RMASSQ=CMQRUN(3)
      ELSEIF(I.EQ.9) THEN
        RMASSQ=CMQRUN(6)
      ENDIF
*
      IF(IFLAGS(IFMISD).EQ.0) THEN
        RQCDV=ZPAR(7+MAX(0,2*INDQ(I+1)-1))
        RQCDA=ZPAR(7+MAX(0,2*INDQ(I+1)  ))
        VSNGAA=RATG*ZPAR(20)
        VSNGZA=RATI*ZPAR(20)
        VSNGZZ=RATZ*ZPAR(20)
        RQCDAS=ZPAR(21)*RMASSQ/SQRT(S)
      ELSE
        RQCDV=QCDCCR(MAX(0,2*INDQ(I+1)-1))
        RQCDA=QCDCCR(MAX(0,2*INDQ(I+1)  ))
        VSNGAA=RATG*QCDCCR(13)
        VSNGZA=RATI*QCDCCR(13)
        VSNGZZ=RATZ*QCDCCR(13)
        RQCDAS=QCDCCR(14)*RMASSQ/SQRT(S)
      ENDIF
*
      IF(IFINAL.GE.1) THEN
* QUARKS, STILL, IF(IFINAL.EQ.0) THE ALQEDZ IS USED
        RQCDV=RQCDV-.75D0*ALQEDZ/PI*QFH2
        RQCDA=RQCDA-.75D0*ALQEDZ/PI*QFH2
      ENDIF
*
      IF(IFINAL.EQ.-1) THEN
        RQCDV =1D0
        RQCDA =1D0
        VSNGAA=0D0
        VSNGZA=0D0
        VSNGZZ=0D0
        RQCDAS=0D0
      ENDIF
*
      CORINT=COLOR*RQCDV
*
      J=I-3
*
      XXVEFI(J)=QEM*(QFHM*XVEF*RQCDV+XVE*1D0/5*VSNGZA)*XRO*COLOR   *FZA
      AVEFA (J)=(QE2*QFH2*RQCDV +1D0/5*VSNGAA)            *COLOR   *FAA
*
      VZ1 = (1+XVE*DCONJG(XVE))*(RQCDA+1D0/5*VSNGZZ)
     &    + (XVF*DCONJG(XVF)+XVEF*DCONJG(XVEF))*RQCDV
      VZ10= (1+XVE*DCONJG(XVE))
     &    + (XVF*DCONJG(XVF)+XVEF*DCONJG(XVEF))
*
      ANORM=GMU*AMZ**3/(24D0*PI*SQRT(2D0))
      AVEEZ(J)=((GAMLE+PARTZI(1)/1D3)/ANORM/CORQEE-WIDTHS(1)/1D3/ANORM
     &/CORQEE + (1D0+XVE*DCONJG(XVE))*RO2)*RQCDA*COLOR             *FZZ
      GAMFII=WIDTHS(I)/1D3
      FACNOR=(GAMFI+GAMHAI)/GAMHAD
      AVEFZ(J) =((GAMLE+PARTZI(1)/1D3)/ANORM/CORQEE*GAMFII/ANORM
     &                 -WIDTHS(1)/1D3 /ANORM/CORQEE*GAMFII/ANORM
     &     + VZ1*RO2*COLOR)*FACNOR                                 *FZZ
      AVEFZ0(J)=((GAMLE+PARTZI(1)/1D3)/ANORM/CORQEE*GAMFII/ANORM
     &                 -WIDTHS(1)/1D3 /ANORM/CORQEE*GAMFII/ANORM
     &     + VZ10*RO2*COLOR)*FACNOR                                *FZZ
*
* summed couplings are needed for BOXINT(IFI=1) QFH takes QF from ALQEF
*
      IF(IFLAGS(IFMISD).EQ.1) THEN
        IF(IFLAGS(IFASCR).EQ.0) THEN
         AZ1= 2*(DREAL(XVE)*DREAL(XVF)+DREAL(XVEF))*(1D0+RQCDAS)
        ELSE
         AZ1= 2*DREAL(XVE*DCONJG(XVF)+XVEF)*(1D0+RQCDAS) 
        ENDIF
        AAEFZ(J)=AZ1*RO2*COLOR                                     *FZZ
        XXAEFI(J)=ABS(QFH)*XRO*(1D0+RQCDAS)         *COLOR         *FZA
        AAEFA(J)=0D0
        VEFA=VEFA  +AVEFA(J) *QFH
        VEFZ=VEFZ  +AVEFZ(J) *QFH
        VEEZ=VEEZ  +AVEEZ(J) *QFH
        VEFZ0=VEFZ0+AVEFZ0(J)*QFH
        XVEFI=XVEFI+XXVEFI(J)*QFH
        AEFA=AEFA  +AAEFA(J) *QFH
        XAEFI=XAEFI+XXAEFI(J)*QFH
        AEFZ=AEFZ  +AAEFZ(J) *QFH
      ENDIF
*
  10  CONTINUE
      RETURN
*
*  INTRF = 3  CROSS SECTION AND ASYMMETRY: Effective Couplings Languages
*
300   CONTINUE
*
      IF(IFLAGS(IFMISC).EQ.0) THEN
        ROE1=ROEC
        ROF1=ROFC
        GVE1=GVEC
        GVF1=GVFC
        GAE1=GAEC
        GAF1=GAFC
        ROFZ1=AROTFZ(1)
        ROFZI=AROTFZ(INDF)
      ELSE
        ROE1=ROEC/RENFAC(1)
        ROF1=ROFC/RENFAC(INDF)
        GVE1=GVEC/SRENFC(1)
        GVF1=GVFC/SRENFC(INDF)
        GAE1=GAEC/SRENFC(1)      
        GAF1=GAFC/SRENFC(INDF)   
        ROFZ1=ARROFZ(1)/RENFAC(1)
        ROFZI=ARROFZ(INDF)/RENFAC(INDF)
      ENDIF
      GAE2=GAE1**2
      GAF2=GAF1**2
*
      XRO = SQRT(ROE1*ROF1)-SQRT(ROFZ1*ROFZI)+XFZ(1)
      RO2 = DREAL(XRO)**2 + DIMAG(XRO)**2
      XVE = GVE1     -ARVEFZ(1)             +XVEZ
      XVF = GVF1     -          ARVEFZ(INDF)+XVFZ
      XVEF= GVE1*GVF1-ARVEFZ(1)*ARVEFZ(INDF)+XVEFZ
*
       VZ1 = (GAE2+XVE*DCONJG(XVE))*GAF2*(RQCDA+1D0/5*VSNGZZ)
     &     + (GAE2*XVF*DCONJG(XVF)+XVEF*DCONJG(XVEF))*RQCDV
       VZ2 = 2*DREAL(GAF2*GAE1*XVE*RQCDA+GAE1*XVF*DCONJG(XVEF)*RQCDV)
*
       VZ10= (GAE2+XVE*DCONJG(XVE))*GAF2        
     &     + (GAE2*XVF*DCONJG(XVF)+XVEF*DCONJG(XVEF))       
       VZ20= 2*DREAL(GAF2*GAE1*XVE      +GAE1*XVF*DCONJG(XVEF)      )
       IF(IFLAGS(IFASCR).EQ.0) THEN 
        AZ1=2*(DREAL(XVE)*DREAL(XVF)+DREAL(XVEF))*(1D0+RQCDAS)*GAE1*GAF1
       ELSE
        AZ1=2*DREAL(XVE*DCONJG(XVF)+XVEF)*GAE1*GAF1*(1D0+RQCDAS)
       ENDIF 
       AZ2 =2*DREAL(GAE2*GAF1*XVF+GAF1*XVE*DCONJG(XVEF))
*
       VEEZ=COMB1*HOMB1*
     &     (GAE2+XVE*DCONJG(XVE))*COLOR*RQCDA*GAF2*RO2              *FZZ
       VEFZ=(COMB1*HOMB1*VZ1+COMB2*HOMB1*VZ2+COMB1*HOMB2*AZ2+
     & COMB2*HOMB2*AZ1)*RO2*COLOR                                   *FZZ
       VEFZ0=(COMB1*HOMB1*VZ10+COMB2*HOMB1*VZ20+COMB1*HOMB2*AZ2+
     & COMB2*HOMB2*AZ1)*RO2*COLOR                                   *FZZ
       AEFZ=(COMB1*HOMB1*AZ1+COMB2*HOMB1*AZ2+COMB1*HOMB2*VZ2+
     & COMB2*HOMB2*VZ1)*RO2*COLOR                                   *FZZ
*
      XVEFI=(COMB1*HOMB1*(QEFM*XVEF*RQCDV+QEM*XVE*1D0/5*VSNGZA)
     &     +QEFM*(COMB2*HOMB1*XVF*GAE1+COMB1*HOMB2*XVE*GAF1+
     &         GAE1*GAF1*COMB2*HOMB2)*RQCDV)*XRO*COLOR              *FZA
      XAEFI=QEFM*
     &(COMB1*HOMB1*GAE1*GAF1+COMB2*HOMB1*XVE*GAF1+COMB1*HOMB2*XVF*GAE1
     &+COMB2*HOMB2*XVEF)*XRO*(1D0+RQCDAS)          *COLOR           *FZA
       VEFA=COMB1*HOMB1*(QEF**2*RQCDV+1D0/5*VSNGAA)*COLOR           *FAA
       AEFA=QEF**2*COMB2*HOMB2*                     COLOR           *FAA
*
      RETURN
*
*  INTRF = 4  CROSS SECTION AND ASYMMETRY, LEPTONS ONLY, GA**2 AND GV**2
*             APPLICABLE ONLY FOR LEPTONS, LEPTONS UNIVERSALITY IS ASSUMED
*
400   CONTINUE
*
      IF(IFLAGS(IFMISC).EQ.0) THEN      
        SGVL2=SQRT(GVL2)
        SGVL2=SQRT(GVL2)
        ROL2R=ROL2
        GAL2R=GAL2
        GVL2R=GVL2
        ROFZ1=AROTFZ(1)
        ROFZI=AROTFZ(INDF)
      ELSE
        SGVL2=SQRT(GVL2)/SRENFC(1)
        SGVL2=SQRT(GVL2)/SRENFC(1)
        ROL2R=ROL2/RENFAC(1)**2
        GAL2R=GAL2/RENFAC(1)
        GVL2R=GVL2/RENFAC(1)
        ROFZ1=ARROFZ(1)/RENFAC(1)
        ROFZI=ARROFZ(INDF)/RENFAC(INDF)
      ENDIF
*
      XRO = SQRT(ROL2R)-SQRT(ROFZ1*ROFZI)+XFZ(1)
      RO2 = DREAL(XRO)**2 + DIMAG(XRO)**2
      XVE = SGVL2-ARVEFZ(1)             +XVEZ
      XVF = SGVL2          -ARVEFZ(INDF)+XVFZ
      XVEF= GVL2R-ARVEFZ(1)*ARVEFZ(INDF)+XVEFZ
*
       VEEZ= (GAL2R+XVE*DCONJG(XVE))*GAL2R*RO2                      *FZZ
       VEFZ=((GAL2R+XVE*DCONJG(XVE))*GAL2R        
     &      +(GAL2R*XVF*DCONJG(XVF)+XVEF*DCONJG(XVEF)) )
     &         *RO2*COLOR                                           *FZZ
       VEFZ0=VEFZ
*
       IF(IFLAGS(IFASCR).EQ.0) THEN
       AEFZ=2*GAL2R*(DREAL(XVE)*DREAL(XVE)+DREAL(XVEF))*RO2         *FZZ
       ELSE
       AEFZ=2*GAL2R*DREAL(XVE*DCONJG(XVE)+XVEF)*RO2                 *FZZ
       ENDIF
      XVEFI=QEFM*XVEF*           XRO                                *FZA
      XAEFI=QEFM*GAL2R*          XRO                                *FZA
       VEFA=QEF**2*                                                 *FAA
       AEFA=0D0                                                     *FAA
*
      RETURN
*
*  INTRF = 5  PARTIAL WIDTHS AND FORWARD-BACKWARD ASYMMETRY (FOUR)
*             APPLICABLE ONLY FOR LEPTONS
*
500   CONTINUE
*
      IF(IFLAGS(IFMISC).EQ.0) THEN    
        VAE2R=VAE2
        VAF2R=VAF2
        FOURR=FOUR
        ROFZ1=AROTFZ(1)
        ROFZI=AROTFZ(INDF)        
      ELSE
        VAE2R=VAE2/RENFAC(1)
        VAF2R=VAF2/RENFAC(INDF)
        FOURR=FOUR/RENFAC(1)/RENFAC(INDF)
        ROFZ1=ARROFZ(1)/RENFAC(1)
        ROFZI=ARROFZ(INDF)/RENFAC(INDF)
      ENDIF      
       VEEZ=(VAE2R-(1+ARVEFZ(1)**2)*ROFZ1
     &            +(1+XVEZ*DCONJG(XVEZ))*RO2)                       *FZZ
       VEFZ=(VAE2R*VAF2R
     &     -(1+ARVEFZ(1)**2)*(1+ARVEFZ(INDF)**2)*ROFZ1*ROFZI
     &     +(1+XVEZ*DCONJG(XVEZ)+XVFZ*DCONJG(XVFZ)+XVEFZ*DCONJG(XVEFZ))
     &                                                 *RO2)        *FZZ
       VEFZ0=VEFZ
       IF(IFLAGS(IFASCR).EQ.0) THEN 
       AEFZ=4*(FOURR-ARVEFZ(1)*ARVEFZ(INDF)*ROFZ1*ROFZI)
     &     +2*(DREAL(XVEZ)*DREAL(XVFZ)+DREAL(XVEFZ))*RO2            *FZZ
       ELSE
       AEFZ=4*(FOURR-ARVEFZ(1)*ARVEFZ(INDF)*ROFZ1*ROFZI)
     &     +2*DREAL(XVEZ*DCONJG(XVFZ)+XVEFZ)*RO2                    *FZZ
       ENDIF
      XVEFI=QEFM*XVEFZ*XRO                                          *FZZ
      XAEFI=QEFM*XRO                                                *FZA
       VEFA=QEF**2                                                  *FAA
       AEFA=0D0                                                     *FAA
*
      RETURN
*
*  INTRF = 6, EMPTY
*
600   CONTINUE
*                                                             END EWCOUP
      END

      SUBROUTINE BORN(IFINLA,R1,R2,SBORN,ABORN,SBORNS,ABORNS)
*     ========== ============================================
      IMPLICIT COMPLEX*16(X)
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
*
      COMMON / SVAR/ S
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /CHARGZ/ QE,QF,QEM,QFM,QEF,QEFM,QE2,QF2
      COMMON /MASSZ / AME,AMF,AME2,AMF2
      COMMON /COUPL / VEFA,XVEFI,VEFZ,AEFA,XAEFI,AEFZ,VEEZ,XVPOL,VPOL2
      COMMON /FORCHI/ XKAPP,XKAPPC,XMZ2,XMZ2C
      COMMON /PCONST/ ALFAI,AL1PI,ALQE2,ALQF2,ALQEF,GMU,CSIGNB
      COMMON /SMATRS/ AMZS,GAMZS,RESR,RESI,RES0,RES1,RES2,RESG,ISMA
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
      COMMON /INTRFS/ INTRF
      COMMON /CDZRUN/ CMQRUN(8)
      COMMON /PSCONS/ SW2,AMZ,GAMZ
*
* For total hadronic cross-section
*
      COMMON /INDFIT/ IND,INDF
      COMMON /ZFCHMS/ ALLCH(0:11),ALLMS(0:11)
      COMMON /HADRON/ XXVEFI(6),XXAEFI(6),AVEFA(6),AAEFA(6),
     &                AVEEZ(6),AVEFZ(6),AAEFZ(6)
      COMMON /COUPL0/ VEFZ0,AVEFZ0(6)
      COMMON /CZAKCO/ CZAKFF
      COMMON /CDAL5H/ DAL5H
      COMMON /CALQED/ ALQEDZ,ALQEDS
      COMMON /KAPPAC/ AKAPPA
      COMMON /FORSPR/ ECUT,ACUT
*
* M. JACK, 18/03/1999 13:00 ADDED COMMON/INDFIN/:
*
      COMMON /INDFIN/ IFUNFIN
*
* flags
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
      S1=S*R1
      S2=S*R2
      AMZ2 = AMZ*AMZ
      GAMZ2=GAMZ*GAMZ
*
* Running electromagnetic coupling
*
* Here XVPOLS, VPOL2S and ALQEDS are: IF ALEM=0,1 at M^2_Z; 
*                                     IF ALEM=2,3 at S.
*
      XVPOLS=XVPOL
      VPOL2S=VPOL2
*
      IF(IFLAGS(IFCONV).GE.1.AND.IFLAGS(IFALEM).GE.2) THEN
        IF(S1.LT.1D-2) THEN
          XFOT = DCMPLX(1D0,0D0)
        ELSE
          IF(IFLAGS(IFALE2).EQ.0) THEN
            XFOT = 1D0+AL1PI/4D0*XFOTF3
     &      (IFLAGS(IFALEM),             1,IFLAGS(IFVPOL),1,1,DAL5H,-S1)
          ELSE
            XFOT = 1D0+AL1PI/4D0*XFOTF3
     &      (IFLAGS(IFALEM),IFLAGS(IFALE2),IFLAGS(IFVPOL),1,1,DAL5H,-S1)
          ENDIF
        ENDIF
*
* Here XVPOLS, VPOL2S and ALQEDS are at S'
*
        XVPOLS=1D0/(2D0-XFOT)
        VPOL2S=DREAL(XVPOLS)**2+DIMAG(XVPOLS)**2
        ALQEDS=1D0/ALFAI/(2D0-DREAL(XFOT))
      ELSEIF(IFLAGS(IFCONV).EQ.-1) THEN
        XVPOLS=DCMPLX(1D0,0D0)
        VPOL2S=1D0
      ENDIF
*
* Running EW couplings
*
* MG's fix
      IF(IFLAGS(IFCONV).EQ.2) CALL EWCOUP(INTRF,INDF,MAX(1D2,S1))
*
* Czarnecki-Kuehn corrections
*
       IF(IFLAGS(IFCZAK).EQ.0.OR.IFLAGS(IFCZAK).EQ.2
     &                       .OR.IFINLA.EQ.-1) THEN
         CZAKUE=0D0
       ELSEIF(IFLAGS(IFCZAK).EQ.1) THEN
         CZAKUE=CZAKFF
       ENDIF
*
      IF(ISMA .EQ. 0) THEN
*
        XCHI1=XKAPP *S1/(S1-XMZ2 )
        XCHI2=XKAPPC*S2/(S2-XMZ2C)
        XCHI =XCHI1+XKAPP*S2/(S2-XMZ2)
*
* to have the same expression for fit and analytics
*
        IF(INDF.NE.10.OR.IND.EQ.0) THEN
*
* Regular chain for all INDF and INTRF
*
          SBORN=0D0
          ABORN=0D0
*
          IF(INDF.NE.6.AND.INDF.NE.9) THEN
            AMFH2=(ALLMS(INDF))**2
          ELSE
            AMFH2=(CMQRUN(INDF-3))**2
          ENDIF           
          QFH2=(ALLCH(INDF))**2
          IF(4D0*AMFH2.GE.S1) RETURN
*
* Final State QED corrections, governed by IFINLA=-1 or IFINAL
*
          SFIN=1D0
          AFIN=1D0
          IF(IFINLA.EQ.0) THEN
            IF(INDF.GE.1.AND.INDF.LE.3) THEN
              SFIN=1D0+3D0/4*ALQEDS/PI*QFH2 
            ELSE
              SFIN=1D0
            ENDIF
          ELSEIF(IFINLA.EQ.1.AND.INDF.NE.0) THEN
            CALL FUNFIN(S,AMFH2,QFH2,R1,SFIN,AFIN)
          ENDIF       
*
          AMF2S1=0D0
          IF(INDF.LE.3) THEN
            IF(IFLAGS(IFPOWR).EQ.1) AMF2S1=AMF2/S1
          ELSE
            IF(IFLAGS(IFFINR).EQ.-1.AND.IFLAGS(IFPOWR).EQ.1)
     &                              AMF2S1=(CMQRUN(INDF-3))**2/S1
          ENDIF
          THRESH=SQRT(MAX(1D0-4D0*AMF2S1,0D0))
          CORF2 =1D0+2D0*AMF2S1
          CORF3 =   -6D0*AMF2S1
*
*         BORN CROSS-SECTION
*
          SBORN=THRESH*(
     &    SFIN*CORF2*(VEFA*VPOL2S
     &                 +DREAL(XVEFI*XCHI*DCONJG(XVPOLS))
     &               )
     &   +(CORF2*(VEFZ+(SFIN-1D0+CZAKUE)*VEFZ0)+CORF3*VEEZ)
     &                                          *DREAL(XCHI1*XCHI2)
     &                 )/R1
*
*         print *,'INDF=',INDF,CORF2,CORF3,VEFZ
*
*         BORN ASYMMETRY
*
          ABORN=THRESH**2*AFIN*(AEFA*VPOL2S
     &   +DREAL(XAEFI*XCHI*DCONJG(XVPOLS))+AEFZ*DREAL(XCHI1*XCHI2))/R1
*
        ELSE
*
* Special chain for the total hadronic cross-section for INTRF=2 and INDF=10
*
          SBORN=0D0
          ABORN=0D0
*
          DO 1 I=4,9
            IF(I.EQ.8) GO TO 1
*
            IF(I.NE.6.AND.I.NE.9) THEN
              AMFH2=(ALLMS(I))**2
            ELSE
              AMFH2=(CMQRUN(I-3))**2
            ENDIF           
            QFH2=(ALLCH(I))**2
            IF(4D0*AMFH2.GE.S1) GOTO 1
*
* Final State QED corrections
*
            SFIN=1D0
            IF(IFINLA.EQ.0) THEN
              IF(I.GE.1.AND.I.LE.3) THEN
                SFIN=1D0+3D0/4*ALQEDS/PI*QFH2 
              ELSE
                SFIN=1D0
              ENDIF
            ELSEIF(IFINLA.EQ.1) THEN
              CALL FUNFIN(S,AMFH2,QFH2,R1,SFIN,AFIN)       
            ENDIF
*
            AMQ2=0D0
            IF(IFLAGS(IFFINR).EQ.-1.AND.IFLAGS(IFPOWR).EQ.1)
     &                              AMQ2=ALLMS(I)**2
*
            THRESH=SQRT(MAX(1D0-4D0*AMQ2/S1,0.D0))
            CORF2 =(1D0+2D0*AMQ2/S1)
            CORF3 =(   -6D0*AMQ2/S1)
*
*           BORN CROSS-SECTION
*
            J=I-3
            SBORN=SBORN+THRESH*(SFIN*CORF2*(AVEFA(J)*VPOL2S
     &           +DREAL(XXVEFI(J)*XCHI*DCONJG(XVPOLS))
     &                                     )
     &           +(CORF2*(AVEFZ(J)+(SFIN-1D0+CZAKUE)*AVEFZ0(J))
     &           +CORF3*AVEEZ(J)
     &            )*DREAL(XCHI1*XCHI2)
     &                         )/R1
*
*           BORN ASYMMETRY needed for IFI
*
            ABORN=ABORN+THRESH**2*AFIN*(AAEFA(J)*VPOL2S
     &           +DREAL(XXAEFI(J)*XCHI*DCONJG(XVPOLS))
     &           +AAEFZ(J)*DREAL(XCHI1*XCHI2))/R1
 1          CONTINUE
*
        ENDIF
*
        SBORNS=0D0
        ABORNS=0D0
        IF(INDF.EQ.10.AND.INTRF.EQ.2) THEN
*
* New addition for INTRF=2 and INDF=10 
*
*         BORN CROSS-SECTION
*
          SBORNS=(VEFA*VPOL2S+DREAL(XVEFI*XCHI*DCONJG(XVPOLS))
     &   +VEFZ*DREAL(XCHI1*XCHI2))/R1
*
*         BORN ASYMMETRY
*
          ABORNS=(AEFA*VPOL2S
     &   +DREAL(XAEFI*XCHI*DCONJG(XVPOLS))+AEFZ*DREAL(XCHI1*XCHI2))/R1
*
        ENDIF
*
      ELSE
*      
      print *,'USING S-MATRIX APPROACH VIA SMATASY'
      STOP
*
      ENDIF
*                                                          END BORN 
      END

      SUBROUTINE BORNN(R1,R2,C,SIGBRN)
*     ================================
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      COMMON / SVAR/ S
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /CHARGZ/ QE,QF,QEM,QFM,QEF,QEFM,QE2,QF2
      COMMON /MASSZ / AME,AMF,AME2,AMF2
      COMMON /COUPL / VEFA,XVEFI,VEFZ,AEFA,XAEFI,AEFZ,VEEZ,XVPOL,VPOL2
      COMMON /COUPL2/ XROZ,XVEZ
      COMMON /FORCHI/ XKAPP,XKAPPC,XMZ2,XMZ2C
      COMMON /PCONST/ ALFAI,AL1PI,ALQE2,ALQF2,ALQEF,GMU,CSIGNB
      COMMON /PSCONS/ SW2,AMZ,GAMZ
*
* flags
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
      S1=S*R1
      S2=S*R2
*
      IF(IFLAGS(IFCONV).EQ.2) CALL EWCOUP(INTRF,INDF,MAX(1D2,S1))
*
      XCHI1=XKAPP *S1/(S1-XMZ2 )
      XCHI2=XKAPPC*S2/(S2-XMZ2C)
*
*     BORN S-CHANNEL CROSS-SECTION
*
      SBORN0=VEFZ*DREAL(XCHI1*XCHI2)/R1
      ABORN0=AEFZ*DREAL(XCHI1*XCHI2)/R1
*
* New branch for IBA in ee\to\nue\nue 
*
* S -channel Born
      SIGS=3D0/4*CSIGNB*S*(1D0/2*(SBORN0*(1D0+C**2)+ABORN0*2D0*C)/S**2)
* T -channel Born. Attention! TT and UU are usual Mandelstamm variables. 
      BETAE=SQRT(1D0-4D0*AME2/S)
      TT=AME2-S/2D0*(1D0-BETAE*C)
      UU=AME2-S/2D0*(1D0+BETAE*C)
      AMW2=AMZ**2*(1D0-SW2)
      AKAP=GMU/(SQRT(2D0)*8D0*PI)*ALFAI
      IF(IFLAGS(IFWEAK).EQ.0.OR.IFLAGS(IFWEAK).EQ.-1) THEN
        ROW=1D0
      ELSEIF(IFLAGS(IFWEAK).EQ.1.OR.IFLAGS(IFWEAK).EQ.2) THEN
        ALSMW=LOG(S/AMW2) 
        CALL RHOCC(UU,-TT,S,0D0,-1D0,-1D0,0D0,ROW)
        ROW=ROW+AL1PI*(-3D0/4*ALSMW+1D0/4*(LOG(-TT/S))**2-2D0*F1+1D0)
        IF(ROW.GT.1.30D0) ROW=1.30D0
        IF(ROW.LT.0.70D0) ROW=0.70D0
      ELSE
      ENDIF
      XROW=DCMPLX(ROW,0D0)
      ROW2=ROW**2
      GAMW=0D0
      XCHIW=AKAP*AMW2/(AMW2-DCMPLX(0D0,1D0)*GAMW-TT)
      CHIW2=DREAL(XCHIW*DCONJG(XCHIW))
      SIGT =3D0/4*CSIGNB*S*(8D0*CHIW2*ROW2*(1D0+C)**2)
* ST-interference
      XCHIZ=XKAPP/(XMZ2-S) ! Tord, XCHIZ changed the sign, right?
      SIGST=3D0/4*CSIGNB*S*(-4D0*DREAL(XROZ*XCHIZ*(1D0+XVEZ)
     &          *DCONJG(XROW*XCHIW))*(1D0+C)**2)
* Three ENUE variants
      IF(IFLAGS(IFENUE).EQ.-1) THEN
        SIGBRN=SIGS
      ELSEIF(IFLAGS(IFENUE).EQ.0) THEN
        SIGBRN=SIGS+SIGT
      ELSEIF(IFLAGS(IFENUE).EQ.1) THEN
        SIGBRN=SIGS+SIGT+SIGST
      ENDIF
* Three SIGBRN (S-channel) are equal
**     CHIZ2=DREAL(XKAPP/(S-XMZ2)*XKAPPC/(S-XMZ2C))
**     SIGBRN=3D0/8*CSIGNB*S*CHIZ2*(VEFZ*COPL2+2*C*AEFZ)
**     SIGBRN=3D0/8*CSIGNB*S*DREAL(XCHIZ*DCONJG(XCHIZ))
**   &       *(VEFZ*COPL2+2*C*AEFZ)
** CSIGNB=4D0/3*PI/ALFAI**2*CONHC
** VEFZ =2*|rho|^2*(1+|v_e|^2)
** AEFZ =4*|rho|^2*Re(v_e)
** XFEFI=v_e in the case of ee->\nu\nu
** AKAPPA=GMU*AMZ**2/(SQRT(2D0)*8D0*PI)*ALFAI
*
      END

      FUNCTION FCROS(X)
*     ======== ========
*  
* Modified by M. Jack on 03/03/99 05:00pm
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /PCONST/ ALFAI,AL1PI,ALQE2,ALQF2,ALQEF,GMU,CSIGNB
      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
      COMMON /BORN0I/ SBORN0,ABORN0,SBORNI,ABORNI
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
      COMMON /CORINT/ CORINT
      COMMON /SOFTPR/ SOFTPR,SFPRFB
      COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT
      COMMON /RMINV / EPSHR,EPSSR
*
      COMMON /SVAR  / S
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /DELVAR/ DEL,DELM,DIF,DIF3,ALD,ALDM,ALDPL,ALDMI
     +               ,RD,RDI,RDI2,RDI3,RDM,RDMI,RDMI2,RDMI3
     +               ,ALRD,ALRDM
      COMMON /MASSZ / AME,AMF,AME2,AMF2
*
      COMMON /INDCUT/ INDC
      COMMON /INDEXA/ INDA,INDM
      COMMON /INDFIT/ IND,INDF
      COMMON /INTRFS/ INTRF
      COMMON /IMISDC/ IMISD
*-----------------------------------------------------------------------
*
      R=1D0-X
      SFACTI=X**(1D0-BETTAE)/BETTAE
*========================================================
* Choose different acceptance/acollinearity cut options:
*
      IF(INDC.LE.1) THEN
*
*************************************************
* 1. s' and M^2 cuts (original coding) :
*************************************************
*
        IF(IPHOT2.EQ.5) THEN
          SINI=4D0/3*FYFS(X)
        ELSE
          H0=-(1D0+R)*TE
          IF (IPHOT2.GE.0) H0=H0+SH2(R,1)
          SINI=4D0/3D0*((1D0+SOFTER)+ALQE2*H0*SFACTI)
        ENDIF
*
        ksoft=1
        IF(INTERF.GE.1.AND.R.LT.1D0-EPSHR) THEN
          H4=-4D0*(1D0+R)*R
          SHINTF=H4
        ENDIF
*
* MJ additions act here:
*
      ELSEIF(INDC.EQ.2) THEN
*
****************************************************
* 2. Acollinearity cut, but no acceptance cut (c=1):
****************************************************
*
        CALL SETR
        IF(R.GE.RECUT2) THEN
          ksoft=1
          INDA=1
          INDM=0
          CALL SHACOL (H0,H4,H5)
          SHINI =H0
          SHINTF=H4
          SHFIN =H5
        ELSE
*       ENERGY CUTS
          ksoft=0
          DEL=(RECUT2-R)*R1MI
          CALL SETDEL
          INDA=0
          CALL SHACOL (H0,H4,H5)
          SHINI =H0
          SHINTF=H4
          SHFIN =H5
        ENDIF
*     
        IF(R.LT.RACUT) THEN
*       ACOLLINEARITY CUT
          ksoft=0
          DEL=1D0/2D0*(1D0-SQRT(ABS(1D0-SINAC2
     &        *(1D0+2D0*R+R2)*R1MI2)/COSAC2))
          CALL SETDEL
          INDA=0
          CALL SHACOL (H0,H4,H5)
          SHINI =SHINI -H0
          SHINTF=SHINTF-H4
          SHFIN =SHFIN -H5
        ENDIF
*
        SINI=4D0/3*(1D0+SOFTER)*TET(R-RCUT)+ALQE2*SHINI*SFACTI
        SHINTF=SHINTF*X
*
      ELSE
        WRITE(*,*) 'WRONG VALUE FOR INDEX INDC !'
      ENDIF
*
      CALL BORN(IFINAL,R,R,SBORN,DUMMY,SBORNS,ABORNS)
      FCROS=SBORN*SINI
*
      IF (INTERF.GE.1.AND.R.LT.1D0-EPSHR) THEN
        CALL BORN(-1,R,1D0,DUMMY,ABORN,DUMMYS,ABORNS)
        IF(INDF.EQ.10.AND.INTRF.EQ.2.AND.IMISD.EQ.1) THEN
          ABORN=ABORNS
        ENDIF
        FCROS=FCROS+ALQEF*(SHINTF*ABORN+8*ABORNI*ksoft)/X*SFACTI
      ENDIF
*
      RETURN
      END
 
      FUNCTION FASYM(X)
*     ======== ========
*  
* Modified by M. Jack on 03/03/99 05:00pm
* Pairs are added (see SOFTER + SFPRFB) 15.10.99
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /PCONST/ ALFAI,AL1PI,ALQE2,ALQF2,ALQEF,GMU,CSIGNB
      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
      COMMON /BORN0I/ SBORN0,ABORN0,SBORNI,ABORNI
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
      COMMON /CORINT/ CORINT
      COMMON /SOFTPR/ SOFTPR,SFPRFB
      COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT
      COMMON /RMINV / EPSHR,EPSSR
*
      COMMON /SVAR  / S
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /DELVAR/ DEL,DELM,DIF,DIF3,ALD,ALDM,ALDPL,ALDMI
     +               ,RD,RDI,RDI2,RDI3,RDM,RDMI,RDMI2,RDMI3
     +               ,ALRD,ALRDM
      COMMON /MASSZ / AME,AMF,AME2,AMF2
*
      COMMON /INDCUT/ INDC
      COMMON /INDEXA/ INDA,INDM
      COMMON /INDFIT/ IND,INDF
      COMMON /INTRFS/ INTRF
      COMMON /IMISDC/ IMISD
*
* flags
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
*-----------------------------------------------------------------------
*
      R=1D0-X
      SFACTI=X**(1D0-BETTAE)/BETTAE
*========================================================
* Choose different acceptance/acollinearity cut options:
*
      IF (INDC.LE.1) THEN
*
****************************************************
* 1. s' and M^2 cuts (original coding) :
****************************************************
*
        R1P=1D0+R
        R2=R*R
        ALR=LOG(R)
        ALR1P=LOG(R1P)
        IF (R.GT.1D0-EPSHR) THEN
          H3=0D0
          SFACTI=0D0
        ELSE
          BR=4.D0*R/R1P**2
          ALBR=2.D0*AL2+ALR-2.D0*ALR1P
          H3=((1.D0+R2)*BR*(TE-LOG(BR))-2.D0*TE)/X
          IF (IPHOT2.GE.0) THEN
           IF(IFLAGS(IFFBHO).EQ.0) THEN
            H3=H3+AH2(R)
           ELSE
            YY= - (1D0-R)/(4D0*R)*( 4D0*LOG(1D0-R) + 3D0 )
            H3=H3+(AH2(R)+YY)*4D0*R/(1D0+R)**2
           ENDIF
          ENDIF
        ENDIF
*
        AINI=(1D0+SOFTER+SFPRFB)+ALQE2*H3*SFACTI
*
        ksoft=1
        IF(INTERF.GE.1.AND.R.LT.1D0-EPSHR) THEN
          H1=2D0/3*(2D0*R*(R2-4.D0*R+1.D0)/R1P+R*(5.D0*R2+3.D0)*ALR
     &    +(2.D0*R-5D0*(1.D0+R2))*R1P*ALR1P)
          AHINTF=H1
        ENDIF
*
* MJ additions act here:
*
      ELSEIF (INDC.EQ.2) THEN
*
****************************************************
* 2. Acollinearity cut, but no acceptance cut (c=1):
****************************************************
*
        CALL SETR
        IF(R.GE.RECUT2) THEN
          ksoft=1
          INDA=1
          INDM=0
          CALL AHACOL (H3,H1,H6)
          AHINI =H3
          AHINTF=H1
          AHFIN =H6
        ELSE
*       ENERGY CUTS
          ksoft=0
          DEL=(RECUT2-R)*R1MI
          CALL SETDEL
          INDA=0
          CALL AHACOL (H3,H1,H6)
          AHINI =H3
          AHINTF=H1
          AHFIN =H6
        ENDIF
*     
        IF (R.LT.RACUT) THEN
*       ACOLLINEARITY CUT
          ksoft=0
          DEL=1D0/2D0*(1D0-SQRT(ABS(1D0-SINAC2
     &    *(1D0+2*R+R2)*R1MI2)/COSAC2))
          CALL SETDEL
          INDA=0
          CALL AHACOL (H3,H1,H6)
          AHINI =AHINI -H3
          AHINTF=AHINTF-H1
          AHFIN =AHFIN -H6
        ENDIF
*
        AINI=(1D0+SOFTER+SFPRFB)*TET(R-RCUT)+ALQE2*AHINI*SFACTI
        AHINTF=AHINTF*X
*
      ELSE
        WRITE(*,*) 'WRONG VALUE FOR INDEX INDC !'
      ENDIF
*
      CALL BORN(IFINAL,R,R,DUMMY,ABORN,SBORNS,ABORNS)
      FASYM=AINI*ABORN
*
      IF (INTERF.GE.1.AND.R.LT.1D0-EPSHR) THEN
        CALL BORN(-1,R,1D0,SBORN,DUMMY,SBORNS,DUMMYS)
        IF(INDF.EQ.10.AND.INTRF.EQ.2.AND.IMISD.EQ.1) THEN
          SBORN=SBORNS
        ENDIF
        FASYM=FASYM+ALQEF*(AHINTF*SBORN
     &       -4D0/3*(-8D0*AL2-1D0)*SBORNI*ksoft)/X*SFACTI
      ENDIF
*
      END
 
      FUNCTION SSOFT(X)
*     ======== ========
      IMPLICIT REAL*8(A-H,O-Z)
*
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
*
      R=1D0-X
      CALL BORN(IFINAL,R,R,SBORN,DUMMY,SBORNS,ABORNS)
      IF(IPHOT2.EQ.5) THEN
        YFS=FYFS(X)
      ELSE
        YFS=1D0
      ENDIF
      SSOFT=SBORN*YFS
*
      END
 
      FUNCTION ASOFT(X)
*     ======== ========
      IMPLICIT REAL*8(A-H,O-Z)
*
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
*
      R=1D0-X
      CALL BORN(IFINAL,R,R,DUMMY,ABORN,SBORNS,ABORNS)
      ASOFT=ABORN
*
      END
 
      FUNCTION FACT(X)
*     ======== =======
      IMPLICIT REAL*8(A-H,O-Z)
*
      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
*     SOFT PHOTON FACTOR BETTAE*(1-R)**(BETTAE-1)
*     IS EXTRACTED (FUN. FACT)
*     TO IMPROVE SIMPSON INTEGRATION (SUB. FDSIMP)
      FACT=X**BETTAE
      END
 
      FUNCTION FACINV(F)
*     ======== =========
      IMPLICIT REAL*8(A-H,O-Z)
*
      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
      FACINV=F**(1/BETTAE)
      END

      FUNCTION SHARD(RR)
*     ======== =========
*
* Modified by M. Jack on 03/03/99 05:00pm
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      COMMON /PCONST/ ALFAI,AL1PI,ALQE2,ALQF2,ALQEF,GMU,CSIGNB
      COMMON /BORN0I/ SBORN0,ABORN0,SBORNI,ABORNI
      COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /DELVAR/ DEL,DELM,DIF,DIF3,ALD,ALDM,ALDPL,ALDMI
     +               ,RD,RDI,RDI2,RDI3,RDM,RDMI,RDMI2,RDMI3
     +               ,ALRD,ALRDM
      COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
      COMMON /CORINT/ CORINT
*
      COMMON /SVAR  / S
      COMMON /MASSZ / AME,AMF,AME2,AMF2
*
      COMMON /INDCUT/ INDC
      COMMON /INDEXA/ INDA,INDM
      COMMON /INDEX/  MF,IIF,JF,KF
      COMMON /INDFIT/ IND,INDF
      COMMON /INTRFS/ INTRF
      COMMON /IMISDC/ IMISD
*-----------------------------------------------------------------------
* 
      R=RR
      CALL SETR
*========================================================
* Choose different acceptance/acollinearity cut options:
*
      IF(INDC.LE.1) THEN
*
*************************************************
* 1. Acceptance cut and acollinearity cut 
*    (original coding by M.Bilenky) :
*************************************************
*
        IF(R.GE.RECUT2) THEN
          ksoft=1
          CALL SHFULL (H0,H4)
          SHINI =H0
          SHINTF=H4
        ELSE
* ENERGY CUTS
          ksoft=0
          DEL=(RECUT2-R)*R1MI
          CALL SETDEL
          CALL SHCUT (H0,H4)
          SHINI =H0
          SHINTF=H4
        ENDIF
*
        IF(R.LT.RACUT) THEN
*     ACOLLINEARITY CUT
          ksoft=0
          DEL=1D0/2*(1-SQRT(ABS(1-SINAC2
     &       *(1+2*R+R2)*R1MI2)/COSAC2))
          CALL SETDEL
          CALL SHCUT (H0,H4)
          SHINI =SHINI -H0
          SHINTF=SHINTF-H4
        ENDIF
*
* MJ additions act here:
*
      ELSEIF (INDC.EQ.3) THEN
c
c*************************************************
c* 3. Acceptance cut and acollinearity cut 
c*    (new coding by M.Jack) :
c*************************************************
c
        IF(R.GE.RECUT2) THEN
          ksoft=1
          CALL SHFULL (H0,H4)
          SHFIN =0D0
cc          INDA=1
cc          INDM=0
cc          DEL=1D0/2D0*(1D0-SQRT(1D0-4D0*AMF2/S/R))
cc          CALL SETDEL
cc          CALL SAVAR
cc          CALL PHASEREGC
cc          CALL SHCUTACOL(H0,H4,H5)
cc          SHFIN =H5
          SHINI =H0
          SHINTF=H4
        ELSE
*       ENERGY CUTS
          ksoft=0
          INDA=0
          DEL=(RECUT2-R)*R1MI
          CALL SETDEL
          CALL SAVAR
          CALL PHASEREGC
*========================================================
cc          CALL SHCUT (H0,DUMMY)
cc          CALL SHCUTACOL (DUMMY,H4,H5)
*========================================================
          CALL SHCUTACOL(H0,H4,H5)
          SHINI =H0
          SHINTF=H4
          SHFIN =H5
        ENDIF
*     
        IF(R.LT.RACUT) THEN
*       ACOLLINEARITY CUT
          ksoft=0
          INDA=0
          DEL=1D0/2D0*(1D0-SQRT(ABS(1D0-SINAC2
     &       *(1D0+2D0*R+R2)*R1MI2)/COSAC2))
          CALL SETDEL
          CALL SAVAR
          CALL PHASEREGC
*========================================================
cc          CALL SHCUT (H0,DUMMY)
cc          CALL SHCUTACOL (DUMMY,H4,H5)
*========================================================
          CALL SHCUTACOL(H0,H4,H5)
          SHINI =SHINI -H0
          SHINTF=SHINTF-H4
          SHFIN =SHFIN -H5
        ENDIF
*
      ELSE
        WRITE(*,*) 'WRONG VALUE FOR INDEX INDC !'
      ENDIF
*
cc      WRITE(*,*) 'INDC,INDA,R =',INDC,INDA,R
cc      WRITE(*,*) 'H0,H4,H5 =',H0,H4,H5
cc      WRITE(*,*) '--------------------------------------------'
c
*========================================================
*
      CALL BORN(IFINAL,R,R,SBORN,DUMMY,SBORNS,ABORNS)
*
      IF(IPHOT2.EQ.5) THEN
        SHARD=0D0
      ELSE
        SHARD=ALQE2*SHINI*SBORN
      ENDIF
*
      IF (INTERF.GE.1) THEN
        CALL BORN(-1,R,1D0,DUMMY,ABORN,DUMMYS,ABORNS)
        IF(INDF.EQ.10.AND.INTRF.EQ.2.AND.IMISD.EQ.1) THEN
          ABORN=ABORNS
        ENDIF
        SHARD=SHARD+ALQEF*(SHINTF*ABORN+4*(ALMI*(C2-1)+2*C)*ABORNI*R1MI
     &                                                   *ksoft)
      ENDIF
*
      END

      SUBROUTINE SHFULL(H0,H4)
*     ========== =============
      IMPLICIT REAL*8(A-H,O-Z)
*
      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM
      COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
*
      R1P=1D0+R
      R1PM2=1D0/R1P**2
      R1M=1D0-R
      R1M2=R1M*R1M
      R1M4=R1M2*R1M2
      R1PMI=R1P*R1MI
      R2P=1D0+R2
      R3P=1D0+R3
      R4P=1D0+R4
      TER=TE-ALR
*
      CC1=CM+R*CP
      CC1I3=1D0/CC1**3
      ALC1=LOG(CC1)
*
      CC2=CP+R*CM
      CC2I3=1D0/CC2**3
      ALC2=LOG(CC2)
*
      CC12I3=CC1I3*CC2I3
*
      H0=
     &+4D0/3*R2P*(ALC1*CC1I3*(R3*CP3-CM3)+ALC2*CC2I3*(CP3-R3*CM3))
     &+C*CC12I3*
     & (TER*R2P*(R3*4D0/3*(1-CPM)+R1M2*(2D0*R2*CPM+2*R2P*R*CPM2))
     &+R1M2*((4*R2P*R2-4*R*R4P-88D0/3*R3)*CPM2+16D0/3*R3*CPM)
     &+R1M4*(20D0/3*R2P*R-8*R2-4D0/3*R4P)*CPM3)
*
      IF (R.GE.RCUT) THEN
       H0=(H0-2*TE*COPL3)*R1MI
      IF (IPHOT2.GE.0) H0=H0+SH2(R,1)*COPL3
      ELSE
       H0=H0*R1MI
      END IF
*
      H4=0D0
      IF (INTERF.EQ.0) RETURN
      H4=4*(-C*R1PMI*R+ALR*C*R+ALMI*R1PMI*CPM*R2P
     &+(ALC2-ALC1)*(C2*R-R1M2*CPM))
*
      END
 
      SUBROUTINE SHCUT(H0,H4)
*     ========== ============
**  05/30/90 07:42pm
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CALL SETRC ( 1)
      CALL SHCUTM (H0M,H4M,I)
      H0=H0M
      H4=H4M
      CALL SETRC (-1)
      CALL SHCUTM (H0M,H4M,J)
      H0=H0-H0M
      H4=H4-H4M
*
      IF (I.EQ.J) RETURN
      CALL SJUMP (H0M,H4M,I,J)
      H0=H0+H0M
      H4=H4+H4M
*
      END
 
      SUBROUTINE SHCUTM(H0M,H4M,I)
*     ========== =================
**  05/31/90 07:13pm
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /DELVAR/ DEL,DELM,DIF,DIF3,ALD,ALDM,ALDPL,ALDMI
     +               ,RD,RDI,RDI2,RDI3,RDM,RDMI,RDMI2,RDMI3
     +               ,ALRD,ALRDM
      COMMON /RCVAR / C,C2,C3,CC,CCI,CCI2,CCI3,CCI4,CDEL,CDELM
     +               ,ALCC,ALCD,ALCDM,ALMI,ALPL
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
*
      ALCDMI=ALCDM-ALCD
      ALCDPL=ALCDM+ALCD
*
      H4M=0D0
      IF(CDEL.GT.0.D0.AND.CDELM.GT.0.D0) THEN
* REGION ( +,+ )
        I=1
        DH0M=2*DIF*C*R*(R4-4*R3+6*R2-4*R+1)
      IF (INTERF.GE.1) H4M=-2*C*R*(DIF*(R-1)+ALDMI*R)
      ELSE IF(CDELM.LT.0.D0.AND.CDEL.LT.0.D0) THEN
* REGION ( -,- )
        I=-1
        DH0M=-2*DIF*C*(R4-4*R3+6*R2-4*R+1)
      IF (INTERF.GE.1) H4M=2*C*(-DIF*(R-1)+ALDMI*R)
      ELSE
* REGION ( -,+ )
        I=0
        ALPM=2*ALCC+ALCDPL-2*ALR-ALPL+TE+1
        H0M=-(2*ALPM*(-1D0/3*(R5+R3+R2+1)+CCI*(R5+2*R3+R)
     . -CCI2*(R5+R4+R3+R2)+CCI3*2D0/3*(R5+R3))
     . +2D0/3*(ALCDMI*DIF3-ALMI)*(R5-3*R4+4*R3-4*R2+3*R-1)
     . +4*DEL*DELM*(C*(R4-3*R3+3*R2-R)-CCI*1D0/3*(R5-2*R4+2*R3-2*R2+R))
     . +2*DEL*C*(R5-5*R4+10*R3-10*R2+5*R-1)
     . +2D0/3*(C*(-R5+6*R4-11*R3+11*R2-6*R+1)
     . +CCI*(-3*R5-24*R4-14*R3-24*R2-3*R)
     . +CCI2*(9*R5+19*R4+19*R3+9*R2)+2*CCI3*(-3*R5-2*R4-3*R3)) )
      IF (INTERF.GE.1) THEN
        H4M=-(2*ALR*C*(R2-R)
     +   +ALCC*(C2*(R3+R2-R-1)-R3+3*R2-3*R+1)
     +   +1D0/2*ALMI*(C2*(R3+R2+R+1)-R3-R2-R-1)
     +   +C*(DIF*(-R3+R2+R-1)+R3+R2+R+1))
        H4M=H4M*R1MI
      END IF
        H0M=H0M*R1MI4
      RETURN
      END IF
*
      GENH0=2*ALCDMI*( 1D0/3*(R5+R3+R2+1)-CCI*(R5+2*R3+R)
     . +CCI2*(R5+R4+R3+R2)-CCI3*2D0/3*(R5+R3) )
     . +2D0/3*DIF3*( ALCDPL*(-R5+3*R4-4*R3+4*R2-3*R+1)
     . +C*(-R5+3*R4-2*R3-2*R2+3*R-1) )
     . +4D0/3*DIF*( CCI*(-R5+R)+CCI2*(R5-R4+R3-R2))
*
      H0M=I*GENH0+DH0M
      H0M=H0M*R1MI4
*
      END
 
      SUBROUTINE SJUMP(H0M,H4M,I,J)
*     =============================
**  05/30/90 07:12pm
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /DELVAR/ DEL,DELM,DIF,DIF3,ALD,ALDM,ALDPL,ALDMI
     +               ,RD,RDI,RDI2,RDI3,RDM,RDMI,RDMI2,RDMI3
     +               ,ALRD,ALRDM
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
*
      ALDPLE=ALDPL-TE-1
*
      H0MSUM=4D0/3*((ALDPLE+ALR)*DIF3+ALDMI)*(R5-3*R4+4*R3-4*R2+3*R-1)
     + +4D0/3*DIF3*(R5+R4-8*R3+8*R2-R-1)
     + +1D0/3*DIF*(-6*R5+22*R4-20*R3+20*R2-22*R+6)
      H0M=H0MSUM*R1MI4*ISIGN(1,I-J)
*
      H4M=0D0
      IF (INTERF.GE.1) THEN
       R6=R*R5
       H4MSUM=-2*(2*((ALRDM-ALRD)*2*(-R4+R3-R2+R)
     + +(ALRDM*RDMI-ALRD*RDI)*(R5+2*R4-2*R2-R)
     + +(ALRDM*RDMI2-ALRD*RDI2)*(-R5-R4+R3+R2))
     + +ALR*((RDMI-RDI)*(-3*R5-8*R4+2*R3+R)
     + +(RDMI2-RDI2)*(3*R5+3*R4-R3-R2))
     + +ALDMI*(4*(R4+R)+(RDMI+RDI)*(-R5-4*R4+2*R3-4*R2-R)
     + +(RDMI2+RDI2)*(R5+R4+R3+R2))
     + +(RDM-RD)*(R4-2*R3+2*R-1)+2*(RDMI-RDI)*(-R4+R2))
       H4M=H4MSUM*R1MI3*ISIGN(1,I-J)
      END IF
*
      IF (ABS(I-J).EQ.2) RETURN
*
       H0MDIF=4D0/3*ALR*(R5-3*R4+4*R3-4*R2+3*R-1)
     + +4D0/3*DEL*DELM*(R5+R4-2*R3-2*R2+R+1)
     + -2D0/3*(R5-15*R4-14*R3-14*R2-15*R+1)
       H0MDIF=H0MDIF*R1MI4*ISIGN(1,I-J)
*
      H4MDIF=0D0
      IF (INTERF.GE.1) THEN
       H4MDIF=-2*(2*((ALRDM+ALRD)*2*(-R4+R3-R2+R)
     + +(ALRDM*RDMI+ALRD*RDI)*(R5+2*R4-2*R2-R)
     + +(ALRDM*RDMI2+ALRD*RDI2)*(-R5-R4+R3+R2))
     + +ALR*(4*(3*R4-2*R3+2*R2-R)+(RDMI+RDI)*(-3*R5-8*R4+2*R3+R)
     + +(RDMI2+RDI2)*(3*R5+3*R4-R3-R2))
     + +ALDMI*((RDMI-RDI)*(-R5-4*R4+2*R3-4*R2-R)
     + +(RDMI2-RDI2)*(R5+R4+R3+R2))
     + +(RDM+RD)*(R4-2*R3+2*R-1)+2*(-R4+7*R3-7*R2+R)
     + +2*(RDMI+RDI)*(-R4+R2))
       H4MDIF=H4MDIF*R1MI3*ISIGN(1,I-J)
      END IF
*
      IF ((I+J).EQ. 1) THEN
       H0M=(H0M+H0MDIF)/2
      IF (INTERF.GE.1) H4M=(H4M+H4MDIF)/2
      ELSE
       H0M=(H0M-H0MDIF)/2
      IF (INTERF.GE.1) H4M=(H4M-H4MDIF)/2
      END IF
*
      END
 
      FUNCTION SH2(R,ISH2)
*     ======== ===========
*-----------------------------------------------------------------------
*     ALPHA**2*L**2 AND ALPHA**2*L PHOTONIC CORRECTIONS
*-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
*
      COMMON /PCONST/ ALFAI,AL1PI,ALQE2,ALQF2,ALQEF,GMU,CSIGNB
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
      COMMON /FLAGZP/ ISRPPR,IFSPPR,IFUNAN
      COMMON /CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
      COMMON /CPAIRS/ FPI2,RINF,RRR0,RRR1,PAIRDL,CORFAC
      COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT
      COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM
      COMMON /PRECIS/ NPREC
      EXTERNAL FUNAN2,FUNAN3
*
      SH2=0D0
      R1M=1D0-R
      IF(R1M.GT.1D-14) THEN
       R2=R*R
       R1P=1D0+R
       R2P=1D0+R2
       ALR=DLOG(R)
       R1MI=1D0/R1M
       ALR1M=DLOG(R1M)
       D2R1M=DDILOG(R1M)
       TL=TE+1D0
       IF(IPHOT2.GE.0) THEN
        IF(IFUNAN.EQ.0) THEN
         SH2A= ALQE2*
     +   TL**2*(-R2P*R1MI*ALR+R1P*(1D0/2D0*(ALR-1D0)-2D0*ALR1M)-2D0)
        ELSE
         SH2A = 0D0
         REPS1=1D-03/NPREC
CAC      REPS1=1D-05/NPREC
         RAPS1= REPS1*REPS1
         DELA = 1D-3
CAC      DELA = 1D-6
         ALDA = DLOG(DELA)
         Z2N  = R/(1D0-DELA)
         Z2X  = 1D0 - DELA
         Z2ST = 0.25D0*(Z2X-Z2N)
         IF(Z2ST.LE.0D0) THEN
          SH2A= ALQE2*
     +    TL**2*(-R2P*R1MI*ALR+R1P*(1D0/2D0*(ALR-1D0)-2D0*ALR1M)-2D0)
          GOTO 1
         ENDIF
         IF(ISH2.EQ.1) THEN
          IF(C.EQ.0D0) THEN
           SH2A= ALQE2*
     +     TL**2*(-R2P*R1MI*ALR+R1P*(1D0/2D0*(ALR-1D0)-2D0*ALR1M)-2D0)
           GOTO 1
          ENDIF
          CALL SIMPT(Z2N,Z2X,Z2ST,REPS1,RAPS1,FUNAN2,TT,AN2,RES2,RES3)
          ADD2= ( 0.5D0*P2(R) + P1(R)*P10(ALDA) )
     &         *( FUNANG(R,1D0,C) + FUNANG(1D0,R,C) )
     &        + P1(R)*FUNANG(R,1D0,C)*2D0*( - ALDA + DLOG(1D0-R) )
     &        + P1(R)*FUNANG(1D0,R,C)*2D0*( - ALDA + DLOG((1D0-R)/R)
     &          + 1D0/R - 1D0 )
          SH2A= ALQE2*TL**2/4D0*( ADD2 + AN2 )
     &          /2D0/(-COPL3)
     &        - ALQE2*TL**2*( R2P*R1MI*(2D0*ALR1M + 3D0/2D0) + 1.5D0*R1P
     &                      + 2D0*R1P*ALR1M ) 
         ELSEIF(ISH2.EQ.2) THEN
          CALL SIMPT(Z2N,Z2X,Z2ST,REPS1,RAPS1,FUNAN3,TT,AN3,RES2,RES3)
          ADD3= ( 0.5D0*P2(R) + P1(R)*P10(ALDA) )
     &         *( FUNANC(R,1D0,C) + FUNANC(1D0,R,C) )
     &        + P1(R)*FUNANC(R,1D0,C)*2D0*( - ALDA + DLOG(1D0-R) )
     &        + P1(R)*FUNANC(1D0,R,C)*2D0*( - ALDA + DLOG((1D0-R)/R)
     &          + 1D0/R - 1D0 )
          SH2A= ALQE2*TL**2/4D0*( ADD3 + AN3 )/(COPL2)
     &        - ALQE2*TL**2*( R2P*R1MI*(2D0*ALR1M + 3D0/2D0) + 1.5D0*R1P
     &                      + 2D0*R1P*ALR1M ) 
         ELSE
          PRINT *,'ISH2=',ISH2,'  IS NOT PROPERLY DEFINED FOR SH2' 
          STOP
         ENDIF
 1       CONTINUE
        ENDIF
        SH2=SH2 + SH2A
       ENDIF
       IF(IPHOT2.GE.1) THEN
        SH2=SH2+ALQE2*(
     +  +TL  *( R2P*R1MI*(D2R1M+ALR*(ALR1M+7D0/2D0-1D0/2D0*ALR))
     +        +R1P*(1D0/4D0*ALR**2+4D0*ALR1M-2D0*F1)-ALR+7D0+1D0/2D0*R))
       ENDIF
       IF(IPHOT2.GE.2) THEN
        D3R1M=TRILOG(R1M)
        S1R1M=   S12(R1M)
        SH2=SH2+
     +         ALQE2*(-2D0*R1P*(ALR1M-F1+1D0)+R2P*R1MI*(-1D0/6D0*ALR**3
     +         +.5D0*ALR*D2R1M+.5D0*ALR**2*ALR1M-1.5D0*D2R1M
     +         -1.5D0*ALR*ALR1M+F1*ALR-17D0/6D0*ALR-ALR**2)
     +         +R1P*(1.5D0*D3R1M-2D0*S1R1M-ALR1M*D2R1M-.5D0)
     +         -.25D0*(1D0-5D0*R)*ALR1M**2+.5D0*(1D0-7D0*R)*ALR*ALR1M
     +         -25D0/6D0*R*D2R1M+(-1D0+13D0/3D0*R)*F1+(1.5D0-R)*ALR1M
     +         +1D0/6D0*(11D0+10D0*R)*ALR+2D0*R1MI**2*ALR**2
     +         -25D0/11D0*R*ALR**2-2D0/3D0*R*R1MI*(1D0+R1MI*ALR)**2 )
       ENDIF
       IF(IPHOT2.EQ.3) THEN
         SH2=SH2+ALQE2**2/6D0*TE**3*(-27D0/2+15D0/4*R1M
     &      +4D0*(1D0-R1M/2)*(6D0*F1-6D0*ALR1M**2+3D0*D2R1M)
     &      +3D0*ALR*(7D0-6D0/R1M-3D0/2*R1M)
     &      +ALR**2*(-7D0+4D0/R1M+7D0/2*R1M)
     &      -6D0*ALR1M*(6D0-R1M)
     &      +6D0*ALR*ALR1M*(6D0-4D0/R1M-3D0*R1M))
       ELSEIF(IPHOT2.EQ.4) THEN
         CEILER=.5772157D0
         FEXPGL=EXP(-BETTAE*CEILER)/DGAMMA(1D0+BETTAE)
         DELVS1=ALQE2*(1.5D0*TL+2D0*F1-2D0)
         DELVS2=ALQE2**2*(TL**2*(9D0/8-2D0*F1)
     &   +TL*(-45D0/16+11D0/2*F1+3D0*ZET3) 
     &   -6D0/5*F1**2-9D0/2*ZET3-6*F1*AL2+3D0/8*F1+19D0/4)
         DELVS =1D0+DELVS1+DELVS2
         DELVSD=1D0+ALQE2**2*(TL*(3D0/16-1.5D0*F1+3D0*ZET3)
     &   -16D0/5*F1**2-4.5D0*ZET3-6D0*F1*ALR+51D0/8*F1+11D0/4)
         SH2=SH2+TE*(1D0-R)**(BETTAE-1D0)
     &      *(R2P*FEXPGL*EXP(DELVS1)*DELVSD-2D0*DELVS)
     &      +TE*R1P*(1D0+ALQE2*2D0*TE*ALR1M+DELVS1)
       ENDIF
*
* YFS is treated separately
* 
       IF(IPHOT2.EQ.5) THEN
         SH2=SH2
       ENDIF
*
* PAIRS are decoupled from photons if ISRPPR=1 
* 
* Below is the old treatment of PAIRS, ZFITTERs versions before 5.20
*
       IF(ISRPPR.EQ.-1) THEN
*
* Now changed back to RCUT:
*
         IF(R.GE.RCUT)  THEN
         TPIH=LOG(SCOM*R1M**2/FPI2/R)
         SH2=SH2
     &         +ALQE2/3*(1D0/2*R2P/R1M*TEE**2
     &      +(R2P/R1M*(2*ALR1M-ALR-5D0/3)-2*R1M)*TEE
     & +R2P/R1M*(1D0/2*(2*ALR1M-ALR)**2-5D0/3*(2*ALR1M-ALR)-2*F1+28D0/9)
     &  -R1M*(2*(2*ALR1M-ALR)-19D0/3)-R2/R1M*(ALR**2/2+D2R1M)-ALR)
     &      *TET((1-SQRT(R))-PAIRDL)
     &         +ALQE2/3*(1D0/2*R2P/R1M*TMU**2
     &      +(R2P/R1M*(2*ALR1M-ALR-5D0/3)-2*R1M)*TMU
     & +R2P/R1M*(1D0/2*(2*ALR1M-ALR)**2-5D0/3*(2*ALR1M-ALR)-2*F1+28D0/9)
     &  -R1M*(2*(2*ALR1M-ALR)-19D0/3)-R2/R1M*(ALR**2/2+D2R1M)-ALR)
     &      *TET((1-SQRT(R))-PAIRDL)*CORFAC
     &         +ALQE2/3*(1D0/2*R2P/R1M*TTU**2
     &      +(R2P/R1M*(2*ALR1M-ALR-5D0/3)-2*R1M)*TTU
     & +R2P/R1M*(1D0/2*(2*ALR1M-ALR)**2-5D0/3*(2*ALR1M-ALR)-2*F1+28D0/9)
     &  -R1M*(2*(2*ALR1M-ALR)-19D0/3)-R2/R1M*(ALR**2/2+D2R1M)-ALR)
     &      *TET((1-SQRT(R))-PAIRDL)
     &       +ALQE2/3*( R2P/R1M*(RINF*(1D0/2*TPIH**2-F1)+RRR0*TPIH+RRR1)
     &      -(1-R)*(RINF*(2*TPIH-3)+2*RRR0)
     &      -RINF*(R**2/R1M*(1D0/2*ALR**2+D2R1M)+ALR))
     &      *TET((1-SQRT(R))-PAIRDL)
        ENDIF
       ENDIF
      ENDIF
*
      END
 
      FUNCTION PH2(X)
*     ======== ======
*-----------------------------------------------------------------------
*     ALPHA**2 PAIR CORRECTIONS (HIGHER ORDERS CAN BE CALLED ALSO)
*-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
*
* flags
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
      COMMON /PCONST/ ALFAI,AL1PI,ALQE2,ALQF2,ALQEF,GMU,CSIGNB
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
      COMMON /FLAGZP/ ISRPPR,IFSPPR,IFUNAN
      COMMON /CPAIRS/ FPI2,RINF,RRR0,RRR1,PAIRDL,CORFAC
      COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT
      COMMON /CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
*
      R=1D0-X
      R2 =R*R
      R1P=1+R
      R2P=1+R2
      ALR=LOG(R)
      ALX=LOG(X)
      D2X=DDILOG(X)
*
      CALL BORN(IFINAL,R,R,SBORN,ABORN,SBORNS,ABORNS)
*
      TPI=DLOG(SCOM*X**2/FPI2/R)
*
      ALEADE = 1D0/2.*R2P/X*TEE**2
      ANLLAE = (R2P/X*(2.*ALX-ALR-5D0/3.)-2.*X)*TEE
     & +R2P/X*(1D0/2*(2*ALX-ALR)**2-5D0/3*(2*ALX-ALR)-2.*F1+28D0/9.)
     &   -X*(2*(2*ALX-ALR)-19D0/3)-R2/X*(ALR**2/2.+D2X)-ALR
      ALEADM = 1D0/2.*R2P/X*TMU**2
      ANLLAM = (R2P/X*(2.*ALX-ALR-5D0/3.)-2.*X)*TMU
     & +R2P/X*(1D0/2*(2*ALX-ALR)**2-5D0/3*(2*ALX-ALR)-2.*F1+28D0/9.)
     &   -X*(2*(2*ALX-ALR)-19D0/3)-R2/X*(ALR**2/2.+D2X)-ALR
      ALEADT = 1D0/2.*R2P/X*TTU**2
      ANLLAT = (R2P/X*(2*ALX-ALR-5D0/3)-2*X)*TTU
     & +R2P/X*(1D0/2*(2*ALX-ALR)**2-5D0/3*(2*ALX-ALR)-2.*F1+28D0/9.)
     &   -X*(2*(2*ALX-ALR)-19D0/3)-R2/X*(ALR**2/2.+D2X)-ALR
      ALEADH = R2P/X*(RINF*(1D0/2.*TPI**2))
      ANLLAH = R2P/X*(RINF*(-F1)+RRR0*TPI+RRR1)
     &   -X*(RINF*(2*TPI-3.)+2.*RRR0)
     &      -RINF*(R**2/X*(1D0/2*ALR**2+D2X)+ALR)

       IF(IFLAGS(IFISPP).LE.2) THEN
        IF(IFLAGS(IFIPFC).EQ.1) ALEAD = ALEADE
        IF(IFLAGS(IFIPFC).EQ.2) ALEAD = ALEADM
        IF(IFLAGS(IFIPFC).EQ.3) ALEAD = ALEADT
        IF(IFLAGS(IFIPFC).EQ.4) ALEAD = ALEADH
        IF(IFLAGS(IFIPFC).EQ.5) ALEAD = ALEADE+ALEADM+ALEADT+ALEADH
        IF(IFLAGS(IFIPFC).EQ.6) ALEAD = ALEADE+ALEADM+ALEADT
        IF(IFLAGS(IFIPFC).EQ.1) ANLLA = ANLLAE
        IF(IFLAGS(IFIPFC).EQ.2) ANLLA = ANLLAM
        IF(IFLAGS(IFIPFC).EQ.3) ANLLA = ANLLAT
        IF(IFLAGS(IFIPFC).EQ.4) ANLLA = ANLLAH
        IF(IFLAGS(IFIPFC).EQ.5) ANLLA = ANLLAE+ANLLAM+ANLLAT+ANLLAH
        IF(IFLAGS(IFIPFC).EQ.6) ANLLA = ANLLAE+ANLLAM+ANLLAT
       ELSEIF(IFLAGS(IFISPP).GE.3) THEN
        ALEAD = 0D0
        IF(IFLAGS(IFIPFC).EQ.4.OR.IFLAGS(IFIPFC).EQ.5) ALEAD = ALEADH
        ANLLA = 0D0
        IF(IFLAGS(IFIPFC).EQ.4.OR.IFLAGS(IFIPFC).EQ.5) ANLLA = ANLLAH
       ENDIF

       PH2=SBORN*ALQE2**2/3D0*( ALEAD + ANLLA )
*
      IF(ISRPPR.EQ.2.OR.ISRPPR.EQ.4) PH2 = PH2 + PH2ADD(X)*SBORN
      END

      FUNCTION FYFS(X)
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /PCONST/ ALFAI,AL1PI,ALQE2,ALQF2,ALQEF,GMU,CSIGNB
      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
*
      R=1D0-X
      R2=R**2
      R1P=1+R
      R2P=1+R2
      ALR=LOG(R)
      CEILER=.5772157D0
      FEXPGL=EXP((3D0/4-CEILER)*BETTAE+ALQE2*(2D0*F1-1D0/2))
     &      /DGAMMA(1D0+BETTAE)
      YFS0 =R2P/2D0
      YFSAB=ALQE2*BETTAE*(3D0/32-PI**2/8+3D0/2*ZET3)
*
      IF(X.GT.1D-20) THEN
       ALX=LOG(X)
       D2X=DDILOG(X)
       YFSA1=ALQE2/8*(-(1D0+3D0*R2)*ALR**2+4D0*R2P*(D2X+ALR*ALX)
     &                +2D0*X*(3D0-2D0*R)+2D0*(3D0+2D0*R+R2)*ALR)
       YFSB1=BETTAE/4*(-(1D0+3D0*R2)/2*ALR-X**2)
       YFSB2=BETTAE**2/8*(1D0/2*X*(1D0-3D0*R)*ALR
     &                +1D0/12*(1D0+7D0*R2)*ALR**2+X*R1P*D2X+X**2)
       FYFS=FEXPGL*(YFS0+YFSA1+YFSB1+YFSAB+YFSB2)
      ELSE
       FYFS=FEXPGL*(YFS0+YFSAB)
      ENDIF
*
      END

      FUNCTION AHARD(RR)
*     ======== =========
*
* Modified by M. Jack on 03/03/99 05:00pm
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /PCONST/ ALFAI,AL1PI,ALQE2,ALQF2,ALQEF,GMU,CSIGNB
      COMMON /BORN0I/ SBORN0,ABORN0,SBORNI,ABORNI
      COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /DELVAR/ DEL,DELM,DIF,DIF3,ALD,ALDM,ALDPL,ALDMI
     +               ,RD,RDI,RDI2,RDI3,RDM,RDMI,RDMI2,RDMI3
     +               ,ALRD,ALRDM
      COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
      COMMON /CORINT/ CORINT
*
      COMMON /SVAR  / S
      COMMON /MASSZ / AME,AMF,AME2,AMF2
*
      COMMON /INDEXA/ INDA,INDM
      COMMON /INDCUT/ INDC
      COMMON /INDEX/  MF,IIF,JF,KF
      COMMON /INDFIT/ IND,INDF
      COMMON /INTRFS/ INTRF
      COMMON /IMISDC/ IMISD
*-----------------------------------------------------------------------
*
      R=RR
      CALL SETR
*
*========================================================
* Choose different acceptance/acollinearity cut options:
*
      IF(INDC.LE.1) THEN
*
****************************************************
* 1. Acceptance cut and acollinearity cut:
****************************************************
*
        IF(R.GE.RECUT2) THEN
          ksoft=1
          CALL AHFULL (H3,H1)
          AHINI =H3
          AHINTF=H1
        ELSE
*       ENERGY CUTS
          ksoft=0
          DEL=(RECUT2-R)*R1MI
          CALL SETDEL
          CALL AHCUT (H3,H1)
          AHINI =H3
          AHINTF=H1
        ENDIF
*     
        IF(R.LT.RACUT) THEN
*       ACOLLINEARITY CUT
          ksoft=0
          DEL=1D0/2D0*(1D0-SQRT(ABS(1D0-SINAC2
     &    *(1D0+2D0*R+R2)*R1MI2)/COSAC2))
          CALL SETDEL
          CALL AHCUT (H3,H1)
          AHINI =AHINI -H3
          AHINTF=AHINTF-H1
        ENDIF
*
      ELSEIF (INDC.EQ.3) THEN
*    
*************************************************
* 3. Acceptance cut and acollinearity cut 
*    (new coding by M.Jack) :
*************************************************
*
        IF(R.GE.RECUT2) THEN
          ksoft=1
          CALL AHFULL (H3,H1)
          AHFIN = 0D0
cc          INDA=1
cc          INDM=0
cc          DEL=1D0/2D0*(1D0-SQRT(1D0-4D0*AMF2/S/R))
cc          CALL SETDEL
cc          CALL SAVAR
cc          CALL PHASEREGC
cc          CALL AHCUTACOL (H3,H1,H6)
          AHINI =H3
          AHINTF=H1
        ELSE
*       ENERGY CUTS
          ksoft=0
          INDA=0
          DEL=(RECUT2-R)*R1MI
          CALL SETDEL
          CALL SAVAR
          CALL PHASEREGC
*========================================================
cc          CALL AHCUT (H3,DUMMY)
cc          CALL AHCUTACOL (DUMMY,H1,H6)
*========================================================
          CALL AHCUTACOL(H3,H1,H6)
          AHINI =H3
          AHINTF=H1
          AHFIN =H6
        ENDIF
*     
        IF(R.LT.RACUT) THEN
*       ACOLLINEARITY CUT
          ksoft=0
          INDA=0
          DEL=1D0/2*(1-SQRT(ABS(1-SINAC2
     &       *(1+2*R+R2)*R1MI2)/COSAC2))
          CALL SETDEL
          CALL SAVAR
          CALL PHASEREGC
*========================================================
cc          CALL AHCUT (H3,DUMMY)
cc          CALL AHCUTACOL (DUMMY,H1,H6)
*========================================================
          CALL AHCUTACOL(H3,H1,H6)
          AHINI =AHINI -H3
          AHINTF=AHINTF-H1
          AHFIN =AHFIN -H6
        ENDIF
*
      ELSE
        WRITE(*,*) 'WRONG VALUE FOR INDEX INDC !'
      ENDIF
*
c      WRITE(*,*) 'INDC,INDA,R =',INDC,INDA,R
c      WRITE(*,*) 'H3,H1,H6 =',H3,H1,H6
c      WRITE(*,*) '--------------------------------------------'
*========================================================
*
      CALL BORN(IFINAL,R,R,DUMMY,ABORN,SBORNS,ABORNS)
*
      AHARD=ALQE2*AHINI*ABORN
*
      IF (INTERF.GE.1) THEN
        CALL BORN(-1,R,1D0,SBORN,DUMMY,SBORNS,DUMMYS)
        IF(INDF.EQ.10.AND.INTRF.EQ.2.AND.IMISD.EQ.1) THEN
          SBORN=SBORNS
        ENDIF
       AHARD=AHARD+ALQEF*(AHINTF*SBORN
     & -4d0/3*(-8*AL2-4*ALPL-ALMI*(C3+3*C)-C2)*R1MI*SBORNI*ksoft)
      ENDIF
*
      END

      SUBROUTINE AHFULL(H3,H1)
*     ========== =============
**  06/04/90 10:53pm
      IMPLICIT REAL*8(A-H,O-Z)
*
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM
      COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
*
* flags
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
      R1P=1D0+R
      R1PM2=1D0/R1P**2
      R1M=1D0-R
      R1M2=R1M*R1M
      R1M4=R1M2*R1M2
      R1PMI=R1P*R1MI
      R2P=1D0+R2
      R3P=1D0+R3
      R4P=1D0+R4
      TER=TE-ALR
*
      ALR1P=LOG(R1P)
      ALR2=ALR1P-AL2
*
      CC1=CM+R*CP
      CC12=CC1*CC1
      CC1M2=1D0/CC12
      ALC1=LOG(CC1)
*
      CC2=CP+R*CM
      CC22=CC2*CC2
      CC2M2=1D0/CC22
      ALC2=LOG(CC2)
*
      CC12M2=CC1M2*CC2M2
*
      H3=8*R2P*R1PM2*ALR2
     &+CC12M2*(TER*(4*R2P*R2*R1PM2-2*R2P**2*CPM
     &+8*R1M2*R1PM2*R2P*(R*CPM+R2P*CPM2))
     &+8*R1M2*(R*CPM-CPM2*R1P**2)+2*R1M4*CPM)
     &-4*R2P*CPM*(ALC1*CC1M2+ALC2*CC2M2)
      H3=R*H3
*
      IF (R.GE.RCUT) THEN
       H3=(H3-TE*2*C2)*R1MI
       IF (IPHOT2.GE.0) THEN
        IF(IFLAGS(IFFBHO).EQ.0) THEN
         H3=H3+AH2(R)*C2
        ELSE
         YY= - R1M/(4D0*R)*( 4D0*LOG(1D0-R) + 3D0 )
         H3=H3+(AH2(R)+YY)*4D0*R*R1PM2*C2
        ENDIF
       ENDIF
      ELSE
       H3=H3*R1MI
      END IF
*
      H1=0D0
      IF (INTERF.EQ.0) RETURN
       H1=C2*R*(2D0/3*R1M2-4D0/3*R1MI+4*R/R1P
     &+2D0/3*C2*R1M*R1P)/(CC1*CC2)
     &-ALR*2*C2*R1M*R
     &-ALPL*4*R1PMI*CPM*R2P
     &-16D0/3*R1MI*R3P*(CM3*ALCM+CP3*ALCP)
     &-ALC1*(8*C*CM2-4*CC1*(R1M+4D0/3*R1MI*CC12-4*CM2+CC1))
     &+ALC2*(8*C*CP2+4*CC2*(R1M+4D0/3*R1MI*CC22-4*CP2+CC2))
     &+ALR1P*(4D0/3*R-10D0/3*R2P)*R1PMI
*
      END
 
      SUBROUTINE AHCUT(H3,H1)
*     ========== ============
**  05/30/90 07:43pm
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CALL SETRC ( 1)
      CALL AHCUTM (H3M,H1M,I)
      H3=H3M
      H1=H1M
      CALL SETRC (-1)
      CALL AHCUTM (H3M,H1M,J)
      H3=H3+H3M
      H1=H1+H1M
      CALL SETRC ( 0)
      CALL AHCUTM (H3M,H1M,K)
      H3=H3-2*H3M
      H1=H1-2*H1M
*
      IF (I.NE.K) THEN
       CALL AJUMP (H3M,H1M,I,K)
       H3=H3+H3M
       H1=H1+H1M
      END IF
      IF (J.NE.K) THEN
       CALL AJUMP (H3M,H1M,J,K)
       H3=H3+H3M
       H1=H1+H1M
      END IF
*
      END
 
      SUBROUTINE AHCUTM(H3M,H1M,I)
*     ========== =================
**  05/31/90 02:10pm
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /DELVAR/ DEL,DELM,DIF,DIF3,ALD,ALDM,ALDPL,ALDMI
     +               ,RD,RDI,RDI2,RDI3,RDM,RDMI,RDMI2,RDMI3
     +               ,ALRD,ALRDM
      COMMON /RCVAR / C,C2,C3,CC,CCI,CCI2,CCI3,CCI4,CDEL,CDELM
     +                 ,ALCC,ALCD,ALCDM,ALMI,ALPL
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
*
      ALCDMI=ALCDM-ALCD
      ALCDPL=ALCDM+ALCD
*
      H1M=0D0
      IF(CDEL.GT.0.D0.AND.CDELM.GT.0.D0) THEN
* REGION ( +,+ )
       I=1
      IF (INTERF.GE.1) H1M=C2*((ALRD-ALRDM)*R+DIF*(R2-R)-ALDMI*R2)
      ELSE IF(CDELM.LT.0.D0.AND.CDEL.LT.0.D0) THEN
* REGION ( -,- )
       I=-1
      IF (INTERF.GE.1) H1M=C2*((ALRD-ALRDM)*R+DIF*(R-1)+ALDMI)
      ELSE
* REGION ( +,- )
       I=0
       ALPM=2*ALCC+ALCDPL-2*ALR-ALPL+TE+1
       H3M=2*ALPM*(-CCI*(R4+R3+R2+R)+CCI2*(R4+R2)+R3+R)
     + +2*ALCDPL*DEL*DELM*(R4-2*R3+2*R2-2*R+1)
     + +4*CCI*(R4+3*R3+3*R2+R)
     + -4*CCI2*(R4+R3+R2)
       H3M=H3M*R1MI3
      IF (INTERF.GE.1) THEN
       H1M=ALR*C2*(R4-2*R3+R2)
     +   +ALCC*(2D0/3*C3*(R4-R3-R+1)+C2*(-R4+2*R3-2*R+1)
     +   +2*C*(R4-3*R3+4*R2-3*R+1)+1D0/3*(-5*R4+2*R3-2*R+5))
     +   +ALMI*(1D0/3*C3+C)*(R4-R3+R-1)
     +   +ALPL*(1D0/2*C2*(R4-2*R3+2*R-1)+1D0/6*(5*R4-2*R3+2*R-5))
     +   +(ALRD-ALRDM)*C2*(R3-2*R2+R)
     +   +ALDM*C2*(-R4+2*R3-2*R+1)
     +   +DIF*C2*1D0/2*(R4-2*R3+2*R-1)
     +   +C2*1D0/6*(-3*R4+14*R3-14*R+3)
     +   +4*CCI*(R4+R2)*R1MI
       H1M=H1M*R1MI2
      END IF
       RETURN
      END IF
*
       H3M=I*(2*ALCDMI*DEL*DELM*(R4-2*R3+2*R2-2*R+1)
     . +2*ALCDMI*R*(R2+1)
     . -(2*CCI*ALCDMI*R)*(R3+R2+R+1)
     . +2*CCI2*ALCDMI*R2*(R2+1)
     . +2*CCI*DIF*R*(-R3+R2-R+1))
 
      H3M=H3M*R1MI3
*
      END
 
      SUBROUTINE AJUMP(H3M,H1M,I,J)
*     ========== ==================
**  05/30/90 07:43pm
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /DELVAR/ DEL,DELM,DIF,DIF3,ALD,ALDM,ALDPL,ALDMI
     +               ,RD,RDI,RDI2,RDI3,RDM,RDMI,RDMI2,RDMI3
     +               ,ALRD,ALRDM
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
*
      ALDPLE=ALDPL-TE-1
*
      H3MSUM=2*DIF*(R4+2*R3-2*R-1)
      H3M=H3MSUM*R1MI3*ISIGN(1,I-J)
*
      H1M=0D0
      R6 =0D0
      FM =0D0
      FM2=0D0
      FM3=0D0
      FP =0D0
      FP2=0D0
      FP3=0D0
      IF (INTERF.GE.1) THEN
       R6=R*R5
       T0=ALRDM-ALRD
       T1=ALRDM*RDMI-ALRD*RDI
       T2=ALRDM*RDMI2-ALRD*RDI2
       T3=ALRDM*RDMI3-ALRD*RDI3
       FM =RDMI-RDI
       FM2=RDMI2-RDI2
       FM3=RDMI3-RDI3
       FP =RDMI+RDI
       FP2=RDMI2+RDI2
       FP3=RDMI3+RDI3
       H1MSUM=4*(T0*1D0/3*(-2*R6+3*R5-6*R4+6*R2-3*R+2)
     + +2*T1*(R6-R5+4*R4-4*R3+R2-R)+2*T2*(-R6-R5+R3+R2)
     + +4D0/3*T3*(R6-R3))
     + +ALR*(+2*FM*(-7*R6+5*R5-18*R4+14*R3-3*R2+R)
     + +2*FM2*(7*R6+4*R5+2*R4-4*R3-R2)
     + +8D0/3*FM3*(-3*R6+R3))
     + +ALDMI*(1D0/3*(11*R6-3*R4+16*R3-3*R2+11)
     + +2*FP*(-3*R6+R5-2*R4-2*R3+R2-3*R)
     + +2*FP2*(3*R6+2*R4+3*R2)-8D0/3*FP3*(R6+R3))
     + +4D0/3*(FM*(-7*R5+5*R4-5*R3+7*R2)
     + +4*FM2*(R5-R3))+(RDM-RD)*(-R5-3*R4+6*R3-6*R2+3*R+1)
       H1M=H1MSUM*R1MI4*ISIGN(1,I-J)
      END IF
*
      IF (ABS(I-J).EQ.2) RETURN
*
      H3MDIF=4*(ALDPLE+ALR)*DEL*DELM*(R4-2*R3+2*R2-2*R+1)
     + +8*DEL*DELM*(R3-2*R2+R)+2*(R4+4*R3+10*R2+4*R+1)
      H3MDIF= H3MDIF*R1MI3*ISIGN(1,I-J)
*
      H1MDIF=0D0
      IF (INTERF.GE.1) THEN
       T0=ALRDM+ALRD
       T1=ALRDM*RDMI+ALRD*RDI
       T2=ALRDM*RDMI2+ALRD*RDI2
       T3=ALRDM*RDMI3+ALRD*RDI3
       H1MDIF=ALD*2D0/3*(5*R6-12*R5+9*R4-9*R2+12*R-5)
     + +4*(T0*1D0/3*(-2*R6+3*R5-6*R4+6*R2-3*R+2)
     + +2*T1*(R6-R5+4*R4-4*R3+R2-R)+2*T2*(-R6-R5+R3+R2)
     + +4D0/3*T3*(R6-R3))
     + +ALR*(2D0/3*(11*R6-6*R5+18*R4+8*R3-21*R2+6*R)
     + +2*FP*(-7*R6+5*R5-18*R4+14*R3-3*R2+R)
     + +2*FP2*(7*R6+4*R5+2*R4-4*R3-R2)+8D0/3*FP3*(-3*R6+R3))
     + +ALDMI*(1D0/3*(5*R6-12*R5+9*R4-9*R2+12*R-5)
     + +2*FM*(-3*R6+R5-2*R4-2*R3+R2-3*R)+2*FM2*(3*R6+2*R4+3*R2)
     + -8D0/3*FM3*(R6+R3))+4D0/3*(FP*(-7*R5+5*R4-5*R3+7*R2)+4*FP2
     + *(R5-R3)+8*(R5-R4+R2-R))+(RDM+RD)*(-R5-3*R4+6*R3-6*R2+3*R+1)
       H1MDIF= H1MDIF*R1MI4*ISIGN(1,I-J)
      END IF
*
      IF ((I+J).EQ. 1) THEN
       H3M=(H3M+H3MDIF)/2
      IF (INTERF.GE.1) H1M=(H1M+H1MDIF)/2
      ELSE
       H3M=(H3M-H3MDIF)/2
      IF (INTERF.GE.1) H1M=(H1M-H1MDIF)/2
      END IF
*
      END
 
      FUNCTION AH2(R)
*     ======== ======
*-----------------------------------------------------------------------
*     ALPHA**2*L**2 AND ALPHA**2*L PHOTONIC CORRECTIONS
*     NEW: ALPHA**2*L**2 PAIR CORRECTIONS (ABA 21/09/99)
*-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
*
* flags
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
      COMMON /PCONST/ ALFAI,AL1PI,ALQE2,ALQF2,ALQEF,GMU,CSIGNB
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
      COMMON /CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
      COMMON /CPAIRS/ FPI2,RINF,RRR0,RRR1,PAIRDL,CORFAC
*
      AH2=0D0
      R1M=1D0-R
      R2 =R**2
      IF(R1M.GT.1.D-14) THEN
       R1MI=1D0/R1M
       R1P=1D0+R
       R2P=1D0+R2
       ALR1M=LOG(R1M)
       D2R1M=DDILOG(R1M)
       SQR=SQRT(R)
       ALR=LOG(R)
       TL=TE+1D0
       IF(IPHOT2.GE.0) THEN
        AH2=AH2+ALQE2*
     +  TL**2*(-R2P*R1MI*ALR+R1P*(1D0/2D0*(ALR-1D0)-2D0*ALR1M)-2D0)
     +  +ALQE2*1D0/4*TL**2*(R1M**3/(2*R)+R1M**2/SQR*(ATAN(1D0/SQR)
     +  -ATAN(SQR))-R1P*ALR+2*R1M)
       ENDIF
       IF(IPHOT2.GE.1) THEN
        AH2=AH2+ALQE2*(
     +  +TL  *( R2P*R1MI*(D2R1M+ALR*(ALR1M+7D0/2D0-1D0/2D0*ALR))
     +        +R1P*(1D0/4D0*ALR**2+4D0*ALR1M-2D0*F1)-ALR+7D0+1D0/2D0*R))
       ENDIF
       IF(IPHOT2.GE.2) THEN
        D3R1M=TRILOG(R1M)
        S1R1M=   S12(R1M)
        AH2=AH2+
     +         ALQE2*(-2D0*R1P*(ALR1M-F1+1D0)+R2P*R1MI*(-1D0/6D0*ALR**3
     +         +.5D0*ALR*D2R1M+.5D0*ALR**2*ALR1M-1.5D0*D2R1M
     +         -1.5D0*ALR*ALR1M+F1*ALR-17D0/6D0*ALR-ALR**2)
     +         +R1P*(1.5D0*D3R1M-2D0*S1R1M-ALR1M*D2R1M-.5D0)
     +         -.25D0*(1D0-5D0*R)*ALR1M**2+.5D0*(1D0-7D0*R)*ALR*ALR1M
     +         -25D0/6D0*R*D2R1M+(-1D0+13D0/3D0*R)*F1+(1.5D0-R)*ALR1M
     +         +1D0/6D0*(11D0+10D0*R)*ALR+2D0*R1MI**2*ALR**2
     +         -25D0/11D0*R*ALR**2-2D0/3D0*R*R1MI*(1D0+R1MI*ALR)**2 )
       ENDIF
      ENDIF
*
      IF(IFLAGS(IFISPP).GE.2.AND.IFLAGS(IFFBHO).NE.0) THEN
       AH2PR = 0D0
       TPI = DLOG(SCOM/FPI2)
       YY1 = (2D0*R**2 - R - 1D0)/3D0
       IF(IFLAGS(IFIPFC).EQ.1) AH2PR = - ALQE2*TEE**2*YY1
       IF(IFLAGS(IFIPFC).EQ.2) AH2PR = - ALQE2*TMU**2*YY1
       IF(IFLAGS(IFIPFC).EQ.3) AH2PR = - ALQE2*TTU**2*YY1
       IF(IFLAGS(IFIPFC).EQ.4) AH2PR = - ALQE2*TPI**2*YY1*RINF
       IF(IFLAGS(IFIPFC).EQ.5) AH2PR = - ALQE2*TEE**2*YY1
     &                                 - ALQE2*TMU**2*YY1
     &                                 - ALQE2*TTU**2*YY1
     &                                 - ALQE2*TPI**2*YY1*RINF
       IF(IFLAGS(IFIPFC).EQ.6) AH2PR = - ALQE2*TEE**2*YY1
     &                                 - ALQE2*TMU**2*YY1
     &                                 - ALQE2*TTU**2*YY1
       IF(IFLAGS(IFIPSC).GE.1) AH2PR = AH2PR + ALQE2/4D0*TEE**2
     &   *( R1M/(3D0*R)*(4D0 + 7D0*R + 4D0*R**2) + 2D0*R1P*ALR )
       AH2 = AH2 + AH2PR
      ENDIF
      END
 
      SUBROUTINE SETR
*     ========== ====
      IMPLICIT DOUBLE PRECISION (A-H,O-W,Y-Z)
*
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
*
      R2=R*R
      R3=R2*R
      R4=R3*R
      R5=R4*R
      R1MI=1D0/(1-R)
      R1MI2=R1MI*R1MI
      R1MI3=R1MI2*R1MI
      R1MI4=R1MI3*R1MI
      ALR=LOG(R)
*
      END
 
      SUBROUTINE SETDEL
*     ========== ======
      IMPLICIT DOUBLE PRECISION (A-H,O-W,Y-Z)
*
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /DELVAR/ DEL,DELM,DIF,DIF3,ALD,ALDM,ALDPL,ALDMI
     +               ,RD,RDI,RDI2,RDI3,RDM,RDMI,RDMI2,RDMI3
     +               ,ALRD,ALRDM
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
*
      DELM=1-DEL
      DIF=DELM-DEL
      DIF3=DELM**3-DEL**3
      ALD =LOG(DEL)
      ALDM=LOG(DELM)
      ALDPL=ALDM+ALD
      ALDMI=ALDM-ALD
*
      IF (INTERF.GE.1) THEN
       RD=R*DEL+DELM
       RDI =1D0/RD
       RDI2=RDI*RDI
       RDI3=RDI*RDI2
       ALRD=LOG(RD)
       RDM=R*DELM+DEL
       RDMI =1D0/RDM
       RDMI2=RDMI*RDMI
       RDMI3=RDMI*RDMI2
       ALRDM=LOG(RDM)
      END IF
*
      END
 
      SUBROUTINE SETRC(J)
*     ========== ========
      IMPLICIT DOUBLE PRECISION (A-H,O-W,Y-Z)
*
      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /DELVAR/ DEL,DELM,DIF,DIF3,ALD,ALDM,ALDPL,ALDMI
     +               ,RD,RDI,RDI2,RDI3,RDM,RDMI,RDMI2,RDMI3
     +               ,ALRD,ALRDM
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM
      COMMON /RCVAR / CT,CT2,CT3,CC,CCI,CCI2,CCI3,CCI4,CDEL,CDELM
     +                  ,ALCC,ALCD,ALCDM,ALMIT,ALPLT
*
      IF (J.EQ.1) THEN
       CT= C
       CT2= C2
       CT3= C3
       CC=CP+R*CM
       CDEL =CC*DEL -CM*R
       CDELM=CC*DELM-CM*R
       ALPLT= ALPL
       ALMIT= ALMI
      ELSE IF (J.EQ.-1) THEN
       CT=-C
       CT2= C2
       CT3=-C3
       CC=CM+R*CP
       CDEL =CC*DEL -CP*R
       CDELM=CC*DELM-CP*R
       ALPLT= ALPL
       ALMIT=-ALMI
      ELSE
       CT=0D0
       CT2= 0D0
       CT3= 0D0
       CC=(1D0+R)/2
       CDEL =CC*DEL -1D0/2*R
       CDELM=CC*DELM-1D0/2*R
       ALPLT= -2*AL2
       ALMIT= 0D0
      END IF
*
      CCI =1D0/CC
      CCI2=CCI*CCI
      CCI3=CCI*CCI2
      CCI4=CCI*CCI3
      ALCC=-LOG(CCI)
      ALMIN=1d0/2*(ALPLT-TE-1)+AL2+ALR-ALCC
      ALCD =MAX( ALMIN, LOG(ABS(CDEL)) )
      ALCDM=MAX( ALMIN, LOG(ABS(CDELM)) )
*
      END
 
      SUBROUTINE SETFIN
*     ========== ======
      IMPLICIT REAL*8(A-H,O-Z)
*
      PARAMETER (NP=20)
      COMMON /FORINT/ ALE1(NP),ALE2(NP),ALA1(NP),ALA2(NP),NPOINT
      COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT
      COMMON /FORFAL/ KREG
*
      EXTERNAL FAL1,FAL2
*
      NPOINT=NP
      REPS1 = 1D-5
      AEPS1 = 1D-5
      KREG=2
      ASTEP=(RACUT-RECUT1)/(NP-1)
      S1=0D0
      S2=0D0
      H=ASTEP/5
      ALA1(NP)=0D0
      ALA2(NP)=0D0
      A=RACUT
      DO 100 I=1,NP-1
        B=A
        A=A-ASTEP
        CALL SIMPS(A,B,H,REPS1,AEPS1,FAL1,X,SI,AIHH,AIHA)
        S1=S1+SI
        ALA1(NP-I)=S1
        CALL SIMPS(A,B,H,REPS1,AEPS1,FAL2,X,SI,AIHH,AIHA)
        S2=S2+SI
        ALA2(NP-I)=S2
100   CONTINUE
*
      KREG=3
      ASTEP=(RECUT2-RECUT1)/(NP-1)
      S1=0D0
      S2=0D0
      ALE1(NP)=0D0
      ALE2(NP)=0D0
      H=ASTEP/5
      A=RECUT2
      DO 200 I=1,NP-1
        B=A
        A=A-ASTEP
        CALL SIMPS(A,B,H,REPS1,AEPS1,FAL1,X,SI,AIHH,AIHA)
        S1=S1+SI
        ALE1(NP-I)=S1
        CALL SIMPS(A,B,H,REPS1,AEPS1,FAL2,X,SI,AIHH,AIHA)
        S2=S2+SI
        ALE2(NP-I)=S2
200   CONTINUE
*
      END
 
      FUNCTION FAL1(Z)
*     ======== =======
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT
      COMMON /FORFAL/ KREG
*
      Z1MI=1D0/(1-Z)
      IF(KREG.EQ.2) THEN
       RATZ2=((1+Z)*Z1MI)**2
       ARG=SQRT(MAX(0D0,((1-SINAC2*RATZ2)/COSAC2)))
      ELSE
       ARG=(1-2*RECUT2+Z)*Z1MI
      END IF
      FAL1=-(1+Z**2)*Z1MI*LOG( MAX(1.D-32,ABS((1-ARG)/(1+ARG))) )
*
      END
 
      FUNCTION FAL2(Z)
*     ======== =======
* T. RIEMANN, 04/03/1999 17:00:
* A TYPO CORRECTED
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT
      COMMON /FORFAL/ KREG
*
      COMMON /INDFIN/ IFUNFIN
*
      Z1MI=1D0/(1-Z)
* T. RIEMANN 04/03/99 ADDED:
      RATZ1=(1+Z)*Z1MI
****************************
      RATZ2=((1+Z)*Z1MI)**2
      IF(KREG.EQ.2) THEN
       ARG=SQRT(MAX(0D0,((1-SINAC2*RATZ2)/COSAC2)))
      ELSE
       ARG=(1-2*RECUT2+Z)*Z1MI
      END IF
*
* M. JACK 18/03/99 MODIFIED :
      IF (IFUNFIN.EQ.0) THEN
        FAL2=(1+Z)*LOG( MAX(1.D-32,ABS((RATZ2-ARG)/(RATZ2+ARG))) )
      ELSE
* T. RIEMANN 04/03/99 CORRECTED (see notes 27/01/99-1):
        FAL2=(1+Z)*LOG( MAX(1.D-32,ABS((RATZ1-ARG)/(RATZ1+ARG))) )
cc M. JACK 10/03/99 MODIFIED
cc       Z1PI=1D0/(1D0+Z)       
cc       ARG0 = (1D0-Z)*Z1PI
cc       FAL2=(1+Z)*LOG( MAX(1.D-32,ABS((1D0-ARG*ARG0)/(1D0+ARG*ARG0))) )
      ENDIF
*
      END

      SUBROUTINE FUNFIN(S,AMF2,QF2,R,SFIN,AFIN)
*     ========== ==============================
* T. RIEMANN, 04/03/1999 17:00:
* A TYPO CORRECTED
*-----------------------------------------------------------------------
* ROUTINE CALCULATES FINAL COS-EVEN AND COS-ODD FACTORS
*-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /PCONST/ ALFAI,AL1PI,ALQE2,ALQF2,ALQEF,GMU,CSIGNB
      COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT
*
      PARAMETER (NP=20)
      COMMON /FORINT/ ALE1(NP),ALE2(NP),ALA1(NP),ALA2(NP),NPOINT
      COMMON /CALQED/ ALQEDZ,ALQEDS
*
      COMMON /INDFIN/ IFUNFIN
*
      ALQF2H=AL1PI*QF2
      IF(IFLAGS(IFFSRS).EQ.0) THEN
        CORFAC=1d0
      ELSE
        CORFAC=ALQEDS*ALFAI
      ENDIF
      ALQF2Z=ALQF2H*CORFAC
*
      TF=LOG(S/AMF2)-1
      BETTAF=TF*2D0*ALQF2H
      SOFTFR=3D0/4*BETTAF+ALQF2H*(2D0*F1-1D0/2)
*
      ALIM=RECUT1/R
      IF(ALIM.GE.RECUT2) THEN
        A=ALIM                     
        SFIN=0D0
        AFIN=0D0
      ELSE
        CALL INTERP(NPOINT,RECUT1,RECUT2,ALE1,ALIM,ALE1I)
        CALL INTERP(NPOINT,RECUT1,RECUT2,ALE2,ALIM,ALE2I)
        D1=ALIM-RECUT2
        D2=RECUT2-1D0
*
* M. JACK 18/03/99 MODIFIED :
        IF (IFUNFIN.EQ.0) THEN
          SFIN=ALQF2Z*(ALE1I+D1*(D1/2+D2))
        ELSE 
* T. RIEMANN 04/03/99 CORRECTED (see notes 03/02/99-1):
          SFIN=ALQF2Z*(ALE1I+D1*(D1/2-D2))
        ENDIF
*
        AFIN=ALQF2Z*(ALE1I+ALE2I)
        A=RECUT2
      ENDIF
*
      ALR=LOG(R)
      ADD=ALQF2Z*ALR
      BETF=BETTAF*CORFAC+2*ADD
      AM1=1D0-A
      A2=A*A
      ALNA=LOG(A)
      IF(AM1.LE.1D-16) THEN
        FINS=0D0
      ELSE
        SFT=1D0+SOFTFR*CORFAC+3D0/2*ADD
        FINS=AM1**BETF*SFT
      ENDIF
*
      FIN=FINS+BETF*(A2+2D0*A-3D0)/4-ALQF2Z*2D0*DDILOG(AM1)
      SFIN=SFIN+FIN+ALQF2Z*(ALNA*(A2/2+A)-1D0/4*A2-A+5D0/4)
      AFIN=AFIN+FIN+ALQF2Z*(A2/2-A+1D0/2)
*
      IF(ALIM.LT.RACUT) THEN
        CALL INTERP(NPOINT,RECUT1,RACUT,ALA1,ALIM,ALA1I)
        CALL INTERP(NPOINT,RECUT1,RACUT,ALA2,ALIM,ALA2I)
        B=-(1D0+SINAC2)/COSAC2
        SQR=SQRT(MAX(0D0,(ALIM**2+2D0*B*ALIM+1D0)))
        SFIN=SFIN-ALQF2Z*(ALA1I+1D0/2
     &      *((ALIM+B)*SQR-(1-B*B)*LOG(ABS((RACUT+B)/(SQR+ALIM+B)))) )
        AFIN=AFIN-ALQF2Z*(ALA1I+ALA2I)
      ENDIF
*
      END
 
      SUBROUTINE INTERP(NPOINT,XMIN,XMAX,Y,XI,YI)
*     ========== ================================
*-----------------------------------------------------------------------
*     ROUTINE MAKES QUADRATIC INTERPOLATION FOR ARRAY "Y" OF "NPOINT"
*     EQUIDISTANT POINTS IN THE INTERVAL (XMIN,XMAX)
*-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
*
      PARAMETER (NR=20)
      DIMENSION Y(NR)
*
      IF(NR.LT.NPOINT) STOP 'WRONG NUMBER OF POINTS'
      XSTEP=(XMAX-XMIN)/(NPOINT-1)
      N=INT((XI-XMIN)/XSTEP+0.5D0)+1
      IF(N.LT.2) THEN
       N=2
      ELSE IF(N.GT.NPOINT-1) THEN
       N=NPOINT-1
      END IF
*
      XCLOSE=XMIN+XSTEP*(N-1)
      F1=(XI-XCLOSE+XSTEP)
      F2=(XI-XCLOSE)
      F3=(XI-XCLOSE-XSTEP)
*
      YI=(0.5D0*F2*F3*Y(N-1)-F3*F1*Y(N)+0.5D0*F1*F2*Y(N+1))/XSTEP**2
      END
 
      SUBROUTINE BOXINT(BOXIS,BOXIA)
*     ========== ===================
      IMPLICIT COMPLEX*16(X)
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
*
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /PSCONS/ SW2,AMZ,GAMZ
      COMMON /PCONST/ ALFAI,AL1PI,ALQE2,ALQF2,ALQEF,GMU,CSIGNB
      COMMON /FORCHI/ XKAPP,XKAPPC,XMZ2,XMZ2C
      COMMON /SVAR  / S
      COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM
      COMMON /COUPL / VEFA,XVEFI,VEFZ,AEFA,XAEFI,AEFZ,VEEZ,XVPOL,VPOL2
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
*
      XR=XMZ2/S
      XR2=XR*XR
      XR3=XR2*XR
      XLR=LOG(XR)
      XL1R=LOG(1-1D0/XR)
      XDLR=XSPENZ(1-1D0/(2*XR))
      XDL1R=XSPENZ(1-1D0/XR)
      XDLCPR=XSPENZ(1-CP/XR)
      XDLCMR=XSPENZ(1-CM/XR)
      XDLPLR=XDLCPR+XDLCMR
      XDLMIR=XDLCPR-XDLCMR
*
      REF4BX=-ALPL*ALMI*1D0/2*(C2-1)-2*ALPL*C+ALMI*(C2-3)+6*C
      AMF4BX=ALMI*PI*(C2-1)
      XF4BX=DCMPLX ( REF4BX, AMF4BX )
*
      XH4BX=4*XL1R*ALPL*C*(XR2-XR)
     & -2*XL1R*ALMI*(-2*XR2+XR*(C2+1)+C2-1)
     & -4*XL1R*C*(XR2+XR)-8*XDL1R*C*(XR2-XR)
     & +4*XLR*C-ALPL*ALMI*(C2-1)
     & -4*ALPL*C-2*ALMI*(XR*(C2-1)-C2+3)
     & +4*XDLPLR*C*(XR2-XR)+2*XDLMIR*(2*XR2-XR*(C2+1))
     & -4*XR*C+12*C
*
      XG4BX=1D0/2*(DCONJG(XF4BX)+XH4BX)
*
      XCHI=XKAPP*S/(S-XMZ2)
      CHI2=DREAL(XCHI)**2+DIMAG(XCHI)**2
      BOXIS=-ALQEF*
     & (AEFA*VPOL2*REF4BX
     &+2*DREAL(XAEFI*DCONJG(XVPOL)*XCHI*XG4BX)+AEFZ*CHI2*DREAL(XH4BX))
*
      IF ((IAFB.EQ.0).AND.(ISYM.EQ.1).AND.INTERF.EQ.0) RETURN
*
      REF1BX=-(ALCP**2+ALCM**2)*1D0/2*(C2-1)
     & -AL2**2-6*AL2+ALPL*(C2-3)-2*ALMI*C-C2
      AMF1BX=-PI*(10D0/3*AL2+ALPL*(C2+5D0/3)+ALMI*2*COPL3
     & -4D0/3*C2)
      XF1BX=DCMPLX ( REF1BX, AMF1BX )
*
      XH1BX=ALPL*ALMI*COPL3+4D0/3*(ALCP**2+ALCM**2)
     & +XL1R*AL2*(8*XR2-4*XR+20D0/3)
     & +XL1R*ALPL*(4*XR2-2*XR*(C2+1)+2*C2+10D0/3)
     & +XL1R*ALMI*4d0*(C*(XR2-XR)+COPL3)
     & +2*XL1R*C2*(-XR2+3*XR-4D0/3)+4*XDL1R*C2*(XR-1)
     & +XLR*C2*(8D0/3*XR-10D0/3)+XDLR*(32D0/3*XR3-8*XR2+4*XR-4D0/3)
     & -8D0/3*AL2**2+AL2*(16D0/3*XR2+16D0/3*XR-12)
     & +ALPL*(8D0/3*XR2+XR*(-4D0/3*C2+8D0/3)+2*(C2-3))
     & +ALMI*C*(8D0/3*XR2+4D0/3*XR-4)
     & +XDLPLR*(-16D0/3*XR3+4*XR2-2*XR*(C2+1)+2*C2+2D0/3)
     & +XDLMIR*(4*C*(XR2-XR)+2*COPL3)
     & +(4D0/3*XR-2)*C2
*
      XG1BX=1D0/2*(DCONJG(XF1BX)+XH1BX)
*
      BOXIA=ALQEF*(VEFA*VPOL2*REF1BX
     & +2*DREAL(XVEFI*DCONJG(XVPOL)*XCHI*XG1BX)+VEFZ*CHI2*DREAL(XH1BX))
*
      END
 
      SUBROUTINE SFTINT(SFTIS,SFTIA)
*     ========== ===================
      IMPLICIT REAL*8(A-H,O-Z)
*
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /SVAR  / S
      COMMON /BORN0I/ SBORN0,ABORN0,SBORNI,ABORNI
      COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM
      COMMON /PCONST/ ALFAI,AL1PI,ALQE2,ALQF2,ALQEF,GMU,CSIGNB
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
      COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT
*
      ALDEL=LOG(1-RCUT)
*
      SFTI4=4*ALDEL*(ALMI*(C2-1)+2*C)
     & -2*(DLCP-DLCM)*(C2-1)+ALPL*ALMI*(C2-1)+4*ALPL*C+4*ALMI-8*C
      SFTIS=-SFTI4*ABORNI*ALQEF
*
      IF ((IAFB.EQ.0).AND.(ISYM.EQ.1).AND.INTERF.EQ.0) RETURN
*
      SFTI1=4D0/3*ALDEL*(-8*AL2-4*ALPL-ALMI*3*COPL3-C2)
     + +(DLCP-DLCM)*2*COPL3+(DLCP+DLCM)*8D0/3
     + -ALPL*ALMI*COPL3-(ALCM**2+ALCP**2)*4D0/3
     + +16D0/3*AL2**2+4D0/3*AL2-2D0/3*ALPL*(C2-1)
     + +2D0/3*C2-8D0/3*F1
      SFTIA=SFTI1*SBORNI*ALQEF
*
      END
 
      FUNCTION SIGINT(NCOS,SIGC,CC1,CC2)
*     ======== =========================
*-----------------------------------------------------------------------
*      ROUTINE USES "NCOS" VALUES OF FUNCTION "SIGC"
*      IN EQUIDISTANT POINTS FROM -1 TO 1 FOR QUADRATIC INTERPOLATION
*      AND RETURNS AN INTEGRAL FROM  CC1  TO  CC2
*      NCOS > 2
*-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      DIMENSION SIGC(*)
      PART(N,D)=2*D*(2*D**2*(SIGC(N-1)-2*SIGC(N)+SIGC(N+1))
     + -3*D*(SIGC(N-1)-SIGC(N+1))+12*SIGC(N))
*
      SIGINT=0.D0
      RNCOS2=NCOS/2
      IF(NCOS.LT.3) RETURN
      SCALE=2.D0/24/NCOS
*
      IF(CC2.GE.CC1) THEN
        IDER=1
        C1=CC1
        C2=CC2
      ELSE
        IDER=-1
        C2=CC1
        C1=CC2
      END IF
      N1=(C1+1)*RNCOS2+1
      N2=(C2+1)*RNCOS2+1
      IF(N1.EQ.1) N1=2
      IF(N1.GE.NCOS) N1=NCOS-1
      D1=(C1+1)*RNCOS2-N1+0.5D0
      P1=PART(N1,D1)-PART(N1,-0.5D0)
      IF(N2.EQ.1) N2=2
      IF(N2.GE.NCOS) N2=NCOS-1
      D2=(C2+1)*RNCOS2-N2+0.5D0
      P2=PART(N2,D2)-PART(N2, 0.5D0)
*
      DO 1 K=N1,N2
      SIGINT=SIGINT+SIGC(K-1)+22*SIGC(K)+SIGC(K+1)
1     CONTINUE
      SIGINT=IDER*(SIGINT-P1+P2)*SCALE
*
      END
*
* Below are differential functions
*
      SUBROUTINE SOFTIN(SFTIN)
*     ========== =============
**  06/08/90 12:09pm
      IMPLICIT REAL*8(A-H,O-Z)
*
      COMMON /BORN0I/ SBORN0,ABORN0,SBORNI,ABORNI
      COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM
      COMMON /PCONST/ ALFAI,AL1PI,ALQE2,ALQF2,ALQEF,GMU,CSIGNB
      COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT
*
      ALE=LOG(1-RCUT)
*
      SFTI=2*(2*ALE*(ALCM-ALCP)+DLCP-DLCM-1D0/2*(ALCP**2-ALCM**2))
      SFTIN=((1+C2)*SBORNI+2*C*ABORNI)*SFTI*ALQEF
*
      END
 
      FUNCTION SOFT(X)
*     ======== =======
      IMPLICIT REAL*8 (A-H,O-Z)
*
      COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
      COMMON /INDFOS/ INDF
*
      R=1D0-X
      IF(INDF.NE.-1) THEN
        CALL BORN(IFINAL,R,R,SBORN,ABORN,SBORNS,ABORNS)
        IF(IPHOT2.EQ.5) THEN
          YFS=FYFS(X)
        ELSE
          YFS=1D0
        ENDIF
        SOFT=(SBORN*(1D0+C2)+ABORN*2D0*C)*YFS
      ELSE
        IF(IPHOT2.EQ.5) THEN
          YFS=FYFS(X)
        ELSE
          PRINT *,'Wrong setting of FOT2 flag for INDF=-1'
          STOP 
        ENDIF        
        CALL BORNN(R,R,C,SBORNN)
        SOFT=SBORNN*YFS
      ENDIF
*
      END
 
      FUNCTION HARD(RR)
*     ======== ========
*======================================================================*
*         Modified by M. JACK on 26/06/99 07:00pm                      *
*======================================================================*
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      COMMON /PCONST/ ALFAI,AL1PI,ALQE2,ALQF2,ALQEF,GMU,CSIGNB
      COMMON /BORN0I/ SBORN0,ABORN0,SBORNI,ABORNI
      COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /DELVAR/ DEL,DELM,DIF,DIF3,ALD,ALDM,ALDPL,ALDMI
     +               ,RD,RDI,RDI2,RDI3,RDM,RDMI,RDMI2,RDMI3
     +               ,ALRD,ALRDM
      COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
      COMMON /CORINT/ CORINT
      COMMON /SVAR  / S
      COMMON /MASSZ / AME,AMF,AME2,AMF2
      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
      COMMON /INDCUT/ INDC
      COMMON /INDEXA/ INDA,INDM
      COMMON /INDEX/  MF,IIF,JF,KF
      COMMON /INDFIT/ IND,INDF
      COMMON /INTRFS/ INTRF
      COMMON /IMISDC/ IMISD
*
* flags
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
*-----------------------------------------------------------------------
*
      R=RR
      CALL SETR
*========================================================
* Choose different acceptance/acollinearity cut options:
*
      IF(INDC.LE.0) THEN
*
*************************************************
* 1. Acollinearity cut (ICUT=0) 
*    (original coding by M.Bilenky) :
*************************************************
*
      IF (R.GE.RECUT2) THEN
        ksoft=1
        DEL=1D0/2D0*(1D0-DSQRT(1D0-4D0*AMF2/(S*R)))
C--------------------------------------------------------------- 
CBARDIN!!!
C 03/02/1997 db detected a problem in MB's logics: 
C everything gets destroyed in the angular distribution 
C if AMF2 is large enough (tau, c, b).
C I checked, that for large AMF2, H0 and H3 returned by HCUT
C do not compensate each other as they do for small AMF2.
C This is due to the fact that two CALLs of HINIM followed by
C two CALLs of SETRC in HCUT returns NON-compensating
C contributions for H0 and H3. 
C Db's solution of this problem is the following IF,
C which passes the code to (-,+) - region of HINIM 
C (see comments therein).
C--------------------------------------------------------------- 
        IF(IRCUT.NE.0) DEL=1D-14
        CALL SETDEL
        CALL HCUT(H0,H4,H3,H1)
        IF (R.GE.RCUT) THEN
          H0=H0-2*TE*COPL2*R1MI
          H3=H3-TE*4D0*C*R1MI
          IF (IPHOT2.GE.0) THEN
           H0=H0+SH2(R,2)*COPL2
           IF(IFLAGS(IFFBHO).EQ.0) THEN
            H3=H3+AH2(R)*2D0*C
           ELSE
            YY= - (1D0-R)/(4D0*R)*( 4D0*LOG(1D0-R) + 3D0 )
            H3=H3+(AH2(R)+YY)*8D0*R/(1D0+R)**2*C
           ENDIF
          END IF
        END IF
*
        SHINI =H0
        SHINTF=H4
        AHINI =H3
        AHINTF=H1
      ELSE
* ENERGY CUTS
        ksoft=0
        DEL=(RECUT2-R)*R1MI
        CALL SETDEL
        CALL HCUT(H0,H4,H3,H1)
        SHINI =H0
        SHINTF=H4
        AHINI =H3
        AHINTF=H1
      END IF
*
      IF (R.LT.RACUT) THEN
* ACOLLINEARITY CUT
        ksoft=0
        DEL=1D0/2*(1-SQRT(ABS(1-SINAC2*(1+2*R+R2)*R1MI2)/COSAC2))
        CALL SETDEL
        CALL HCUT(H0,H4,H3,H1)
        SHINI =SHINI-H0
        SHINTF=SHINTF-H4
        AHINI =AHINI-H3
        AHINTF=AHINTF-H1
      ENDIF
*
*------------------------------------------------------------------
*
      ELSE
*
*************************************************
* 2. Acollinearity cut (ICUT=2,3) or s' cut (ICUT=1) 
*    (new coding by M.Jack) :
*************************************************
*
      IF (R.GE.RECUT2) THEN
        ksoft=1
        DEL=1D0/2D0*(1D0-DSQRT(1D0-4D0*AMF2/(S*R)))
C---------------------------------------------------------------  
C See comments above for INDC<=0 by BARDIN from 03/02/1997 !!! 
C--------------------------------------------------------------- 
C        IF(IRCUT.NE.0) DEL=1D-14
C============================================== 
C Modification by M.Jack 26/06/1999:
*
        IF((IRCUT.NE.0).AND.(IRCUT.NE.2)
     &  .AND.(IRCUT.NE.3)) DEL=1D-14
C============================================== 
*
        CALL SETDEL
*
        IF (INDC.EQ.2) THEN
*
        INDA=1
        CALL SAVAR        
        CALL HCUTC(H0,H4,H3,H1)
*
        ELSE 
*
        INDA=1
        CALL SAVAR
        CALL PHASEREGC
        CALL HCUTACOL(H0,H4,H3,H1)
*
cc        IF ((IIF.NE.0).OR.(JF.NE.0)) THEN
cc         WRITE(*,*) 'INDC = ',INDC
cc         WRITE(*,*) 'R,C = ',R,C
cc         WRITE(*,*) 'I,J = ',IIF,JF
cc        ENDIF
*
        ENDIF
*
        IF (R.GE.RCUT) THEN
          H0=H0-2*TE*COPL2*R1MI
          H3=H3-TE*4D0*C*R1MI
          IF (IPHOT2.GE.0) THEN
           H0=H0+SH2(R,2)*COPL2
           IF(IFLAGS(IFFBHO).EQ.0) THEN
            H3=H3+AH2(R)*2D0*C
           ELSE
            YY= - (1D0-R)/(4D0*R)*( 4D0*LOG(1D0-R) + 3D0 )
            H3=H3+(AH2(R)+YY)*8D0*R/(1D0+R)**2*C
           ENDIF
          END IF
        END IF
*
        SHINI =H0
        SHINTF=H4
        AHINI =H3
        AHINTF=H1
      ELSE
* ENERGY CUTS
        ksoft=0
        DEL=(RECUT2-R)*R1MI
        CALL SETDEL
        INDA=0
        CALL SAVAR
        CALL PHASEREGC
        CALL HCUTACOL(H0,H4,H3,H1)
        SHINI =H0
        SHINTF=H4
        AHINI =H3
        AHINTF=H1
      END IF
*
      IF (R.LT.RACUT) THEN
* ACOLLINEARITY CUT
        ksoft=0
        DEL=1D0/2*(1-SQRT(ABS(1-SINAC2*(1+2*R+R2)*R1MI2)/COSAC2))
        CALL SETDEL
        INDA=0
        CALL SAVAR
        CALL PHASEREGC
        CALL HCUTACOL(H0,H4,H3,H1)
        SHINI =SHINI-H0
        SHINTF=SHINTF-H4
        AHINI =AHINI-H3
        AHINTF=AHINTF-H1
      ENDIF
*
      ENDIF
*
*------------------------------------------------------------------
*
      CALL BORN(IFINAL,R,R,SBORN,ABORN,SBORNS,ABORNS)
*
      SHARD=ALQE2*SHINI*SBORN
      AHARD=ALQE2*AHINI*ABORN
*
      IF (INTERF.GE.1) THEN
        CALL BORN(-1,R,1D0,SBORN,ABORN,SBORNS,ABORNS)
        IF(INDF.EQ.10.AND.INTRF.EQ.2.AND.IMISD.EQ.1) THEN
          SBORN=SBORNS
          ABORN=ABORNS
        ENDIF
        SHARD=SHARD
     &       +ALQEF*(SHINTF*ABORN+8*C*R1MI*ALMI*ABORNI*ksoft)
        AHARD=AHARD
     &       +ALQEF*(AHINTF*SBORN+4*COPL2*R1MI*ALMI*SBORNI*ksoft)
      ENDIF
      HARD=SHARD+AHARD
*
      END

      SUBROUTINE HCUT(H0,H4,H3,H1)
*     ========== =================
**  06/06/90 08:18pm
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
*
      H4=0.D0
      H1=0.D0
      CALL SETRC( 1)
      CALL HINIM(H0M,H3M)
      H0=H0M
      H3=H3M
      IF (INTERF.GE.1) THEN
        CALL HINTFM(H4M,H1M)
        H4=H4M
        H1=H1M
      END IF
      CALL SETRC(-1)
      CALL HINIM(H0M,H3M)
      H0=H0+H0M
      H3=H3-H3M
      IF (INTERF.GE.1) THEN
        CALL HINTFM(H4M,H1M)
        H4=H4+H4M
        H1=H1-H1M
      END IF
*
      END
 
      SUBROUTINE HINIM(H0M,H3M) 
*     =========================
**  06/06/90 11:22pm
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /RCVAR / C,C2,C3,CC,CCI,CCI2,CCI3,CCI4,CDEL,CDELM
     +                 ,ALCC,ALCD,ALCDM,ALMI,ALPL
      COMMON /DELVAR/ DEL,DELM,DIF,DIF3,ALD,ALDM,ALDPL,ALDMI
     +               ,RD,RDI,RDI2,RDI3,RDM,RDMI,RDMI2,RDMI3
     +               ,ALRD,ALRDM
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
*
      IF(CDEL.GT.0.D0.AND.CDELM.GT.0.D0) THEN
* REGION ( +,+ )
       ALCDMI=ALCDM-ALCD
       H0M=DIF3*2D0/3*(R4-2*R3+2*R-1)
     . +2*DIF*(-R4+3*R3-3*R2+R)
     . +CCI*DIF3*2D0/3*(-R5+3*R4-4*R3+4*R2-3*R+1)
     . -CCI2*ALCDMI*(R5+2*R3+R)
     . +CCI2*DIF*(R5-R)
     . +2*CCI3*ALCDMI*(R5+R4+R3+R2)
     . +2*CCI3*DIF*(-R5+R4-R3+R2)
     . -2*CCI4*ALCDMI*(R5+R3)
*
       H3M=-CCI2*ALCDMI*(R4+R3+R2+R)
     . +2*CCI2*DIF*(R4-R3+R2-R)
     . +2*CCI3*ALCDMI*(R4+R2)
*
      ELSEIF(CDELM.LT.0.D0.AND.CDEL.LT.0.D0) THEN
* REGION ( -,- )
       ALCDMI=ALCDM-ALCD
       H0M=-(2D0/3*DIF3*(R4-2*R3+2*R-1)
     . +2*DIF*(-R3+3*R2-3*R+1)
     . +CCI*DIF3*2D0/3*(-R5+3*R4-4*R3+4*R2-3*R+1)
     . -CCI2*ALCDMI*(R5+2*R3+R)
     . +CCI2*DIF*(R5-R)
     . +2*CCI3*ALCDMI*(R5+R4+R3+R2)
     . +2*CCI3*DIF*(-R5+R4-R3+R2)
     . -2*CCI4*ALCDMI*(R5+R3))
*
       H3M=-(-CCI2*ALCDMI*(R4+R3+R2+R)
     . +2*CCI2*DIF*(R4-R3+R2-R)
     . +2*CCI3*ALCDMI*(R4+R2))
      ELSE
* REGION ( -,+ )
* 
       ALPM=2*ALCC+ALCD+ALCDM-2*ALR-ALPL+TE+1
       H0M=CCI2*ALPM*(R5+2*R3+R)
     . -2*CCI3*ALPM*(R5+R4+R3+R2)
     . +2*CCI4*ALPM*(R5+R3)
     . +4*DEL*DELM*(R3-2*R2+R)
     . +2*CCI2*DEL*DELM*(-R5+2*R4-2*R3+2*R2-R)
     . +2*DEL*(R4-4*R3+6*R2-4*R+1)
     . +2D0/3*(-R4+5*R3-6*R2+5*R-1)
     . +CCI*2D0/3*(R5+R3+R2+1)
     . +CCI2*R*1D0/3*(-9*R4-24*R3-26*R2-24*R-9)
     . +CCI3*4D0/3*(6*R5+11*R4+11*R3+6*R2)
     . +CCI4*2D0/3*(-11*R5-6*R4-11*R3)
*
       H3M=CCI2*ALPM*(R4+R3+R2+R)
     . -2*CCI3*ALPM*(R4+R2)
     . +2*CCI*DEL*DELM*(R4-2*R3+2*R2-2*R+1)
     . +2*CCI*(R3+R)
     . -4*CCI2*(R4+2*R3+2*R2+R)
     . +2*CCI3*(3*R4+2*R3+3*R2)
*
      ENDIF
*
      H0M=H0M*R1MI3
      H3M=H3M*R1MI2
*
      END
 
      SUBROUTINE HINTFM(H4M,H1M)
*     ========== ===============
**  06/07/90 06:44pm
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /RCVAR / C,C2,C3,CC,CCI,CCI2,CCI3,CCI4,CDEL,CDELM
     +                 ,ALCC,ALCD,ALCDM,ALMI,ALPL
      COMMON /DELVAR/ DEL,DELM,DIF,DIF3,ALD,ALDM,ALDPL,ALDMI
     +               ,RD,RDI,RDI2,RDI3,RDM,RDMI,RDMI2,RDMI3
     +               ,ALRD,ALRDM
*
      ALDELR=ALRD-ALRDM
      IF(CDEL.GT.0.D0.AND.CDELM.GT.0.D0) THEN
* REGION ( +,+ )
       H1M=-2*C*R*(2*DEL*(R-1)-ALDELR+R*(ALDM-ALD)-R+1)
       H4M=-2*R*(-2*DEL*(R-1)+ALDM-ALD+R-1)
      ELSE IF(CDELM.LT.0.D0.AND.CDEL.LT.0.D0) THEN
* REGION ( -,- )
       H4M=2*(2*DEL*(R-1)+ALDM*R-ALD*R-R+1)
       H1M=2*C*(2*DEL*(-R+1)+ALDELR*R+ALDM-ALD+R-1)
      ELSE
* REGION ( -,+ )
       H4M=2*DEL*(R2-1)+2*ALCC*C*(R2+2*R+1)+4*CCI*R1MI*(R3+R)
     . -ALMI*R1MI*C*(R3+R2+R+1)+2*ALR*R-8*R1MI+4*R+8
*
       H1M=-2*DEL*C*(R2-1)+2*ALDELR*C*R
     . +2*ALCC*(C2*(R2+R+1)-C*(R2-1)+R2-R+1)
     . -4*CCI*R1MI2*(R3+R2)-2*CCI2*R1MI2*(R4+R2)
     . +ALPL*C*(R2-1)-ALMI*R1MI*(C2+1)*(R3+1)
     . -2*ALDM*C*(R2-1)+2*ALR*C*R2-8*R1MI*C+4*C*(R+2)
      END IF
      END
 
      SUBROUTINE BOXIN(BOXI)
*     ========== ===========
      IMPLICIT COMPLEX*16(X)
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
*
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /PSCONS/ SW2,AMZ,GAMZ
      COMMON /PCONST/ ALFAI,AL1PI,ALQE2,ALQF2,ALQEF,GMU,CSIGNB
      COMMON /FORCHI/ XKAPP,XKAPPC,XMZ2,XMZ2C
      COMMON /SVAR  / S
      COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     +                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM
      COMMON /COUPL / VEFA,XVEFI,VEFZ,AEFA,XAEFI,AEFZ,VEEZ,XVPOL,VPOL2
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
*
      XR=XMZ2/S
      XR2=XR*XR
      XLR=LOG(XR)
      XL1R=LOG(1-1D0/XR)
      XDL1R=XSPENZ(1-1D0/XR)
      XDLCPR=XSPENZ(1-CP/XR)
      XDLCMR=XSPENZ(1-CM/XR)
      XDLPLR=XDLCPR+XDLCMR
      XDLMIR=XDLCPR-XDLCMR
*
      REF4BX=-(ALPL*ALMI*C+ALPL-ALMI*C)
      AMF4BX=2*PI*(ALMI*C-1)
      XF4BX=DCMPLX ( REF4BX, AMF4BX )
*
      XH4BX=-2*(ALPL*ALMI*C-2*ALPL*XL1R*(XR2-XR)-ALPL*(XR-1)
     . +2*ALMI*XL1R*C*(XR+1)+ALMI*C*(XR-1)-2*XDLPLR*(XR2-XR)
     . +2*XDLMIR*C*XR+2*XLR*(XR-1)-2*XL1R*(XR2-2*XR+1)
     . +4*XDL1R*(XR2-XR))
*
      XG4BX=1D0/2*(XF4BX+XH4BX)
*
      XCHI=XKAPP*S/(S-XMZ2)
      CHI2=DREAL(XCHI)**2+DIMAG(XCHI)**2
      BOXIS=-ALQEF*
     & (AEFA*VPOL2*REF4BX
     &+2*DREAL(XAEFI*DCONJG(XVPOL)*XCHI*XG4BX)+AEFZ*CHI2*DREAL(XH4BX))
*
      REF1BX=-(ALCP**2*C+ALCM**2*C-ALPL*C+ALMI)
      AMF1BX=-2*PI*(ALPL*C+ALMI*COPL2-C)
      XF1BX=DCMPLX ( REF1BX, AMF1BX )
*
      XH1BX=ALPL*ALMI*COPL2-4*ALPL*XL1R*C*(XR-1)-2*ALPL*C*(XR-1)
     . +4*ALMI*XL1R*(COPL2+XR2-XR)+2*ALMI*(XR-1)-4*XDLPLR*C*(XR-1)
     . +2*XDLMIR*(COPL2+2*(XR2-XR))+4*XLR*C*(XR-1)
     . -4*XL1R*C*(XR2-2*XR+1)+8*XDL1R*C*(XR-1)
*
      XG1BX=1D0/2*(XF1BX+XH1BX)
*
      BOXIA=ALQEF*(VEFA*VPOL2*REF1BX
     & +2*DREAL(XVEFI*DCONJG(XVPOL)*XCHI*XG1BX)+VEFZ*CHI2*DREAL(XH1BX))
*
      BOXI=BOXIA+BOXIS
*
      END
*
* Below are SUBROUTINEs for BHANG
*
      SUBROUTINE BHAFFF(IFIT,S,T)
*     ===========================
*----------------------------------------------------------------------
* SMALL ROUTINE FOR BHANG ( WEAK FORM FACTORS TRANSFER)
*----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
* FROM...ZINITF
      COMMON/XFORMZ/XFZ(4),XFZT,XFZTA
      COMMON/FRINIT/ NPAR(30),ZPAR(31)
       COMMON/FORBHA/ SW2
* FROM...BHANG
       COMMON /BHAPHY/ ALFAI,AL1PI,DUM(3)
* FOR...BHANG
       COMMON /BHASM / XFFS(4),XFFT(4),SW2SIR
*
       SW2SIR=SW2
       XFFS(1)=1D0/(2D0-XFZT)
       XFFS(2)=XFZ(1)
       XFFS(3)=XFZ(2)
       XFFS(4)=XFZ(4)
       IF (IFIT.EQ.-1) THEN
         XFFT(1)=(1D0,0D0)
       ELSE
         IHVP=NPAR(2)
         IQCD=NPAR(3)
         Q2=-T
         XFOTT=1.D0+0.25D0*AL1PI*XFOTF1(3,IHVP,IQCD,1,Q2)
         XFFT(1)=1D0/(2D0-XFOTT)
       END IF
       DO 10 I=2,4
         XFFT(I)=(1D0,0D0)
 10    CONTINUE
*
      END
 
      SUBROUTINE BHANG(INTRF,S,AMZ,GAMZ,WIDTHS,SW2,NPAR,ZPAR
     +                                     ,SIGBRN,SIGQED,AFBBRN,AFBQED)
*     ==================================================================
*-----------------------------------------------------------------------
*     Author M.Bilenky
*-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
*
      COMMON/XFORMZ/XFZ(4),XFZT,XFZTA
      COMMON/EWFORM/XALLCH(5,4),XFOTF
      COMMON/GAMFIT/GAMEEC,GAMFIC
      COMMON/ASSFIT/ROEC,GVEC,GAEC,ROFC,GVFC,GAFC
      COMMON/AS2FIT/ROL2C,GVL2C,GAL2C
      COMMON/FORBHA/SW2SIR
      COMMON/CORINT/CORINT
      COMMON/CDZBOX/QE,QF,ALAM1,ALAM2,HELI1,HELI2,SC,VE,VF,CRFAC
      COMMON/CDZIBF/IBFLA
*
      DIMENSION NPAR(30),ZPAR(31),WIDTHS(0:11)
*
      EXTERNAL DZEWBX
*
* Setting of IBFLA for DZEWBX
*
      IBFLA=0
*   
      SQS=SQRT(S)
      SW2SIR=SW2
      IBOX = NPAR(4)
      INTERF=NPAR(8)
      IFINAL=NPAR(9)
      IPHOT2=NPAR(10)
      IPART =NPAR(17)
      IFIT  =INTRF-1
      CALL BHAINI(IPART,IFINAL,INTERF,IPHOT2)
      QE=ZPAR(1)
      QF=QE
      THEMAX=ZPAR(29)
      THEMIN=ZPAR(30)
      ZPRIPP=ZPAR(31)
      ACOL=ZPAR(27)
      EMIN=ZPAR(28)
      XFZ(1)=XALLCH(2,1)
      XFZ(2)=XALLCH(2,2)
      XFZ(3)=XALLCH(2,3)
      XFZ(4)=XALLCH(2,4)
      XFZT  =XFOTF
*
      IF(INTRF.EQ.1) THEN
        Q2=S/2D0
        U2=S/2D0
* IF(NPAR(1).EQ.-1) the same as 1 and 2 here
        IF(NPAR(1).GE.1.OR.NPAR(1).EQ.-1) 
     &  CALL ROKANC(0,0,U2,-S,-Q2,QE,QF,XFZ,XFZT,XFZTA)
        IAFB=NPAR(13)
        CALL BHATOT(IAFB,IFIT
     &             ,SQS,THEMIN,THEMAX,ACOL,EMIN,EMIN
     &             ,AMZ,GAMZ,GAMEE,RO,GV,GA,RO2,GV2,GA2
     &             ,CROSB,CROSQ,ASYMB,ASYMQ)
*
      ELSEIF(INTRF.EQ.2) THEN
*     BHANGF(S,AMZ,GAMZ,WIDTHS,SW2,NPAR,ZPAR,SBORN,STOT)
*
        IAFB=0
        GAMEE=GAMEEC
        CALL BHATOT(IAFB,IFIT
     &             ,SQS,THEMIN,THEMAX,ACOL,EMIN,EMIN
     &             ,AMZ,GAMZ,GAMEE,RO,GV,GA,RO2,GV2,GA2
     &             ,CROSB,CROSQ,ASYMB,ASYMQ)
*
      ELSEIF(INTRF.EQ.3) THEN
*     BHANGA(S,AMZ,GAMZ,WIDTHS,SW2,NPAR,ZPAR,SBORN,STOT,ABORN,ATOT)
*
        IAFB=1
        RO=ROEC
        GV=GVEC
        GA=GAEC
        CALL BHATOT(IAFB,IFIT
     &             ,SQS,THEMIN,THEMAX,ACOL,EMIN,EMIN
     &             ,AMZ,GAMZ,GAMEE,RO,GV,GA,RO2,GV2,GA2
     &             ,CROSB,CROSQ,ASYMB,ASYMQ)
      ELSEIF(INTRF.EQ.4) THEN
*
        IAFB=1
        RO2=ROL2C
        GV2=GVL2C
        GA2=GAL2C
        CALL BHATOT(IAFB,IFIT
     &             ,SQS,THEMIN,THEMAX,ACOL,EMIN,EMIN
     &             ,AMZ,GAMZ,GAMEE,RO,GV,GA,RO2,GV2,GA2
     &             ,CROSB,CROSQ,ASYMB,ASYMQ)
      ENDIF
      SIGBRN=CROSB
      AFBBRN=ASYMB
      SIGQED=CROSQ
      AFBQED=ASYMQ
*
      IBOX = NPAR(4)
      SBOX = 0.D0
      ABOX = 0.D0
      IF(IBOX.NE.1) GOTO 90
*
      QE    =ZPAR(1)
      QF    =ZPAR(2)
      AMF   =ZPAR(4)
      ALAMP= ZPAR(22)
      ALAME= ZPAR(23)
      HELPL= ZPAR(24)
      HELMI= ZPAR(25)
      ALAM1=1D0-ALAMP*ALAME
      ALAM2=ALAMP-ALAME
      HELI1=(1D0-HELPL*HELMI)/4D0
      HELI2=(HELPL-HELMI)/4D0
      QFM=DABS(QF)
      QEM=DABS(QE)
      VE=1D0-4D0*SW2*QEM
      VF=1D0-4D0*SW2*QFM
      SC=S
      FMUSQ= DSQRT(1D0-4D0*AMF*AMF/SC)
      CNANOB=.38937966D6
      CRFAC=CNANOB*FMUSQ
*
      ANG1=ZPAR(29)
      ANG2=ZPAR(30)
      ZPRIPP=ZPAR(31)
      PI=4D0*ATAN(1D0)
      c1=cos(ang1*pi/180d0)
      c2=cos(ang2*pi/180d0)
      eps = 1.d-3
      amin= c1+1.d-8
      amax= c2-1.d-8
      st  = (amax-amin)/2
      call simps(amin,0d0,st,eps,1.d-18,DZEWBX,ac,resb,res2,res3)
      call simps(0d0,amax,st,eps,1.d-18,DZEWBX,ac,resf,res2,res3)
      sbox = resf+resb
      abox = resf-resb
      afbqed=afbqed*sigqed
      sigqed=sigqed + sbox*corint
      afbqed=afbqed + abox*corint
      if(sigqed.ne.0d0) afbqed=afbqed/sigqed
90    continue
*
      END
