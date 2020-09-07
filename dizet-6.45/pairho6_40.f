************************************************************************
*     higher order and singlet contributions of paris                  *
*     1) according to ABA                                              *
*     2) according to JMS                                              *
************************************************************************
      DOUBLE PRECISION FUNCTION PH2ADD(X)
C     -----------------------------------
C.... ADDENDUM TO RADIATOR FOR LIGHT PAIR CONTRIBUTION 
C     TO E^+E^- ANNIHILATION
C     HARD PAIRS (SEE SOFT+VIRT ONES IN PAIRS(X))
C     CALLED FROM FUNCTION PH2
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
      COMMON /CPAIRS/ FPI2,RINF,RRR0,RRR1,PAIRDL,CORFAC
      PAIRH= 0D0
      Z    = 1D0 - X

      DLZ  = DLOG(Z)
      DL1Z = DLOG(X)
      DL1PZ= DLOG(1D0+Z)
      SP1Z = SPENCE(X)
      SPMZ = SPENCE(-Z)
      TR1Z = TRILOG(X)
      TPI  = DLOG(SCOM*X**2/FPI2/Z)
C
      IF(IFLAGS(IFISPP).EQ.4) GOTO 77
C
C.... SINGLET ELECTRON PAIRS ACCORDING TO BNB:     
      IF(IFLAGS(IFIPFC).EQ.1.OR.IFLAGS(IFIPFC).EQ.5) THEN
      IF(IFLAGS(IFIPSC).GE.1) 
     & PAIRH= PAIRH
     &      + TEE**2*( 0.5D0*(1D0+Z)*DLZ + 1D0/(3D0*Z) + 0.25D0
     &      - 0.25D0*Z - Z**2/3D0 ) 
      IF(IFLAGS(IFIPSC).GE.2) 
     & PAIRH= PAIRH
     &      + TEE*( (1D0+Z)*(2D0*DLZ*DL1Z - DLZ**2
     &      + 2D0*SP1Z ) + ( 4D0/(3D0*Z) + 1D0 - Z 
     &      - 4D0/3D0*Z**2 )*DL1Z - ( 2D0/(3D0*Z) + 1D0
     &      - Z/2D0 - 4D0/3D0*Z**2 )*DLZ - 8D0/(9D0*Z) - 8D0/3D0
     &      + 8D0/3D0*Z + 8D0/9D0*Z**2 )
     &      + (1D0+Z)*( 2D0*DLZ*DL1Z**2 
     &      - 2D0*DLZ**2*DL1Z + 4D0*DL1Z*SP1Z
     &      - 5D0*SPMZ - 5D0*DLZ*DL1PZ )
     &      + ( 4D0/(3D0*Z) + 1D0 - Z - 4D0/3D0*Z**2 )*DL1Z**2
     &      - ( 4D0/(3D0*Z) + 2D0 - Z - 8D0/3D0*Z**2 )*DLZ*DL1Z
     &      - ( 16D0/(9D0*Z) + 16D0/3D0 - 16D0/3D0*Z
     &      - 16D0/9D0*Z**2 )*DL1Z 
     &      + ( - 4D0/Z + 2D0*Z )*S12(X)
     &      + ( - 2D0/Z - 1D0 - 11D0/2D0*Z )*TR1Z
     &      + ( 10D0/Z - 10D0 + 3D0*Z )*(0.25D0*TRILOG(Z**2)-TRILOG(Z))
     &      + ( 2D0/Z + 2D0 + Z )*( 6D0*( S12(Z/(1D0+Z))
     &      + DL1PZ**3/6D0 ) + 6D0*DL1PZ*SPMZ
     &      + 3D0*DLZ*DL1PZ**2
     &      - 5D0/2D0*DLZ**2*DL1PZ + 3D0*F1*DL1PZ )
     &      + ( - 5D0 + Z/2D0 )*DLZ*SP1Z
     &      + ( - 10D0/Z - 4D0*Z )*DLZ*SPMZ 
     &      + ( 2D0/3D0 + Z )*DLZ**3
     &      + ( 2D0/(3D0*Z) - 1.5D0 + 2D0/3D0*Z**2 )*SP1Z
     &      + ( 1D0/(3D0*Z) + 1.5D0 + 9D0/4D0*Z - Z**2 )*DLZ**2
     &      + ( 8D0/(9D0*Z) + 51D0/4D0 + Z/2D0 - 16D0/9D0*Z**2
     &      - 1.5D0*Z*F1 - 6D0*F1 )*DLZ
     &      + ( - 2D0/(3D0*Z) - 3D0 - 2D0*Z + 2D0/3D0*Z**2 )*F1
     &      + ( 6D0/Z - 9D0 + 1.5D0*Z )*ZET3
     &      + ( 44D0/(27D0*Z) + 1477D0/72D0 - 385D0/18D0*Z 
     &      - 163D0/216D0*Z**2 ) 
C.... INTERFERENCE OF SINGLET AND NON-SINGLET PAIRS (BNB):
     &      + TEE*( (1D0+Z**2)/X*( - SP1Z
     &      - DLZ**2/2D0 - 3D0/4D0*DLZ ) - 7D0/4D0*(1D0+Z)*DLZ
     &      - 4D0 + 7D0/2D0*Z )
     &      + 2D0*DL1Z*( (1D0+Z**2)/X*( - SP1Z
     &      - DLZ**2/2D0 - 3D0/4D0*DLZ ) - 7D0/4D0*(1D0+Z)*DLZ 
     &      - 4D0 + 7D0/2D0*Z ) 
     &      + (1D0+Z**2)/X*( 3D0*TR1Z - 8D0*S12(X)
     &      + 4D0*( 0.25D0*TRILOG(Z**2) - TRILOG(Z) ) + 3D0*ZET3 
     &      - 2D0*DLZ*SPMZ - 6D0*DLZ*SP1Z
     &      + 3D0/2D0*F1*DLZ - DLZ**3/4D0
     &      + 35D0/12D0*DLZ - 5D0/4D0*SP1Z + DLZ**2/2D0 )
     &      + (1D0+Z)*( - TR1Z + DLZ**3/12D0 - 0.5D0*DLZ*SP1Z )
     &      + ( 2D0 + 1D0/X - 2D0/X**2 )*( F1
     &      + 2D0*SPMZ + 2D0*DLZ*DL1PZ )
     &      + 1D0/(1D0+Z)**2*( - 18D0*F1 - 6D0*DLZ**2
     &      - 12D0*SP1Z - 24D0*DLZ + 12D0 )
     &      + 1D0/(1D0+Z)*( 36D0*F1 + 12D0*DLZ**2 
     &      + 24D0*SP1Z + 7D0*DLZ - 18D0 ) + 12D0/(1D0+Z)**3*DLZ 
     &      + SP1Z*( - 5D0/4D0*Z - 73D0/4D0 )
     &      - 33D0/2D0*F1 + ( 5D0/2D0*Z - 8D0 )*DLZ**2
     &      + DLZ**2/X**2 
     &      + ( 83D0/12D0 - 13D0/12D0*Z )*DLZ
     &      - 2D0*DLZ/(X)**2 - 2D0/(X) 
     &      + 1495D0/72D0 - 905D0/72D0*Z - Z**2/9D0 
      ENDIF
C.... ADD IF REQUIRED THE THIRD ORDER LLA (NEW):
      IF(IFLAGS(IFIPTO).GE.1) THEN 
       RSZ  = (1D0-Z)/(3D0*Z)*( 4D0 + 7D0*Z + 4D0*Z**2 ) 
     &      + 2D0*(1D0+Z)*DLZ
       P1Z  = (1D0+Z**2)/(X)
       P2Z  = 2D0*( (1D0+Z**2)/(X)*( 2D0*DL1Z - DLZ
     &      + 1.5D0 ) + (1D0+Z)/2D0*DLZ - 1D0 + Z )
       P1RS = ( 1.5D0 + 2D0*DL1Z )*RSZ + (1D0+Z)*( - DLZ**2
     &      + 4D0*SP1Z ) + 1D0/3D0*( - 9D0 - 3D0*Z
     &      + 8D0*Z**2 )*DLZ + 2D0/3D0*( - 3D0/Z - 8D0 + 8D0*Z 
     &      + 3D0*Z**2 )
       IF(IFLAGS(IFIPFC).EQ.1) FCON3 = PPH(TEE,Z)
       IF(IFLAGS(IFIPFC).EQ.2) FCON3 = PPH(TMU,Z)
       IF(IFLAGS(IFIPFC).EQ.3) FCON3 = PPH(TTU,Z)
       IF(IFLAGS(IFIPFC).EQ.4) FCON3 = PHA(TPI,Z)
       IF(IFLAGS(IFIPFC).EQ.5) FCON3 = PPH(TEE,Z) 
     &                               + PPH(TMU,Z)
     &                               + PHA(TPI,Z)
     &                               + PPH(TTU,Z)
       IF(IFLAGS(IFIPFC).EQ.6) FCON3 = PPH(TEE,Z) 
     &                               + PPH(TMU,Z)
     &                               + PPH(TTU,Z)
       PAIRH= PAIRH + ALQE2*(TEE-1D0)*FCON3
       IF(IFLAGS(IFIPTO).GE.2) THEN
        IF(IFLAGS(IFIPFC).EQ.1.OR.IFLAGS(IFIPFC).EQ.5) PAIRH= PAIRH 
     &      + ALQE2*TEE**3*P1Z/27D0
       ENDIF

       IF(IFLAGS(IFIPSC).EQ.3) THEN
        IF(IFLAGS(IFIPFC).EQ.1.OR.IFLAGS(IFIPFC).EQ.5) PAIRH= PAIRH 
     &      + ALQE2*(TEE-1D0)**3*( 5D0/24D0*P1RS - RSZ/36D0 )
       ENDIF
      ENDIF
C.... NEW! THE FOURTH ORDER
      IF(IFLAGS(IFIPTO).EQ.3) THEN
       IF(IFLAGS(IFIPFC).EQ.1.OR.IFLAGS(IFIPFC).EQ.5) THEN 
       P3Z  = 
     &        48D0*( 0.5D0*(1D0+Z**2)/(X)*( 9D0/32D0 - PI**2/12D0
     &      + 3D0/4D0*DL1Z - 3D0/8D0*DLZ + 0.5D0*DL1Z**2
     &      + 1D0/12D0*DLZ**2 - 0.5D0*DLZ*DL1Z )
     &      + 1D0/8D0*(1D0+Z)*DLZ*DL1Z
     &      - 1D0/4D0*(X)*DL1Z + (5D0-3D0*Z)/32D0*DLZ
     &      - (1D0-Z)/16D0 - (1D0+Z)/32D0*DLZ**2
     &      + (1D0+Z)/8D0*SP1Z )
       FCOE4= P3Z/12D0 + P2Z*22D0/(16D0*27D0) + P1Z/(4D0*27D0)
       PAIRH= PAIRH + ALQE2**2*(TEE-1D0)**4*FCOE4
       ENDIF
      ENDIF
C.... RESTORE THE OVERALL FACTOR:
      PH2ADD= PAIRH*ALQE2**2
      RETURN
 77   CONTINUE
C.... PURE HADRONIC THIRD ORDER CONTRIBUTION TO BE ADDED TO JMS
      FCON3 = 0D0
      IF(IFLAGS(IFIPFC).EQ.4.OR.IFLAGS(IFIPFC).EQ.5) FCON3 = PHA(TPI,Z)
      PH2ADD= ALQE2**3*(TEE-1D0)*FCON3
      RETURN
      END

      DOUBLE PRECISION FUNCTION PPH(TXX,Z)
C     ------------------------------------
C.... P*PH2 CONVOLUTED FUNCTION (THETA PART)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /NCONST/ PI,F1,AL2,ZET3
      DZ2  = F1
      DZ3  = ZET3
      ALS  = TXX
 
      G    = (1D0+Z**2)/(1D0-Z)
      DLZ  = DLOG(Z)
      DL1Z = DLOG(1D0-Z)
      SP1Z = SPENCE(1D0-Z)
      TRIZ = TRILOG(1D0-Z)
      S12Z = S12(1D0-Z)
 
      PPH =  G*ALS**2 * ( 1./2. - 1./3.*DLZ + 2./3.*DL1Z )
 
      PPH = PPH + G*ALS * (  - 9./4. - 4./3.*dz2 + 1./3.*DLZ**2 -2.*DLZ
     +    *DL1Z + 11./18.*DLZ + 2.*DL1Z**2 - 11./9.*DL1Z )
 
      PPH = PPH + G * ( 1073./162. + 8./3.*dz2*DLZ - 16./3.*dz2*DL1Z - 
     +    4./3.*dz2 + 4.*dz3 - 2./3.*SP1Z*DLZ - 1./3.*SP1Z*DL1Z -1./4.*
     +    SP1Z - 1./18.*DLZ**3 + 5./6.*DLZ**2*DL1Z - 31./72.*DLZ**2 - 8.
     +    /3.*DLZ*DL1Z**2 + 7./3.*DLZ*DL1Z - 67./54.*DLZ + 16./9.*
     +    DL1Z**3 - 7./3.*DL1Z**2 + 67./27.*DL1Z - 5./3.*S12Z )
 
      PPH = PPH + ALS**2 * (  - 1./3. + 1./6.*z*DLZ + 1./3.*z + 1./6.*
     +    DLZ )
 
      PPH = PPH + ALS * ( 19./9. + 2./3.*z*SP1Z - 1./6.*z*DLZ**2 + 2./3.
     +    *z*DLZ*DL1Z - 20./9.*z*DLZ + 8./3.*z*DL1Z - 19./9.*z + 2./3.*
     +    SP1Z - 1./6.*DLZ**2 + 2./3.*DLZ*DL1Z + 4./9.*DLZ - 8./3.*DL1Z
     +     )
 
      PPH = PPH - 319./54. -4.*z*dz2 - 2./3.*z*dz3 - 2./3.*z*SP1Z*DLZ
     +     + 1./3.*z*SP1Z*DL1Z - 67./36.*z*SP1Z - 1./2.*z*DLZ**2*DL1Z
     +     + 77./72.*z*DLZ**2 - 52./9.*z*DLZ*DL1Z + 184./27.*z*DLZ +4.*
     +    z*DL1Z**2 - 76./9.*z*DL1Z - 1./3.*z*TRIZ + 319./54.*z +4.*dz2
     +     - 2./3.*dz3 - 2./3.*SP1Z*DLZ + 1./3.*SP1Z*DL1Z - 67./36.*
     +    SP1Z - 1./2.*DLZ**2*DL1Z + 5./72.*DLZ**2 + 14./9.*DLZ*DL1Z - 
     +    97./54.*DLZ -4.*DL1Z**2 + 76./9.*DL1Z - 1./3.*TRIZ + 4./3.*
     +    S12(Z) + 4./3.*S12(Z)*z 
     +    - 2./3.*TRILOG(Z) - 2./3.*TRILOG(Z)*z
      RETURN
      END

      DOUBLE PRECISION FUNCTION PPH0(TXX,DELL)
C     ----------------------------------------
C.... P*PH2 CONVOLUTED FUNCTION (DELTA PART)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /NCONST/ PI,F1,AL2,ZET3
      DZ2  = F1
      DZ3  = ZET3
      ALS  = TXX
      ALD  = DELL/2D0
      PPH0 = 
     +  + ALS*ALD * (  - 9./2. - 8./3.*dz2 )
      PPH0 = PPH0 + ALS*ALD**2 * (  - 11./9. )
      PPH0 = PPH0 + ALS*ALD**3 * ( 4./3. )
      PPH0 = PPH0 + ALS * (  - 17./8. )
      PPH0 = PPH0 + ALS**2*ALD * ( 1. )
      PPH0 = PPH0 + ALS**2*ALD**2 * ( 2./3. )
      PPH0 = PPH0 + ALS**2 * ( 3./8. - 2./3.*dz2 )
      PPH0 = PPH0 + ALD * ( 1073./81. - 8./3.*dz2 + 8*dz3 )
      PPH0 = PPH0 + ALD**2 * ( 67./27. - 16./3.*dz2 )
      PPH0 = PPH0 + ALD**3 * (  - 14./9. )
      PPH0 = PPH0 + ALD**4 * ( 8./9. )
      PPH0 = PPH0 + 821./108. - 23./6.*dz2 + 2.*dz3
      RETURN
      END

      DOUBLE PRECISION FUNCTION PHA(TPI,Z)
C     ------------------------------------
C.... P*HADRONS CONVOLUTED FUNCTION (THETA PART)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /CPAIRS/ FPI2,RINF,RRR0,RRR1,PAIRDL,CORFAC
      RI   = RINF
      R0   = RRR0
      R1   = RRR1
 
      DZ2  = F1
      DZ3  = ZET3
      G    = (1D0+Z**2)/(1D0-Z)
      DLZ  = DLOG(Z)
      DL1Z = DLOG(1D0-Z)
      SP1Z = SPENCE(1D0-Z)
      TRIZ = TRILOG(1D0-Z)
      S12Z = S12(1D0-Z)
      ALH  = TPI + 2D0*AL2 + DLZ - 2D0*DL1Z
 
      PHA =
     +  + G*ALH * (  - 4./3.*dz2*Ri - 2*DLZ*DL1Z*Ri - 1./2.*DLZ*Ri - 2./
     +    3.*DLZ*R0 + 1./3.*DLZ**2*Ri + DL1Z*Ri + 4./3.*DL1Z*R0 + 2*
     +    DL1Z**2*Ri - 7./12.*Ri + R0 + 4./3.*AL2*DLZ*Ri - 8./3.*AL2
     +    *DL1Z*Ri - 2*AL2*Ri )
      PHA = PHA + G*ALH**2 * (  - 1./3.*DLZ*Ri + 2./3.*DL1Z*Ri + 1./2.*
     +    Ri )
      PHA = PHA + G * ( 2*dz2*DLZ*Ri - 4*dz2*DL1Z*Ri - dz2*Ri - 4./3.*
     +    dz2*R0 + 10./3.*dz3*Ri - 2./3.*SP1Z*DLZ*Ri - 1./3.*SP1Z*DL1Z*
     +    Ri - 1./4.*SP1Z*Ri - DLZ*DL1Z*Ri - 2*DLZ*DL1Z*R0 - 8./3.*DLZ*
     +    DL1Z**2*Ri - 1./2.*DLZ*R0 - 2./3.*DLZ*R1 + 5./6.*DLZ**2*DL1Z*
     +    Ri + 1./8.*DLZ**2*Ri + 1./3.*DLZ**2*R0 - 1./18.*DLZ**3*Ri + 
     +    DL1Z*R0 + 4./3.*DL1Z*R1 + DL1Z**2*Ri + 2*DL1Z**2*R0 + 16./9.*
     +    DL1Z**3*Ri - 5./3.*S12Z*Ri + 5./8.*Ri - 7./12.*R0 + R1 + 8./3.
     +    *AL2*dz2*Ri + 4*AL2*DLZ*DL1Z*Ri + AL2*DLZ*Ri + 4./3.*
     +    AL2*DLZ*R0 - 2./3.*AL2*DLZ**2*Ri - 2*AL2*DL1Z*Ri - 8./3.
     +    *AL2*DL1Z*R0 - 4*AL2*DL1Z**2*Ri + 7./6.*AL2*Ri - 2.*AL2
     +    *R0 - 4./3.*AL2*AL2*DLZ*Ri + 8./3.*AL2*AL2*DL1Z*Ri
     +     + 2*AL2*AL2*Ri )
      PHA = PHA + ALH * ( 2./3.*z*SP1Z*Ri + 2./3.*z*DLZ*DL1Z*Ri - 5./3.
     +    *z*DLZ*Ri + 1./3.*z*DLZ*R0 - 1./6.*z*DLZ**2*Ri + 8./3.*z*DL1Z
     +    *Ri - z*Ri + 2./3.*z*R0 + 2./3.*SP1Z*Ri + 2./3.*DLZ*DL1Z*Ri
     +     + DLZ*Ri + 1./3.*DLZ*R0 - 1./6.*DLZ**2*Ri - 8./3.*DL1Z*Ri + 
     +    Ri - 2./3.*R0 - 2./3.*AL2*z*DLZ*Ri - 4./3.*AL2*z*Ri - 2./
     +    3.*AL2*DLZ*Ri + 4./3.*AL2*Ri )
      PHA = PHA + ALH**2 * ( 1./6.*z*DLZ*Ri + 1./3.*z*Ri + 1./6.*DLZ*Ri
     +     - 1./3.*Ri )
      PHA = PHA + 1./3.*z*dz2*DLZ*Ri - 10./3.*z*dz2*Ri - 2./3.*z*dz3*Ri
     +     - 2./3.*z*SP1Z*DLZ*Ri + 1./3.*z*SP1Z*DL1Z*Ri - 3./4.*z*SP1Z*
     +    Ri + 2./3.*z*SP1Z*R0 - 14./3.*z*DLZ*DL1Z*Ri + 2./3.*z*DLZ*
     +    DL1Z*R0 + 3*z*DLZ*Ri - 5./3.*z*DLZ*R0 + 1./3.*z*DLZ*R1 - 1./2.
     +    *z*DLZ**2*DL1Z*Ri + 19./24.*z*DLZ**2*Ri - 1./6.*z*DLZ**2*R0
     +     - 4*z*DL1Z*Ri + 8./3.*z*DL1Z*R0 + 4*z*DL1Z**2*Ri - 1./3.*z*
     +    TRIZ*Ri + 13./6.*z*Ri - z*R0 + 2./3.*z*R1 + 1./3.*dz2*DLZ*Ri
     +     + 10./3.*dz2*Ri - 2./3.*dz3*Ri - 2./3.*SP1Z*DLZ*Ri + 1./3.*
     +    SP1Z*DL1Z*Ri - 3./4.*SP1Z*Ri + 2./3.*SP1Z*R0 + 8./3.*DLZ*DL1Z
     +    *Ri + 2./3.*DLZ*DL1Z*R0 - 7./6.*DLZ*Ri + DLZ*R0 + 1./3.*DLZ*
     +    R1 - 1./2.*DLZ**2*DL1Z*Ri - 5./24.*DLZ**2*Ri - 1./6.*DLZ**2*
     +    R0 + 4*DL1Z*Ri - 8./3.*DL1Z*R0 - 4*DL1Z**2*Ri - 1./3.*TRIZ*Ri
     +     - 13./6.*Ri + R0 - 2./3.*R1 - 4./3.*AL2*z*SP1Z*Ri - 4./3.*
     +    AL2*z*DLZ*DL1Z*Ri + 10./3.*AL2*z*DLZ*Ri - 2./3.*AL2*z*
     +    DLZ*R0 + 1./3.*AL2*z*DLZ**2*Ri - 16./3.*AL2*z*DL1Z*Ri + 2
     +    *AL2*z*Ri
      PHA = PHA - 4./3.*AL2*z*R0 - 4./3.*AL2*SP1Z*Ri - 4./3.*AL2*
     +    DLZ*DL1Z*Ri - 2*AL2*DLZ*Ri - 2./3.*AL2*DLZ*R0 + 1./3.*AL2
     +    *DLZ**2*Ri + 16./3.*AL2*DL1Z*Ri - 2*AL2*Ri + 4./3.*AL2
     +    *R0 + 2./3.*AL2*AL2*z*DLZ*Ri + 4./3.*AL2*AL2*z*Ri + 
     +    2./3.*AL2*AL2*DLZ*Ri - 4./3.*AL2*AL2*Ri + 4./3.*S12(Z
     +    )*z*Ri + 4./3.*S12(Z)*Ri 
     +    - 2./3.*TRILOG(Z)*z*Ri - 2./3.*TRILOG(Z)*
     +    Ri
      HADH = PHA

      RETURN
      END

      DOUBLE PRECISION FUNCTION PHA0(TPI,DELL)
C     ----------------------------------------
C.... P*HADRONS CONVOLUTED FUNCTION (DELTA PART)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /CPAIRS/ FPI2,RINF,RRR0,RRR1,PAIRDL,CORFAC
      RI   = RINF
      R0   = RRR0
      R1   = RRR1
      DZ3  = ZET3
      ALD  = DELL/2D0
      ALH  = TPI + 2D0*AL2
 
      PHA0 = 
     +  + ALD * (  - 8./3.*ALH*F1*Ri - 7./6.*ALH*Ri + 2.*ALH*R0 + 
     +    ALH**2*Ri - 2.*F1*Ri - 8./3.*F1*R0 + 20./3.*dz3*Ri + 5./4.*
     +    Ri - 7./6.*R0 + 2.*R1 - 4.*AL2*ALH*Ri + 16./3.*AL2*F1*Ri
     +     + 7./3.*AL2*Ri - 4*AL2*R0 + 4*AL2*AL2*Ri )
      PHA0 = PHA0 + ALD**2 * ( ALH*Ri + 4./3.*ALH*R0 + 2./3.*ALH**2*Ri
     +     - 4.*F1*Ri + R0 + 4./3.*R1 - 8./3.*AL2*ALH*Ri - 2.*AL2*Ri
     +     - 8./3.*AL2*R0 + 8./3.*AL2*AL2*Ri )
      PHA0 = PHA0 + ALD**3 * ( 4./3.*ALH*Ri + 2./3.*Ri + 4./3.*R0 - 8./
     +    3.*AL2*Ri )
      PHA0 = PHA0 + ALD**4 * ( 8./9.*Ri )
      PHA0 = PHA0 - 7./8.*ALH*Ri + 3./4.*ALH*R0 + 3./8.*ALH**2*Ri - 3./
     +    4.*F1*Ri + dz3*Ri + 15./16.*Ri - 7./8.*R0 + 3./4.*R1 - 3./2.
     +    *AL2*ALH*Ri + 7./4.*AL2*Ri - 3./2.*AL2*R0 + 3./2.*AL2
     +    *AL2*Ri
      RETURN
      END

      DOUBLE PRECISION FUNCTION AJMSB(T)
C     ----------------------------------
C.... CONTRIBUTION FROM PAIRS AND PHOTONS WITH EXPONENTIATION
C     ACCORDING TO JMS 
C     VARIABLE EXCHANGED 
C     PHOTONS (THE "PRAGMATIC" FORMULA) ARE EXTRACTED
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

      BETE = 2D0*ALQE2*(TEE-1D0)
      A    = 1D0/BETE
      Z    = 1D0 - T**A
      X    = 1D0 - Z
      GEXP = A
C     GEXP = (1D0-Z)**(BETE-1D0)
C     DZ   = A*T**(A-1D0)

      DLZ  = DLOG(Z)
      IF(T.LE.0D0) THEN
       DL1Z= 0D0
      ELSE
       DL1Z= A*DLOG(T)
      ENDIF
      SP1Z = SPENCE(1D0-Z) 
 
      ALSE = TEE
      ALSF = TEE
      DZ2  = F1
      DZ3  = ZET3
      ALPI = ALQE2
 
      ALFE = DL1Z + TEE/2D0 - 5D0/6D0
      ALFM = DL1Z + TMU/2D0 - 5D0/6D0
      ALFT = DL1Z + TTU/2D0 - 5D0/6D0
 
      BETBE= BETE + 4D0/3D0*ALPI**2*( ALFE**2 + 31D0/36D0 )
      BETBM= BETE + 4D0/3D0*ALPI**2*( ALFM**2 + 31D0/36D0 )
      BETBT= BETE + 4D0/3D0*ALPI**2*( ALFT**2 + 31D0/36D0 )
 
      IF(IFLAGS(IFIPFC).EQ.1) BETB = BETBE
      IF(IFLAGS(IFIPFC).EQ.2) BETB = BETBM
      IF(IFLAGS(IFIPFC).EQ.3) BETB = BETBT
      IF(IFLAGS(IFIPFC).EQ.4) BETB = BETE
      IF(IFLAGS(IFIPFC).EQ.5.OR.IFLAGS(IFIPFC).EQ.6) BETB = BETE
     &     + 4D0/3D0*ALPI**2*( ALFE**2 + 31D0/36D0 )
     &     + 4D0/3D0*ALPI**2*( ALFM**2 + 31D0/36D0 )
     &     + 4D0/3D0*ALPI**2*( ALFT**2 + 31D0/36D0 )
 
      IF(IFLAGS(IFIPFC).EQ.1) BLLR = 2D0*ALPI*ALSE
     &     + 1D0/3D0*ALPI**2*TEE**2
      IF(IFLAGS(IFIPFC).EQ.2) BLLR = 2D0*ALPI*ALSE
     &     + 1D0/3D0*ALPI**2*TMU**2
      IF(IFLAGS(IFIPFC).EQ.3) BLLR = 2D0*ALPI*ALSE
     &     + 1D0/3D0*ALPI**2*TTU**2
      IF(IFLAGS(IFIPFC).EQ.4) BLLR = 2D0*ALPI*ALSE
      IF(IFLAGS(IFIPFC).EQ.5.OR.IFLAGS(IFIPFC).EQ.6) BLLR=2D0*ALPI*ALSE
     &     + 1D0/3D0*ALPI**2*TEE**2
     &     + 1D0/3D0*ALPI**2*TMU**2
     &     + 1D0/3D0*ALPI**2*TTU**2

      FBB  = 1D0 - BETB**2/2D0*DZ2 + BETB**3/3D0*DZ3
 
      IF(IFLAGS(IFIPFC).EQ.1) PHI = 
     &        DEXP(4D0/9D0*ALPI**2*PI**2/2D0*ALFE)
     &       *FBB*( BETB - 4D0/9D0*ALPI**2
     &       *( PI**2/2D0 + (3D0/2D0*PI**2*ALFE - 8D0*DZ3)*BETBE ) )
      IF(IFLAGS(IFIPFC).EQ.2) PHI = 
     &        DEXP(4D0/9D0*ALPI**2*PI**2/2D0*ALFM)
     &       *FBB*( BETB - 4D0/9D0*ALPI**2
     &       *( PI**2/2D0 + (3D0/2D0*PI**2*ALFM - 8D0*DZ3)*BETBM ) )
      IF(IFLAGS(IFIPFC).EQ.3) PHI = 
     &        DEXP(4D0/9D0*ALPI**2*PI**2/2D0*ALFT)
     &       *FBB*( BETB - 4D0/9D0*ALPI**2
     &       *( PI**2/2D0 + (3D0/2D0*PI**2*ALFT - 8D0*DZ3)*BETBT ) )
      IF(IFLAGS(IFIPFC).EQ.4) PHI = FBB*BETB
      IF(IFLAGS(IFIPFC).EQ.5.OR.IFLAGS(IFIPFC).EQ.6) PHI 
     &       = DEXP(4D0/9D0*ALPI**2*PI**2/2D0*(ALFE+ALFM+ALFT))
     &       *FBB*( BETB - 4D0/9D0*ALPI**2*( 
     &         ( PI**2/2D0 + (3D0/2D0*PI**2*ALFE - 8D0*DZ3)*BETB ) 
     &       + ( PI**2/2D0 + (3D0/2D0*PI**2*ALFM - 8D0*DZ3)*BETB ) 
     &       + ( PI**2/2D0 + (3D0/2D0*PI**2*ALFT - 8D0*DZ3)*BETB ) ) )
 
      DLL  = (1D0+Z**2)/2D0 + 0.25D0*BLLR*( - (1D0-Z)**2
     &     - 0.5D0*(1D0+3D0*Z**2)*DLZ )
     &     + 1D0/8D0*BLLR**2*( (3D0*Z**2-4D0*Z+1D0)/2D0*DLZ
     &     + (1D0+7D0*Z**2)/12D0*DLZ**2 + (1D0-Z)**2
     &     + (1D0-Z**2)*SP1Z ) 
      DNLL = ALPI*BETE*( 3D0/32D0 - PI**2/8D0 + 3D0/2D0*DZ3 )
     &     + ALPI/8D0*( 4D0*(1D0+Z**2)*( SP1Z
     &     + DLZ*DL1Z ) - (1D0+3D0*Z**2)*DLZ**2 
     &     + 4D0*( 2D0 + Z + 2D0*Z**2 )*DLZ 
     &     + 2D0*(1D0-Z)*(5D0-4D0*Z) )
 
      IF(IFLAGS(IFIPFC).EQ.1) ARGU = ALPI**2*( 4D0/9D0*ALFE**3
     &     + 2D0/3D0*( 31D0/18D0 - 2D0*DZ2 )*ALFE
     &     - TEE**3/18D0 + 19D0/36D0*TEE**2
     &     + ( PI**2/9D0 - 265D0/108D0 )*TEE ) 
      IF(IFLAGS(IFIPFC).EQ.2) ARGU = ALPI**2*( 4D0/9D0*ALFM**3
     &     + 2D0/3D0*( 31D0/18D0 - 2D0*DZ2 )*ALFM
     &     - TMU**3/18D0 + 19D0/36D0*TMU**2
     &     + ( PI**2/9D0 - 265D0/108D0 )*TMU ) 
      IF(IFLAGS(IFIPFC).EQ.3) ARGU = ALPI**2*( 4D0/9D0*ALFT**3
     &     + 2D0/3D0*( 31D0/18D0 - 2D0*DZ2 )*ALFT
     &     - TTU**3/18D0 + 19D0/36D0*TTU**2
     &     + ( PI**2/9D0 - 265D0/108D0 )*TTU ) 
      IF(IFLAGS(IFIPFC).EQ.4) ARGU = 1D-10
      IF(IFLAGS(IFIPFC).EQ.5.OR.IFLAGS(IFIPFC).EQ.6) ARGU = ALPI**2*
     &     ( 4D0/9D0*ALFE**3
     &     + 2D0/3D0*( 31D0/18D0 - 2D0*DZ2 )*ALFE
     &     - TEE**3/18D0 + 19D0/36D0*TEE**2
     &     + ( PI**2/9D0 - 265D0/108D0 )*TEE
     &     + 4D0/9D0*ALFM**3
     &     + 2D0/3D0*( 31D0/18D0 - 2D0*DZ2 )*ALFM
     &     - TMU**3/18D0 + 19D0/36D0*TMU**2
     &     + ( PI**2/9D0 - 265D0/108D0 )*TMU
     &     + 4D0/9D0*ALFT**3
     &     + 2D0/3D0*( 31D0/18D0 - 2D0*DZ2 )*ALFT
     &     - TTU**3/18D0 + 19D0/36D0*TTU**2
     &     + ( PI**2/9D0 - 265D0/108D0 )*TTU
     &     )  
 
      AJMST= GEXP
     &       *DEXP( 3D0/4D0*BETE 
     &            + ALPI*( 2D0*DZ2 - 0.5D0 ) 
     &            )
     &       *PHI
     &       *DEXP(ARGU)
     &       *( DLL + DNLL )
C.... PRAGMATIC FORMULA FOR PHOTONIC RC:
      FBBE = 1D0 - BETE**2/2D0*DZ2 + BETE**3/3D0*DZ3
      PRAGM= BETE*GEXP*FBBE
     &       *DEXP( 3D0/4D0*BETE + ALPI*(2D0*DZ2-0.5D0) )
     &       *( (1D0+Z**2)/2D0 + ALPI*BETE*( 3D0/32D0 - 3D0/4D0*DZ2
     &        + 3D0/2D0*DZ3 ) + BETE/8D0*( - (1D0+3D0*Z**2)*DLZ
     &        - 2D0*(1D0-Z)**2 ) + ALPI/8D0*( 4D0*(1D0+Z**2)*( SP1Z
     &        + DL1Z*DLZ ) - (1D0+3D0*Z**2)*DLZ**2 + 2D0*( 3D0 + 2D0*Z
     &        + Z**2 )*DLZ + 2D0*(1D0-Z)*(3D0-2D0*Z) )
     &        + BETE**2/8D0*( 1D0/12D0*(1D0+7D0*Z**2)*DLZ**2
     &        + (1D0-Z)/2D0*(1D0-3D0*Z)*DLZ + (1D0-Z)**2
     &        + (1D0-Z**2)*SP1Z ) )
*
      CALL BORN(IFINAL,Z,Z,SBORN,DUMMY,SBORNS,ABORNS)
*
      AJMSB= (AJMST-PRAGM)*SBORN
      RETURN
      END

      DOUBLE PRECISION FUNCTION RRR(S)
C     --------------------------------
C.... PARAMETRIZATION OF R-FUNCTION FOR ANNIHILATION INTO HADRONS
C     BETA-VERSION
C     AUTHOR: A. ARBUZOV /arbuzov@to.infn.it/
C     S = (E_CM)^2   [GeV^2]
C     R = SIGMA(E+E- --> HADRONS)/SIGMA(E+E- --> MUONS)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NDA=98,AMPI=139.56995D-03)
      DIMENSION EP(NDA),RP(NDA)
      DATA EP/ 
     &  0.280, 0.300, 0.350, 0.400, 0.450, 0.500,
     &  0.550, 0.600, 0.650, 0.700, 0.720, 0.740,
     &  0.760, 0.780, 0.800, 0.810, 0.830, 0.850,
CAC  &  0.280, 0.300, 0.350, 
C    &  0.360, 0.380, 0.410, 0.430, 0.438, 0.470,
C    &  0.540, 0.580, 0.620, 0.640, 0.660, 0.700,
C    &  0.740, 0.760, 0.770, 0.774, 0.778, 0.780,
C    &  0.782, 0.786, 0.790, 0.794, 0.800, 0.820,
CAC  &                              0.830, 0.850,
     &  0.870, 0.890, 0.910, 0.930, 0.950, 0.970,
     &  0.990, 1.050, 1.070, 1.090, 1.110, 1.130,
     &  1.150, 1.170, 1.190, 1.210, 1.230, 1.250,
     &  1.270, 1.290, 1.310, 1.330, 1.350, 1.370,
     &  1.390, 1.450, 1.470, 1.500, 1.540, 1.570,
     &  1.600, 1.630, 1.650, 1.670, 1.700, 1.800,
     &  1.900, 2.000, 2.180, 2.380, 2.600, 2.700,
     &  2.900, 3.000, 3.200, 3.400, 3.600, 3.750,
     &  3.750, 3.800, 3.850, 3.900, 3.920, 3.950,
     &  4.000, 4.050, 4.100, 4.200, 4.300, 4.350,
     &  4.400, 4.450, 4.500, 4.550, 4.600, 4.700,
     &  4.800, 6.000, 6.600, 7.000, 7.200,10.000,
     & 12.500,25.000,35.000,40.000,45.000,50.000,
     & 55.000,60.000 /
C
      DATA RP/ 
     &  0.163, 0.179, 0.277, 0.405, 0.617, 0.752,
     &  1.107, 1.697, 2.530, 5.590, 7.057, 8.734,
     &  9.945,10.686, 6.428, 5.545, 4.204, 3.219,
CAC  &  0.163, 0.179, 0.277,
C    &  0.118, 0.166, 0.232, 0.269, 0.317, 0.362,
C    &  0.699, 0.949, 1.538, 1.798, 2.343, 4.892,
C    &  9.028, 8.226, 8.895, 8.808, 9.554, 9.162,
C    &  6.689, 5.698, 6.251, 5.908, 5.563, 4.033,
CAC  &                              4.204, 3.219,
     &  2.516, 1.958, 1.617, 1.355, 1.134, 1.018,
     &  0.930, 1.492, 1.245, 1.118, 1.056, 0.998,
     &  0.995, 1.090, 1.140, 1.021, 1.264, 1.320,
     &  1.257, 1.424, 1.486, 1.622, 1.774, 1.789,
     &  1.877, 1.470, 1.590, 2.000, 2.300, 2.200,
     &  2.000, 2.100, 1.800, 2.180, 2.050, 1.700,
     &  1.600, 1.650, 2.000, 2.600, 2.800, 2.900,
     &  2.230, 2.450, 2.800, 2.400, 2.500, 2.600,
     &  3.600, 4.600, 2.600, 3.000, 3.200, 4.000,
     &  4.200, 5.600, 5.000, 4.000, 3.800, 4.000,
     &  5.000, 4.800, 4.000, 3.800, 3.600, 4.000,
     &  4.000, 4.200, 4.600, 4.350, 4.400, 4.000,
     &  4.000, 3.800, 3.750, 3.850, 4.000, 4.400,
     &  4.600, 5.000 /
CAC  &  4.800, 6.500 /

      RRR  = RP(NDA)
      SS   = DSQRT(S)

CAC   DO I3= 1,98
C      PRINT 10,EP(I3),RP(I3)
C 10   FORMAT(2F12.3)
C     END DO
CAC   STOP

      IF(SS.GE.EP(98)) RETURN
      IF(SS.LT.2D0*AMPI) THEN
       RRR = 0D0
       RETURN
      ENDIF
      IF(SS.LE.EP(1)) THEN
       RRR = RP(1)
       RETURN
      ENDIF
      I4   = 1
CAC   DO WHILE(SS.GE.EP(I4))
C      I4  = I4 + 1
C     ENDDO
CAC   DO I4 = 1,109
      DO I4 = 1,98
       IF(SS.LT.EP(I4)) GOTO 11
      ENDDO
 11   CONTINUE  
 12   RRR  = RP(I4-1) + (RP(I4)-RP(I4-1))*(SS-EP(I4-1))
     &      /(EP(I4)-EP(I4-1))
CAC
      IF(I4.LE.15) RRR = RRR*1.2D0
CAC     
      RETURN
      END

      FUNCTION FSRPR(S,ZMIN,AMPRC)
C     ----------------------------
C.... FINAL STATE PAIRS IN E^+E^- ANNIHILATION
C     ZMIN*S  CUT ON THE INVARIANT MASS OF THE PRIMARY PAIR
C     AMPRC*S CUT ON THE INVARIANT MASS OF THE SECONDARY PAIR
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
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /PCONST/ ALFAI,AL1PI,ALQE2,ALQF2,ALQEF,GMU,CSIGNB
      COMMON /MASSZ / AME,AMF,AME2,AMF2
      COMMON /CHARGZ/ QE,QF,QEM,QFM,QEF,QEFM,QE2,QF2
      COMMON /CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
      COMMON /VARY / YY,AM1S,AM2S,AM2PR,SS,ALM1
      EXTERNAL HKT,VHP,RHP
      PARAMETER (REPS=1D-5,RAPS=1D-10,AMPI = 139.56995D-03)
      FSRPR = 0D0
C
C----------------------------------------------------      
C.... PREPARATIONS
      SS    = S
C.... PRIMARY PAIR:
      AM1   = AMF
      AM1S  = AM1*AM1/S
      ALM1  = DLOG(AM1S)
      AM2PR = AMPRC
 
C.... ELECTRON SECONDARY PAIR
      CORRE = 0D0
      VIRTE = 0D0
      IF(IFLAGS(IFIPFC).EQ.1.OR.IFLAGS(IFIPFC).GE.5) THEN
       AM2  = AME
       AM2S = AM2*AM2/S
       YN   = MAX(4D0*AM1S,ZMIN)
       YX   = (1D0-2D0*AM2/SQRT(S))**2
       YST  = 0.25D0*(YX-YN)
       RES  = 0D0
       IF(YST.GT.0D0.AND.IFLAGS(IFIPTO).NE.-1) 
     &   CALL SIMPS(YN,YX,YST,REPS,RAPS,HKT,TTT,RES,RES2,RES3)
       CORRE= 4D0/3D0*RES
       VIRTE= VHKT(AM2S)
      ENDIF
C.... MUON SECONDARY PAIR
      CORRM = 0D0
      VIRTM = 0D0
      IF(IFLAGS(IFIPFC).EQ.2.OR.IFLAGS(IFIPFC).GE.5) THEN
       AM2  = AML(4)
       AM2S = AM2*AM2/S
       YN   = MAX(4D0*AM1S,ZMIN)
        YX   = (1D0-2D0*AM2/SQRT(S))**2
       YST  = 0.25D0*(YX-YN)
       RES  = 0D0
       IF(YST.GT.0D0.AND.IFLAGS(IFIPTO).NE.-1) 
     &   CALL SIMPS(YN,YX,YST,REPS,RAPS,HKT,TTT,RES,RES2,RES3)
       CORRM= 4D0/3D0*RES
       VIRTM= VHKT(AM2S)
      ENDIF
C.... HADRONIC SECONDARY PAIR
      CORRH = 0D0
      VIRTH = 0D0
      IF(IFLAGS(IFIPFC).EQ.4.OR.IFLAGS(IFIPFC).EQ.5) THEN
       XN   = ZMIN
       XX   = (1D0 - 2D0*AMPI/SQRT(S))**2
       XST  = 0.0025D0*(XX-XN)
       RESX = 0D0
       IF(XST.GT.0D0.AND.IFLAGS(IFIPTO).NE.-1) 
     &   CALL SIMPS(XN,XX,XST,REPS,RAPS,RHP,TTT,RESX,RES2,RES3)
       CORRH= RESX*S
       Q2N  = 4D0*AMPI**2
       Q2X  = 1D8
       Q2ST = 0.00025D0*(Q2X-Q2N)
       CALL SIMPS(Q2N,Q2X,Q2ST,REPS,RAPS,VHP,TTT,RES,RES2,RES3)
       VIRTH= RES*2D0
      ENDIF

C----------------------------------------------------      
      FSRPR = (CORRE + VIRTE + CORRM + VIRTM + CORRH + VIRTH)
     &        *ALQE2**2*QF2
CAC
C     PRINT *,'EE',CORRE,VIRTE/2D0
C     PRINT *,'EE',CORRE*ALQE2**2,VIRTE*ALQE2**2
C     PRINT *,'MU',CORRM,VIRTM/2D0
C     PRINT *,'MU',CORRM*ALQE2**2,VIRTM*ALQE2**2
C     PRINT *,'HA',CORRH,VIRTH/2D0
C     PRINT *,'HA',CORRH*ALQE2**2,VIRTH*ALQE2**2
C     STOP
CAC
      RETURN
      END

      DOUBLE PRECISION FUNCTION HKT(Y)
C     --------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /VARY / YY,AM1S,AM2S,AM2PR,SS,ALM1
      EXTERNAL HKT1
      PARAMETER (REP=1D-6,RAP=1D-15)
      HKT   = 0D0
      YY    = Y 
       
      ZN    = 4D0*AM2S
      ZX    = MIN((1D0-SQRT(Y))**2,AM2PR)
      ZST   = 0.25D0*(ZX-ZN)
      IF(ZST.LE.0D0) RETURN
      CALL SIMPT(ZN,ZX,ZST,REP,RAP,HKT1,TTT,RES1,RES2,RES3)
      HKT   = RES1
      RETURN
      END

      DOUBLE PRECISION FUNCTION HKT1(Z)
C     ---------------------------------
C.... DIFFERENTIAL DDISTRIBUTION IN INVARIANT MASS
C     ACCORDING TO: A.H.HOANG J.H.KUHN T.TEUBNER
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /VARY / YY,AM1S,AM2S,AM2PR,SS,ALM1
      HKT1  = 0D0
      Y     = YY

      ALAM  = 1D0 + Y**2 + Z**2 - 2D0*(Y+Z+Y*Z)
      IF(ALAM.LE.0D0) RETURN
      SALA  = SQRT(ALAM)
      DEN1  = 1D0 - Y + Z
      IF(DEN1.GE.-1D-12.AND.DEN1.LE.0D0) DEN1 = - 1D-12
      IF(DEN1.LE. 1D-12.AND.DEN1.GT.0D0) DEN1 =   1D-12
      ARG1  = ( 1D0 - Y + Z - SQRT(1D0-4D0*AM1S/Y)*SALA )
     &      / ( 1D0 - Y + Z + SQRT(1D0-4D0*AM1S/Y)*SALA )
      HKT1  = 1D0/Z*(1D0+2D0*AM2S/Z)*SQRT(1D0-4D0*AM2S/Z)*(
     &        ( 2D0*AM2S**2 + AM1S*(1D0-Y+Z) - 0.25D0*(1D0-Y+Z)**2
     &         - 0.5D0*(1D0+Z)*Y )/DEN1*LOG(ABS(ARG1))
     &      - SQRT(1D0-4D0*AM1S/Y)*SALA*( 0.25D0  + ( 2D0*AM1S 
     &        + 4D0*AM1S**2 + (1D0+2D0*AM1S)*Z )/( (1D0-Y+Z)**2
     &        - (1D0-4D0*AM1S/Y)*ALAM ) ) 
     &                                                    )
      RETURN
      END

      DOUBLE PRECISION FUNCTION VHKT(X)
C     ---------------------------------
C.... VIRTUAL PAIR RC (AM1S,AM2S << 1) (KEEP ONLY F1)
C     ACCORDING TO: A.H.HOANG J.H.KUHN T.TEUBNER, EQ.(33)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /VARY / YY,AM1S,AM2S,AM2PR,SS,ALM1
      VHKT  = 0D0
      ALX   = DLOG(X)
      IF(DABS(X-AM1S)/X.LE.1D-5) THEN
        VHKT= ALX**3/18D0 + 19D0/36D0*ALX**2 + 2D0/3D0*( 265D0/72D0
     &        - F1 )*ALX - 11D0/3D0*F1 + 383D0/54D0
      ELSEIF(X.LT.AM1S*0.9999D0) THEN
        VHKT= ALX**2/6D0*(ALM1+1D0) + ALX*( - ALM1**2/6D0 
     &        + 13D0/18D0*ALM1 - 4D0/3D0*F1 + 11D0/9D0 )
     &        + ALM1**3/18D0 - 13D0/36D0*ALM1**2 
     &        + ( 133D0/108D0 + 2D0/3D0*F1 )*ALM1
     &        - 2D0/3D0*ZET3 - 32D0/9D0*F1 + 67D0/27D0
     &        + 3D0/4D0*X/AM1S*PI**2 
      ELSEIF(X.GT.AM1S*1.0001D0) THEN
        VHKT= ALX**3/18D0 + 19D0/36D0*ALX**2 + 2D0/3D0*( 265D0/72D0
     &        - F1 )*ALX - 2D0/3D0*ZET3 - 19D0/9D0*F1 + 3355D0/648D0
      ELSE
        PRINT *,'PLEASE CHECK MASSES:'
        PRINT *,'X*S=  ',X*SS
        PRINT *,'AM1=  ',AM1S*SS
        STOP
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION VHP(Q2)
C     ---------------------------------
C.... FUNCTION UNDER INTEGRAL
      IMPLICIT REAL*8(A-H,O-Y)
      IMPLICIT COMPLEX*16(Z)
      COMMON /VARY / YY,AM1S,AM2S,AM2PR,SS,ALM1
      ARG   = SS/Q2
      VHP   = 1D0/3D0/Q2*RRR(Q2)*DREAL(ZRHO(ARG))
      RETURN
      END

      FUNCTION ZRHO(U)
C     ----------------
C.... RHO-FUNCTION FOR VERTEX CORRECTION
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 ZI,ZLMU,ZRHO
      COMMON /NCONST/ PI,F1,AL2,ZET3
      ZRHO  = DCMPLX(0D0,0D0)
      IF(U.EQ.0D0) RETURN
      ZI    = DCMPLX(0D0,1D0)
      ZLMU  = DCMPLX(DLOG(DABS(U)),0D0)
      IF(U.GT.0D0) ZLMU = ZLMU - ZI*PI
      ZRHO  = - 7D0/8D0 - 1D0/2D0/U + (3D0/4D0+1D0/2D0/U)*ZLMU
     &      - 0.5D0*(1D0+1D0/U)**2*( SPENCE(-U) + ZLMU*LOG(1D0+U) )
      END

      DOUBLE PRECISION FUNCTION RHP(X)
C     --------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /VARY / X0,AM1S,AM2S,AM2PR,SS,ALM1
      PARAMETER (REP=1D-6,RAP=1D-15,AMPI=139.56995D-03)
      EXTERNAL RHP1
      RHP   = 0D0
      X0    = X

      Q2N   = 4D0*AMPI**2
      Q2X   = MIN(SS*(1D0-SQRT(X))**2,AM2PR*SS)
      Q2ST  = 0.0025D0*(Q2X-Q2N)
      IF(Q2ST.LE.0D0) RETURN
      CALL SIMPT(Q2N,Q2X,Q2ST,REP,RAP,RHP1,TT,RESH,RES2,RES3)
      RHP   = RESH
      RETURN
      END

      DOUBLE PRECISION FUNCTION RHP1(Q2)
C     ----------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /NCONST/ PI,F1,AL2,ZET3
      COMMON /VARY / X0,AM1S,AM2S,AM2PR,SS,ALM1
      RHP1  = 0D0
      P2    = X0*SS
  
      ARG   = SS**2 + P2**2 + Q2**2 - 2D0*SS*P2
     &           - 2D0*P2*Q2 - 2D0*Q2*SS 
      IF(ARG.LE.0D0) RETURN
      SQAL  = SQRT(ARG)
      RHP1  = RRR(Q2)/3D0/SS**2/Q2*(
     &        ( SS**2 + (P2+Q2)**2 )/( SS - P2 - Q2)
     &        *DLOG( (SS - P2 - Q2 + SQAL)/(SS - P2 - Q2 - SQAL) )
     &        - 2D0*SQAL )
      RETURN
      END
