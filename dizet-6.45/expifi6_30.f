      SUBROUTINE SXPIFI(C,SIFI)
C     -------------------------
C.... EXPONENTIATION OF INITIAL-FINAL PHOTON RC
C     RC TO INTEGRATED X-SECTION
C     TEST VERSION
C     A.ARBUZOV 20.03.2000
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL SIFI1
      SIFI  = 0D0
      IF(C.LE.0D0) RETURN
      C0    = C
      IF(C.GE.1D0-1D-4) C0 = 1D0 - 1D-4
      REPS  = 1D-6
      RAPS  = REPS*REPS
      CN    = - C0
      CX    = + C0
      CST   = 0.25D0*(CX-CN)
      CALL SIMPS(CN,CX,CST,REPS,RAPS,SIFI1,UU,RESS,AIHH,AIHA) 
      SIFI  = RESS/(C+C**3/3D0)
      RETURN
      END
      DOUBLE PRECISION FUNCTION SIFI1(CC)
C     -----------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /PCONST/ ALFAI,AL1PI,ALQE2,ALQF2,ALQEF,GMU,CSIGNB
      COMMON /CUTVAR/ SINAC2,COSAC2,RCUT,RACUT,RECUT1,RECUT2,RECUTA,PCUT
      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR 
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
      COMMON /BETAC / BETACC,BETAEF,CCC
      COMMON /PSCONS/ SW2,AMZ,GAMZ
      DATA EPSIFI/1.0D-7/
      EXTERNAL SIFI2
      REPS  = 1D-7
      RAPS  = 1D-15

      CCC   = CC

      ALEPSI= DLOG(EPSIFI)
      BETACC= BETTAE + 4D0*ALQEF*DLOG((1D0-CC)/(1D0+CC))
      BETAEF= BETTAE
      IF(IFINAL.EQ.1) THEN
       BETAEF = BETAEF + BETTAF
       BETACC = BETACC + BETTAF
      ENDIF

      YN    = ALEPSI
      YX    = DLOG(1D0-RCUT)
      YST   = 0.25D0*(YX-YN)
      IF(YST.LE.0D0) THEN
       PRINT *,'RCUT=',RCUT,' IS NOT GOOD FOR EXPIFI.F'
       STOP
      ENDIF
      CALL SIMPT(YN,YX,YST,REPS,RAPS,SIFI2,RR,RES4,AIHH,AIHA)

      CALL BORN(IFINAL,1D0,1D0,SBORN,ABORN,SBORNS,ABORNS)
      RES2  = DEXP(BETACC*ALEPSI)
      RES3  = DEXP(BETAEF*ALEPSI)
      RES4D = ( BETACC - BETAEF )*ALEPSI
      SIFI1 = ( RES2 - RES3 - RES4D )
     &     *0.5D0*( SBORN*(1D0+CC**2) + ABORN*2D0*CC )
     &      + RES4

      RETURN
      END
      DOUBLE PRECISION FUNCTION SIFI2(Y)
C     ----------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /SFTVAR/ SCOM,TE,TEE,TMU,TTU,TF,BETTAE,BETTAF,SOFTER,SOFTFR
      COMMON /PCONST/ ALFAI,AL1PI,ALQE2,ALQF2,ALQEF,GMU,CSIGNB
      COMMON /FLAGZ / IAFB,IBORN,IRCUT,IFINAL,INTERF,IWEAK,IPHOT2,ISYM
      COMMON /BETAC / BETACC,BETAEF,CCC
      COMMON /PSCONS/ SW2,AMZ,GAMZ
      R     = 1D0 - DEXP(Y)

CAC   CALL BORN(IFINAL,1D0,1D0,SBORN,ABORN,SBORNS,ABORNS)
      CALL BORN(IFINAL,R,1D0,SBORN,ABORN,SBORNS,ABORNS)
  
      SIFI2 = 0.5D0*( SBORN*(1D0+CCC**2) + ABORN*2D0*CCC )
     &      *( 0D0
     &       + BETACC*DEXP(BETACC*Y)
     &       - BETAEF*DEXP(BETAEF*Y)
     &       - BETACC + BETAEF 
     &       ) 
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE AXPIFI(C,AIFI)
C     -------------------------
C.... EXPONENTIATION OF INITIAL-FINAL PHOTON RC
C     RC TO INTEGRATED ASYMMETRY
C     A.ARBUZOV 20.03.2000
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL SIFI1
      AIFI  = 0D0
      IF(C.LE.0D0) RETURN
      C0    = C
      IF(C.GE.+1D0-1D-4) C0 =+1D0 - 1D-4
      IF(C.LE.-1D0+1D-4) C0 =-1D0 + 1D-4
      REPS  = 1D-6
      RAPS  = REPS*REPS
      CN    = - C0
      CX    = + C0
      CST   = 0.125D0*(CX-CN)
      CALL SIMPS(0D0,CX,CST,REPS,RAPS,SIFI1,UU,RESF,AIHH,AIHA)
      CALL SIMPS(CN,0D0,CST,REPS,RAPS,SIFI1,UU,RESB,AIHH,AIHA)
      AIFI  = (RESF-RESB)/C**2
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE CXPIFI(C,CIFI)
C     -------------------------
C.... EXPONENTIATION OF INITIAL-FINAL PHOTON RC
C     RC TO DIFFERENTIAL X-SECTION
C     TEST VERSION
C     A.ARBUZOV 20.03.2000
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL SIFI1
      C0    = C 
      IF(C.GE.+1D0-1D-4) C0 =+1D0 - 1D-4
      IF(C.LE.-1D0+1D-4) C0 =-1D0 + 1D-4
      CIFI  = 2D0*SIFI1(C0)
C.... THE FACTOR 2 ABOVE IS DUE TO 1/2 IN SIFI1,2
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------