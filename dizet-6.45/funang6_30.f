      DOUBLE PRECISION FUNCTION FUNANC(Z1,Z2,C)
C     -----------------------------------------
C.... PART OF CORRECTED ANGULAR DISTRIBUTION OF SECOND ORDER
C     LLA PHOTONIC CORRECTIONS
C     DIFFERENTIAL IN C
      IMPLICIT REAL*8(A-H,O-Z)
      A      = Z1 + Z2 - C*(Z1-Z2)
      CTI    = ( Z2 - Z1 + C*(Z1+Z2) )/A
      DCTI   = 4D0*Z1*Z2/A**2
      FUNANC = (1D0+CTI**2)*DCTI
      RETURN
      END

      DOUBLE PRECISION FUNCTION FUNANG(Z1,Z2,C)
C     -----------------------------------------
C.... PART OF CORRECTED ANGULAR DISTRIBUTION OF SECOND ORDER
C     LLA PHOTONIC CORRECTIONS
C     (INTEGRATED OVER C)
      IMPLICIT REAL*8(A-H,O-Z)
      AX   = Z1 + Z2 + C*(Z1-Z2)
      AN   = Z1 + Z2 - C*(Z1-Z2)
      Z1PZ2= Z1 + Z2
      FUNANG =
     +    8*z1/an - 8*z1/ax - 24*z1**2*c/an*
     +    z1pz2**(-1) - 24*z1**2*c/ax/z1pz2 - 16.*z1**2*
     +    an**(-2) + 16*z1**2/ax**2 + 64*z1**3*c/an**2*
     +    z1pz2**(-1) + 64*z1**3*c/an/z1pz2**2 + 64.*z1**3*c*
     +    ax**(-2)/z1pz2 + 64*z1**3*c/ax/z1pz2**2 - 64.*
     +    z1**3*c/z1pz2**3 + 32*z1**3*c**2/an/z1pz2**2 - 32.
     +    *z1**3*c**2/ax/z1pz2**2 + 64./3.*z1**3/an**3 - 64.
     +    /3.*z1**3/ax**3 - 64*z1**4*c/an**3/z1pz2 - 64*
     +    z1**4*c/an**2/z1pz2**2 - 64*z1**4*c/an*
     +    z1pz2**(-3) - 64*z1**4*c/ax**3/z1pz2 - 64*z1**4*c*
     +    ax**(-2)/z1pz2**2 - 64*z1**4*c/ax/z1pz2**3 + 160*
     +    z1**4*c/z1pz2**4 - 80*z1**4*c**2/an**2/z1pz2**2 - 
     +    160*z1**4*c**2/an/z1pz2**3 + 80*z1**4*c**2/ax**2*
     +    z1pz2**(-2) + 160*z1**4*c**2/ax/z1pz2**3 - 16*z1**4*
     +    c**3/an/z1pz2**3 - 16*z1**4*c**3/ax*
     +    z1pz2**(-3) 
     +    + 64*z1**5*c**2/an**3/z1pz2**2 + 128*z1**5*
     +    c**2/an**2/z1pz2**3 + 192*z1**5*c**2/an*
     +    z1pz2**(-4) - 64*z1**5*c**2/ax**3/z1pz2**2 - 128*z1**5*
     +    c**2/ax**2/z1pz2**3 - 192*z1**5*c**2/ax*
     +    z1pz2**(-4) + 32*z1**5*c**3/an**2/z1pz2**3 + 96*z1**5*
     +    c**3/an/z1pz2**4 + 32*z1**5*c**3/ax**2*
     +    z1pz2**(-3) + 96*z1**5*c**3/ax/z1pz2**4 - 64./3.*
     +    z1**6*c**3/an**3/z1pz2**3 - 64*z1**6*c**3/an**2*
     +    z1pz2**(-4) - 128*z1**6*c**3/an*z1pz2**(-5) - 64./3.*
     +    z1**6*c**3/ax**3/z1pz2**3 - 64*z1**6*c**3/ax**2*
     +    z1pz2**(-4) - 128*z1**6*c**3/ax*z1pz2**(-5) 
      RETURN
      END

      FUNCTION P1(X)
C     --------------
C.... P1 SPLITTING FUNCTION (THETA-PART)
      IMPLICIT REAL*8(A-H,O-Z)
      P1   = ( 1D0 + X**2 )/( 1D0 - X )
      RETURN
      END

      FUNCTION P1Y(Y)
C     ---------------
C.... P1 SPLITTING FUNCTION (THETA-PART)
C     Y = LOG(1-X)    <-- EXCHANGED VARIABLES
      IMPLICIT REAL*8(A-H,O-Z)
      X    = 1D0 - DEXP(Y)
      P1Y  = 1D0 + X**2
      RETURN
      END

      FUNCTION P10(ALD)
C     -----------------
C.... P1 SPLITTING FUNCTION (DELTA-PART)
      IMPLICIT REAL*8(A-H,O-Z)
      P10  = 2D0*ALD + 1.5D0
      RETURN
      END

      FUNCTION P2(X)
C     --------------
C.... P2 SPLITTING FUNCTION (THETA-PART)
      IMPLICIT REAL*8(A-H,O-Z)
      P2   = 2D0     
     &       *( (1D0+X**2)/(1D0-X)*( 2D0*DLOG(1D0-X) - DLOG(X)
     &       + 3D0/2D0 ) + 1D0/2D0*(1D0+X)*DLOG(X) - 1D0 + X )
      RETURN
      END

      FUNCTION P20(ALD)
C     -----------------
C.... P2 SPLITTING FUNCTION (DELTA-PART)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /NCONST/ PI,F1,AL2,ZET3
      P20  = ( 2D0*ALD + 1.5D0 )**2 - 4D0*F1
      RETURN
      END

      DOUBLE PRECISION FUNCTION FUNAN2(Z2)
C     ------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     &                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM
      Z1   = R/Z2
      FUNAN2 = FUNANG(Z1,Z2,C)*P1(Z2)*P1(Z1)/Z2
     &     - 2D0/(1D0-Z2)*P1(R)*FUNANG(R,1D0,C)
     &     - 2D0/R/(1D0-R/Z2)*P1(R)*FUNANG(1D0,R,C)
      RETURN
      END

      DOUBLE PRECISION FUNCTION FUNAN3(Z2)
C     ------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /RVAR  / R,R2,R3,R4,R5,R1MI,R1MI2,R1MI3,R1MI4,ALR
      COMMON /CVAR  / C,COPL3,C2,C3,CP,CM,CP2,CM2,CP3,CM3,CPM,CPM2,CPM3
     &                 ,COPL2,ALCP,ALCM,ALPL,ALMI,DLCP,DLCM
      Z1   = R/Z2
      FUNAN3 = FUNANC(Z1,Z2,C)*P1(Z2)*P1(Z1)/Z2
     &     - 2D0/(1D0-Z2)*P1(R)*FUNANC(R,1D0,C)
     &     - 2D0/R/(1D0-R/Z2)*P1(R)*FUNANC(1D0,R,C)
      RETURN
      END
