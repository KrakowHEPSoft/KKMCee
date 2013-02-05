      SUBROUTINE formkuehn(s,FPI)
*-----------------------------------------------------------------------
* Here is the subroutine which returns the Kuehn-Santamaria-form-factor
* F_\pi(s) for a given \pi\pi invariant mass s (Z. Phys. C45 (1990) 445):
*-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 s
      REAL*8 alpha,beta,gamma,mrho,mrho1,mrho2,Gammarho,Gammarho1
      REAL*8 Gammarho2,momega,Gammaomega
      COMPLEX*16 BWrho,BWrho1,BWrho2,BWomega
      COMPLEX*16 FPI

      alpha = 1.85d-3
      beta  = -0.145d0
      gamma = 0.d0

      mrho  = 773.d-3
      mrho1 = 1370.d-3
      mrho2 = 1700.d-3

      Gammarho  = 145.d-3
      Gammarho1 = 510.d-3
      Gammarho2 = 235.d-3
 
      momega   = 0.78194d0
      Gammaomega = 8.43d-3

      BWrho  = mrho**2/(mrho**2-s-(0.d0,1.d0)*sqrt(s)*Gammarho)
      BWrho1 = mrho1**2/(mrho1**2-s-(0.d0,1.d0)*sqrt(s)*Gammarho1)
      BWrho2 = mrho2**2/(mrho2**2-s-(0.d0,1.d0)*sqrt(s)*Gammarho2)
      BWomega = momega**2/(momega**2-s
     &         - (0.d0,1.d0)*sqrt(s)*Gammaomega)

      FPI = (BWrho*(1.d0+alpha*BWomega)/(1.d0+alpha)
     &     + beta*BWrho1+gamma*BWrho2)/(1.d0+beta+gamma)
      END
  
      DOUBLE PRECISION  FUNCTION F_Pi_Kuehn_SQ(s)
      IMPLICIT NONE
      REAL*8  s
      COMPLEX*16 F_pi
*
      CALL formkuehn(s,F_pi)
      F_Pi_Kuehn_SQ = F_pi * CONJG(F_pi)
      END
