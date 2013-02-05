      SUBROUTINE formkuehn02(s,FPI)
c KLOE 2002 data - best fit (G. Venanzoni) unpublished? 
      IMPLICIT NONE
      REAL*8 s 
      REAL*8 alpha,beta,gamma,mrho,mrho1,mrho2,Gammarho,Gammarho1
      REAL*8 Gammarho2,momega,Gammaomega
      COMPLEX*16 BWrho,BWrho1,BWrho2,BWomega
      COMPLEX*16 FPI
c old 1999
c old: alpha = 1.85d-3
c old: beta  = -0.145d0 
c old: mrho  = 773.d-3
c old: mrho1 = 1370.d-3
c old: Gammarho  = 145.d-3 
c old: Gammarho1 = 510.d-3
c old: momega   = 0.78194d0
c old: Gammaomega = 8.43d-3
c [masses and widths in GeV]
      alpha = 1.48d-3           ! alpha = (1.48 +- 0.12) * 10^-3
      beta = -0.1473d0          ! beta =  -0.1473 +- 0.002
      gamma = 0.d0              ! gamma put to 0 in parameterization as before
      mrho = 772.6d-3           ! mrho = 772.6 +- 0.5 MeV   
      mrho1 = 1460.d-3          ! mrho1 = 1.46 GeV
      mrho2 = 1700.d-3          ! old value - but anyway multiplied by gamma=0
      Gammarho = 143.7d-3       ! Gammarho = 143.7 +- 0.7 MeV 
      Gammarho1 = 310.d-3       ! Gammarho1 = 0.31 GeV
      Gammarho2 = 235.d-3       ! old value - but anyway multiplied by gamma=0
      momega = 782.78d-3        ! momega = 0.78278 GeV
      Gammaomega = 8.68d-3      ! Gammaomega = 8.68 * 10^-3 GeV 
      BWrho  = mrho**2/(mrho**2-s-(0.d0,1.d0)*sqrt(s)*Gammarho)
      BWrho1 = mrho1**2/(mrho1**2-s-(0.d0,1.d0)*sqrt(s)*Gammarho1)  
      BWrho2 = mrho2**2/(mrho2**2-s-(0.d0,1.d0)*sqrt(s)*Gammarho2)
      BWomega = momega**2/(momega**2-s
     &         - (0.d0,1.d0)*sqrt(s)*Gammaomega)
      FPI = (BWrho*(1.d0+alpha*BWomega)/(1.d0+alpha)
     &     + beta*BWrho1+gamma*BWrho2)/(1.d0+beta+gamma)
      END
