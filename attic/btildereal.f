
      function Btilde(p1,p2,am1,am2,aKmax,amgam)
*     ******************************************
!----------------------------------------------------------------------!
! This function provides a value of YFS real photon IR function        !
! B-tilde for any pair of charged particles.                           !
! INPUT: p1,p2   - particles 4-momenta;                                !
!        am1,am2 - particles masses;                                   !
!        amgam   - "photon mass" (IR regulator)                        !
!        aKmax   - maximum soft photon energy [GeV]                    !
!----------------------------------------------------------------------!
! Written by:  Wieslaw Placzek                Knoxville, January  1996 !  
! Last update: 30.01.1996                by: W.P.                      !
!----------------------------------------------------------------------! 
      IMPLICIT REAL*8(a-h,o-z)
      PARAMETER ( pi = 3.1415926535897932D0, alfinv = 137.0359895d0 )
      PARAMETER ( alfpi = 1/alfinv/pi )
      REAL*8 p1(4),p2(4)
!
      Btilde = 0
      E1 = p1(4)
      E2 = p2(4)
      am12 = am1*am2
      p1p2 = p1(4)*p2(4) - p1(3)*p2(3) - p1(2)*p2(2) - p1(1)*p2(1)
      IF (p1p2-am12.lt.1d-10) RETURN
      xlam = SQRT( (p1p2 - am12)*(p1p2 + am12) )
! Function A(p1,p2)
      A  = 1/xlam *LOG( (p1p2 + xlam)/am12 )
      bet1 = SQRT(1-am1**2/E1**2)
      bet2 = SQRT(1-am2**2/E2**2)
      b1ln = 2*LOG( (1+bet1)*E1/am1 )
      b2ln = 2*LOG( (1+bet2)*E2/am2 )
! Function A4(p1,p2)
      A4 = A4anal(p1p2,E1,E2,am1,am2)
! B-tilde(p1,p2;aKmax,amgam)
      Btian = (p1p2*A - 1) *2*LOG(2*aKmax/amgam)
     &      + 0.5/bet1*b1ln + 0.5/bet2*b2ln + p1p2*A4
      Btilde = alfpi*Btian
      END

      function A4anal(p1p2,En1,En2,xm1,xm2)
*     *************************************
!----------------------------------------------------------------------!
! This function provides an analytical result for the integral         !
! A4(p1,p2) being a part of the YFS IR function B-tilde.               !
! Note: This is a general case without any approximation!              !
! INPUT: p1p2    - scalar product of the 4-momenta p1 and p2;          !
!        E1,E2   - particles energies;                                 !
!        xm1,xm2 - particles masses;                                   !
!----------------------------------------------------------------------!
! Written by:  Wieslaw Placzek                Knoxville, January  1996 !  
! Last update: 30.01.1996                by: W.P.                      !
!----------------------------------------------------------------------! 
      IMPLICIT REAL*8(a-h,o-z)
! Statement function
      etaln(x1,x2,x3,x4,z) = LOG(ABS( (z-x1)*(z-x2)/(z-x3)/(z-x4) ))
! Some auxiliary variables
      E1 = En1
      E2 = En2
      am1 = xm1
      am2 = xm2
      p1s = E1**2 - am1**2
      p2s = E2**2 - am2**2
      IF (p1s.lt.p2s) THEN
        am1 = xm2
        am2 = xm1
        E1 = En2
        E2 = En1
      ENDIF
      Ep  = E1 + E2
      Em  = E1 - E2
      sm  = am1 + am2 
      dm  = am1 - am2
      Q2  = 2*p1p2 - am1**2 - am2**2
      xl  = SQRT( (Q2 + sm**2)*(Q2 + dm**2) )
      xq  = SQRT(Q2 + Em**2)
      qp = xq + Em
      qm = xq - Em
      et0 = SQRT(E2**2 - am2**2)
      IF (p1p2.gt.E1*E2) et0 = -et0
      et1 = SQRT(E1**2 - am1**2) + xq
      y1  = 0.5*( (xq - Ep) + (sm*dm + xl)/qp )
      y2  = y1 - xl/qp
      y3  = 0.5*( (xq + Ep) + (sm*dm + xl)/qm )
      y4  = y3 - xl/qm       
! Some auxiliary functions
      IF (ABS(Em).gt.1d-10) THEN
        Eln = LOG(ABS(qm/qp))*( etaln(y1,y4,y2,y3,et1) 
     &                        - etaln(y1,y4,y2,y3,et0) )
      ELSE
        Eln = 0
      ENDIF
      Vet0 = Yijeta(y1,y4,et0) + Yijeta(y2,y1,et0)
     &     + Yijeta(y3,y2,et0) - Yijeta(y3,y4,et0)
     &     + 0.5*etaln(y1,y2,y3,y4,et0)*etaln(y2,y3,y1,y4,et0)
      Vet1 = Yijeta(y1,y4,et1) + Yijeta(y2,y1,et1)
     &     + Yijeta(y3,y2,et1) - Yijeta(y3,y4,et1)
     &     + 0.5*etaln(y1,y2,y3,y4,et1)*etaln(y2,y3,y1,y4,et1)
! Function A4(p1,p2) 
      A4anal = 1/xl*(Eln + Vet1 - Vet0 )
      END

      function Yijeta(yi,yj,eta)
*     **************************
!----------------------------------------------------------------------!
! Some auxiliary function (combination of Logs and Dilogs) used in     !
! the function A4anal for A4(p1,p2).                                   !
!----------------------------------------------------------------------!
! Written by:  Wieslaw Placzek                Knoxville, January  1996 !  
! Last update: 30.01.1996                by: W.P.                      !
!----------------------------------------------------------------------! 
      IMPLICIT REAL*8(a-h,o-z)
!
      Yijeta = 2*DILOGY( (yj-yi)/(eta-yi) ) 
     &       + 0.5*LOG(ABS( (eta-yi)/(eta-yj) ))**2
      END
 

