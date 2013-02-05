      SUBROUTINE for_init(nout,fname)
!     **********************************
! -------------------------------------------------------
      CHARACTER fname*(*)
      CHARACTER*100 filename

      filename = fname
      nout2 = nout
!     WRITE(6,'(A,A20,A)') '======>',filename,'<========='
      OPEN(nout,file=filename)

      END

      SUBROUTINE for_fini(nout)
!     *************************
! -------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 

      nout2 = nout
      CLOSE(nout2)

      END

      SUBROUTINE koralw_out(qq1,qq2,nnphot,ssphot,pp1,pp2,pp3,pp4)
!     ***************************************************************
! Initialization of internal histogramming package GLIBK
! -------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DIMENSION qq1(4),qq2(4),ssphot(4,100)
      DIMENSION pp1(4),pp2(4),pp3(4),pp4(4)

      COMMON / momset / qeff1(4),qeff2(4),sphum(4),sphot(100,4),nphot 
      COMMON / momdec / q1(4),q2(4),p1(4),p2(4),p3(4),p4(4)
      COMMON / decays / iflav(4), amdec(4), br(2), brel

      nnphot = nphot
      DO iph = 1,nphot
         DO k=1,4
            ssphot(k,iph) = sphot(iph,k)
         ENDDO
      ENDDO
         DO k=1,4
            qq1(k) = q1(k)
            qq2(k) = q2(k)
            pp1(k) = p1(k)
            pp2(k) = p2(k)
            pp3(k) = p3(k)
            pp4(k) = p4(k)
         ENDDO
      END

      SUBROUTINE lujets_out(n2,k2,p2,v2)
!     **********************************
      DIMENSION k2(5,4000),p2(5,4000),v2(5,4000)
!...Purpose: to store one parton/particle in commonblock LUJETS. 
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)

      n2 = n
      DO i=1,4000
         DO j = 1,5
            k2(j,i) = k(i,j)
            p2(j,i) = p(i,j)
            v2(j,i) = v(i,j)
         ENDDO
      ENDDO
      END

      SUBROUTINE ludat1_out(mu,pu,mj,pj)
!     **********************************
      DIMENSION mu(200),pu(200),mj(200),pj(200)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 

      DO i=1,200
         mu(i)=MSTU(i)
         pu(i)=PARU(i)
         mj(i)=MSTJ(i)
         pj(i)=PARJ(i)
      ENDDO
      END

      SUBROUTINE mylugive(chin)
!     *************************
      CHARACTER chin*(*)
      CHARACTER*100 ch

      ch = chin
      WRITE(6,'(A,A20,A)') '======>',chin,'<========='
      WRITE(6,'(A,A20,A)') '======>',ch,'<========='

!      CALL lugive('MSTU(11)=9')
      CALL lugive(ch)

      END
