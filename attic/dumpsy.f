
      SUBROUTINE dumps(nout)
*     **********************
* This prints out ALL four momenta on unit no. nout
*     **********************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON / momset / qf1(4),qf2(4),sphum(4),sphot(100,4),nphot
      SAVE   / momset /
      DIMENSION sum(4)
      WRITE(nout,*) '=============ALL====MOMSET===================='
      WRITE(nout,3100) 'qf1',(  qf1(  k),k=1,4)
      WRITE(nout,3100) 'qf2',(  qf2(  k),k=1,4)
      DO i=1,nphot
         WRITE(nout,3100) 'pho',(sphot(i,j),j=1,4)
      ENDDO
      DO k=1,4
         sum(k)=qf1(k)+qf2(k)
      ENDDO
      DO i=1,nphot
         DO k=1,4
            sum(k)=sum(k)+sphot(i,k)
         ENDDO
      ENDDO
      WRITE(nout,3100) 'sum',(  sum(  j),j=1,4)
 3100 FORMAT(1x,a3,1x,5f20.14)
      END

      SUBROUTINE dumpf(nout)
*     **********************
* This prints out four momenta of FSR photons on unit no. nout
*     **********************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON / momfin / qf1(4),qf2(4),sphum(4),sphot(100,4),nphot
      SAVE   / momfin /
      DIMENSION sum(4)
      WRITE(nout,*) '============FSR=====MOMFIN===================='
      WRITE(nout,3100) 'qf1',(  qf1(  k),k=1,4)
      WRITE(nout,3100) 'qf2',(  qf2(  k),k=1,4)
      DO i=1,nphot
         WRITE(nout,3100) 'pho',(sphot(i,k),k=1,4)
      ENDDO
      DO k=1,4
         sum(k)=qf1(k)+qf2(k)
      ENDDO
      DO i=1,nphot
         DO k=1,4
            sum(k)=sum(k)+sphot(i,k)
         ENDDO
      ENDDO
      WRITE(nout,3100) 'sum',(  sum(  k),k=1,4)
 3100 FORMAT(1x,a3,1x,5f20.14)
      END

      SUBROUTINE dumpi(nout)
*     **********************
* This prints out four momenta of ISR photons on unit no. nout
*     **********************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON / momini / qf1(4),qf2(4),sphum(4),sphot(100,4),nphot
      SAVE   / momini /
      DIMENSION sum(4)
      WRITE(nout,*) '============ISR=====MOMINI===================='
      WRITE(nout,3100) 'qf1',(  qf1(  k),k=1,4)
      WRITE(nout,3100) 'qf2',(  qf2(  k),k=1,4)
      DO i=1,nphot
         WRITE(nout,3100) 'pho',(sphot(i,k),k=1,4)
      ENDDO
      DO k=1,4
         sum(k)=qf1(k)+qf2(k)
      ENDDO
      DO i=1,nphot
         DO k=1,4
            sum(k)=sum(k)+sphot(i,k)
         ENDDO
      ENDDO
      WRITE(nout,3100) 'sum',(  sum(  k),k=1,4)
 3100 FORMAT(1x,a3,1x,5f20.14)
      END

      SUBROUTINE dumpn(nout,iev)
*     **************************
* this prints out four momenta of final state
* and the serial number of event iev on unit nout
*     **********************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON / momini / qf1(4),qf2(4),sphum(4),sphot(100,4),nphot
      SAVE   / momini /
      DIMENSION sum(4)
      WRITE(nout,*) '======ISR=====MOMINI====================>',iev
      WRITE(nout,3100) 'qf1',(  qf1(  k),k=1,4)
      WRITE(nout,3100) 'qf2',(  qf2(  k),k=1,4)
      DO i=1,nphot
         WRITE(nout,3100) 'pho',(sphot(i,k),k=1,4)
      ENDDO
      DO k=1,4
         sum(k)=qf1(k)+qf2(k)
      ENDDO
      DO i=1,nphot
         DO k=1,4
            sum(k)=sum(k)+sphot(i,k)
         ENDDO
      ENDDO
      WRITE(nout,3100) 'sum',(  sum(  k),k=1,4)
 3100 FORMAT(1x,a3,1x,5f20.14)
      END
