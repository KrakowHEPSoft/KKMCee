*----------------------------------------------------------------
* 	gmake bornex
*----------------------------------------------------------------
      PROGRAM MAIN
*     ************
      IMPLICIT NONE
      INTEGER   ninp,nout
      CHARACTER*60      Tesnam, Dname
      DOUBLE PRECISION  BornPb,QCDcor,Born
      DOUBLE PRECISION  QCDcorR(10)
*------------------------------------------------------------------------------
      INTEGER    imax,i,j,k
      PARAMETER (imax=10000)
      REAL*8     xpar(imax)
      REAL*8     ymin,ymax
*--------------------------------------------------------------------
      ninp=  5
      nout= 16
      Tesnam    = 'bornex'
      OPEN( nout, file='output-'//Tesnam)
      CALL GLK_SetNout(nout)

      Dname  = './bornex.input'
      CALL KK2f_ReaDataX('../../.KK2f_defaults', 1,imax,xpar)  ! Read data, the same as in MC run
      CALL KK2f_ReaDataX(                 Dname, 0,imax,xpar)  ! Read user input
      CALL KK2f_Initialize( xpar)                    ! Initialize generator with the production data
      CALL KKsem_Initialize(xpar)                    ! Initialize semianalytical package
*=================================================================================================
      WRITE(   *,*) '-------------------------------------------------------------' !
      CALL KKsem_GetBorn(Born)  ! This re-sets QCDcor
      BornPb = Born*1000d0
      WRITE(   *,'(a,f15.7)') 'Bornex:: Born [Pb], QCD on = ',BornPb
      WRITE(nout,'(a,f15.7)') 'Bornex:: Born [Pb], QCD on = ',BornPb
      CALL BornV_GetQCDcor(QCDcor)
      WRITE(   *,'(a,f15.7)') 'Bornex:: QCDcor= ',QCDcor
      WRITE(nout,'(a,f15.7)') 'Bornex:: QCDcor= ',QCDcor
      CALL BornV_SetKeyQCD(0)
      CALL KKsem_GetBorn(Born)
      BornPb = Born*1000d0
      WRITE(   *,'(a,f15.7)') 'Bornex:: Born [Pb], QCD off= ',BornPb
      WRITE(nout,'(a,f15.7)') 'Bornex:: Born [Pb], QCD off= ',BornPb
      WRITE(   *,*) '-------------------------------------------------------------' !
      CALL KKsem_GetBorn(Born)  ! This re-sets QCDcor
      CALL BornV_GetQCDcorR(QCDcorR)
      WRITE(   *,'(a,9f15.7)') 'Bornex:: QCDcorR= ',(QCDcorR(i),i=1,4)
      WRITE(nout,'(a,9f15.7)') 'Bornex:: QCDcorR= ',(QCDcorR(i),i=1,4)
      WRITE(   *,*) '-------------------------------------------------------------' !
      CLOSE(nout)
      END
