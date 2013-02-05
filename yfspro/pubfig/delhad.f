      PROGRAM MAIN
*****************************
*     gmake delhad-data
*****************************

      IMPLICIT NONE

*------------------------------------------------------------------------------
      INTEGER    IdGenYFS3
      PARAMETER (IdGenYFS3 = 6)
*------------------------------------------------------------------------------
      INTEGER isigG1            !GPS
      INTEGER isigG1nin         !GPS int OFF
      INTEGER isigG1int         !GPS Interf. only
*------------------------------------------------------------------------------
      INTEGER idi,idc,mdc,idb,mdb,nbiv
      INTEGER nout,kdb,ibin
*------------------------------------------------------------------------------
      PARAMETER(
     $  idc =  10000,
     $  mdc = 110000,
     $  idb =  50000,
     $  mdb = 150000,
     $  idi = 300000,
     $  kdb = 450000,
*
     $  isigG1          = mdb+ 202,        !sig(vmax) GPS O1
     $  isigG1nin       = mdb+ 252,        !sig(vmax) GPS O1 int OFF
     $  isigG1int       = mdb+ 282         !sig(vmax) GPS O1 interf. only
     $ )
*------------------------------------------------------------------------------
      DOUBLE PRECISION  GLK_hi,GLK_hie
      DOUBLE PRECISION  SigG1, dSigG1, SigG1nin, dSigG1nin
      DOUBLE PRECISION  SigG1Int, dSigG1Int
      DOUBLE PRECISION  MZ,CMSene,dCMSene
      INTEGER loop
      CHARACTER*60  Hname
      INTEGER      i,j,k, Nb, npoint
      REAL*8       xbin(400),xerr(400),ybin(400),yerr(400),xdel(400),ydel(400)
*------------------------------------------------------------------------------
      nout= 16
      OPEN( nout, file='output-delhad')
      MZ =91.187d0

      OPEN( 13,File='delhad.data')
      WRITE( *,*) ' Table of delta-sigma-total [pb] due to ISR*FSR interference '
      WRITE( *,*) '  CMSene-MZ,     SigG1,    SigG1Int,  dSigG1Int'
      WRITE(13,*) ' Table of delta-sigma-total due to ISR*FSR interference '
      WRITE(13,*) '  CMSene-MZ,     SigG1,    SigG1Int,  dSigG1Int'

      npoint = 5                !! number of energy data points
      WRITE(13,*) npoint, '   points'
      k=npoint/2

      DO loop=-k,k
         SigG1     = -1d9
         dSigG1    =  0d0
         SigG1nin  = -1d9
         dSigG1nin =  0d0
         SigG1Int  = -1d9
         dSigG1Int =  0d0
         IF(loop.EQ. -2) THEN
            dCMSene = -3d0
            Hname  = '../E91GeV/H91-3GeV.hst'
         ELSEIF(loop.EQ. -1) THEN
            dCMSene = -1.8d0
            Hname  = '../E91GeV/H91-1.8GeV.hst'
            GOTO 100 !!!!!<<--- no data
         ELSEIF(loop.EQ. 0) THEN
            dCMSene = + 0d0
            Hname  = '../E91GeV/H91+0GeV.hst'
         ELSEIF(loop.EQ. +1) THEN
            dCMSene = + 1.8d0
            Hname  = '../E91GeV/H91+1.8GeV.hst'
            GOTO 100 !!!!!<<--- no data
         ELSEIF(loop.EQ. +2) THEN
            dCMSene = + 3d0
            Hname  = '../E91GeV/H91+3GeV.hst'
         ENDIF
******************** Extract MC data ********************
         CALL GLK_ReadFile(Hname)
         CALL GLK_CumHis(     IdGenYFS3, idb+ 202, isigG1) ! GPS O(alf1)
         CALL GLK_CumHis(     IdGenYFS3, idb+ 252, isigG1nin) ! GPS O(alf1) int OFF
*     CALL GLK_Print(isigG1)
*     CALL GLK_Print(isigG1nin)
         CALL GLK_CumHis(     IdGenYFS3, idb+ 282, isigG1int) ! GPS O(alf1)
* vmax = 0.99 ==> ib=99
         ibin =99
         SigG1     = GLK_hi(  isigG1   ,ibin)*1d3
         dSigG1    = GLK_hie( isigG1   ,ibin)*1d3
         SigG1nin  = GLK_hi(  isigG1nin,ibin)*1d3
         dSigG1nin = GLK_hie( isigG1nin,ibin)*1d3
         SigG1Int  = GLK_hi(  isigG1int,ibin)*1d3
         dSigG1Int = GLK_hie( isigG1int,ibin)*1d3
         CALL GLK_Flush
 100     CONTINUE
*********************************************************
         WRITE(13,'(20g12.4)') dCMSene, SigG1, SigG1Int,       dSigG1Int
         WRITE( *,'(20g12.4)') dCMSene, SigG1, SigG1Int,       dSigG1Int
         WRITE( *,'(20g12.4)') dCMSene, SigG1, SigG1Int/SigG1, dSigG1Int/SigG1
*********************************************************
      ENDDO
      CLOSE(13)

      END
