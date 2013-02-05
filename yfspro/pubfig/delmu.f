      PROGRAM MAIN
*****************************
* This prepared file delmu.data for Delmu-plot
*     gmake delmu-data
*     gmake delta-ps
*****************************

      IMPLICIT NONE

*------------------------------------------------------------------------------
      INTEGER    IdGenYFS3
      PARAMETER (IdGenYFS3 = 6)
*------------------------------------------------------------------------------
      INTEGER isigG1            !GPS
      INTEGER isigG1nin         !GPS int OFF
      INTEGER isigG1int         !GPS Interf. only
      INTEGER isigG0,isigG0nin,isigG0int  ! GPS zero order version
*------------------------------------------------------------------------------
      INTEGER nout,nbiv
*------------------------------------------------------------------------------
      INTEGER idi,idc,mdc,idb,mdb,kdb,ibin,mds,ids
      PARAMETER(
     $  idc =  10000,
     $  mdc = 110000,
     $  idb =  50000,
     $  mdb = 150000,
     $  idi = 300000,
     $  kdb = 450000,
     $  ids  = 40000,
     $  mds  =140000)
*------------------------------------------------------------------------------
*  SEMIREALISTIC  v-like-distributios, all cumulative
      INTEGER         isigS1, isigS1nin, isigS1int
      PARAMETER(
     $  isigS1          = mds+ 202,        !sig(vmax) CEEX O1 full
     $  isigS1nin       = mds+ 252,        !sig(vmax) CEEX O1 int OFF
     $  isigS1int       = mds+ 282)        !sig(vmax) CEEX O1 interf. only
*
      DOUBLE PRECISION  SigG1, dSigG1, SigG1nin, dSigG1nin
      PARAMETER(
     $  isigG1          = mdb+ 202,        !sig(vmax) GPS O1
     $  isigG1nin       = mdb+ 252,        !sig(vmax) GPS O1 int OFF
     $  isigG1int       = mdb+ 282,        !sig(vmax) GPS O1 interf. only
     $  isigG0          = mdb+ 201,        !sig(vmax) GPS O0
     $  isigG0nin       = mdb+ 251,        !sig(vmax) GPS O0 int OFF
     $  isigG0int       = mdb+ 281)        !sig(vmax) GPS O0 interf. only
*------------------------------------------------------------------------------
      DOUBLE PRECISION  GLK_hi,GLK_hie
      DOUBLE PRECISION  SigG1Int, dSigG1Int
      DOUBLE PRECISION  MZ,CMSene,dCMSene
      INTEGER loop
      CHARACTER*60  Hname,Zname,Zname2,Kname
      INTEGER      i,j,k, Nb, npoint
      REAL*8       xbin(400),xerr(400),ybin(400),yerr(400),xdel(400),ydel(400)
      REAL*8       SigZf,SigZfn,SigZfInt
      REAL*8       SigKmInt,dSigKmInt
*------------------------------------------------------------------------------
      nout= 16
      OPEN( nout, file='output-delmu')
      MZ =91.187d0

      OPEN( 13,File='delmu.data')
      WRITE( *,*) ' Table of delta-sigma-total [pb] due to ISR*FSR interference '
      WRITE( *,*) 
     $ '  CMSene-MZ,     SigG1,    SigZfInt,   SigG1Int,  dSigG1Int,  SigKmInt,  dSigKmInt'
      WRITE(13,*) ' Table of delta-sigma-total due to ISR*FSR interference '
      WRITE(13,*) 
     $ '  CMSene-MZ,     SigG1,    SigZfInt,   SigG1Int,  dSigG1Int,  SigKmInt,  dSigKmInt'

      npoint = 5                !! number of energy data points
      WRITE(13,*) npoint, '   points'
      k=npoint/2

      DO loop=-k,k
         IF(    loop.EQ.-2) THEN
            dCMSene = -3d0
            Hname  = '../E91GeV/E91-3GeV.hst'
            Zname  = 'ZfitterMZ-3'
            Zname2 = 'ZfitterMZ-3noint'
            Kname  = 'MustraMZ-3'
         ELSEIF(loop.EQ.-1) THEN
            dCMSene = -1.8d0
            Hname  = '../E91GeV/E91-1.8GeV.hst'
            Zname  = 'ZfitterMZ-1.8'
            Zname2 = 'ZfitterMZ-1.8noint'
            Kname  = 'MustraMZ-1.8'
         ELSEIF(loop.EQ.0) THEN
            dCMSene = + 0d0
            Hname  = '../E91GeV/E91+0GeV.hst'
            Zname  = 'ZfitterMZ'
            Zname2 = 'ZfitterMZnoint'
            Kname  = 'MustraMZ'
         ELSEIF(loop.EQ.+1) THEN
            dCMSene = + 1.8d0
            Hname  = '../E91GeV/E91+1.8GeV.hst'
            Zname  = 'ZfitterMZ+1.8'
            Zname2 = 'ZfitterMZ+1.8noint'
            Kname  = 'MustraMZ+1.8'
         ELSEIF(loop.EQ.+2) THEN
            dCMSene = + 3d0
            Hname  = '../E91GeV/E91+3GeV.hst'
            Zname  = 'ZfitterMZ+3'
            Zname2 = 'ZfitterMZ+3noint'
            Kname  = 'MustraMZ+3'
         ENDIF
******************** Extract MC data ********************
         CALL GLK_ReadFile(Hname)
*
         CALL GLK_CumHis(     IdGenYFS3, idb+ 201, isigG0)    ! GPS O(alf0)
         CALL GLK_CumHis(     IdGenYFS3, idb+ 251, isigG0nin) ! GPS O(alf0) int OFF
         CALL GLK_CumHis(     IdGenYFS3, idb+ 281, isigG0int) ! GPS O(alf0)
*
         CALL GLK_CumHis(     IdGenYFS3, idb+ 202, isigG1)    ! GPS O(alf1)
         CALL GLK_CumHis(     IdGenYFS3, idb+ 252, isigG1nin) ! GPS O(alf1) int OFF
         CALL GLK_CumHis(     IdGenYFS3, idb+ 282, isigG1int) ! GPS O(alf1)
* SEMIREALISTIC
         CALL GLK_CumHis(     IdGenYFS3, ids+ 202, isigS1)    ! CEEX O(alf0)
         CALL GLK_CumHis(     IdGenYFS3, ids+ 252, isigS1nin) ! CEEX O(alf0) int OFF
         CALL GLK_CumHis(     IdGenYFS3, idb+ 282, isigS1int) ! GPS O(alf0)
* vmax = 0.99 ==> ibin=99
         ibin =90
         ibin =10
         ibin =01
         ibin =40
         ibin =50
         ibin =30
* Standard, with loose cut
         ibin =99
         SigG1     = GLK_hi(  isigG1   ,ibin)*1d3
         dSigG1    = GLK_hie( isigG1   ,ibin)*1d3
         SigG1nin  = GLK_hi(  isigG1nin,ibin)*1d3
         dSigG1nin = GLK_hie( isigG1nin,ibin)*1d3
         SigG1Int  = GLK_hi(  isigG1int,ibin)*1d3
         dSigG1Int = GLK_hie( isigG1int,ibin)*1d3
* Semirealistic, for test
* s'>0.64s, vmax=0.36 ==> ibin=36, ALEPH
         ibin =30
         SigG1     = GLK_hi(  isigS1   ,ibin)*1d3
         dSigG1    = GLK_hie( isigS1   ,ibin)*1d3
         SigG1nin  = GLK_hi(  isigS1nin,ibin)*1d3
         dSigG1nin = GLK_hie( isigS1nin,ibin)*1d3
         SigG1Int  = GLK_hi(  isigS1int,ibin)*1d3
         dSigG1Int = GLK_hie( isigS1int,ibin)*1d3
* test: O(alf0) instead O(alf1)
c         SigG1     = GLK_hi(  isigG0   ,ibin)*1d3
c         dSigG1    = GLK_hie( isigG0   ,ibin)*1d3
c         SigG1nin  = GLK_hi(  isigG0nin,ibin)*1d3
c         dSigG1nin = GLK_hie( isigG0nin,ibin)*1d3
c         SigG1Int  = GLK_hi(  isigG0int,ibin)*1d3
c         dSigG1Int = GLK_hie( isigG0int,ibin)*1d3
*********************************************************
         CALL GLK_GetNb(isigG1,Nb)
         CALL RData1(Zname ,Nb, xbin, xerr, ybin, yerr)
         SigZf   = xbin(ibin)
         CALL RData1(Zname2,Nb, xbin, xerr, ybin, yerr)
         SigZfn  = xbin(ibin)
         SigZfInt = SigZf-SigZfn
         CALL  RData4(Kname,Nb, xbin, xerr, ybin, yerr, xdel,ydel)
         SigKmInt  = xdel(ibin)*1d3
         dSigKmInt = xerr(ibin)*1d3
*********************************************************
         WRITE(13,'(20f16.6)') dCMSene, SigG1, SigZfInt, SigG1Int, dSigG1Int, SigKmInt, dSigKmInt!
         WRITE( *,'(20f16.6)') dCMSene, SigG1, SigZfInt, SigG1Int, dSigG1Int, SigKmInt, dSigKmInt!
         WRITE( *,'(20g16.6)') dCMSene, SigG1, SigZfInt, SigG1Int, dSigG1Int, SigKmInt, dSigKmInt!
         WRITE( *,'(20g16.6)') dCMSene, SigG1, SigZfInt/SigG1, SigG1Int/SigG1, dSigG1Int/SigG1!
*********************************************************
         CALL GLK_Flush
      ENDDO
      CLOSE(13)
      END

      SUBROUTINE   RData1(DataFile1,Nb,xbin,xerr,ybin,yerr)
*/////////////////////////////////////////////////////////////////////////////
*// reading sigma and  from diskfile, error is zero                         //
*// data from Zfiter                                                        //
*/////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      CHARACTER*(*)     DataFile1
      DOUBLE PRECISION  xbin(*),xerr(*),ybin(*),yerr(*)
      INTEGER           Nb,k,i,npoint,ibin
      DOUBLE PRECISION  vmax,sig,afb
      CHARACTER*4       ch

***      WRITE(*,*) '======== RData1 will read diskfile ',DataFile1
      DO k=1,Nb
         xbin(  k) =  1d9
         xerr(  k) =  0d0
         ybin(  k) =  1d9
         yerr(  k) =  0d0
      ENDDO
      OPEN(20,File=DataFile1)
      READ(20,*) ch
      READ(20,*) ch
      READ(20,*) ch
      READ(20,*) ch
      READ(20,*) ch
      READ(20,*) npoint
      DO i=1,npoint
         READ(20,*) vmax,sig,afb
         ibin = vmax*100d0
         IF(ibin.LT.0 .OR. ibin.GT.Nb) GOTO 900
         xbin(ibin) = sig
         ybin(ibin) = afb
***         WRITE(*,*) ibin,xbin(ibin),ybin(ibin)
      ENDDO
      CLOSE(20)
      RETURN
 900  WRITE(*,*) ' STOP in  RData1 !!!!'
      STOP
      END

      SUBROUTINE   RData4(DataFile1,Nb,xbin,xerr,ybin,yerr,xdel,ydel)
*/////////////////////////////////////////////////////////////////////////////
*// reading sigma and AFB from diskfile, data from KoralZ/Mustraal          //
*/////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      CHARACTER*(*)     DataFile1
      DOUBLE PRECISION  xbin(*),xerr(*),ybin(*),yerr(*),xdel(*),ydel(*)
      INTEGER           Nb,k,i,npoint,ibin
      DOUBLE PRECISION  vmax,sig,dsig,afb,dafb,sigint,afbint
      CHARACTER*26      ch

      WRITE(*,*) '======== RData4 will read diskfile ',DataFile1
      DO k=1,Nb
         xbin(  k) =  1d9
         xerr(  k) =  0d0
         ybin(  k) =  1d9
         yerr(  k) =  0d0
         xdel(  k) =  1d9
         ydel(  k) =  1d9
      ENDDO
      OPEN(20,File=DataFile1)
      READ(20,*) ch
      READ(20,*) ch
      READ(20,*) ch
      READ(20,*) ch
      READ(20,*) ch
      READ(20,*) npoint
      DO i=1,npoint
         READ(20,*) vmax,  sig,sigint,dsig,  afb,afbint,dafb
         ibin = NINT(vmax*100d0)
         IF(ibin.LT.0 .OR. ibin.GT.Nb) GOTO 900
         xbin(ibin) =  sig
         xerr(ibin) = dsig
         ybin(ibin) =  afb
         yerr(ibin) = dafb
         xdel(ibin) = sigint
         ydel(ibin) = afbint
         WRITE(*,*) vmax,sig,dsig,afb,dafb,sigint,afbint,ibin
      ENDDO
      CLOSE(20)
      RETURN
 900  WRITE(*,*) ' STOP in  RData1 !!!!'
      STOP
      END

