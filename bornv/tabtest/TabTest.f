      PROGRAM tabtest
*     ***************
* Testing BornV class which provides Born(s,theta) from tables
* >     make tabtest
* >     make tabtest-debug
*----------------------------------------------------------------------------
      IMPLICIT none

      INTEGER KFdown, KFup, KFstran, KFcharm, KFbotom, KFtop
      PARAMETER( 
     $     KFdown  = 1,   KFup    = 2,
     $     KFstran = 3,   KFcharm = 4,
     $     KFbotom = 5,   KFtop   = 6)
      INTEGER KFel,KFelnu,KFmu,KFmunu,KFtau,KFtaunu
      PARAMETER(
     $     KFel    = 11,  KFelnu  = 12,
     $     KFmu    = 13,  KFmunu  = 14,
     $     KFtau   = 15,  KFtaunu = 16)

      DOUBLE PRECISION    ene

      DOUBLE PRECISION    seps1,seps2
      DOUBLE PRECISION    BornV_Make
      DOUBLE PRECISION    BornV_Differential

      DOUBLE PRECISION    ta,tb,CMSene,CosTheta,bor1,bor2
      INTEGER  ninp, ninp2, nout
      INTEGER  i,imax
      DOUBLE PRECISION    born0,bornr,svar

      CHARACTER*40 fmt,fmt1,fmt2
      PARAMETER ( fmt  ='(   i5, 2f5.1, 9g19.13)')
      PARAMETER ( fmt1 ='(   i2, 2f5.1, 9g19.13)')
      PARAMETER ( fmt2 ='(i2,i5, 2f5.1, 9g19.13)')

      CHARACTER*1 chr1
      INTEGER mode
      CHARACTER*40 TableFile, BenchFile

      INTEGER  KFini,KFfin
      DOUBLE PRECISION   amz,amh,amtop

      INTEGER jmax,Itest1,Itest2,iloop
      PARAMETER (jmax = 1000)
      DOUBLE PRECISION      xpar(jmax)
*----------------------------------------------------------------------------
      Itest1 =1
      Itest2 =1

      ninp  =7
      ninp2 =9
      nout =16

      CALL fort_OPEN(  ninp,'./TabTest.input')
      CALL ReaDataX(   ninp,xpar,jmax)
      CALL fort_CLOSE( ninp)

      OPEN(nout,FILE='./TabTest.output')

* beam params
      CMSene =  200d0
      ene    =  CMSene/2d0
      KFini  =  KFel

      amz   = 91.187d0
      amh   = 100d0
      amtop = 175d0

      seps1=0d0
      seps2=0d0

* Initializations
      CALL BornV_Initialize(xpar)
*---------------------- test1 ------------------------------------
      IF(Itest1.EQ.1) THEN
* Here we compare differ. xsections from Koralz stored in *.benchmark
* files with the ones provided by BornV_Differential
* It is done separately for each flavour.
* EW corrections are off and on.

      WRITE(*,*) ' -------------- test1 --------------'
      DO iloop =1,5
* note that present KORALZ benchmarks were done with negative KF's
* i.e. positron along z-axis
      IF(iloop.EQ.1) KFfin  = KFdown
      IF(iloop.EQ.2) KFfin  = KFup
      IF(iloop.EQ.3) KFfin  = KFbotom
      IF(iloop.EQ.4) KFfin  = KFtau
      IF(iloop.EQ.5) KFfin  = KFmu
*
      IF(    ABS(KFfin) .EQ. KFmu) THEN
         BenchFile= './benchmark.mu'
      ELSEIF(    ABS(KFfin) .EQ. KFtau) THEN
         BenchFile= './benchmark.tau'
      ELSEIF( ABS(KFfin) .EQ. KFdown) THEN
         BenchFile= './benchmark.down'
      ELSEIF( ABS(KFfin) .EQ. KFup) THEN
         BenchFile= './benchmark.up'
      ELSEIF( ABS(KFfin) .EQ. KFbotom) THEN
         BenchFile= './benchmark.botom'
      ELSE
         WRITE(*,*) 'TabTest: Wrong KFfin = ',KFfin
         STOP
      ENDIF

      OPEN(ninp2,FILE=BenchFile)
      READ(ninp2,*) chr1
      READ(ninp2,*) chr1

      WRITE(*,*) ' '
      WRITE(*,*) 'Comparison with results in: ', BenchFile

      imax=25
      DO i=1,imax
* read benchmark from koralz
         READ(ninp2,fmt) mode,ta,tb,CMSene,CosTheta,bor1,bor2
         svar = CMSene**2
* calculate and check x-section at several svar and CosTheta
* true KFini,KFfin
***         bornr= BornV_Make( mode,KFini,KFfin,svar,CosTheta,seps1,seps2,ta,tb)
* reversed KFini,KFfin, result is the same!!!
***         bornr= BornV_Make( mode,-KFini,-KFfin,svar,CosTheta,seps1,seps2,-ta,-tb)
* KFini defined in input data ! TabTest.input
         CALL BornV_SetKeyElw(1) ! reset KeyElw
         bornr= BornV_Differential( mode,KFfin,svar,CosTheta,seps1,seps2,ta,tb)
         WRITE( *,fmt1)  mode,ta,tb,CMSene,CosTheta, bor1,bor2, bor1/bornr, bor2/bornr
      ENDDO
      CLOSE(ninp2)
      ENDDO
      ENDIF ! Itest1

*---------------------- test2 ------------------------------------
* This test is done inclusively: in ./inclusive.benchmark we merged
* entries from previous test1 with various KFfin (KF is added in first column)
* EW corrections are OFF. Final state helicities ta,tb are zero.
* We compare Born without EW as implemented in old way without spin
* amplitudes and new way with spin amplitudes, the difference is
* only (?) due to different treatment of mass terms and threshold bevavior.
      IF(Itest2.EQ.1) THEN
      WRITE(*,*) ' -------------- test2 --------------'
      OPEN(ninp2,FILE='./benchmark.inclusive')
      READ(ninp2,*) chr1
      READ(ninp2,*) chr1
      imax=32
      DO i=1,imax
         READ(ninp2,fmt2) KFfin,mode,ta,tb,CMSene,CosTheta,bor1,bor2
         svar = CMSene**2
* Old Born0 as emulated in new one
         CALL BornV_SetKeyElw(0) ! KeyElw=0 is ancient formula
         born0 = BornV_Differential( 0,KFfin,svar,CosTheta,0d0,0d0,0d0,0d0)
* BornV
         CALL BornV_SetKeyElw(1) ! KeyElw=1 is new Born from spin ampls.
         Mode = 1 ! for Mode = 1 we check how big are EW corrections
         bornr = BornV_Differential( mode,KFfin,svar,CosTheta,seps1,seps2,ta,tb)
* Born0 versus BornV
***         WRITE(  *,fmt2) KFfin,mode,ta,tb,CMSene,CosTheta,bornr,born0,bornr/born0
***         WRITE(  *,fmt2) KFfin,mode,ta,tb,CMSene,CosTheta,bor1, born0, bor1/born0
* BornV versus KORALZ benchmark
         WRITE(  *,fmt2) KFfin,mode,ta,tb,CMSene,CosTheta,bor1,bornr,bor1/bornr
      ENDDO
      CLOSE(ninp2)
      ENDIF
*-----------------------------------------------------------------
      END

