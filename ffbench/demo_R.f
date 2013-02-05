      PROGRAM main
*//////////////////////////////////////////////////////////////////////
*//                                                                  //
*//  Very Simple demonstration program with loop over events         //
*//                                                                  //
*//  To execute:                                                     //
*//      make demo-start                                             //
*//                                                                  //
*//////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
* input/output files
      INTEGER         ninp,nout
      COMMON /c_demo/ ninp,nout
      CHARACTER*4 semaph
      INTEGER     ntot2n, ntotin, ninp2, ninph, ijklin, ninp3
      SAVE
*---------------
* general output for everybody including glibk
      ninp =5
      nout =16
      OPEN(nout,FILE='./demo.output')
      REWIND(nout)
      CALL GLK_SetNout(nout)
*---------------
* READ semaphore flag
      CALL givsem(semaph)
*---------------
      IF(semaph  .EQ.  'STAR') THEN
         WRITE(6,*) ' ------- starting from the scratch ----------'
* READ initial (root) random number seed
         ninp3=3
         OPEN(ninp3,FILE='./iniseed')
         READ(ninp3,'(i10)') ijklin
         READ(ninp3,'(i10)') ntotin
         READ(ninp3,'(i10)') ntot2n
         CALL PseuMar_Initialize(ijklin,ntotin,ntot2n)
      ELSEIF(semaph  .EQ.  'CONT') THEN
         WRITE(6,*) ' ------- restoring from the disk   ----------'
* restore histograms from the disk
         ninph=10
         OPEN(ninph,FILE='./'//'pro.hst')
         CALL GLK_hrfile(ninph,' ',' ')      !transfer FILE number
         CALL GLK_hrin(   0,9999,0)          !READ from the disk
         CALL GLK_hrend(' ')                 !CLOSE FILE
* READ random number seed stored in semaphore FILE
         ninp2=2
         OPEN(ninp2,FILE='./semaphore')
         READ(ninp2,'(a4)') semaph
         READ(ninp2,'(i10)') ijklin
         READ(ninp2,'(i10)') ntotin
         READ(ninp2,'(i10)') ntot2n
         CALL PseuMar_Initialize(ijklin,ntotin,ntot2n)
         CLOSE(ninp2)
      ELSEIF(semaph  .EQ.  'STOP') THEN
         WRITE(6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++'
         WRITE(6,*) '++++ STOP: please change semaph to cont or start !'
         WRITE(6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++'
         STOP
      ELSE
         WRITE(6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++'
         WRITE(6,*) '++++ STOP: wrong key semaph = ', semaph
         WRITE(6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++'
         STOP
      ENDIF
* *******************
      CALL splot
* *******************
      END


      SUBROUTINE splot
*//////////////////////////////////////////////////////////////////////
*//                                                                  //
*//  Very Simple demonstration program with loop over events         //
*//                                                                  //
*//////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER         ninp,nout
      COMMON /c_demo/ ninp,nout
      SAVE
      INTEGER    imax
      PARAMETER (imax = 10000)               ! length ox xpar
      DOUBLE PRECISION       xpar(imax)
      INTEGER      igroup, loop, iev, nevt, ngroup
      DOUBLE PRECISION       xSecPb, xErrPb
      CHARACTER*4  semaph,chdum
      SAVE
*----------------------------------------------------------------------
      WRITE(nout,*) '   '
      WRITE(nout,*) '=============================================='
      WRITE(nout,*) '=============================================='
      WRITE(nout,*) '============ Demo for KK MC =================='
      WRITE(nout,*) '=============================================='
      WRITE(nout,*) '=============================================='
      WRITE(nout,*) '   '
* Read data for main program
      
      OPEN(13,FILE='./pro.input')
      READ(13,'(/,i10)') nevt
      CLOSE(13)
      WRITE(   6,*)   nevt,' requested events '
      WRITE(nout,*)   nevt,' requested events '
*
      CALL KK2f_ReaDataX('../../.KK2f_defaults', 1,imax,xpar)  ! reading general defaults
      CALL KK2f_ReaDataX(         './pro.input', 0,imax,xpar)  ! reading user input
*
C (M.B.) KK2f_Initialize -> KK2f_Initialize_R
      CALL KK2f_Initialize(xpar)                  ! initialize generator
C end (M.B.)
      CALL robol1(-1,xpar)


*-------------------------------------------------------!
*                 main MC loop                          !
*-------------------------------------------------------!
      ngroup = 5000
      iev=0
      DO loop=1,10000000
        DO igroup =1,ngroup
          iev=iev+1
          IF(MOD(iev, ngroup)   .EQ.   1) WRITE( 6,*)  'iev= ',iev
C (M.B.) KK2f_Make -> KK2f_Make_R
          CALL KK2f_Make                        ! make single event
C end (M.B.)
          IF(iev .LE. 10) THEN
             CALL pygive('MSTU(11)=16')
             CALL pylist(2)
             CALL pygive('MSTU(11)=6')
             CALL pylist(2)
          ENDIF
          CALL robol1( 0,xpar)  ! histograming
          IF(iev  .EQ.  nevt)  GOTO 300
        ENDDO
        CALL givsem(semaph)     ! check on semaphore flag
        IF(semaph  .EQ.  'STOP') GOTO 300
        CALL  dumpeh(iev)  !dump partial results on the disk after every ngroup
      ENDDO
 300  CONTINUE
*-------------------------------------------------------!
*                 End  MC loop                          !
*-------------------------------------------------------!
      WRITE(6,*) ' generation finished '
C (M.B.) KK2f_Finalize -> KK2f_Finalize_R
      CALL KK2f_Finalize                        ! final bookkeping, printouts etc.
C end (M.B.)
      CALL KK2f_GetXsecMC(xSecPb, xErrPb)         ! get MC x-section

      CALL robol1(1,xpar)
      CALL  dumpeh(iev)  ! dump final results on the disk

      END

      SUBROUTINE givsem(semaph)
*//////////////////////////////////////////////////////////////////////
*//                                                                  //
*//               READ semaphore flag                                //
*//                                                                  //
*//////////////////////////////////////////////////////////////////////
      IMPLICIT    NONE
      CHARACTER*4 semaph
      INTEGER     ninp2
* ------------------------------------------------------
      ninp2=2
      OPEN(ninp2,FILE='./semaphore')
      READ(ninp2,'(a4)') semaph
      CLOSE(ninp2)
      END

      SUBROUTINE dumpeh(nev)
*//////////////////////////////////////////////////////////////////////
*//                                                                  //
*//           WRITE histos on the disk                               //
*//                                                                  //
*//////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER  nev
      INTEGER  ijklin,ntotin,ntot2n,icy,nouth, ninp2
* ------------------------------------------------------
      nouth=11
      OPEN(nouth,FILE='./pro.hst')
      CALL GLK_hrfile(nouth,' ','n')   !transfer FILE number
      CALL GLK_hrout( 0,icy,' ')       !WRITE on the disk
      CALL GLK_hrend(' ')              !CLOSE FILE
* ------------------------------------------------------
* overWRITE semaphore FILE flag
      ninp2=2
      OPEN(ninp2,FILE='./semaphore')
      WRITE(ninp2,'(a4)') 'CONT'
* append semaphore FILE with new random number seed in
      CALL PseuMar_Out(ijklin,ntotin,ntot2n)
      WRITE(ninp2,'(i10,a)') ijklin, ' = ijklin '
      WRITE(ninp2,'(i10,a)') ntotin, ' = ntotin '
      WRITE(ninp2,'(i10,a)') ntot2n, ' = ntot2n '
      WRITE(ninp2,'(i10,a)') nev,    ' =    nev '
      CLOSE(ninp2)
* ------------------------------------------------------
      END

      SUBROUTINE Robol1(mode,xpar)
*//////////////////////////////////////////////////////////////////////
*//                                                                  //
*//  histogramming is done here                                      //
*//  this subprogram is stil incomplete                              //
*//                                                                  //
*//////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION  xpar(*)
      INTEGER mode
      DOUBLE PRECISION    p1(4),p2(4),p3(4),p4(4)
      INTEGER             k, KFfin
      DOUBLE PRECISION    ss2, vv, ss
C (M.B)
      INTEGER ikf
      DOUBLE PRECISION wt, wtcrud, rs, ene
C end (M.B.)
*=======================================================================
      IF(mode .EQ. -1 ) THEN
C (M.B.) book histograms
c ... general
         CALL KK2f_GetCMSene(ene)
         CALL GLK_Book1(1003,'sprime distribution  $',100,0d0,ene) ! sprime, all flavours
c ... full spectrum per flavor
         CALL GLK_Book1(1005,'sprime gen wtmain u,d$',100,0d0,ene) ! sprime, u+d
         CALL GLK_Book1(1006,'sprime gen wtmain s  $',100,0d0,ene) ! sprime, s
         CALL GLK_Book1(1007,'sprime gen wtmain c  $',100,0d0,ene) ! sprime, c
         CALL GLK_Book1(1008,'sprime gen wtmain b  $',100,0d0,ene) ! sprime, b
c ... zoom resonances
         CALL GLK_Book1(1011,'ZOOM rho, phi$',100,0.2d0,1.2d0)     ! u,d,s
         CALL GLK_Book1(1012,'ZOOM J/Psi   $', 50,3.0964d0,3.0974d0) ! c
         CALL GLK_Book1(1013,'ZOOM Upsilon $', 40,9.4602d0,9.4606d0) ! b
C end (M.B.)
*=======================================================================
      ELSEIF(mode .EQ. 0 ) THEN
* some simple histograming using fermion momenta only
         CALL HepEvt_GetBeams(p1,p2)
         CALL HepEvt_GetFfins(p3,p4)
* Flavour distribution (KF code)
         CALL HepEvt_GetKFfin(KFfin)
         ss  =  (p1(4)+p2(4))**2 -(p1(3)+p2(3))**2 -(p1(2)+p2(2))**2 -(p1(1)+p2(1))**2
         ss2 =  (p3(4)+p4(4))**2 -(p3(3)+p4(3))**2 -(p3(2)+p4(2))**2 -(p3(1)+p4(1))**2
         vv  = 1d0 -ss2/ss
C (M.B.) fill histos
         CALL KK2f_GetWt(wt,wtcrud)
         rs = dsqrt(1d0-vv)*ene
         CALL GLK_Fil1(1003,rs,wt)
         CALL MBrA_GetKF(ikf)
c ... full spectrum per flavor
         IF(ikf.LE.2) CALL GLK_Fil1(1005,rs,wt)
         IF(ikf.EQ.3) CALL GLK_Fil1(1006,rs,wt)
         IF(ikf.EQ.4) CALL GLK_Fil1(1007,rs,wt)
         IF(ikf.EQ.5) CALL GLK_Fil1(1008,rs,wt)
c ... zoom resonances
         IF(ikf.LE.3) CALL GLK_Fil1(1011,rs,wt)
         IF(ikf.EQ.4) CALL GLK_Fil1(1012,rs,wt)
         IF(ikf.EQ.5) CALL GLK_Fil1(1013,rs,wt)
C end (M.B.)
*=======================================================================
      ELSE
*     finalization, printouts
C (M.B.) and plot
         CALL GLK_Print(1003)
         CALL GLK_Print(1005)
         CALL GLK_Print(1006)
         CALL GLK_Print(1007)
         CALL GLK_Print(1008)
         CALL GLK_Print(1011)
         CALL GLK_Print(1012)
         CALL GLK_Print(1013)
C end (M.B.)
      ENDIF
      END
