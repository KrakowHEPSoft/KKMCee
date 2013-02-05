*///////////////////////////////////////////////////////////////////////////////
*//
*//   make Inclusive-start
*//   make Tau-start
*//   make Mu-start
*//   make Down-start
*//   make Beast-start          <-- Beamsstrahlung
*//
*///////////////////////////////////////////////////////////////////////////////
      PROGRAM Main
*///////////////////////////////////////////////////////////////////////////////
*//                   main program for MC production                          //
*///////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
*
      INTEGER            ninp, nout
      COMMON /c_MainPro/ ninp, nout
      CHARACTER*4 semaph
      INTEGER     ninp3, ntot2n, ninph, ninp2, ntotin, ijklin
      SAVE
*--------------------------------------------------------------------------------
      ninp =5   ! standard input
      nout =16  ! general output for everybody including glibk
      OPEN(nout,FILE='./pro.output')
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
      CALL ProdMC
* *******************
      END


      SUBROUTINE ProdMC
*///////////////////////////////////////////////////////////////////////////////
*//                   MC production                                           //
*///////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
*
      INTEGER            ninp, nout
      COMMON /c_MainPro/ ninp, nout
      INTEGER    imax
*
      PARAMETER (imax = 10000)               ! length ox xpar
      DOUBLE PRECISION        xpar(imax)
*
      CHARACTER*4   semaph,chdum
      INTEGER       kat1,kat2,kat3,kat4
      INTEGER       igroup, ngroup, nevt, loop, iev
      DOUBLE PRECISION        xSecPb, xErrPb
      SAVE
*-------------------------------------------------------------------------------
* Read data for main program
      OPEN( ninp,FILE='./pro.input')
      READ( ninp,'(4i2)') kat1,kat2,kat3,kat4
      WRITE(nout,'(4a6/4i6)')
     $     'kat1','kat2','kat3','kat4',
     $      kat1 , kat2 , kat3 , kat4 
      READ(ninp,'(i10)') nevt
      CLOSE(ninp)
      WRITE(   6,*)   nevt,' requested events '
      WRITE(nout,*)   nevt,' requested events '
*
      CALL KK2f_ReaDataX('../../.KK2f_defaults', 1,imax,xpar)  ! reading general defaults
      CALL KK2f_ReaDataX(         './pro.input', 0,imax,xpar)  ! reading user input
*
      CALL KK2f_Initialize(xpar)                  ! initialize generator
*
      IF(kat1 .EQ. 1) CALL Robol1(-1,xpar)
*-------------------------------------------------------!
*                 main MC loop                          !
*-------------------------------------------------------!
      ngroup = 500
      iev=0
      DO loop=1,10000000
        DO igroup =1,ngroup
          iev=iev+1
          IF(MOD(iev, ngroup) .EQ. 1) WRITE( 6,*)  'iev= ',iev
          CALL KK2f_Make                          ! make single event
*   Control printouts
*          CALL momprt(' YFSPRO ', 6,iev,1,10,pf1,pf2,qf1,qf2,nphot,sphot,KFfin)
*          CALL dumpri('*momini*', 6,iev,1,10,xf1,xf2,nphox,xphot)
          IF(iev .LE. 10) THEN
             CALL PYgive('MSTU(11)=16')
             CALL PYlist(1)
             CALL PYgive('MSTU(11)=6')
             CALL PYlist(1)
          ENDIF

          IF(kat1 .EQ. 1) CALL Robol1( 0,xpar)    ! histograming
          IF(iev  .EQ.  nevt)     GOTO 300
        ENDDO
        CALL givsem(semaph)  ! check on semaphore flag
        IF(semaph  .EQ.  'STOP') GOTO 300
        CALL  dumpeh(iev)  ! dump partial results on the disk
      ENDDO
 300  CONTINUE
*-------------------------------------------------------!
*            End  main MC loop                          !
*-------------------------------------------------------!
      WRITE(6,*) ' generation finished '
      CALL KK2f_Finalize                          ! final bookkeping, printouts etc.
      CALL KK2f_GetXsecMC(xSecPb, xErrPb)         ! get MC x-section
*
      IF(kat1 .EQ. 1) CALL Robol1(1,xpar)
      CALL  dumpeh(iev)    ! dump final results on the disk
      END

      SUBROUTINE givsem(semaph)
*///////////////////////////////////////////////////////////////////////////////
*//             READ semaphore flag                                           //
*///////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      CHARACTER*4 semaph
      INTEGER     ninp2
* ------------------------------------------------------
      ninp2=2
      OPEN(ninp2,FILE='./semaphore')
      READ(ninp2,'(a4)') semaph
      CLOSE(ninp2)
      END

      SUBROUTINE dumpeh(nev)
*///////////////////////////////////////////////////////////////////////////////
*//             WRITE histos on the disk                                      //
*///////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER    ijklin, ntotin, ntot2n, ninp2, nev, nouth, icy
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
*////////////////////////////////////////////////////////////////////////////////////////
*//                                                                                    //
*//   Example of histogramming, with acces to events through getters                   //
*//                                                                                    //
*//                                                                                    //
*////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION    xpar(*)
      INTEGER   mode
*
      DOUBLE PRECISION    p1(4),p2(4),p3(4),p4(4)
      INTEGER             NphAll,NphIni,NphFin,NphBst
      DOUBLE PRECISION    PhoAll(100,4),PhoIni(100,4),PhoFin(100,4),PhoBst(100,4)
      DOUBLE PRECISION    WtSet(1000), WtMain, WtCrud
      DOUBLE PRECISION    vv, vmin, vmax, ss, ss2, CMSene
      DOUBLE PRECISION    EneTot, EneIni, EneFin, EneBst
      DOUBLE PRECISION    xlTot, xlIni, xlFin, xlBst
      DOUBLE PRECISION    xlmin, xlmax, cosThe, wtNOint
      INTEGER             i, nbv, idv, KFfin, nbl
      SAVE
*----------------------------------------------------------------------------------------
      IF(mode .EQ. -1 ) THEN
*     book histograms
         nbv  = 400
         vmin = 0d0
         vmax = 1d0
         idv  = 50000
         CALL GLK_Book1(idv+1,'d_sigma/dv(v)      $',nbv,vmin,vmax)
         nbl   =  60
         xlmin = -6d0
         xlmax =  1d0
         CALL GLK_Book1(idv +2,'ln10(Ene/Ebeam) total ALL      $',nbl,xlmin,xlmax)
         CALL GLK_Book1(idv +3,'ln10(Ene/Ebeam) total ISR      $',nbl,xlmin,xlmax)
         CALL GLK_Book1(idv +4,'ln10(Ene/Ebeam) total FSR      $',nbl,xlmin,xlmax)
         CALL GLK_Book1(idv +5,'ln10(Ene/Ebeam) beamstrahlung  $',nbl,xlmin,xlmax)
         CALL GLK_Book1(idv +7,'NphAll             $',16,0d0,16d0)
         CALL GLK_Book1(idv +9,'KFcode             $',16,1d0,17d0)
         CALL GLK_Book1(idv+20,'2f mass [GeV]      $',30,88d0,94d0)
         CALL GLK_Book1(idv+51,'dsig/dcos(theta), wt=1, IFI on       $',8, -1d0,1d0)
         CALL GLK_Book1(idv+52,'dsig/dcos(theta), wted, IFI off      $',8, -1d0,1d0)
         CALL GLK_Book1(idv+53,'dsig/dcos(theta), ratio IFIon/IFIoff $',8, -1d0,1d0)
         CALL GLK_Book1(idv+59,'IFIoff/IFIon weight $',50, 0d0, 10d0)
*=======================================================================
      ELSEIF(mode .EQ. 0 ) THEN
* Histograming
         CALL HepEvt_GetBeams(p1,p2)
         CALL HepEvt_GetFfins(p3,p4)
         CALL HepEvt_GetPhotAll(NphAll,PhoAll)
         CALL HepEvt_GetPhotIni(NphIni,PhoIni)
         CALL HepEvt_GetPhotFin(NphFin,PhoFin)
         CALL HepEvt_GetPhotBst(NphBst,PhoBst)

         CALL KK2f_GetWtList(WtMain,WtSet)
         CALL KK2f_GetWtCrud(WtCrud)
         
         ss =  (p1(4)+p2(4))**2 -(p1(3)+p2(3))**2
     $        -(p1(2)+p2(2))**2 -(p1(1)+p2(1))**2
         ss2=  (p3(4)+p4(4))**2 -(p3(3)+p4(3))**2
     $        -(p3(2)+p4(2))**2 -(p3(1)+p4(1))**2
         vv=1d0-ss2/ss
         CMSene = p1(4)+p2(4)
         cosThe = p3(3)/DSQRT(p3(1)**2 +p3(2)**2+p3(3)**2)
* Energy loss
         CALL GLK_Fil1(idv+1, vv, WtMain)
* Sum of energies of ALL photons
         EneTot = 0d0
         DO i=1,NphAll
            EneTot=EneTot+PhoAll(i,4)
         ENDDO
         CALL GLK_Fil1(idv+2, EneTot, WtMain)
         EneTot = 2d0*EneTot/CMSene
         xlTot = -1000d0
         IF( EneTot.GT.0d0) xlTot = LOG10(EneTot)
* Z radiative return mass
         CALL GLK_Fil1(idv+20, DSQRT(ss2), WtMain)
* Sum of energies of ISR photons
         EneIni = 0d0
         DO i=1,NphIni
            EneIni=EneIni+PhoIni(i,4)
         ENDDO
         EneIni = 2d0*EneIni/CMSene
         xlIni = -1000d0
         IF( EneIni.GT.0d0) xlIni = LOG10(EneIni)
         CALL GLK_Fil1(idv+3, xlIni, WtMain)
* Sum of energies of FSR photons
         EneFin = 0d0
         DO i=1,NphFin
            EneFin=EneFin+PhoFin(i,4)
         ENDDO
         EneFin = 2d0*EneFin/CMSene
         xlFin = -1000d0
         IF( EneFin.GT.0d0) xlFin = LOG10(EneFin)
         CALL GLK_Fil1(idv+4, xlFin, WtMain)
* Sum of energies of beamsstrahlung photons
         EneBst = 0d0
         DO i=1,NphBst
            EneBst=EneBst+PhoBst(i,4)
         ENDDO
         EneBst = 2d0*EneBst/CMSene
         xlBst = -1000d0
         IF( EneBst.GT.0d0) xlBst = LOG10(EneBst)
         CALL GLK_Fil1(idv+5, xlBst, WtMain)
* Photon multiplicity (total) including beamsstrahlung
         CALL GLK_Fil1(idv+7, NphAll*1.00001d0, WtMain)
* Flavour distribution (KF code)
         CALL HepEvt_GetKFfin(KFfin)
         CALL GLK_Fil1(idv+9, KFfin*1.00001d0, WtMain)
* Angular distribution with and without IFI
         wtNOint = WtSet(253) ! WtBest=WtSet(203)
         CALL GLK_Fil1(idv+59, wtNOint, 1d0)
         IF( vv .LT. 0.1d0 ) THEN
            CALL GLK_Fil1(idv+51, cosThe, WtMain)
            CALL GLK_Fil1(idv+52, cosThe, wtNOint)
         ENDIF
*=======================================================================
      ELSE
*     finalization, printouts
         CALL GLK_Print(idv+1)
         CALL GLK_Print(idv+2)
         CALL GLK_Print(idv+3)
         CALL GLK_Print(idv+4)
         CALL GLK_Print(idv+5)
         CALL GLK_Print(idv+7)
         CALL GLK_Print(idv+9)
         CALL GLK_Print(idv+20)
         CALL GLK_idopt(idv+51,'ERRO')
         CALL GLK_idopt(idv+52,'ERRO')
         CALL GLK_idopt(idv+53,'ERRO')
         CALL GLK_Print(idv+51)
         CALL GLK_Print(idv+52)
         CALL GLK_Operat(idv+51, '/',idv+52, idv+53, 1d0, 1d0)
         CALL GLK_Print(idv+53)
         CALL GLK_Print(idv+59)
      ENDIF
      END

