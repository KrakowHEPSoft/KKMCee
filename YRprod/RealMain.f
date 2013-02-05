*///////////////////////////////////////////////////////////////////////////////
*//
*//   make Inclusive-start
*//   make Tau-start
*//   make Mu-start
*//   make DOwn-start
*//   make Beast-start          <-- Beamsstrahlung
*//
*///////////////////////////////////////////////////////////////////////////////
      PROGRAM Main
*///////////////////////////////////////////////////////////////////////////////
*//                   main PROGRAM for MC production                          //
*/////////////////ifphot=ifphot .AND. ifpart//////////////////////////////////////////////////////////////
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
      IF(semaph   .EQ.   'STAR') THEN
         WRITE(6,*) ' ------- starting from the scratch ----------'
* READ initial (root) random number seed
         ninp3=3
         OPEN(ninp3,FILE='./iniseed')
         READ(ninp3,'(i10)') ijklin
         READ(ninp3,'(i10)') ntotin
         READ(ninp3,'(i10)') ntot2n
         CALL PseuMar_Initialize(ijklin,ntotin,ntot2n)
      ELSEIF(semaph   .EQ.   'CONT') THEN
         WRITE(6,*) ' ------- restoring from the disk   ----------'
* restore histograms from the disk
         ninph=10
         OPEN(ninph,FILE='./'//'pro.hst')
         CALL GLK_hrFILE(ninph,' ',' ')      !transfer FILE number
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
      ELSEIF(semaph   .EQ.   'STOP') THEN
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
      INTEGER       kat1,kat2,kat3,kat4,k
      INTEGER       igroup, ngroup, nevt, loop, iev
      DOUBLE PRECISION        xSecPb, xErrPb,dum,dum2,duma,duma2,cmsene
      DOUBLE PRECISION p1(4),p2(4)
      SAVE
*-------------------------------------------------------------------------------
* READ DATA for main PROGRAM
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
      CALL KK2f_READATAX('../../.KK2f_defaults', 1,imax,xpar)  ! READing general defaults
      CALL KK2f_READATAX(         './pro.input', 0,imax,xpar)  ! READing user input
*
      CALL KK2f_Initialize(xpar)                  ! initialize generator
*
      IF(kat1  .EQ.  1) CALL Robol1(-1,xpar)
      k=xpar(10)
      CALL WSHOP2000(-1,k)
*-------------------------------------------------------!
*                 main MC loop                          !
*-------------------------------------------------------!
      ngroup = 5000
      ngroup = 500
      iev=0
      DO loop=1,10000000
        DO igroup =1,ngroup
          iev=iev+1
          IF(MOD(iev, ngroup)  .EQ.  1) WRITE( 6,*)  'iev= ',iev
          CALL KK2f_Make                          ! make single event
*   Control printouts
*          CALL momprt(' YFSPRO ', 6,iev,1,10,pf1,pf2,qf1,qf2,nphot,sphot,KFfin)
*          CALL dumpri('*momini*', 6,iev,1,10,xf1,xf2,nphox,xphot)
          IF(iev  .LE.  10) THEN
             CALL pygive('MSTU(11)=16')
             CALL pylist(1)
             CALL pygive('MSTU(11)=6')
             CALL pylist(1)
          ENDIF

          IF(kat1  .EQ.  1) CALL Robol1( 0,xpar)    ! histograming
          CALL WSHOP2000(0,k)
          IF(iev   .EQ.   nevt)     GOTO 300
        ENDDO
        CALL givsem(semaph)  ! check on semaphore flag
        IF(semaph   .EQ.   'STOP') GOTO 300
        CALL  dumpeh(iev)  ! dump partial results on the disk
      ENDDO
 300  CONTINUE
*-------------------------------------------------------!
*            END  main MC loop                          !
*-------------------------------------------------------!
      WRITE(6,*) ' generation finished '
      CALL KK2f_Finalize                          ! final bookkeping, printouts etc.
      CALL KK2f_GetXsecMC(xSecPb, xErrPb)         ! get MC x-section
*
      IF(kat1  .EQ.  1) CALL Robol1(1,xpar)
      CALL WSHOP2000(1,k)
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
      CALL GLK_hrFILE(nouth,' ','n')   !transfer FILE number
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
      DOUBLE PRECISION    xlmin, xlmax
      INTEGER             i, nbv, idv, KFfin, nbl
      SAVE
*----------------------------------------------------------------------------------------
      IF(mode  .EQ.  -1 ) THEN
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
*=======================================================================
      ELSEIF(mode  .EQ.  0 ) THEN
* Histograming
         CALL HepEvt_GetBeams(p1,p2)
         CALL HepEvt_GetFfins(p3,p4)
         CALL HepEvt_GetPhotAll(NphAll,PhoAll)
         CALL HepEvt_GetPhotIni(NphIni,PhoIni)
         CALL HepEvt_GetPhotFin(NphFin,PhoFin)
         CALL HepEvt_GetPhotBst(NphBst,PhoBst)
         CALL KK2f_GetWtAll(WtMain,WtCrud,WtSet)
         
         ss =  (p1(4)+p2(4))**2 -(p1(3)+p2(3))**2
     $        -(p1(2)+p2(2))**2 -(p1(1)+p2(1))**2
         ss2=  (p3(4)+p4(4))**2 -(p3(3)+p4(3))**2
     $        -(p3(2)+p4(2))**2 -(p3(1)+p4(1))**2
         vv=1d0-ss2/ss
         CALL GLK_Fil1(idv+1, vv, WtMain)

         CMSene = p1(4)+p2(4)
* Sum of energies of ALL photons
         EneTot = 0d0
         DO i=1,NphAll
            EneTot=EneTot+PhoAll(i,4)
         ENDDO
         CALL GLK_Fil1(idv+2, EneTot, WtMain)
         EneTot = 2d0*EneTot/CMSene
         xlTot = -1000d0
         IF( EneTot .GT. 0d0) xlTot = LOG10(EneTot)

* Z radiative RETURN mass
         CALL GLK_Fil1(idv+20, DSQRT(ss2), WtMain)
* Sum of energies of ISR photons
         EneIni = 0d0
         DO i=1,NphIni
            EneIni=EneIni+PhoIni(i,4)
         ENDDO
         EneIni = 2d0*EneIni/CMSene
         xlIni = -1000d0
         IF( EneIni .GT. 0d0) xlIni = LOG10(EneIni)
         CALL GLK_Fil1(idv+3, xlIni, WtMain)

* Sum of energies of FSR photons
         EneFin = 0d0
         DO i=1,NphFin
            EneFin=EneFin+PhoFin(i,4)
         ENDDO
         EneFin = 2d0*EneFin/CMSene
         xlFin = -1000d0
         IF( EneFin .GT. 0d0) xlFin = LOG10(EneFin)
         CALL GLK_Fil1(idv+4, xlFin, WtMain)

* Sum of energies of beamsstrahlung photons
         EneBst = 0d0
         DO i=1,NphBst
            EneBst=EneBst+PhoBst(i,4)
         ENDDO
         EneBst = 2d0*EneBst/CMSene
         xlBst = -1000d0
         IF( EneBst .GT. 0d0) xlBst = LOG10(EneBst)
         CALL GLK_Fil1(idv+5, xlBst, WtMain)

* Photon multiplicity (total) including beamsstrahlung
         CALL GLK_Fil1(idv+7, NphAll*1.00001d0, WtMain)

* Flavour distribution (KF code)
         CALL HepEvt_GetKFfin(KFfin)
         CALL GLK_Fil1(idv+9, KFfin*1.00001d0, WtMain)
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
      ENDIF
      END

      SUBROUTINE WSHOP2000(MODE,ikey)
C ikey is the keywgt to know how normalize. however FUNCTION 
C calculating weight must be adjusted by hand
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION P1(4),P2(4)
      COMMON /adresnik/ IDV

      IF     (MODE .EQ. -1) THEN
        IDV=3219
        CALL memorizer4(-1,ikey,dum,dum2,duma,duma2)
        IDV=3319
        CALL memorizer4(-1,ikey,dum,dum2,duma,duma2)
        IDV=3419
        CALL memorizer4(-1,ikey,dum,dum2,duma,duma2)
        IDV=3519
        CALL memorizer4(-1,ikey,dum,dum2,duma,duma2)
      ELSEIF (MODE .EQ.  0) THEN
        IDV=3219
        CALL buker4(IDV)
        IDV=3319
        CALL buker4(IDV)
        IDV=3419
        CALL buker4(IDV)
        IDV=3519
        CALL buker4(IDV)
      ELSE
        CALL KK2f_GetXsecMC(xSecPb, xErrPb)         ! get MC x-section
        CALL HepEvt_GetBeams(p1,p2)
        CMSene = p1(4)+p2(4)
        IDV=3219
        CALL memorizer4(1,ikey,xSecPb, xErrPb,cmsene,duma2)
        IDV=3319
        CALL memorizer4(1,ikey,xSecPb, xErrPb,cmsene,duma2)
        IDV=3419
        CALL memorizer4(1,ikey,xSecPb, xErrPb,cmsene,duma2)
        IDV=3519
        CALL memorizer4(1,ikey,xSecPb, xErrPb,cmsene,duma2)
      ENDIF
      END

      SUBROUTINE memorizer4(mode,id,x,x2,y,y2)
C this routine is prepared to store and print info for observables.
C at the moment inFORMATion is stored in matrices, but it is 
C straightforward to use GLK also.
C mode denotes: -1 initialization all other PARAMETERs are THEN dummy
C                  id transmit in keywgt though
C                0 collecting active input are mode,id,x,y
C                  id -observable index, x-weight, y-weight*sign for asymetry
C                1 printout  active input are mode,id,x,x2,y which denote:
C                  ID-no of generated events, x,x2-sig and its stat err
C                  y-center of mass energy
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON / INOUT / INUT,IOUT
*/////////////////////////////////////////////////////////////////////////////////////////////
*//                                                                                         //
*//                                                                                         //
*//                               Pseudo-CLASS  tabutil                                     //
*//                                                                                         //
*//             DATA on table entries and options on tab makin.                             //
*//                                                                                         //
*//                                                                                         //
*/////////////////////////////////////////////////////////////////////////////////////////////

      CHARACTER*13 label
      CHARACTER*14 progname
      CHARACTER*1  sora
      REAL         obs,errstat,errsyst
      CHARACTER*9 unitobs,uniterr
      CHARACTER*26 comment
 
      PARAMETER(NMAX=1000)
      CHARACTER*13 m_label(NMAX) 
      CHARACTER*14 m_progname(NMAX)
      CHARACTER*1  m_sora(NMAX)
      REAL         m_cmsene(NMAX),m_obs(NMAX),m_errstat(NMAX),m_errsyst(NMAX)
      CHARACTER*9 m_unitobs(NMAX),m_uniterr(NMAX)
      CHARACTER*26 m_comment(NMAX) 
      COMMON /c_store/ m_lev_tab,m_lines_max,m_lines(NMAX,2)
     $                ,Nfilmax,m_label,m_progname,m_sora,m_cmsene,m_obs,m_errstat,m_errsyst,m_unitobs,m_uniterr,m_comment
      SAVE /c_store/
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//                      END of CLASS  tabutil                               //
*//////////////////////////////////////////////////////////////////////////////
      REAL*8 wto(30),wto2(30),wta(30),wta2(30) 
      CHARACTER*1 sorb
      CHARACTER*9 unitobsa,uniterra
      COMMON /adresnik/ IDV
      SAVE

      IF (MODE .EQ. -1) THEN
        keywgt=id
        WRITE(*,*) 'keywgt=',keywgt
        CALL GLK_Book1(idv+1,'wto      $',30,0.5D0,30.5D0)
        CALL GLK_Book1(idv+2,'wto2     $',30,0.5D0,30.5D0)
        CALL GLK_Book1(idv+3,'wta      $',30,0.5D0,30.5D0)
        CALL GLK_Book1(idv+4,'wta2     $',30,0.5D0,30.5D0)
        DO k=1,29
          wto(k)=0d0
          wto2(k)=0d0
          wta(k)=0
          wta2(k)=0
        ENDDO
      ELSEIF (MODE .EQ. 0) THEN
          wto(ID) =wto (ID)+x
          wto2(ID)=wto2(ID)+x**2
          wta(ID) =wta (ID)+y
          wta2(ID)=wta2(ID)+y**2
          CALL GLK_Fil1(idv+1,1D0*ID,x )
          CALL GLK_Fil1(idv+2,1D0*ID,x**2 )
          CALL GLK_Fil1(idv+3,1D0*ID,y )
          CALL GLK_Fil1(idv+4,1D0*ID,y**2 )
      ELSEIF (MODE .EQ. 1) THEN
C taking info on the run from main PROGRAM 
        sig=x
        err=x2
        cmsene=y

        CALL GLK_UnPak(idv+1,wto,'enry',idum)
        CALL GLK_UnPak(idv+2,wto2,'enry',idum)
        CALL GLK_UnPak(idv+3,wta,'enry',idum)
        CALL GLK_UnPak(idv+4,wta2,'enry',idum)
        nevtes=wto(30)
        IF (keywgt .EQ. 1) THEN
          LevelPrint=0
          CALL KarLud_Finalize(LevelPrint,sig,err)
          sig0pb =  BornV_Sig0nb(CMSene)*1000
          sig=sig*sig0pb
          err=err*sig0pb
         ENDIF  
        WRITE(iout,*) 'nevtes:',nevtes
        WRITE(iout,*) 'sig and err ',sig,err
C setting text variables etc for printouts.
                         progname='KKMC 4.13  '
        IF (IDV .EQ. 3319) progname='KKMC4.13 noIFI'
        IF (IDV .EQ. 3419) progname='KKMC4.13 EEX2 '
        IF (IDV .EQ. 3519) progname='KKMC4.13 EEX3 '
        unitobs='pb'
        uniterr='pb'
        unitobsa='  '
        uniterra='  '
        comment='preliminary'
        comment='25 05 tests'
        sora='S'
        sorb='A'
        m_label(1)='ALEPH-12'
        m_label(2)='ALEPH-15'
        m_label(3)='DELPHI9'
        m_label(4)='DELPHI12'
        m_label(5)='LT15'
        m_label(6)='LT18'
        m_label(7)='Opal14'
        m_label(8)='Opal15'
        m_label(9)='Opal20'
        m_label(10)='Opal21'
        m_errsyst(1)=0.02
        m_errsyst(2)=0.05
        m_errsyst(3)=0.02
        m_errsyst(4)=0.05
        m_errsyst(5)=0.02
        m_errsyst(6)=0.05
        m_errsyst(7)=0.02
        m_errsyst(8)=0.02
        m_errsyst(9)=0.05
        m_errsyst(10)=0.05

        m_label(11)='Aleph5'
        m_label(12)='Aleph6'
        m_label(13)='Delphi4'
        m_label(14)='Delphi5'
        m_label(15)='LT9'
        m_label(16)='LT10'
        m_label(17)='Opal6'
        m_label(18)='Opal7'
        m_errsyst(11)=0.02
        m_errsyst(12)=0.02
        m_errsyst(13)=0.02
        m_errsyst(14)=0.02
        m_errsyst(15)=0.02
        m_errsyst(16)=0.02
        m_errsyst(17)=0.02
        m_errsyst(18)=0.02

        m_label(19)='IAleph5'
        m_label(20)='IAleph6'
        m_label(21)='IDelphi5'
        m_label(22)='IDelphi6'
        m_label(23)='ILT9'
        m_label(24)='ILT10'
        m_label(25)='ILT11'
        m_label(26)='IOpal6'
        m_label(27)='IOpal7'
        m_label(28)='IOpal8'
        m_label(29)='IOpal9'
        m_errsyst(19)=0.02
        m_errsyst(20)=0.02
        m_errsyst(21)=0.02
        m_errsyst(22)=0.02
        m_errsyst(23)=0.02
        m_errsyst(24)=0.02
        m_errsyst(25)=0.02
        m_errsyst(26)=0.02
        m_errsyst(27)=0.02
        m_errsyst(28)=0.02
        m_errsyst(29)=0.02

C calculating and printing entries for table -------
        DO k=1,29
          wts=wto(k)/nevtes
          erwts=SQRT((wto2(k)/nevtes-wts**2)/nevtes)
          obs=sig*wts
          errstat=SQRT(wts**2*err**2+sig**2*erwts**2)
          errsyst=obs*m_errsyst(k)
          WRITE(  16,10) m_label(k),progname,sora,cmsene,obs,errstat,errsyst,unitobs,uniterr,comment
          IF(wta(k) .GT. 0d0) THEN
            ass=wta(k)/wto(k)
            errstat=1d0/SQRT(wto(k))*SQRT(1d0-ass**2) !surely not complete
            errsyst=m_errsyst(k)
            WRITE(  16,10) m_label(k),progname,sorb,cmsene,ass,errstat,errsyst,unitobsa,uniterra,comment
          ENDIF
        ENDDO
 10   FORMAT(1H=,a13,1H=,a14,1H=,a1,1H=,f7.3,1H=,e14.8,1H=,e14.8,1H=,e14.8,1H=,a9,1H=,a9,1H=,a26,1H=)
      ELSE
         WRITE(*,*) 'memorizer4 wrong mode=',mode
         STOP
      ENDIF
      END
      
      FUNCTION FHISTWT(ID)
C this FUNCTION can be used for weighted events.
C In buker4 it is CALLed with ID=0, it is prepared
C that actual weight option will be set by memorizer4

      IMPLICIT REAL*8 (A-H,O-Z)
      DOUBLE PRECISION    WtSet(1000), WtMain, WtCrud
      FHISTWT=1D0
      CALL KK2f_GetWtAll(WtMain,WtCrud,WtSet)
      FHISTWT=WtMain
                      FHISTWT=WTCrud*WtSet(203)
      IF (ID .EQ. 3319) FHISTWT=WTCrud*WtSet(253)
      IF (ID .EQ. 3419) FHISTWT=WTCrud*WtSet(73)
      IF (ID .EQ. 3519) FHISTWT=WTCrud*WtSet(74)
      END

      SUBROUTINE buker4(ID)
! observables for muons
! routine calculates whether for observable of the given index KK
! it is in or out the cuts.
! INPUT: HEPEVT COMMON block
! 
! OUTPUT: for every observable numbered KK CALL on storing routine is executed
!         (if in cuts) or is not executed ELSEwhere

      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMXHEP=2000)
      COMMON/d_HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      REAL*8 phep,vhep,pol,sig,err
       SAVE  /d_HEPEVT/
      COMMON / INOUT / INUT,IOUT
      COMMON /kalinout/ wtkal(6)
      REAL*8 wtkal,wtsum(6),wtsum2(6),wto,wto2,wta,wta2,wt
      REAL*8  phsum(4),PH(4),PI,PHTR(4),qp(4),qm(4),xmp(10),xmm(10),ph1(4)
      LOGICAL IFacc,ifphot,ifpart,ifphot1
      CHARACTER*1 sorb
      CHARACTER*9 unitobsa,uniterra
      DATA pi /3.141592653589793238462643D0/
         cmsene=phep(4,1)*2
         wt=FHISTWT(ID)
          DO L=1,4
           QP(l)=PHEP(l,3)
           QM(l)=PHEP(l,4)
          ENDDO

         NEVTES=NEVTES+1
         CALL memorizer4(0,30,1D0,0D0,0D0,0D0)
C===============================================================
C     'ALEPH-12'
      kk=1
        NPH=0
        ITR=0
        IFacc=(xCOS(qp) .LT. 0.95 .AND. xCOS(qm) .LT. 0.95)
        IFacc=ifacc .AND. (enrg(qp) .GT. 5.0 .OR. enrg(qm) .GT. 5.0)
        IFacc=ifacc .AND. (xacol(qp,qm,3) .GT. 30.0)
        IFphot=.false.
        DO K=5,NHEP
          DO L=1,4
           PH(l)=PHEP(l,k)
          ENDDO
          IFpart=(XCOS(PH) .LT. 0.95)
          IFpart=ifpart .AND. (ENRG(PH) .GT. 15D0)
          IFpart=ifpart .AND. (xacol(PH,qm,3) .GT. 10.0) .AND. (xacol(PH,qp,3) .GT. 10.0)

          IFphot=ifphot .OR. ifpart
        ENDDO
          IFacc=ifacc .AND. ifphot .AND. ((qp(4)+QM(4)) .GT. 0.6*cmsene)
        IF (ifacc) THEN
          CALL memorizer4(0,KK,wt,0D0,0D0,0D0)
        ENDIF

C===============================================================
C     'ALEPH-15'
       kk=2
        NPH=0
        ITR=0
        IFacc=(xCOS(qp) .LT. 0.95 .AND. xCOS(qm) .LT. 0.95)
        IFacc=ifacc .AND. (enrg(qp) .GT. 5.0 .OR. enrg(qm) .GT. 5.0)
        IFacc=ifacc .AND. (xacol(qp,qm,3) .GT. 30.0)
        DO K=5,NHEP
          DO L=1,4
           PH(l)=PHEP(l,k)
          ENDDO
          IFpart=(XCOS(PH) .LT. 0.95)
          IFpart=ifpart .AND. (ENRG(PH) .GT. 15D0)
          IFpart=ifpart .AND. (xacol(PH,qm,3) .GT. 20.0) .AND. (xacol(PH,qp,3) .GT. 20.0)
          IF (ifpart) THEN
             nph=nph+1
             CALL pMsum(QP,PH,PHSUM)
             xmp(nph)=XMAS(phsum)
             CALL pMsum(QM,PH,PHSUM)
             xmm(nph)=XMAS(phsum)
          ENDIF
        ENDDO
        IFacc=ifacc .AND. (nph .GE. 2)
        IFpart=.false.
        DO k=1,nph
          DO l=1,nph
            IF (k .NE. l) THEN
              IFpart=ifpart .OR. (ABS(xmp(l)-xmm(k)) .LT. 5D0)
            ENDIF
          ENDDO
        ENDDO

        IF (ifacc .AND. ifpart) THEN
          CALL memorizer4(0,KK,wt,0D0,0D0,0D0)
        ENDIF


C===============================================================
C     'DELPHI9'
      kk=3
        NPH=0
        ITR=0
        esum=0
        ejetmax=0
        IFacc=(xCOS(qp) .LT. 0.9063 .AND. (enrg(qp) .GT. 5.0))
        IF (ifacc)                                   ejetmax=qp(4)

        IFacc=ifacc .OR. (xCOS(qm) .LT. 0.9063 .AND. (enrg(qm) .GT. 5.0))
        IF (xCOS(qm) .LT. 0.9063 .AND. (enrg(qm) .GT. 5.0)) ejetmax=max(ejetmax,qm(4))

        IF (xCOS(qp) .LT. 0.9397) esum=esum+qp(4)
        IF (xCOS(qm) .LT. 0.9397) esum=esum+qm(4)
        IFphot=.false.
        DO K=5,NHEP
          DO L=1,4
           PH(l)=PHEP(l,k)
          ENDDO
          IFpart=(XCOS(PH) .LT. 0.9848)
          IFpart=ifpart .AND. (ENRG(PH) .GT. 5D0)
          IF (ifpart) ejetmax=max(ejetmax,ph(4)) 
          IF (xCOS(PH) .LT. 0.9397) esum=esum+PH(4)
          IF(xacol(ph,qp,3) .LT. 15 .AND. qp(4) .GT. 1d0) IFpart=.false.
          IF(xacol(ph,qm,3) .LT. 15 .AND. qm(4) .GT. 1d0) IFpart=.false.

          erec=0
          DO i=5,NHEP
            DO j=1,4
             ph1(j)=phep(j,i)
            ENDDO
            IF (xacol(ph,ph1,3) .LT. 15 .AND. k .NE. i) erec=erec+ph1(4)
          ENDDO
          IFpart=ifpart .AND. erec .LT. 2

          IFphot=ifphot .OR. ifpart
        ENDDO
          IFacc=ifacc .AND. ifphot .AND. esum .GT. 0.2*cmsene .AND. ejetmax .GT. 10D0
        IF (ifacc) THEN
          CALL memorizer4(0,KK,wt,0D0,0D0,0D0)
        ENDIF


C===============================================================
C     'DELPHI12'
       kk=4
        NPH=0
        ITR=0
        esum=0
        IF (xCOS(qp) .LT. 0.9397) esum=esum+qp(4)
        IF (xCOS(qm) .LT. 0.9397) esum=esum+qm(4)
        ejetmax=0
        IFacc=(xCOS(qp) .LT. 0.9063 .AND. (enrg(qp) .GT. 5.0))
        IFacc=ifacc .AND. (xCOS(qm) .LT. 0.9063 .AND. (enrg(qm) .GT. 5.0))
        IFphot=.false.
        DO K=5,NHEP
          DO L=1,4
           PH(l)=PHEP(l,k)
          ENDDO
          IFpart=(XCOS(PH) .LT. 0.9848)
          IFpart=ifpart .AND. (ENRG(PH) .GT. 5D0)
          IF (ifpart) ejetmax=max(ejetmax,ph(4)) 
          IF (xCOS(PH) .LT. 0.9397) esum=esum+PH(4)
          IF(xacol(ph,qp,3) .LT. 15 .AND. qp(4) .GT. 1d0) IFpart=.false.
          IF(xacol(ph,qm,3) .LT. 15 .AND. qm(4) .GT. 1d0) IFpart=.false.

          erec=0
          DO i=5,NHEP
            DO j=1,4
             ph1(j)=phep(j,i)
            ENDDO
            IF (xacol(ph,ph1,3) .LT. 15 .AND. k .NE. i) erec=erec+ph1(4)
          ENDDO
          IFpart=ifpart .AND. erec .LT. 2

          IF (ifpart) THEN
             nph=nph+1
             CALL pMsum(QP,PH,PHSUM)
             xmp(nph)=XMAS(phsum)
             CALL pMsum(QM,PH,PHSUM)
             xmm(nph)=XMAS(phsum)
          ENDIF
        ENDDO

        IFacc=ifacc .AND. (nph .GE. 2)

        IFpart=.false.
        DO k=1,min(2,nph)
          DO l=1,min(2,nph)
            IF (k .NE. l) THEN
              IFpart=ifpart .OR. (ABS(xmp(l)-xmm(k)) .LT. 10D0)
            ENDIF
          ENDDO
        ENDDO

        IF (ifacc .AND. ifpart .AND. esum .GT. 0.2*cmsene) THEN
          CALL memorizer4(0,KK,wt,0D0,0D0,0D0)
        ENDIF
C===============================================================
C     'LT15'
      kk=5
        NPH=0
        ITR=0
        esum=0
        ejetmax=0
        IFacc=(xCOS(qp) .LT. 0.94)
        IFacc=ifacc .OR. (xCOS(qm) .LT. 0.94)
        IFphot=.false.
        DO K=5,NHEP
          DO L=1,4
           PH(l)=PHEP(l,k)
          ENDDO
          IFpart=(XCOS(PH) .LT. 0.97)
          IFpart=ifpart .AND. (ENRG(PH) .GT. 15D0)
          IF (ifpart) THEN
             nph=nph+1
             CALL pMsum(QP,PH,PHSUM)
             xmp(nph)=XMAS(phsum)
             IF(xCOS(qp) .GT. 0.94) xmp(nph)=0
             CALL pMsum(QM,PH,PHSUM)
             xmm(nph)=XMAS(phsum)
            IF(xCOS(qm) .GT. 0.94) xmm(nph)=0
          ENDIF
          IFphot=ifphot .OR. (ENRG(PH) .GT. 20D0 .AND. XCOS(PH) .LT. 0.75 .AND. max(xmp(nph),xmm(nph)) .GT. 70D0)

!          IFphot=ifphot .AND. ifpart
        ENDDO
          IFacc=ifacc .AND. ifphot .AND. nph .LE. 2
        IF (ifacc) THEN
          CALL memorizer4(0,KK,wt,0D0,0D0,0D0)
        ENDIF
c################################################################
C===============================================================
C     'LT18'
      kk=6
        NPH=0
        ITR=0
        esum=0
        ejetmax=0
        IFacc=(xCOS(qp) .LT. 0.94)
        IFacc=ifacc .AND. (xCOS(qm) .LT. 0.94)
        IFphot=.false.
        IFphot1=.true.
        DO K=5,NHEP
          DO L=1,4
           PH(l)=PHEP(l,k)
          ENDDO
          IFpart=(XCOS(PH) .LT. 0.97)
          IFpart=ifpart .AND. (ENRG(PH) .GT. 15D0)
          IF (ifpart) THEN
             nph=nph+1
             CALL pMsum(QP,PH,PHSUM)
             xmp(nph)=XMAS(phsum)
             CALL pMsum(QM,PH,PHSUM)
             xmm(nph)=XMAS(phsum)
           ENDIF

          IFphot=ifphot .OR. (ENRG(PH) .GT. 15D0 .AND. XCOS(PH) .LT. 0.75)
!          IFphot1=ifphot1 .AND. (ENRG(PH) .GT. 15D0 .AND. XCOS(PH) .LT. 0.94)
        ENDDO
          del1=ABS(xmp(1)-xmm(2))
          del2=ABS(xmm(1)-xmp(2))
          sum1=xmp(1)+xmm(2)
          sum2=xmp(2)+xmm(1)
          IFacc=ifacc .AND. ifphot .AND. ifphot1 .AND. nph .EQ. 2
          IFacc=ifacc .AND. ((del2 .LT. 10d0 .AND. sum2 .GT. 100) .OR. (del1 .LT. 10d0 .AND. sum1 .GT. 100))
        IF (ifacc) THEN
          CALL memorizer4(0,KK,wt,0D0,0D0,0D0)
        ENDIF
C===============================================================
C     'Opal14'
      kk=7
        NPH=0
        ITR=0
        evis=0
        egam=0
        IFacc=(xCOS(qp) .LT. 0.95 .AND. xCOS(qm) .LT. 0.95)
        IFacc=ifacc .AND. (pete(qp) .GT. 1.0 .AND. pete(qm) .GT. 1.0)
        IFacc=ifacc .AND. (xacol(qp,qm,3) .GT. 20.0)
        evis=qp(4)+qm(4)
        IFphot=.false.
        DO K=5,NHEP
          DO L=1,4
           PH(l)=PHEP(l,k)
          ENDDO
          IFpart=(XCOS(PH) .LT. 0.95)
          IFpart=ifpart .AND. (ENRG(PH) .GT. 0.05*cmsene/2)
          IFpart=ifpart .AND. (xacol(PH,qm,3) .GT. 20.0) .AND. (xacol(PH,qp,3) .GT. 20.0)
          IF (ifpart) egam=max(egam,ph(4))
          IFphot=ifphot .OR. ifpart
        ENDDO
          evis=evis+egam
          IFacc=ifacc .AND. ifphot .AND. (evis .GT. 1.6*cmsene/2)
        IF (ifacc) THEN
          CALL memorizer4(0,KK,wt,0D0,0D0,0D0)
        ENDIF
C===============================================================
C     'Opal15'
      kk=8
        NPH=0
        ITR=0
        evis=0
        egam=0
        IFacc=(xCOS(qp) .LT. 0.95 .AND. xCOS(qm) .LT. 0.95)
        IFacc=ifacc .AND. (pete(qp) .GT. 1.0 .AND. pete(qm) .GT. 1.0)
        IFacc=ifacc .AND. (xacol(qp,qm,3) .GT. 20.0)
        evis=qp(4)+qm(4)
        IFphot=.false.
        DO K=5,NHEP
          DO L=1,4
           PH(l)=PHEP(l,k)
          ENDDO
          IFpart=(XCOS(PH) .LT. 0.95)
          IFpart=ifpart .AND. (ENRG(PH) .GT. 0.05*cmsene/2)
          IFpart=ifpart .AND. (xacol(PH,qm,3) .GT. 20.0) .AND. (xacol(PH,qp,3) .GT. 20.0)
          IF (ifpart) egam=max(egam,ph(4))
          IFphot=ifphot .OR. ifpart
        ENDDO
          evis=evis+egam
          CALL pmsum(qp,qm,phsum)
          IFacc=ifacc .AND. (85D0 .GT. xmas(phsum) .OR. xmas(phsum) .GT. 95d0)
          IFacc=ifacc .AND. ifphot .AND. (evis .GT. 1.6*cmsene/2)
        IF (ifacc) THEN
          CALL memorizer4(0,KK,wt,0D0,0D0,0D0)
        ENDIF
C===============================================================
C     'Opal20'
      kk=9
        NPH=0
        ITR=0
        evis=0
        egam=0
        egam1=0
        IFacc=(xCOS(qp) .LT. 0.95 .AND. xCOS(qm) .LT. 0.95)
        IFacc=ifacc .AND. (pete(qp) .GT. 1.0 .AND. pete(qm) .GT. 1.0)
        IFacc=ifacc .AND. (xacol(qp,qm,3) .GT. 20.0)
        evis=qp(4)+qm(4)
        IFphot=.false.
        DO K=5,NHEP
          DO L=1,4
           PH(l)=PHEP(l,k)
          ENDDO
          IFpart=(XCOS(PH) .LT. 0.95)
          IFpart=ifpart .AND. (ENRG(PH) .GT. 0.05*cmsene/2)
          IFpart=ifpart .AND. (xacol(PH,qm,3) .GT. 20.0) .AND. (xacol(PH,qp,3) .GT. 20.0)
          IF (ifpart) egam=max(egam,ph(4))
          IF (ifpart) nph=nph+1
        ENDDO
        DO K=5,NHEP
          DO L=1,4
           PH(l)=PHEP(l,k)
          ENDDO
          IFpart=(XCOS(PH) .LT. 0.95)
          IFpart=ifpart .AND. (ENRG(PH) .GT. 0.05*cmsene/2) .AND. ph(4) .LT. egam
          IFpart=ifpart .AND. (xacol(PH,qm,3) .GT. 20.0) .AND. (xacol(PH,qp,3) .GT. 20.0)
          IF (ifpart) egam1=max(egam1,ph(4))
        ENDDO

          evis=evis+egam+egam1
          CALL pmsum(qp,qm,phsum)
          IFacc=ifacc .AND. nph .GT. 1 .AND. (evis .GT. 1.6*cmsene/2)
        IF (ifacc) THEN
          CALL memorizer4(0,KK,wt,0D0,0D0,0D0)
        ENDIF
C===============================================================
C     'Opal21'
      kk=10
        NPH=0
        ITR=0
        evis=0
        egam=0
        egam1=0
        IFacc=(xCOS(qp) .LT. 0.95 .AND. xCOS(qm) .LT. 0.95)
        IFacc=ifacc .AND. (pete(qp) .GT. 1.0 .AND. pete(qm) .GT. 1.0)
        IFacc=ifacc .AND. (xacol(qp,qm,3) .GT. 20.0)
        evis=qp(4)+qm(4)
        IFphot=.false.
        DO K=5,NHEP
          DO L=1,4
           PH(l)=PHEP(l,k)
          ENDDO
          IFpart=(XCOS(PH) .LT. 0.95)
          IFpart=ifpart .AND. (ENRG(PH) .GT. 0.05*cmsene/2)
          IFpart=ifpart .AND. (xacol(PH,qm,3) .GT. 20.0) .AND. (xacol(PH,qp,3) .GT. 20.0)
          IF (ifpart) egam=max(egam,ph(4))
          IF (ifpart) nph=nph+1
        ENDDO
        DO K=5,NHEP
          DO L=1,4
           PH(l)=PHEP(l,k)
          ENDDO
          IFpart=(XCOS(PH) .LT. 0.95)
          IFpart=ifpart .AND. (ENRG(PH) .GT. 0.05*cmsene/2) .AND. ph(4) .LT. egam
          IFpart=ifpart .AND. (xacol(PH,qm,3) .GT. 20.0) .AND. (xacol(PH,qp,3) .GT. 20.0)
          IF (ifpart) egam1=max(egam1,ph(4))
        ENDDO

          evis=evis+egam+egam1
          CALL pmsum(qp,qm,phsum)
          IFacc=ifacc .AND. (85D0 .GT. xmas(phsum) .OR. xmas(phsum) .GT. 95d0)
          IFacc=ifacc .AND. nph .GT. 1 .AND. (evis .GT. 1.6*cmsene/2)
        IF (ifacc) THEN
          CALL memorizer4(0,KK,wt,0D0,0D0,0D0)
        ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C general variables for observables of the mumu (eventually gamma)
          CALL pmsum(qp,qm,phsum)
          xminv=xmas(phsum)
          CALL xprop_view(1,phsum)
          xmprop=xmas(phsum)
          xmprop=xmas1(NHEP,PHEP)
          xmang=xmasang(cmsene,qp,qm)
C          xmang1=xmasang1(cmsene,qp,qm,PH)
          cosmu=xCOS(qm)
          sign=1
          IF(qm(3) .GT. 0d0) sign=-1


C----------------------------------------------------
C Aleph5  !! need clarifications before checking what DOes it mean photon
        KK=11
        IFacc=(xCOS(qp) .LT. 0.95 .AND. xCOS(qm) .LT. 0.95)
        IFacc=ifacc .AND. (enrg(qp) .GT. 6.0 .AND. enrg(qm) .GT. 6.0)
        IFacc=ifacc .AND. (enrg(qp)+enrg(qm) .GT. 60D0)
        xmangu=xmang
        phmax=0
        DO K=5,NHEP
          DO L=1,4
           PH(l)=PHEP(l,k)
          ENDDO
           IF (xCOS(ph) .LT. 0.97 .AND. enrg(ph) .GT. 5d0 .AND. enrg(pH) .GT. phmax) THEN
             phmax=enrg(ph)
             xmangu=xmasang1(cmsene,qp,qm,PH)
           ENDIF
        ENDDO
        IF (ifacc .AND. xmangu .GT. 0.1*cmsene) THEN
          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
        ENDIF


C----------------------------------------------------
C Aleph6 !! need clarification
        KK=12
        IFacc=(xCOS(qp) .LT. 0.95 .AND. xCOS(qm) .LT. 0.95)
        IFacc=ifacc .AND. (enrg(qp) .GT. 6.0 .AND. enrg(qm) .GT. 6.0)
        IFacc=ifacc .AND. (enrg(qp)+enrg(qm) .GT. 60D0)
        xmangu=xmang
        phmax=0
        DO K=5,NHEP
          DO L=1,4
           PH(l)=PHEP(l,k)
          ENDDO
           IF (xCOS(ph) .LT. 0.97 .AND. enrg(ph) .GT. 5d0 .AND. enrg(pH) .GT. phmax) THEN
             phmax=enrg(ph)
             xmangu=xmasang1(cmsene,qp,qm,PH)
           ENDIF
        ENDDO
        IF (ifacc .AND. xmangu .GT. 0.9*cmsene .AND. xminv .GT. 0.74*cmsene) THEN
          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
        ENDIF


C----------------------------------------------------
C Delphi4
        KK=13
        IFacc=(xCOS(qp) .LT. COS(20./180.*pi) .AND. xCOS(qm) .LT. COS(20./180.*pi))
        IFacc=ifacc .AND. (enrg(qp) .GT. 30.0 .OR. enrg(qm) .GT. 30.0)
        xmangu=xmang
        phmax=0
        IF (ifacc .AND. xminv .GT. 75) THEN
          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
        ENDIF



C----------------------------------------------------
C Delphi5
        KK=14
        IFacc=(xCOS(qp) .LT. COS(20./180.*pi) .AND. xCOS(qm) .LT. COS(20./180.*pi))
        IFacc=ifacc .AND. (enrg(qp) .GT. 30.0 .OR. enrg(qm) .GT. 30.0)
        xmangu=xmang
        phmax=0
        IF (ifacc .AND. xminv .GT. 0.85*cmsene) THEN
          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
        ENDIF


C----------------------------------------------------
C LT9
        KK=15
        IFacc=(xCOS(qp) .LT. COS(20./180.*pi) .AND. xCOS(qm) .LT. COS(20./180.*pi))
        IFacc=ifacc .AND. (enrg(qp) .GT. 35.0 .OR. enrg(qm) .GT. 35.0)
        xmangu=xmang
        phmax=0
        DO K=5,NHEP
          DO L=1,4
           PH(l)=PHEP(l,k)
          ENDDO
           IF (xCOS(ph) .LT. 0.985 .AND. enrg(ph) .GT. 15d0 .AND. enrg(pH) .GT. phmax
     $          .AND. (xacol(PH,qm,3) .GT. 10.0) .AND. (xacol(PH,qp,3) .GT. 10.)) THEN
             phmax=enrg(ph)
             xmangu=xmasang1(cmsene,qp,qm,PH)
           ENDIF
        ENDDO
        IF (ifacc .AND. xmangu .GT. 75D0) THEN
          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
        ENDIF


C----------------------------------------------------
C LT10
        KK=16
        IFacc=(xCOS(qp) .LT. COS(20./180.*pi) .AND. xCOS(qm) .LT. COS(20./180.*pi))
        IFacc=ifacc .AND. (enrg(qp) .GT. 35.0 .OR. enrg(qm) .GT. 35.0)
        xmangu=xmang
        phmax=0
        DO K=5,NHEP
          DO L=1,4
           PH(l)=PHEP(l,k)
          ENDDO
           IF (xCOS(ph) .LT. 0.985 .AND. enrg(ph) .GT. 15d0 .AND. enrg(pH) .GT. phmax
     $          .AND. (xacol(PH,qm,3) .GT. 10.0) .AND. (xacol(PH,qp,3) .GT. 10.)) THEN
             phmax=enrg(ph)
             xmangu=xmasang1(cmsene,qp,qm,PH)
           ENDIF
        ENDDO
        IF (ifacc .AND. xmangu .GT. 0.85*cmsene) THEN
          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
        ENDIF


C----------------------------------------------------
C Opal6
        KK=17
        IFacc=(xCOS(qp) .LT. 0.95 .AND. xCOS(qm) .LT. 0.95)
        IFacc=ifacc .AND. (enrg(qp) .GT. 6.0 .AND. enrg(qm) .GT. 6.0)
        IFacc=ifacc .AND. (xacol(qp,qm,2) .GT. 0.32/pi*180)
        xmangu=xmang
        phmax=0
        esum=qp(4)+qm(4)
        DO K=5,NHEP
          DO L=1,4
           PH(l)=PHEP(l,k)
          ENDDO
           IF (xCOS(ph) .LT. 0.985 .AND. enrg(ph) .GT. 0.8d0) THEN
!!!!!!!!!!!!!!!!             esum=esum+ph(4)
             IF      (enrg(pH) .GT. phmax
     $           .AND. (xacol(PH,qm,3) .GT. 0.2/pi*180)
     $           .AND. (xacol(PH,qp,3) .GT. 0.2/pi*180)) THEN
               xmangu=xmasang1(cmsene,qp,qm,PH)
             ENDIF
             IF (enrg(pH) .GT. phmax) THEN
               phmax=enrg(ph)
             ENDIF

           ENDIF
        ENDDO
        esum=esum+phmax
        IFacc=ifacc .AND. xmangu .GT. 0.1*cmsene
        IF (ifacc .AND. (
     $  (esum .GT. (0.35*cmsene+0.5*91.17**2/cmsene) .AND. esum .LT. (0.75*cmsene+0.5*91.17**2/cmsene) .AND. xminv .GT. 70D0)) 
     $   .OR. esum .LT. (0.35*cmsene+0.5*91.17**2/cmsene) .OR. esum .GT. (0.75*cmsene+0.5*91.17**2/cmsene)) THEN
          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
        ENDIF


C----------------------------------------------------
C Opal7
        KK=18
        IFacc=(xCOS(qp) .LT. 0.95 .AND. xCOS(qm) .LT. 0.95)
        IFacc=ifacc .AND. (enrg(qp) .GT. 6.0 .AND. enrg(qm) .GT. 6.0)
        IFacc=ifacc .AND. (xacol(qp,qm,2) .GT. 0.32/pi*180)
        xmangu=xmang
        phmax=0
        esum=qp(4)+qm(4)
        DO K=5,NHEP
          DO L=1,4
           PH(l)=PHEP(l,k)
          ENDDO
           IF (xCOS(ph) .LT. 0.985 .AND. enrg(ph) .GT. 0.8d0) THEN
!!!!!!!!!!!!!!!             esum=esum+ph(4)
             IF      (enrg(pH) .GT. phmax
     $           .AND. (xacol(PH,qm,3) .GT. 0.2/pi*180)
     $           .AND. (xacol(PH,qp,3) .GT. 0.2/pi*180)) THEN
               xmangu=xmasang1(cmsene,qp,qm,PH)
             ENDIF
             IF (enrg(pH) .GT. phmax) THEN
               phmax=enrg(ph)
             ENDIF
           ENDIF
        ENDDO
        esum=esum+phmax
c        WRITE(*,*) '>>>>>>>>>>>>>>>>>>>>>>'
c        DO K=3,NHEP
c         WRITE(*,1000) PHEP(1,k),PHEP(2,k),PHEP(3,k),PHEP(4,k)
c        ENDDO
c 1000    FORMAT(1x,4(f10.4,2x))
C        WRITE(*,*) cmsene,'>>',xminv, xmang, xmangu
        IF (ifacc .AND. esum .GT. (0.35*cmsene+0.5*91.17**2/cmsene) 
     $     .AND. xmangu .GT. 0.85*cmsene .AND. xminv .GT. SQRT(0.1*cmsene**2+91.17**2) ) THEN
          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
        ENDIF



C idealized observables !!!!!!!!!!!!!!!!!
C----------------------------------------------------
C IAleph5
        KK=19

        IF (cosmu .LT. 0.95D0 .AND. xminv .GT. 0.1*cmsene) THEN
          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
        ENDIF

C----------------------------------------------------
C IAleph6
        KK=20

        IF (cosmu .LT. 0.95D0 .AND. xminv .GT. 0.9*cmsene) THEN
          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
        ENDIF

C----------------------------------------------------
C IDelphi5
        KK=21

        IF (cosmu .LT. 0.95D0 .AND. xminv .GT. 75D0) THEN
          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
        ENDIF

C----------------------------------------------------
C IDelphi6
        KK=22

        IF (cosmu .LT. 0.95D0 .AND. xminv .GT. 0.85*cmsene) THEN
          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
        ENDIF

C----------------------------------------------------
C ILT9
        KK=23

        IF (cosmu .LT. 0.9D0 .AND. xminv .GT. 75) THEN
          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
        ENDIF

C----------------------------------------------------
C ILT10
        KK=24

        IF (cosmu .LT. 0.9D0 .AND. xminv .GT. 0.85*cmsene) THEN
          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
        ENDIF

C----------------------------------------------------
C ILT11
        KK=25

        IF (cosmu .LT. 1D0 .AND. xminv .GT. 0.85*cmsene) THEN
          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
        ENDIF

C----------------------------------------------------
C IOpal6
        KK=26

        IF (cosmu .LT. 0.95D0 .AND. xmprop .GT. 0.1*cmsene) THEN
          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
        ENDIF

C----------------------------------------------------
C IOpal7
        KK=27

        IF (cosmu .LT. 1D0 .AND. xmprop .GT. 0.1*cmsene) THEN
          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
        ENDIF

C----------------------------------------------------
C IOpal8
        KK=28

        IF (cosmu .LT. 0.95D0 .AND. xmprop .GT. 0.85*cmsene) THEN
          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
        ENDIF

C----------------------------------------------------
C IOpal9
        KK=29

        IF (cosmu .LT. 1D0 .AND. xmprop .GT. 0.85*cmsene) THEN
          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
        ENDIF
c################################################################
      END


C ================================================================
C ================================================================
C universal library for calculating observables for any generator
C storing final states in HEPEVT COMMON block
C in HEPEVT: FILEds 1,2 must be beams
C            fields 3,4 must be fs fermions.
C ================================================================
C ================================================================

      FUNCTION ENRG(PH)
c takes energy of the 4-vector
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PH(4)
      ENRG=PH(4)
      END

      FUNCTION XCOS(PH)
c for 4-vector calculates module of the cosine with beam 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PH(4)
      XCOS=ABS(PH(3))/SQRT(PH(1)**2+PH(2)**2+PH(3)**2)
      END

      FUNCTION PETE(PH)
c for 4-vector calculates p_t
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PH(4)
      PETE=SQRT(PH(1)**2+PH(2)**2)
      END

      FUNCTION XMAS(PH)
c for 4-vector calculates mass
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PH(4)
      XMAS=SQRT(ABS(PH(4)**2-PH(3)**2-PH(2)**2-PH(1)**2))
      XMAS1=XMAS
      END
      FUNCTION XMAS1(N,PHEP)
c calculate mass of virtual Z assuming leading log approx.
c event must have beams first fs fermions later and photons at END
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMXHEP=2000)
      REAL*8 PHEP(5,NMXHEP)
      DIMENSION qp(4),qm(4),psum(4),pp(4),pm(4),ph(4)
      AMZ=91.17
      GAMZ=2.5
      DO K=1,4
       PP(k)=PHEP(K,1)
       PM(k)=PHEP(K,2)
       QP(k)=PHEP(K,3)
       QM(k)=PHEP(K,4)
       PSUM(K)=QP(k)+QM(K)
      ENDDO
      DO L=5,N
        DO K=1,4
          PH(k)=PHEP(K,1)
        ENDDO
        XM1=XINV(PP,PH)
        XM1=MAX(XM1,XINV(PM,PH))
        XM2=XINV(QP,PH)
        XM2=MAX(XM2,XINV(QM,PH))
        XMV1=(XMAS(PSUM))**2
        XMV2=XINV(PSUM,PH)
        XM1=1D0/XM1/((XMV1-AMZ**2)**2+GAMZ**2*AMZ**2)
        XM2=1D0/XM2/((XMV2-AMZ**2)**2+GAMZ**2*AMZ**2)
        IF (XM2 .GT. XM1) THEN
          DO K=1,4
           PSUM(K)=PSUM(K)+PH(K)
          ENDDO
        ENDIF
      ENDDO
      XMAS1=XMAS(PSUM)
      END

      FUNCTION XINV(qp,qm)
C invariant calculated from 2 four-momenta
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION qp(4),qm(4)
      XINV=(qp(4)+qm(4))**2-(qp(3)+qm(3))**2-(qp(2)+qm(2))**2-(qp(1)+qm(1))**2
      XINV=ABS(XINV)
      END

      FUNCTION XMASANG(cmsene,qp,qm)
C mass calculated from directions of qp qm assuming there is just 
c one extra collinear with beam ph.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION qp(4),qm(4)
      c1=qp(3)/SQRT(qp(1)**2+qp(2)**2+qp(3)**2)
      c2=qm(3)/SQRT(qm(1)**2+qm(2)**2+qm(3)**2)
      s1=SQRT(1-c1**2)
      s2=SQRT(1-c2**2)
      t1=s1/c1
      t2=s2/c2
      tt=ABS(t1+t2)
      xk=2*tt/(tt+SQRT(1+1/t1**2)+SQRT(1+1/t2**2))
      R=ABS(c1 +c2*s1/s2)/(1+s1/s2)
      xk=2*r/(1+r)
      xmasang=cmsene*SQRT(1-xk)
      END
      FUNCTION XMASANG1(cmsene,qp,qm,PH)
C mass calculated from directions of qp qm assuming there is just one extra visible ph.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION qp(4),qm(4),ph(4)
      c1=(qp(1)*PH(1)+qp(2)*PH(2)+qp(3)*PH(3))
     $   /SQRT(qp(1)**2+qp(2)**2+qp(3)**2)/SQRT(PH(1)**2+PH(2)**2+PH(3)**2)
      c2=(qm(1)*PH(1)+qm(2)*PH(2)+qm(3)*PH(3))
     $   /SQRT(qm(1)**2+qm(2)**2+qm(3)**2)/SQRT(PH(1)**2+PH(2)**2+PH(3)**2)
      s1=SQRT(1-c1**2)
      s2=SQRT(1-c2**2)
      t1=s1/c1
      t2=s2/c2
      tt=ABS(t1+t2)
      xk=2*tt/(tt+SQRT(1+1/t1**2)+SQRT(1+1/t2**2))
      R=ABS(c1 +c2*s1/s2)/(1+s1/s2)
      xk=2*r/(1+r)
      xmasang1=cmsene*SQRT(ABS(1-xk))
      END

      SUBROUTINE PMSUM(P,Q,PH)
c adds two 4-vectors
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PH(4),P(4),Q(4)
      DO K=1,4
       PH(k)=P(k)+Q(K)
      ENDDO
      END

      FUNCTION XACOL(X,Y,N)
C     ********************
c calculates acollinearity 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8    X(*),Y(*)
      DIMENSION X1(4),Y1(4)
      DATA PI /3.1415926535897932D0/
      S=0.D0
      X2=0.D0
      Y2=0.D0
      DO 9  I=1,N
      X1(I)=X(I)
    9 Y1(I)=Y(I)
      DO 10 I=1,N
      S=S+X1(I)*Y1(I)
      X2=X2+X1(I)**2
   10 Y2=Y2+Y1(I)**2
      XACOL=ACOS(S/SQRT(X2*Y2))*180.D0/PI
      RETURN
      END
      SUBROUTINE xprop_view(mode,PP)
      IMPLICIT REAL*8 (A-H,O-Z)
c this routine is digging out of yfs3 non-obervable 4-momentum of Z.
c but can also get it from hepevt, in LL approximation
c event must have beams first fs fermions later and photons at END
c IF you DO not use guts of yfs subr. karlud comment out 1 line.
      PARAMETER (NMXHEP=2000)
      COMMON/d_HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      REAL*8 phep,vhep
      REAL*8 pol,sig,err
      SAVE  /HEPEVT/
      REAL*8 PP(4),PS(4)
      SAVE ps
      DIMENSION qp(4),qm(4),psum(4),pm(4),ph(4)
      AMZ=91.17
      GAMZ=2.5


      IF (mode .EQ. 0) THEN
        DO K=1,4
           ps(k)=pp(k)
        ENDDO
        RETURN
      ENDIF

      DO K=1,4
       PP(k)=PHEP(K,1)
       PM(k)=PHEP(K,2)
       QP(k)=PHEP(K,3)
       QM(k)=PHEP(K,4)
       PSUM(K)=QP(k)+QM(K)
       PP(k)=PSUM(K)
      ENDDO
      DO L=5,NHEP
        DO K=1,4
          PH(k)=PHEP(K,L)
        ENDDO
        XM1=XINV(PP,PH)
        XM1=MIN(XM1,XINV(PM,PH))
        XM2=XINV(QP,PH)
        XM2=MIN(XM2,XINV(QM,PH))
        XMV1=(XMAS(PSUM))**2
        XMV2=XINV(PSUM,PH)
        XM1=1D0/XM1/((XMV1-AMZ**2)**2+GAMZ**2*AMZ**2)
        XM2=1D0/XM2/((XMV2-AMZ**2)**2+GAMZ**2*AMZ**2)
        IF (XM2 .GT. XM1) THEN
          DO K=1,4
           PSUM(K)=PSUM(K)+PH(K)
           PP(k)=PSUM(K)
          ENDDO
        ENDIF
      ENDDO

C case when it could be directly taken from host PROGRAM (in mode 0)
C overWRITEs calculated in LL pp. IF it is not available
C next 3 lines should be commented out.
      DO k=1,4
        pp(k)=ps(k)
      ENDDO

      END
