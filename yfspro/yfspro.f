*//////////////////////////////////////////////////////////////////////////////
*//                  main program for MC production                          //
*//////////////////////////////////////////////////////////////////////////////
*//   makeE189GeV-start
*//
*//   make mix200-start
*//   make ini200-start
*//   make fin140-start
*//
*//   make figini
*//   make figfin
*//   make figmix
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      PROGRAM main
*     ************
      IMPLICIT NONE
!     ***********
      CALL yfspro
!     ***********
      END


      SUBROUTINE yfspro
*     **********************************
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
* input/output files
      COMMON / inout  / ninp,nout
      CHARACTER*4 semaph
!
!---------------
! general output for everybody including glibk
      ninp =5
      nout =16
      OPEN(nout,FILE='./pro.output')
      REWIND(nout)
      CALL GLK_SetNout(nout)
!
      WRITE(nout,*) '   '
      WRITE(nout,*) '=============================================='
      WRITE(nout,*) '==========*********************==============='
      WRITE(nout,*) '==========***    YFSPRO     ***==============='
      WRITE(nout,*) '==========*********************==============='
      WRITE(nout,*) '=============================================='
      WRITE(nout,*) '   '
!---------------

! READ semaphore flag
      CALL givsem(semaph)

!---------------
      IF(semaph  .EQ.  'STAR') THEN
         WRITE(6,*) ' ------- starting from the scratch ----------'
! READ initial (root) random number seed
         ninp3=3
         OPEN(ninp3,FILE='./iniseed')
         READ(ninp3,'(i10)') ijklin
         READ(ninp3,'(i10)') ntotin
         READ(ninp3,'(i10)') ntot2n
         CALL PseuMar_Initialize(ijklin,ntotin,ntot2n)
      ELSEIF(semaph  .EQ.  'CONT') THEN
         WRITE(6,*) ' ------- restoring from the disk   ----------'
! restore histograms from the disk
         ninph=10
         OPEN(ninph,FILE='./'//'pro.hst')
         CALL GLK_hrfile(ninph,' ',' ')      !Set FILE number
         CALL GLK_hrin(   0,9999,0)          !READ from the disk
         CALL GLK_hrend(' ')                 !CLOSE FILE
! READ random number seed stored in semaphore FILE
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

! *******************************************
      CALL splot
! ********************************************

      END


      SUBROUTINE splot
*     ****************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
************************************************************************
* All photons and fermions
      COMMON / momset /  pf1(4),pf2(4),qf1(4),qf2(4),sphum(4),sphot(100,4),nphot,KFfin
* ISR photons, note that xf1,xf2 are not filled
      COMMON / momini / xf1(4),xf2(4),xphum(4),xphot(100,4),nphox
************************************************************************
* ----------------------------------------------------------------------
      INTEGER nmxhep         ! maximum number of particles
      PARAMETER (nmxhep=2000)
      REAL*8  phep, vhep
      INTEGER nevhep, nhep, isthep, idhep, jmohep, jdahep
      COMMON /d_hepevt/
     $     nevhep,           ! serial number
     $     nhep,             ! number of particles
     $     isthep(nmxhep),   ! status code
     $     idhep(nmxhep),    ! particle ident KF
     $     jmohep(2,nmxhep), ! parent particles
     $     jdahep(2,nmxhep), ! childreen particles
     $     phep(5,nmxhep),   ! four-momentum, mass [GeV]
     $     vhep(4,nmxhep)    ! vertex [mm]
* ----------------------------------------------------------------------
      COMMON / inout  / ninp,nout
      SAVE

      PARAMETER( imax = 10000)
      REAL*8   xpar(imax)

      CHARACTER*4 semaph
      CHARACTER*80 DiskFile
*=============================================================================
      WRITE(nout,*) '   '
      WRITE(nout,*) '=============================================='
      WRITE(nout,*) '==========*********************==============='
      WRITE(nout,*) '==========***  Y F S P R O  ***==============='
      WRITE(nout,*) '==========*********************==============='
      WRITE(nout,*) '=============================================='
      WRITE(nout,*) '   '
*=============================================================================
* Read data for main program
      OPEN( ninp,FILE='./pro.input')
      READ(ninp,'(8i2)') kat1,kat2,kat3,kat4,kat5,kat6,kat7,kat8
      WRITE(nout,'(8a6/8i6)')
     $ 'kat1','kat2','kat3','kat4','kat5','kat6','kat7','kat8',
     $  kat1 , kat2 , kat3 , kat4 , kat5 , kat6 , kat7 , kat8
      READ(ninp,'(i10)') nevt
      CLOSE(ninp)
      WRITE(   6,*)   nevt,' requested events '
      WRITE(nout,*)   nevt,' requested events '
*=============================================================================
      CALL KK2f_ReaDataX('../../.KK2f_defaults', 1,imax,xpar)  ! reading general defaults
      CALL KK2f_ReaDataX(         './pro.input', 0,imax,xpar)  ! reading user input
      CALL KK2f_Initialize(xpar)                               ! initialize generator
*=============================================================================
      IF(kat1   .EQ.   1) CALL robol1(-1,xpar)
      IF(kat2   .EQ.   1) CALL robol2(-1,xpar)
      IF(kat3   .EQ.   1) CALL robol3(-1,xpar)
      IF(kat4   .EQ.   1) CALL robol4(-1,xpar)
      IF(kat5   .EQ.   1) CALL robol5(-1,xpar)
!-------------------------------------------------------!
!                 main MC loop                          !
!-------------------------------------------------------!
      ngroup = 100000
      ngroup = 10000
      iev=0
      DO loop=1,10000000
        DO igroup =1,ngroup
          iev=iev+1
          IF(MOD(iev, ngroup)   .EQ.   1) WRITE( 6,*)  'iev= ',iev
!         *************************
          CALL KK2f_Make
!         ****************************
*---------------------------------------------------------------
* Fill-in /momset/ and /momini/ common blocks for the local use
* fermions
          DO k=1,4
             pf1(k) =phep(k,1)
             pf2(k) =phep(k,2)
             qf1(k) =phep(k,3)
             qf2(k) =phep(k,4)
             sphum(k)=0d0
             xphum(k)=0d0
          ENDDO
          KFfin  = idhep(3)
*     and photons
          nphot=0
          nphox=0
          DO j=1,100
             ih = 4+j
             kf = idhep(ih)
*     stop if /hepevt/ ended or non-photon found
             IF(ih .GT. nhep) GOTO 110
             IF(kf .NE. 22)   GOTO 110
*     All photons here
             nphot=nphot+1
             DO  k=1,4
                sphot(nphot,k) =phep(k,ih)
                sphum(k) =sphum(k)+sphot(nphot,k)
             ENDDO
*     ISR photons separately
             IF(jmohep(1,ih) .EQ. 1) THEN
                nphox=nphox+1
                DO  k=1,4
                   xphot(nphox,k) =phep(k,ih)
                   xphum(k) =xphum(k)+xphot(nphox,k)
                ENDDO
             ENDIF
          ENDDO
 110      CONTINUE
*---------------------------------------------------------------
*   Control printout
          CALL momprt(' YFSPRO ', 6,iev,1,20,pf1,pf2,qf1,qf2,nphot,sphot,KFfin)
*         CALL momprt('*YFSPRO*',16,iev,1,20,pf1,pf2,qf1,qf2,nphot,sphot,KFfin)
*         CALL dumpri('*momini*', 6,iev,1,10,xf1,xf2,nphox,xphot)
*         CALL dumpri('*momini*',16,iev,1,10,xf1,xf2,nphox,xphot)

          IF(iev .LE. 10) THEN
             CALL pygive("MSTU(11)=16")
             CALL pylist(1)
          ENDIF
          IF(iev .LE. 10) THEN
             CALL pygive("MSTU(11)=6")
             CALL pylist(1)
          ENDIF

!         ============================================
!         histograming
          IF(kat1  .EQ.  1) CALL robol1( 0,xpar)
          IF(kat2  .EQ.  1) CALL robol2( 0,xpar)
          IF(kat3  .EQ.  1) CALL robol3( 0,xpar)
          IF(kat4  .EQ.  1) CALL robol4( 0,xpar)
          IF(kat5  .EQ.  1) CALL robol5( 0,xpar)
!         ============================================
!         check on requested no. of events
          IF(iev  .EQ.  nevt)     GOTO 300
        ENDDO
!       check on semaphore flag
        CALL givsem(semaph)
        IF(semaph  .EQ.  'STOP') GOTO 300
!       dump partial results on the disk after every ngroup
        CALL  dumpeh(iev)
      ENDDO
 300  CONTINUE
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      WRITE(6,*) ' generation finished '
      CALL KK2f_Finalize
*------------------------------------------------
      IF(kat1  .EQ.  1) CALL robol1(1,xpar)
      IF(kat2  .EQ.  1) CALL robol2(1,xpar)
      IF(kat3  .EQ.  1) CALL robol3(1,xpar)
      IF(kat4  .EQ.  1) CALL robol4(1,xpar)
      IF(kat5  .EQ.  1) CALL robol5(1,xpar)
*
      CALL  dumpeh(iev)
!     +++++++++++++++++
      END


      SUBROUTINE givsem(semaph)
*     ************************
      IMPLICIT REAL*8(a-h,o-z)
      CHARACTER*4 semaph
! ------------------------------------------------------
! READ semaphore flag
! ------------------------------------------------------
      ninp2=2
      OPEN(ninp2,FILE='./semaphore')
      READ(ninp2,'(a4)') semaph
      CLOSE(ninp2)
      END

      SUBROUTINE dumpeh(nev)
*     ************************
      IMPLICIT REAL*8(a-h,o-z)
! ------------------------------------------------------
! WRITE histos on the disk
! ------------------------------------------------------
      nouth=11
      OPEN(nouth,FILE='./pro.hst')
      CALL GLK_hrfile(nouth,' ','n')   !transfer FILE number
      CALL GLK_hrout( 0,icy,' ')       !WRITE on the disk
      CALL GLK_hrend(' ')              !CLOSE FILE
! ------------------------------------------------------
! overWRITE semaphore FILE flag
      ninp2=2
      OPEN(ninp2,FILE='./semaphore')
      WRITE(ninp2,'(a4)') 'CONT'
! append semaphore FILE with new random number seed in
      CALL PseuMar_Out(ijklin,ntotin,ntot2n)
      WRITE(ninp2,'(i10,a)') ijklin, ' = ijklin '
      WRITE(ninp2,'(i10,a)') ntotin, ' = ntotin '
      WRITE(ninp2,'(i10,a)') ntot2n, ' = ntot2n '
      WRITE(ninp2,'(i10,a)') nev,    ' =    nev '
      CLOSE(ninp2)
! ------------------------------------------------------
      END


      SUBROUTINE robol1(mode,xpar)
*/////////////////////////////////////////////////////////////////////////////////////////
*//                                                                                     //
*//   ISR study of v=1-s'/s distribution,                                               //
*//   total MC result and beta contributions in O(alf0,1,2,3)                           //
*//   Provides input for figini                                                         //
*//                                                                                     //
*/////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      REAL*8   xpar(*)
      REAL*8 pi
      PARAMETER( pi=3.1415926535897932d0)
      COMMON / inout  / ninp,nout
      COMMON / momset / pf1(4),pf2(4),qf1(4),qf2(4),sphum(4),sphot(100,4),nphot,KFfin
      REAL*8  WtSet(1000),WtList(1000),WtList2(1000)
      REAL*8  pk1(4), pk2(4), pk(4)
      SAVE
*
      IF(mode  .EQ.  -1) THEN
*     =======================
* Z peak in final mass
      idm  = 10000
      nbt  = 60
      AMZ  = 91.187D0
      tmax = AMZ +3d0
      tmin = AMZ -3d0
* Rapidity
      ida  = 30000
      nba  = 60
      rmin= 1.57d0 -0.075
      rmax= 1.57d0 +0.075
* linear scale
      idb  = 50000
      nbv  = 100
      vmin = 0
      vmax = 1
* new log scale
      idz  = 70000
      nbz  =   100
      zmin =  -3d0
      zmin =  -10d0
      zmax =   7d0
      DO k= 1, 8
         CALL GLK_Book1(ida+k,'rapidity                  $',nba,rmin,rmax)
         CALL GLK_Book1(idm+k,'ff mass                   $',nbt,tmin,tmax)
      ENDDO
      DO k= 1, 7
         CALL GLK_Book1(idb+k,'rho(log(v)) totals, diffs $',nbv,vmin,vmax)
         CALL GLK_Book1(idz+k,'rho(log(v)) totals, diffs $',nbz,zmin,zmax)
      ENDDO
      DO k=10,11
         CALL GLK_Book1(idb+k,'rho(log(v)) Betas O(alf1) $',nbv,vmin,vmax)
         CALL GLK_Book1(idz+k,'rho(log(v)) Betas O(alf1) $',nbz,zmin,zmax)
      ENDDO
      DO k=20,22
         CALL GLK_Book1(idb+k,'rho(log(v)) Betas O(alf2) $',nbv,vmin,vmax)
         CALL GLK_Book1(idz+k,'rho(log(v)) Betas O(alf2) $',nbz,zmin,zmax)
      ENDDO
      DO k=30,33
         CALL GLK_Book1(idb+k,'rho(log(v)) Betas O(alf3) $',nbv,vmin,vmax)
         CALL GLK_Book1(idz+k,'rho(log(v)) Betas O(alf3) $',nbz,zmin,zmax)
      ENDDO
      DO k=201,202
         CALL GLK_Book1(idb+k,'GPS new exponentiation $',nbv,vmin,vmax)
      ENDDO
*============================
* technical tests
      DO k=50,55
         CALL GLK_Book1(idb+k,'rho(z) beta0, technical  $',nbv,vmin,vmax)
         CALL GLK_Book1(idz+k,'rho(z) beta0, technical  $',nbz,zmin,zmax)
      ENDDO
* UNEXP
      DO k=160,166
         CALL GLK_Book1(idz+k,'rho(z) beta0, technical  $',nbz,zmin,zmax)
      ENDDO

      ELSEIF(mode  .EQ.  0) THEN
*     ==========================
      CMSene = pf1(4)+pf2(4)

      CALL KK2f_GetWtAll(WtMain,WtCrud,WtSet)
      DO k=1,1000
         WtList(k)=WtSet(k)*WtCrud
      ENDDO
*
      ffm  = 0d0
c[[[[[
      CALL BornV_GetVV(vv)
      CALL vv2zz(vv,z10)
c[[[[[      vv   = 0d0
      thet1 = acos( qf1(3)/qf1(4))
      thet2 = acos(-qf2(3)/qf2(4))
      acol =0d0
      rapt =0d0
      IF( (nphot.NE.0) .AND. (wtcrud .NE. 0d0) ) THEN
         s1=   (qf1(4)+qf2(4))**2 -(qf1(3)+qf2(3))**2
     $        -(qf1(2)+qf2(2))**2 -(qf1(1)+qf2(1))**2
         ffm= sqrt(s1)
c[[[[[        vv = 1d0 -s1/cmsene**2
c[[[[[        CALL vv2zz(vv,z10)
c[[[[[        write(*,*) ' =====robol1==>>>> vv= ',vv,nphot,WtCrud,sqrt(s1),s1/cmsene**2
         acol = abs(thet1-thet2)
         rapt = abs( log( (qf2(4)-qf2(3)) / (qf1(4)+qf1(3)) ) )
      ENDIF
*//////////////////////////////////////////////////////////////////////////////
*//      GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS     //
      DO k=201,202
         CALL GLK_Fil1(idb+k,  vv,WtList(k))
      ENDDO
*//////////////////////////////////////////////////////////////////////////////
* Total
      DO k=1,4
         CALL GLK_Fil1(ida+k,rapt,WtList(k)) !!!new
         CALL GLK_Fil1(idm+k, ffm,WtList(k))
         CALL GLK_Fil1(idb+k,  vv,WtList(k))
         CALL GLK_Fil1(idz+k, z10,WtList(k))
      ENDDO
* Total, differences O(alf^(j))-O(alf^(j-1))
      DO k=2,4
         CALL GLK_Fil1(ida+k+3,rapt,(WtList(k)-WtList(k-1))) !!!new
         CALL GLK_Fil1(idm+k+3, ffm,(WtList(k)-WtList(k-1)))
         CALL GLK_Fil1(idb+k+3,  vv,(WtList(k)-WtList(k-1)))
         CALL GLK_Fil1(idz+k+3, z10,(WtList(k)-WtList(k-1)))
      ENDDO
* Betas O(alf1)
      DO k=10,11
         CALL GLK_Fil1(idb+k, vv,WtList(k))
         CALL GLK_Fil1(idz+k,z10,WtList(k))
      ENDDO     
* Betas O(alf2)
      DO k=20,22
         CALL GLK_Fil1(idb+k, vv,WtList(k))
         CALL GLK_Fil1(idz+k,z10,WtList(k))
      ENDDO     
* Betas O(alf3)
cc      DO k=30,33
cc         CALL GLK_Fil1(idb+k, vv,WtList(k))
cc         CALL GLK_Fil1(idz+k,z10,WtList(k))
cc      ENDDO
* Betas O(alf3)-O(alf2)
      DO k=30,33
         CALL GLK_Fil1(idb+k, vv,(WtList(k)-WtList(k-10)))
         CALL GLK_Fil1(idz+k,z10,(WtList(k)-WtList(k-10)))
      ENDDO

*******************************************************
********* Technical tests at beta0 level **************
*******************************************************
      CALL BornV_GetVV(vv)

      CALL vv2zz(vv,z10)

      KeyISR = xpar(20)
      CALL GLK_Fil1(idz+50,z10,WtList(1)) ! true beta00
      CALL GLK_Fil1(idz+51,z10,wtcrud)   ! crude, without YFS formfactor!
      CALL GLK_Fil1(idb+50,vv ,WtList(1)) ! true beta00
      CALL GLK_Fil1(idb+51,vv ,wtcrud)   ! crude, without YFS formfactor!
****   UNEXP UNEXP UNEXP NEXP ****
      DO k=160,162
         CALL GLK_Fil1(idz+k,z10,WtList(k))
      ENDDO
*=======================================================================
*=======  SALAMI SALAMI SALAMI SALAMI SALAMI SALAMI SALAMI   ===========
*=======================================================================
* Notes:
*     It seems that at present (contrary to september) wtcrud
*     calculated localy is different from the one supplied
*     by generator, it requires debugging!!!
*     Found: wt_isr is now divided my wtmax in karlud
*-----------------------------------------------------------------------
      IF( KeyISR .EQ. 1) THEN
*! Guide on ISR weights:
*     ! wtcrud = wt_ISR = wtves*wtini
*     ! wtini  = wt_mas*wt_dil*wt_cut
*     ! WtList(1) = forini*beta00/discru
* weights from inside of ISR generator are available through ypar:
         CALL KK2f_GetOneY(250,wtves)
         CALL KK2f_GetOneY(251,wtini)
         CALL KK2f_GetOneY(252,wt_mas)
         CALL KK2f_GetOneY(253,wt_dil)
         CALL KK2f_GetOneY(254,wt_cut)
         CALL KK2f_GetOneY(255,wt_KF)
         CALL KK2f_GetOneY(263,wt_dil0)
         CALL KK2f_GetOneY(264,wt_cut0)
* this is pure vesko, only dilatation factor removed
* wt_KF because crude normalization multiplied by WtMax
         wtves = wtves*wt_KF
c[[[         CALL GLK_Fil1(idz+52,z10,wtves*wt_dil0)
c[[[         CALL GLK_Fil1(idb+52,vv ,wtves*wt_dil0)
         CALL GLK_Fil1(idz+52,z10,wtves)
         CALL GLK_Fil1(idb+52,vv ,wtves)
* mass weight and dilatation factor
         CALL GLK_Fil1(idz+53,z10,wtves*wt_dil0*wt_mas)
         CALL GLK_Fil1(idb+53,vv ,wtves*wt_dil0*wt_mas)
* mass weight, dilatation factor and IR-cutoff effect
* Mock wtcrud = wt1
         wt1 = wtves*wt_mas *wt_dil0*wt_cut0
         CALL GLK_Fil1(idz+54,z10,wt1)
         CALL GLK_Fil1(idb+54,vv ,wt1)
* Standard wtcrud = wt2 (as in [51])
         wt2 = wtves*wt_mas *wt_dil*wt_cut
         wtdif1 = wt2-wt1
***         wtdif1 = wtves*wt_mas *(wt_dil*wt_cut0 -wt_dil0*wt_cut0) !ok
***         wtdif1 = wtves*wt_mas *(wt_dil*wt_cut  -wt_dil0*wt_cut)
         CALL GLK_Fil1(idz+55,z10,wtdif1)
         CALL GLK_Fil1(idb+55,vv ,wtdif1)
*------------------------------------------------------
*             ***   UNEXP   ***
*------------------------------------------------------
***      wtdif1 = wtves*wt_mas *(wt_dil*wt_cut0 -wt_dil0*wt_cut0) !ok
         wtdif1 = wtves*wt_mas *(wt_dil*wt_cut  -wt_dil0*wt_cut0)
         CALL GLK_Fil1(idz+163,z10,wtdif1*wtset(162))
c[[[[[[[[[[[[[[[[[[
         wt_dilX = 1d0
         IF(nphot .EQ. 1) THEN
            djac0 = (1d0+1d0/sqrt(1d0-vv))/2d0
            wt_dilX = 1/djac0
         ELSEIF(nphot .EQ. 2) THEN
            DO k=1,4
               pk1(k)=sphot(1,k)
               pk2(k)=sphot(2,k)
               pk(k) = pk1(k)+pk2(k)
            ENDDO
            ppdpp = cmsene**2
***         pkdpk =pk(4)**2-pk(3)**2-pk(2)**2-pk(1)**2  ! original
***         pkdpk =
***  $     2*(pk1(4)*pk2(4)-pk1(3)*pk2(3)-pk1(2)*pk2(2)-pk1(1)*pk2(1))
            pkdpk = 2*(pk1(4)*pk2(4)-pk1(3)*pk2(3) ) !trans. mom. neglected
            ppdpk = cmsene*pk(4)
            aa    = ppdpp*pkdpk/(ppdpk)**2
            djac  = (1d0+1d0/sqrt(1d0-vv*aa))/2d0
            djac0 = (1d0+1d0/sqrt(1d0-vv))/2d0
            wt_dilX = djac/djac0
         ENDIF
****     if(vv.gt. 0.99999 .and. nphot .le. 2) 
****  $     write(*,*) nphot,wt_dilX,wt_dil,wt_dilX/wt_dil
         wtdif2 = wtves*wt_mas *(wt_dilX*wt_cut  -wt_dil0*wt_cut)
         CALL GLK_Fil1(idz+164,z10,wtdif2*wtset(162))
c]]]]]]]]]]]]]]]]]]
      ELSEIF( KeyISR .EQ. 2) THEN
c[[[[[[[[[[[[[[[[[[
* very special tests for forward-forward, forward-backward photon emission
         IF(nphot .EQ. 0) THEN
            vvnll = 0d0
            vvref = 0d0
         ELSEIF(nphot .EQ. 1) THEN
            vvnll = 2*sphot(1,4)/cmsene
            vvref = 2*sphot(1,4)/cmsene
         ELSEIF(nphot .EQ. 2) THEN
            DO k=1,4
               pk1(k)=sphot(1,k)
               pk2(k)=sphot(2,k)
               pk(k) = pk1(k)+pk2(k)
            ENDDO
            ppdpp = cmsene**2
            pkdpk = 2*(pk1(4)*pk2(4)-pk1(3)*pk2(3) ) !trans. mom. neglected
            ppdpk = cmsene*pk(4)
            vvnll = (2*ppdpk - pkdpk)/ppdpp
            vvref = (2*ppdpk)/ppdpp
         ENDIF
         CALL vv2zz(vvnll,znll)
         CALL vv2zz(vvref,zref)
         CALL GLK_Fil1(idz+164,znll,wtset(162))
         wta = wtset(162)
         CALL GLK_Fil1diff(idz+165,z10,wta,zref,wta)
         CALL GLK_Fil1diff(idz+166,z10,wta,znll,wta)
c]]]]]]]]]]]]]]]]]]
      ENDIF
*******************************************************
****END of Technical tests at beta0 level *************
*******************************************************
*------------------------------------------------------
* Virtual Variation of MZ
      CALL BornV_GetMZ(aMZ)
* redefine MZ
      bMZ = aMZ + 0.010d0
      CALL BornV_SetMZ(bMZ)
* recalculate weights
      IF(WtCrud .NE. 0d0) CALL QED3_Make
      CALL QED3_GetWtSet(  WtBest,WtSet )
      DO k=1,1000
         WtList2(k)=WtSet(k)*WtCrud
      ENDDO
      CALL GLK_Fil1(ida+8,rapt,(WtList2(2)-WtList(2))) !!!O(alf1)
      CALL GLK_Fil1(idm+8, ffm,(WtList2(2)-WtList(2))) !!!O(alf1)
      CALL BornV_SetMZ(aMZ)     !reset to original value
*------------------------------------------------------
      ELSE
*     ===========
         WRITE(nout,*) '==================================='
         WRITE(nout,*) '============ robol1 ==============='
         WRITE(nout,*) '==================================='
*----------------------
cccc         CALL GLK_PrintAll
*----------------------
      ENDIF
*     =====
      END

      SUBROUTINE vv2zz(vv,zz)
*     *********************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
*
      IF( vv .LE. 0d0 ) THEN
         zz  = -10000000.
      ELSEIF( vv .GT. 0d0 .AND. vv .LT. 1d0) THEN
         zz= log10(vv/(1-vv))
      ELSE
         zz  = +10000000.         
      ENDIF
      END

      SUBROUTINE robol2(mode,xpar)
*/////////////////////////////////////////////////////////////////////////////////////////
*//                                                                                     //
*//  FINAL STATE study of v=1-s'/s distribution                                         //
*//  v*rho(v) distribution as a function of log(v)                                      //
*//  total mc result and beta contributions in o(alf0,1,2)                              //
*//                                                                                     //
*/////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      REAL*8   xpar(*)
      COMMON / inout  / ninp,nout
      COMMON / momset / pf1(4),pf2(4),qf1(4),qf2(4),sphum(4),sphot(100,4),nphot,KFfin
      REAL*8  WtSet(1000),WtList(1000)
      SAVE
*
      IF(mode  .EQ.  -1) THEN
*     ===================
      tmax = 0.0
      tmin = -4.0
      nbt  = 80
      ida  = 20000
      CALL GLK_Book1(ida+300,'v*rho(v) o(alf0)         $',nbt,tmin,tmax)
      CALL GLK_Book1(ida+301,'v*rho(v) o(alf1)         $',nbt,tmin,tmax)
      CALL GLK_Book1(ida+302,'v*rho(v) o(alf2)         $',nbt,tmin,tmax)
      CALL GLK_Book1(ida+305,'v*rho(v) o(alf1)-o(alf0) $',nbt,tmin,tmax)
      CALL GLK_Book1(ida+306,'v*rho(v) o(alf2)-o(alf1) $',nbt,tmin,tmax)
      CALL GLK_Book1(ida+310,'v*rho(v)  beta0  o(alf1) $',nbt,tmin,tmax)
      CALL GLK_Book1(ida+311,'v*rho(v)  beta1  o(alf1) $',nbt,tmin,tmax)
      CALL GLK_Book1(ida+320,'v*rho(v)  beta0  o(alf2) $',nbt,tmin,tmax)
      CALL GLK_Book1(ida+321,'v*rho(v)  beta1  o(alf2) $',nbt,tmin,tmax)
      CALL GLK_Book1(ida+322,'v*rho(v)  beat2  o(alf2) $',nbt,tmin,tmax)
      nbv = 50
      vmin= 0
      vmax= 1
      idb  = 50000
      DO k=201,202
         CALL GLK_Book1(idb+k,'GPS new exponentiation $',nbv,vmin,vmax)
      ENDDO
      CALL GLK_Book1(idb+300,'rho(v) ini O(0)        $',nbv,vmin,vmax)
      CALL GLK_Book1(idb+301,'rho(v) ini O(1)        $',nbv,vmin,vmax)
      CALL GLK_Book1(idb+302,'rho(v) ini O(2)        $',nbv,vmin,vmax)
      CALL GLK_Book1(idb+305,'rho(v) ini O(1-0)      $',nbv,vmin,vmax)
      CALL GLK_Book1(idb+306,'rho(v) ini O(2-1)      $',nbv,vmin,vmax)
      CALL GLK_Book1(idb+310,'rho(v)  beta0  o(alf1) $',nbv,vmin,vmax)
      CALL GLK_Book1(idb+311,'rho(v)  beta1  o(alf1) $',nbv,vmin,vmax)
      CALL GLK_Book1(idb+320,'rho(v)  beta0  o(alf2) $',nbv,vmin,vmax)
      CALL GLK_Book1(idb+321,'rho(v)  beta1  o(alf2) $',nbv,vmin,vmax)
      CALL GLK_Book1(idb+322,'rho(v)  beat2  o(alf2) $',nbv,vmin,vmax)

      ELSEIF(mode  .EQ.  0) THEN
*     ======================
      CMSene = pf1(4)+pf2(4)

      CALL KK2f_GetWtAll(WtMain,WtCrud,WtSet)
      DO k=1,1000
         WtList(k)=WtSet(k)*WtCrud
      ENDDO
      wt1  = WtList(71)
      wt2  = WtList(72)
      wt3  = WtList(73)
      wt10 = WtList(80)
      wt11 = WtList(81)
      wt20 = WtList(90)
      wt21 = WtList(91)
      wt22 = WtList(92)
*---------//////----->>>>
      wt=wtMOD
      vv = 0d0
      t10= -100000000
      IF(nphot.ne.0) THEN
        s1=   (qf1(4)+qf2(4))**2 -(qf1(3)+qf2(3))**2
     $     -(qf1(2)+qf2(2))**2 -(qf1(1)+qf2(1))**2
        vv = 1d0 -s1/cmsene**2
        t10= log10(vv)
      ENDIF
      CALL GLK_Fil1(ida+300,t10,wt1)
      CALL GLK_Fil1(ida+301,t10,wt2)
      CALL GLK_Fil1(ida+302,t10,wt3)
      CALL GLK_Fil1(ida+305,t10,wt2-wt1)
      CALL GLK_Fil1(ida+306,t10,wt3-wt2)
      CALL GLK_Fil1(ida+310,t10,wt10)
      CALL GLK_Fil1(ida+311,t10,wt11)
      CALL GLK_Fil1(ida+320,t10,wt20)
      CALL GLK_Fil1(ida+321,t10,wt21)
      CALL GLK_Fil1(ida+322,t10,wt22)

*//////////////////////////////////////////////////////////////////////////////
*//      GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS GPS     //
      DO k=201,202
         CALL GLK_Fil1(idb+k,  vv,WtList(k))
      ENDDO
c[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(*,'(a,9g20.10)') 
c     $ '------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%------------'
c      write(*,'(a,9g20.10)') 
c     $ ' WtList(201),WtList(300) =', WtList(201),WtList(71), WtList(201)/WtList(71)
c]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
*//////////////////////////////////////////////////////////////////////////////
      CALL GLK_Fil1(idb+300, vv,wt1)
      CALL GLK_Fil1(idb+301, vv,wt2)
      CALL GLK_Fil1(idb+302, vv,wt3)
      CALL GLK_Fil1(idb+305, vv,wt2-wt1)
      CALL GLK_Fil1(idb+306, vv,wt3-wt2)
      CALL GLK_Fil1(idb+310, vv,wt10)
      CALL GLK_Fil1(idb+311, vv,wt11)
      CALL GLK_Fil1(idb+320, vv,wt20)
      CALL GLK_Fil1(idb+321, vv,wt21)
      CALL GLK_Fil1(idb+322, vv,wt22)

      ELSE
*     ===========
      WRITE(nout,*) '==================================='
      WRITE(nout,*) '============ robol2 ==============='
      WRITE(nout,*) '==================================='
*----------------------
***         CALL GLK_PrintAll
*----------------------
      ENDIF
*     =====
      END


      SUBROUTINE robol3(mode,xpar)
*/////////////////////////////////////////////////////////////////////////////////////////
*//                                                                                     //
*//  sussex related plots. initial state.                                               //
*//  dependence on energy cut on second photon,                                         //
*//  study on energy/pt of the hardest 3 photons, initial state.                        //
*//                                                                                     //
*/////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      REAL*8   xpar(*)
      COMMON / inout  / ninp,nout
      COMMON / momset / pf1(4),pf2(4),qf1(4),qf2(4),sphum(4),sphot(100,4),nphot,KFfin
      REAL*8  WtSet(1000),WtList(1000)
      SAVE
      REAL*8   soflim(20)
      REAL*8   pt(100),ind(100)
      DATA lsof, soflim / 10,
     $  1d0,  0.100d0,  0.030d0,  0.010d0,  0.003d0, 0.001d0,
     $ 0.0003d0, 0.0001d0, 0.00003d0, 0.00001,    10*0d0 /
      SAVE lsof,soflim,ide,idt,nbx,xmin,xmax
*
      IF(mode  .EQ.  -1) THEN
*     ===================
      ide = 10050
      idt = 10060
      CALL GLK_Book1(10010,'photon multiplicity all      $',30,0d0,30d0)
      nbx =    70
      xmin= -14.0
      xmax=   0.0
      DO j=1,3
        CALL GLK_Book1(ide+j,'log10(e/eeam) $',nbx,xmin,xmax)
        CALL GLK_Book1(idt+j,'log10(pt/eeam)$',nbx,xmin,xmax)
      ENDDO
* two-dim. controll scatteregram of sudakov variables
*      zz= -20d0
*      CALL GLK_Book2(10090,' scatt. log10 $',40,zz,0d0, 80,zz,0d0)

      ELSEIF(mode  .EQ.  0) THEN
*     ======================
      CMSene = pf1(4)+pf2(4)
      CALL KK2f_GetWtAll(WtMain,WtCrud,WtSet)
      DO k=1,1000
         WtList(k)=WtSet(k)*WtCrud
      ENDDO
      wt=WtMain
*...ordering pt
      CALL ordpt(pt,ind)
      DO j=1,3
         ze = -1000
         zpt= -1000
         IF(nphot   .GE.   j) THEN
*     ...hardest photons  (no need to order!?)
            ze=log10(2*sphot(j,4)/cmsene)
*     ...photons with highest pt
            zpt=log10(2*pt(j)/cmsene)
         ENDIF
         CALL GLK_Fil1(ide+j, ze, wt)
         CALL GLK_Fil1(idt+j,zpt, wt)
      ENDDO
*
* controll scattergram
*      DO 27 i=1,nphot
*      xl=  log10((sphot(i,4)+sphot(i,3))/cmsene)
*      yl=  log10((sphot(i,4)-sphot(i,3))/cmsene)
*      CALL hfill(10090,xl,yl, wt)
   27 CONTINUE
      ELSE
*     ===========
      WRITE(nout,*) '==================================='
      WRITE(nout,*) '============ robol3 ==============='
      WRITE(nout,*) '==================================='
*----------------------
***         CALL GLK_PrintAll
*----------------------
      ENDIF
*     =====
      END

      SUBROUTINE ordpt(pt,ind)
*     ************************
* ordering photons according to pt, pt is list of ordered pt
* and ind is the adress list of photons (indices in sphot)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      REAL*8  pt(*),ind(*)
      COMMON / momset / pf1(4),pf2(4),qf1(4),qf2(4),sphum(4),sphot(100,4),nphot,KFfin
      SAVE
      DO 10 i=1,nphot
      pt(i)=sqrt(sphot(i,1)**2+sphot(i,2)**2)
   10 ind(i)=i
      IF(nphot  .LE.  1) RETURN
      DO 30 i=2,nphot
      DO 30 j=nphot,i,-1
      IF(pt(j)  .GT.  pt(j-1)) THEN
        x=pt(j)
        pt(j)=pt(j-1)
        pt(j-1)=x
        l=ind(j)
        ind(j)=ind(j-1)
        ind(j-1)=l
      ENDIF
   30 CONTINUE
      END

      SUBROUTINE robol4(Mode,xpar)
*/////////////////////////////////////////////////////////////////////////////////////////
*//                                                                                     //
*//   Study of ISR+FSR                                                                  //
*//   MC input for figmix                                                               //
*//                                                                                     //
*/////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION   xpar(*)
      INTEGER            Mode
*------------------------------------------------------------------------------------
      INTEGER           ninp,nout
      COMMON / inout  / ninp,nout
*------------------------------------------------------------------------------------
      INTEGER           nphox
      DOUBLE PRECISION  xf1,   xf2,   xphum,   xphot
      COMMON / momini / xf1(4),xf2(4),xphum(4),xphot(100,4),nphox
*------------------------------------------------------------------------------------
      INTEGER           nphot,KFfin
      DOUBLE PRECISION  pf1,   pf2,   qf1,   qf2,   sphum,   sphot
      COMMON / momset / pf1(4),pf2(4),qf1(4),qf2(4),sphum(4),sphot(100,4),nphot,KFfin
*------------------------------------------------------------------------------------
      DOUBLE PRECISION  WtSet(1000),WtList(1000),WtList2(1000)
      DOUBLE PRECISION  xxf(4)
*------------------------------------------------------------------------------------
      DOUBLE PRECISION  CMSene,WtMain,WtCrud,Flip,aMZ,bMZ
      DOUBLE PRECISION  s1,costh,rapt,WtBest
      INTEGER           iev,k,i,j
*------------------------------------------------------------------------------------
      DOUBLE PRECISION  cmin,cmax,cost1,cost2,yy1,yy2,sint1,sint2
      INTEGER           idc,nbc
      DOUBLE PRECISION  vmin,vmax,t10,z10,vvB
      INTEGER           idb,nbv
      DOUBLE PRECISION  rmin,rmax
      INTEGER           ida,nba
      DOUBLE PRECISION  zmin,zmax
      INTEGER           idz,nbz
      DOUBLE PRECISION  smin,smax,vvAL,theta1,theta2
      INTEGER           ids,nbs
      DOUBLE PRECISION  vvPR
* energy cut parameters
      DOUBLE PRECISION  vvALmax,vvBcut
      PARAMETER (vvBcut  = 0.9d0,  vvALmax = 0.2775d0   ) ! ALEPH exper. cut
      DOUBLE PRECISION  vvWGmax1,vvWGmax2,vvWGmax9
      PARAMETER (vvWGmax1 = 0.1900d0 ) ! LEP2 WG cut: 1-0.85^2
      PARAMETER (vvWGmax2 = 0.2775d0 ) ! LEP2 WG cut: 1-0.90^2
      PARAMETER (vvWGmax9 = 0.9000d0 ) ! LEP2 WG cut, ZRR included
*---- robol4, table of KF-Cut
      INTEGER           nCut, iCut, nKF, iKF, idf, iBin
      PARAMETER (       nCut = 11, nKF = 16 )
      DOUBLE PRECISION  VmaxTab(nCut)
      DATA VmaxTab / 0.01d0, 0.10d0, 0.20d0, 0.30d0, 0.40d0, 0.50d0,
     $                       0.60d0, 0.70d0, 0.80d0, 0.90d0, 0.99d0/
      DOUBLE PRECISION  GLK_hi,GLK_hie
      INTEGER           KeyWgt
      SAVE
*------------------------------------------------------------------------------------
      IF(Mode  .EQ.  -1) THEN
*     =======================
* Flavour/cut spreadsheet
      idf  = 20000
* Costheta
      idc  = 10000
      nbc  = 40
      cmin= -1d0
      cmax=  1d0
* Rapidity
      ida  = 30000
      nba  = 60
**********************!for ZRR
      rmin= 1.57d0 -0.075
      rmax= 1.57d0 +0.075
**********************!for ZRR
* vvAL
      ids  = 40000
      nbs  = 100
      smin = 0d0
      smax = 1d0
*for ZRR
      nbs  = 60     !for ZRR
      nbs  = 30     !for ZRR
      smin= .775d0  !for ZRR, 0.792 at 200GeV
      smax= .815d0  !for ZRR, 0.792 at 200GeV
* linear scale
      idb  = 50000
      nbv  =   100
      vmin =     0
      vmax =     1
* new log scale
      idz  = 70000
      nbz  =   120
      zmin =  -5d0
      zmax =   7d0
* Flavour-cut table
ccc         CALL GLK_Book1(idf,' table (KFcode,vmax)  $', 400, 1d0, 401d0)
         CALL GLK_Book1(idf,' table (KFcode,vmax)  $', nKF*nCut, 1d0, nKF*nCut+1d0) !!!
* Totals and differences SEMIREALISTIC
      DO k=202,203
         CALL GLK_Book1(ids+k,     'Sig(v) CEEX, interf. ON, $',nbs,smin,smax)
         CALL GLK_Book1(ids+k+2000,'AFB(v) CEEX, interf. ON, $',nbs,smin,smax)
         CALL GLK_Book1(ids+k  +50,'Sig(v) CEEX, interf. OFF $',nbs,smin,smax)
         CALL GLK_Book1(ids+k+2050,'AFB(v) CEEX, interf. OFF $',nbs,smin,smax)
         CALL GLK_Book1(ids+k  +80,'Sig(v) CEEX, int. ON-OFF $',nbs,smin,smax)
         CALL GLK_Book1(ids+k+2080,'AFB(v) CEEX, int. ON-OFF $',nbs,smin,smax)
      ENDDO
      DO k=71,90
         CALL GLK_Book1(ids+k,'spec-kat aleph $',nbs,smin,smax)
      ENDDO
      DO k=71,74
         CALL GLK_Book1(ida+k,'rapidity   $',nba,rmin,rmax)
      ENDDO
* EEX all
      DO k=71,77
         CALL GLK_Book1(idb+k,' rho(v) EEX $',nbv,vmin,vmax)
         CALL GLK_Book1(idz+k,' rho(z) EEX $',nbz,zmin,zmax)
      ENDDO
* EEX beta contributions first order
      DO k=80,84
         CALL GLK_Book1(idb+k,' rho(v) EEX   $',nbv,vmin,vmax)
         CALL GLK_Book1(idz+k,' rho(z)  EEX   $',nbz,zmin,zmax)
      ENDDO
* EEX beta contributions second order
      DO k=90,97
         CALL GLK_Book1(idb+k,' rho(v)  EEX   $',nbv,vmin,vmax)
         CALL GLK_Book1(idz+k,' rho(z)  EEX   $',nbz,zmin,zmax)
      ENDDO
* EEX beta contributions third order
      DO k=100,111
         CALL GLK_Book1(idb+k,' rho(v)  EEX   $',nbv,vmin,vmax)
         CALL GLK_Book1(idz+k,' rho(z)  EEX   $',nbz,zmin,zmax)
      ENDDO
      DO k=71,74
         CALL GLK_Book1(idc +1000+ k  ,'CosTheta, v<1d-4  EEX$',nbc,cmin,cmax)
         CALL GLK_Book1(idc +2000+ k  ,'CosTheta, v<1d-2  EEX$',nbc,cmin,cmax)
         CALL GLK_Book1(idc +3000+ k  ,'CosTheta, v<0.1   EEX$',nbc,cmin,cmax)
         CALL GLK_Book1(idc +4000+ k  ,'CosTheta, v<0.9   EEX$',nbc,cmin,cmax)
         CALL GLK_Book1(idc +5000+ k  ,'CosTheta, v<0.9   EEX$',nbc,cmin,cmax)
      ENDDO
*//////////////////////////////////////////////////////////////////////////////
*//                 CEEX tests  dSigma/dCosTheta                             //
      DO k=201,203
         CALL GLK_Book1(idc +1000+ k ,'CosTheta, v<1d-4 CEEX$',nbc,cmin,cmax) !! CEEX
         CALL GLK_Book1(idc +2000+ k ,'CosTheta, v<1d-2 CEEX$',nbc,cmin,cmax) !! CEEX
         CALL GLK_Book1(idc +3000+ k ,'CosTheta, v<0.1  CEEX$',nbc,cmin,cmax) !! CEEX
         CALL GLK_Book1(idc +4000+ k ,'CosTheta, v<0.9  CEEX$',nbc,cmin,cmax) !! CEEX
         CALL GLK_Book1(idc +5000+ k ,'CosTheta, v<0.9  CEEX$',nbc,cmin,cmax) !! CEEX
      ENDDO
      DO k=251,253
         CALL GLK_Book1(idc +1000+ k ,'CosTheta, v<1d-4 CEEX intOFF$',nbc,cmin,cmax) !intOFF
         CALL GLK_Book1(idc +2000+ k ,'CosTheta, v<1d-2 CEEX intOFF$',nbc,cmin,cmax) !intOFF
         CALL GLK_Book1(idc +3000+ k ,'CosTheta, v<0.1  CEEX intOFF$',nbc,cmin,cmax) !intOFF
         CALL GLK_Book1(idc +4000+ k ,'CosTheta, v<0.9  CEEX intOFF$',nbc,cmin,cmax) !intOFF
         CALL GLK_Book1(idc +5000+ k ,'CosTheta, v<0.9  CEEX intOFF$',nbc,cmin,cmax) !intOFF
      ENDDO
*/////////////////////////////////////////////////////////////////////////////
*//                   EEX tests v-dsirtibutions                             //
      DO k=71,74
         CALL GLK_Book1(idb+k+2000,'rho(v) AFB EEX $',nbv,vmin,vmax)
         CALL GLK_Book1(idz+k+2000,'rho(v) AFB EEX $',nbz,zmin,zmax)
      ENDDO
*//////////////////////////////////////////////////////////////////////////////
*//                   CEEX tests v-distributions                             //
      DO k=201,203
         CALL GLK_Book1(idb+k,     'Sig(v) CEEX, interf. ON  $',nbv,vmin,vmax)
         CALL GLK_Book1(idb+k+2000,'AFB(v) CEEX, interf. ON  $',nbv,vmin,vmax)
         CALL GLK_Book1(idb+k  +50,'Sig(v) CEEX, interf. OFF $',nbv,vmin,vmax)
         CALL GLK_Book1(idb+k+2050,'AFB(v) CEEX, interf. OFF $',nbv,vmin,vmax)
         CALL GLK_Book1(idb+k  +80,'Sig(v) CEEX, int. ON-OFF $',nbv,vmin,vmax)
         CALL GLK_Book1(idb+k+2080,'AFB(v) CEEX, int. ON-OFF $',nbv,vmin,vmax)
         CALL GLK_Book1(idb+k  +70,'Sig(v) diff. CEEX-EEX int. OFF $',nbv,vmin,vmax)
         CALL GLK_Book1(idb+k+2070,'AFB(v) diff. CEEX-EEX int. OFF $',nbv,vmin,vmax)
      ENDDO
*//////////////////////////////////////////////////////////////////////////////
*//                  CEEX Differences between two orders                     //
      DO k=202,203
         CALL GLK_Book1(idb+k+1000,'intOFF, CEEX O(alf1-alf0) sig(vmax)$',nbv,vmin,vmax)
         CALL GLK_Book1(idb+k+3000,'intOFF, CEEX O(alf1-alf0) AFB      $',nbv,vmin,vmax)
      ENDDO
      DO k=252,253
         CALL GLK_Book1(idb+k+1000,'intOFF, CEEX O(alf1-alf0) sig(vmax)$',nbv,vmin,vmax)
         CALL GLK_Book1(idb+k+3000,'intOFF, CEEX O(alf1-alf0) AFB      $',nbv,vmin,vmax)
      ENDDO
*//////////////////////////////////////////////////////////////////////////////
*//                   CEEX ISR*FSR Interference effects                      //
      DO k=201,203
         CALL GLK_Book1(idz+k,     'Sig(v) CEEX, interf. ON, $',nbz,zmin,zmax)
         CALL GLK_Book1(idz+k+2000,'AFB(v) CEEX, interf. ON, $',nbz,zmin,zmax)
         CALL GLK_Book1(idz+k  +50,'Sig(v) CEEX, interf. OFF $',nbz,zmin,zmax)
         CALL GLK_Book1(idz+k+2050,'AFB(v) CEEX, interf. OFF $',nbz,zmin,zmax)
         CALL GLK_Book1(idz+k  +80,'Sig(v) CEEX, int. ON-OFF $',nbz,zmin,zmax)
         CALL GLK_Book1(idz+k+2080,'AFB(v) CEEX, int. ON-OFF $',nbz,zmin,zmax)
         CALL GLK_Book1(idz+k  +70,'AFB(v) diff. CEEX-OLD, int. OFF $',nbz,zmin,zmax)
         CALL GLK_Book1(idz+k+2070,'AFB(v) diff. CEEX-OLD, int. OFF $',nbz,zmin,zmax)
      ENDDO
*//////////////////////////////////////////////////////////////////////////////
* Totals ISR ONLY
      DO k=300,302
         CALL GLK_Book1(idb+k,'rho(v)      $',nbv,vmin,vmax)
         CALL GLK_Book1(idz+k,'rho(v)      $',nbz,zmin,zmax)
      ENDDO
      DO k=305,306
         CALL GLK_Book1(idb+k,'rho(v)      $',nbv,vmin,vmax)
         CALL GLK_Book1(idz+k,'rho(v)      $',nbz,zmin,zmax)
      ENDDO
*//////////////////////////////////////////////////////////////////////////////
      iev=0
      ELSEIF(Mode  .EQ.  0) THEN
*     ==========================
      iev=iev+1
      CMSene = pf1(4)+pf2(4)
***   CALL KK2f_GetKeyWgt(KeyWgt)
      CALL KK2f_GetWtAll( WtMain,WtCrud,WtSet)
      CALL KK2f_GetWtList(WtMain,WtList)
*//////////////////////////////////////
*//   Scattering angle    variables  //
*//////////////////////////////////////
      Flip  = 1d0               ! avoid undefined value
      costh = 0d0
      rapt  = 0d0
      vvAL  = 0d0
      IF( wtcrud .NE. 0d0 ) THEN
*-------------------------------
         cost1 = qf1(3) /SQRT( qf1(1)**2 +qf1(2)**2 +qf1(3)**2)
         cost2 = qf2(3) /SQRT( qf2(1)**2 +qf2(2)**2 +qf2(3)**2)
** definition of P.L. B219, 103 (1989)
**         costh = ( qf1(4)*cost1 -qf2(4)*cost2)/(qf1(4)+qf2(4))
* definition of P.R. D41, 1425 (1990)
         sint1 = DSQRT(DABS((1d0-cost1)*(1d0+cost1)))
         sint2 = DSQRT(DABS((1d0-cost2)*(1d0+cost2)))
         yy1 = sint2/(sint1+sint2)
         yy2 = sint1/(sint1+sint2)
         costh = yy1*cost1 - yy2*cost2
*-------------------------------
* LL formula for s'/s from angles according to ALEPH note 1996
         theta1= ACOS( cost1 )
         theta2= ACOS( cost2 )
         vvAL = (SIN(theta1)+SIN(theta2) -ABS(SIN(theta1+theta2)))
     $         /(SIN(theta1)+SIN(theta2) +ABS(SIN(theta1+theta2)))
         vvAL = 1d0 -vvAL
         rapt = ABS( log( (qf2(4)-qf2(3)) / (qf1(4)+qf1(3)) ) )
*****    write(*,'(a,9g20.10)') '////////////  vvAL,vv = ',vvAL,vvB,vvAL/vvB
*-------------------------------
* for comparison with Zfiter
*         costh = cost1
*-------------------------------
         Flip  = 1d0
         IF( costh .LT. 0d0 ) Flip = -1d0
*-------------------------------
      ENDIF
*//////////////////////////////////////
*//   Total photon Energy variables  //
*//////////////////////////////////////
      t10  = -10000.
      z10  = -10000.
      vvB   = 0d0
      IF( (nphot .NE. 0) .AND. (wtcrud .NE. 0d0) ) THEN
* vvB from outgoing charged particles
         s1=   (qf1(4)+qf2(4))**2 -(qf1(3)+qf2(3))**2
     $        -(qf1(2)+qf2(2))**2 -(qf1(1)+qf2(1))**2
         vvB = 1d0 -s1/cmsene**2
         t10= LOG10(vvB)
         z10= LOG10(vvB/(1-vvB))
      ENDIF
* vvP propagator abstract, from inside of MC generator, makes sense only for IFI OFF!!!
      vvPR = 0d0
      IF( (nphot .NE. 0) .AND. (wtcrud .NE. 0d0) ) THEN
         CALL BornV_GetVV(vvPR)
      ENDIF
*//////////////////////////////////////////////////////////////////////////////
*//                   Flavour distributions                                  //
         CALL HepEvt_GetKFfin(iKF)
         DO iCut=1,nCut
            ibin = nCut*(iKF-1) +iCut
            IF( vvB .LT. VmaxTab(iCut) ) CALL GLK_Fil1(idf, ibin*1.000001d0, WtMain) !
         ENDDO
*//////////////////////////////////////////////////////////////////////////////
*//                   CEEX tests v-distributions                             //
      DO k=201,203
         CALL GLK_Fil1(idb+k,      vvB, WtList(k)         )           !!! Interference ON
         CALL GLK_Fil1(idb+k+2000, vvB, WtList(k)*Flip    )           !!! Interference ON
         CALL GLK_Fil1(idb+k+  50, vvB, WtList(k+50)      )           !!! Interference OFF
         CALL GLK_Fil1(idb+k+2050, vvB, WtList(k+50)*Flip )           !!! Interference OFF
         CALL GLK_Fil1(idb+k+  80, vvB,(WtList(k)-WtList(k+50)) )     !!! delta = ON - OFF
         CALL GLK_Fil1(idb+k+2080, vvB,(WtList(k)-WtList(k+50))*Flip) !!! delta = ON - OFF
      ENDDO
*//////////////////////////////////////////////////////////////////////////////
*//                          CEEX-EEX                                        //
      DO k=1,3
         CALL GLK_Fil1(idb+ 270+k, vvB,(WtList(250+k)-WtList(k+70)))      !!! Sig: CEEXk-EEXk
***      CALL GLK_Fil1(idb+ 270+k, vvB,(WtList(250+k)-WtList(4+70)))      !!! Sig: CEEXk-EEX3ref
         CALL GLK_Fil1(idb+2270+k, vvB,(WtList(250+k)-WtList(k+70))*Flip) !!! AFB: CEEXk-EEXk
      ENDDO
*//////////////////////////////////////////////////////////////////////////////
*//                  CEEX Differences between two orders                     //
      DO k=202,203
         CALL GLK_Fil1(idb+k+1000, vvB,(WtList(k)-WtList(k-1)))       !!! Sig: CEEXi-CEEX(i-1)
         CALL GLK_Fil1(idb+k+3000, vvB,(WtList(k)-WtList(k-1))*Flip)  !!! AFB: CEEXi-CEEX(i-1)
      ENDDO
      DO k=252,253
         CALL GLK_Fil1(idb+k+1000, vvB,(WtList(k)-WtList(k-1)))       !!! Sig: CEEXi-CEEX(i-1)
         CALL GLK_Fil1(idb+k+3000, vvB,(WtList(k)-WtList(k-1))*Flip)  !!! AFB: CEEXi-CEEX(i-1)
      ENDDO
*//////////////////////////////////////////////////////////////////////////////
*//                   CEEX ISR*FSR Interference effects                      //
      DO k=201,203
         CALL GLK_Fil1(idz+k,      z10, WtList(k)         )           !!! Interference ON
         CALL GLK_Fil1(idz+k+2000, z10, WtList(k)*Flip    )           !!! Interference ON
         CALL GLK_Fil1(idz+k+  50, z10, WtList(k+50)      )           !!! Interference OFF
         CALL GLK_Fil1(idz+k+2050, z10, WtList(k+50)*Flip )           !!! Interference OFF
         CALL GLK_Fil1(idz+k+  80, z10,(WtList(k)-WtList(k+50)) )     !!! delta = ON - OFF
         CALL GLK_Fil1(idz+k+2080, z10,(WtList(k)-WtList(k+50))*Flip) !!! delta = ON - OFF
      ENDDO
*//////////////////////////////////////////////////////////////////////////////
*//                                 CEEX-EEX                                 //
      DO k=1,3
         CALL GLK_Fil1(idz+ 270+k, z10,(WtList(250+k)-WtList(k+70)))      !!! Sig: CEEX-EEX
         CALL GLK_Fil1(idz+2270+k, z10,(WtList(250+k)-WtList(k+70))*Flip) !!! AFB: CEEX-EEX
      ENDDO
*//////////////////////////////// s' study ////////////////////////////////////
      IF(vvB  .LT. vvBcut) THEN
         CALL GLK_Fil1(ids +90, vvAL,   WtList(203)-WtList(202)  ) !!!O(alf2-alf1)
         DO k=  202,203
            CALL GLK_Fil1(ids     +k, vvAL,   WtList(   k)  ) !!!Interference ON
            CALL GLK_Fil1(ids+2000+k, vvAL,   WtList(   k)*Flip  ) !!!Interference ON
            CALL GLK_Fil1(ids+  50+k, vvAL,   WtList(50+k)  ) !!!Interference OFF
            CALL GLK_Fil1(ids+2050+k, vvAL,   WtList(50+k)*Flip  ) !!!Interference OFF
            CALL GLK_Fil1(ids+  80+k, vvAL,  (WtList(   k)-WtList(50+k))  ) !!!delta = ON - OFF
            CALL GLK_Fil1(ids+2080+k, vvAL,  (WtList(   k)-WtList(50+k))*Flip ) !!!delta = ON - OFF
         ENDDO
      ENDIF
*//////////////////////////////////////////////////////////////////////////////
*//                       Angular distributions                              //
      DO k=71,74
         IF((vvAL.LT.vvALmax).AND.(vvB .LT.vvBcut)) THEN  !!! ALEPH energy cut
            CALL GLK_Fil1(idc+1000+k, costh,WtList(k))    !!! EEX
         ENDIF
         IF( vvB .LT. vvWGmax1 ) THEN                     !!! LEP2 WG
            CALL GLK_Fil1(idc+2000+k, costh,WtList(k))    !!! EEX
         ENDIF
         IF( vvB .LT. vvWGmax2 ) THEN                     !!! LEP2 WG
            CALL GLK_Fil1(idc+3000+k, costh,WtList(k))    !!! EEX
         ENDIF
         IF( vvB .LT. vvWGmax9 ) THEN                     !!! private
            CALL GLK_Fil1(idc+4000+k, costh,WtList(k))    !!! EEX
         ENDIF
         IF( vvPR .LT. vvWGmax2 ) THEN                    !!! s'-PRopagator
            CALL GLK_Fil1(idc+5000+k, costh,WtList(k))    !!! CEEX
         ENDIF
      ENDDO
      DO k=201,203
         IF((vvAL.LT.vvALmax).AND.(vvB .LT.vvBcut)) THEN  !!! ALEPH energy cut(= 0.2775d0 )
            CALL GLK_Fil1(idc+1000+k, costh,WtList(k))    !!! CEEX
            CALL GLK_Fil1(idc+1050+k, costh,WtList(k+50)) !!! CEEX intOFF
         ENDIF
         IF( vvB .LT. vvWGmax1 ) THEN                     !!! LEP2 WG (= 0.1900d0 )
            CALL GLK_Fil1(idc+2000+k, costh,WtList(k))    !!! CEEX
            CALL GLK_Fil1(idc+2050+k, costh,WtList(k+50)) !!! CEEX intOFF
         ENDIF
         IF( vvB .LT. vvWGmax2 ) THEN                     !!! LEP2 WG (= 0.2775d0 )
            CALL GLK_Fil1(idc+3000+k, costh,WtList(k))    !!! CEEX
            CALL GLK_Fil1(idc+3050+k, costh,WtList(k+50)) !!! CEEX intOFF
         ENDIF
         IF( vvB .LT. vvWGmax9 ) THEN                     !!! private (= 0.9000d0 )
            CALL GLK_Fil1(idc+4000+k, costh,WtList(k))    !!! CEEX
            CALL GLK_Fil1(idc+4050+k, costh,WtList(k+50)) !!! CEEX intOFF
         ENDIF
         IF( vvPR .LT. vvWGmax2 ) THEN                    !!! s'-PRopagator
            CALL GLK_Fil1(idc+5000+k, costh,WtList(k))    !!! CEEX
            CALL GLK_Fil1(idc+5050+k, costh,WtList(k+50)) !!! CEEX intOFF <---NONSENSE!!!!
         ENDIF
      ENDDO
********************
* EEX Total xsections
      DO k=71,74
         CALL GLK_Fil1(idb+k,       vvB,WtList(k))
         CALL GLK_Fil1(idb+k+2000,  vvB,WtList(k)*Flip)
         CALL GLK_Fil1(idz+k,      z10,WtList(k))
         CALL GLK_Fil1(idz+k+2000, z10,WtList(k)*Flip)
         CALL GLK_Fil1(ida+k,     rapt,WtList(k))
         IF(vvB  .LT. vvBcut) THEN
            CALL GLK_Fil1(ids+k,     vvAL,WtList(k))
         ENDIF
      ENDDO
* EEX Differences
      DO k=72,74
         CALL GLK_Fil1(idb+k+3,  vvB,(WtList(k)-WtList(k-1)))
         CALL GLK_Fil1(idz+k+3, z10,(WtList(k)-WtList(k-1)))
         CALL GLK_Fil1(ida+k+3,rapt,(WtList(k)-WtList(k-1)))
         IF(vvB  .LT. vvBcut) THEN
            CALL GLK_Fil1(ids+k+3,vvAL,(WtList(k)-WtList(k-1)))
         ENDIF
      ENDDO
* EEX beta contributions first order
      DO k=80,84
         CALL GLK_Fil1(idb+k, vvB ,WtList(k))
         CALL GLK_Fil1(idz+k, z10,WtList(k)) ! new
      ENDDO
* EEX beta contributions second order
      DO k=90,97
         CALL GLK_Fil1(idb+k, vvB ,WtList(k))
         CALL GLK_Fil1(idz+k, z10,WtList(k)) ! new
      ENDDO
* beta contributions third order
      DO k=100,111
         CALL GLK_Fil1(idb+k, vvB ,WtList(k))
         CALL GLK_Fil1(idz+k, z10,WtList(k)) ! new
      ENDDO
***************************************
**  INITIAL ONLY  (nphox, xf1, xf2)  **
***************************************
      t10  = -10000.
      vvB   = 0d0
* Should we use wtcru1 here????
      IF( (nphox .NE. 0) .AND. (wtcrud .NE. 0d0) ) THEN
         DO k=1,4
            xxf(k) = pf1(k)+pf2(k)-xphum(k)
         ENDDO
         s1= xxf(4)**2-xxf(3)**2-xxf(2)**2-xxf(1)**2
         vvB = 1d0 -s1/cmsene**2
         t10= log10(vvB)
      ENDIF
* Totals 300-302
      DO k=1,3
         CALL GLK_Fil1(idb+ 299+k, vvB,WtList(k))
         CALL GLK_Fil1(idz+ 299+k,z10,WtList(k))
      ENDDO
* Differences 305-306
      DO k=2,3
         CALL GLK_Fil1(idb+ 303+k, vvB,(WtList(k)-WtList(k-1)))
         CALL GLK_Fil1(idz+ 303+k,z10,(WtList(k)-WtList(k-1)))
      ENDDO
***************************************
*------------------------------------------------------
* Virtual Variation of MZ
      CALL BornV_GetMZ(aMZ)
* redefine MZ
      bMZ = aMZ + 0.010d0
      CALL BornV_SetMZ(bMZ)
* recalculate weights
      IF(WtCrud .NE. 0d0) CALL QED3_Make
      CALL QED3_GetWtSet(  WtBest,WtSet )
      DO k=1,1000
         WtList2(k)=WtSet(k)*WtCrud
      ENDDO
      CALL BornV_SetMZ(aMZ)     !!!!!<--RESET TO ORIGINAL VALUE!!!!
*
      IF(vvB  .LT. vvBcut) THEN
         CALL GLK_Fil1(ids +78, vvAL,(WtList2(72) -WtList(72))) !!!O(alf1)
      ENDIF
      CALL GLK_Fil1(ida +78, rapt,(WtList2(72) -WtList(72)))   !!!O(alf1)
*------------------------------------------------------
      ELSE
*     ===========
         WRITE(nout,*) '==================================='
         WRITE(nout,*) '============ robol4 ==============='
         WRITE(nout,*) '==================================='
*----------------------
******         CALL GLK_PrintAll
* CEEX related histograms
         DO k=202,203
            CALL GLK_Print(ids+k     )
            CALL GLK_Print(ids+k+2000)
            CALL GLK_Print(ids+k  +50)
            CALL GLK_Print(ids+k+2050)
            CALL GLK_Print(ids+k  +80)
            CALL GLK_Print(ids+k+2080)
         ENDDO
         DO k=201,203
            DO i=1,5 
               CALL GLK_Print(idc +1000*i    + k) !CosTheta distr CEEX
               CALL GLK_Print(idc +1000*i+50 + k) !CosTheta distr CEEX IntOff
            ENDDO
         ENDDO
         DO k=201,203
            CALL GLK_Print(idb+k     )
            CALL GLK_Print(idb+k+2000)
            CALL GLK_Print(idb+k  +50)
            CALL GLK_Print(idb+k+2050)
            CALL GLK_Print(idb+k  +80)
            CALL GLK_Print(idb+k+2080)
            CALL GLK_Print(idb+k  +70)
            CALL GLK_Print(idb+k+2070)
         ENDDO
         DO k=202,203
            CALL GLK_Print(idb+k+1000)
            CALL GLK_Print(idb+k+3000)
         ENDDO
         DO k=252,253
            CALL GLK_Print(idb+k+1000)
            CALL GLK_Print(idb+k+3000)
         ENDDO
         CALL GLK_Print(idf)
         WRITE(nout,'(3a5,5a20)') 'ibin','iCut','KF' ,'bin', 'err'!
         DO iKF=1,5
            DO iCut=1,nCut
               ibin = nCut*(iKF-1) +iCut
               WRITE(nout,'(3i5,5f20.9)') ibin,iCut,iKF ,GLK_hi( idf,ibin) ,GLK_hie( idf,ibin) !
            ENDDO
         ENDDO
*----------------------
      ENDIF
      END                       ! Robol4


      SUBROUTINE robol5(Mode,xpar)
*/////////////////////////////////////////////////////////////////////////////////////////
*//                                                                                     //
*//  study on experimental cut-offs                                                     //
*//                                                                                     //
*/////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      REAL*8   xpar(*)
      COMMON / inout  / ninp,nout
      COMMON / momset / pf1(4),pf2(4),qf1(4),qf2(4),sphum(4),sphot(100,4),nphot,KFfin
      REAL*8  WtSet(1000),WtList(1000)
      SAVE

      IF(Mode  .EQ.  -1) THEN
*     ===================
      CALL avryg( -1,20,' mupair all cuts $',dum1,0d0,1d0)
      CALL avryg( -1,21,' minimum energy   muons $',dum1,0d0,1d0)
      CALL avryg( -1,22,' minimum costheta muons $',dum1,0d0,1d0)
      CALL avryg( -1,23,' acollinearity $',dum1,0d0,1d0)
      CALL GLK_Book1(  70,'log dn/dv     -all           $',40, 0d0,1d0)
      CALL GLK_Book1(  71,'log dn/dv    acoll cut       $',40, 0d0,1d0)
      CALL GLK_idopt(70,'logy')
      CALL GLK_idopt(71,'logy')
**     CALL GLK_Book2(200,' acol/v-var.  $',40,0d0 ,40d0, 50,0d0,1.00d0)

      ELSEIF(Mode  .EQ.  0) THEN
*     ======================
      CMSene = pf1(4)+pf2(4)
      CALL KK2f_GetWtAll(WtMain,WtCrud,WtSet)
      DO k=1,1000
         WtList(k)=WtSet(k)*WtCrud
      ENDDO
      wt=WtMain

* define mupair cut-off
      CALL defmup(w1,w2,w3)
      ww=wt*w1*w2*w3
      CALL avryg(  0, 20,' ',ww,0d0,0d0)
      CALL avryg(  0, 21,' ',w1,0d0,0d0)
      CALL avryg(  0, 22,' ',w2,0d0,0d0)
      CALL avryg(  0, 23,' ',w3,0d0,0d0)
      s1=   (qf1(4)+qf2(4))**2 -(qf1(3)+qf2(3))**2
     $     -(qf1(2)+qf2(2))**2 -(qf1(1)+qf2(1))**2
      vv = 1d0 -s1/cmsene**2
      CALL GLK_Fil1(  70, vv, wt)
      CALL GLK_Fil1(  71, vv, w3)
      acol=0d0
      IF(nphot.ne.0) acol=
     $ acos(-(qf1(1)*qf2(1)+qf1(2)*qf2(2)+qf1(3)*qf2(3))
     $          /sqrt((qf1(1)**2 +qf1(2)**2 +qf1(3)**2)*
     $                (qf2(1)**2 +qf2(2)**2 +qf2(3)**2)))*180/3.141593
*c      CALL hfill( 200, acol,vv,  wt)
      ELSE
*     ====
      WRITE(nout,*) '==================================='
      WRITE(nout,*) '============ robol5 ==============='
      WRITE(nout,*) '==================================='
      CALL avryg(110, 20,' ', xav,xer,xnev)
      CALL avryg(110, 21,' ', xav,xer,xnev)
      CALL avryg(110, 22,' ', xav,xer,xnev)
      CALL avryg(110, 23,' ', xav,xer,xnev)
*----------------------
***         CALL GLK_PrintAll
*----------------------
      ENDIF
      END
      SUBROUTINE defmup(w1,w2,w3)
*     ***************************
*-------------------------------------------------------------------
*  MARKII-like cuts from bennie:
*  visible muon is
*        |cos(theta-mu)|     < 0.80;
*        e-mu       > 2.0 gev;
*  visible gammas are
*        |cos(theta-gamma)|  < 0.95;
*        e-gamma             > 0.2 gev;
*  event is accepted IF there are two visible muons and
*        e-visible           > 0.1*sqrt(s),
*  where e-visible includes the visible mu-mu-bar and gamma energy.
*-------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON / momset / pf1(4),pf2(4),qf1(4),qf2(4),sphum(4),sphot(100,4),nphot,KFfin
      SAVE

      CMSene = pf1(4)+pf2(4)
      w1=1d0
      w2=1d0
      w3=1d0
      cmulim=0.80d0
      emulim=5.0d0
      cgalim=0.95d0
      egalim=0.2d0
      xvslim=0.3d0
      acolim=20.d0
* two visible muons
      IF(qf1(4)  .LT.  emulim) w1=0d0
      IF(qf2(4)  .LT.  emulim) w1=0d0
      cost1= qf1(3)/sqrt(qf1(1)**2+qf1(2)**2+qf1(3)**2)
      IF(abs(cost1)  .GT.  cmulim) w2=0d0
      cost2= qf2(3)/sqrt(qf2(1)**2+qf2(2)**2+qf2(3)**2)
      IF(abs(cost2)  .GT.  cmulim) w2=0d0
* minimum visible energy
*     evis= qf1(4)+qf2(4)
*     DO 10 i=1,nphot
*     costg= sphot(i,3)/sqrt(sphot(i,1)**2+sphot(i,2)**2+sphot(i,3)**2)
*     eg= sphot(i,4)
*     IF(abs(costg)  .LT.  cgalim  .AND.  eg  .GT.  egalim)  evis=evis+eg
*  10 CONTINUE
*     IF(evis  .LT.  xvslim*cmsene) w3=0d0
* acollinearity
      acol=0d0
      IF(nphot.ne.0) THEN
         dq12= qf1(1)*qf2(1) +qf1(2)*qf2(2) +qf1(3)*qf2(3)
         qq1 = qf1(1)**2 +qf1(2)**2 +qf1(3)**2
         qq2 = qf2(1)**2 +qf2(2)**2 +qf2(3)**2
         acol= acos(-dq12/sqrt(qq1*qq2))*180/3.141593
      ENDIF
      IF(acol  .GT.  acolim) w3=0d0
      END

      SUBROUTINE avryg(Mode,id,title, x,xmin,xmax)
*     ********************************************
* Utility program for calculating up to idmax averages.
* IF(Mode  .EQ.  -1) THEN
*          initialization of entry id,
*          title is up to 32 CHARACTER title, last CHARACTER
*                 must necessarily be a $
*          x ignored, xmin and xmax are limits on x.
* ELSEIF(Mode  .EQ.  0) THEN
*          summing up x value for a given entry id,
*          (x is summed up in average even IF outside xmin,xmax),
*          title, xmin, xmax ignored.
* ELSEIF(Mode  .GE.  1) THEN
*          x    = average <x>
*          xmin = absolute error
*          xmax = number of all sampled events
*          printout depending on the value of the Mode
*               Mode=1       -none
*               Mode=2       -title  only
*               Mode=10,100  -partial
*               Mode=110     -maximal
* ENDIF
*     ************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(idmax=80,lentit=32)
      CHARACTER*32 title
      COMMON / cvryg  / averx,errelx,nevtot,nevacc,nevund,nevove,nevzer
      COMMON / inout  / ninp,nout
      SAVE / cvryg  /,/ inout  /
      CHARACTER*32 mtitle(idmax)
      INTEGER   ntot(idmax),   nacc(idmax)
      INTEGER   nund(idmax),   nove(idmax),   nzer(idmax)
      REAL*8     swt(idmax),  sswt(idmax)
      REAL*8   bbmin(idmax), bbmax(idmax)
      REAL*8   xxmin(idmax), xxmax(idmax)
      LOGICAL met
      DATA  ntot    /idmax*-1/    swt   /idmax*0d0/  sswt /idmax*0d0/
      DATA xxmin /idmax*-1d20/  xxmax /idmax*-1d20/
      DATA bbmin /idmax*-1d20/  bbmax /idmax*-1d20/
      SAVE mtitle,ntot,nacc,nund,nove,nzer,swt,sswt
      SAVE bbmin,bbmax,xxmin,xxmax

      IF(id  .LE.  0  .OR.  id  .GT.  idmax) THEN
           WRITE(nout,*) ' =====avryg: wrong id'
           STOP
      ENDIF
* initialization
      IF(Mode  .EQ.  -1) THEN
           ntot(id)=0
           nacc(id)=0
           nund(id)=0
           nove(id)=0
           nzer(id)=0
           swt(id)   =0d0
           sswt(id)  =0d0
           xxmax(id)  = -1d+30
           xxmin(id)  = +1d+30
           bbmax(id)  = xmax
           bbmin(id)  = xmin
*        storing the title
           met = .false.
           DO 10 i=1,lentit
           IF( title(i:i)  .EQ.  '$'   .OR.   met )   THEN
             mtitle(id)(i:i)=  '='
             met=.true.
           ELSE
             mtitle(id)(i:i)=title(i:i)
           ENDIF
  10       CONTINUE
* suming up the contributions
      ELSEIF(Mode  .EQ.  0) THEN
           IF(ntot(id)  .LT.  0) THEN
              WRITE(nout,*) ' ==== warning from avryg: '
              WRITE(nout,*) ' lack of initialization, id=',id
           ENDIF
           ntot(id)  =ntot(id) +1
           swt(id)   =swt(id)  +x
           sswt(id)  =sswt(id) +x*x
           xxmax(id)= max(xxmax(id), x)
           xxmin(id)= min(xxmin(id), x)
           IF( x  .EQ.  0d0)   nzer(id)=nzer(id)+1
           IF(     x  .LT.  bbmin(id))  THEN
             nund(id) = nund(id)+1
           ELSEIF( x  .GT.  bbmax(id))  THEN
             nove(id) = nove(id)+1
           ELSE
             nacc(id) = nacc(id)+1
           ENDIF
* calculating averages - final report
      ELSEIF( Mode  .GE.  1) THEN
           IF(ntot(id)  .LE.  0  .OR.  swt(id)  .EQ.  0d0)  THEN
              averx  =0d0
              errelx =0d0
           ELSE
              averx =swt(id)/float(ntot(id))
              varia =sswt(id)/ntot(id)-(swt(id)/ntot(id))**2
              errox =sqrt(abs(varia)/ntot(id))
           ENDIF
           nevtot=ntot(id)
           nevacc=nacc(id)
           nevund=nund(id)
           nevzer=nzer(id)
           nevove=nove(id)
           x    =averx
           xmin =errox
           xmax =nevtot
           k1= MOD((Mode/10 ),10)
           k2= MOD((Mode/100),10)
 1000 FORMAT(1x,a,i2,a,a,a)
           IF(Mode  .GT.  1) WRITE(nout,1000)
     $      '========[', id, ']============' ,mtitle(id),
     $      '======================'
 1003 FORMAT(1x,a,4a17,a / 1x,a,4e17.9,a)
           IF(k1  .EQ.  1) WRITE(nout,1003)
     $      '====','< x >','abs_err', 'min(x)'  , 'max(x)'  ,'  ===='
     $     ,'====', averx ,  errox  ,  xxmin(id),  xxmax(id),'  ===='
 1004 FORMAT(1x,a,5(a,i6),a)
           IF(k2  .EQ.  1) WRITE(nout,1004)
     $     '====   ','  ntot=',ntot(id),'  nacc=',nacc(id)
     $              ,'  nund=',nund(id),'  nove=',nove(id)
     $              ,'  nzer=',nzer(id),'  ===='
      ELSE
           WRITE(nout,*) ' =====avryg:  wrong Mode'
           STOP
      ENDIF
      END
*=======================================================================
*=======================================================================
*     subprograms for control printouts of events
*=======================================================================
*=======================================================================

      SUBROUTINE  momprt(txt,nout,iev,ie1,ie2,pf1,pf2,qf1,qf2,nphot,sphot,KFfin)
*     ***********************************************************
* Prints out four momenta of Beams and Final state particles,
* and the serial number of event iev on unit nout
*     **********************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION  pf1(4),pf2(4),qf1(4),qf2(4),sphum(4),sphot(100,4)
      CHARACTER*8 txt
      DIMENSION sum(4)

      IF( (iev .GE. ie1) .AND. (iev .LE. ie2) ) THEN
         WRITE(nout,*) 
     $        '=========== ',txt,' ======================>',iev
*
         amf1 = pf1(4)**2-pf1(3)**2-pf1(2)**2-pf1(1)**2
         amf1 = sqrt(abs(amf1))
         amf2 = pf2(4)**2-pf2(3)**2-pf2(2)**2-pf2(1)**2
         amf2 = sqrt(abs(amf2))
         WRITE(nout,3100) 'pf1',(  pf1(  k),k=1,4),amf1
         WRITE(nout,3100) 'pf2',(  pf2(  k),k=1,4),amf2
*
         amf1 = qf1(4)**2-qf1(3)**2-qf1(2)**2-qf1(1)**2
         amf1 = sqrt(abs(amf1))
         amf2 = qf2(4)**2-qf2(3)**2-qf2(2)**2-qf2(1)**2
         amf2 = sqrt(abs(amf2))
         WRITE(nout,3100) 'qf1',(  qf1(  k),k=1,4),amf1,KFfin
         WRITE(nout,3100) 'qf2',(  qf2(  k),k=1,4),amf2,KFfin
*
         DO i=1,nphot
            amph = sphot(i,4)**2-sphot(i,3)**2
     $            -sphot(i,2)**2-sphot(i,1)**2
            amph = sqrt(abs(amph))
            WRITE(nout,3100) 'pho',(sphot(i,k),k=1,4),amph
         ENDDO
         DO k=1,4
            sum(k)=qf1(k)+qf2(k)
         ENDDO
         DO i=1,nphot
            DO k=1,4
               sum(k)=sum(k)+sphot(i,k)
            ENDDO
         ENDDO
         ams = sum(4)**2-sum(3)**2-sum(2)**2-sum(1)**2
         ams = sqrt(abs(ams))
         WRITE(nout,3100) 'sum',(  sum(  k),k=1,4),ams
      ENDIF

 3100 FORMAT(1x,a3,1x,5f20.14,i5)
      END

      SUBROUTINE dumpri(txt,nout,iev,ie1,ie2,qf1,qf2,nphot,sphot)
*     ***********************************************************
* Prints out four momenta of FINAL state
* and the serial number of event iev on unit nout
*     **********************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION  qf1(4),qf2(4),sphum(4),sphot(100,4)
      CHARACTER*8 txt
      DIMENSION sum(4)

      IF( (iev .GE. ie1) .AND. (iev .LE. ie2) ) THEN
         WRITE(nout,*) 
     $        '=========== ',txt,' ======================>',iev
         amf1 = qf1(4)**2-qf1(3)**2-qf1(2)**2-qf1(1)**2
         amf1 = sqrt(abs(amf1))
         amf2 = qf2(4)**2-qf2(3)**2-qf2(2)**2-qf2(1)**2
         amf2 = sqrt(abs(amf2))
         WRITE(nout,3100) 'qf1',(  qf1(  k),k=1,4),amf1
         WRITE(nout,3100) 'qf2',(  qf2(  k),k=1,4),amf2
         DO i=1,nphot
            amph = sphot(i,4)**2-sphot(i,3)**2
     $            -sphot(i,2)**2-sphot(i,1)**2
            amph = sqrt(abs(amph))
            WRITE(nout,3100) 'pho',(sphot(i,k),k=1,4),amph
         ENDDO
         DO k=1,4
            sum(k)=qf1(k)+qf2(k)
         ENDDO
         DO i=1,nphot
            DO k=1,4
               sum(k)=sum(k)+sphot(i,k)
            ENDDO
         ENDDO
         ams = sum(4)**2-sum(3)**2-sum(2)**2-sum(1)**2
         ams = sqrt(abs(ams))
         WRITE(nout,3100) 'sum',(  sum(  k),k=1,4),ams
      ENDIF

 3100 FORMAT(1x,a3,1x,5f20.14)
      END

