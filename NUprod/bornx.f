*----------------------------------------------------------------
*     make ew-mu2
*     make ew-numu
*     make ew-nuel
*     make bornx
*----------------------------------------------------------------
* Two Born dSigma/dCosTheta are available:
*  (1)   KKsem_BornV(svar,costhe) in KKsem/KKsem.f
*  (2)   BornV_Dizet(Mode,KFi,KFf,svar,CosThe,eps1,eps2,ta,tb) in bornv/BorbV.f
* EW corrections are in (2), search for  "DelW" for example
* It is KeyFob which manages which Born dSigma/dCosTheta is used
* KeyElw is also involved, as usual.
*----------------------------------------------------------------
      PROGRAM MAIN
*     ************
      IMPLICIT NONE
      INTEGER   ninp,nout
      CHARACTER*60      Tesnam, Dname
*------------------------------------------------------------------------------
      INTEGER    imax
      PARAMETER (imax=10000)
      DOUBLE PRECISION     xpar(imax)
      DOUBLE PRECISION     ymin,ymax
      DOUBLE PRECISION     CMSene
      DOUBLE PRECISION  BornPb,Born,Born2,BornNb,BorFor,BorBac,BorAll,AFB
      DOUBLE PRECISION  Cmax
      INTEGER    nene
      PARAMETER (nene=12)
      DOUBLE PRECISION  CMSlist(nene)
      DATA              CMSlist /
cc     $     19.87d0, 26.73d0, 32.73d0, 37.80d0, 42.26d0, 59.76d0,
     $     10.00d0, 20.00d0, 30.00d0, 40.00d0, 50.00d0, 60.00d0,
     $     91.19d0,   100d0,   140d0,   189d0,   200d0,   206d0/
      INTEGER           i,k,j,jmax
      INTEGER           KeyA,KeyF,KeyB,KeyFoB,KeyELW
      DOUBLE PRECISION  Sw2,MZ,GammZ,MW,Gmu,AlfInv
*--------------------------------------------------------------------
      ninp=  5
      nout= 16
      Tesnam    = 'bornx'
      OPEN( nout, file='output-'//Tesnam)
      CALL GLK_SetNout(nout)

      Dname  = './bornx.input'
      CALL KK2f_ReaDataX('../../.KK2f_defaults', 1,imax,xpar)  ! Read data, the same as in MC run
      CALL KK2f_ReaDataX(                 Dname, 0,imax,xpar)  ! Read user input
      CALL KK2f_Initialize( xpar)                    ! Initialize generator with the production data
      CALL KKsem_Initialize(xpar)                    ! Initialize semianalytical package
*=================================================================================================
      CMSene = xpar(1)
      KeyELW = xpar(12)
      Cmax= 1d0
      Cmax= 0.9999999d0
      Cmax= 0.99999d0
      CALL KKsem_SetCmax(Cmax)
      KeyFoB= -100
      CALL KKsem_SetKeyFoB(KeyFoB)
***      CALL KKsem_GetBorn(Born2)  !
      IF(KeyELW .EQ. 0) THEN    ! for KeyELW=0
         KeyF=   1              ! Forward
         KeyB=  -1              ! Backward
         KeyA=  10              ! All
      ELSE                      ! for KeyELW=1
         KeyF=   1              ! Forward
         KeyB=  -1              ! Backward
         KeyA=  10              ! All
      ENDIF
*------------------------------
*     EW paratemetrs taken from BornV
      CALL  BornV_GetMZ(MZ)
      CALL  BornV_GetGammZ(GammZ)
      CALL  BornV_GetSwsq(Sw2)
      CALL  BornV_GetMW(MW)
      CALL  BornV_GetGmu(Gmu)
      CALL  BornV_GetAlfInv(AlfInv)
      WRITE(   *,'(a)') '|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||'
      WRITE(   *,'(3(a,f20.14))') ' MZ = ',MZ, '  GammZ= ',GammZ
      WRITE(   *,'(3(a,f20.14))') ' MW = ',MW, '    Sw2= ',Sw2
      WRITE(   *,'(3(a,f20.14))') ' Gmu= ',Gmu,' AlfInv= ',AlfInv
      WRITE(   *,'(a)') '-------------------------------------------------------------'
      WRITE(   *,'(a10,4a20)') 'CMSene','Gauss BornNb','AFB','3/8*Dsig(0)'
      WRITE(   *,'(a)') '-------------------------------------------------------------'
      DO j=1,nene
         CMSene=CMSlist(j)
         CALL KKsem_SetKeyFoB(KeyF)
         CALL KKsem_MakeBorn(CMSene**2,BorFor)
         CALL KKsem_SetKeyFoB(KeyB)
         CALL KKsem_MakeBorn(CMSene**2,BorBac)
         CALL KKsem_SetKeyFoB(KeyA)
         CALL KKsem_MakeBorn(CMSene**2,BorAll)
         BornNb=BorFor+BorBac
         AFB = (BorFor-BorBac)/BornNb
         WRITE(   *,'(f10.2,4f20.14)') CMSene,BornNb,AFB,BorAll
         WRITE(nout,'(f10.2,4f20.14)') CMSene,BornNb,AFB,BorAll
      ENDDO
      WRITE(   *,'(a)') '-------------------------------------------------------------' !
      CLOSE(nout)
      END
