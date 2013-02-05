*///////////////////////////////////////////////////////////////////////////////
*//                   main program for MC production                          //
*///////////////////////////////////////////////////////////////////////////////
*//
*//   make Tau-start
*//
*///////////////////////////////////////////////////////////////////////////////
      PROGRAM main
*     **********************************
      IMPLICIT NONE
* input/output files
      INTEGER               ninp,nout
      COMMON / c_MainPro  / ninp,nout
      CHARACTER*4 semaph
      INTEGER ntot2n, ninph, ninp2, ninp3, ijklin, ntotin
      SAVE

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
* *******************************************
      CALL splot
* ********************************************
      END


      SUBROUTINE splot
*     ****************
      IMPLICIT NONE
*     ***********************************
      INTEGER               ninp,nout
      COMMON / c_MainPro  / ninp,nout
      SAVE
      INTEGER    imax
      PARAMETER (imax = 10000)               ! length ox xpar
      REAL*8     xpar(imax)
      CHARACTER*4   semaph,chdum
      CHARACTER*80  DiskFile
      INTEGER       kat1,kat2,kat3,kat4
      INTEGER       igroup, ngroup, nevt, loop, iev
      SAVE
*======================================================================
      WRITE(nout,*) '   '
      WRITE(nout,*) '=============================================='
      WRITE(nout,*) '==========*********************==============='
      WRITE(nout,*) '==========***    MainPro    ***==============='
      WRITE(nout,*) '==========*********************==============='
      WRITE(nout,*) '=============================================='
      WRITE(nout,*) '   '
*======================================================================
* Read data for main program
      OPEN( ninp,FILE='./pro.input')
      READ( ninp,'(4i2)') kat1,kat2,kat3,kat4
      WRITE(nout,'(4a6/4i6)')
     $ 'kat1','kat2','kat3','kat4',
     $  kat1 , kat2 , kat3 , kat4 
      READ(ninp,'(i10)') nevt
      CLOSE(ninp)
      WRITE(   6,*)   nevt,' requested events '
      WRITE(nout,*)   nevt,' requested events '
*
      CALL KK2f_ReaDataX('../../.KK2f_defaults', 1,imax,xpar)  ! reading general defaults
      CALL KK2f_ReaDataX(         './pro.input', 0,imax,xpar)  ! reading user input
      CALL KK2f_Initialize(xpar)
*
      IF(kat1 .EQ. 1) CALL Robol1(-1,xpar)
      IF(kat2 .EQ. 1) CALL Robol2(-1)
*-------------------------------------------------------!
*                 main MC loop                          !
*-------------------------------------------------------!
      ngroup = 5000
      iev=0
      DO loop=1,10000000
        DO igroup =1,ngroup
          iev=iev+1
          IF(MOD(iev, ngroup) .EQ. 1) WRITE( 6,*)  'iev= ',iev
*         **************
          CALL KK2f_Make
*         **************
*   Control printouts
*          CALL momprt(' YFSPRO ', 6,iev,1,10,pf1,pf2,qf1,qf2,nphot,sphot,KFfin)
*          CALL dumpri('*momini*', 6,iev,1,10,xf1,xf2,nphox,xphot)
          IF(iev .LE. 10) THEN
             CALL lugive("MSTU(11)=16")
             CALL lulist(1)
          ENDIF
          IF(iev .LE. 10) THEN
             CALL lugive("MSTU(11)=6")
             CALL lulist(1)
          ENDIF
*         ============================================
*         histograming
          IF(kat1 .EQ. 1) CALL Robol1( 0,xpar)
          IF(kat2 .EQ. 1) CALL Robol2( 0)
*         ============================================
*         check on requested no. of events
          IF(iev  .EQ.  nevt)     GOTO 300
        ENDDO
*       check on semaphore flag
        CALL givsem(semaph)
        IF(semaph  .EQ.  'STOP') GOTO 300
*       dump partial results on the disk after every ngroup
        CALL  dumpeh(iev)
      ENDDO
 300  CONTINUE
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
      WRITE(6,*) ' generation finished '
      CALL KK2f_Finalize
*------------------------------------------------
      IF(kat1 .EQ. 1) CALL Robol1(1,xpar)
      IF(kat2 .EQ. 1) CALL Robol2(1)
*
      CALL  dumpeh(iev)
*     +++++++++++++++++
      END

      SUBROUTINE givsem(semaph)
*     ************************
      IMPLICIT NONE
      CHARACTER*4 semaph
      INTEGER     ninp2
* ------------------------------------------------------
* READ semaphore flag
* ------------------------------------------------------
      ninp2=2
      OPEN(ninp2,FILE='./semaphore')
      READ(ninp2,'(a4)') semaph
      CLOSE(ninp2)
      END

      SUBROUTINE dumpeh(nev)
*     ************************
      IMPLICIT NONE
      INTEGER    ijklin, ntotin, ntot2n, ninp2, nev, nouth, icy
* ------------------------------------------------------
* WRITE histos on the disk
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
*//   Example of histogramming, with acces to output through getters                   //
*//                                                                                    //
*//                                                                                    //
*////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      REAL*8 xpar(*)
*
      INTEGER mode
      REAL*8  pf1(4),pf2(4),qf1(4),qf2(4)
      INTEGER NphAll,NphIni,NphFin
      REAL*8  PhoAll(100,4),PhoIni(100,4),PhoFin(100,4)
      REAL*8  WtSet(1000)
      INTEGER  i, idv, nbv, KFfin
      REAL*8   enetot, vv, cmsene, xkf, eneini, enefin, vmax, vmin, ss, ss2, wtmain, wtcrud
      SAVE
*=======================================================================
      IF(mode .EQ. -1 ) THEN
*     book histograms
         nbv = 400
         vmin= 0
         vmax= 1
         idv  = 50000
         CALL GLK_Book1(idv+1,'d_sigma/dv(v)  $',nbv,vmin,vmax)
         CALL GLK_Book1(idv+2,'Ene total ALL  $',nbv,vmin,vmax)
         CALL GLK_Book1(idv+3,'Ene total ISR  $',nbv,vmin,vmax)
         CALL GLK_Book1(idv+4,'Ene total FSR  $',nbv,vmin,vmax)
         CALL GLK_Book1(idv+7,'NphAll         $',16,0d0,16d0)
         CALL GLK_Book1(idv+9,'KFcode         $',16,1d0,17d0)
*=======================================================================
      ELSEIF(mode .EQ. 0 ) THEN
* Histograming
         CALL HepEvt_GetBeams(pf1,pf2)
         CALL HepEvt_GetFfins(qf1,qf2)
         CALL HepEvt_GetPhotAll(NphAll,PhoAll)
         CALL HepEvt_GetPhotIni(NphIni,PhoIni)
         CALL HepEvt_GetPhotFin(NphFin,PhoFin)
         CALL KK2f_GetWtAll(WtMain,WtCrud,WtSet)
         
         ss =  (pf1(4)+pf2(4))**2 -(pf1(3)+pf2(3))**2
     $        -(pf1(2)+pf2(2))**2 -(pf1(1)+pf2(1))**2
         ss2=  (qf1(4)+qf2(4))**2 -(qf1(3)+qf2(3))**2
     $        -(qf1(2)+qf2(2))**2 -(qf1(1)+qf2(1))**2
         vv=1d0-ss2/ss
         CALL GLK_Fil1(idv+1, vv, WtMain)

         CMSene = pf1(4)+pf2(4)
* Sum of energies of ALL photons
         EneTot = 0d0
         DO i=1,NphAll
            EneTot=EneTot+PhoAll(i,4)
         ENDDO
         EneTot = 2d0*EneTot/CMSene
         CALL GLK_Fil1(idv+2, EneTot, WtMain)

* Sum of energies of ISR photons
         EneIni = 0d0
         DO i=1,NphIni
            EneIni=EneIni+PhoIni(i,4)
         ENDDO
         EneIni = 2d0*EneIni/CMSene
         CALL GLK_Fil1(idv+3, EneIni, WtMain)

* Sum of energies of FSR photons
         EneFin = 0d0
         DO i=1,NphFin
            EneFin=EneFin+PhoFin(i,4)
         ENDDO
         EneFin = 2d0*EneFin/CMSene
         CALL GLK_Fil1(idv+4, EneFin, WtMain)

* Photon multiplicity (total)
         CALL GLK_Fil1(idv+7, NphAll*1.00001d0, WtMain)

* Flavour distribution (KF code)
         CALL HepEvt_GetKFfin(KFfin)
         xKF = KFfin+1d-6
         CALL GLK_Fil1(idv+9, xKF, WtMain)
*=======================================================================
      ELSE
*     finalization, printouts
         CALL GLK_Print(idv+1)
         CALL GLK_Print(idv+2)
         CALL GLK_Print(idv+3)
         CALL GLK_Print(idv+4)
         CALL GLK_Print(idv+7)
         CALL GLK_Print(idv+9)
      ENDIF
      END

      SUBROUTINE Robol2(mode)
*////////////////////////////////////////////////////////////////////////////////////////
*//                                                                                    //
*//   This is for tau chanel only                                                      //
*//   Old good tests on tau polarizations, double tau decay into pions                 //
*//                                                                                    //
*////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT REAL*8(a-h,o-z)
      COMMON / c_MainPro  / ninp,nout
*
      REAL*8  pf1(4),pf2(4),qf1(4),qf2(4)
      REAL*4 pp1,pp2,pbeam,vec1,vec2
      INTEGER NphAll,NphIni,NphFin
      REAL*8  PhoAll(100,4),PhoIni(100,4),PhoFin(100,4)
      REAL*8  WtSet(1000)
      REAL*8 Piminus(4),Piplus(4)
      DIMENSION p1(4),p2(4), v1(3),v2(3), pdif(3)
      DIMENSION b1(3),b2(3),db(3), xn(3),yn(3),qp(4),qm(4),ph(4)
      DIMENSION pp1(4),pp2(4),pbeam(4),vec1(4),vec2(4)
      REAL*4 rrr(2)
      LOGICAL ltrig
      DATA PI /3.1415926535897932D0/
      SAVE
*=======================================================================
      IF(mode .EQ. -1 ) THEN
*     book histograms
         nbv = 40
         vmin= 0
         vmax= 1
         idv  = 60000
         CALL GLK_Book1(idv+1,'pi+ energy     $',nbv,vmin,vmax)
         CALL GLK_Book1(idv+2,'pi- energy     $',nbv,vmin,vmax)
         CALL GLK_Book1(idv+3,'pi fast/slow     $',4,vmin,4*vmax)
         pbeam(1)=0.0
         pbeam(2)=0.0
         pbeam(3)=1.0
         pbeam(4)=1.0
         write(6,*) " Robol2: initialization"
         CALL  GLK_Book1(idv+100,'energy pions    $',50,0D0,1D0)
         CALL  GLK_Book1(idv+200,'flight distance $',50,0D0,1D0)
         CALL  GLK_Book1(idv+301,'impact parm 1,2 $',50,-0.03D0,0.03D0)
         CALL  GLK_Book1(idv+302,'impact parm 1-2 $',50,-0.03D0,0.03D0)
         CALL  GLK_Book1(idv+400,'impact parm 1-2 $',50, 0.D0,0.03D0)
         CALL  GLK_Book1(idv+600,'delta fi        $',100, -1.D0,1.D0)
         CALL  GLK_Book1(idv+601,'x*y             $',100, -.001D0,.001D0)
         CALL  GLK_Book1(idv+602,'x**2            $',100, -.001D0,.001D0)
         CALL  GLK_WtMon(-1,idv+10,0d0,1d0,20d0) !20D0 --> 1d0
         CALL  GLK_WtMon(-1,idv+11,0d0,1d0,1d0)
         CALL  GLK_WtMon(-1,idv+12,0d0,1d0,1d0)
         xlam= 0.2d0
         iev=0
         CALL GLK_Book1(idv+1001,'x-x            $',100, .00D0,10.001D0)
         CALL GLK_Book1(idv+1002,'x-y            $',100, .00D0,10.001D0)
         CALL GLK_Book1(idv+1003,'x-z            $',100, .00D0,10.001D0)
         CALL GLK_Book1(idv+1004,'y-x            $',100, .00D0,10.001D0)
         CALL GLK_Book1(idv+1005,'y-y            $',100, .00D0,10.001D0)
         CALL GLK_Book1(idv+1006,'y-z            $',100, .00D0,10.001D0)
         CALL GLK_Book1(idv+1007,'z-x            $',100, .00D0,10.001D0)
         CALL GLK_Book1(idv+1008,'z-y            $',100, .00D0,10.001D0)
         CALL GLK_Book1(idv+1009,'z-z            $',100, .00D0,10.001D0)
*=======================================================================
      ELSEIF(mode .EQ. 0 ) THEN
* Histograming
         CALL HepEvt_GetBeams(pf1,pf2)
         CALL HepEvt_GetFfins(qf1,qf2)
         CALL HepEvt_GetPhotAll(NphAll,PhoAll)
         CALL HepEvt_GetPhotIni(NphIni,PhoIni)
         CALL HepEvt_GetPhotFin(NphFin,PhoFin)
         CALL KK2f_GetWtAll(WtMain,WtCrud,WtSet)
         IF( WtCrud.EQ.0d0 ) RETURN
         CALL HepEvt_GetParticle(-211,1,Iadress,piminus)        
         CALL HepEvt_GetParticle( 211,1,Iadress,piplus)        
         xpiplus=piplus(4)/pf1(4)
         xpiminus=piminus(4)/pf1(4)
         a=0.
         b=0.
         if(xpiplus.gt.0.5d0) a=a+1
         if(xpiminus.gt.0.5d0) b=b+2
         CALL GLK_Fil1(idv+1, xpiplus         , 1d0)
         CALL GLK_Fil1(idv+2, xpiminus        , 1d0)
         CALL GLK_Fil1(idv+3, 0.5d0+a+b       , 1d0)
         ene=pf1(4)

c[[[         DO  K=1,4
c[[[            pp1(K)=piminus(k)
c[[[            pp2(K)=piplus(K)
c[[[            qm(k)=qf1(k)
c[[[            qp(k)=qf2(k)
c[[[            ph(k)=0d0
c[[[            DO l=1,nphini
c[[[               ph(k)=ph(k)+phoini(l,k)
c[[[            ENDDO
c[[[         ENDDO

         DO  K=1,4
            sgn=1d0
            if(k.eq.2) sgn=-1d0
            if(k.eq.3) sgn=-1d0
! qp pp1 are moments of tau+ and pi+ in this test 
! we rotate  by angle pi to the system of KORALB
! we take only events `semi-born'
           pp2(K)=sgn*piminus(k)
           pp1(K)=sgn*piplus(K)
            qm(k)=sgn*qf1(k)
            qp(k)=sgn*qf2(k)
            ph(k)=0d0
            DO l=1,nphini
               ph(k)=ph(k)+phoini(l,k)
            ENDDO
         ENDDO

         IF (ph(4).gt.1d-4) return
* asymmetry A_trans from APP B15 1984 1151
         CALL MULVT(PP1,PP2,VEC1)
         CALL MULVT(PBEAM,PP1,VEC2)
         ANGLE= XAKOL(VEC1,VEC2,3)
         IF (ANGLE.GT.90D0) ANGLE=180D0-ANGLE
         ANGLE=ANGLE-45D0
         CALL AASYM(1,ANGLE)
         x1=PP1(4)/ene
         x2=PP2(4)/ene
         IF((ABS(x1-0.5).LT.0.2D0).AND.(ABS(x2-0.5).LT.0.2D0)) THEN
            CALL AASYM(2,ANGLE)
         ENDIF
         IF ((QP(3).GT.0D0)
     $ .AND.((QP(2)).lT.0)
     $ .AND.((QP(1)).lT.0)
*     $ .AND.(ABS(QP(1)).LT.0.4*ABS(QP(2)))
*     $ .AND.(PH(4).GT.0D0).AND.(PH(4).LT.01D0
     $        )  THEN  
            DO I=1,3
               DO J=1,4
                  II=(1000+I+(J-1)*3)
                  xx=(pp1(I)*pp2(J))
                  call GLK_Fil1(idv+II,xx,1d0)
               ENDDO
            ENDDO
         ENDIF
         iev=iev+1
***** vertex generation
         CALL PseuMar_MakeVec(rrr,2)
         xlen1=-xlam*log(rrr(1))
         xlen2=-xlam*log(rrr(2))
         DO k=1,3
            v1(k)=qp(k)/qp(4)*xlen1
            v2(k)=qm(k)/qm(4)*xlen2
         ENDDO
         CALL GLK_Fil1(idv+200,xlen1,1d0)
         CALL GLK_Fil1(idv+200,xlen2,1d0)
*     
         DO k=1,4
            p1(k)=pp1(k)
            p2(k)=pp2(k)
         ENDDO
         c1 = p1(3)/p1(4)
         c2 = p2(3)/p2(4)
         z1 = p1(4)/ene
         z2 = p2(4)/ene
         ltrig = abs(c1).lt.0.8.and.abs(c2).lt.0.8
         ltrig = ltrig.and. z1.gt.0.2d0 .and. z2.gt.0.2d0 
         IF(ltrig) THEN
            CALL GLK_Fil1(idv+100,z1,1d0)
            CALL GLK_Fil1(idv+100,z2,1d0)
*     two-vector impact params
            a1 = (v1(1)*p1(1)+v1(2)*p1(2))/(p1(1)**2+p1(2)**2)
            a2 = (v2(1)*p2(1)+v2(2)*p2(2))/(p2(1)**2+p2(2)**2)
            DO k=1,2
               b1(k) = v1(k)-a1*p1(k)
               b2(k) = v2(k)-a2*p2(k)
               db(k) = b1(k)-b2(k)
               xn(k) = p1(k)-p2(k)
            ENDDO
*     scalar impact params
            xn(1)=xn(1)/sqrt(xn(1)**2+xn(2)**2)
            xn(2)=xn(2)/sqrt(xn(1)**2+xn(2)**2)
            yn(1)=  xn(2)
            yn(2)= -xn(1)
            bp1= yn(1)*b1(1)+yn(2)*b1(2)
            bp2= yn(1)*b2(1)+yn(2)*b2(2)
*     alternative scalar impact params
            bs1= sqrt(b1(1)**2+b1(2)**2)
            bs2= (b1(1)*b2(1)+b1(2)*b2(2))/bs1
            IF(iev.lt.20) WRITE(6,*) "Robol2: ***bp1,bp2=",bp1,bp2,bp1-bp2
            IF(iev.lt.20) WRITE(6,*) "Robol2:    bs1,bs2=",bs1,bs2,bs1-bs2
            CALL GLK_Fil1(idv+301,bp1,1d0)
            CALL GLK_Fil1(idv+301,bp2,1d0)
            CALL GLK_Fil1(idv+ 302,bp1-bp2,1d0)
* vector relative impact param
            rb = sqrt(db(1)**2+db(2)**2)
            CALL GLK_Fil1(idv+400,rb,1d0)
* azimuthal angles for tau-track and decay tracks
            phi1 = angfi(p1(1),p1(2))
            phi2 = angfi(p2(1),p2(2))
            psi1 = angfi(qp(1),qp(2))
            psi2 = angfi(qm(1),qm(2))
            delf1 = phi1-psi1
            delf2 = phi2-psi2
            IF(delf1.gt.pi)  delf1=delf1-2*pi
            IF(delf1.lt.-pi) delf1=delf1+2*pi
            IF(delf2.gt.pi)  delf2=delf2-2*pi
            IF(delf2.lt.-pi) delf2=delf2+2*pi
            IF(iev.lt.20) WRITE(6,*) "Robol2:  delf1,delf2=", delf1,delf2
            acc = 0d0
            IF(delf1*delf2.gt.0d0) acc=1d0
            CALL GLK_WtMon(0,idv+10,acc,1d0,0d0)
* aleph method, x-y variables  ===========
            DO k =1,3
               pdif(k) = p1(k)-p2(k)
            ENDDO
            sinthe = sqrt(pdif(1)**2 + pdif(2)**2)
     $           /sqrt(pdif(1)**2 + pdif(2)**2 +pdif(3)**2)
            delphi = phi1-phi2 +pi
            IF(delphi.gt.pi)  delphi=delphi-2*pi
            IF(delphi.lt.-pi) delphi=delphi+2*pi
            x = delphi * sinthe
            bs1 = sqrt( b1(1)**2 + b1(2)**2 )
            bs1 = bs1* sign(1d0,(p1(1)*b1(2)-p1(2)*b1(1)) )
            bs2 = sqrt( b2(1)**2 + b2(2)**2 )
            bs2 = bs2* sign(1d0,(p2(1)*b2(2)-p2(2)*b2(1)) )
            y   = bs2-bs1 
            CALL GLK_WtMon(0,idv+11, x*y,1d0,0d0)
            CALL GLK_WtMon(0,idv+12, x*x,1d0,0d0)
            CALL GLK_Fil1(idv+600,delphi,1D0)
            CALL GLK_Fil1(idv+601,x*y,1D0)
            CALL GLK_Fil1(idv+602,x**2,1D0)
         ENDIF
*=======================================================================
      ELSE
*     finalization, printouts
         CALL GLK_Yminim(idv+1,0d0)
         CALL GLK_Yminim(idv+2,0d0)
         CALL GLK_Yminim(idv+3,0d0)
         CALL GLK_Print(idv+1)
         CALL GLK_Print(idv+2)
         CALL GLK_Print(idv+3)
         
         WRITE(6,*) " Robol2: postgeneration"
         DO j=1,9
            II=1000+j
            CALL GLK_Print(idv+II)
         ENDDO
         CALL GLK_Yminim(idv+100,0d0)
         CALL GLK_Yminim(idv+200,0d0)
         CALL GLK_Yminim(idv+301,0d0)
         CALL GLK_Yminim(idv+302,0d0)
         CALL GLK_Yminim(idv+400,0d0)
         CALL GLK_Yminim(idv+600,0d0)
         CALL GLK_Yminim(idv+601,0d0)
         CALL GLK_Yminim(idv+602,0d0)
         CALL GLK_Ymaxim(idv+601,4293D0)
         CALL GLK_Print(idv+100)
         CALL GLK_Print(idv+200)
         CALL GLK_Print(idv+301)
         CALL GLK_Print(idv+302)
         CALL GLK_Print(idv+400)
         CALL GLK_Print(idv+600)
         CALL GLK_Print(idv+601)
         CALL GLK_Print(idv+602)
         CALL GLK_WtMon(3,idv+10,aver,erela,dum1)
*
         bgam = ene/amf2
         CALL GLK_WtMon(1,idv+11,avexy,davexy,dum1)
         CALL GLK_WtMon(1,idv+12,avexx,davexx,dum1)
         xltau = avexy/avexx
         dltau = sqrt( davexy**2 + davexx**2)
         WRITE(   6,"(a,f10.6,a,f10.6)") "xltau= ", xltau, " +-",xltau*dltau
         WRITE(nout,"(a,f10.6,a,f10.6)") "xltau= ", xltau, " +-",xltau*dltau
         CALL aasym(-1,x)
         CALL aasym(-2,x)
      ENDIF
      END

      DOUBLE PRECISION FUNCTION XAKOL(X,Y,N)
*     **************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      REAL*4    X(*),Y(*)
      DIMENSION X1(4),Y1(4)
      DATA PI /3.1415926535897932D0/
      S=0.D0
      X2=0.D0
      Y2=0.D0
      DO I=1,N
         X1(I)=X(I)
         Y1(I)=Y(I)
      ENDDO
      DO I=1,N
         S=S+X1(I)*Y1(I)
         X2=X2+X1(I)**2
         Y2=Y2+Y1(I)**2
      ENDDO
      XAKOL=ACOS(-S/SQRT(X2*Y2))*180.D0/PI
      END

      SUBROUTINE MULVT(X,Y,R)
*     **********************
      DIMENSION  X(4),Y(4),R(4)
      R(1)=X(2)*Y(3)-X(3)*Y(2)
      R(2)=X(3)*Y(1)-X(1)*Y(3)
      R(3)=X(1)*Y(2)-X(2)*Y(1)
      R(4)=0.0
      END

      SUBROUTINE AASYM(N,ZZ)
*     **********************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      COMMON / c_MainPro  / ninp,nout
      DIMENSION NENT(5),NSUM(5)
      DATA NENT /0,0,0,0,0/
      DATA NSUM /0,0,0,0,0/

      IF( (IABS(N).GT.5) .OR. (N.EQ.0) ) RETURN
      IF(N .GT. 0) THEN
         Z=ZZ
         NENT(N)=NENT(N)+1
         IF(Z.GT.0.) NSUM(N)=NSUM(N)+1
      ELSE
         J=IABS(N)
         A1=NSUM(J)
         A2=NENT(J)-A1
         AA=A1+A2+1.
         ASYM=(A1-A2)/AA *100.
         DASYM=2.*SQRT(A1*A2/AA**3)*100.
         WRITE(nout, 1000) IABS(N),ASYM,DASYM
         WRITE(   *, 1000) IABS(N),ASYM,DASYM
      ENDIF
 1000 FORMAT(1H0,5X,4HASYM,I1,F10.1,4X,2H+-,F4.1)
      END
