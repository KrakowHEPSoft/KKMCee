      PROGRAM ZWIDTHTB
*
* this is to run 5_13 to reproduce Tables from Note with Giampiero
* it should not reproduce CERN 95-03 numbers
* f77 zwidthtb.f dizet5_13.f m2tcor5_11.f bcqcdl.f bkqcdl.f bhang4_61.f 
* zf512_aux.f -o zwidthtb.exe
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
      DIMENSION XFZ(4)
*
      PARAMETER(NOBS=35,NAMT=1,NAMH=1,NALS=2)
      PARAMETER(ALFAI=137.0359895D0)
*
      COMMON/CDZRKZ/ARROFZ(0:10),ARKAFZ(0:10),ARVEFZ(0:10),ARSEFZ(0:10)
     &             ,AROTFZ(0:10),AIROFZ(0:10),AIKAFZ(0:10),AIVEFZ(0:10)
      COMMON/CDZPHM/CMASS,BMASS,TMASS
*
      DIMENSION NPAR(21),ZPAR(31),PARTZ(0:11),PARTW(3)
      DIMENSION ARMT(NAMT),ARMH(NAMH),ARAL(NALS)
      DIMENSION TABLE (NOBS,NALS,NAMH,NAMT)
      DIMENSION TABLEA(NOBS,NALS,NAMH,NAMT)
      DIMENSION TABLEX(NOBS,NALS,NAMH,NAMT)
      DIMENSION TABLEN(NOBS,NALS,NAMH,NAMT)
*
cbard DATA ARMT/150D0,160D0,170D0,180D0,190D0,200D0,210D0/
cbard DATA ARMH/65D0,100D0,140D0,300D0,1000D0/
cbard DATA ARMT/175.6D0/
cbard DATA ARMT/175D0/
cbard DATA ARMT/168.8D0/
      DATA ARMT/173.8D0/
cbard DATA ARMT/178.8D0/
cbard DATA ARMH/300D0/
      DATA ARMH/100D0/
cbard DATA ARMH/76D0,100D0,300D0/
cbard DATA ARMH/10D0,30D0,76D0,100D0,300D0,1000D0/
cbard DATA ARAL/.00D0,.117D0,.123D0,.129D0/
      DATA ARAL/0D0,.1190D0/
cbard DATA ARAL/0D0,.1200D0/
cbard DATA ARAL/.1194D0/
cbard DATA ARAL/.1200D0/
cbard DATA ARAL/.1250D0/
*
* Default setting of integer flags
*
      IHVP =1
      IAMT4=4
      IQCDD=3
      IMOMS=1
      IMASS=0
      ISCRE=0
      IALEM=3
      IMASK=0
      ISCAL=0
      IBARB=2
      IFTJR=1
      IFACR=0
      IFACT=0
      IHIGS=0
      IAFMT=1
      IEWLC=1
      ICZAK=1
      IHIG2=0
      IALE2=3
      IGFER=2
      IDDZZ=1
*
      PRINT 1,IHVP,IAMT4,IQCDD,ISCRE,IFACT,ISCAL,IALEM,IBARB,IFTJR
 1    FORMAT(1X,'LEP-I PCWG TABLES PRODUCTION RUN',/,1X,'FLAGS:'
     &  ,1X,'IHVP,IAMT4,IQCD,ISCRE,IFACT,ISCAL,IALEM,IBARB,IFTJR =',9I2)
* Preferrred setting of parameter vector NPAR
      NPAR(1)=IHVP
      NPAR(2)=IAMT4
      NPAR(3)=IQCDD
      NPAR(4)=IMOMS
      NPAR(5)=IMASS
      NPAR(6)=ISCRE
      NPAR(7)=IALEM
      NPAR(8)=IMASK
      NPAR(9)=ISCAL
      NPAR(10)=IBARB
      NPAR(11)=IFTJR
      NPAR(12)=IFACR
      NPAR(13)=IFACT
      NPAR(14)=IHIGS
      NPAR(15)=IAFMT
      NPAR(16)=IEWLC
      NPAR(17)=ICZAK
      NPAR(18)=IHIG2
      NPAR(19)=IALE2
      NPAR(20)=IGFER
      NPAR(21)=IDDZZ
*
      AMZ=91.1867D0
cbard AMZ=91.1888D0
      AMW=80.000D0
*
      IOPT=0
*
* \xi -effect
      DO 1000 iii=1,3,2
* in 1998 only iii=3 is allowed
* since AFMT is known and confirmed 
cbard iii=1
      IF(iii.EQ.1) then
      NPAR(9)=0
      NPAR(15)=1 
      endif
      IF(iii.EQ.2) then
      NPAR(9)=1
      NPAR(15)=0
      endif
      IF(iii.EQ.3) then
      NPAR(9)=4
      NPAR(15)=0
      endif
*
* amt4-effect
* 1994
* only IAMT4=2 is excluded,
* since it gives abnormal remainder in \rho and \kappa, see run #1
* 1998
* only IAMT4=4
cbard DO 1000 IAMT4=4,1,-1
      IAMT4=4
      NPAR(2)=IAMT4
* Scale of the remainder terms both \Delta r and \rho, \kappa
* ISCRE=2: R(1-\Delta\alpha) scale in \Delta r
      DO 1000 ISCRE=0,2,1
cbard ISCRE=0
      NPAR(6)=ISCRE
* Factorized/Expanded \Delta r
* IFACR=2,3 are excluded since they contradict to Sirlin's study
* and produce artifacts of positive M_W shift at low m_t
      DO 1000 IFACR=0,2,1
cbard IFACR=0
      NPAR(12)=IFACR
* Factorized/Expanded \rho and \kappa
* IFACT=1,2,3,4 are excluded,
* since they produce more narrow bands then IFACT=5, fully expanded
      DO 1000 IFACT=0,2,1
cbard IFACT=0
      NPAR(13)=IFACT
* Non-resummed/resummed (1) leading higgs contribution: M_H>M_W*e^(5/12)
      DO 1000 IHIGS=0,1
cb    IHIGS=0
      NPAR(14)=IHIGS
      DO 1000 IHIG2=0,1
cb    IHIG2=0
      NPAR(18)=IHIG2
cbard DO 1000 IALE2=3,0,-3
      IALE2=3
      NPAR(19)=IALE2
*
      IOPT=IOPT+1
*
* In this run IA is used for DD-ZZ game
cbard DO 100 IA=1,NALS
cbard IF(IA.EQ.1) THEN
cbard   IQCD=0
cbard     ELSE
      ia=2
        IQCD=IQCDD
        IDDZZ=IA-1
cbard ENDIF
      NPAR(3) =IQCD
      NPAR(21)=IDDZZ
      ALPHST=ARAL(2)
      DO 100 IH=1,NAMH
      AMH=ARMH(IH)
      DO 100 IT=1,NAMT
      AMT=ARMT(IT)
cbard AMT=115D0+(IT-1)*.1D0
      ARMT(IT)=AMT
*
* value for IALE2=0
*
      DAL5H=0.0279801D0
      DAL5H=0.0280 D0
c     DAL5H=0.027963 D0
c     DAL5H=0.028039808929734D0
c     DAL5H=0.0280D0
*
      CALL DIZET
     &(NPAR,AMW,AMZ,AMT,AMH,DAL5H,ALQED,ALPHST,ALPHTT,ZPAR,PARTZ,PARTW)
      DR   =ZPAR(1)
      DRREM=ZPAR(2)
      SW2  =ZPAR(3)
*
      S  = AMZ**2
      Q2 = S/2D0
      U2 = Q2-S
      CALL ROKANC(0,0,U2,-S,-Q2,-1D0,-1D0,XFZ,XFZT,XFZTA)
      IF(IOPT.EQ.1.AND.IA.EQ.1.AND.IT.EQ.1.AND.IH.EQ.1) THEN
      PRINT 2,AMZ,ALFAI*DREAL(2D0-XFZTA)
     &           ,ALFAI*DREAL(2D0-XFZT ),CMASS,BMASS
 2    FORMAT(1X,'INPUT:',1X,'MZ=',F7.4,
     &2X,'ALFAI(5)=',F7.3,2X,'ALFAI(6)=',F7.3,
     &2X,'MC=',F4.2,2X,'MB=',F4.2)
*
      print *,'dal5h=',dal5h
      print 21,amt,alphst,armh
 21   format(1x,'AMT=',f6.2,2x,'ALST=',f6.4,2x,'AMH=',6f9.3)
c     do iz=16,30
c     print 22,iz,zpar(iz)
c22   format(1x,'iz=',i2,2x,f8.6)
c     enddo
      ENDIF
*
* observables
*
      VLAL=1-4D0  *ARSEFZ(2)
      VBAB=1-4D0/3*ARSEFZ(9)
      VCAC=1-8D0/3*ARSEFZ(6)
      AFBL=3*VLAL**2/(VLAL**2+1)**2
*
* PO's conventions game
*
c      print *,'VLAL=',VLAL
c      print *,'VECL=',ARVEFZ(2)
c      print *,'ILAL=',-4D0*SW2*AIKAFZ(2)
c      print *,'IVEC=',AIVEFZ(2)
c      print *,'AFBL=',AFBL
      AFBLi=3*VLAL**2/(VLAL**2+(AIVEFZ(2))**2+1)**2
c     print *,'AFBL,AFBLi=',AFBL,AFBLi
c      stop
*
      AFBB=3*VLAL*VBAB/(VLAL**2+1)/(VBAB**2+1)
      AFBC=3*VLAL*VCAC/(VLAL**2+1)/(VCAC**2+1)
      ALRL=2*VLAL/(VLAL**2+1)
      ALRB=2*VBAB/(VBAB**2+1)
      ALRC=2*VCAC/(VCAC**2+1)
      ZMASS=AMZ
      GAMEE=PARTZ(1)
      GAMUU=PARTZ(4)
      GAMDD=PARTZ(5)
      GAMCC=PARTZ(6)
      GAMBB=PARTZ(9)
      GZHAD=PARTZ(10)
      GZTOT=PARTZ(11)
      RHAD=GZHAD/GAMEE
      CNANOB=.389386D6
      PI=4D0*ATAN(1D0)
      CQEDEE=1D0+3D0/4*ALQED
      CQEDUU=1D0+3D0/4*ALQED*4D0/9
      CQEDDD=1D0+3D0/4*ALQED*1D0/9
      SIGH=12*PI*GZHAD       *GAMEE       /GZTOT**2/ZMASS**2*CNANOB
      SIEE=12*PI*GAMEE/CQEDEE*GAMEE/CQEDEE/GZTOT**2/ZMASS**2*CNANOB
      SIUU=12*PI*GAMEE/CQEDEE*GAMUU/CQEDUU/GZTOT**2/ZMASS**2*CNANOB
      SIDD=12*PI*GAMEE/CQEDEE*GAMDD/CQEDDD/GZTOT**2/ZMASS**2*CNANOB
      SICC=12*PI*GAMEE/CQEDEE*GAMCC/CQEDUU/GZTOT**2/ZMASS**2*CNANOB
      SIBB=12*PI*GAMEE/CQEDEE*GAMBB/CQEDDD/GZTOT**2/ZMASS**2*CNANOB
* T6
      TABLE( 1,IA,IH,IT)=AMW
      TABLE( 2,IA,IH,IT)=SIGH
      TABLE( 3,IA,IH,IT)=PARTZ(0)/1D3
      TABLE( 4,IA,IH,IT)=PARTZ(1)/1D3
      TABLE( 5,IA,IH,IT)=PARTZ(2)/1D3
      TABLE( 6,IA,IH,IT)=PARTZ(3)/1D3
      TABLE( 7,IA,IH,IT)=PARTZ(4)/1D3
      TABLE( 8,IA,IH,IT)=PARTZ(5)/1D3
      TABLE( 9,IA,IH,IT)=PARTZ(6)/1D3
      TABLE(10,IA,IH,IT)=PARTZ(9)/1D3
      TABLE(11,IA,IH,IT)=GZHAD/1D3
      TABLE(12,IA,IH,IT)=3*PARTZ(0)/1D3
      TABLE(13,IA,IH,IT)=GZTOT/1D3
* T7
      TABLE(14,IA,IH,IT)=RHAD
      TABLE(15,IA,IH,IT)=GAMBB/GZHAD
      TABLE(16,IA,IH,IT)=GAMCC/GZHAD
      TABLE(17,IA,IH,IT)=ARSEFZ(2)
      TABLE(18,IA,IH,IT)=ARSEFZ(9)
      TABLE(19,IA,IH,IT)=ARSEFZ(6)
      TABLE(20,IA,IH,IT)=AFBL
      TABLE(21,IA,IH,IT)=AFBB
      TABLE(22,IA,IH,IT)=AFBC
      TABLE(23,IA,IH,IT)=ALRL
      TABLE(24,IA,IH,IT)=ALRB
      TABLE(25,IA,IH,IT)=ALRC
      TABLE(26,IA,IH,IT)=AROTFZ(2)
      TABLE(27,IA,IH,IT)=AROTFZ(9)
      TABLE(28,IA,IH,IT)=AROTFZ(6)
*Additional
      TABLE(29,IA,IH,IT)=SW2
      TABLE(30,IA,IH,IT)=SIGH/RHAD
      TABLE(31,IA,IH,IT)=12*PI*GAMEE*GAMEE/GZTOT**2/ZMASS**2*CNANOB
      TABLE(32,IA,IH,IT)=12*PI*GAMEE*GAMUU/GZTOT**2/ZMASS**2*CNANOB
      TABLE(33,IA,IH,IT)=12*PI*GAMEE*GAMDD/GZTOT**2/ZMASS**2*CNANOB
      TABLE(34,IA,IH,IT)=12*PI*GAMEE*GAMCC/GZTOT**2/ZMASS**2*CNANOB
      TABLE(35,IA,IH,IT)=12*PI*GAMEE*GAMBB/GZTOT**2/ZMASS**2*CNANOB

  100 CONTINUE
*
      IF(IOPT.EQ.1) THEN
       DO 110 IO=1,NOBS
       DO 110 IA=1,NALS
       DO 110 IH=1,NAMH
       DO 110 IT=1,NAMT
        TABLEA(IO,IA,IH,IT)=TABLE(IO,IA,IH,IT)
        TABLEX(IO,IA,IH,IT)=TABLE(IO,IA,IH,IT)
        TABLEN(IO,IA,IH,IT)=TABLE(IO,IA,IH,IT)
 110   CONTINUE
      ELSE
       DO 120 IO=1,NOBS
       DO 120 IA=1,NALS
       DO 120 IH=1,NAMH
       DO 120 IT=1,NAMT
        IF(TABLE(IO,IA,IH,IT).LT.TABLEN(IO,IA,IH,IT)) THEN
         TABLEN(IO,IA,IH,IT)=TABLE(IO,IA,IH,IT)
        ENDIF
        IF(TABLE(IO,IA,IH,IT).GT.TABLEX(IO,IA,IH,IT)) THEN
         TABLEX(IO,IA,IH,IT)=TABLE(IO,IA,IH,IT)
        ENDIF
  120  CONTINUE
      ENDIF
 1000 CONTINUE
*
*-----------------------------------------------------------------------
*
* PRINTING BLOCK FOR PLOTS, THREE SETS OF TABLES: PREFERRED, A, miN, maX
*
      DO 60 IOBZ=1,NOBS
      DO 50 IALS=1,NALS
      DO 40 IAMT=1,NAMT
      PRINT 30,IOBZ,ARAL(IALS),ARMT(IAMT),
     &                        (TABLEA(IOBZ,IALS,IAMH,IAMT),IAMH=1,NAMH)
  30  FORMAT(1X,I2,1X,F5.3,1X,F5.1,1X,6(D13.6))
  40  CONTINUE
  50  CONTINUE
  60  CONTINUE
*
* Construction of DD-ZZ analogs
*
      DO IH=1,NAMH
        AMH=ARMH(IH)
        GAMEE0=TABLE( 5,1,IH,1)*1D3
        GAMUU0=TABLE( 7,1,IH,1)*1D3
        GAMDD0=TABLE( 8,1,IH,1)*1D3
        GAMCC0=TABLE( 9,1,IH,1)*1D3
        GAMBB0=TABLE(10,1,IH,1)*1D3
        GAMZT =TABLE(13,2,IH,1)*1D3
        CQEDEE=1D0
        CQEDUU=1D0
        CQEDDD=1D0
        SIEE=12*PI*GAMEE0/CQEDEE*GAMEE0/CQEDEE/GAMZT**2/ZMASS**2*CNANOB
        SIUU=12*PI*GAMEE0/CQEDEE*GAMUU0/CQEDUU/GAMZT**2/ZMASS**2*CNANOB
        SIDD=12*PI*GAMEE0/CQEDEE*GAMDD0/CQEDDD/GAMZT**2/ZMASS**2*CNANOB
        SICC=12*PI*GAMEE0/CQEDEE*GAMCC0/CQEDUU/GAMZT**2/ZMASS**2*CNANOB
        SIBB=12*PI*GAMEE0/CQEDEE*GAMBB0/CQEDDD/GAMZT**2/ZMASS**2*CNANOB
        PRINT 70,AMH,SIEE,SIUU,SIDD,SICC,SIBB
 70     FORMAT(1X,'MH=',F6.1,2X,'DD-ZZ=',5F10.6)
      ENDDO
*
       DO 61 IOBZ=1,NOBS
       DO 51 IALS=1,NALS
       DO 41 IAMT=1,NAMT
       PRINT 31,IOBZ,ARAL(IALS),ARMT(IAMT),
     &                        (TABLEN(IOBZ,IALS,IAMH,IAMT),IAMH=1,NAMH)
   31  FORMAT(1X,I2,1X,F5.3,1X,F5.1,1X,6(D13.6))
   41  CONTINUE
   51  CONTINUE
   61  CONTINUE
*
       DO 62 IOBZ=1,NOBS
       DO 52 IALS=1,NALS
       DO 42 IAMT=1,NAMT
       PRINT 32,IOBZ,ARAL(IALS),ARMT(IAMT),
     &                        (TABLEX(IOBZ,IALS,IAMH,IAMT),IAMH=1,NAMH)
   32  FORMAT(1X,I2,1X,F5.3,1X,F5.1,1X,6(D13.6))
   42  CONTINUE
   52  CONTINUE
   62  CONTINUE
*
* END OF PRINTING BLOCK FOR PLOTS
*
*-----------------------------------------------------------------------
*
      END


