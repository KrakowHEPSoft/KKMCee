*      CALL ZFTEST(0)
*  CALL ZFLEP2(IBOX,IAD,NV)
*      IBOX=0,1,2 as usual, NO, ADD, EWFF
*       IAD=1 - study of angular dependence of weak boxes, NV-irrelevant
*       IAD=0 - no boxes study, RO's in accordance with NV's
      CALL ZFLTMN(1,0,-3)
*
      END
 
      SUBROUTINE ZFLTMN(IBOX,IAD,NV)
*     ========== ===================
************************************************************************
*                                                                      *
*     Main used for PCP Winter 98 project                              *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 VECARG,VECRE2,VECRE4,VECRE5,VECRE6,VECRE9,VECR10
      COMPLEX*16 XVPOL,XFOTF3
      COMPLEX*16 XALLCH,XFOTF
*
* DIMENSIONs of variable dimension
*
* LEP1
cb    PARAMETER(NECMP=8 ,NCOS=7)
* LEP2
cb    PARAMETER(NECMP=16,NCOS=7)
* LEP2 reduced
      PARAMETER(NECMP=11,NCOS=7,NMH=7,NMCUT=7)
*
      DIMENSION VECARG(190)
     *         ,VECRE2(190)
     *         ,VECRE4(190)
     *         ,VECRE5(190)
     *         ,VECRE6(190)
     *         ,VECRE9(190)
     *         ,VECR10(190)
      DIMENSION RSARR(NECMP),RCOS(NCOS),ARMH(NMH),RCUTA(NMCUT)
      DIMENSION INDC(7),DXSEC(5,NECMP,NCOS)
      DIMENSION REOBSV(17,NECMP)
      DIMENSION REOBMI(17,NECMP),REOBAV(17,NECMP),REOBMA(17,NECMP)
      DIMENSION ANGLEL(5)
*
* ZFITTER constants
*
      PARAMETER(ZMASS=91.1867D0,TMASS=175D0,ALFAS=.120D0)
*
      PARAMETER(RSMN=87.D0,DRS=1.D0,NRS=9)
*
* ZFITTER common blocks
*
      COMMON /ZUPARS/QDF,QCDCOR(0:14),ALPHST,SIN2TW,S2TEFF(0:11),
     & WIDTHS(0:11)
      COMMON /CDZRKZ/ARROFZ(0:10),ARKAFZ(0:10),ARVEFZ(0:10),ARSEFZ(0:10)
      COMMON /EWFORM/XALLCH(5,4),XFOTF
*
      COMMON/FORFUZ/RS,ZMASSC,TMASSC,HMASSC,ALQED,ALFASC,INDF
      COMMON /ZUDATA/ ACOLIN(0:11),EF_MIN(0:11),SPRIME(0:11),
     + ANGMIN(0:11),ANGMAX(0:11),IRCUTS(0:11),IRFAST(0:11)
*
      EXTERNAL FUZAN 
*
* Array of LEP2 energies
*
      DATA RSARR/5*91d0,100d0,140d0,175d0,183d0,189d0,200d0/
      DATA RCUTA/0.01D0,0.1D0,0.3D0,0.5D0,0.7D0,0.9D0,0.99D0/
*
* Array of COS of scattering angles to see the angular distribution
*
      DATA RCOS /-.99D0,-.9D0,-.5D0,0D0,+.5D0,+.9D0,+.99D0/
*
* Array of Higgs energies
*
      DATA ARMH/76D0,10D0,30D0,76D0,100D0,300D0,1000D0/
*
* Selecled channel indices: 2-muon, 4-up, 5-down, 6-cc, 9-bb, 10-hadrons
*
      DATA INDC/2,3,4,5,6,9,10/
      DATA ANGLEL/0d0,5d0,10d0,20d0,40d0/
*-----------------------------------------------------------------------
c      OPEN (UNIT=21,FILE='for_lk2_a.dat')
c      OPEN (UNIT=22,FILE='for_lk2_2.dat')
c      OPEN (UNIT=24,FILE='for_lk2_4.dat')
c      OPEN (UNIT=25,FILE='for_lk2_5.dat')
c      OPEN (UNIT=26,FILE='for_lk2_6.dat')
c      OPEN (UNIT=29,FILE='for_lk2_9.dat')
c      OPEN (UNIT=30,FILE='for_lk2_10.dat')
*
* Initialization
*
      ZMASSC=ZMASS
      TMASSC=TMASS
      HMASSC=HMASS
      ALFASC=ALFAS
*
* Limits in \cos of scattering angle
*
      COTMIN=-.99999D0
      COTMAX=+.99999D0
*
* This job with five nines ran about 1 hour. The results differ within
* integration error, 1e-4, from those with four nines. 
* The latter ran about 45'.
* The code fell in a infinite do loop (was killed after 17 hours of running)
* with SIX nines. Although, the results with five nines look better than
* with four nines, I got a feeling (while analyzing them and given the fact
* of increasing of computational time from 4 to 5 9's) that the running with
* 5 9's is `at the edge' of numerical stability. Note, the program checked
* that the differential cross section was always positive. 
*
      CALL ZUINIT
*
* Set ZFITTER flags and print flag values
*
      PI=4*ATAN(1D0)
      ALPHAI=1/137.035 989 5d0
*
*
* setup, C3ps
*
       CALL ZUFLAG('PRNT',0) !
       CALL ZUFLAG('BORN',0) !
       CALL ZUFLAG('DIAG',1) !
       CALL ZUFLAG('POWR',1) !
       CALL ZUFLAG('BOXD',1) !
       CALL ZUFLAG('PREC',10)!
       CALL ZUFLAG('INTF',1) 
       CALL ZUFLAG('FINR',1)   
       CALL ZUFLAG('FOT2',3) 
       CALL ZUFLAG('ISPP',0) 
       CALL ZUFLAG('FSRS',0)
*
      CALL ZUINFO(0)
*
      ICUT=1
      RATS=1d-2
cb      RATS=0.8d0
      Acol=10d0
cb    Acol=25d0
      DAL5H=2.8039808929734D-02
*
c      PRINT *
c      PRINT *,'Deconvoluted differential/total:'
c      PRINT *
*
c      PRINT *
c      PRINT 1110,IBOX
c      PRINT * 
      PRINT *
      PRINT 1100,ICUT
      PRINT *  
*   
      DO IMH=5,5
      HMASS=ARMH(IMH)
      PRINT *
      PRINT *
      PRINT *,'M_H=',HMASS
      PRINT *
*
* Static quantities
*
       CALL ZUWEAK(ZMASS,TMASS,HMASS,DAL5H,ALFAS)
*
       XVPOL=ALPHAI/(1-XFOTF3(3,3,1,0,0,DAL5H,-ZMASS**2)*ALPHAI/4/PI)
       ALQED=DREAL(XVPOL)
       print *,DREAL(1D0/XVPOL)
*
       INDMAX=7
       LAMAX =1
       IF(ICUT.EQ.0.OR.ICUT.EQ.2.OR.ICUT.EQ.3) INDMAX=1 
       IF(ICUT.EQ.0.OR.ICUT.EQ.3) LAMAX=5
*
       do la=1,lamax
       ang0=0d0  +anglel(la)
       ang1=180d0-anglel(la)
       print *,'theta=',anglel(la)
*
*
       DO IS = 10,10
*
* LEP1-array
*
         IF(IS.EQ.1) RS = ZMASS-3D0
         IF(IS.EQ.2)    THEN
           IF(IMH.EQ.1) THEN 
                     RS = ZMASS-1D0
           ELSE
                     RS = ZMASS-1.8D0
           ENDIF
         ENDIF
         IF(IS.EQ.3) RS = ZMASS
         IF(IS.EQ.4)    THEN 
           IF(IMH.EQ.1) THEN
                     RS = ZMASS+1D0
           ELSE
                     RS = ZMASS+1.8D0
           ENDIF
         ENDIF        
         IF(IS.EQ.5) RS = ZMASS+3D0
*
* LEP2-array
*
         IF(IS.GT.5) RS = RSARR(IS)
*
         S=RS**2
*
        DO IRC=1,NMCUT
           SPR=RCUTA(IRC)*S
           PPR=SPR
           print *,'CUT=',1D0-RCUTA(IRC)
        DO IND=1,INDMAX
          INDF=INDC(IND)
          PRINT *
          PRINT 1111,INDF
          PRINT *
*
        IF(ICUT.EQ.-1.OR.ICUT.EQ.1) THEN
         CALL ZUCUTS(INDF,ICUT,0D0,0D0,SPR,0D0,180D0,PPR)
        ELSE
         CALL ZUCUTS(2,ICUT,Acol,1d0,0d0,ang0,ang1,0.25*S)
        ENDIF
*
* Angular distribution
*
         IF(IAD.EQ.1) THEN
          DO ICOS=1,NCOS
           CSA=RCOS(ICOS)
           CALL ZUATSM(INDF,RS,ZMASS,TMASS,HMASS,DAL5H,ALFAS,CSA,DXS)
           DXSEC(IND,IS,ICOS)=DXS
           PRINT 1112,RS,CSA,DXS
          ENDDO
         ENDIF
*
         IF(IBOX.NE.2) THEN
          CALL ZUTHSM(INDF,RS,ZMASS,TMASS,HMASS,DAL5H,ALFAS,TCS,AFB)
          PRINT 1113,RS,TCS,AFB
          vecarg(IS)   = RS
          IF(INDF.EQ.2) vecre2(IS)=TCS
          IF(INDF.EQ.4) vecre4(IS)=TCS
          IF(INDF.EQ.5) vecre5(IS)=TCS
          IF(INDF.EQ.6) vecre6(IS)=TCS
          IF(INDF.EQ.9) vecre9(IS)=TCS
          IF(INDF.EQ.10)vecr10(IS)=TCS
         ELSE
          STECOT=.25D0*(COTMAX-COTMIN)
          CALL 
     &    SIMPV(COTMIN,COTMAX,STECOT,1D-4,1D-20,FUZAN,ARG,TCS,TC2,TC3)
          PRINT 1114,RS,TCS
         ENDIF
        ENDDO
       ENDDO
      ENDDO
      ENDDO
      ENDDO
*cbardin
cb      write(21,FMT=998) vecarg
cb      write(22,FMT=999) vecre2
cb      write(24,FMT=999) vecre4
cb      write(25,FMT=999) vecre5
cb      write(26,FMT=999) vecre6
cb      write(29,FMT=999) vecre9
cb      write(30,FMT=999) vecr10
cb 998  format(2x,f10.5)
cb 999  format(1x,e15.5)
*
cb      CLOSE (UNIT=21)
cb      CLOSE (UNIT=22)
cb      CLOSE (UNIT=24)
cb      CLOSE (UNIT=25)
cb      CLOSE (UNIT=26)
cb      CLOSE (UNIT=29)
cb      CLOSE (UNIT=30)
*
 1100 FORMAT(1X,'ICUT=',I2,2X,'RATS=',D8.1)
 1110 FORMAT(1X,'IBOX=',I2)
 1111 FORMAT(1X,'INDF=',I2)
 1112 FORMAT(' Ecm=',F8.4,'   COS(THETA)=',F6.2,'   XS=',D12.5 )
 1113 FORMAT(' Ecm=',F8.4,'   INTEGR. XS=',D15.7,'  AFB=',D12.5)
 1114 FORMAT(' Ecm=',F8.4,'   INTEGR. XS=',D15.7)
*
      STOP
*
      END

      FUNCTION FUZAN(CSA)
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      COMMON/FORFUZ/RS,ZMASS,TMASS,HMASS,ALQED,ALFAS,INDF
*
      CALL ZUATSM(INDF,RS,ZMASS,TMASS,HMASS,ALQED,ALFAS,CSA,DXS)
*
      FUZAN=ABS(DXS)
*     
      if(dxs.lt.0d0) then
        print *,'CSA,DXS=',CSA,DXS
      endif
*
      END
