*
* MAIN to call DIZET
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION NPARD(25),ZPARD(30)
      real*8 T_MASS,Z_MASS,W_MASS,H_MASS,DAL5H,V_TB,ALQED,ALST,dm(5)
      DATA dm/-2.,-1.,0.,1.,2./
*
      COMMON /CDZRKZ/ARROFZ(0:10),ARKAFZ(0:10),ARVEFZ(0:10),ARSEFZ(0:10)
     &              ,AROTFZ(0:10),AIROFZ(0:10),AIKAFZ(0:10),AIVEFZ(0:10)

*
      COMMON/PARTZW/PARTZ(0:11),PARTW(3)
      COMMON /ZUPARS/QDF,QCDCOR(0:14),ALPHST,SIN2TW,S2TEFF(0:11),
     & WIDTHS(0:11)
      CHARACTER FCHAN(0:11)*10,FCHANW(1:3)*10
      DATA FCHAN/'nu,nubar','e+,e-','mu+,mu-','tau+,tau-','u,ubar',
     + 'd,dbar','c,cbar','s,sbar','t,tbar','b,bbar','hadron','total'/
      DATA FCHANW/'lept,nubar','down,ubar','total'/
      
      PARAMETER(ALFAI=137.0359895D0,ALFA=1.D0/ALFAI,CONS=1.D0)
      PARAMETER(ZMASS=91.1876D0,TMASS=172.9D0,HMASS=125.1D0)
      PARAMETER(ALFAS=.1181D0)
*     1 - old default, 4,5 - new options for dal5h      
      NPARD(1) = 1
*     4 - old default, 6,7,8 - new option for SIN2_EFF
      NPARD(2) = 4
      NPARD(3) = 3
      NPARD(4) = 1
      NPARD(5) = 0
      NPARD(6) = 0
      NPARD(7) = 3
      NPARD(8) = 0
      NPARD(9) = 0
      NPARD(10)= 2
      NPARD(11)= 1
      NPARD(12)= 0
      NPARD(13)= 0
      NPARD(14)= 0
      NPARD(15)= 3
      NPARD(16)= 1
      NPARD(17)= 1
      NPARD(18)= 0
      NPARD(19)= 3
      NPARD(20)= 2
      NPARD(21)= 1
      NPARD(22)= 0
      NPARD(23)= 1
      NPARD(24)= 0
      NPARD(25)= 0
      
      
     
      W_MASS = 0d0
      Z_MASS = ZMASS
      H_MASS = HMASS
      T_MASS = TMASS
      DAL5H = 0d0
      V_TB = 1d0
      
    
      CALL DIZET (NPARD,W_MASS,Z_MASS,T_MASS,H_MASS,DAL5H,V_TB,ALFAS
     &           ,ALQED,ALST,ZPARD,PARTZ,PARTW)
     
     
      SIN2TW = ZPARD(3)
      ALPHST = ZPARD(15)
      DO IQCDC=0,14
      QCDCOR(IQCDC) = ZPARD(16+IQCDC)
      ENDDO
     
     
        PRINT *
        PRINT *,'DIZET input parameters:'
        PRINT *,'ZMASS',Z_MASS,'TMASS ',T_MASS
        PRINT *,'HMASS',H_MASS,'WMASS ',0d0
        PRINT *,'DAL5H',0d0,'ALQED5',ALFAI
        PRINT *,'ALFAS ',ALFAS
        PRINT *
        PRINT *,'DIZET results:'
        WMASS=ZMASS*SQRT(1D0-SIN2TW)
        PRINT *,'SIN2TW   ',SIN2TW
        PRINT *,'WMASSsin',WMASS
        PRINT *,'WMASS   ',W_MASS
        PRINT *,'DAL5H   ',DAL5H
        PRINT *,'ALQED5  ',1d0/ALQED
        PRINT *,'ALST     ',ALST
        PRINT *,'ALPAS    ',ALPHST
        PRINT *
        PRINT *,'CHANNEL         WIDTH         RHO_F_R        RHO_F_T '
     + ,'       SIN2_EFF'
        PRINT *,'-------        -------       --------       --------'
     + ,'       --------'
        DO I=0,9
          PRINT '(1X,A10,2X,F10.3,7X,F8.6,7X,F8.6,7X,F8.6)',
     +    FCHAN(I),PARTZ(I),ARROFZ(I),AROTFZ(I),ZPARD(I+5)
        ENDDO
        DO I=10,11
          PRINT '(1X,A10,2X,F10.3)',FCHAN(I),PARTZ(I)
        ENDDO
        PRINT *
        PRINT *,'W-widths'
        DO I=1,3
          PRINT '(1X,A10,2X,F10.3)',FCHANW(I),PARTW(I)
        ENDDO
        PRINT *,'****************************************************'

     
      END
