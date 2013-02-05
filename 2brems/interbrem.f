C  this part of the code may find its place in GPS class
C  the example routine is placed in FILE cellar.f just after
C  key-providing routine  GPS_brem_klu which is on top of that FILE





      SUBROUTINE GPS_HffAdd(CNorm,KFi,KFf,PX,pA,mA,pB,mB,pC,mC,pD,mD,Hel1,ph1,Hel2,ph2,mph,AmpWork) !
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   Genuine IR-finite non 1-photon amplitudes for FSR-FSR are added to AmpWork    //
*//   Photon helicity imported from the calling program.                            //
*//                                                                                 //
*//   For Y_IR=1 IR-finite 1-photon contributions are calculated here as well.      //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
*
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                                    FSR                                            //
*      //                                                                                   //
*      //                            1                  2                                   //
*      //                            |                  |                                   //
*      //             c              |                  |          d                        //
*      //    u  ------<------OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO----<----- v                 //
*      //                                     |                                             //
*      //                                     |X                                            //
*      //                                     |                                             //
*      //                                                                                   //
*      ///////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE "GPS.h"
      INTEGER            KFi,KFf
      DOUBLE PRECISION   PX(4),pA(4),pB(4),pC(4),pD(4),ph1(4),ph2(4),ph(4)
      DOUBLE PRECISION   mA,mB,mC,mD,mph
      DOUBLE COMPLEX     CNorm,sProd
      INTEGER            Hel1, Hel2
      DOUBLE COMPLEX     BornAB1D(2,2,2,2), BornAB2D(2,2,2,2), BornABCD(2,2,2,2)
      DOUBLE COMPLEX     BornABC1(2,2,2,2), BornABC2(2,2,2,2)
      DOUBLE COMPLEX     BornAB12(2,2,2,2), BornAB21(2,2,2,2)
      DOUBLE COMPLEX     AmpWork(2,2,2,2)
      DOUBLE COMPLEX     Uc11(2,2),   V11d(2,2),   U122(2,2),  V221(2,2)
      DOUBLE COMPLEX     Uc22(2,2),   V22d(2,2),   U211(2,2),  V112(2,2)
      DOUBLE COMPLEX     U21c(2,2),   Vd12(2,2),   Uc12(2,2),  V21d(2,2)
      DOUBLE COMPLEX     U12c(2,2),   Vd21(2,2),   Uc21(2,2),  V12d(2,2)
      DOUBLE COMPLEX     U121(2,2),   V121(2,2),   U212(2,2),  V212(2,2)
      DOUBLE COMPLEX     Su1,Su2,Su3
      INTEGER            j,j1,j2,j3,j4,k,Sig,l
      DOUBLE PRECISION   Fprop1, Fprop2
      DOUBLE PRECISION   prC1, prD1, prC2, prD2, prC12, prD12
      DOUBLE PRECISION   BornV_GetCharge, ChaIni,ChaFin
      DOUBLE COMPLEX     GPS_Sof1,GPS_Sof1b
      DOUBLE COMPLEX     sC(2,2),sD(2,2)
      DOUBLE COMPLEX     gF
      INTEGER            Y_IR, N_IR
      DOUBLE PRECISION   PP12(4),PP1(4),PP2(4),QQ(4),SvarQ,SvarX12,SvarX1,SvarX2
*------------------------------------------------------------------------------------------
      Y_IR=1                  ! YES, IR included
      Y_IR=0                  ! No,  IR not included
      N_IR=1-Y_IR
*--------------------
      ChaFin =  BornV_GetCharge( KFf)
      gF = DCMPLX(ChaFin *m_e_QED)
      IF( m_KeyArb  .EQ.  0 ) THEN
         sC(1,1)  =  gF *GPS_Sof1( 1,ph1,pC)
         sC(2,1)  =  gF *GPS_Sof1( 1,ph2,pC)
         sD(1,1)  = -gF *GPS_Sof1( 1,ph1,pD)
         sD(2,1)  = -gF *GPS_Sof1( 1,ph2,pD)
      ELSE
         sC(1,1)  =  gF *GPS_Sof1b( 1,ph1,pC,mC)
         sC(2,1)  =  gF *GPS_Sof1b( 1,ph2,pC,mC)
         sD(1,1)  = -gF *GPS_Sof1b( 1,ph1,pD,mD)
         sD(2,1)  = -gF *GPS_Sof1b( 1,ph2,pD,mD)
      ENDIF
      sC(1,2) = -DCONJG(sC(1,1))
      sC(2,2) = -DCONJG(sC(2,1))
      sD(1,2) = -DCONJG(sD(1,1))
      sD(2,2) = -DCONJG(sD(2,1))
* Calculate Born spin amplitudes, also with substitutions
      CALL GPS_Born(KFi,KFf,PX, pA,mA,  pB,-mB,  pC,   mC,   pD,   -mD, BornABCD) ! Standard
      CALL GPS_Born(KFi,KFf,PX, pA,mA,  pB,-mB,  ph1, mph,   pD,   -mD, BornAB1D) ! C->1
      CALL GPS_Born(KFi,KFf,PX, pA,mA,  pB,-mB,  pC,   mC,   ph1, -mph, BornABC1) ! D->1
      CALL GPS_Born(KFi,KFf,PX, pA,mA,  pB,-mB,  ph2, mph,   pD,   -mD, BornAB2D) ! C->2
      CALL GPS_Born(KFi,KFf,PX, pA,mA,  pB,-mB,  pC,   mC,   ph2, -mph, BornABC2) ! D->2
      CALL GPS_Born(KFi,KFf,PX, pA,mA,  pB,-mB,  ph1, mph,   ph2, -mph, BornAB12) ! C->1,D->2
      CALL GPS_Born(KFi,KFf,PX, pA,mA,  pB,-mB,  ph2, mph,   ph1, -mph, BornAB21) ! C->2,D->1
      DO k=1,4
         PP12(k) = pC(k)+pD(k)+ph1(k)+ph2(k)
         PP1 (k) = pC(k)+pD(k)+ph1(k)
         PP2 (k) = pC(k)+pD(k)+ph2(k)
         QQ(k)   = pC(k)+pD(k)
      ENDDO
      svarX12 = PP12(4)**2 -PP12(3)**2 -PP12(2)**2 -PP12(1)**2
      svarX1  =  PP1(4)**2  -PP1(3)**2  -PP1(2)**2  -PP1(1)**2
      svarX2  =  PP2(4)**2  -PP2(3)**2  -PP2(2)**2  -PP2(1)**2
      svarQ   =   QQ(4)**2   -QQ(3)**2   -QQ(2)**2   -QQ(1)**2
* Fermion propagarotors 1
      prC1=  1d0/(pC(4)*ph1(4)-pC(3)*ph1(3)-pC(2)*ph1(2)-pC(1)*ph1(1))/2d0
      prD1= -1d0/(pD(4)*ph1(4)-pD(3)*ph1(3)-pD(2)*ph1(2)-pD(1)*ph1(1))/2d0
* Fermion propagarotors 2
      prC2=  1d0/(pC(4)*ph2(4)-pC(3)*ph2(3)-pC(2)*ph2(2)-pC(1)*ph2(1))/2d0
      prD2= -1d0/(pD(4)*ph2(4)-pD(3)*ph2(3)-pD(2)*ph2(2)-pD(1)*ph2(1))/2d0
* Double propagators
      prC12= 1d0/( pC(4)*ph1(4)- pC(3)*ph1(3)- pC(2)*ph1(2)- pC(1)*ph1(1)
     $           + pC(4)*ph2(4)- pC(3)*ph2(3)- pC(2)*ph2(2)- pC(1)*ph2(1)
     $           +ph1(4)*ph2(4)-ph1(3)*ph2(3)-ph1(2)*ph2(2)-ph1(1)*ph2(1))/2d0
      prD12=-1d0/( pD(4)*ph1(4)- pD(3)*ph1(3)- pD(2)*ph1(2)- pD(1)*ph1(1)
     $            +pD(4)*ph2(4)- pD(3)*ph2(3)- pD(2)*ph2(2)- pD(1)*ph2(1)
     $           +ph1(4)*ph2(4)-ph1(3)*ph2(3)-ph1(2)*ph2(2)-ph1(1)*ph2(1))/2d0
      Fprop1= (1d0/prC1+1d0/prC2)*prC12 -1d0
      Fprop2= (1d0/prD1+1d0/prD2)*prD12 -1d0
      Sig = 3-2*Hel1
      IF( m_KeyArb  .EQ.  0 ) THEN
* end of line
         CALL GPS_MatrV( gF, ph1,Sig,  ph1,mph, pD,mD,     V11d) ! <1|{1}|D>
         CALL GPS_MatrU( gF, ph1,Sig,  pC,mC,   ph1,mph,   Uc11) ! <C|[1]|1>
* false second
         CALL GPS_MatrV( gF, ph1,Sig,  ph2,mph, pD,mD,     V21d) ! <2|{1}|D>
         CALL GPS_MatrU( gF, ph1,Sig,  pC,mC,   ph2,mph,   Uc12) ! <C|[1]|2>
* reverse order
         CALL GPS_MatrV( gF, ph1,Sig,  pD,mD,   ph2,mph,   Vd12) ! <D|{1}|2>
         CALL GPS_MatrU( gF, ph1,Sig,  ph2,mph, pC,mC,     U21c) ! <2|[1]|C>
* xk-xk term case, ph2 first
         CALL GPS_MatrV( gF, ph1,Sig,  ph1,mph, ph2,mph,   V112) ! <1|{1}|2>
         CALL GPS_MatrU( gF, ph1,Sig,  ph2,mph, ph1,mph,   U211) ! <2|[1]|1>
* xk-xk term case ph2 first
         CALL GPS_MatrV( gF, ph1,Sig,  ph2,mph, ph2,mph,   V212) ! <2|{1}|2>
         CALL GPS_MatrU( gF, ph1,Sig,  ph2,mph, ph2,mph,   U212) ! <2|[1]|2>
      ELSE
         CALL GPS_MatrVb(gF, ph1,Sig,  ph1,mph, pD,mD,     V11d)
         CALL GPS_MatrUb(gF, ph1,Sig,  pC,mC,   ph1,mph,   Uc11)
* falSe second
         CALL GPS_MatrVb(gF, ph1,Sig,  ph2,mph, pD,mD,     V21d)
         CALL GPS_MatrUb(gF, ph1,Sig,  pC,mC,   ph2,mph,   Uc12)
* reverse order
         CALL GPS_MatrVb(gF, ph1,Sig,  pD,mD,   ph2,mph,   Vd12)
         CALL GPS_MatrUb(gF, ph1,Sig,  ph2,mph, pC,mC,     U21c)
* for the case when there was ph2 first xk-xk term
         CALL GPS_MatrVb(gF, ph1,Sig,  ph1,mph, ph2,mph,   V112)
         CALL GPS_MatrUb(gF, ph1,Sig,  ph2,mph, ph1,mph,   U211)
* for the case when there was ph2 first xk-xk term
         CALL GPS_MatrVb(gF, ph1,Sig,  ph2,mph, ph2,mph,   V212)
         CALL GPS_MatrUb(gF, ph1,Sig,  ph2,mph, ph2,mph,   U212)
      ENDIF
      Sig = 3-2*Hel2
      IF( m_KeyArb  .EQ.  0 ) THEN
         CALL GPS_MatrV( gF, ph2,Sig,  ph2,mph,  pD,mD,    V22d) ! <2|{2}|D>
         CALL GPS_MatrU( gF, ph2,Sig,  pC,mC,   ph2,mph,   Uc22) ! <C|[2]|2>
* falSe second
         CALL GPS_MatrV( gF, ph2,Sig,  ph1,mph,  pD,mD,    V12d) ! <1|{2}|D>
         CALL GPS_MatrU( gF, ph2,Sig,  pC,mC,   ph1,mph,   Uc21) ! <C|[2]|1>
* reverse order
         CALL GPS_MatrV( gF, ph2,Sig,  pD,mD,   ph1,mph,   Vd21) ! <D|{2}|1>
         CALL GPS_MatrU( gF, ph2,Sig,  ph1,mph, pC,mC,     U12c) ! <1|[2]|C>
* xk-xk term, ph1 first
         CALL GPS_MatrV( gF, ph2,Sig,  ph2,mph, ph1,mph,   V221) ! <2|{2}|1>
         CALL GPS_MatrU( gF, ph2,Sig,  ph1,mph, ph2,mph,   U122) ! <1|[2]|2>
* xk-xk term, ph1 first 
         CALL GPS_MatrV( gF, ph2,Sig,  ph1,mph, ph1,mph,   V121) ! <1|{2}|1>
         CALL GPS_MatrU( gF, ph2,Sig,  ph1,mph, ph1,mph,   U121) ! <1|[2]|1>
      ELSE
         CALL GPS_MatrVb(gF, ph2,Sig,  ph2,mph,  pD,mD,    V22d)
         CALL GPS_MatrUb(gF, ph2,Sig,  pC,mC,   ph2,mph,   Uc22)
* falSe second
         CALL GPS_MatrVb(gF, ph2,Sig,  ph1,mph,  pD,mD,    V12d)
         CALL GPS_MatrUb(gF, ph2,Sig,  pC,mC,   ph1,mph,   Uc21)
* reverse order
         CALL GPS_MatrVb(gF, ph2,Sig,  pD,mD,   ph1,mph,   Vd21)
         CALL GPS_MatrUb(gF, ph2,Sig,  ph1,mph, pC,mC,     U12c)
* for the case when there was ph1 first xk-xk term
         CALL GPS_MatrVb(gF, ph2,Sig,  ph2,mph, ph1,mph,   V221)
         CALL GPS_MatrUb(gF, ph2,Sig,  ph1,mph, ph2,mph,   U122)
* for the case when there was ph1 first xk-xk term
         CALL GPS_MatrVb(gF, ph2,Sig,  ph1,mph, ph1,mph,   V121)
         CALL GPS_MatrUb(gF, ph2,Sig,  ph1,mph, ph1,mph,   U121)
      ENDIF
      DO j1=1,2
         DO j2=1,2
            DO j3=1,2
               DO j4=1,2
                  Su1 = DCMPLX(0d0,0d0)
                  DO j=1,2
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                     1|            2|                                              //
*      //               c      |    c+m+1    |    c+m+1+2          -d                       //
*      //       u  -----<------S-----<-------U------<-------O-------<------ v               //
*      //                                                   |X                              //
*      ///////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +(sC(1,Hel1)*(prC12-prC2*N_IR))*Uc22(j3,j)*BornAB2D(j1,j2,j,j4) !<c|(1)c[2]2|X|d>
                     Su1=Su1  +sC(1,Hel1)* prC12            *Uc21(j3,j)*BornAB1D(j1,j2,j,j4) !<c|(1)c[2]1|X|d>
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                     2|            1|                                              //
*      //               c      |    c+m+2    |    c+m+1+2          -d                       //
*      //       u  -----<------S-----<-------U------<-------O-------<------ v               //
*      //                                                   |X                              //
*      ///////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1+(sC(2,Hel2)*(prC12-prC1*N_IR))*Uc11(j3,j)*BornAB1D(j1,j2,j,j4) !<c|(2)c[1]1|X|d>
                     Su1=Su1 +sC(2,Hel2)* prC12            *Uc12(j3,j)*BornAB2D(j1,j2,j,j4) !<c|(2)c[1]2|X|d>
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                                    |1             |2                              //
*      //               c          -d+m-1-2  |   -d+m-2     |       -d                      //
*      //       u  -----<------O-----<-------V-----<--------S--------<----- v               //
*      //                     X|                                                            //
*      ///////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +BornABC1(j1,j2,j3,j)*V11d(j,j4)*( sD(2,Hel2)*(prD12-prD1*N_IR))!<c|X|1{1}d(2)|d>
                     Su1=Su1 +BornABC2(j1,j2,j3,j)*V21d(j,j4)*  sD(2,Hel2)* prD12            !<c|X|2{1}d(2)|d>
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                                    |2             |1                              //
*      //               c          -d+m-1-2  |   -d+m-1     |       -d                      //
*      //       u  -----<------O-----<-------V-----<--------S--------<----- v               //
*      //                     X|                                                            //
*      ///////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +BornABC2(j1,j2,j3,j)*V22d(j,j4)*( sD(1,Hel1)*(prD12-prD2*N_IR))!<c|X|2{2}d(1)|d>
                     Su1=Su1 +BornABC1(j1,j2,j3,j)*V12d(j,j4)  *sD(1,Hel1)* prD12            !<c|X|1{2}d(1)|d>
                  ENDDO
                  Su3 = DCMPLX(0d0,0d0)
                  DO j=1,2
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                     2|                            |1                              //
*      //               c      |    c+m+1         c+m+1     |      -d                       //
*      //       u  -----<------U-----<-------O------<-------S-------<------ v               //
*      //                                    |X                                             //
*      ///////////////////////////////////////////////////////////////////////////////////////
                     Su3=Su3 +Uc22(j3,j)*prC2  *BornAB2D(j1,j2,j,j4) *sD(1,Hel1) *Y_IR !<c|[2]c|X|1(1)|d>
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                     1|                            |2                              //
*      //               c      |    c+m+1         c+m+2     |      -d                       //
*      //       u  -----<------U-----<-------O------<-------S-------<------ v               //
*      //                                    |X                                             //
*      ///////////////////////////////////////////////////////////////////////////////////////
                     Su3=Su3 +Uc11(j3,j)*prC1 *BornAB1D(j1,j2,j,j4)  *sD(2,Hel2) *Y_IR !<c|[1]c|X|2(2)|d>
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                     2|                            |1                              //
*      //               c      |    c+m+2         c+m+1     |      -d                       //
*      //       u  -----<------S-----<-------O------<-------U-------<------ v               //
*      //                                    |X                                             //
*      ///////////////////////////////////////////////////////////////////////////////////////
                     Su3=Su3 +sC(2,Hel2) *BornABC1(j1,j2,j3,j) *prD1 *V11d(j,j4) *Y_IR !<c|(2)c|X|1{1}|d>
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                     1|                            |2                              //
*      //               c      |    c+m+1         c+m+2     |      -d                       //
*      //       u  -----<------S-----<-------O------<-------U-------<------ v               //
*      //                                    |X                                             //
*      ///////////////////////////////////////////////////////////////////////////////////////
                     Su3=Su3 +sC(1,Hel1) *BornABC2(j1,j2,j3,j) *prD2 *V22d(j,j4) *Y_IR !<c|(1)c|X|2{2}|d>
                  ENDDO
                  Su2 = DCMPLX(0d0,0d0)
                  DO j=1,2
                     DO l=1,2
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                      |2                           |1                              //
*      //               c      |    c+m+2        -d+m-1     |       -d                      //
*      //       u  -----<------U-----<-------O-----<--------V--------<----- v               //
*      //                                    |X                                             //
*      ///////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Uc22( j3,l)*prC2 *BornAB21(j1,j2,l,j ) *V11d( j,j4)*prD1 !<c|[2]2|X|1{1}|d>
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                      |1                           |2                              //
*      //               c      |    c+m+1        -d+m-2     |       -d                      //
*      //       u  -----<------U-----<-------O-----<--------V--------<----- v               //
*      //                                    |X                                             //
*      ///////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Uc11( j3,l)*prC1 *BornAB12(j1,j2,l,j ) *V22d( j,j4)*prD2 !<c|[1]1|X|2{2}|d>
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                     1|            2|                                              //
*      //               c      |    c+m+1    |    c+m+1+2          -d                       //
*      //       u  -----<------U-----<-------O------<-------V-------<------ v               //
*      //                                                  X|                               //
*      ///////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Uc11( j3,l)*prC1  *U122(l,j)*prC12  *BornAB2D(j1,j2,j ,j4) ! <c|[1]1[2]2|X|d>
                        Su2=Su2 +Uc11( j3,l)*prC1  *U121(l,j)*prC12  *BornAB1D(j1,j2,j ,j4) ! <c|[1]1[2]1|X|d>
                        Su2=Su2 +Uc11( j3,l)*prC1  *U12c(l,j)*prC12  *BornABCD(j1,j2,j ,j4) ! <c|[1]1[2]c|X|d>
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                     2|            1|                                              //
*      //               c      |    c+m+2    |    c+m+1+2          -d                       //
*      //       u  -----<------U-----<-------U------<-------O-------<------ v               //
*      //                                                  X|                               //
*      ///////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Uc22( j3,l)*prC2  *U211(l,j)*prC12  *BornAB1D(j1,j2,j ,j4) ! <c|[2]2[1]1|X|d>
                        Su2=Su2 +Uc22( j3,l)*prC2  *U212(l,j)*prC12  *BornAB2D(j1,j2,j ,j4) ! <c|[2]2[1]2|X|d>
                        Su2=Su2 +Uc22( j3,l)*prC2  *U21c(l,j)*prC12  *BornABCD(j1,j2,j ,j4) ! <c|[2]2[1]c|X|d>
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                                    |2             |1                              //
*      //               c          -d+m-1-2  |   -d+m-1     |       -d                      //
*      //       u  -----<------O-----<-------V-----<--------V--------<----- v               //
*      //                      |X                                                           //
*      ///////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +BornABC2(j1,j2,j3,j ) *V221(j,l )*prD12 *V11d( l,j4)*prD1  ! <c|X|2{2}1{1}|d>
                        Su2=Su2 +BornABC1(j1,j2,j3,j ) *V121(j,l )*prD12 *V11d( l,j4)*prD1  ! <c|X|1{2}1{1}|d>
                        Su2=Su2 +BornABCD(j1,j2,j3,j ) *Vd21(j,l )*prD12 *V11d( l,j4)*prD1  ! <c|X|d{2}1{1}|d>
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                                    |1             |2                              //
*      //               c          -d+m-1-2  |   -d+m-2     |       -d                      //
*      //       u  -----<------O-----<-------V-----<--------V--------<----- v               //
*      //                      |X                                                           //
*      ///////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +BornABC1(j1,j2,j3,j ) *V112(j,l )*prD12 *V22d( l,j4)*prD2  ! <c|X|1{1}2{2}|d>
                        Su2=Su2 +BornABC2(j1,j2,j3,j ) *V212(j,l )*prD12 *V22d( l,j4)*prD2  ! <c|X|2{1}2{2}|d>
                        Su2=Su2 +BornABCD(j1,j2,j3,j ) *Vd12(j,l )*prD12 *V22d( l,j4)*prD2  ! <c|X|d{1}2{2}|d>
                     ENDDO
                  ENDDO
                  sProd = (sC(1,Hel1)+sD(1,Hel1)) *( sC(2,Hel2)+sD(2,Hel2)) !
                  AmpWork(j1,j2,j3,j4) = AmpWork(j1,j2,j3,j4)
     $                 +CNorm*( Su1 +Su2 +Su3)
     $                 +CNorm*BornABCD(j1,j2,j3,j4)*( sC(1,Hel1)*sC(2,Hel2)*Fprop1   !
     $                                               +sD(1,Hel1)*sD(2,Hel2)*Fprop2 ) !
     $                 +CNorm*BornABCD(j1,j2,j3,j4)*sProd *(1d0 -svarX12/svarQ)  *N_IR !
     $                 +CNorm*BornABCD(j1,j2,j3,j4)*sProd *( svarX1/svarQ -1d0)  *N_IR !
     $                 +CNorm*BornABCD(j1,j2,j3,j4)*sProd *( svarX2/svarQ -1d0)  *N_IR !
     $                 +CNorm*BornABCD(j1,j2,j3,j4)*sProd                        *Y_IR !
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      END                       ! GPS_HffAdd



      SUBROUTINE GPS_HiiAdd(CNorm,KFi,KFf,PX,pA,mA,pB,mB,pC,mC,pD,mD,Hel1,ph1,Hel2,ph2,mph,AmpWork) !
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   Genuine IR-finite non 1-photon amplitudes for ISR-ISR are added to AmpWork    //
*//   That is for dip-switch Y_IR=0.                                                //
*//   For Y_IR=1 IR-finite 1-photon contributions are calculated here as well       //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
*//                                        |                                        //
*//                              1         |          2                             //
*//                              |         |X         |                             //
*//      _       -b              |         |          |          a                  //
*//      v  ------<------OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO----<----- u           //
*//                                                                                 //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE "GPS.h"
      INTEGER           KFi,KFf
      DOUBLE PRECISION  PX(4),pA(4),pB(4),pC(4),pD(4),ph1(4),ph2(4),ph(4)
      DOUBLE PRECISION  mA,mB,mC,mD,mph
      DOUBLE COMPLEX    CNorm
      INTEGER           Hel1, Hel2
      DOUBLE COMPLEX    Born1BCD(2,2,2,2), Born2BCD(2,2,2,2), BornABCD(2,2,2,2)
      DOUBLE COMPLEX    BornA1CD(2,2,2,2), BornA2CD(2,2,2,2)
      DOUBLE COMPLEX    Born12CD(2,2,2,2), Born21CD(2,2,2,2)
      DOUBLE COMPLEX    AmpWork(2,2,2,2)
      DOUBLE COMPLEX    U11a(2,2),Vb11(2,2),U221(2,2),V122(2,2)
      DOUBLE COMPLEX    U22a(2,2),Vb22(2,2),U112(2,2),V211(2,2)
      DOUBLE COMPLEX    Ua12(2,2),V21b(2,2),U21a(2,2),Vb12(2,2)
      DOUBLE COMPLEX    Ua21(2,2),V12b(2,2),U12a(2,2),Vb21(2,2)
      DOUBLE COMPLEX    U121(2,2),V121(2,2),U212(2,2),V212(2,2)
      DOUBLE COMPLEX    Su1,Su2,sProd
      INTEGER           j,j1,j2,j3,j4,k,Sig,l
      DOUBLE PRECISION  prA1,prB1,prA2,prB2,prA12,prB12
      DOUBLE PRECISION  BornV_GetCharge,ChaIni
      DOUBLE PRECISION  Fprop1,Fprop2
      DOUBLE COMPLEX    gI
      DOUBLE COMPLEX    GPS_Sof1,GPS_Sof1b
      DOUBLE COMPLEX    sA(2,2),sB(2,2)
      INTEGER           Y_IR,N_IR
*------------------------------------------------------------------------------------------
      Y_IR=1                  ! YES, IR included
      Y_IR=0                  ! No,  IR not included
      N_IR=1-Y_IR
*--------------------
      ChaIni =  BornV_GetCharge( KFi)
      gI = DCMPLX(ChaIni*m_e_QED)
      IF( m_KeyArb  .EQ.  0 ) THEN
         sA(1,1)  = -gI*GPS_Sof1( 1,ph1,pA)
         sA(2,1)  = -gI*GPS_Sof1( 1,ph2,pA)
         sB(1,1)  =  gI*GPS_Sof1( 1,ph1,pB)
         sB(2,1)  =  gI*GPS_Sof1( 1,ph2,pB)
      ELSE
         sA(1,1)  = -gI*GPS_Sof1b( 1,ph1,pA,mA)
         sA(2,1)  = -gI*GPS_Sof1b( 1,ph2,pA,mA)
         sB(1,1)  =  gI*GPS_Sof1b( 1,ph1,pB,mB)
         sB(2,1)  =  gI*GPS_Sof1b( 1,ph2,pB,mB)
      ENDIF
      sA(1,2) = -DCONJG(sA(1,1))
      sA(2,2) = -DCONJG(sA(2,1))
      sB(1,2) = -DCONJG(sB(1,1))
      sB(2,2) = -DCONJG(sB(2,1))
* Calculate Born spin amplitudes
      CALL GPS_Born(KFi,KFf,PX, pA,mA,    pB,-mB,      pC,MC,   pD,-mD,   BornABCD) ! Standard
      CALL GPS_Born(KFi,KFf,PX, ph1,mph,  pB,-mB,      pC,mC,   pD,-mD,   Born1BCD) ! A->1
      CALL GPS_Born(KFi,KFf,PX, pA,mA,    ph1,-mph,    pC,mC,   pD,-mD,   BornA1CD) ! B->1
      CALL GPS_Born(KFi,KFf,PX, ph2,mph,  pB,-mB,      pC,mC,   pD,-mD,   Born2BCD) ! A->2
      CALL GPS_Born(KFi,KFf,PX, pA,mA,    ph2,-mph,    pC,mC,   pD,-mD,   BornA2CD) ! B->2
      CALL GPS_Born(KFi,KFf,PX, ph1,mph,  ph2,-mph,    pC,mC,   pD,-mD,   Born12CD) ! A->1,B->2
      CALL GPS_Born(KFi,KFf,PX, ph2,mph,  ph1,-mph,    pC,mC,   pD,-mD,   Born21CD) ! A->2,B->1
* Fermion propagarotors ini1
      prA1= 1d0/(pA(4)*ph1(4)-pA(3)*ph1(3)-pA(2)*ph1(2)-pA(1)*ph1(1))/2d0
      prB1=-1d0/(pB(4)*ph1(4)-pB(3)*ph1(3)-pB(2)*ph1(2)-pB(1)*ph1(1))/2d0
* Fermion propagarotors ini2
      prA2= 1d0/(pA(4)*ph2(4)-pA(3)*ph2(3)-pA(2)*ph2(2)-pA(1)*ph2(1))/2d0
      prB2=-1d0/(pB(4)*ph2(4)-pB(3)*ph2(3)-pB(2)*ph2(2)-pB(1)*ph2(1))/2d0
* DOUBLE propagators
      prA12= 1d0/( pA(4)*ph1(4)- pA(3)*ph1(3)- pA(2)*ph1(2)- pA(1)*ph1(1)
     $           + pA(4)*ph2(4)- pA(3)*ph2(3)- pA(2)*ph2(2)- pA(1)*ph2(1)
     $           -ph1(4)*ph2(4)+ph1(3)*ph2(3)+ph1(2)*ph2(2)+ph1(1)*ph2(1)
     $           )/2d0
      prB12=-1d0/( pB(4)*ph1(4)- pB(3)*ph1(3)- pB(2)*ph1(2)- pB(1)*ph1(1)
     $            +pB(4)*ph2(4)- pB(3)*ph2(3)- pB(2)*ph2(2)- pB(1)*ph2(1)
     $           -ph1(4)*ph2(4)+ph1(3)*ph2(3)+ph1(2)*ph2(2)+ph1(1)*ph2(1)
     $           )/2d0
      Fprop1=(1d0/prA1+1d0/prA2)*prA12-1d0
      Fprop2=(1d0/prB1+1d0/prB2)*prB12-1d0
      Sig = 3-2*Hel1
      IF( m_KeyArb  .EQ.  0 ) THEN
         CALL GPS_MatrU( gI, ph1,Sig,  ph1,mph, pA,mA,     U11a) ! <1|[1]|a>
         CALL GPS_MatrV( gI, ph1,Sig,  pB,mB,   ph1,mph,   Vb11) ! <b|{1}|1>
* falSe second
         CALL GPS_MatrU( gI, ph1,Sig,  ph2,mph, pA,mA,     U21a) ! <2|[1]|a>
         CALL GPS_MatrV( gI, ph1,Sig,  pB,mB,   ph2,mph,   Vb12) ! <b|{1}|2>
* reverse order
         CALL GPS_MatrU( gI, ph1,Sig,  pA,mA,   ph2,mph,   Ua12) ! <a|[1]|2>
         CALL GPS_MatrV( gI, ph1,Sig,  ph2,mph, pB,mB,     V21b) ! <2|{1}|b>
* for the case when there was ph2 first xk-xk term
         CALL GPS_MatrU( gI, ph1,Sig,  ph1,mph, ph2,mph,   U112) ! <1|[1]|2>
         CALL GPS_MatrV( gI, ph1,Sig,  ph2,mph, ph1,mph,   V211) ! <2|{1}|1>
* for the case when there was ph2 first xk-xk term
         CALL GPS_MatrU( gI, ph1,Sig,  ph2,mph, ph2,mph,   U212) ! <2|[1]|2>
         CALL GPS_MatrV( gI, ph1,Sig,  ph2,mph, ph2,mph,   V212) ! <2|{1}|2>
      ELSE
         CALL GPS_MatrUb(gI, ph1,Sig,  ph1,mph, pA,mA,     U11a)
         CALL GPS_MatrVb(gI, ph1,Sig,  pB,mB,   ph1,mph,   Vb11)
* falSe second
         CALL GPS_MatrUb(gI, ph1,Sig,  ph2,mph, pA,mA,     U21a)
         CALL GPS_MatrVb(gI, ph1,Sig,  pB,mB,   ph2,mph,   Vb12)
* reverse order
         CALL GPS_MatrUb(gI, ph1,Sig,  pA,mA,   ph2,mph,   Ua12)
         CALL GPS_MatrVb(gI, ph1,Sig,  ph2,mph, pB,mB,     V21b)
* for the case when there was ph2 first xk-xk term
         CALL GPS_MatrUb(gI, ph1,Sig,  ph1,mph, ph2,mph,   U112)
         CALL GPS_MatrVb(gI, ph1,Sig,  ph2,mph, ph1,mph,   V211)
* for the case when there was ph2 first xk-xk term
         CALL GPS_MatrUb(gI, ph1,Sig,  ph2,mph, ph2,mph,   U212)
         CALL GPS_MatrVb(gI, ph1,Sig,  ph2,mph, ph2,mph,   V212)
      ENDIF
      Sig = 3-2*Hel2
      IF( m_KeyArb  .EQ.  0 ) THEN
         CALL GPS_MatrU( gI, ph2,Sig,  ph2,mph,  pA,mA,    U22a) ! <2|[2]|a>
         CALL GPS_MatrV( gI, ph2,Sig,  pB,mB,   ph2,mph,   Vb22) ! <b|{2}|2>
* falSe second
         CALL GPS_MatrU( gI, ph2,Sig,  ph1,mph,  pA,mA,    U12a) ! <1|[2]|a>
         CALL GPS_MatrV( gI, ph2,Sig,  pB,mB,   ph1,mph,   Vb21) ! <b|{2}|1>
* reverse order
         CALL GPS_MatrU( gI, ph2,Sig,  pA,mA,   ph1,mph,   Ua21) ! <a|[2]|1>
         CALL GPS_MatrV( gI, ph2,Sig,  ph1,mph, pB,mB,     V12b) ! <1|{2}|b>
* for the case when there was ph1 first xk-xk term
         CALL GPS_MatrU( gI, ph2,Sig,  ph2,mph, ph1,mph,   U221) ! <2|[2]|1>
         CALL GPS_MatrV( gI, ph2,Sig,  ph1,mph, ph2,mph,   V122) ! <1|{2}|2>
c for the case when there was ph1 first xk-xk term
         CALL GPS_MatrU( gI, ph2,Sig,  ph1,mph, ph1,mph,   U121) ! <1|[2]|1>
         CALL GPS_MatrV( gI, ph2,Sig,  ph1,mph, ph1,mph,   V121) ! <1|{2}|1>
      ELSE
         CALL GPS_MatrUb(gI, ph2,Sig,  ph2,mph,  pA,mA,    U22a)
         CALL GPS_MatrVb(gI, ph2,Sig,  pB,mB,   ph2,mph,   Vb22)
c falSe second
         CALL GPS_MatrUb(gI, ph2,Sig,  ph1,mph,  pA,mA,    U12a)
         CALL GPS_MatrVb(gI, ph2,Sig,  pB,mB,   ph1,mph,   Vb21)
c reverse order
         CALL GPS_MatrUb(gI, ph2,Sig,  pA,mA,   ph1,mph,   Ua21)
         CALL GPS_MatrVb(gI, ph2,Sig,  ph1,mph, pB,mB,     V12b)
c for the case when there was ph1 first xk-xk term
         CALL GPS_MatrUb(gI, ph2,Sig,  ph2,mph, ph1,mph,   U221)
         CALL GPS_MatrVb(gI, ph2,Sig,  ph1,mph, ph2,mph,   V122)
c for the case when there was ph1 first xk-xk term
         CALL GPS_MatrUb(gI, ph2,Sig,  ph1,mph, ph1,mph,   U121)
         CALL GPS_MatrVb(gI, ph2,Sig,  ph1,mph, ph1,mph,   V121)
      ENDIF
      DO j1=1,2
         DO j2=1,2
            DO j3=1,2
               DO j4=1,2
                  Su1 = DCMPLX(0d0,0d0)
                  DO j=1,2
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    |                  1               2                         //
*        //                    |X                 |               |                         //
*        //      _       -b    |      b+m-1-2     |      a+m-2    |      a                  //
*        //      v  ------<----O--------<---------U--------<------S------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1+Born1BCD(j,j2,j3,j4) *U11a(j,j1)*(prA12-prA1*N_IR)*sA(2,Hel2) !<b|X|1[1]a(2)|a>
                     Su1=Su1+Born2BCD(j,j2,j3,j4) *U21a(j,j1)*   prA12         *sA(2,Hel2) !<b|X|2[1]a(2)|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    |                  2               1                         //
*        //                    |X                 |               |                         //
*        //      _       -b    |      b+m-1-2     |      a+m-1    |      a                  //
*        //      v  ------<----O--------<---------U-------<-------S------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1+Born2BCD(j,j2,j3,j4) *U22a(j,j1)*(prA12-prA2*N_IR)*sA(1,Hel1) !<b|X|2[2]a(1)|a>
                     Su1=Su1+Born1BCD(j,j2,j3,j4) *U12a(j,j1)* prA12           *sA(1,Hel1) !<b|X|1[2]a(1)|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    2                  1               |                         //
*        //                    |                  |               |X                        //
*        //      _       -b    |     -b+m+2       |    -b+m+1+2   |      a                  //
*        //      v  ------<----S--------<---------V--------<------O------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1+( sB(2,Hel2)*(-prB1*N_IR+prB12))*Vb11(j2,j)*BornA1CD(j1,j,j3,j4)!<b|(2)b[1]1|X|a>
                     Su1=Su1+( sB(2,Hel2))           *prB12  *Vb12(j2,j)*BornA2CD(j1,j,j3,j4)!<b|(2)b[1]2|X|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    1                  2               |                         //
*        //                    |                  |               |X                        //
*        //      _       -b    |     -b+m+1       |    -b+m+1+2   |      a                  //
*        //      v  ------<----S--------<---------V--------<------O------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +(sB(1,Hel1)*(-prB2*N_IR+prB12))*Vb22(j2,j)*BornA2CD(j1,j,j3,j4)!<b|(1)b[2]2|X|a>
                     Su1=Su1 +(sB(1,Hel1))           *prB12  *Vb21(j2,j)*BornA1CD(j1,j,j3,j4)!<b|(1)b[2]1|X|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    2                  |               1                         //
*        //                    |                  |X              |                         //
*        //      _       -b    |     -b+m+2       |     a+m-1     |      a                  //
*        //      v  ------<----S--------<---------O--------<------U------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +sB(2,Hel2)    *Born1BCD(j,j2,j3,j4) *prA1*U11a(j,j1) *Y_IR !<b|(2)2|X|1[1]|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    1                  |               2                         //
*        //                    |                  |X              |                         //
*        //      _       -b    |     -b+m+1       |     a+m-2     |      a                  //
*        //      v  ------<----S--------<---------O--------<------U------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +sB(1,Hel1)    *Born2BCD(j,j2,j3,j4) *prA2*U22a(j,j1) *Y_IR !<b|(1)1|X|2[2]|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    1                  |               2                         //
*        //                    |                  |X              |                         //
*        //      _       -b    |     -b+m+1       |     a+m-2     |      a                  //
*        //      v  ------<----V--------<---------O--------<------S------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +Vb11(j2,j)*prB1 *BornA1CD(j1,j,j3,j4) *sA(2,Hel2)    *Y_IR !<b|[1]1|X|a[2]|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    2                  |               1                         //
*        //                    |                  |X              |                         //
*        //      _       -b    |     -b+m+2       |     a+m-1     |      a                  //
*        //      v  ------<----V--------<---------O--------<------S------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +Vb22(j2,j)*prB2 *BornA2CD(j1,j,j3,j4) *sA(1,Hel1)    *Y_IR !<b|[2]2|X|a(1)|a>
                  ENDDO
                  Su2 = DCMPLX(0d0,0d0)
                  DO j=1,2
                     DO l=1,2
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    2                  |               1                         //
*        //                    |                  |X              |                         //
*        //      _       -b    |     -b+m+2       |     a+m-1     |      a                  //
*        //      v  ------<----*--------<---------O--------<------O------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2  +Vb22(j2,l)*prB2  *Born12CD(j,l,j3,j4 )  *U11a(j,j1)*prA1 ! <b|[2]2|X|1[1]|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    1                  |               2                         //
*        //                    |                  |X              |                         //
*        //      _       -b    |     -b+m+1       |     a+m-2     |      a                  //
*        //      v  ------<----*--------<---------O--------<------O------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2  +Vb11(j2,l)*prB1  *Born21CD(j,l,j3,j4 )  *U22a(j,j1)*prA2 ! <b|[1]1|X|2[2]|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    |                  2               1                         //
*        //                    |X                 |               |                         //
*        //      _       -b    |      b+m-1-2     |      a+m-1    |      a                  //
*        //      v  ------<----O--------<---------O--------<------*------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Born1BCD(j,j2,j3,j4) *U121(j,l)*prA12 *U11a(l,j1)*prA1  ! <b|X|1[2]1(1)|a>
                        Su2=Su2 +Born2BCD(j,j2,j3,j4) *U221(j,l)*prA12 *U11a(l,j1)*prA1  ! <b|X|2[2]1(1)|a>
                        Su2=Su2 -BornABCD(j,j2,j3,j4) *Ua21(j,l)*prA12 *U11a(l,j1)*prA1  ! <b|X|a[2]1(1)|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    |                  1               2                         //
*        //                    |X                 |               |                         //
*        //      _       -b    |      b+m-1-2     |      a+m-2    |      a                  //
*        //      v  ------<----O--------<---------O--------<------*------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Born2BCD(j,j2,j3,j4) *U212(j,l)*prA12  *U22a(l,j1)*prA2 ! <b|X|2[1]2(2)|a>
                        Su2=Su2 +Born1BCD(j,j2,j3,j4) *U112(j,l)*prA12  *U22a(l,j1)*prA2 ! <b|X|1[1]2(2)|a>
                        Su2=Su2 -BornABCD(j,j2,j3,j4) *Ua12(j,l)*prA12  *U22a(l,j1)*prA2 ! <b|X|a[1]2(2)|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    1                  2               |                         //
*        //                    |                  |               |X                        //
*        //      _       -b    |     -b+m+1       |    -b+m+1+2   |      a                  //
*        //      v  ------<----*--------<---------O--------<------O------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Vb11(j2,l)*prB1 *V121(l,j)*prB12 *BornA1CD(j1,j,j3,j4) ! <b|(1)1[2]1|X|a>
                        Su2=Su2 +Vb11(j2,l)*prB1 *V122(l,j)*prB12 *BornA2CD(j1,j,j3,j4) ! <b|(1)1[2]2|X|a>
                        Su2=Su2 -Vb11(j2,l)*prB1 *V12b(l,j)*prB12 *BornABCD(j1,j,j3,j4) ! <b|(1)1[2]b|X|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    2                  1               |                         //
*        //                    |                  |               |X                        //
*        //      _       -b    |     -b+m+2       |    -b+m+1+2   |      a                  //
*        //      v  ------<----*--------<---------O--------<------O------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Vb22(j2,l)*prB2 *V212(l,j)*prB12  *BornA2CD(j1,j,j3,j4) ! <b|(2)2[1]2|X|a>
                        Su2=Su2 +Vb22(j2,l)*prB2 *V211(l,j)*prB12  *BornA1CD(j1,j,j3,j4) ! <b|(2)2[1]1|X|a>
                        Su2=Su2 -Vb22(j2,l)*prB2 *V21b(l,j)*prB12  *BornABCD(j1,j,j3,j4) ! <b|(2)2[2]b|X|a>
                     ENDDO
                  ENDDO
                  sProd = ( sA(1,Hel1)+sB(1,Hel1))*( sA(2,Hel2)+sB(2,Hel2))
                  AmpWork(j1,j2,j3,j4) =AmpWork(j1,j2,j3,j4)
     $                 +CNorm*( Su1+Su2 )
     $                 +CNorm*BornABCD(j1,j2,j3,j4)*( sA(1,Hel1)*sA(2,Hel2)*Fprop1  !
     $                                               +sB(1,Hel1)*sB(2,Hel2)*Fprop2) !
     $                 +CNorm*BornABCD(j1,j2,j3,j4)* sProd                   *Y_IR  !
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      END                       ! GPS_HiiAdd    

      SUBROUTINE GPS_HifAdd(CNorm,KFi,KFf,PX,pA,mA,pB,mB,pC,mC,pD,mD,Hel1,ph1,Hel2,ph2,mph,AmpWork) !
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   Genuine IR-finite non 1-photon amplitudes for ISR-FSR are added to AmpWork    //
*//                                                                                 //
*//   For Y_IR=1 IR-finite 1-photon contributions are calculated here as well       //
*//                                                                                 //
*//   1-st photon in Initials state,  symmetrisation 1<-->2 required                //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//                              1                2                                 //
*//                              |                |                                 //
*//                              |                |                                 //
*//               c              |    OOOOOOOOOOOOOOOOOO          d                 //
*//     u  -------<------------- | ---OOOOOOOOOOOOOOOOOO----------<----- v          //
*//                              |         |                                        //
*//                              |         |X                                       //
*//                              |         |                                        //
*//      _       -b          OOOOOOOOOOOOOOOOOOOOO                a                 //
*//      v  ------<----------OOOOOOOOOOOOOOOOOOOOO----------------<----- u          //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE "GPS.h"
      INTEGER           KFi,KFf
      DOUBLE PRECISION  PX(4),pA(4),pB(4),pC(4),pD(4),ph1(4),ph2(4),ph(4)
      DOUBLE PRECISION  mA,mB,mC,mD,mph
      DOUBLE COMPLEX    CNorm
      INTEGER           Hel1,Hel2
      DOUBLE COMPLEX    Born1BCD(2,2,2,2),BornAB2D(2,2,2,2),BornABCD(2,2,2,2)
      DOUBLE COMPLEX    BornA1CD(2,2,2,2),BornABC2(2,2,2,2)
      DOUBLE COMPLEX    Born1B2D(2,2,2,2),Born1BC2(2,2,2,2)
      DOUBLE COMPLEX    BornA12D(2,2,2,2),BornA1C2(2,2,2,2)
      DOUBLE COMPLEX    AmpWork(2,2,2,2)
      DOUBLE COMPLEX    U11a(2,2),Vb11(2,2),Uc22(2,2),V22d(2,2)
      DOUBLE COMPLEX    Su1,Su2
      DOUBLE COMPLEX    GPS_soft,GPS_softb
      DOUBLE COMPLEX    Sini(2),Sfin(2)
      INTEGER           j,j1,j2,j3,j4,k,Sig,l
      DOUBLE PRECISION  prA1,prB1,prC2,prD2
      DOUBLE PRECISION  BornV_GetCharge,ChaIni,ChaFin
      DOUBLE COMPLEX    gI,gF,sProd
      INTEGER           Y_IR,N_IR
      DOUBLE PRECISION  PP2(4),QQ(4),SvarQ,SvarX2
*----------------------------------------
      Y_IR=1                  ! YES, IR included
      Y_IR=0                  ! No,  IR not included
      N_IR=1-Y_IR
*--------------------
      ChaIni =  BornV_GetCharge( KFi)
      ChaFin =  BornV_GetCharge( KFf)
      gI = DCMPLX(ChaIni *m_e_QED)
      gF = DCMPLX(ChaFin *m_e_QED)
      IF( m_KeyArb  .EQ.  0 ) THEN
         Sini(1)  =  gI *GPS_soft(  1,ph1,pA,pB)
         Sfin(1)  = -gF *GPS_soft(  1,ph2,pC,pD)
      ELSE
         Sini(1)  =  gI *GPS_softb( 1,ph1,pA,mA,pB,mB)
         Sfin(1)  = -gF *GPS_softb( 1,ph2,pC,mC,pD,mD)
      ENDIF
      Sini(2) = -DCONJG(Sini(1))
      Sfin(2) = -DCONJG(Sfin(1))
* Calculate Born spin amplitudes
      CALL GPS_Born(KFi,KFf,PX, pA,mA,    pB,-mB,    pC,MC,    pD,-md,   BornABCD) ! standard
      CALL GPS_Born(KFi,KFf,PX, ph1,mph,  pB,-mB,    pC,mC,    pD,-mD,   Born1BCD) ! A->1
      CALL GPS_Born(KFi,KFf,PX, pA,mA,    ph1,-mph,  pC,mC,    pD,-mD,   BornA1CD) ! B->1
      CALL GPS_Born(KFi,KFf,PX, pA,mA,    pB,-mB,    ph2,mph,  pD,-mD,   BornAB2D) ! C->2
      CALL GPS_Born(KFi,KFf,PX, pA,mA,    pB,-mB,    pC,mC,    ph2,-mph, BornABC2) ! D->2
      CALL GPS_Born(KFi,KFf,PX, ph1,mph,  pB,-mB,    ph2,mph,  pD,-mD,   Born1B2D) ! A->1,C->2
      CALL GPS_Born(KFi,KFf,PX, ph1,mph,  pB,-mB,    pC,mC,    ph2,-mph, Born1BC2) ! A->1,D->2
      CALL GPS_Born(KFi,KFf,PX, pA,mA,    ph1,-mph,  ph2,mph,  pD,-mD,   BornA12D) ! B->1,C->2
      CALL GPS_Born(KFi,KFf,PX, pA,mA,    ph1,-mph,  pC,mC,    ph2,-mph, BornA1C2) ! B->1,D->2
      DO k=1,4
         PP2 (k) = pC(k)+pD(k)+ph2(k)
         QQ(k)   = pC(k)+pD(k)
      ENDDO
      svarX2  =  PP2(4)**2  -PP2(3)**2  -PP2(2)**2  -PP2(1)**2
      svarQ   =   QQ(4)**2   -QQ(3)**2   -QQ(2)**2   -QQ(1)**2
* Fermion propagarotors ini
      prA1=  1d0/(pA(4)*ph1(4)-pA(3)*ph1(3)-pA(2)*ph1(2)-pA(1)*ph1(1))/2d0
      prB1= -1d0/(pB(4)*ph1(4)-pB(3)*ph1(3)-pB(2)*ph1(2)-pB(1)*ph1(1))/2d0
* Fermion propagarotors fin
      prC2=  1d0/(pC(4)*ph2(4)-pC(3)*ph2(3)-pC(2)*ph2(2)-pC(1)*ph2(1))/2d0
      prD2= -1d0/(pD(4)*ph2(4)-pD(3)*ph2(3)-pD(2)*ph2(2)-pD(1)*ph2(1))/2d0
      Sig = 3-2*Hel1
      IF( m_KeyArb  .EQ.  0 ) THEN
         CALL GPS_MatrU( gI, ph1,Sig,  ph1,mph,  pA,mA,    U11a) ! <1|[1]|a>
         CALL GPS_MatrV( gI, ph1,Sig,  pB,mB,   ph1,mph,   Vb11) ! <b|{1}|1>
      ELSE
         CALL GPS_MatrUb(gI, ph1,Sig,  ph1,mph,  pA,mA,    U11a)
         CALL GPS_MatrVb(gI, ph1,Sig,  pB,mB,   ph1,mph,   Vb11)
      ENDIF
      Sig = 3-2*Hel2
      IF( m_KeyArb  .EQ.  0 ) THEN
         CALL GPS_MatrU( gF, ph2,Sig,  pC,mC,   ph2,mph,   Uc22) ! <c|[2]|2>
         CALL GPS_MatrV( gF, ph2,Sig,  ph2,mph,  pD,mD,    V22d) ! <2|{2}|d>
      ELSE
         CALL GPS_MatrUb(gF, ph2,Sig,  pC,mC,   ph2,mph,   Uc22) ! 
         CALL GPS_MatrVb(gF, ph2,Sig, ph2,mph,  pD,mD,     V22d) ! 
      ENDIF
      DO j1=1,2
         DO j2=1,2
            DO j3=1,2
               DO j4=1,2
                  Su1= DCMPLX(0d0,0d0)
                  DO j=1,2
*      /////////////////////////////////////////////////////////////////////////////////////
*      //               c                |2                             d                 //
*      //      u  ------<-------SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS-----<----- v          //
*      //                                        |                                        //
*      //                                        |X             |1                        //
*      //      _       -b                        |     a+m-1    |       a                 //
*      //      v  ------<------------------------O--------------O-------<----- u          //
*      /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +Born1BCD( j,j2,j3,j4) *U11a(j,j1)*prA1 *Sfin(Hel2)*Y_IR !<b|X|1[1]|a><c|(+2)|X|(+2)|d>
*      /////////////////////////////////////////////////////////////////////////////////////
*      //               c                |2                            -d                 //
*      //      u  ------<-------SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS-----<----- v          //
*      //                                        |                                        //
*      //                          |1            |X                                       //
*      //      _       -b          |   -b+m+1    |                      a                 //
*      //      v  ------<----------O-------------O----------------------<----- u          //
*      /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +Vb11(j2,j)*prB1 *BornA1CD(j1,j,j3,j4) *Sfin(Hel2)*Y_IR !<b|[1]1|X|a><c|(+2)|X|(+2)|d>
*      /////////////////////////////////////////////////////////////////////////////////////
*      //                         |2                                                      //
*      //               c         |   c+m+2                            -d                 //
*      //      u  ------<---------O--------------O----------------------<----- v          //
*      //                                        |X                                       //
*      //      _       -b                        |       |1             a                 //
*      //      v  ------<------SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS-----<----- u          //
*      /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +Uc22(j3,j)*prC2 *BornAB2D(j1,j2,j,j4) *Sini(Hel1)*Y_IR !<b|(+1)|X|(+1)|a><c|[2]2|X|d>
*      /////////////////////////////////////////////////////////////////////////////////////
*      //                                                       |2                        //
*      //               c                            -d+m+2     |      -d                 //
*      //      u  ------<------------------------O----------------------<----- v          //
*      //                                        |X                                       //
*      //      _       -b                        |       |1             a                 //
*      //      v  ------<------SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS-----<----- u          //
*      /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +BornABC2(j1,j2,j3,j) *V22d(j,j4)*prD2  *Sini(Hel1)*Y_IR !<b|(+1)|X|(+1)|a><c|X|2[2]d>
                  ENDDO
                  Su2= DCMPLX(0d0,0d0)
                  DO j=1,2
                     DO l=1,2
*      /////////////////////////////////////////////////////////////////////////////////////
*      //                         |2                                                      //
*      //               c         |   c+m+2                             d                 //
*      //      u  ------<---------O--------------O----------------------<----- v          //
*      //                                        |                                        //
*      //                                        |X             |1                        //
*      //      _       -b                        |     a+m-1    |       a                 //
*      //      v  ------<------------------------O--------------O-------<----- u          //
*      /////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Uc22(j3,l )*prC2 *Born1B2D(j,j2,l,j4) *U11a(j,j1)*prA1 !<b|X|1[1]|a><c|[2]2|X|d>
*      /////////////////////////////////////////////////////////////////////////////////////
*      //                                                       |2                        //
*      //               c                             -d+m-2    |       d                 //
*      //      u  ------<------------------------O--------------O-------<----- v          //
*      //                                        |                                        //
*      //                                        |X             |1                        //
*      //      _       -b                        |     a+m-1    |       a                 //
*      //      v  ------<------------------------O--------------O-------<----- u          //
*      /////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Born1BC2(j,j2,j3,l) *U11a(j,j1)*prA1  *V22d(l,j4)*prD2 !<b|X|1[1]|a><c|X|2{2}|d>
*      /////////////////////////////////////////////////////////////////////////////////////
*      //                       |2                                                        //
*      //               c       |    c+m+2                              d                 //
*      //      u  ------<-------O----------------O----------------------<----- v          //
*      //                                        |                                        //
*      //                       |1               |X                                       //
*      //      _       -b       |   -b+m+1       |                      a                 //
*      //      v  ------<-------O----------------O----------------------<----- u          //
*      /////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Vb11(j2,j)*prB1  *Uc22(j3,l )*prC2 *BornA12D(j1,j,l,j4)!<b|{1}1|X|a><c|[2]2|X|d>
*      /////////////////////////////////////////////////////////////////////////////////////
*      //                                                       |2                        //
*      //               c                             -d+m-2    |       d                 //
*      //      u  ------<------------------------O--------------O-------<----- v          //
*      //                                        |                                        //
*      //                       |1               |X                                       //
*      //      _       -b       |   -b+m+1       |                      a                 //
*      //      v  ------<-------O----------------O----------------------<----- u          //
*      /////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Vb11(j2,j)*prB1  *BornA1C2(j1,j,j3,l) *V22d(l,j4)*prD2 !<b|{1}1|X|a><c|X|2[2]|d>
                     ENDDO
                  ENDDO
                  sProd = Sini(Hel1)* Sfin(Hel2)
                  AmpWork(j1,j2,j3,j4) =AmpWork(j1,j2,j3,j4) 
     $                 +CNorm*( Su1+Su2 )
     $                 +CNorm*BornABCD(j1,j2,j3,j4)*sProd                      *Y_IR !
     $                 +CNorm*BornABCD(j1,j2,j3,j4)*sProd *(svarX2/svarQ-1d0)  *N_IR !
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      END                       ! GPS_HifAdd


      SUBROUTINE GPS_HiniAdd(KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4,kPhel,ph,mph,Sactu,SProd,AmpBorn,AmpExpo1) !
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   IR-finite part od 1-photon amplitudes for ISR  (equiv. to GPS_Hini)           //
*//   Photon helicity imported from the CALLing PROGRAM                             //
*//                                                                                 //
*//   AmpExpo1 is working space  (not m_AmpExpo1 and that is the diff with HiniPlus)//
*//   m_AmpBorn   is hidden INPUT (not yet used)                                    //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE "GPS.h"
*
      INTEGER                KFi,KFf
      DOUBLE PRECISION       PX(4),p1(4),p2(4),p3(4),p4(4),ph(4)
      DOUBLE PRECISION       m1,m2,m3,m4,mph
      DOUBLE COMPLEX         Sactu,SProd
      INTEGER                kPhel
      DOUBLE COMPLEX         AmpBorn( 2,2,2,2),AmpExpo1(2,2,2,2)
      DOUBLE COMPLEX         AmpBornV(2,2,2,2),AmpBornU(2,2,2,2)
      DOUBLE COMPLEX         Csum1,Csum2,U(2,2),V(2,2)
      INTEGER                j,j1,j2,j3,j4,k,Sig
      DOUBLE PRECISION       pr1,pr2,Fleps
      DOUBLE PRECISION       BornV_GetCharge,ChaIni
*----------------------------------------
      Fleps =  1d-100
      ChaIni =  BornV_GetCharge( KFi)
* ISR non-infrared two parts: (1) p1 -> photon, contracted with U-matrix
*                             (2) p2 -> photon, contracted with V-matrix
* Calculate Born spin amplitudes
****> CALL GPS_Born(KFi,KFf,PX, p1,Mbeam,  p2,-Mbeam,  p3,Massf,p4,-Massf,m_AmpBorn) !!!!<****
      CALL GPS_Born(KFi,KFf,PX, ph,mph,    p2,-Fleps,  p3,m3,   p4,-m4,   AmpBornU)
      CALL GPS_Born(KFi,KFf,PX, p1,Fleps,  ph,-mph,    p3,m3,   p4,-m4,   AmpBornV)
* Fermion propagarotors
      pr1 = 1d0/(p1(4)*ph(4)-p1(3)*ph(3)-p1(2)*ph(2)-p1(1)*ph(1))/2d0
      pr2 =-1d0/(p2(4)*ph(4)-p2(3)*ph(3)-p2(2)*ph(2)-p2(1)*ph(1))/2d0
      Sig = 3-2*kPhel
      IF( m_KeyArb  .EQ.  0 ) THEN
         CALL GPS_MakeU(ph,Sig,  ph,mph,  p1,m1,    U)
         CALL GPS_MakeV(ph,Sig,  p2,m2,   ph,mph,   V)
      ELSE
         CALL GPS_MakeUb(ph,Sig, ph,mph,  p1,m1,    U)
         CALL GPS_MakeVb(ph,Sig, p2,m2,   ph,mph,   V)
      ENDIF
      DO j1=1,2
         DO j2=1,2
            DO j3=1,2
               DO j4=1,2
                  Csum1=DCMPLX(0d0,0d0)
                  Csum2=DCMPLX(0d0,0d0)
                  DO j=1,2
                     Csum1=Csum1 +DCMPLX(ChaIni *m_e_QED) *U(j,j1)*pr1 *AmpBornU( j,j2,j3,j4) !
                     Csum2=Csum2 +DCMPLX(ChaIni *m_e_QED) *V(j2,j)*pr2 *AmpBornV(j1, j,j3,j4) !
                  ENDDO
                  AmpExpo1(j1,j2,j3,j4) =AmpExpo1(j1,j2,j3,j4) +Sprod/Sactu*(Csum1+Csum2) !
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      END                       !!! GPS_HiniPlus


      SUBROUTINE GPS_HfinAdd(KFi,KFf,PX, p1,m1,p2,m2,p3,m3,p4,m4,kPhel,ph,mph,Sactu,SProd,AmpBorn,AmpExpo1) !
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   IR-finite part od 1-photon amplitudes for FSR  (equiv. to GPS_HfinPlus)       //
*//   Photon helicity is give by the CALLing PROGRAM                                //
*//                                                                                 //
*//   Missing contribution in FSR non-IR part due to  svarX/svarQ                   //
*//   Contribution -svarX/svarQ from HERE cancels exactly with svarX/svarQ in beta0 //
*//                                                                                 //
*//   AmpExpo1  is working space (not m_AmpExpo1 and that is the diff with HiniPlus)// 
*//   m_AmpBorn   is hidden INPUT (not yet used)                                    //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE "GPS.h"
*
      INTEGER                KFi,KFf
      DOUBLE PRECISION       PX(4),p1(4),p2(4),p3(4),p4(4),ph(4)
      DOUBLE PRECISION       m1,m2,m3,m4,mph
*
      DOUBLE COMPLEX         Sactu,SProd
      INTEGER                kPhel
      DOUBLE COMPLEX         AmpBorn( 2,2,2,2),AmpExpo1(2,2,2,2)
      DOUBLE COMPLEX         AmpBornV(2,2,2,2),AmpBornU(2,2,2,2)
      DOUBLE COMPLEX         Csum1,Csum2,U(2,2),V(2,2)
      INTEGER                j,j1,j2,j3,j4,k,Sig
      DOUBLE PRECISION       pr1,pr2,Fleps
      DOUBLE PRECISION       BornV_GetCharge,ChaFin
      DOUBLE PRECISION       svarX1,svarQ,PP1(4),QQ(4),CKine1
*----------------------------------------
      Fleps =  1d-100
      ChaFin =  BornV_GetCharge( KFf)
* FSR non-infrared two parts: (1) p1 -> photon, contracted with U-matrix
*                             (2) p2 -> photon, contracted with V-matrix
****> CALL GPS_Born(KFi,KFf,PX, p1,Mbeam, p2,-Mbeam,  p3,Massf, p4,-Massf,m_AmpBorn) !!!!<****
      CALL GPS_Born(KFi,KFf,PX, p1,Fleps, p2,-Fleps,  ph,mph,   p4,-m4,   AmpBornU)
      CALL GPS_Born(KFi,KFf,PX, p1,Fleps, p2,-Fleps,  p3,m3,    ph,-mph,  AmpBornV)
* Flux correction CKine, due to ONE photon
      DO k=1,4
         QQ(k)  = p3(k)+p4(k)
         PP1(k) = p3(k)+p4(k)+ph(k)
      ENDDO
      svarX1= PP1(4)**2 -PP1(3)**2 -PP1(2)**2 -PP1(1)**2
      svarQ =  QQ(4)**2  -QQ(3)**2  -QQ(2)**2  -QQ(1)**2
      CKine1 = svarX1/svarQ
* Fermion propagarotors
      pr1 = 1d0/(p3(4)*ph(4)-p3(3)*ph(3)-p3(2)*ph(2)-p3(1)*ph(1))/2d0
      pr2 =-1d0/(p4(4)*ph(4)-p4(3)*ph(3)-p4(2)*ph(2)-p4(1)*ph(1))/2d0
      Sig = 3-2*kPhel
      IF( m_KeyArb  .EQ.  0 ) THEN
         CALL GPS_MakeU(ph,Sig,    p3,m3,  ph,mph,   U)
         CALL GPS_MakeV(ph,Sig,    ph,mph, p4,m4,    V)
      ELSE
         CALL GPS_MakeUb(ph,Sig,   p3,m3,  ph,mph,   U)
         CALL GPS_MakeVb(ph,Sig,   ph,mph, p4,m4,    V)
      ENDIF
      DO j1=1,2
         DO j2=1,2
            DO j3=1,2
               DO j4=1,2
                  Csum1=DCMPLX(0d0,0d0)
                  Csum2=DCMPLX(0d0,0d0)
                  DO j=1,2
                     Csum1=Csum1 +DCMPLX(ChaFin *m_e_QED) *U(j3,j)*pr1* AmpBornU(j1,j2, j,j4) !
                     Csum2=Csum2 +DCMPLX(ChaFin *m_e_QED) *V(j,j4)*pr2* AmpBornV(j1,j2,j3, j) !
                  ENDDO
                  AmpExpo1(j1,j2,j3,j4) =AmpExpo1(j1,j2,j3,j4) +Sprod/Sactu*(Csum1+Csum2) !
     $                 +sProd*(1d0-CKine1)*AmpBorn(j1,j2,j3,j4)     ! compenstion for (svarX/svarQ)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      END                       !!! GPS_HfinPlus !!!

      SUBROUTINE GPS_BornAdd(Sprod,AmpBorn,AmpWork) !
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE "GPS.h"
*
      DOUBLE COMPLEX  Sprod
      DOUBLE COMPLEX  AmpBorn(2,2,2,2),AmpWork(2,2,2,2)
      INTEGER    j1,j2,j3,j4
*
      DO j1=1,2
         DO j2=1,2
            DO j3=1,2
               DO j4=1,2
                  AmpWork(j1,j2,j3,j4) =AmpWork(j1,j2,j3,j4) +Sprod*AmpBorn(j1,j2,j3,j4) !
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      END                       !!! GPS_BornAdd !!!
