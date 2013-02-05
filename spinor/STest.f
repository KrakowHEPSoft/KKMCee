*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//                                                                                 //
*//                         Pseudo-CLASS  STest                                     //
*//                                                                                 //
*//       Purpose:  Tests of spinor techniques as programmed in GPS                 //
*//                                                                                 //
*//       Notes:                                                                    //
*//                                                                                 //
*//   GPS class should not be littered with tests any more, all tests from GPS      //
*//   should be mooved here.                                                        //
*//   If some funcionality of GPS lacks then it shoud be copied here.               //
*//                                                                                 //
*//   See README  for notes on the tests.                                           //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////

*///////////////////////////////////////////////////////////////////////////////
*//   make Test_SfacISR
*//   make Test_SfacFSR
*//   make Test_Rmat
*//   make Test_UV
*//   make Test_BornXsec
*//   make Test_IsrSingle
*//   make Test_FsrSingle
*//   make Test_NewSingle
*//   make Test_MultiPhot
*//   make Test_Virtual
*//   make Test_Rules
*//   make Test_DsigOverDtau
*//-----------------------------------------------------------------------------
*//   make Tau-start
*//   make Tau-start-debug
*//
*///////////////////////////////////////////////////////////////////////////////

      SUBROUTINE STest_Initialize
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//    Class initialization                                                         //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
*
      INCLUDE "BXformat.h"
      INCLUDE "STest.h"
*
      INTEGER i,j,k,m,j1,j2
*
      m_out = 16
********************************
*     Define Pauli matrices
      DO k = 0,3
         DO j1 = 1,2
            DO j2 = 1,2
               m_Pauli( k,j1,j2) = DCMPLX(0d0,0d0)
            ENDDO
         ENDDO
      ENDDO
* Sigma0
      m_Pauli( 0,1,1) = DCMPLX( 1d0, 0d0)
      m_Pauli( 0,2,2) = DCMPLX( 1d0, 0d0)
* SigmaX
      m_Pauli( 1,1,2) = DCMPLX( 1d0, 0d0)
      m_Pauli( 1,2,1) = DCMPLX( 1d0, 0d0)
* SigmaY
      m_Pauli( 2,1,2) = DCMPLX( 0d0,-1d0)
      m_Pauli( 2,2,1) = DCMPLX( 0d0, 1d0)
* SigmaZ
      m_Pauli( 3,1,1) = DCMPLX( 1d0, 0d0)
      m_Pauli( 3,2,2) = DCMPLX(-1d0, 0d0)
*
* The other notation for 4-vector index
      DO k = 1,3
         DO j1 = 1,2
            DO j2 = 1,2
               m_Pauli4( k,j1,j2) = m_Pauli( k,j1,j2)
            ENDDO
         ENDDO
      ENDDO
      DO j1 = 1,2
         DO j2 = 1,2
            m_Pauli4( 4,j1,j2) = m_Pauli( 0,j1,j2)
         ENDDO
      ENDDO
********************************
      WRITE(m_out,bxope)
      WRITE(m_out,bxtxt) '  STest   Initializator     '
      WRITE(m_out,bxclo)
      END


      SUBROUTINE STest_TestTemplate(nout,nmax)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   Template of the test routine                                                  //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER nout,nmax
      REAL*8    x,y,z
      INTEGER i,j,k,icont
      SAVE icont
      DATA icont/0/
*-------------
      IF(icont .GE. nmax) RETURN
*     ****************
*     Body of the test
*     ****************
*     Printouts
      WRITE(nout,'(a)') '  '
      WRITE(nout,'(3a)') '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',
     $                              ' STest_TestTemplate ',
     $                   '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      icont=icont+1
      END                       !!!!! STest_TestTemplate


      SUBROUTINE STest_SfacISR(mout,nmax)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   make  Test_SfacISR                                                            //
*//   test of GPS_bfact on ISR photons                                              //
*//   INPUT:                                                                        //
*//   KeyISR=1, KeyFSR=0                                                            //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER mout,nmax
      REAL*8  pf1(4),pf2(4),pp(4)
      INTEGER NphIni
      REAL*8  PhoIni(100,4),yini(100),zini(100),ph(4)

      REAL*8  x,y,z,SfaSud,delp,svar,amel
      INTEGER i,j,k,jph,nout
      INTEGER KFbeam
      REAL*8  BornV_GetMass
      REAL*8 pfd,pk1,pk2,SfaBhl,SfaCur,xk
      REAL*8 jmu(4)

      COMPLEX*16 Sm,sp1,sp2
      COMPLEX*16 GPS_bfact,GPS_soft
      REAL*8  SfaSpi
*     *********************************
      REAL*8  wm0,a,b,del
      wm0(del,a,b)= 1d0 -2d0*del -del*(a/b+b/a)
*     *********************************
      INTEGER icont
      SAVE icont
      DATA icont/0/
*-------------
      nout = mout
ccc      nout = 6
      IF(icont .GE. nmax) RETURN

      KFbeam = 11                      ! KF=11 is electron
      amel   = BornV_GetMass(KFbeam)

      CALL KarLud_GetBeams(pf1,pf2)
      CALL KarLud_GetSudakov(NphIni,yini,zini)
      CALL KarLud_GetPhotons(NphIni,PhoIni)

      DO k=1,4
         pp(k)= pf1(k)+ pf2(k)
      ENDDO
      svar  = pp(4)**2-pp(3)**2-pp(2)**2-pp(1)**2
      delp=  amel**2/svar

      pfd = pf1(4)*pf2(4)-pf1(3)*pf2(3)-pf1(2)*pf2(2)-pf1(1)*pf2(1)

      DO jph=1,NphIni
         DO k=1,4
            ph(k)=PhoIni(jph,k)
         ENDDO
* dot-products k*p using memorized information from MC generation
         pk1  = yini(jph)*svar/2
         pk2  = zini(jph)*svar/2
         SfaBhl = 2*pfd/(2*pk1*2*pk2) -amel**2/(2*pk1)**2 -amel**2/(2*pk2)**2
* clasical current
         DO k=1,4
            jmu(k)= pf1(k)/(2*pk1) - pf2(k)/(2*pk2)
         ENDDO
         SfaCur = -(jmu(4)**2 -jmu(3)**2-jmu(2)**2-jmu(1)**2 )
* as in YFS3 and BHLUMI
         y = yini(jph)
         z = zini(jph)
         SfaSud  =  2d0/(y*z) *wm0(delp,y,z) / (2*svar)
* Spinor version
* Ordinary dot-products k*p, may be replaced by something better
cc         pk1 = pf1(4)*ph(4)-pf1(3)*ph(3)-pf1(2)*ph(2)-pf1(1)*ph(1)
cc         pk2 = pf2(4)*ph(4)-pf2(3)*ph(3)-pf2(2)*ph(2)-pf2(1)*ph(1)
cc         sp1 = GPS_bfact(1,ph,pf1)
cc         sp2 = GPS_bfact(1,ph,pf2)
cc         Sm  = sp1/(2*pk1) - sp2/(2*pk2)
         Sm  = GPS_soft(1,ph,pf1,pf2)
         SfaSpi = Sm*DCONJG(Sm)/2d0

*##################################################################################
         IF( abs(SfaSpi/SfaSud-1).GT. 1d-3  .OR. abs(SfaCur/SfaSud-1).GT. 1d-3) THEN
            Icont=Icont+1
            WRITE(nout,'(a)') 
     $           '------------------------------ISR-------------------------------'
            WRITE(nout,'(a,1f18.12,3e20.7)') 
     $           'y+z, y/(z*delp), z/(y*delp) = ',y+z, y/(z*delp), z/(y*delp)
            WRITE(nout,'(a,6f20.14)') 'wm0(delp,y,z)= ', wm0(delp,y,z)
            WRITE(nout,'(a,6e20.14)') 
     $           'SfaBhl/SfaSud, SfaCur/SfaSud = ', SfaBhl/SfaSud, SfaCur/SfaSud
            WRITE(nout,'(a,6e20.14)') 
     $           'SfaSpi/SfaSud                = ', SfaSpi/SfaSud
         ENDIF
*##################################################################################
      ENDDO
      END



      SUBROUTINE STest_SfacFSR(mout,nmax)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   test of GPS_bfact on FSR photons                                              //
*//   INPUT:                                                                        //
*//   KeyISR=0, KeyFSR=1                                                            //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER mout,nmax
      REAL*8  qf1(4),qf2(4),pp(4)
      INTEGER NphFin
      REAL*8  PhoFin(100,4),yFin(100),zFin(100),ph(4)

      REAL*8  x,y,z,SfaSud,delq,svar,amfin
      INTEGER i,j,k,jph,nout
      INTEGER KFfin
      REAL*8  BornV_GetMass
      REAL*8 qfd,pk1,pk2,SfaBhl,SfaCur,xk
      REAL*8 jmu(4)

      COMPLEX*16 Sm,sp1,sp2
      COMPLEX*16 GPS_bfact,GPS_soft
      REAL*8  SfaSpi
*     *********************************
      REAL*8  wm0,a,b,del
      wm0(del,a,b)= 1d0 -2d0*del -del*(a/b+b/a)
*     *********************************
      INTEGER icont
      SAVE icont
      DATA icont/0/
*-------------
      nout = mout
cc      nout = 6
      IF(icont .GE. nmax) RETURN

* Actual KFcode of final fermion
      CALL KarLud_GetKFfin(KFfin)
      amfin  = BornV_GetMass(KFfin)

      CALL KarFin_GetFermions(qf1,qf2)
      CALL KarFin_GetSudakov(NphFin,yFin,zFin)
      CALL KarFin_GetPhotons(NphFin,PhoFin)

      DO k=1,4
         pp(k)= qf1(k)+ qf2(k)
      ENDDO
      svar  = pp(4)**2-pp(3)**2-pp(2)**2-pp(1)**2
      delq=  amfin**2/svar

      qfd = qf1(4)*qf2(4)-qf1(3)*qf2(3)-qf1(2)*qf2(2)-qf1(1)*qf2(1)

      DO jph=1,NphFin
         DO k=1,4
            ph(k)=PhoFin(jph,k)
         ENDDO
*     dot-products k*p using memorized information from MC generation
         pk1  = yFin(jph)*svar/2
         pk2  = zFin(jph)*svar/2
         SfaBhl = 2*qfd/(2*pk1*2*pk2) -amfin**2/(2*pk1)**2 -amfin**2/(2*pk2)**2
*     clasical current
         DO k=1,4
            jmu(k)= qf1(k)/(2*pk1) - qf2(k)/(2*pk2)
         ENDDO
         SfaCur = -(jmu(4)**2 -jmu(3)**2-jmu(2)**2-jmu(1)**2 )
*     as in YFS3 and BHLUMI
         y = yFin(jph)
         z = zFin(jph)
         SfaSud  =  2d0/(y*z) *wm0(delq,y,z) / (2*svar)
*     Spinor version
*     Ordinary dot-products k*p, may be replaced by something better
cc         pk1 = qf1(4)*ph(4)-qf1(3)*ph(3)-qf1(2)*ph(2)-qf1(1)*ph(1)
cc         pk2 = qf2(4)*ph(4)-qf2(3)*ph(3)-qf2(2)*ph(2)-qf2(1)*ph(1)
cc         sp1 = GPS_bfact(1,ph,qf1)
cc         sp2 = GPS_bfact(1,ph,qf2)
cc         Sm  = sp1/(2*pk1) - sp2/(2*pk2)
         Sm  = GPS_soft(1,ph,qf1,qf2)
         SfaSpi = Sm*DCONJG(Sm) /2d0
         xk = y+z
*##################################################################################
         IF( ABS(SfaSpi/SfaSud-1).GT. 1d-3  .OR. ABS(SfaCur/SfaSud-1).GT. 1d-12) THEN
            Icont=Icont+1
            WRITE(nout,'(a)') 
     $           '------------------------------FSR-------------------------------'
            WRITE(nout,'(a,1f18.12,3e20.7)') 
     $           'y+z, y/(z*delp), z/(y*delq) = ',y+z, y/(z*delq), z/(y*delq)
            WRITE(nout,'(a,6f20.14)') 'wm0(delq,y,z)= ', wm0(delq,y,z)
            WRITE(nout,'(a,6e20.14)') 
     $           'SfaBhl/SfaSud, SfaCur/SfaSud = ', SfaBhl/SfaSud, SfaCur/SfaSud
            WRITE(nout,'(a,6e20.14)') 
     $           'SfaSpi/SfaSud                = ', SfaSpi/SfaSud
         ENDIF
*##################################################################################
      ENDDO
      END




      SUBROUTINE STest_BornRmat(mout,nmax)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   Self-test of Born Rmat tensor                                                 //
*//   Compares Rmat spin correlation tensors from calculeted from spin amplitudes   //
*//   in  GPS_MakeBorn with KORALB analytical formula.                              //
*//   Compensating rotation due to difference in quantization axes is applied.      //
*//   Two Tralor-type routines are combined in order to find out compens. rotation. //
*//                                                                                 //
*//   Validity range:                                                               //
*//      - Z switched off                                                           //
*//      - energy fron ffbar threshold to infinity                                  //
*//   Notes:                                                                        //
*//      In order to have 12-digit agreement please put electron mass to zero       //
*//      in GPS_MakeBorn.                                                           //
*//      Example input for tau:                                                     //
*//      CMSene=4d0, KeyZet=0,  KeyELW=0 (for TH preprint first event )             //
*//      May try CMSene=3.5542d0                                                    //
*//      KeyISR=0, KeyFSR=1  to get more events with no photons                     //
*//                                                                                 //
*//   Resulting printout of the difference of 2 R-matrices should contain zero's!   //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER mout,nmax
      REAL*8 e1(4),e2(4),e3(4)
      REAL*8 f1(4),f2(4),f3(4)
      REAL*8 g1(4),g2(4),g3(4)
      REAL*8 Rot1(4,4),Rot2(4,4)
      REAL*8 PX(4),pf1(4),pf2(4),qf1(4),qf2(4)
      REAL*8 BornX
      INTEGER nout,i,j,k,icont
      INTEGER KFini,KFfin
      SAVE icont
      DATA icont/0/
*-------------
      nout = mout
cccc      nout = 6
      IF(icont .GE. nmax) RETURN
      CALL KarLud_GetBeams(    pf1,pf2)
      CALL KarFin_GetFermions( qf1,qf2)
* check on energy conservation
      IF( ABS((pf1(4)+pf2(4)) - (qf1(4)+qf2(4))) .GT. 1d-10 ) RETURN
      icont=icont+1

      WRITE(nout,'(a)')  '  '
      WRITE(nout,'(3a)') '===============================================',
     $                              ' STest_BornRmTest ',
     $                   '==============================================='
*
      KFini = 11                ! KF=11 is electron
      CALL KarLud_GetKFfin(KFfin)     ! Actual KFcode of final fermion
* Spin amplitudes are defined at this point 
      CALL KarLud_GetPX(PX)
      CALL GPS_MakeBorn(KFini,KFfin,PX,pf1,pf2,qf1,qf2,BornX)
      CALL GPS_BornPrint(nout,0)

      CALL GPS_RmatMake
      CALL GPS_RmatPrint(nout,'GPS frame')
* Initialize GPS tralor
      CALL GPS_TralorPrepare(qf1,1)
      CALL GPS_TralorPrepare(qf2,2)
* We shall transform each versor using Tralor because in general
* we may not know 4x4 Lorenz trasf. matrix corresponding to Tralor
      DO k=1,4
         e1(k) = 0d0
         e2(k) = 0d0
         e3(k) = 0d0
      ENDDO
      e1(1) = 1d0
      e2(2) = 1d0
      e3(3) = 1d0
*----------
* TraJacobWick is to be replaced by any kind of external TRALOR
      CALL GPS_TraJacobWick(-1,qf1,e1,f1) ! to CMS
      CALL GPS_TraJacobWick(-1,qf1,e2,f2) ! to CMS
      CALL GPS_TraJacobWick(-1,qf1,e3,f3) ! to CMS
      CALL GPS_TralorUnDo(   1,f1,g1)     ! to ferm_rest GPS
      CALL GPS_TralorUnDo(   1,f2,g2)     ! to ferm_rest GPS
      CALL GPS_TralorUnDo(   1,f3,g3)     ! to ferm_rest GPS
      CALL KinLib_RotColumn(g1,g2,g3,Rot1)
      CALL KinLib_RotTranspose(Rot1)
***      CALL KinLib_LorPrint(nout,'Rot1    ',Rot1)
*---------
      DO k=1,4
         e1(k) = 0d0
         e2(k) = 0d0
         e3(k) = 0d0
      ENDDO
      e1(1) = 1d0
      e2(2) = 1d0
      e3(3) = 1d0
      CALL GPS_TraJacobWick(-1,qf2,e1,f1) ! to CMS
      CALL GPS_TraJacobWick(-1,qf2,e2,f2) ! to CMS
      CALL GPS_TraJacobWick(-1,qf2,e3,f3) ! to CMS
      CALL GPS_TralorUnDo(   2,f1,g1)     ! to ferm_rest GPS
      CALL GPS_TralorUnDo(   2,f2,g2)     ! to ferm_rest GPS
      CALL GPS_TralorUnDo(   2,f3,g3)     ! to ferm_rest GPS
      CALL KinLib_RotColumn(g1,g2,g3,Rot2)
      CALL KinLib_RotTranspose(Rot2)
***      CALL KinLib_LorPrint(nout,'Rot2    ',Rot2)
*----------
* Transform tensor Rmat from GPS to JacobWick (or any other frame defined by other Tralor)
      CALL GPS_RmatMult(Rot1,Rot2)
* Control printout of spinor result
      CALL GPS_RmatPrint(nout,'JacWick   ')
* Control printout of KORALB style result
* Z-exchange should be off in order to agree!!!
      CALL STest_RmatMock(nout,pf1,pf2,qf1,qf2)
      END


      SUBROUTINE STest_RmatMock(nout,pf1,pf2,qf1,qf2)
*///////////////////////////////////////////////////////////////////////////////
*//                                                                           //
*//   This is MOCK Rmat correlation tensor calculated directly from momenta.  //
*//   using eq. (2.6) in Acta. Phys. Pol. B15 (1984).                         //
*//   It is valid for s-chanel photon exchange only (no Z).                   //
*//   Final state fermion is massive however.                                 //
*//   Quantication axes are as in GPS_TraJacobWick routine, i.e. slightly     //
*//   different from original KORALB convention, z-axix along flight,         //
*//   x-axis in-reaction plane, the same for two fermions                     //
*//   y-axix perpendicular to reaction palne, opposite for two fermions       //
*//                                                                           //
*///////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER      nout
      REAL*8 pf1(4),pf2(4)
      REAL*8 qf1(4),qf2(4)
      INTEGER    k,l
      REAL*8 svar,CosTheta,SinTheta,amfin,Mfin,c2,s2,M2,xn
      REAL*8 Rmat(0:3,0:3),RmatG(0:3,0:3)
*----------------
      svar = (pf1(4) +pf2(4))**2
      CosTheta=        (pf1(3)*qf1(3)+pf1(2)*qf1(2)+pf1(1)*qf1(1))
     $          /SQRT(pf1(3)*pf1(3)+pf1(2)*pf1(2)+pf1(1)*pf1(1))
     $          /SQRT(qf1(3)*qf1(3)+qf1(2)*qf1(2)+qf1(1)*qf1(1))
      amfin=SQRT(qf1(4)**2-qf1(3)**2-qf1(2)**2-qf1(1)**2)
      Mfin = amfin*2d0/sqrt(svar)
      IF(abs(CosTheta) .GT. 1d0) THEN
         WRITE(nout,*) 'STest_RmatMock: CosTheta=',CosTheta
         STOP
      ENDIF
      SinTheta = SQRT(1d0-CosTheta**2)

      DO k=0,3
         DO l=0,3
            Rmat(k,l)=0d0
         ENDDO
      ENDDO
* permutations and reflections with respect to APP are done "by hand"
      c2 = CosTheta**2
      s2 = SinTheta**2
      M2 = Mfin**2
      Rmat(0,0)=    1 +c2  +M2*s2  ! unpolarized
      Rmat(2,2)= -(-1 +c2  +M2*s2) ! +(1-M2)*s2 ! perpendicular to reaction plane
      Rmat(1,1)=    1 -c2  +M2*s2  ! +(1+M2)*s2 ! in-reaction plane
      Rmat(1,3)=  2*Mfin*CosTheta*SinTheta
      Rmat(3,1)= -2*Mfin*CosTheta*SinTheta
      Rmat(3,3)= -(1 +c2  -M2*s2)  ! flight direction (helicity)
      xn = Rmat(0,0)
      DO k=0,3
         DO l=0,3
            Rmat(k,l)=Rmat(k,l)/xn
         ENDDO
      ENDDO
*----------------
      WRITE(nout,*) ' '
      WRITE(nout,'(4a)') '*******************',
     $                   '   R analytical    ',
     $                   ' *******************'
      DO k=0,3
         WRITE(nout,'(4(a,f18.12,a))')  ('[',Rmat(k,l),'] ', l=0,3)
      ENDDO
*----------------
      WRITE(nout,*) ' '
      WRITE(nout,'(4a)') '*******************',
     $                   '   R Difference    ',
     $                   ' *******************'
      CALL GPS_GetRmat(RmatG)
      DO k=0,3
         WRITE(nout,'(4(a,f18.12,a))')  
     $        ('[',RmatG(k,l)-Rmat(k,l),'] ', l=0,3)
      ENDDO
*----------------
      WRITE(nout,'(a,8f18.12)') 'cos(theta),sin(theta)= ',CosTheta,SinTheta
      WRITE(nout,'(a,8f18.12)') 'amfin,sqrt(s)= ',amfin,pf1(4)+pf2(4)
      END




      SUBROUTINE STest_UV(mout,nmax)
*/////////////////////////////////////////////////////////////////////////////////////
*//   make Test_UV                                                                  //
*//                                                                                 //
*//   tests of U and V transition matrices                                          //
*//   INPUT  FSR type                                                               //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER mout,nmax
*
      REAL*8     ph(4),p1(4),p2(4),p3(4),p4(4),pk(4)
      COMPLEX*16 U(2,2),Ub(2,2),UE(2,2)
      COMPLEX*16 V(2,2),Vb(2,2),VE(2,2)
      REAL*8     xk,amel,amfin,Xdist
      REAL*8     BornV_GetMass
      INTEGER    KFi,KFf,sigma
      INTEGER    nphot,nphoy,i,j,k,l,icont,nout
      COMPLEX*16 bf1,bf2
      COMPLEX*16 GPS_bfact, GPS_bfacb
      COMPLEX*16 sm1,sm2,smb1,smb2
      COMPLEX*16 GPS_soft,GPS_softb
      REAL*8     y,z,yy,zz,Mpk
      DATA sigma /1/
*-------------
      SAVE icont
      DATA icont/0/
*-------------
      nout = mout
ccc      nout = 6
      IF(icont .GE. nmax) RETURN
      CALL KK2f_GetNphot(      Nphot)
      CALL KarFin_GetNphot(    nphoy) 
*/////////////////////////////////
*//   select 1-phot pure FSR    //
      IF( nphot .NE. 1) RETURN
      IF( nphoy .NE. 1) RETURN
*/////////////////////////////////
      CALL KarLud_GetBeams(    p1,p2)
      CALL KarFin_GetFermions( p3,p4)
      CALL KarFin_GetPhoton1(   1,ph)
      CALL KarFin_GetSudakov1( 1,yy,zz)
      y  = yy/(1 +yy+zz)
      z  = zz/(1 +yy+zz)
      DO i=1,4
         pk(i) = p4(i)+ph(i)
      ENDDO
      Mpk  = DSQRT(pk(4)**2-pk(3)**2-pk(2)**2-pk(1)**2)

      KFi = 11                ! KF=11 is electron
      CALL KarLud_GetKFfin(KFf)     ! Actual KFcode of final fermion
      amel  =  BornV_GetMass(KFi)
      amfin =  BornV_GetMass(KFf)

      xk= ph(4)/p1(4)
***   IF( (xk .LT. 0.3d0) .OR. (xk .GT. 0.7d0) ) RETURN    !   0.3<vv<0.7
***   IF( (xk  .LT. 0.1d0) )  RETURN      ! minimum v=0.1
***   IF( (xk  .LT. 0.9d0) )  RETURN      ! minimum v=0.9
      IF( (y*z .LT. 0.01d0) ) RETURN      ! minimum pT**2/s=0.01

*//////////////////////////////////////////////////////////////////////////////
*//  Here is test in case of diagonal transition p_1=p_2, m_1=m_2            //
*//////////////////////////////////////////////////////////////////////////////
      sigma= -sigma

      WRITE(nout,'(a)') '  '
      WRITE(nout,'(3a)') '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',
     $                           ' STest_UTest diag ',
     $                   '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      CALL KinLib_VecPrint(nout,'ph=      ',ph)
      WRITE(nout,'(a,i4)') ' Nphot = ', Nphot
****  CALL GPS_Setb2            ! <-- other than default b1, (not-paralel to xi)
      CALL GPS_Setb3            ! <-- beta almost paralel to xi
*///////////////////////
*//     U,Ub,UE       //
      WRITE(nout,'(a,i4)') '|||||| sigma= ',sigma
      CALL GPS_MakeU(  ph,sigma,     p3,amfin,  p3,amfin, U)
      CALL GPS_MakeUb( ph,sigma,     p3,amfin,  p3,amfin, Ub)
      CALL GPS_MakeUE( ph,sigma,     p3,amfin,  p3,amfin, UE)
      CALL GPS_UPrint(nout,'U       ',U)
      CALL GPS_UPrint(nout,'Ub      ',Ub)
      CALL GPS_UPrint(nout,'UE      ',UE)
*///////////////////////
*//      V,Vb,VE      //
      CALL GPS_MakeV(  ph,sigma,     p3,amfin,  p3,amfin, V)
      CALL GPS_MakeVb( ph,sigma,     p3,amfin,  p3,amfin, Vb)
      CALL GPS_MakeVE( ph,sigma,     p3,amfin,  p3,amfin, VE)
      CALL GPS_UPrint(nout,'V       ',V)
      CALL GPS_UPrint(nout,'Vb      ',Vb)
      CALL GPS_UPrint(nout,'VE      ',VE)
*///////////////////////
*//    Sfac, bfac     //
      WRITE(nout,'(a,8g20.13)') ' Sfac, bfac '
* Here you see a miracle, bf1,bf2 are different...
      bf1 = GPS_bfact( 1,ph,p3)
      bf2 = GPS_bfact(-1,ph,p3)
      WRITE(nout,'(a,4g20.13)') 'bfact +1,-1= ',bf1,bf2
      bf1 = GPS_bfacb( 1,ph,p3,amfin)
      bf2 = GPS_bfacb(-1,ph,p3,amfin)
      WRITE(nout,'(a,4g20.13)') 'bfacb +1,-1= ',bf1,bf2
* But soft-factor amplitudes are the same!!!
      sm1 = GPS_soft( 1,ph,p3,p4)
      sm2 = GPS_soft(-1,ph,p3,p4)
      WRITE(nout,'(a,8g20.13)') 'sm1,sm2    = ',sm1,sm2
*/////////////////////
      smb1= GPS_softb( 1,ph,p3,amfin,p4,amfin)
      smb2= GPS_softb(-1,ph,p3,amfin,p4,amfin)
      WRITE(nout,'(a,8g20.13)') 'smb1,smb2  = ',smb1,smb2
      WRITE(nout,'(a,8g20.13)') '|smb1/sm1|, |sm1|= ',CDABS(smb1)**2/CDABS(sm1)**2,CDABS(sm1)**2
*//             End of Diagonal case                                         //
*//////////////////////////////////////////////////////////////////////////////

      WRITE(nout,'(3a)') '==========================================',
     $                           ' STest_UTest gen ',
     $                   '=========================================='
      CALL KinLib_VecPrint(nout,'pk=      ',pk)
*///////////////////////
*//     U,Ub,UE       //
      WRITE(nout,'(a,i4)') '%%%%%% sigma= ',sigma
      CALL GPS_MakeU(  ph,sigma,     p3,amfin,  pk,Mpk, U)
      CALL GPS_MakeUb( ph,sigma,     p3,amfin,  pk,Mpk, Ub)
      CALL GPS_MakeUE( ph,sigma,     p3,amfin,  pk,Mpk, UE)
      CALL GPS_UPrint(nout,'U       ',U)
      CALL GPS_UPrint(nout,'Ub      ',Ub)
      CALL GPS_UPrint(nout,'UE      ',UE)
*
      CALL GPS_MakeV(  ph,sigma,     p3,amfin,  pk,Mpk, V)
      CALL GPS_MakeVb( ph,sigma,     p3,amfin,  pk,Mpk, Vb)
      CALL GPS_MakeVE( ph,sigma,     p3,amfin,  pk,Mpk, VE)
      CALL GPS_UPrint(nout,'V       ',V)
      CALL GPS_UPrint(nout,'Vb      ',Vb)
      CALL GPS_UPrint(nout,'VE      ',VE)

      icont=icont+1
      END                       !!!!! STest_UTest


      SUBROUTINE STest_BornXsec(mout,nmax)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   make Test_BornXsec                                                            //
*//                                                                                 //
*//   Self-test of Born xross section                                               //
*//   Two xsections are compared, one from spinor methods and one calculated        //
*//   in "pedestrial" way using GPS_BornSimple                                      //
*//                                                                                 //
*//   Validity range:                                                               //
*//       KeyISR=1, KeyFSR=1  is OK.                                                //
*//       Z may be ON                                                               //
*//       Energies should be far from the threshold                                 //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER mout,nmax
      REAL*8    p1(4),p2(4),p3(4),p4(4),PX(4)
      REAL*8    BornX,Born0,Born1,Born2,Born3,Born4,Born5,Born6
      REAL*8    Svar,CosTheta
      REAL*8    ct11,ct12,ct21,ct22
      REAL*8    Born11,Born12,Born21,Born22
      REAL*8    tt,uu,t1,u1
      REAL*8    SvarX,SvarQ,vv,vQ
      REAL*8    BornV_Simple
      REAL*8    BornV_GetMass, amfin, betaf
      INTEGER   i,j,k,icont,nout
      INTEGER   KFi,KFf
      SAVE icont
      DATA icont/0/
*-------------
      nout = mout
cccc      nout = 6
      IF(icont .GE. nmax) RETURN
      CALL KarLud_GetBeams(    p1,p2)
      CALL KarFin_GetFermions( p3,p4)
      CALL KarLud_GetPX(PX)
      svarX = PX(4)**2 - PX(3)**2 - PX(2)**2 - PX(1)**2
      Svar  = (p1(4)+p2(4))**2-(p1(3)+p2(3))**2-(p1(2)+p2(2))**2-(p1(1)+p2(1))**2
      SvarQ = (p3(4)+p4(4))**2-(p3(3)+p4(3))**2-(p3(2)+p4(2))**2-(p3(1)+p4(1))**2
      tt    = 2d0*( p1(4)*p3(4)-p1(3)*p3(3)-p1(2)*p3(2)-p1(1)*p3(1) )
      uu    = 2d0*( p1(4)*p4(4)-p1(3)*p4(3)-p1(2)*p4(2)-p1(1)*p4(1) )
      t1    = 2d0*( p2(4)*p4(4)-p2(3)*p4(3)-p2(2)*p4(2)-p2(1)*p4(1) )
      u1    = 2d0*( p2(4)*p3(4)-p2(3)*p3(3)-p2(2)*p3(2)-p2(1)*p3(1) )
      vv= 1d0-svarX/svar
      vQ= 1d0-svarQ/svar
* check on energy conservation
c[[      IF( vQ .GT. 1d-6 ) RETURN
*-------------
      KFi = 11                ! KF=11 is electron
      CALL KarLud_GetKFfin(KFf)     ! Actual KFcode of final fermion
      amfin   = BornV_GetMass(KFf)
      betaf   = SQRT(1-4*amfin**2/svarQ)
      CALL KinLib_ThetaR(PX,p1,p2,p3,p4,ct11,ct12,ct21,ct22)
      CALL GPS_BornSimple(KFi,KFf,  SvarX, ct11, Born11)
      CALL GPS_BornSimple(KFi,KFf,  SvarX, ct12, Born12)
      CALL GPS_BornSimple(KFi,KFf,  SvarX, ct21, Born21)
      CALL GPS_BornSimple(KFi,KFf,  SvarX, ct22, Born22)
      Born1 = (Born11+Born12+Born21+Born22)/4d0
      Born2 = (BornV_Simple( KFi,KFf,svarX, ct11  )
     $        +BornV_Simple( KFi,KFf,svarX, ct12  )
     $        +BornV_Simple( KFi,KFf,svarX, ct21  )
     $        +BornV_Simple( KFi,KFf,svarX, ct22  ))/4d0
*     /////////////////////////////////////////////////
      Born0 = BornV_Simple( KFi,KFf,svarX, 0d0  )
*     /////////////////////////////////////////////////
      Born3 = (tt**2 + uu**2 + t1**2 +u1**2)/svar/svarQ *(svarX/svar) !!<- Usefull ONLY for Z off
      Born4 = 2*(tt*t1 + uu*u1)/(svar*svarQ) *Born0                   !!<- Usefull ONLY for Z off
      Born5 = (tt**2 + uu**2 + t1**2 +u1**2)/svar**2
c[[[      Born6 = 2*(tt*t1 + uu*u1)/(svar**2)
      Born6 = 2*(tt*t1 + uu*u1)/svar/svarQ
*     /////////////////////////////////////////////////
      CALL GPS_MakeBorn(  KFi,KFf,PX,p1,p2,p3,p4,BornX)
*     /////////////////////////////////////////////////
****  BornX = BornX*(svar/svarQ) *(svarX/svar)**2
      BornX = BornX  *svarX**2 /(svar *svarQ)
*******************************************************************************************
*******************************************************************************************
* For high v values, hunt for special cases
*      IF( Born4 .LT. 20d0 ) RETURN
      IF( SQRT(svarQ) . GT. 10d0 ) RETURN
*******************************************************************************************
*******************************************************************************************
      WRITE(nout,'(a)') '  '
      WRITE(nout,'(3a)') '/////////////////////////////////////////',
     $                   ' STest_BornXsTest ',
     $                   '/////////////////////////////////////////'
***   CALL KarLud_Print1(6)
***   CALL KK2f_Print1(6)
      WRITE(nout,'(a,6f25.14)') 'vv,vQ, betaf=  ',           vv,vQ,betaf
      WRITE(nout,'(a,6f25.14)') 'Born4 =        ',           Born4
      WRITE(nout,'(a,6f25.14)') 'Born5, Born6 = ',           Born5,Born6
      WRITE(nout,'(a,6f25.14)') 'Born0, BornX = ',           Born0, BornX
      WRITE(nout,'(a,6f25.14)') 'BornX/Born0, BornX/Born4',  BornX/Born0, BornX/Born4
      WRITE(nout,'(a,6f25.14)') 'Born1/Born2, Born3/Born0',  Born1/Born2, Born3/Born0
      WRITE(nout,'(a,6f25.14)') 
     $                ' GPS_MakeBorn, BornSimple, ratio= ',  BornX,Born1,BornX/Born1
*
      icont=icont+1
      END                       !STest_BornXsTest


      SUBROUTINE STest_IsrSingle(mout,nmax)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   make Test_IsrSingle                                                           //
*//                                                                                 //
*//   test of ISR single bremss matrix element                                      //
*//   compared are also ERW subprograms                                             //
*//                                                                                 //
*//   Validity range:                                                               //
*//       Z may be ON                                                               //
*//       Energies should be far from the threshold                                 //
*//       KeyISR=1 !!!!                                                             //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER  mout,nmax
      REAL*8   p1(4),p2(4),p3(4),p4(4),ph(4)
      REAL*8   y,z,xk
      COMPLEX*16 Amp1isr(2,2,2,2,2)
      INTEGER  i,j,k,icont,idv,nout
      INTEGER  KFi,KFf,nphot
      REAL*8   svar,svar1,amel,amfin,delp,ratio,vv
      REAL*8   KLsect,Xsect3,Dig1,Xdist,Xdist0
      REAL*8   BornV_GetMass
      REAL*8   wm0,wmd,wm1,wt
*-------------
      SAVE icont
      DATA icont/0/
*-------------
      nout = mout
ccc      nout = 6
      CALL KarLud_GetBeams(    p1,p2)
      CALL KarFin_GetFermions( p3,p4)
      CALL KarLud_GetNphot(nphot)     
* check photon multiplicity
*---------------------------------------
      IF( nphot .NE. 1 ) RETURN
*---------------------------------------
      KFi = 11                ! KF=11 is electron
      CALL KarLud_GetKFfin(KFf)     ! Actual KFcode of final fermion
      CALL KarLud_GetSudakov1( 1,y,z)
      CALL KarLud_GetPhoton1(  1,ph)
      amel   = BornV_GetMass(KFi)
      svar  = (p1(4)+p2(4))**2
      svar1 = (p3(4)+p4(4))**2-(p3(3)+p4(3))**2-(p3(2)+p4(2))**2-(p3(1)+p4(1))**2
      vv    = 1-svar1/svar
      delp=  amel**2/svar
      wm0= 1d0 -2d0*delp -delp*(y/z+z/y)
      wmd= 1d0 + delp*( y/z +z/y ) *(y**2+z**2) / ((1-y)**2+(1-z)**2)
      wm1= 1d0 - delp*( y/z+z/y ) *(1d0-y)*(1d0-z) *2d0/((1d0-y)**2+(1d0-z)**2)
********************************
c      CALL GPS_Setb2
c      CALL GPS_Setb3
c      CALL  GPS_SetKeyArb(1)

      CALL STest_Make1isr(KFi,KFf,p1,p2,p3,p4,ph,Amp1isr,Xdist,Xdist0)
**      CALL GPS_Amp1Print(nout,'Amp1isr',Amp1isr)
**
      wt = Xdist/Xdist0*wm0
      idv  = 50000
      CALL GLK_Fil1(idv+19, wt, 1d0)
**
      CALL QED3_GetDig1(Dig1)
      ratio = Xdist/Dig1*(svar/2d0)
      CALL ERW_initialize
      CALL ERW_SKLINI(KLsect)
      CALL ERW_SingIni(Xsect3)
********************************
*----------------------------------------------------
      IF(icont .GE. nmax) RETURN
*----------------------------------------------------
*-------------------------------------------------------
cc      IF( vv .LT. 0.9d0 ) RETURN
      IF( (vv .LT. 0.2d0) .OR. (vv .GT. 0.9d0) ) RETURN
cc      IF( ABS(ratio-1d0)  .LT. 1d-5)  RETURN
cc      IF( ABS(ratio-1d0)  .LT. 1d-4)  RETURN
*-------------------------------------------------------
      WRITE(nout,'(a)') '  '
      WRITE(nout,'(3a)') '///////////////////////////////////////////////',
     $                                  ' STest_isrSingle ',
     $                   '///////////////////////////////////////////////'
      WRITE(nout,'(a,3f20.14)') 'wm0,wmd,wm1 = ',wm0,wmd,wm1
c      CALL KinLib_VecPrint(nout,'p1=      ',p1)
c      CALL KinLib_VecPrint(nout,'p2=      ',p2)
c      CALL KinLib_VecPrint(nout,'p3=      ',p3)
c      CALL KinLib_VecPrint(nout,'p4=      ',p4)
      CALL KinLib_VecPrint(nout,'ph=      ',ph)
      WRITE(nout,'(a,1f18.12,3e20.7)') 
     $     'y+z, y/(z*delp), z/(y*delp) = ',y+z, y/(z*delp), z/(y*delp)
      WRITE(nout,'(a,5g20.14)') 'Make1isr:  Xdist,Xdist0,Xdist/Xdist0= ', Xdist,Xdist0,Xdist/Xdist0
* svar factor due to normalization of s-factors
      WRITE(nout,'(a,5g20.14)') ' ratio = Xdist/Dig1*(svar/2d0)            = ', ratio
      WRITE(nout,'(a,5g20.14)') ' KLsect = ', KLsect, 2*KLsect/Xdist
      WRITE(nout,'(a,5g20.14)') ' Xsect3 = ', Xsect3,   Xsect3/Xdist
*------------------
      icont=icont+1
      END                       !STest_isrTest


      SUBROUTINE STest_FsrSingle(mout,nmax)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   make Test_FsrSingle                                                           //
*//                                                                                 //
*//   test of FSR single bremss matrix element                                      //
*//   compared are also ERW subprograms                                             //
*//                                                                                 //
*//   Validity range:                                                               //
*//       Z may be ON --- still some problem for keyZet=1 with ERW subprograms.     //
*//       Energies should be far from the threshold                                 //
*//       KeyFSR=1!!!                                                               //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER  mout,nmax
      REAL*8   p1(4),p2(4),p3(4),p4(4),ph(4)
      REAL*8   yy,zz,y,z,xk,vv
      COMPLEX*16 Amp1fsr(2,2,2,2,2)
      INTEGER  i,j,k,icont,idv,nout
      INTEGER  KFi,KFf,nphot
      REAL*8   svar,svar1,amel,amfin,delq,ratio
      REAL*8   KLsect,Xsect3,Dig1,Xdist,Xdist0
      REAL*8   BornV_GetMass
      REAL*8   wm0,wmd,wm1
      REAL*8   wt,WtMain,WtCrud
*-------------
      SAVE icont
      DATA icont/0/
*----------------------------------------------------
      IF(icont .GE. nmax) RETURN
*----------------------------------------------------
      nout = mout
ccc      nout = 6
      CALL KarLud_GetBeams(    p1,p2)
      CALL KarFin_GetFermions( p3,p4)
      CALL KarFin_GetNphot(nphot)     
* check photon multiplicity
*---------------------------------------
      IF( nphot .NE. 1 ) RETURN
*---------------------------------------
      KFi = 11                ! KF=11 is electron
      CALL KarLud_GetKFfin(KFf)     ! Actual KFcode of final fermion
      CALL KarFin_GetSudakov1( 1,yy,zz)
      y  = yy/(1 +yy+zz)
      z  = zz/(1 +yy+zz)
      CALL KarFin_GetPhoton1(  1,ph)
      amfin   = BornV_GetMass(KFf)
      svar  = (p1(4)+p2(4))**2
      svar1 = (p3(4)+p4(4))**2-(p3(3)+p4(3))**2-(p3(2)+p4(2))**2-(p3(1)+p4(1))**2
      vv = 1-svar1/svar
      delq=  amfin**2/svar1
      wm0= 1d0 -2d0*delq -delq*(y/z+z/y)
      wmd= 1d0 + delq*( y/z +z/y ) *(y**2+z**2) / ((1-y)**2+(1-z)**2)
      wm1= 1d0 - delq*( y/z+z/y ) *(1d0-y)*(1d0-z) *2d0/((1d0-y)**2+(1d0-z)**2)
*------------------------------------------------------------------
cc      IF( (vv .LT. 0.2d0) .OR. (vv .GT. 0.9d0) ) RETURN
      IF( vv .LT. 0.8d0 ) RETURN
cc      IF( vv .LT. 0.1d0 ) RETURN
*------------------------------------------------------------------
********************************
c      CALL  GPS_SetKeyArb(1)
c      CALL GPS_Setb2
c      CALL GPS_Setb3

      CALL STest_Make1fsr(KFi,KFf,p1,p2,p3,p4,ph,Amp1fsr,Xdist,Xdist0)
***   CALL GPS_Amp1Print(nout,'Amp1fsr',Amp1fsr)

      wt = Xdist/Xdist0
**
**      CALL KK2f_GetWt(WtMain,WtCrud) !! events shoud be weighted !!
      idv  = 50000
      CALL GLK_Fil1(idv+19, wt*wm0, 1d0)

      CALL QED3_GetDig1(Dig1)
      ratio = Xdist/Dig1*(svar/2d0)
      CALL ERW_initialize       !!!
      CALL ERW_SklFin(KLsect)
      CALL ERW_SingFin(Xsect3)  !!!
********************************
cc      IF( ABS(ratio-1d0)  .LT. 1d-4)  RETURN
*------------------------------------------------------------------
      WRITE(nout,'(a)') '  '
      WRITE(nout,'(3a)') '///////////////////////////////////////////////',
     $                                  ' STest_fsrSingle ',
     $                   '///////////////////////////////////////////////'
      WRITE(nout,'(a,3f20.14)') 'wm0,wmd,wm1 = ',wm0,wmd,wm1
c      CALL KinLib_VecPrint(nout,'p1=      ',p1)
c      CALL KinLib_VecPrint(nout,'p2=      ',p2)
c      CALL KinLib_VecPrint(nout,'p3=      ',p3)
c      CALL KinLib_VecPrint(nout,'p4=      ',p4)
      CALL KinLib_VecPrint(nout,'ph=      ',ph)
      WRITE(nout,'(a,1f18.12,3e20.7)') 
     $     'y+z, y/(z*delq), z/(y*delq) = ',y+z, y/(z*delq), z/(y*delq)
      WRITE(nout,'(a,5g20.14)') ' Dig1=                ', Dig1
      WRITE(nout,'(a,5g20.14)') ' Make1fsr: Xdist,Xdist0,X/X0/2= ', Xdist,Xdist0,Xdist/Xdist0/2
* svar factor due to normalization of s-factors
      WRITE(nout,'(a,5g20.14)') ' ratio = Xdist/Dig1*(svar/2d0)            = ', ratio
      WRITE(nout,'(a,5g20.14)') ' Xsect3, Xdist/Xsect3    = ',  Xsect3, Xdist/Xsect3
      WRITE(nout,'(a,5g20.14)') '2*KLsect, Xdist/(2*KLsect)=',2*KLsect, Xdist/(2*KLsect)
*------------------
      icont=icont+1
      END                       !STest_fsrSingle



      SUBROUTINE STest_MultiPhot(mout,nmax)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   make Test_MultiPhot                                                           //
*//                                                                                 //
*//   Test of ISR+FSR multiphoton exponentiated matrix element                      //
*//                                                                                 //
*//   Weight distribution is examined, note that mass weight has to be included     //
*//   in order to judge properly the tail of the weight!                            //
*//                                                                                 //
*//                                                                                 //
*//   Validity range:                                                               //
*//       KeyISR=1 and KeyFSR=1                                                     //
*//                                                                                 //
*//                                                                                 //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER    mout,nmax
*
      REAL*8     x,y,z
      INTEGER    i,j,k,icont
      REAL*8     m1,m2,m3,m4,mph,Fleps,betaf
      INTEGER    KFi,KFf
      REAL*8     amfin,amel
      REAL*8     BornV_GetMass
      REAL*8     p1(4),p2(4),p3(4),p4(4),ph(4)
      REAL*8     PhoAll(100,4)
      REAL*8     svar,svar1,delp,delq,vv
      REAL*8     Xdist,Xdist0
      INTEGER    nphot,nphox,nphoy,idv,nout
      COMPLEX*16 Amp1pho( 2,2,2,2,2)
      REAL*8     wt,WtMain,WtCrud,WtMass
      REAL*8     WtSet(1000), WtExp0, WtExp1, WtExp2
      REAL*8     WtExp0max, WtExp1max, WtExp2max
*--------------------------------------------------------------
      SAVE icont
      DATA icont/0/
      DATA WtExp0max, WtExp1max /0d0,0d0/
*--------------------------------------------------------------
      nout = mout
cc      nout = 6
ccccccccccccccc      CALL GPS_Initialize
      KFi = 11                ! KF=11 is electron
      CALL KarLud_GetKFfin(KFf)     ! Actual KFcode of final fermion
      amel  =  BornV_GetMass(KFi)
      amfin =  BornV_GetMass(KFf)
      Fleps =  1d-100
      m1  = amel
      m2  = amel
      m3  = amfin
      m4  = amfin
      mph = Fleps
*
      CALL KarLud_GetBeams(    p1,p2)
      CALL KarFin_GetFermions( p3,p4)

      svar  = (p1(4)+p2(4))**2
      svar1 = (p3(4)+p4(4))**2-(p3(3)+p4(3))**2-(p3(2)+p4(2))**2-(p3(1)+p4(1))**2
      vv = 1-svar1/svar
      delp=  amel**2/svar
      delq=  amfin**2/svar1

      CALL KK2f_GetPhotAll(Nphot,PhoAll) ! ordered in energy

      CALL KarLud_GetNphot(nphox)
      CALL KarFin_GetNphot(nphoy)

      WtMass = 1d0
      DO j=1,nphox
         CALL KarLud_GetSudakov1( j,y,z)
         WtMass= WtMass*(1d0 -2d0*delp -delp*(y/z+z/y))
      ENDDO
      DO j=1,nphoy
         CALL KarFin_GetSudakov1( j,y,z)
         WtMass= WtMass*(1d0 -2d0*delq -delq*(y/z+z/y))
      ENDDO

* Exponentiation weight
      CALL KK2f_GetWtAll(WtMain,WtCrud,WtSet)
***   WtExp0 = WtSet(251)*WtCrud  ! Interf OFF
***   WtExp1 = WtSet(252)*WtCrud  ! Interf OFF
      WtExp0 = WtSet(201)*WtCrud
      WtExp1 = WtSet(202)*WtCrud
      WtExp2 = WtSet(203)*WtCrud
******************************************************************************
c      CALL KK2f_GetWt(WtMain,WtCrud) !! events shoud be weighted !!
c      DO k=1,4
c         ph(k) = PhoAll(1,k)
c      ENDDO
c      CALL GPS_Make1pho(KFi,KFf,p1,p2,p3,p4,ph,Amp1pho,Xdist,Xdist0)
c      WtExp1 = Xdist/Xdist0*WtMass
c      idv  = 50000
c      CALL GLK_Fil1(idv+19, WtExp1*WtMass , 1d0)
c      CALL GLK_Fil1(idv+20, WtExp0, 1d0)
******************************************************************************
*************************@@@@@@@@@@@******************************************
***      IF( Nphot .LT. 2 ) RETURN
***      IF( vv .GT. 0.999d0 ) RETURN          !<-- simultaneous cut of ISR and FSR
ccc      IF( WtExp0 .LT. WtExp0max ) RETURN
ccc      WtExp0max = WtExp0
ccc      IF( WtExp1 .LT. WtExp1max ) RETURN
ccc      WtExp1max = WtExp1
      IF( WtExp2 .LT. WtExp2max ) RETURN
      WtExp2max = WtExp2
******************************************************************************
*=============================================================================
      IF(icont .GE. nmax) RETURN
*=============================================================================
******************************************************************************
      WRITE(nout,'(a)') '  '
      WRITE(nout,'(3a)') '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',
     $                              ' STest_MultiPhot ',
     $                   '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

      WRITE(nout,'(a,i10,5g20.9)') '--> ISR+FSR nphot, vvQ, massQ  = ',  nphot,vv, SQRT(svar1)
      WRITE(   *,'(a,i10,5g20.9)') '--> ISR+FSR nphot, vvQ, massQ  = ',  nphot,vv, SQRT(svar1)

      betaf = SQRT(1-4*amfin**2/svar1)
      WRITE(   *,'(a,5g20.9)') '--> betaf  = ',  betaf

      CALL KK2f_Print1(nout)
c[[      CALL KinLib_VecPrint(nout,'p1=      ',p1)
c[[      CALL KinLib_VecPrint(nout,'p2=      ',p2)
      CALL KinLib_VecPrint(nout,'p3=      ',p3)
      CALL KinLib_VecPrint(nout,'p4=      ',p4)
      DO j=1,nphot
         DO k=1,4
            ph(k) = PhoAll(j,k)
         ENDDO
         CALL KinLib_VecPrint(nout,'ph=      ',Ph)
      ENDDO
*-------
      WtMass = 1d0
      DO j=1,nphox
         CALL KarLud_GetSudakov1( j,y,z)
         WtMass= WtMass*(1d0 -2d0*delp -delp*(y/z+z/y))
         WRITE(nout,'(a,1f18.12,3e20.7)') 
     $     'ISR: y+z, y/(z*delp), z/(y*delp) = ',y+z, y/(z*delp), z/(y*delp),WtMass
      ENDDO
      DO j=1,nphoy
         CALL KarFin_GetSudakov1( j,y,z)
         WtMass= WtMass*(1d0 -2d0*delq -delq*(y/z+z/y))
         WRITE(nout,'(a,1f18.12,3e20.7)') 
     $     'FSR: y+z, y/(z*delq), z/(y*delq) = ',y+z, y/(z*delq), z/(y*delq),WtMass
      ENDDO

      WRITE(nout,'(a,5g20.14)') ' WtExp0,1,2 WtMass=    ',WtExp0,WtExp1,WtExp2, WtMass
      WRITE(   *,'(a,5g20.14)') ' WtExp0,1,2 WtMass=    ',WtExp0,WtExp1,WtExp2, WtMass

      icont=icont+1
      END                       !!!!! STest_MultiPhot


      SUBROUTINE STest_Virtual(mout,nmax)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//                                                                                 //
*//   make Test_Virtual                                                             //
*//                                                                                 //
*//   Formfactors and Virtual corrections                                           //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE "STest.h"

      INTEGER mout,nmax
      REAL*8     x,y,z
      INTEGER    i,j,k,icont
      REAL*8     alfpi,Emin,AlfInv
      REAL*8     p1(4),p2(4),p3(4),p4(4),m1,m2,m3,m4,mBeam,Massf,MasPhot
      INTEGER    KFi,KFf,Nphot,nout
      REAL*8     Yint
      REAL*8     BVR_SForFac,BVR_TForFac
      REAL*8     BornV_GetCharge
      REAL*8     BornV_GetMass
      REAL*8     BVR_TBvirtExact,BVR_TBvirt
      REAL*8     BVR_SBvirt, BVR_SBvirt2
      REAL*8     BVR_TBvirt2
      REAL*8     SBvirt,TBvirt,TBvirt2
      REAL*8     p1p2,p3p4, p1p3,p2p4, p1p4,p2p3, ss,s1,tt,t1,uu,u1
      REAL*8     result1,result2
      REAL*8     BVR_Btilda,BVR_Dilog
      COMPLEX*16 BVR_CBoxGG,BVR_CBoxGZ,BVR_IntIR
      COMPLEX*16 CBoxGG,CBoxGZ
      REAL*8     MassZ,GammZ
      REAL*8     ReBox,ImBox,ReIR,DelBox
      COMPLEX*16 s31,s24,s14,s32
      COMPLEX*16 GPS_iProd2
      COMPLEX*16 BVR_CnuA
      REAL*8     Xborn,Xboxy
      REAL*8     Yfsr
      REAL*8 TBvirtDif1, TBvirtDif2
      REAL*8 TBrealDif1, TBrealDif2
*---------------------------
      REAL*8     Ymax
      DATA       Ymax/0d0/
*--------------------------
      SAVE icont
      DATA icont/0/
*--------------------------------------------------
      IF(icont .GE. nmax) RETURN
*--------------------------------------------------
      nout =  mout
ccc      nout = 6
      CALL KK2f_GetNphot(Nphot)
      CALL GPS_GetDebg(1,XBorn)
      CALL GPS_GetDebg(2,XBoxy)

      CALL KK2f_GetKFini( KFi)     ! Normaly for beam KFi=11 is electron
      CALL KarLud_GetKFfin(       KFf)     ! Actual KFcode of the final fermion
      Mbeam  =  BornV_GetMass(   KFi)
      Massf  =  BornV_GetMass(   KFf)
      m1  = Mbeam
      m2  = Mbeam
      m3  = Massf
      m4  = Massf
      CALL KK2f_GetBeams(    p1,p2)
      CALL KK2f_GetFermions( p3,p4)
      CALL KK2f_GetMasPhot(MasPhot)
      CALL KK2f_GetEmin( Emin)
      CALL BornV_GetAlfInv(AlfInv)
      alfpi =1/AlfInv/m_pi

      p1p3  =  p1(4)*p3(4) -p1(3)*p3(3) -p1(2)*p3(2) -p1(1)*p3(1)
      p2p4  =  p2(4)*p4(4) -p2(3)*p4(3) -p2(2)*p4(2) -p2(1)*p4(1)
      p1p4  =  p1(4)*p4(4) -p1(3)*p4(3) -p1(2)*p4(2) -p1(1)*p4(1)
      p2p3  =  p2(4)*p3(4) -p2(3)*p3(3) -p2(2)*p3(2) -p2(1)*p3(1)
      p1p2  =  p1(4)*p2(4) -p1(3)*p2(3) -p1(2)*p2(2) -p1(1)*p2(1)
      p3p4  =  p3(4)*p4(4) -p3(3)*p4(3) -p3(2)*p4(2) -p3(1)*p4(1)
      ss = m1**2 +m2**2 +2d0*p1p2
      s1 = m3**2 +m4**2 +2d0*p3p4
      tt = m1**2 +m3**2 -2d0*p1p3
      t1 = m2**2 +m4**2 -2d0*p2p4
      uu = m1**2 +m4**2 -2d0*p1p4
      u1 = m2**2 +m3**2 -2d0*p2p3
      s31 = GPS_iProd2(  1, p3, m3,   -1, p1, m1) ! t
      s24 = GPS_iProd2(  1, p2, m2,   -1, p4, m4) ! t1
      s14 = GPS_iProd2( -1, p1,-m1,   -1, p4, m4) ! u
      s32 = GPS_iProd2(  1, p3, m3,    1, p2,-m2) ! u1

      TBvirtDif1 = alfpi*( DLOG(tt/uu)*DLOG(MasPhot**2/DSQRT(tt*uu)) +0.5d0*DLOG(tt/uu))
      TBvirtDif2 = BVR_TBvirt(alfpi,p1p3,m1,m3,MasPhot)
     $            -BVR_TBvirt(alfpi,p1p4,m1,m4,MasPhot)

      TBrealDif1 = alfpi*( DLOG(tt/uu)*DLOG(4*Emin**2/MasPhot**2)
     $                     +0.5d0*DLOG(-tt/ss)**2   -0.5d0*DLOG(-uu/ss)**2
     $                      -BVR_Dilog(-tt/ss)       +BVR_Dilog(-uu/ss)
     $     )
      TBrealDif2 = 
     $     +BVR_Btilda(alfpi, p1p3, p1(4),p3(4), m1, m3,  Emin, MasPhot) !! Exact
     $     -BVR_Btilda(alfpi, p1p4, p1(4),p4(4), m1, m4,  Emin, MasPhot) !! Exact

      CALL BornV_GetMZ(   MassZ)
      CALL BornV_GetGammZ(GammZ)
      CBoxGG = BVR_CBoxGG(MasPhot,ss,tt,uu)
     $        -BVR_IntIR( MasPhot,ss,tt,uu)
**********
****  MassZ = MasPhot
****  GammZ = MasPhot*1d-8
**********
      CBoxGZ = BVR_CBoxGZ(MasPhot,MassZ,GammZ,ss,tt,uu)
     $        -BVR_IntIR( MasPhot,            ss,tt,uu)
      CALL BVR_RBoxGG(MasPhot,ss,tt,uu,ReBox,ImBox,ReIR,DelBox)

      TBvirt    = BVR_TBvirt(      alfpi, p1p3, m1, m3, MasPhot)
****  TBvirt2   = BVR_TBvirt2( alfpi, p1p3, m1, m3, MasPhot)
      TBvirt2   = BVR_TBvirtExact( alfpi, p1p3, m1, m3, MasPhot)

      SBvirt    = BVR_SBvirt( alfpi,p3p4,m3,m4,MasPhot)

      Yfsr = BVR_SForFac( alfpi, p3,Massf, p4,Massf, Emin, MasPhot)

*****************************************************************************
*****************************************************************************
*****************************************************************************
cccc      IF( ABS(XBoxy/Xborn) .LT. 0.10d0)    RETURN
cccc      IF( Nphot .NE. 1 ) RETURN
cccc The test with Yfsr is not finished, pi**2/beta problem needs to be solved
      IF( Yfsr .LT. Ymax )  RETURN
      Ymax=Yfsr
*****************************************************************************
*****************************************************************************
*****************************************************************************
      WRITE(nout,'(a)') '  '
      WRITE(nout,'(3a)') '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',
     $                              ' STest_VIRTUAL ',
     $                   '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      WRITE(nout,'(a,5g20.12)') 'Nphot =', Nphot

      CALL KK2f_Print1(nout)
      WRITE(nout,'(a,5g20.10)') '#######>>> Boxy/Born = ',Xboxy/XBorn

      WRITE(nout,'(a,5g20.12)') 'MasPhot =', MasPhot
*----------------------------------------------------------------------------
c[[      Yint= BVR_TForFac(-alfpi, p1,Mbeam, p3,Massf, Emin, MasPhot)
*----------------------------------------------------------------------------
      WRITE(nout,'(a,9g19.12)') 'ss ...u1  = ', ss,s1,tt,t1,uu,u1
      WRITE(nout,'(a,9g19.12)') ' |s31|**2/tt ... = ', 
     $     -CDABS(s31)**2/tt, -CDABS(s24)**2/t1, -CDABS(s14)**2/uu, -CDABS(s32)**2/u1
      WRITE(nout,'(a,9g19.12)') 'SQRT(ss*s1),SQRT(tt*tt1),SQRT(uu*u1)= ', 
     $     SQRT(ss*s1), SQRT(tt*t1), SQRT(uu*u1)
*----------------------------------------------------------------------------
      WRITE(nout,'(a,5g20.12)') 'Virtual B(t)-B(u) = ',TBvirtDif1,TBvirtDif2,TBvirtDif1/TBvirtDif2
*----------------------------------------------------------------------------
      WRITE(nout,'(a,5g20.12)') 'Real B(t)-B(u) = '  , TBrealDif1,TBrealDif2,TBrealDif1/TBrealDif2
*----------------------------------------------------------------------------
      WRITE(nout,'(a,9g19.12)') ' CBoxGG,CBoxGZ = ', CBoxGG,CBoxGZ
****     $                  , DREAL(CBoxGG)/DREAL(CBoxGZ), DIMAG(CBoxGG)/DIMAG(CBoxGZ)
      WRITE(nout,'(a,9g20.12)') ' CBoxGG again ',DREAL(CBoxGG)/(ReBox-ReIR), DIMAG(CBoxGG)/ImBox
      WRITE(nout,'(a,9g20.12)') ' alfpi*DelBox    =  ', alfpi*DelBox
*----------------------------------------------------------------------------
      WRITE(nout,'(a,5g20.12)') 'TBvirt,TBvirt2 = ',   TBvirt, TBvirt2, TBvirt/TBvirt2
*----------------------------------------------------------------------------
      WRITE(nout,'(a,5g20.12)') 'A funct= ', BVR_CnuA(ss,m1,m2), sqrt(ss)
*----------------------------------------------------------------------------
*----------------------------------------------------------------------------
      WRITE(nout,'(a,5g20.12)')   'SBvirt= ',          SBvirt
      WRITE(nout,'(a,5g20.12)')   'SForFac= ',          Yfsr
*----------------------------------------------------------------------------
*----------------------------------------------------------------------------
      WRITE(nout,'(3a)') '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',
     $                              '===============',
     $                   '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      icont=icont+1
      END


      SUBROUTINE STest_Rules(mout,nmax)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   test of GPS rules for Weyl spinors                                            //
*//   does not use events                                                           //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER    mout,nmax
      REAL*8     x,y,z
      INTEGER    i,j,k,icont
      REAL*8     pi
      PARAMETER( pi=3.1415926535897932d0)

      COMPLEX*16 uup(4),u1(4),u2(4)
      COMPLEX*16 Uzplus(4), Uzminu(4),Umasiv(4), Cfactor
      REAL*8     vec(4),zeta(4),veta(4),pmom(4),phat(4)
      REAL*8     thet,phi,chi,eta,exe
      REAL*8     KinLib_AngPhi, beta, mass, tpz
      INTEGER    nout
*-------------*----------*----------------*-------------
      SAVE icont
      DATA icont/0/
*-------------*----------*----------------*-------------
      IF(icont .GE. nmax) RETURN
      nout = mout
cccc  nout = 6
*     ****************
*     Body of the test
*     ****************
*     Printouts
      WRITE(nout,'(a)') '  '
      WRITE(nout,'(3a)') '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',
     $                              ' STest_Rules ',
     $                   '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      icont=icont+1

      CALL Weyl_Initialize

      DO i=1,4
         uup(i)    = DCMPLX(0d0)
         Uzplus(i) = DCMPLX(0d0)
         Uzminu(i) = DCMPLX(0d0)
         Umasiv(i) = DCMPLX(0d0)
         zeta(i)   =0d0
         veta(i)   =0d0
         pmom(i)   =0d0
      ENDDO
      uup(1)    = DCMPLX(1d0)
      Uzplus(2) = DCMPLX(1d0)
      Uzminu(3) = DCMPLX(1d0)
*
      zeta(3)   = -1d0
      zeta(4)   =  1d0
      veta(1)   =  1d0
*
      mass      =  5d0
      pmom(4)   =  mass

      CALL  Weyl_Print(nout,' uup=   ',uup)
      CALL  Weyl_Print(nout,'Uzplus= ',Uzplus)
      CALL  Weyl_Print(nout,'Uzminu= ',Uzminu)


      CALL KinLib_VecPrint(nout,' zeta=  ',zeta)
      CALL KinLib_VecPrint(nout,' veta=  ',veta)
      CALL KinLib_VecPrint(nout,' pmom=  ',pmom)

      chi = 1.43d0
      exe = EXP(chi)
      thet=     pi*0.393d0
      phi = 2d0*pi*0.221d0
      CALL KinLib_Boost(3,   exe,  pmom,pmom)
      CALL KinLib_Rotor(3,1, thet, pmom,pmom)
      CALL KinLib_Rotor(1,2, phi,  pmom,pmom)
      CALL KinLib_VecPrint(nout,' pmom=  ',pmom)

      tpz = pmom(4)+pmom(3)
      Cfactor = DCMPLX(1d0/DSQRT(2*tpz))
* massive u(p,plus)
      Umasiv(1) = Cfactor *DCMPLX(tpz)
      Umasiv(2) = Cfactor *DCMPLX(pmom(1),pmom(2))
      Umasiv(3) = Cfactor *DCMPLX(mass)
      Umasiv(4) = Cfactor *DCMPLX(0d0)
* massive u(p,minus)
c      Umasiv(1) = Cfactor *DCMPLX(0d0)
c      Umasiv(2) = Cfactor *DCMPLX(mass)
c      Umasiv(3) = Cfactor *DCMPLX(-pmom(1),pmom(2))
c      Umasiv(4) = Cfactor *DCMPLX(tpz)
* massive v(p,plus)
c      Umasiv(1) = Cfactor *DCMPLX(0d0)
c      Umasiv(2) = Cfactor *DCMPLX(-mass)
c      Umasiv(3) = Cfactor *DCMPLX(-pmom(1),pmom(2))
c      Umasiv(4) = Cfactor *DCMPLX(tpz)
* massive v(p,minus)
c      Umasiv(1) = Cfactor *DCMPLX(tpz)
c      Umasiv(2) = Cfactor *DCMPLX(pmom(1),pmom(2))
c      Umasiv(3) = Cfactor *DCMPLX(-mass)
      Umasiv(4) = Cfactor *DCMPLX(0d0)

      CALL  Weyl_Print(nout,'Umasiv= ',Umasiv)

cc      phi = pi
cc      CALL  Weyl_Rotor(2, phi, uup, u1)
cc      CALL  Weyl_Print(nout,' u1=    ',u1)

      WRITE(nout,*) ' start '
*----------------------------------------------------------------------------------
      beta = KinLib_AngPhi( pmom(1), pmom(2))
      CALL KinLib_Rotor(1,2, -beta,  pmom,pmom) ! z-rotation
      CALL KinLib_Rotor(1,2, -beta,  zeta,zeta)
      CALL KinLib_Rotor(1,2, -beta,  veta,veta)
      CALL Weyl_Rotor(3, -beta, Uzplus, Uzplus)
      CALL Weyl_Rotor(3, -beta, Uzminu, Uzminu)
      CALL Weyl_Rotor(3, -beta, Umasiv, Umasiv)
              WRITE(nout,*) '/////////////////////////////// R3 puts pmom in x-z plane'
              CALL KinLib_VecPrint(nout,' pmom=  ',pmom)
              CALL KinLib_VecPrint(nout,' zeta=  ',zeta)
              CALL Weyl_Print(nout,'Uzplus= ',Uzplus)
              CALL Weyl_Print(nout,'Uzminu= ',Uzminu)
              CALL Weyl_Print(nout,'Umasiv= ',Umasiv)
*-----------------------------------------------------------------------------------
      beta = KinLib_AngPhi( pmom(3), pmom(1))
      CALL KinLib_Rotor(3,1, -beta, pmom,pmom) ! y-rotation
      CALL KinLib_Rotor(3,1, -beta, zeta,zeta)
      CALL KinLib_Rotor(3,1, -beta, veta,veta)
      CALL Weyl_Rotor(2, -beta, Uzplus, Uzplus)
      CALL Weyl_Rotor(2, -beta, Uzminu, Uzminu)
      CALL Weyl_Rotor(2, -beta, Umasiv, Umasiv)
              WRITE(nout,*) '/////////////////////////////// R2 puts pmom on z-axis'
              CALL KinLib_VecPrint(nout,' pmom=  ',pmom)
              CALL KinLib_VecPrint(nout,' zeta=  ',zeta)
              CALL Weyl_Print(nout,'Uzplus= ',Uzplus)
              CALL Weyl_Print(nout,'Uzminu= ',Uzminu)
              CALL Weyl_Print(nout,'Umasiv= ',Umasiv)
*----------------------------------------------------------------------------------
      beta = SQRT((pmom(4)+pmom(3))/(pmom(4)-pmom(3)))
      CALL KinLib_Boost(3, 1d0/beta,  pmom,pmom) ! z-boost
      CALL KinLib_Boost(3, 1d0/beta,  zeta,zeta)
      CALL KinLib_Boost(3, 1d0/beta,  veta,veta)
      beta = LOG(beta)
      CALL Weyl_Boost(3, -beta, Uzplus, Uzplus)
      CALL Weyl_Boost(3, -beta, Uzminu, Uzminu)
      CALL Weyl_Boost(3, -beta, Umasiv, Umasiv)
              WRITE(nout,*) '////////////////////////////// B3 puts pmom in rest frame'
              CALL KinLib_VecPrint(nout,' pmom=  ',pmom)
              CALL KinLib_VecPrint(nout,' zeta=  ',zeta)
              CALL Weyl_Print(nout,'Uzplus= ',Uzplus)
              CALL Weyl_Print(nout,'Uzminu= ',Uzminu)
              CALL Weyl_Print(nout,'Umasiv= ',Umasiv)
*----------------------------------------------------------------------------------
      beta = KinLib_AngPhi( zeta(3), zeta(1))
      CALL KinLib_Rotor(3,1, -pi-beta, pmom,pmom) ! y-rotation
      CALL KinLib_Rotor(3,1, -pi-beta, zeta,zeta)
      CALL KinLib_Rotor(3,1, pi-beta, veta,veta)
      CALL Weyl_Rotor(2, -pi-beta, Uzplus, Uzplus)
      CALL Weyl_Rotor(2, -pi-beta, Uzminu, Uzminu)
      CALL Weyl_Rotor(2, -pi-beta, Umasiv, Umasiv)
              WRITE(nout,*) '///////////////////////////// R2 puts zeta on z-axis'
              CALL KinLib_VecPrint(nout,' pmom=  ',pmom)
              CALL KinLib_VecPrint(nout,' zeta=  ',zeta)
              CALL KinLib_VecPrint(nout,' veta=  ',veta)
              CALL Weyl_Print(nout,'Uzplus= ',Uzplus)
              CALL Weyl_Print(nout,'Uzminu= ',Uzminu)
              CALL Weyl_Print(nout,'Umasiv= ',Umasiv)
*---------------------------------------------------------------------------------
      beta = KinLib_AngPhi( veta(1), veta(2))
      CALL KinLib_Rotor(1,2, -beta,  zeta,zeta)
      CALL KinLib_Rotor(1,2, -beta,  veta,veta)
      CALL Weyl_Rotor(3, -beta, Uzplus, Uzplus)
      CALL Weyl_Rotor(3, -beta, Uzminu, Uzminu)
      CALL Weyl_Rotor(3, -beta, Umasiv, Umasiv)
              WRITE(nout,*) '///////////////////////////// R3 puts veta in x-z plane'
              CALL KinLib_VecPrint(nout,' zeta=  ',zeta)
              CALL KinLib_VecPrint(nout,' veta=  ',veta)
              CALL Weyl_Print(nout,'Uzplus= ',Uzplus)
              CALL Weyl_Print(nout,'Uzminu= ',Uzminu)
              CALL Weyl_Print(nout,'Umasiv= ',Umasiv)
      WRITE(nout,'(a,g20.13)') ' last 1-beta/2pi = ',1-beta/2/pi
      END                       !!!!! STest_Rules



      SUBROUTINE STest_NewSingle(nout,nmax)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   See README  for notes on the tests.                                           //
*//                                                                                 //
*//   test of FSR single bremss matrix element                                      //
*//   compared are also ERW subprograms                                             //
*//                                                                                 //
*//   Validity range:                                                               //
*//       Z may be ON --- still some problem for keyZet=1 with ERW subprograms.     //
*//       Energies should be far from the threshold                                 //
*//       KeyFSR=1!!!                                                               //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER  nout,nmax
c[[[      REAL*8   p1(4),p2(4),p3(4),p4(4),ph(4),h1(3),h2(3)
      REAL*8   p1(4),p2(4),p3(4),p4(4),ph(4),h1(4),h2(4)
      REAL*8   yy,zz,y,z,xk,vv
      COMPLEX*16 Amp1fsr(2,2,2,2,2),tensor(2,2,2,2,2,2,2,2)
      COMPLEX*16 Amp1Phot(2,2,2,2,2),Amp1isr(2,2,2,2,2)
      INTEGER  i,j,k,icont,idv
      INTEGER  KFi,KFf,nphot
      REAL*8   svar,svar1,amel,amfin,delq,ratio
      REAL*8   KLsect,Xsect3,Dig1,Xdist,Xdist0
      REAL*8   BornV_GetMass
      REAL*8   wm0,wmd,wm1
      REAL*8   wt,WtMain,WtCrud
      REAL*8 oldist(4), arra,arra1,arra2,fako,cykus
*-------------
      SAVE icont,arra
      DATA icont/0/
      DATA arra/0d0/
*----------------------------------------------------
*      IF(icont .GE. nmax) RETURN
*----------------------------------------------------
      CALL KarLud_GetBeams(    p1,p2)
      CALL KarFin_GetFermions( p3,p4)
      CALL KarFin_GetNphot(nphot)     
* check photon multiplicity
*---------------------------------------
      IF( nphot .NE. 1 ) RETURN
*---------------------------------------
      KFi = 11                ! KF=11 is electron
      CALL KarLud_GetKFfin(KFf)     ! Actual KFcode of final fermion
      CALL KarFin_GetSudakov1( 1,yy,zz)
      y  = yy/(1 +yy+zz)
      z  = zz/(1 +yy+zz)
      CALL KarFin_GetPhoton1(  1,ph)
      amfin   = BornV_GetMass(KFf)
      svar  = (p1(4)+p2(4))**2
      svar1 = (p3(4)+p4(4))**2-(p3(3)+p4(3))**2-(p3(2)+p4(2))**2-(p3(1)+p4(1))**2
      vv = 1-svar1/svar
      delq=  amfin**2/svar1
      wm0= 1d0 -2d0*delq -delq*(y/z+z/y)
      wmd= 1d0 + delq*( y/z +z/y ) *(y**2+z**2) / ((1-y)**2+(1-z)**2)
      wm1= 1d0 - delq*( y/z+z/y ) *(1d0-y)*(1d0-z) *2d0/((1d0-y)**2+(1d0-z)**2)
*------------------------------------------------------------------
cc      IF( (vv .LT. 0.2d0) .OR. (vv .GT. 0.9d0) ) RETURN
CC      IF( vv .LT. 0.8d0 ) RETURN
        IF( (ph(1)**2+ph(2)**2).lt.0.0001) return
        IF( vv .LT. 0.01d0 ) RETURN
*------------------------------------------------------------------
********************************
c      CALL  GPS_SetKeyArb(1)
c      CALL GPS_Setb2
c      CALL GPS_Setb3

      CALL Brem_Porownanie(KFi,KFf,p1,p2,p3,p4,ph,oldist)
*     ==================================================

      CALL STest_Make1isr(KFi,KFf,p1,p2,p3,p4,ph,Amp1isr,Xdist,Xdist0) ! OK
      CALL STest_Make1fsr(KFi,KFf,p1,p2,p3,p4,ph,Amp1fsr,Xdist,Xdist0) ! OK

      CALL STest_MakeTensor(Amp1isr,Amp1fsr,tensor) ! OK

      CALL STest_Make2(KFi,KFf,p1,p2,p3,p4,tensor,Xdist) ! OK
***   CALL GPS_Amp1Print(nout,'Amp1fsr',Amp1fsr)

      wt = Xdist/Xdist0
**
**      CALL KK2f_GetWt(WtMain,WtCrud) !! events shoud be weighted !!
      idv  = 50000
      CALL GLK_Fil1(idv+19, wt*wm0, 1d0)

      CALL QED3_GetDig1(Dig1)
      ratio = Xdist/Dig1*(svar/2d0)
      CALL ERW_initialize       !!!
      CALL ERW_SklFin(KLsect)
      CALL ERW_SingFin(Xsect3)  !!!
********************************
cc      IF( ABS(ratio-1d0)  .LT. 1d-4)  RETURN
*------------------------------------------------------------------
      fako=(Xdist/oldist(2)-1)**2! +(Xdist/Xsect3-1d0)**2!+
      if (fako.gt.arra) then
      arra=max(fako,arra)
!      arra=0.01**2
      arra1=Xdist/oldist(2)
      arra2=Xdist/oldist(3)
      WRITE(nout,'(a)') '  '
      WRITE(nout,'(3a)') '///////////////////////////////////////////////',
     $                                  ' STest_NewSingle ',
     $                   '///////////////////////////////////////////////'
      WRITE(nout,'(a,3f20.14)') 'wm0,wmd,wm1 = ',wm0,wmd,wm1

c[[[[[
cc      CALL Brem_SpinStore(1,h1,h2)  ! original code
c]]]]]
      CALL STest_GetHvectors(h1,h2)

      WRITE(nout,*) 'h1= ',h1
      WRITE(nout,*) 'h2= ',h2

      CALL KinLib_VecPrint(nout,'p1=      ',p1)
      CALL KinLib_VecPrint(nout,'p2=      ',p2)
      CALL KinLib_VecPrint(nout,'p3=      ',p3)
      CALL KinLib_VecPrint(nout,'p4=      ',p4)
      CALL KinLib_VecPrint(nout,'ph=      ',ph)
      WRITE(nout,*) 
     $     'p3+p4=',sqrt((p3(4)+p4(4))**2-(p3(3)+p4(3))**2-(p3(2)+p4(2))**2-(p3(1)+p4(1))**2)
      WRITE(nout,*) 
     $     'ph+p4=',sqrt((ph(4)+p4(4))**2-(ph(3)+p4(3))**2-(ph(2)+p4(2))**2-(ph(1)+p4(1))**2)
      WRITE(nout,*) 
     $     'p3+ph=',sqrt((p3(4)+ph(4))**2-(p3(3)+ph(3))**2-(p3(2)+ph(2))**2-(p3(1)+ph(1))**2)
      WRITE(nout,'(a,1f18.12,3e20.7)') 
     $     'y+z, y/(z*delq), z/(y*delq) = ',y+z, y/(z*delq), z/(y*delq)
      WRITE(nout,'(a,5g20.14)') ' Dig1=                ', Dig1
      WRITE(nout,'(a,5g20.14)') ' Make1fsr: Xdist,Xdist0,X/X0/2= ', Xdist,Xdist0,Xdist/Xdist0/2
* svar factor due to normalization of s-factors
      WRITE(nout,'(a,5g20.14)') ' ratio = Xdist/Dig1*(svar/2d0)            = ', ratio
      WRITE(nout,'(a,5g20.14)') ' Xsect3, Xdist/Xsect3    = ',  Xsect3, Xdist/Xsect3
      WRITE(nout,'(a,5g20.14)') '2*KLsect, Xdist/(2*KLsect)=',2*KLsect, Xdist/(2*KLsect)
      WRITE(nout,'(a,5g20.14)') 'oldbrem ratio 1,2/Xdist   =',xdist/oldist(1),xdist/oldist(2)
      WRITE(nout,'(a,5g20.14)') 'oldbrem ratio 3,4/Xdist   =',xdist/oldist(3),xdist/oldist(4)

      ENDIF
*------------------
      icont=icont+1

      END                       !STest_NewSingle


      SUBROUTINE STest_MakeTensor(Amp1isr,Amp1fsr,tensor)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//              TESTS ONLY                                                         //
*//                                                                                 //
*//   Input:                                                                        //
*//   Amp1isr   = 32 spin amplitudes                                                //
*//   Amp1fsr   = 32 spin amplitudes                                                //
*//                                                                                 //
*//   Output:                                                                       //
*//   Spin density tensor for (ini+fin) summed/averaged over photon spin            //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
*
      COMPLEX*16 Amp1fsr( 2,2,2,2,2)
      COMPLEX*16 Amp1isr( 2,2,2,2,2),tensor(2,2,2,2,2,2,2,2)
*
      INTEGER    j1,j2,j3,j4,i1,i2,i3,i4
*--------------------------------------
         DO j1=1,2
            DO i1=1,2
               DO j2=1,2
                  DO i2=1,2
                     DO j3=1,2
                        DO i3=1,2
                           DO j4=1,2
                              DO i4=1,2
                              tensor(i1,j1,i2,j2,i3,j3,i4,j4)=
     $                              (Amp1fsr(i1,i2,i3,i4,1)+Amp1isr(i1,i2,i3,i4,1))*
     $                        DCONJG(Amp1fsr(j1,j2,j3,j4,1)+Amp1isr(j1,j2,j3,j4,1))+
     $                              (Amp1fsr(i1,i2,i3,i4,2)+Amp1isr(i1,i2,i3,i4,2))*
     $                        DCONJG(Amp1fsr(j1,j2,j3,j4,2)+Amp1isr(j1,j2,j3,j4,2))
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
      END


      SUBROUTINE STest_Make2(KFi,KFf,p1,p2,p3,p4,tensor,Xdist)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   !!!!!!!!! initial state density matrix should be transposed !!!!!!            //
*//   !!!!!!!!! initial state density matrix should be transposed !!!!!!            //
*//   !!!!!!!!! initial state density matrix should be transposed !!!!!!            //
*//                                                                                 //
*//          FOR TESTS ONLY                                                         //
*//                                                                                 //
*//   calculates x-section from spin tensor using spinor techniques.                //
*//   Complete spin effects of fermions are taken into account                      //
*//                                                                                 //
*//   Input:                                                                        //
*//   KFi, Kff = beam and final fermion flavour codes (to define charges)           //
*//   pi,mi    are for spinors, not for gamma and Z propagators                     //
*//   p1,m1    =fermion momentum and mass (beam)                                    //
*//   p2,m2    =fermion momentum and mass (beam)                                    //
*//   p3,m3    =fermion momentum and mass final state                               //
*//   p4,m4    =fermion momentum and mass final state                               //
*//   ph       =photon  momentum                                                    //
*//                                                                                 //
*//   Output:                                                                       //
*//   Xdist                                                                         //
*//                                                                                 //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
*
      INTEGER    KFi,KFf
      REAL*8     p1(4),p2(4),p3(4),p4(4),Xdist
      COMPLEX*16 tensor(2,2,2,2,2,2,2,2)
      COMPLEX*16 DensityMe1(2,2), DensityMe2(2,2), DensityMf1(2,2), DensityMf2(2,2)
*
      INTEGER    j1,j2,j3,j4,i1,i2,i3,i4
      REAL*8     Sum
      INTEGER    BornV_GetColor
      COMPLEX*16 ccsum
*--------------------------------------

      CALL STest_SpinGiveBeam(1,DensityMe1)
      CALL STest_SpinGiveBeam(2,DensityMe2)

      CALL GPS_tralorPrepare(p3,1) !!!
      CALL GPS_tralorPrepare(p4,2) !!!

      CALL STest_SpinGive(1,DensityMf1)
      CALL STest_SpinGive(2,DensityMf2)

      ccSum=0d0
      DO j1=1,2
         DO i1=1,2
            DO j2=1,2
               DO i2=1,2
                  DO j3=1,2
                     DO i3=1,2
                        DO j4=1,2
                           DO i4=1,2
                              ccSum=ccSum+ tensor(i1,j1,i2,j2,i3,j3,i4,j4)
     $                                            *DensityMe1(j1,i1) !!!<-- should be (i1,j1)
     $                                            *DensityMe2(j2,i2) !!!<-- should be (i2,j2)
     $                                            *DensityMf1(j3,i3) ! OK.
     $                                            *DensityMf2(j4,i4) ! OK.
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      sum=ccsum
* Color factor
      Xdist = Sum *BornV_GetColor(KFf)
      END




      SUBROUTINE STest_SpinGive(ID,DensityM)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//  calculates 2x2 density matrix for a single fs fermion.                         //
*//  polarization vector is SET by  CALL  STest_SetHvectors                         //
*//  and the appropriate boosts rotations are performed                             //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE "STest.h"
*
      INTEGER      ID,K,i,j
      REAL*8       polar(4)
      COMPLEX*16   DensityM(2,2)
*-------------------------------------------------------------------------------------

      IF(ID.EQ.1) THEN
         polar(1) = m_HvecFer1(1)
         polar(2) = m_HvecFer1(2)
         polar(3) = m_HvecFer1(3)
         polar(4) = m_HvecFer1(4)
      ELSE
         polar(1) = m_HvecFer2(1)
         polar(2) = m_HvecFer2(2)
         polar(3) = m_HvecFer2(3)
         polar(4) = m_HvecFer2(4)
      ENDIF

c{{{{{{ TEST TEST TEST TEST TEST TEST
c* Brem_Tralor of KORALB is performing internally alignement of the spin quatization
c* Axes of tau+ and tau- . Rotation of the angle pi.
c* It is Wigner rotation from JW2 (of KORALB) to GPS frame!
      polar(4) =0d0                ! <--  Apparently Brem_Tralor requires this
      CALL  Brem_Tralor(   id,polar)        ! JW2-->CMS
      CALL  GPS_TralorUnDo(id,polar,polar)  ! CMS-->GPS
      polar(4) =1d0
c}}}}}}

      DO i=1,2
         DO j=1,2
            DensityM(i,j)=0D0
            DO k=1,4
               DensityM(i,j)=DensityM(i,j)+m_Pauli4( k,i,j)*polar(k)
            ENDDO
         ENDDO
      ENDDO
      END



      SUBROUTINE STest_SpinGiveBeam(ID,DensityM)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              //
*//   Note that Wigner rotation for beam polarization can be done with              //
*//   KK2f_WignerIni(KFbeam,CMSene,PolBeam1,PolBeam2, Polar1,Polar2)                //
*//   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              //
*//                                                                                 //
*//  calculates 2x2 density matrix for a single fermion.                            //
*//  polarization vector is SET by  CALL  STest_SetPolBeams                         //
*//  and the appropriate  rotations/sign arrangements are performed                 //
*//  note sign difference withy respect to fs and rotation of e1 around 2-nd axis   //
*//  the two operations are not done in an explicit way                             //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE "STest.h"
      INTEGER ID,K,i,j
      REAL*8     polar(4)
      COMPLEX*16 DensityM(2,2)
*--------------------------------------------------------------------------------------
* ????????????????????????????????????????????????????
* This looks like transformation to left-handed frame
* ????????????????????????????????????????????????????
      IF(ID.EQ.1) THEN
         polar(1)  =  m_PolBeam1(1)
         polar(2)  = -m_PolBeam1(2)
         polar(3)  =  m_PolBeam1(3)
      ELSE
         polar(1)  = -m_PolBeam2(1)
         polar(2)  = -m_PolBeam2(2)
         polar(3)  = -m_PolBeam2(3)
      ENDIF
      polar(4)=1D0              ! Just to be sure

      DO i=1,2
         DO j=1,2
            DensityM(i,j)=0D0
            DO k=1,4
               DensityM(i,j)=DensityM(i,j)+m_Pauli4( k,i,j)*polar(k)
            ENDDO
         ENDDO
      ENDDO
      END


*/////////////////////////////////////////////////////////////////////////////////////
*/////////////////////////////////////////////////////////////////////////////////////
*/////////////////////////////////////////////////////////////////////////////////////
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//                                                                                 //
*//           Test subprograms moved here from GPS class                            //
*//                                                                                 //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
*/////////////////////////////////////////////////////////////////////////////////////
*/////////////////////////////////////////////////////////////////////////////////////
*/////////////////////////////////////////////////////////////////////////////////////


      SUBROUTINE STest_Make1isr(KFi,KFf,p1,p2,p3,p4,ph,Amp1isr,Xdist,Xdist0)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   ISR 1-photon spin amplitudes using spinor techniques.                         //
*//                                                                                 //
*//   Input:                                                                        //
*//   KFi, Kff = beam and final fermion flavour codes (to define charges)           //
*//   pi,mi    are for spinors, not for gamma and Z propagators                     //
*//   p1,m1    =fermion momentum and mass (beam)                                    //
*//   p2,m2    =fermion momentum and mass (beam)                                    //
*//   p3,m3    =fermion momentum and mass final state                               //
*//   p4,m4    =fermion momentum and mass final state                               //
*//   ph       =photon  momentum                                                    //
*//                                                                                 //
*//   Output:                                                                       //
*//   Amp1isr   = 32 spin amplitudes                                                //
*//   Xdist     = spin summed O(alf1)                                               //
*//   Xdist0    = spin summed O(alf0) crude, for test                               //
*//                                                                                 //
*//   Notes:                                                                        //
*//   X-checked at the level of unpolarized differential distribution               //
*//   Clear split into infrared part and the rest                                   //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
*
      INTEGER    KFi,KFf
      REAL*8     p1(4),p2(4),p3(4),p4(4),ph(4),Xdist,Xdist0
      REAL*8     PP(4),PX(4)
      COMPLEX*16 Amp1isr( 2,2,2,2,2)
*
      COMPLEX*16 AmBorIsr(2,2,2,2)
      COMPLEX*16 AmHarIsr(2,2,2,2,2)
*
      INTEGER    j,j1,j2,j3,j4,k,Sig,KeyArb
      REAL*8     Massf,Mbeam,Fleps,BetaFin,SfacIni,Sum,Sum0
      REAL*8     svar,svar1
      INTEGER    BornV_GetColor
      REAL*8     BornV_GetMass
      COMPLEX*16 GPS_soft,GPS_softb
      COMPLEX*16 Sini(2)
      REAL*8     m1,m2,m3,m4,mph
*--------------------------------------
      CALL GPS_Initialize
      Mbeam  =  BornV_GetMass(KFi)
      Massf =  BornV_GetMass(KFf)
      Fleps =  1d-100
      m1  = Mbeam
      m2  = Mbeam
      m3  = Massf
      m4  = Massf
      mph = Fleps
*
      DO k=1,4
         PP(k)=p1(k)+p2(k)
         PX(k)=p1(k)+p2(k)-ph(k)
      ENDDO
      svar  = PP(4)**2-PP(3)**2-PP(2)**2-PP(1)**2
      svar1 = PX(4)**2-PX(3)**2-PX(2)**2-PX(1)**2
      IF(svar1 .LE. 4*Massf**2) GOTO 900
* Soft factor for two photon helicities
      CALL GPS_GetKeyArb(KeyArb)
      SfacIni=0d0
      DO k=1,2
         Sig = 3-2*k
         IF( KeyArb .EQ. 0 ) THEN
            Sini(k)  = GPS_soft(  Sig,ph,p1,p2)
         ELSE
            Sini(k)  = GPS_softb( Sig,ph,p1,m1,p2,m2)
         ENDIF
         SfacIni =SfacIni +CDABS(Sini(k))**2
      ENDDO
***   WRITE(*,'(a,5g20.14)') 'Make1Phot: SfacIni*svar= ',SfacIni*svar/2d0
*----------------------------------------
* Calculate normal Born spin amplitudes,
* (However, fermion 4-momenta not restricted to within 2-body phase space!)
* Electron mass set to infinitesimaly small (sign is however significant!!!)
* The Born is associated with ISR s-factor, therefore PX (not PP).
      CALL GPS_Born(KFi,KFf, PX, p1,Fleps, p2,-Fleps, p3,m3, p4,-m4, AmBorIsr)
****  CALL GPS_BPrint(6,'AmBorIsr',AmBorIsr)
*----------------------------------------
* ISR non-infrared part
      CALL GPS_Hini(KFi,KFf,PX, p1,m1,p2,m2,p3,m3,p4,m4,ph,mph, AmHarIsr)
*----------------------------------
* Grand total is sum of soft and hard:
      Sum  =0d0
      Sum0 =0d0
      DO k=1,2
         DO j1=1,2
            DO j2=1,2
               DO j3=1,2
                  DO j4=1,2
                     Amp1isr(j1,j2,j3,j4,k) = AmBorIsr( j1,j2,j3,j4)*Sini(k)
     $                                       +AmHarIsr( j1,j2,j3,j4, k)
                     Sum =Sum+  CDABS( Amp1isr(j1,j2,j3,j4,k) )**2
                     Sum0=Sum0+ CDABS( AmBorIsr( j1,j2,j3,j4)*Sini(k) )**2
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
***   CALL GPS_Amp1Print(6,'AmHarIsr  ',AmHarIsr)
***   CALL GPS_Amp1Print(6,'Amp1isr   ',Amp1isr)
* Color factor
      Xdist  = Sum  *BornV_GetColor(KFf)
      Xdist0 = Sum0 *BornV_GetColor(KFf)
* What about phase space factor???
      RETURN
 900  Xdist  = 0d0
      Xdist0 = 0d0
      END                       !!! GPS_Make1isr

      SUBROUTINE STest_Make1fsr(KFi,KFf,p1,p2,p3,p4,ph,Amp1fsr,Xdist,Xdist0)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   FSR 1-photon spin amplitudes using spinor techniques.                         //
*//                                                                                 //
*//   Input:                                                                        //
*//   KFi, Kff = beam and final fermion flavour codes (to define charges)           //
*//   pi,mi    are for spinors, not for gamma and Z propagators                     //
*//   p1,m1    =fermion momentum and mass (beam)                                    //
*//   p2,m2    =fermion momentum and mass (beam)                                    //
*//   p3,m3    =fermion momentum and mass final state                               //
*//   p4,m4    =fermion momentum and mass final state                               //
*//   ph       =photon  momentum                                                    //
*//                                                                                 //
*//   Output:                                                                       //
*//   Amp1fsr   = 32 spin amplitudes                                                //
*//   Xdist     = spin summed O(alf1) form Amp1fsr                                  //
*//   Xdist0    = spin summed O(alf0) for crude MC, for tests of the weight         //
*//                                                                                 //
*//   X-checked at the level of unpolarized differential distribution               //
*//   Clear split into infrared part and the rest                                   //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
cccc      INCLUDE "STest.h"
*
      INTEGER    KFi,KFf
      REAL*8     p1(4),p2(4),p3(4),p4(4),ph(4),Xdist,Xdist0
      REAL*8     PP(4),PX(4)
      COMPLEX*16 Amp1fsr( 2,2,2,2,2)
*
      COMPLEX*16 AmBorFsr(2,2,2,2)
      COMPLEX*16 AmHarFsr(2,2,2,2,2)
      COMPLEX*16 AmpBornU(2,2,2,2)
      COMPLEX*16 AmpBornV(2,2,2,2)
*
      INTEGER    j,j1,j2,j3,j4,k,Sig,KeyArb
      REAL*8     Massf,Mbeam,Fleps,BetaFin
      REAL*8     svar,svar1
      INTEGER    BornV_GetColor
      REAL*8     BornV_GetMass
      REAL*8     pr1,pr2, SfacFin
      COMPLEX*16 GPS_soft,GPS_softb
      COMPLEX*16 Csum1,Csum2,U(2,2),V(2,2)
      COMPLEX*16 Sfin(2)
      REAL*8     m1,m2,m3,m4,mph
      REAL*8     Sum,Sum0
*--------------------------------------
      CALL GPS_Initialize
      Mbeam  =  BornV_GetMass(KFi)
      Massf =  BornV_GetMass(KFf)
      Fleps =  1d-100
      m1  = Mbeam
      m2  = Mbeam
      m3  = Massf
      m4  = Massf
      mph = Fleps
*
      DO k=1,4
         PP(k)=p1(k)+p2(k)
         PX(k)=p1(k)+p2(k)-ph(k)
      ENDDO
      svar  = PP(4)**2-PP(3)**2-PP(2)**2-PP(1)**2
      svar1 = PX(4)**2-PX(3)**2-PX(2)**2-PX(1)**2
      IF(svar1 .LE. 4*Massf**2) GOTO 900
*
* Soft part is easy...
      CALL GPS_GetKeyArb(KeyArb)
      SfacFin=0d0
      DO k=1,2
         Sig = 3-2*k
         IF( KeyArb .EQ. 0 ) THEN
            Sfin(k) = -GPS_soft(  Sig,ph,p3,p4)
         ELSE
            Sfin(k) = -GPS_softb( Sig,ph,p3,m3,p4,m4)
         ENDIF
         SfacFin = SfacFin +CDABS(Sfin(k))**2
      ENDDO
***   WRITE(*,'(a,5g20.14)') 'Make1Phot: SfacFin*svar= ',SfacFin*svar/2d0
*----------------------------------------
* Calculate normal Born spin amplitudes,
* (However, fermion 4-momenta not restricted to within 2-body phase space!)
* Electron mass set to infinitesimaly small (sign is however significant!!!)
* The Born is multiplied by FSR s-factor, therefore PP (not PX).
      CALL GPS_Born(KFi,KFf, PP, p1,Fleps, p2,-Fleps, p3,m3, p4,-m4, AmBorFsr)
*----------------------------------------
      CALL GPS_Hfin(KFi,KFf,PP, p1,m1,p2,m2,p3,m3,p4,m4,ph,mph, AmHarFsr)
*----------------------------------
* Grand total is sum of soft and hard:
      Sum =0d0
      Sum0=0d0
      DO k=1,2
         DO j1=1,2
            DO j2=1,2
               DO j3=1,2
                  DO j4=1,2
                     Amp1fsr(j1,j2,j3,j4,k) = 
     $                    +AmBorFsr(j1,j2,j3,j4)*Sfin(k)
     $                    +AmHarFsr(j1,j2,j3,j4, k)
                     Sum =Sum + CDABS( Amp1fsr( j1,j2,j3,j4,k) )**2
                     Sum0=Sum0+ CDABS( AmBorFsr(j1,j2,j3,j4)*Sfin(k)*svar/svar1 )**2
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
***   CALL GPS_BPrint(   6,'AmBorFsr',AmBorFsr)
***   CALL GPS_Amp1Print(6,'AmHarFsr',AmHarFsr)
***   CALL GPS_Amp1Print(6,'Amp1fsr ',Amp1fsr)
* Color factor
      Xdist  = Sum  *BornV_GetColor(KFf)
      Xdist0 = Sum0 *BornV_GetColor(KFf)
* What about phase space factor
      RETURN
 900  Xdist  = 0d0
      Xdist0 = 0d0
      END                       !!!! GPS_Make1fsr



      SUBROUTINE STest_Make1pho(KFi,KFf,p1,p2,p3,p4,ph,Amp1pho,Xdist,Xdist0)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   ISR+FSR 1-photon spin amplitudes using spinor techniques.                     //
*//                                                                                 //
*//   Input:                                                                        //
*//   KFi, Kff = beam and final fermion flavour codes (to define charges)           //
*//   pi,mi    are for spinors, not for gamma and Z propagators                     //
*//   p1,m1    =fermion momentum and mass (beam)                                    //
*//   p2,m2    =fermion momentum and mass (beam)                                    //
*//   p3,m3    =fermion momentum and mass final state                               //
*//   p4,m4    =fermion momentum and mass final state                               //
*//   ph       =photon  momentum                                                    //
*//                                                                                 //
*//   Output:                                                                       //
*//   Amp1pho   = 32 spin amplitudes                                                //
*//                                                                                 //
*//   Notes:                                                                        //
*//                                                                                 //
*//   UNTESTED?????                                                                 //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
ccc      INCLUDE "STest.h"
*
      INTEGER    KFi,KFf
      REAL*8     p1(4),p2(4),p3(4),p4(4),ph(4),Xdist,Xdist0
      REAL*8     PP(4),PX(4)
      COMPLEX*16 Amp1pho( 2,2,2,2,2)
*
      COMPLEX*16 AmBorIsr(2,2,2,2),   AmBorFsr(2,2,2,2)
      COMPLEX*16 AmHarIsr(2,2,2,2,2), AmHarFsr(2,2,2,2,2)
*
      INTEGER    j,j1,j2,j3,j4,k,Sig,KeyArb
      REAL*8     Massf,Mbeam,Fleps,BetaFin,Sfac,Sum,Sum0
      REAL*8     svar,svar1
      INTEGER    BornV_GetColor
      REAL*8     BornV_GetMass
      COMPLEX*16 GPS_soft,GPS_softb
      COMPLEX*16 Sini(2),Sfin(2)
      REAL*8     m1,m2,m3,m4,mph
*--------------------------------------
      CALL GPS_Initialize
      Mbeam  =  BornV_GetMass(KFi)
      Massf =  BornV_GetMass(KFf)
      Fleps =  1d-100
      m1  = Mbeam
      m2  = Mbeam
      m3  = Massf
      m4  = Massf
      mph = Fleps
*
      DO k=1,4
         PP(k)=p1(k)+p2(k)
         PX(k)=p1(k)+p2(k)-ph(k)
      ENDDO
      svar  = PP(4)**2-PP(3)**2-PP(2)**2-PP(1)**2
      svar1 = PX(4)**2-PX(3)**2-PX(2)**2-PX(1)**2
      IF(svar1 .LE. 4*Massf**2) GOTO 900
* Soft factor for two photon helicities
      CALL GPS_GetKeyArb(KeyArb)
      Sfac=0d0
      DO k=1,2
         Sig = 3-2*k
         IF( KeyArb .EQ. 0 ) THEN
            Sini(k)  = GPS_soft(  Sig,ph,p1,p2)
            Sfin(k) = -GPS_soft(  Sig,ph,p3,p4)
         ELSE
            Sini(k)  = GPS_softb( Sig,ph,p1,m1,p2,m2)
            Sfin(k) = -GPS_softb( Sig,ph,p3,m3,p4,m4)
         ENDIF
         Sfac =Sfac +CDABS(Sini(k))**2
      ENDDO
*----------------------------------------
* Calculate normal Born spin amplitudes,
* (However, fermion 4-momenta not restricted to within 2-body phase space!)
* Electron mass set to infinitesimaly small (sign is however significant!!!)
* The Born is associated with ISR has PX, with FSR has PP.
      CALL GPS_Born(KFi,KFf, PX, p1,Fleps, p2,-Fleps, p3,m3, p4,-m4, AmBorIsr)
      CALL GPS_Born(KFi,KFf, PP, p1,Fleps, p2,-Fleps, p3,m3, p4,-m4, AmBorFsr)
*----------------------------------------
* ISR and FSR non-infrared parts
      CALL GPS_Hini(KFi,KFf,PX, p1,m1,p2,m2,p3,m3,p4,m4,ph,mph, AmHarIsr)
      CALL GPS_Hfin(KFi,KFf,PP, p1,m1,p2,m2,p3,m3,p4,m4,ph,mph, AmHarFsr)
*----------------------------------
* Grand total is sum of soft and hard:
      Sum  = 0d0
      Sum0 = 0d0
      DO k=1,2
         DO j1=1,2
            DO j2=1,2
               DO j3=1,2
                  DO j4=1,2
                     Amp1pho(j1,j2,j3,j4,k) = AmBorIsr( j1,j2,j3,j4)*Sini(k)
     $                                       +AmBorFsr( j1,j2,j3,j4)*Sfin(k)
     $                                       +AmHarIsr( j1,j2,j3,j4, k)
     $                                       +AmHarFsr( j1,j2,j3,j4, k)
                     Sum  =Sum + CDABS(Amp1pho(j1,j2,j3,j4,k))**2
                     Sum0 =Sum0+ CDABS( AmBorIsr( j1,j2,j3,j4)*Sini(k)
     $                                 +AmBorFsr( j1,j2,j3,j4)*Sfin(k)*svar/svar1 
     $                                )**2
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
****  WRITE(*,'(a,5g20.14)') 'Make1Phot: Sfac*svar= ',Sfac*svar/2d0
****  CALL GPS_BPrint(6,'AmBorIsr',AmBorIsr)
****  CALL GPS_BPrint(6,'AmBorFsr',AmBorFsr)
****  CALL GPS_Amp1Print(6,'AmHarIsr  ',AmHarIsr)
****  CALL GPS_Amp1Print(6,'AmHarFsr  ',AmHarFsr)
****  CALL GPS_Amp1Print(6,'Amp1pho   ',Amp1pho)
* Color factor
      Xdist  = Sum  *BornV_GetColor(KFf)
      Xdist0 = Sum0 *BornV_GetColor(KFf)
* What about phase space factor???
      RETURN
 900  Xdist  = 0d0
      Xdist0 = 0d0
      END                       !!! GPS_Make1pho

*/////////////////////////////////////////////////////////////////////////////////////
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//           Setters and Getters of CLASS  STest                                   //
*//           Setters and Getters of CLASS  STest                                   //
*//           Setters and Getters of CLASS  STest                                   //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
*/////////////////////////////////////////////////////////////////////////////////////



      SUBROUTINE STest_SetPolBeams(PolBeam1,PolBeam2)
*/////////////////////////////////////////////////////////////////////////////////////
*//   Seting beam POLARIZATION vectors                                              //
*//   Dont forget Wigner rotation to GPS frame!!!!                                  //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE "STest.h"
      INTEGER i,j,k
      REAL*8 PolBeam1(4),PolBeam2(4)
*------------------------------------------------------------------------------------
      DO k=1,4
         m_PolBeam1( k) = PolBeam1(k)
         m_PolBeam2( k) = PolBeam2(k)
      ENDDO
* Define spin density matriced
      DO i=1,2
         DO j=1,2
            m_SDMat1(i,j)=0D0
            m_SDMat2(i,j)=0D0
            DO k=1,4
               m_SDMat1(i,j)=m_SDMat1(i,j)+m_Pauli4( k,i,j) *PolBeam1(k)
               m_SDMat2(i,j)=m_SDMat2(i,j)+m_Pauli4( k,i,j) *PolBeam2(k)
            ENDDO
         ENDDO
      ENDDO
      END

      SUBROUTINE STest_SetHvectors(HvecFer1,HvecFer2)
*/////////////////////////////////////////////////////////////////////////////////////
*//   Seting final fermion POLARIMETER vectors                                      //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE "STest.h"
      INTEGER i,j,k
      REAL*8 HvecFer1(4),HvecFer2(4)
*-------------------------------------------------------------------------------------
      DO k=1,4
         m_HvecFer1( k) = HvecFer1(k)
         m_HvecFer2( k) = HvecFer2(k)
      ENDDO
* Define immediately polarimeter density matriced
      DO i=1,2
         DO j=1,2
            m_SDMat3(i,j)=0D0
            m_SDMat4(i,j)=0D0
            DO k=1,4
               m_SDMat3(i,j)=m_SDMat3(i,j)+m_Pauli4( k,i,j) *HvecFer1(k)
               m_SDMat4(i,j)=m_SDMat4(i,j)+m_Pauli4( k,i,j) *HvecFer2(k)
            ENDDO
         ENDDO
      ENDDO
      END


      SUBROUTINE STest_GetPolBeams(PolBeam1,PolBeam2)
*/////////////////////////////////////////////////////////////////////////////////////
*//   Geting beam POLARIZATION vectors                                              //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE "STest.h"
      INTEGER k
      REAL*8 PolBeam1(4),PolBeam2(4)
*
      DO k=1,4
         PolBeam1( k) = m_PolBeam1(k)
         PolBeam2( k) = m_PolBeam2(k)
      ENDDO
      END

      SUBROUTINE STest_GetHvectors(HvecFer1,HvecFer2)
*/////////////////////////////////////////////////////////////////////////////////////
*//   Geting final fermion POLARIMETER vectors                                      //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE "STest.h"
      INTEGER k
      REAL*8 HvecFer1(4),HvecFer2(4)
*
      DO k=1,4
         HvecFer1( k) = m_HvecFer1(k)
         HvecFer2( k) = m_HvecFer2(k)
      ENDDO
      END


*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//                           End of CLASS  STest                                   //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////

