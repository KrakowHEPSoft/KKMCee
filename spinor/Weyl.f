*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//                                                                                 //
*//                         Pseudo-CLASS  Weyl                                      //
*//                                                                                 //
*//       Purpose:  Operations on Weyl spinors                                      //
*//                                                                                 //
*//       Notes:    At the moment only for tests of GPS rules                       //
*//                                                                                 //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////

      SUBROUTINE Weyl_Initialize
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//    Class initialization                                                         //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
*
      INCLUDE "BXformat.h"
      INCLUDE "Weyl.h"

      INTEGER  j1,j2,k
*------------------------------------
      INTEGER init
      SAVE    init
      DATA init/0/
*------------------------------------
      IF(init .EQ. 1) RETURN
      init = 1
*------------------------------------
      m_out    = 16

      WRITE(m_out,bxope)
      WRITE(m_out,bxtxt) '     Initialization of Weyl class    '
      WRITE(m_out,bxclo)
*/////////////////////////////////////////////////////////////
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
*/////////////////////////////////////////////////////////////

      END                       !!!end of Weyl_Initialize!!!


      SUBROUTINE Weyl_Rotor(kaxis,phi1,u1,u2)
*//////////////////////////////////////////////////////////////////////////////////
*//                                                                              //
*//  This is for rotation around any of the 3 axes, kaxis =1,2,3                 //
*//                                                                              //
*//////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE "Weyl.h"
*
      INTEGER    kaxis
      REAL*8     phi1
      COMPLEX*16 u1(4),u2(4)
      INTEGER    i,j,k
*
      COMPLEX*16 Umat(2,2)      !  2x2 blocks in the transformation matrix
      COMPLEX*16 u(4)
      REAL*8     phi
*///////////////////////////////////////////////////
      phi = phi1/2d0
      DO i=1,4
         u(i) = u1(i)
      ENDDO
      DO i=1,2
         DO j=1,2
            Umat(i,j) = m_Pauli(0,     i, j) *DCMPLX( DCOS(phi),       0d0 )
     $                 +m_Pauli(kaxis, i, j) *DCMPLX(       0d0, -DSIN(phi))
         ENDDO
      ENDDO
      u2(1) = Umat(1,1)*u(1)  + Umat(1,2)*u(2)
      u2(2) = Umat(2,1)*u(1)  + Umat(2,2)*u(2)
      u2(3) = Umat(1,1)*u(3)  + Umat(1,2)*u(4)
      u2(4) = Umat(2,1)*u(3)  + Umat(2,2)*u(4)
*
      END                       !!! Weyl_Rotor !!!

      SUBROUTINE Weyl_Boost(kaxis,eta1,u1,u2)
*//////////////////////////////////////////////////////////////////////////////////
*//                                                                              //
*//  This is for rotation around any of the 3 axes, kaxis =1,2,3                 //
*//                                                                              //
*//////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE "Weyl.h"
*
      INTEGER kaxis
      REAL*8  eta1
      COMPLEX*16 u1(4),u2(4)
      INTEGER i,j,k
*
      COMPLEX*16 Umat(2,2),Lmat(2,2)      !  2x2 blocks in the transformation matrix
      COMPLEX*16 u(4)
      REAL*8  eta
*///////////////////////////////////////////////////
      eta = eta1/2d0
      DO i=1,4
         u(i) = u1(i)
      ENDDO
      DO i=1,2
         DO j=1,2
            Umat(i,j) = m_Pauli(0,     i, j) *DCMPLX( DCOSH(eta),     0d0 )
     $                 +m_Pauli(kaxis, i, j) *DCMPLX( +SINH(eta),     0d0 )
            Lmat(i,j) = m_Pauli(0,     i, j) *DCMPLX( DCOSH(eta),     0d0 )
     $                 +m_Pauli(kaxis, i, j) *DCMPLX( -SINH(eta),     0d0 )
         ENDDO
      ENDDO
      u2(1) = Umat(1,1)*u(1)  + Umat(1,2)*u(2)
      u2(2) = Umat(2,1)*u(1)  + Umat(2,2)*u(2)
      u2(3) = Lmat(1,1)*u(3)  + Lmat(1,2)*u(4)
      u2(4) = Lmat(2,1)*u(3)  + Lmat(2,2)*u(4)
*
      END                       !!! Weyl_Boost !!!

      SUBROUTINE Weyl_Print(nunit,word,u)
*///////////////////////////////////////////////////////////////////////////////
*//                                                                           //
*//   prints single Weyl spinor u on unit "nunit" with comment "word"         //
*//                                                                           //
*///////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER      nunit
      CHARACTER*8  word
      COMPLEX*16   u(4)
      INTEGER      i
*----
      WRITE(nunit,'(a8,4(/,a,(2f20.13),a))' ) word,   (' [',u(i),' ] ',  i=1,4)
      END


*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//                           End of CLASS  Weyl                                    //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
