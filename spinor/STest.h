*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//                                                                          //
*//                     Pseudo-CLASS  STest                                  //
*//                                                                          //
*//   Purpose:  Tests on spin amplitudes calculated using spinor methods     //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////

      REAL*8     m_pi
      PARAMETER( m_pi=3.1415926535897932d0)
*
      INTEGER    m_out
      COMPLEX*16 m_Pauli,    m_Pauli4
      REAL*8     m_PolBeam1,         m_PolBeam2
      REAL*8     m_HvecFer1,         m_HvecFer2
      COMPLEX*16 m_SDMat1,       m_SDMat2,       m_SDMat3,       m_SDMat4
*
      COMMON /c_STest/
     $    m_Pauli( 0:3, 1:2, 1:2), ! Pauli matrices
     $    m_Pauli4(1:4, 1:2, 1:2), ! Pauli matrices, other vector index numbering
     $    m_PolBeam1(4),           ! POLARIZATION vector 1-st beam
     $    m_PolBeam2(4),           ! POLARIZATION vector 2-nd beam
     $    m_HvecFer1(4),           ! POLARIMETER  vector 1-st fin.ferm.
     $    m_HvecFer2(4),           ! POLARIMETER  vector 2-nd fin.ferm.
     $    m_SDMat1(2,2),           ! Spin Density Matrix  1-st beam 
     $    m_SDMat2(2,2),           ! Spin Density Matrix  2-nd beam 
     $    m_SDMat3(2,2),           ! Polarimeter Density Matrix  1-st fin.ferm
     $    m_SDMat4(2,2),           ! Polarimeter Density Matrix  2-nd fin.ferm
     $    m_out                    ! output unit number

      SAVE /c_STest/
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//                                                                          //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////

