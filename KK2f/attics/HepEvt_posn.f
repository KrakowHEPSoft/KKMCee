* from John Holt, not used

      SUBROUTINE HepEvt_PosnSetPhotStart(Posn)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*// Purpose:  Set start position of photons in HepEvt common                 //
*//           First photon is at StartPosn+1                                 //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE

      INCLUDE 'HepEvt_posn.h'

      INTEGER Posn

      m_PhotStart=Posn

      RETURN
      END

      SUBROUTINE HepEvt_PosnGetPhotStart(Posn)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*// Purpose:  Get start position of photons in HepEvt common                 //
*//           First photon is at StartPosn+1                                 //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE

      INCLUDE 'HepEvt_posn.h'

      INTEGER Posn

      Posn=m_PhotStart

      RETURN
      END

      SUBROUTINE HepEvt_PosnSetPhotEnd(Posn)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*// Purpose:  Get end position of photons in HepEvt common                   //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE

      INCLUDE 'HepEvt_posn.h'

      INTEGER Posn

      m_PhotEnd=Posn

      RETURN
      END



      SUBROUTINE HepEvt_PosnGetPhotEnd(Posn)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*// Purpose:  Get end position of photons in HepEvt common                   //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE

      INCLUDE 'HepEvt_posn.h'

      INTEGER Posn

      Posn=m_PhotEnd

      RETURN
      END



      SUBROUTINE HepEvt_PosnGetPhotNumb(Numb)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*// Purpose:  Get number of photons in HepEvt common                         //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE

      INCLUDE 'HepEvt_posn.h'

      INTEGER Numb

      Numb=m_PhotEnd-m_PhotStart

      RETURN
      END


      SUBROUTINE HepEvt_PosnSetF(Posn)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*// Purpose:  Set position of final state fermion in HepEvt common           //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE

      INCLUDE 'HepEvt_posn.h'

      INTEGER Posn

      m_PosnF=Posn

      RETURN
      END

      SUBROUTINE HepEvt_PosnGetF(Posn)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*// Purpose:  Set position of final state fermion in HepEvt common           //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE

      INCLUDE 'HepEvt_posn.h'

      INTEGER Posn

      Posn=m_PosnF

      RETURN
      END



      SUBROUTINE HepEvt_PosnSetFbar(Posn)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*// Purpose:  Set position of final state anti-fermion in HepEvt common      //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE

      INCLUDE 'HepEvt_posn.h'

      INTEGER Posn

      m_PosnFbar=Posn

      RETURN
      END

      SUBROUTINE HepEvt_PosnGetFbar(Posn)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*// Purpose:  Set position of final state anti-fermion in HepEvt common      //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE

      INCLUDE 'HepEvt_posn.h'

      INTEGER Posn

      Posn=m_PosnFbar

      RETURN
      END


