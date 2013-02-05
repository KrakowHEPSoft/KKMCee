* from John Holt, not used 
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//                     Pseudo-CLASS  HepEvt_Posn                            //
*//                                                                          //
*//  Purpose:  keep and serve info on positions of fermions and KKphotons    //
*//           in HEPEVT records                                              //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////

* ----------------------------------------------------------------------
      INTEGER m_PhotStart ! start posn of photons (1st phot at m_PhotStart+1)
      INTEGER m_PhotEnd   ! end posn photons
      INTEGER m_PosnF     ! posn of final state fermion
      INTEGER m_PosnFbar  ! posn of final state anti-fermion

      COMMON/HepEvt_posn/m_PhotStart,m_PhotEnd,m_PosnF,m_PosnFbar

      SAVE/HepEvt_posn/

*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//                      End of CLASS  HepEvt_Posn                           //
*//////////////////////////////////////////////////////////////////////////////
