*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//                     Pseudo-CLASS  TabCheck                                      //
*//                                                                                 //
*//  Purpose:                                                                       //
*//      Xcheck implementation of electroweak corrections                           //
*//                                                                                 //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
*
*/////////////////////////////////////////////////////////////////////////////////////
*//       Energy limits in the EW grid, w=sqrt(s) in GeV units.                     //
*/////////////////////////////////////////////////////////////////////////////////////
      DOUBLE PRECISION  m_WminZ, m_WmaxZ                             ! Z pole interval
      DOUBLE PRECISION  m_WminLEP1,          m_WmaxLEP1              ! LEP1 basic interval
      PARAMETER(        m_WminLEP1=0.010d0,  m_WmaxLEP1=120.001d0 )  ! LEP1 basic interval
      DOUBLE PRECISION  m_WmaxLEP2                                   ! LEP2 interval
      PARAMETER(        m_WmaxLEP2  =240.001d0 )                     ! LEP2 interval
      DOUBLE PRECISION  m_WmaxNLC                                    ! NLC interval
      PARAMETER(        m_WmaxNLC  =1040.001d0 )                     ! NLC interval
      INTEGER           m_poin1, m_poin2, m_poin3, m_poin4, m_poinT
*//////////////////////////////////////////////////////////////////////////////
* 340-point grid, only 80pt for NLC to be improved/tested in future
      PARAMETER( m_poin1 = 120 ) ! LEP1 range        (m_WminLEP1,m_WmaxLEP1)
      PARAMETER( m_poin2 =  20 ) ! near Z pole range (amz +- 2*gammz)
      PARAMETER( m_poin3 = 120 ) ! LEP2 energy range (m_WmaxLEP1,m_WmaxLEP2)
      PARAMETER( m_poin4 =  80 ) ! NLC energy range  (m_WmaxLEP2,m_WmaxNLC)
      PARAMETER( m_poinT =  14 ) ! cost(heta) grid
*/////////////////////////////////////////////////////////////////////////////////////
      DOUBLE COMPLEX    m_GSW
      DOUBLE PRECISION  m_CMSene,  m_MZ,     m_GammZ, m_amh,   m_amtop,  m_QCDcor
      DOUBLE PRECISION  m_seps1,   m_seps2,  m_ta,    m_tb
      INTEGER           m_KFini,   m_KFfin,  m_Mode,  m_out
      COMMON  /c_TabCheck/
     $  m_GSW(100),                 ! Electroweak formfactors
     $  m_QCDcor,                   ! QCD correction
     $  m_CMSene,                   ! CMS energy
     $  m_MZ,                       ! Z mass
     $  m_GammZ,                    ! Z width
     $  m_amh,                      ! Higgs mass
     $  m_amtop,                    ! Top mass
     $  m_KFini,                    ! KF code beam
     $  m_KFfin,                    ! KF code final
     $  m_seps1,                    ! helicity init  KORALZ convention
     $  m_seps2,                    ! helicity init  KORALZ convention
     $  m_ta,                       ! helicity final KORALZ convention
     $  m_tb,                       ! helicity final KORALZ convention
     $  m_Mode,                     ! mode for BornV_Differential
     $  m_out                       ! output file number
      SAVE /c_TabCheck/
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////

