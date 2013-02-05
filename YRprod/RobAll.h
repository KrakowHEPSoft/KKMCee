*///////////////////////////////////////////////////////////////////////////////////////
*//                                                                                   //
*//          Pseudoclass RobAll                                                       //
*//                                                                                   //
*///////////////////////////////////////////////////////////////////////////////////////
*  Sigma+AFB     All flavours           50 selections           F+B and F-B           10 models
      INTEGER    m_iKFf1,  m_iKFf2,     m_iSel1,  m_iSel2,      m_iThe1,  m_iThe2,    m_iMod1,  m_iMod2     !
      PARAMETER( m_iKFf1=1,m_iKFf2=16,  m_iSel1=1,m_iSel2=50,   m_iThe1=1,m_iThe2=2,  m_iMod1=1,m_iMod2=10) !
*  dSig/dCos     Muon only              25 selections           40 bins in theta      10 models
      INTEGER    m_jKFf1,   m_jKFf2,    m_jSel1,  m_jSel2,      m_jThe1,  m_jThe2,    m_jMod1,  m_jMod2     !
      PARAMETER( m_jKFf1=13,m_jKFf2=13, m_jSel1=1,m_jSel2=25,   m_jThe1=1,m_jThe2=40, m_jMod1=1,m_jMod2=10) !
*  Dimensions of big storage histograms
      INTEGER    m_nBinSAB, m_nBinANG
      PARAMETER( m_nBinSAB= (m_iKFf2-m_iKFf1+1)*(m_iSel2-m_iSel1+1)*(m_iThe2-m_iThe1+1)*(m_iMod2-m_iMod1+1)) !
      PARAMETER( m_nBinANG= (m_jKFf2-m_jKFf1+1)*(m_jSel2-m_jSel1+1)*(m_jThe2-m_jThe1+1)*(m_jMod2-m_jMod1+1)) !
      INTEGER    m_kSAB,        m_kANG,           m_kSABren,        m_kANGren           !
      PARAMETER( m_kSAB=100000, m_kANG=200000,    m_kSABren=300000, m_kANGren=400000  ) !
      INTEGER    m_kSABsem
      PARAMETER( m_kSABsem=500000)

      INTEGER      m_nVmax
      PARAMETER  ( m_nVmax = 21)
      DOUBLE PRECISION  m_VmaxList(m_nVmax)
*  This is original LEP2 wshop assignement
      DATA m_VmaxList /    0.01d0, 0.05d0,                 ! 1,2
     $     0.10d0, 0.15d0, 0.20d0, 0.25d0, 0.30d0, 0.35d0, ! 3-8
     $     0.40d0, 0.45d0, 0.50d0, 0.55d0, 0.60d0, 0.65d0, ! 9-14
     $     0.70d0, 0.75d0, 0.80d0, 0.85d0, 0.90d0, 0.95d0, 0.99d0/ ! 15-21
*  This is for nunu exercises
c$$$      DATA m_VmaxList /    0.01d0, 0.05d0,         ! 1,2
c$$$     $     0.10d0, 0.15d0, 0.20d0, 0.30d0,         ! 3-6
c$$$     $     0.40d0, 0.50d0, 0.55d0, 0.60d0, 0.65d0, ! 7-11
c$$$     $     0.70d0, 0.75d0, 0.80d0, 0.85d0, 0.90d0, ! 12-16
c$$$     $     0.95d0, 0.96d0, 0.97d0, 0.98d0, 0.99d0/ ! 17-21
      SAVE m_VmaxList
* Ofset poniters for YR selections
      INTEGER      m_nYR
***      PARAMETER  ( m_nYR  = m_nVmax+3) ! OLDER one, production before 6-th June
      PARAMETER  ( m_nYR  = m_nVmax) ! NEW one after 6-th June
* Angular range
      INTEGER      m_nBinThe
      PARAMETER  ( m_nBinThe  =  m_jThe2-m_jThe1+1 )

      DOUBLE PRECISION    m_WtSel, m_WtThe, m_WtMod, m_sum1(100), m_sum2(100) !
      COMMON /c_RobAll/   
     $ m_WtSel(200),                ! Mask Selections
     $ m_WtThe(m_nBinThe),          ! +1,-1 for F,B or mask funcion for single bin
     $ m_WtMod(m_iMod1:m_iMod2),    ! Model weights
     $ m_sum1, m_sum2               ! for tests

*///////////////////////////////////////////////////////////////////////////////////////
*//                                                                                   //
*//          End of Pseudoclass Robol1                                                //
*//                                                                                   //
*///////////////////////////////////////////////////////////////////////////////////////
