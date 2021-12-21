*/////////////////////////////////////////////////////////////////////////////
*//                                                                         //
*// Actual version rewritten by S. Jadach, Nov 1997.                        //
*//                                                                         //
*// Universal random number generator proposed by MARSAGLIA and ZAMAN       //
*// in report FSU-SCRI-87-50                                                //
*//        modified by F. James, 1988 and 1989, to generate a vector        //
*//        of pseudorandom numbers rvec of length lenv, and to put in       //
*//        the COMMON block everything needed to specify currrent state,    //
*//        and to add input and output entry points rmarin, rmarut.         //
*// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   //
*// ++  CALLing sequences for ranmar:                                  ++   //
*// ++      CALL ranmar (rvec, len)   returns a vector rvec of len     ++   //
*// ++                   32-bit random floating point numbers between  ++   //
*// ++                   zero and one.                                 ++   //
*// ++      CALL rmarin(i1,n1,n2)   initializes the generator from one ++   //
*// ++                  32-bit integer i1, and number counts n1,n2     ++   //
*// ++                  (for initializing, set n1=n2=0, but to restart ++   //
*// ++                    a previously generated sequence, use values  ++   //
*// ++                    output by rmarut)                            ++   //
*// ++      CALL rmarut(i1,n1,n2)   outputs the value of the original  ++   //
*// ++                  seed and the two number counts, to be used     ++   //
*// ++                  for restarting by initializing to i1 and       ++   //
*// ++                  skipping n2*100000000+n1 numbers.              ++   //
*// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   //
*//                                                                         //
*//         Initializing routine for ranmar, may be called before           //
*//         generating pseudorandom numbers with ranmar. the input          //
*//         values should be in the ranges:  0<=ijklin<=900 000 000         //
*//                                          0<=ntotin<=999 999 999         //
*//                                          0<=ntot2n<<999 999 999!        //
*// to get the standard values in MARSAGLIA's paper, ijklin=54217137        //
*//                                            ntotin,ntot2n=0              //
*//                                                                         //
*/////////////////////////////////////////////////////////////////////////////


      SUBROUTINE PseuMar_Initialize(ijkl_new, ntot_new,ntot2_new)
*/////////////////////////////////////////////////////////////////////////////
*//                                                                         //
*/////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PseuMar.h'
      INTEGER ijkl_new, ntot_new, ntot2_new
*
      REAL    t,uni,s
      INTEGER m,i24,jj,idum,loop2,now
      INTEGER i,j,ij,k,l,ii,kl

      m_iwarm = 9000009   ! special mark on initialization

      m_ijkl = ijkl_new
      m_ntot = max(ntot_new,0)
      m_ntot2= max(ntot2_new,0)
      
      ij = m_ijkl/30082
      kl = m_ijkl - 30082*ij
      i = MOD(ij/177, 177) + 2
      j = MOD(ij, 177)     + 2
      k = MOD(kl/169, 178) + 1
      l = MOD(kl, 169)
      WRITE(6,'(a,5i10)')
     $     'ranmar initialized: ij,kl,ijkl,ntot,ntot2=',
     $     ij,kl,m_ijkl,m_ntot,m_ntot2
      DO ii= 1, 97
         s = 0.
         t = .5
         DO jj= 1, 24
            m = MOD(MOD(i*j,179)*k, 179)
            i = j
            j = k
            k = m
            l = MOD(53*l+1, 169)
            IF (MOD(l*m,64)  .GE.  32)  s = s+t
            t = 0.5*t
         ENDDO
         m_U(ii) = s
      ENDDO
      m_twom24 = 1.0
      DO i24= 1, 24
         m_twom24 = 0.5*m_twom24
      ENDDO
      m_C  =   362436.*m_twom24
      m_CD =  7654321.*m_twom24
      m_CM = 16777213.*m_twom24
      m_i97 = 97
      m_j97 = 33
*       complete initialization by skipping
*            (ntot2*modcns + ntot) random numbers
      DO loop2= 1, m_ntot2+1
         now = modcns
         IF (loop2  .EQ.  m_ntot2+1)  now=m_ntot
         IF (now  .GT.  0)  THEN
            WRITE(6,'(a,i15)') ' rmarin skipping over ',now
            DO idum = 1, m_ntot
               uni = m_U(m_i97)-m_U(m_j97)
               IF (uni  .LT.  0.)  uni=uni+1.
               m_U(m_i97) = uni
               m_i97 = m_i97-1
               IF (m_i97  .EQ.  0)  m_i97=97
               m_j97 = m_j97-1
               IF (m_j97  .EQ.  0)  m_j97=97
               m_C = m_C - m_CD
               IF (m_C  .LT.  0.)   m_C=m_C+m_CM
            ENDDO
         ENDIF
      ENDDO
      END


      SUBROUTINE PseuMar_MakeVec(rvec,lenv)
*/////////////////////////////////////////////////////////////////////////////
*//                                                                         //
*/////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PseuMar.h'
*
      REAL              rvec(*)
      INTEGER lenv
      INTEGER ijkl_new, ntot_new,ntot2_new
      REAL              zuni,uni
      INTEGER ivec, icont
      DATA icont /0/
*-------------------------------------------------------------------------
      icont = icont+1
      IF (m_iwarm  .NE.  9000009)  THEN
* Default initialization. User has called ranmar without rmarin.
         ijkl_new  = 54217137
         ntot_new  = 0
         ntot2_new = 0
         CALL PseuMar_Initialize(ijkl_new, ntot_new,ntot2_new)
      ENDIF

* Normal entry to generate lenv random numbers
      DO ivec= 1, lenv
         uni = m_U(m_i97)-m_U(m_j97)
         IF (uni  .LT.  0.)  uni=uni+1.
         m_U(m_i97) = uni
         m_i97 = m_i97-1
         IF (m_i97  .EQ.  0)  m_i97=97
         m_j97 = m_j97-1
         IF (m_j97  .EQ.  0)  m_j97=97
         m_C = m_C - m_CD
         IF (m_C  .LT.  0.)   m_C=m_C+m_CM
         uni = uni-m_C
         IF (uni  .LT.  0.) uni=uni+1.
         rvec(ivec) = uni
* Replace exact zeros by uniform distr. *2**-24
         IF (uni  .EQ.  0.)  THEN
            zuni = m_twom24*m_U(2)
*     an exact zero here is very unlikely, but let's be safe.
            IF (zuni  .EQ.  0.) zuni= m_twom24*m_twom24
            rvec(ivec) = zuni
         ENDIF
      ENDDO
c[[[[[[[[[[[[[[[[[[[[[
c      write(16,*) '###PseuMar_MakeVec icont, lenv, rvec=',icont,lenv, (rvec(ivec), ivec=1,lenv)
c]]]]]]]]]]]]]]]]]]]]]
      m_ntot  = m_ntot + lenv
      IF (m_ntot  .GE.  modcns)  THEN
         m_ntot2  =  m_ntot2 + 1
         m_ntot   =  m_ntot - modcns
      ENDIF
      END


      SUBROUTINE PseuMar_Out(ijkl_out, ntot_out, ntot2_out)
*/////////////////////////////////////////////////////////////////////////////
*//                                                                         //
*/////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PseuMar.h'
      INTEGER ijkl_out, ntot_out, ntot2_out
*
      ijkl_out  = m_ijkl
      ntot_out  = m_ntot
      ntot2_out = m_ntot2
      END

