* from John Holt, not used
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//                     Pseudo-CLASS  HepEvt                                 //
*//                                                                          //
*//  Purpose:  keep and serve event in HEPEVT format                         //
*//                                                                          //
*//  Output of KK2f is encoded in double precission /d_hepevt/               //
*//  which is double precision version of /hepevt/                           //
*//  It was necessary to rename /hepevt/                                     //
*//  because older Jetset uses single precision version of /hepevt/          //
*//                                                                          //
*//  We introduce luhepc->luhepcd which is                                   //
*//  Double Precision version of luhepc. It translates                       //
*//  from /hepevtd/ to (from) Jetset-common-blocks                           //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////


      SUBROUTINE HepEvt_Fill
*//////////////////////////////////////////////////////////////////////////////////////
*//                                                                                  //
*//  Filling HEP COMMON block using HepEvt_Fil1 routine                              //
*//                                                                                  //
*//  INPUT:                                                                          //
*//     KFfin          final particle code       (KarLud)                            //
*//     pf1,pf2        beam (electron) momenta   (KarLud)                            //
*//     nphox          photon multiplicity       (KarLud)                            //
*//     xphot          photon momenta            (KarLud)                            //
*//     qf1,qf2        final momenta             (Karfin)                            //
*//     nphoy          photon multiplicity       (KarFin)                            //
*//     yphot          photon momenta            (KarFin)                            //
*//                                                                                  //
*//  !!!!!! To be disscussed with the users. !!!!!!                                  //
*//  In present version we add beamstrahlung photons to the record                   //
*//  along with the ISR photons. It can be missleading!                              //
*//  May be they should be rather listed                                             //
*//  as the initial state particles, together with beams?                            //
*//                                                                                  //
*//////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION  pf1(4),pf2(4),qf1(4),qf2(4),xphot(100,4),yphot(100,4)
      DOUBLE PRECISION  aph(5), bph(5),pph(5)
      DOUBLE PRECISION  amel,amfin,BornV_GetMass
      DOUBLE PRECISION  Psum(4), Etot
      DOUBLE PRECISION  WtMain,WtCrud
      INTEGER           KFfin,ip,kstat,j,i,nphoy,nphox,kfbeam

      INTEGER           nbeams,nphot,nfinal,nboson
      INTEGER           ida1,ida2
      INTEGER           imo1,imo2
CD
      double precision  pbos(4),xmtmp
*------------------------------------------------------------------------------
* Actual KFcode of final fermion
      CALL KarLud_GetKFfin(KFfin)
* Beams and ISR photons
      CALL KarLud_GetPhotons(nphox,xphot)
      CALL KarLud_GetBeams(pf1,pf2)
* Beamstrahlung photons
      CALL KarLud_GetBeasts( aph, bph)
* Final fermions and FSR photons
      CALL KarFin_GetPhotons(nphoy,yphot)
      CALL KarFin_GetFermions(qf1,qf2)
      DO j=1,4
         Psum(j)=qf1(j)+qf2(j)
      ENDDO
*
* Count total number of photons
      nphot=nphox+nphoy
      IF( aph(4).NE.0d0) nphot=nphot+1
      IF( bph(4).NE.0d0) nphot=nphot+1

      KFbeam = 11    ! KF=11 denotes electron
      amel   =BornV_GetMass(KFbeam)

* Put beams into record
* NB HepEvt_Fil1 will update daughters if decay products added
      ip=1
      CALL HepEvt_Fil1(ip,3, 11, 0,0,0,0, pf1,amel,.FALSE.)
      ip=ip+1
      CALL HepEvt_Fil1(ip,3,-11, 0,0,0,0, pf2,amel,.FALSE.)

* Final state fermions
      IF( KFfin .EQ. 0  ) RETURN
      amfin = BornV_GetMass(KFfin)
      KStat = 1

* Put final state fermions into record
* originally only particle 1 was identified as mother
c      ip=ip+1
c      CALL HepEvt_Fil1(ip,KStat, KFfin, 1,2,0,0, qf1,amfin,.FALSE.)
c      CALL HepEvt_PosnSetF(ip)
c      ip=ip+1
c      CALL HepEvt_Fil1(ip,KStat,-KFfin, 1,2,0,0, qf2,amfin,.FALSE.)
c      CALL HepEvt_PosnSetFbar(ip)

CD: ffbar mother added in event record 
* 1,2             - beams
* 3               - Z/gamma*
* 4...4+nphot     - photons 
* 4+nphot,5+nphot - final fermions   

C Add boson to list
      do j=1,4
         pbos(j)=pf1(j)+pf2(j)
         do i=1,nphox
            pbos(j)=pbos(j)-xphot(i,j)
         enddo
cc       do i=1,nphoy
cc        pbos(j)=pbos(j)-yphot(i,j)
cc       enddo
cc       if( aph(4).NE.0d0) then
cc        pbos(j)=pbos(j)-aph(j)
cc       endif
cc       if( bph(4).NE.0d0) then
cc        pbos(j)=pbos(j)-bph(j)
cc       endif
      enddo
      xmtmp=pbos(4)**2-(pbos(1)**2+pbos(2)**2+pbos(3)**2)
      if(xmtmp.lt.0.)xmtmp=0.d0
      xmtmp=dsqrt(xmtmp)
      ip=ip+1
      CALL HepEvt_Fil1(ip,2,23,1,2,0,0,pbos,xmtmp,.FALSE.)


* Radiative photons (6 ... 5+nphot) (pdg-code for gamma is 22)
* Note that both ISR and FSR photons have beams as parents.
* This is because (a) JETSET does not tolerate FSR photons attached to quarks
*                 (b) photons for CEEX model cannot be split into ISR and FSR
      CALL HepEvt_PosnSetPhotStart(ip)
      IF(nphox .NE. 0) THEN
         DO i=1,nphox
            DO j=1,4
               pph(j) = xphot(i,j)
               Psum(j)= Psum(j)+xphot(i,j)
            END DO
            ip=ip+1
            CALL HepEvt_Fil1(ip,1,22, 1,2,0,0, pph,0d0,.FALSE.) ! ISR
         END DO
      END IF
      IF(nphoy .NE. 0) THEN
         DO i=1,nphoy
            DO j=1,4
               pph(j) = yphot(i,j)
               Psum(j)= Psum(j)+yphot(i,j)
            END DO
            ip=ip+1
c- PJH originally final state photons where filled with mothers 2,1
c- change to proper HepEvt order 1,2
cc            CALL HepEvt_Fil1(ip,1,22, 2,1,0,0, pph,0d0,.FALSE.) ! FSR
            CALL HepEvt_Fil1(ip,1,22, 1,2,0,0, pph,0d0,.FALSE.) ! FSR
         END DO
      END IF
* add beamstrahlung photons to the record along with ISR photons, can be missleading!!
      IF( aph(4).NE.0d0) THEN
         ip=ip+1
         CALL HepEvt_Fil1(ip,1,22, 1,2,0,0, aph,0d0,.FALSE.)
         DO j=1,4
            Psum(j)= Psum(j)+aph(j)
         ENDDO
      ENDIF
      IF( bph(4).NE.0d0) THEN
         ip=ip+1
         CALL HepEvt_Fil1(ip,1,22, 1,2,0,0, bph,0d0,.FALSE.)
         DO j=1,4
            Psum(j)= Psum(j)+bph(j)
         ENDDO
      ENDIF
      CALL HepEvt_PosnSetPhotEnd(ip)



* add final state particles to list (assuming come from boson)
      ip=ip+1
      CALL HepEvt_Fil1(ip,KStat, KFfin, 3,3,0,0, qf1,amfin,.FALSE.)
      CALL HepEvt_PosnSetF(ip)
      ip=ip+1
      CALL HepEvt_Fil1(ip,KStat,-KFfin, 3,3,0,0, qf2,amfin,.FALSE.)
      CALL HepEvt_PosnSetFbar(ip)



      Etot= SQRT(ABS(Psum(4)**2 -Psum(3)**2 -Psum(2)**2 -Psum(1)**2))
* Check on total 4-momentum conservation
      IF( ABS(Etot/(pf1(4)+pf2(4)+aph(4)+bph(4))-1d0) .GT.1d-4) THEN
         WRITE(*,*) '++++ HepEvt_Fill: something wrong with Etot=',Etot
      ENDIF
      IF( ABS(Psum(4)/(pf1(4)+pf2(4)+aph(4)+bph(4))-1d0) .GT.1d-4) THEN
         WRITE(*,*) '++++ HepEvt_Fill: something wrong with Psum(4)=',Psum(4)
      ENDIF

* Finaly fill also LUND common block
      CALL HepEvt_LuHepc(2)
      END


      SUBROUTINE HepEvt_Hadronize(HadMin)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//          Aranging jets and hadronization                                 //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'HepEvt.h'
*
      INTEGER  ijoin(2)
      DOUBLE PRECISION    HadMin, sqrs1

      INTEGER  ih1,ih2
* ----------------------------------------------------------------------

* Quarks only, KeyHad=1 required
cc      ih1=4     ! fermion line
cc      ih2=5     ! antifermion line
      Call HepEvt_PosnGetF(ih1)
      Call HepEvt_PosnGetFbar(ih2)

      IF ( ABS(idhep(ih1)) .LT. 10 ) THEN 
** Explicit string arangement:
         ijoin(1) = ih1
         ijoin(2) = ih2
** q-qbar effective mass
         sqrs1=(phep(4,ih1)+phep(4,ih2))**2
     $        -(phep(3,ih1)+phep(3,ih2))**2
     $        -(phep(2,ih1)+phep(2,ih2))**2
     $        -(phep(1,ih1)+phep(1,ih2))**2
         sqrs1=sqrt(abs(sqrs1))
* Showering < HadMas cut-off value (this also deals with WT=0 events)
* see also bornv
         IF( sqrs1 .GT. HadMin**2) THEN
****      IF ( ABS(idhep(ih1)) .LT. 10 ) THEN 
            CALL PYjoin(2,ijoin)
****      ENDIF	
            CALL PYshow(ih1,ih2,sqrs1)
            CALL PYexec
         ENDIF
      ENDIF
      
      END


      SUBROUTINE HepEvt_Fil1( n,ist,id,jmo1,jmo2,jda1,jda2,p4,pinv,phflag)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*// This subroutine fills one entry into the hepevt common                   //
*// and updates the information for affected mother entries                  //
*// WRITTEN by Martin W. Gruenewald (91/01/28)                               //
*// Re-Edited by S. Jadach, 6 july 97                                        //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'HepEvt.h'
* input
      INTEGER           n,ist,id,jmo1,jmo2,jda1,jda2
      DOUBLE PRECISION  p4(4),pinv
      LOGICAL           phflag
* locals
      INTEGER           ip,i,ihep
      SAVE
* ----------------------------------------------------------------------
*
* check address mode
      IF (n .EQ. 0) THEN
* append mode
        ihep=nhep+1
      ELSE IF (n .GT. 0) THEN
* absolute position
        ihep=n
      ELSE
* relative position
        ihep=nhep+n
      END IF
*
* check on ihep
      IF ((ihep .LE. 0) .OR. (ihep .GT. nmxhep)) RETURN
*
* add entry
      nhep=ihep
      isthep(ihep)   =ist
      idhep(ihep)    =id
      jmohep(1,ihep) =jmo1
      IF(jmo1 .LT. 0)  jmohep(1,ihep)=jmohep(1,ihep)+ihep
      jmohep(2,ihep) =jmo2
      IF(jmo2 .LT. 0)  jmohep(2,ihep)=jmohep(2,ihep)+ihep
      jdahep(1,ihep) =jda1
      jdahep(2,ihep) =jda2
*
      DO i=1,4
         phep(i,ihep)=p4(i)
* KORAL-B and KORAL-Z do not provide vertex and/or lifetime informations
         vhep(i,ihep)=0d0
      END DO
      phep(5,ihep)=pinv
* flag for photos...
      qedrad(ihep)=phflag
* update process:
      DO ip=jmohep(1,ihep),jmohep(2,ihep)
         IF(ip .GT. 0)THEN
* IF there is a daughter at ihep, mother entry at ip has decayed
            IF(isthep(ip) .EQ. 1)isthep(ip)=2
* and daughter pointers of mother entry must be updated
            IF(jdahep(1,ip) .EQ. 0)THEN
               jdahep(1,ip)=ihep
               jdahep(2,ip)=ihep
            ELSE
               jdahep(2,ip)=max(ihep,jdahep(2,ip))
            END IF
         END IF
      END DO
*
      END

      SUBROUTINE HepEvt_LuHepc(mconv)
*     *******************************
* This is double precision version of luhepc from jetset
* It translates between DOUBLE PRECISION  /c_HepEvt/ and 
*****************************************************REAL
* Lund commons.
* Corrections marked with !!! STJ
*     ******************************************************************
*...purpose: to convert jetset event record contents to or from 
*...the standard event record common block. 
      IMPLICIT DOUBLE PRECISION (a-h,o-z)  !!! STJ
* ----------------------------------------------------------------------
      INCLUDE 'HepEvt.h'
* ----------------------------------------------------------------------
*new: LU -> PY, real -> DP
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE
      INTEGER PYCOMP	
* ---------------------------------------------------------------------- 
*...conversion from jetset to standard, the easy part. 
      IF(mconv .EQ. 1) THEN 
        nevhep=0 
        IF(n .GT. nmxhep) CALL PYerrm(8, 
     &  '(luhepc:) no more space in /hepevtd/') 
        nhep=min(n,nmxhep) 
        DO 140 i=1,nhep 
        isthep(i)=0 
        IF(k(i,1) .GE. 1  .AND. k(i,1) .LE. 10) isthep(i)=1 
        IF(k(i,1) .GE. 11 .AND. k(i,1) .LE. 20) isthep(i)=2 
        IF(k(i,1) .GE. 21 .AND. k(i,1) .LE. 30) isthep(i)=3 
        IF(k(i,1) .GE. 31 .AND. k(i,1) .LE. 100) isthep(i)=k(i,1) 
        idhep(i)=k(i,2) 
        jmohep(1,i)=k(i,3) 
        jmohep(2,i)=0 
        IF(k(i,1) .NE. 3 .AND. k(i,1) .NE. 13 .AND. k(i,1) .NE. 14) THEN
          jdahep(1,i)=k(i,4) 
          jdahep(2,i)=k(i,5) 
        ELSE 
          jdahep(1,i)=0 
          jdahep(2,i)=0 
        ENDIF 
        DO 100 j=1,5 
        phep(j,i)=p(i,j) 
  100   CONTINUE 
        DO 110 j=1,4 
        vhep(j,i)=v(i,j) 
  110   CONTINUE 
 
*...check IF new event (from pileup). 
        IF(i .EQ. 1) THEN 
          inew=1 
        ELSE 
          IF(k(i,1) .EQ. 21 .AND. k(i-1,1) .NE. 21) inew=i 
        ENDIF 
 
*...fill in missing mother inFORMATion. 
        IF(i .GE. inew+2 .AND. k(i,1) .EQ. 21 .AND. k(i,3) .EQ. 0) THEN 
           imo1=i-2 
          IF(i .GE. inew+3 .AND. k(i-1,1) .EQ. 21 .AND. k(i-1,3) .EQ. 0)
     &         imo1=imo1-1 
          jmohep(1,i)=imo1 
          jmohep(2,i)=imo1+1 
        ELSEIF(k(i,2) .GE. 91 .AND. k(i,2) .LE. 93) THEN 
          i1=k(i,3)-1 
  120     i1=i1+1 
          IF(i1 .GE. i) CALL PYerrm(8, 
     &    '(luhepc:) translation of inconsistent event history') 
          IF(i1 .LT. i .AND. k(i1,1) .NE. 1 .AND. k(i1,1) .NE. 11) 
     $         GOTO 120 
          kc=PYcomp(k(i1,2)) 
          IF(i1 .LT. i .AND. kc .EQ. 0) GOTO 120 
          IF(i1 .LT. i .AND. kchg(kc,2) .EQ. 0) GOTO 120 
          jmohep(2,i)=i1 
        ELSEIF(k(i,2) .EQ. 94) THEN 
          njet=2 
          IF(nhep .GE. i+3 .AND. k(i+3,3) .LE. i) njet=3 
          IF(nhep .GE. i+4 .AND. k(i+4,3) .LE. i) njet=4 
          jmohep(2,i)=MOD(k(i+njet,4)/mstu(5),mstu(5)) 
          IF(jmohep(2,i) .EQ. jmohep(1,i)) jmohep(2,i)= 
     &    MOD(k(i+1,4)/mstu(5),mstu(5)) 
        ENDIF 
 
*...fill in missing daughter inFORMATion. 
        IF(k(i,2) .EQ. 94 .AND. mstu(16) .NE. 2) THEN 
          DO 130 i1=jdahep(1,i),jdahep(2,i) 
          i2=MOD(k(i1,4)/mstu(5),mstu(5)) 
          jdahep(1,i2)=i 
  130     CONTINUE 
        ENDIF 
        IF(k(i,2) .GE. 91 .AND. k(i,2) .LE. 94) GOTO 140 
        i1=jmohep(1,i) 
        IF(i1 .LE. 0 .OR. i1 .GT. nhep) GOTO 140 
        IF(k(i1,1) .NE. 13 .AND. k(i1,1) .NE. 14) GOTO 140 
        IF(jdahep(1,i1) .EQ. 0) THEN 
          jdahep(1,i1)=i 
        ELSE 
          jdahep(2,i1)=i 
        ENDIF 
  140   CONTINUE 
        DO 150 i=1,nhep 
        IF(k(i,1) .NE. 13 .AND. k(i,1) .NE. 14) GOTO 150 
        IF(jdahep(2,i) .EQ. 0) jdahep(2,i)=jdahep(1,i) 
  150   CONTINUE 
 
*...conversion from standard to jetset, the easy part. 
      ELSE 
        IF(nhep .GT. mstu(4)) CALL PYerrm(8, 
     &  '(luhepc:) no more space in /PYjets/') 
        n=min(nhep,mstu(4)) 
        nkq=0 
        kqsum=0 
        DO 180 i=1,n 
        k(i,1)=0 
        IF(isthep(i) .EQ. 1) k(i,1)=1 
        IF(isthep(i) .EQ. 2) k(i,1)=11 
        IF(isthep(i) .EQ. 3) k(i,1)=21 
        k(i,2)=idhep(i) 
        k(i,3)=jmohep(1,i) 
        k(i,4)=jdahep(1,i) 
        k(i,5)=jdahep(2,i) 
        DO 160 j=1,5 
        p(i,j)=phep(j,i) 
  160   CONTINUE 
        DO 170 j=1,4 
        v(i,j)=vhep(j,i) 
  170   CONTINUE 
        v(i,5)=0d0 
        IF(isthep(i) .EQ. 2 .AND. phep(4,i) .GT. phep(5,i)) THEN 
          i1=jdahep(1,i) 
          IF(i1 .GT. 0 .AND. i1 .LE. nhep) 
     $         v(i,5)=(vhep(4,i1)-vhep(4,i))*phep(5,i)/phep(4,i) 
        ENDIF 

*...fill in missing information on colour connection in jet systems. 
        IF(isthep(i) .EQ. 1) THEN 
          kc=PYcomp(k(i,2)) 
          kq=0 
          IF(kc .NE. 0) kq=kchg(kc,2)*isign(1,k(i,2)) 
          IF(kq .NE. 0) nkq=nkq+1 
          IF(kq .NE. 2) kqsum=kqsum+kq 
          IF(kq .NE. 0 .AND. kqsum .NE. 0) THEN 
            k(i,1)=2 
          ELSEIF(kq .EQ. 2 .AND. i .LT. n) THEN 
            IF(k(i+1,2) .EQ. 21) k(i,1)=2 
          ENDIF 
        ENDIF 
  180   CONTINUE 
        IF(nkq .EQ. 1 .OR. kqsum .NE. 0) CALL PYerrm(8, 
     &  '(luhepc:) input parton configuration not colour singlet') 
      ENDIF 
      END 

      SUBROUTINE HepEvt_GetKFfin(KFfin)
*//////////////////////////////////////////////////////////////////////////////
*// Purpose:  Get KFcode of final fermion out of /hepevt/                    //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'HepEvt.h'
      INTEGER KFfin

      INTEGER ih1
* ----------------------------------------------------------------------
cc      ih1=4	
      CALL HepEvt_PosnGetF(ih1)

      KFfin = idhep(ih1)
      END ! GetKFfin

      SUBROUTINE HepEvt_GetBeams(pf1,pf2)
*//////////////////////////////////////////////////////////////////////////////
*// Purpose:  Get Beam momenta out of /hepevt/                               //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'HepEvt.h'
      DOUBLE PRECISION  pf1(4),pf2(4)
      INTEGER k
* ----------------------------------------------------------------------
      DO k=1,4
         pf1(k) =phep(k,1)
         pf2(k) =phep(k,2)
      ENDDO
      END ! GetBeams
      SUBROUTINE HepEvt_GetParticle(Id,Istart,Iadress,pf1)
*//////////////////////////////////////////////////////////////////////////////
*// Purpose:  Get Id-partcle momenta out of /hepevt/                         //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'HepEvt.h'
      DOUBLE PRECISION  pf1(4)
      INTEGER k,l,Id,Istart,Iadress
* ----------------------------------------------------------------------
      Iadress=-1
      DO l=Istart,nhep
        IF(Idhep(l).eq.Id) then
          Iadress=l
          DO k=1,4
             pf1(k) =phep(k,l)
          ENDDO
          RETURN
        ENDIF
      ENDDO
      END ! GetParticle

      SUBROUTINE HepEvt_GetFfins(pf1,pf2)
*//////////////////////////////////////////////////////////////////////////////
*// Purpose:  Get final fermion momenta out of /hepevt/                      //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'HepEvt.h'

      DOUBLE PRECISION  pf1(4),pf2(4)
      INTEGER k

      integer ih1,ih2
* ----------------------------------------------------------------------
cc      ih1=4
cc      ih2=5
      CALL HepEvt_PosnGetF(ih1)
      CALL HepEvt_PosnGetFbar(ih1)

      DO k=1,4
         pf1(k) =phep(k,ih1)
         pf2(k) =phep(k,ih2)
      ENDDO
      END ! GetFfins

      SUBROUTINE HepEvt_GetNPhot(nphot)
*//////////////////////////////////////////////////////////////////////////////
*// Purpose:  Get number of bremsstrahlung photons  out of /hepevt/          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'HepEvt.h'

      INTEGER nphot,kf,j,ih
      integer ifot,ih1
*
cc      ih1=4
cc      ifot=5	! photons start after eeZff	

      CALL HepEvt_PosnGetF(ih1)
      CALL HepEvt_PosnGetPhotStart(ifot)

      nphot=0
      DO j=1,100
         ih = ifot+j
         kf = idhep(ih)
* STOP if /hepevt/ ended or non-photon found
         IF(ih .GT. nhep) GOTO 110
         IF(kf .NE. 22)   GOTO 110
* Bremsstrahlung photons have by convention 1-st parent being 1-st beam (ISR)
* at ih=1 or first final fermion at ih1=4 (FSR)
         IF(jmohep(1,ih).EQ.1 .OR. jmohep(1,ih).EQ.ih1) nphot=nphot+1
      ENDDO
      RETURN
 110  CONTINUE
      END


      SUBROUTINE HepEvt_GetPhotAll(nphot,phot)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*// Purpose:  Get ALL photons out of /hepevt/                                //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'HepEvt.h'
      DOUBLE PRECISION  phot(100,4)
      INTEGER nphot,j,k,ih,kf

      integer ifot,ih1
* ----------------------------------------------------------------------
cc      ifot=5	! photons start after eeZff	
cc      ih1=4
      CALL HepEvt_PosnGetF(ih1)
      CALL HepEvt_PosnGetPhotStart(ifot)

      nphot=0
      DO j=1,100
         ih = ifot+j
         kf = idhep(ih)
* STOP if /hepevt/ ended or non-photon found
         IF(ih .GT. nhep) GOTO 110
         IF(kf .NE. 22)   GOTO 110
* Bremsstrahlung photons have by convention 1-st parent being 1-st beam (ISR)
* at ih=1 or first final fermion at ih1=4 (FSR)
         IF(jmohep(1,ih).EQ.1 .OR. jmohep(1,ih).EQ.ih1) THEN
            nphot=nphot+1
            DO  k=1,4
               phot(nphot,k) =phep(k,ih)
            ENDDO
         ENDIF
      ENDDO
      GOTO 900
 110  CONTINUE
      RETURN
 900  WRITE(*,*) '++++ HepEvt_GetPhotAll: To many photons!!!'
      END ! HepEvt_GetPhotAll

      SUBROUTINE HepEvt_GetPhotIni(nphot,phot)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*// Purpose:  Get ISR photons out of /hepevt/                                //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'HepEvt.h'
      DOUBLE PRECISION  phot(100,4)
      INTEGER nphot,j,k,ih,kf
      integer ifot
*
cc      ifot=5	! photons start after eeZff	
      CALL HepEvt_PosnGetPhotStart(ifot)
* ----------------------------------------------------------------------
      nphot=0
      DO j=1,100
         ih = ifot+j               
         kf = idhep(ih)
* STOP if /hepevt/ ended or non-photon found
         IF(ih .GT. nhep) GOTO 110
         IF(kf .NE. 22)   GOTO 110
* ISR photons have by convention 1-st parent being 1-st beam (convention)
         IF(jmohep(1,ih) .EQ. 1) THEN
            nphot=nphot+1
            DO  k=1,4
               phot(nphot,k) =phep(k,ih)
            ENDDO
         ENDIF
      ENDDO
      GOTO 900
 110  CONTINUE
      RETURN
 900  WRITE(*,*) '++++ HepEvt_GetPhotIni: To many photons!!!'
      END ! HepEvt_GetPhotIni

      SUBROUTINE HepEvt_GetPhotFin(nphot,phot)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*// Purpose:  Get FSR photons out of /hepevt/                                //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'HepEvt.h'
      DOUBLE PRECISION  phot(100,4)
      INTEGER nphot,j,k,ih,kf
      integer ifot
*
cc      ifot=5	! photons start after eeZff	
      CALL HepEvt_PosnGetPhotStart(ifot)
* ----------------------------------------------------------------------
      nphot=0
      DO j=1,100
         ih = ifot+j               
         kf = idhep(ih)
* STOP if /hepevt/ ended or non-photon found
         IF(ih .GT. nhep) GOTO 110
         IF(kf .NE. 22)   GOTO 110
* FSR photons have 1-st parent being 2-nd beam (convention)
         IF(jmohep(1,ih) .EQ. 2) THEN
            nphot=nphot+1
            DO  k=1,4
               phot(nphot,k) =phep(k,ih)
            ENDDO
         ENDIF
      ENDDO
      GOTO 900
 110  CONTINUE
      RETURN
 900  WRITE(*,*) '++++ HepEvt_GetPhotFin: To many photons!!!'
      END ! HepEvt_GetPhotFin

      SUBROUTINE HepEvt_GetPhotBst(nphot,phot)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*// Purpose:  Get Beamstrahlung photons out of /hepevt/                      //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'HepEvt.h'
      DOUBLE PRECISION  phot(100,4)
      INTEGER nphot,j,k,ih,kf
      integer ifot
*
cc      ifot=5	! photons start after eeZff	
      CALL HepEvt_PosnGetPhotStart(ifot)
* ----------------------------------------------------------------------
      nphot=0
      DO j=1,100
         ih = ifot+j               
         kf = idhep(ih)
* STOP if /hepevt/ ended or non-photon found
         IF(ih .GT. nhep) GOTO 110
         IF(kf .NE. 22)   GOTO 110
* Beamsstrahlung photons have pT exactly zero
         IF( (phep(1,ih)**2 +phep(2,ih)) .EQ. 0d0) THEN
            nphot=nphot+1
            DO  k=1,4
               phot(nphot,k) =phep(k,ih)
            ENDDO
         ENDIF
      ENDDO
      GOTO 900
 110  CONTINUE
      RETURN
 900  WRITE(*,*) '++++ HepEvt_GetPhotBst: To many photons!!!'
      END ! HepEvt_GetPhotBst

*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//                      End of CLASS  HepEvt                                //
*//////////////////////////////////////////////////////////////////////////////
