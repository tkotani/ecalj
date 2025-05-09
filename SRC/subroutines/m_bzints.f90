module m_bzints
  public bzints,bzints2x
  private
contains
  subroutine bzints(n1n2n3,ep,wp,nq,nband,nbmx,nspx,nsp,emin,emax, dos,nr,ef,job,ntet,idtet,sumev,sumwp,spinweightsoc)!- BZ integrations by linear method
    use m_lgunit,only:stdo
    use m_lmfinit,only: lso!,nspx
    use m_ftox
    !i   nq    :no. of irr. k-points
    !i   ep    :energy bands
    !i   nband :number of bands
    !i   nbmx  :dimensions ep
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   emin  :(job=1) lower bound for energy
    !i   emax  :(job=1) upper bound for energy
    !i   nr    :(job=1) number of points for idos
    !i   job   :(job=1) makes idos;
    !i         :(job=2) makes Bloechl weights.
    !i         :sign of job<0:  spin up, spin down bands are coupled.
    !i         :In this case, nsp must be unity!
    !i   ntet  :No. of different tetrahedra
    !i   idtet :idtet(0,i) = number of tetrahedra of the i'th kind
    !i         :idtet(1-4,i) points to the 4 irr. k-points defining tetr.
    ! o Inputs/Outputs
    ! o  ef    :Fermi energy
    ! o        :job=1: output
    ! o        :job=2: input
    !o Outputs
    !o   sumev :sum of eigenvalues (job = 2)
    !o   dos   :Integrated density of states (idos) (job = 1)
    !o   wp    :Bloechl quadrature weights (job = 2)
    !u Updates
    !u   17 Jan 05 Returns sumwp
    ! ----------------------------------------------------------------------
    implicit none
    integer :: n1n2n3,nq,nband,nbmx,nsp,idtet(0:4,*),nr,job,ntet,nspx
    real(8) :: dos(nr,nsp),wp(nband,nsp,nq),emin,emax,ef,sumev,sumwp,ep(nbmx,nspx,nq) 
    integer :: ib,iq,iq1,iq2,iq3,iq4,isp,itet,jjob,ispx
    integer :: ipr,iqq(4)
    real(8) :: ec(4),wc(4,2),ebot,etop,sev1,sev2,sumwm, volwgt,wt
    character(10):: i2char
    real(8),optional:: spinweightsoc(:,:,:)
    real(8),allocatable::epp(:,:,:)
    call getpr(ipr)
    jjob = iabs(job)
    if(job<0.AND.nsp==2.OR.jjob/=1.AND.jjob /= 2) &
         call rx('BZINTS: job='//trim(i2char(job))//' and nsp='//trim(i2char(nsp))//' not allowed')
    if (jjob == 1) dos=0d0
    if (jjob == 2) wp=0d0
    sev1 = 0d0
    sev2 = 0d0
    volwgt = dble(3d0-nsp)/(n1n2n3*6d0)
    if(job<0) volwgt = volwgt/2d0
    allocate( epp,source=reshape(ep,[nbmx,merge(nspx,nsp,lso==1),nq]))
    isploop: do  40  isp = 1, nsp 
       LoopOverTetrahedra: do  201  itet = 1, ntet
          iqq = idtet(1:4,itet)
          ibandloop: do  20  ib = 1, nband
             ispx = merge(1,isp,lso==1)
             ec(1:4) = epp(ib,ispx,iqq(1:4)) ! --- Set up energies at 4 corners of tetrahedron ---
             etop = maxval(ec)
             ebot = minval(ec)
             if (jjob == 1) then !so=1 for 2024-5-10
                wt = 1d0
                if(lso==1.and.present(spinweightsoc)) wt = sum(spinweightsoc(ib,isp,idtet(1:4,itet)))/4d0
                if( ebot < emax ) call slinz(wt*volwgt*idtet(0,itet),ec,emin,emax,dos(1,isp),nr)
             else
                if( ef >= ebot ) then
                   call fswgts(volwgt*idtet(0,itet),ec,ef,etop,wc)
                   sev1 = sev1 + sum(ec*wc(:,1)) !wc(1,1)*ec(1) + wc(2,1)*ec(2) + wc(3,1)*ec(3) + wc(4,1)*ec(4)
                   sev2 = sev2 + sum(ec*wc(:,2)) !wc(1,2)*ec(1) + wc(2,2)*ec(2) + wc(3,2)*ec(3) + wc(4,2)*ec(4)
                   wp(ib,isp,iqq(1:4)) = wp(ib,isp,iqq(1:4)) + wc(1:4,1) + wc(1:4,2)
                endif
             endif
20        enddo ibandloop
201    enddo LoopOverTetrahedra
40  enddo isploop
    deallocate(epp)
    if(jjob == 2) then
       sumev = sev1 + sev2
       sumwp = sum(wp(1:nband,1:nsp,1:nq))
       if(nsp==2) sumwm=sum(wp(1:nband,1,1:nq))-sum(wp(1:nband,2,1:nq))
       printout:if(ipr>= 10 ) then ! when bands are coupled, moment from bands makes no sense
          if (nsp == 1) then
             write(stdo,922) ef,sumwp,sev1+sev2,sev2
          elseif (nsp == 2 .AND. job > 0) then
             write(stdo,923) ef,sumwp,sumwm,sev1+sev2,sev2
          endif
          if (dabs(sumwp-nint(sumwp)) > 1d-6) write(stdo,924)
       endif printout
    endif
922 format(' BZINTS: Fermi energy:',f14.6,';',F11.6,' electrons'/9x,'Sum occ. bands:',f12.6, ', incl. Bloechl correction:',f10.6)
923 format(' BZINTS: Fermi energy:',f14.6,';',F11.6,' electrons','  amom=',f8.4/9x,'Sum occ. bands:',f12.6,'&
         , incl. Bloechl correction:',f10.6)
924 format(' (warning): non-integral number of electrons ---',' possible band crossing at E_f')
  end subroutine bzints
  ! subroutine fswgts(volwgt,e,ef,etop,w)
  !   implicit none
  !   real(8):: e(4),ef,volwgt,etop,w(4,2),wx(4,2),efm,efp,kbt
  !   call fswgts_(volwgt,e,ef,etop,w)
  !   return
  !   ! kbt=0.003d0 !room temperature? (I tested but little make congergence smoother)
  !   ! call fswgts_(volwgt,e,ef-2*kbt,etop,wx); w=1d0/10d0*wx
  !   ! call fswgts_(volwgt,e,ef-kbt,etop,wx);   w=w+2d0/10d0*wx
  !   ! call fswgts_(volwgt,e,ef,etop,wx);       w=w+4d0/10d0*wx
  !   ! call fswgts_(volwgt,e,ef+kbt,etop,wx);   w=w+2d0/10d0*wx
  !   ! call fswgts_(volwgt,e,ef+2*kbt,etop,wx); w=w+1d0/10d0*wx
  ! end subroutine fswgts
  subroutine fswgts(volwgt,e,ef,etop,w)    !- Makes weights for integration up to Ef for one tetrahedron.
    ! ----------------------------------------------------------------
    !i Inputs
    !i
    !o Outputs
    !o   w
    !r Remarks
    !r   w(i,1): normal weights.  w(i,2): Bloechl-correction.
    ! ----------------------------------------------------------------
    implicit none
    real(8):: e(4),ef,volwgt,etop,w(4,2)
    real(8)::ex(4),w1(4),w2(4),a14,a21,a24,a31,a32,a34,a41,a42,df=0d0,v1,v2,v3,vw4,xxx
    real(8),target:: ec(4)
    real(8),pointer:: e1,e2,e3,e4
    integer :: i,i00,j,m,n,isort(4)
    !  complex(8):: img=(0d0,1d0)
    !  real(8),parameter:: eps=1d-3
    vw4 = volwgt/4d0
    w=0d0
    if(ef >= etop) then
       w(:,1)=vw4
       return
    endif
    ! --- Sort energies into ec ---
    ex = e
    do  3  i = 1, 4
       i00 = 1
       do j = 2, 4
          if (ex(j) < ex(i00)) i00 = j
       enddo
       ec(i) = ex(i00)
       isort(i) = i00
       ex(i00) = etop + 1d0
3   enddo
    e1 => ec(1)
    e2 => ec(2)
    e3 => ec(3)
    e4 => ec(4)
    ! --- Case Ef between e2,e3 ---
    if (e2 < ef .AND. ef <= e3) then
       a31 = (ef-e1)/(e3-e1)
       a41 = (ef-e1)/(e4-e1)
       a32 = (ef-e2)/(e3-e2)
       a42 = (ef-e2)/(e4-e2)
       v1 = a31*a41
       v2 = a31*a42*(1-a41)
       v3 = a42*a32*(1-a31)
       w1(1) = (v1*(3-a31-a41) + v2*(2-a31-a41) + v3*(1-a31))*vw4
       w1(2) = (v1 + v2*(2-a42) + v3*(3-a32-a42))*vw4
       w1(3) = (v1*a31 + v2*a31 + v3*(a31+a32))*vw4
       w1(4) = (v1*a41 + v2*(a41+a42) + v3*a42)*vw4
       df = ((e1+e2-e3-e4)*a32*a42 + 2*ef - e1 - e2)/((e3-e1)*(e4-e1))
       df = 3*volwgt*df
       ! --- Case Ef between e1,e2 ---
    else if (e1 < ef .AND. ef <= e2) then
       a21 = (ef-e1)/(e2-e1)
       a31 = (ef-e1)/(e3-e1)
       a41 = (ef-e1)/(e4-e1)
       xxx = a21*a31*a41*vw4
       w1(1) = xxx*(4-a21-a31-a41)
       w1(2) = xxx*a21
       w1(3) = xxx*a31
       w1(4) = xxx*a41
       df = 3*volwgt*a31*a41/(e2-e1)
       ! --- Case Ef between e3,e4 ---
    else if (e3 < ef .AND. ef <= e4) then
       a14 = (ef-e4)/(e1-e4)
       a24 = (ef-e4)/(e2-e4)
       a34 = (ef-e4)/(e3-e4)
       xxx = a14*a24*a34*vw4
       w1(1) = vw4 - xxx*a14
       w1(2) = vw4 - xxx*a24
       w1(3) = vw4 - xxx*a34
       w1(4) = vw4 - xxx*(4-a14-a24-a34)
       df = -3*volwgt*a14*a24/(e3-e4)
    endif
    ! --- Bloechl correction ---
    w2=0d0
    do     m = 1, 4
       do  n = 1, 4
          w2(m) = w2(m) + (ec(n)-ec(m))*df/40
       enddo
    enddo
    ! ----------------------------------------------
    do  i = 1, 4
       j = isort(i)
       w(j,1) = w1(i)
       w(j,2) = w2(i)
    enddo
  end subroutine fswgts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE BZINTS2x(volwgt,EP,WP,NQ,nband,NB,NSP,EMIN,EMAX,DOS,NR,EF,JOB,NTET,IDTET)!-  Bz integrations by linear method.
    IMPLICIT real(8) (A-H,O-Z)
    implicit integer (i-n)
    DIMENSION EP(nb,nsp,nq),DOS(NR),EC(4),WC(4,2),WP(nband,nsp,nq), idtet(0:4,*)
    integer::iqq(4)
    IF (JOB /= 1 .AND. JOB /= 2) STOP '*** BAD JOB IN BZINTS2x'
    IF (JOB == 1) dos=0d0 !call dinit(dos,nr)
    IF (JOB == 2) wp=0d0  !call dinit(wp,nband*nsp*nq)
    SEV1 = 0.D0
    SEV2 = 0.D0     !      volwgt = (3.d0 - nsp) / (n1*n2*n3*6.d0)
    do  40  isp = 1, nsp
       DO  201  ITET = 1, NTET
          iqq = idtet(1:4,itet)
          DO  20  IB = 1, nb !nband
             ec  = ep(ib,isp,iqq)
             etop = maxval(ec)
             ebot = minval(ec)
             IF (JOB == 1) THEN
                if ( ebot < emax ) CALL SLINZ2(volwgt*idtet(0,ITET),EC,EMIN,EMAX,DOS,NR)
             ELSE
                if ( ef >= ebot ) then
                   CALL FSWGTS(volwgt*idtet(0,ITET),EC,EF,ETOP,WC)
                   SEV1 = SEV1 + sum(wc(:,1)*ec) 
                   SEV2 = SEV2 + sum(wc(:,2)*ec) 
                   WP(ib,isp,iqq) = WP(ib,isp,iqq) + WC(1:4,1) + WC(1:4,2)
                endif
             ENDIF
20        enddo
201    enddo
40  enddo
  end SUBROUTINE BZINTS2x
  SUBROUTINE SLINZ2(VOLWGT,EC,EMIN,EMAX,DOSI,NR)    !- Adds to number-of-states for one tetrahedron
    !i Inputs
    !i   volwgt, weight on tetrahedron; ec energies at corners of tethdn.;
    !i   emin, emax, energy window; nr, number of bins + 1
    !o Outputs
    !o   dosi(k), integrated density in kth bin from tethdn.
    IMPLICIT real(8) (A-H,O-Z)
    implicit integer (i-n)
    DIMENSION EC(4),DOSI(NR)
    DO   I = 1, 3
       DO   J = 1, 4-I
          IF (EC(J) > EC(J+1)) THEN
             E=EC(J)
             EC(J) = EC(J+1)
             EC(J+1) = E
          ENDIF
       enddo
    enddo
    E1 = EC(1)
    E2 = EC(2)
    E3 = EC(3)
    E4 = EC(4)
    if (e4 < emin) then
       i4=1
       go to 26
    endif
    DE = (EMAX-EMIN)/(NR-1)
    d2 = 2.0d0*( 1.0d0-d1mach(3) )
    !x takao ---------------
    !x This is a correction in order to get the very accurate Fermi energy.
    !x
    E1x=(E1-EMIN)/DE
    if( e1x >= dble(nr+10) ) e1x = nr + 10
    if( e1x <= -10.0d0)      e1x =  - 10
    I01 = e1x + d2

    E2x=(E2-EMIN)/DE
    if( e2x >= dble(nr+10) ) e2x = nr + 10
    if( e2x <= -10.0d0)      e2x =  - 10
    I02 = e2x + d2

    E3x=(E3-EMIN)/DE
    if( e3x >= dble(nr+10) ) e3x = nr + 10
    if( e3x <= -10.0d0)      e3x =  - 10
    I03 = e3x + d2

    E4x=(E4-EMIN)/DE
    if( e4x >= dble(nr+10) ) e4x = nr + 10
    if( e4x <= -10.0d0)      e4x =  - 10
    I04 = e4x + d2
    !x this was sometimes poor and cause a problem when you want to get very accurate Ef.
    !      write(6,"('EEEE=',4d13.6)") E1, E2, E3, E4
    !      print *, I01, I02, I03, I04
    !      I01   = (E1   -EMIN)/DE + 1.9999999D0
    !      I02   = (E2   -EMIN)/DE + 1.9999999D0
    !      I03   = (E3   -EMIN)/DE + 1.9999999D0
    !      I04   = (E4   -EMIN)/DE + 1.9999999D0

    ! --------------------------------
    I1 = MAX0(I01  ,1)
    I2 = MIN0(I02-1, NR)
    IF (I1 <= I2) THEN
       CC = VOLWGT/((E2-E1)*(E3-E1)*(E4-E1))
       DO  20  I = I1, I2
          X = EMIN - E1 + (I-1)*DE
          DOSI(I) = DOSI(I) + CC*X**3
20     enddo
    ENDIF
    I2 = MAX0(I02  ,1)
    I3 = MIN0(I03-1,NR)
    IF (I2 <= I3) THEN
       C3 = VOLWGT*(E1+E2-E3-E4)/((E3-E1)*(E4-E1)*(E3-E2)*(E4-E2))
       C2 = VOLWGT*3.D0/((E3-E1)*(E4-E1))
       C1 = C2*(E2-E1)
       C0 = C1*(E2-E1)/3.D0
       DO  21  I = I2, I3
          X = EMIN - E2 + (I-1)*DE
          DOSI(I) = DOSI(I) + C0 + X*(C1 + X*(C2 + X*C3))
21     enddo
    ENDIF
    I3 = MAX0(I03  ,1)
    I4 = MIN0(I04  -1,NR)
    IF (I3 <= I4) THEN
       CC = VOLWGT/((E3-E4)*(E2-E4)*(E1-E4))
       DO  22  I = I3, I4
          X = EMIN - E4 + (I-1)*DE
          DOSI(I) = DOSI(I) + VOLWGT - CC*X**3
22     enddo
    ENDIF
    I4 = MAX0(I04  ,1)
26  continue
    DO  25  I = I4, NR
       DOSI(I) = DOSI(I) + VOLWGT
25  enddo
  END SUBROUTINE SLINZ2
endmodule m_bzints

  ! !--------------------------------------------------
  ! subroutine getvaln(konfig,z,nl,natom,iclass,nclass, valn)
  !   ! - Get valn
  !   !o valn    = number of valence electron.

  !   implicit none
  !   integer(4):: nclass,natom,nl,ia,ic,l
  !   real(8)   :: valn,ef, z(nclass)
  !   integer(4):: iclass(natom),konfig(0:nl-1,nclass)
  !   valn    = 0d0
  !   do ia   = 1,natom
  !      ic    = iclass(ia)
  !      valn  = valn + z(ic)
  !      do    l = 0,nl-1
  !         valn  = valn - (konfig(l,ic)-l-1) *( 2*l +1)*2
  !      end do
  !   end do
  !   print *,' getvaln: valn=',valn
  ! end subroutine getvaln
