subroutine bzints(n1n2n3,ep,wp,nq,nband,nbmx,nsp,emin,emax,dos,nr,ef,job,ntet,idtet,sumev,sumwp)!- BZ integrations by linear method
  use m_lgunit,only:stdo
  !i Inputs
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
  integer :: n1n2n3,nq,nband,nbmx,nsp,idtet(0:4,*),nr,job,ntet
  double precision :: ep(nbmx,nsp,nq),dos(nr,nsp),wp(nband,nsp,nq),emin,emax,ef,sumev,sumwp
  integer :: ib,iq,iq1,iq2,iq3,iq4,isp,itet,jjob
  integer :: ipr
  double precision :: ec(4),wc(4,2),ebot,etop,sev1,sev2,sumwm, volwgt
  character(10):: i2char
  call getpr(ipr)
  jjob = iabs(job)
  if(job < 0 .AND. nsp == 2 .OR. jjob /= 1 .AND. jjob /= 2) &
       call rx('Exit -1 BZINTS: job='//trim(i2char(job))//' and nsp='//trim(i2char(nsp))//' not allowed')
  if (jjob == 1) dos=0d0
  if (jjob == 2) wp=0d0
  sev1 = 0d0
  sev2 = 0d0
  volwgt = dble(3d0-nsp)/(n1n2n3*6d0)
  if (job < 0) volwgt = volwgt/2d0
  do  40  isp = 1, nsp
     ! --- Loop over tetrahedra ---
     do  201  itet = 1, ntet
        iq1 = idtet(1,itet)
        iq2 = idtet(2,itet)
        iq3 = idtet(3,itet)
        iq4 = idtet(4,itet)
        do  20  ib = 1, nband! --- Set up energies at 4 corners of tetrahedron ---
           ec(1:4) = ep(ib,isp,[iq1,iq2,iq3,iq4])
           etop = dmax1(ec(1),ec(2),ec(3),ec(4))
           ebot = dmin1(ec(1),ec(2),ec(3),ec(4))
           if (jjob == 1) then
              if( ebot < emax ) call slinz(volwgt*idtet(0,itet),ec,emin,emax,dos(1,isp),nr)
           else
              if( ef >= ebot ) then
                 call fswgts(volwgt*idtet(0,itet),ec,ef,etop,wc)
                 sev1 = sev1 + wc(1,1)*ec(1) + wc(2,1)*ec(2) + wc(3,1)*ec(3) + wc(4,1)*ec(4)
                 sev2 = sev2 + wc(1,2)*ec(1) + wc(2,2)*ec(2) + wc(3,2)*ec(3) + wc(4,2)*ec(4)
                 wp(ib,isp,iq1) = wp(ib,isp,iq1) + wc(1,1) + wc(1,2)
                 wp(ib,isp,iq2) = wp(ib,isp,iq2) + wc(2,1) + wc(2,2)
                 wp(ib,isp,iq3) = wp(ib,isp,iq3) + wc(3,1) + wc(3,2)
                 wp(ib,isp,iq4) = wp(ib,isp,iq4) + wc(4,1) + wc(4,2)
              endif
           endif
20      enddo
201  enddo
40 enddo
  if (jjob == 2) then
     sumev = sev1 + sev2
     sumwp = 0d0
     do   isp = 1, nsp
        sumwm = sum(wp(1:nband,isp,1:nq))
        sumwp = sumwp + sum(wp(1:nband,isp,1:nq))
     enddo
     if ( ipr >= 10 ) then ! when bands are coupled, moment from bands makes no sense
        if (nsp == 1) then
           write(stdo,922) ef,sumwp,sev1+sev2,sev2
        elseif (nsp == 2 .AND. job > 0) then
           sumwm = sumwp-2*sumwm
           write(stdo,923) ef,sumwp,sumwm,sev1+sev2,sev2
        endif
        if (dabs(sumwp-nint(sumwp)) > 1d-6) write(stdo,924)
     endif
  endif
922 format(1x,'BZINTS: Fermi energy:',f14.6,';',F11.6,' electrons'/ &
       9x,'Sum occ. bands:',f12.6, ', incl. Bloechl correction:',f10.6)
923 format(' BZINTS: Fermi energy:',f14.6,';',F11.6,' electrons', &
       '  amom=',f8.4/9x,'Sum occ. bands:',f12.6,', incl. Bloechl correction:',f10.6)
924 format(' (warning): non-integral number of electrons ---',' possible band crossing at E_f')
end subroutine bzints
subroutine fswgts(volwgt,e,ef,etop,w)
  implicit none
  real(8):: e(4),ef,volwgt,etop,w(4,2),wx(4,2),efm,efp,kbt
  call fswgts_(volwgt,e,ef,etop,w)
  return
  ! kbt=0.003d0 !room temperature? (I tested but little make congergence smoother)
  ! call fswgts_(volwgt,e,ef-2*kbt,etop,wx); w=1d0/10d0*wx
  ! call fswgts_(volwgt,e,ef-kbt,etop,wx);   w=w+2d0/10d0*wx
  ! call fswgts_(volwgt,e,ef,etop,wx);       w=w+4d0/10d0*wx
  ! call fswgts_(volwgt,e,ef+kbt,etop,wx);   w=w+2d0/10d0*wx
  ! call fswgts_(volwgt,e,ef+2*kbt,etop,wx); w=w+1d0/10d0*wx
end subroutine fswgts
subroutine fswgts_(volwgt,e,ef,etop,w)
  !- Makes weights for integration up to Ef for one tetrahedron.
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
3 enddo
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
end subroutine fswgts_



