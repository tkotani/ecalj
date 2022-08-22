subroutine bzints(n1,n2,n3,ep,wp,nq,nband,nbmx,nsp,emin,emax,dos,nr,ef,job,ntet,idtet,sumev,sumwp)
  use m_lgunit,only:stdo
  !- BZ integrations by linear method
  ! ----------------------------------------------------------------------
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
  integer :: n1,n2,n3,nq,nband,nbmx,nsp,idtet(0:4,*),nr,job,ntet
  double precision :: ep(nbmx,nsp,nq),dos(nr,nsp),wp(nband,nsp,nq),emin,emax,ef,sumev,sumwp
  integer :: ib,iq,iq1,iq2,iq3,iq4,isp,itet,jjob
  integer :: ipr
  double precision :: ec(4),wc(4,2),ebot,etop,sev1,sev2,sumwm, volwgt
  character(10):: i2char
  call getpr(ipr)
  jjob = iabs(job)
  if(job < 0 .AND. nsp == 2 .OR. jjob /= 1 .AND. jjob /= 2) &
       call rx('Exit -1 BZINTS: job='//trim(i2char(job))//' and nsp='//trim(i2char(nsp))//' not allowed')
  if (jjob == 1) call dpzero(dos,nsp*nr)
  if (jjob == 2) call dpzero(wp,nband*nsp*nq)
  sev1 = 0d0
  sev2 = 0d0
  volwgt = dble(3-nsp)/(n1*n2*n3*6)
  if (job < 0) volwgt = volwgt/2
  do  40  isp = 1, nsp
     ! --- Loop over tetrahedra ---
     do  201  itet = 1, ntet
        iq1 = idtet(1,itet)
        iq2 = idtet(2,itet)
        iq3 = idtet(3,itet)
        iq4 = idtet(4,itet)
        do  20  ib = 1, nband
           ! --- Set up energies at 4 corners of tetrahedron ---
           ec(1) = ep(ib,isp,iq1)
           ec(2) = ep(ib,isp,iq2)
           ec(3) = ep(ib,isp,iq3)
           ec(4) = ep(ib,isp,iq4)
           etop = dmax1(ec(1),ec(2),ec(3),ec(4))
           ebot = dmin1(ec(1),ec(2),ec(3),ec(4))
           if (jjob == 1) then
              if ( ebot < emax ) &
                   call slinz(volwgt*idtet(0,itet),ec,emin,emax,dos(1,isp),nr)
           else
              if ( ef >= ebot ) then
                 call fswgts(volwgt*idtet(0,itet),ec,ef,etop,wc)
                 sev1 = sev1 + wc(1,1)*ec(1) + wc(2,1)*ec(2) + &
                      wc(3,1)*ec(3) + wc(4,1)*ec(4)
                 sev2 = sev2 + wc(1,2)*ec(1) + wc(2,2)*ec(2) + &
                      wc(3,2)*ec(3) + wc(4,2)*ec(4)
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
     do    isp = 1, nsp
        sumwm = 0d0
        do    ib = 1, nband
           do    iq = 1, nq
              sumwm = sumwm + wp(ib,isp,iq)
              sumwp = sumwp + wp(ib,isp,iq)
           enddo
        enddo
     enddo
     if ( ipr >= 10 ) then
        !     ... when bands are coupled, moment from bands makes no sense
        if (nsp == 1) then
           write(stdo,922) ef,sumwp,sev1+sev2,sev2
           !           write(stdl,922) ef,sumwp,sev1+sev2, sev2
        elseif (nsp == 2 .AND. job > 0) then
           sumwm = sumwp-2*sumwm
           write(stdo,923) ef,sumwp,sumwm,sev1+sev2,sev2
        endif
        if (dabs(sumwp-nint(sumwp)) > 1d-6) then
           write(stdo,924)
        endif
     endif
     !      if (nsp .eq. 1) call awrit4('bzi tetra ef %,6;6d  q %,6;6d  '//
     !     .    'sev %,6;6d  Blo %,6;6d',' ',80,stdl,ef,sumwp,sev1+sev2,sev2)
     !      if (nsp .eq. 2) call awrit5('bzi tetra ef %,6;6d  q %,6;6d  mom'//
     !     .  ' %,6;6d  sev %,6;6d  Blo %,6;6d',' ',80,stdl,ef,sumwp,sumwm,
     !     .  sev1+sev2,sev2)
  endif
922 format(1x,'BZINTS: Fermi energy:',f14.6,';',F11.6,' electrons'/ &
       9x,'Sum occ. bands:',f12.6, ', incl. Bloechl correction:',f10.6)
923 format(' BZINTS: Fermi energy:',f14.6,';',F11.6,' electrons', &
       '  amom=',f8.4/9x,'Sum occ. bands:',f12.6,', incl. Bloechl correction:',f10.6)
924 format(' (warning): non-integral number of electrons ---', &
       ' possible band crossing at E_f')
end subroutine bzints


