subroutine checkeq(i,j)
  integer::i,j
  if(i/=j) call rx( " checkeq in hsfp0: dim of WVR and WVI not compatible")
end subroutine checkeq
!--------------------------------------------------------------------
subroutine rsexx2 (nspin, itq, q, ntq,nq,ginv, symgg,ng, vxco)
  implicit real*8 (a-h,o-z)
  implicit integer (i-n)
  dimension vxco(ntq,nq,nspin),q(3,nq),itq(ntq) !itq is not dependent on q, right?
  real(8),allocatable :: qqq(:,:),vxcfpx(:,:,:)
  logical ::nocore,lfind
  real(8)::  rydberg,tolq=1d-5,qx(3),ginv(3,3),qr(3),symgg(3,3,ng),sym(3,3)
  write(6,*)' OPEN VXCFP '
  open(newunit=ifvxcfp,file='VXCFP',form='unformatted')
  read(ifvxcfp) ldim,nqbz
  write(6,*)' rsexx ldim,nqbz',ldim,nqbz
  allocate(qqq(3,nqbz),vxcfpx(ldim,nqbz,nspin))
  do ikp = 1,nqbz
     read(ifvxcfp) qqq(1:3,ikp),vxcfpx(1:ldim,ikp,1:nspin)
     write(6,"(i5,100d13.5)") ikp,qqq(1:3,ikp)
  enddo
  close(ifvxcfp)
  do iq=1,nq
     do ikp=1,nqbz
        do ig = 1,ng
           sym = symgg(:,:,ig)
           qr=matmul(sym,qqq(1:3,ikp))
           lfind=.false.
           if(sum( (qr-q(1:3,iq))**2) <tolq) then
              lfind=.true.
           else
              call rangedq( matmul(ginv,q(1:3,iq)-qr), qx)
              if(sum(abs(qx))< tolq) lfind= .TRUE. 
           endif
           if(lfind) then
              ikpx=ikp
              goto 100
           endif
        enddo
     enddo
     call rx( ' rsexx: not find ikp')
100  continue
     vxco(1:ntq,iq,1:nspin)=rydberg()*vxcfpx(itq(1:ntq),ikpx,1:nspin)
  enddo
end subroutine rsexx2
!------------------------------------------------------------------------
subroutine rsexx (nspin, itq, q, ntq,nq,ginv, &
     vxco)
  implicit real*8 (a-h,o-z)
  implicit integer (i-n)
  dimension vxco(ntq,nq,nspin),q(3,nq),itq(ntq) !itq is not dependent on q, right?
  real(8),allocatable :: qqq(:,:),vxcfpx(:,:,:)
  logical ::nocore,lfind
  real(8)::  rydberg,tolq=1d-5,qx(3),ginv(3,3)
  write(6,*)' OPEN VXCFP '
  open(newunit=ifvxcfp,file='VXCFP',form='unformatted')
  read(ifvxcfp) ldim,nqbz
  write(6,*)' rsexx ldim,nqbz',ldim,nqbz
  allocate(qqq(3,nqbz),vxcfpx(ldim,nqbz,nspin))
  do ikp = 1,nqbz
     read(ifvxcfp) qqq(1:3,ikp),vxcfpx(1:ldim,ikp,1:nspin)
     write(6,"(i5,100d13.5)") ikp,qqq(1:3,ikp)
  enddo
  close(ifvxcfp)
  do iq=1,nq
     do ikp=1,nqbz
        lfind=.false.
        if(sum( (qqq(1:3,ikp)-q(1:3,iq))**2) <tolq) then
           lfind=.true.
        else
           call rangedq( matmul(ginv,q(1:3,iq)-qqq(:,ikp)), qx)
           if(sum(abs(qx))< tolq) lfind= .TRUE. 
        endif
        if(lfind) then
           ikpx=ikp
           goto 100
        endif
     enddo
     call rx( ' rsexx: not find ikp')
100  continue
     vxco(1:ntq,iq,1:nspin)=rydberg()*vxcfpx(itq(1:ntq),ikpx,1:nspin)
  enddo
end subroutine rsexx
!--------------------------------------------------------------------
double precision  function egex (q,ef)
  ! 92.02.22
  ! electron gas bare exchange
  ! SEx(q) = -(kf/pi) { 1 + ((kf^2-q^2)/2q*kf) ln|q+kf/q-kf|
  ! ef     = kf^2
  ! q      = magnitude of q-vector
  ! ef     = fermi energy in Rydberg
  ! egex   = exchange energy in Hartree
  implicit real*8 (a-h,o-z)
  implicit integer (i-n)
  pi         = 4.d0*datan(1.d0)
  qf         = dsqrt(ef)
  if(q==0d0) then
     fac =  1d0
  else
     fac =  ((qf*qf-q*q)/(2.d0*q*qf) )*dlog(dabs( (q+qf)/(q-qf) ))
  endif
  egex  = -(qf/pi)*(1d0 + fac)
  return
END function egex

subroutine winfo(ifi,nspin,nq,ntq,is,nbloch,ngp,ngc,nqbz,nqibz,ef &
     ,deltaw,alat,esmr)
  implicit none
  integer(4) :: ifi,nspin,nq,ntq,is,nbloch,ngp,ngc,nqbz,nqibz
  real(8) :: ef,deltaw,alat,esmr
  write (ifi,*)' ***'
  write (ifi,6700) nspin,nq,ntq
  write (ifi,6501) is,nbloch,ngp,ngc,nqbz,nqibz,ef &
       ,deltaw,alat,ef,esmr
6501 format (' spin =',i2,'   nbloch ngp ngc=',3i4 &
       ,'  nqbz =',i6,'  nqibz =', &
       i6,'   ef=', f10.4,' Rydberg' &
       ,/,d23.16,' <= deltaw(Hartree)' &
       ,/,d23.16,' <= alat' &
       ,/,d23.16,' <= ef ' &
       ,/,d23.16,' <= esmr')
6700 format (1x,3i4,'  nspin  nq  ntq')
end subroutine winfo
!-----------------------------------------------------------------
subroutine readxx(ifil)
  character(72) :: rchar
  integer(4):: n=1000,ifil,i
  do 10 i = 1,n
     read(ifil,'(a72)')rchar
     rchar=trim(adjustl(rchar))
     if(rchar(1:5) == '*****') return
     if(rchar(1:5) == '#####') return
10 enddo
  call rx( 'readx: cannot find the string (gwsrc/wse.f/readxx)')
end subroutine readxx
!-------------------------------------------------------------------
subroutine winfo2(ifi,nspin,nq,ntq,is,nbloch,ngp,ngc,nqbz,nqibz,ef &
     ,ef2,deltaw,alat,esmr,esmr2)
  implicit none
  integer(4) :: ifi,nspin,nq,ntq,is,nbloch,ngp,ngc,nqbz,nqibz
  real(8) :: ef,ef2,deltaw,alat,esmr,esmr2
  write (ifi,*)' ***'
  write (ifi,6700) nspin,nq,ntq
  write (ifi,6501) is,nbloch,ngp,ngc,nqbz,nqibz,ef &
       ,deltaw,alat,ef,ef2,esmr,esmr2
6501 format (' spin =',i2,'   nbloch ngp ngc=',3i4 &
       ,'  nqbz =',i6,'  nqibz =', &
       i6,'   ef=', f10.4,' Rydberg' &
       ,/,d23.16,' <= deltaw(Hartree)' &
       ,/,d23.16,' <= alat' &
       ,/,2d23.15,' <= ef ' &
       ,/,2d23.15,' <= esmr')
6700 format (1x,3i4,'  nspin  nq  ntq')
end subroutine winfo2
!! ---------------------------------------------------------
subroutine q0iwgt3(allq0i,symops,ngrp,wqt,q0i,nq0i, &
     wgt0)
  !! Get weight for each k-points near 0.
  !! wgt0(irreducible-k, irotation)
  implicit none
  logical :: allq0i !if true --> all rot. incl., even equivivalent
  ! f false -> include only nonequivivalent rot.
  integer(4) :: ixx,ix,i, ngrp,ig,nq0i
  real(8)     :: q0(3,6),q0i(3,6),symops(3,3,ngrp),sym(3,3), &
       sym2(3,3),qt(3), q0in(3,ngrp*nq0i), wgt0(nq0i,ngrp), wqt(nq0i)
  wgt0 = 1d0
  do i = 1,nq0i
     ixx = 0
     qt = q0i(:,i)
     ! equivalence check
     do ig = 1,ngrp
        sym = symops(:,:,ig)
        if (allq0i) then
           wgt0 (i,ig)  = wqt(i) /dble(ngrp)
           if(wgt0(i,ig)/=0d0) write(6,'(a,2i3,f12.5)')'allq0i=T ',i,ig, wgt0(i,ig)
        else
           do ix = 1,ig-1
              sym2=symops(:,:,ix)
              if(sum(abs(  matmul(sym2,qt)- matmul(sym,qt)  ))<1d-6) then
                 wgt0(i,ig) = 0d0
                 goto 1111
              endif
           enddo
           ixx = ixx+1
        endif
1111    continue
     enddo !ig
     if ( .NOT. allq0i) then
        do ig = 1, ngrp
           wgt0(i,ig) = wgt0(i,ig)*wqt(i)/dble(ixx)
           if(wgt0(i,ig)/=0d0) write(6,'(a,2i3, f12.5)')'allq0i=F ',i,ig, wgt0(i,ig)
        enddo !ig
     endif
  enddo !i
end subroutine q0iwgt3
!---------------------------------------------------
integer*4 function ivsumxxx(ia,n)
  integer(4) ::ia(n),n,i
  ivsumxxx=0
  do i=1,n
     if(ia(i)/=0) ivsumxxx=ivsumxxx+1
  enddo
END function ivsumxxx
