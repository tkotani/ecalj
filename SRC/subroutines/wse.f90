!utils for self-energy things.
subroutine winfo(ifi,nspin,nq,ntq,is,nbloch,ngp,ngc,nqbz,nqibz,ef,deltaw,alat,esmr)
  implicit none
  integer(4) :: ifi,nspin,nq,ntq,is,nbloch,ngp,ngc,nqbz,nqibz
  real(8) :: ef,deltaw,alat,esmr
  write (ifi,*)' ***'
  write (ifi,6700) nspin,nq,ntq
  write (ifi,6501) is,nbloch,ngp,ngc,nqbz,nqibz,ef,deltaw,alat,ef,esmr
6501 format (' spin =',i2,'   nbloch ngp ngc=',3i4,'  nqbz =',i6,'  nqibz =', &
       i6,'   ef=', f10.4,' Rydberg',/,d23.16,' <= deltaw(Hartree)',/,d23.16,' <= alat' &
       ,/,d23.16,' <= ef ',/,d23.16,' <= esmr')
6700 format (1x,3i4,'  nspin  nq  ntq')
end subroutine winfo
!-------------------------------------------------------------------
subroutine winfo2(ifi,nspin,nq,ntq,is,nbloch,ngp,ngc,nqbz,nqibz,ef,ef2,deltaw,alat,esmr,esmr2)
  implicit none
  integer(4) :: ifi,nspin,nq,ntq,is,nbloch,ngp,ngc,nqbz,nqibz
  real(8) :: ef,ef2,deltaw,alat,esmr,esmr2
  write (ifi,*)' ***'
  write (ifi,6700) nspin,nq,ntq
  write (ifi,6501) is,nbloch,ngp,ngc,nqbz,nqibz,ef,deltaw,alat,ef,ef2,esmr,esmr2
6501 format (' spin =',i2,'   nbloch ngp ngc=',3i4 &
       ,'  nqbz =',i6,'  nqibz =',i6,'   ef=', f10.4,' Rydberg',/,d23.16,' <= deltaw(Hartree)' &
       ,/,d23.16,' <= alat',/,2d23.15,' <= ef ',/,2d23.15,' <= esmr')
6700 format (1x,3i4,'  nspin  nq  ntq')
end subroutine winfo2
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
!! ---------------------------------------------------------
subroutine q0iwgt3(allq0i,symops,ngrp,wqt,q0i,nq0i, wgt0)
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
