subroutine praldm(ifi,ipr1,ipr2,sharm,nbas,nsp,lmaxu,lldau,sspec, &
     ssite,strn,dmatu)
  use m_struc_def  !Cgetarg
  !- Writes out a site density-matrix-like object for all sites
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ifi   :if zero, write to stdo, in screen style format
  !i         :else, write to file ifi in high-precision format
  !i   ipr1  :if verbosity ipr>ipr1, print header
  !i   ipr2  :if verbosity ipr>ipr2, print contents of dmats
  !i   sharm :0 if in real harmonics, 1 if in spherical harmonics
  !i   nbas  :size of basis
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   lmaxu :dimensioning parameter for U matrix
  !i   lldau :lldau(ib)=0 => no U on this site otherwise
  !i          U on site ib with dmat in dmats(*,lldau(ib))
  !i   sspec :struct for species-specific information; see routine uspec
  !i     Elts read: lmxa idu
  !i     Stored:
  !i     Passed to:
  !i   ssite :struct for site-specific information; see routine usite
  !i     Elts read: spec
  !i     Stored:
  !i     Passed to:
  !i   strn  :string put into header
  !i   dmatu :density matrix for LDA+U
  !o Outputs
  !o   dmatu is written to file ifi
  !l Local variables
  !l         :
  !r Remarks
  !r
  !u Updates
  !u   27 Jan 06 First created
  ! ----------------------------------------------------------------------
  implicit none
  integer :: nbas,nsp,lldau(nbas),ifi,lmaxu,ipr1,ipr2,sharm,i_copy_size
  type(s_spec)::sspec(*)
  type(s_site)::ssite(*)
  real(8)::   dmatu(2,-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,*)
  character strn*(*)
  integer :: iblu,ib,is,igetss,lmxa,idu(4),l
  iblu = 0
  do  ib = 1, nbas
     if (lldau(ib) /= 0) then
        is = int(ssite(ib)%spec)
        lmxa=sspec(is)%lmxa
        idu =sspec(is)%idu
        do  l = 0, min(lmxa,3)
           if (idu(l+1) /= 0) then
              iblu = iblu+1
              call prdmts(ifi,ipr1,ipr2,sharm,strn,ib,l,lmaxu,iblu,dmatu, &
                   nsp,1)
           endif
        enddo
     endif
  enddo
end subroutine praldm

subroutine prdmts(ifi,ipr1,ipr2,sharm,strn,ib,l,lmaxu,iblu,dmats, &
     nsp,nspc)
  use m_ftox
  use m_lmfinit,only: stdo
  !- Writes out a site density-matrix-like object for a single l
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ifi   :if zero, write to stdo, in screen style format
  !i         :else, write to file ifi in high-precision format
  !i   ipr1  :if verbosity ipr>ipr1, print header
  !i   ipr2  :if verbosity ipr>ipr2, print contents of dmats
  !i   sharm :0 if in real harmonics, 1 if in spherical harmonics
  !i   strn  :string put into header
  !i   ib    :site index (ib=0 suppresses printout)
  !i   l     :dmats defined for l block
  !i   lmaxu :dimensioning parameter for dmats
  !i   iblu  :index to current block
  !i   dmats :site density matrix
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
  !o Outputs
  !o  header and dmats are printed to stdo
  !l Local variables
  !l         :
  !r Remarks
  !r
  !u Updates
  !u   09 Nov 05 dmats changed to a complex matrix
  !u   02 Jun 05 First created
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: ifi,ipr1,ipr2,l,ib,lmaxu,nsp,nspc,iblu,sharm
  double precision :: dmats(2,-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,iblu)
  character strn*(*)
  integer :: isp,ipr,m1,m2,mpipid,nlm
  character strnl*120,strn1*30,lll*60
  if (nspc == 2) call rx('prdmts not ready for nspc=2')
  if (mpipid(1) /= 0) return
  call getpr(ipr)
  nlm = 2*l+1
  if (ifi /= 0) then
     if(sharm==0) lll='complex rharm'
     if(sharm/=0) lll='complex sharm'
     do  isp = 1, nsp
        if (ipr >= ipr1) &
             write(ifi,ftox)'% rows ',nlm,' cols ',nlm,' ', &
             trim(lll),trim(strn),'l=',l,'site',ib,'spin',isp
        if (ipr < ipr2) return
        do  m1 = -l, l
           write(ifi,'(7(f12.7,2x))') (dmats(1,m1,m2,isp,iblu),m2=-l,l)
        enddo
        write(ifi,'(1x)')
        do  m1 = -l, l
           write(ifi,'(7(f12.7,2x))') (dmats(2,m1,m2,isp,iblu),m2=-l,l)
        enddo
     enddo
  else
     if(sharm == 0) then
        strnl = strn // ' real harmonics'
     else
        strnl = strn // ' spherical harmonics'
     endif
     !       Header: printout l, ib (if ib>0), spin (if nsp=2, ipr2>=ipr)
     do  isp = 1, nsp
        if (ipr < ipr2) return
        write(stdo,ftox) trim(strnl),'l=',l,'ib',ib,'isp',isp
        do  m1 = -l, l
           write(stdo,'(7(f9.5,2x))')(dmats(1,m1,m2,isp,iblu),m2=-l,l)
        enddo
        write(stdo,'(1x)')
        do  m1 = -l, l
           write(stdo,'(7(f9.5,2x))')(dmats(2,m1,m2,isp,iblu),m2=-l,l)
        enddo
     enddo
  endif
end subroutine prdmts


