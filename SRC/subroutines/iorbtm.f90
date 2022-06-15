subroutine iorbtm()
  use m_lmfinit,only: ssite=>v_ssite,sspec=>v_sspec,nl,nsp,nbas,slabl
  use m_struc_def           !Cgetarg
  use m_bandcal,only: orbtm=>orbtm_rv
  use m_lgunit,only:stdo
  !- Printout of orbital moments
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   sspec :struct containing species-specific information
  !i     Elts read: name
  !i   ics   :species table: class ic belongs to species ics(ic)
  !i   nl    :(global maximum l) + 1
  !i   nclass:number of inequivalent classes
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   orbtm :orbital moments
  !o Outputs
  !u Updates
  !u   09 Aug 04 (A. Chantis) Correct sign of orbl
  !u   08 Dec 00 First implementation
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  !      integer ics(nbas)
  !      real(8):: orbtm(nl,nsp,nbas) !nclass)
  !      type(s_spec)::sspec(*)

  ! ... Local parameters
  integer :: isp,l,im,lm,m,l1,ipr,is,ibas
  double precision :: amom,orbl(10)
  !character(8) :: slabl
  ! ... External calls
  ! ino      external dpzero,getpr,spacks
  external dpzero,getpr
  call getpr(ipr)
  if (ipr < 20) return
  !      stdo = lgunit(1)
  write(stdo,332)
332 format(/' IORBTM:  orbital moments :'/ &
       ' ibas  Spec        spin   Moment decomposed by l ...')
  do  ibas = 1,nbas
     is = ssite(ibas)%spec
     !        icyy = ibas
     !slabl=sspec(is)%name
     amom = 0
     do  isp = 1, nsp
        call dpzero(orbl,nl)
        lm = 0
        do  l = 0, nl-1
           l1 = l+1
           im = l
           if (nl == nl) im = 0
           do  m = -im, im
              lm = lm+1
              !              print *, l,m,isp,ibas,orbtm(lm,isp,ibas)
              orbl(l1) = orbl(l1) + orbtm(lm,isp,ibas)
              amom = amom + orbtm(lm,isp,ibas)
           enddo
        enddo
        write(stdo,333) ibas,slabl(is),isp,(orbl(l1),l1=1,nl)
333     format(i5,4x,a8,i6,8f12.6)
     enddo
     write(stdo,334) ibas, amom
334  format(' total orbital moment',i4,':',f12.6)
  enddo

end subroutine iorbtm


