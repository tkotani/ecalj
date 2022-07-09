subroutine iorbtm()
  use m_lmfinit,only: ispec,sspec=>v_sspec,nl,nsp,nbas,slabl
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
  implicit none
  integer :: isp,l,im,lm,m,l1,ipr,is,ibas
  double precision :: amom,orbl(10)
  call getpr(ipr)
  if (ipr < 20) return
  write(stdo,332)
332 format(/'IORBTM:  orbital moments :'/' ibas  Spec        spin   Moment decomposed by l ...')
  do  ibas = 1,nbas
     is = ispec(ibas)
     amom = 0
     do  isp = 1, nsp
        orbl=0d0
        lm = 0
        do  l = 0, nl-1
           l1 = l+1
           im = l
           if (nl == nl) im = 0
           do  m = -im, im
              lm = lm+1
              orbl(l1) = orbl(l1) + orbtm(lm,isp,ibas)
              amom = amom + orbtm(lm,isp,ibas)
           enddo
        enddo
        write(stdo,"(i5,4x,a8,i6,8f12.6)") ibas,slabl(is),isp,(orbl(l1),l1=1,nl)
     enddo
     write(stdo,"(' total orbital moment',i4,':',f12.6)") ibas, amom
  enddo
end subroutine iorbtm


