subroutine tbhsi(nspec,nermx,net,et,ipet,nrt,rt,iprt,ltop)
  use m_uspecb,only:uspecb
  use m_lmfinit,only: nkaphh,lhh
  !- Table of orbital envelope energies and smoothing radii
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   sspec :struct containing species-specific information
  !i     Passed to: uspecb
  !i   nspec  :number of species
  !i   nermx :maximum allowed size of et,rt
  !o Outputs
  !o   net   :size of table et
  !o   et    :table of all inequivalent energies
  !o   ipet  :index to which entry in et a given orbital belongs
  !o   nrt   :size of table rt
  !o   rt    :table of all inequivalent smoothing radii
  !o   iprt  :index to which entry in rt a given orbital belongs
  !o   ltop  :largest l at any site
  !r Remarks
  !r   An orbital (l,ik,is) has smoothing radius index ir=iprt(l,ik,is)
  !r   and energy index ie=ipet(l,ik,is).  Its smoothing radius and
  !r   energy are rt(ir) and et(ie), respectively.
  !u Updates
  !u   12 Aug 04 First implementation of extended local orbitals
  !u   10 Apr 02 Redimensionsed ipet,iprt to accomodate larger lmax
  !u   18 May 00 Adapted from nfp tb_hsi.f
  ! ----------------------------------------------------------------------
  implicit none
  integer :: nspec,nermx,net,nrt,ltop,n0,nkap0
  parameter (n0=10,nkap0=3)
  integer :: ipet(n0,nkap0,nspec),iprt(n0,nkap0,nspec)
  double precision :: et(nermx),rt(nermx)
  integer :: is,j,k
  integer :: lh(nkap0),nkape,ik,l
  double precision :: x1,x2,rsmh(n0,nkap0),eh(n0,nkap0)
  logical :: dcmpre
  dcmpre(x1,x2) = dabs(x1-x2) .lt. 1d-8
  net = 0
  nrt = 0
  ltop = -1
  ! --- Loop over orbitals (is,io) ---
  do  is = 1, nspec
     call uspecb(is,rsmh,eh)
     do  ik = 1, nkaphh(is) !nkape
        ltop = max0(ltop,lhh(ik,is))
        do  l = 0, lhh(ik,is)
           if (rsmh(l+1,ik) > 0) then
              !     ... Find eh, or add it to list
              j = 0
              do  k = 1, net
                 if (dcmpre(eh(l+1,ik),et(k))) then
                    j = k
                    goto 31
                 endif
              enddo
31            continue
              if (j > 0) then
                 ipet(l+1,ik,is) = j
              else
                 net = net+1
                 if (net > nermx) call rx('tbhsi: nermx exceeded for et')
                 et(net) = eh(l+1,ik)
                 ipet(l+1,ik,is) = net
              endif
              !     ... Find rmsh, or add it to list
              j = 0
              do  k = 1, nrt
                 if (dcmpre(rsmh(l+1,ik),rt(k))) then
                    j = k
                    goto 32
                 endif
              enddo
32            continue
              if (j > 0) then
                 iprt(l+1,ik,is) = j
              else
                 nrt = nrt+1
                 if (nrt > nermx) call rx('tbhsi: nermx exceeded for rt')
                 rt(nrt) = rsmh(l+1,ik)
                 iprt(l+1,ik,is) = nrt
              endif
           endif
        enddo
     enddo
  enddo
end subroutine tbhsi

