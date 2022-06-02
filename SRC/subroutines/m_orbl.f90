module m_orbl
  public orblib,orblib1,orblib2,orblinit
  integer,public,pointer:: ltab(:),ktab(:),offl(:),ndim,norb
  integer,public,pointer:: ltab1(:),ktab1(:),offl1(:),ndim1,norb1
  integer,public,pointer:: ltab2(:),ktab2(:),offl2(:),ndim2,norb2

  integer,parameter,public :: n0=10,nkap0=3,n00=n0*nkap0
  integer,allocatable,public,target:: ltabx(:,:),ktabx(:,:),offlx(:,:),ndimx(:),norbx(:)

  private
  logical :: init=.true.
contains
  subroutine orblib(ia)
    integer::ia
    if(init) call orblinit()
    ktab=>ktabx(:,ia)
    ltab=>ltabx(:,ia)
    offl=>offlx(:,ia)
    ndim=>ndimx(ia)
    norb=>norbx(ia)
  end subroutine orblib

  subroutine orblib1(ia)
    integer::ia
    if(init) call orblinit()
    ktab1=>ktabx(:,ia)
    ltab1=>ltabx(:,ia)
    offl1=>offlx(:,ia)
    ndim1=>ndimx(ia)
    norb1=>norbx(ia)
  end subroutine orblib1

  subroutine orblib2(ia)
    integer::ia
    if(init) call orblinit()
    ktab2=>ktabx(:,ia)
    ltab2=>ltabx(:,ia)
    offl2=>offlx(:,ia)
    ndim2=>ndimx(ia)
    norb2=>norbx(ia)
  end subroutine orblib2

  subroutine orblinit() !norb,ltab,ktab,off,offl,ndim)
    use m_lmfinit,only: nl,nkaph,iprmb,nlmto,nbas
    !- Extract a list of l's in given hamiltonian block for one site
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   ib    :site for which to extract list of l's and k's
    !o Outputs
    !o   norb  :number of orbital types for ib; see Remarks
    !o   ltab  :table of l-quantum numbers for each type
    !o   ktab  :table of energy index for each type
    !o   offl  :offl(norb) offset in h to this block of orbitals
    !o   ndim  :dimension of hamiltonian for this site
    !r Remarks
    !r   Each orbital type is label by a 'l' and a 'k' index
    !r   Each orbital corresponds to a unique radial wave function at the
    !r   site where the orbit is centered.  There can be multiple 'k'
    !r   indices (radial wave function shapes) for a particular l.
    !r
    implicit none
    integer:: ib,ik,l,lmr,ia
    if(init) then
       allocate(ltabx(n00,nbas),ktabx(n00,nbas),offlx(n00,nbas),ndimx(nbas),norbx(nbas))
       norbx=0
       ndimx=0
       do ib=1,nbas
          lmr = nl*nl*nkaph*(ib-1)
          do  ik = 1, nkaph
             do  l = 0, nl-1
                offlx(norbx(ib)+1,ib) = -1
                if (iprmb(lmr+1) >0 .AND. iprmb(lmr+1) <= nlmto) then
                   offlx(norbx(ib)+1,ib) = iprmb(lmr+1) - 1
                   norbx(ib) = norbx(ib) +1
                   ndimx(ib) = ndimx(ib) + 2*l+1
                   ltabx(norbx(ib),ib) = l
                   ktabx(norbx(ib),ib) = ik
                endif
                lmr = lmr + 2*l+1
             enddo
          enddo
          if (norbx(ib) > n00) call rx('orbl: norb> n00')
       enddo
       init=.false.
    endif
  end subroutine orblinit

end module m_orbl
