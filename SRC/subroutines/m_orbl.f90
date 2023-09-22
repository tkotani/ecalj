!> orbital index system 
module m_orbl 
  use m_lmfinit,only: ltabx,ktabx,offlx,ndimxx,norbx
  public orblib,orblib1,orblib2!,gtbsl4!,orblinit
  integer,public,pointer:: ltab(:),ktab(:),offl(:),ndim,norb
  integer,public,pointer:: ltab1(:),ktab1(:),offl1(:),ndim1,norb1
  integer,public,pointer:: ltab2(:),ktab2(:),offl2(:),ndim2,norb2
  integer,allocatable,public:: ntab(:),blks(:)
  private
contains
  subroutine orblib(ia)
    integer::ia
    ktab=>ktabx(:,ia)
    ltab=>ltabx(:,ia)
    offl=>offlx(:,ia)
    ndim=>ndimxx(ia)
    norb=>norbx(ia)
    call gtbsl4(ia)
  end subroutine orblib

  subroutine orblib1(ia)
    integer::ia
    ktab1=>ktabx(:,ia)
    ltab1=>ltabx(:,ia)
    offl1=>offlx(:,ia)
    ndim1=>ndimxx(ia)
    norb1=>norbx(ia)
  end subroutine orblib1

  subroutine orblib2(ia)
    integer::ia
    ktab2=>ktabx(:,ia)
    ltab2=>ltabx(:,ia)
    offl2=>offlx(:,ia)
    ndim2=>ndimxx(ia)
    norb2=>norbx(ia)
  end subroutine orblib2
  
  subroutine gtbsl4(ia)
    !- Marks blocks of contiguous l for which rsm and e are unchanged
    ! requires l be consecutive and kappa index be constant in contiguous block
    !i   norb  :number of orbital types, i.e. kinds of radial functions 
    !i   ltab  :table of l-quantum numbers for each type (orbl.f)
    !i   ktab  :table of energy index for each type (orbl.f)
    !o Outputs
    !o   ntab  :table of upper ranges for each orbital block
    !o   blks sum of 2*l+1
    ! ----------------------------------------------------------------------
    implicit none
    integer :: ia,iorb,jorb,io,iorbe
    logical:: agree
    iorb=1
    if(allocated(ntab)) deallocate(ntab,blks)
    allocate(ntab(norb),blks(norb))
    do while(iorb<=norb)
       do jorb=iorb,norb
          agree =  ltabx(jorb,ia)==ltabx(iorb,ia)+(jorb-iorb) .and. ktabx(jorb,ia)==ktabx(iorb,ia)
          if(agree) iorbe=jorb
          if(.not.agree) exit
       enddo
       ntab(iorb)= iorbe
       ntab(iorb+1:iorbe)=0
       iorb=iorbe+1
    enddo
    blks =[(sum( [ (2*ltabx(io,ia)+1, io=iorb,ntab(iorb)) ] ),iorb=1,norb)]
  end subroutine gtbsl4
  
end module m_orbl
