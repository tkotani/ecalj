module m_orbl
  use m_lmfinit,only: ltabx,ktabx,offlx,ndimxx,norbx
  public orblib,orblib1,orblib2!,orblinit
  integer,public,pointer:: ltab(:),ktab(:),offl(:),ndim,norb
  integer,public,pointer:: ltab1(:),ktab1(:),offl1(:),ndim1,norb1
  integer,public,pointer:: ltab2(:),ktab2(:),offl2(:),ndim2,norb2
  private
contains
  subroutine orblib(ia)
    integer::ia
    ktab=>ktabx(:,ia)
    ltab=>ltabx(:,ia)
    offl=>offlx(:,ia)
    ndim=>ndimxx(ia)
    norb=>norbx(ia)
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

end module m_orbl
