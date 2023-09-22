!>generate index table for product basis.
module m_pbindex 
  Use NaNum,only: NaN
  implicit none

  integer,protected:: &
       norbt=NaN,  &!  number of PB block.
  max_ibas_tbl=NaN, max_l_tbl=NaN, max_k_tbl=NaN, max_offset_tbl=NaN
  integer,allocatable,protected:: &
       ibas_tbl(:),l_tbl(:),k_tbl(:),offset_tbl(:),offset_rev_tbl(:,:,:)

  logical,private::init=.true.
contains !-----------------------------------------------
  subroutine PBindex(natom,lx,l2nl,nx)
    intent(in)::       natom,lx,l2nl,nx
    integer::natom,l2nl,ibas,lb,nb,iorbt,kn,lx(natom),nx(0:l2nl,natom),mb
    if( .NOT. init) return
    init=.false.
    iorbt=0
    do ibas=1, natom
       do lb  = 0, lx (ibas)
          do nb  = 1, nx (lb,ibas)
             iorbt=iorbt+1
          enddo
       enddo
    enddo
    norbt = iorbt !number of product basis block
    write(6,*)'npbindex: norbt=',norbt
    allocate( ibas_tbl(norbt),l_tbl(norbt),k_tbl(norbt),offset_tbl(norbt) )
    iorbt=0
    offset_tbl(1)=0
    do ibas=1, natom
       do lb = 0, lx (ibas)
          do kn = 1, nx (lb,ibas)
             iorbt=iorbt+1
             ibas_tbl(iorbt)=ibas
             l_tbl(iorbt)=lb
             k_tbl(iorbt)=kn
             if(iorbt<norbt) offset_tbl(iorbt+1)=offset_tbl(iorbt)+2*lb+1
             !        do mb  = -lb, lb
             !          i = i+1  !The number of product basis is  =(i at the end of loop).
             !          write(6,"(' === product basis index: iorbt ibas l n m',10i4)")iorbt,ic,lb,kn
             !        enddo
          enddo
       enddo
    enddo
    max_ibas_tbl=natom
    max_l_tbl=maxval(l_tbl(1:norbt))
    max_k_tbl=maxval(k_tbl(1:norbt))
    allocate(offset_rev_tbl(max_ibas_tbl,0:max_l_tbl,max_k_tbl))
    do iorbt=1,norbt
       ibas=ibas_tbl(iorbt)
       lb= l_tbl(iorbt)
       kn= k_tbl(iorbt)
       offset_rev_tbl(ibas,lb,kn)=offset_tbl(iorbt)
       !        write(6,"(' === product basis index: iorbt ibas l n offset',10i4)")iorbt,ibas,lb,kn,offset_rev_tbl(ibas,lb,kn)
    enddo
  end subroutine PBindex
end module m_pbindex
