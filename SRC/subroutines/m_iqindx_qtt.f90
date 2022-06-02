module m_iqindx_qtt
  !! To accelarate inverse mapping, q to 'integer index of q',
  !! We prepare integer index of q, its resolution is given by epsd.
  use m_hamindex,only: qtt,nqtt
  use m_read_bzdata,only: ginv
  implicit none

  public:: Iqindx2_, Init_iqindx_qtt

  private
  integer,allocatable,private:: key(:,:),kk1(:),kk2(:),kk3(:),iqkkk(:,:,:)
  integer,private:: nkey(3)
  !      real(8),private:: ginv(3,3)
  real(8),private:: epsd=1d-7 !key parameter to map to integer index

contains
  !!
  subroutine iqindx2_(q, iqindx,qu)
    intent(in)::        q
    intent(out)::          iqindx,qu ! qu(i) = q(i) + matmul(qlat(i,:)* nxx(:))
    !> Find index as q=qq(:,iq) with modulo of premitive vector.===
    !! ginv is the inverse of plat (premitive translation vector).
    !! Use kk1,kk2,kk3,nkey(1:3),iqkkk to get iqindx.
    real(8) :: q(3)
    integer :: iqindx
    real(8) :: qu(3)
    integer:: i_out, iq,iqx ,kkk3(3),ik1,ik2,ik3
    real(8):: qx(3),qzz(3)
    call rangedq(matmul(ginv,q), qzz)
    !! we generate qzz integer index for qzz
    kkk3 = (qzz+0.5d0*epsd)/epsd
    call tabkk(kkk3(1), kk1,nkey(1), ik1)
    call tabkk(kkk3(2), kk2,nkey(2), ik2)
    call tabkk(kkk3(3), kk3,nkey(3), ik3)
    iqindx = iqkkk(ik1,ik2,ik3)
    qu =qtt(:,iqindx)
  end subroutine iqindx2_
  !!----
  subroutine init_iqindx_qtt()
    !! === mapping of qtt ===
    !! nkey, kk1,kk2,kk3, iqkkk are to get iqindx.
    !!  q --> call rangedq(matmul(ginv,q), qx) ---> n= (qx+0.5*epsd)/epsd
    !!       --->  ik1,ik2,ik3= tabkk(kkk,iqk,nkey) ---> iqkkk(ik1,ik2,ik3)
    real(8):: qzz(3)
    real(8),allocatable:: qxx(:,:)
    integer:: isig,i,ix,kkk,kkk3(3),ik1,ik2,ik3,iq,ik
    integer,allocatable:: ieord(:)
    logical::debug=.false.
    !      ginv=ginv_
    allocate(ieord(nqtt))
    if(debug) write(6,"(a,2i5,20f9.4)")' iiiiii nqtt=',nqtt,size(qtt),ginv(1:3,1:3)
    allocate(key(3,0:nqtt),qxx(3,nqtt))
    key=-99999
    do iq=1,nqtt
       call rangedq(matmul(ginv,qtt(:,iq)), qxx(:,iq))
    enddo
    !! get key and nkey for each ix.
    key(:,0)=0 !dummy
    do ix =1,3
       call sortea(qxx(ix,:),ieord,nqtt,isig)
       ik=0
       do i=1,nqtt
          kkk=(qxx(ix,ieord(i))+0.5d0*epsd)/epsd  !kkk is digitized by 1/epsd
          if(i==1 .OR. key(ix,ik)<kkk) then
             ik=ik+1
             key(ix,ik) = kkk
          elseif (key(ix,ik)>kkk) then
             write(6,*)ix, ik,i, key(ix,ik), qxx(ix,ieord(i))
             call rx( 'iqindx: bug not sorted well')
          endif
       enddo
       nkey(ix)=ik
    enddo
    deallocate(ieord)
    !!  key is reallocated. inverse mattping, iqkkk
    allocate( kk1(nkey(1)),kk2(nkey(2)),kk3(nkey(3)) )
    kk1(:) = key(1,1:nkey(1))
    kk2(:) = key(2,1:nkey(2))
    kk3(:) = key(3,1:nkey(3))
    deallocate(key)
    allocate( iqkkk(nkey(1),nkey(2),nkey(3)) )
    iqkkk=-99999
    do i=1,nqtt
       kkk3= (qxx(:,i)+0.5d0*epsd)/epsd !kkk is digitized by 1/epsd
       call tabkk(kkk3(1), kk1,nkey(1), ik1)
       call tabkk(kkk3(2), kk2,nkey(2), ik2)
       call tabkk(kkk3(3), kk3,nkey(3), ik3)
       iqkkk(ik1,ik2,ik3)=i
    enddo
    deallocate(qxx)
  end subroutine init_iqindx_qtt

  !! ---
  subroutine tabkk(kkin, kktable,n, nout)
    intent(in)::     kkin, kktable,n
    intent(out)::                     nout
    integer:: nout,n, kkin, kktable(n),i,mm,i1,i2
    i1=1
    i2=n
    if(kkin==kktable(1)) then
       nout=1
       return
    elseif(kkin==kktable(n)) then
       nout=n
       return
    endif
    do i=1,n
       mm=(i1+i2)/2
       if(kkin==kktable(mm)) then
          nout=mm
          return
       elseif(kkin>kktable(mm)) then
          i1=mm
       else
          i2=mm
       endif
    enddo
    write(6,*) 'xxxxx takk ', i1,i2,kkin
    write(6,*) 'xxxxx takk ',kktable(i1),kktable(i2)
    call rx( 'takk: error')
  end subroutine tabkk
end module m_iqindx_qtt
