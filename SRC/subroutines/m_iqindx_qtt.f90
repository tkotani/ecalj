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
  real(8),private:: epsd=1d-7 !key parameter to map to integer index
contains
  subroutine iqindx2_(q, iqindx,qu) !Find index for q=qq(:,iqindx).Modulo of premitive vector.
    intent(in)::      q
    intent(out)::        iqindx,qu ! qu(i) = q(i) + matmul(qlat(i,:)* nxx(:))
    !! ginv is the inverse of plat (premitive translation vector).
    real(8) :: q(3),qu(3),qx(3),qzz(3)
    integer :: iqindx, kkk3(3), ik1(1),ik2(1),ik3(1)
    call rangedq(matmul(ginv,q), qzz) ! we generate qzz integer index for qzz
    kkk3 = (qzz+0.5d0*epsd)/epsd
    ik1= findloc(kkk3(1)-kk1,value=0)
    ik2= findloc(kkk3(2)-kk2,value=0)
    ik3= findloc(kkk3(3)-kk3,value=0)
    iqindx = iqkkk(ik1(1),ik2(1),ik3(1))
    qu =qtt(:,iqindx)
  end subroutine iqindx2_
  subroutine init_iqindx_qtt()
    !! === mapping of qtt ===
    !! nkey, kk1,kk2,kk3, iqkkk are to get iqindx.
    !!  q --> call rangedq(matmul(ginv,q), qx) ---> n= (qx+0.5*epsd)/epsd
    real(8):: qzz(3)
    integer:: isig,i,ix,kkk,kkk3(3),ik1(1),ik2(1),ik3(1),iq,ik
    real(8),allocatable:: qxx(:,:)
    integer,allocatable:: ieord(:)
    allocate(ieord(nqtt),key(3,nqtt),qxx(3,nqtt))
    do iq=1,nqtt
       call rangedq(matmul(ginv,qtt(:,iq)), qxx(:,iq))
    enddo
    !! get key and nkey for each ix.
    call getqkey(qxx(1,:),nqtt,epsd, nkey(1),key(1,:))
    call getqkey(qxx(2,:),nqtt,epsd, nkey(2),key(2,:))
    call getqkey(qxx(3,:),nqtt,epsd, nkey(3),key(3,:))
    !!  key is reallocated. inverse mapping, iqkkk
    allocate( kk1(nkey(1)),kk2(nkey(2)),kk3(nkey(3)) )
    kk1(:) = key(1,1:nkey(1))
    kk2(:) = key(2,1:nkey(2))
    kk3(:) = key(3,1:nkey(3))
    deallocate(ieord,key)
    allocate( iqkkk(nkey(1),nkey(2),nkey(3)) )
    iqkkk=-99999
    do i=1,nqtt
       kkk3= (qxx(:,i)+0.5d0*epsd)/epsd 
       ik1= findloc(kkk3(1)-kk1,value=0)
       ik2= findloc(kkk3(2)-kk2,value=0)
       ik3= findloc(kkk3(3)-kk3,value=0)
       iqkkk(ik1(1),ik2(1),ik3(1))=i
    enddo
    deallocate(qxx)
  end subroutine init_iqindx_qtt
end module m_iqindx_qtt
