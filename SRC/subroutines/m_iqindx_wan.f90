! inverse mapping, q to 'integer index' for wan (okumura)
module m_iqindx_wan
  implicit none
  !
  !! we prepare integer index of q, its resolution is given by epsd.
  !!      use m_hamindex, only: qtt,nqtt
  !      use m_iqindx_qtt, only: init_iqindx_qtt
  public:: Iqindx2_wan

  private
  integer,allocatable,private:: key(:,:),kk1(:),kk2(:),kk3(:),iqkkk(:,:,:)
  real(8),private:: epsd=1d-7 !key parameter to map to integer index
  integer,private:: nkey(3)
  real(8),private:: ginv(3,3)
  real(8),private,allocatable::qtt(:,:)
  integer,private::nqtt
contains
  !------------------------------------------------------
  subroutine iqindx2_wan(q,  iqindx,qu)
    intent(in) ::        q
    intent(out) ::           iqindx,qu
    !! ginv is the inverse of plat (premitive translation vector).
    !! Use kk1,kk2,kk3,nkey(1:3),iqkkk to get iqindx.
    real(8):: q(3),qu(3)
    integer:: iqindx
    integer:: i_out, iq,iqx ,kkk3(3),ik1,ik2,ik3
    real(8):: qx(3),qzz(3)
    logical::debug=.true., init=.true.
    if (init) then
       call wan_getqbz()      !!! get nqtt, qtt(3,nqtt),ginv
       call init_iqindx_wan() !ginv)
       init=.false.
    endif
    ! ccccccccccccccccccc
    !debug=.false.
    !if(abs(q(1)+0.1d0)+abs(q(2)+0.1d0)<1d-3) then
    !   debug=.true.
    !endif
    call rangedq(matmul(ginv,q), qzz)
    if(debug) write(6,"(' iqindx2_: q=',3f20.15)") qzz
    kkk3 = (qzz+0.5d0*epsd)/epsd
    ik1= findloc(kk1,value=kkk3(1),dim=1)
    ik2= findloc(kk2,value=kkk3(2),dim=1)
    ik3= findloc(kk3,value=kkk3(3),dim=1)
    iqindx = iqkkk(ik1,ik2,ik3)
!    if(debug) then
!       do iqx=1,nqtt
!          write(6,"(i5,3f13.5)")iqx,qtt(:,iqx)
!       enddo
!    endif
    qu =qtt(:,iqindx)
    if(debug) write(6,*) iqindx,qu
  end subroutine iqindx2_wan
  !!-------------------
  subroutine wan_getqbz()
    implicit none
    integer:: ifwqb
    open(newunit=ifwqb,file="wanqbz",form='unformatted')
    read(ifwqb) ginv
    read(ifwqb) nqtt
    allocate(qtt(3,nqtt))
    read(ifwqb) qtt
    close(ifwqb)
  end subroutine wan_getqbz
  !--------------------------------
  subroutine init_iqindx_wan() ! For magnon+wannier.  mapping of qtt 
    real(8):: qzz(3)
    real(8),allocatable:: qxx(:,:)
    integer:: isig,i,ix,kkk,kkk3(3),ik1,ik2,ik3,iq,ik
    integer,allocatable:: ieord(:)
    logical::debug=.false.
    allocate(ieord(nqtt))
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
       ik1= findloc(kk1,value=kkk3(1),dim=1)
       ik2= findloc(kk2,value=kkk3(2),dim=1)
       ik3= findloc(kk3,value=kkk3(3),dim=1)
       iqkkk(ik1,ik2,ik3)=i
    enddo
    deallocate(qxx)
  end subroutine init_iqindx_wan
end module m_iqindx_wan

