!> Return QGcou and QGpsi ===
module m_readQG
  use m_read_bzdata,only: ginv
  use NaNum,only:NaN
  use m_nvfortran,only:findloc
  use m_keyvalue,only: getkeyvalue
  implicit none
  !--------------------------------------------
  public:: readqg, readqg0, readngmx, readngmx2
#ifdef __GPU
  public:: readqg_d
#endif
  integer,protected,public:: ngpmx=NaN, ngcmx=NaN!,nblochpmx
  !--------------------------------------------

  private
  real(8),allocatable,private,target:: qc(:,:),qp(:,:)
  logical,private:: init(2)=.true.
  real(8),private:: QpGcut_cou, QpGcut_psi
  integer,private,target::   nqnumc,nqnump
  integer,allocatable,private:: ngvecp(:,:,:),ngp(:),ngvecc(:,:,:),ngc(:)
  integer,pointer,private::nqtt
  real(8),pointer,private::qtt(:,:)
  real(8),private:: epsd=1d-7
  integer,private,pointer:: nkey(:),kk1(:),kk2(:),kk3(:),iqkkk(:,:,:)
  integer,target,private :: nkeyp(3),nkeyc(3)
  integer,target,allocatable,private:: keyp(:,:),kk1p(:),kk2p(:),kk3p(:),iqkkkp(:,:,:)
  integer,target,allocatable,private:: keyc(:,:),kk1c(:),kk2c(:),kk3c(:),iqkkkc(:,:,:)
  logical,private:: keepqg
  integer,private:: imx, imxc
  !      real(8),private:: ginv(3,3)
contains
  !----------------------------------
  subroutine readngmx2()
    integer:: ngmx,ifiqg
    integer :: idummy
    open(newunit=ifiqg, file='QGpsi',form='unformatted')
    read(ifiqg) nqnump, ngpmx, QpGcut_psi, idummy, idummy, imx
    close(ifiqg)
    open(newunit=ifiqg, file='QGcou',form='unformatted')
    read(ifiqg) nqnumc, ngcmx, QpGcut_cou, idummy, idummy, imxc
    close(ifiqg)
  end subroutine readngmx2

  subroutine readngmx(key,ngmx)
    intent(in) ::        key
    intent(out)::           ngmx
    !- get ngcmx or mgpmx
    ! integer:: ngmx,ifiqg,ngcmx,ngpmx
    integer:: ngmx,ifiqg
    integer :: idummy
    character*(*) key
    if    (key=='QGpsi') then
       open(newunit=ifiqg, file='QGpsi',form='unformatted')
       read(ifiqg) nqnump, ngpmx, QpGcut_psi, idummy, idummy, imx
       ngmx=ngpmx
    elseif(key=='QGcou') then
       open(newunit=ifiqg, file='QGcou',form='unformatted')
       read(ifiqg) nqnumc, ngcmx, QpGcut_cou, idummy, idummy, imxc
       ngmx=ngcmx
    else
       call rx( "readngmx: key is not QGpsi QGcou")
    endif
    close(ifiqg)
  end subroutine readngmx

  !> Get ngv and ngvec(3,ngv) for given qin(3)
  !! key=='QGcou' or 'QGpsi'
  subroutine readqg(key,qin,  qu,ngv,ngvec)
    intent(in) ::     key,qin
    intent(out)::                    qu,ngv,ngvec
    character*(*) :: key
    real(8) :: qin(3)!,ginv(3,3)
    real(8) :: qu(3)
    integer :: ngv, ngvec(3,*)
    integer:: ifi=-999999, iq,verbose
    integer:: ifiqg
    if    (key=='QGpsi') then
       ifi=1
       if(verbose()>=80) write (6,"(' readqg psi: qin=',3f8.3,i5)") qin
    elseif(key=='QGcou') then
       ifi=2
       if(verbose()>=80) write (6,"(' readqg cou: qin=',3f8.3,i5)") qin
    else
       call rx( "readqg: wrongkey")
    endif
    if(init(ifi)) then
       call init_readqg(ifi)
       init(ifi)=.false.
    endif
    if(verbose()>=40) write(6,*)'end of init_readqg'
    call iqindx2qg(qin,ifi, iq,qu)
    if(ifi==1) then
       ngv  = ngp(iq)
       if(keepqg) then
         ngvec(1:3,1:ngv) = ngvecp(1:3,1:ngv,iq)
       else
         open(newunit=ifiqg, file='QGpsi_rec',form='unformatted', access='direct', &
              recl=4*(3*ngpmx+(imx*2+1)**3), status='old')
         read(ifiqg,rec=iq) ngvec(1:3,1:ngv)
         close(ifiqg)
       endif
       return
    elseif(ifi==2) then
       ngv  = ngc(iq)
       if(keepqg) then
         ngvec(1:3,1:ngv) = ngvecc(1:3,1:ngv,iq)
       else
         open(newunit=ifiqg, file='QGcou_rec',form='unformatted', access='direct', &
              recl=4*(3*ngcmx+(2*imxc+1)**3), status='old')
         read(ifiqg,rec=iq) ngvec(1:3,1:ngv)
         close(ifiqg)
       endif
       return
    endif
    call rx( "readqg: can not find QGpsi or QPcou for given q")
  end subroutine readqg

  !> Get ngv and ngvec(3,ngv) for given qin(3)
  !! key=='QGcou' or 'QGpsi'
  subroutine readqg_d(key, qin, qu, ngv, ngvec_d)
    implicit none
    character*(*), intent(in) :: key
    real(8), intent(in) :: qin(3)
    real(8), intent(out) :: qu(3)
    integer, intent(out) :: ngv, ngvec_d(3,*)
    integer :: ifi, iq, ifiqg
    integer, allocatable :: ngvec_h(:,:)
#ifdef __GPU
    attributes(device) :: ngvec_d
#endif
    ifi = -999999
    if (key=='QGpsi') then
      ifi=1
    elseif(key=='QGcou') then
      ifi=2
    else
      call rxi( "readqg: wrongkey ifi", ifi)
    endif

    if(init(ifi)) then
       call init_readqg(ifi)
       init(ifi)=.false.
    endif
    call iqindx2qg(qin,ifi, iq,qu)
    if(ifi==1) then
       ngv  = ngp(iq)
       if(keepqg) then
         !$acc kernels present(ngvecp)
         ngvec_d(1:3,1:ngv) = ngvecp(1:3,1:ngv,iq)
         !$acc end kernels
       else
         allocate(ngvec_h(3,ngv))
         open(newunit=ifiqg, file='QGpsi_rec',form='unformatted', access='direct', &
              recl=4*(3*ngpmx+(imx*2+1)**3), status='old')
         read(ifiqg,rec=iq) ngvec_h(1:3,1:ngv)
         close(ifiqg)
         ngvec_d(1:3,1:ngv) = ngvec_h(1:3,1:ngv)
         deallocate(ngvec_h)
       endif
       return
    elseif(ifi==2) then
       ngv  = ngc(iq)
       if(keepqg) then
         !$acc kernels present(ngvecc)
         ngvec_d(1:3,1:ngv) = ngvecc(1:3,1:ngv,iq)
         !$acc end kernels
       else
         allocate(ngvec_h(3,ngv))
         open(newunit=ifiqg, file='QGcou_rec',form='unformatted', access='direct', &
              recl=4*(3*ngcmx+(2*imxc+1)**3), status='old')
         read(ifiqg,rec=iq) ngvec_h(1:3,1:ngv)
         close(ifiqg)
         ngvec_d(1:3,1:ngv) = ngvec_h(1:3,1:ngv)
         deallocate(ngvec_h)
       endif
       return
    endif
    call rx( "readqg: can not find QGpsi or QPcou for given q")
  end subroutine readqg_d

  !> Get ngv
  !! key=='QGcou' or 'QGpsi'
  subroutine readqg0(key,qin,qu,ngv)
    intent(in)::       key,qin
    intent(out)::               qu,ngv
    character*(*) :: key
    integer :: ngv
    real(8):: qin(3),ginv(3,3)
    real(8):: qu(3)
    integer:: ifi, iq,verbose
    if    (key=='QGpsi') then
       ifi=1
       if(verbose()>=80) write (6,"('readqg0 psi: qin=',3f8.3,i5)") qin
    elseif(key=='QGcou') then
       ifi=2
       if(verbose()>=80) write (6,"('readqg0 cou: qin=',3f8.3,i5)") qin
    else
       call rx( "readqg: wrongkey")
    endif
    if(init(ifi)) then
       call init_readqg(ifi)
       init(ifi)=.false.
    endif
    call iqindx2qg(qin,ifi, iq,qu)
    if(ifi==1) then
       ngv  = ngp(iq)
       if(verbose()>=80) write(6,*)'ngp=',ngv
    elseif(ifi==2) then
       ngv  = ngc(iq)
       if(verbose()>=80) write(6,*)'ngc=',ngv
    endif
    return
    call rx( "readqg0: can not find QGpsi or QPcou for given q")
  end subroutine readqg0

  !> initialization. readin QGpsi or QGcou.
  subroutine init_readqg(ifi)
    integer, intent(in) :: ifi
    ! integer:: ifiqg,iq,verbose,ngcmx,ngpmx
    integer:: ifiqg,iq,verbose
    real(8)::qq(3)
    real(8),allocatable:: qxx(:,:)
    integer:: isig,i,ix,kkk,kkk3(3),ik1(1),ik2(1),ik3(1),ik
    integer,allocatable:: ieord(:),key(:,:)
    write(6,*)' init_readqg ifi=',ifi
    call getkeyvalue("GWinput","KeepQG",keepqg,default=.true.)
    if(.not.keepqg) write(6,*) 'keepQG = .false. in readqg'
    if(ifi==1) then
       open(newunit=ifiqg, file='QGpsi',form='unformatted')
       read(ifiqg) nqnump, ngpmx, QpGcut_psi
       if(verbose()>49) write(6,"('init_readqg ngnumc ngcmx QpGcut_psi=',2i5,f8.3)") &
            nqnump, ngpmx, QpGcut_psi
       ! allocate(ngvecp(3,ngpmx,nqnump),qp(3,nqnump),ngp(nqnump))
       allocate(qp(3,nqnump),ngp(nqnump))
       if(keepqg) allocate(ngvecp(3,ngpmx,nqnump))
       do iq=1, nqnump
          read (ifiqg) qp(1:3,iq), ngp(iq)
          if(keepqg) then
            read (ifiqg) ngvecp(1:3,1:ngp(iq),iq)
          else
            read (ifiqg)
          endif
          if(verbose()>40) write(6,"('init_readqg psi qp ngp =',3f8.3,i5)") qp(1:3,iq),ngp(iq)
       enddo
       if(keepqg) then
          !$acc enter data copyin(ngvecp)
       endif
    elseif(ifi==2) then
       open(newunit=ifiqg, file='QGcou',form='unformatted')
       read(ifiqg) nqnumc, ngcmx, QpGcut_cou
       ! allocate(ngvecc(3,ngcmx,nqnumc),qc(3,nqnumc),ngc(nqnumc))
       allocate(qc(3,nqnumc),ngc(nqnumc))
       if(keepqg) allocate(ngvecc(3,ngcmx,nqnumc))
       do iq=1, nqnumc
          read(ifiqg) qc(1:3,iq), ngc(iq)
          write (6,"('init_readqg cou  qc ngc =',3f8.3,i5)") qc(1:3,iq), ngc(iq)
          if(keepqg) then
            read (ifiqg) ngvecc(1:3,1:ngc(iq),iq)
          else
            read (ifiqg)
          endif
       enddo
       if(keepqg) then
          !$acc enter data copyin(ngvecc)
       endif
    endif
    close(ifiqg)
    !! === mapping of qtt ===
    !! nkey, kk1,kk2,kk3, iqkkk are to get iqindx.
    !!  q --> call rangedq(matmul(ginv,q), qx) ---> n= (qx+0.5*epsd)/epsd
    if(ifi==1) then
       nqtt => nqnump
       qtt  => qp
       nkey => nkeyp
    elseif(ifi==2) then
       nqtt => nqnumc
       qtt  => qc
       nkey => nkeyc
    endif
    allocate(ieord(nqtt))
    allocate(key(3,nqtt),qxx(3,nqtt))
    key=-99999
    do iq=1,nqtt
       call rangedq(matmul(ginv,qtt(:,iq)), qxx(:,iq))
    enddo
    call getqkey(qxx(1,:),nqtt,epsd, nkey(1),key(1,:))
    call getqkey(qxx(2,:),nqtt,epsd, nkey(2),key(2,:))
    call getqkey(qxx(3,:),nqtt,epsd, nkey(3),key(3,:))
    !!  key is reallocated. inverse mattping, iqkkk
    if(ifi==1) then
       allocate( kk1p(nkey(1)),kk2p(nkey(2)),kk3p(nkey(3)) )
       allocate( iqkkkp(nkey(1),nkey(2),nkey(3)) )
       iqkkk => iqkkkp
       kk1 =>kk1p
       kk2 =>kk2p
       kk3 =>kk3p
    elseif(ifi==2) then
       allocate( kk1c(nkey(1)),kk2c(nkey(2)),kk3c(nkey(3)) )
       allocate( iqkkkc(nkey(1),nkey(2),nkey(3)) )
       iqkkk => iqkkkc
       kk1 =>kk1c
       kk2 =>kk2c
       kk3 =>kk3c
    endif
    kk1(:) = key(1,1:nkey(1))
    kk2(:) = key(2,1:nkey(2))
    kk3(:) = key(3,1:nkey(3))
    deallocate(key)
    do i=1,nqtt
       kkk3= (qxx(:,i)+0.5*epsd)/epsd !kkk is digitized by 1/epsd
       ik1= findloc(kk1,value=kkk3(1))
       ik2= findloc(kk2,value=kkk3(2))
       ik3= findloc(kk3,value=kkk3(3))
       iqkkk(ik1(1),ik2(1),ik3(1))=i
    enddo
    deallocate(qxx)
  end subroutine init_readqg
  subroutine iqindx2qg(q,ifi, iqindx,qu)! Find index as q=qq(:,iq) with modulo of premitive vec.
    integer, intent(in):: ifi
    integer, intent(out):: iqindx
    real(8), intent(in) :: q(3)
    real(8), intent(out) :: qu(3)
    integer:: i_out, iq,iqx ,kkk3(3),ik1(1),ik2(1),ik3(1)
    real(8):: qx(3),qzz(3)
    logical::debug=.false.
    if(ifi==1) then       !    nqtt => nqnump
       qtt  => qp
       nkey => nkeyp
       iqkkk => iqkkkp
       kk1 =>kk1p
       kk2 =>kk2p
       kk3 =>kk3p
    elseif(ifi==2) then   !    nqtt => nqnumc
       qtt  => qc
       nkey => nkeyc
       iqkkk => iqkkkc
       kk1 =>kk1c
       kk2 =>kk2c
       kk3 =>kk3c
    endif
    call rangedq(matmul(ginv,q), qzz)
    kkk3 = (qzz+0.5*epsd)/epsd
    ik1= findloc(kk1,value=kkk3(1))
    ik2= findloc(kk2,value=kkk3(2))
    ik3= findloc(kk3,value=kkk3(3))
    iqindx = iqkkk(ik1(1),ik2(1),ik3(1))
    qu = qtt(:,iqindx)
  end subroutine iqindx2qg
end module m_readQG
