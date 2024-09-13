!>Write W-v. Gamma-cell averaged W-v
module m_llw
  use m_rdpp,only: nblochpmx
  use m_genallcf_v3,only: nclass,natom,nspin,nl,nn, ndima,nlnmx, nctot,alat, deltaw,clabl,iclass, plat, pos, ecore, tpioa
  use m_freq,only: frhis,freq_r,freq_i, nwhis,nw_i,nw,npm,niw !output of getfreq
  use m_qbze,only: Setqbze, nqbze,nqibze,qbze,qibze
  use m_read_bzdata,only: Read_bzdata, ngrp2=>ngrp,nqbz,nqibz,n1,n2,n3,ginv,dq_,qbz,wbz,qibz,wibz, ntetf,idtetf,ib1bz
  use m_read_bzdata,only: qbzw,nqbzw, q0i,nq0i ,nq0iadd,ixyz
  use m_readVcoud,only: vcousq,zcousq,ngb
  use m_rdpp,only: nbloch,mrecl
  use m_x0kf,only: zxq,zxqi
  use m_mpi, only: MPI__GatherXqw, mpi__root_k, mpi__root_q
  use m_kind,only: kp => kindrcxq
  use m_stopwatch
  use m_blas, only: zminv
  implicit none
  public:: WVRllwR,WVIllwI,  MPI__sendllw,MPI__sendllw2
  complex(8),allocatable,protected,public:: llw(:,:), llwI(:,:)
  complex(8),allocatable,protected,public:: wmuk(:,:)
  logical,protected,public:: w4pmode
  integer,protected,public:: ngbq0
  private
  complex(8),allocatable:: zw0(:,:)
  complex(kind=kp), allocatable :: zw(:,:)
  complex(8),allocatable:: epstinv(:,:),epstilde(:,:)
  real(8),parameter:: pi=4d0*datan(1d0),fourpi = 4d0*pi
contains
  subroutine WVRllwR(q,iq,nmbas1,nmbas2)
    use m_readqg,only: Readqg0
    intent(in)::       q,iq,    nmbas1,nmbas2 !zxq can be twiced when nspin=2
    integer:: iq,iq0,nwmax,nwmin,iw,imode,ix,igb1,igb2,ifllw
    integer:: nmbas1,nmbas2,ngc0,ifw4p,ifrcw,mreclx
    real(8):: frr,q(3),vcou1,quu(3),eee
    logical::  localfieldcorrectionllw,cmdopt0,emptyrun
    complex(8), allocatable :: zxqw(:,:)
    logical,save:: init=.true.
    type(stopwatch) :: t_sw_matinv, t_sw_x_gather
    integer :: istat
!    complex(8):: zxq(nmbas1,nmbas2,nw_i:nw)
    character(10):: i2char
    mreclx=mrecl
    emptyrun=.false. !cmdopt0('--emptyrun')
    if(init) then !initialization related to w4pmode, zw, tpioa...
       allocate( llw(nw_i:nw,nq0i),source=(1d99,0d0) )
       if(sum(ixyz)/=0) w4pmode= .TRUE. 
       if(w4pmode) allocate( wmuk(2:nblochpmx,3),source=(1d99,0d0))
       allocate( zw(nblochpmx,nblochpmx) )
       init=.false.
       call stopwatch_init(t_sw_matinv, 'matinv')
       call stopwatch_init(t_sw_x_gather, 'gather')
    endif
    call readqg0('QGcou', (/0d0,0d0,0d0/),  quu,ngc0) ! ngb is q-dependent. released at the end of WVIllwi
    ngbq0 = nbloch+ngc0
    allocate( zw0(ngb,ngb), epstilde(ngb,ngb), epstinv(ngb,ngb))
    allocate (zxqw(ngb, ngb))
    flush(6)
    if(nspin == 1) zxq = 2d0*zxq !if paramagnetic, multiply x0 by 2
    nwmax = nw
    nwmin = nw_i
    write(6, *)" === trace check for W-V === nwmin nwmax=",nwmin,nwmax, 'qqqqqxx=',iq,q
    if(iq<=nqibz) then        !for mmmw
      if(mpi__root_q) then
        open(newunit=ifrcw, file='WVR.'//i2char(iq),form='unformatted',access='direct',recl=mreclx)
      endif
      iwloop: do 1015 iw  = nwmin,nwmax
        if(emptyrun) exit
        frr= dsign(freq_r(abs(iw)),dble(iw))
        if(iq==1) then
          ix=1
          zw0(:,1)=0d0
          zw0(1,:)=0d0
        else
          ix=0
        endif
        call stopwatch_start(t_sw_x_gather)
        call MPI__GatherXqw(zxq(:,:,iw), zxqw, nmbas1, nmbas2)
        call stopwatch_pause(t_sw_x_gather)
        if(mpi__root_q) then
        do igb1=ix+1,ngb !  Eqs.(37),(38) in PRB81 125102 (Friedlich)
          do igb2=ix+1,ngb
            ! epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
            epstilde(igb1,igb2)= -vcousq(igb1)*zxqw(igb1,igb2)*vcousq(igb2)
            if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
          enddo
        enddo
        epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
        call stopwatch_start(t_sw_matinv)
        ! call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
        !$acc data copy(epstinv)
        !$acc host_data use_device(epstinv)
        istat = zminv(epstinv(ix+1,ix+1), n=ngb-ix, lda=ngb)
        !$acc end host_data
        !$acc end data
        call stopwatch_pause(t_sw_matinv)
        !  w4p writing eps
        if(iw==0 .AND. w4pmode) then ! static epstinv is saved. For q=0 epstilde (mu=1 skipped). For q/=0 full matrix inversion. ix=1 is set for q=0)
          open(newunit=ifw4p,file='W4PHONON.'//i2char(iq),form='unformatted')
          write(ifw4p) iq,q,ngb,ix !ix=0, or ix=1 for q=0 (iq=1)
          write(ifw4p) epstinv(ix+1:ngb,ix+1:ngb)
          close(ifw4p)
        endif
        do igb1=1+ix,ngb
          do igb2=1+ix,ngb
            zw0(igb1,igb2)= vcousq(igb1)*epstinv(igb1,igb2)*vcousq(igb2)
            if(igb1==igb2) zw0(igb1,igb2)= zw0(igb1,igb2)-vcousq(igb1)*vcousq(igb2)
          enddo
        enddo
        zw(1:ngb,1:ngb) = cmplx(zw0(1:ngb,1:ngb),kind=kp)
1012    continue
        write(ifrcw, rec= iw-nw_i+1 ) zw !  WP = vsc-v
        call tr_chkwrite("freq_r iq iw realomg trwv=", zw, iw, frr,nblochpmx, nbloch,ngb,iq)
        endif
1015  enddo iwloop
      if(mpi__root_q) close(ifrcw)
      if(mpi__root_q) then 
         call stopwatch_show(t_sw_x_gather)
         call stopwatch_show(t_sw_matinv)
      endif
    else  ! llw, Wing elements of W. See PRB81 125102
      iq0 = iq - nqibz
      vcou1 = fourpi/sum(q**2*tpioa**2) ! --> vcousq(1)**2!  !fourpi/sum(q**2*tpioa**2-eee)
      do 1115 iw  = nwmin,nwmax
        if(emptyrun) exit
        frr= dsign(freq_r(abs(iw)),dble(iw))
        !! Full inversion to calculalte eps with LFC.
        !if(localfieldcorrectionllw()) then
        call MPI__GatherXqw(zxq(:,:,iw), zxqw, nmbas1, nmbas2)
        if(mpi__root_q) then
        ix=0
        eee=0d0
        do igb1=ix+1,ngb
          do igb2=ix+1,ngb
            if(igb1==1 .AND. igb2==1) then
              ! epstilde(igb1,igb2)= 1d0 - vcou1*zxq(1,1,iw)
              epstilde(igb1,igb2)= 1d0 - vcou1*zxqw(1,1)
              cycle
            endif
            ! epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
            epstilde(igb1,igb2)= -vcousq(igb1)*zxqw(igb1,igb2)*vcousq(igb2)
            if(igb1==igb2) then
              epstilde(igb1,igb2)=1d0 + epstilde(igb1,igb2)
            endif
          enddo
        enddo
        epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
        call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
        if(iq0<=nq0i) llw(iw,iq0)= 1d0/epstinv(1,1)
        !     ! Wing elements calculation july2016    ! We need check nqb is the same as that of q=0
        if(ixyz(iq0)/=0 .AND. iw==0) then
          if(ngb/=ngbq0) then
            write(6,*)q,iq0,ngb,ngbq0
            call rx('hx0p0_sc: ngb/=ngbq0')
          endif
          wmuk(2:ngb,ixyz(iq0))=epstinv(1,2:ngb)/epstinv(1,1) ! this is dot(q(:)*w_mu(:,igb)). See PRB125102(2016) eq.(36)
        endif
        !else
        !   if(iq0<=nq0i) llw(iw,iq0)= 1d0 - vcou1*zxq(1,1,iw)
        !endif
        if(iq0<=nq0i) write(6,"('epsWVR: iq iw_R omg(iw) eps(wFC) eps(woLFC) ', &
             2i5,x,10(d13.6,2x,d13.6,x,d13.6,2x,d13.6,x,d13.6))") &
             iq,iw,freq_r(iw),llw(iw,iq0),1d0-vcou1*zxqw(1,1)
             ! iq,iw,freq_r(iw),llw(iw,iq0),1d0-vcou1*zxq(1,1,iw)
        continue               !iw
        endif
1115  enddo
    endif
    deallocate(zxqw)
  end subroutine WVRllwR
  subroutine WVIllwi(q,iq,nmbas1,nmbas2)
    intent(in)::       q,iq,     nmbas1,nmbas2 !zxqi can be twiced when nspin=2
    integer:: nmbas1,nmbas2,mreclx
    integer:: iq,iq0,nwmax,nwmin,iw,imode,ix,igb1,igb2,ifllwi,ifrcwi
    real(8):: frr,q(3),vcou1
    logical::  localfieldcorrectionllw,cmdopt0,emptyrun
    logical,save:: init=.true.
    complex(8), allocatable :: zxqw(:,:)
!    complex(8):: zxqi(nmbas1,nmbas2,niw)
    character(10):: i2char
    integer :: istat
    allocate (zxqw(ngb, ngb))
    mreclx=mrecl
    emptyrun=.false. !cmdopt0('--emptyrun')
    if(init) then
       allocate( llwI(niw,nq0i) )
       init=.false.
       llwI= 1d99
    endif
    write(6,*)'WVRllwI: init'
    if (nspin == 1) zxqi = 2d0*zxqi ! if paramagnetic, multiply x0 by 2
    if( iq<=nqibz ) then
       if(mpi__root_q) then
         open(newunit=ifrcwi,file='WVI.'//i2char(iq),form='unformatted',access='direct',recl=mreclx)
       endif
       do 1016 iw  = 1,niw
          if(emptyrun) exit
          !!  Eqs.(37),(38) in PRB81 125102
          if(iq==1) then
             ix=1
             zw0(:,1)=0d0
             zw0(1,:)=0d0
          else
             ix=0
          endif
          call MPI__GatherXqw(zxqi(:,:,iw), zxqw, nmbas1, nmbas2)
          if(mpi__root_q) then
          do igb1=ix+1,ngb
             do igb2=ix+1,ngb
                ! epstilde(igb1,igb2)= -vcousq(igb1)*zxqi(igb1,igb2,iw)*vcousq(igb2)
                epstilde(igb1,igb2)= -vcousq(igb1)*zxqw(igb1,igb2)*vcousq(igb2)
                if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
             enddo
          enddo
          epstinv=epstilde
          !$acc data copy(epstinv)
          !$acc host_data use_device(epstinv)
          istat = zminv(epstinv(ix+1,ix+1), n=ngb-ix, lda=ngb)
          !$acc end host_data
          !$acc end data
          ! call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
          do igb1=ix+1,ngb
             do igb2=ix+1,ngb
                zw0(igb1,igb2)= vcousq(igb1)*epstinv(igb1,igb2)*vcousq(igb2)
                if(igb1==igb2) zw0(igb1,igb2)= zw0(igb1,igb2)-vcousq(igb1)*vcousq(igb2)
             enddo
          enddo
1014      continue
          zw(1:ngb,1:ngb) = cmplx(zw0(1:ngb,1:ngb),kind=kp)
          write(ifrcwi, rec= iw)  zw !  WP = vsc-v
          call tr_chkwrite("freq_i iq iw imgomg trwv=",zw,iw,freq_i(iw),nblochpmx,nbloch,ngb,iq)
          endif
1016   enddo
       if(mpi__root_q) close(ifrcwi)
    else
       !! Full inversion to calculalte eps with LFC.
       iq0 = iq - nqibz
       vcou1 = fourpi/sum(q**2*tpioa**2) ! --> vcousq(1)**2!  !fourpi/sum(q**2*tpioa**2-eee)
       do 1116 iw  = 1,niw
          if(emptyrun) exit
          !if(localfieldcorrectionllw()) then
          call MPI__GatherXqw(zxqi(:,:,iw), zxqw, nmbas1, nmbas2)
          if(mpi__root_q) then
             ix=0
             do igb1=ix+1,ngb
                do igb2=ix+1,ngb
                   if(igb1==1 .AND. igb2==1) then
                      ! epstilde(igb1,igb2)= 1d0 - vcou1*zxqi(1,1,iw)
                      epstilde(igb1,igb2)= 1d0 - vcou1*zxqw(1,1)
                      cycle
                   endif
                   ! epstilde(igb1,igb2)= -vcousq(igb1)*zxqi(igb1,igb2,iw)*vcousq(igb2)
                   epstilde(igb1,igb2)= -vcousq(igb1)*zxqw(igb1,igb2)*vcousq(igb2)
                   if(igb1==igb2) then
                      epstilde(igb1,igb2)=1d0 + epstilde(igb1,igb2)
                   endif
                enddo
             enddo
             epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
             call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
             if(iq0<=nq0i) llwI(iw,iq0)= 1d0/epstinv(1,1)
          !else
          !   if(iq0<=nq0i) llwI(iw,iq0)=  1d0 -vcou1*zxqi(1,1,iw)
          !endif
          ! if(iq0<=nq0i) write(6,"('iq iw_img eps(wLFC) eps(noLFC)',i4,i4,2f10.4,2x,2f10.4)") &
          !      iq,iw,llwI(iw,iq0),1d0-vcou1*zxqi(1,1,iw)
             if(iq0<=nq0i) write(6,"('iq iw_img eps(wLFC) eps(noLFC)',i4,i4,2f10.4,2x,2f10.4)") &
                  iq,iw,llwI(iw,iq0),1d0-vcou1*zxqw(1,1)
          endif
1116   enddo
    endif
    deallocate(epstinv,epstilde,zw0)
    deallocate(zxqw)
  end subroutine WVILLWI
  subroutine MPI__sendllw2(iqxend,MPI__ranktab) !for hx0fp0
    use m_mpi,only: MPI__root,MPI__DbleCOMPLEXsend,MPI__DbleCOMPLEXrecv,MPI__rank,MPI__size
    intent(in)::             iqxend
    integer:: iq0,dest,src,iq,iqxend,MPI__ranktab(:)
    !! === Recieve llw and llwI at node 0, where q=0(iq=1) is calculated. ===
    if(MPI__size==1) return
    do iq=nqibz+1,iqxend
      iq0 = iq - nqibz
      if(MPI__ranktab(iq)==0) cycle 
      if(MPI__ranktab(iq) == MPI__rank) then
        dest=0
        call MPI__DbleCOMPLEXsend(llw(nw_i,iq0),(nw-nw_i+1),dest)
        call MPI__DbleCOMPLEXsend(llwI(1,iq0),niw,dest)
      elseif(MPI__root) then
        src=MPI__ranktab(iq)
        call MPI__DbleCOMPLEXrecv(llw(nw_i,iq0),(nw-nw_i+1),src)
        call MPI__DbleCOMPLEXrecv(llwI(1,iq0),niw,src)
      endif
    enddo
  end subroutine MPI__sendllw2
  subroutine MPI__sendllw(iqxend,MPI__Qranktab) !for hx0fp0_sc
    use m_mpi,only: MPI__DbleCOMPLEXsendQ,MPI__DbleCOMPLEXrecvQ,MPI__size,MPI__rank,MPI__root
    ! === Recieve llw and llwI at node 0, where q=0(iq=1) is calculated. ===
    intent(in)::            iqxend
    integer:: iq0,dest,src,iq,iqxend,MPI__Qranktab(:)
    if(MPI__size==1) return
    do iq=nqibz+1,iqxend
      iq0 = iq - nqibz
      if(MPI__Qranktab(iq)==0) cycle
      if(MPI__Qranktab(iq) == MPI__rank) then
        dest=0
        if(iq0<=nq0i) then
          call MPI__DbleCOMPLEXsendQ(llw(nw_i,iq0),(nw-nw_i+1),dest)
          call MPI__DbleCOMPLEXsendQ(llwI(1,iq0),niw,dest)
        endif
        if(ixyz(iq0)/=0) call MPI__DbleCOMPLEXsendQ(wmuk(2:ngbq0,ixyz(iq0)),ngbq0-1,dest)
      elseif(MPI__root) then
        src=MPI__Qranktab(iq)
        if(iq0<=nq0i) then
          call MPI__DbleCOMPLEXrecvQ(llw(nw_i,iq0),(nw-nw_i+1),src)
          call MPI__DbleCOMPLEXrecvQ(llwI(1,iq0),niw,src)
        endif
        if(ixyz(iq0)/=0) call MPI__DbleCOMPLEXrecvQ(wmuk(2:ngbq0,ixyz(iq0)),ngbq0-1,src)
      endif
    enddo
  end subroutine MPI__sendllw
end module m_llw
!===================================================================
subroutine tr_chkwrite(tagname,zw,iw,freqq,nblochpmx,nbloch,ngb,iq)
  use m_kind,only: kp => kindrcxq
  implicit none
  integer:: nblochpmx,nbloch,ngb,iw,i,iq
  complex(kind=kp):: zw(nblochpmx,nblochpmx)
  complex(8):: trwv,trwv2
  real(8):: freqq
  character*(*)::tagname
  trwv=0d0
  do i = 1,nbloch
     trwv = trwv + zw(i,i)
  enddo
  trwv2 = 0d0
  do i = 1,ngb
     trwv2 = trwv2 + zw(i,i)
  enddo  !  write(6,'(" realomg trwv=",2i6,4d22.14)') iq,iw,trwv(iw),trwv2(iw)
  write(6,'(a,f10.4,2i5,4d22.14)')tagname,freqq,iq,iw,trwv,trwv2
end subroutine tr_chkwrite
