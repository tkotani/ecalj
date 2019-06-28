      subroutine readppovl0(q,ngc,ppovl)
      implicit none
      integer, intent(in) :: ngc
      complex(8), intent(out) :: ppovl(ngc,ngc)
      real(8), intent(in) :: q(3)
      integer:: ngc_r,ippovl0,ifile_handle
      real(8):: qx(3),tolq=1d-8
      ippovl0=ifile_handle()
      open(ippovl0,file='PPOVL0',form='unformatted')
      do
        read(ippovl0) qx,ngc_r
        if(sum(abs(qx-q))<tolq) then
          if(ngc_r/=ngc) call rx( 'readin ppovl: ngc_r/=ngc')
          read(ippovl0) ppovl
          exit
        endif
      enddo
      close(ippovl0)
      end

!> Return QGcou and QGpsi ===
      module m_readQG
      implicit none
      real(8),allocatable,private,target:: qc(:,:),qp(:,:)
      logical,private:: init(2)=.true.
      real(8),private:: QpGcut_cou, QpGcut_psi
      integer(4),private,target::   nqnumc,nqnump,ngcmx,ngpmx
      integer(4),allocatable,private:: ngvecp(:,:,:),ngp(:),ngvecc(:,:,:),ngc(:)
      integer,pointer,private::nqtt
      real(8),pointer,private::qtt(:,:)
      real(8),private:: epsd=1d-7
      integer,private,pointer:: nkey(:),kk1(:),kk2(:),kk3(:),iqkkk(:,:,:)
      integer,target,private :: nkeyp(3),nkeyc(3)
      integer,target,allocatable,private:: keyp(:,:),kk1p(:),kk2p(:),kk3p(:),iqkkkp(:,:,:)
      integer,target,allocatable,private:: keyc(:,:),kk1c(:),kk2c(:),kk3c(:),iqkkkc(:,:,:)
      real(8),private:: ginv_(3,3)
      contains
c----------------------------------
      subroutine readngmx(key,ngmx)
c- get ngcmx or mgpmx
      implicit none
      integer(4):: ngmx,ifile_handle,ifiqg
      character*(*) key
      ifiqg=ifile_handle()
      if    (key=='QGpsi') then
        open(ifiqg, file='QGpsi',form='unformatted')
        read(ifiqg) nqnump, ngpmx, QpGcut_psi
        ngmx=ngpmx
      elseif(key=='QGcou') then
        open(ifiqg, file='QGcou',form='unformatted')
        read(ifiqg) nqnumc, ngcmx, QpGcut_cou
        ngmx=ngcmx
      else
        call rx( "readngmx: key is not QGpsi QGcou")
      endif
      close(ifiqg)
      end subroutine

!> Get ngv and ngvec(3,ngv) for given qin(3)
!! key=='QGcou' or 'QGpsi'
      subroutine readqg(key,qin,ginv,  qu,ngv,ngvec)
      implicit none
      character*(*), intent(in) :: key
      real(8), intent(in) :: qin(3),ginv(3,3)
      real(8), intent(out) :: qu(3)
      integer(4), intent(out) :: ngv, ngvec(3,*)

      integer(4):: ifi, iq,verbose
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
        call init_readqg(ifi,ginv)
        init(ifi)=.false.
      endif
      if(verbose()>=40) write(6,*)'end of init_readqg'
      call iqindx2qg(qin,ifi, iq,qu)
      if(ifi==1) then
        ngv  = ngp(iq)
        ngvec(1:3,1:ngv) = ngvecp(1:3,1:ngv,iq)
        return
      elseif(ifi==2) then
        ngv  = ngc(iq)
        ngvec(1:3,1:ngv) = ngvecc(1:3,1:ngv,iq)
        return
      endif
      call rx( "readqg: can not find QGpsi or QPcou for given q")
      end subroutine readqg

!> Get ngv
!! key=='QGcou' or 'QGpsi'
      subroutine readqg0(key,qin,ginv,  qu,ngv)
      implicit none
      character*(*), intent(in) :: key
      integer(4), intent(out) :: ngv
      real(8), intent(in):: qin(3),ginv(3,3)
      real(8), intent(out):: qu(3)

      integer(4):: ifi, iq,verbose
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
        call init_readqg(ifi,ginv)
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
      end subroutine

!> initialization. readin QGpsi or QGcou.
      subroutine init_readqg(ifi,ginv)
      implicit none
      integer(4), intent(in) :: ifi
      real(8), intent(in) :: ginv(3,3)

      integer(4):: ifiqg,iq,verbose,ifile_handle
      real(8)::qq(3)
      real(8),allocatable:: qxx(:,:)
      integer:: isig,i,ix,kkk,kkk3(3),ik1,ik2,ik3,ik
      integer,allocatable:: ieord(:),key(:,:)
      ginv_=ginv
      write(6,*)' init_readqg ifi=',ifi
      ifiqg=ifile_handle()
      if(ifi==1) then
        open(ifiqg, file='QGpsi',form='unformatted')
        read(ifiqg) nqnump, ngpmx, QpGcut_psi
        if(verbose()>49) 
     &     write (6,"('init_readqg ngnumc ngcmx QpGcut_psi=',2i5,f8.3)") 
     &     nqnump, ngpmx, QpGcut_psi
        allocate(ngvecp(3,ngpmx,nqnump),qp(3,nqnump),ngp(nqnump))
        do iq=1, nqnump
          read (ifiqg) qp(1:3,iq), ngp(iq)
          read (ifiqg) ngvecp(1:3,1:ngp(iq),iq)
          if(verbose()>40) 
     &     write (6,"('init_readqg psi qp ngp =',3f8.3,i5)") qp(1:3,iq),ngp(iq)
        enddo
      elseif(ifi==2) then
        open(ifiqg, file='QGcou',form='unformatted')
        read(ifiqg) nqnumc, ngcmx, QpGcut_cou
c         write (6,"('init_readqg ngnumc ngcmx QpGcut_cou=',2i5,f8.3)")
c     &     nqnumc, ngcmx, QpGcut_cou
        allocate(ngvecc(3,ngcmx,nqnumc),qc(3,nqnumc),ngc(nqnumc))
        do iq=1, nqnumc
          read(ifiqg) qc(1:3,iq), ngc(iq)
c          if(verbose()>40)  write (6,"('init_readqg cou  qc ngc =',3f8.3,i5)") qc(1:3,iq), ngc(iq)
          write (6,"('init_readqg cou  qc ngc =',3f8.3,i5)") qc(1:3,iq), ngc(iq)
          read (ifiqg) ngvecc(1:3,1:ngc(iq),iq)
        enddo
      endif
      close(ifiqg)

!! === mapping of qtt ===
!! nkey, kk1,kk2,kk3, iqkkk are to get iqindx.
!!  q --> call rangedq(matmul(ginv,q), qx) ---> n= (qx+0.5*epsd)/epsd
!!       --->  ik1,ik2,ik3= tabkk(kkk,iqk,nkey) ---> iqkkk(ik1,ik2,ik3)
      if(ifi==1) then
         nqtt => nqnump
         qtt  => qp
         nkey => nkeyp
      elseif(ifi==2) then
         nqtt => nqnumc
         qtt  => qc
         nkey => nkeyc
      endif
!! followings are the same as codes in readeigen.F
      allocate(ieord(nqtt))
      allocate(key(3,0:nqtt),qxx(3,nqtt))
      key(:,0)=0 !dummy
      key=-99999
c      print *,'ginv_=',ginv_
      do iq=1,nqtt
         call rangedq(matmul(ginv_,qtt(:,iq)), qxx(:,iq))
ccccccccccccc
c        ix=1
c        print *,' xxxx ix qxx=',ix,iq,qtt(ix,iq),matmul(ginv_,qtt(:,iq))
ccccccccccccc
      enddo
!! get key and nkey for each ix.
      do ix =1,3
ccccccccccccc
c        do i=1,nqtt
c        print *,' ix qxx=',ix,i,qtt(ix,i),qxx(ix,i)
c        enddo
ccccccccccccc
         call sortea(qxx(ix,:),ieord,nqtt,isig)
ccccccccccccc
c        do i=1,nqtt
c        print *,' ix qxx=',ix,i,qxx(ix,ieord(i))
c        enddo
ccccccccccccc
         ik=0
         do i=1,nqtt
            kkk=(qxx(ix,ieord(i))+0.5d0*epsd)/epsd  !kkk is digitized by 1/epsd
            if(i==1.or.key(ix,ik)<kkk) then
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
c      write(6,*)' ifi init_qqq nqtt=',ifi,nqtt
c      write(6,*)'kkk3=',kkk3
c      write(6,*)'nkey=',nkey
c      write(6,*)'kk1=',kk1
      do i=1,nqtt
         kkk3= (qxx(:,i)+0.5*epsd)/epsd !kkk is digitized by 1/epsd
         call tabkk(kkk3(1), kk1,nkey(1), ik1)
         call tabkk(kkk3(2), kk2,nkey(2), ik2)
         call tabkk(kkk3(3), kk3,nkey(3), ik3)
         iqkkk(ik1,ik2,ik3)=i
cccccccccccccccc
c      write(6,*)' ik1,ik2,ik3 i=',ik1,ik2,ik3,i
cccccccccccccccc
      enddo
      deallocate(qxx)
      end subroutine init_readqg
!! ---
      subroutine tabkk(kkin, kktable,n, nout)
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
      write(6,*) i1,i2,kkin
      write(6,*) kktable(i1),kktable(i2)
      call rx( 'takk: error')
      end subroutine

c$$$c--- release to save memory area.
c$$$      subroutine releaseqg_notusednow(key)
c$$$      implicit none
c$$$      character*(*) key
c$$$      integer(4):: ifi
c$$$      if    (key=='QGpsi') then
c$$$        ifi=1
c$$$        deallocate(qp,ngvecp)
c$$$      elseif(key=='QGcou') then
c$$$        ifi=2
c$$$        deallocate(qc,ngvecc)
c$$$      else
c$$$        stop "releaseqg: in readQGcou"
c$$$      endif
c$$$      init(ifi)=.false.
c$$$      end subroutine
!!---------------------------------------------------------

!> Find index as q=qq(:,iq) with modulo of premitive vector.
!! ginv is the inverse of plat (premitive translation vector).
!! Use kk1,kk2,kk3,nkey(1:3),iqkkk to get iqindx.
      subroutine iqindx2qg(q,ifi, iqindx,qu)
      implicit none
      integer, intent(in):: ifi
      integer, intent(out):: iqindx
      real(8), intent(in) :: q(3)
      real(8), intent(out) :: qu(3)

      integer:: i_out, iq,iqx ,kkk3(3),ik1,ik2,ik3
      real(8):: qx(3),qzz(3)
      logical::debug=.false.
      if(ifi==1) then
c         nqtt => nqnump
         qtt  => qp
         nkey => nkeyp
         iqkkk => iqkkkp
         kk1 =>kk1p
         kk2 =>kk2p
         kk3 =>kk3p
      elseif(ifi==2) then
c         nqtt => nqnumc
         qtt  => qc
         nkey => nkeyc
         iqkkk => iqkkkc
         kk1 =>kk1c
         kk2 =>kk2c
         kk3 =>kk3c
      endif
c      if(debug) write(*,"(' iqindx2_: q=',3f20.15)") q
c      write(*,"(' iqindx2_: q=',3f20.15)") q
      call rangedq(matmul(ginv_,q), qzz)
      kkk3 = (qzz+0.5*epsd)/epsd
c      write(6,*)'kkk3=',kkk3
c      write(6,*)'kk1,nkey1',kk1,nkey(1)
c      write(6,*)'kk2,nkey2',kk2,nkey(2)
c      write(6,*)'kk3,nkey3',kk3,nkey(3)
      call tabkk(kkk3(1), kk1,nkey(1), ik1)
      call tabkk(kkk3(2), kk2,nkey(2), ik2)
      call tabkk(kkk3(3), kk3,nkey(3), ik3)
c      write(6,*)' ik1ik2ik3=',ik1,ik2,ik3
      iqindx = iqkkk(ik1,ik2,ik3)
c      write(6,*)'iqindx=',iqindx
      qu = qtt(:,iqindx)
c      write(6,*)'iqindx=',iqindx
c      write(6,*)'qu=',qu
      end subroutine

!> mini-sort routine.
      subroutine sortea(ea,ieaord,n,isig)
      real(8), intent(in) :: ea(n)
      integer(4), intent(inout) :: ieaord(n)
      integer, intent(in) :: n
      integer, intent(out) :: isig
      integer :: ix,i
      isig = 1
      do i = 1,n
        ieaord(i) = i
      enddo
      do ix= 2,n
        do i=ix,2,-1
          if( ea(ieaord(i-1)) >ea(ieaord(i) ) ) then
            call iswap (ieaord(i-1),ieaord(i))
            isig= -isig
            cycle
          endif
          exit
        enddo
      enddo
      end subroutine
      subroutine iswap(i,j)
      implicit none
      integer,intent(inout) :: i, j
      integer:: iwork
      iwork = j
      j = i
      i = iwork
      end subroutine
      end module m_readQG
