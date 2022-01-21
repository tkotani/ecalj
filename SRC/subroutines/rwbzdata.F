!! read BZDATA
      module m_read_bzdata
      implicit none

      public:: Read_bzdata
!! We set following data when you call read_BZDATA() 
      integer,protected,public :: n1,n2,n3,ngrp,nqbz,nqibz,nqbzw,nteti,ntetf ,itet
      integer,allocatable,protected,public :: 
     &    idtetf(:,:),ib1bz(:),idteti(:,:),
     &    nstar(:),irk(:,:),nstbz(:) !,index_qbz(:,:,:)
      real(8),allocatable,protected,public:: dq_(:) !qbas(:,:),qbasinv(:,:)
     &    , qbz(:,:),wbz(:),qibz(:,:) 
     &    , wibz(:),qbzw(:,:)
      real(8),protected,public :: qlat(3,3),ginv(3,3)
      logical,protected,public:: done_read_bzdata=.false.
      private
      
!! =====================================================
      contains
      subroutine read_BZDATA(hx0)
      intent(in)::           hx0
!! After you call this, you can access Brillowin Zone datas above ----
      integer :: intq(3),iqbz,ifbz,n,verbose,i 
      real(8) :: qout(3),deltaq(3)
      logical,optional:: hx0
      logical:: qbzreg
      write(6,*)' ### readin BZDATA ###'
c     allocate(,dq_(3),qbasinv(3,3)) !,qbasmc(3,3)
      open(newunit=ifbz, file='BZDATA',form='unformatted')
      read(ifbz) nqbz,nqibz, nqbzw, ntetf, nteti, ngrp, n1,n2,n3,qlat,ginv
      allocate( qibz(1:3,1:nqibz),wibz(1:nqibz),nstar(1:nqibz),irk(1:nqibz,1:ngrp))
      read(ifbz)qibz(1:3,1:nqibz),wibz(1:nqibz),nstar(1:nqibz),irk(1:nqibz,1:ngrp)
      allocate( qbz(1:3,1:nqbz),wbz(1:nqbz),nstbz(1:nqbz))
      read(ifbz)qbz(1:3,1:nqbz),wbz(1:nqbz),nstbz(1:nqbz)
      allocate( idtetf(0:3,1:ntetf))
      read(ifbz)idtetf(0:3,1:ntetf)
      allocate( ib1bz(1:nqbzw), qbzw(1:3,1:nqbzw))
      read(ifbz)ib1bz(1:nqbzw), qbzw(1:3,1:nqbzw)
      allocate( idteti(0:4,1:nteti))
      read(ifbz)idteti(0:4,1:nteti)
      read(ifbz) dq_
      close(ifbz)
c      allocate(n1,n2,n3,ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,nqibz_r)
c      call rbzdata(ifbz, !1,
c     &    ngrp,qbas,ginv, !qbasmc,
c     i    qbz, wbz,nstbz,           nqbz, 
c     i    qibz,wibz, nstar,irk,  nqibz,
c     i    idtetf, ntetf, qbzw,ib1bz, nqbzw,
c     i    idteti, nteti,dq_) !, qibz_r, nqibz_r )
c      write(6,"(a,9f9.4)")'read_BZDATA:ginv=',ginv
c      close(ifbz)
c      call minv33(qbas,qlatinv)
      if(present(hx0) .and. (.not.qbzreg())) then ! set off-gamma mesh
         deltaq = qlat(:,1)/n1 + qlat(:,2)/n2 +qlat(:,3)/n3
         do i=1,nqbz
            qbz(:,i) = qbz(:,i) - deltaq/2d0
            write(6,"('i qbz=',i3,3f8.4)") i,qbz(:,i)
         enddo
      endif
      if(abs(sum(wibz(1:nqibz))-2d0)>1d-10) then
        print *, 'sum (wibz)=', sum(wibz(1:nqibz))
        call Rx( 'read_BZDATA  sum (wibz) is not 2.')
      endif
      done_read_bzdata=.true.
      write(6,*)' end of read_BZdata'
      end subroutine
      end module

!! --------------------------------------------------------------------------
      subroutine rbzdata(ifbz, !job,
     x  ngrp,qlat,ginv, !qlatmc,
     x  qbz, wbz  ,nstbz, nqbz, 
     x  qibz,wibz, nstar,irk,  nqibz,
     x  idtetf, ntetf, qbzw,ib1bz, nqbzw,
     x  idteti, nteti,dq_) !, qibz_r,nqibz_r )
!! Read  BZ mesh data reuired for GW 
      implicit none
      integer(4):: nqbz,n1q,n2q,n3q,ntetf,nteti,nqbzw,iqbz,ifbz
     & ,nqibz,iqibz,itet,ngrp,job !,nqibz_r
      real(8)   :: plat(3,3),vol !,qlatmc(3,3)
      real(8):: qbz(3,nqbz),wbz(nqbz),qibz(3,nqibz),wibz(nqibz)
     &         ,qbzw(3,nqbzw),qlat(3,3),ginv(3,3),dq_(3) !, qibz_r(3,nqibz_r)
      integer,allocatable:: ipq(:),iw1(:)
      integer:: idtetf(0:3,ntetf),ib1bz(nqbzw),idteti(0:4,nteti)
     &       ,irk(nqibz,ngrp),nstar(nqibz),nstbz(nqbz),bzcase
      logical tetraf,tetrai
c      do iqibz = 1,nqibz
c        read(ifbz,"(4d24.16,i9)") 
c     &   qibz(1:3,iqibz),wibz(iqibz),nstar(iqibz)
c        read(ifbz,"(100i8)") irk(iqibz,1:ngrp)
c      enddo
c      do iqbz = 1,nqbz
c        read(ifbz, "(4d24.16,i10)") qbz(1:3,iqbz),wbz(iqbz),nstbz(iqbz)
c! bug fix (kino); May 2004 it was write(ifbz,"(4d24.16)") read(ifbz,"(4d24.16)")
c! But it was working in ifc
c!   (the strange number in ifbz file was correctly recoverd to the original values!)
c      enddo
c      if(ntetf>0) then
c        read(ifbz,"(4i10)") (idtetf(0:3,itet),itet=1,ntetf)
c        read(ifbz,"(i9,3d24.16)") 
c     &   (ib1bz(iqbz), qbzw(1:3,iqbz),iqbz=1,nqbzw)
c      endif
c      if(nteti>0) read(ifbz,"(5i10)") (idteti(0:4,itet),itet=1,nteti)
c      read (ifbz,"(3d24.16)") dq_
      end
!!
      subroutine checkrange(intq,n1,n2)
      integer:: intq,n1,n2
      if(intq<n1 .or. intq>n2) then
        print *,'checkrange: intq n1 n2= ',intq,n1,n2
        call rx( 'checkrange: stop ')
      endif
      end


