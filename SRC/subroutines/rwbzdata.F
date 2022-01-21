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
      integer,protected,public:: nq0i,nq0itrue,nq0iadd,iq0pin,nq0ix,neps       ! Number of Q0P
      real(8),allocatable,protected,public:: q0i(:,:),wt(:)  ! Q0P and its weight.
      integer,allocatable,public:: epslgroup(:)              !EPSwklm
      integer,public:: lxklm                  !EPSwklm
      real(8),allocatable,public:: epinv(:,:,:),wklm(:), dmlx(:,:), epinvq0i(:,:) !EPSwklm
      integer,protected,allocatable,public:: ixyz(:) ! ixyz(1:nq0i+nq0iadd) q0i for x,y,z directions
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
      read(ifbz) iq0pin,nq0i,iq0pin,nq0iadd
      allocate(  wt(1:nq0i+nq0iadd),q0i(1:3,1:nq0i+nq0iadd))
      read(ifbz) wt(1:nq0i+nq0iadd),q0i(1:3,1:nq0i+nq0iadd)
      if(iq0pin==1) then
         allocate(  ixyz(1:nq0i+nq0iadd))
         read(ifbz) ixyz(1:nq0i+nq0iadd)
         allocate(  dmlx(1:nq0i,1:9), epinv(1:3,1:3,1:nq0i), epinvq0i(1:nq0i,1:nq0i))
         read(ifbz) lxklm,dmlx(1:nq0i,1:9), epinv(1:3,1:3,1:nq0i), epinvq0i(1:nq0i,1:nq0i)
         allocate(  wklm(1:(lxklm+1)**2))
         read(ifbz) wklm(1:(lxklm+1)**2)
      elseif(iq0pin==2) then
         allocate(  epslgroup(1:nq0i))
         read(ifbz) epslgroup(1:nq0i)
      endif
      close(ifbz)
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
!      
      nq0ix = nq0i
      do i=1,nq0i
        if(wt(i)==0d0 ) then
           nq0ix = i-1
            exit
         endif
      enddo
      neps=nq0i-nq0ix           ! number of zero weight q0p which are used for ixc=2 or 3 mode.
      
      done_read_bzdata=.true.
      write(6,*)' end of read_BZdata'
      end subroutine
      end module
      
      subroutine checkrange(intq,n1,n2)
      integer:: intq,n1,n2
      if(intq<n1 .or. intq>n2) then
        print *,'checkrange: intq n1 n2= ',intq,n1,n2
        call rx( 'checkrange: stop ')
      endif
      end


