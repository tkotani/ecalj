c$$$      module m_readQ0P
c$$$      implicit none
c$$$      real(8),allocatable,protected:: wqt(:), wgt0(:,:),q0i(:,:) !,nx(:,:),nblocha(:)
c$$$      integer,protected:: nq0i,nq0iadd,nq0ix,neps
c$$$      integer,protected,allocatable:: ixyz(:)
c$$$
c$$$      contains
c$$$      subroutine readq0p()
c$$$      integer:: ifiq0p,ifile_handle,i,iq0pin
c$$$      logical:: debug=.false.
c$$$c      write(6,*) 'reading QOP'
c$$$c      ifiq0p=ifile_handle()
c$$$      open (newunit=ifiq0p,file='Q0P')
c$$$      read (ifiq0p,*) nq0i,iq0pin,nq0iadd
c$$$      allocate( wqt(1:nq0i),q0i(1:3,1:nq0i+nq0iadd),ixyz(nq0i+nq0iadd) )
c$$$      do i=1,nq0i+nq0iadd
c$$$         read (ifiq0p, * ) wqt(i),q0i(1:3,i),ixyz(i)
c$$$c         write (*, * ) wqt(i),q0i(1:3,i),ixyz(i)
c$$$      enddo
c$$$      nq0ix = nq0i
c$$$      do i=1,nq0i
c$$$         if(wqt(i)==0d0 ) then
c$$$            nq0ix = i-1
c$$$            exit
c$$$         endif
c$$$      enddo
c$$$      neps=nq0i-nq0ix ! number of zero weight q0p which are used for ixc=2 or 3 mode.
c$$$      write(6,*) ' num of zero weight q0p=',neps
c$$$      write(6,"(i3,f14.6,2x, 3f14.6)" )(i, wqt(i),q0i(1:3,i),i=1,nq0i+nq0iadd)
c$$$      close(ifiq0p)
c$$$      end subroutine
c$$$      end module
c$$$
