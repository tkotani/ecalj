      module m_readQ0P
      implicit none
      real(8),allocatable,protected:: wqt(:), wgt0(:,:),q0i(:,:) !,nx(:,:),nblocha(:)
      integer,protected:: nq0i,nq0iadd,nq0ix,neps
      integer,protected,allocatable:: ixyz(:)

      contains
      subroutine readq0p()
      integer:: ifiq0p,ifile_handle,i,iq0pin
      logical:: debug=.false.
c      write(6,*) 'reading QOP'
c      ifiq0p=ifile_handle()
      open (newunit=ifiq0p,file='Q0P')
      read (ifiq0p,*) nq0i,iq0pin,nq0iadd
      allocate( wqt(1:nq0i),q0i(1:3,1:nq0i+nq0iadd),ixyz(nq0i+nq0iadd) )
      do i=1,nq0i+nq0iadd
         read (ifiq0p, * ) wqt(i),q0i(1:3,i),ixyz(i)
c         write (*, * ) wqt(i),q0i(1:3,i),ixyz(i)
      enddo
      nq0ix = nq0i
      do i=1,nq0i
         if(wqt(i)==0d0 ) then
            nq0ix = i-1
            exit
         endif
      enddo
      neps=nq0i-nq0ix ! number of zero weight q0p which are used for ixc=2 or 3 mode.
      write(6,*) ' num of zero weight q0p=',neps
      write(6,"(i3,f14.6,2x, 3f14.6)" )(i, wqt(i),q0i(1:3,i),i=1,nq0i+nq0iadd)
      close(ifiq0p)
      end subroutine
      end module

