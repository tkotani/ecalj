      module m_w0w0i
!! The algolizm is a little too compolcated.
!! See the Friedlich's paper. 
!! output
      complex(8),allocatable,protected :: w0(:),w0i(:),llmat(:,:)

      contains

!! all are input
      subroutine w0w0i(llw,llwI,nw_i,nw,nq0i,niw,q0i)
      use m_read_bzdata,only: 
     &  wbz
      real(8),intent(in):: q0i(1:3,1:nq0i) 
      complex(8),intent(in):: llw(nw_i:nw,nq0i), llwI(niw,nq0i)
      integer:: nw_i,nw,nq0i,niw,ifidmlx,ifile_handle
      real(8),allocatable:: wklm(:)
      real(8),allocatable:: dmlx(:,:),epinvq0i(:,:),epinv(:,:,:) !,epinvq0i_m1(:,:),qeibz(:,:,:)
      logical:: readw0w0itest
      complex(8):: llmat_dummy(3,3)
      write(6,*) 'w0w0i:'
      write(6,*)' ==== newaniso mode W(0) divergent part ==== '

!! == W(0) divergent part ==
!! getw0 routine...
c        ifidmlx = iopen('EPSwklm',0,0,0)
      ifidmlx=ifile_handle()
      open(ifidmlx,file='EPSwklm',form='unformatted')
      read(ifidmlx) nq0i,lxklm
      allocate( dmlx(nq0i,9))
      allocate( epinvq0i(nq0i,nq0i),epinv(3,3,nq0i))
      nlxklm=(lxklm+1)**2
      allocate( wklm(nlxklm))
      read(ifidmlx) dmlx, epinv,epinvq0i
      read(ifidmlx) wklm
      close(ifidmlx) ! = iclose('EPSwklm') !jan2013 iclose(ifidmlx)
!! starting from llw(iw,iq0),llwI(iw,iq0)
!! === <e|L|e> (eq.36 in Friedrich paper) is expanded in YL -->stored in llwyl. ===
      allocate(w0(nw_i:nw),w0i(niw),llmat(3,3))
c       write(6,*)' goto getw0 nq0i epinvq0i=',nq0i,epinvq0i
!! wbz(1) is the weight for q=0 = 1/(n1*n2*n3)
!! llmat is added. July2016. llw is calculated at iw=0, when nw_i<=0
      call getw0(llw, nw_i,nw,nq0i,dmlx,epinvq0i,wklm,wbz(1), lxklm,  q0i,epinv,w0,  llmat)
      call getw0(llwI,1,niw  ,nq0i,dmlx,epinvq0i,wklm,wbz(1), lxklm,  q0i,epinv,w0i, llmat_dummy)

!! test mode
      if(ixc/=1011) then
      ifw0w0i = ifile_handle()
      open(ifw0w0i,file='W0W0I',form='unformatted')
      write(ifw0w0i) nw_i,nw,niw,nq0i
      write(ifw0w0i) llw(nw_i:nw,1:nq0i)
      write(ifw0w0i) llwI(1:niw,1:nq0i)
      write(ifw0w0i) w0(nw_i:nw)
      write(ifw0w0i) w0i(1:niw)
      close(ifw0w0i)
      endif
      do i=nw_i,nw
        write(6,"('w0 =',i4,2f13.4)")i,w0(i)
      enddo   
      do i=1,niw
        write(6,"('w0i=',i4,2f13.4)")i,w0i(i)
      enddo   

      end subroutine w0w0i

      end module m_w0w0i
