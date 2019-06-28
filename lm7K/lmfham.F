      module m_readqplist
      integer:: ndat
      real(8),allocatable:: xdat(:),qplistsy(:,:)
      real(8):: eferm
      contains
      subroutine readqplistsy()
      implicit none
      integer:: ifqplistsy,ifile_handle,nnn,ix
      ifqplistsy = ifile_handle()
      open(ifqplistsy,file='qplist.dat')
      nnn=1000
      allocate(xdat(nnn),qplistsy(3,nnn))
      read(ifqplistsy,*) eferm
      ix=0
      do
        ix=ix+1
        read(ifqplistsy,*,end=1011) xdat(ix),qplistsy(:,ix)
        print *,'xxxxxx ix=',ix
      enddo
 1011 continue
      ndat=ix-1
      end subroutine
      end module
      
! Started at Jan2018 at KAIST.
! Get Real space Hamiltonian for MTO-only case.
      module m_readHamiltonianMTO
      real(8),allocatable,protected:: plat(:,:),pos(:,:),qplist(:,:),qlat(:,:)
      integer,allocatable,protected:: nlat(:,:,:,:),npair(:,:),
     &  ib_table(:),l_table(:),k_table(:),nqwgt(:,:,:)
      integer,protected:: nkk1,nkk2,nkk3,nbas,nkp,npairmx,ldim,jsp,lso,ndimhx,nsp,nspx
      real(8),protected:: epsovl,alat
      complex(8),allocatable,protected:: ovlmr(:,:,:,:),hammr(:,:,:,:)

      contains
      subroutine readHamiltonianMTO()
      implicit none
      integer:: ififft,ifile_handle,ikpd,ikp,ib1,ib2,ifih,it,iq,nev,nmx,ifig=-999,i,j,ndimhx2,lold,m
      complex(8):: img=(0d0,1d0),aaaa,phase
      character*4:: cccx
      real(8)::qp(3),pi=4d0*atan(1d0)
      complex(8),allocatable:: hamm(:,:),ovlm(:,:),t_zv(:,:)
      logical:: lprint=.true.,savez=.false.,getz=.false.,skipdiagtest=.true.
      real(8),allocatable:: evl(:)

      ififft = ifile_handle() !return unused file handle
      open(ififft,file='HamiltonianMTOInfo',form='unformatted')
      allocate(plat(3,3),qlat(3,3))
      read(ififft) plat,nkk1,nkk2,nkk3,nbas,qlat
      nkp = nkk1*nkk2*nkk3
      allocate(qplist(3,nkp))
      allocate(pos(3,nbas))
      read(ififft) pos,alat   !atomic positions
      read(ififft) qplist !qpoint list. all in the BZ
      allocate(npair(nbas,nbas))
      read(ififft) npair,npairmx
      allocate( nlat(3,npairmx,nbas,nbas), nqwgt(npairmx,nbas,nbas) )
      read(ififft) nlat,nqwgt
      read(ififft) ldim,lso,nsp ! size of Hamiltonian: MTO part
      allocate(ib_table(ldim),l_table(ldim),k_table(ldim))
      read(ififft)ib_table,l_table,k_table
      close(ififft)
!!      
!     !
      write(6,"('MHAM: --- MTO Hamiltonian index (See real-harmonics table in job_pdos script) --- ')")
      write(6,'("MHAM: Ham dimension=",i5)') ldim
      lold=-999
      do i = 1,ldim
         if(l_table(i)/= lold) then !reset m of lm
            m=-l_table(i)
            lold=l_table(i)
         else
            m=m+1
         endif
         write(6,"('MHAM: i ib(atom) l m k(EH,EH2,PZ)=',5i3)")
     &        i, ib_table(i),l_table(i),m,k_table(i)
      enddo
      
c$$$!! delta fun check for FFT: k --> T --> k 
c$$$!!    \delta_{kk'} = \sum_{T \in T(i,j)} W_T exp( i (k-k') T)
c$$$      ikpd=7
c$$$      write(6,*)'test for ikpd=',ikpd
c$$$      do ikp=1,nkp
c$$$        qp = qplist(:,ikp) - qplist(:,ikpd)
c$$$        do ib1=1,nbas
c$$$          do ib2=1,nbas
c$$$            aaaa=0d0
c$$$            do it = 1,npair(ib1,ib2) 
c$$$              aaaa =  aaaa + 1d0/(nkp*nqwgt(it,ib1,ib2))*exp(img*2d0*pi* sum(qp*matmul(plat,nlat(:,it,ib1,ib2))))
c$$$            enddo
c$$$            cccx=''
c$$$            if(ikp==ikpd) cccx=' <--'
c$$$            write(6,"('\delta-fun test',i4,3f10.4,2i3,2f23.15,a)") ikp, qplist(:,ikp),ib1,ib2,aaaa,cccx
c$$$          enddo
c$$$        enddo
c$$$      enddo
!
      ifih=ifile_handle()
      open(ifih,file='HamiltonianMTO',form='unformatted')
      write(6,*)'Reaing HamiltonianMTO...'
      ndimhx=ldim
      if(lso==1) ldim=ldim*2 !L.S mode
      nspx=nsp
      if(lso==1) nspx=1
      allocate(ovlmr(1:ndimhx,1:ndimhx,npairmx,nspx), hammr(1:ndimhx,1:ndimhx,npairmx,nspx))
      print *,'ndimhx　ldim lso=',ndimhx,ldim,lso
      hammr=0d0
      ovlmr=0d0
      iq=0
      do 
         read(ifih,end=2019) qp,ndimhx2,lso,epsovl,jsp
!     jsp=isp in the collinear case; jsp=1 in the noncollinear
         if(jsp==1) iq=iq+1
         write(6,"('=== iq spin q=',2i4,3f9.5)") iq,jsp,qp
         if(ndimhx2/=ldim.and.(lso==0.or.lso==2)) call rx('lmfham:   ndimhx/=ldim')
         if(ndimhx2/=2*ldim.and.lso==1)           call rx('lmfham: 2*ndimhx/=ldim') ! L.S mode or not
         allocate(ovlm(1:ndimhx,1:ndimhx),hamm(1:ndimhx,1:ndimhx))
         read(ifih) ovlm(1:ndimhx,1:ndimhx)
         read(ifih) hamm(1:ndimhx,1:ndimhx)
!     ! H(k) ->  H(T) FourierTransformation
         do i=1,ndimhx
            do j=1,ndimhx
               ib1 = mod(ib_table(i),ldim)
               ib2 = mod(ib_table(j),ldim)
               do it =1,npair(ib1,ib2) ! hammr_ij (T)= \sum_k hamm(k) exp(ikT).   it is the index for T
                  phase = 1d0/dble(nkp)* exp(img*2d0*pi* sum(qp*matmul(plat,nlat(:,it,ib1,ib2))))
                  hammr(i,j,it,jsp)= hammr(i,j,it,jsp)+ hamm(i,j)*phase
                  ovlmr(i,j,it,jsp)= ovlmr(i,j,it,jsp)+ ovlm(i,j)*phase
               enddo
            enddo
         enddo
         
!! skip diagonalization test or not
         if(skipdiagtest) then
            deallocate(ovlm,hamm)
            cycle
         else
!!    diagonalization test (H-eO)z=0
            allocate(t_zv(ndimhx,ndimhx),evl(ndimhx))
            nmx = ndimhx
            if(epsovl<1.000001d-14) then
               call zhev_tk2( ndimhx , hamm , ovlm , nmx , nev,
     .              evl, t_zv, lprint,savez,getz,ifig) 
            else                    
               call zhev_tk3( ndimhx , hamm , ovlm , nmx , nev,
     .              evl, t_zv, lprint,savez,getz,ifig,epsovl) 
            endif
            do i=1,nev
               if(jsp==1) write(6,"('eigen111_1 ',3i4,f15.5)") iq,jsp,i,evl(i)
               if(jsp==2) write(6,"('eigen111_2 ',3i4,f15.5)") iq,jsp,i,evl(i)
            enddo
            deallocate(ovlm,hamm,t_zv,evl)
         endif  
      enddo
 2019 continue

      write(6,*)'read:total # of q =',iq
      close(ifih)
      end subroutine readHamiltonianMTO
      end module m_readHamiltonianMTO


!!-------------------------------------------
      program lmfham
!! Read Hamiltonian files for MTO only case.
      use m_readHamiltonianMTO, only: 
     &  plat, npair,nlat,nqwgt, ldim,ndimhx, hammr,ovlmr, nkp,qplist, epsovl, ib_table,nspx,alat,
     &     ReadHamiltonianMTO
      use m_readqplist,only:
     &  readqplistsy,eferm,qplistsy,ndat,xdat
      implicit none
      integer:: i,j,ikp,ib1,ib2,it,nmx,nev,jsp
      complex(8)::img=(0d0,1d0),phase
      complex(8),allocatable:: hamm(:,:),ovlm(:,:),t_zv(:,:)
      real(8),allocatable:: evl(:)
      real(8)::qp(3),pi=4d0*atan(1d0)
      logical:: lprint=.true.,savez=.false.,getz=.false. !dummy
      integer:: ifig=-999       !dummy
      integer:: ndatx,ifsy1,ifsy2,ifile_handle,ifsy

!! syml case. Read q points for test from qplist.dat file generated by lmf (band plot mode).
      logical:: symlcase=.true.

!! Show git id
      call show_programinfo(6)

!! Read Hamiltonian files and Get H(T) (real space Hamiltonian).
      call readHamiltonianMTO() 

!! If symlcase=T, we read qplist.dat (q point list, see bndfp.F), and use q points in it to compare
!! bands by original ecalj (by job_band), and TB hamiltonian read by readHamiltonianMTO.
      if(symlcase) then
         call readqplistsy()
         ifsy1 = ifile_handle()
         ifsy2 = ifile_handle()
         open(ifsy1,file='band_lmfham1.dat')
         open(ifsy2,file='band_lmfham2.dat')
         write(6,*)'ndat=',ndat
      endif
!     ! Get H(k) from H(T), and diagonalize
      allocate(ovlm(1:ndimhx,1:ndimhx),hamm(1:ndimhx,1:ndimhx))
      allocate(t_zv(ndimhx,ndimhx),evl(ndimhx))
      nmx = ndimhx
      ndatx = nkp
      if(symlcase) ndatx=ndat
      do ikp=1,ndatx
        if(symlcase) then
           qp= qplistsy(:,ikp)
        else
           qp = qplist(:,ikp) 
        endif
        write(6,"(' ikp, q=',i5,3f9.4)")ikp,qp
!     ! true q(i)= 2pi/alat * qplist(i,ikp)
      do jsp=1,nspx !nsp is the number of spin.  When lso=1(Lz.Sz), nspx=1
        ovlm = 0d0
        hamm = 0d0
        do i=1,ndimhx
        do j=1,ndimhx
          ib1 = mod(ib_table(i),ldim) !atomic-site index in the primitive cell
          ib2 = mod(ib_table(j),ldim)
          do it =1,npair(ib1,ib2)
            phase = 1d0/nqwgt(it,ib1,ib2)*exp(-img*2d0*pi* sum(qp*matmul(plat,nlat(:,it,ib1,ib2))))
            hamm(i,j)= hamm(i,j)+ hammr(i,j,it,jsp)*phase
            ovlm(i,j)= ovlm(i,j)+ ovlmr(i,j,it,jsp)*phase
          enddo
        enddo
        enddo

!! Diagonalization
        if(epsovl<1.000001d-14) then
          call zhev_tk2( ndimhx , hamm , ovlm , nmx , nev,
     .     evl, t_zv, lprint,savez,getz,ifig) !stock eigenfunctions z
        else  
          call zhev_tk3( ndimhx , hamm , ovlm , nmx , nev,
     .     evl, t_zv, lprint,savez,getz,ifig,epsovl) !stock eigenfunctions z
        endif

        do i=1,nev
          if(jsp==1) write(6,"('eigen222_1 ',3i4,f15.5)") ikp,jsp,i,evl(i)
          if(jsp==2) write(6,"('eigen222_2 ',3i4,f15.5)") ikp,jsp,i,evl(i)
        enddo
!!
        if(symlcase) then
        if(jsp==1) ifsy=ifsy1
        if(jsp==2) ifsy=ifsy2
        do i=1,nev
           write(ifsy,"(f15.5,f15.5,2i4)") xdat(ikp),evl(i),jsp,i
        enddo
        endif

      enddo    
      enddo
      close(ifsy1)
      close(ifsy2)
      deallocate(ovlm,hamm)
      end program lmfham
