!! Jan2018 at KAIST. May2019 Ameslab, Nov2019tottori

!!---  qplist.dat reader -----------------------------------
      module m_readqplist     
      integer:: ndat
      real(8),allocatable,protected:: xdat(:),qplistsy(:,:)
      real(8),protected:: eferm !<-- temporary use. just as a memo.
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
c        write(6,'(" qplist ix xdata q=",i5,f9.4,x,3f9.4)')ix,xdat(ix),qplistsy(:,ix)
      enddo
 1011 continue
      ndat=ix-1
      end subroutine
      end module m_readqplist

!! -- Read HamiltionanPMTinfo and HamiltonianPMT. Then convert HamPMT to HamRsMTO  ------
      module m_HamPMT
      real(8),allocatable,protected:: plat(:,:),pos(:,:),qplist(:,:),qlat(:,:)
      integer,allocatable,protected:: nlat(:,:,:,:),npair(:,:),
     &  ib_table(:),l_table(:),k_table(:),nqwgt(:,:,:)
      integer,protected:: nkk1,nkk2,nkk3,nbas,nkp,npairmx,ldim,jsp,lso,ndimMTO,nsp,nspx
      real(8),protected:: epsovl,alat
      complex(8),allocatable,protected:: ovlmr(:,:,:,:),hammr(:,:,:,:)
      contains

      subroutine ReadHamPMTInfo()
!!  read information for crystal strucre, k points, neighbor pairs.
      implicit none
      integer:: ififft,ifile_handle,i,lold,m
      character*4:: cccx
      ififft = ifile_handle() !return unused file handle
      open(ififft,file='HamiltonianPMTInfo',form='unformatted')
      allocate(plat(3,3),qlat(3,3)) !plat primitive vectors, qlat:primitive in reciprocal space
      read(ififft) plat,nkk1,nkk2,nkk3,nbas,qlat
      nkp = nkk1*nkk2*nkk3
      allocate(qplist(3,nkp))
      allocate(pos(3,nbas))
      read(ififft) pos,alat  !atomic positions, unit of the scale. 
      read(ififft) qplist    !qpoint list. all mesh points in the BZ mesh
      allocate(npair(nbas,nbas)) ! pair table of atoms corresponding to the mesh points
      read(ififft) npair,npairmx
      allocate( nlat(3,npairmx,nbas,nbas), nqwgt(npairmx,nbas,nbas) )
      read(ififft) nlat,nqwgt
      read(ififft) ldim,lso,nsp ! size of Hamiltonian: PMT part
      allocate(ib_table(ldim),l_table(ldim),k_table(ldim))
      read(ififft)ib_table,l_table,k_table
      close(ififft)
      write(6,"('MHAM: --- MTO part of Hamiltonian index (real-harmonics table is in job_pdos script) --- ')")
      write(6,'("MHAM: MTO block dim=",i5)') ldim
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
      end subroutine ReadHamPMTInfo

c$$$  !! delta fun check for FFT: k --> T --> k 
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

      subroutine HamPMTtoHamRsMTO() !Convert HamPMT(k mesh) to HamRsMTO(real space)
      implicit none
      integer:: ifihmto
      integer:: ifile_handle,ikpd,ikp,ib1,ib2,ifih,it,iq,nev,nmx,ifig=-999,i,j,ndimPMT,lold,m
      complex(8),allocatable:: hamm(:,:),ovlm(:,:),t_zv(:,:)
      logical:: lprint=.true.,savez=.false.,getz=.false.,skipdiagtest=.true.
      real(8),allocatable:: evl(:)
      complex(8):: img=(0d0,1d0),aaaa,phase
      real(8)::qp(3),pi=4d0*atan(1d0)
      ifih=ifile_handle()
      open(ifih,file='HamiltonianPMT',form='unformatted')
      write(6,*)'Reaing HamiltonianPMT...'
      ndimMTO=ldim
      if(lso==1) ldim=ldim*2 !L.S mode
      nspx=nsp
      if(lso==1) nspx=1
      allocate(ovlmr(1:ndimMTO,1:ndimMTO,npairmx,nspx), hammr(1:ndimMTO,1:ndimMTO,npairmx,nspx))
      write(6,"('ndimMTO ldim lso=',i6,4i3)") ndimMTO,ldim,lso
      hammr=0d0
      ovlmr=0d0
      iq=0
      do 
         read(ifih,end=2019) qp,ndimPMT,lso,epsovl,jsp
!     jsp=isp in the collinear case; jsp=1 in the noncollinear
         if(jsp==1) iq=iq+1
         write(6,"('=== Reading Ham for iq,spin,q=',2i4,3f9.5)") iq,jsp,qp
c         if(ndimPMT/=ldim.and.(lso==0.or.lso==2)) call rx('lmfham:   ndimMTO/=ldim')
c         if(ndimPMT/=2*ldim.and.lso==1)           call rx('lmfham: 2*ndimMTO/=ldim') ! L.S mode or not
         allocate(ovlm(1:ndimPMT,1:ndimPMT),hamm(1:ndimPMT,1:ndimPMT))
         read(ifih) ovlm(1:ndimPMT,1:ndimPMT)
         read(ifih) hamm(1:ndimPMT,1:ndimPMT)
!! H(k) ->  H(T) FourierTransformation to real space
!!       only MTO part ndimMTO (ndimPMT = ndimMTO + ndimAPW)
         do i=1,ndimMTO
            do j=1,ndimMTO
               ib1 = mod(ib_table(i),ldim)
               ib2 = mod(ib_table(j),ldim)
               do it =1,npair(ib1,ib2) ! hammr_ij (T)= \sum_k hamm(k) exp(ikT).   it is the index for T
                  phase = 1d0/dble(nkp)* exp(img*2d0*pi* sum(qp*matmul(plat,nlat(:,it,ib1,ib2))))
                  hammr(i,j,it,jsp)= hammr(i,j,it,jsp)+ hamm(i,j)*phase
                  ovlmr(i,j,it,jsp)= ovlmr(i,j,it,jsp)+ ovlm(i,j)*phase
               enddo
            enddo
         enddo
!! skip diagonalization test or not ---> daigonalization test should reproduce original enery bands.
         if(skipdiagtest) then
            deallocate(ovlm,hamm)
            cycle
         else
!! Diagonalization test (H-eO)z=0
!!    These eigenvalues must generate the same eigenvalues as
!!    we perfermoed lmf-MPIK --writeham mode to generete HamiltonianPMT
            allocate(t_zv(ndimPMT,ndimPMT),evl(ndimPMT))
            nmx = ndimPMT
            if(epsovl<1.000001d-14) then !epsovl is the trancation to remove poor linear-dependency basis
               call zhev_tk2( ndimPMT , hamm(1:ndimPMT,1:ndimPMT) , ovlm(1:ndimPMT,1:ndimPMT) , nmx , nev,
     .              evl, t_zv, lprint,savez,getz,ifig) 
            else                    
               call zhev_tk3( ndimPMT , hamm(1:ndimPMT,1:ndimPMT) , ovlm(1:ndimPMT,1:ndimPMT) , nmx , nev,
     .              evl, t_zv, lprint,savez,getz,ifig,epsovl) 
            endif
            do i=1,nev
               if(jsp==1) write(6,"('eigenPMT_spin1 ',3i4,f15.5)") iq,jsp,i,evl(i)
               if(jsp==2) write(6,"('eigenPMT_spin2 ',3i4,f15.5)") iq,jsp,i,evl(i)
            enddo
            deallocate(ovlm,hamm,t_zv,evl)
         endif  
      enddo
 2019 continue
      write(6,*)'Read: total # of q for Ham=',iq
      close(ifih)
!! write RealSpace MTO Hamiltonian
      write(6,*)' Writing HamRsMTO...'
      ifihmto = ifile_handle()
      open(ifihmto,file='HamRsMTO',form='unformatted')
      write(ifihmto) ndimMTO,npairmx,nspx
      write(ifihmto) hammr(1:ndimMTO,1:ndimMTO,1:npairmx,1:nspx),
     &               ovlmr(1:ndimMTO,1:ndimMTO,1:npairmx,1:nspx)
      flush(ifihmto)
      close(ifihmto)
c      deallocate(hammr,ovlmr)
      print*,' OK: Wrote HamRsMTO file!'
      end subroutine HamPMTtoHamRsMTO
      end module m_HamPMT

!!-------------------------------------------------------------
      module m_HamRsMTO
!! read real-space MTO Hamiltonian
      integer,private:: ndimMTO,npairmx,nspx
      complex(8),allocatable,protected:: ovlmr(:,:,:,:),hammr(:,:,:,:)
      contains
!! read RealSpace MTO Hamiltonian
      subroutine ReadHamRsMTO()
      integer:: ifihmto,ifile_handle
      ifihmto=ifile_handle()
      open(ifihmto,file='HamRsMTO',form='unformatted')
      read(ifihmto) ndimMTO,npairmx,nspx
      write(6,*)'ndimMTO,npairmx,nspx=',ndimMTO,npairmx,nspx
      allocate(ovlmr(1:ndimMTO,1:ndimMTO,npairmx,nspx), hammr(1:ndimMTO,1:ndimMTO,npairmx,nspx))
      read(ifihmto) hammr(1:ndimMTO,1:ndimMTO,1:npairmx,1:nspx),
     &              ovlmr(1:ndimMTO,1:ndimMTO,1:npairmx,1:nspx)
      close(ifihmto)
      print*,' OK: Read HamRsMTO file!'
      end subroutine ReadHamRsMTO
      end module m_HamRsMTO



!! -----------------------------------------------------------------------------------      
!! -----------------------------------------------------------------------------------      
      program lmfham1
!! Read HamiltonianPMT and generates MTO-only Hamiltonian
      use m_HamPMT, only: 
     &  plat, npair,nlat,nqwgt, ldim, ndimMTO, nkp,qplist, epsovl, ib_table,alat,npairmx,nspx, 
     &  ReadHamPMTInfo, HamPMTtoHamRsMTO
! note.  HamPMTtoHamRsMTO do not change variables. Only generate HamRsMTO file.
      use m_HamRsMTO,  only: hammr,ovlmr,  ReadHamRsMTO
      use m_readqplist,only: eferm,qplistsy,ndat,xdat, Readqplistsy
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

!! Read infomation for Hamiltonian (lattice structures and index of basis).
      call ReadHamPMTInfo()
!! Read Hamiltonian files on k mesh point, and Get H(T) (real space Hamiltonian).
!! HamRsMTO is generated.
      call HamPMTtoHamRsMTO() 

!! HamRsMTO is readin
!! If HamRSMTO exist, you can skip 'call HamPMTtoHamRsMTO()'.
      call ReadHamRsMTO()

!! From here on, we get Hamitonian for given k points and get eigenvalues by diagonalizaztion.
!!      
!! If symlcase=T, we read qplist.dat (q point list, see bndfp.F), and use q points in it to compare
!! bands by original ecalj (by job_band), and TB hamiltonian read by ReadHamiltonianPMTInfo.
      if(symlcase) then
         call readqplistsy()
         ifsy1 = ifile_handle()
         open(ifsy1,file='band_lmfham_spin1.dat')
         if(nspx==2) ifsy2 = ifile_handle()
         if(nspx==2) open(ifsy2,file='band_lmfham_spin2.dat')
         write(6,*)  'ndat =',ndat
      endif

      write(6,*)  'ndimMTO=',ndimMTO
!! Get H(k) from H(T), and diagonalize
      allocate(ovlm(1:ndimMTO,1:ndimMTO),hamm(1:ndimMTO,1:ndimMTO))
      allocate(t_zv(ndimMTO,ndimMTO),evl(ndimMTO))
      nmx = ndimMTO
      ndatx = nkp
      if(symlcase) ndatx=ndat
      do ikp=1,ndatx
        if(symlcase) then
           qp= qplistsy(:,ikp)
        else
           qp = qplist(:,ikp) 
        endif
        write(6,"(' ikp along qplist, q=',i5,3f9.4)")ikp,qp
!     ! true q(i)= 2pi/alat * qplist(i,ikp)
      do jsp=1,nspx !nsp is the number of spin.  When lso=1(Lz.Sz), nspx=1
        ovlm = 0d0
        hamm = 0d0
        do i=1,ndimMTO
        do j=1,ndimMTO
          ib1 = mod(ib_table(i),ldim) !atomic-site index in the primitive cell
          ib2 = mod(ib_table(j),ldim)
          do it =1,npair(ib1,ib2)
            phase = 1d0/nqwgt(it,ib1,ib2)*exp(-img*2d0*pi* sum(qp*matmul(plat,nlat(:,it,ib1,ib2))))
            hamm(i,j)= hamm(i,j)+ hammr(i,j,it,jsp)*phase
            ovlm(i,j)= ovlm(i,j)+ ovlmr(i,j,it,jsp)*phase
          enddo
        enddo
        enddo
!!
!! Diagonalization. savez,getz,ifig are dummy.
!! Both zhev_tk2,3 works in the same manner except epsovl mechanism (linear dependency trancation of basis).
!! evl are eigenvalues, t_zv are eigenvuncitons. nmx is how many eigenvalues functions you get.
!! Usually obtaied # is nev, which is the same as nmx.        
        if(epsovl<1.000001d-14) then
          call zhev_tk2( ndimMTO , hamm , ovlm , nmx , nev,
     .     evl, t_zv, lprint,savez,getz,ifig) 
        else  
          call zhev_tk3( ndimMTO , hamm , ovlm , nmx , nev,
     .     evl, t_zv, lprint,savez,getz,ifig,epsovl) 
        endif

c        do i=1,nev
c          if(jsp==1) write(6,"('Eigen of MTO ham spin1 ',3i4,f15.5)") ikp,jsp,i,evl(i)
c          if(jsp==2) write(6,"('Eigen of MTO ham spin2 ',3i4,f15.5)") ikp,jsp,i,evl(i)
c        enddo
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
      if(nspx==2) close(ifsy2)
      write(6,"(a)")'OK! band_lmfham_spin*.dat has generated!'
      write(6,"(a)")'NOTE: We need to implement Fermi energy for MTO Hamiltonian.'
      write(6,"(a)")'  For a while, you may use Fermi energy when you plot bands (shown at L1:qplist.dat)'
      write(6,"(a)")'  See README_MATERIALS.org for how to make a plot for band_lmfham_spin*.dat'
      write(6,"(a)")'    For example, gnuplot scrpt can be'
      write(6,"(a)")'    ef=0.2239816400'
      write(6,"(a)")'    plot \'
      write(6,"(a)")'    "bnd001.spin1" u ($2):($3) lt 1 pt 1 w lp,\'
      write(6,"(a)")'    "bnd002.spin1" u ($2):($3) lt 1 pt 1 w lp,\'
      write(6,"(a)")'    "bnd003.spin1" u ($2):($3) lt 1 pt 1 w lp,\'
      write(6,"(a)")'    "bnd004.spin1" u ($2):($3) lt 1 pt 1 w lp,\'
      write(6,"(a)")'    "bnd005.spin1" u ($2):($3) lt 1 pt 1 w lp,\'
      write(6,"(a)")'    "bnd006.spin1" u ($2):($3) lt 1 pt 1 w lp,\'
      write(6,"(a)")'    "band_lmfham_spin1.dat" u ($1):(13.605*($2-ef)) pt 2'
      write(6,"(a)")'    pause -1'
      deallocate(ovlm,hamm)
      end program lmfham1
