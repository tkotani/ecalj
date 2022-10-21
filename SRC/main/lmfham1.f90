!! -- Read HamiltionanPMTinfo and HamiltonianPMT. Then convert HamPMT to HamRsMTO  ------
module m_HamPMT
  real(8),allocatable,protected:: plat(:,:),pos(:,:),qplist(:,:),qlat(:,:)
  integer,allocatable,protected:: nlat(:,:,:,:),npair(:,:),&
       ib_table(:),l_table(:),k_table(:),nqwgt(:,:,:)
  integer,protected:: nkk1,nkk2,nkk3,nbas,nkp,npairmx,ldim,jsp,lso,ndimMTO,nsp,nspx
  real(8),protected:: epsovl,alat
  complex(8),allocatable,protected:: ovlmr(:,:,:,:),hammr(:,:,:,:)
contains
  subroutine ReadHamPMTInfo() ! read information for crystal strucre, k points, neighbor pairs.
    implicit none
    integer:: ififft,ifile_handle,i,lold,m
    character*4:: cccx
    !ififft = ifile_handle() !return unused file handle
    open(newunit=ififft,file='HamiltonianPMTInfo',form='unformatted')
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
       write(6,"('MHAM: i ib(atom) l m k(EH,EH2,PZ)=',5i3)")i,ib_table(i),l_table(i),m,k_table(i)
    enddo
  end subroutine ReadHamPMTInfo

  !c$$$  !! delta fun check for FFT: k --> T --> k 
  !c$$$!!    \delta_{kk'} = \sum_{T \in T(i,j)} W_T exp( i (k-k') T)
  !c$$$      ikpd=7
  !c$$$      write(6,*)'test for ikpd=',ikpd
  !c$$$      do ikp=1,nkp
  !c$$$        qp = qplist(:,ikp) - qplist(:,ikpd)
  !c$$$        do ib1=1,nbas
  !c$$$          do ib2=1,nbas
  !c$$$            aaaa=0d0
  !c$$$            do it = 1,npair(ib1,ib2) 
  !c$$$              aaaa =  aaaa + 1d0/(nkp*nqwgt(it,ib1,ib2))*exp(img*2d0*pi* sum(qp*matmul(plat,nlat(:,it,ib1,ib2))))
  !c$$$            enddo
  !c$$$            cccx=''
  !c$$$            if(ikp==ikpd) cccx=' <--'
  !c$$$            write(6,"('\delta-fun test',i4,3f10.4,2i3,2f23.15,a)") ikp, qplist(:,ikp),ib1,ib2,aaaa,cccx
  !c$$$          enddo
  !c$$$        enddo
  !c$$$      enddo

  subroutine HamPMTtoHamRsMTO() !Convert HamPMT(k mesh) to HamRsMTO(real space)
    use m_ftox
    implicit none
    integer:: ifihmto
    integer:: ifile_handle,ikpd,ikp,ib1,ib2,ifih,it,iq,nev,nmx,ifig=-999,i,j,ndimPMT,lold,m
    complex(8),allocatable:: hamm(:,:),ovlm(:,:)!,t_zv(:,:)!,ovlmx(:,:)
    logical:: lprint=.true.,savez=.false.,getz=.false.,skipdiagtest=.false.
    real(8),allocatable:: evl(:)
    complex(8):: img=(0d0,1d0),aaaa,phase
    real(8)::qp(3),pi=4d0*atan(1d0)
    !ifih=ifile_handle()
    open(newunit=ifih,file='HamiltonianPMT',form='unformatted')
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
    qploop: do 
       read(ifih,end=2019) qp,ndimPMT,lso,epsovl,jsp ! jsp=isp in the collinear case; jsp=1 in the noncollinear
       if(jsp==1) iq=iq+1
       write(6,"('=== Reading Ham for iq,spin,q=',2i4,3f9.5)") iq,jsp,qp
       !c if(ndimPMT/=ldim.and.(lso==0.or.lso==2)) call rx('lmfham:   ndimMTO/=ldim')
       !c if(ndimPMT/=2*ldim.and.lso==1)        call rx('lmfham: 2*ndimMTO/=ldim') ! L.S mode or not
       allocate(ovlm(1:ndimPMT,1:ndimPMT),hamm(1:ndimPMT,1:ndimPMT))
       read(ifih) ovlm(1:ndimPMT,1:ndimPMT)
       read(ifih) hamm(1:ndimPMT,1:ndimPMT)
!       if(.not.skipdiagtest) then
       !! Diagonalization test (H-eO)z=0
       !!    These eigenvalues must generate the same eigenvalues as
       !!    we perfermoed lmf-MPIK --writeham mode to generete HamiltonianPMT
       epsovl=0d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
       hammwgt:block
         integer:: nx
         real(8)::beta,emu,val,wgt(ndimPMT),evlmto(ndimMTO),evl(ndimPMT),evlx(ndimMTO)
         complex(8):: oz(ndimPMT,ndimPMT),wnj(ndimPMT,ndimMTO),wnm(ndimPMT,ndimMTO),wnn(ndimMTO,ndimMTO),wnmbk(ndimPMT,ndimMTO)
         complex(8):: evecmto(ndimMTO,ndimMTO),evecpmt(ndimPMT,ndimPMT)
         complex(8):: ovlmx(ndimPMT,ndimPMT),hammx(ndimPMT,ndimPMT),fac(ndimPMT,ndimMTO),ddd(ndimMTO,ndimMTO)
         ovlmx= ovlm
         hammx= hamm
         nmx = ndimMTO
         call zhev_tk4(ndimMTO,hamm(1:ndimMTO,1:ndimMTO),ovlm(1:ndimMTO,1:ndimMTO), &
              nmx,nev, evlmto, evecmto, epsovl)
         ovlm= ovlmx
         hamm= hammx
         nmx = ndimPMT
         call zhev_tk4(ndimPMT,hamm(1:ndimPMT,1:ndimPMT),ovlm(1:ndimPMT,1:ndimPMT), &
              nmx,nev, evl,    evecpmt, epsovl)
!         do i=1,6! nev
!            if(jsp==1) write(6,"('eigenPMT_spin1 pmt',3i4,f15.5)") iq,jsp,i,evl(i)
!            if(jsp==1) write(6,"('eigenPMT_spin1 mto',3i4,f15.5)") iq,jsp,i,evlmto(i)
!            if(jsp==2) write(6,"('eigenPMT_spin2 pmt',3i4,f15.5)") iq,jsp,i,evl(i)
!            if(jsp==2) write(6,"('eigenPMT_spin2 mto',3i4,f15.5)") iq,jsp,i,evlmto(i)
!         enddo
         !write(6,*)" ndimPMT=",ndimPMT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         beta= 1d0 ! 1/Ry
         !val= ndimMTO !*nspx/nspx !test for si
         !call detemu(beta,val,evl,ndimPMT, emu) !determine emu for nocc=val
         !wgt = 1d0 /(exp(beta*(evl(:)-emu))+1d0)
         !oz = matmul(ovlmx,t_zv)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
         ovlm=ovlmx
         wnm=0d0
         do j=1,ndimMTO
            do i=1,ndimPMT
               fac(i,j)= sum(dconjg(evecpmt(:,i))*matmul(ovlm(:,1:ndimMTO),evecmto(1:ndimMTO,j)))
               wnm(i,j)= fac(i,j) *abs(fac(i,j)) 
            enddo
         enddo
         call GramSchmidt(ndimPMT,ndimMTO,wnm)
! unitatiry check         
         ddd=matmul(dconjg(transpose(wnm)),wnm)
         do i=1,ndimMTO
            do j=1,ndimMTO
               if(i==j.and. abs(ddd(i,j)-1d0)>1d-7) then
                  write(6,ftox)'eeeeeeee i i',i,i,ddd(i,j)
                  call rx('xxxxxxxxxxxx')
               elseif(i/=j.and. abs(ddd(i,j))>1d-7) then
                  write(6,ftox)'eeeeeeee i j',i,j,ddd(i,j)
                  call rx('xxxxxxxxxxxx')
               endif
            enddo
         enddo
         ! Mapping operator
         wnj = matmul(wnm,matmul(transpose(dconjg(evecmto(:,:))),ovlm(1:ndimMTO,1:ndimMTO)))
         nx=ndimPMT          
         do i=1,ndimMTO
            do j=1,ndimMTO
               hamm(i,j)= sum( dconjg(wnj(1:nx,i))*evl(1:nx)*wnj(1:nx,j)) !+ 0.02*hammx(i,j)
               ovlm(i,j)= sum( dconjg(wnj(1:nx,i))*wnj(1:nx,j) )
            enddo
         enddo
         !! Hamiltonian modified. Psi_MLO are eigenfunctions
         !wnn = matmul(transpose(dconjg(evecmto(:,:))),ovlm(1:ndimMTO,1:ndimMTO))
         !nx=ndimMTO
         !do i=1,ndimMTO
         !   evlx(i) = evl(i) !sum(dconjg(wnm(1:ndimPMT,i))*evl(1:ndimPMT)*wnm(1:ndimPMT,i))
         !enddo   
         !do i=1,ndimMTO
         !   do j=1,ndimMTO
         !      hamm(i,j)= sum( dconjg(wnn(1:nx,i))*evlx(1:nx)*wnn(1:nx,j)) 
         !      ovlm(i,j)= sum( dconjg(wnn(1:nx,i))*wnn(1:nx,j) )
         !   enddo
         !enddo
       !! Real space Hamiltonian. H(k) ->  H(T) FourierTransformation to real space
       !!       Only MTO part ndimMTO (ndimPMT = ndimMTO + ndimAPW)
       do i=1,ndimMTO
          do j=1,ndimMTO
             ib1 = mod(ib_table(i),ldim)
             ib2 = mod(ib_table(j),ldim)
             do it =1,npair(ib1,ib2)! hammr_ij (T)= \sum_k hamm(k) exp(ikT). it is the index for T
                phase = 1d0/dble(nkp)* exp(img*2d0*pi* sum(qp*matmul(plat,nlat(:,it,ib1,ib2))))
                hammr(i,j,it,jsp)= hammr(i,j,it,jsp)+ hamm(i,j)*phase
                ovlmr(i,j,it,jsp)= ovlmr(i,j,it,jsp)+ ovlm(i,j)*phase
!                hammr(i,j,it,jsp)= hammr(i,j,it,jsp)+ hammx(i,j)*phase !original MTO only 
!                ovlmr(i,j,it,jsp)= ovlmr(i,j,it,jsp)+ ovlmx(i,j)*phase
             enddo
          enddo
       enddo
       endblock hammwgt
       deallocate(ovlm,hamm)
       !! skip diagonalization test or not ---> daigonalization test should reproduce original enery bands.
    enddo qploop
2019 continue
    write(6,*)'Read: total # of q for Ham=',iq
    close(ifih)
    !! write RealSpace MTO Hamiltonian
    write(6,*)' Writing HamRsMTO...'
    !ifihmto = ifile_handle()
    open(newunit=ifihmto,file='HamRsMTO',form='unformatted')
    write(ifihmto) ndimMTO,npairmx,nspx
    write(ifihmto) hammr(1:ndimMTO,1:ndimMTO,1:npairmx,1:nspx),&
         ovlmr(1:ndimMTO,1:ndimMTO,1:npairmx,1:nspx)
    close(ifihmto)
    print*,' OK: Wrote HamRsMTO file!'
  end subroutine HamPMTtoHamRsMTO

  subroutine detemu(beta,val,evl,nev, emu)
    use m_ftox
    integer:: nev,ix
    real(8):: evl(nev),val,emu,beta,demu,valn
    emu=-1d0
    demu=1d0
    valn=-9999d0
    ix=0
!    write(6,ftox)evl
    do while(abs(valn-val)>1d-10)
       valn= sum(1d0/(exp(beta*(evl(:)-emu))+1d0))
       if(valn<val) emu=emu+demu
       if(valn>val) then
          emu = emu-demu
          demu= 0.5d0*demu
       endif
       ix=ix+1
!       write(6,*)'ix=',ix,emu
    enddo
!    stop 'vvvvvvvvvvvvxxxxxxxx'
  end subroutine detemu
    
  ! ssssssssssssssssssssssssssssssssss
  subroutine GramSchmidt(nv,n,zmel) 
    integer:: igb=1,it,itt,n,nv
    complex(8):: ov(n),vec(nv),dnorm2(nv),zmel(nv,n)
    real(8):: dnorm
    do it = 1,n        
       vec(:)= zmel(:,it)
       do itt = 1,it-1
          ov(itt) = sum( dconjg(zmel(:,itt))*vec(:))
       enddo
       vec = vec - matmul(zmel(:,1:it-1),ov(1:it-1))
       zmel(:,it) = vec/sum(dconjg(vec)*vec)**.5
    enddo
  end subroutine GramSchmidt
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
    !ifihmto=ifile_handle()
    open(newunit=ifihmto,file='HamRsMTO',form='unformatted')
    read(ifihmto) ndimMTO,npairmx,nspx
    write(6,*)'ndimMTO,npairmx,nspx=',ndimMTO,npairmx,nspx
    allocate(ovlmr(1:ndimMTO,1:ndimMTO,npairmx,nspx), hammr(1:ndimMTO,1:ndimMTO,npairmx,nspx))
    read(ifihmto) hammr(1:ndimMTO,1:ndimMTO,1:npairmx,1:nspx), ovlmr(1:ndimMTO,1:ndimMTO,1:npairmx,1:nspx)
    close(ifihmto)
    print*,' OK: Read HamRsMTO file!'
  end subroutine ReadHamRsMTO
end module m_HamRsMTO

!! -----------------------------------------------------------------------------------      
program lmfham1
  !! Read HamiltonianPMT and generates MTO-only Hamiltonian
  use m_HamPMT, only: plat, npair,nlat,nqwgt, ldim, ndimMTO, nkp,qplist, epsovl,&
       ib_table,alat,npairmx,nspx, ReadHamPMTInfo, HamPMTtoHamRsMTO
  ! note.  HamPMTtoHamRsMTO do not change variables. Only generate HamRsMTO file.
  use m_HamRsMTO,  only: hammr,ovlmr,  ReadHamRsMTO
  use m_readqplist,only: eferm,qplistsy,ndat,xdat, Readqplistsy
  use m_mpi,only: MPI__hx0fp0_rankdivider2Q, MPI__Qtask, &
       MPI__Initialize, MPI__Finalize,MPI__root, &
       MPI__Broadcast, MPI__rank,MPI__size, MPI__consoleout,MPI__barrier
  implicit none
  integer:: i,j,ikp,ib1,ib2,it,nmx,nev,jsp
  complex(8)::img=(0d0,1d0),phase
  complex(8),allocatable:: hamm(:,:),ovlm(:,:),t_zv(:,:)
  real(8),allocatable:: evl(:)
  real(8)::qp(3),pi=4d0*atan(1d0)
  logical:: lprint=.true.,savez=.false.,getz=.false. !dummy
  integer:: ifig=-999       !dummy
  integer:: ndatx,ifsy1,ifsy2,ifile_handle,ifsy
  logical:: symlcase=.true.
  call MPI__Initialize()
  call ReadHamPMTInfo()! Read infomation for Hamiltonian (lattice structures and index of basis).
                       
  call HamPMTtoHamRsMTO() ! HamRsMTO (real-space MTO based Hamiltonian) is generated.
                          ! HamPMT is for k mesh points, and Get H(T) (real space Hamiltonian).  
  call ReadHamRsMTO() !Read HamRsMTO. ! If HamRSMTO exist, you can skip 'call HamPMTtoHamRsMTO()'.

  GetEigenvaluesForSYML: block!Get Hamitonian at k points for hammr,ovlmr(realspace), then diagnalize.
  ! bands by original ecalj (by job_band), and TB hamiltonian read by ReadHamiltonianPMTInfo.
  if(symlcase) then ! When symlcase=T, read qplist.dat (q points list, see bndfp.F). 
     call readqplistsy()
     !ifsy1 = ifile_handle()
     open(newunit=ifsy1,file='band_lmfham_spin1.dat')
     !if(nspx==2) ifsy2 = ifile_handle()
     if(nspx==2) open(newunit=ifsy2,file='band_lmfham_spin2.dat')
     write(6,*)  'ndat =',ndat
  endif
  write(6,*)  'ndimMTO=',ndimMTO 
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
     write(6,"(' ikp along qplist, q=',i5,3f9.4)")ikp,qp! true q(i)= 2pi/alat * qplist(i,ikp)
     do jsp=1,nspx !nsp is the number of spin.  When lso=1(Lz.Sz), nspx=1
        ovlm = 0d0
        hamm = 0d0
        do i=1,ndimMTO
           do j=1,ndimMTO
              ib1 = mod(ib_table(i),ldim) !atomic-site index in the primitive cell
              ib2 = mod(ib_table(j),ldim)
              do it =1,npair(ib1,ib2)
                 phase=1d0/nqwgt(it,ib1,ib2)*exp(-img*2d0*pi* sum(qp*matmul(plat,nlat(:,it,ib1,ib2))))
                 hamm(i,j)= hamm(i,j)+ hammr(i,j,it,jsp)*phase
                 ovlm(i,j)= ovlm(i,j)+ ovlmr(i,j,it,jsp)*phase
              enddo
           enddo
        enddo
        call zhev_tk4(ndimMTO,hamm,ovlm,nmx,nev, evl,t_zv, epsovl)!Diangonale (hamm- evl ovlm) z=0
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
  write(6,"(a)")'!!! OK! band_lmfham_spin*.dat has generated! MTO only band plot !!!!!!!!!!!!'
  write(6,"(a)")'NOTE: We need to implement Fermi energy for MTO Hamiltonian.'
  write(6,"(a)")'For a while, you may use Fermi energy when you plot bands (shown at L1:qplist.dat)'
  write(6,"(a)")'See README_MATERIALS.org for how to make a plot for band_lmfham_spin*.dat'
  write(6,"(a)")'  For example, gnuplot scrpt can be'
  write(6,"(a)")'  ef=0.2239816400 (take it from efermi.lmf)'
  write(6,"(a)")'  plot \'
  write(6,"(a)")'  "bnd001.spin1" u ($2):($3) lt 1 pt 1 w lp,\'
  write(6,"(a)")'  "bnd002.spin1" u ($2):($3) lt 1 pt 1 w lp,\'
  write(6,"(a)")'  "bnd003.spin1" u ($2):($3) lt 1 pt 1 w lp,\'
  write(6,"(a)")'  "bnd004.spin1" u ($2):($3) lt 1 pt 1 w lp,\'
  write(6,"(a)")'  "bnd005.spin1" u ($2):($3) lt 1 pt 1 w lp,\'
  write(6,"(a)")'  "bnd006.spin1" u ($2):($3) lt 1 pt 1 w lp,\'
  write(6,"(a)")'  "band_lmfham_spin1.dat" u ($1):(13.605*($2-ef)) pt 2'
  write(6,"(a)")'  pause -1'
  deallocate(ovlm,hamm)
  end block GetEigenvaluesForSYML
endprogram lmfham1
