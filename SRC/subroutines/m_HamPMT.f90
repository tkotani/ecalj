!> Read HamiltionanPMTinfo and HamiltonianPMT. Then convert HamPMT to HamRsMTO  
module m_HamPMT 
  use m_MPItk,only: procid, master_mpi, nsize,master
  use m_lgunit,only:stdo
  use m_ftox
  use m_lmfinit,only: oveps
  real(8),allocatable,protected:: plat(:,:),pos(:,:),qplist(:,:),qlat(:,:),symops(:,:,:)
  integer,allocatable,protected:: nlat(:,:,:,:),npair(:,:),ib_table(:),l_table(:),k_table(:),ispec_table(:),nqwgt(:,:,:)
  character(8),allocatable,protected:: slabl_table(:)
  integer,protected:: nkk1,nkk2,nkk3,nbas,nkp,npairmx,ldim,jsp,lso,nsp,nspx,ngrp !ldim is number of MTOs
  real(8),protected:: alat
  complex(8),allocatable,protected:: ovlmr(:,:,:,:),hammr(:,:,:,:)
  integer:: ndimMTO 
contains
  subroutine ReadHamPMTInfo() ! read information for crystal strucre, k points, neighbor pairs.
    implicit none
    integer:: ififft,i,lold,m,ibold,ioff
    character*4:: cccx
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
    read(ififft) ngrp
    allocate(symops(3,3,ngrp))
    read(ififft) ldim,lso,nsp,symops ! size of Hamiltonian: PMT part
    allocate(ib_table(ldim),l_table(ldim),k_table(ldim),ispec_table(ldim),slabl_table(ldim))
    read(ififft)ib_table,l_table,k_table,ispec_table,slabl_table
    close(ififft)
    if(master_mpi) write(stdo,"('MHAM: --- MTO part of PMT Hamiltonian index (real-harmonics table is in job_pdos script) --- ')")
    if(master_mpi) write(stdo,'("MHAM: MTO block dim=",i5)') ldim
    lold=-999
    ibold=-999
    ioff=0
    do i = 1,ldim
       if(l_table(i)/= lold) then !reset m of lm
          m=-l_table(i)
          lold=l_table(i)
       else
          m=m+1
       endif
       if(ib_table(i)/=ibold) then
          ioff=i-1
          ibold=ib_table(i)
       endif
       if(master_mpi) write(stdo,"('MHAM: i i-ioffib ib(atom) l m k(1:EH,2:EH2,3:PZ)=',i4,5i3)")&
            i,i-ioff,ib_table(i),l_table(i),m,k_table(i)
    enddo
  end subroutine ReadHamPMTInfo
  !c$$$  !! delta fun check for FFT: k --> T --> k 
  !c$$$!!    \delta_{kk'} = \sum_{T \in T(i,j)} W_T exp( i (k-k') T)
  !c$$$      ikpd=7
  !c$$$      write(stdo,*)'test for ikpd=',ikpd
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
  !c$$$            write(stdo,"('\delta-fun test',i4,3f10.4,2i3,2f23.15,a)") ikp, qplist(:,ikp),ib1,ib2,aaaa,cccx
  !c$$$          enddo
  !c$$$        enddo
  !c$$$      enddo
  subroutine HamPMTtoHamRsMPO(facw,ecutw,eww) !Convert HamPMT(k mesh) to HamRsMPO(real space)
    use m_zhev,only:zhev_tk4
    use m_readqplist,only: eferm
!    use m_mpi,only: MPI__reduceSum
    implicit none
    integer:: ifihmto
    integer::ikpd,ikp,ib1,ib2,ifih,it,iq,nev,nmx,ifig=-999,i,j,ndimPMT,lold,m
    complex(8),allocatable:: hamm(:,:),ovlm(:,:)!,t_zv(:,:)!,ovlmx(:,:)
    logical:: lprint=.true.,savez=.false.,getz=.false.,skipdiagtest=.false.
    complex(8):: img=(0d0,1d0),aaaa,phase
    real(8)::qp(3),pi=4d0*atan(1d0),fff,ef,fff1=2,fff2=2,fff3=0 ,facw,ecutw,eww,xxx
    integer:: nn,ib,k,l,ix5,imin,ixx,j2,j1,j3,nx,ix(ldim),iqini,iqend,ndiv
    integer:: ndimMTO !ndimMTO<ldim if we throw away f MTOs, for example. 
    integer:: ib_tableM(ldim),k_tableM(ldim),l_tableM(ldim),ierr
    nn=0
    do i=1,ldim  !only MTOs. Further restrictions.
       !if(l_table(i)>=3) cycle !only spd. skip f orbitals.     !if(k_table(i)==2.and.l_table(i)>=2) cycle ! throw away EH2 for d
       !if( (ib_table(i)==3.or.ib_table(i)==4).and.l_table(i)>=2) cycle !for Oxygen of NiO. skip 3d
       !   write(stdo,*) 'ham1 index', i,ib_table(i),l_table(i),k_table(i)
 !      if( k_table(i)>2) cycle !skip PZ orbitals !
       nn=nn+1
       ix(nn)=i
       ib_tableM(nn)=ib_table(i)
       k_tableM(nn) =k_table(i)
       l_tableM(nn) =l_table(i)
    enddo
    ndimMTO=nn
    open(newunit=ifih,file='HamiltonianPMT',form='unformatted')
    if(lso==1) ldim=ldim*2 !L.S mode
    nspx=nsp
    if(lso==1) nspx=1
    if(master_mpi) write(stdo,ftox)'Reading HamiltonianPMT: ndimMTO ldim lso=',ndimMTO,ldim,lso
    allocate(ovlmr(1:ndimMTO,1:ndimMTO,npairmx,nspx), hammr(1:ndimMTO,1:ndimMTO,npairmx,nspx))
    hammr=0d0
    ovlmr=0d0
    ndiv= nkp/nsize 
    if(nkp>ndiv*nsize) ndiv=ndiv+1    !MPI division
    iqini =         ndiv*procid+1     !initial for each procid
    iqend = min(nkp,ndiv*procid+ndiv) !end  for each procid
    write(stdo,ftox)'nsize procid iqini iqend=',nsize,procid,iqini,iqend
    iq=0
    qploop: do 
       read(ifih,end=2019) qp,ndimPMT,lso,xxx,jsp !jsp for isp; if so=1, jsp=1 only
       if(jsp==1) iq=iq+1
       if(iq<iqini.or.iqend<iq) cycle 
       if(master_mpi) write(stdo,"('=== Reading Ham for iq,spin,q=',2i4,3f9.5)") iq,jsp,qp
       allocate(ovlm(1:ndimPMT,1:ndimPMT),hamm(1:ndimPMT,1:ndimPMT))
       read(ifih) ovlm(1:ndimPMT,1:ndimPMT)
       read(ifih) hamm(1:ndimPMT,1:ndimPMT)
       !epsovl=0d0 !1d-8
       GETham_ndimMTO: block
         real(8):: evlmlo(ndimMTO)
         call Hreduction(.false.,facw,ecutw,eww,ndimPMT,hamm(1:ndimPMT,1:ndimPMT),ovlm(1:ndimPMT,1:ndimPMT), & !Get reduced Hamitonian for ndimMTO
              ndimMTO,ix,fff1,      hamm(1:ndimMTO,1:ndimMTO),ovlm(1:ndimMTO,1:ndimMTO))
         
         if(iq==3) then !            if(sum([qp(2),qp(3)]**2)<1d-3) then !check write
         Checkfinaleigen: block
           real(8):: rydberg
           complex(8):: evec(ndimMTO**2), hh(ndimMTO,ndimMTO),oo(ndimMTO,ndimMTO)
           hh = hamm(1:ndimMTO,1:ndimMTO) 
           oo = ovlm(1:ndimMTO,1:ndimMTO)
           nmx= ndimMTO
           !do i=1,ndimMTO
           !   write(stdo,ftox)i,ftof(abs(hh(1:ndimMTO,i)))
           !   write(stdo,ftox)i,ftof(abs(oo(1:ndimMTO,i)))
           !enddo   
           write(stdo,ftox) ' checkfinaleigen zhev_tk4',ftof(qp,3)
           call zhev_tk4(ndimMTO,hh,oo, nmx,nev,evlmlo, evec, oveps) !epsovl)
           write(stdo,ftox) '  evlmto=',nev,ftof(evlmlo(1:10)*rydberg())
           write(stdo,ftox) '  evlmto=',nev,ftof(evlmlo(11:20)*rydberg())
         endblock Checkfinaleigen
         endif
         !            endif
         
       endblock GETham_ndimMTO
       GETrealspaceHamiltonian: block ! H(k) ->  H(T) FourierTransformation to real space
         do i=1,ndimMTO
            do j=1,ndimMTO
               ib1 = ib_tableM(i) 
               ib2 = ib_tableM(j) 
               do it =1,npair(ib1,ib2)! hammr_ij (T)= \sum_k hamm(k) exp(ikT). it is the index for T
                  phase = 1d0/dble(nkp)* exp(img*2d0*pi* sum(qp*matmul(plat,nlat(:,it,ib1,ib2))))
                  hammr(i,j,it,jsp)= hammr(i,j,it,jsp)+ hamm(i,j)*phase
                  ovlmr(i,j,it,jsp)= ovlmr(i,j,it,jsp)+ ovlm(i,j)*phase
               enddo
            enddo
         enddo
       endblock GETrealspaceHamiltonian
       deallocate(ovlm,hamm)
    enddo qploop
2019 continue
    call mpibc2_complex(hammr,size(hammr),'m_HamPMT_hammr') !to master
    call mpibc2_complex(ovlmr,size(ovlmr),'m_HamPMT_ovlmr') !to master
    if(master_mpi) then
       write(stdo,*)'Read: total # of q for Ham=',iq
       close(ifih)
       !! write RealSpace MTO Hamiltonian
       !ix(1:ndimMTO)=ix1(1:ndimMTO) !for atom idex
       write(stdo,*)' Writing HamRsMTO... ndimMTO=',ndimMTO
       open(newunit=ifihmto,file='HamRsMTO',form='unformatted')
       write(ifihmto) ndimMTO,npairmx,nspx
       write(ifihmto) hammr(1:ndimMTO,1:ndimMTO,1:npairmx,1:nspx),ovlmr(1:ndimMTO,1:ndimMTO,1:npairmx,1:nspx) !,ix(1:ndimMTO)
       write(ifihmto) ib_tableM(1:ndimMTO),k_tableM(1:ndimMTO),l_tableM(1:ndimMTO)
       close(ifihmto)
    endif   
    if(master_mpi) write(stdo,*)" Wrote HamRsMTO file! End of lmfham1"
  end subroutine HamPMTtoHamRsMPO
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
       zmel(:,it) = vec/sum(dconjg(vec)*vec)**.5d0
    enddo
  end subroutine GramSchmidt
end module m_HamPMT
module m_HamRsMPO ! read real-space MPO Hamiltonian
  use m_lgunit,only:stdo
  use m_ftox
  integer,protected:: ndimMTO,npairmx,nspx  !ndimMTO<ldim if we throw away f MTOs, for example. 
  integer,allocatable,protected:: ib_tableM(:),l_tableM(:),k_tableM(:)
  complex(8),allocatable,protected:: ovlmr(:,:,:,:),hammr(:,:,:,:)
contains
  subroutine ReadHamRsMPO()! read RealSpace MTO Hamiltonian
    use m_MPItk,only: master_mpi
    integer:: ifihmto
    open(newunit=ifihmto,file='HamRsMTO',form='unformatted')
    read(ifihmto) ndimMTO,npairmx,nspx
!    allocate(ix(ndimMTO))
    if(master_mpi) write(stdo,ftox)'MTOHamiltonian: ndimMTO,npairmx,nspx=',ndimMTO,npairmx,nspx
    allocate(ovlmr(1:ndimMTO,1:ndimMTO,npairmx,nspx), hammr(1:ndimMTO,1:ndimMTO,npairmx,nspx))
    read(ifihmto)hammr(1:ndimMTO,1:ndimMTO,1:npairmx,1:nspx),ovlmr(1:ndimMTO,1:ndimMTO,1:npairmx,1:nspx) !,ix(1:ndimMTO)
    allocate(ib_tableM(1:ndimMTO),k_tableM(1:ndimMTO),l_tableM(1:ndimMTO))
    read(ifihmto) ib_tableM(1:ndimMTO),k_tableM(1:ndimMTO),l_tableM(1:ndimMTO)
    close(ifihmto)
    if(master_mpi) write(stdo,*)'OK: Read HamRsMTO file! Use i-ioffib for setting <Worb>'
  end subroutine ReadHamRsMPO
end module m_HamRsMPO
subroutine Hreduction(iprx,facw,ecutw,eww,ndimPMT,hamm,ovlm,ndimMTO,ix,fff1, hammout,ovlmout)!Reduce H(ndimPMT) to H(ndimMTO)
  use m_zhev,only:zhev_tk4
  use m_readqplist,only: eferm
  use m_HamPMT,only: GramSchmidt!,epsovl
  use m_lgunit,only:stdo
  use m_lmfinit,only:oveps
  implicit none
  integer::i,j,ndimPMT,ndimMTO,nx,nmx,ix(ndimMTO),nev,nxx,jj,ndimPMTx
  real(8)::beta,emu,val,wgt(ndimPMT),evlmto(ndimMTO),evl(ndimPMT),evlx(ndimPMT),facw,ecutw,eww
  complex(8):: evecmto(ndimMTO,ndimMTO),evecpmt(ndimPMT,ndimPMT)
  complex(8):: ovlmx(ndimPMT,ndimPMT),hammx(ndimPMT,ndimPMT),fac(ndimPMT,ndimMTO),ddd(ndimMTO,ndimMTO)
  complex(8):: hamm(ndimPMT,ndimPMT),ovlm(ndimPMT,ndimPMT)
  complex(8):: hammout(ndimMTO,ndimMTO),ovlmout(ndimMTO,ndimMTO)
  complex(8),allocatable :: wnj(:,:),wnm(:,:)
  real(8):: fff1,fff !epsovl=1d-8 epsovlm=0d0 ,
  logical:: iprx
  ovlmx= ovlm
  hammx= hamm
  nmx = ndimMTO !  write(stdo,*)'Start Hreduction: 111'
  call zhev_tk4(ndimMTO,hamm(ix(1:ndimMTO),ix(1:ndimMTO)),ovlm(ix(1:ndimMTO),ix(1:ndimMTO)), nmx,nev, evlmto, evecmto, oveps)
  if(nev/=ndimMTO) call rx('Hreduction: nev/=ndimMTO We didnot get eigenfuncitons of ndimMTO. Linear dependency problem?')
  ovlm= ovlmx
  hamm= hammx
  nmx = ndimPMT !  write(stdo,*)'Start Hreduction: 222'
  call zhev_tk4(ndimPMT,hamm(1:ndimPMT,1:ndimPMT),ovlm(1:ndimPMT,1:ndimPMT), nmx,nev, evl,evecpmt, oveps) !PMT
  ovlm=ovlmx
  ndimPMTx=nev
  do j=1,ndimMTO !wnm is corrected matrix element of <psi_PMT|psi_MTO>
     do i=1,ndimPMTx
        fac(i,j)= sum(dconjg(evecpmt(:,i))*matmul(ovlmx(:,ix(1:ndimMTO)),evecmto(1:ndimMTO,j)))
     enddo
  enddo
  ModifyMatrixElements :block
    integer:: ie,nidxevlmto,nidxevl,ibx,jx,idxevlmto(ndimMTO),idxevl(ndimPMT),jbx
    real(8):: eee,fffx,ecut,xxx,ewcutf,rydberg,facww
    real(8),allocatable::mulfac(:,:),mulfacw(:,:)
    !    complex(8):: wnj(ndimPMTx,ndimMTO),wnm(ndimPMTx,ndimMTO)
    allocate(wnj(ndimPMTx,ndimMTO),wnm(ndimPMTx,ndimMTO))!this is to avoid bug in ifort18.0.5
    ewcutf = ecutw+eferm
    do j=1,ndimMTO
       do i=1,ndimPMTx
          facww = facw*fermidist((evlmto(j)-ewcutf)/eww) 
!          wnm(i,j) = fac(i,j)*fac(i,j)**facww
          wnm(i,j) = fac(i,j)*abs(fac(i,j))**facww !2023-12-5 abs(fac) needed with PWMODE=11 to keep symmetry
       enddo
       !do i=1,ndimPMTx; if(j<4.and.abs(fac(i,j))>.1d0) write(6,*)' j=',j,i,' fac=',abs(fac(i,j)); enddo
    enddo
    call GramSchmidt(ndimPMTx,ndimMTO,wnm)
    if(iprx) then
       do j=1,ndimMTO !wnm is corrected matrix element of <psi_PMT|psi_MTO>
          do i=1,ndimPMTx
             if(abs(wnm(i,j))>.1) write(stdo,*)'wnm matrix ',j,i,abs(wnm(i,j))**2
          enddo
       enddo
    endif
    ! Mapping operator wnm*<psi_MTO|F_i>, where F_i is MTO basis.
    nx=ndimPMTx
    wnj = matmul(wnm(1:ndimPMTx,1:ndimMTO),matmul(transpose(dconjg(evecmto(:,:))),ovlmx(ix(1:ndimMTO),ix(1:ndimMTO))))
    do i=1,ndimMTO
       do j=1,ndimMTO
          hammout(i,j)= sum( dconjg(wnj(1:nx,i))*evl(1:nx)*wnj(1:nx,j))
          ovlmout(i,j)= sum( dconjg(wnj(1:nx,i))*wnj(1:nx,j) ) !<MLO|MLO>
       enddo
    enddo
  endblock ModifyMatrixElements
  return
contains
  real(8) function fermidist(x)
    real(8),intent(in) :: x
    if(x>100d0) then
       fermidist=0d0
    elseif(x<-100d0) then
       fermidist=1d0
    else
       fermidist=1d0/(exp(x)+1)
    endif
  end function fermidist
end subroutine Hreduction

