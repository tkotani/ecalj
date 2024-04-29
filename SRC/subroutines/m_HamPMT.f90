!> Read HamiltionanPMTinfo and HamiltonianPMT. Then convert HamPMT to HamRsMPO
module m_HamPMT
   use m_MPItk,only: procid, master_mpi, nsize,master
   use m_lgunit,only:stdo
   use m_ftox
   use m_lmfinit,only: oveps
   real(8),allocatable,protected:: plat(:,:),pos(:,:),qlat(:,:),symops(:,:,:)
   real(8),allocatable,protected,target:: qplist(:,:)
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
      use m_setqibz_lmfham,only: qibz,irotq,irotg,ndiff,iqbzrep,qbzii,igiqibz,nqibz,iqii,wiqibz,ngx,igx
      use m_zhev,only:zhev_tk4
      use m_readqplist,only: eferm
      use m_rotwave,only:  rotmatMTO!,rotmatPMT
      implicit none
      integer:: ifihmto,nqbz
      integer::ikpd,ikp,ib1,ib2,ifih,it,iq,nev,nmx,ifig=-999,i,j,ndimPMT,lold,m,ndimPMTx
      complex(8),allocatable:: hamm(:,:),ovlm(:,:), cmpo(:,:) 
      logical:: lprint=.true.,savez=.false.,getz=.false.,skipdiagtest=.false.
      complex(8):: img=(0d0,1d0),aaaa,phase
      real(8)::qp(3),pi=4d0*atan(1d0),fff,ef,fff1=2,fff2=2,fff3=0 ,facw,ecutw,eww,xxx,posd(3)
      integer:: nn,ib,k,l,ix5,imin,ixx,j2,j1,j3,nx,ix(ldim),iqini,iqend,ndiv
      integer:: ndimMTO !ndimMTO<ldim if we throw away f MTOs, for example.
      integer:: ib_tableM(ldim),k_tableM(ldim),l_tableM(ldim),ierr,ificmpo,iqibz,iqbz,igg
      logical:: cmdopt0
      real(8),pointer::qbz(:,:)
      complex(8),allocatable::ovlmi(:,:,:,:),hammi(:,:,:,:),rotmat(:,:)
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
      ndiv= nqibz/nsize
      if(nqibz>ndiv*nsize) ndiv=ndiv+1    !MPI division
      iqini =            ndiv*procid+1     !initial for each procid
      iqend = min(nqibz,ndiv*procid+ndiv) !end  for each procid
      write(stdo,ftox)'nsize procid iqini iqend=',nsize,procid,iqini,iqend
! Readin Hamiltonian only at iqibz
      allocate(ovlmi(1:ndimMTO,1:ndimMTO,nqibz,nspx),hammi(1:ndimMTO,1:ndimMTO,nqibz,nspx),source=(0d0,0d0))
      allocate(rotmat(ndimMTO,ndimMTO))
      HreductionIqibz: block
        logical:: cmdopt0
        integer:: i,iqxx,jspxx
        real(8)::eps=1d-8
        nqbz=nkp
        qbz=>qplist      
        iq=0
        iqiloop: do iqxx=1,nqibz !xx=1,nqibz !iqini,iqend !iqxx=1,nqibz 
           do jspxx=1,nspx
              read(ifih,end=2029) qp,ndimPMT,lso,xxx,jsp !jsp for isp; if so=1, jsp=1 only
              if(iqxx<iqini.or.iqend<iqxx) then
                 read(ifih)
                 read(ifih)
                 cycle
              endif   
              do i=iqibz,nqibz !Find iqibz and jsp here
                 if(sum(abs(qibz(:,i)-qp))<eps) then
                    iqibz=i
                 endif
              enddo
              write(stdo,"('=== Reading Ham for iq,iqibz,spin,q=',3i4,3f9.5)") iq,i,jsp,qp
              allocate(ovlm(1:ndimPMT,1:ndimPMT),hamm(1:ndimPMT,1:ndimPMT),cmpo(ndimPMT,ndimMTO))      !write(stdo,ftox)'qp=',ftof(qp)
              read(ifih) ovlm(1:ndimPMT,1:ndimPMT)
              read(ifih) hamm(1:ndimPMT,1:ndimPMT)
              call Hreduction(.false.,facw,ecutw,eww,ndimPMT,hamm(1:ndimPMT,1:ndimPMT),ovlm(1:ndimPMT,1:ndimPMT), &!Get reduced Hamitonian for ndimMTO
                   ndimMTO,ix,fff1, hamm(1:ndimMTO,1:ndimMTO),ovlm(1:ndimMTO,1:ndimMTO),cmpo)
              if(cmdopt0('--cmlo')) then  !at qibz only
                 writempo: block
                   character(8):: xt
                   iqibz=irotq(iq)
                   open(newunit=ificmpo, file='Cmpo' //trim(xt(i))//trim(xt(jsp)),form='unformatted')
                   write(ificmpo) ndimPMT,ndimMTO
                   write(ificmpo) cmpo
                   close(ificmpo)
                 endblock writempo
              endif
              hammi(:,:,iqibz,jsp)=0d0
              ovlmi(:,:,iqibz,jsp)=0d0
              do igg=1,ngx(iqibz) !symmetrized for rotations keeping qibz
                 call rotmatMTO(igg=igx(igg,iqibz),q=qibz(:,iqibz),qtarget=qibz(:,iqibz),ndimh=ndimMTO,rotmat=rotmat)
                 ovlmi(:,:,iqibz,jsp)=ovlmi(:,:,iqibz,jsp)&
                      +matmul(rotmat,matmul(ovlm(1:ndimMTO,1:ndimMTO),dconjg(transpose(rotmat))))
                 hammi(:,:,iqibz,jsp)=hammi(:,:,iqibz,jsp)&
                      +matmul(rotmat,matmul(hamm(1:ndimMTO,1:ndimMTO),dconjg(transpose(rotmat))))
              enddo
              hammi(:,:,iqibz,jsp)=hammi(:,:,iqibz,jsp)/ngx(iqibz)  
              ovlmi(:,:,iqibz,jsp)=ovlmi(:,:,iqibz,jsp)/ngx(iqibz)  
              deallocate(ovlm,hamm,cmpo)
           enddo
        enddo iqiloop
2029    continue
      endblock HreductionIqibz
      close(ifih)
      call mpibc2_complex(hammi,size(hammi),'m_HamPMT_hammi') 
      call mpibc2_complex(ovlmi,size(ovlmi),'m_HamPMT_ovlmi') 
      
      write(stdo,ftox)' --- goto qploop------'
      ndiv= nqbz/nsize
      if(nqbz>ndiv*nsize) ndiv=ndiv+1    !MPI division
      iqini =           ndiv*procid+1    !initial for each procid
      iqend = min(nqbz,ndiv*procid+ndiv) !end  for each procid
      write(stdo,ftox)'nqbz: nsize procid iqini iqend=',nsize,procid,iqini,iqend
      !iq=0
      allocate(ovlm(1:ndimMTO,1:ndimMTO),hamm(1:ndimMTO,1:ndimMTO))
      qploop: do iqbz=1,nqbz
         jsploop:do jsp=1,nspx
         !read(ifih,end=2019) qp,ndimPMT,lso,xxx,jsp !jsp for isp; if so=1, jsp=1 only
         !if(jsp==1) iq=iq+1
         if(iqbz<iqini.or.iqend<iqbz) cycle
         qp=qbz(:,iqbz)
         write(stdo,"('=== Ham for procid,iq,spin,q=',3i4,' ', 3f9.5)") procid,iqbz,jsp, qp
         !read(ifih) !ovlm(1:ndimPMT,1:ndimPMT)
         !read(ifih) !hamm(1:ndimPMT,1:ndimPMT)
         GETham_ndimMTO: block
            character(8):: xt
            real(8):: evlmlo(ndimMTO),cmptmto(ndimPMT,ndimMTO)
!            call Hreduction(.false.,facw,ecutw,eww,ndimPMT,hamm(1:ndimPMT,1:ndimPMT),ovlm(1:ndimPMT,1:ndimPMT), & !Get reduced Hamitonian for ndimMTO
!                 ndimMTO,ix,fff1, hamm(1:ndimMTO,1:ndimMTO),ovlm(1:ndimMTO,1:ndimMTO),cmpo)
!            iqbz = iq
            iqibz = irotq(iqbz)
            call rotmatMTO(igg=irotg(iqbz),q=qibz(:,iqibz),qtarget=qp+matmul(qlat,ndiff(:,iqibz)),ndimh=ndimMTO,rotmat=rotmat)
            ovlm = matmul(rotmat,matmul(ovlmi(:,:,iqibz,jsp),dconjg(transpose(rotmat))))
            hamm = matmul(rotmat,matmul(hammi(:,:,iqibz,jsp),dconjg(transpose(rotmat))))
            ! if(iq==3) then !            if(sum([qp(2),qp(3)]**2)<1d-3) then !check write
            !    Checkfinaleigen: block
            !       real(8):: rydberg
            !       complex(8):: evec(ndimMTO**2), hh(ndimMTO,ndimMTO),oo(ndimMTO,ndimMTO)
            !       hh = hamm(1:ndimMTO,1:ndimMTO)
            !       oo = ovlm(1:ndimMTO,1:ndimMTO)
            !       nmx= ndimMTO
            !       write(stdo,ftox) ' checkfinaleigen zhev_tk4',ftof(qp,3)
            !       call zhev_tk4(ndimMTO,hh,oo, nmx,nev,evlmlo, evec, oveps) !epsovl)
            !       write(stdo,ftox) '  evlmto=',nev,ftof(evlmlo(1:10)*rydberg())
            !       write(stdo,ftox) '  evlmto=',nev,ftof(evlmlo(11:20)*rydberg())
            !    endblock Checkfinaleigen
            ! endif
         endblock GETham_ndimMTO
         GETrealspaceHamiltonian: block ! H(k) ->  H(T) FourierTransformation to real space
            do i=1,ndimMTO
               do j=1,ndimMTO
                  ib1 = ib_tableM(i)
                  ib2 = ib_tableM(j)
                  do it =1,npair(ib1,ib2)! hammr_ij (T)= \sum_k hamm(k) exp(ikT). it is the index for T
                     phase = 1d0/dble(nkp)* exp(img*2d0*pi* sum(qp*(matmul(plat,nlat(:,it,ib1,ib2)))))
                     hammr(i,j,it,jsp)= hammr(i,j,it,jsp) + hamm(i,j) *phase
                     ovlmr(i,j,it,jsp)= ovlmr(i,j,it,jsp) + ovlm(i,j) *phase
                  enddo
               enddo
            enddo
          endblock GETrealspaceHamiltonian
       enddo jsploop
      enddo qploop
      deallocate(ovlm,hamm)!,cmpo)
!2019  continue
      call mpibc2_complex(hammr,size(hammr),'m_HamPMT_hammr') !to master
      call mpibc2_complex(ovlmr,size(ovlmr),'m_HamPMT_ovlmr') !to master
      if(master_mpi) then
         write(stdo,*)'Read: total # of q for Ham=',iq
         !! write RealSpace MTO Hamiltonian
         !ix(1:ndimMTO)=ix1(1:ndimMTO) !for atom idex
         write(stdo,*)' Writing HamRsMPO... ndimMTO=',ndimMTO
         open(newunit=ifihmto,file='HamRsMPO',form='unformatted')
         write(ifihmto) ndimMTO,npairmx,nspx
         write(ifihmto) hammr(1:ndimMTO,1:ndimMTO,1:npairmx,1:nspx),ovlmr(1:ndimMTO,1:ndimMTO,1:npairmx,1:nspx) !,ix(1:ndimMTO)
         write(ifihmto) ib_tableM(1:ndimMTO),k_tableM(1:ndimMTO),l_tableM(1:ndimMTO)
         close(ifihmto)
      endif
      if(master_mpi) write(stdo,*)" Wrote HamRsMPO file! End of lmfham1"
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
      open(newunit=ifihmto,file='HamRsMPO',form='unformatted')
      read(ifihmto) ndimMTO,npairmx,nspx
!    allocate(ix(ndimMTO))
      if(master_mpi) write(stdo,ftox)'MTOHamiltonian: ndimMTO,npairmx,nspx=',ndimMTO,npairmx,nspx
      allocate(ovlmr(1:ndimMTO,1:ndimMTO,npairmx,nspx), hammr(1:ndimMTO,1:ndimMTO,npairmx,nspx))
      read(ifihmto)hammr(1:ndimMTO,1:ndimMTO,1:npairmx,1:nspx),ovlmr(1:ndimMTO,1:ndimMTO,1:npairmx,1:nspx) !,ix(1:ndimMTO)
      allocate(ib_tableM(1:ndimMTO),k_tableM(1:ndimMTO),l_tableM(1:ndimMTO))
      read(ifihmto) ib_tableM(1:ndimMTO),k_tableM(1:ndimMTO),l_tableM(1:ndimMTO)
      close(ifihmto)
      if(master_mpi) write(stdo,*)'OK: Read HamRsMPO file! Use i-ioffib for setting <Worb>'
   end subroutine ReadHamRsMPO
end module m_HamRsMPO
subroutine Hreduction(iprx,facw,ecutw,eww,ndimPMT,hamm,ovlm,ndimMTO,ix,fff1, hammout,ovlmout, cmpo) !> Reduce H(ndimPMT) to H(ndimMTO)
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
   complex(8):: hammout(ndimMTO,ndimMTO),ovlmout(ndimMTO,ndimMTO) , cmpo(ndimPMT,ndimMTO)
   complex(8),allocatable :: wnm(:,:)
   real(8):: fff1,fff !epsovl=1d-8 epsovlm=0d0 ,
   logical:: iprx
   logical:: cmdopt0
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
   ndimPMTx=nev !obtained. oveps may reduce ndimPMT to be ndimPMTx
   do j=1,ndimMTO !wnm is corrected matrix element of <psi_PMT|psi_MTO>
      do i=1,nev
         fac(i,j)= sum(dconjg(evecpmt(:,i))*matmul(ovlmx(:,ix(1:ndimMTO)),evecmto(1:ndimMTO,j)))
      enddo
   enddo
   ModifyMatrixElements :block
      integer:: ie,nidxevlmto,nidxevl,ibx,jx,idxevlmto(ndimMTO),idxevl(ndimPMT),jbx
      real(8):: eee,fffx,ecut,xxx,ewcutf,rydberg,facww
      real(8),allocatable::mulfac(:,:),mulfacw(:,:)
      allocate(wnm(ndimPMTx,ndimMTO))!this is to avoid bug in ifort18.0.5
      ewcutf = ecutw+eferm
      do j=1,ndimMTO
         do i=1,ndimPMTx
            facww = facw*fermidist((evlmto(j)-ewcutf)/eww)
            wnm(i,j) = fac(i,j)*abs(fac(i,j))**facww !2023-12-5 abs(fac) needed with PWMODE=11 to keep symmetry
            !  wnm(i,j) = fac(i,j)*fac(i,j)**facww
         enddo
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
      cmpo(ndimPMTx+1:ndimPMT,1:ndimMTO)=0d0
      cmpo(1:ndimPMTx,1:ndimMTO) = matmul(wnm(1:ndimPMTx,1:ndimMTO),&
         matmul(transpose(dconjg(evecmto(:,:))),ovlmx(ix(1:ndimMTO),ix(1:ndimMTO))))
      do i=1,ndimMTO
         do j=1,ndimMTO
            hammout(i,j)= sum( dconjg(cmpo(1:nx,i))*evl(1:nx)*cmpo(1:nx,j)) !|FMPO_i>=|PsiPMT_j> cmpo(j,i)
            ovlmout(i,j)= sum( dconjg(cmpo(1:nx,i))*cmpo(1:nx,j) ) !<MPO|MPO>
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

