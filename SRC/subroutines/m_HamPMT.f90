!> Read HamiltionanPMTinfo and HamiltonianPMT. Then convert HamPMT to HamRsMLO
module m_HamPMT
   use m_MPItk,only: procid, master_mpi, nsize,master,strprocid
   use m_lgunit,only:stdo
   use m_ftox
   use m_lmfinit,only: oveps
   use m_keyvalue,only: getkeyvalue
   use m_hreduction,only: hreduction
   real(8),external::tolq !eps=1d-8
   real(8),allocatable,protected:: plat(:,:),pos(:,:),qlat(:,:),symops(:,:,:)
   real(8),allocatable,protected,target:: qplist(:,:)
   integer,allocatable,protected:: nlat(:,:,:,:),npair(:,:),ib_table(:),l_table(:),k_table(:),ispec_table(:),nqwgt(:,:,:),m_table(:)
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
      write(stdo,ftox)'readhamPMTinfo nkp=',nkp
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
      write(stdo,"('MHAM: --- MTO part of PMT Hamiltonian index (real-harmonics table is in job_pdos script) --- ')")
      write(stdo,'("MHAM: MTO block dim=",i5)') ldim
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
         write(stdo,"('MHAM: i i-ioffib ib(atom) l k(1:EH,2:EH2,3:PZ)=',i4,5i3)")&
            i,i-ioff,ib_table(i),l_table(i),k_table(i)
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
   subroutine HamPMTtoHamRsMLO()!eww) !Convert HamPMT(k mesh) to HamRsMLO(real space)
      use m_setqibz_lmfham,only: qibz,irotq,irotg,ndiff,iqbzrep,qbzii,igiqibz,nqibz,iqii,wiqibz,ngx,igx
      use m_zhev,only:zhev_tk4
      use m_readqplist,only: eferm
      use m_rotwave,only:  rotmatMTO!,rotmatPMT
      implicit none
      integer:: ifihmto,nqbz
      integer::ikpd,ikp,ib1,ib2,ifih,it,iq,nev,nmx,ifig=-999,i,j,ndimPMT,lold,m,ndimPMTmx
      complex(8),allocatable:: hamm(:,:),ovlm(:,:)
      logical:: lprint=.true.,savez=.false.,getz=.false.,skipdiagtest=.false.
      complex(8):: img=(0d0,1d0),aaaa,phase
      real(8)::qp(3),pi=4d0*atan(1d0),fff,ef,fff1=2,fff2=2,fff3=0 ,xxx,posd(3) !,ecutw,eww
      integer:: nn,ib,k,l,ix5,imin,ixx,j2,j1,j3,nx,ix(ldim),iqini,iqend,ndiv,ifihh,niqisp,nqirr,numprocs_in
      integer:: ndimMTO !ndimMTO<ldim if we throw away f MTOs, for example.

      integer:: ib_tableM(ldim),k_tableM(ldim),l_tableM(ldim),ierr,iqibz,iqbz,igg,nMTO,mlomethod,nskip !,procid_in,numprocs_in
      logical:: cmdopt0
      integer, allocatable:: ib_tableI(:)
      real(8),pointer::qbz(:,:)
      complex(8),allocatable::ovlmi(:,:,:,:),hammi(:,:,:,:),rotmat(:,:)
      integer,allocatable::ndimPMTq(:), iqproc(:),isproc(:)
      logical,allocatable:: lqibz(:)
      logical::debug=.true.
      
      ReadInfoFromGWinput: block ! Input orbital index for MLO, stored into idmto (s,p,d=1,2,3,4,5,6,7,8,9)
        use m_nvfortran,only : findloc
        integer::lmindex(16,nbas),ifmloc,ret,lm
        character(256):: labl,aaa
        call getkeyvalue("GWinput","mlo_method",mlomethod,default=0)
!        mlomethod=-999
        call getkeyvalue("GWinput","<Worb>",unit=ifmloc,status=ret)
        do 
          read(ifmloc,"(a)") aaa
          if(aaa(1:1) == '!') then
            read(aaa,*)
            cycle
          endif
          aaa=trim(aaa)//repeat(' -999 ',16)
          read(aaa,*,end=1201,err=1201) ib,labl,lmindex(1:16,ib)
        enddo
1201    continue
        close(ifmloc)
        nn=0
        lold=-999
!        nskip=0
        do i=1,ldim  !Only MTOs for EH 
          if( k_table(i)==2) cycle  !skip 2nd
          if( k_table(i)==3) then
!            nskip=nskip+1 
            cycle  !skip local orbita l! we assume nskip is for local orbtail to skip 
          endif  
          ib=ib_table(i)
          if(lold/=l_table(i)) then
            m= -l_table(i)
          else
            m=m+1
          endif
          lm = l_table(i)**2 + l_table(i) + m +1
          if(.not. any(lm==lmindex(1:16,ib))) cycle
          !     if(cmdopt0('--skipf') .and.   l_table(i)>=3) cycle ! skip f orbitals. !if(k_table(i)==2.and.l_table(i)>=2) cycle ! throw away EH2 for d
          !     if(cmdopt0('--skipd') .and.   l_table(i)>=2) cycle ! throw away EH2 for d
          !     if(cmdopt0('--skip2nd') .and. k_table(i)==2) cycle 
          !     if(cmdopt0('--skip2ndp').and. k_table(i)==2.and.l_table(i)==1) cycle 
          !     if(cmdopt0('--skip2nds').and. k_table(i)==2.and.l_table(i)==0) cycle  !skip EH2 
          !     if(cmdopt0('--skiplo') .and. k_table(i)==3) cycle !skip local orbital
          !     if(cmdopt0('--skip2ndd') .and. k_table(i)>=2.and.l_table(i)>=2) cycle ! throw away EH2 for d
          !     if(cmdopt0('--skip1d') .and.  k_table(i)==1 .and. l_table(i)==2) cycle ! throw away EH2 for d
          nn=nn+1
          ix(nn)=i
          ib_tableM(nn)= ib_table(i)
          k_tableM(nn) = k_table(i)
          l_tableM(nn) = l_table(i)
        enddo
      endblock ReadInfoFromGWinput
      ndimMTO=nn
      if(lso==1) ndimMTO=nn*2 !L.S mode
      nMTO=ldim
      nspx=nsp
      if(lso==1) nspx=1
      ! Readin Hamiltonian only at iqibz
      ib_tableI = pack(ib_tableM(1:ndimMTO), [(all(ib_tableM(:i-1)/=ib_tableM(i)), i=1,ndimMTO)])
      allocate(ovlmi(1:ndimMTO,1:ndimMTO,nqibz,nspx),hammi(1:ndimMTO,1:ndimMTO,nqibz,nspx),source=(0d0,0d0))
      allocate(rotmat(nMTO,nMTO))
      allocate(ndimPMTq(nqibz),source=0)

!2026-1-27      
      cmlo4GWinput: if(cmdopt0('--mlo')) then !from __Hamiltoniangw to __cmlo.data, __cmlo.info
        HreductionIqibzGWinput: block
          use m_mpiio,only: openm,writem,closem
          integer:: i,iqxx,jspxx,idat,ifizz,isp,mrecbb,ndble,ifi,nbandmx,nqbzgw
          complex(8)::rotmatt(ndimMTO,ndimMTO)
          real(8),allocatable:: qplistgw(:,:)
          !for sugw output for GWinput to get zcplz for q point in qg4gw.
          open(newunit=ifihh,file='__Hamiltoniangw.'//trim(strprocid),form='unformatted') !sugw for GWinput
          read(ifihh) niqisp,nqirr,nbandmx,numprocs_in,nqbzgw !nbandmx=ndimPMTmx: max dim of PMT Hamiltonian for all q
          allocate(  qplistgw(1:3,1:nqirr),iqproc(1:niqisp),isproc(1:niqisp))
          read(ifihh)qplistgw(1:3,1:nqirr),iqproc(1:niqisp),isproc(1:niqisp)
          ndble = 8
          mrecbb = 2*nbandmx*ndimMTO* ndble !byte size  !Use -assume byterecl for ifort, so that ifort recognizes the recored in the unit of bytes.
          i = openm(newunit=ifizz, file='__cmlo.data',recl=mrecbb)
          if(numprocs_in/=nsize) call rxii('m_HamPMT: nsize for lmf and lmfham1 should be the same', numprocs_in,nsize)
          iqibzloops: do idat=1,niqisp
            iq  = iqproc(idat) ! iq index
            isp = isproc(idat) ! spin index: Note isp=1:nspx, where nspx=nsp/nspc. See sugw.f90
            qp  = qplistgw(:,iq) ! q vector containing nqirr
            read(ifihh) ndimPMT
            write(stdo,ftox)' iqibzloops: m_HamPMT for GWinput=',procid,' iq isp=',iq,isp,' q=',ftof(qp)
            !          if(.not.lqibz(iq) ) cycle ! if qp is not qibz in GWinput
            write(stdo,ftox)'=== Reading Ham for iqibz spin procid q= ', iq,jsp,procid,ftof(qp)
            readingovlmp: block
              integer:: ificmlo,iqqisp
              character(8):: xt
              logical:: cmdopt0
              complex(8):: ovlmp(1:ndimPMT,1:ndimPMT),hammp(1:ndimPMT,1:ndimPMT),cmlo(nbandmx,ndimMTO),&
                   ovlm(1:ndimMTO,1:ndimMTO),hamm(1:ndimMTO,1:ndimMTO)
              read(ifihh) ovlmp
              read(ifihh) hammp
              cmlo=0d0 !zero padding for 1:nbandmx in advance
              call Hreduction(mlomethod,.false.,ndimPMT,hammp,ovlmp, ndimMTO,ix,fff1, hamm,ovlm,qp,cmlo(1:ndimPMT,1:ndimMTO))
              !                                                                          Get reduced Hamitonian for ndimMTO
              if(cmdopt0('--mlo')) then  !at qibz only
                iqqisp= isp + nspx*(iq-1)
!                write(*,*)'cccccccccc cmlowrite',isp,iq,iqqisp, sum(abs(cmlo))
                i = writem(ifizz,rec=iqqisp,data=cmlo) 
              endif
            endblock readingovlmp
          enddo iqibzloops
2019      continue
          i=closem(ifizz)
          close(ifihh)
          if(master_mpi) then
            open(newunit=ifi,file='__cmlo.info',form='unformatted') 
            write(ifi) ndimMTO,nqbzgw,nqirr,nMTO,mrecbb
            write(stdo,ftox)'nnnnn ndimMTO nqirr=',ndimMTO,nqirr
            write(ifi) ix(1:ndimMTO),qplistgw(1:3,1:nqirr)
            close(ifi)
          endif
        endblock HreductionIqibzGWinput
      endif cmlo4GWinput
    
! --- base line for ctrl.foobar
      HreductionIqibz: block
        integer:: i,iqxx,jspxx,idat,isp,ndble
        complex(8)::rotmatt(ndimMTO,ndimMTO)
        open(newunit=ifih, file='HamiltonianPMT.'//trim(strprocid),form='unformatted')
        read(ifih)  !procid_in,numprocs_in
        iqiloop: do iqxx=1,nqibz !nqibz !xx=1,nqibz !iqini,iqend !iqxx=1,nqibz 
           if(debug)write(6,*)' start iqiloop=',iqxx,nqibz
           do jspxx=1,nspx
              if(debug)write(6,*)' jspxx=',jspxx
              read(ifih,end=2029) qp,ndimPMT,lso,xxx,jsp !jsp for isp; if so=1, jsp=1 only
              iqibz = findloc( [(sum(abs(qibz(:,i)-qp))<tolq(),i=1,nqibz)],value=.true.,dim=1)
              write(stdo,ftox)'=== Reading Ham for iqibz spin procid q= ', iqibz,jsp,procid,ftof(qp)
              allocate(ovlm(1:ndimMTO,1:ndimMTO),hamm(1:ndimMTO,1:ndimMTO))
              readingovlmp2: block
                character(8):: xt
                logical:: cmdopt0
                complex(8):: ovlmp(1:ndimPMT,1:ndimPMT),hammp(1:ndimPMT,1:ndimPMT),cmlo(ndimPMT,ndimMTO)
                read(ifih) ovlmp
                read(ifih) hammp
                call Hreduction(mlomethod,.false.,ndimPMT,hammp,ovlmp, ndimMTO,ix,fff1, hamm,ovlm,qp,cmlo) !Get reduced Hamitonian for ndimMTO
              endblock readingovlmp2
              ndimPMTq(iqibz)=ndimPMT
              do igg=1,ngx(iqibz) !symmetrized for rotations keeping qibz
                 call rotmatMTO(igg=igx(igg,iqibz),q=qibz(:,iqibz),qtarget=qibz(:,iqibz),ndimh=nMTO,rotmat=rotmat)
                 !associate( rotmatt=>rotmat(ix(1:ndimMTO),ix(1:ndimMTO)))
                   forall(i=1:ndimMTO,j=1:ndimMTO) rotmatt(i,j)=rotmat(ix(i),ix(j))
                   ovlmi(:,:,iqibz,jsp)=ovlmi(:,:,iqibz,jsp) +matmul(rotmatt,matmul(ovlm,dconjg(transpose(rotmatt))))
                   hammi(:,:,iqibz,jsp)=hammi(:,:,iqibz,jsp) +matmul(rotmatt,matmul(hamm,dconjg(transpose(rotmatt))))
                 !endassociate
              enddo
              hammi(:,:,iqibz,jsp)=hammi(:,:,iqibz,jsp) /ngx(iqibz)  
              ovlmi(:,:,iqibz,jsp)=ovlmi(:,:,iqibz,jsp) /ngx(iqibz)
              deallocate(ovlm,hamm)
           enddo
           if(debug)write(6,*)' end of iqiloop=',iqxx,nqibz
        enddo iqiloop
2029    continue
        close(ifih)
      endblock HreductionIqibz
      call mpibc2_complex(hammi,size(hammi),'m_HamPMT_hammi') 
      call mpibc2_complex(ovlmi,size(ovlmi),'m_HamPMT_ovlmi') 
      call mpibc2_int(ndimPMTq,size(ndimPMTq),'m_HamPMT_ndimPMTq')
! hammr ovlmr     
      ! allocate(ovlmr(1:ndimMTO,1:ndimMTO,npairmx,nspx), hammr(1:ndimMTO,1:ndimMTO,npairmx,nspx),source=(0d0,0d0))
      allocate(ovlmr(npairmx,ndimMTO,ndimMTO,nspx), source=(0d0,0d0))
      allocate(hammr(npairmx,ndimMTO,ndimMTO,nspx), source=(0d0,0d0))
      nqbz=nkp
      qbz=>qplist      
      ndiv= nqbz/nsize
      if(nqbz>ndiv*nsize) ndiv=ndiv+1  !MPI division
      iqini =      ndiv*procid+1       !initial for each procid
      iqend =      ndiv*procid+ndiv    !end  for each procid
      if(iqini>nqbz) then
         iqini=0
         iqend=-1
      elseif(iqend>nqbz) then
         iqend=nqbz
      endif
      !iqend = min(nqbz,ndiv*procid+ndiv) !end  for each procid
      write(stdo,ftox)'nnnn nsize procid iqini iqend=',nsize,procid,iqini,iqend,'  ',ndiv
      qploop: do iqbz=iqini,iqend
        qp     = qbz(:,iqbz)
        iqibz  = irotq(iqbz)
        ndimPMT= ndimPMTq(iqibz)
        hammovlm: block
          complex(8):: ovlm(1:ndimPMT,1:ndimPMT),hamm(1:ndimPMT,1:ndimPMT),rotmatt(ndimMTO,ndimMTO)
          do jsp=1,nspx
            if(master_mpi) write(stdo,ftox)'=== Rotate Ham from iqibz to iqbz; iqibz iqbz isp q=',iqibz,iqbz,jsp,'q ig=',ftof(qp,4),irotg(iqbz)
            call rotmatMTO(igg=irotg(iqbz),q=qibz(:,iqibz),qtarget=qp+matmul(qlat,ndiff(:,iqibz)),ndimh=nMTO,rotmat=rotmat)
            forall(i=1:ndimMTO,j=1:ndimMTO) rotmatt(i,j)=rotmat(ix(i),ix(j))
            ovlm(1:ndimMTO,1:ndimMTO) = matmul(rotmatt,matmul(ovlmi(:,:,iqibz,jsp),dconjg(transpose(rotmatt))))
            hamm(1:ndimMTO,1:ndimMTO) = matmul(rotmatt,matmul(hammi(:,:,iqibz,jsp),dconjg(transpose(rotmatt))))
            GETrealspaceHamiltonian:block !optimized version
              integer :: ibt1, ibt2, np, ii, jj
              integer, allocatable :: idims(:), jdims(:)
              complex(8) :: phases(npairmx)
              do ibt1 = 1, size(ib_tableI)
                do ibt2 = 1, size(ib_tableI)
                  ib1 = ib_tableI(ibt1)
                  ib2 = ib_tableI(ibt2)
                  np = npair(ib1,ib2)
                  idims = pack([(i,i=1,ndimMTO)], mask=(ib_tableM(:)==ib1))
                  jdims = pack([(j,j=1,ndimMTO)], mask=(ib_tableM(:)==ib2))
                  phases(1:np) = [(1d0/dble(nqbz)* exp(img*2d0*pi* sum(qp*(matmul(plat,nlat(:,it,ib1,ib2))))),it=1,np)]
                  do concurrent(it=1:np, ii=1:size(idims), jj=1:size(jdims))
                    hammr(it,idims(ii),jdims(jj),jsp) = hammr(it,idims(ii),jdims(jj),jsp) + hamm(idims(ii),jdims(jj))*phases(it)
                    ovlmr(it,idims(ii),jdims(jj),jsp) = ovlmr(it,idims(ii),jdims(jj),jsp) + ovlm(idims(ii),jdims(jj))*phases(it)
                 enddo
                enddo
              enddo
            endblock GetrealspaceHamiltonian
            ! GETrealspaceHamiltonian: block ! H(k) ->  H(T) FourierTransformation to real space
            !   do i=1,ndimMTO; do j=1,ndimMTO
            !     ib1 = ib_tableM(i)
            !     ib2 = ib_tableM(j)
            !     do it =1,npair(ib1,ib2)! hammr_ij (T)= \sum_k hamm(k) exp(ikT). it is the index for T
            !       phase = 1d0/dble(nqbz)* exp(img*2d0*pi* sum(qp*(matmul(plat,nlat(:,it,ib1,ib2)))))
            !       hammr(i,j,it,jsp)= hammr(i,j,it,jsp)+ hamm(i,j)*phase
            !       ovlmr(i,j,it,jsp)= ovlmr(i,j,it,jsp)+ ovlm(i,j)*phase
            !     enddo
            !   enddo; enddo
            ! endblock GETrealspaceHamiltonian
          enddo
        endblock hammovlm
      enddo qploop
      call mpibc2_complex(hammr,size(hammr),'m_HamPMT_hammr') !to master
      call mpibc2_complex(ovlmr,size(ovlmr),'m_HamPMT_ovlmr') !to master
      if(master_mpi) then ! write RealSpace MTO Hamiltonian          !ix(1:ndimMTO)=ix1(1:ndimMTO) !for atom idex
        write(stdo,*)' Writing HamRsMLO... ndimMTO=',ndimMTO
        open(newunit=ifihmto,file='HamRsMLO',form='unformatted')
        write(ifihmto) ndimMTO,npairmx,nspx
        ! write(ifihmto) hammr(1:ndimMTO,1:ndimMTO,1:npairmx,1:nspx)
        ! write(ifihmto) ovlmr(1:ndimMTO,1:ndimMTO,1:npairmx,1:nspx) !,ix(1:ndimMTO)
        write(ifihmto) hammr(1:npairmx,1:ndimMTO,1:ndimMTO,1:nspx)
        write(ifihmto) ovlmr(1:npairmx,1:ndimMTO,1:ndimMTO,1:nspx) !,ix(1:ndimMTO)
        write(ifihmto) ib_tableM(1:ndimMTO),k_tableM(1:ndimMTO),l_tableM(1:ndimMTO)
        close(ifihmto)
        write(stdo,*)" Wrote HamRsMLO file! End of lmfham1"
      endif
   end subroutine HamPMTtoHamRsMLO
end module m_HamPMT


module m_HamRsMLO ! read real-space MLO Hamiltonian
   use m_lgunit,only:stdo
   use m_ftox
   integer,protected:: ndimMTO,npairmx,nspx  !ndimMTO<ldim if we throw away f MTOs, for example.
   integer,allocatable,protected:: ib_tableM(:),l_tableM(:),k_tableM(:)
   complex(8),allocatable,protected:: ovlmr(:,:,:,:),hammr(:,:,:,:)
contains
   subroutine ReadHamRsMLO()! read RealSpace MTO Hamiltonian
      use m_MPItk,only: master_mpi
      integer:: ifihmto
      open(newunit=ifihmto,file='HamRsMLO',form='unformatted')
      read(ifihmto) ndimMTO,npairmx,nspx !    allocate(ix(ndimMTO))
      if(master_mpi) write(stdo,ftox)'MTOHamiltonian: ndimMTO,npairmx,nspx=',ndimMTO,npairmx,nspx
      allocate(ovlmr(1:ndimMTO,1:ndimMTO,npairmx,nspx), hammr(1:ndimMTO,1:ndimMTO,npairmx,nspx))
      ! read(ifihmto) hammr(1:ndimMTO,1:ndimMTO,1:npairmx,1:nspx)
      ! read(ifihmto) ovlmr(1:ndimMTO,1:ndimMTO,1:npairmx,1:nspx) !,ix(1:ndimMTO)
      read(ifihmto) hammr(1:npairmx,1:ndimMTO,1:ndimMTO,1:nspx)
      read(ifihmto) ovlmr(1:npairmx,1:ndimMTO,1:ndimMTO,1:nspx) !,ix(1:ndimMTO)
      allocate(ib_tableM(1:ndimMTO),k_tableM(1:ndimMTO),l_tableM(1:ndimMTO))
      read(ifihmto) ib_tableM(1:ndimMTO),k_tableM(1:ndimMTO),l_tableM(1:ndimMTO)
      close(ifihmto)
      if(master_mpi) write(stdo,*)'OK: Read HamRsMLO file! Use i-ioffib for setting <Worb>'
   end subroutine ReadHamRsMLO
end module m_HamRsMLO
