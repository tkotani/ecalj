!! -- Read HamiltionanPMTinfo and HamiltonianPMT. Then convert HamPMT to HamRsMTO  ------
module m_HamPMT
  real(8),allocatable,protected:: plat(:,:),pos(:,:),qplist(:,:),qlat(:,:)
  integer,allocatable,protected:: nlat(:,:,:,:),npair(:,:),&
       ib_table(:),l_table(:),k_table(:),nqwgt(:,:,:)
  integer,protected:: nkk1,nkk2,nkk3,nbas,nkp,npairmx,ldim,jsp,lso,nsp,nspx
  real(8),protected:: epsovl,alat
  complex(8),allocatable,protected:: ovlmr(:,:,:,:),hammr(:,:,:,:)
  integer:: ndimMTO
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
    use m_readqplist,only: eferm
    implicit none
    integer:: ifihmto
    integer:: ifile_handle,ikpd,ikp,ib1,ib2,ifih,it,iq,nev,nmx,ifig=-999,i,j,ndimPMT,lold,m
    complex(8),allocatable:: hamm(:,:),ovlm(:,:)!,t_zv(:,:)!,ovlmx(:,:)
    logical:: lprint=.true.,savez=.false.,getz=.false.,skipdiagtest=.false.
!    real(8),allocatable:: evl(:)
    complex(8):: img=(0d0,1d0),aaaa,phase
    real(8)::qp(3),pi=4d0*atan(1d0),fff,ef,fff1=8,fff2=4 !,fff1=2,fff2=2
    integer:: ix(ldim),ixm(ldim),iix(ldim),nnn,ib,k,l,ix5,ix21,imin,ixx,ndimMTO2,j2
    integer:: nx
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1    
    nnn=0
    do i=1,ldim
       if(l_table(i)>=3) cycle
!       if(l_table(i)>=2.and.k_table(i)==2) cycle
       if(k_table(i)==2) cycle
       write(6,*) i,ib_table(i),l_table(i),k_table(i)
       nnn=nnn+1
       ixm(nnn)=i
    enddo
    ndimMTO2=nnn
!
    nnn=0
    do i=1,ldim !MLO dimention 
       if(l_table(i)>=3) cycle
       !if(l_table(i)>=2.and.k_table(i)>=2) cycle
!       if(k_table(i)==2) cycle
       write(6,*) i,ib_table(i),l_table(i),k_table(i)
       nnn=nnn+1
       ix(nnn)=i
    enddo
    ndimMTO=nnn
!
    do j2=1,ndimMTO2
       do j=1,ndimMTO
          if(ixm(j2)==ix(j)) then
             iix(j2)=j        !index for MLO ndimMTO
             !exit
          endif
       enddo
    enddo
    
    open(newunit=ifih,file='HamiltonianPMT',form='unformatted')
    write(6,*)'Reaing HamiltonianPMT...'
    if(lso==1) ldim=ldim*2 !L.S mode
    
    nspx=nsp
    if(lso==1) nspx=1
    allocate(ovlmr(1:ndimMTO,1:ndimMTO,npairmx,nspx), hammr(1:ndimMTO,1:ndimMTO,npairmx,nspx))
    write(6,"('ndimMTO ldim lso=',i6,4i3)") ndimMTO,ldim,lso
    
    hammr=0d0
    ovlmr=0d0
    iq=0
    qploop: do 
       read(ifih,end=2019) qp,ndimPMT,lso,epsovl,jsp !jsp for isp; if so=1, jsp=1 only
       if(jsp==1) iq=iq+1
       write(6,"('=== Reading Ham for iq,spin,q=',2i4,3f9.5)") iq,jsp,qp
       !c if(ndimPMT/=ldim.and.(lso==0.or.lso==2)) call rx('lmfham:   ndimMTO/=ldim')
       !c if(ndimPMT/=2*ldim.and.lso==1)        call rx('lmfham: 2*ndimMTO/=ldim') ! L.S mode or not
       allocate(ovlm(1:ndimPMT,1:ndimPMT),hamm(1:ndimPMT,1:ndimPMT))
       read(ifih) ovlm(1:ndimPMT,1:ndimPMT)
       read(ifih) hamm(1:ndimPMT,1:ndimPMT)
       epsovl=1d-8

       eigen: block
         real(8):: evlmlo(ndimMTO),evl(ndimPMT) !,evlmlo2(ndimMTO2),evlmlo2x(ndimMTO2)
         ! <PsiPMT|PsiMLO> 
         call Hreduction(ndimPMT,hamm(1:ndimPMT,1:ndimPMT),ovlm(1:ndimPMT,1:ndimPMT), &
              ndimMTO,ix,fff1,   evl,hamm(1:ndimMTO,1:ndimMTO),ovlm(1:ndimMTO,1:ndimMTO))
         ! <PsiMLO|PsiMLO2> 
         call Hreduction(ndimMTO,hamm(1:ndimMTO,1:ndimMTO),ovlm(1:ndimMTO,1:ndimMTO),&
              ndimMTO2,iix,fff2, evlmlo,hamm(1:ndimMTO2,1:ndimMTO2),ovlm(1:ndimMTO2,1:ndimMTO2))

         echeck: block
           complex(8):: evecmlo2(ndimMTO2,ndimMTO2)
           real(8):: evlmlo2(ndimMTO2),evlmlo2x(ndimMTO2)
           complex(8):: hh(1:ndimMTO2,1:ndimMTO2),oo(1:ndimMTO2,1:ndimMTO2)
           if(sum([qp(2),qp(3)]**2)<1d-3) then
              hh=hamm(1:ndimMTO2,1:ndimMTO2)
              oo=ovlm(1:ndimMTO2,1:ndimMTO2)
              nmx = ndimMTO2
              call zhev_tk4(ndimMTO2,hh,oo, nmx,nev,evlmlo2x, evecmlo2, epsovl) !PsiMLO
              do i=1,ndimMTO2
                 write(6,"(a,i5,3f12.4)")'mmmmm', i,evl(i),evlmlo(i),evlmlo2x(i)-evlmlo(i)
              enddo
              !if(abs(qp(1)-1)<1d-4) stop 'vvvvvvvvqp'
           endif
         endblock echeck
       endblock eigen
       !! Real space Hamiltonian. H(k) ->  H(T) FourierTransformation to real space
       !!       Only MTO part ndimMTO (ndimPMT = ndimMTO + ndimAPW)
       do i=1,ndimMTO2
          do j=1,ndimMTO2
             ib1 = ib_table(ixm(i)) 
             ib2 = ib_table(ixm(j)) 
             do it =1,npair(ib1,ib2)! hammr_ij (T)= \sum_k hamm(k) exp(ikT). it is the index for T
                phase = 1d0/dble(nkp)* exp(img*2d0*pi* sum(qp*matmul(plat,nlat(:,it,ib1,ib2))))
                hammr(i,j,it,jsp)= hammr(i,j,it,jsp)+ hamm(i,j)*phase
                ovlmr(i,j,it,jsp)= ovlmr(i,j,it,jsp)+ ovlm(i,j)*phase
             enddo
          enddo
       enddo
       deallocate(ovlm,hamm)
    enddo qploop
2019 continue
    write(6,*)'Read: total # of q for Ham=',iq
    close(ifih)
   
    ndimMTO=ndimMTO2
    ix(1:ndimMTO)=ixm(1:ndimMTO) !for atom idex
   
    !! write RealSpace MTO Hamiltonian
    write(6,*)' Writing HamRsMTO... ndimMTO=',ndimMTO
    open(newunit=ifihmto,file='HamRsMTO',form='unformatted')
    write(ifihmto) ndimMTO,npairmx,nspx
    write(ifihmto) hammr(1:ndimMTO,1:ndimMTO,1:npairmx,1:nspx),&
         ovlmr(1:ndimMTO,1:ndimMTO,1:npairmx,1:nspx),ix(1:ndimMTO)
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
    do while(abs(valn-val)>1d-10)
       valn= sum(1d0/(exp(beta*(evl(:)-emu))+1d0))
       if(valn<val) emu=emu+demu
       if(valn>val) then
          emu = emu-demu
          demu= 0.5d0*demu
       endif
       ix=ix+1
    enddo
  end subroutine detemu
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
  integer,protected:: ndimMTO,npairmx,nspx
  integer,allocatable,protected:: ix(:)
  complex(8),allocatable,protected:: ovlmr(:,:,:,:),hammr(:,:,:,:)
contains
  !! read RealSpace MTO Hamiltonian
  subroutine ReadHamRsMTO()
    integer:: ifihmto,ifile_handle
    open(newunit=ifihmto,file='HamRsMTO',form='unformatted')
    read(ifihmto) ndimMTO,npairmx,nspx
    allocate(ix(ndimMTO))
    write(6,*)'ndimMTO,npairmx,nspx=',ndimMTO,npairmx,nspx
    allocate(ovlmr(1:ndimMTO,1:ndimMTO,npairmx,nspx), hammr(1:ndimMTO,1:ndimMTO,npairmx,nspx))
    read(ifihmto)hammr(1:ndimMTO,1:ndimMTO,1:npairmx,1:nspx),ovlmr(1:ndimMTO,1:ndimMTO,1:npairmx,1:nspx),ix(1:ndimMTO)
    close(ifihmto)
    print*,' OK: Read HamRsMTO file!'
  end subroutine ReadHamRsMTO
end module m_HamRsMTO
!! -----------------------------------------------------------------------------------      
program lmfham1
  !! Read HamiltonianPMT and generates MTO-only Hamiltonian
  use m_HamPMT, only: plat, npair,nlat,nqwgt, ldim, nkp,qplist,&
       ib_table,alat, ReadHamPMTInfo, HamPMTtoHamRsMTO
  ! note.  HamPMTtoHamRsMTO do not change variables. Only generate HamRsMTO file.
  use m_HamRsMTO,  only: hammr,ovlmr,ndimMTO,  ReadHamRsMTO,ix,npairmx,nspx
  use m_readqplist,only: eferm,qplistsy,ndat,xdat, Readqplistsy
  use m_mpi,only: MPI__hx0fp0_rankdivider2Q, MPI__Qtask, &
       MPI__Initialize, MPI__Finalize,MPI__root, &
       MPI__Broadcast, MPI__rank,MPI__size, MPI__consoleout,MPI__barrier
  implicit none
  integer:: i,j,ikp,ib1,ib2,it,nmx,nev,jsp
  complex(8)::img=(0d0,1d0),phase
  complex(8),allocatable:: hamm(:,:),ovlm(:,:),t_zv(:,:)
  real(8),allocatable:: evl(:)
  real(8)::qp(3),pi=4d0*atan(1d0),epsovl=1d-8
  logical:: lprint=.true.,savez=.false.,getz=.false. !dummy
  integer:: ifig=-999       !dummy
  integer:: ndatx,ifsy1,ifsy2,ifile_handle,ifsy,iix(36)
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
       open(newunit=ifsy1,file='band_lmfham_spin1.dat')
       !if(nspx==2) ifsy2 = ifile_handle()
       if(nspx==2) open(newunit=ifsy2,file='band_lmfham_spin2.dat')
       write(6,*)  'ndat =',ndat
    endif
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  ndimMTO=8
    !  iix(1:4)=[1,2,3,4]+9 !,6,7,8,9]
    !  iix(4+1:4+4)=[1,2,3,4]+18+9 !,5,6,7,8,9]+9+9
    !  iix(18+1:18+9)=[1,2,3,4,5,6,7,8,9]+18
    !  iix(27+1:27+9)=[1,2,3,4,5,6,7,8,9]+27

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
                ib1 = ib_table(ix(i))!,ldim) !atomic-site index in the primitive cell
                ib2 = ib_table(ix(j))!,ldim)
                !              ib1 = ib_table(ix(iix(i)))!,ldim) !atomic-site index in the primitive cell
                !              ib2 = ib_table(ix(iix(j)))!,ldim)
                do it =1,npair(ib1,ib2)
                   phase=1d0/nqwgt(it,ib1,ib2)*exp(-img*2d0*pi* sum(qp*matmul(plat,nlat(:,it,ib1,ib2))))
                   !                 hamm(i,j)= hamm(i,j)+ hammr(iix(i),iix(j),it,jsp)*phase
                   !                 ovlm(i,j)= ovlm(i,j)+ ovlmr(iix(i),iix(j),it,jsp)*phase
                   hamm(i,j)= hamm(i,j)+ hammr(i,j,it,jsp)*phase
                   ovlm(i,j)= ovlm(i,j)+ ovlmr(i,j,it,jsp)*phase
                enddo
             enddo
          enddo
          call zhev_tk4(ndimMTO,hamm,ovlm,nmx,nev, evl,t_zv, epsovl)!Diangonale (hamm- evl ovlm) z=0
          !do i=1,12
          !   write(6,"(i2,f9.3,2x,9f5.1,2x,9f5.1)")i,evl(i),(abs(t_zv(j,i)),j=1,ndimMTO/2)
          !enddo   
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
    write(6,"(a)")'!!! OK! band_lmfham_spin*.dat has generated! MTO only band plot'
    write(6,"(a)")'See README_MATERIALS.org for how to make a plot for band_lmfham_spin*.dat'
    write(6,"(a)")'  For example, gnuplot scrpt can be'
    write(6,"(a)")' plot \'
    write(6,"(a)")' "bnd001.spin1" u ($2):($3) lt 1 pt 1 w lp,\'
    write(6,"(a)")' "bnd002.spin1" u ($2):($3) lt 1 pt 1 w lp,\'
    write(6,"(a)")' "bnd003.spin1" u ($2):($3) lt 1 pt 1 w lp,\'
    write(6,"(a)")' "bnd004.spin1" u ($2):($3) lt 1 pt 1 w lp,\'
    write(6,"(a)")' "bnd005.spin1" u ($2):($3) lt 1 pt 1 w lp,\'
    write(6,"(a)")' "bnd006.spin1" u ($2):($3) lt 1 pt 1 w lp,\'
    write(6,"(a)")' "band_lmfham_spin1.dat" u ($1):(13.605*($2-ef)) pt 2'
    write(6,"(a)")' pause -1'
    deallocate(ovlm,hamm)
  end block GetEigenvaluesForSYML
endprogram


subroutine Hreduction(ndimPMT,hamm,ovlm,ndimMTO,ix,fff1, evl,hammout,ovlmout)
  !       call Hreduction(ndimPMT,hamm,ovlm, ndimMTO,ix,fff1, hamm,ovlm)
  !       firstreduction:block
  use m_readqplist,only: eferm
  use m_HamPMT,only: GramSchmidt
  implicit none
  integer::i,j,ndimPMT,ndimMTO,nx,nmx,ix(ndimMTO),nev
  real(8)::beta,emu,val,wgt(ndimPMT),evlmto(ndimMTO),evl(ndimPMT),evlx(ndimPMT)
  complex(8):: oz(ndimPMT,ndimPMT),wnj(ndimPMT,ndimMTO),wnm(ndimPMT,ndimMTO),wnn(ndimMTO,ndimMTO)
  complex(8):: evecmto(ndimMTO,ndimMTO),evecpmt(ndimPMT,ndimPMT)
  complex(8):: ovlmx(ndimPMT,ndimPMT),hammx(ndimPMT,ndimPMT),fac(ndimPMT,ndimMTO),ddd(ndimMTO,ndimMTO)
  complex(8):: hamm(ndimPMT,ndimPMT),ovlm(ndimPMT,ndimPMT)
  complex(8):: hammout(ndimMTO,ndimMTO),ovlmout(ndimMTO,ndimMTO)
  real(8):: epsovl=1d-8,fff1,fff
  ovlmx= ovlm
  hammx= hamm
  nmx = ndimMTO
  call zhev_tk4(ndimMTO,hamm(ix(1:ndimMTO),ix(1:ndimMTO)),ovlm(ix(1:ndimMTO),ix(1:ndimMTO)), &
       nmx,nev, evlmto, evecmto, epsovl) !MTO
  ovlm= ovlmx
  hamm= hammx
  nmx = ndimPMT
  call zhev_tk4(ndimPMT,hamm(1:ndimPMT,1:ndimPMT),ovlm(1:ndimPMT,1:ndimPMT), &
       nmx,nev, evl,    evecpmt, epsovl) !PMT
  ovlm=ovlmx
  do j=1,ndimMTO !wnm is corrected matrix element of <psi_PMT|psi_MTO>
     do i=1,ndimPMT
        fac(i,j)= sum(dconjg(evecpmt(:,i))*matmul(ovlmx(:,ix(1:ndimMTO)),evecmto(1:ndimMTO,j)))
        !fff= 1.5d0/(exp( (evl(i)-(ef+5d0))/1d0 )+ 1d0)
        fff= 1d0/(exp( (evl(i)-(eferm+5d0))/1d0 )+ 1d0)
        wnm(i,j)= fac(i,j) *abs(fac(i,j)) **fff1 !**fff ! Without abs(fac(i,j)), simple projection.
     enddo
  enddo
  call GramSchmidt(ndimPMT,ndimMTO,wnm)
  ! do i=1,ndimPMT
  !    do j=1,ndimMTO !wnm is corrected matrix element of <psi_PMT|psi_MTO>
  !       if(abs(wnm(i,j))>.1) write(6,*)'wnm matrix ',i,j,abs(wnm(i,j))
  !    enddo
  ! enddo
  ! Mapping operator wnm*<psi_MTO|F_i>, where F_i is MTO basis.
  wnj = matmul(wnm(1:ndimPMT,1:ndimMTO),matmul(transpose(dconjg(evecmto(:,:))),&
       ovlmx(ix(1:ndimMTO),ix(1:ndimMTO))))
  nx=ndimPMT          
  do i=1,ndimMTO
     do j=1,ndimMTO
        hammout(i,j)= sum( dconjg(wnj(1:nx,i))*evl(1:nx)*wnj(1:nx,j))
        ovlmout(i,j)= sum( dconjg(wnj(1:nx,i))*wnj(1:nx,j) ) !<MLO|MLO>
     enddo
  enddo
  !       endblock firstreduction
end subroutine Hreduction
