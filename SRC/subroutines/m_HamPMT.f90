module m_HamPMT ! -- Read HamiltionanPMTinfo and HamiltonianPMT. Then convert HamPMT to HamRsMTO  ------
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
    integer:: ififft,i,lold,m
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
    use m_zhev,only:zhev_tk4
    use m_ftox
    use m_readqplist,only: eferm
    implicit none
    integer:: ifihmto
    integer::ikpd,ikp,ib1,ib2,ifih,it,iq,nev,nmx,ifig=-999,i,j,ndimPMT,lold,m
    complex(8),allocatable:: hamm(:,:),ovlm(:,:)!,t_zv(:,:)!,ovlmx(:,:)
    logical:: lprint=.true.,savez=.false.,getz=.false.,skipdiagtest=.false.
!    real(8),allocatable:: evl(:)
    complex(8):: img=(0d0,1d0),aaaa,phase
    real(8)::qp(3),pi=4d0*atan(1d0),fff,ef,fff1=2,fff2=2,fff3=0 
    integer:: ix(ldim),ix1(ldim),ix2(ldim),ix3(ldim),ix12(ldim),ix23(ldim),nnn,ib,k,l,ix5,imin,ixx &
         ,ndimMTO2,j2,j1,j3,ndimMTO3
    integer:: nx

!    [fff1,fff2,fff3]=[8d0,8d0,8d0]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1    
!     nnn=0
!     do i=1,ldim !MLO3 dimension
!        if(l_table(i)>=3) cycle 
!        if(k_table(i)==2) cycle !.and.l_table(i)>=1) cycle
! !       if( (ib_table(i)==3.or.ib_table(i)==4).and.l_table(i)>=2) cycle
! !       if( (ib_table(i)==3.or.ib_table(i)==4).and.l_table(i)==0) cycle
! !       if(k_table(i)==2.and.l_table(i)>=2) cycle
!        write(6,*) 'ham3 index', i,ib_table(i),l_table(i),k_table(i)
!        nnn=nnn+1
!        ix3(nnn)=i
!     enddo
!     ndimMTO3=nnn
!    

    ! nnn=0
    !  do i=1,ldim !MLO2 dimension
    !     if(l_table(i)>=3) cycle 
    !     if(k_table(i)==2) cycle 
    !     nnn=nnn+1
    !     ix2(nnn)=i
    !  enddo
    !  ndimMTO2=nnn
!
    nnn=0
    do i=1,ldim 
       if(l_table(i)>=3) cycle !only spd
       !if(k_table(i)==2) cycle !.and.l_table(i)>=2) cycle ! throw away EH2
       !if(k_table(i)==2) cycle ! throw away EH2

       
       ! for Si LDA 9=(
       ! for Si QSGW 13+13
       !MLO dimension !for NiO LDA 9+9+4+4
       !if( (ib_table(i)==3.or.ib_table(i)==4).and.l_table(i)>=2) cycle !for Oxygen of NiO. skip 3d
       write(6,*) 'ham1 index', i,ib_table(i),l_table(i),k_table(i)
       nnn=nnn+1
       ix1(nnn)=i
    enddo
    ndimMTO=nnn
!
    do j2=1,ndimMTO2 !Map from MLO2->MLO
        do j1=1,ndimMTO
           if(ix2(j2)==ix1(j1)) ix12(j2)=j1  
        enddo
    enddo
!    do j3=1,ndimMTO3 !Map from MLO3->MLO2
!       do j2=1,ndimMTO2
!          if(ix3(j3)==ix2(j2)) ix23(j3)=j2  
!       enddo
!    enddo
    
    open(newunit=ifih,file='HamiltonianPMT',form='unformatted')
    write(6,*)'Reaing HamiltonianPMT ...'
    if(lso==1) ldim=ldim*2 !L.S mode
    nspx=nsp
    if(lso==1) nspx=1
    write(6,"('ndimMTO ldim lso=',i6,4i3)") ndimMTO,ldim,lso
    
    allocate(ovlmr(1:ndimMTO,1:ndimMTO,npairmx,nspx), hammr(1:ndimMTO,1:ndimMTO,npairmx,nspx))
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
       epsovl=0d0 !1d-8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
       ndimMTO3=ndimMTO
       ix3(1:ndimMTO3)=ix1(1:ndimMTO)
!       ndimMTO3=ndimMTO2
!       ix3(1:ndimMTO2)=ix2(1:ndimMTO2)
       eigen: block
         real(8):: evlmlo(ndimMTO),evl(ndimPMT),evlmlo2(ndimMTO2),evlmlo3(ndimMTO3)
         complex(8)::hh(ndimMTO3,ndimMTO3),oo(ndimMTO3,ndimMTO3)
         ! <PsiPMT|PsiMLO> 
         call Hreduction(.false.,ndimPMT,hamm(1:ndimPMT,1:ndimPMT),ovlm(1:ndimPMT,1:ndimPMT), &
              ndimMTO,ix1,fff1,   evl,hamm(1:ndimMTO,1:ndimMTO),ovlm(1:ndimMTO,1:ndimMTO))
!
! <PsiMLO|PsiMLO2> 
!         call Hreduction(.false.,ndimMTO,hamm(1:ndimMTO,1:ndimMTO),ovlm(1:ndimMTO,1:ndimMTO),&
!              ndimMTO2,ix12,fff2, evlmlo,hamm(1:ndimMTO2,1:ndimMTO2),ovlm(1:ndimMTO2,1:ndimMTO2))
! <PsiMLO2|PsiMLO3> 
!         call Hreduction(.false.,ndimMTO2,hamm(1:ndimMTO2,1:ndimMTO2),ovlm(1:ndimMTO2,1:ndimMTO2),&
         !              ndimMTO3,ix23,fff3, evlmlo2,hamm(1:ndimMTO3,1:ndimMTO3),ovlm(1:ndimMTO3,1:ndimMTO3))
!         print *,'goto finaleigen'
         finaleigen: block
           complex(8):: evec(ndimMTO3**2), hh(ndimMTO3,ndimMTO3),oo(ndimMTO3,ndimMTO3)
           if(sum([qp(2),qp(3)]**2)<1d-3) then
              hh=hamm(1:ndimMTO3,1:ndimMTO3)
              oo=ovlm(1:ndimMTO3,1:ndimMTO3)
              nmx = ndimMTO3
              call zhev_tk4(ndimMTO3,hh,oo, nmx,nev,evlmlo3, evec, epsovl)
              !do i=1,ndimMTO3 
                 !if(i<=4.and.abs(evlmlo3(i)-evl(i))>1d-8) call rx('eeeeeeeee dif')
                 ! write(6,"(a,i5,4f12.4)")'mmm', i,evl(i),evlmlo2(i)-evl(i),evlmlo3(i)-evl(i)
                 !if(sum([qp(2),qp(3)]**2)<1d-3)  write(6,"(a,i5,4f12.4)")'mmmx', i,evl(i),evlmlo3(i)-evl(i)
              !   write(6,"(a,i5,4f12.4)")'mmm', i,evl(i),evlmlo3(i)-evl(i)
              !enddo
           endif
         endblock finaleigen
       endblock eigen
       !! Real space Hamiltonian. H(k) ->  H(T) FourierTransformation to real space
       do i=1,ndimMTO3
          do j=1,ndimMTO3
             ib1 = ib_table(ix3(i)) 
             ib2 = ib_table(ix3(j)) 
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
    !! write RealSpace MTO Hamiltonian
    ndimMTO=ndimMTO3
    ix(1:ndimMTO)=ix3(1:ndimMTO) !for atom idex
    write(6,*)' Writing HamRsMTO... ndimMTO=',ndimMTO
    open(newunit=ifihmto,file='HamRsMTO',form='unformatted')
    write(ifihmto) ndimMTO,npairmx,nspx
    write(ifihmto) hammr(1:ndimMTO,1:ndimMTO,1:npairmx,1:nspx),ovlmr(1:ndimMTO,1:ndimMTO,1:npairmx,1:nspx),ix(1:ndimMTO)
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
       zmel(:,it) = vec/sum(dconjg(vec)*vec)**.5d0
    enddo
  end subroutine GramSchmidt
end module m_HamPMT
module m_HamRsMTO
  !! read real-space MTO Hamiltonian
  integer,protected:: ndimMTO,npairmx,nspx
  integer,allocatable,protected:: ix(:)
  complex(8),allocatable,protected:: ovlmr(:,:,:,:),hammr(:,:,:,:)
contains
  !! read RealSpace MTO Hamiltonian
  subroutine ReadHamRsMTO()
    integer:: ifihmto
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
subroutine Hreduction(iprx,ndimPMT,hamm,ovlm,ndimMTO,ix,fff1, evl,hammout,ovlmout)!Reduce H(ndimPMT) to H(ndimMTO)
  use m_zhev,only:zhev_tk4
  use m_readqplist,only: eferm
  use m_HamPMT,only: GramSchmidt
  use m_ftox
  implicit none
  integer::i,j,ndimPMT,ndimMTO,nx,nmx,ix(ndimMTO),nev,nxx,jj
  real(8)::beta,emu,val,wgt(ndimPMT),evlmto(ndimMTO),evl(ndimPMT),evlx(ndimPMT)
  complex(8):: oz(ndimPMT,ndimPMT),wnj(ndimPMT,ndimMTO),wnm(ndimPMT,ndimMTO),wnn(ndimMTO,ndimMTO)
  complex(8):: evecmto(ndimMTO,ndimMTO),evecpmt(ndimPMT,ndimPMT)
  complex(8):: ovlmx(ndimPMT,ndimPMT),hammx(ndimPMT,ndimPMT),fac(ndimPMT,ndimMTO),ddd(ndimMTO,ndimMTO)
  complex(8):: hamm(ndimPMT,ndimPMT),ovlm(ndimPMT,ndimPMT)
  complex(8):: hammout(ndimMTO,ndimMTO),ovlmout(ndimMTO,ndimMTO)
  real(8):: epsovl=0d0 ,fff1,fff !epsovl=1d-8
  logical:: iprx
  ovlmx= ovlm
  hammx= hamm
  nmx = ndimMTO
!  write(6,*)'Hreduction: eferm=',eferm
  call zhev_tk4(ndimMTO,hamm(ix(1:ndimMTO),ix(1:ndimMTO)),ovlm(ix(1:ndimMTO),ix(1:ndimMTO)), &
       nmx,nev, evlmto, evecmto, epsovl) !MTO
!  write(6,ftox)'evecmto=',ftof(evlmto)
  ovlm= ovlmx
  hamm= hammx
  nmx = ndimPMT
  call zhev_tk4(ndimPMT,hamm(1:ndimPMT,1:ndimPMT),ovlm(1:ndimPMT,1:ndimPMT), &
       nmx,nev, evl,    evecpmt, epsovl) !PMT
  ovlm=ovlmx
!  write(6,ftox)'evevl=',ftof(evl)
  
  do j=1,ndimMTO !wnm is corrected matrix element of <psi_PMT|psi_MTO>
     do i=1,ndimPMT
        fac(i,j)= sum(dconjg(evecpmt(:,i))*matmul(ovlmx(:,ix(1:ndimMTO)),evecmto(1:ndimMTO,j)))
!        if(fff1>=0d0) then
!           fff=fff1
!        else
!           fff= 4d0
           !fff= 0.5d0 !1d0 !2d0/(exp( (evl(i)-(eferm+5d0))/1d0 )+ 1d0)
           !fff= 1d0/(exp( (evl(i)-(eferm+5d0))*2d0 )+ 1d0)
           !if(evl(i)< eferm+2d0) fff=3d0
!        endif
!        wnm(i,j)= fac(i,j) *abs(fac(i,j)) **fff  !1 !**fff ! Without abs(fac(i,j)), simple projection.
     enddo
  enddo
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  mulfac :block
    integer:: ie,nidxevlmto,nidxevl,ibx,jx,idxevlmto(ndimMTO),idxevl(ndimPMT),jbx
    real(8):: eee,fffx,ecut,www,xxx
    real(8),allocatable::mulfac(:,:),mulfacw(:,:)
    
    ! ie=0
    ! eee=-1d10
    ! do j=1,ndimMTO
    !    if(evlmto(j)-eee>1d-3) then
    !       ie=ie+1
    !       eee=evlmto(j)
    !    endif   
    !    idxevlmto(j)=ie !degeneracy index
    ! enddo
    ! ie=0
    ! eee=-1d10
    ! do i=1,ndimPMT
    !    if(evl(i)-eee>1d-3) then
    !       ie=ie+1
    !       eee=evl(i)
    !    endif   
    !    idxevl(i)=ie
    ! enddo
    ! nidxevlmto = idxevlmto(ndimMTO) !number of bands degenerated
    ! nidxevl    = idxevl(ndimPMT)
    ! allocate( mulfac(nidxevl,nidxevlmto),mulfacw(nidxevl,nidxevlmto))
    ! if(fff1>=0d0) then
    !    fff= fff1
    ! else
    !    fff= 1d0
    ! endif
    ! ecut=1d0
    ! mulfac=0d0
    ! mulfacw=0d0
    ! do j=1,ndimMTO 
    !    do i=1,ndimPMT
    !       if(evlmto(j)<eferm+ecut) then
    !          www = abs(fac(i,j))**2
    !       else
    !          www=1d0
    !       endif
    !       ibx= idxevl(i)
    !       jbx= idxevlmto(j)
    !       mulfac (ibx,jbx) = mulfac (ibx,jbx)+ www
    !       mulfacw(ibx,jbx) = mulfacw(ibx,jbx)+ 1d0
    !    enddo
    ! enddo
    ! mulfac=mulfac/mulfacw

    
    do j=1,ndimMTO 
       do i=1,ndimPMT
          !          wnm(i,j) = fac(i,j) * mulfac(idxevl(i),idxevlmto(j))
          
          !wnm(i,j) = fac(i,j)
          !if(evlmto(j)<eferm+ecut) wnm(i,j) = fac(i,j) *abs(fac(i,j))
          
          !ecut=1d0
          !eee= evlmto(j)-(eferm+ecut)
          !xxx= 4d0/(exp(eee*100d0)+1d0)
          
          !xxx=0d0
          !if(eee<0d0) xxx=1d0
          !www = abs(fac(i,j))**xxx

          !ecut=1d0
          !www=1d0
          !if(evlmto(j)<eferm+ecut)
          www=abs(fac(i,j))
          wnm(i,j) = fac(i,j) * www
       enddo
    enddo   
  endblock mulfac

  !wnm=fac
  
  ! if(fff1<0) then
  !    do j=1,8 !ndimMTO !wnm is corrected matrix element of <psi_PMT|psi_MTO>
  !       do i=1,ndimPMT !occupied bands
  !          if(abs(wnm(i,j))>0.5) then
  !             wnm(i,j)=wnm(i,j)/abs(wnm(i,j))
  !          else   
  !             wnm(i,j)=0d0
  !          endif
  !          if(abs(wnm(i,j))>0.5) write(6,*)'wwwwnm=',i,j,wnm(i,j)
  !       enddo
  !    enddo
  ! endif
  
  !if(fff1<0) then
  !do i=1,20
  !   write(6,ftox)i,(ftof(abs(wnm(i,jj)),2),jj=1,20)
  !enddo
  !endif

  call GramSchmidt(ndimPMT,ndimMTO,wnm)
  
  if(iprx) then
     do j=1,ndimMTO !wnm is corrected matrix element of <psi_PMT|psi_MTO>
        do i=1,ndimPMT
           if(abs(wnm(i,j))>.1) write(6,*)'wnm matrix ',j,i,abs(wnm(i,j))**2
        enddo
     enddo
  endif
  ! Mapping operator wnm*<psi_MTO|F_i>, where F_i is MTO basis.
  wnj = matmul(wnm(1:ndimPMT,1:ndimMTO),matmul(transpose(dconjg(evecmto(:,:))),ovlmx(ix(1:ndimMTO),ix(1:ndimMTO))))
  nx=ndimPMT
!  if(fff1<0) then
!     nxx=8
!     evl(1:nxx)=evl(1:nxx)
!     evl(nxx+1:nx)=eferm+1.5d0
  !  endif
  
  !do i=1,nx
  !   if(evl(i)>eferm+2d0) evl(i)=eferm+10d0 !evl(ndimMTO)
  !enddo
  !nx=ndimMTO
  do i=1,ndimMTO
     do j=1,ndimMTO
!        hammout(i,j)= sum( dconjg(wnj(1:nx,i))*evl(1:nx)*wnj(1:nx,j))
        hammout(i,j)= sum( dconjg(wnj(1:nx,i))*evl(1:nx)*wnj(1:nx,j))
        ovlmout(i,j)= sum( dconjg(wnj(1:nx,i))*wnj(1:nx,j) ) !<MLO|MLO>
!        hammout(i,j)= hammx(ix(i),ix(j))
!        ovlmout(i,j)= ovlmx(ix(i),ix(j))
     enddo
  enddo
end subroutine Hreduction

