program lmfham
  !! Read HamiltonianMTO and generates MTO-only Hamiltonian
  use m_hamMTO,only: ReadHamMTOInfo,plat,pos,nlat,npair,ib_table,nqwgt,lso &
       ,nsp,nkp,npairmx,NMTO,startprint
  implicit none
  character(200) inname,oname1,oname2,mtoname
  character(8) tgs
  logical:: debug=.true.,lprint=.true.,savez=.false.,getz=.false.
  integer:: ifqplistsy,i,j,n,m,s,t,l,access,ifoh1,ifoh2,idum,jsp,NM,NS,NP,Ndum,ifomto
  integer:: ifih,iq,nevP,nevS,nevC,ifig=-999,nspx,lsodum,ilog,ifrmat,ifrmat2,di1,di2,iftxt,iftxt2
  real(8)::qp(3),pi=4d0*atan(1d0),eferm,epsovl,el1,el2,expn
  real(8),parameter ::  ryd=13.6058d0
  real(8),allocatable ::  evl_S(:),evl_P(:),evl_C(:)
  complex(8),allocatable,dimension(:,:) :: hamP,hamI,hamS,ovlP,ovlI,ovlS,hamC,ovlC,rdmat &
       ,c_S,c_C,C_P,inpro,dmat,aea,lft,rgt,prj,emat,rext

  call startprint(0)

  open(newunit=idum,file='qplist.dat')
  read(idum,*) eferm
  close(idum)

  open(newunit=idum,file="mtoctrl")
  do
     read(idum,"(a8)") tgs
     if(tgs=="<lmfham>") exit
  end do
  read(idum,*) expn
  close(idum)

  call ReadHamMTOInfo(1)    !initmlo.Info
  NM=NMTO
  open(newunit=idum,file="ir_list")
  read(idum,*) NS
  close(idum)

  nspx=nsp  !nspx=1: non-mag, nspx=2: mag
  if(lso==1) nspx=1 !lso==0: without SO
  !     if(lso==1) NS=NS*2 !L.S mode
  open(newunit=ifih,file ="HamI",form='unformatted') !input
  open(newunit=ifomto,file ="HamM",form='unformatted') !output1: full-MTO
  open(newunit=ifoh1,file="HamS",form='unformatted')   !output2: trancated without projection
  open(newunit=ifoh2,file="HamC",form='unformatted') !output3: MLO(projected)
  if(debug) open(newunit=idum,file="lmfham.debug")
  open(newunit=ilog,file="sizelog")

  allocate(hamS(1:NS,1:NS),ovlS(1:NS,1:NS),hamC(1:NS,1:NS),c_S(1:NS,1:NS),c_C(1:NS,1:NS))
  allocate(evl_S(1:NS),ovlC(1:NS,1:NS),rdmat(1:NM,1:NM),evl_C(1:NS))
  open(newunit=ifrmat,file="rotmat")
  do i=1,NM
     do j=1,NM
        read(ifrmat,*) di1,di2,el1,el2
        rdmat(di1,di2)=dcmplx(el1,el2)
     end do
  end do
  close(ifrmat)

  iq=0
  if(nspx==1)write(6,"(a)") "**************** Reading Ham for iq(non-magnetic)***************"
  if(nspx==2)write(6,"(a)") "**************** Reading Ham for iq(magnetic)***************"
  write(6,"(a,i5)") "proceeded iq are"
  do
     read(ifih,end=2019) qp,NP,lsodum,epsovl,jsp
     if(lsodum /= lso) stop "inconsistent lso"
     !     jsp=isp in the collinear case; jsp=1 in the noncollinear
     if(jsp==1)then
        iq=iq+1
        write(6,"(i4)",advance="no") iq
        if(mod(iq,10)==0) write(6,"(a,i4,a)") "  : ",idnint((dble(iq)/dble(nkp)*1d2))," %"
     end if
     !     making Ham_MTO

     allocate(rext(1:NP,1:NP),dmat(1:NP,1:NP),ovlI(1:NP,1:NP),hamI(1:NP,1:NP)) !dummy: will be deallocated by line 126
     rext=0d0
     do i=1,NP
        rext(i,i)=1d0
     end do
     rext(1:NM,1:NM)=rdmat(1:NM,1:NM)
     read(ifih) ovlI(1:NP,1:NP) !ovl(1:NP,1:NP) !<-- NP respects the size of the original file
     read(ifih) hamI(1:NP,1:NP) !ham(1:NP,1:NP)

     allocate(ovlP(1:NP,1:NP),hamP(1:NP,1:NP))
     dmat=matmul(transpose(rext),matmul(ovlI,rext))
     ovlS(1:NS,1:NS)=dmat(1:NS,1:NS)
     ovlP(1:NP,1:NP)=dmat(1:NP,1:NP)
     dmat=matmul(transpose(rext),matmul(hamI,rext))
     hamS(1:NS,1:NS)=dmat(1:NS,1:NS)
     hamP(1:NP,1:NP)=dmat(1:NP,1:NP)
     deallocate(rext,dmat,ovlI,hamI)
     !     Diagonalize
     allocate(evl_P(1:NP),c_P(1:NP,1:NP))

     call zhev_tk4(NP,hamP,ovlP,NP,nevP,evl_P,c_P,epsovl)
     call zhev_tk4(NS,hamS,ovlS,NS,nevS,evl_S,c_S,epsovl)
     if(nevS /= NS) stop "error! nevS"
     write(ilog,*) "NS,nevS,NM,NP,nevP",NS,nevS,NM,NP,nevP

!!!!!!!making modified-Hamiltonian taking APW effect into account
     allocate(prj(1:nevP,1:nevS),inpro(1:nevP,1:nevS))
     inpro=0d0
     inpro=dconjg(matmul(transpose(dconjg(c_P(1:NP,1:nevP))),matmul(ovlP(1:NP,1:NS),c_S)))

     prj=0d0
     do j=1,nevP
        do s=1,nevS
           prj(j,s)=inpro(j,s) *abs(inpro(j,s))**expn
           !             if(j==s)prj(j,s)=1d0
        end do
     end do

     call gsortho(prj,nevP,nevS) !Gram-Schmid orthogonalization
!!!!!!!!!!!!!!!!!!!!
     allocate(emat(1:nevP,1:nevP),lft(1:NS,1:nevS),rgt(1:nevS,1:NS),aea(1:NS,1:NS))
     emat=0d0
     do n=1,nevP
        emat(n,n)=evl_P(n)
     end do
     aea=matmul(matmul(transpose(dconjg(prj)),emat),prj)

     lft   = matmul( ovlS(1:NS,1:NS),c_S)
     rgt   = transpose(dconjg(lft))
     hamC  = matmul( lft,matmul(aea,rgt))
     ovlC  = ovlS
     !        ovlC  = matmul( lft,rgt)

     call zhev_tk4(NS,hamC,ovlC,NS,nevC,evl_C,c_C,epsovl)

     if(debug)then
        do i=1,nevC
           write(idum,"(2i4,4f13.7)") iq,jsp,evl_C(i),evl_S(i),evl_P(i),evl_C(i)-evl_P(i)
        end do
        write(idum,*) " "
     end if

!!!!!!!!!!!!!! output hamiltonian data !!!!!!!!!!
     write(ifoh1) qp,NS,lso,epsovl,jsp
     write(ifoh1) ovlS(1:NS,1:NS)
     write(ifoh1) hamS(1:NS,1:NS)
     write(ifoh2) qp,NS,lso,epsovl,jsp
     write(ifoh2) ovlC(1:NS,1:NS)
     write(ifoh2) hamC(1:NS,1:NS)
     write(ifomto) qp,NM,lso,epsovl,jsp !HamM
     write(ifomto) ovlP(1:NM,1:NM) !HamM excludes only PW bais
     write(ifomto) hamP(1:NM,1:NM)
     ! ccccccccccend of reduction step cccccccccccccccccccc
     deallocate(lft,rgt,emat,aea,inpro,prj,ovlP,hamP,c_P,evl_P)
  enddo
2019 if(iq /= nkp) stop "error! something wrong between iq and nkp"
  write(6,*) " "
  write(6,*) "**************** end Reading Ham for iq **********"
  stop "exiting lmfham  (return 0)"
end program lmfham
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

subroutine gsortho(vecIO,L,S) !Gram-Schmid(GS) orthogonalization
  implicit none
  logical(4) :: debug=.false.
  integer(4),intent(in) :: L,S !L=large dim. S=Small dim
  integer(4) i,j,n,m
  complex(8) ,intent(inout):: vecIO(L,S)
  real(8)     norm
  complex(8) prefc
  norm=0.0d0                !renormalization of the first vector
  do n=1,L
     norm = norm + dconjg(vecIO(n,1))*vecIO(n,1)
  end do
  norm=dsqrt(norm)

  vecIO(:,1)=vecIO(:,1)/norm !normalization!
  do m=2,S                  ! m vectors containing nevMTO components
     do i=1,m-1
        prefc=0.0d0
        do n=1,L
           prefc = prefc + dconjg(vecIO(n,i))*vecIO(n,m)
        end do
        vecIO(:,m)=vecIO(:,m)-prefc*vecIO(:,i)
     end do
     !     !renormalization!
     norm=0.0d0
     do n=1,L
        norm = norm + dconjg(vecIO(n,m))*vecIO(n,m)
     end do
     norm=dsqrt(norm)
     vecIO(:,m)=vecIO(:,m)/norm
     !     !end renorm
  end do

  if(debug)then
     write(*,*) "GS-ortho, "
     do i=1,S
        do n=1,L
           if(i==n) write(*,"(f5.2,2x,f5.2,3x)",advance="no") dble(vecIO(i,n)),dimag(vecIO(i,n))
        end do
        write(*,*)
     end do
  end if
  return
end subroutine gsortho
!!-------------------------------------------------------------
