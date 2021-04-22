      subroutine gsortho(vecIO,L,S) !Gram-Schmid(GS) orthogonalization
      implicit none
      integer(4),intent(in) :: L,S !L=large dim. S=Small dim
      integer(4) i,j,n,m
      complex(8) ,intent(inout):: vecIO(S,L)
      real(8)     norm
      complex(8) prefc,dum(S,L)
      norm=0.0d0                !renormalization of the first vector
      do n=1,L
         norm = norm + dconjg(vecIO(1,n))*vecIO(1,n)
      end do
      norm=dsqrt(norm)      
      
      vecIO(1,:)=vecIO(1,:)/norm !normalization!         
      do m=2,S                  ! m vectors containing nevMTO components
         do i=1,m-1
            prefc=0.0d0
            do n=1,L
               prefc = prefc + dconjg(vecIO(i,n))*vecIO(m,n)
            end do                  
            vecIO(m,:)=vecIO(m,:)-prefc*vecIO(i,:)
         end do           
c     !renormalization!
         norm=0.0d0
         do n=1,L
            norm = norm + dconjg(vecIO(m,n))*vecIO(m,n)
         end do
         norm=dsqrt(norm)
         vecIO(m,:)=vecIO(m,:)/norm
c     !end renorm
      end do
      
c         write(*,*) "GS-ortho, "
c         do i=1,S
c            do n=1,L
c               if(i==n) write(*,"(f5.2,2x,f5.2,3x)",advance="no") dble(vecIO(i,n)),dimag(vecIO(i,n))              
c            end do
c            write(*,*) 
c         end do
         return
      end subroutine gsortho
!!-------------------------------------------------------------
      program lmfham
!! Read HamiltonianMTO and generates MTO-only Hamiltonian
      use m_hamMTO,only: ReadHamMTOInfo,plat,pos,nlat,npair,ib_table,nqwgt,lso
     &   ,nsp,nkp,npairmx,NMTO,startprint
      implicit none
      character(200) inname,oname1,oname2,mtoname
      logical:: lprint=.true.,savez=.false.,getz=.false.,skipdiagtest=.true.,explf=.false.
      integer:: ifqplistsy,i,j,n,m,s,t,l,access,ifoh1,ifoh2,idum,jsp,NM,NS,NP,ifdmat,Ndum,ifomto
      integer:: ifih,iq,nev,nevS,nevC,ifig=-999,nspx,lsodum,ilog,ifrmat,ifrmat2,di1,di2
      real(8)::qp(3),pi=4d0*atan(1d0),eferm,epsovl,el1,el2
      real(8),parameter ::  ryd=13.6058d0
      real(8),allocatable ::  evl_S(:),evl(:),evl_C(:)
      complex(8),allocatable,dimension(:,:) :: ham,hamS,ovl,ovlS,hamC,ovlC,rdmat,prjdum,aeadum
     .     ,amat,ovlrect,gsvec,inpro,dmat,c_S,c,aea,lft,rgt,prj,emat,rext
      complex(8),allocatable,dimension(:,:,:) :: rmat

      call startprint(0)

c      open(newunit=idum,file="input_lmfham")
c      read(idum,*) inname,mtoname,oname1,oname2
c      close(idum)
      open(newunit=idum,file='qplist.dat')
      read(idum,*) eferm
      close(idum)
      
c      call ReadHamMTOInfo(inname)
      call ReadHamMTOInfo(1)    !initmlo.Info      
      NM=NMTO
      open(newunit=idum,file="ir_list")
      read(idum,*) NS
      close(idum)     

      nspx=nsp  !nspx=1: non-mag, nspx=2: mag
      if(lso==1) nspx=1 !lso==0: without SO
c      if(lso==1) NS=NS*2 !L.S mode      
c      open(newunit=ifih,file =trim(adjustl(inname))//".K",form='unformatted')
c      open(newunit=ifomto,file =trim(adjustl(mtoname))//".K",form='unformatted')      
c      open(newunit=ifoh1,file=trim(adjustl(oname1))//".K",form='unformatted')
c      open(newunit=ifoh2,file=trim(adjustl(oname2))//".K",form='unformatted')
      open(newunit=ifih,file ="HamR.K",form='unformatted')
      open(newunit=ifomto,file ="HamM.K",form='unformatted')      
      open(newunit=ifoh1,file="HamS.K",form='unformatted')
      open(newunit=ifoh2,file="HamC.K",form='unformatted')
      open(newunit=ilog,file="sizelog")
      
      allocate(hamS(1:NS,1:NS),ovlS(1:NS,1:NS),hamC(1:NS,1:NS),c_S(1:NS,1:NS),evl_S(1:NS),ovlC(1:NS,1:NS),rdmat(1:NM,1:NM))
c     open(newunit=ifrmat,file="rotmatbin1",form='unformatted')
      open(newunit=ifrmat,file="rotmat1")
      do i=1,NM         
         do j=1,NM
            read(ifrmat,*) di1,di2,el1,el2
            rdmat(di1,di2)=dcmplx(el1,el2)
         end do
      end do            
c     read(ifrmat) rmat(1:NM,1:NM,1)
c     if(nsp==2)open(newunit=ifrmat2,file="rotmatbin2",form='unformatted')
c     if(nsp==2)read(ifrmat2) rmat(1:NM,1:NM,2)
         
      iq=0      
      if(nspx==1)write(6,"(a)") "**************** Reading Ham for iq(non-magnetic)***************"
      if(nspx==2)write(6,"(a)") "**************** Reading Ham for iq(magnetic)***************"
      write(6,"(a,i5)") "proceeded iq are"
      do 
         read(ifih,end=2019) qp,NP,lsodum,epsovl,jsp         
         if(lsodum.ne.lso) stop "inconsistent lso"
         write(ilog,*) "size of Hamiltonian (Large,Small)=",NS,NM,NP
!     jsp=isp in the collinear case; jsp=1 in the noncollinear
         if(jsp==1)then
            iq=iq+1
            write(6,"(i4)",advance="no") iq
            if(mod(iq,10)==0) write(6,"(a,i4,a)") "  : ",idnint((dble(iq)/dble(nkp)*1d2))," %"
         end if
!     making Ham_MTO
         allocate(ovl(1:NP,1:NP),ham(1:NP,1:NP),rext(1:NP,1:NP),dmat(1:NP,1:NP))
         rext=0d0
         do i=1,NP
            rext(i,i)=1d0
         end do
         rext(1:NM,1:NM)=rdmat(1:NM,1:NM)                           
         read(ifih) ovl(1:NP,1:NP) !ovl(1:NP,1:NP) !<-- NP respects the size of the original file                 
         dmat=matmul(transpose(rext),matmul(ovl,rext))         
         ovlS(1:NS,1:NS)=dmat(1:NS,1:NS)
         read(ifih) ham(1:NP,1:NP) !ham(1:NP,1:NP)
         dmat=matmul(transpose(rext),matmul(ham,rext))
         hamS(1:NS,1:NS)=dmat(1:NS,1:NS)         
         deallocate(rext,dmat)         
!     Diagonalize
         allocate(evl(1:NP),c(1:NP,1:NP))
         
         if(epsovl<1.000001d-14) then !epsovl is the trancation to remove poor linear-dependency basis
            call zhev_tk2(NP,ham,ovl,NP,nev,evl,c,lprint,savez,getz,ifig)
            call zhev_tk2(NS,hamS,ovlS,NS,nevS,evl_S,c_S,lprint,savez,getz,ifig) 
         else                    
            call zhev_tk3(NP,ham,ovl,NP,nev,evl,c,lprint,savez,getz,ifig,epsovl)
            call zhev_tk3(NS,hamS,ovlS,NS,nevS,evl_S,c_S,lprint,savez,getz,ifig,epsovl)
         endif
         
!!!!!!!making modified-Hamiltonian taking APW effect into account
         allocate(prj(1:nevS,1:nev),inpro(1:nevS,1:nev))
         inpro=transpose(matmul(transpose(dconjg(c)),matmul(ovl(1:NP,1:NS),c_S)))
         
         do j=1,nevS
            do s=1,nev
               prj(j,s)=inpro(j,s) *abs(inpro(j,s))**4.0d0
            end do
         end do        
         
         call gsortho(prj,nev,nevS) !Gram-Schmid orthogonalization         
!!!!!!!!!!!!!!!!!!!!
         allocate(emat(1:nev,1:nev),lft(1:NS,1:NS),rgt(1:NS,1:NS),aea(1:NS,1:NS))
         emat=0d0
         do n=1,nev
            emat(n,n)=evl(n)
         end do         
         aea=matmul(matmul(dconjg(prj),emat),transpose(prj))

         lft   = matmul( ovlS(1:NS,1:NS),c_S)
         rgt   = transpose(dconjg(lft))
         hamC  = matmul( lft,matmul(aea,rgt))
c     ovlC  = ovlS
         ovlC  = matmul( lft,rgt)                           
!!!!!!!!!!!!!! output hamiltonian data !!!!!!!!!!         
         write(ifoh1) qp,NS,lso,epsovl,jsp 
         write(ifoh1) ovlS(1:NS,1:NS)
         write(ifoh1) hamS(1:NS,1:NS) 
         write(ifoh2) qp,NS,lso,epsovl,jsp 
         write(ifoh2) ovlC(1:NS,1:NS)
         write(ifoh2) hamC(1:NS,1:NS)
         write(ifomto) qp,NM,lso,epsovl,jsp 
         write(ifomto) ovl(1:NM,1:NM)
         write(ifomto) ham(1:NM,1:NM)
ccccccccccccend of reduction step cccccccccccccccccccc
         deallocate(lft,rgt,emat,aea,inpro,prj,ovl,ham,c,evl)
      enddo
 2019 if(iq.ne.nkp) stop "error! something wrong between iq and nkp"
      write(6,*) " "
      write(6,*) "**************** end Reading Ham for iq **********"    
      stop "exiting lmfham  (return 0)"
      end program lmfham
