!!---  qplist.dat reader -----------------------------------
      module m_readqplist     
      integer,protected:: ndat
      real(8),allocatable,protected:: xdat(:),qplistsy(:,:)
      real(8),protected:: eferm !<-- temporary use. just as a memo.
      contains
      subroutine readqplistsy()
      implicit none
      integer:: ifqplistsy,nnn,ix
      open(newunit=ifqplistsy,file='qplist.dat')
      nnn=10000
      if(allocated(xdat)) deallocate(xdat)
      if(allocated(qplistsy)) deallocate(qplistsy)
      allocate(xdat(nnn),qplistsy(1:3,nnn))
      read(ifqplistsy,*) eferm
      ix=0
      do
        ix=ix+1
        read(ifqplistsy,*,end=1011) xdat(ix),qplistsy(1:3,ix)
c        write(6,'(" qplist ix xdata q=",i5,f9.4,x,3f9.4)')ix,xdat(ix),qplistsy(:,ix)
      enddo
 1011 continue
      ndat=ix-1
      close(ifqplistsy)
      return
      end subroutine
      end module m_readqplist      
!! -----------------------------MAIN ROUTINES---------------------------------
      program plotmodel
      use m_HamMTO,only: plat,npair,nlat,nqwgt,nkp,qplist,ib_table,alat,ReadHamMTOInfo,npairmx,nsp,
     . hamr,ovlr,startprint,mkRsMTO,epsovl,NS,ib,resetMTOInfo
      use m_readqplist,only: eferm,qplistsy,ndat,xdat, Readqplistsy
      implicit none
      integer:: i,j,n,ikp,ib1,ib2,it,nev,jsp,access,ii,jj,idum
      complex(8)::img=(0d0,1d0),phase
      complex(8),dimension(:,:)  ,allocatable:: ham,ovl,t_zv,t_zv2
      real(8)   ,dimension(:,:),allocatable:: evl
      real(8)::qp(3),pi=4d0*atan(1d0),ryd,iprd
      logical:: debug=.true.,lprint=.true.,savez=.false.,getz=.false.
      integer:: ifig=-999,ifi,ifbnd 
      character(200) fname
      ryd=13.605693d0
      call startprint(3)      
!     ! Reading infomation for Hamiltonian (lattice structures and index of basis).
      call readqplistsy()
            
      open(newunit=ifi,file="input_plotmodel")
      call ReadHamMTOInfo(1)
      do
         read(ifi,*,end=1919) fname,jsp !jsp=1(up) =2(down)                  
         write(*,"(a23,a20)") "Hamiltonian's name is ",trim(adjustl(fname))
         if(jsp==1)open(newunit=ifbnd,file=trim(adjustl(fname))//".bnd.spin1")
         if(jsp==2)open(newunit=ifbnd,file=trim(adjustl(fname))//".bnd.spin2")
         if(jsp==1.and.debug)open(newunit=idum,file=trim(adjustl(fname))//".debug")
         call mkRsMTO(fname,jsp) !making Real-space Hamiltonian
         
         if(allocated(ovl))  deallocate(ovl)
         if(allocated(ham))  deallocate(ham)
         if(allocated(t_zv)) deallocate(t_zv)
         if(allocated(evl))  deallocate(evl)
         allocate(ovl(1:NS,1:NS),ham(1:NS,1:NS),t_zv(NS,NS),evl(ndat,NS))
         write(6,*) " "
         write(6,*) "*************** ikp along qplist *************************"
         
         do ikp=1,ndat
            qp(1:3)= qplistsy(1:3,ikp)
            write(6,"(i4)",advance="no") ikp
            if(mod(ikp,10)==0) write(6,"(a,i4,a)") "  : ",idnint((dble(ikp)/dble(ndat)*1d2))," %"
            ham  = 0d0
            ovl  = 0d0
            do i=1,NS
               do j=1,NS
                  do it =1,npair(ib(i),ib(j))
                     iprd=sum(qp*matmul(plat,nlat(:,it,ib(i),ib(j))))
                     phase=1d0/nqwgt(it,ib(i),ib(j))*exp(-img*2d0*pi*iprd)
                     ham(i,j)    = ham(i,j)    + hamr(i,j,it)*phase
                     ovl(i,j)    = ovl(i,j)    + ovlr(i,j,it)*phase
                  end do
               end do
            end do

            if(epsovl<1.000001d-14) then
               call zhev_tk2(NS,ham,ovl,NS,nev,evl(ikp,1:NS),t_zv,lprint,savez,getz,ifig) 
            else  
               call zhev_tk3(NS,ham,ovl,NS,nev,evl(ikp,1:NS),t_zv,lprint,savez,getz,ifig,epsovl) 
            endif
            if(jsp==1.and.debug)then
               do n=1,nev
                  write(idum, "(2i4,f13.7)") i,jsp,evl(ikp,n)
c                  write(ifbnd,   "(f15.5,f15.5,2i4)")   xdat(ikp),(evl(ikp,n)-eferm)*ryd,jsp,n
c                  write(ifbnd,*) " "
               end do
               write(idum,*) " "
            end if
         enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         write(6,*) " "
         write(6,*) "**********************************************************"      
         
         do n=1,nev
            do ikp=1,ndat
               write(ifbnd,   "(f15.5,f15.5,2i4)")   xdat(ikp),(evl(ikp,n)-eferm)*ryd,jsp,n
            enddo               
            write(ifbnd,*) " "
         enddo
         close(ifbnd)
         if(jsp==1.and.debug) close(idum)
         deallocate(ham,ovl,t_zv,evl)
      end do      
 1919 continue
cccccccccccccccccccccccc finalization ccccccccccccccccccccccccccccccc
      write(6,*) " "
      stop "program plotmodel exit(return 0)"
      end program plotmodel

