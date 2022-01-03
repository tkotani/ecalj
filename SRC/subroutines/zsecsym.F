      subroutine zsecsym(zsec,ntq,nq,nband,nbandmx,nspinmx,nspin, eibzsym,ngrp,tiii,q,is)
!! --- symmetrize zsec for eibz4sig mode. -----------------
!! Read a file lmfgw_kdivider, which contains info for vxc and evec (they are in separated files in MPI)
!!
c      use m_mpi,only: MPI__AllreduceSum
      use m_readeigen,only: READEVAL
      implicit none 
      complex(8),intent(inout)::zsec(ntq,ntq,nq)
      integer,intent(in)::ntq,nq,nspinmx,nband,nbandmx(nq,nspinmx),is,nspin
      integer,intent(in):: ngrp,eibzsym(ngrp,-1:1,nq)
      logical,intent(in):: tiii !time reversal switch
      real(8),intent(in):: q(3,nq)

      complex(8),allocatable::zsect(:,:)
      integer:: ifile_handle,iqq
      integer:: procid,nrankv,ifvxc_,ifevec_,ifiproc,iqqxx,
     & isp,ixx,ixxx,nqixx,nspxx,ispxx,iqbz,i,igrp,iq
      character*256:: extn,ext
      character*256,allocatable:: extp(:)
      integer,allocatable:: ifevec__(:),ifvxc__(:),iprocq(:,:)

      integer:: nsym,nhdim,it,nblk,iband,napw,ldim,ierr,ispx,nbsize,nbsizemx
     & ,iblk1,iblk2,ii1,ii2,ie1,ie2,ne1,ne2,iqxx, ndimhx, nspx,nnnx
      integer,allocatable::iblki(:),iblke(:)
      complex(8),allocatable:: evec(:,:),evec_inv(:,:),evecrot(:,:),rmatjj(:,:,:)
      real(8),allocatable::evaliq(:)
      real(8)::tolry=1d-4,qqqx(3),qtarget(3),tolq=1d-8
      complex(8),allocatable:: ovl(:,:)
      integer::nev,j
      write(6,*)'zsecsym:'
      allocate( zsect(ntq,ntq))
!! === readin lmfgw_kdivider, and get extensions === apr2013
      ifiproc=ifile_handle()
      open(unit=ifiproc,file='lmfgw_kdivider',status='old')
      read(ifiproc,*) ext
      read(ifiproc,*) nqixx, nspxx, nrankv
      if(allocated(iprocq)) deallocate(iprocq)
      allocate(iprocq(nqixx,nspxx))
      do isp=1,nspxx
        do iqq=1,nqixx
          read(ifiproc,*) iqqxx, ispxx, ixxx
          if(iqqxx/=iqq) call rx( 'iqqxx/=iqq')
          if(ispxx/=isp) call rx( 'ispxx/=isp')
          iprocq(iqq,isp) = ixxx
          write(6,"('iqq isp irank=',i8,i2,i6)") iqq,isp, iprocq(iqq,isp)
        enddo
      enddo
      close(ifiproc)
!! for multiple files.
c      if(allocated(extp)) deallocate(extp,ifvxc__,ifevec__)
      allocate(extp(0:nrankv-1),ifvxc__(0:nrankv-1),ifevec__(0:nrankv-1))
      extp(0) = trim(ext)
      write(6,"('  0 ext= ',a,a)") trim(extp(0)),' ----------'
      do procid=1,nrankv-1
        write(extn,"(i10)") procid
        extp(procid)=trim(adjustl(ext))//'_'//trim(adjustl(extn))
        write(6,"(i3,' ext= ',a,a)") procid,trim(extp(procid)),' ----------'
      enddo
      do procid=0,nrankv-1
        ifvxc__(procid) = ifile_handle()
        open( ifvxc__(procid), file='vxc'//extp(procid),form='unformatted')
        ifevec__(procid)= ifile_handle()
        open( ifevec__(procid), file='evec'//extp(procid),form='unformatted')
      enddo
      ifvxc_ = ifvxc__(0)       !0 is root
      ifevec_= ifevec__(0)
      read(ifevec_)  ndimhx, nspx,nnnx 
      read(ifvxc_)              !skip ndimh, nsp,nnn
      allocate(evaliq(nband),iblki(nband),iblke(nband))
      iqq=0                     !iqq is to read multiple vxc.* evec.*
      do 3020 iq=1,nq           !nq means iq for which we will calculate sigma
        iqq=iqq+1
        do 3030 ispx=1,nspin !nspinmx  bugfix!!! at 2021-11-01 !ispx loop is to find isx=is
          ifvxc_  = ifvxc__ (iprocq(iqq,ispx)) 
          ifevec_ = ifevec__(iprocq(iqq,ispx)) 
          if(ispx==is) then 
!this if-block is due to evec and v_xc file-->they shall be divieded into spin files.
            read(ifvxc_)  nhdim,ldim
            read(ifvxc_)
            allocate( evec(nhdim,nhdim),evecrot(nhdim,nhdim))
            read(ifevec_) qqqx(1:3), evec(1:nhdim,1:nhdim),nev !nev number of true bands nov2015
            zsect = 0d0
          else                  !skip isx/=is. Need to get access sequential files evec and v_xc.
            read(ifvxc_)  
            read(ifvxc_)  
            read(ifevec_)
            cycle
          endif  
          do i=1,nnnx           !nq     !qqqx from evec v_xc.
            if(sum(abs(qqqx-q(:,i)))<tolq) then
              iqxx=i
              goto 3011
            endif  
          enddo  
          deallocate(evec,evecrot)
          call rx( 'hsfp0_sc: bug:qqqx can not find ...')
 3011     continue
          if(tiii) call rx( 'timereversal is not yet implemented')

!! evec_inv(ib1,iww)= \sum_ib2  ovlinv(ib1,ib2)*dconjg(evec(iww,ib2))  nov2015, we introduce nev. iww is for PMT basis. ib for band index.
!! This is for converting rotated evec (=evecrot(ib)) in the representation of original evec(ib).
          allocate(ovl(nev,nev)) 
c          print *,'nnnnnnnnn zsecsym: nband=',nhdim,nband,nev
          do i=1,nev
          do j=1,nev
c            write(6,*)'evec orth=',i,j,sum(dconjg(evec(:,i)*evec(:,j)))
            ovl(i,j)=sum(dconjg(evec(:,i))*evec(:,j))
          enddo
          enddo
          call matcinv(nev,ovl) !ovl --> ovlinv
          allocate(evec_inv(nev,nhdim))
          evec_inv = matmul(ovl(1:nev,1:nev),dconjg(transpose(evec(:,1:nev)))) !note ovl means ovlinv
          deallocate(ovl)
          evaliq = READEVAL(q(:,iqxx), is)
          nsym = sum(eibzsym(:,:,iqxx))
          do it=1,1             !no-time reversal yet !it=1,-1,-2 !c.f. x0kf_v4h
            do igrp=1,ngrp      !A-rotator
              if( eibzsym(igrp,it,iqxx)==0) cycle
              nblk=0
              iblki=0
              iblke=0
              iblki(1)=1
!!  degeneracy divider for evaliq. See How to apply EIBZ to 
!! Is this procedure really make speed up so much?
              tolry= 0.2d0      !Degeneracy tol. if tolry is large, 
!! larger tolry is safer, although a little inefficient.
!! If tolry is too small to divide degenerated values to different blocks --> then we have wrong results.
!(NOTE that Hamiltonian can be not so symmetric in some reasons) 
              nbsizemx=0
              do iband=2,nbandmx(iqxx,is) 
! nbandmx is the number of bands for which we calculate self-energy.
! We assume nbandmx(iqxx,is) is well separated for degeneracy.
                if(evaliq(iband) > evaliq(iband-1)+tolry 
     &           .or.iband==nbandmx(iqxx,is)) then
                  nblk=nblk+1
                  if(nblk>=2) iblki(nblk)=iblke(nblk-1)+1
                  if(iband==nbandmx(iqxx,is)) then
                    iblke(nblk)=iband
                  else
                    iblke(nblk)=iband-1
                  endif  
                  nbsize = iblke(nblk)- iblki(nblk)+1
                  if( nbsize>nbsizemx ) nbsizemx = nbsize
                endif
              enddo             ! iband
!! rotation of evec. Generate evecrot. (Within degenerated block, evec are mapped).e
              allocate(rmatjj(nbsizemx,nbsizemx,nblk))
              napw=nhdim-ldim
              do iblk1=1,nblk  
                ii1=iblki(iblk1)
                ie1=iblke(iblk1) 
                ne1=ie1-ii1+1
                call rotwvigg(igrp,q(:,iqxx),q(:,iqxx),nhdim,
     &           napw,ne1,evec(:,ii1:ie1),evecrot(:,ii1:ie1),ierr )
                rmatjj(1:ne1,1:ne1,iblk1) = 
     &           matmul(evec_inv(ii1:ie1,:),evecrot(:,ii1:ie1))
              enddo             ! iblk1
              do iblk1=1,nblk  
                do iblk2=1,nblk  
                  ii1=iblki(iblk1)
                  ie1=iblke(iblk1) 
                  ne1=ie1-ii1+1
                  ii2=iblki(iblk2)
                  ie2=iblke(iblk2) 
                  ne2=ie2-ii2+1
                  zsect(ii1:ie1,ii2:ie2)= zsect(ii1:ie1,ii2:ie2) 
     &             + matmul( dconjg(transpose(rmatjj(1:ne1,1:ne1,iblk1))), 
     &             matmul(zsec(ii1:ie1,ii2:ie2,iqxx),
     &             rmatjj(1:ne2,1:ne2,iblk2)) )
                enddo           ! iblk2
              enddo             ! iblk1
              deallocate(rmatjj)
            enddo               ! igrp
          enddo                 ! it
          deallocate(evec, evec_inv, evecrot)
          zsec(:,:,iqxx) = zsect(:,:)/dble(nsym)
c          call MPI__AllreduceSum( zsec(:,:,iqxx),ntq*ntq ) ! MIZUHO-IR
 3030   continue                ! ispx
 3020 continue                  ! iq
      do procid=0,nrankv-1
        close(ifvxc__(procid) )
        close(ifevec__(procid))
      enddo
      deallocate(iblki,iblke,evaliq)
      deallocate(zsect,extp,ifevec__,ifvxc__,iprocq)
      end subroutine zsecsym
