       program hdipole
!! == dipole and Hso matrix based on Wannier functions
c   use m_keyvalue,only: getkeyvalue
c   use m_readqg,only: readqg,readngmx
c   use m_readeigen,only: init_readeigen,init_readeigen2,readeval,lowesteval,readcphif,readgeigf,readcphifq
c      use m_DATA4GW,only: read_data4gw,set_mnla,iclass,nclass,zz,alat,nbas,nsp,plat,ldim2,bas
c     use m_QG,only: read_qg,ngp
c      use m_readhbe,only:Readhbe,nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
      use m_genallcf_v3,only: genallcf_v3,plat,nsp=>nspin
      use m_read_bzdata,only: read_bzdata,qibz, nqbz,nqibz,irk !,nqbzw,nteti,ntetf
     &     ,n1,n2,n3,ginv,qbz,wbz !wibz,qbzw,idtetf,ib1bz,idteti !qbasmc,
c     &     ,irk !,nstbz!,nstar,ngrp !2=>ngrp !,qibz_r,nqibz_r
      use m_hamindex,only: readhamindex,symops,ngrp,napwmx,ndham
      use m_ftox
      implicit none
      integer(4) :: nwf,iko_ix,iko_fx,nband_wfn,ik1(2),ik2(2),ik1m,ik2m
      real(8) :: q(3)
      real(8),allocatable:: qreg(:,:),ddd(:)
      complex(8),allocatable :: dnk(:,:,:,:)
      integer :: iq,ib,  ngpmx
      integer:: i,id,j,ifi,iqbz2,ifile_handle,incwfin
      integer :: isp,iqbz,iwf,iopen!,nprecb,mrecb,mrece,nlmtot,nqbzt,ifhbed, mrecg,nband 
      character*4::fname
      logical:: debug=.false.
      integer::nqbzx,ia,ifievec_,ifiproc,ifvxc_,ig,iqq,iqqxx,ispxx,ixxx !,iqibz
     &     ,is,iwfs,ixx,ndimh,ndimhx,nhq,nkp,iwfd,nmto,nnn,nnnx,nqixx,nrank,nspx,nspxx,
     &     ntq,nz,procid,nz1,ipmt,iz,izz,nwfo,iz1,iz2,iz1z,iz2z,ioffatom,ldim,napw
      complex(8),allocatable:: evecw(:,:,:),dipole(:,:,:,:),dipoig(:,:,:),dipoig1(:,:,:),evecr(:,:),
     &    evec(:,:,:),dipo(:,:,:,:),hammso(:,:,:),hso(:,:,:),hammsoig1(:,:,:),ovl(:,:,:),orth(:,:),orthcheck(:,:,:)
c      integer,allocatable:: ifevec__(:),ifvxc__(:),iprocq(:,:)
      character*256:: extn,ext
      character*256,allocatable:: extp(:)
      real(8)::qqq(3),qqqx(3),qtarget(3),tolq=1d-8   ,diff,ispblock,norm,ddd2
      integer(4),allocatable:: nev(:,:), nhqx(:,:)
      character(5):: charnum5
      complex(8):: ovi(2,2),aa,bb,cc,dd
      logical:: socmatrix
      character*8 xt
c!! MPI dummy
      include 'mpif.h'
      integer:: ierr
      call mpi_init(ierr)
!!      
      incwfin= -1   
      call GENALLCF_V3(incwfin)    ! readin basic data
      call readhamindex()
      call read_bzdata()
      
c$$$!! === readin lmfgw_kdivider, and get extensions
c$$$      open(newunit=ifiproc,file='lmfgw_kdivider',status='old')
c$$$      read(ifiproc,*) ext
c$$$      read(ifiproc,*) nqixx, nsp, nrank
c$$$      allocate(iprocq(nqixx,nsp))
c$$$      do isp=1,nsp
c$$$      do iqq=1,nqixx
c$$$        read(ifiproc,*) iqqxx, ispxx, ixxx
c$$$        if(iqqxx/=iqq) call rx( 'iqqxx/=iqq')
c$$$        if(ispxx/=isp) call rx( 'ispxx/=isp')
c$$$        iprocq(iqq,isp) = ixxx
c$$$        write(6,"('iqq isp irank=',i8,i2,i6)") iqq,isp, iprocq(iqq,isp)
c$$$      enddo
c$$$      enddo
c$$$      close(ifiproc)
      
c$$$!! Open vxc* evec*
c$$$      allocate(extp(0:nrank-1),ifvxc__(0:nrank-1),ifevec__(0:nrank-1))
c$$$      extp(0) = trim(ext)
c$$$      write(6,"('  0 ext= ',a,a)") trim(extp(0)),' ----------'
c$$$      do procid=1,nrank-1
c$$$        write(extn,"(i10)") procid
c$$$        extp(procid)=trim(adjustl(ext))//'_'//trim(adjustl(extn))
c$$$        write(6,"(i3,' ext= ',a,a)") procid,trim(extp(procid)),' ----------'
c$$$      enddo
c      do procid=0,nrank-1
c        open(newunit=ifvxc__(procid), file='vxc'//extp(procid),form='unformatted')
c        open(newunit=ifevec__(procid),file='evec'//extp(procid),form='unformatted')
c      enddo
      
!! Read dnk
      do isp=1,nsp
         if (isp.eq.1) fname='MLWU'
         if (isp.eq.2) fname='MLWD'
         open(newunit=ifi,file=fname,form='unformatted',status='old', action='read')
         read(ifi)nqbzx,nwf,iko_ix,iko_fx
         close(ifi)
         ik1(isp)=iko_ix
         ik2(isp)=iko_fx
         if(isp==2.and. nwfo/=nwf) call rx('nwf is spin-dependent')
         nwfo=nwf
      enddo          
      ik1m= minval(ik1(1:nsp))
      ik2m= maxval(ik2(1:nsp))
      allocate(dnk(ik1m:ik2m,nwf,nqbz,nsp) ,qreg(3,nqbz))

      do isp=1,nsp
         if (isp.eq.1) fname='MLWU'
         if (isp.eq.2) fname='MLWD'
         open(newunit=ifi,file=fname,form='unformatted',status='old', action='read')
         read(ifi)nqbzx,nwf,iko_ix,iko_fx
         write(6,'("nqbz isp iko=",4i5)')nqbz,isp,iko_ix,iko_fx
         if(nqbz/=nqbzx) then
            call rx('wanplot:nqbz/=nqbzx')
         endif   
         do iqbz = 1,nqbz
            read(ifi)iqbz2,q(1:3)
            if(debug)write(6,"(i5,3f13.5)") iqbz,q(:)
            qreg(:,iqbz) = q
            read(ifi)dnk(iko_ix:iko_fx,1:nwf,iqbz,isp)
         enddo
         close(ifi)
      enddo          
      write(6,*)'read end of MLWU/D ...'
!! readin vxc and evec ----------------------------------------------------
c      ifvxc_ = ifvxc__(0)
c      ifevec_= ifevec__(0)
c      read(ifvxc_)  ndimh,  nspxx,nnn 
c      read(ifevec_) ndimhx, nspx,nnnx
c      if (ndimhx /= ndimh)call rx( 'hdipole: ndimhx/=ndimx')
c      if (nsp .ne. nspx)  call rx( 'hdipole: nsp/=nspx')
c      if (nsp .ne. nspxx)  call rx( 'hdipole: nsp/=nspxx')
c      if (nnn .ne. nnnx)  call rx( 'hdipole: nnn/=nnnx')
c      if(nnn/=nqibz) call rx( ' nnn/nqibz ')
      nhq=ndham+napwmx !max of ndimhamiltonian
c      allocate(nhqx(nnn,nsp))
      allocate(evec(nhq,nhq,nsp),evecr(nhq,nhq))
      allocate(dipo(nhq,nhq,3,nsp),hso(nhq,nhq,1:3),ovl(nhq,nhq,nsp))
!! evec,dnk, dipo,hamso      
      allocate(evecw(nhq,nwf,nsp),dipole(3,nwf,nwf,nsp),dipoig(3,nwf,nwf),dipoig1(3,nwf,nwf))
      allocate(hammso(nwf,nwf,3),hammsoig1(nwf,nwf,3),orth(nwf,nwf),orthcheck(nwf,nwf,nsp),ddd(nhq))
      write(6,*)' ndimh nsp nnn =',nhq,nsp,nnn
      
!! -------------------------     
      hammso=0d0
      dipole=0d0
      orthcheck=0d0
      socmatrix=.true.
      do iqbz = 1,nqbz          !nqibz !index for q in irrecucible BZ.
         write(6,"('iqbz qbz(:,iqbz)=',i6,3f9.4)")iqbz,qbz(:,iqbz)
      do isp = 1,nsp !isp parallel for dipole, isp=nsp only for hammso
c          ifvxc_  = ifvxc__ (iprocq(iqbz,isp)) 
c          ifevec_ = ifevec__(iprocq(iqbz,isp)) 
          open(newunit=ifievec_,   file='evec'//trim(xt(iqbz))//trim(xt(isp)),form='unformatted')
          open(newunit=ifvxc_,      file='vxc'//trim(xt(iqbz))//trim(xt(isp)),form='unformatted')
          read(ifvxc_) nz!,ldim  ! Hamiltonian dimension, dependent of iqbz
          read(ifievec_) nz!,ldim  ! Hamiltonian dimension, dependent of iqbz
          if(isp==1) nz1=nz
          if(isp==2.and.nz1/=nz) call rx('hdipole: Hamiltonian dimension is spin dependent') 
          read(ifvxc_) 
          read(ifvxc_) dipo(1:nz,1:nz,1:3,isp),ovl(1:nz,1:nz,isp) !dipole matrix for pmt basis.
          if(isp==1.and.socmatrix) read(ifvxc_) hso(1:nz,1:nz,1:3) !see sugw.F for writing vxc file.
c          write(6,*) 'evec'//trim(xt(iqbz))//trim(xt(isp)),iqbz,isp,nz
          read(ifievec_) qqq,evec(1:nz,1:nz,isp) !evec(ipmt,iband,isp) for pmt basis.
c          write(6,"('iiiii iqbz isp qqq nz=',2i5,3f9.4,i5)") iqbz,isp,qqq,nz
c          write(6,ftox)' v_xc dipo hso: iq is nz=', iqbz,isp,nz,ftof(qqq,3)
c          if (sum(abs( qqq-qqqx ))        > tolq) call rx( 'hdipole: qqq/=qqqx')
          if (sum(abs( qqq-qbz(:,iqbz) )) > tolq) call rx( 'hdipole: qqq/=qbz')
          iko_ix=ik1(isp)
          iko_fx=ik2(isp)
          
c     do ig  = 1,ngrp  ! space group operation 
c     if(ig==1.and.irk(iqibz,ig)==0) call rx('sanity check: irk(iqibz,ig)==0 for ig=1')
c     if(irk(iqibz,ig)==0) cycle ! irk(iqibz,ig)=0 if iq is given by another rotation from iqibz.
c     iqbz = irk(iqibz,ig) !iqbz=irk(iqibz,ig) runs all iqbz in the 1st BZ.
!     q consistency check
c     q=      qbz(:,iqbz)-matmul(symops(:,:,ig),qqq)
c     &         -nint(qbz(:,iqbz)-matmul(symops(:,:,ig),qqq))
c     diff = sum(abs(  matmul(q,plat)))
c     if(diff>tolq) write(6,*) 'diff111=',iqbz, diff
c     q=      qreg(:,iqbz)-matmul(symops(:,:,ig),qqq)
c     &         -nint(qreg(:,iqbz)-matmul(symops(:,:,ig),qqq))
c     diff = sum(abs(  matmul(q,plat)))
c     if(diff>tolq) write(6,*) 'diff222=',iqbz, diff
c     ! rotate wave function coefficient  
c     napw= nz-ldim
c     c             write(6,"('rotwvigg qqq q=',3f9.4,x,3f9.4)") qqq,qbz(:,iqbz)
c     call rotwvigg(ig,qqq,qbz(:,iqbz),nz,
c     &            napw,nz,evec(1:nz,1:nz,isp),evecr(1:nz,1:nz),ierr)
c     if(ig==1) then !ig=0 is the identical rotation of space group
c     if(irk(iqibz,ig)==0) call rx('bug: irk(iqibz,ig)==0 for ig=1')
          
          do iwf = 1,nwf        ! iko_ix:iko_fx is the range of bands to construct MlocWannier.
             evecw(1:nz,iwf,isp)=  matmul(
     &            evec(1:nz,iko_ix:iko_fx,isp),dnk(iko_ix:iko_fx,iwf,iqbz,isp)) !iqbz component of Wannier based on PMT basis.
          enddo                 !evecw is the component of Wannier for qbz(:,iqbz)

          
!! component check at Gamma point (iqbz=1)
          if(iqbz==1) then !.and.ig==1) then
             do iwf= 1,nwf
                norm = sum(dconjg(evecw(1:nz,iwf,isp))*  matmul(ovl(1:nz,1:nz,isp),evecw(1:nz,iwf,isp)))
                ioffatom=0 
                do iz = 1,5     !nz
                   iz1= ioffatom + 1+ 3+ iz !PMT index for 3d EH,    EH=-1 -1 -1 -1 =1+3+5+7
                   iz2= ioffatom + 1+ 3+ 5+ 7+ 1+ 3+ iz !PMT index for 3d EH2,   EH2=-2 -2 -2    =1+3+5
                   aa= ovl(iz1,iz1,isp)
                   bb= ovl(iz1,iz2,isp)
                   cc= ovl(iz2,iz1,isp)
                   dd= ovl(iz2,iz2,isp)
                   ovi(1,1)=  dd/(aa*dd-bb*cc)
                   ovi(1,2)= -bb/(aa*dd-bb*cc)
                   ovi(2,1)= -cc/(aa*dd-bb*cc)
                   ovi(2,2)=  aa/(aa*dd-bb*cc)
                   ddd(iz)=
     & sum(dconjg(evecw(1:nz,iwf,isp))*ovl(1:nz,iz1,isp))*ovi(1,1)*sum(ovl(iz1,1:nz,isp)*evecw(1:nz,iwf,isp))
     &+sum(dconjg(evecw(1:nz,iwf,isp))*ovl(1:nz,iz1,isp))*ovi(1,2)*sum(ovl(iz2,1:nz,isp)*evecw(1:nz,iwf,isp))
     &+sum(dconjg(evecw(1:nz,iwf,isp))*ovl(1:nz,iz2,isp))*ovi(2,1)*sum(ovl(iz1,1:nz,isp)*evecw(1:nz,iwf,isp))
     &+sum(dconjg(evecw(1:nz,iwf,isp))*ovl(1:nz,iz2,isp))*ovi(2,2)*sum(ovl(iz2,1:nz,isp)*evecw(1:nz,iwf,isp))
                enddo   
      write(6,"('dddd 3d comp at Gamma:isp iwf norm evecw(3d) =',i2,i3,f6.3,2x,5f7.4)") isp,iwf,norm,ddd(1:5)
             enddo
          endif
!!          
          do iwf  = 1,nwf
          do iwfd = 1,nwf
                do ia=1,3
c                   dipoig1(ia,iwf,iwfd)=
c     &                  sum(dconjg(evecw(1:nz,iwf,isp))*matmul(dipo(1:nz,1:nz,ia,isp),evecw(1:nz,iwfd,isp)))
                   orth(iwf,iwfd)= !Orthgonalization check of Wannier functions                     
     &                  sum(dconjg(evecw(1:nz,iwf,isp))*matmul(ovl(1:nz,1:nz,isp),evecw(1:nz,iwfd,isp)))
                enddo          
                if(isp==nsp) then
                   hammsoig1(iwf,iwfd,1)=
     &                  sum(dconjg(evecw(1:nz,iwf,1))*matmul(hso(1:nz,1:nz,1),evecw(1:nz,iwfd,1)))
                   hammsoig1(iwf,iwfd,2)=
     &                  sum(dconjg(evecw(1:nz,iwf,2))*matmul(hso(1:nz,1:nz,2),evecw(1:nz,iwfd,2)))
                   hammsoig1(iwf,iwfd,3)= ! Hso(isp=1,isp=2) block elements
     &                  sum(dconjg(evecw(1:nz,iwf,1))*matmul(hso(1:nz,1:nz,3),evecw(1:nz,iwfd,2)))
                endif  
          enddo          
          enddo
 ! accumulate hammso and dipole for q points in the 1st BZ
          if(isp==nsp) hammso(:,:,:) = hammso(:,:,:) + hammsoig1
c          dipole(:,:,:,isp) = dipole(:,:,:,isp) + dipoig1
          orthcheck(:,:,isp)= orthcheck(:,:,isp) + orth
      enddo 
      enddo
!! normalization. averaged in       
      orthcheck=orthcheck/dble(nqbz)
c      dipole=dipole/dble(nqbz)
      hammso=hammso/dble(nqbz)
      
!! check orthgonalization of Wannier
      do isp=1,nsp
      do iwf  = 1,nwf           !dipoig is the rotation of dipoig1 by symops
         iwfd = iwf
         write(6,ftox)'orthcheck diag=',iwf,iwfd,isp,ftod(orthcheck(iwf,iwfd,isp))
      enddo          
      do iwf  = 1,nwf           !dipoig is the rotation of dipoig1 by symops
      do iwfd = 1,nwf
         if(iwf==iwfd) cycle   
         if(abs(orthcheck(iwf,iwfd,isp))>1d-10)
     &      write(6,"('WARN! orthcheck offd=',2i3,i2,3(2d11.3,x))")iwf,iwfd,isp, orthcheck(iwf,iwfd,isp)
      enddo          
      enddo
      enddo
!!      
c      do isp=1,nsp
c      do iwf  = 1,nwf           !dipoig is the rotation of dipoig1 by symops
c      do iwfd = 1,nwf
c         write(6,ftox)'dipole real=',iwf,iwfd,isp,ftof(dipole(:,iwf,iwfd,isp),6)
c      enddo          
c      enddo
c      enddo
      do iwf  = 1,nwf        
      do iwfd = 1,nwf
         write(6,ftox)'dddd hammso real =',iwf,iwfd,ftof(hammso(iwf,iwfd,1:3),8)
      enddo          
      enddo
      print *,'hdipole ok'
      print *,'2021-5-11: something wrong. diagonal dipole should point at the wannier center but not...'
      call rx0s('hdipole: ok')
      end program hdipole

