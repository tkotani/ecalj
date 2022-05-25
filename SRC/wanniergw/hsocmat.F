       program hsocmat
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
     &     ntq,nz,procid,nz1,ipmt,iz,izz,nwfo,iz1,iz2,iz1z,iz2z,ioffatom,ldim,napw,iwf2
      complex(8),allocatable:: evecw(:,:,:),dipole(:,:,:,:),dipoig(:,:,:),dipoig1(:,:,:),evecr(:,:),
     &    evec(:,:,:),dipo(:,:,:,:),hammso(:,:,:),hso2(:,:,:),hso(:,:,:),hammsoig1(:,:,:),ovl(:,:,:),orth(:,:),orthcheck(:,:,:)
     &     ,hammv(:,:,:),psig(:,:,:,:)
c      integer,allocatable:: ifevec__(:),ifvxc__(:),iprocq(:,:)
      character*256:: extn,ext
      character*256,allocatable:: extp(:)
      real(8)::qqq(3),qqqx(3),qtarget(3),tolq=1d-8   ,diff,ispblock,norm,ddd2,sss
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

c$$$cccccccccccccccccc
c$$$      allocate(psig(iko_ix:iko_fx,nwf,nqbz,nsp))
c$$$      do isp=1,nsp
c$$$         if (isp.eq.1) fname='U'
c$$$         if (isp.eq.2) fname='D'
c$$$         open(newunit=ifi,file='PSIG'//trim(fname),form='unformatted',status='old', action='read')
c$$$         block
c$$$         integer :: nqbz2,iti2,itf2,nwf2,iqbz2,j1,j2,iti,itf
c$$$         iti=iko_ix
c$$$         itf=iko_fx
c$$$         read(ifi)nqbz2,iti2,itf2,nwf2
c$$$         block
c$$$         complex(8):: psigread(iti2:itf2,nwf2,nqbz2)
c$$$         do iqbz = 1,nqbz
c$$$            read(ifi)iqbz2
c$$$            if(iqbz2 .ne. iqbz) call rx('readpsig: iqbz error')
c$$$            read(ifi)((psigread(j1,j2,iqbz),j1=iti2,itf2),j2=1,nwf)
c$$$         enddo
c$$$         psig(iti:itf,:,:,isp)=psigread(iti:itf,:,:)
c$$$         endblock
c$$$         endblock
c$$$         close(ifi)
c$$$      enddo
c$$$cccccccccccccccccccc

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
ccccccccccccccccccccccc
c      allocate(hso2(nhq,nhq,1:3))
cccccccccccccccccccccccc
!! -------------------------     
      hammso=0d0
      dipole=0d0
c      orthcheck=0d0
      socmatrix=.true.
      do 1012 iqbz = 1,nqbz          !nqibz !index for q in irrecucible BZ.
         write(6,"('iqbz qbz(:,iqbz)=',i6,3f9.4)")iqbz,qbz(:,iqbz)
      do 1014 isp = 1,nsp !isp parallel for dipole, isp=nsp only for hammso
c          ifvxc_  = ifvxc__ (iprocq(iqbz,isp)) 
c          ifevec_ = ifevec__(iprocq(iqbz,isp)) 
          open(newunit=ifievec_,   file='evec'//trim(xt(iqbz))//trim(xt(isp)),form='unformatted')
          open(newunit=ifvxc_,      file='vxc'//trim(xt(iqbz))//trim(xt(isp)),form='unformatted')
          read(ifvxc_) nz!,ldim  ! Hamiltonian dimension, dependent of iqbz
          read(ifievec_) nz     !,ldim  ! Hamiltonian dimension, dependent of iqbz
c          print *,'iiiiqbz isp nz',iqbz,isp,nz
c          if(isp==1) nz1=nz
c          if(isp==2.and.nz1/=nz) call rx('hdipole: Hamiltonian dimension is spin dependent') 
          read(ifvxc_) 
c          read(ifvxc_) dipo(1:nz,1:nz,1:3,isp),ovl(1:nz,1:nz,isp) !dipole matrix for pmt basis.
          if(isp==1.and.socmatrix) read(ifvxc_) hso(1:nz,1:nz,1:3) !see sugw.F for writing vxc file.
          close(ifvxc_)
c          if(isp==1.and.socmatrix) write(6,ftox)'hhhhso',iqbz, sum(abs( hso(1:nz,1:nz,1)))
c     &         , sum(abs( hso(1:nz,1:nz,2))), sum(abs( hso(1:nz,1:nz,3)))
ccccccccccccccccccccccccccccccccccccccccccccc
c          open(newunit=ifvxc_, file='uuu/vxc'//trim(xt(iqbz))//trim(xt(isp)),form='unformatted')
c          read(ifvxc_) nz!,ldim  ! Hamiltonian dimension, dependent of iqbz
c          read(ifvxc_) 
c          if(isp==1.and.socmatrix) read(ifvxc_) hso2(1:nz,1:nz,1:3) !see sugw.F for writing vxc file.
c          close(ifvxc_)
c          if(isp==1.and.socmatrix) then
c             write(6,ftox)'hhh222hso',iqbz,(sum(abs(hso(1:nz,1:nz,i)-hso2(1:nz,1:nz,i))),i=1,3),(sum(abs(hso(1:nz,1:nz,i))),i=1,3)
c          endif

c         write(6,*) 'evec'//trim(xt(iqbz))//trim(xt(isp)),iqbz,isp,nz
          read(ifievec_) qqq,evec(1:nz,1:nz,isp) !evec(ipmt,iband,isp) for pmt basis.
c          write(6,"('iiiii iqbz isp qqq nz=',2i5,3f9.4,i5)") iqbz,isp,qqq,nz
c          write(6,ftox)' v_xc dipo hso: iq is nz=', iqbz,isp,nz,ftof(qqq,3)
c          if (sum(abs( qqq-qqqx ))        > tolq) call rx( 'hdipole: qqq/=qqqx')
          if (sum(abs( qqq-qbz(:,iqbz) )) > tolq) call rx( 'hdipole: qqq/=qbz')
          iko_ix=ik1(isp)
          iko_fx=ik2(isp)
c          if(allocated(hammv)) deallocate(hammv)
c          allocate(hammv(iko_ix:iko_fx,iko_ix:iko_fx,3))
          
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
c             evecw(1:nz,iwf,isp)=  matmul(
c     &            evec(1:nz,iko_ix:iko_fx,isp),psig(iko_ix:iko_fx,iwf,iqbz,isp)) !iqbz component of Wannier based on PMT basis.
          enddo                 !evecw is the component of Wannier for qbz(:,iqbz)

c$$$ccccccccccccccccccccccccccccccccc
c$$$          write(6,*) '--------- eeev1 eeev2 eeev3 start',iqbz,isp
c$$$          iwf=1
c$$$          write(6,ftox)'eeev2',isp,iqbz,ftof(evecw(10:16,iwf,isp),3)
c$$$          if(isp==2.and.iqbz==5) then
c$$$             write(6,ftox)' dnk abs ',isp,iqbz, ftod(abs(dnk(iko_ix:iko_fx,1,iqbz,isp)),2)
c$$$             write(6,ftox)' dnk phas',isp,iqbz, ftod(dimag(log(dnk(iko_ix:iko_fx,1,iqbz,isp))),2)
c$$$             do iz=10,16
c$$$c                write(6,ftox)'evec abs  ',isp,iqbz,iz,ftod(abs(evec(iz,iko_ix:iko_fx,isp)),2)
c$$$c                write(6,ftox)'evec phase',isp,iqbz,iz
c$$$c     &               ,ftof(dimag(log(evec(iz,iko_ix:iko_fx,isp)/evec(1,iko_ix:iko_fx,isp))),2)
c$$$                write(6,ftox)'     ab3  ',isp,iqbz,iz,ftof(abs(evec(iz,iko_ix:iko_fx,isp)*
c$$$     &                                                dnk(iko_ix:iko_fx,1,iqbz,isp)),3)
c$$$                write(6,ftox)'     ph3   ',isp,iqbz,iz,ftof(dimag(log(evec(iz,iko_ix:iko_fx,isp)*
c$$$     &                                                dnk(iko_ix:iko_fx,1,iqbz,isp))),3)
c$$$             enddo
c$$$          endif
c$$$cccccccccccccccccccccccccccccc

          
c$$$          
c$$$!! component check at Gamma point (iqbz=1)
c$$$          if(iqbz==1) then !.and.ig==1) then
c$$$             do iwf= 1,nwf
c$$$                norm = sum(dconjg(evecw(1:nz,iwf,isp))*  matmul(ovl(1:nz,1:nz,isp),evecw(1:nz,iwf,isp)))
c$$$                ioffatom=0 
c$$$                do iz = 1,5     !nz
c$$$                   iz1= ioffatom + 1+ 3+ iz !PMT index for 3d EH,    EH=-1 -1 -1 -1 =1+3+5+7
c$$$                   iz2= ioffatom + 1+ 3+ 5+ 7+ 1+ 3+ iz !PMT index for 3d EH2,   EH2=-2 -2 -2    =1+3+5
c$$$                   aa= ovl(iz1,iz1,isp)
c$$$                   bb= ovl(iz1,iz2,isp)
c$$$                   cc= ovl(iz2,iz1,isp)
c$$$                   dd= ovl(iz2,iz2,isp)
c$$$                   ovi(1,1)=  dd/(aa*dd-bb*cc)
c$$$                   ovi(1,2)= -bb/(aa*dd-bb*cc)
c$$$                   ovi(2,1)= -cc/(aa*dd-bb*cc)
c$$$                   ovi(2,2)=  aa/(aa*dd-bb*cc)
c$$$                   ddd(iz)=
c$$$     & sum(dconjg(evecw(1:nz,iwf,isp))*ovl(1:nz,iz1,isp))*ovi(1,1)*sum(ovl(iz1,1:nz,isp)*evecw(1:nz,iwf,isp))
c$$$     &+sum(dconjg(evecw(1:nz,iwf,isp))*ovl(1:nz,iz1,isp))*ovi(1,2)*sum(ovl(iz2,1:nz,isp)*evecw(1:nz,iwf,isp))
c$$$     &+sum(dconjg(evecw(1:nz,iwf,isp))*ovl(1:nz,iz2,isp))*ovi(2,1)*sum(ovl(iz1,1:nz,isp)*evecw(1:nz,iwf,isp))
c$$$     &+sum(dconjg(evecw(1:nz,iwf,isp))*ovl(1:nz,iz2,isp))*ovi(2,2)*sum(ovl(iz2,1:nz,isp)*evecw(1:nz,iwf,isp))
c$$$                enddo   
c$$$c      write(6,"('dddd 3d comp at Gamma:isp iwf norm evecw(3d) =',i2,i3,f6.3,2x,5f7.4)") isp,iwf,norm,ddd(1:5)
c$$$c             enddo
c$$$          endif
c$$$!!          
          do iwf  = 1,nwf
          do iwfd = 1,nwf
c                do ia=1,3
c                   dipoig1(ia,iwf,iwfd)=
c     &                  sum(dconjg(evecw(1:nz,iwf,isp))*matmul(dipo(1:nz,1:nz,ia,isp),evecw(1:nz,iwfd,isp)))
c                   orth(iwf,iwfd)= !Orthgonalization check of Wannier functions                     
c     &                  sum(dconjg(evecw(1:nz,iwf,isp))*matmul(ovl(1:nz,1:nz,isp),evecw(1:nz,iwfd,isp)))
c                enddo          

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
c          if(isp==nsp) then
c          do iwf  = iko_ix,iko_fx 
c          do iwfd = iko_ix,iko_fx 
c                   hammv(iwf,iwfd,1)=
c     &                  sum(dconjg(evec(1:nz,iwf,1))*matmul(hso(1:nz,1:nz,1),evec(1:nz,iwfd,1)))
c                   hammv(iwf,iwfd,2)=
c     &                  sum(dconjg(evec(1:nz,iwf,2))*matmul(hso(1:nz,1:nz,2),evec(1:nz,iwfd,2)))
c                   hammv(iwf,iwfd,3)= ! Hso(isp=1,isp=2) block elements
c     &                  sum(dconjg(evec(1:nz,iwf,1))*matmul(hso(1:nz,1:nz,3),evec(1:nz,iwfd,2)))
c          enddo
c          enddo
c          endif
 ! accumulate hammso and dipole for q points in the 1st BZ
          if(isp==nsp) hammso(:,:,:) = hammso(:,:,:) + hammsoig1
c          if(isp==nsp) then
c          do iwf=iko_ix,iko_fx 
c             write(6,ftox) 'vvvhammv1=',iqbz,iwf, ftof([(hammv(iwf2,iwf,1),iwf2=iko_ix,iko_fx) ],3)
c             write(6,ftox) 'vvvhammv2=',iqbz,iwf, ftof([(hammv(iwf2,iwf,2),iwf2=iko_ix,iko_fx) ],3)
c             write(6,ftox) 'vvvhammv3=',iqbz,iwf, ftof([(hammv(iwf2,iwf,3),iwf2=iko_ix,iko_fx) ],3)
c          enddo   
c          endif
cc          dipole(:,:,:,isp) = dipole(:,:,:,isp) + dipoig1
cc          orthcheck(:,:,isp)= orthcheck(:,:,isp) + orth
 1014  enddo 
 1012 enddo

!! normalization. averaged in       
cc      orthcheck=orthcheck/dble(nqbz)
c      dipole=dipole/dble(nqbz)
      hammso=hammso/dble(nqbz)
      
!! check orthgonalization of Wannier
      do isp=1,nsp
c      do iwf  = 1,nwf           !dipoig is the rotation of dipoig1 by symops
c         iwfd = iwf
c         write(6,ftox)'orthcheck diag=',iwf,iwfd,isp,ftod(orthcheck(iwf,iwfd,isp))
c      enddo          
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
!
      block
      real(8):: hammso2(nwf,nwf),sss=0d0
      hammso2 = matmul(hammso(:,:,1),hammso(:,:,1)) + matmul(hammso(:,:,2),hammso(:,:,2))
     & + matmul(hammso(:,:,3),transpose(hammso(:,:,3)))+ matmul(transpose(hammso(:,:,3)),hammso(:,:,3))
      do iwf  = 1,nwf        
         sss=sss+ hammso2(iwf,iwf)
      enddo
      write(6,*)"hhhddd hammso**2 diagonal=",sss
      endblock

      print *,'hsocmat ok'
      call rx0s('hsocmat: ok')
      end program hsocmat

