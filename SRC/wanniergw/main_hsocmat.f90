subroutine hsocmat()
  !! == dipole and Hso matrix based on Wannier functions
  !   use m_keyvalue,only: getkeyvalue
  !   use m_readqg,only: readqg,readngmx
  !   use m_readeigen,only: init_readeigen,init_readeigen2,readeval,lowesteval,readcphif,readgeigf,readcphifq
  !      use m_DATA4GW,only: read_data4gw,set_mnla,iclass,nclass,zz,alat,nbas,nsp,plat,ldim2,bas
  !     use m_QG,only: read_qg,ngp
  !      use m_readhbe,only:Readhbe,nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
  use m_genallcf_v3,only: genallcf_v3,plat,nsp=>nspin
  use m_read_bzdata,only: read_bzdata,qibz, nqbz,nqibz,irk,n1,n2,n3,ginv,qbz,wbz
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
  integer::nqbzx,ia,ifievec_,ifiproc,ifvxc_,ig,iqq,iqqxx,ispxx,ixxx  &
       ,is,iwfs,ixx,ndimh,ndimhx,nhq,nkp,iwfd,nmto,nnn,nnnx,nqixx,nrank,nspx,nspxx, &
       ntq,nz,procid,nz1,ipmt,iz,izz,nwfo,iz1,iz2,iz1z,iz2z,ioffatom,ldim,napw,iwf2
  complex(8),allocatable:: evecw(:,:,:),dipole(:,:,:,:),dipoig(:,:,:),dipoig1(:,:,:),evecr(:,:), &
       evec(:,:,:),dipo(:,:,:,:),hammso(:,:,:),hso2(:,:,:),hso(:,:,:),hammsoig1(:,:,:),ovl(:,:,:),orth(:,:),orthcheck(:,:,:) &
       ,hammv(:,:,:),psig(:,:,:,:)
  !      integer,allocatable:: ifevec__(:),ifvxc__(:),iprocq(:,:)
  character*256:: extn,ext
  character*256,allocatable:: extp(:)
  real(8)::qqq(3),qqqx(3),qtarget(3),tolq=1d-8   ,diff,ispblock,norm,ddd2
  integer(4),allocatable:: nev(:,:), nhqx(:,:)
  character(5):: charnum5
  complex(8):: ovi(2,2),aa,bb,cc,dd
  logical:: socmatrix
  character(8) :: xt
!!! MPI dummy
  include 'mpif.h'
  integer:: ierr
  call mpi_init(ierr)
  !!
  incwfin= -1
  call GENALLCF_V3(incwfin)    ! readin basic data
  call readhamindex()
  call read_bzdata()

  !! Read dnk
  do isp=1,nsp
     if (isp == 1) fname='MLWU'
     if (isp == 2) fname='MLWD'
     open(newunit=ifi,file=fname,form='unformatted',status='old', action='read')
     read(ifi)nqbzx,nwf,iko_ix,iko_fx
     close(ifi)
     ik1(isp)=iko_ix
     ik2(isp)=iko_fx
     if(isp==2 .AND. nwfo/=nwf) call rx('nwf is spin-dependent')
     nwfo=nwf
  enddo
  ik1m= minval(ik1(1:nsp))
  ik2m= maxval(ik2(1:nsp))
  allocate(dnk(ik1m:ik2m,nwf,nqbz,nsp) ,qreg(3,nqbz))

  do isp=1,nsp
     if (isp == 1) fname='MLWU'
     if (isp == 2) fname='MLWD'
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
  nhq=ndham+napwmx !max of ndimhamiltonian
  allocate(evec(nhq,nhq,nsp),evecr(nhq,nhq))
  allocate(dipo(nhq,nhq,3,nsp),hso(nhq,nhq,1:3),ovl(nhq,nhq,nsp))
  !! evec,dnk, dipo,hamso
  allocate(evecw(nhq,nwf,nsp),dipole(3,nwf,nwf,nsp),dipoig(3,nwf,nwf),dipoig1(3,nwf,nwf))
  allocate(hammso(nwf,nwf,3),hammsoig1(nwf,nwf,3),orth(nwf,nwf),orthcheck(nwf,nwf,nsp),ddd(nhq))
  write(6,*)' ndimh nsp nnn =',nhq,nsp,nnn
  !! -------------------------
  hammso=0d0
  dipole=0d0
  !      orthcheck=0d0
  socmatrix=.true.
  do 1012 iqbz = 1,nqbz          !nqibz !index for q in irrecucible BZ.
     write(6,"('iqbz qbz(:,iqbz)=',i6,3f9.4)")iqbz,qbz(:,iqbz)
     do 1014 isp = 1,nsp !isp parallel for dipole, isp=nsp only for hammso
        open(newunit=ifievec_,   file='evec'//trim(xt(iqbz))//trim(xt(isp)),form='unformatted')
        open(newunit=ifvxc_,      file='vxc'//trim(xt(iqbz))//trim(xt(isp)),form='unformatted')
        read(ifvxc_) nz!,ldim  ! Hamiltonian dimension, dependent of iqbz
        read(ifievec_) nz     !,ldim  ! Hamiltonian dimension, dependent of iqbz
        read(ifvxc_)
        !          read(ifvxc_) dipo(1:nz,1:nz,1:3,isp),ovl(1:nz,1:nz,isp) !dipole matrix for pmt basis.
        if(isp==1 .AND. socmatrix) read(ifvxc_) hso(1:nz,1:nz,1:3) !see sugw.F for writing vxc file.
        close(ifvxc_)
        read(ifievec_) qqq,evec(1:nz,1:nz,isp) !evec(ipmt,iband,isp) for pmt basis.
        if (sum(abs( qqq-qbz(:,iqbz) )) > tolq) call rx( 'hdipole: qqq/=qbz')
        iko_ix=ik1(isp)
        iko_fx=ik2(isp)
        do iwf = 1,nwf        ! iko_ix:iko_fx is the range of bands to construct MlocWannier.
           evecw(1:nz,iwf,isp)=  matmul( &
                evec(1:nz,iko_ix:iko_fx,isp),dnk(iko_ix:iko_fx,iwf,iqbz,isp)) !iqbz component of Wannier based on PMT basis.
           !             evecw(1:nz,iwf,isp)=  matmul(
           !     &            evec(1:nz,iko_ix:iko_fx,isp),psig(iko_ix:iko_fx,iwf,iqbz,isp)) !iqbz component of Wannier based on PMT basis.
        enddo                 !evecw is the component of Wannier for qbz(:,iqbz)

        !$$$ccccccccccccccccccccccccccccccccc
        !$$$          write(6,*) '--------- eeev1 eeev2 eeev3 start',iqbz,isp
        !$$$          iwf=1
        !$$$          write(6,ftox)'eeev2',isp,iqbz,ftof(evecw(10:16,iwf,isp),3)
        !$$$          if(isp==2.and.iqbz==5) then
        !$$$             write(6,ftox)' dnk abs ',isp,iqbz, ftod(abs(dnk(iko_ix:iko_fx,1,iqbz,isp)),2)
        !$$$             write(6,ftox)' dnk phas',isp,iqbz, ftod(dimag(log(dnk(iko_ix:iko_fx,1,iqbz,isp))),2)
        !$$$             do iz=10,16
        !$$$c                write(6,ftox)'evec abs  ',isp,iqbz,iz,ftod(abs(evec(iz,iko_ix:iko_fx,isp)),2)
        !$$$c                write(6,ftox)'evec phase',isp,iqbz,iz
        !$$$c     &               ,ftof(dimag(log(evec(iz,iko_ix:iko_fx,isp)/evec(1,iko_ix:iko_fx,isp))),2)
        !$$$                write(6,ftox)'     ab3  ',isp,iqbz,iz,ftof(abs(evec(iz,iko_ix:iko_fx,isp)*
        !$$$     &                                                dnk(iko_ix:iko_fx,1,iqbz,isp)),3)
        !$$$                write(6,ftox)'     ph3   ',isp,iqbz,iz,ftof(dimag(log(evec(iz,iko_ix:iko_fx,isp)*
        !$$$     &                                                dnk(iko_ix:iko_fx,1,iqbz,isp))),3)
        !$$$             enddo
        !$$$          endif
        !$$$cccccccccccccccccccccccccccccc


        !$$$
        !$$$!! component check at Gamma point (iqbz=1)
        !$$$          if(iqbz==1) then !.and.ig==1) then
        !$$$             do iwf= 1,nwf
        !$$$                norm = sum(dconjg(evecw(1:nz,iwf,isp))*  matmul(ovl(1:nz,1:nz,isp),evecw(1:nz,iwf,isp)))
        !$$$                ioffatom=0
        !$$$                do iz = 1,5     !nz
        !$$$                   iz1= ioffatom + 1+ 3+ iz !PMT index for 3d EH,    EH=-1 -1 -1 -1 =1+3+5+7
        !$$$                   iz2= ioffatom + 1+ 3+ 5+ 7+ 1+ 3+ iz !PMT index for 3d EH2,   EH2=-2 -2 -2    =1+3+5
        !$$$                   aa= ovl(iz1,iz1,isp)
        !$$$                   bb= ovl(iz1,iz2,isp)
        !$$$                   cc= ovl(iz2,iz1,isp)
        !$$$                   dd= ovl(iz2,iz2,isp)
        !$$$                   ovi(1,1)=  dd/(aa*dd-bb*cc)
        !$$$                   ovi(1,2)= -bb/(aa*dd-bb*cc)
        !$$$                   ovi(2,1)= -cc/(aa*dd-bb*cc)
        !$$$                   ovi(2,2)=  aa/(aa*dd-bb*cc)
        !$$$                   ddd(iz)=
        !$$$     & sum(dconjg(evecw(1:nz,iwf,isp))*ovl(1:nz,iz1,isp))*ovi(1,1)*sum(ovl(iz1,1:nz,isp)*evecw(1:nz,iwf,isp))
        !$$$     &+sum(dconjg(evecw(1:nz,iwf,isp))*ovl(1:nz,iz1,isp))*ovi(1,2)*sum(ovl(iz2,1:nz,isp)*evecw(1:nz,iwf,isp))
        !$$$     &+sum(dconjg(evecw(1:nz,iwf,isp))*ovl(1:nz,iz2,isp))*ovi(2,1)*sum(ovl(iz1,1:nz,isp)*evecw(1:nz,iwf,isp))
        !$$$     &+sum(dconjg(evecw(1:nz,iwf,isp))*ovl(1:nz,iz2,isp))*ovi(2,2)*sum(ovl(iz2,1:nz,isp)*evecw(1:nz,iwf,isp))
        !$$$                enddo
        !$$$c      write(6,"('dddd 3d comp at Gamma:isp iwf norm evecw(3d) =',i2,i3,f6.3,2x,5f7.4)") isp,iwf,norm,ddd(1:5)
        !$$$c             enddo
        !$$$          endif
        !$$$!!
        do iwf  = 1,nwf
           do iwfd = 1,nwf
              !                do ia=1,3
              !                   dipoig1(ia,iwf,iwfd)=
              !     &                  sum(dconjg(evecw(1:nz,iwf,isp))*matmul(dipo(1:nz,1:nz,ia,isp),evecw(1:nz,iwfd,isp)))
              !                   orth(iwf,iwfd)= !Orthgonalization check of Wannier functions
              !     &                  sum(dconjg(evecw(1:nz,iwf,isp))*matmul(ovl(1:nz,1:nz,isp),evecw(1:nz,iwfd,isp)))
              !                enddo

              if(isp==nsp) then
                 hammsoig1(iwf,iwfd,1)= &
                      sum(dconjg(evecw(1:nz,iwf,1))*matmul(hso(1:nz,1:nz,1),evecw(1:nz,iwfd,1)))
                 hammsoig1(iwf,iwfd,2)= &
                      sum(dconjg(evecw(1:nz,iwf,2))*matmul(hso(1:nz,1:nz,2),evecw(1:nz,iwfd,2)))
                 hammsoig1(iwf,iwfd,3)= & ! &  Hso(isp=1,isp=2) block elements
                 sum(dconjg(evecw(1:nz,iwf,1))*matmul(hso(1:nz,1:nz,3),evecw(1:nz,iwfd,2)))
              endif
           enddo
        enddo
        !          if(isp==nsp) then
        !          do iwf  = iko_ix,iko_fx
        !          do iwfd = iko_ix,iko_fx
        !                   hammv(iwf,iwfd,1)=
        !     &                  sum(dconjg(evec(1:nz,iwf,1))*matmul(hso(1:nz,1:nz,1),evec(1:nz,iwfd,1)))
        !                   hammv(iwf,iwfd,2)=
        !     &                  sum(dconjg(evec(1:nz,iwf,2))*matmul(hso(1:nz,1:nz,2),evec(1:nz,iwfd,2)))
        !                   hammv(iwf,iwfd,3)= ! Hso(isp=1,isp=2) block elements
        !     &                  sum(dconjg(evec(1:nz,iwf,1))*matmul(hso(1:nz,1:nz,3),evec(1:nz,iwfd,2)))
        !          enddo
        !          enddo
        !          endif
        ! accumulate hammso and dipole for q points in the 1st BZ
        if(isp==nsp) hammso(:,:,:) = hammso(:,:,:) + hammsoig1
        !          if(isp==nsp) then
        !          do iwf=iko_ix,iko_fx
        !             write(6,ftox) 'vvvhammv1=',iqbz,iwf, ftof([(hammv(iwf2,iwf,1),iwf2=iko_ix,iko_fx) ],3)
        !             write(6,ftox) 'vvvhammv2=',iqbz,iwf, ftof([(hammv(iwf2,iwf,2),iwf2=iko_ix,iko_fx) ],3)
        !             write(6,ftox) 'vvvhammv3=',iqbz,iwf, ftof([(hammv(iwf2,iwf,3),iwf2=iko_ix,iko_fx) ],3)
        !          enddo
        !          endif
        !c          dipole(:,:,:,isp) = dipole(:,:,:,isp) + dipoig1
        !c          orthcheck(:,:,isp)= orthcheck(:,:,isp) + orth
        close(ifievec_) !suzuki 2023feb
1014 enddo
1012 enddo

!! normalization. averaged in
!c      orthcheck=orthcheck/dble(nqbz)
!      dipole=dipole/dble(nqbz)
hammso=hammso/dble(nqbz)

!! check orthgonalization of Wannier
do isp=1,nsp
   !      do iwf  = 1,nwf           !dipoig is the rotation of dipoig1 by symops
   !         iwfd = iwf
   !         write(6,ftox)'orthcheck diag=',iwf,iwfd,isp,ftod(orthcheck(iwf,iwfd,isp))
   !      enddo
do iwf  = 1,nwf           !dipoig is the rotation of dipoig1 by symops
  do iwfd = 1,nwf
     if(iwf==iwfd) cycle
     if(abs(orthcheck(iwf,iwfd,isp))>1d-10) &
          write(6,"('WARN! orthcheck offd=',2i3,i2,3(2d11.3,x))")iwf,iwfd,isp, orthcheck(iwf,iwfd,isp)
  enddo
enddo
enddo
!!
!      do isp=1,nsp
!      do iwf  = 1,nwf           !dipoig is the rotation of dipoig1 by symops
!      do iwfd = 1,nwf
!         write(6,ftox)'dipole real=',iwf,iwfd,isp,ftof(dipole(:,iwf,iwfd,isp),6)
!      enddo
!      enddo
!      enddo
do iwf  = 1,nwf
do iwfd = 1,nwf
  write(6,ftox)'dddd hammso real =',iwf,iwfd,ftof(hammso(iwf,iwfd,1:3),8)
enddo
enddo

block
real(8):: hammso2(nwf,nwf),sss=0d0
hammso2 = matmul(hammso(:,:,1),hammso(:,:,1)) + matmul(hammso(:,:,2),hammso(:,:,2)) &
   + matmul(hammso(:,:,3),transpose(hammso(:,:,3)))+ matmul(transpose(hammso(:,:,3)),hammso(:,:,3))
do iwf  = 1,nwf
 sss=sss+ hammso2(iwf,iwf)
enddo
write(6,*)"hhhddd hammso**2 diagonal=",sss
endblock

print *,'hsocmat ok'
call rx0s('hsocmat: ok')
end subroutine hsocmat
