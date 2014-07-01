      program read_gw0

c-------------------------------------------------
c read  <psi(k,m) | HLDA | psi(k,n) >
c       <psi(k,m) | Vxc | psi(k,n)>
c       <psi(k,m) | SEex | psi(k,n)>
c       <psi(k,m) | SEcor(w) | psi(k,n)>
c |psi(k,m)> : Bloch fn. in the maxloc Wannier gauge
c |phi_i > : the i-th maxloc Wannier orbital centered at the cell R=0
c measured from the Fermi level 
c energies in eV (HLDA,SE,V,Wc)
c    
c Takashi Miyake, Sep.22, 2008, calc. interpolated bands
c Takashi Miyake, Jul.02, 2008, cutoff high energy part
c Takashi Miyake, Jan.22, 2008, interpolation to a fine k mesh
c Takashi Miyake, Oct.17, 2007, from read_gwt
c Takashi Miyake, Jun.12, 2007, added nel1 and nel2
c Takashi Miyake, Oct.28, 2006
c-------------------------------------------------
      implicit none
      complex(8),parameter:: img=(0d0,1d0)
      complex(8):: ctmp,ztmp,cegw
      real(8) :: eps ! broadening for A(w)
      real(8),parameter :: rydberg=13.6058d0
      real(8):: ef,alat,aa(3,3),bb(3,3),nel1,nel2,vn,pi,dwplot,dwev,
     &      rdel,cdel,omg,efgw,eflda,eqp,ds0,nel,nsum,
     &      nocc,nall,ec1,ec2,c1,c2,domg,omg0,beta,zfac,rtmp,rtmp2,
     &      rtmp3,rtmp4,xk,dk(3),egw,dsde,dwz,dwz2,zfac2
      integer(4):: nwf,nqbz,nqibz,nq1,nq2,nq3,ntau,nspin,iq,is,ik,ik2,
     &             iwini,iwend,ifin,niw,nrws,ii,jj,kk,nkbz,nk1,nk2,nk3,
     &             iw,iwef,iwcut,inef,ikef,nkef1,nkef2,nkef3,nbin,nloop,
     &             iwc1,iwc2,ibin,niv,iv,iflda,ifgw,ifegw,im,in,iwz2,
     &             iwa,iwb,iwgw,ifgwc
      complex(8),allocatable :: zw(:,:),zw0(:,:),zw2(:,:),zw3(:,:,:),
     &     zlda(:,:,:,:),zsc(:,:,:,:,:),zsx(:,:,:,:),zvxc(:,:,:,:)
     &    ,zlda2(:,:,:,:),zsc2(:,:,:,:),zsx2(:,:,:,:),zvxc2(:,:,:,:)
     &    ,zlda3(:,:),zsc3(:,:),zsx3(:,:),zvxc3(:,:)
     &    ,zse0(:,:),zwr(:,:,:),zwq(:,:,:),evecc(:,:)
     &    ,zs0(:,:),zsw(:,:,:),zgw(:,:,:,:),zeqp(:,:,:),zvec(:,:,:,:)
     &    ,zscd(:,:),zscd2(:,:,:),zgwc(:,:,:,:)
      real(8),allocatable :: qbz(:,:),kbz(:,:),tau(:),rws(:,:),drws(:)
     &    ,eval(:),aw(:,:),aw0(:,:),rw(:),znk(:,:,:)
      integer(4),allocatable :: irws(:),ikq(:),iqstar(:)
      logical :: lwin,lband,ldebug,lrshift
      integer(4):: nlinex
      parameter (nlinex=100)
      integer(4)::nline,np(nlinex)
      real(8):: qi(3,nlinex),qf(3,nlinex),p(3)
c-------------------------------------------------
c debug: >
!      open(99,file='kkrot',status='old')
!      read(99,*)nqbz,nqibz
!      allocate(ikq(nqbz),iqstar(nqibz))
!      iqstar=0
!      do ik=1,nqbz
!         read(99,*)ik2,ikq(ik)
!         if (ik2.ne.ik) stop 'read_sew2: kkrot error'
!      enddo   
!      close(99)

      pi = 4d0*datan(1d0)
      ifin = 901
      iflda = 902
      ifgw = 903
      ifegw= 904
      ifgwc = 905
cc debug:
!      open(ifin,file='GW0.diag',form='formatted')
      open(ifin,file='GW0',form='formatted')
      open(iflda,file='HLDA',form='formatted')
      open(ifgw,file='HGW',form='formatted')
      open(ifegw,file='EGW',form='formatted')
      open(ifgwc,file='HGWC',form='formatted')
! nwf: # "wannier" bands
! nqbz: # q points in FBZ
! n1,n2,n3: q-mesh (nqbz = n1*n2*n3)
! ef: LDA Fermi energy (in Ry) (irrelevant)
! alat: lattice constant in Bohr
! dwplot: frequency mesh in Ry
! omega(iw) = dwplot*13.6058*iw (eV), (iw=iwini,...iwend, measured from Ef)
      read(ifin,*)nwf,nqbz,nq1,nq2,nq3,iwini,iwend,dwplot,nspin,ef,alat
      dwev = dwplot*rydberg
      allocate(zw(nwf,nwf),zw0(nwf,nwf),zw2(nwf,nwf),qbz(3,nqbz),
     &     zlda(nwf,nwf,nqbz,nspin),zsc(nwf,nwf,nqbz,iwini:iwend,nspin),
     &     zsx(nwf,nwf,nqbz,nspin),zvxc(nwf,nwf,nqbz,nspin),
     &     zgw(nwf,nwf,nqbz,nspin),
     &     zgwc(nwf,nwf,nqbz,nspin),znk(nwf,nqbz,nspin),
     &     zeqp(nwf,nqbz,nspin),zvec(nwf,nwf,nqbz,nspin))
! aa: unit cell vectors in units of alat
      read(ifin,*)aa(1:3,1)
      read(ifin,*)aa(1:3,2)
      read(ifin,*)aa(1:3,3)
! bb: reciprocal vectors in units of 2pi/alat
      read(ifin,*)bb(1:3,1)
      read(ifin,*)bb(1:3,2)
      read(ifin,*)bb(1:3,3)
      read(ifin,*)qbz(:,:) ! in units of 2pi/alat
! Wigner-Seitz supercell
      allocate(irws(nq1*nq2*nq3*8),rws(3,nq1*nq2*nq3*8),
     &         drws(nq1*nq2*nq3*8))
      call wigner_seitz(alat,aa,nq1,nq2,nq3,nrws,rws,irws,drws)
! standard input
      read(*,*)nk1,nk2,nk3 ! fine k-mesh
      nkbz = nk1*nk2*nk3
      allocate(kbz(3,nkbz))
      call fine_mesh(nk1,nk2,nk3,bb,kbz)
      read(*,*)nkef1,nkef2,nkef3 ! k-mesh for LDA Fermi level
      read(*,*)nel ! number of electrons
      read(*,*)eps ! broadening
      read(*,*)ec1,ec2 ! energy window
      read(*,*)nbin
      read(*,*)beta ! inverse temperature in eV^{-1}
      read(*,*)niv ! number of Matsubara frequencies
      read(*,*)dwz ! energy mesh for Z
      iwz2 = int(0.5d0*dwz/dwev) + 1
      dwz2 = dble(2d0*iwz2)*dwev      
!
      allocate(zlda2(nwf,nwf,nkbz,nspin),
     &         zsc2(nwf,nwf,nkbz,nspin),
     &         zsx2(nwf,nwf,nkbz,nspin),zvxc2(nwf,nwf,nkbz,nspin),
     &         zlda3(nwf,nwf),zsc3(nwf,nwf),
     &         zsx3(nwf,nwf),zvxc3(nwf,nwf),
     &         zwr(nwf,nwf,nrws),zwq(nwf,nwf,nqbz),
     &         aw(iwini:iwend,nspin),aw0(iwini:iwend,nspin))
      allocate(eval(nwf),evecc(nwf,nwf))
      zlda2 = 0d0
      zvxc2 = 0d0
      zsx2 = 0d0
      zsc2 = 0d0
      aw  = 0d0
      aw0  = 0d0

      do is = 1,nspin
         do ik = 1,nqbz
            read(ifin,*)zlda(:,:,ik,is) ! HLDA (eV)
!            call diag_hm(zlda(:,:,ik,is),nwf,eval,evecc)
            read(ifin,*)zvxc(:,:,ik,is) ! Vxc (eV)
c Wannier gauge to Hamiltonian gauge: Vxc
!            call wan2ham(zvxc(:,:,ik,is),evecc,nwf,zw)
            read(ifin,*)zsx(:,:,ik,is) ! SEex (eV)
c Wannier gauge to Hamiltonian gauge: Sx
!            call wan2ham(zsx(:,:,ik,is),evecc,nwf,zw)
            do iw = iwini,iwend
              read(ifin,*)zsc(:,:,ik,iw,is) ! SEcor (eV)
            enddo ! iw
cc debug:
            zsc(:,:,ik,iw,is) = zsc(:,:,ik,iw,is) - 
     &                img* imag(zsc(:,:,ik,0,is))
c Wannier gauge to Hamiltonian gauge: Sc(w)
!            if (ik.eq.1) then
!            do iw=iwini,iwend
!              call wan2ham(zsc(:,:,ik,iw,is),evecc,nwf,zw)
!              write(11,*)iw,real(zw(22,22)),imag(zw(22,22))
!            enddo
!            endif
         enddo ! ik
       enddo ! is  
! end of reading data

! construct GW Hamiltonian
       zgw = 0d0
       zgwc = 0d0
       zeqp = 0d0
       zvec = 0d0
       do is = 1,nspin
! determine Ef(GW) = E(LDA)+ SE(w=0) for the highest occupied LDA level
         if (is.ne.1) stop 'spin polarized case not implemented'

         inquire(file='ReadShift',exist=lrshift)
         if (lrshift) then
            write(*,*)'read ReadShift'
            open(98,file='ReadShift')
            read(98,*)eflda,efgw,ds0
            close(98)
         else ! lrshift
         call calc_ef(zlda(:,:,:,is),zvxc(:,:,:,is),zsx(:,:,:,is),
     &             zsc(:,:,:,:,is),nel,
     &             irws,rws,qbz,bb,dwev,
     &             nkef1,nkef2,nkef3,
     &             nwf,nqbz,nrws,iwini,iwend,
     &             eflda,efgw,ds0)
         write(*,*)'eflda=',eflda
         write(*,*)'efgw=',efgw
         write(*,*)'ds0=',ds0
         open(98,file='ReadShift.tmp')
         write(98,*)eflda,efgw,ds0
         close(98)
         endif ! lrshift
! shift the Fermi level
         do ii = 1,nwf
            zlda(ii,ii,:,:) = zlda(ii,ii,:,:) - eflda ! shift HLDA so that EF(LDA) = 0
            zsx(ii,ii,:,:)  = zsx(ii,ii,:,:)  - ds0   ! shift SE(exchange) so that EF(GW) = 0
         enddo   

         allocate(zs0(nwf,nwf),zsw(nwf,nwf,4),zscd(nwf,nwf))
         do ik = 1,nqbz
           call diag_hm(zlda(:,:,ik,is),nwf,eval,evecc)
           zvec(:,:,ik,is) = evecc
           zw = zsx(:,:,ik,is) - zvxc(:,:,ik,is)
           call wan2ham(zw,evecc,nwf,zs0)
           write(99,*)
! check Self-energy
           if (ik.eq.1) then
           do iw=iwini,iwend
             call wan2ham(zsc(:,:,ik,iw,is),evecc,nwf,zsw(:,:,1))
             zsw(:,:,1) = zsw(:,:,1) + zs0
             omg=dwev*dble(iw)
             write(71,"(7f12.6)")omg,(real(zsw(ii,ii,1)),ii=1,6)
             write(72,"(7f12.6)")omg,(imag(zsw(ii,ii,1)),ii=1,6)
           enddo   
           endif

           do kk = 1,nwf
!! quasiparticle energy
              iw = int(eval(kk)/dwev)
              if (eval(kk).le.0d0) iw = iw - 1
              call wan2ham(zsc(:,:,ik,iw,is),evecc,nwf,zsw(:,:,1))
              call wan2ham(zsc(:,:,ik,iw+1,is),evecc,nwf,zsw(:,:,2))
              call wan2ham(zsc(:,:,ik,0,is),evecc,nwf,zsw(:,:,3))
! real part
!              zfac = (real(zsw(kk,kk,2))-real(zsw(kk,kk,1)))/dwev
              iwa = iw-iwz2
              if (iwa.lt.iwini) stop 'iwa too low'
              iwb = iw+iwz2
              if (iwb.gt.iwend) stop 'iwb too high'
              zscd = zsc(:,:,ik,iwb,is)-zsc(:,:,ik,iwa,is)
              call wan2ham(zscd,evecc,nwf,zsw(:,:,4))
              zfac = (real(zsw(kk,kk,4)))/(dble(iwb-iwa)*dwev)
              zfac = 1.d0/(1.0d0-zfac)
              rtmp = (eval(kk) - dwev*dble(iw))/dwev
              rtmp2 = real(zsw(kk,kk,2))*rtmp
     &              + real(zsw(kk,kk,1))*(1.0d0-rtmp)
              egw = eval(kk) + zfac*(real(zs0(kk,kk))+rtmp2)
              znk(kk,ik,is) = zfac
!              egw = eval(kk) ! for debug
! imaginary part
              iw = int(egw/dwev)
!              iwa = iw-iwz2
!              if (iwa.lt.iwini) stop 'iwa too low'
!              iwb = iw+iwz2
!              if (iwb.gt.iwend) stop 'iwb too high'
!              zscd = zsc(:,:,ik,iwb,is)-zsc(:,:,ik,iwa,is)
              call wan2ham_d(zsc(:,:,ik,iw,is),evecc,nwf,zsw(:,:,1))
              call wan2ham_d(zsc(:,:,ik,iw+1,is),evecc,nwf,zsw(:,:,2))
              call wan2ham_d(zsc(:,:,ik,0,is),evecc,nwf,zsw(:,:,3))
!              call wan2ham_d(zscd,evecc,nwf,zsw(:,:,4))
!              dsde = (imag(zsw(kk,kk,2))-imag(zsw(kk,kk,1)))/dwev
!              dsde = (imag(zsw(kk,kk,4)))/(dble(iwb-iwa)*dwev)
!              zfac2 = 1.d0/(1.0d0-dsde)
              rtmp = (egw - dwev*dble(iw))/dwev
              rtmp2 = imag(zsw(kk,kk,2))*rtmp
     &              + imag(zsw(kk,kk,1))*(1.0d0-rtmp)
     &              - imag(zsw(kk,kk,3))
c              rtmp3 = dsde*(egw-eval(kk))
c              rtmp4 = rtmp2 + rtmp3
c              cegw = egw + img*rtmp4
              cegw = egw + img*rtmp2
              zeqp(kk,ik,is) = cegw
! imaginary part, z factor at E(LDA)
              iw = int(eval(kk)/dwev)
              iwa = iw-iwz2
              if (iwa.lt.iwini) stop 'iwa too low'
              iwb = iw+iwz2
              if (iwb.gt.iwend) stop 'iwb too high'
              zscd = zsc(:,:,ik,iwb,is)-zsc(:,:,ik,iwa,is)
              call wan2ham_d(zscd,evecc,nwf,zsw(:,:,4))
!              dsde = (imag(zsw(kk,kk,2))-imag(zsw(kk,kk,1)))/dwev
              dsde = (imag(zsw(kk,kk,4)))/(dble(iwb-iwa)*dwev)
c debug: >>
c              write(81,"(6f12.6)")eval(kk),egw,zfac,rtmp2,rtmp3,rtmp4
              write(81,"(4f12.6)")eval(kk),egw,zfac,rtmp2
              write(82,"(4f12.6)")eval(kk),egw,zfac,-dsde
c debug: <<
              do im = 1,nwf
              do in = 1,nwf
                zgw(im,in,ik,is) = zgw(im,in,ik,is)
     &        + evecc(im,kk)*egw*dconjg(evecc(in,kk))
                zgwc(im,in,ik,is) = zgwc(im,in,ik,is)
     &        + evecc(im,kk)*cegw*dconjg(evecc(in,kk))
              enddo ! in
              enddo ! im
           enddo ! kk
         enddo ! ik
         enddo ! is
         close(99)
         deallocate(zs0,zsw,zscd)

! output Hamiltonian on mesh points
         write(iflda,*)nspin,nqbz,nwf
         write(iflda,*)qbz
         write(iflda,*)zlda
         write(ifgw,*)nspin,nqbz,nwf
         write(ifgw,*)qbz
         write(ifgw,*)zgw
         write(ifegw,*)nspin,nqbz,nwf
         write(ifegw,*)qbz
         write(ifegw,*)zeqp
         write(ifegw,*)zvec
         write(ifegw,*)znk
         write(ifgwc,*)nspin,nqbz,nwf
         write(ifgwc,*)qbz
         write(ifgwc,*)zgwc
         write(ifgwc,*)znk
!         stop 'output HLDA and HGW done'

! bands
       inquire(file='SYML',exist=lband)
       if (lband) then
         open(99,file='SYML',status='old')
         do ii = 1,nlinex
           read(99,*)np(ii),qi(1:3,ii),qf(1:3,ii)
           if (np(ii).eq.0) exit
           nline = ii
         enddo ! ii
         if (nline.eq.nlinex) stop 'read_gw0: too many lines in SYML'
         close(99)
       endif

       do is = 1,nspin
! determine Ef(GW) = E(LDA)+ SE(w=0) for the highest occupied LDA level
         if (is.ne.1) stop 'spin polarized case not implemented'
!!         call calc_ef(zlda(:,:,:,is),zvxc(:,:,:,is),zsx(:,:,:,is),
!!     &             zsc(:,:,:,:,is),nel,
!!     &             irws,rws,qbz,bb,dwev,
!!     &             nkef1,nkef2,nkef3,
!!     &             nwf,nqbz,nrws,iwini,iwend,
!!     &             eflda,efgw,ds0)
!!         write(*,*)'eflda=',eflda
!!         write(*,*)'efgw=',efgw
!!         write(*,*)'ds0=',ds0
!!! shift HLDA so that EF(LDA) = 0
!!         do ii = 1,nwf
!!            zlda(ii,ii,:,:) = zlda(ii,ii,:,:) - eflda
!!         enddo   
!!! shift SE(exchange) so that EF(GW) = 0
!!         do ii = 1,nwf
!!            zsx(ii,ii,:,:) = zsx(ii,ii,:,:) - ds0
!!         enddo   

! band (linearization, diagonal only)
         open(99,file='qpbnd')
         allocate(zs0(nwf,nwf),zsw(nwf,nwf,4),zscd2(nwf,nwf,nqbz))
         xk = 0d0
         do ii = 1,nline
         do jj = 0,np(ii)   
           dk(:) = (qf(:,ii)-qi(:,ii))/dble(np(ii))
           if (jj.gt.0) xk = xk + dsqrt(sum(dk(:)**2))
           p(:) = qi(:,ii) + dk(:)*dble(jj)
           write(*,"('k=',3f12.6)")p
           call q2p(irws,rws,qbz,p,nrws,nqbz,nwf,zlda(:,:,:,is),zlda3)
           call q2p(irws,rws,qbz,p,nrws,nqbz,nwf,zvxc(:,:,:,is),zvxc3)
           call q2p(irws,rws,qbz,p,nrws,nqbz,nwf,zsx(:,:,:,is),zsx3)
           call diag_hm(zlda3,nwf,eval,evecc)
           zw = zsx3 - zvxc3
           call wan2ham(zw,evecc,nwf,zs0)
           write(99,*)
           do kk = 1,nwf
              iw = int(eval(kk)/dwev)
              if (eval(kk).le.0d0) iw = iw - 1
              call q2p(irws,rws,qbz,p,nrws,nqbz,nwf,zsc(:,:,:,iw,is),zw)
              call wan2ham(zw,evecc,nwf,zsw(:,:,1))
              call q2p(irws,rws,qbz,p,nrws,nqbz,nwf,zsc(:,:,:,iw+1,is),zw)
              call wan2ham(zw,evecc,nwf,zsw(:,:,2))
              call q2p(irws,rws,qbz,p,nrws,nqbz,nwf,zsc(:,:,:,0,is),zw)
              call wan2ham(zw,evecc,nwf,zsw(:,:,3))
! real part
!              zfac = (real(zsw(kk,kk,2))-real(zsw(kk,kk,1)))/dwev
              iwa = iw-iwz2
              if (iwa.lt.iwini) stop 'iwa too low'
              iwb = iw+iwz2
              if (iwb.gt.iwend) stop 'iwb too high'
              zscd2 = zsc(:,:,:,iwb,is)-zsc(:,:,:,iwa,is)
              call q2p(irws,rws,qbz,p,nrws,nqbz,nwf,zscd2,zw)
              call wan2ham(zw,evecc,nwf,zsw(:,:,4))
              zfac = (real(zsw(kk,kk,4)))/(dble(iwb-iwa)*dwev)
              zfac = 1.d0/(1.0d0-zfac)
              rtmp = (eval(kk) - dwev*dble(iw))/dwev
              rtmp2 = real(zsw(kk,kk,2))*rtmp
     &              + real(zsw(kk,kk,1))*(1.0d0-rtmp)
              egw = eval(kk) + zfac*(real(zs0(kk,kk))+rtmp2)
!              egw = eval(kk) + (real(zs0(kk,kk)))
              write(99,"(f12.6,i5,3f12.6)")xk,kk,eval(kk),egw,zfac
           enddo ! kk
         enddo ! jj
         enddo ! ii
         close(99)
         deallocate(zs0,zsw,zscd2)

         stop 'debug: before fine-mesh'

! interpolation to a fine k mesh
         call q2k(irws,rws,qbz,kbz,nrws,nqbz,nkbz,nwf,zlda(:,:,:,is),
     &             zlda2(:,:,:,is))
         call q2k(irws,rws,qbz,kbz,nrws,nqbz,nkbz,nwf,zvxc(:,:,:,is),
     &             zvxc2(:,:,:,is))
         call q2k(irws,rws,qbz,kbz,nrws,nqbz,nkbz,nwf,zsx(:,:,:,is),
     &             zsx2(:,:,:,is))

! calc A(w) and nsum(w) = I[:w] A(w)
         open(777,file='AW')
         nsum = 0.0d0
         iwc1 = ec1 / dwev - 1
         iwc2 = ec2 / dwev + 1
         do iw = iwini,iwend
c         do iw = 0,iwini,-1
! check if inside energy window
           lwin = .false.
           nloop = 1
           if (iw.ge.iwc1.and.iw.le.iwc2) then
              lwin = .true.
              nloop = nbin
           endif   
           do ibin = 1,nloop
              zwq = zsc(:,:,:,iw,is)
              omg0 = dwev*dble(iw)
              omg = omg0
              if (lwin) then
                 domg = -0.5d0*(1.0d0+1.0d0/dble(nbin))
     &                         + dble(ibin)/dble(nbin)
                 omg = omg0 + dwev*domg
                 if (domg.lt.0.0d0) then
                    c1 = abs(domg)*2.0d0
                    c2 = 1.0d0 - c1
                    zwq = c1*zsc(:,:,:,iw-1,is) + c2*zsc(:,:,:,iw,is)
                 else
                    c1 = abs(domg)*2.0d0
                    c2 = 1.0d0 - c1
                    zwq = c1*zsc(:,:,:,iw+1,is) + c2*zsc(:,:,:,iw,is)
                 endif   
              endif   
! -(H+SE)
           call q2k(irws,rws,qbz,kbz,nrws,nqbz,nkbz,nwf,zwq
     &             ,zsc2(:,:,:,is))
           do ik = 1,nkbz
              zw0 = - zlda2(:,:,ik,is)
              zw  = - zlda2(:,:,ik,is)
     &        - (zsx2(:,:,ik,is)+zsc2(:,:,ik,is)-zvxc2(:,:,ik,is))
! omega-(H+SE)
             if (omg.le.0.d0) then
               do ii=1,nwf
                 zw0(ii,ii) = zw0(ii,ii) + omg - img*eps
                 zw(ii,ii) = zw(ii,ii) + omg - img*eps
               enddo ! ii
             else
               do ii=1,nwf
                 zw0(ii,ii) = zw0(ii,ii) + omg + img*eps
                 zw(ii,ii) = zw(ii,ii) + omg + img*eps
               enddo ! ii
             endif   ! omg
! spectral function
             call inv_mat(zw0,nwf,zw2)
             do ii=1,nwf
              aw0(iw,is) = aw0(iw,is)+abs(imag(zw2(ii,ii)))/dble(nloop)
             enddo
             call inv_mat(zw,nwf,zw2)
             do ii=1,nwf
               aw(iw,is) = aw(iw,is) + abs(imag(zw2(ii,ii)))/dble(nloop)
             enddo
           enddo ! ik   
           enddo ! ibin
           aw0(iw,is)=aw0(iw,is)/(pi*dble(nkbz))
           aw(iw,is)=aw(iw,is)/(pi*dble(nkbz))
           if (nspin.eq.1) aw0(iw,is) = aw0(iw,is)*2.0d0
           if (nspin.eq.1) aw(iw,is) = aw(iw,is)*2.0d0
           nsum = nsum + aw(iw,is)*dwev
           write(777,"(i5,4f16.7)")iw,omg0,aw0(iw,is),aw(iw,is),nsum
         enddo ! iw
         close(777)
         stop 'debug: after A(w)'

c real w to Matsubara iv
         allocate(zw3(nwf,nwf,nkbz))
         do iv = -niv,niv-1
            vn = dble(2*iv+1) * pi / beta
            zw3 = 0.0d0
            do iw = iwc1,iwc2
            do ibin = 1,nbin
              omg0 = dwev*dble(iw)
              domg = -0.5d0*(1.0d0+1.0d0/dble(nbin))
     &                      + dble(ibin)/dble(nbin)
              omg = omg0 + dwev*domg
              if (domg.lt.0.0d0) then
                 c1 = abs(domg)*2.0d0
                 c2 = 1.0d0 - c1
                 zwq = c1*zsc(:,:,:,iw-1,is) + c2*zsc(:,:,:,iw,is)
              else
                 c1 = abs(domg)*2.0d0
                 c2 = 1.0d0 - c1
                 zwq = c1*zsc(:,:,:,iw+1,is) + c2*zsc(:,:,:,iw,is)
              endif
              call q2k(irws,rws,qbz,kbz,nrws,nqbz,nkbz,nwf,zwq
     &             ,zsc2(:,:,:,is))
              ctmp = 1.0d0/(img*vn-omg)*(dwev/dble(nbin))/pi
              do ik = 1,nkbz
              do jj = 1,nwf
              do ii = 1,nwf
!                 zw3(ii,jj,ik) = zw3(ii,jj,ik) + ctmp*abs(zwq(ii,jj,ik))
                 zw3(ii,jj,ik) = zw3(ii,jj,ik) + ctmp*imag(zsc2(ii,jj,ik,is))
              enddo ! ii
              enddo ! jj
              enddo ! ik
            enddo ! ibin
            enddo ! iw   
         enddo ! iv   

c divide into low-energy and high-energy parts
         do iw = iwini,iwend
           omg = dwev*dble(iw)
         enddo ! iw
      enddo ! is

      end
c--------------------------------------------------------------------
      subroutine q2k(irws,rws,qbz,kbz,nrws,nqbz,nkbz,nwf,zin,zout)
      implicit none

      integer(4):: irws(nrws),nrws,nqbz,nkbz,nwf,iq,ik,ir,ii,jj
      real(8):: rws(3,nrws),qbz(3,nqbz),kbz(3,nkbz),pi,rk
      complex(8):: zin(nwf,nwf,nqbz),zout(nwf,nwf,nkbz),zw(nwf,nwf,nrws)
     &            ,ci,ceikr

      ci = (0d0,1d0)
      pi = 4d0*datan(1.0d0)

      zw = 0d0
      do iq = 1,nqbz
      do ir = 1,nrws
         rk = 2d0*pi*sum(rws(:,ir)*qbz(:,iq))
         ceikr = exp(-ci*rk) / dble(nqbz)
         zw(:,:,ir) = zw(:,:,ir) + ceikr * zin(:,:,iq)
      enddo ! ir
      enddo ! iq

      zout = 0d0
      do ik = 1,nkbz
      do ir = 1,nrws
         rk = 2d0*pi*sum(rws(:,ir)*kbz(:,ik))
         ceikr = exp(ci*rk) / dble(irws(ir))
         zout(:,:,ik) = zout(:,:,ik) + ceikr * zw(:,:,ir)
      enddo ! ir
      enddo ! ik

      return
      end
c-----------------------------------------------------------------------
      subroutine q2p(irws,rws,qbz,p,nrws,nqbz,nwf,zin,zout)
      implicit none

      integer(4):: irws(nrws),nrws,nqbz,nwf,iq,ik,ir,ii,jj
      real(8):: rws(3,nrws),qbz(3,nqbz),p(3),pi,rk
      complex(8):: zin(nwf,nwf,nqbz),zout(nwf,nwf),zw(nwf,nwf,nrws)
     &            ,ci,ceikr

      ci = (0d0,1d0)
      pi = 4d0*datan(1.0d0)

      zw = 0d0
      do iq = 1,nqbz
      do ir = 1,nrws
         rk = 2d0*pi*sum(rws(:,ir)*qbz(:,iq))
         ceikr = exp(-ci*rk) / dble(nqbz)
         zw(:,:,ir) = zw(:,:,ir) + ceikr * zin(:,:,iq)
      enddo ! ir
      enddo ! iq

      zout = 0d0
      do ir = 1,nrws
         rk = 2d0*pi*sum(rws(:,ir)*p(:))
         ceikr = exp(ci*rk) / dble(irws(ir))
         zout(:,:) = zout(:,:) + ceikr * zw(:,:,ir)
      enddo ! ir

      return
      end
c-----------------------------------------------------------------------
      subroutine q2r(irws,rws,qbz,nrws,nqbz,nwf,zin,zout)
      implicit none

      integer(4):: irws(nrws),nrws,nqbz,nwf,iq,ik,ir,ii,jj
      real(8):: rws(3,nrws),qbz(3,nqbz),pi,rk
      complex(8):: zin(nwf,nwf,nqbz),zout(nwf,nwf,nrws)
     &            ,ci,ceikr

      ci = (0d0,1d0)
      pi = 4d0*datan(1.0d0)

      zout = 0d0
      do iq = 1,nqbz
      do ir = 1,nrws
         rk = 2d0*pi*sum(rws(:,ir)*qbz(:,iq))
         ceikr = exp(-ci*rk) / dble(nqbz)
         zout(:,:,ir) = zout(:,:,ir) + ceikr * zin(:,:,iq)
      enddo ! ir
      enddo ! iq

      return
      end
c-----------------------------------------------------------------------
      subroutine r2k(irws,rws,kbz,nrws,nwf,zin,zout)
      implicit none

      integer(4):: irws(nrws),nrws,nwf,iq,ik,ir,ii,jj
      real(8):: rws(3,nrws),kbz(3),pi,rk
      complex(8):: zin(nwf,nwf,nrws),zout(nwf,nwf)
     &            ,ci,ceikr

      ci = (0d0,1d0)
      pi = 4d0*datan(1.0d0)

      zout = 0d0
      do ir = 1,nrws
         rk = 2d0*pi*sum(rws(:,ir)*kbz(:))
         ceikr = exp(ci*rk) / dble(irws(ir))
         zout(:,:) = zout(:,:) + ceikr * zin(:,:,ir)
      enddo ! ir

      return
      end
c-----------------------------------------------------------------------
      subroutine wigner_seitz(alat,plat,n1,n2,n3,nrws,rws,irws,drws)
      implicit real*8(a-h,o-z)

      integer(4):: n1,n2,n3,nrws
      real(8):: alat,plat(3,3)
      integer(4) :: irws(n1*n2*n3*8)
      real(8) :: rws(3,n1*n2*n3*8),drws(n1*n2*n3*8)

      integer(4):: ii0(3,8),isort(8),
     &             iwork1(n1*n2*n3*8),iwork2(n1*n2*n3*8)
      real(8):: rr(3,8),dd(8)
      parameter (tol=1.d-6)


      nrws = 0
      do i1=0,n1-1
      do i2=0,n2-1
      do i3=0,n3-1

         n = 0
         do j1=0,1
         do j2=0,1
         do j3=0,1
            n = n+1
            ii0(1,n) = i1 - j1*n1
            ii0(2,n) = i2 - j2*n2
            ii0(3,n) = i3 - j3*n3
         enddo ! j3
         enddo ! j2
         enddo ! j1

         do n=1,8
            rr(1:3,n) =  ( plat(1:3,1)*dble(ii0(1,n))
     &                 +   plat(1:3,2)*dble(ii0(2,n))
     &                 +   plat(1:3,3)*dble(ii0(3,n)) )
         enddo
         call sortvec2(8,rr,dd,isort)
         ndegen = 1
         do n=2,8
            if ((dd(n)-dd(1)).le.tol) ndegen = n
         enddo

         do n=1,ndegen
            nrws = nrws + 1
            rws(1:3,nrws) = rr(1:3,n)
            drws(nrws) = dd(n)
            irws(nrws) = ndegen
         enddo   

      enddo ! i3
      enddo ! i2
      enddo ! i1

      call sortvec2(nrws,rws,drws,iwork1)
      iwork2(1:nrws) = irws(1:nrws)
      do n=1,nrws
         irws(n) = iwork2(iwork1(n))
      enddo   

      return
      end
c-----------------------------------------------------------------------
      subroutine sortr(a,n,isort)
      implicit real*8(a-h,o-z)
      real(8) :: eps=1.0d-6
      real (8) :: a(n),b(n)
      integer (4) :: isort(n)

      b = a
      do i = 1,n
         isort(i) = i
      enddo   

      do j = 2,n
         c = b(j)
         do i = j-1,1,-1
            if (b(i)<=c) goto 999
            b(i+1) = b(i)
            isort(i+1) = isort(i)
         enddo
         i = 0
 999     b(i+1) = c
         isort(i+1) = j
      enddo

      do i = 1,n-1
         if (b(i).gt.b(i+1)) stop 'sortr: sorting error!'
      enddo   
      do i = 1,n
         if (abs(b(i)-a(isort(i))).gt.eps) stop 'sortr: sorting error!'
      enddo   

      return
      end
c-----------------------------------------------------------------------
      subroutine  sortvec2(ndat,vec,dist,idat)

      implicit real*8(a-h,o-z)
      real (8) :: vec(3,ndat),vtmp(3,ndat),dist(ndat)
      integer (4) :: idat(ndat)

      vtmp = vec
      do i = 1,ndat
         dist(i) = dsqrt(sum(vtmp(:,i)**2))
         idat(i) = i
      enddo

      do j = 2,ndat
         d = dist(j)
         do i = j-1,1,-1
            if (dist(i)<=d) goto 999
            dist(i+1) = dist(i)
            idat(i+1) = idat(i)
         enddo
         i = 0
 999     dist(i+1) = d
         idat(i+1) = j
      enddo

      do i = 1,ndat
         vec(1:3,i) = vtmp(1:3,idat(i))
      enddo   

      do i = 1,ndat-1
         d1 = dsqrt(sum(vec(:,i)**2))
         d2 = dsqrt(sum(vec(:,i+1)**2))
         if (d1 .gt. d2) stop 'sortvec: sorting error!'
         if (abs(d1-dist(i)).gt.1.d-4)
     &   stop 'sortvec: sorting error in d!'
      enddo   

      return
      end
c-----------------------------------------------------------------------
      subroutine diag_hm(zmat,ndim,eval,evecc)
      implicit real*8(a-h,o-z)

      complex(8),allocatable :: zmat2(:,:),ovlpc(:,:)
      complex(8):: zmat(ndim,ndim),evecc(ndim,ndim)
      real(8):: eval(ndim),wk(ndim,11)
      integer iwk(ndim)

      allocate(zmat2(ndim,ndim),ovlpc(ndim,ndim))

      nev  = ndim
      nmx  = ndim

      zmat2 = zmat

      ovlpc = (0d0,0d0)
      do i=1,ndim
         ovlpc(i,i) = (1d0,0d0)
      enddo

      evecc = (0d0,0d0)
      eval = 0d0

!      call diagno(ndim,zmat2,ovlpc,wk,iwk,evecc,eval)
      call diagcv(ovlpc,zmat2, evecc, ndim, eval, nmx, 1d99, nev)

      deallocate(zmat2,ovlpc)

      return
      end
c-----------------------------------------------------------------------
      subroutine wan2ham(zin,evecc,ndim,zout)
      implicit real*8(a-h,o-z)

      integer ii,jj
      complex(8) :: zin(ndim,ndim),evecc(ndim,ndim),zout(ndim,ndim),
     &              zw(ndim,ndim)

      zw = matmul(zin,evecc)
      do jj=1,ndim
      do ii=1,ndim
         zout(ii,jj)=sum(dconjg(evecc(:,ii))*zw(:,jj))
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine wan2ham_d(zin,evecc,ndim,zout)
      implicit real*8(a-h,o-z)

      integer ii,jj
      complex(8) :: zin(ndim,ndim),evecc(ndim,ndim),zout(ndim,ndim),
     &              zw(ndim,ndim),zz(ndim,ndim)

      zz = zin
      do jj=1,ndim
      do ii=1,ndim
         if (ii.ne.jj)zz(ii,jj) = 0d0
      enddo   
      enddo   

      zw = matmul(zz,evecc)
      do jj=1,ndim
      do ii=1,ndim
         zout(ii,jj)=sum(dconjg(evecc(:,ii))*zw(:,jj))
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine inv_mat(zin,ndim,zout)
      implicit real*8(a-h,o-z)

      complex(8) :: zin(ndim,ndim),zout(ndim,ndim)
      complex(8) :: zw(ndim,ndim),work(ndim)
      integer :: ipiv(ndim)

      zw = zin
      call zgetrf(ndim,ndim,zw,ndim,ipiv,info)
      if (info.ne.0) stop 'zgetrf error'
      call zgetri(ndim,zw,ndim,ipiv,work,ndim,info)
      if (info.ne.0) stop 'zgetri error'
      zout = zw

      return
      end
c-----------------------------------------------------------------------
      subroutine fine_mesh(nk1,nk2,nk3,bb,kbz)
      implicit real*8(a-h,o-z)

      real(8):: bb(3,3),kbz(3,nk1*nk2*nk3)
      real(8):: k1,k2,k3

      ik = 0
      do i3 = 1,nk3
      do i2 = 1,nk2
      do i1 = 1,nk1
         ik = ik + 1
         k1 = dble(i1-1) / dble(nk1)
         k2 = dble(i2-1) / dble(nk2)
         k3 = dble(i3-1) / dble(nk3)
         do ii = 1,3
           kbz(ii,ik) = bb(ii,1)*k1 + bb(ii,2)*k2 + bb(ii,3)*k3
         enddo ! ii
      enddo ! i1
      enddo ! i2
      enddo ! i3
      if (ik.ne.nk1*nk2*nk3) stop 'read_gwt2: nkbz error'

      return
      end
c-----------------------------------------------------------------------
      subroutine calc_ef(zlda,zvxc,zsx,zsc,nel,
     &             irws,rws,qbz,bb,dwev,
     &             n1,n2,n3,
     &             nwf,nqbz,nrws,iwini,iwend,
     &             eflda,efgw,ds0)
      implicit real*8(a-h,o-z)

      integer(4):: n1,n2,n3
      complex(8):: zlda(nwf,nwf,nqbz),zvxc(nwf,nwf,nqbz),
     &             zsx(nwf,nwf,nqbz),zsc(nwf,nwf,nqbz,iwini:iwend)
      real(8):: rws(3,nrws),qbz(3,nqbz),bb(3,3),dwev,nel
      integer(4):: irws(nrws)
c local
      complex(8):: zhr(nwf,nwf,nrws),zvr(nwf,nwf,nrws),
     &             zxr(nwf,nwf,nrws),zcr(nwf,nwf,nrws),
     &             evecc(nwf,nwf)
      complex(8):: zhkf(nwf,nwf),zvkf(nwf,nwf),
     &             zxkf(nwf,nwf),zckf(nwf,nwf,3)
      complex(8):: zw1(nwf,nwf),zw2(nwf,nwf),zw3(nwf,nwf)
      real(8):: k1,k2,k3,q(3),kf(3),eval(nwf),eig(nwf,n1*n2*n3),
     &        nocc,nall
      integer(4):: isort(nwf*n1*n2*n3)

! compute Kohn-Sham levels for the n1*n2*n3 mesh
      call q2r(irws,rws,qbz,nrws,nqbz,nwf,zlda,zhr)
      ik = 0
      do i3 = 1,n3
      do i2 = 1,n2
      do i1 = 1,n1
         ik = ik + 1
         k1 = dble(i1-1) / dble(n1)
         k2 = dble(i2-1) / dble(n2)
         k3 = dble(i3-1) / dble(n3)
         do ii = 1,3
           q(ii) = bb(ii,1)*k1 + bb(ii,2)*k2 + bb(ii,3)*k3
         enddo ! ii
         call r2k(irws,rws,q,nrws,nwf,zhr,zw1)
         call diag_hm(zw1,nwf,eval,evecc)
         eig(:,ik) = eval
      enddo ! i1
      enddo ! i2
      enddo ! i3

! the highest LDA level
      call sortr(eig,nwf*n1*n2*n3,isort)
      ink = nint(nel*n1*n2*n3/2.0d0) ! assuming spin non-polarized
      j = isort(ink)
      iwf = mod(j,nwf)
      if (iwf==0) iwf = nwf
      ik = (j-iwf) / nwf + 1
      i1 = mod(ik,n1)
      if (i1==0) i1 = n1
      i2 = mod((ik-i1)/n1,n2) + 1
      i3 = (ik-i1-(i2-1)*n1)/(n1*n2) + 1
      eflda = eig(iwf,ik)

! check occupation number
      nocc = 0.0d0
      do i3 = 1,n3
      do i2 = 1,n2
      do i1 = 1,n1
         ik = ik + 1
         k1 = dble(i1-1) / dble(n1)
         k2 = dble(i2-1) / dble(n2)
         k3 = dble(i3-1) / dble(n3)
         do ii = 1,3
           q(ii) = bb(ii,1)*k1 + bb(ii,2)*k2 + bb(ii,3)*k3
         enddo ! ii
         call r2k(irws,rws,q,nrws,nwf,zhr,zw1)
         call diag_hm(zw1,nwf,eval,evecc)
         do ii = 1,nwf
            if (eval(ii).le.eflda) nocc = nocc + 1.0d0
         enddo
      enddo ! i1
      enddo ! i2
      enddo ! i3
      write(*,*)'nocc =',2.0d0*nocc/(n1*n2*n3) ! spin non-polarizaed

! Fermi wave vector
      k1 = dble(i1-1) / dble(n1)
      k2 = dble(i2-1) / dble(n2)
      k3 = dble(i3-1) / dble(n3)
      do ii = 1,3
        kf(ii) = bb(ii,1)*k1 + bb(ii,2)*k2 + bb(ii,3)*k3
      enddo ! ii
! HLDA,SEx,SEc,Vxc at kf
      call r2k(irws,rws,kf,nrws,nwf,zhr,zhkf) ! HLDA(kf)
      call q2r(irws,rws,qbz,nrws,nqbz,nwf,zsx,zxr)
      call r2k(irws,rws,kf,nrws,nwf,zxr,zxkf) ! SEx(kf)
      call q2r(irws,rws,qbz,nrws,nqbz,nwf,zvxc,zvr)
      call r2k(irws,rws,kf,nrws,nwf,zvr,zvkf) ! Vxc(kf)
!      iw = eflda / dwev
      iw = 0
      if (iw-1.lt.iwini) stop 'iwini too high'
      if (iw+1.gt.iwend) stop 'iwend too low'
      call q2r(irws,rws,qbz,nrws,nqbz,nwf,zsc(:,:,:,iw-1),zcr)
      call r2k(irws,rws,kf,nrws,nwf,zcr,zckf(:,:,1)) ! SEc(kf,Ef-de)
      call q2r(irws,rws,qbz,nrws,nqbz,nwf,zsc(:,:,:,iw),zcr)
      call r2k(irws,rws,kf,nrws,nwf,zcr,zckf(:,:,2)) ! SEc(kf,Ef)
      call q2r(irws,rws,qbz,nrws,nqbz,nwf,zsc(:,:,:,iw+1),zcr)
      call r2k(irws,rws,kf,nrws,nwf,zcr,zckf(:,:,3)) ! SEc(kf,Ef+de)
      call diag_hm(zhkf,nwf,eval,evecc)
! Self-energy correction
      zw1 = zxkf+zckf(:,:,2)-zvkf
      call wan2ham(zw1,evecc,nwf,zw2) ! Hamiltonian gauge 
      ds = dreal(zw2(iwf,iwf))
! Z factor
      zw1 = zckf(:,:,1)
      call wan2ham(zw1,evecc,nwf,zw2) ! Hamiltonian gauge 
      dsm = dreal(zw2(iwf,iwf))
      zw1 = zckf(:,:,3)
      call wan2ham(zw1,evecc,nwf,zw2) ! Hamiltonian gauge 
      dsp = dreal(zw2(iwf,iwf))
      zfac = 1.0d0 / (1.0d0 - (dsp - dsm)/(2d0*dwev))
! Ef(GW)
      efgw = eflda + zfac * ds
!      efgw = eflda + ds
      ds0 = ds
!
      write(*,*)'Ef(LDA)  =',eflda,'eV'
      write(*,*)'Ef(GW)   =',efgw,'eV'
      write(*,*)'Z factor =',zfac
! 
      return
      end
c-----------------------------------------------------------------------
