      program hqpband_ent

c Nov,2008 Takashi Miyake, from hband_ent_MPI.m.F
c Nov,2008 Takashi Miyake, from hband_ent_single_MPI.m.F
c Nov,2008 Takashi Miyake, from hqpe_mlw.m.F
c calculates disentangled bands
c outputs
c d space
c   psi(k,d) in the Wannier gauge, expanded by KS
c   Ham(k;d,d') in the Wannier gauge
c   psi(k,d) in the Hamiltonian gauge, expanded by Wannier
c   E(k,d) in the KS gauge and Wannier gauge
c r space
c   psi(k,r) in the Wannier gauge, expanded by KS
c   Ham(k;r,r') in the Wannier gauge
c   psi(k,r) in the Hamiltonian gauge, expanded by Wannier
c   E(k,r) in the KS gauge and Wannier gauge
      use m_readeigen
      use m_read_bzdata,ngrp2=>ngrp
      use m_genallcf_v3
      use keyvalue
c RS: MPI module
      use rsmpi
      implicit real*8 (a-h,o-z)

c parameter
      integer(4):: nqbandx
      parameter (nqbandx = 1000)
      real(8):: qband(3,nqbandx)

c local data
      logical laf,lbnds

      integer(4):: ret
      integer(4):: iko_ix,iko_fx,nwf,mrece,ixc,idummy,nz
      integer(4),allocatable :: iqk4(:,:),itq(:),irws(:)
      real(8),allocatable :: elda(:,:),elda_bz(:,:),rws(:,:),drws(:),
     &    wqk4(:,:),qbz2(:,:),eval(:)
      real(8):: q(3),rydberg,pi,eps,aa(3,3),bb(3,3),dk(3),xk
      complex(8),allocatable:: umn_bz(:,:,:),umn(:,:,:),uuqk4(:,:,:,:),
     &     zw(:,:),zw0(:,:),zlda(:,:,:,:),zsc(:,:,:,:,:),
     &     zsx(:,:,:,:),zvxc(:,:,:,:),evecc(:,:),
     &     zlda2(:,:,:),zsc2(:,:,:,:),zsx2(:,:,:),zvxc2(:,:,:),
     &     zlda3(:,:),zsc3(:,:,:),zsx3(:,:),zvxc3(:,:),
     &     zs0(:,:),zsw(:,:,:)

c RS: program begins here
      call RSMPI_Init()
      eps = 1.d-4
      pi = 4d0*datan(1d0)

c mode switch
      call readin5(ixc,nz,idummy)
      write(6,*) ' ixc=',ixc
      if(ixc==0) stop ' --- ixc=0 --- Choose computational mode!'
      if (Is_IO_Root_RSMPI()) then
      if (ixc==1) then
         write(*,*)'read GW0'
      elseif (ixc==2) then
         write(*,*)'read GW0.wan'
      endif ! ixc   
      endif ! IO   
      call read_BZDATA()
      ifi = 99
      open(ifi,file='hbe.d',status='old')
      read (ifi,*) nprecb,mrecb,mrece,nlmtot,nqbzt, nband,mrecg
      close(ifi)
      call genallcf_v3(-999,-999d0,-1) !in module m_genallcf_v3
      call init_readeigen(ginv,nspin,nband,mrece) !initialization of readEigen 
      call init_readeigen2(mrecb,nlmto,mrecg) !initialize m_readeigen     

c Fermi level
      lbnds=.false.
      inquire(file='BNDS',exist=lbnds)      
      if (lbnds) then
        write(*,*)'Read EF from BNDS'
        open(99,file='BNDS',status='old')
        read(99,*)ntmp,ef
        close(99)
      else ! lbnds 
        stop' hqpband_ent: Cannot find BNDS'
      endif ! lbnds  

c read q points from SYML
      call read_syml(qband,nqbandx,nqband)
      if (nqband.gt.nqbandx) stop 'hqpband_ent: too many q points'

c WS supercell
      allocate(irws(n1*n2*n3*8),rws(3,n1*n2*n3*8),drws(n1*n2*n3*8))
      call wigner_seitz(alat,plat,n1,n2,n3,nrws,rws,irws,drws)

c read GW0
      open(ifi,file='GW0',form='formatted')
      read(ifi,*)nwf,nqbz2,nq1,nq2,nq3,iwini,iwend,dwplot,nspin2,
     &  ef2,alat2
      if (nqbz2.ne.nqbz) stop 'hqpband_ent: nqbz error in GW0'
      if (nq1.ne.n1) stop 'hqpband_ent: nq1 error in GW0'
      if (nq2.ne.n2) stop 'hqpband_ent: nq2 error in GW0'
      if (nq3.ne.n3) stop 'hqpband_ent: nq3 error in GW0'
      if (nspin2.ne.nspin) stop 'hqpband_ent: nspin error in GW0'
!      if (abs(ef-ef2).gt.eps) stop 'hqpband_ent: ef error in GW0'
      if (abs(alat-alat2).gt.eps) stop 'hqpband_ent: alat error in GW0'
      dwev = dwplot*rydberg()
      allocate(zw(nwf,nwf),zw0(nwf,nwf),qbz2(3,nqbz),
     &     zlda(nwf,nwf,nqbz,nspin),zsc(nwf,nwf,nqbz,iwini:iwend,nspin),
     &     zsx(nwf,nwf,nqbz,nspin),zvxc(nwf,nwf,nqbz,nspin),
     &     eval(nwf),evecc(nwf,nwf))
! aa: unit cell vectors in units of alat
      read(ifi,*)aa(1:3,1)
      read(ifi,*)aa(1:3,2)
      read(ifi,*)aa(1:3,3)
! bb: reciprocal vectors in units of 2pi/alat 
      read(ifi,*)bb(1:3,1)
      read(ifi,*)bb(1:3,2)
      read(ifi,*)bb(1:3,3)
      read(ifi,*)qbz2(:,:) ! in units of 2pi/alat
      do ik = 1,nqbz
      do ii = 1,3
         if (abs(qbz2(ii,ik)-qbz(ii,ik)).gt.eps)
     &   stop 'hqpband_ent: qbz error in GW0'
      enddo   
      enddo   
      do is = 1,nspin
      do ik = 1,nqbz
         read(ifi,*)zlda(:,:,ik,is) ! HLDA (eV)
!         call diag_hm(zlda(:,:,ik,is),nwf,eval,evecc)
         read(ifi,*)zvxc(:,:,ik,is) ! Vxc (eV)
c Wannier gauge to Hamiltonian gauge: Vxc
!         call wan2ham(zvxc(:,:,ik,is),evecc,nwf,zw)
         read(ifi,*)zsx(:,:,ik,is) ! SEex (eV)
c Wannier gauge to Hamiltonian gauge: Sx 
!         call wan2ham(zsx(:,:,ik,is),evecc,nwf,zw)
         do iw = iwini,iwend
           read(ifi,*)zsc(:,:,ik,iw,is) ! SEcor (eV)
         enddo ! iw 
c Wannier gauge to Hamiltonian gauge: Sc(w)
!         if (ik.eq.1) then 
!         do iw=iwini,iwend
!           call wan2ham(zsc(:,:,ik,iw,is),evecc,nwf,zw) 
!           write(11,*)iw,real(zw(22,22)),imag(zw(22,22))
!         enddo
!         endif
       enddo ! ik
       enddo ! is
       close(ifi)
! end of reading GW0

c loop over spin
      do      is = 1,nspin
      print *, ' --- is=',is

c> output file
      ifks  = 96
      ifd   = 95
      ifr   = 94
      ifd2  = 93
      if (is.eq.1) then
        open(ifks,file='qp.ks.up')
        open(ifd,file='qp.d.up')
        open(ifr,file='qp.r.up')
!        open(ifd2,file='qp.d2.up')
      else
        open(ifks,file='qp.ks.dn')
        open(ifd,file='qp.d.dn')
        open(ifr,file='qp.r.dn')
!        open(ifd2,file='qp.d2.dn')
      endif

c read Umn_bz
      if (is==1) then
        open(ifi,file='MLWU',form='unformatted',status='old')
      else
        open(ifi,file='MLWD',form='unformatted',status='old')
      endif ! is   
      read(ifi)nqbz2,nwf2,iko_ix,iko_fx
      if (nwf2 .ne. nwf) stop 'hqpband_ent: nwf error in MLWU/D'
      if (nqbz2 .ne. nqbz) stop 'hqpband_ent: nqbz error in MLWU/D'
      allocate(umn_bz(iko_ix:iko_fx,nwf,nqbz))
      do iqbz = 1,nqbz
          read(ifi)iqbz2,q
          if (iqbz2 .ne. iqbz) stop 'hqpband_ent: iqbz error'
          read(ifi)umn_bz(iko_ix:iko_fx,1:nwf,iqbz)
      enddo ! iqbz
      close(ifi)
!      call chk_umn(umn_bz,iko_ix,iko_fx,nwf,nqbz)

c read uumat(q,k)
c q from SYML, k on the BZ mesh
      allocate(uuqk4(nband,iko_ix:iko_fx,4,nqband),wqk4(4,nqband)
     &        ,iqk4(4,nqband))
      call readuu6(is,iko_ix,iko_fx,wqk4,iqk4,nqband,nqbz,nband,uuqk4)
!      call chk_uuqk(uuqk,iko_ix,iko_fx,nband,nqband)

c read Wannier orbital index to be included
      call getkeyvalue("GWinput","wan_ntq",ntq,default=-999)
      if (ntq.le.0.or.ntq.gt.nwf) stop 'hqpband_ent: ntq error'
      allocate(itq(ntq))
      call getkeyvalue("GWinput","wan_itq",itq,ntq)
      write(*,*)'ntq =',ntq
      write(*,*)'itq =',(itq(i),i=1,ntq)

      allocate(zlda2(ntq,ntq,nqbz),zsc2(ntq,ntq,nqbz,iwini:iwend),
     &     zsx2(ntq,ntq,nqbz),zvxc2(ntq,ntq,nqbz))
      do j = 1,ntq
      do i = 1,ntq
         zlda2(i,j,1:nqbz) = zlda(itq(i),itq(j),1:nqbz,is)
         zvxc2(i,j,1:nqbz) = zvxc(itq(i),itq(j),1:nqbz,is)
         zsx2(i,j,1:nqbz) = zsx(itq(i),itq(j),1:nqbz,is)
         zsc2(i,j,1:nqbz,iwini:iwend)
     & = zsc(itq(i),itq(j),1:nqbz,iwini:iwend,is)
      enddo   
      enddo   

c Umn_bz(qbz) to Umn(q)
c Umn(i,j,q) = <psi(q,i)|psi^(W)(q,j)> (i=1...nband,j=1...ntq)
c            = <psi(q,i)|psi^(w)(k,j)>
c            ~ <psi(q,i)|psi(k,l)><psi(k,l)|psi^(W)(k,j)>
c            = S[ii=1,4]wqk4(ii,q)*uuqk4(i,l,ii,q) * uum_bz(l,j,k)
      allocate(umn(nband,ntq,nqband) )
      umn = 0d0
      do iq  = 1,nqband
      do ii = 1,4   
         ik = iqk4(ii,iq)
         do ib = 1,nband
         do i = 1,ntq
         do it = iko_ix,iko_fx
           umn(ib,i,iq) = umn(ib,i,iq)
     &        + uuqk4(ib,it,ii,iq)*umn_bz(it,itq(i),ik)*wqk4(ii,iq)
         enddo ! it
         enddo ! i
         enddo ! ib
      enddo ! ii
      enddo ! iq
!      call chk_umn(umn,1,nband,ntq,nqband)

c read eks,nbks
      allocate(elda(nband,nqband))
      do iq = 1,nqband
         q(:) = qband(:,iq)
         call readeval(q,is,elda(:,iq))
      enddo ! iq
      elda = (elda-ef)*rydberg()
      allocate(elda_bz(nband,nqbz))
      do iq = 1,nqbz
         q(:) = qbz(:,iq)
         call readeval(q,is,elda_bz(:,iq))
      enddo ! iq
      elda_bz = (elda_bz-ef)*rydberg()

c disentangled band
      call disent_qpband(ixc,ifks,ifd,ifr,ifd2,
     &                 zlda2,zvxc2,zsx2,zsc2,
     &                 elda,elda_bz,umn_bz,umn,
     &                 dwev,iwini,iwend,
     &                 itq,iko_ix,iko_fx,iqk,qbz,qband,
     &                 rws,irws,nband,ntq,nwf,nqbz,nqband,nrws)

      deallocate(umn_bz,umn,uuqk4,wqk4,iqk4,elda,elda_bz,itq,
     &   zlda2,zvxc2,zsx2,zsc2)

      close(ifks)
      close(ifd)
      close(ifr)
!      close(ifd2)

      end do ! is

      call RSMPI_Finalize()
!      write(ifile_rsmpi,*) "RSMPI_Finalize is called."

      write(*,*)' OK! hqpband_ent '  

      end
c-----------------------------------------------------------------------
      subroutine readuu6(is,iti,itf,wqk4,iqk4,nqband,nqbz,nband,uuqk4)
      implicit real*8(a-h,o-z)

      integer(4) :: is,iti,itf,nqband,nband,iqk4(4,nqband)
      real(8) :: wqk4(4,nqband)
      complex(8) :: uuqk4(nband,iti:itf,4,nqband)

      nbb = 4
      ifi = 99
      if (is==1) then
        open(ifi,file='UUqk4U',form='unformatted',status='old')
      else 
        open(ifi,file='UUqk4D',form='unformatted',status='old')
      endif  

      read(ifi)
      read(ifi)nqband2,nqbz2,nbb2,iti2,itf2,nband2
      if (nqband2 .ne. nqband) stop 'readuu6: nqband error in UUqk4U/D'
      if (nqbz2 .ne. nqbz) stop 'readuu6: nqbz error in UUqk4U/D'
      if (nbb2 .ne. nbb) stop 'readuu6: nbb error in UUqk4U/D'
      if (nband2 .ne. nband) stop 'readuu6: nband error in UUqk4U/D'
      if (iti2 .ne. iti) stop 'readuu6: iti error'
      if (itf2 .ne. itf) stop 'readuu6: itf error'

      do iq = 1,nqband
      do ibb = 1,nbb
         read(ifi)iflg
         if (iflg.ne.-10) stop 'readuu6: iflg error'
            read(ifi) iq2,ibb2,iqk4(ibb,iq),wqk4(ibb,iq)
            read(ifi)((uuqk4(j1,j2,ibb,iq),j1=1,nband),j2=iti,itf)
            if (iq2 .ne. iq) stop 'readuu6: iq error'
            if (ibb2 .ne. ibb) stop 'readuu6: ibb error'
      enddo ! ibb
      enddo ! iq

      close(ifi)

      return
      end
c-----------------------------------------------------------------------
      subroutine chk_umn(umn,iti,itf,nwf,nqbz)
      implicit real*8(a-h,o-z)

      integer(4) :: iti,itf,nwf,nqbz
      complex(8) :: umn(iti:itf,nwf,nqbz)

      write(*,*)'*** check umn'
      do iq = 1,nqbz
         do iwf = 1,nwf
            a = sum(umn(:,iwf,iq)*conjg(umn(:,iwf,iq)))
            write(*,*)iwf,iq,a
         enddo ! iwf
      enddo ! iq   
      write(*,*)'*** end of chk_umn'

      return
      end
c-----------------------------------------------------------------------
      subroutine chk_uuqk(uuqk,iti,itf,nband,nq)
      implicit real*8(a-h,o-z)

      integer(4) :: iti,itf,nband,nq
      complex(8) :: uuqk(nband,iti:itf,nq)

      write(*,*)'*** check uuqk(1:nband,:)'
      do iq = 1,nq
         do it = iti,itf
            a = sum(uuqk(:,it,iq)*conjg(uuqk(:,it,iq)))
            write(*,*)it,iq,a
         enddo ! it
      enddo ! iq   
      write(*,*)'*** check uuqk(iti:itf,:)'
      do iq = 1,nq
         do it = iti,itf
            a = sum(uuqk(iti:itf,it,iq)*conjg(uuqk(iti:itf,it,iq)))
            write(*,*)it,iq,a
         enddo ! it
      enddo ! iq   
      write(*,*)'*** end of chk_uuqk'

      return
      end
c-----------------------------------------------------------------------
      subroutine disent_qpband(ixc,ifks,ifd,ifr,ifd2,
     &                 zlda,zvxc,zsx,zsc,
     &                 elda,elda_bz,umn_bz,umn,
     &                 dwev,iwini,iwend,
     &                 itq,iko_ix,iko_fx,iqk,qbz,qband,
     &                 rws,irws,nband,ntq,nwf,nqbz,nqband,nrws)
      use m_readeigen
      implicit real*8(a-h,o-z)

      integer(4) :: nband,ntq,nqbz,nqband,nrws,irws(nrws),itq(ntq)
     &             ,iqk(nqband),ixc
      real(8) :: elda(nband,nqband),elda_bz(nband,nqbz),
     &           qbz(3,nqbz),qband(3,nqband),
     &           rws(3,nrws),x(nqband),eval(ntq),dq(3),q(3),
     &           ed(ntq,nqband),ed2(ntq,nqband),er(nband,nqband),
     &   egw(nband,nqband,5),edgw(ntq,nqband,5),ergw(nband,nqband,5),
     &   dk(3)
      complex(8) :: ci,ceikr,
     &       umn_bz(iko_ix:iko_fx,nwf,nqbz),umn(nband,ntq,nqband),
     &       hk(ntq,ntq,nqbz),hr(ntq,ntq,nrws),hq(ntq,ntq),
     &       evecc(ntq,ntq),umnh(nband,ntq),
     &       zlda(ntq,ntq,nqbz),zvxc(ntq,ntq,nqbz),
     &       zsx(ntq,ntq,nqbz),zsc(ntq,ntq,nqbz,iwini:iwend),
     &       zse0(ntq,ntq,nqbz),
     &       se0r(ntq,ntq,nrws),sewr(ntq,ntq,nrws,iwini:iwend),
     &       se0q(ntq,ntq),sewq(ntq,ntq),zw(ntq,ntq,2)
      integer(4),allocatable :: iw1(:)
      real(8),allocatable :: dd(:),evald(:),evalr1(:),rw1(:),edr1r2(:)
      complex(8),allocatable :: aa(:,:),zz(:,:),vv(:,:),
     &  evecd(:,:),evecr1(:,:),zw1(:,:),zw2(:,:),zw3(:,:),
     &  zw4(:,:),zw5(:,:)
     & ,zlda3(:,:),zvxc3(:,:),zsx3(:,:),zs0(:,:),zsw(:,:,:)
     & ,zks(:,:,:),zham(:,:,:)

      pi = 4d0*datan(1.0d0)
      ci = (0.0d0,1.0d0)
      etol = 1.0d-3

! check the data
!       allocate(zlda3(ntq,ntq),zvxc3(ntq,ntq),zsx3(ntq,ntq),
!     &          zs0(ntq,ntq),zsw(ntq,ntq,2),zw1(ntq,ntq))
!         open(99,file='test2')
!         do iq = 1,nqband
!           if (iq .gt. 1) then
!             dk = qband(:,iq) - qband(:,iq-1)
!             xk = xk + dsqrt(sum(dk(:)**2))
!           else
!             xk = 0.0d0
!           endif
!           q(:) = qband(:,iq)
!           write(*,"('k=',3f12.6)")q
!           call q2p(irws,rws,qbz,q,nrws,nqbz,ntq,zlda,zlda3)
!           call q2p(irws,rws,qbz,q,nrws,nqbz,ntq,zvxc,zvxc3)
!           call q2p(irws,rws,qbz,q,nrws,nqbz,ntq,zsx,zsx3)
!           call diag_hm(zlda3,ntq,eval,evecc)
!           zw1 = zsx3 - zvxc3
!           call wan2ham(zw1,evecc,ntq,zs0)
!           write(99,*)
!           do kk = 1,ntq
!              iw = int(eval(kk)/dwev)
!              if (eval(kk).le.0d0) iw = iw - 1
!              call q2p(irws,rws,qbz,q,nrws,nqbz,ntq,
!     &                 zsc(:,:,iq,iw),zw1)
!              call wan2ham(zw1,evecc,ntq,zsw(:,:,1))
!              call q2p(irws,rws,qbz,q,nrws,nqbz,ntq,
!     &                 zsc(:,:,iq,iw+1),zw1)
!              call wan2ham(zw1,evecc,ntq,zsw(:,:,2))
!              zfac = (real(zsw(kk,kk,2))-real(zsw(kk,kk,1)))/dwev
!              zfac = 1.d0/(1.0d0-zfac)
!              rtmp = (eval(kk) - dwev*dble(iw))/dwev
!              rtmp2 = real(zsw(kk,kk,2))*rtmp
!     &              + real(zsw(kk,kk,1))*(1.0d0-rtmp)
!              etmp = eval(kk) + (real(zs0(kk,kk)))
!              write(99,"(f12.6,i5,3f12.6)")xk,kk,eval(kk),etmp,zfac
!           enddo ! kk
!         enddo ! iq
!         close(99)
!       deallocate(zlda3,zvxc3,zsx3,zs0,zsw,zw1)
!       stop 'debug: ok'

c Hamiltonian in the Wannier gauge
      hk = 0d0
      do iq = 1,nqbz
         do in = 1,ntq
         do im = 1,ntq
            ii = itq(im)
            ij = itq(in)
         do ib = iko_ix,iko_fx
            hk(im,in,iq) = hk(im,in,iq) + 
     &        dconjg(umn_bz(ib,ii,iq))*elda_bz(ib,iq)*umn_bz(ib,ij,iq)
         enddo ! ib   
         enddo ! im
         enddo ! in
      enddo ! iq   

c k to R
      call k2r(irws,rws,qbz,nrws,nqbz,ntq,hk,hr)
      zse0 = zsx-zvxc
      call k2r(irws,rws,qbz,nrws,nqbz,ntq,zse0,se0r)
      do iw = iwini,iwend
         call k2r(irws,rws,qbz,nrws,nqbz,ntq,zsc(:,:,:,iw),sewr(:,:,:,iw))
      enddo   

c loop over qband
      do iq = 1,nqband         
c x axis
        if (iq .gt. 1) then
          dq = qband(:,iq) - qband(:,iq-1)
          x(iq) = x(iq-1) + dsqrt(sum(dq(:)**2))
        else
          x(iq) = 0.0d0
        endif

c d band via Fourier transform
        call r2q(irws,rws,qband(:,iq),nrws,ntq,hr,hq)
        call diag_hm(hq,ntq,ed(:,iq),evecc)

c d band from another interpolation
c switch to the Hamiltonian gauge
        umnh = 0d0
!        iti = 1
!        itf = nband
!        iti = iko_ix
!        itf = iko_fx
        iti = max0(iko_ix-5,1)
        itf = min0(iko_fx+5,nband)
        ndr1 = itf - iti + 1
        do ib = iti,itf
        do iw = 1,ntq
        do ih = 1,ntq
          umnh(ib,ih) = umnh(ib,ih) + umn(ib,iw,iq)*evecc(iw,ih)
        enddo
        enddo
        enddo
c orthonormalize umnh
        allocate (aa(ndr1,ntq),zz(ndr1,ndr1),vv(ntq,ntq),dd(ntq))
        aa(1:ndr1,1:ntq) = umnh(iti:itf,1:ntq)
        call zgesvd_mn(ndr1,ntq,aa,dd,zz,vv)
        umnh = 0d0
        do ij = iti,itf
          jj = ij - iti + 1
          do ii = 1,ntq
          do kk = 1,ntq
            umnh(ij,ii) = umnh(ij,ii) + zz(jj,kk)*vv(kk,ii)
          enddo ! kk
          enddo ! ii
        enddo ! ij
        deallocate(aa,zz,vv,dd)
        call chk_orth(umnh(iti:itf,:),ndr1,ntq)
        call chk_norm(umnh(iti:itf,:),ndr1,ntq)

c d band from the second approach
        ed2(:,iq) = 0d0
        do ih = 1,ntq
        do ib = iti,itf
           ed2(ih,iq) = ed2(ih,iq)
     &      + dconjg(umnh(ib,ih))*elda(ib,iq)*umnh(ib,ih)
        enddo ! ib   
        enddo ! ih

c r1 space
        nd = ntq
        nr1 = ndr1 - nd
        allocate(evald(nd),evecd(iti:itf,nd),
     &     evalr1(nr1),evecr1(iti:itf,nr1),
     &     zw1(iti:itf,iti:itf),
     &     zw2(iti:itf,nr1),zw3(iti:itf,nr1),
     &     zw4(nr1,nr1),zw5(nr1,nr1),rw1(iti:itf),
     &     iw1(ndr1),edr1r2(nband))

! projection
        evecd(iti:itf,:) = umnh(iti:itf,:)
        zw1 = 0d0
        do ib = iti,itf
           zw1(ib,ib) = 1d0
        enddo ! ib   
        do ib = iti,itf
         do ib2 = iti,itf
          do iwf = 1,nd
           zw1(ib2,ib) = zw1(ib2,ib) - 
     &      evecd(ib2,iwf)*dconjg(evecd(ib,iwf))
          enddo ! iwf
         enddo ! ib2   
! norm
         rw1(ib) = sum(conjg(zw1(iti:itf,ib))
     &                      *zw1(iti:itf,ib) )
        enddo ! ib   
! sort the norm
        call sort_rarray(rw1(iti:itf),ndr1,iw1)
        zw2 = 0d0
        do ib = 1,nr1
          itmp = ndr1+1-ib
          itmp2 = iw1(itmp)+iti-1
          zw2(iti:itf,ib) = zw1(iti:itf,itmp2)
        enddo ! ib
! orthonormalize
        allocate (aa(ndr1,nr1),zz(ndr1,ndr1),vv(nr1,nr1),dd(nr1))
        aa(1:ndr1,1:nr1) = zw2(iti:itf,1:nr1)
        call zgesvd_mn(ndr1,nr1,aa,dd,zz,vv)
        zw3 = 0d0
        do ij = iti,itf
           jj = ij - iti + 1
           do ii = 1,nr1
           do kk = 1,nr1
            zw3(ij,ii) = zw3(ij,ii) + zz(jj,kk)*vv(kk,ii)
           enddo ! kk
           enddo ! ii
        enddo ! ij
        deallocate(aa,zz,vv,dd)
        call chk_orth(zw3(iti:itf,:),ndr1,nr1)
        call chk_norm(zw3(iti:itf,:),ndr1,nr1)
! Hamiltonian in the r1 space
        zw4 = 0d0
        do ii = 1,nr1
        do kk = 1,nr1
          do ij = iti,itf
           zw4(ii,kk) = zw4(ii,kk) + 
     &      conjg(zw3(ij,ii))*elda(ij,iq)*zw3(ij,kk)
          enddo ! ij 
        enddo ! kk
        enddo ! ii   
! diagonalize
        call diag_hermite(zw4,nr1,evalr1,zw5)
        evecr1 = 0d0
        do ij = iti,itf
        do ii = 1,nr1
        do kk = 1,nr1
          evecr1(ij,kk) = evecr1(ij,kk) + zw3(ij,ii)*zw5(ii,kk)
        enddo ! kk
        enddo ! ii 
        enddo ! ij 
! end of r1 space

        edr1r2 = 0.0d0
!        do ib = 1,nd
!          edr1r2(ib) = evald(ib)
!        enddo ! ib
        if (nr1.ge.1) then
         do ib = 1,nr1
!          edr1r2(ib+nd) = evalr1(ib) 
          edr1r2(ib) = evalr1(ib) 
         enddo ! ib
        endif ! nr1
        if (iti.ge.2) then
          do ib = 1,iti-1
!            edr1r2(ib+ndr1) = elda(ib,iq)
            edr1r2(ib+nr1) = elda(ib,iq)
          enddo ! ib
        endif ! iti
        if (itf.lt.nband) then
          do ib = itf+1,nband
!            edr1r2(ib) = elda(ib,iq)
            edr1r2(ib-nd) = elda(ib,iq)
          enddo ! ib
        endif ! itf

        er(:,iq) = edr1r2

c self-energy correction
!!!!! Kohn-Sham + GW
        allocate(zks(iti:itf,iti:itf,2))
        do ii = 1,5
          egw(:,iq,ii) = elda(:,iq)
        enddo  
        call r2q(irws,rws,qband(:,iq),nrws,ntq,se0r,se0q)
        do ib = iti,itf
           iw = int(elda(ib,iq)/dwev)
           if (elda(ib,iq).le.0.0d0) iw = iw-1
           if (iw.lt.iwini) iw = iwini
           if (iw.ge.iwend) iw = iwend-1
           rtmp = (elda(ib,iq) - dwev*dble(iw))/dwev
           call r2q(irws,rws,qband(:,iq),nrws,ntq,sewr(:,:,:,iw),
     &              zw(:,:,1))
           call r2q(irws,rws,qband(:,iq),nrws,ntq,sewr(:,:,:,iw+1),
     &              zw(:,:,2))
! full matrix in the Wannier gauge
           sewq = zw(:,:,1) + se0q
           call wan2ks(sewq,umn(iti:itf,1:ntq,iq),itf-iti+1,ntq,
     &                 zks(iti:itf,iti:itf,1))
           sewq = zw(:,:,2) + se0q
           call wan2ks(sewq,umn(iti:itf,1:ntq,iq),itf-iti+1,ntq,
     &                 zks(iti:itf,iti:itf,2))
           zfac = (real(zks(ib,ib,2))-real(zks(ib,ib,1)))/dwev
           zfac = 1.d0/(1.0d0-zfac)
           rtmp = (elda(ib,iq) - dwev*dble(iw))/dwev
           rtmp2 = real(zks(ib,ib,2))*rtmp
     &           + real(zks(ib,ib,1))*(1.0d0-rtmp)
           egw(ib,iq,2) = egw(ib,iq,2)+zfac*rtmp2
           do ibp = iti,itf
              ediff = elda(ib,iq) - elda(ibp,iq)
              if (ibp.ne.ib.and.abs(ediff).gt.etol) then
               rtmp2 = abs(zks(ibp,ib,2)*rtmp
     &           +         zks(ibp,ib,1)*(1.0d0-rtmp))
               egw(ib,iq,3) = egw(ib,iq,2)+rtmp2**2/ediff
              endif 
           enddo ! ibp   
! diagonal terms in the Wannier gauge
           sewq = zw(:,:,1) + se0q
           call wan2ks_d(sewq,umn(iti:itf,1:ntq,iq),itf-iti+1,ntq,
     &                 zks(iti:itf,iti:itf,1))
           sewq = zw(:,:,2) + se0q
           call wan2ks_d(sewq,umn(iti:itf,1:ntq,iq),itf-iti+1,ntq,
     &                 zks(iti:itf,iti:itf,2))
           zfac = (real(zks(ib,ib,2))-real(zks(ib,ib,1)))/dwev
           zfac = 1.d0/(1.0d0-zfac)
           rtmp = (elda(ib,iq) - dwev*dble(iw))/dwev
           rtmp2 = real(zks(ib,ib,2))*rtmp
     &           + real(zks(ib,ib,1))*(1.0d0-rtmp)
           egw(ib,iq,4) = egw(ib,iq,4)+zfac*rtmp2
           do ibp = iti,itf
              ediff = elda(ib,iq) - elda(ibp,iq)
              if (ibp.ne.ib.and.abs(ediff).gt.etol) then
               rtmp2 = abs(zks(ibp,ib,2)*rtmp
     &           +         zks(ibp,ib,1)*(1.0d0-rtmp))
               egw(ib,iq,5) = egw(ib,iq,4)+rtmp2**2/ediff
              endif 
           enddo ! ibp   
        enddo ! ib   
        deallocate(zks)

!!!!! d-band + GW
c d band via Fourier transform
        allocate(zham(ntq,ntq,2))
!        call r2q(irws,rws,qband(:,iq),nrws,ntq,hr,hq)
!        call diag_hm(hq,ntq,ed(:,iq),evecc)
        do ii = 1,5
          edgw(:,iq,ii) = ed(:,iq)
        enddo  
        do ib = 1,ntq
           iw = int(ed(ib,iq)/dwev)
           if (ed(ib,iq).le.0.0d0) iw = iw-1
           if (iw.lt.iwini) iw = iwini
           if (iw.ge.iwend) iw = iwend-1
           rtmp = (ed(ib,iq) - dwev*dble(iw))/dwev
           call r2q(irws,rws,qband(:,iq),nrws,ntq,sewr(:,:,:,iw),
     &              zw(:,:,1))
           call r2q(irws,rws,qband(:,iq),nrws,ntq,sewr(:,:,:,iw+1),
     &              zw(:,:,2))
! full matrix in the Wannier gauge
           sewq = zw(:,:,1) + se0q
           call wan2ham(sewq,evecc,ntq,zham(:,:,1))
           sewq = zw(:,:,2) + se0q
           call wan2ham(sewq,evecc,ntq,zham(:,:,2))
           zfac = (real(zham(ib,ib,2))-real(zham(ib,ib,1)))/dwev
           zfac = 1.d0/(1.0d0-zfac)
           rtmp = (ed(ib,iq) - dwev*dble(iw))/dwev
           rtmp2 = real(zham(ib,ib,2))*rtmp
     &           + real(zham(ib,ib,1))*(1.0d0-rtmp)
           edgw(ib,iq,2) = edgw(ib,iq,2)+zfac*rtmp2
           do ibp = 1,ntq
              ediff = ed(ib,iq) - ed(ibp,iq)
              if (ibp.ne.ib.and.abs(ediff).gt.etol) then
               rtmp2 = abs(zham(ibp,ib,2)*rtmp
     &           +         zham(ibp,ib,1)*(1.0d0-rtmp))
               edgw(ib,iq,3) = edgw(ib,iq,2)+rtmp2**2/ediff
              endif 
           enddo ! ibp   
! diagonal terms in the Wannier gauge
           sewq = zw(:,:,1) + se0q
           call wan2ham_d(sewq,evecc,ntq,zham(:,:,1))
           sewq = zw(:,:,2) + se0q
           call wan2ham_d(sewq,evecc,ntq,zham(:,:,2))
           zfac = (real(zham(ib,ib,2))-real(zham(ib,ib,1)))/dwev
           zfac = 1.d0/(1.0d0-zfac)
           rtmp = (ed(ib,iq) - dwev*dble(iw))/dwev
           rtmp2 = real(zham(ib,ib,2))*rtmp
     &           + real(zham(ib,ib,1))*(1.0d0-rtmp)
           edgw(ib,iq,4) = edgw(ib,iq,4)+zfac*rtmp2
           do ibp = 1,ntq
              ediff = ed(ib,iq) - ed(ibp,iq)
              if (ibp.ne.ib.and.abs(ediff).gt.etol) then
               rtmp2 = abs(zham(ibp,ib,2)*rtmp
     &           +         zham(ibp,ib,1)*(1.0d0-rtmp))
               edgw(ib,iq,5) = edgw(ib,iq,4)+rtmp2**2/ediff
              endif 
           enddo ! ibp   
        enddo ! ib   
        deallocate(zham)
!
        deallocate(evald,evecd,evalr1,evecr1,
     &     zw1,zw2,zw3,zw4,zw5,rw1,iw1,edr1r2)


      enddo ! iq

c output
      nr = nband - nd
      call wqpband(ifks,x,egw,nband,nqband)
      call wqpband(ifd,x,edgw,ntq,nqband)
!      call wband(ifd2,x,ed2,ntq,nqband)
      call wband(ifr,x,er(1:nr,1:nqband),nr,nqband)

      return
      end
c-----------------------------------------------------------------------
      subroutine wband(ifo,x,y,nb,nq)
      implicit real*8(a-h,o-z)

      real(8) :: x(nq),y(nb,nq)

      do ib  = 1,nb
      do iq = 1,nq
         write (ifo,"(1x,2f11.5)") x(iq),y(ib,iq)
      enddo
      write(ifo,*)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine wqpband(ifo,x,y,nb,nq)
      implicit real*8(a-h,o-z)

      real(8) :: x(nq),y(nb,nq,5)

      do ib  = 1,nb
      do iq = 1,nq
         write (ifo,"(1x,6f11.5)") x(iq),(y(ib,iq,ii),ii=1,5)
      enddo
      write(ifo,*)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine k2r(irws,rws,qbz,nrws,nqbz,nwf,zk,zr)
      implicit none

      integer(4):: irws(nrws),nrws,nqbz,nwf,ik,ir
      real(8):: rws(3,nrws),qbz(3,nqbz),pi,rk
      complex(8):: zr(nwf,nwf,nrws),zk(nwf,nwf,nqbz),ci,ceikr

      ci = (0d0,1d0)
      pi = 4d0*datan(1.0d0)

      zr = 0d0
      do ir = 1,nrws
      do ik = 1,nqbz
         rk = 2d0*pi*sum(rws(:,ir)*qbz(:,ik))
         ceikr = exp(-ci*rk)
         zr(:,:,ir) = zr(:,:,ir) + ceikr * zk(:,:,ik) / dble(nqbz)
      enddo ! ik
      enddo ! ir

      return
      end
c-----------------------------------------------------------------------
      subroutine r2q(irws,rws,p,nrws,nwf,zr,zq)
      implicit none

      integer(4):: irws(nrws),nrws,nwf,iq,ik,ir,ii,jj
      real(8):: rws(3,nrws),p(3),pi,rk
      complex(8):: zr(nwf,nwf,nrws),zq(nwf,nwf),ci,ceikr

      ci = (0d0,1d0)
      pi = 4d0*datan(1.0d0)

      zq = 0d0
      do ir = 1,nrws
         rk = 2d0*pi*sum(rws(:,ir)*p(:))
         ceikr = exp(ci*rk) / dble(irws(ir))
         zq(:,:) = zq(:,:) + ceikr * zr(:,:,ir)
      enddo ! ir 

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
      subroutine wan2ham_d(zwan,evecc,nwan,zham)
      implicit real*8(a-h,o-z)

      integer ii,jj
      complex(8) :: zwan(nwan,nwan),evecc(nwan,nwan),zham(nwan,nwan)

      zwan = 0.0d0
      do jj=1,nwan
      do ii=1,nwan
        do kk = 1,nwan
         zwan(ii,jj)=zwan(ii,jj)
     &   + dconjg(evecc(kk,ii))*zwan(kk,kk)*evecc(kk,jj)
        enddo 
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine wan2ks(zwan,evecc,nks,nwan,zks)
      implicit real*8(a-h,o-z)

      integer ii,jj
      complex(8) :: zwan(nwan,nwan),evecc(nks,nwan),zks(nks,nks),
     &              zw(nks,nwan)

      zw = matmul(evecc,zwan)
      do jj=1,nks
      do ii=1,nks
         zks(ii,jj)=sum(zw(ii,:)*dconjg(evecc(jj,:)))
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine wan2ks_d(zwan,evecc,nks,nwan,zks)
      implicit real*8(a-h,o-z)

      integer ii,jj
      complex(8) :: zwan(nwan,nwan),evecc(nks,nwan),zks(nks,nks)

      zks = 0.0d0
      do jj=1,nks
      do ii=1,nks
        do kk = 1,nwan
         zks(ii,jj)=zks(ii,jj)
     &   + evecc(ii,kk)*zwan(kk,kk)*dconjg(evecc(jj,kk))
        enddo 
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
