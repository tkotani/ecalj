      program hband_ent

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
      integer(4):: iko_ix,iko_fx,nwf,mrece
      integer(4),allocatable :: iqk(:),itq(:),irws(:)
      real(8),allocatable :: elda(:,:),elda_bz(:,:),rws(:,:),drws(:)
      real(8):: q(3),rydberg
      complex(8),allocatable:: umn_bz(:,:,:),umn(:,:,:),uuqk(:,:,:)

c RS: program begins here
      call RSMPI_Init()

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
        stop' hband_ent: Cannot find BNDS'
      endif ! lbnds  

c open files
c> input files

c read q points from SYML
      call read_syml(qband,nqbandx,nqband)
      if (nqband.gt.nqbandx) stop 'hband_ent: too many q points'

c WS supercell
      allocate(irws(n1*n2*n3*8),rws(3,n1*n2*n3*8),drws(n1*n2*n3*8))
      call wigner_seitz(alat,plat,n1,n2,n3,nrws,rws,irws,drws)

c loop over spin
      do      is = 1,nspin
      print *, ' --- is=',is

c> output file
!      ifbnd = 98
!      ifpsi = 97
      ifks  = 96
      ifd   = 95
      ifr   = 94
      ifd2  = 93
      if (is.eq.1) then
        open(ifks,file='bnds.ks.up')
        open(ifd,file='bnds.d.up')
        open(ifr,file='bnds.r.up')
        open(ifd2,file='bnds.d2.up')
!        open(ifpsi,file='wf.disent.up')
      else
        open(ifks,file='bnds.ks.dn')
        open(ifd,file='bnds.d.dn')
        open(ifr,file='bnds.r.dn')
        open(ifd2,file='bnds.d2.dn')
!        open(ifpsi,file='wf.disent.dn')
      endif

c read Umn_bz
      if (is==1) then
        open(ifi,file='MLWU',form='unformatted',status='old')
      else
        open(ifi,file='MLWD',form='unformatted',status='old')
      endif ! is   
      read(ifi)nqbz2,nwf,iko_ix,iko_fx
      if (nqbz2 .ne. nqbz) stop 'hband_ent: nqbz error in MLWU/D'
      allocate(umn_bz(iko_ix:iko_fx,nwf,nqbz))
      do iqbz = 1,nqbz
          read(ifi)iqbz2,q
          if (iqbz2 .ne. iqbz) stop 'hband_ent: iqbz error'
          read(ifi)umn_bz(iko_ix:iko_fx,1:nwf,iqbz)
      enddo ! iqbz
      close(ifi)
!      call chk_umn(umn_bz,iko_ix,iko_fx,nwf,nqbz)

c read uumat(q,k)
c q from SYML, k on the BZ mesh
      allocate(uuqk(nband,iko_ix:iko_fx,nqband),iqk(nqband))
      call readuu5(is,iko_ix,iko_fx,iqk,nqband,nqbz,nband,uuqk)
check q2
!      do iq = 1,nqband
!        q(:) = qband(:,iq)
!        iq2 = iqk(iq)
!        call findk(q,qbz,plat,qbas,nqbz,iq2)
!        if (iq2.ne.iqk(iq)) then 
!          write(*,*)'*** iqk error',iq,iq2,iqk(iq)
!          write(*,"(3f12.6)")q
!          write(*,"(3f12.6)")qbz(:,iq2)
!          write(*,"(3f12.6)")qbz(:,iqk(iq))
!        endif  
!        write(60,*)iq,iq2,q,qbz(:,iq2)
!      enddo  
!      stop
c
!      call chk_uuqk(uuqk,iko_ix,iko_fx,nband,nqband)

c read Wannier orbital index to be included
      call getkeyvalue("GWinput","wan_ntq",ntq,default=-999)
      if (ntq.le.0.or.ntq.gt.nwf) stop 'hband_ent: ntq error'
      allocate(itq(ntq))
      call getkeyvalue("GWinput","wan_itq",itq,ntq)
      write(*,*)'ntq =',ntq
      write(*,*)'itq =',(itq(i),i=1,ntq)

c Umn_bz(qbz) to Umn(q)
c Umn(i,j,q) = <psi(q,i)|psi^(W)(q,j)> (i=1...nband,j=1...ntq)
c            = <psi(q,i)|psi^(w)(k,j)>
c            ~ <psi(q,i)|psi(k,l)><psi(k,l)|psi^(W)(k,j)>
c            = uuqk(i,l,q) * uum_bz(l,j,k)
      allocate(umn(nband,ntq,nqband) )
      umn = 0d0
      do iq  = 1,nqband
         ik = iqk(iq)
         do ib = 1,nband
         do i = 1,ntq
         do it = iko_ix,iko_fx
           umn(ib,i,iq) = umn(ib,i,iq)
     &        + uuqk(ib,it,iq)*umn_bz(it,itq(i),ik)
         enddo ! it
         enddo ! i
         enddo ! ib
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

c Kohn-Sham band (for debug)
!      ifi = 99
!      open(ifi,file='BNDS2')
!         write(ifi,*)nband,ef,nqband
!      do iq = 1,nqband
!         write(ifi,"('*** q =',3f12.6)")qband(:,iq)
!         write(ifi,"(10f8.4)")elda(1:10,iq)
!      enddo ! iq
!      close(ifi)

c disentangled band
      call disent_band(ifks,ifd,ifr,ifd2,
     &                 elda,elda_bz,umn_bz,umn,
     &                 itq,iko_ix,iko_fx,iqk,qbz,qband,
     &                 rws,irws,nband,ntq,nwf,nqbz,nqband,nrws)

      deallocate(umn_bz,umn,uuqk,iqk,elda,elda_bz,itq)

      close(ifks)
      close(ifd)
      close(ifr)
      close(ifd2)
!      close(ifbnd)
!      close(ifpsi)

      end do ! is

      call RSMPI_Finalize()
!      write(ifile_rsmpi,*) "RSMPI_Finalize is called."

      write(*,*)' OK! hband_ent '  

      end
c-----------------------------------------------------------------------
      subroutine readuu5(is,iti,itf,iqk,nqband,nqbz,nband,uuqk)
      implicit real*8(a-h,o-z)

      integer(4) :: is,iti,itf,nqband,nband,iqk(nqband)
      complex(8) :: uuqk(nband,iti:itf,nqband)

      nbb = 1
      ifi = 99
      if (is==1) then
        open(ifi,file='UUqkU',form='unformatted',status='old')
      else 
        open(ifi,file='UUqkD',form='unformatted',status='old')
      endif  

      read(ifi)
      read(ifi)nqband2,nqbz2,nbb2,iti2,itf2,nband2
      if (nqband2 .ne. nqband) stop 'readuu5: nqband error in UUqkU/D'
      if (nqbz2 .ne. nqbz) stop 'readuu5: nqbz error in UUqkU/D'
      if (nbb2 .ne. nbb) stop 'readuu5: nbb error in UUqkU/D'
      if (nband2 .ne. nband) stop 'readuu5: nband error in UUqkU/D'
      if (iti2 .ne. iti) stop 'readuu5: iti error'
      if (itf2 .ne. itf) stop 'readuu5: itf error'

      do iq = 1,nqband
      do ibb = 1,nbb
         read(ifi)iflg
         if (iflg.ne.-10) stop 'readuu5: iflg error'
            read(ifi) iq2,iqk(iq)
            read(ifi)((uuqk(j1,j2,iq),j1=1,nband),j2=iti,itf)
            if (iq2 .ne. iq) stop 'readuu5: iq error'
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
      subroutine disent_band(ifks,ifd,ifr,ifd2,
     &                 elda,elda_bz,umn_bz,umn,
     &                 itq,iko_ix,iko_fx,iqk,qbz,qband,
     &                 rws,irws,nband,ntq,nwf,nqbz,nqband,nrws)
      use m_readeigen
      implicit real*8(a-h,o-z)

      integer(4) :: nband,ntq,nqbz,nqband,nrws,irws(nrws),itq(ntq)
     &             ,iqk(nqband)
      real(8) :: elda(nband,nqband),elda_bz(nband,nqbz),
     &           qbz(3,nqbz),qband(3,nqband),
     &           rws(3,nrws),x(nqband),eval(ntq),dq(3),
     &           ed(ntq,nqband),ed2(ntq,nqband),er(nband,nqband)
      complex(8) :: ci,ceikr,
     &       umn_bz(iko_ix:iko_fx,nwf,nqbz),umn(nband,ntq,nqband),
     &       hk(ntq,ntq,nqbz),hr(ntq,ntq,nrws),hq(ntq,ntq),
     &       evecc(ntq,ntq),umnh(nband,ntq)
      integer(4),allocatable :: iw1(:)
      real(8),allocatable :: dd(:),evald(:),evalr1(:),rw1(:),edr1r2(:)
      complex(8),allocatable :: aa(:,:),zz(:,:),vv(:,:),
     &  evecd(:,:),evecr1(:,:),zw1(:,:),zw2(:,:),zw3(:,:),
     &  zw4(:,:),zw5(:,:)

      pi = 4d0*datan(1.0d0)
      ci = (0.0d0,1.0d0)

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
      hr = 0d0
      do ir = 1,nrws
      do iq = 1,nqbz
         rk = 2d0*pi*sum(rws(:,ir)*qbz(:,iq))
         ceikr = exp(-ci*rk)
         do im = 1,ntq
         do in = 1,ntq
            hr(im,in,ir) = hr(im,in,ir) + 
     &      ceikr * hk(im,in,iq) / dble(nqbz)
         enddo ! in
         enddo ! im
      enddo ! iq
      enddo ! ir

c loop over qband
      do iq = 1,nqband
         
c x axis
        if (iq .gt. 1) then
          dq = qband(:,iq) - qband(:,iq-1)
          x(iq) = x(iq-1) + dsqrt(sum(dq(:)**2))
        else
          x(iq) = 0.0d0
        endif

c debug:
c        write(60,"(f8.4,2i5,6f8.4)")
c     &   x(iq),iq,iqk(iq),qband(:,iq),qbz(:,iqk(iq))

c d band via Fourier transform
        hq = 0d0
        do ir = 1,nrws
         rk = 2d0*pi*sum(rws(:,ir)*qband(:,iq))
         ceikr = exp(ci*rk) / dble(irws(ir))
         do im = 1,ntq
         do in = 1,ntq
            hq(im,in) = hq(im,in) + ceikr * hr(im,in,ir)
         enddo
         enddo
        enddo
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
c debug:
!          write(*,*)'***',iqbz,ikp
!          do jj = 1,nr1
!          do ii = 1,nd
!            ctmp = sum(conjg(zw3(:,jj))*evecd(:,ii))
!            rnorm = real(ctmp)
!            cnorm = imag(ctmp)
!            if (abs(ctmp).gt.1.d-3)write(*,"(2i5,2f12.6)")ii,jj,ctmp
!          enddo ! ii
!          enddo ! jj
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
            edr1r2(ib-nd) = elda(ib,iq)
          enddo ! ib
        endif ! itf

        er(:,iq) = edr1r2

c debug:
c        write(80,*)'*** iq=',iq,nr1,iti,itf
c        write(80,*)'ndim',nband,nd,nr1,ndr1
c        write(80,"(6f8.4)")evalr1(1:6)
c        write(80,"(6f8.4)")er(1:6,iq)
c        write(80,"(10f8.4)")er(7:16,iq)
c        write(80,"(10f8.4)")er(17:26,iq)
c        write(80,"(8f8.4)")er(27:34,iq)

        deallocate(evald,evecd,evalr1,evecr1,
     &     zw1,zw2,zw3,zw4,zw5,rw1,iw1,edr1r2)

      enddo ! iq

c output
      nr = nband - nd
      call wband(ifks,x,elda,nband,nqband)
      call wband(ifd,x,ed,ntq,nqband)
      call wband(ifd2,x,ed2,ntq,nqband)
      call wband(ifr,x,er(1:nr,1:nqband),nr,nqband)

      return
      end
c-----------------------------------------------------------------------
      subroutine wband(ifo,x,y,nb,nq)
      implicit real*8(a-h,o-z)

      real(8) :: x(nq),y(nb,nq)

      ndiv = nb/5
      if (ndiv*5.ne.nb) ndiv = ndiv + 1
      do i  = 1,ndiv
         i1 = 1 + (i-1)*5
         i2 = min0 (i1+4,nb)
         do iq = 1,nq
            write (ifo,"(1x,6f11.5)") x(iq),(y(ib,iq),ib=i1,i2)
         enddo
         write(ifo,*)
      enddo

      return
      end
c-----------------------------------------------------------------------
