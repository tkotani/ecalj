!>  make sigm.* file which contains Sigma-Vxc 
!     This uses amix in order to guess a better Sigma-Vxc from previous iterations.
!     * iSigma_en==5 is for diagonal-only Sigma-Vxc in LDA basis set.
!      (Then you need evec0, which contains eigenvector of LDA).
      program hqpe_sc
c------------------------------------------------------------------
c     calculates quasiparticle energies
c     E(k,t) = e(k,t) + Z [SEx(k,t) + SEc(k,t) - xcLDA(k,t)]
c     e(k,t) = LDA eigenvalue
c     Z      = [1 - dSEc(e(k,t))/dw]^(-1)
c     SEx(k,t)   = <psi(k,t)| SEx |psi(k,t)>
c     SEc(k,t)   = <psi(k,t)| SEc |psi(k,t)>, SEc = GWc
c     xcLDA(k,t) = <psi(k,t)| vxc |psi(k,t)>
c     SEx and xcLDA are in file SEX
c     SEc is in file SEC
c-----------------------------------------------------------------------
      use m_keyvalue,only: Getkeyvalue
      use m_read_bzdata, only: Read_bzdata, nstar, nqibz2=>nqibz
!! laf is used (antiferro or not) after anfcond.
      use m_anf,only: anfcond,laf
      use m_hamindex,only: Readhamindex, nhq=>ndham
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
c     local data
      logical sigma_mixing,lsi, lsigin
      dimension ifsex(2),ifsexcore(2),ifxc(2),ifsec(2),ifqpe(2)
     &     ,iftote(2),iftote2(2),ifsex2(2),ifsexcore2(2),ifsec2(2) !sf..3June
      integer(4),allocatable :: itxc(:),itc(:),itx(:)
      real(8),allocatable :: qxc(:,:,:),eldaxc(:,:),vxc(:,:),
     &     qc(:,:,:),eldac(:,:),sex(:,:),sexcore(:,:),
     &     qx(:,:,:),eldax(:,:),rsec(:,:),csec(:,:) !sf..13May
     &     ,qqq(:,:,:)          !sf..3June
     &     ,qqqx_m(:,:,:)
      complex(8), allocatable :: sex2(:,:,:),sexcore2(:,:,:),
     &     sec2(:,:,:), se(:,:,:),work(:),evec_inv(:,:),evec_invt(:,:), !sf.beg
     &     se_ev(:,:), ev_se_ev(:,:), sen(:,:),sen2(:,:),se_in(:,:),           
     &     v_xc(:,:,:,:),evec(:,:,:,:),evec0(:,:),sigma_m(:,:,:,:),
     &     sigin(:,:,:,:),
     &     evec00(:,:,:,:) ,evec00inv(:,:)
      real(8) qx2(3) ,qqqx(3) ,del ,mix_fac
      integer, allocatable :: ipiv(:) !sf.end
      real(8) ::  rydberg,hartree, qqqx0(3)
      integer(4) ::  n1,n2,n3, ifse_out, ifse_in
      logical :: AddDummySig,evec0ex=.false.

      complex(8),allocatable:: sigma_m_out(:,:,:,:),pmat(:,:)
     &     ,pmatd(:,:),sigmv(:,:,:)

c----------------------------------------------
      character(2):: soflag
      real(8) :: wex !,  eseavrs(2)
      logical :: exonly,mtosigmaonly,nexist
      integer(4):: ret,ifile_handle
      character(3):: iaaa
      real(8),allocatable:: eseavr(:,:) !,eseavr_in(:,:)
      integer,allocatable:: nhqx(:,:)
      integer:: nz,ntqmin, nmto,ndimsig,ndimsig2
      integer:: procid,nrank,ifigwb_,ifigwx1_,ifigwx2_,ifvxc_,ifevec_, ipx
      character*256:: extn,ext
      character*256,allocatable:: extp(:)
      integer,allocatable:: ifevec__(:),ifvxc__(:),iprocq(:,:)

      integer,allocatable:: ntqxx(:)
      real(8):: eseavrmean,eseadd,tolq=1d-8
      integer,allocatable:: nev(:,:)
      complex(8),allocatable:: ovl(:,:)
      character(8):: xt
      hartree=2d0*rydberg()

!! Shift quasiparticle energies (eV)
      write (*,*)' q+band index for zero?'
      read (*,*)jin
c      call headver('hqpe_sc',jin)
      call getkeyvalue("GWinput","EXonly",wex,default=0d0,status=ret)
      if(wex==0d0) then
        exonly=.false.
      else
        exonly=.true.
        write(6,*)' exonly=T wex=',wex
      endif

!! antiferro or not. Only calculate up spin
!! For AF case, we have laf=.true. and we have data set for 'call anfsig', stored in m_anf.
      call anfcond() 
      call readhamindex()
c$$$      inquire(file='ANFcond',EXIST=laf)
c$$$      if(laf) then
c$$$        ifianf = 211
c$$$        open(ifianf,file='ANFcond')
c$$$        read(ifianf,*) soflag
c$$$        if(soflag=='SO') then
c$$$          write(6,*)' ANFcond found. but SO case is not implimented. Not using ANF symm'
c$$$          laf=.false.
c$$$        else
c$$$          write(6,*) 'Find ANFcond--- up only'
c$$$        endif
c$$$        close(ifianf)
c$$$      endif

      ifsex(1) = iopenxx ('SEXU')
      ifxc(1)  = iopenxx ('XCU')
      ifsex2u = ifile_handle()
      open(ifsex2u, file='SEX2U',form='UNFORMATTED', status='OLD') !sf.beg
      ifsex2(1)=ifsex2u

      if(.not.exonly) then
        ifsec(1) = iopenxx ('SECU')
        ifsec2u = ifile_handle() 
        open(ifsec2u, file='SEC2U',form='UNFORMATTED', status='OLD')
        ifsec2(1)=ifsec2u     !sf.end
      endif
      ifsexcore(1) = iopenxx ('SEXcoreU')
      ifSEXcore2U=ifile_handle()
      open(ifSEXcore2U,file='SEXcore2U',form='UNFORMATTED',status='OLD')
      ifsexcore2(1)=ifSEXcore2U

      call readx (ifsex(1),50)
      read (ifsex(1),*) nspin,nq,ntq
      if(nspin .eq. 2 .and. .not. laf) then
        ifsex(2)   = iopenxx ('SEXD')
        ifxc(2)    = iopenxx ('XCD')
        ifSEX2D=ifile_handle()
        open(ifSEX2D, file='SEX2D',form='UNFORMATTED', status='OLD') !sf.beg
        ifsex2(2)=ifSEX2D
        if(.not.exonly) then
          ifsec(2)   = iopenxx ('SECD')
          ifSEC2D=ifile_handle()
          open(ifSEC2D, file='SEC2D',form='UNFORMATTED', status='OLD')
          ifsec2(2)=ifSEC2D      !sf.end
        endif
        ifsexcore(2)   = iopenxx ('SEXcoreD')
        ifSEXcore2D=ifile_handle()
        open(ifSEXcore2D,file='SEXcore2D',form='UNFORMATTED',status='OLD')
        ifsexcore2(2)=ifSEXcore2D
      endif

      rewind (ifsex(1))

c     june209
      write(6,*)'nq nspin=',nq,nspin
      allocate(eseavr(nq,nspin)) !,eseavr_in(nq,nspin))


c     > output file
      ifqpe(1)   = iopenxx ('QPU')
      iftote(1)  = iopenxx ('TOTE.UP')
      iftote2(1) = iopenxx ('TOTE2.UP')



C     3 Feb 2004 read sigma and pass to mixing routine
      INQUIRE (FILE='sigm', EXIST = lsigin)
      ifsigm=ifile_handle()
      open(ifsigm, file='sigm',form='UNFORMATTED') !sf
      ifse_out=ifsigm
      if (lsigin) then
        rewind ifse_out
        write(6,*) ' ... reading input sigma from file sigm'
c        read(ifse_out) nspin,ndimsigin,n1,n2,n3,nqx
        read(ifse_out) nspin,ndimsigin,n1,n2,n3,nqx
        if (nqx .ne. nq) then
          print 368, nqx,nq
 368      format (6x,' (warning) file mismatch : file nq=',i4,
     .           ' but expected',i4)
          lsigin = .false.
        else
          rewind ifse_out
          allocate(sigin(ndimsigin,ndimsigin,nq,nspin),qqqx_m(3,nq,nspin))
          call rwsigma(ifse_out,sigin,qqqx_m,nspin,ndimsigin,n1,n2,n3,nq) !,eseavr_in)
          deallocate(qqqx_m)
        endif
      endif
      if(laf) nspin=2

      rewind ifse_out
      !sf
      if (nspin == 2) then
        ifqpe(2)   = iopenxx ('QPD')
        iftote(2)  = iopenxx ('TOTE.DN')
        iftote2(2) = iopenxx ('TOTE2.DN')
      endif

      if(jin == -101) goto 9998

c-----------reading self-energy mode parameter from file 'GWIN_V2'
c      call readd_iSigma_en(ifinin,iSigma_en) !sf..21May02
      call getkeyvalue("GWinput","iSigMode",iSigma_en )
      if(isigma_en==5) then     !.or.core3ptest) then
c     inquire(file='evec0',exist=evec0ex)
cccccccccccccccccccccccccc
        evec0ex=.false.    !true before 12Aug2006 ---> but it caused a problem maybe because of degeneracy.
cccccccccccccccccccccccccc
        if(evec0ex) then
          ifevec0 =ifile_handle()
          open(ifevec0,file='evec0',form='UNFORMATTED',status='OLD')
          ifevecchk=ifile_handle()
          if(isigma_en==5) open(ifevecchk,file='evecfix.chk')
        endif
      endif

c$$$!! === readin lmfgw_kdivider, and get extensions === apr2013
c$$$      ifiproc=ifile_handle()
c$$$      open(unit=ifiproc,file='lmfgw_kdivider',status='old')
c$$$      read(ifiproc,*) ext
c$$$      read(ifiproc,*) nqixx, nspxx, nrank
c$$$      allocate(iprocq(nqixx,nspxx))
c$$$      do isp=1,nspxx
c$$$      do iqq=1,nqixx
c$$$        read(ifiproc,*) iqqxx, ispxx, ixxx
c$$$        if(iqqxx/=iqq) call rx( 'iqqxx/=iqq')
c$$$        if(ispxx/=isp) call rx( 'ispxx/=isp')
c$$$        iprocq(iqq,isp) = ixxx
c$$$        write(6,"('iqq isp irank=',i8,i2,i6)") iqq,isp, iprocq(iqq,isp)
c$$$      enddo
c$$$      enddo
c$$$      close(ifiproc)
c$$$!! for multiple files
c$$$      allocate(extp(0:nrank-1),ifvxc__(0:nrank-1),ifevec__(0:nrank-1))
c$$$      extp(0) = trim(ext)
c$$$      write(6,"('  0 ext= ',a,a)") trim(extp(0)),' ----------'
c$$$      do procid=1,nrank-1
c$$$        write(extn,"(i10)") procid
c$$$        extp(procid)=trim(adjustl(ext))//'_'//trim(adjustl(extn))
c$$$        write(6,"(i3,' ext= ',a,a)") procid,trim(extp(procid)),' ----------'
c$$$      enddo
c$$$      do procid=0,nrank-1
c$$$        ixx=ifile_handle()
c$$$c        open(ixx,file='vxc'//extp(procid),form='unformatted')
c$$$        ifvxc__(procid)=ixx
c$$$        ixx=ifile_handle()
c$$$c        open(ixx,file='evec'//extp(procid),form='unformatted')
c$$$        ifevec__(procid)=ixx
c$$$      enddo
      
      call read_BZDATA()
      write(6,*)' read from bzdata nqibz2; nqibz nq nhq=',nqibz2,nq,nhq
!!---------------------------------------------------------------------
c      ifvxc_ = ifvxc__(0)
c      ifevec_= ifevec__(0)
c      read(ifvxc_) ndimh, nsp,nnn !sf.beg
c      if(evec0ex) read(ifevec0) ndimhx0, nspx0,nnnx0
c      read(ifevec_) ndimhx, nspx,nnnx
c      nhq=ndimhx
      nsp=nspin
      nnn=nqibz2
      write(6,*)' ndimh ntq nsp nnn =',nhq,ntq,nsp,nnn !NOTE this ndimh is the maximum dimention of Hamiltonian.
      !In PMT, ndimh is q-dependent. See lm*/gwd/sugw.F  june2009 takao
c      if (nspin .ne. nsp)   call rx( 'hqpe: wrong nspin v_xc')
c      if (nspin .ne. nspx)  call rx( 'hqpe: wrong nspin evec')
c      if (ndimhx /= ndimh)   call rx( 'hqpe: SC ndimhx/= ndimx')
c      if (nnn .ne. nnnx) call rx( 'hqpe: wrong nnn evec or v_xc')
      allocate(qqq(3,nnn,nspin)) !,se_in(ntq,ntq))                 !sf.beg
      allocate(v_xc(nhq,nhq,nnn,nspin),evec(nhq,nhq,nnn,nspin),nev(nnn,nspin))
      if(evec0ex) allocate(evec00(nhq,nhq,nnn,nspin))
      allocate(nhqx(nnn,nspin))
c     write(6,*)' nnn =',nnn

      iqq=0
      do iq=1,nnn               !now nnn is not necessary to be nqbz !nnn=nqbz
        iqq=iqq+1
        do is=1,nspin
!!
c          ifvxc_  = ifvxc__ (iprocq(iqq,is)) 
c          ifevec_ = ifevec__(iprocq(iqq,is)) 
          open(newunit=ifvxc_,    file='vxc'//trim(xt(iq))//trim(xt(is)),form='unformatted')
          open(newunit=ifevec_,  file='evec'//trim(xt(iq))//trim(xt(is)),form='unformatted')
          read(ifvxc_ ) ndimh, nsp,nnn ,nnnx,nmto
          nz=ndimh
          read(ifevec_) ndimhx, nspx,nnnx
c          read(ifvxc_) nz,nmto    ! june209 takao nhqx is true dimension which is dependent on q. 
                                  ! nmto(not iq,is is dependent) is MTO dimension.!2012sep
          write(6,*) ' reading v_xc ... iq is nz=',iq,is,nz
          nhqx(iq,is) = nz   !nz is introduced instead of nhq
          read(ifvxc_) v_xc(1:nz,1:nz,iq,is)
          read(ifevec_) qqq(1:3,iq,is),evec(1:nz,1:nz,iq,is), nev(iq,is) !nev number of true bands nov2015
          if(evec0ex) read(ifevec0) qqqx0(1:3), evec00(1:nz,1:nz,iq,is)
c          if (sum(abs( qqq(:,iq,is )-qqqx(:) )) .gt. tolq)
c     &         call rx( 'hqpe: qqq/=qqqx')
          close(ifvxc_)
          close(ifevec_)
         enddo
      enddo                     !sf.end
      print *,' end of reading vxc evec'
c      allocate(sigmv(nhq,nhq,nq))
      if(mtosigmaonly()) then
           ndimsig = nmto
      else
           ndimsig = nhq
      endif
      if (lsigin .and. ndimsigin .ne. ndimsig) then
          deallocate(sigin)
          write(6,*) '... input sigma dimension mismatch ... discarding'
          lsigin = .false.
      endif
      allocate(sigmv(ndimsig,ndimsig,nq))

      if_eseavr = iopenxx ('ESEAVR')
      
c----------------------------------------
 9998 continue
c----------------------------------------

!! ------------------------------------------------------
c meanvalue of eseavr june2009.
c---  readin BZDATA. See gwsrc/rwbzdata.f
c--------readin data set when you call read_BZDATA ---------------
c     integer(4)::ngrp,nqbz,nqibz,nqbzw,nteti,ntetf
ccccc ! &   ,n_index_qbz
c     integer(4):: n1,n2,n3
c     real(8):: qbas(3,3),ginv(3,3),qbasmc(3,3)
c     real(8),allocatable:: qbz(:,:),wbz(:),qibz(:,:)
c     &    ,wibz(:),qbzw(:,:)
c     integer(4),allocatable:: idtetf(:,:),ib1bz(:),idteti(:,:)
c     &    ,nstar(:),irk(:,:),nstbz(:)          !,index_qbz(:,:,:)
c-----------------------------------------------------------------
      nstarsum= sum(nstar(:))
      write(6,*)' nstarsum=',nstarsum
Cstop2rx 2013.08.09 kino      if( nstarsum/= nqbz ) stop ' nstarsum/= nqbz'
      allocate(ntqxx(nq))


!!    loop over spin -------------------------------
      do 1001  is = 1,nspin
        write(6,*) ' --- is=',is
        call readx   (ifsex(is),50)
        read (ifsex(is),*) nspinx,nqx,ntqx
        read (ifsex(is),*)
        read (ifsex(is),*) deltaw
        read (ifsex(is),*) alat
        read (ifsex(is),*) ef
        read (ifsex(is),*) esmr
        write(6,*)' xxxxxxxxxxx111111111'
        call readx   (ifxc(is),50)
        read (ifxc(is),*) nspinxc,nqxc,ntqxc
        read(ifsex2(is))  nspinx2, nqx2, ntqx2, nqbz ,  nqibz,  n1,n2,n3 !sf.beg

        if( nstarsum/= nqbz ) call rx( ' nstarsum/= nqbz')

        if(.not.exonly)  then
          call readx(ifsec(is),50)
          read (ifsec(is),*) nspinc,nqc,ntqc
          read(ifsec2(is))  nspinc2, nqc2, ntqc2, nqbzc2, nqibzc2,n1,n2,n3 !sf.end
        endif
        read(ifsexcore2(is)) nspinxc2,nqxc2,ntqxc2,nqbzxc2,nqibzxc2
        if (nspin .ne. nspinx)   call rx( 'hqpe: wrong nspin SEx')
        if (nspin .ne. nspinxc)  call rx( 'hqpe: wrong nspin vxc')
        if (nspin .ne. nspinx2)  call rx( 'hqpe: wrong nspin SEx2')
        if (nq .ne. nqx)         call rx( 'hqpe: wrong nq SEx')
        if (nq .ne. nqxc)        call rx( 'hqpe: wrong nq vxc')
        if (nq .ne. nqx2)        call rx( 'hqpe: wrong nq SEx2')
        if (ntq .ne. ntqx)       call rx( 'hqpe: wrong ntq SEx')
        if (ntq .ne. ntqxc)      call rx( 'hqpe: wrong ntq vxc')
        if (ntq .ne. ntqx2)      call rx( 'hqpe: wrong ntq SEx2')
        if(.not.exonly) then
          if (nqbz .ne.  nqbzxc2)  call rx( 'hqpe: wrong nqbzx2')
          if (nqibz .ne. nqibzxc2) call rx( 'hqpe: wrong nqibzxc2')
          if (nspin .ne. nspinxc2) call rx( 'hqpe: wrong nspin SExcore2')
          if (nspin .ne. nspinc)   call rx( 'hqpe: wrong nspin SEc')
          if (nspin .ne. nspinc2)  call rx( 'hqpe: wrong nspin SEc2')
          if (nq .ne. nqc)         call rx( 'hqpe: wrong nq SEc')
          if (nq .ne. nqc2)        call rx( 'hqpe: wrong nq SEc2')
          if (nqbz .ne. nqbzc2)    call rx( 'hqpe: wrong nqbzxc2')
          if (ntq .ne. ntqc2)      call rx( 'hqpe: wrong ntq SEc2')
          if (nq .ne. nqxc2)       call rx( 'hqpe: wrong nq SExcore2')
          if (ntq .ne. ntqxc2)     call rx( 'hqpe: wrong ntq SExcore2')
        endif
        if(is==1) write(6,*)' ###  readin XCU'
        if(is==2) write(6,*)' ###  readin XCD'
        allocate( itxc(ntq),qxc(3,ntq,nq),eldaxc(ntq,nq),vxc(ntq,nq) )
        call readx (ifxc(is),50)
        read(ifxc(is),*)
        do ip = 1,nq
          do i  = 1,ntq
            read(ifxc(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16)")
     &              itxc(i),ipxx,isxxx, qxc(1:3,i,ip), eldaxc(i,ip), 
     &              vxc(i,ip)
          enddo
        enddo
c
        if(is==1) write(6,*)' ###  readin SEXU'
        if(is==2) write(6,*)' ###  readin SEXD'
        allocate( itx(ntq), qx (3,ntq,nq),eldax (ntq,nq),sex(ntq,nq) )
        allocate( sex2(ntq,ntq,nq)) !sf..3June
        call readx   (ifsex(is),50)
        read(ifsex(is),*)
        do ip = 1,nq
          do i  = 1,ntq
            read(ifsex(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16)")
     &              itx(i),ipxx,isxxx, qx(1:3,i,ip), eldax(i,ip), 
     &              sex(i,ip)
          enddo
          read(ifsex2(is)) isx,qx2,sex2(1:ntq,1:ntq,ip) !sf..3June
        enddo
c
        allocate( sexcore(ntq,nq),sexcore2(ntq,ntq,nq) ) !sf..3June
c     if(exonly) then !zero for exonly case
c     write(6,*)' set sexcore=0 for exonly case'
c     sexcore=0d0; sexcore2=0d0
        if(is==1) write(6,*)' ###  readin SEXcoreU'
        if(is==2) write(6,*)' ###  readin SEXcoreD'
        call readx   (ifsexcore(is),50)
        call readx   (ifsexcore(is),50)
        read(ifsexcore(is),*)
        do ip = 1,nq
          do i  = 1,ntq
            read(ifsexcore(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16)")
     &              ixx1,ixx2,ixx3, qxxx1,qxxx2,qxxx3, exxx, sexcore(i,ip)
          enddo
          read(ifsexcore2(is)) isx,qx2,sexcore2(1:ntq,1:ntq,ip) !sf..3June
        enddo
c
        allocate( itc(ntq), qc (3,ntq,nq),eldac (ntq,nq)
     &        ,rsec(ntq,nq),csec(ntq,nq),sec2(ntq,ntq,nq)) !sf..3June
        if(exonly) then        !zero for exonly case
          write(6,*)' set sec=0 for exonly case'
          itc=0d0;qc=0d0;eldac=0d0;rsec=0d0;csec=0d0;sec2=0d0
        else
          if(is==1) write(6,*)' ###  readin SECU'
          if(is==2) write(6,*)' ###  readin SECD'
          call readx   (ifsec(is),50)
          read(ifsec(is),*)
          do ip = 1,nq
            do i  = 1,ntq
              read(ifsec(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16, 
     &                 3x,d24.16)") !sf13May2002     ,3x,d24.16)")
     &                 itc(i),ipxxx,isxxx, qc(1:3,i,ip), eldac(i,ip), 
     &                 rsec(i,ip),csec(i,ip) !sf..13May2002   ,zfac(i,ip)
            enddo
            read(ifsec2(is)) isx,qx2,sec2(1:ntq,1:ntq,ip) !sf..3June
          enddo
        endif

c     check that LDA eigenvalues are the same
c     mod 93.06.29
c     chev       = compvv (w(ieldax),w(ieldac),ntq*nq,tol,'eLDA x,c')
c     chev       = compvv (w(ieldax),w(ieldaxc),ntq*nq,tol,'eLDA x,xc')
c     chq        = compvv (w(iqx),w(iqc),3*nq,tol,'q-vec x,c')
c     chq        = compvv (w(iqx),w(iqxc),3*nq,tol,'q-vec x,xc')
c-----
c     itxc       = icompvv (w(iitx),w(iitc),ntq,'state-label x,c')
c     itxc       = icompvv (w(iitx),w(iitxc),ntq,'state-label x,xc')
        if(.not.exonly)itxcxxx=icompvv2(itx, itc,ntq,' state-label x c')
        if(.not.exonly)itxcxxx=icompvv2(itx,itxc,ntq,' state-label x xc')

c     calculate quasiparticle energies and write them on file QPE
c     call qpe     (ifqpe(is),iftote,w(iitc),w(iqc),
c     i              w(ieldac),w(ivxc),w(isex),
c     i              w(irsec),w(icsec),w(izfac),eshift,
c     d              ntq,nq,is,
c     o              w(ieqp),w(iwqp))
        call qpe1_sc  (ifqpe(is),iftote(is),iftote2(is),itx,qx, !itc,qc
     i        eldax,vxc,sex,sexcore,
     i        rsec,csec,jin,deltaw,alat,ef, !sf..13May |it was zfac after csec
     d        ntq,nq,is,
     o        eshift0,eshift02,eshtlda)
c     o              w(ieqp),w(iwqp))
        write(6,*)' end of qpe1'
        close(ifqpe(is))
        close(iftote(is))
        close(iftote2(is))
        if(jin==-101) cycle

!------------------------------------------------------------------------------
ccccccbegin making SE_ij-VXC_ij, ij are 'basis functions' indexes !sf..3June
!------------------------------------------------------------------------------
        allocate(se(ntq,ntq,nq),ipiv(nhq),work(nhq*nhq)
     &        ,evec_inv(nhq,nhq) ,evec_invt(nhq,nhq))
        allocate(ev_se_ev(ndimsig,ndimsig)) !  sep2012     allocate(ev_se_ev(nhq,nhq)) 
        if(evec0ex) allocate(evec00inv(nhq,nhq))
        do ip=1,nq
          do itp=1,ntq
            do itpp=1,ntq    !make Sigma hermitean
              if(.not.exonly) then
                se(itpp,itp,ip)=sex2(itpp,itp,ip)+sexcore2(itpp,itp,ip)
     &                    +.5d0*(sec2(itpp,itp,ip)+dconjg(sec2(itp,itpp,ip)))
              else
                se(itpp,itp,ip)= wex*sex2(itpp,itp,ip)+sexcore2(itpp,itp,ip)
ccc   se(itpp,itp,ip)= wex*sex2(itpp,itp,ip)
              endif
cccccccccccccccccccccccc
c     c takao test ex only test
c     se(itpp,itp,ip)= 0.3*sex2(itpp,itp,ip)+ sexcore2(itpp,itp,ip)
ccccccccccccccc
c     c takao test2---reduce screen
c     se(itpp,itp,ip)=sex2(itpp,itp,ip)+sexcore2(itpp,itp,ip)+
c     &        .1d0*(sec2(itpp,itp,ip)+dconjg(sec2(itp,itpp,ip)))
cccccccccccccc
            enddo
          enddo
        enddo
c     write(6,*)' xxx1'
        do ip=1,nq
          do itp=1,ntq
            do itpp=1,ntq    !make Sigma hermitian
              if(abs(se(itpp,itp,ip)-dconjg(se(itp,itpp,ip)) ) .gt. 1d-8) then
                write(*,*)'diff=', se(itpp,itp,ip)-dconjg(se(itp,itpp,ip))
                write(*,*)'se   ', se(itpp,itp,ip), se(itp,itpp,ip)
                write(*,*)'sex2 ',sex2(itpp,itp,ip),sex2(itp,itpp,ip)
                write(*,*)'sexc2',sexcore2(itpp,itp,ip),sexcore2(itp,itpp,ip)
                call rx( "hqpe: Sigma_nn' is not hermitian")
              endif
            enddo
          enddo
        enddo
        if (is==1.and.laf) then
           write(ifse_out) 1,ndimsig,n1,n2,n3,nq,0,0,0
        elseif (is == 1) then
           write(ifse_out) nspin,ndimsig,n1,n2,n3,nq,0,0,0
        endif
        do 2001 ip=1,nq
c$$$c sep2013t ccccccccccccccccccccccccccccc
c$$$           if(is==1) if_ntqxxd=iopenxx('NTQXXD.up.'//charnum5(ip))
c$$$           if(is==2) if_ntqxxd=iopenxx('NTQXXD.dn.'//charnum5(ip))
c$$$           read(if_ntqxxd,*) qttt(1:3),ntqxxd
c$$$           close(if_ntqxxd)
c$$$           if( sum((qx(1:3,1,ip)-qttt(1:3))**2 ) < 1d-10 ) then
c$$$              continue
c$$$           else
c$$$              call rx('hqpe_sc:xxxx:sum((qx(1:3,1,ip)-qttt(1:3))**2')
c$$$           endif
c$$$           write(6,*)' ntqxxd=',ip,qttt(1:3),ntqxxd
c$$$cccccccccccccccccccccccccccccccccccccccc
          do ikp=1,nnn        !nqbz
            if (sum ( (qqq(1:3,ikp,is)-qx(1:3,1,ip))**2 ) < tolq ) then
              ikpx=ikp      !qc(:,i,:) does not depents on band index i=1:ntq
              goto 100
            endif
          enddo               !ikp
          call rx( 'hqpe.sc: not find ikp 100')
 100      continue
          nz = nhqx(ikpx,is) !june 2009 takao. we use nz instead of nhq
c$$$          evec_inv(1:nz,1:nz)=evec(1:nz,1:nz,ikpx,is)
c$$$          !evec(i,n), i-basis, n - band argumants
c$$$c---  make inverse evec_inv(n,i) matrix \psi_n=sum_i evec(i,n)\phi_i,
c$$$c---  where \psi is eigenfunction and \phi is basis function
c$$$          call zgetrf(nz,nz,evec_inv,nhq,ipiv,info)  !nhq is leading dimension of evec_inv
c$$$          if (info /= 0)  call rx( 'hqpe: zgetrf info/=0')
c$$$          !inverse matrix evec_inv(n,i), i-basis, n - band argumants
c$$$          call zgetri(nz,evec_inv,nhq,ipiv,work,nhq**2,info)
c$$$          if (info /= 0)  call rx( 'hqpe: zgetri info/=0')

C---  takao
          if(evec0ex.and.iSigma_en==5) then
            call rx('Not support evec0ex.and.iSigma_en==5 now... sep2013')
c$$$            evec00inv(1:nz,1:nz)=evec00(1:nz,1:nz,ikpx,is)
c$$$            !evec(i,n), i-basis, n - band argumants
c$$$            call zgetrf(nz,nz,evec00inv,nhq,ipiv,info)
c$$$Cstop2rx 2013.08.09 kino            if (info /= 0)  stop 'hqpe: zgetrf 00 info/=0'
c$$$            if (info /= 0)  call rx( 'hqpe: zgetrf 00 info/=0')
c$$$            call zgetri(nz,evec00inv,nhq,ipiv,work,nhq**2,info)
c$$$            !inverse matrix evec_inv(n,i), i-basis, n - band argumants
c$$$Cstop2rx 2013.08.09 kino            if (info /= 0)  stop 'hqpe: zgetri 00 info/=0'
c$$$            if (info /= 0)  call rx( 'hqpe: zgetri 00 info/=0')
c$$$            itppxs=0
c$$$            do itp=1,nz
c$$$              sum0=0d0
c$$$              sum1=0d0
c$$$              do itpp=1,nz
c$$$c     write(6,*)evec_inv(itp,itpp)
c$$$                aaa=abs(sum(evec_inv(itp,1:nz)*evec00(1:nz,itpp,ikpx,is)))**2
c$$$                if(aaa>sum1) then
c$$$                  sum1 =aaa
c$$$                  itppx=itpp
c$$$                endif
c$$$                sum0=sum0+aaa
c$$$              enddo
c$$$              evec_inv(itp,1:nz)    =  evec00inv(itppx,1:nz)
c$$$              evec(1:nz,itp,ikpx,is)=  evec00(1:nz,itppx,ikpx,is)
c$$$              itppxs = itppxs+itppx
c$$$              write(ifevecchk,"(3i5,2d13.5)")itp,ikpx,is, sum0, sum1
c$$$c     if(abs(sum0-1d0)>0.1) stop "!!! sum0 is not unity ==="
c$$$            enddo
c$$$            if(nz*(nz+1)/2 /= itppxs) then
c$$$              write(6,*)'nz_sum itppxs=', nz*(nz+1)/2,itppxs
c$$$Cstop2rx 2013.08.09 kino              stop "hqpe: nz*(nz+1)/2 /= itppxs"
c$$$              call rx( "hqpe: nz*(nz+1)/2 /= itppxs")
c$$$            endif
          endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     takao test only diagonal
c     call onlydiag(se(1:ntq,1:ntq,ip),ntq)
c     do itp=1,ntq
c     se(itp,itp,ip)=se(itp,itp,ip) - 0.5*vxc(itp,ip)/rydberg()
c     enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          ntqxx(ip)=ntq
          do itp=ntq,1,-1
            if(se(itp,itp,ip)/=0d0) then
              ntqxx(ip) = itp
              exit
            endif
          enddo
          write(6,*)
          write(6,"( ' ip ntq ntqxx=',3i5)")ip,ntq,ntqxx(ip)
!!
          if(ip==1) then
            ntqxxmin=ntqxx(ip)
          else
            if(ntqxx(ip)<ntqxxmin) ntqxxmin=ntqxx(ip)
          endif

          do itp = 1,ntqxx(ip)
            do itpp= 1,ntqxx(ip)
              if(itp/=itpp.and.isigma_en==5) then
                se(itp,itpp,ip)= 0d0
              else
c     Aug2004
                if(.not.exonly) then
                  se(itp,itpp,ip)= se(itp,itpp,ip)
     &                       -.5d0* sum(dconjg(evec(1:nz,itp,ikpx,is))*
     &                       matmul(v_xc(1:nz,1:nz,ikpx,is),evec(1:nz,itpp,ikpx,is)))
                else
                  se(itp,itpp,ip)= se(itp,itpp,ip)
     &                       -wex*.5d0* sum(dconjg(evec(1:nz,itp,ikpx,is))*
     &                       matmul(v_xc(1:nz,1:nz,ikpx,is),evec(1:nz,itpp,ikpx,is)))
                endif
!!
              endif
            enddo
          enddo
 2001   continue

        do 2002 ip=1,nq
!     e-weighted average
!     eavr  : average of eigenvalues within threshold (itp<=ntqxx)
!     eavr2  : square average of eigenvalues within threshold (itp<=ntqxx)
!     eseavr: average of se*eigenvalue
           eavr  = 0d0
           eavr2  = 0d0
           eseavr0= 0d0
           eseavr02 =0d0
           iix=0
           do itp=1,ntqxx(ip)
              eee = eldax(itp,ip) - rydberg()*ef
              if( eee > 1d-2 ) then
                 eavr   = eavr   + eee
                 eavr2   = eavr2   + eee**2
                 eseavr0 = eseavr0 + eee* se(itp,itp,ip)
                 eseavr02 = eseavr02 + eee**2* se(itp,itp,ip)
                 iix=iix+1
              endif
           enddo
!! jun2014
           if(iix==0) then !if iix=0
             eavr2=1d0                              
             eseavr02=0d0                           
           endif                                    
           eseavr(ip,is) = eseavr02/eavr2 !now eseavr is
c          eseavr(ip,is) = eseavr0/eavr !now eseavr is
cc           do itp = ntqxx+1,ntq
cc              se(itp,itp,ip) = eseavr(ip,is)
cc           enddo
           eseavr(ip,is) = 2d0*eseavr(ip,is) !in Ry.
           write(6,*)"### A correction takao2009June: find this in hqpe.se.m.F"
           write(6,*)"###   constant is added to sigm above threshold."
           write(6,*)"###   the constant (ESEAVR=e-weighted average Ry)= ",is,ip,eseavr(ip,is)
 2002   continue
        eseavrmean=0d0
        do ip=1,nq
          eseavrmean = eseavrmean + nstar(ip)*eseavr(ip,is)
        enddo
        eseavrmean = eseavrmean/nqbz
        call getkeyvalue("GWinput","AddToESEAVR",eseadd,default=0d0,status=ret)
        eseavrmean=  eseavrmean+eseadd
cccccccccccccccccccccccccccccccc
c        eseavrmean = 0d0
cccccccccccccccccccccccccccccccc
        write(if_eseavr,"(d23.15,i3,i8)") eseavrmean,is ,ntqxxmin
c       write(if_eseavr,"(d23.15,i3)") 0d0,is
        write(6,"(' ESEAVRmean (used bands above emax_sigm) isp=',d13.6,i2)")eseavrmean,is
!!-----------------------------------------------------------------
        do 2003 ip=1,nq
          do ikp=1,nnn        !nqbz
            if (sum ((qqq(1:3,ikp,is)-qx(1:3,1,ip))**2 ) < tolq ) then
              ikpx=ikp      !qc(:,i,:) does not depents on band index i=1:ntq
              goto 102
            endif
          enddo               !ikp
          call rx( 'hqpe.sc: not find ikp 102')
 102      continue
          nz = nhqx(ikpx,is) !june 2009 takao. we use nz instead of nhq


c!!  Make inverse evec_inv(n,i) matrix \psi_n=sum_i evec(i,n)\phi_i,
c!!  where \psi is eigenfunction and \phi is basis function
!! nov2015 
!! evec_inv(ib1,iww)= \sum_ib2  ovlinv(ib1,ib2)*dconjg(evec(iww,ib2))  nov2015, we introduce nev. iww is for PMT basis. ib for band index.
!! This is for converting rotated evec (=evecrot(ib)) in the representation of original evec(ib).
          nevv=nev(ikpx,is)
          allocate(ovl(nevv,nevv)) 
c          print *,'nnnnnnnnn zsecsym: nband=',nhdim,nband,nev
          do i=1,nevv
          do j=1,nevv
c            write(6,*)'evec orth=',i,j,sum(dconjg(evec(:,i)*evec(:,j)))
            ovl(i,j)=sum(dconjg(evec(1:nz,i,ikpx,is))*evec(1:nz,j,ikpx,is))
          enddo
          enddo
          call matcinv(nevv,ovl) !ovl --> ovlinv
          evec_inv(1:nevv,1:nz) = matmul(ovl(1:nevv,1:nevv),dconjg(transpose(evec(1:nz,1:nevv,ikpx,is)))) !note ovl means ovlinv
          deallocate(ovl)
c!! ccccccccccccccccccccccccc before nov2015
c          evec_inv(1:nz,1:nz) = evec(1:nz,1:nz,ikpx,is) !evec(i,n), i-basis, n - band argumants
c          call zgetrf(nz,nz,evec_inv,nhq,ipiv,info)   !nhq is leading dimension of evec_inv
c          if (info /= 0)  call rx( 'hqpe: zgetrf info/=0')
c          !inverse matrix evec_inv(n,i), i-basis, n - band argumants
c          call zgetri(nz,evec_inv,nhq,ipiv,work,nhq**2,info)
c          if (info /= 0)  call rx( 'hqpe: zgetri info/=0')


!!
          evec_invt(1:nz,1:nevv) = transpose(dconjg(evec_inv(1:nevv,1:nz)))
          if(mtosigmaonly()) then
             ndimsig2= nmto
          else
             ndimsig2 = nz
          endif   
          ev_se_ev(1:ndimsig2,1:ndimsig2) = matmul(evec_invt(1:ndimsig2,1:ntqxx(ip))
     &    ,matmul(se(1:ntqxx(ip),1:ntqxx(ip),ip),evec_inv(1:ntqxx(ip),1:ndimsig2)))

!! sep2013 exprapolation of se by eseavrmean
          do itp =1,ndimsig2
          do itpp=1,ndimsig2
            ev_se_ev(itp,itpp) = ev_se_ev(itp,itpp) +
     &      sum(evec_invt(itp,ntqxx(ip)+1:nevv)*evec_inv(ntqxx(ip)+1:nevv,itpp))*eseavrmean/2d0 ! in Hartree.
          enddo
          enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc11 sep2013 exprapolation of se by eseavrmean                                               
c          do itp =1,ndimsig2                                                                
c          do itpp=1,ndimsig2                                                                
c            ev_se_ev(itp,itpp) = ev_se_ev(itp,itpp) +                                       
c     &      sum(evec_invt(itp,ntqxx(ip)+1:ntq)*evec_inv(ntqxx(ip)+1:ntq,itpp))*eseavrmean/2d0
c          enddo                                                                             
c          enddo                                                                             
cc12 sep2013 exprapolation of se by eseavrmean                         
c          do itp =1,ndimsig2                                          
c          do itpp=1,ndimsig2                                          
c            ev_se_ev(itp,itpp) = ev_se_ev(itp,itpp) +                 
c     &      sum(evec_invt(itp,ntq+1:nz)*evec_inv(ntq+1:nz,itpp))* 5d0 
c          enddo                                                       
c          enddo                                                       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!! - write SE_ij-Vxc_ij where ij are basis function indices
          do itp=1,ndimsig2
            do itpp=1,ndimsig2
              if( abs(ev_se_ev(itpp,itp) - dconjg(ev_se_ev(itp,itpp)) )> 1d-8  ) then
                write(6,*)itp,itpp
                write(6,*)ev_se_ev(itpp,itp)
                write(6,*)ev_se_ev(itp,itpp)
                call rx( 'hqpe: Sigma_ij is not hermitian')
              endif
              if(abs(v_xc(itp,itpp,ikpx,is)-dconjg(v_xc(itpp,itp,ikpx,is)))> 1d-8) 
     &          call rx( 'hqpe: v_xc is not hermitean')
            enddo
          enddo
          write(ifse_out) qqq(1:3,ikpx,is),is !,eseavr(ip,is) !,ntqxxmin ! ntqxx(ip)  sep2013
          sigmv(:,:,ip) = 1d20
          sigmv(1:ndimsig2,1:ndimsig2,ip) =  2d0*ev_se_ev(1:ndimsig2,1:ndimsig2) !in Ry.
          write(6,*) 'ssssss1=', sum(abs(sigmv(:,:,ip)))
          write(ifse_out) sigmv(:,:,ip)

!!      2*ev_se_ev bacause v_xc in sugw.f was in rydberg while SE was in hartree
c     endif !mixing_factor
c$$$c--   check diagonal elements (SE-Vxc)_nn, where n is band index
c$$$          allocate(sen(nhq,nhq),sen2(nhq,nhq),evec0(nhq,nhq))
c$$$          evec0(1:ndimsig,1:ndimsig)=evec(1:ndimsig,1:ndimsig,ikpx,is)
c$$$          do itp =1,ndimsig
c$$$            do itpp=1,ndimsig
c$$$              sen(itpp,itp)= ev_se_ev(itpp,itp)*hartree
c$$$            enddo
c$$$          enddo
c$$$          call zgemm('n','n',ndimsig,ndimsig,ndimsig,(1d0,0d0),sen,
c$$$     &           nhq,evec0,nhq,(0d0,0d0), sen2 ,nhq)
c$$$          call zgemm('c','n',ndimsig,ndimsig,ndimsig,(1d0,0d0),evec0,
c$$$     &           nhq,sen2,nhq,(0d0,0d0), sen ,nhq)
c$$$          do i=1,ntq
c$$$            del=abs(rsec(i,ip)+sex(i,ip)+sexcore(i,ip)-vxc(i,ip)-
c$$$     &              sen(i,i))
c$$$          enddo               !i
c$$$          deallocate(sen,sen2,evec0)
c$$$c--   end check diagonal elements
 2003   continue                    !ip
        write(6,*)
        deallocate(sex2,sexcore2,sec2,se,ipiv,work)
        deallocate(evec_inv,evec_invt,ev_se_ev) !,se_ev
        if(evec0ex) deallocate(evec00inv)
! - end making SE_ij-VXC_ij cccccccccccccccccccccccccccccccccc  !sf..3June
        deallocate( itxc,qxc,eldaxc,vxc ,itc, qc ,eldac,
     &        sexcore ,rsec,csec, !zfac,  !sf..13May2002 
     &        itx, qx ,eldax,sex)
        if (jin .gt. 0) jin = 999999
        if (laf) exit
 1001 continue !isp

c-------------------------------
      if(jin==-101) goto 9999
      deallocate(v_xc,evec)     !, se_in)                                        !sf..3June
      close(ifse_out)

C     --- mixing sigma -------------------------------------
C     ... Once readin sigma
      ifse_out=ifile_handle()
      open(UNIT=ifse_out, file='sigm',form='UNFORMATTED')
      allocate(sigma_m(ndimsig,ndimsig,nq,nspin),qqqx_m(3,nq,nspin))
      write(6,*)"========= Sigma mixing section using mixsigma ======="
      call rwsigma (ifse_out,sigma_m,qqqx_m, nspin,ndimsig,n1,n2,n3,nq ) !,eseavr_in)
C     write(6,*)"==== Have read sigm and goto mixsigma ===="
C     MvS Feb 04 Allocate sigin if not already done so
      if (lsigin .and. ndimsigin .ne. ndimsig) then
        deallocate(sigin)
        write(6,*) '... input sigma dimension mismatch ... discarding'
        lsigin = .false.
      endif
      if (.not. lsigin) then
        allocate(sigin(1,1,1,1))
        sigin = 0d0
      endif
      call mixsigma(sigma_m, lsigin, sigin, 2*ndimsig**2*nq*nspin)

C     ... write sigm Binary format ---------------------------------------
      write(6,*)
      write(6,*) "=== Write sigma to files (Binary and ASCII) ========"
      rewind ifse_out
      call rwsigma (-ifse_out,sigma_m,qqqx_m, nspin,ndimsig,n1,n2,n3,nq) !,eseavr)
      close(ifse_out)

C     ... Write sigma ASCII format ---------------------------------------
      ifse_out=ifile_handle()
      open(UNIT=ifse_out, file='sigma')
      call rwsigmaa (-ifse_out,sigma_m,qqqx_m, nspin,ndimsig,n1,n2,n3,nq) !,eseavr)
      close(ifse_out)
 9999 continue
      call rx0s( ' OK! hqpe_sc ')
      end

c------------------------------------------------------
      subroutine testfff(a,b,nnn)
      implicit integer (i-n)
      real(8):: a(nnn),b(nnn)
      do i=1,nnn
        if(i/=nnn) then
          b(i)= +  (-1d0)**i*a(i+1) + i/10d0
        else
          b(i)=   a(i)**2 - 1d0
        endif
      enddo
      end

c------------------------------------------------------
      subroutine onlydiag(sigm,ndimh)
      implicit integer (i-n)
      complex(8):: sigm(ndimh,ndimh)
      do ix=1,ndimh
        do iy=1,ndimh
          if(ix/=iy) sigm(ix,iy)=0d0
        enddo
      enddo
      end

C----------------------------------------------------------------------
      subroutine rwsigma(ifs,sigma_m,qqqx_m,nspin,ntq,n1,n2,n3,nq ) !,eseavr)
      implicit integer (i-n)
      complex(8)::sigma_m(ntq,ntq,nq,nspin)
      real(8)::qqqx_m(3,nq,nspin) !,eseavr(nq,nspin)
      ntm1=0;ntm2=0;ntm3=0
      if(ifs>0) write(6,*) " rwsigma: Reading sigm (Binary)==="
      if(ifs<0) write(6,*) " rwsigma  Writing sigm (Binary)==="
      if(ifs>0) read ( ifs) nspin,ntq,n1,n2,n3,nq, ntm1,ntm2,ntm3
      if(ifs<0) write(-ifs) nspin,ntq,n1,n2,n3,nq, ntm1,ntm2,ntm3
      do is=1,nspin
        do ip=1,nq
          if(ifs>0) read (ifs)  qqqx_m(1:3,ip,is)!,isr !,eseavr(ip,is)
          if(ifs>0) read (ifs)  sigma_m(1:ntq,1:ntq,ip,is)
          if(ifs<0) write(-ifs) qqqx_m(1:3,ip,is),is !,eseavr(ip,is)
          if(ifs<0) write(-ifs) sigma_m(1:ntq,1:ntq,ip,is)
          write(6,"('  === ',i5,i3,3f10.5)")ip,is,qqqx_m(1:3,ip,is)
        enddo
      enddo
      write(6,*)' === rwsigma:  sum check of sigma_m=',sum(abs(sigma_m))
      print *
      end
c
      subroutine rwsigmaa(ifs,sigma_m,qqqx_m, nspin,ntq,n1,n2,n3,nq) !,eseavr)
      implicit integer (i-n)
      complex(8)::sigma_m(ntq,ntq,nq,nspin)
      real(8)::qqqx_m(3,nq,nspin) !,eseavr(nq,nspin)
      ntm1=0;ntm2=0;ntm3=0
      if(ifs>0) write(6,*) " rwsigmaa: Reading sigma (ASCII)==="
      if(ifs<0) write(6,*) " rwsigmaa: Writing sigma (ASCII)==="
      if(ifs>0)read ( ifs,*)      nspin,ntq,n1,n2,n3,nq,ntm1,ntm2,ntm3
      if(ifs<0)write(-ifs,"(i10)")nspin,ntq,n1,n2,n3,nq,ntm1,ntm2,ntm3
      do is=1,nspin
        do ip=1,nq
          if(ifs>0) read (ifs,*)  qqqx_m(1:3,ip,is) !,isr!,eseavr(ip,is)
          if(ifs<0) write(-ifs,"(3d24.16,i3)") qqqx_m(1:3,ip,is),is!,eseavr(ip,is)
          do j=1,ntq
            do i=1,ntq
              if(ifs>0) read (ifs,"(2d24.16)")  sigma_m(i,j,ip,is)
              if(ifs<0) write(-ifs,"(2d24.16,2i5)") sigma_m(i,j,ip,is),i,j
            enddo
          enddo
          write(6,"('  === ',i5,i3,3f10.5)")ip,is,qqqx_m(1:3,ip,is)
        enddo
      enddo
      write(6,*)' === rwsigmaa: sum check of sigma_m=',sum(abs(sigma_m))
      print *
      end

C----------------------------------------------------------------------
      subroutine mixsigma(sss, lsigin, sigin, nda)
c  subroutine pqmixa(nda,nmix,mmix,mxsav,beta,rms2,a,tj)
C- Mixing routine for sigma. Modified from pqmixa in subs/pqmix.f
C- Anderson mixing of a vector
Ci  mmix: number of iterates available to mix
Cio nmix: nmix > 0: number of iter to try and mix
Ci        nmix < 0: use mmix instead of nmix.
Co  nmix: (abs)  number of iter actually mixed.
Co        (sign) <0, intended that caller update nmix for next call.
C  MvS Feb 04 use sigin as input sigma if available (lsigin=T)
C             Add mixnit as parameter
      use m_keyvalue,only: getkeyvalue

      implicit none
      logical lsigin
      integer nda,nmix,mmix
      integer(4),parameter:: mxsav=10
      double precision rms2,tj(mxsav),beta
      integer im,imix,jmix,amix,onorm,okpvt,oa
      integer iprintxx,ifi,nitr,ndaf
      real(8)::sss(nda),sigin(nda)
      real(8):: tjmax
      real(8),allocatable::norm(:),a(:,:,:)
      integer(4),allocatable:: kpvt(:)
      integer(4)::ret
      character*8 fff
      logical fexist
      real(8):: acc
      integer(4):: ido,ifile_handle
      iprintxx = 30
      beta=1d0
      call getkeyvalue("GWinput","mixbeta",beta,default=1d0,status=ret)
      print '('' mixsigma: Anderson mixing sigma with mixing beta ='',f12.6)',beta
C' ... reads prior iteration INCLUDING starting sigma for current iteration
      allocate ( a(nda,0:mxsav+1,2) )
      fff="mixsigma"
      INQUIRE (FILE =fff, EXIST = fexist)
      if(fexist)      write(6,*)'... reading file mixsigma'
      if(.not.fexist) write(6,*)'... No file mixsigma'
      ifi=ifile_handle()
      open(ifi,file=fff,form='unformatted')
      if(fexist) then
        read(ifi,err=903,end=903) nitr,ndaf
        if (ndaf .ne. nda) goto 903
        read(ifi,err=903,end=903) a
        goto 902
      endif
      goto 901
  903 continue
      print 368
  368 format(5x,'(warning) file mismatch ... mixing file not read')
  901 continue
      nitr = 0
 902  continue
      a(:,0,1) = sss      !output
C     if input sigma available, use it instead of file a(:,0,2)
      if (lsigin) then
        write(6,*)'... using input sigma read from sigm file'
        a(:,0,2) = sigin  !input
      endif
      write(6,*)'sum sss=',sum(abs(sss))
C     Restrict maximum number of prior iterations
      call getkeyvalue("GWinput","mixpriorit",imix,default=9,status=ret)
      mmix = min(max(nitr-1,0),imix)
      if (mmix > mxsav) mmix = mxsav
C     this information already printed out by amix
C     write(6,*)'mixing parameters for amix are fixed in mixsigma'
C     write(6,*)'   beta       =', beta
C     write(6,*)'   tjmax      =', tjmax
C     write(6,*)'   mmix mxsav =', mmix,mxsav
      call getkeyvalue("GWinput","mixtj",acc,default=0d0,status=ret)
      if(acc/=0d0) then
        write(6,*)' readin mixtj from GWinput: mixtj=',acc
        tjmax=abs(acc)+1d-3
        if(mmix==1) then
          tj(1)=1d0
        else
          tj(1)= acc
          tj(2)= 1-acc
          mmix=2
        endif
        ido=2
      else
        tjmax=5d0
        ido=0
      endif
      allocate(norm(mxsav**2),kpvt(mxsav))
      imix = amix(nda,mmix,mxsav,ido,dabs(beta),iprintxx,tjmax,
     .  norm,kpvt,a,tj,rms2)
      deallocate(norm, kpvt)
C ... Restore PQ array, updating new x
c      call dpscop(a,w(oa),nda,1+nda*(mxsav+2),1+nda*(mxsav+2),1d0)
c      call dcopy(nda*(mxsav+2)*2,w(oa),1,a,1)
c ...
      sss = a(:,0,2)
      rewind(ifi)
      write(ifi) nitr+1,nda
      write(ifi) a
      close(ifi)
      end

c$$$!! ===================================================================
c$$$      subroutine anfsig( natom,ldima,bas, nbloch,
c$$$     i  q,           ! for q+G
c$$$     i  anfvec,iaf,  ! these are antiferro informations.
c$$$     i  zxq)  ! i/o
c$$$C- antiferro part is added to sism
c$$$c We assume that the crystal has a magnetic symmetry described by (translataion + spin flip).
c$$$c The translation is specified by a vector,
c$$$c   AFvector = anfvec(1:3)*alat, which is the true real vector in Cartesian coodinate.
c$$$c   The anfvec(1:3) is given in ANFcond file.
c$$$c
c$$$c Each mixed basis is mapped to the other mixed basis.
c$$$c E.g. the product basis B({\bf r}-{\bf a}) is mapped to
c$$$c  B({\bf r}-{\bf a}-{\bf A}) = B({\bf r}-{\bf a}'-{\bf T}_0),
c$$$c  by the translation specified by AFvec={\bf A}.
c$$$c  Here {\bf T}_0 is some crystal tralslation vector.
c$$$c  In this code you see,
c$$$c      bas(1:3,ia1)+ anfvec  = bas(1:3,iaf(ia1)) + transaf(1:3,ia1)
c$$$c  ==   {\bf a}    +{\bf A}  = {\bf a}'          + {\bf T}_0
c$$$c
c$$$c ---- The corresponding atoms should have the same product basis.
c$$$c
c$$$      implicit none
c$$$      integer natom, nbloch, nclass
c$$$      integer ldima(natom), iaf(natom), iof(natom),
c$$$     &  ix,ia1,ia2,ic1,ic2,ifi,i,j, iaf1,iaf2,igp,im1,im2
c$$$      integer(4):: iam(nbloch),imf(nbloch)
c$$$      real(8) :: qt(natom), rf,cf,q(3),qlat(3,3),bas(3,natom),
c$$$     &         anfvec(3),transaf(3,natom),qg(3)
c$$$      complex(8) :: zxq (nbloch,nbloch),fac(nbloch)
c$$$      complex(8):: imagtwopi ,imag=(0d0,1d0)
c$$$      complex(8),allocatable :: zxqw(:,:)
c$$$      real(8) :: pi=3.1415926535897932d0
c$$$cr
c$$$cr    True_q(1:3)     = 2*pi/alat * q(1:3)
c$$$cr  True G is given by
c$$$cr    True_G(1:3,igp) = 2*pi/alat * matmul(qlat * ngvecc(1:3,igp)) ,igp=1,ngp
c$$$c------------------------
c$$$c      imagtwopi = 2d0*(0d0,3.1415926535897932d0)
c$$$      write(6,"(' q=', 3d23.15)") q
c$$$      do ia1 = 1, natom
c$$$        transaf(1:3,ia1)= bas(1:3,ia1)+ anfvec - bas(1:3,iaf(ia1))
c$$$        write(6,"(' ia1 transaf=',i3, 3d23.15)") ia1,transaf(1:3,ia1)
c$$$      enddo
c$$$      iof(1) = 0
c$$$      do ia1 = 1, natom-1
c$$$        iof(ia1+1)= iof(ia1) + ldima(ia1)
c$$$        iam(iof(ia1)+1:iof(ia1+1)) = ia1
c$$$      enddo
c$$$      iam(iof(natom)+1:nbloch) = natom
c$$$      write(6,*) (ia1, ldima(ia1),ia1 = 1, natom)
c$$$
c$$$c (we had a problem at June18 2007)
c$$$c  this sanitiy check should be removed in fufture...
c$$$      if( nbloch /= iof(natom) +ldima(natom) ) then
c$$$        print *,' nbloch=',nbloch
c$$$        print *,' iof=',iof
c$$$        print *,' ldima=',ldima
c$$$        call rx(' anfsig: nbloch.ne.... not yet for floting orbital case')
c$$$      endif
c$$$c phase shifts
c$$$      do ia1 = 1, natom
c$$$        qt(ia1) = 2d0*pi*sum(q*transaf(1:3,ia1))
c$$$c        write( 6, "(i3,2x,d23.15)") iaf(ia1),sum(q*transaf(1:3,ia1))
c$$$c        write( 6, "(3x,2x,d23.15)") 2d0*pi*sum(q*transaf(1:3,ia1))
c$$$      enddo
c$$$      write(6,*) ' anfsig:  nbloch=',nbloch
c$$$      do im1 = 1,nbloch
c$$$        ia1  = iam(im1)
c$$$        i    = im1 - iof(ia1)
c$$$        imf (im1) = i + iof(iaf(ia1))
c$$$        fac (im1) = exp( imag*qt(ia1) )
c$$$c        write( 6, "(5i6,2d13.5)") ia1,iaf(ia1), i, im1, imf(im1),fac(im1)
c$$$      enddo
c$$$      allocate(zxqw(nbloch,nbloch))
c$$$      zxqw = zxq
c$$$c      write(6,*) ' anfsig: end1 ',sum(abs(zxqw))
c$$$      do im1 = 1,nbloch
c$$$        do im2 = 1,nbloch
c$$$          zxq(im1, im2)= 
c$$$     &               + zxqw(imf(im1),imf(im2))*fac(im1)*dconjg(fac(im2))
c$$$        enddo
c$$$      enddo
c$$$c      write(6,*) ' anfsig: end2 ',sum(abs(zxq))
c$$$      return
c$$$      end
c$$$
c sep2013t
c$$$      character(5) function charnum5(num)
c$$$      integer(4) ::num
c$$$      charnum5 = char(48+mod(num/10000,10))//
c$$$     &           char(48+mod(num/1000,10))//
c$$$     &           char(48+mod(num/100,10))//
c$$$     &           char(48+mod(num/10,10))//
c$$$     &           char(48+mod(num,10))
c$$$      end
