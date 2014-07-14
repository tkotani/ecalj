      program hqpe_mlw

c Oct,2008 Takashi Miyake, from hqpe.m.F
c calculates quasiparticle energies
c E(k,t) = e(k,t) + Z [SEx(k,t) + SEc(k,t) - xcLDA(k,t)]
c e(k,t) = LDA eigenvalue
c Z      = [1 - dSEc(e(k,t))/dw]^(-1)
c SEx(k,t)   = <psi(k,t)| SEx |psi(k,t)>
c SEc(k,t)   = <psi(k,t)| SEc |psi(k,t)>, SEc = GWc
c xcLDA(k,t) = <psi(k,t)| vxc |psi(k,t)>
c SEx and xcLDA are in file SEX
c SEc is in file SEC
      use m_readeigen
      use m_read_bzdata,ngrp2=>ngrp
      use m_genallcf_v3
      use keyvalue
      implicit real*8 (a-h,o-z)

c local data
      logical laf
      dimension ifsex(2),ifsexcore(2),ifxc(2),ifsec(2),ifqpe(2)
     & ,iftote(2),iftote2(2)
      integer(4),allocatable :: itxc(:),itc(:),itx(:)
      real(8),allocatable :: qxc(:,:,:),eldaxc(:,:),vxc(:,:),
     &    qc(:,:,:),eldac(:,:),sex(:,:),sexcore(:,:),
     &    qx(:,:,:),eldax(:,:),rsec(:,:,:),csec(:,:,:),zfac(:,:)

      integer(4):: ret
c      logical:: readgwinput
ctm, 070806 scissors operator
      logical(4):: lscissors
      integer(4):: ifsci,ib1,ib2
      real(8):: deig
ctm, 081021
      integer(4):: iko_ix,iko_fx,nwf,iwini,iwfin,iw,mrece
      integer(4),allocatable :: itq(:)
      real(8),allocatable :: omg(:),qpt(:,:),elda(:,:)
      real(8):: q(3)
      complex(8),allocatable:: umn_bz(:,:,:),umn(:,:,:)

c shift quasiparticle energies (eV)
!      write (*,*)' q+band index for zero?'
!      read (*,*)jin
!      call headver('hqpe',jin)

ctm, 081029
      call read_BZDATA()
      ifhbed = 99
      open(ifhbed,file='hbe.d',status='old')
      read (ifhbed,*) nprecb,mrecb,mrece,nlmtot,nqbzt, nband,mrecg
      close(ifhbed)
      call genallcf_v3(-999,-999d0,-1) !in module m_genallcf_v3
      call init_readeigen(ginv,nspin,nband,mrece) !initialization of readEigen 
      call init_readeigen2(mrecb,nlmto,mrecg) !initialize m_readeigen     

c open files
c> input files
c      ifqpnt     = iopenxx ('QPNT')
c      if(readgwinput()) then
       call getkeyvalue("GWinput","<QPNT>",unit=ifqpnt,status=ret) 
c      else
c        ifqpnt     = iopenxx ('QPNT')
c      endif
      laf        = .false.

      call readx   (ifqpnt,10)
      read (ifqpnt,*) iqall,iaf
      if (iaf .eq. 1) laf = .true.
      call readx   (ifqpnt,10)
      read (ifqpnt,*) ntq
      allocate( itq(ntq) )
      read (ifqpnt,*) (itq(i),i=1,ntq)

      ifsex(1)   = iopenxx ('SEXU.wan')
      ifsexcore(1) = iopenxx ('SEXcoreU.wan')
      ifsec(1)   = iopenxx ('SEComg.UP.wan')
      ifxc(1)    = iopenxx ('XCU.wan')


      call readx   (ifsex(1),50)
      read (ifsex(1),*) nspin2,nq,ntw
      if (nspin2.ne.nspin) stop 'hqpe_mlw: nspin error'
      if (nspin .eq. 2 .and. .not. laf) then
      ifsex(2)   = iopenxx ('SEXD.wan')
      ifsexcore(2)   = iopenxx ('SEXcoreD.wan')
      ifsec(2)   = iopenxx ('SEComg.DN.wan')
      ifxc(2)    = iopenxx ('XCD.wan')
      endif
      rewind (ifsex(1))

c> output file
      ifqpe(1)   = iopenxx ('QPU.wan')
      iftote(1)  = iopenxx ('TOTE.UP.wan')
      iftote2(1) = iopenxx ('TOTE2.UP.wan')
      if (nspin == 2) then
        ifqpe(2)   = iopenxx ('QPD.wan')
        iftote(2)  = iopenxx ('TOTE.DN.wan')
        iftote2(2) = iopenxx ('TOTE2.DN.wan')
      endif


c loop over spin
      do      is = 1,nspin
      print *, ' --- is=',is
c read dimensions
      call readx   (ifsex(is),50)
      read (ifsex(is),*) nspinx,nqx,ntwx
      read (ifsex(is),*)
      read (ifsex(is),*) deltaw
      read (ifsex(is),*) alat
      read (ifsex(is),*) ef
c
      call readx(ifsec(is),50)
      read (ifsec(is),*) iwini,iwfin,ntwc,nqc
c
      call readx   (ifxc(is),50)
      read (ifxc(is),*) nspinxc,nqxc,ntwxc
      if (nspin .ne. nspinx)  stop 'hqpe: wrong nspin SEx'
!      if (nspin .ne. nspinc)  stop 'hqpe: wrong nspin SEc'
      if (nspin .ne. nspinxc) stop 'hqpe: wrong nspin vxc'
      if (nq .ne. nqx)        stop 'hqpe: wrong nq SEx'
      if (nq .ne. nqc)        stop 'hqpe: wrong nq SEc'
      if (nq .ne. nqxc)       stop 'hqpe: wrong nq vxc'
      if (ntw .ne. ntwx)      stop 'hqpe: wrong ntw SEx'
      if (ntw .ne. ntwc)      stop 'hqpe: wrong ntw SEc'
      if (ntw.ne. ntwxc)      stop 'hqpe: wrong ntw vxc'
c
      if(is==1) print *,' ###  readin XCU.wan'
      if(is==2) print *,' ###  readin XCD.wan'
      allocate( itxc(ntw),qxc(3,ntw,nq),eldaxc(ntw,nq),
     &          vxc(ntw,nq) )
      call readx (ifxc(is),50)
	read(ifxc(is),*)
      do ip = 1,nq
      do i  = 1,ntw
        read(ifxc(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16)")
     &      itxc(i),ipxx,isxxx, qxc(1:3,i,ip), eldaxc(i,ip), 
     &      vxc(i,ip)
      enddo
      enddo
c
      if(is==1) print *,' ###  readin SEXU.wan'
      if(is==2) print *,' ###  readin SEXD.wan'
      allocate( itx(ntw), qx (3,ntw,nq),eldax (ntw,nq),sex(ntw,nq) )
      call readx   (ifsex(is),50)
	read(ifsex(is),*)
      do ip = 1,nq
      do i  = 1,ntw
        read(ifsex(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16)")
     &      itx(i),ipxx,isxxx, qx(1:3,i,ip), eldax(i,ip), 
     &      sex(i,ip)
      enddo
      enddo
c
      if(is==1) print *,' ###  readin SEXcoreU.wan'
      if(is==2) print *,' ###  readin SEXcoreD.wan'
      allocate( sexcore(ntw,nq) )
      call readx   (ifsexcore(is),50)
      call readx   (ifsexcore(is),50)
	read(ifsexcore(is),*)
      do ip = 1,nq
      do i  = 1,ntw
        read(ifsexcore(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16)")
     &      ixx1,ixx2,ixx3, qxxx1,qxxx2,qxxx3, exxx, sexcore(i,ip)
      enddo
      enddo
c
      if(is==1) print *,' ###  readin SEComg.UP.wan'
      if(is==2) print *,' ###  readin SEComg.DN.wan'
      allocate( itc(ntw), qc (3,ntw,nq)
     &         ,rsec(iwini:iwfin,ntw,nq),csec(iwini:iwfin,ntw,nq)
     &         ,omg(iwin:iwfin))
      do ip = 1,nq
      do i  = 1,ntw
      do iw = iwini,iwfin
        read(ifsec(is),"(4i5,3f10.6,3x,f10.6,2x,f16.8,x,3f16.8)")
     &     iw2,itc(i),ip2,is2,qc(1:3,i,ip),exxx,
     &     omg(iw),rsec(iw,i,ip),csec(iw,i,ip) 
c      print *," itc=",itc(i)
      enddo
      enddo
      enddo
      
c check that LDA eigenvalues are the same
      itxcxxx     = icompvv2 (itx, itc,ntq,' state-label x c')
      itxcxxx     = icompvv2 (itx,itxc,ntq,'state-label x xc')

ctm, 070806 scissors operator
      inquire(file='SCISSORS',EXIST=lscissors)
      if (lscissors) then
        ifsci = 102
        open (ifsci,file='SCISSORS')
        read(ifsci,*)
        read(ifsci,*)ib1,ib2,deig
        if (ib2.gt.ntw) ib2=ntw
        close(ifsci)
        rsec(:,ib1:ib2,:) = rsec(:,ib1:ib2,:) - deig
      endif  

ctm, 081027
c read Umn
      ifi = 99
      if (is==1) then
        open(ifi,file='MLWU',form='unformatted',status='old')
      else
        open(ifi,file='MLWD',form='unformatted',status='old')
      endif ! is   
      read(ifi)nqbz2,nwf,iko_ix,iko_fx
      if (nqbz2 .ne. nqbz) stop 'hqpe_mlw: nqbz error'
      allocate(umn_bz(iko_ix:iko_fx,nwf,nqbz))
      do iqbz = 1,nqbz
          read(ifmlw)iqbz2,q
          if (iqbz2 .ne. iqbz) stop 'hqpe_mlw: iqbz error'
          read(ifi)umn_bz(iko_ix:iko_fx,1:nwf,iqbz)
      enddo
      close(ifi)
c Umn(qbz) to Umn(q)
      allocate(umn(iko_ix:iko_fx,ntw,nq),qpt(3,nq) )
      do ip  = 1,nq
         qpt(:,ip) = qx(:,1,ip)
         iq = iqindx(qpt(:,ip),ginv,qbz,nqbz)
         do i = 1,ntw
           umn(:,i,ip) = umn_bz(:,itx(i),iq)
         enddo
      enddo ! ip
c read eks,nbks
      allocate(elda(nband,nq))
      do ip = 1,nq
         call readeval(qpt(:,ip),isp,elda(:,ip))
      enddo   

c calculate quasiparticle energies and write them on file QPE
      call qpe1_mlw (ifqpe(is),iftote(is),iftote2(is),itq,qpt,
     i              elda,vxc,sex,sexcore,
     i              rsec,csec,jin,deltaw,alat,ef,omg,
     i              umn,iko_ix,iko_fx,iwini,iwfin,
     d              nband,ntq,ntw,nq,is,
     o              eshift0,eshift02,eshtlda)

      deallocate( itxc,qxc,eldaxc,vxc ,itc, qc ,
     &                 sexcore ,rsec,csec,omg,
     &       itx, qx ,eldax,sex) 
      deallocate(umn_bz,umn,qpt,elda)
      if (laf) exit
!      if (jin .gt. 0) jin = 999999
      end do ! is
      stop ' OK! hqpe '  
      end
