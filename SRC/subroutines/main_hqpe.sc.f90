!>  make sigm.* file which contains Sigma-Vxc
!     This uses amix in order to guess a better Sigma-Vxc from previous iterations.
!     * iSigma_en==5 is for diagonal-only Sigma-Vxc in LDA basis set.
!      (Then you need evec0, which contains eigenvector of LDA).
subroutine hqpe_sc()
  !------------------------------------------------------------------
  !     calculates quasiparticle energies
  !     E(k,t) = e(k,t) + Z [SEx(k,t) + SEc(k,t) - xcLDA(k,t)]
  !     e(k,t) = LDA eigenvalue
  !     Z      = [1 - dSEc(e(k,t))/dw]^(-1)
  !     SEx(k,t)   = <psi(k,t)| SEx |psi(k,t)>
  !     SEc(k,t)   = <psi(k,t)| SEc |psi(k,t)>, SEc = GWc
  !     xcLDA(k,t) = <psi(k,t)| vxc |psi(k,t)>
  !     SEx and xcLDA are in file SEX
  !     SEc is in file SEC
  !-----------------------------------------------------------------------
  use m_keyvalue,only: Getkeyvalue
  use m_read_bzdata, only: Read_bzdata, nstar, nqibz2=>nqibz
  use m_anf,only: anfcond,laf
  use m_hamindex,only: Readhamindex, nhq=>ndham
  use m_mpi,only: MPI__Initialize,MPI__Finalize
  use m_lgunit,only: m_lgunit_init
  implicit real*8 (a-h,o-z)
  implicit integer (i-n)
  logical :: sigma_mixing,lsi, lsigin
  dimension ifsex(2),ifsexcore(2),ifxc(2),ifsec(2),ifqpe(2) &
       ,iftote(2),iftote2(2),ifsex2(2),ifsexcore2(2),ifsec2(2) !sf..3June
  integer,allocatable :: itxc(:),itc(:),itx(:)
  real(8),allocatable :: qxc(:,:,:),eldaxc(:,:),vxc(:,:), &
       qc(:,:,:),eldac(:,:),sex(:,:),sexcore(:,:), &
       qx(:,:,:),eldax(:,:),rsec(:,:),csec(:,:) ,qqq(:,:,:) ,qqqx_m(:,:,:)
  complex(8), allocatable :: sex2(:,:,:),sexcore2(:,:,:), &
       sec2(:,:,:), se(:,:,:),work(:),evec_inv(:,:),evec_invt(:,:), se_ev(:,:), ev_se_ev(:,:),&
       sen(:,:),sen2(:,:),se_in(:,:), v_xc(:,:,:,:),evec(:,:,:,:),evec0(:,:),sigma_m(:,:,:,:), &
       sigin(:,:,:,:), evec00(:,:,:,:) ,evec00inv(:,:)
  real(8) qx2(3) ,qqqx(3) ,del ,mix_fac
  integer, allocatable :: ipiv(:) !sf.end
  real(8) ::  rydberg,hartree, qqqx0(3)
  integer ::  n1,n2,n3, ifse_out, ifse_in
  logical :: AddDummySig,evec0ex=.false.
  complex(8),allocatable:: sigma_m_out(:,:,:,:),pmat(:,:) &
       ,pmatd(:,:),sigmv(:,:,:)
  character(2):: soflag
  real(8) :: wex !,  eseavrs(2)
  logical :: exonly,mtosigmaonly,nexist
  integer:: ret
  character(3):: iaaa
  real(8),allocatable:: eseavr(:,:) !,eseavr_in(:,:)
  integer,allocatable:: nhqx(:,:)
  integer:: nz,ntqmin, nmto,ndimsig,ndimsig2
  integer:: procid,nrank,ifigwb_,ifigwx1_,ifigwx2_,ifvxc_,ifevec_, ipx
  character*256:: extn,ext
  character*256,allocatable:: extp(:)
  integer,allocatable:: ifevec__(:),ifvxc__(:),iprocq(:,:)
  integer,allocatable:: ntqxx(:)
  real(8):: eseavrmean,eseadd,tolq=1d-4
  integer,allocatable:: nev(:,:)
  complex(8),allocatable:: ovl(:,:)
  character(8):: xt
  call MPI__Initialize()    !this is just for exit routine subroutine rx('...') works well.
  call m_lgunit_init()
  hartree=2d0*rydberg()
  !! Shift quasiparticle energies (eV)
  !  write (*,*)' q+band index for zero?'
  !  read (*,*)jin
  jin=0
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
  open(newunit=ifsex(1) ,file='SEXU')
  open(newunit=ifxc(1)  ,file='XCU')
  open(newunit=ifsex2u, file='SEX2U',form='UNFORMATTED', status='OLD') !sf.beg
  ifsex2(1)=ifsex2u
  if( .NOT. exonly) then
     open(newunit=ifsec(1) ,file='SECU')
     open(newunit=ifsec2u, file='SEC2U',form='UNFORMATTED', status='OLD')
     ifsec2(1)=ifsec2u     !sf.end
  endif
  open(newunit=ifsexcore(1) ,file='SEXcoreU')
  open(newunit=ifSEXcore2U,file='SEXcore2U',form='UNFORMATTED',status='OLD')
  ifsexcore2(1)=ifSEXcore2U
  call readx (ifsex(1),50)
  read (ifsex(1),*) nspin,nq,ntq
  if(nspin == 2 .AND. .NOT. laf) then
     open(newunit=ifsex(2)   ,file='SEXD')
     open(newunit=ifxc(2)    ,file='XCD')
     open(newunit=ifSEX2D, file='SEX2D',form='UNFORMATTED', status='OLD') !sf.beg
     ifsex2(2)=ifSEX2D
     if( .NOT. exonly) then
        open(newunit=ifsec(2)   ,file='SECD')
        open(newunit=ifSEC2D, file='SEC2D',form='UNFORMATTED', status='OLD')
        ifsec2(2)=ifSEC2D      !sf.end
     endif
     open(newunit=ifsexcore(2)   ,file='SEXcoreD')
     open(newunit=ifSEXcore2D,file='SEXcore2D',form='UNFORMATTED',status='OLD')
     ifsexcore2(2)=ifSEXcore2D
  endif
  rewind (ifsex(1))
  write(6,*)'nq nspin=',nq,nspin
  allocate(eseavr(nq,nspin)) !,eseavr_in(nq,nspin))
  open(newunit=ifqpe(1)   ,file='QPU')
  open(newunit=iftote(1)  ,file='TOTE.UP')
  open(newunit=iftote2(1) ,file='TOTE2.UP')
  INQUIRE (FILE='sigm', EXIST = lsigin)
  open(newunit=ifsigm, file='sigm',form='UNFORMATTED') !sf
  ifse_out=ifsigm
  if (lsigin) then
     rewind ifse_out
     write(6,*) ' ... reading input sigma from file sigm'
     read(ifse_out) nspin,ndimsigin,n1,n2,n3,nqx
     if (nqx /= nq) then
        print 368, nqx,nq
368     format (6x,' (warning) file mismatch : file nq=',i4, &
             ' but expected',i4)
        lsigin = .false.
     else
        rewind ifse_out
        allocate(sigin(ndimsigin,ndimsigin,nq,nspin),qqqx_m(3,nq,nspin))
        call rwsigma('read',ifse_out,sigin,qqqx_m,nspin,ndimsigin,n1,n2,n3,nq) !,eseavr_in)
        deallocate(qqqx_m)
     endif
  endif
  if(laf) nspin=2
  rewind ifse_out
  if (nspin == 2) then
     open(newunit=ifqpe(2)   ,file='QPD')
     open(newunit=iftote(2)  ,file='TOTE.DN')
     open(newunit=iftote2(2) ,file='TOTE2.DN')
  endif
  if(jin == -101) goto 9998
  call getkeyvalue("GWinput","iSigMode",iSigma_en )
  if(isigma_en==5) then     ! .OR. core3ptest) then
     evec0ex=.false.    !true before 12Aug2006 ---> but it caused a problem maybe because of degeneracy.
     if(evec0ex) then
        open(newunit=ifevec0,file='evec0',form='UNFORMATTED',status='OLD')
        if(isigma_en==5) open(newunit=ifevecchk,file='evecfix.chk')
     endif
  endif
  call read_BZDATA()
  write(6,*)' read from bzdata nqibz2; nqibz nq nhq=',nqibz2,nq,nhq
  nsp=nspin
  nnn=nqibz2
  write(6,*)' ndimh ntq nsp nnn =',nhq,ntq,nsp,nnn !NOTE this ndimh is the maximum dimention of Hamiltonian.
  ! n PMT, ndimh is q-dependent. See lm*/gwd/sugw.F  june2009 takao
  allocate(qqq(3,nnn,nspin)) !,se_in(ntq,ntq))                 !sf.beg
  allocate(v_xc(nhq,nhq,nnn,nspin),evec(nhq,nhq,nnn,nspin),nev(nnn,nspin))
  if(evec0ex) allocate(evec00(nhq,nhq,nnn,nspin))
  allocate(nhqx(nnn,nspin))
  iqq=0
  do iq=1,nnn               !now nnn is not necessary to be nqbz !nnn=nqbz
     iqq=iqq+1
     do is=1,nspin
        open(newunit=ifvxc_,    file='vxc'//trim(xt(iq))//trim(xt(is)),form='unformatted')
        open(newunit=ifevec_,  file='evec'//trim(xt(iq))//trim(xt(is)),form='unformatted')
        read(ifvxc_ ) ndimh, nmto !nsp,nnn ,nnnx,nmto
        nz=ndimh
        read(ifevec_) !ndimhx , nspx,nnnx
        write(6,*) ' reading v_xc ... iq is nz=',iq,is,nz
        nhqx(iq,is) = nz   !nz is introduced instead of nhq
        read(ifvxc_) v_xc(1:nz,1:nz,iq,is)
        read(ifevec_) qqq(1:3,iq,is),evec(1:nz,1:nz,iq,is), nev(iq,is) !nev number of true bands nov2015
        if(evec0ex) read(ifevec0) qqqx0(1:3), evec00(1:nz,1:nz,iq,is)
        close(ifvxc_)
        close(ifevec_)
     enddo
  enddo                     !sf.end
  print *,' end of reading vxc evec'
  if(mtosigmaonly()) then
     ndimsig = nmto
  else
     ndimsig = nhq
  endif
  if (lsigin .AND. ndimsigin /= ndimsig) then
     deallocate(sigin)
     write(6,*) '... input sigma dimension mismatch ... discarding'
     lsigin = .false.
  endif
  allocate(sigmv(ndimsig,ndimsig,nq))
  open(newunit=if_eseavr ,file='ESEAVR')
9998 continue
  !! ------------------------------------------------------
  ! meanvalue of eseavr june2009.
  nstarsum= sum(nstar(:))
  write(6,*)' nstarsum=',nstarsum
  allocate(ntqxx(nq))
  do 1001  is = 1,nspin
     write(6,*) ' --- is=',is
     call readx   (ifsex(is),50)
     read (ifsex(is),*) nspinx,nqx,ntqx
     read (ifsex(is),*)
     read (ifsex(is),*) deltaw
     read (ifsex(is),*) alat
     read (ifsex(is),*) ef
     read (ifsex(is),*) esmr
     call readx   (ifxc(is),50)
     read (ifxc(is),*) nspinxc,nqxc,ntqxc
     read(ifsex2(is))  nspinx2, nqx2, ntqx2, nqbz ,  nqibz,  n1,n2,n3 !sf.beg
     if( nstarsum/= nqbz ) call rx( ' nstarsum/= nqbz')
     if( .NOT. exonly)  then
        call readx(ifsec(is),50)
        read (ifsec(is),*) nspinc,nqc,ntqc
        read(ifsec2(is))  nspinc2, nqc2, ntqc2, nqbzc2, nqibzc2,n1,n2,n3 !sf.end
     endif
     read(ifsexcore2(is)) nspinxc2,nqxc2,ntqxc2,nqbzxc2,nqibzxc2
     if (nspin /= nspinx)   call rx( 'hqpe: wrong nspin SEx')
     if (nspin /= nspinxc)  call rx( 'hqpe: wrong nspin vxc')
     if (nspin /= nspinx2)  call rx( 'hqpe: wrong nspin SEx2')
     if (nq /= nqx)         call rx( 'hqpe: wrong nq SEx')
     if (nq /= nqxc)        call rx( 'hqpe: wrong nq vxc')
     if (nq /= nqx2)        call rx( 'hqpe: wrong nq SEx2')
     if (ntq /= ntqx)       call rx( 'hqpe: wrong ntq SEx')
     if (ntq /= ntqxc)      call rx( 'hqpe: wrong ntq vxc')
     if (ntq /= ntqx2)      call rx( 'hqpe: wrong ntq SEx2')
     if( .NOT. exonly) then
        if (nqbz /=  nqbzxc2)  call rx( 'hqpe: wrong nqbzx2')
        if (nqibz /= nqibzxc2) call rx( 'hqpe: wrong nqibzxc2')
        if (nspin /= nspinxc2) call rx( 'hqpe: wrong nspin SExcore2')
        if (nspin /= nspinc)   call rx( 'hqpe: wrong nspin SEc')
        if (nspin /= nspinc2)  call rx( 'hqpe: wrong nspin SEc2')
        if (nq /= nqc)         call rx( 'hqpe: wrong nq SEc')
        if (nq /= nqc2)        call rx( 'hqpe: wrong nq SEc2')
        if (nqbz /= nqbzc2)    call rx( 'hqpe: wrong nqbzxc2')
        if (ntq /= ntqc2)      call rx( 'hqpe: wrong ntq SEc2')
        if (nq /= nqxc2)       call rx( 'hqpe: wrong nq SExcore2')
        if (ntq /= ntqxc2)     call rx( 'hqpe: wrong ntq SExcore2')
     endif
     if(is==1) write(6,*)' ###  readin XCU'
     if(is==2) write(6,*)' ###  readin XCD'
     allocate( itxc(ntq),qxc(3,ntq,nq),eldaxc(ntq,nq),vxc(ntq,nq) )
     call readx (ifxc(is),50)
     read(ifxc(is),*)
     do ip = 1,nq
        do i  = 1,ntq
           read(ifxc(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16)") &
                itxc(i),ipxx,isxxx, qxc(1:3,i,ip), eldaxc(i,ip), &
                vxc(i,ip)
        enddo
     enddo

     if(is==1) write(6,*)' ###  readin SEXU'
     if(is==2) write(6,*)' ###  readin SEXD'
     allocate( itx(ntq), qx (3,ntq,nq),eldax (ntq,nq),sex(ntq,nq) )
     allocate( sex2(ntq,ntq,nq)) !sf..3June
     call readx   (ifsex(is),50)
     read(ifsex(is),*)
     do ip = 1,nq
        do i  = 1,ntq
           read(ifsex(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16)") &
                itx(i),ipxx,isxxx, qx(1:3,i,ip), eldax(i,ip), &
                sex(i,ip)
        enddo
        read(ifsex2(is)) isx,qx2,sex2(1:ntq,1:ntq,ip) !sf..3June
     enddo

     allocate( sexcore(ntq,nq),sexcore2(ntq,ntq,nq) ) !sf..3June
     if(is==1) write(6,*)' ###  readin SEXcoreU'
     if(is==2) write(6,*)' ###  readin SEXcoreD'
     call readx   (ifsexcore(is),50)
     call readx   (ifsexcore(is),50)
     read(ifsexcore(is),*)
     do ip = 1,nq
        do i  = 1,ntq
           read(ifsexcore(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16)") &
                ixx1,ixx2,ixx3, qxxx1,qxxx2,qxxx3, exxx, sexcore(i,ip)
        enddo
        read(ifsexcore2(is)) isx,qx2,sexcore2(1:ntq,1:ntq,ip) !sf..3June
     enddo

     allocate( itc(ntq), qc (3,ntq,nq),eldac (ntq,nq) &
          ,rsec(ntq,nq),csec(ntq,nq),sec2(ntq,ntq,nq)) !sf..3June
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
              read(ifsec(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16, 3x,d24.16)") &
                   itc(i),ipxxx,isxxx, qc(1:3,i,ip), eldac(i,ip), &
                   rsec(i,ip),csec(i,ip)
           enddo
           read(ifsec2(is)) isx,qx2,sec2(1:ntq,1:ntq,ip)
        enddo
     endif
     if( .NOT. exonly) then
        if (sum(abs(itx(1:ntq)-itc(1:ntq)))/=0)  call rx('hqpe.sc itx /= itc ')
        if (sum(abs(itx(1:ntq)-itxc(1:ntq)))/=0) call rx('hqpe.sc itx /= itxc ')
     endif
     call qpe1_sc  (ifqpe(is),iftote(is),iftote2(is),itx,qx,&
          eldax,vxc,sex,sexcore, rsec,csec,jin,deltaw,alat,ef, &
          ntq,nq,is, eshift0,eshift02,eshtlda)
     write(6,*)' end of qpe1'
     close(ifqpe(is))
     close(iftote(is))
     close(iftote2(is))
     if(jin==-101) cycle
     !------------------------------------------------------------------------------
     ! making SE_ij-VXC_ij, ij are 'basis functions' indexes
     !------------------------------------------------------------------------------
     allocate(se(ntq,ntq,nq),ipiv(nhq),work(nhq*nhq) &
          ,evec_inv(nhq,nhq) ,evec_invt(nhq,nhq))
     allocate(ev_se_ev(ndimsig,ndimsig))
     if(evec0ex) allocate(evec00inv(nhq,nhq))
     do ip=1,nq
        do itp=1,ntq
           do itpp=1,ntq    !make Sigma hermitean
              if( .NOT. exonly) then
                 se(itpp,itp,ip)=sex2(itpp,itp,ip)+sexcore2(itpp,itp,ip) &
                      +.5d0*(sec2(itpp,itp,ip)+dconjg(sec2(itp,itpp,ip)))
              else
                 se(itpp,itp,ip)= wex*sex2(itpp,itp,ip)+sexcore2(itpp,itp,ip)
              endif
           enddo
        enddo
     enddo
     do ip=1,nq
        do itp=1,ntq
           do itpp=1,ntq    !make Sigma hermitian
              if(abs(se(itpp,itp,ip)-dconjg(se(itp,itpp,ip)) ) > tolq) then
                 write(*,*)'diff=', se(itpp,itp,ip)-dconjg(se(itp,itpp,ip))
                 write(*,*)'se   ', se(itpp,itp,ip), se(itp,itpp,ip)
                 write(*,*)'sex2 ',sex2(itpp,itp,ip),sex2(itp,itpp,ip)
                 write(*,*)'sexc2',sexcore2(itpp,itp,ip),sexcore2(itp,itpp,ip)
                 call rx( "hqpe: Sigma_nn' is not hermitian")
              endif
           enddo
        enddo
     enddo
     if (is==1 .AND. laf) then
        write(ifse_out) 1,ndimsig,n1,n2,n3,nq,0,0,0
     elseif (is == 1) then
        write(ifse_out) nspin,ndimsig,n1,n2,n3,nq,0,0,0
     endif
     do 2001 ip=1,nq
        do ikp=1,nnn        !nqbz
           if (sum ( (qqq(1:3,ikp,is)-qx(1:3,1,ip))**2 ) < tolq ) then
              ikpx=ikp      !qc(:,i,:) does not depents on band index i=1:ntq
              goto 100
           endif
        enddo               !ikp
        call rx( 'hqpe.sc: not find ikp 100')
100     continue
        nz = nhqx(ikpx,is)
        if(evec0ex .AND. iSigma_en==5) then
           call rx('Not support evec0ex.and.iSigma_en==5 now... sep2013')
        endif
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
              if(itp/=itpp .AND. isigma_en==5) then
                 se(itp,itpp,ip)= 0d0
              else
                 if( .NOT. exonly) then
                    se(itp,itpp,ip)= se(itp,itpp,ip) &
                         -.5d0* sum(dconjg(evec(1:nz,itp,ikpx,is))* &
                         matmul(v_xc(1:nz,1:nz,ikpx,is),evec(1:nz,itpp,ikpx,is)))
                 else
                    se(itp,itpp,ip)= se(itp,itpp,ip) &
                         -wex*.5d0* sum(dconjg(evec(1:nz,itp,ikpx,is))* &
                         matmul(v_xc(1:nz,1:nz,ikpx,is),evec(1:nz,itpp,ikpx,is)))
                 endif
                 !!
              endif
           enddo
        enddo
2001 enddo
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
        if(iix==0) then
           eavr2=1d0
           eseavr02=0d0
        endif
        eseavr(ip,is) = eseavr02/eavr2 !now eseavr is
        eseavr(ip,is) = 2d0*eseavr(ip,is) !in Ry.
        write(6,*)"### A correction takao2009June: find this in hqpe.se.m.F"
        write(6,*)"###   constant is added to sigm above threshold."
        write(6,*)"###   the constant (ESEAVR=e-weighted average Ry)= ",is,ip,eseavr(ip,is)
2002 enddo
     eseavrmean=0d0
     do ip=1,nq
        eseavrmean = eseavrmean + nstar(ip)*eseavr(ip,is)
     enddo
     eseavrmean = eseavrmean/nqbz
     call getkeyvalue("GWinput","AddToESEAVR",eseadd,default=0d0,status=ret)
     eseavrmean=  eseavrmean+eseadd
     write(if_eseavr,"(d23.15,i3,i8)") eseavrmean,is ,ntqxxmin
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
102     continue
        nz = nhqx(ikpx,is) !june 2009 takao. we use nz instead of nhq
!!!  Make inverse evec_inv(n,i) matrix \psi_n=sum_i evec(i,n)\phi_i,
!!!  where \psi is eigenfunction and \phi is basis function
        !! nov2015
        !! evec_inv(ib1,iww)= \sum_ib2  ovlinv(ib1,ib2)*dconjg(evec(iww,ib2))  nov2015, we introduce nev. iww is for PMT basis. ib for band index.
        !! This is for converting rotated evec (=evecrot(ib)) in the representation of original evec(ib).
        nevv=nev(ikpx,is)
        allocate(ovl(nevv,nevv))
        do i=1,nevv
           do j=1,nevv
              ovl(i,j)=sum(dconjg(evec(1:nz,i,ikpx,is))*evec(1:nz,j,ikpx,is))
           enddo
        enddo
        call matcinv(nevv,ovl) !ovl --> ovlinv
        evec_inv(1:nevv,1:nz) = matmul(ovl(1:nevv,1:nevv),dconjg(transpose(evec(1:nz,1:nevv,ikpx,is)))) !note ovl means ovlinv
        deallocate(ovl)
        evec_invt(1:nz,1:nevv) = transpose(dconjg(evec_inv(1:nevv,1:nz)))
        if(mtosigmaonly()) then
           ndimsig2= nmto
        else
           ndimsig2 = nz
        endif
        ev_se_ev(1:ndimsig2,1:ndimsig2) = matmul(evec_invt(1:ndimsig2,1:ntqxx(ip)) &
             ,matmul(se(1:ntqxx(ip),1:ntqxx(ip),ip),evec_inv(1:ntqxx(ip),1:ndimsig2)))
        !! sep2013 exprapolation of se by eseavrmean
        do itp =1,ndimsig2
           do itpp=1,ndimsig2
              ev_se_ev(itp,itpp) = ev_se_ev(itp,itpp) + &
                   sum(evec_invt(itp,ntqxx(ip)+1:nevv)*evec_inv(ntqxx(ip)+1:nevv,itpp))*eseavrmean/2d0 ! in Hartree.
           enddo
        enddo
        !! - write SE_ij-Vxc_ij where ij are basis function indices
        do itp=1,ndimsig2
           do itpp=1,ndimsig2
              if( abs(ev_se_ev(itpp,itp) - dconjg(ev_se_ev(itp,itpp)) )> tolq  ) then
                 write(6,*)itp,itpp
                 write(6,*)ev_se_ev(itpp,itp)
                 write(6,*)ev_se_ev(itp,itpp)
                 call rx( 'hqpe: Sigma_ij is not hermitian')
              endif
              if(abs(v_xc(itp,itpp,ikpx,is)-dconjg(v_xc(itpp,itp,ikpx,is)))> tolq) &
                   call rx( 'hqpe: v_xc is not hermitean')
           enddo
        enddo
        write(ifse_out) qqq(1:3,ikpx,is),is !,eseavr(ip,is) !,ntqxxmin ! ntqxx(ip)  sep2013
        sigmv(:,:,ip) = 1d20
        sigmv(1:ndimsig2,1:ndimsig2,ip) =  2d0*ev_se_ev(1:ndimsig2,1:ndimsig2) !in Ry.
        write(6,*) 'ssssss1=', sum(abs(sigmv(:,:,ip)))
        write(ifse_out) sigmv(:,:,ip)
        !!      2*ev_se_ev bacause v_xc in sugw.f was in rydberg while SE was in hartree
2003 enddo
     write(6,*)
     deallocate(sex2,sexcore2,sec2,se,ipiv,work)
     deallocate(evec_inv,evec_invt,ev_se_ev) !,se_ev
     if(evec0ex) deallocate(evec00inv)
     ! - end making SE_ij-VXC_ij cccccccccccccccccccccccccccccccccc  !sf..3June
     deallocate( itxc,qxc,eldaxc,vxc ,itc, qc ,eldac, sexcore ,rsec,csec, itx, qx ,eldax,sex)
     if (jin > 0) jin = 999999
     if (laf) exit
1001 enddo
  if(jin==-101) goto 9999
  deallocate(v_xc,evec)
  close(ifse_out)
  ! mixing sigma ----     ...
  open(newUNIT=ifse_out, file='sigm',form='UNFORMATTED') !Once readin sigma
  allocate(sigma_m(ndimsig,ndimsig,nq,nspin),qqqx_m(3,nq,nspin))
  write(6,*)"========= Sigma mixing section using mixsigma ======="
  call rwsigma ('read',ifse_out,sigma_m,qqqx_m, nspin,ndimsig,n1,n2,n3,nq ) !,eseavr_in)
  if (lsigin .AND. ndimsigin /= ndimsig) then
     deallocate(sigin)
     write(6,*) '... input sigma dimension mismatch ... discarding'
     lsigin = .false.
  endif
  if ( .NOT. lsigin) then
     allocate(sigin(1,1,1,1))
     sigin = 0d0
  endif
  call mixsigma(sigma_m, lsigin, sigin, ndimsig**2*nq*nspin)
  write(6,*)
  write(6,*) "=== Write sigm to files ========"
  rewind ifse_out
  call rwsigma ('write',ifse_out,sigma_m,qqqx_m, nspin,ndimsig,n1,n2,n3,nq) !,eseavr)
  close(ifse_out)
9999 continue
  call rx0s( ' OK! hqpe_sc ')
end subroutine hqpe_sc

!------------------------------------------------------
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
end subroutine testfff

!------------------------------------------------------
subroutine onlydiag(sigm,ndimh)
  implicit integer (i-n)
  complex(8):: sigm(ndimh,ndimh)
  do ix=1,ndimh
     do iy=1,ndimh
        if(ix/=iy) sigm(ix,iy)=0d0
     enddo
  enddo
end subroutine onlydiag

!----------------------------------------------------------------------
subroutine rwsigma(rw,ifss,sigma_m,qqqx_m,nspin,ntq,n1,n2,n3,nq ) !,eseavr)
  implicit integer (i-n)
  complex(8)::sigma_m(ntq,ntq,nq,nspin)
  character(*):: rw
  real(8)::qqqx_m(3,nq,nspin) !,eseavr(nq,nspin)
  ntm1=0;ntm2=0;ntm3=0
  if(rw=='read')  ifs= 1
  if(rw=='write') ifs=-1
  if(ifs>0) write(6,*) " rwsigma: Reading sigm (Binary)==="
  if(ifs<0) write(6,*) " rwsigma  Writing sigm (Binary)==="
  if(ifs>0) read ( ifss) nspin,ntq,n1,n2,n3,nq, ntm1,ntm2,ntm3
  if(ifs<0) write( ifss) nspin,ntq,n1,n2,n3,nq, ntm1,ntm2,ntm3
  do is=1,nspin
     do ip=1,nq
        if(ifs>0) read (ifss)  qqqx_m(1:3,ip,is)!,isr !,eseavr(ip,is)
        if(ifs>0) read (ifss)  sigma_m(1:ntq,1:ntq,ip,is)
        if(ifs<0) write(ifss) qqqx_m(1:3,ip,is),is !,eseavr(ip,is)
        if(ifs<0) write(ifss) sigma_m(1:ntq,1:ntq,ip,is)
        write(6,"('  === ',i5,i3,3f10.5)")ip,is,qqqx_m(1:3,ip,is)
     enddo
  enddo
  write(6,*)' === rwsigma:  sum check of sigma_m=',sum(abs(sigma_m))
  print *
end subroutine rwsigma

!----------------------------------------------------------------------
subroutine mixsigma(sss, lsigin, sigin, nda)
  use m_amix,only: amix
  !  subroutine pqmixa(nda,nmix,mmix,mxsav,beta,rms2,a,tj)
  !- Mixing routine for sigma. Modified from pqmixa in subs/pqmix.f
  !- Anderson mixing of a vector
  !i  mmix: number of iterates available to mix
  ! o nmix: nmix > 0: number of iter to try and mix
  !i        nmix < 0: use mmix instead of nmix.
  !o  nmix: (abs)  number of iter actually mixed.
  !o        (sign) <0, intended that caller update nmix for next call.
  !  MvS Feb 04 use sigin as input sigma if available (lsigin=T)
  !             Add mixnit as parameter
  use m_keyvalue,only: getkeyvalue
  implicit none
  logical :: lsigin
  integer :: nda,nmix,mmix
  integer,parameter:: mxsav=10
  double precision :: rms2,tj(mxsav),beta
  integer :: im,imix,jmix,onorm,okpvt,oa !,amix
  integer :: iprintxx,ifi,nitr,ndaf
  complex(8)::sss(nda),sigin(nda),img=(0d0,1d0)
  real(8):: tjmax
  real(8),allocatable::norm(:),a(:,:,:)
  integer,allocatable:: kpvt(:)
  integer::ret
  character(8) :: fff
  logical :: fexist
  real(8):: acc
  integer:: ido
  iprintxx = 30
  beta=1d0
  call getkeyvalue("GWinput","mixbeta",beta,default=1d0,status=ret)
  write(6,*)'('' mixsigma: Anderson mixing sigma with mixing beta ='',f12.6)',beta
  !' ... reads prior iteration INCLUDING starting sigma for current iteration
  allocate ( a(2*nda,0:mxsav+1,2) )
  fff="mixsigma"
  INQUIRE (FILE =fff, EXIST = fexist)
  if(fexist)      write(6,*)'... reading file mixsigma'
  if( .NOT. fexist) write(6,*)'... No file mixsigma'
  open(newunit=ifi,file=fff,form='unformatted')
  if(fexist) then
     read(ifi,err=903,end=903) nitr,ndaf
     if (ndaf /= nda) goto 903
     read(ifi,err=903,end=903) a
     goto 902
  endif
  goto 901
903 continue
  write(6,"(5x,'(warning) file mismatch ... mixing file not read')")
901 continue
  nitr = 0
902 continue
  a(1:nda,0,1)       = dreal(sss)      !output
  a(nda+1:2*nda,0,1) = dimag(sss)      !output
  !     if input sigma available, use it instead of file a(:,0,2)
  if (lsigin) then
     write(6,*)'... using input sigma read from sigm file'
     a(1:nda,0,2)       = dreal(sigin)        !input
     a(nda+1:2*nda,0,2) = dimag(sigin)  !input
  endif
  !     Restrict maximum number of prior iterations
  call getkeyvalue("GWinput","mixpriorit",imix,default=9,status=ret)
  mmix = min(max(nitr-1,0),imix)
  if (mmix > mxsav) mmix = mxsav
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
  imix = amix(nda*2,mmix,mxsav,ido,dabs(beta),iprintxx,tjmax, a,tj,rms2) !norm,kpvt,
  sss= a(1:nda,0,2) + img* a(nda+1:2*nda,0,2) 
  rewind(ifi)
  write(ifi) nitr+1,nda
  write(ifi) a
  close(ifi)
end subroutine mixsigma
