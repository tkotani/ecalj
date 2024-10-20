!>  make sigm.* file which contains Sigma-Vxc
!     This uses amix in order to guess a better Sigma-Vxc from previous iterations.
!     * iSigma_en==5 is for diagonal-only Sigma-Vxc in LDA basis set.
!      (Then you need evec0, which contains eigenvector of LDA).
module m_hqpe_sc
  public hqpe_sc
  private
contains
  subroutine hqpe_sc() bind(C)     ! Summarize  quasiparticle energies
    !     E(k,t) = e(k,t) + Z [SEx(k,t) + SEc(k,t) - xcLDA(k,t)]
    !     e(k,t) = LDA eigenvalue
    !     Z      = [1 - dSEc(e(k,t))/dw]^(-1)
    !     SEx(k,t)   = <psi(k,t)| SEx |psi(k,t)>
    !     SEc(k,t)   = <psi(k,t)| SEc |psi(k,t)>, SEc = GWc
    !     xcLDA(k,t) = <psi(k,t)| vxc |psi(k,t)>
    !     SEx and xcLDA are in file SEX
    !     SEc is in file SEC
    use m_keyvalue,only: Getkeyvalue
    use m_read_bzdata, only: Read_bzdata, nstar, nqibz2=>nqibz, nqbz,nqibz,  n1,n2,n3
    use m_hamindex,only: Readhamindex, nhq=>ndham
    use m_mpi,only: MPI__Initialize, mpi__rank
    use m_lgunit,only: m_lgunit_init,stdo
    use m_genallcf_v3,only: genallcf_v3,laf,nmto=>nlmto 
    use m_ftox
    implicit real*8 (a-h,o-z)
    implicit integer (i-n)
    logical :: sigma_mixing,lsi, lsigin
    dimension ifsex(2),ifsexcore(2),ifxc(2),ifsec(2),ifqpe(2) &
         ,iftote(2),iftote2(2),ifsex2(2),ifsexcore2(2),ifsec2(2)
    integer,allocatable :: itxc(:),itx(:) !,itc(:)
    real(8),allocatable :: qxc(:,:,:),eldaxc(:,:),vxc(:,:), &
         qc(:,:,:),eldac(:,:),sex(:,:),sexcore(:,:), &
         qx(:,:,:),eldax(:,:),rsec(:,:),csec(:,:) ,qqq(:,:,:) ,qqqx_m(:,:,:)
    complex(8), allocatable :: sex2(:,:,:),sexcore2(:,:,:), &
         sec2(:,:,:), se(:,:,:),work(:),evec_inv(:,:),evec_invt(:,:), se_ev(:,:), ev_se_ev(:,:),&
         sen(:,:),sen2(:,:),se_in(:,:), v_xc(:,:,:,:),evec(:,:,:,:),evec0(:,:),sigma_m(:,:,:,:), &
         sigin(:,:,:,:), evec00(:,:,:,:) ,evec00inv(:,:)
    real(8) qx2(3) ,qqqx(3) ,del ,mix_fac
    integer, allocatable :: ipiv(:)
    real(8) ::  rydberg,hartree, qqqx0(3)
    integer :: ifse_out, ifse_in
    logical :: AddDummySig,evec0ex=.false.
    complex(8),allocatable:: sigma_m_out(:,:,:,:),pmat(:,:),pmatd(:,:),sigmv(:,:,:)
    character(2):: soflag
    real(8) :: wex 
    logical :: mtosigmaonly,nexist !exonly,
    integer:: ret
    character(3):: iaaa
    real(8),allocatable:: eseavr(:,:) !,eseavr_in(:,:)
    integer,allocatable:: nhqx(:,:)
    integer:: nz,ntqmin, ndimsig,ndimsig2
    integer:: procid,nrank,ifigwb_,ifigwx1_,ifigwx2_,ifvxcevec,ifevec_, ipx
    character*256:: extn,ext
    character*256,allocatable:: extp(:)
    integer,allocatable:: ifevec__(:),ifvxcevec_(:),iprocq(:,:)
    integer,allocatable:: ntqxx(:)
    real(8):: eseavrmean,eseadd,tolq=1d-4
    integer,allocatable:: nev(:,:)
    complex(8),allocatable:: ovl(:,:)
    character(8):: xt
    logical:: nsp2laf
    ! call getkeyvalue("GWinput","EXonly",wex,default=0d0,status=ret); exonly = .not.(wex==0d0);  if(exonly) write(6,*)' exonly=T wex=',wex
    call MPI__Initialize()    !this is just for exit routine subroutine rx('...') works well.
    call m_lgunit_init()
    call genallcf_v3(0) 
    call readhamindex()
    call read_BZDATA()
    hartree=2d0*rydberg()
! open files    
    open(newunit=ifsex(1) , file='SEXU')
    open(newunit=ifxc(1)  , file='XCU')
    open(newunit=ifsex2(1), file='SEX2U',form='UNFORMATTED', status='OLD') !    open(newunit=ifsec(1),file='SECU')
    open(newunit=ifsec2(1), file='SEC2U',form='UNFORMATTED', status='OLD') !    open(newunit=ifsexcore(1) ,file='SEXcoreU')
    open(newunit=ifsexcore2(1),file='SEXcore2U',form='UNFORMATTED',status='OLD')
    call readx (ifsex(1),50)
    read (ifsex(1),*) nspin,nq,ntq
    rewind (ifsex(1))
    write(stdo,*)'nq nspin=',nq,nspin
    nsp2laf= nspin == 2 .AND. .NOT. laf
    if(nsp2laf)  open(newunit=ifsex(2),     file='SEXD')
    if(nsp2laf)  open(newunit=ifxc(2) ,     file='XCD')
    if(nsp2laf)  open(newunit=ifsex2(2),    file='SEX2D',form='UNFORMATTED', status='OLD')
    if(nsp2laf)  open(newunit=ifsec2(2),    file='SEC2D',form='UNFORMATTED', status='OLD')
    if(nsp2laf)  open(newunit=ifsexcore2(2),file='SEXcore2D',form='UNFORMATTED',status='OLD')
    open(newunit=ifqpe(1)   ,file='QPU')
    open(newunit=iftote(1)  ,file='TOTE.UP')
    open(newunit=iftote2(1) ,file='TOTE2.UP')
    INQUIRE (FILE='sigm', EXIST = lsigin)
    open(newunit=ifsigm, file='sigm',form='UNFORMATTED')
    if(nspin==2) open(newunit=ifqpe(2)   ,file='QPD')
    if(nspin==2) open(newunit=iftote(2)  ,file='TOTE.DN')
    if(nspin==2) open(newunit=iftote2(2) ,file='TOTE2.DN')
!    
    allocate(eseavr(nq,nspin))
    if (lsigin) then !   write(6,*) ' ... reading input sigma from file sigm'
      read(ifsigm) nspin,ndimsigin !,n1,n2,n3,nqx
      if (nqx /= nq) then
        write(stdo,"(' (warning) file mismatch : file nq=',i4,' but expected',i4)") nqx,nq
        lsigin = .false.
      else
        rewind ifsigm
        allocate(sigin(ndimsigin,ndimsigin,nq,nspin),qqqx_m(3,nq,nspin))
        call rwsigma('read',ifsigm,sigin,qqqx_m,nspin,ndimsigin,n1,n2,n3,nq) 
        deallocate(qqqx_m)
      endif
    endif
    if(laf) nspin=2 !right?
    rewind ifsigm
    write(6,*)' read from bzdata nqibz2; nqibz nq nhq=',nqibz2,nq,nhq
    nsp=nspin
    nnn=nqibz2
    write(6,*)' ndimh ntq nsp nnn =',nhq,ntq,nsp,nnn !NOTE this ndimh is the maximum dimention of Hamiltonian.
    allocate(qqq(3,nnn,nspin),v_xc(nhq,nhq,nnn,nspin),evec(nhq,nhq,nnn,nspin),nev(nnn,nspin))
    allocate(nhqx(nnn,nspin))
    iqq=0
    do iq=1,nnn               !now nnn is not necessary to be nqbz !nnn=nqbz
      iqq=iqq+1
      do is=1,nspin
        open(newunit=ifvxcevec, file='vxcevec'//trim(xt(iq))//trim(xt(is)),form='unformatted')
        read(ifvxcevec) qqq(1:3,iq,is),nz, nev(iq,is)
        read(ifvxcevec) v_xc(1:nz,1:nz,iq,is)
        read(ifvxcevec) evec(1:nz,1:nz,iq,is)  !nev number of true bands nov2015
        close(ifvxcevec)
        nhqx(iq,is) = nz   !nz is introduced instead of nhq
        write(6,*) ' reading vxcevec ... iq is nz=',iq,is,nz
      enddo
    enddo
    ndimsig= merge(nmto,nhq, mtosigmaonly())
    if(lsigin .AND. ndimsigin /= ndimsig) then
      deallocate(sigin)
      write(6,*) '... input sigma dimension mismatch ... discarding'
      lsigin = .false.
    endif
    allocate(sigmv(ndimsig,ndimsig,nq))
    nstarsum= sum(nstar(:)) !    open(newunit=if_eseavr ,file='ESEAVR')
    if(nstarsum/= nqbz ) call rx( ' nstarsum/= nqbz')
    allocate(ntqxx(nq))
    isloop: do 1001  is = 1,nspin 
      write(6,*) ' --- is=',is
      call readx(ifsex(is),50)
      read (ifsex(is),*) nspinx,nqx,ntqx
      read (ifsex(is),*)
      read (ifsex(is),*) !deltaw
      read (ifsex(is),*) !alat
      read (ifsex(is),*) ef
      read (ifsex(is),*) !esmr
      
      read(ifsex2(is))  nspinx2, nqx2, ntqx2 !, nqbz ,  nqibz,  n1,n2,n3
      read(ifsec2(is))   nspinc,nqc,ntqc, nqbzc2, nqibzc2!,n1,n2,n3
      read(ifsexcore2(is)) nspinxc2,nqxc2,ntqxc2,nqbzxc2,nqibzxc2
      
      allocate( itxc(ntq),qxc(3,ntq,nq),eldaxc(ntq,nq),vxc(ntq,nq) )
      allocate( itx(ntq), qx (3,ntq,nq),eldax (ntq,nq),sex(ntq,nq) )
      call readx   (ifxc(is),50)
      read (ifxc(is),*) nspinxc,nqxc,ntqxc
      call readx (ifxc(is),50)
      read(ifxc(is),*)
      do ip = 1,nq
        do i  = 1,ntq
          read(ifxc(is),*) itx(i),ipxx,isxxx, qx(1:3,i,ip), eldax(i,ip),vxc(i,ip)
        enddo
      enddo
      allocate( sex2(ntq,ntq,nq))
!      call readx   (ifsex(is),50)
!      read(ifsex(is),*)
      do ip = 1,nq
        read(ifsex2(is)) isx,qx2,sex2(1:ntq,1:ntq,ip)
        do i  = 1,ntq !        read(ifsex(is),*) itx(i) !,ipxx,isxxx, qx(1:3,i,ip)
          sex(i,ip)=sex2(i,i,ip)*hartree !sex in eV
        enddo
      enddo
      allocate( sexcore(ntq,nq),sexcore2(ntq,ntq,nq) )
!      call readx   (ifsexcore(is),50)
!      call readx   (ifsexcore(is),50)
!      read(ifsexcore(is),*)
      do ip = 1,nq
        read(ifsexcore2(is)) isx,qx2,sexcore2(1:ntq,1:ntq,ip) 
        do i  = 1,ntq !     read(ifsexcore(is),*) ixx1,ixx2,ixx3, qxxx1,qxxx2,qxxx3,
          sexcore(i,ip)= sexcore2(i,i,ip)*hartree 
        enddo
      enddo
      allocate( qc (3,ntq,nq) ,rsec(ntq,nq),csec(ntq,nq),sec2(ntq,ntq,nq))  !,eldac (ntq,nq)
      if(is==1) write(6,*)' ###  readin SECU'
      if(is==2) write(6,*)' ###  readin SECD'
!      call readx   (ifsec(is),50)
!      read(ifsec(is),*)
      do ip = 1,nq
        read(ifsec2(is)) isx,qx2,sec2(1:ntq,1:ntq,ip)
        do i  = 1,ntq
!          read(ifsec(is),*) itc(i),ipxxx,isxxx, qc(1:3,i,ip)!, eldac(i,ip)!,rsec(i,ip),csec(i,ip)
          rsec(i,ip)=dreal(sec2(i,i,ip))*hartree
          csec(i,ip)=dimag(sec2(i,i,ip))*hartree
        enddo
      enddo
!      if (sum(abs(itx(1:ntq)-itc(1:ntq)))/=0)  call rx('hqpe.sc itx /= itc ')
!      if (sum(abs(itx(1:ntq)-itxc(1:ntq)))/=0) call rx('hqpe.sc itx /= itxc ')
      qpe1_sc: block !call qpe1_sc(ifqpe(is),iftote(is),iftote2(is),itx,qx,eldax,vxc,sex,sexcore, rsec,csec,deltaw,alat,ef,ntq,nq,is)
        real(8):: eqp(ntq,nq),eqp2(ntq,nq), zfac=1d0 
        character(256):: aaa
        if(is == 1) aaa=' quasiparticle energies MAJORITY'
        if(is == 2) aaa=' quasiparticle energies MINORITY'
        write(iftote(is), *) nq,ntq,ef
        write(iftote2(is),"(2i9,4d24.16)") nq,ntq, ef*rydberg() !, eshtlda!, eshift0, eshift02
        write(ifqpe(is),*) '==============================================================='
        write(ifqpe(is),*) trim(aaa)
        write(ifqpe(is),*) '==============================================================='
        write(ifqpe(is),*) !"('E_shift=', 3d24.16,' eV')")eshtlda,eshift0,eshift02
        write(ifqpe(is),*)
        write(ifqpe(is),"(a)") '           q               state  SEx   SExcore SEc    vxc   ---' &
             // '   dSEnoZ  eQP(starting by lmf)  eHF  Z=1  FWHM=2Z*Simg ReS(elda)'
        do      iq = 1,nq
          do   it = 1,ntq
            eshift      = sex(it,iq)+sexcore(it,iq)+rsec(it,iq)-vxc(it,iq)
            eshift2     = sex(it,iq)+sexcore(it,iq)+rsec(it,iq)-vxc(it,iq)
            eqp(it,iq)  = eldax(it,iq) + eshift  
            eqp2(it,iq) = eldax(it,iq) + eshift2
            fwhm  =  2d0*csec(it,iq) 
            ehf   =  eldax(it,iq) + sex(it,iq)+ sexcore(it,iq) - vxc(it,iq)
            if(eldax(it,iq)<1d20) then
              ehfx = ehf
              dsenoz = eshift2
              if(abs(sex(it,iq))+abs(sexcore(it,iq))+abs(rsec(it,iq))==0d0) then
                dsenoz= 0d0
                ehfx  = 0d0
              endif
              write(ifqpe(is),'(3f9.5,1x,i2,1x,10f7.2,f5.2,f10.5,3x,f10.5)') qx(1:3,it,iq),itx(it),&
                   sex(it,iq),sexcore(it,iq) ,rsec(it,iq),&
                   vxc(it,iq), 0d0, dsenoz, eldax(it,iq), 0d0, 0d0, ehfx,zfac,fwhm, sex(it,iq)+sexcore(it,iq)+rsec(it,iq) 
            endif
            write(iftote(is), "(3f12.7,1x,2i4,1x,4d24.16)") qx(1:3,it,iq),itx(it),iq, eldax(it,iq),eldax(it,iq),eldax(it,iq),zfac
            write(iftote2(is),"(3f12.7,1x,2i4,1x,4d24.16)") qx(1:3,it,iq),itx(it),iq, eldax(it,iq), eqp(it,iq),eqp2(it,iq),  zfac
          end do
          write (ifqpe(is),*)
        end do
      endblock qpe1_sc !  write(6,*)' end of qpe1'
      close(ifqpe(is))
      close(iftote(is))
      close(iftote2(is))
      !------------------------------------------------------------------------------
      ! making SE_ij-VXC_ij, ij are 'basis functions' indexes
      !------------------------------------------------------------------------------
      allocate(se(ntq,ntq,nq),ipiv(nhq),work(nhq*nhq) &
           ,evec_inv(nhq,nhq) ,evec_invt(nhq,nhq))
      allocate(ev_se_ev(ndimsig,ndimsig))
      do ip=1,nq
        do itp=1,ntq
          do itpp=1,ntq    !make Sigma hermitean
            se(itpp,itp,ip)= sex2(itpp,itp,ip)+sexcore2(itpp,itp,ip) +.5d0*(sec2(itpp,itp,ip)+dconjg(sec2(itp,itpp,ip)))
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
        write(ifsigm) 1,ndimsig,n1,n2,n3,nq,0,0,0
      elseif (is == 1) then
        write(ifsigm) nspin,ndimsig,n1,n2,n3,nq,0,0,0
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
!            if( .NOT. exonly) then
              se(itp,itpp,ip)= se(itp,itpp,ip) &
                   -.5d0* sum(dconjg(evec(1:nz,itp,ikpx,is))* matmul(v_xc(1:nz,1:nz,ikpx,is),evec(1:nz,itpp,ikpx,is)))
!            else
!              se(itp,itpp,ip)= se(itp,itpp,ip) &
!                   -wex*.5d0* sum(dconjg(evec(1:nz,itp,ikpx,is))* matmul(v_xc(1:nz,1:nz,ikpx,is),evec(1:nz,itpp,ikpx,is)))
!            endif
          enddo
        enddo
2001  enddo
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
2002  enddo
      eseavrmean = sum(nstar(1:nq)*eseavr(1:nq,is))/nqbz       !call getkeyvalue("GWinput","AddToESEAVR",eseadd,default=0d0,status=ret)
      write(6,"(' ESEAVRmean (used bands above emax_sigm) isp=',d13.6,i2)")eseavrmean,is
      iploop: do 2003 ip=1,nq
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
             ,matmul(se(1:ntqxx(ip),1:ntqxx(ip),ip),evec_inv(1:ntqxx(ip),1:ndimsig2)))   !! sep2013 exprapolation of se by eseavrmean
        do itp =1,ndimsig2
          do itpp=1,ndimsig2
            ev_se_ev(itp,itpp) = ev_se_ev(itp,itpp) + &
                 sum(evec_invt(itp,ntqxx(ip)+1:nevv)*evec_inv(ntqxx(ip)+1:nevv,itpp))*eseavrmean/2d0 ! in Hartree.
          enddo
        enddo 
        do itp=1,ndimsig2 ! - write SE_ij-Vxc_ij where ij are basis function indices
          do itpp=1,ndimsig2
            if( abs(ev_se_ev(itpp,itp) - dconjg(ev_se_ev(itp,itpp)) )> tolq  ) then
              write(6,*)itp,itpp
              write(6,*)ev_se_ev(itpp,itp)
              write(6,*)ev_se_ev(itp,itpp)
              call rx( 'hqpe: Sigma_ij is not hermitian')
            endif
            if(abs(v_xc(itp,itpp,ikpx,is)-dconjg(v_xc(itpp,itp,ikpx,is)))>tolq) call rx( 'hqpe: v_xc is not hermitean')
          enddo
        enddo
        write(ifsigm) qqq(1:3,ikpx,is),is 
        sigmv(:,:,ip) = 1d20
        sigmv(1:ndimsig2,1:ndimsig2,ip)= 2d0*ev_se_ev(1:ndimsig2,1:ndimsig2) !in Ry.
        ! Note 2*ev_se_ev bacause v_xc in sugw.f was in rydberg while SE was in hartree
        write(ifsigm) sigmv(:,:,ip)
2003  enddo iploop
      write(6,*)
      deallocate(sex2,sexcore2,sec2,se,ipiv,work,evec_inv,evec_invt,ev_se_ev)
      deallocate(itxc,qxc,eldaxc,vxc, qc , sexcore ,rsec,csec, itx, qx ,eldax,sex) !,eldac
      if(laf) exit
1001 enddo isloop
    deallocate(v_xc,evec)
    close(ifsigm)
    
    !OUTPUT: Mixing sigm with previous iteration ------------------
    open(newUNIT=ifse_out, file='sigm',form='UNFORMATTED') !Once readin sigma
    allocate(sigma_m(ndimsig,ndimsig,nq,nspin),qqqx_m(3,nq,nspin))
    write(6,*)"========= Sigma mixing section using mixsigma ======="
    call rwsigma ('read',ifse_out,sigma_m,qqqx_m, nspin,ndimsig,n1,n2,n3,nq ) !,eseavr_in)
    if (lsigin .AND. ndimsigin /= ndimsig) then
      deallocate(sigin)
      write(6,*) '... input sigma dimension mismatch ... discarding'
      lsigin = .false.
    endif
    if( .NOT. lsigin) then
      allocate(sigin(1,1,1,1))
      sigin = 0d0
    endif
    call mixsigma(sigma_m, lsigin, sigin, ndimsig**2*nq*nspin)
    write(6,*)
    write(6,*) "=== Write sigm to files ========"
    rewind ifse_out
    call rwsigma ('write',ifse_out,sigma_m,qqqx_m, nspin,ndimsig,n1,n2,n3,nq)
    close(ifse_out)
    if(mpi__rank==0) write(6,ftox) ' OK! hqpe_sc '
    call mpi_finalize(ierr)
  end subroutine hqpe_sc
  subroutine rwsigma(rw,ifss,sigma_m,qqqx_m,nspin,ntq,n1,n2,n3,nq ) !mixing sigm file
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
  subroutine mixsigma(sss, lsigin, sigin, nda) !sigma file mixing
    use m_amix,only: amix
    !- Anderson mixing of a vector
    !i  mmix: number of iterates available to mix
    ! o nmix: nmix > 0: number of iter to try and mix
    !i        nmix < 0: use mmix instead of nmix.
    !o  nmix: (abs)  number of iter actually mixed.
    !o        (sign) <0, intended that caller update nmix for next call.
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
    call getkeyvalue("GWinput","mixpriorit",imix,default=9,status=ret) !  Restrict maximum number of prior iterations
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
end module m_hqpe_sc

! subroutine qpe1_sc (ifqpe,iftote,iftote2,itx,qx, eldax,vxc,sex,sexcore, rsec,csec,deltaw,alat,ef, ntq,nq,is)
!   ! calculates the quasiparticle energies !In QSGW, we have eshift are zero when converged because of self-consistent perturbation.
!   ! E(k,t) = e(k,t) + [SEx(k,t) + SEc(k,t) - xcLDA(k,t)]
!   ! SEx(k,t)   = <psi(k,t)| SEx |psi(k,t)>
!   ! SEc(k,t)   = <psi(k,t)| SEc |psi(k,t)>, SEc = GWc
!   ! xcLDA(k,t) = <psi(k,t)| vxc |psi(k,t)>
!   implicit real*8 (a-h,o-z)
!   integer:: itx(ntq)
!   real(8):: eldax(ntq,nq),vxc(ntq,nq),sex(ntq,nq),sexcore(ntq,nq),rsec(ntq,nq),csec(ntq,nq),qx(3,ntq,nq),&
!        eqp(ntq,nq),eqp2(ntq,nq),wqp(ntq,nq), zfac=1d0 
!   integer:: ifqpe,iftote,iftote2,ntq,jin,nq,is,if1111,ifepq,iq,it, ifeflda,ifefqp,ifefqp1
!   character(256):: aaa
!   if(is == 1) aaa=' quasiparticle energies MAJORITY'
!   if(is == 2) aaa=' quasiparticle energies MINORITY'
!   write(iftote, *) nq,ntq,ef
!   write(iftote2,"(2i9,4d24.16)") nq,ntq, ef*rydberg() !, eshtlda!, eshift0, eshift02
!   write(ifqpe,*) '==============================================================='
!   write(ifqpe,*) trim(aaa)
!   write(ifqpe,*) '==============================================================='
!   write(ifqpe,*) !"('E_shift=', 3d24.16,' eV')")eshtlda,eshift0,eshift02
!   write(ifqpe,*)
!   write(ifqpe,"(a)") '           q               state  SEx   SExcore SEc    vxc   ---' &
!        // '   dSEnoZ  eQP(starting by lmf)  eHF  Z=1  FWHM=2Z*Simg ReS(elda)'
!   do      iq = 1,nq
!      do   it = 1,ntq
!         eshift      = sex(it,iq)+sexcore(it,iq)+rsec(it,iq)-vxc(it,iq)
!         eshift2     = sex(it,iq)+sexcore(it,iq)+rsec(it,iq)-vxc(it,iq)
!         eqp(it,iq)  = eldax(it,iq) + eshift  
!         eqp2(it,iq) = eldax(it,iq) + eshift2
!         fwhm  =  2d0*csec(it,iq) 
!         ehf   =  eldax(it,iq) + sex(it,iq)+ sexcore(it,iq) - vxc(it,iq)
!         if(eldax(it,iq)<1d20) then
!            ehfx = ehf
!            dsenoz = eshift2
!            if(abs(sex(it,iq))+abs(sexcore(it,iq))+abs(rsec(it,iq))==0d0) then
!               dsenoz= 0d0
!               ehfx  = 0d0
!            endif
!            write(ifqpe,'(3f9.5,1x,i2,1x,10f7.2,f5.2,f10.5,3x,f10.5)') qx(1:3,it,iq),itx(it),sex(it,iq),sexcore(it,iq) ,rsec(it,iq),&
!                 vxc(it,iq), 0d0, dsenoz, eldax(it,iq), 0d0, 0d0, ehfx,zfac,fwhm, sex(it,iq)+sexcore(it,iq)+rsec(it,iq) 
!         endif
!         write(iftote, "(3f12.7,1x,2i4,1x,4d24.16)") qx(1:3,it,iq),itx(it),iq, eldax(it,iq),eldax(it,iq),eldax(it,iq), zfac
!         write(iftote2,"(3f12.7,1x,2i4,1x,4d24.16)") qx(1:3,it,iq),itx(it),iq, eldax(it,iq), eqp(it,iq),eqp2(it,iq), zfac
!      end do
!      write (ifqpe,*)
!   end do
! end subroutine qpe1_sc
