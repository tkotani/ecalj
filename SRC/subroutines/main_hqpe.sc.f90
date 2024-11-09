!>  make sigm.* file which contains Sigma-Vxc
!     This uses amix in order to guess a better Sigma-Vxc from previous iterations.
!     * iSigma_en==5 is for diagonal-only Sigma-Vxc in LDA basis set.
!      (Then you need evec0, which contains eigenvector of LDA).
module m_hqpe_sc
  use m_nvfortran
  use m_ftox
  use m_lgunit,only: m_lgunit_init,stdo
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
    use m_genallcf_v3,only: genallcf_v3,laf,nmto=>nlmto ,nspin
    !    use m_readefermi,only: readefermi,ef
    implicit none
    integer:: ifsex(2),ifsexcore(2),ifxc(2),ifsec(2),ifqpe(2),ifsex2(2),ifsexcore2(2),ifsec2(2) !,iftote(2),iftote2(2)
    integer:: i,ierr,ifsigm,iix,ikp,ip,ipxx,iq,iqq,is,isx,j,nevv,isxxx,itp,itpp,nstarsum,nt,ntqx,ntqxxmin
    integer,allocatable :: itx(:) 
    integer,allocatable:: nhqx(:,:)
    integer, allocatable :: ipiv(:)
    real(8):: eavr,eavr2,eee,eseavr0,eseavr02,dsenoz,zfac=1d0
    real(8):: qx2(3) ,qqqx(3) ,del ,mix_fac
    real(8) ::  rydberg,hartree, qqqx0(3)
    logical :: sigma_mixing,lsi, lsigin
    real(8),allocatable :: vxc(:,:), sex(:,:),sexcore(:,:), qx(:,:,:),eldax(:,:),rsec(:,:),csec(:,:) ,&
         qqq(:,:,:) ,qqqx_m(:,:,:),evl(:,:,:),eqp(:,:),sed(:)
    complex(8), allocatable :: sex2(:,:,:),sexcore2(:,:,:), &
         sec2(:,:,:), se(:,:,:),work(:),evec_inv(:,:),evec_invt(:,:), se_ev(:,:), ev_se_ev(:,:),&
         sen(:,:),sen2(:,:),se_in(:,:), v_xc(:,:,:,:),evec(:,:,:,:),evec0(:,:),sigma_m(:,:,:,:), &
         sigin(:,:,:,:), evec00(:,:,:,:) ,evec00inv(:,:)
    integer ::  ifse_out, ifse_in
    logical :: AddDummySig,evec0ex=.false.
    complex(8),allocatable:: sigma_m_out(:,:,:,:),pmat(:,:),pmatd(:,:),sigmv(:,:,:)
    character(2):: soflag
    logical :: mtosigmaonly,nexist !exonly,
    real(8),allocatable:: eseavr(:,:) !,eseavr_in(:,:)
    integer:: nz,ntqmin, ndimsig,ns2,ndimsigin,procid,nrank,ifigwb_,ifigwx1_,ifigwx2_,ifvxcevec,ifevec_, ipx
    character*256:: extn,ext
    integer,allocatable:: ifevec__(:),ifvxcevec_(:),iprocq(:,:),ntqxx(:)
    real(8):: eseavrmean,eseadd,tolq=1d-4
    integer,allocatable:: nev(:,:)
    complex(8),allocatable:: ovl(:,:)
    character(8):: xt
    logical:: nsp2laf
    integer:: ntq,it,n1x,n2x,n3x,nqx,nspinx,nx
    real(8):: ehf,ehfx,eshift,eshift2,fwhm,exx,elow=1d-2
    real(8) :: wex 
    ! call getkeyvalue("GWinput","EXonly",wex,default=0d0,status=ret); exonly = .not.(wex==0d0);  if(exonly) write(6,*)' exonly=T wex=',wex
    call MPI__Initialize()    !this is just for exit routine subroutine rx('...') works well.
    call m_lgunit_init()
    call genallcf_v3(0) 
    call readhamindex()
    call read_BZDATA()
    hartree=2d0*rydberg()
    ndimsig= merge(nmto,nhq, mtosigmaonly())
!!! open files    
    open(newunit=ifxc(1)  , file='XCU')
    open(newunit=ifsex2(1), file='SEX2U',form='UNFORMATTED', status='OLD') !    open(newunit=ifsec(1),file='SECU')
    open(newunit=ifsec2(1), file='SEC2U',form='UNFORMATTED', status='OLD') !    open(newunit=ifsexcore(1) ,file='SEXcoreU')
    open(newunit=ifsexcore2(1),file='SEXcore2U',form='UNFORMATTED',status='OLD')
    nsp2laf= nspin == 2 .AND. .NOT. laf
    if(nsp2laf) open(newunit=ifxc(2) ,     file='XCD')
    if(nsp2laf) open(newunit=ifsex2(2),    file='SEX2D',form='UNFORMATTED', status='OLD')
    if(nsp2laf) open(newunit=ifsec2(2),    file='SEC2D',form='UNFORMATTED', status='OLD')
    if(nsp2laf) open(newunit=ifsexcore2(2),file='SEXcore2D',form='UNFORMATTED',status='OLD')
    open(newunit=ifqpe(1)   ,file='QPU')
    if(nsp2laf) open(newunit=ifqpe(2)   ,file='QPD')
!!! Get ntq. Maximum band index
    call readx(ifxc(1),50)
    read(ifxc(1),*) nspinx,nqx,ntq !readin ntq
    rewind(ifxc(1))
    write(stdo,*)'nqibz nspin=',nqibz,nspin
!!! Read sigm if availabe
    INQUIRE(FILE='sigm', EXIST = lsigin)
    open(newunit=ifsigm, file='sigm',form='UNFORMATTED')
    if(lsigin) then !   write(stdo,ftox) ' ... reading input sigma from file sigm'
      read(ifsigm) nspinx,ndimsigin,n1x,n2x,n3x,nqx
      if((nqx/=nqibz).or.(n1/=n1x).or.(n2/=n2x).or.(n2/=n2x).or.(ndimsigin /= ndimsig)) then
        write(stdo,ftox)' (warning) file mismatch : file nq=',nqx,n1,n2,n3,ndimsig,'  ',nqibz,n1x,n2x,n3x,ndimsigin
        lsigin = .false.
        allocate(sigin(1,1,1,1),source=(0d0,0d0))
      else
        rewind ifsigm
        allocate(sigin(ndimsigin,ndimsigin,nqibz,nspin),qqqx_m(3,nqibz,nspin))
        call rwsigma('read',ifsigm,sigin,qqqx_m,nspinx,ndimsigin,n1x,n2x,n3x,nqx) 
        deallocate(qqqx_m)
      endif
      rewind ifsigm
    endif     !if(laf) nspin=2
!!! Read vxcevec
    write(stdo,ftox)' read from bzdata nqibz; nqibz nq nhq=',nqibz,nqibz,nhq
    write(stdo,ftox)' ndimh ntq nsp nqibz =',nhq,ntq,nspin,nqibz !NOTE this ndimh is the maximum dimention of Hamiltonian.
    allocate(qqq(3,nqibz,nspin),v_xc(nhq,nhq,nqibz,nspin),evec(nhq,nhq,nqibz,nspin),evl(nhq,nqibz,nspin),nev(nqibz,nspin))
    allocate(nhqx(nqibz,nspin))
    Readvxcever: do iq=1,nqibz               !now nqibz is not necessary to be nqbz !nqibz=nqbz
      do is=1,nspin
        open(newunit=ifvxcevec, file='vxcevec'//trim(xt(iq))//trim(xt(is)),form='unformatted')
        read(ifvxcevec) qqq(1:3,iq,is),nz, nev(iq,is)
        read(ifvxcevec) v_xc(1:nz,1:nz,iq,is)
        read(ifvxcevec) evec(1:nz,1:nz,iq,is) !,evl(1:nz,iq,is)
        close(ifvxcevec)
        nhqx(iq,is) = nz   !nz is introduced instead of nhq
        write(stdo,ftox) ' reading vxcevec ... iq is nz=',iq,is,nz
      enddo
    enddo Readvxcever
    if(sum(nstar(:))/= nqbz ) call rx( ' nstarsum/= nqbz')
    HeaderQPU: do is=1,nspin
      write(ifqpe(is),*) '==============================================================='
      write(ifqpe(is),*) ' quasiparticle energies isp=',is
      write(ifqpe(is),*) '==============================================================='
      write(ifqpe(is),*) 
      write(ifqpe(is),*)
      write(ifqpe(is),"(a)") '           q               state   SEx    SExcore   SEc     vxc     ---' &
           // '    dSEnoZ   eQP(starting by lmf)    eHF   Z=1  FWHM=2Z*Simg ReS(elda)'
    enddo HeaderQPU
    nspinx=merge(1,nspin,laf)
    write(ifsigm) nspinx,ndimsig,n1,n2,n3,nqibz,0,0,0
    SETreadinPOINT: do is=1,nspin
      read(ifsex2(is))   
      read(ifsec2(is))   
      read(ifsexcore2(is)) 
      call readx(ifxc(is),50); call readx(ifxc(is),50); read(ifxc(is),*)
    enddo SETreadinPOINT
!!! Main loop for isp
    allocate(sigmv(ndimsig,ndimsig,nqibz), eseavr(nqibz,nspin))
    allocate(vxc(ntq,nqibz), itx(ntq), qx(3,ntq,nqibz),eldax(ntq,nqibz),sex(ntq,nqibz) )
    allocate(sex2(ntq,ntq,nqibz), sexcore(ntq,nqibz),sexcore2(ntq,ntq,nqibz) )
    allocate(rsec(ntq,nqibz),csec(ntq,nqibz),sec2(ntq,ntq,nqibz))
    allocate(eqp(ntq,nqibz), ntqxx(nqibz))!,eqp2(ntq,nqibz))
    allocate(se(ntq,ntq,nqibz),ipiv(nhq),work(nhq*nhq),evec_inv(nhq,nhq) ,evec_invt(nhq,nhq),ev_se_ev(ndimsig,ndimsig),sed(ntq))
    MAINspinloop: do 1001  is = 1,nspin ; write(stdo,ftox) ' --- is=',is
      do 1010 ip = 1,nqibz
        read(ifsex2(is))     isx,qx2,sex2(1:ntq,1:ntq,ip)
        read(ifsexcore2(is)) isx,qx2,sexcore2(1:ntq,1:ntq,ip) 
        read(ifsec2(is))     isx,qx2,sec2(1:ntq,1:ntq,ip)
        do i  = 1,ntq 
          read(ifxc(is),*) itx(i),ipxx,isxxx, qx(1:3,i,ip), eldax(i,ip),vxc(i,ip)
          !  eldax(i,ip)=evl(i,ip,is)*rydberg() -eftrue*rydberg() !We have to use eftrue for G, determined by smearing 
          sex(i,ip)    = sex2(i,i,ip)*hartree           ! in eV
          sexcore(i,ip)= sexcore2(i,i,ip)*hartree 
          rsec(i,ip)   = dreal(sec2(i,i,ip))*hartree
          csec(i,ip)   = dimag(sec2(i,i,ip))*hartree
        enddo
        if(sum((qqq(1:3,ip,is)-qx(1:3,1,ip))**2 ) > tolq ) call rx( 'hqpe.sc: qqq /=qx')
        ntqxx(ip) = findloc([(sexcore2(itp,itp,ip)/=0d0,itp=1,ntq)],back=.true.,value=.true.,dim=1)
        WRITEqpe: do   it = 1,ntq
          eshift     = sex(it,ip)+sexcore(it,ip)+rsec(it,ip)-vxc(it,ip)  !eshift2     = sex(it,ip)+sexcore(it,ip)+rsec(it,ip)-vxc(it,ip)
          eqp(it,ip) = eldax(it,ip) + eshift        !eqp2(it,ip) = eldax(it,ip) + eshift2
          fwhm  =  2d0*csec(it,ip) 
          ehf   =  eldax(it,ip) + sex(it,ip)+ sexcore(it,ip) - vxc(it,ip)
          if(eldax(it,ip)>1d10) cycle ! padding by huge number for it for no data
          ehfx   = merge(0d0,ehf,    abs(sexcore(it,ip))==0d0)
          dsenoz = merge(0d0,eshift, abs(sexcore(it,ip))==0d0)
          write(ifqpe(is),'(3f9.5,1x,i3,1x,10f8.3,f5.2,f10.5,3x,f10.5)') qx(1:3,it,ip),itx(it),&
               sex(it,ip),sexcore(it,ip) ,rsec(it,ip),&
               vxc(it,ip), 0d0, dsenoz, eldax(it,ip), 0d0, 0d0, ehfx,zfac,fwhm, sex(it,ip)+sexcore(it,ip)+rsec(it,ip) 
          !write(iftote2(is),"(3f12.7,1x,2i4,1x,4d24.16)") qx(1:3,it,ip),itx(it),ip, eldax(it,ip), eqp(it,ip),eqp2(it,ip),zfac
        enddo WRITEqpe
        write (ifqpe(is),*)
1010  enddo
      do concurrent(ip=1:nqibz) ! Make SE_ij-VXC_ij, where ij are band index
        nx = ntqxx(ip)
        nz = nhqx(ip,is)
        se(1:nx,1:nx,ip)=&
             sex2(1:nx,1:nx,ip)+ sexcore2(1:nx,1:nx,ip) +.5d0*(sec2(1:nx,1:nx,ip)+dconjg(transpose(sec2(1:nx,1:nx,ip)))) &
             -.5d0* matmul(transpose(dconjg(evec(1:nz,1:nx,ip,is))), matmul(v_xc(1:nz,1:nz,ip,is),evec(1:nz,1:nx,ip,is))) !in Hartree
        nx=ntqxx(ip)
        forall(itp=1:nx) sed(itp)=se(itp,itp,ip)
        eavr2   = sum(eldax(1:nx,ip)**2*sed(1:nx),mask= eldax(1:nx,ip)>elow) &
             /    sum(eldax(1:nx,ip)**2,          mask= eldax(1:nx,ip)>elow)
        eseavr(ip,is) = merge(eavr2,0d0,nx/=0) !in Hartree since sed is in Hartree SquareAverage4extrapolationOFsigma
        write(stdo,ftox)"###  the constant (ESEAVR=e-weighted average Ry)= ",is,ip,2d0*eseavr(ip,is)
      enddo
      eseavrmean = sum(nstar(1:nqibz)*eseavr(1:nqibz,is))/nqbz    
      write(6,"(' ESEAVRmean (exprapolated SE above emax_sigm) isp=',d13.6,i2)")eseavrmean,is
!!! Make inverse evec_inv(n,i) matrix \psi_n=sum_i evec(i,n)\phi_i, where \psi is eigenfunction and \phi is basis function
!!! evec_inv(ib1,iww)= \sum_ib2 ovlinv(ib1,ib2)*dconjg(evec(iww,ib2)), we introduce nev. iww is for PMT basis. ib for band index.
!!! This is for converting rotated evec (=evecrot(ib)) in the representation of original evec(ib).
      SIGMiploop: do 3003 ip=1,nqibz ! Make SE_ij-VXC_ij, where ij are MTO(PMT) basis index
        nz   = nhqx(ip,is) 
        nevv = nev(ip,is)
        ns2  = merge(nmto,nz,mtosigmaonly()) !we allow only mtosigmaonly=T 2024-10-21
        allocate(ovl(nevv,nevv))
        ovl = matmul(dconjg(transpose(evec(1:nz,1:nevv,ip,is))),evec(1:nz,1:nevv,ip,is))
        call matcinv(nevv,ovl) !ovl --> ovlinv
        evec_inv(1:nevv,1:nz) = matmul(ovl(1:nevv,1:nevv),dconjg(transpose(evec(1:nz,1:nevv,ip,is)))) !note ovl means ovlinv
        deallocate(ovl)
        evec_invt(1:nz,1:nevv) = transpose(dconjg(evec_inv(1:nevv,1:nz))) 
        ev_se_ev(1:ns2,1:ns2) = &
             matmul( evec_invt(1:ns2,1:ntqxx(ip)),matmul(se(1:ntqxx(ip),1:ntqxx(ip),ip),evec_inv(1:ntqxx(ip),1:ns2))) & 
             + matmul(evec_invt(1:ns2,ntqxx(ip)+1:nevv),evec_inv(ntqxx(ip)+1:nevv,1:ns2))*eseavrmean ! in Hartree.!extrapo by eseavrmean
        write(ifsigm) qqq(1:3,ip,is),is 
        sigmv(:,:,ip) = 1d20
        sigmv(1:ns2,1:ns2,ip)= 2d0*ev_se_ev(1:ns2,1:ns2) !in Ry.
        ! Note 2*ev_se_ev bacause v_xc in sugw.f was in rydberg while SE was in hartree
        write(ifsigm) sigmv(:,:,ip)
3003  enddo SIGMiploop
      if(laf) exit
      close(ifqpe(is)) !      close(iftote(is))!        close(iftote2(is))
1001 enddo MAINspinloop
    deallocate(v_xc,evec,se,ipiv,work,evec_inv,evec_invt,ev_se_ev,sed)
    close(ifsigm)
!!! OUTPUT: Mixing sigm with previous iteration.  GWinput mixbeta=0.3 should mix new sigm with the weight of 0.3.
    open(newUNIT=ifse_out, file='sigm',form='UNFORMATTED') !Once readin sigma
    allocate(sigma_m(ndimsig,ndimsig,nqibz,nspin),qqqx_m(3,nqibz,nspin))
    write(stdo,ftox)"========= Sigma mixing section using mixsigma ======="
    call rwsigma ('read',ifse_out,sigma_m,qqqx_m, nspin,ndimsig,n1x,n2x,n3x,nqx ) !,eseavr_in)
    call mixsigma(sigma_m, lsigin, sigin, ndimsig**2*nqibz*nspin)
    write(stdo,ftox)
    write(stdo,ftox) "=== Write sigm to files ========"
    rewind ifse_out
    call rwsigma ('write',ifse_out,sigma_m,qqqx_m, nspin,ndimsig,n1,n2,n3,nqibz)
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
    if(ifs>0) write(stdo,ftox) " rwsigma: Reading sigm (Binary)==="
    if(ifs<0) write(stdo,ftox) " rwsigma  Writing sigm (Binary)==="
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
    write(stdo,ftox)' === rwsigma:  sum check of sigma_m=',sum(abs(sigma_m))
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
    write(stdo,ftox)'('' mixsigma: Anderson mixing sigma with mixing beta ='',f12.6)',beta
    allocate ( a(2*nda,0:mxsav+1,2) )
    fff="mixsigma"
    INQUIRE (FILE =fff, EXIST = fexist)
    if(fexist)      write(stdo,ftox)'... reading file mixsigma'
    if( .NOT. fexist) write(stdo,ftox)'... No file mixsigma'
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
      write(stdo,ftox)'... using input sigma read from sigm file'
      a(1:nda,0,2)       = dreal(sigin)        !input
      a(nda+1:2*nda,0,2) = dimag(sigin)  !input
    endif
    call getkeyvalue("GWinput","mixpriorit",imix,default=9,status=ret) !  Restrict maximum number of prior iterations
    mmix = min(max(nitr-1,0),imix)
    if (mmix > mxsav) mmix = mxsav
    call getkeyvalue("GWinput","mixtj",acc,default=0d0,status=ret)
    if(acc/=0d0) then
      write(stdo,ftox)' readin mixtj from GWinput: mixtj=',acc
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
