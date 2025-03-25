!>plot wannier wanplot in m_wanplot is main module but not checked well.
module m_wan_lmf2gw !Note this is now only for wanplot for backward comatibility. Not in gwsc. So may unused varirables.
  !> lmf2gw() set variables to module variables by reading files from following input files.
  integer,allocatable,protected:: nindx(:),lindx(:),ibasindx(:),mnla(:,:) !,iantiferro(:)
  integer,protected :: nbandmx,nphimx,&
       nsp,       &  !=1 or 2, corresponding to para or ferro.
       nbas,      &  !Number of atom in the primitive cell
       nclass,    &  !Number of atomic class (or type) for in the primitive cell
       nrmx,      &  != maxval(nr(1:nclass))  Maximum number of nr
       ncoremx,   &  != maxval(ncore(1:nclass))
       lmxamx,    &  != maxval(lmxa(1:nclass))
       ngpmx,     &  !Maximum number of G vector.
       ldim2    ! = total number of augmentation functions nRlm
  character(8),allocatable,protected:: spid(:)
  integer,allocatable:: &
       iclass(:),lmxa_d  (:),nr(:),konf_d(:,:),ncore_d(:),ibasf(:)
  real(8),allocatable :: &
       zz(:),aa(:),bb(:),ec_d (:,:,:),evl_d(:,:,:), &
       gx_d(:,:,:,:,:),gcore_d(:,:,:,:), bas(:,:)
  real(8):: plat(3,3), alat, efermi,qval
  complex(8),allocatable:: cphi_d(:,:,:,:)
  complex(8),allocatable:: geig_d(:,:,:,:)
  !      real(8),allocatable :: qirr(:,:)
  logical,protected:: laf   !! - laf: antiferro switch
  integer,protected:: ngcmx,nqnum,nqnumc,nqtt,nq0i,iq0pin,nq0iadd,nqbz,nqibz,nqbzx
  real(8):: QpGcut_psi,QpGcut_cou
  real(8),allocatable :: wt(:),q0i(:,:)
  integer,allocatable,target:: ngvecptt(:,:,:),ngvecctt(:,:,:),ngptt(:),ngctt(:)
  real(8),allocatable:: qtt(:,:)
  !  integer,allocatable:: ngplist(:),ndimhall(:),iqindex(:)
  !  real(8),allocatable:: qplist(:,:)
  integer:: nqirr     ! = Number of q points for irr=1 (see m_qplist, output of qg4gw).
  real(8),allocatable,protected :: qibz(:,:)
contains
  subroutine lmf2gw() !read atomic part wanplotatom.dat written in sugw.f90. This is not clean historically.
    !wanplot is expected to be unsuppported 2024-6-18
    use m_keyvalue,only: getkeyvalue
    use m_hamindex0,only: readhamindex0,nclass_in=>nclass,iclass_in=>iclasst
    use m_hamindex0,only: nindx_in=>nindx,lindx_in=>lindx,ibasindx_in=>ibasindx,nphimx_in=>nphimx
    implicit none
    integer:: iq0p
    integer:: ldim,      & ! = sum ( (lmxa(1:nbas)+1)**2 )
         nband    ! Number of bands given by GWIN0
    integer :: icor1,icorex,i,i1,i2,ibas,ibasx,ibx,ic,icore, &
         ifichkv,ifigw0,ifigwa,ifigwb,ifigwx1,ifigwx2,ifigwx3,isp,ispx, &
         ispxx,ix,kkk,kkkdummy,l,ldummy,lmxa,lxx,nclassx,m,n, &
         ncore,ndimh,ngp,nnc,nspdummy,IKP,NR_A !takao feb2012 ngc,ngcmx,
    real(8) ovv(20),ef0,z,a_a,b_a,rofi_anr
    character(120) ::  ext0(256), ext(256)
    real(8):: qqq(3),qxx(3), vvvv(18)
    integer:: ifi,ifefclass,icors(2)
    complex(8),allocatable:: zegf(:,:) ,geig(:,:)
    complex(8),allocatable:: cphi(:,:)
    real(8),allocatable:: evl(:), vvv1(:),vvv2(:),vvv3(:),rofi_A(:),gcore_A(:,:), ec_A(:)
    integer,allocatable:: konf(:,:),nncx(:,:),ngvecp(:,:),lmxaa(:)
    real(8),parameter ::  rydberg=13.6058d0
    ! nocore is obtained by inquire(file='NoCore',exist=nocore) in the upper routine.
    ! If nocore exist. you have to supply
    !  <psi|Vxc(n_valence)|psi>  to  vxclda (nband, nqirr).
    ! If not, <psi|Vxc(n_total)|psi>  to  vxclda.
    !----------------------------------------------
    integer:: ificg
    integer:: procid,nrank,ifigwb_,ifigwx1_,ifigwx2_
    integer:: iq,iqq,iqqx,nxxx,ifibz
    character*256:: extn,aaa,fname
    integer,parameter :: nsize= 1000000
    integer:: ifiproc,nqixx,nspxx,numprocxx,ixxx,ifiqibz
    integer::  id,nsizex,iqqxx,ib,ii,ipqn,nn,nnn(3),ifiqg,ifiqgc,irr,irrq,iqibz
    !! =================================================================
    open(newunit=ifigwa,file='wanplotatom.dat',form='unformatted')
    read (ifigwa) nbas,nsp,ldim2,nbandmx,lmxamx,ncoremx,nrmx,plat,alat!,nqirr
    allocate(bas(3,nbas),lmxaa(nbas))!,qplist(3,nqirr),ngplist(nqirr),ndimhall(nqirr))
    read(ifigwa) bas,lmxaa!,qplist,ngplist,ndimhall,qval
    call readhamindex0()
    nclass=nclass_in
    allocate(nindx(ldim2),lindx(ldim2),ibasindx(ldim2))
    nindx=nindx_in
    lindx=lindx_in
    ibasindx=ibasindx_in
    nphimx=nphimx_in
    allocate( iclass(nbas) )
    iclass=iclass_in
    allocate(lmxa_d(nclass), nr(nclass), ncore_d(nclass), konf_d(0:lmxamx,nclass), zz(nclass),aa(nclass),bb(nclass) )
    lmxa_d(iclass(1:nbas)) = lmxaa(1:nbas)
    !! ATOMIC PART ic = ibas scheme ==,  GET nrxx and ncoremx ----------------------
    !    open(newunit=ifigwa,file='gwa',form='unformatted')
    allocate(nncx(0:lmxamx,nbas),konf(lmxamx+1,nbas),spid(nbas),ec_d(ncoremx, nclass, nsp),&
         gx_d(nrmx,0:lmxamx,nphimx,nclass,nsp), gcore_d(nrmx,ncoremx,nclass,nsp)  )
    do 3001 ibas = 1, nbas
      read(ifigwa) z, nr_A, a_A, b_A, rofi_Anr,lmxa,nspdummy,ncore,spid(ibas)
      allocate(rofi_A(nr_A), gcore_A(nr_A,ncore),ec_A(ncore))
      read(ifigwa) konf(1:lmxa+1,ibas)
      read(ifigwa) rofi_A(1:nr_A)
      write(6,"(' site',i3,'  z=',f5.1,'  rmax=',f8.5,'  lmax=',i1,'  konf=',10i1)")ibas,z,rofi_A(nr_A),lmxa,konf(1:lmxa+1,ibas)
      ic = iclass(ibas)
      zz(ic)= z
      aa(ic)= a_A
      bb(ic)= b_A
      nr(ic)= nr_A
      ncore_d(ic) = ncore/nsp
      konf_d(0:lmxa,ic) = konf(1:lmxa+1,ibas)
      write(6,"('  l    g(rmax)    gp(rmax)',4x,'<g g>',9x,'<gp gp>',9x,'<g gp>')")
      do  l = 0, lmxa
        do  isp = 1, nsp
          read(ifigwa) lxx,ispxx
          if(lxx /= l .OR. isp /=ispxx) call rx('lmf2gw:lxx or isp wrong')
          read(ifigwa) gx_d(1:nr_A,l,1,ic,isp) !phi
          read(ifigwa) gx_d(1:nr_A,l,2,ic,isp) !phidot
          if (konf_d(l,ic) >= 10) read(ifigwa) gx_d(1:nr_A,l,3,ic,isp) !phiz
        enddo
      enddo
      if(ncore/=0) write(6,'(''  l  k isp       ecore      gc(rmax)     <gc gc>'')')! core part
      icore = 0
      icors = 0
      do isp = 1, nsp
        do l = 0, lmxa
          nncx(l,ibas) = mod(konf(l+1,ibas),10)-1 -(l+1) +1
          nnc          = max(nnc,nncx(l,ibas))
          do kkk = l+1, mod(konf(l+1,ibas),10)-1
            icore = icore+1
            icors(isp) = icors(isp) +1
            icor1=icors(isp)
            read(ifigwa) icorex,ldummy,ispx,kkkdummy,ec_A(icore)
            if(icore/=icorex)  call rx('lmf2gw:icore/=icorex')
            read(ifigwa) gcore_A(1:nr_A,icore) ! gcore
            ec_d(icor1, ic, isp) = ec_A(icore)
            gcore_d(1:nr_A,icor1,ic,isp)  = gcore_A(1:nr_A,icore)
          enddo
        enddo
      enddo
      deallocate(rofi_A,gcore_A,ec_A)
3001 enddo
    close(ifigwa)
  end subroutine lmf2gw
end module m_wan_lmf2gw
module m_wan_qg ! Information for plane-wave basis
  implicit none
  integer ::  nqnum, ngpmx_qg, nnnn
  integer ::  nqnumc, ngcmx
  double precision :: QpGcut_psi,QpGcut_Cou
  double precision,allocatable :: qqqa(:,:),qqqb(:,:)
  integer,allocatable :: ngp(:),ngvecp(:,:,:)
  integer,allocatable :: ngc(:),ngvecc(:,:,:)
contains
  subroutine read_QG()
    implicit none
    integer :: ikp
    integer :: ifiqg,ifiqgc
    write(6,*) '--- read_QG ---'
    open(newunit=ifiqg ,file='QGpsi',form='unformatted')
    open(newunit=ifiqgc,file='QGcou',form='unformatted')
    read(ifiqg  ) nqnum , ngpmx_qg, QpGcut_psi,nnnn
    read(ifiqgc ) nqnumc, ngcmx, QpGcut_Cou
    write(6,*) 'nqnum,nqnumc=',nqnum,nqnumc
    write(6,*) 'QpGcut_psi QpGcutCou nnnn=' &
         , QpGcut_psi,QpGcut_Cou ,nnnn
    if (nqnum /= nqnumc) then
      write(6,*) 'Error : nqnum!=nqnumc'
      write(6,*) 'nqnum,nqnumc=',nqnum,nqnumc
      stop 'Error : nqnum!=nqnumc'
    endif
    allocate(qqqa(3,nqnum),qqqb(3,nqnum))
    allocate(ngp(nqnum),ngc(nqnum))
    allocate( ngvecp(3,ngpmx_qg,nqnum),ngvecc(3,ngcmx,nqnum) )
    ngvecp(1:3,1:ngpmx_qg,1:nqnum)=0
    ngvecc(1:3,1:ngcmx,1:nqnum)=0
    do ikp = 1,nqnum
      read (ifiqg)  qqqa(1:3,ikp), ngp(ikp)
      read (ifiqgc) qqqb(1:3,ikp), ngc(ikp)

      !        write(6,"(i5,3f8.4,f10.5)")
      !     &       ikp,qqqa(1:3,ikp),sum(abs(qqqa(1:3,ikp)-qqqb(1:3,ikp)))
      if (sum(abs(qqqa(1:3,ikp)-qqqb(1:3,ikp))) > 1.0d-8) then
        stop 'qqqa!=qqqb'
      endif

      read (ifiqg ) ngvecp(1:3,1:ngp(ikp),ikp)
      read (ifiqgc) ngvecc(1:3,1:ngc(ikp),ikp)
    enddo
    close(ifiqg)
    close(ifiqgc)
  end subroutine read_QG
end module m_wan_qg
module m_wanutil
  public myinv3,mymatvec,calc_ylm,calc_phiall_abc2,b2w,wrt_xsf,expand_mesh,calc_npw
  private
contains
  subroutine calc_npw(nfac, npw)
    use m_wan_QG,only: ngvecp,qqqa,nqnum,ngp
    use m_wan_lmf2gw,only: alat,plat
    implicit none
    ! input
    integer :: nfac
    ! output
    integer :: npw(3)
    ! local
    integer :: iq,ig,id,itmp(3),ntmp(3)
    double precision :: pi,gtmp(3),gcutmax,gcuttmp,at(3,3),g(3,3)
    logical :: debug=.false.
    write(*,"(a)") '--- calc_npw ---'
    !call mytranspose(plat,At,3,3)
    At=transpose(plat)
    call myinv3(At,G)
    pi=4.0d0*atan(1.0d0)
    ntmp(1:3)=0
    do iq=1,nqnum
      gcutmax=-1.0d0
      do ig=1,ngp(iq)
        call mymatvec(G,dble(ngvecp(1:3,ig,iq)),gtmp,3,3)
        gtmp(1:3)=gtmp(1:3)+qqqa(1:3,iq)
        gtmp(1:3)=gtmp(1:3)*2.0d0*pi/alat
        gcuttmp=sqrt(sum(gtmp(1:3)**2))
        if (gcutmax < gcuttmp) gcutmax=gcuttmp
        do id=1,3
          itmp(id)=abs(ngvecp(id,ig,iq))
          if (ntmp(id) < itmp(id)) ntmp(id)=itmp(id)
        enddo
      enddo
      if(debug) write(*,"(a,2i5,f10.5)") '# iq ngp gcutmax= ',iq,ngp(iq),gcutmax
    enddo
    !      npw(1:3)=2*ntmp(1:3)+2
    npw(1:3)=nfac*ntmp(1:3)+2
    write(*,"(a,3i6)") '# npw(1:3)=',npw(1:3)
  end subroutine calc_npw

  !-----------------------------------------------------
  ! Linear interpolation of gx/r
  double precision function calc_gxr(r,l,n,ic,isp)
    !      use m_LMTO
    use m_wan_lmf2gw,only: bb,nr,aa,alat,gx=>gx_d
    implicit none
    ! input
    double precision :: r
    integer :: l,n,ic,isp
    ! local
    double precision :: r1,r2
    integer :: ir
    ir=1+int(log(r/bb(ic)+1.0d0)/aa(ic))
    if (ir < 1) stop 'ir < 1'
    if (ir > nr(ic)-1) stop 'ir > nr(ic)-1'
    r1=bb(ic)*(exp((ir-1)*aa(ic))-1d0)
    r2=bb(ic)*(exp((ir  )*aa(ic))-1d0)
    if (r1 > r) stop 'r1 > r'
    if (r2 <= r) stop 'r2 <= r'
    calc_gxr=(r-r2)/(r1-r2)*gx(ir,l,n,ic,isp) &
         + (r-r1)/(r2-r1)*gx(ir+1,l,n,ic,isp)
    calc_gxr=calc_gxr/(r+1.0d-20)
  end function calc_gxr

  !-----------------------------------------------------
  subroutine b2w(nq_wfn,nband_wfn,q_wfn,bindx_wfn,tlat,npw, &
       phi,wan)
    !! Make Wannier functions from Bloch functions in real space representation.
    use m_wan_lmf2gw,only: bb,nr,aa,alat,nsp,plat
    implicit none
    integer :: nq_wfn,nband_wfn,npw(3),bindx_wfn(nband_wfn),tlat(3)
    double precision :: q_wfn(3,nq_wfn),tvec(3),phase,pi,rtmp(3)
    double complex :: &
         phi(npw(1)+1,npw(2)+1,npw(3)+1,nband_wfn,nq_wfn,nsp), &
         wan(npw(1)+1,npw(2)+1,npw(3)+1,nband_wfn,nsp) &
         ,ephase
    integer :: iq,isp
    ! debug:
    !      wan(:,:,:,:,1) = phi(:,:,:,:,2,1)
    !      return
    pi = 4.0d0*atan(1.d0)
    rtmp(:) = dble(tlat(:))
    call mymatvec(plat,rtmp,tvec,3,3)
    tvec(1:3)=alat*tvec(1:3)
    wan = (0.0d0,0.0d0)
    do isp = 1,nsp
      do iq = 1,nq_wfn
        phase=2.0d0*pi/alat*sum(q_wfn(1:3,iq)*tvec(1:3))
        ephase=dcmplx(cos(phase),-sin(phase))
        wan(:,:,:,:,isp) = wan(:,:,:,:,isp) + phi(:,:,:,:,iq,isp)*ephase
      enddo ! iq
    enddo ! isp
    wan = wan / dble(nq_wfn)
  end subroutine b2w
  !      module m_wfrho_abc
  !      contains
  subroutine calc_phiall_abc2(nq_wfn,nband_wfn,q_wfn,bindx_wfn, &
       npw,mesh,nsp,nband,ldim2,ngpmx, &
       geig,cphi,nwf, &
       phipw,phiaug,phitot)
    !      use m_readeigen,only: readcphif,readgeig
    ! cccccccccccccccccccccccccccccccccccc
    use m_wan_qg,only:ngp
    ! cccccccccccccccccccccccccccccccccccc
    !      use m_LMTO
    !      use m_FFT3D
    implicit none
    ! inputs
    integer :: nq_wfn,nband_wfn,bindx_wfn(nband_wfn),nsp
    double precision :: q_wfn(3,nq_wfn)
    integer :: npw(3),mesh(3),ldim2,nband,ngpmx
    ! outputs
    double complex :: &
         phipw(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp), &
         phiaug(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp), &
         phitot(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp)
    integer :: isp,iq_wfn,ib,i1,i2,i3
    integer :: augregion(4,mesh(1)+1,mesh(2)+1,mesh(3)+1)
    double complex :: phipwtmp(mesh(1)+1,mesh(2)+1,mesh(3)+1), &
         phiaugtmp(mesh(1)+1,mesh(2)+1,mesh(3)+1)
    double complex :: eikr(mesh(1)+1,mesh(2)+1,mesh(3)+1), &
         eikT(mesh(1)+1,mesh(2)+1,mesh(3)+1)
    real(8):: q(3),quu(3)
    logical :: debug=.false.
    integer:: nwf
    complex(8):: geig(ngpmx,nwf,nq_wfn,nsp),cphi(ldim2,nwf,nq_wfn,nsp)

    !------------------------------------------------------------------
    if(debug) write(*,"(a)") '--- calc_phiall_abc2 ---'

    !      print *,'eeee dim =', mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp

    !      call fft_init(npw,'B')
    call calc_augregion_abc(mesh(1),mesh(2),mesh(3),augregion)
    !      allocate(geig2(ngpmx,nband))
    !      allocate(cphi2(ldim2,nband))
    ! omp parallel do private(iq, eikr,eikT, phipwtmp,phiaugtmp )
    do iq_wfn=1,nq_wfn
      q = q_wfn(1:3,iq_wfn)
      !         write(6,"('ffff iq q=',i5,3f10.4)") iq_wfn,q
      call calc_eikreikT_abc(q,mesh,augregion,eikr,eikT)
      do isp=1,nsp
        !            call readgeig(q,ngpmx,isp, quu, geig2)
        !            if(sum(abs(q-quu))>1d-6) stop 'mmlf111eeeee'
        !            call readcphi(q,ldim2,isp, quu, cphi2)
        !            if(sum(abs(q-quu))>1d-6) stop 'mmlf222eeeee'
        do ib=1,nband_wfn
          ! omp critical
          !     write(*,"(a,i2,2i5,3f10.4,i5)")
          !     &         '# isp,iq_wfn,iq,q,ib=',isp,iq_wfn,iq,qtt(1:3,iq),ib
          !               write(*,"(a,i2,2i5,3f10.4,i5)")
          !     &         '# isp,iq_wfn,iq,q,ib=',isp,iq_wfn,iq_wfn,q_wfn(1:3,iq_wfn),ib
          ! ccccccccccccccccccccccccc
          !      if(iq_wfn>=63) write(6,*)'bbbbb222     ',ib,iq_wfn,isp,ngp(iq_wfn)
          !      if(iq_wfn>=63) write(6,*)'bbbbb222 sum ',sum(geig(1:ngp(iq_wfn),ib,iq_wfn,isp))
          ! ccccccccccccccccccccccccc

          ! omp end critical
          call calc_phi_abc2(geig(1,ib,iq_wfn,isp),ngpmx,cphi(1,ib,iq_wfn,isp),ldim2, &
               iq_wfn,npw,mesh,augregion, &
               phipwtmp,phiaugtmp)
          !               write(6,*)'sumccc=',sum(abs(phipwtmp)),sum(abs(phiaugtmp))
          do i3=1,mesh(3)+1
            do i2=1,mesh(2)+1
              do i1=1,mesh(1)+1 !   bloch function
                ! cccccccccccccccccccccccc
                !         write(6,"('eeeee iq =',100i5)") iq_wfn,nq_wfn,i1,i2,i3,ib,isp
                ! cccccccccccccccccccccccc
                phipw(i1,i2,i3,ib,iq_wfn,isp) =  eikr(i1,i2,i3) *phipwtmp(i1,i2,i3)
                phiaug(i1,i2,i3,ib,iq_wfn,isp)=  eikT(i1,i2,i3) *phiaugtmp(i1,i2,i3)
                phitot(i1,i2,i3,ib,iq_wfn,isp)= &
                     phipw(i1,i2,i3,ib,iq_wfn,isp)+ &
                     phiaug(i1,i2,i3,ib,iq_wfn,isp)
                ! cccccccccccccccccccccccc
                !         write(6,"('ggggg iq =',100i5)") iq_wfn,nq_wfn,i1,i2,i3,ib,isp
                ! cccccccccccccccccccccccc
              enddo
            enddo
          enddo
          !       write(*,'(6f10.5)') phitot(:,:,:,ib,iq_wfn,isp)
        enddo               !ib
      enddo                  !isp

    enddo                     !iq
  end subroutine calc_phiall_abc2
  subroutine calc_augregion_abc(n1,n2,n3,augregion)
    use m_wan_lmf2gw,only: bb,nr,aa,alat,iclass,nclass,bas,nbas,plat
    !      use m_genallcf_v3,only: nbas=>natom, bas=>pos,plat
    implicit none
    ! input
    integer :: n1,n2,n3
    ! output
    integer :: augregion(4,n1+1,n2+1,n3+1)
    ! local
    integer :: nshell
    parameter (nshell=4)
    integer :: i1,i2,i3,j1,j2,j3,ibas,ic
    double precision :: rmax,ratom(3),r(3),rtmp(3),dr
    logical:: debug=.false.

    write(*,*) '--- calc_augregion ---',nclass
    augregion(:,:,:,:)=0

    do ibas=1,nbas
      ic=iclass(ibas)
      rmax = bb(ic)*(exp((nr(ic)-1)*aa(ic))-1d0)
      write(6,*)'ibas, rmax=',ibas,ic,rmax
      do j1=-nshell,nshell
        do j2=-nshell,nshell
          do j3=-nshell,nshell
            rtmp(1)=j1
            rtmp(2)=j2
            rtmp(3)=j3
            call mymatvec(plat,rtmp,ratom,3,3)
            ratom(1:3)=alat*(ratom(1:3)+bas(1:3,ibas))

            do i3=1,n3+1
              do i2=1,n2+1
                do i1=1,n1+1

                  rtmp(1)=(i1-1)/dble(n1)
                  rtmp(2)=(i2-1)/dble(n2)
                  rtmp(3)=(i3-1)/dble(n3)
                  !            call mymatvec(plat,rtmp,r,3,3)
                  !            r(:) = rini(:) + (rfin(:)-rini(:))*rtmp(:)
                  r(:) = plat(:,1)*rtmp(1)+plat(:,2)*rtmp(2)+plat(:,3)*rtmp(3)
                  r(1:3)=alat*r(1:3)
                  dr=sqrt(sum((r(1:3)-ratom(1:3))**2))
                  if (dr < rmax) then
                    if (augregion(4,i1,i2,i3) /= 0) then
                      stop 'calc_augregion_abc: Overlap in augmented region!'
                    endif
                    augregion(1,i1,i2,i3)=j1
                    augregion(2,i1,i2,i3)=j2
                    augregion(3,i1,i2,i3)=j3
                    augregion(4,i1,i2,i3)=ibas
                  endif
                enddo !i1
              enddo !i2
            enddo !i3
          enddo !j3
        enddo !j2
      enddo !j1
    enddo !ibas
  end subroutine calc_augregion_abc
  subroutine calc_eikreikT_abc(kvec,mesh, &
       augregion,eikr,eikT)
    !      use m_LMTO
    use m_wan_lmf2gw,only: bb,nr,aa,alat,nsp,plat
    implicit none
    ! input
    double precision :: kvec(3)
    integer :: mesh(3),augregion(4,mesh(1)+1,mesh(2)+1,mesh(3)+1)
    ! output
    double complex ::  eikr(mesh(1)+1,mesh(2)+1,mesh(3)+1), &
         eikT(mesh(1)+1,mesh(2)+1,mesh(3)+1)
    ! local
    integer :: i1,i2,i3
    double precision :: rtmp(3),r(3),tvec(3)
    double precision :: phase,pi
    pi=4.0d0*atan(1.0d0)
    if( .FALSE. ) write(*,*) 'kvec=',kvec (1:3)
    ! Calculate e^{ikr}
    do i3=1,mesh(3)+1
      do i2=1,mesh(2)+1
        do i1=1,mesh(1)+1
          rtmp(1)=(i1-1)/dble(mesh(1))
          rtmp(2)=(i2-1)/dble(mesh(2))
          rtmp(3)=(i3-1)/dble(mesh(3))
          !        call mymatvec(plat,rtmp,r,3,3)
          !        r(:) = rini(:) + (rfin(:)-rini(:))*rtmp(:)
          r(:)= plat(:,1)*rtmp(1)+plat(:,2)*rtmp(2)+plat(:,3)*rtmp(3)
          r(1:3)=alat*r(1:3)
          phase=2.0d0*pi/alat*sum(kvec(1:3)*r(1:3))
          eikr(i1,i2,i3)=dcmplx(cos(phase),sin(phase))
        enddo
      enddo
    enddo

    ! Calculate e^{ikT}
    do i3=1,mesh(3)+1
      do i2=1,mesh(2)+1
        do i1=1,mesh(1)+1

          if (augregion(4,i1,i2,i3) /= 0) then
            rtmp(1:3)=augregion(1:3,i1,i2,i3)
            !             tvec(i) =plat(i,j)*rtmp(j)
            call mymatvec(plat,rtmp,tvec,3,3)
            tvec(1:3)=alat*tvec(1:3)
            !  2 pi  k(i)*tvec(i)
            phase=2.0d0*pi/alat*sum(kvec(1:3)*tvec(1:3))
            eikT(i1,i2,i3)=dcmplx(cos(phase),sin(phase))
          else
            eikT(i1,i2,i3)=0.0d0
          endif
        enddo
      enddo
    enddo

  end subroutine calc_eikreikT_abc
  subroutine calc_phi_abc2(geig,ngpmx,cphi,ldim2,iq, npw,mesh, &
                                !!-- Plane wave expansion of an eigenfunciton (geig,cphi).
       augregion, &
       phipwtmp,phiaugtmp)
    use m_wan_QG,only:ngvecp,ngp
    use m_wan_lmf2gw,only:mnla,iclass,nbas,bas,plat,alat
    !      use m_genallcf_v3,only: nbas=>natom, bas=>pos,plat,alat
    implicit none
    integer :: isp,iq,npw(3),mesh(3),ngpmx,ldim2 !,iband
    integer :: augregion(4,mesh(1)+1,mesh(2)+1,mesh(3)+1)
    double precision :: qlat(3,3)
    double complex :: eigr,ci = (0.0d0,1.0d0)
    double complex :: &
         phipwtmp(mesh(1)+1,mesh(2)+1,mesh(3)+1), &
         phiaugtmp(mesh(1)+1,mesh(2)+1,mesh(3)+1)
    integer :: itmp(3),ig,id,i1,i2,i3,j1,j2,j3,ii
    double precision :: rtmp(3),r(3),r0(3) !points to plot
    double precision :: ratom(3) ! atomic points
    double precision :: dr(3)
    complex(8):: geig(ngpmx),cphi(ldim2)
    integer,parameter :: lmax=6
    double complex :: Y(2*lmax+1,lmax+1)
    double precision :: Yreal(2*lmax+1,lmax+1)
    !  double precision :: calc_gxr
    double precision :: drlength,theta,pphi,sintheta
    integer :: idim,il,mtmp,ntmp,ltmp, ibas
    double complex :: phia
    double precision :: pi=4.0d0*atan(1.0d0),tpi = 8d0*atan(1.0d0)
    logical :: debug=.false.
    if(debug) write(6,*)' --- calc_phi_abc2 ---'
    ! ccccccccccccccccccccccccc
    !      if(iq>=63) write(6,*)'aaaaaa222     ',ngp(iq)
    !      if(iq>=63) write(6,*)'aaaaaa222 sum ',sum(geig(1:ngp(iq)))
    ! ccccccccccccccccccccccccc

    call minv33tp(plat,qlat)
    !      call chkinv33(plat,qlat)
    phipwtmp = 0d0
    !      write(6,*)'aaaaaa1',mesh
    do i3=1,mesh(3)+1
      do i2=1,mesh(2)+1
        do i1=1,mesh(1)+1
          !         write(6,*)'aaaaaa',i1,i2,i3
          rtmp(1)=(i1-1)/dble(mesh(1))
          rtmp(2)=(i2-1)/dble(mesh(2))
          rtmp(3)=(i3-1)/dble(mesh(3))
          !         r(:) = rini(:) + (rfin(:)-rini(:))*rtmp(:)
          r(:) = plat(:,1)*rtmp(1)+plat(:,2)*rtmp(2)+plat(:,3)*rtmp(3)
          !         r0(:) = matmul(qlat,r)
          do ii=1,3
            r0(ii) = sum(qlat(:,ii)*r(:))
          enddo ! ii
          !   r0(i)=G0(j,i)*r(j)*
          !   G(i)= G0(j,i)*nG(i)
          !   exp (i 2 pi  G(i)*r(i) )
          !         if(iq>=63) write(6,*)'aaaaaa222     ',i1,i2,i3,iq,ngp(iq)
          !         if(iq>=63) write(6,*)'aaaaaa222 sum ',sum(geig(1:ngp(iq)))

          !, size(ngvecp(:,:,iq)), size(ngvecp(1,:,iq))
          do ig=1,ngp(iq)
            !           if(iq>=63) write(6,*)'ig ngpmx',ig
            eigr=exp(ci*tpi*sum(r0(:)*dble(ngvecp(:,ig,iq))))
            !           if(iq>=63) write(6,*)'eigr',eigr
            !           if(iq>=63) write(6,*)'geig',geig(ig)
            phipwtmp(i1,i2,i3) = phipwtmp(i1,i2,i3) &
                 + eigr*geig(ig) !,iband) !,iq,isp)
            !           if(iq>=63) write(6,*)'end phip'
            !        phipwtmp(i1,i2,i3)=out_fft(mod(i1-1,npw(1))+1,
            !     &       mod(i2-1,npw(2))+1,mod(i3-1,npw(3))+1)
          enddo ! ig
          !         if(iq>=63) write(6,*)'aaaaaa3333',i1,i2,i3,ngp(iq)
        enddo ! i1
      enddo ! i2
    enddo ! i3
    ! Augmented part
    if(debug) write(6,*)' ----- goto augmented part ------------------'
    phiaugtmp(:,:,:)=0.0d0
    do i3=1,mesh(3)+1
      do i2=1,mesh(2)+1
        do i1=1,mesh(1)+1
          if (augregion(4,i1,i2,i3) /= 0) then
            !          write(6,*)i1,i2,i3,mesh(1)+1,mesh(2)+1,mesh(3)+1
            ! set plane-wave part to zero
            !          phiaugtmp(i1,i2,i3)=0.0d0
            rtmp(1)=(i1-1)/dble(mesh(1))
            rtmp(2)=(i2-1)/dble(mesh(2))
            rtmp(3)=(i3-1)/dble(mesh(3))
            !          call mymatvec(plat,rtmp,r,3,3)
            !          r(1:3)=alat*r(1:3)
            !          r(:) = rini(:) + (rfin(:)-rini(:))*rtmp(:)
            r(:) = plat(:,1)*rtmp(1)+plat(:,2)*rtmp(2)+plat(:,3)*rtmp(3)
            r(1:3)=alat*r(1:3)
            rtmp(1:3)=augregion(1:3,i1,i2,i3)
            call mymatvec(plat,rtmp,ratom,3,3)
            ratom(1:3)=alat*(ratom(1:3)+bas(1:3,augregion(4,i1,i2,i3)))
            dr(1:3)=r(1:3)-ratom(1:3)
            drlength=sqrt(sum(dr(1:3)**2))
            !---
            !c          call calc_phiaug(dr,augregion(4,i1,i2,i3),
            !c     &         phiaugtmp(i1,i2,i3),isp,iq,iband)
            ! x=r*sin(theta)*cos(pphi)
            ! y=r*sin(theta)*sin(pphi)
            ! z=r*cos(theta)
            theta    = acos(dr(3)/(drlength+1.0d-15))
            sintheta = sqrt(1.0d0-cos(theta)**2)
            pphi     = acos(dr(1)/(drlength*sintheta+1.0d-15))
            if (dr(2) < 0.0d0) pphi=2*pi-pphi
            do il=0,lmax
              call calc_Ylm(il,theta,pphi, &
                   Y(1:2*il+1,il+1), &
                   Yreal(1:2*il+1,il+1))
            enddo
            !          phia=0.0d0
            do idim=1,ldim2
              if (mnla(4,idim) == ibas) then
                mtmp=mnla(1,idim)
                ntmp=mnla(2,idim)
                ltmp=mnla(3,idim)
                if (ltmp > lmax) then
                  stop 'ltmp.gt.lmax!'
                endif
                phiaugtmp(i1,i2,i3)=phiaugtmp(i1,i2,i3) +cphi(idim) &
                     *calc_gxr(drlength,ltmp,ntmp,iclass(ibas),isp) &
                     *Yreal(mtmp+ltmp+1,ltmp+1)
              endif
            enddo
          endif
        enddo
      enddo
    enddo
    if(debug) write(6,*)'---- end of calc_phi_abc2 ----------'
  end subroutine calc_phi_abc2
  subroutine expand_mesh(plat, &
       nq_wfn,nband_wfn,q_wfn,nsp, &
       rini,rfin, &
       mesh0,phipw0,phiaug0,phitot0, &
       mesh, phipw,phiaug,phitot )

    implicit none
    integer,intent(in):: nq_wfn,nband_wfn,nsp
    real(8),intent(in)::  q_wfn(3,nq_wfn)
    integer:: rini(3), rfin(3)
    integer:: mesh0(3),mesh(3)
    double complex:: phipw0(mesh0(1)+1,mesh0(2)+1,mesh0(3)+1, &
         nband_wfn,nq_wfn,nsp)
    double complex:: phiaug0(mesh0(1)+1,mesh0(2)+1,mesh0(3)+1, &
         nband_wfn,nq_wfn,nsp)
    double complex:: phitot0(mesh0(1)+1,mesh0(2)+1,mesh0(3)+1, &
         nband_wfn,nq_wfn,nsp)

    double complex::phipw(mesh(1)+1,mesh(2)+1,mesh(3)+1, &
         nband_wfn,nq_wfn,nsp)
    double complex::phiaug(mesh(1)+1,mesh(2)+1,mesh(3)+1, &
         nband_wfn,nq_wfn,nsp)
    double complex::phitot(mesh(1)+1,mesh(2)+1,mesh(3)+1, &
         nband_wfn,nq_wfn,nsp)

    real(8):: pi ,cell_l(3),phase
    integer:: i,j,k,iq
    complex(8):: ci=(0.0d0,1.0d0), expikt
    real(8):: sss(3) ,plat(3,3)!
    pi=2.0d0*asin(1.0d0)

    !  phase factor, pi*sum_i(q_wfn(i)*cell_lbound(i))

    iqloop: do iq=1,nq_wfn

      do i=rini(1),rfin(1)-1
        do j=rini(2),rfin(2)-1
          do k=rini(3),rfin(3)-1
            !        cell_l=(/i,j,k/)
            sss=plat(:,1)*i+plat(:,2)*j+plat(:,3)*k !
            ! similar to exp (ikT)
            phase=sum(q_wfn(1:3,iq)*sss) !bugfix jul2014 cell_l(1:3) )
            expikt=exp(ci*2.0d0*pi*phase )  ! phase of the cell
            ! ugfix july2014 mesh0* 1+(k-rini(1))*mesh0(1) ---> 1+(k-rini(3))*mesh0(3)
            phipw(1+(i-rini(1))*mesh0(1):1+(1+i-rini(1))*mesh0(1), &
                 1+(j-rini(2))*mesh0(2):1+(1+j-rini(2))*mesh0(2), &
                 1+(k-rini(3))*mesh0(3):1+(1+k-rini(3))*mesh0(3),:,iq,:) = &
                 expikt*  phipw0(1:1+mesh0(1), &
                 1:1+mesh0(2), &
                 1:1+mesh0(3),:,iq,:)
            phiaug(1+(i-rini(1))*mesh0(1):1+(1+i-rini(1))*mesh0(1), &
                 1+(j-rini(2))*mesh0(2):1+(1+j-rini(2))*mesh0(2), &
                 1+(k-rini(3))*mesh0(3):1+(1+k-rini(3))*mesh0(3),:,iq,:) = &
                 expikt*  phiaug0(1:1+mesh0(1), &
                 1:1+mesh0(2), &
                 1:1+mesh0(3),:,iq,:)
            phitot(1+(i-rini(1))*mesh0(1):1+(1+i-rini(1))*mesh0(1), &
                 1+(j-rini(2))*mesh0(2):1+(1+j-rini(2))*mesh0(2), &
                 1+(k-rini(3))*mesh0(3):1+(1+k-rini(3))*mesh0(3),:,iq,:) = &
                 expikt*  phitot0(1:1+mesh0(1), &
                 1:1+mesh0(2), &
                 1:1+mesh0(3),:,iq,:)
          enddo
        enddo
      enddo
    enddo iqloop
  end subroutine expand_mesh
  subroutine wrt_cube( &
       basename, &
       alat,plat,nsp,nq_wfn,nband_wfn,q_wfn,bindx_wfn, &
       mesh,rini,rfin,phipw,phiaug,phitot, &
       natom,apos,nclass,iclass,zz )
    implicit none
    ! input
    character(*),intent(in):: basename
    double precision,intent(in) :: alat,plat(3,3),rini(3),rfin(3)
    integer,intent(in) :: nsp,nq_wfn,nband_wfn,bindx_wfn(nband_wfn),mesh(3)
    double precision,intent(in) :: q_wfn(3,nq_wfn)
    double complex,intent(in) :: &
         phipw(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp), &
         phiaug(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp), &
         phitot(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp)

    integer,intent(in) :: natom,nclass,iclass(natom)
    double precision,intent(in) :: apos(3,natom),zz(nclass)

    character(200):: filename
    integer:: ifile !150

    integer :: isp,iq,ib,iband,i1,i2,i3,natomall
    double precision :: rtmp(3),r(3)
    double precision,parameter:: zero=0.0d0

    double precision :: Z(natom)
    integer:: ic,ia,nrange(3),idim, iimg
    write(filename,'(a,a)')  basename(:len_trim(basename)), '.xyz'
    open(newunit=ifile,file=filename)
    do ia=1,natom
      ic=iclass(ia)
      Z(ia)=zz(ic)
    enddo

    do ia=1,3
      nrange(ia)= max( -floor(rini(ia)), ceiling(rfin(ia)) )
    enddo
    write(*,*) 'nrange=',nrange(1:3)

    call wrt_pos_xyz('query natom',ifile,alat,plat, rini,rfin,natom,apos, Z, nrange, &
         natomall )
    !      call wrt_pos_2('query natom',ifile,alat,plat, rini,rfin,natom,apos, Z, nrange,
    !     o  natomall )


    do iimg=1,2
      do isp=1,nsp
        do iq=1,nq_wfn
          do ib=1,nband_wfn
            iband=bindx_wfn(ib)
            if (iimg == 1) then
              write(filename,"(a,a1,i1,a1,i4.4,a1,i4.4,a6)") &
                   basename(:len_trim(basename)), 's', &
                   isp,'q',iq,'b',iband,'r.cube'
            else
              write(filename,"(a,a1,i1,a1,i4.4,a1,i4.4,a6)") &
                   basename(:len_trim(basename)), 's', &
                   isp,'q',iq,'b',iband,'i.cube'
            endif
            write(*,*) 'open ',filename
            open(newunit=ifile,file=filename,status='unknown')
            write(ifile,'(a)') 'wavefunction'
            write(ifile,'(a,4I5)') 'isp,iq,ib,iband=',isp,iq,ib,iband

            r(:)= rini(:) *alat
            write(ifile,100) natomall, r(1:3)

            r = (rfin-rini)*alat
            idim=1
            write(ifile,100) mesh(idim)+1,r(idim)/mesh(idim),zero,zero
            idim=2
            write(ifile,100) mesh(idim)+1,zero,r(idim)/mesh(idim),zero
            idim=3
            write(ifile,100) mesh(idim)+1,zero,zero,r(idim)/mesh(idim)

            call wrt_pos_xyz('write',ifile,alat,plat, rini,rfin,natom,apos, Z, nrange, &
                 natomall )
            !      call wrt_pos_2('write',ifile,alat,plat, rini,rfin,natom,apos, Z, nrange,
            !     o  natomall )

            do i1=1,mesh(1)+1
              do i2=1,mesh(2)+1
                if (iimg == 1) then
                  write(ifile,200) &
                       (  real(phitot(i1,i2,i3,ib,iq,isp)), i3=1,mesh(3)+1)
                else
                  write(ifile,200) &
                       (  imag(phitot(i1,i2,i3,ib,iq,isp)), i3=1,mesh(3)+1)
                endif
              enddo
            enddo

            close(ifile)

          enddo ! ib
        enddo ! iq
      enddo ! isp

    enddo ! iimg

100 format(i6,4f20.10)
200 format(6f20.10)
  end  subroutine wrt_cube
  subroutine wrt_pos_xyz( job, ifile,alat,plat, rini,rfin,natom,apos, Z, nrange, natomall )
    !! atoms in range rini:rfin is given. These are in fractional coordinate.
    !! nrange is allowed range in the unit of plat
    implicit none
    integer:: ifile
    character(*):: job
    double precision,intent(in) :: alat,plat(3,3),rini(3),rfin(3), apos(3,natom),Z(natom)
    integer,intent(in) :: natom,nrange(3)
    integer,intent(out) :: natomall
    integer:: natomx
    integer :: i,i1,i2,i3
    double precision :: v1(3),v2(3),aini(3),afin(3),eps
    double precision ,allocatable :: rall(:,:),zall(:)
    integer::ix
    real(8)::qlat(3,3)
    call minv33tp(plat,qlat)
    eps = 0.05d0
    !      aini = alat*(rini-eps)
    !      afin = alat*(rfin+eps)
    write(6,*)' xxxxxxxxx nrange=',nrange
    write(6,*)' xxxxxxxxx rini=',rini
    write(6,*)' xxxxxxxxx rfin=',rfin
    natomx=natom*(2*nrange(1)+1)*(2*nrange(2)+1)*(2*nrange(3)+1)
    allocate(rall(3,natomx),zall(natomx))
    natomall = 0
    do i=1,natom
      do i1=-nrange(1),nrange(1)
        do i2=-nrange(2),nrange(2)
          do i3=-nrange(3),nrange(3)
            v1(1)=dble(i1)
            v1(2)=dble(i2)
            v1(3)=dble(i3)
            !          call mymatvec(plat,v1,v2,3,3)
            !          v2(1:3)=alat*(v2(1:3)+apos(1:3,i)) !
            do ix=1,3
              v2(ix)= v1(ix) + sum( apos(1:3,i)*qlat(1:3,ix) ) !atomic position in fractional coordinate
            enddo
            if(  rini(1)-eps<v2(1) .AND. v2(1)<rfin(1)+eps .AND. &
                 rini(2)-eps<v2(2) .AND. v2(2)<rfin(2)+eps .AND. &
                 rini(3)-eps<v2(3) .AND. v2(3)<rfin(3)+eps ) then
              !          if ( (v2(1).ge.aini(1).and.v2(1).le.afin(1))
              !     &    .and.(v2(2).ge.aini(2).and.v2(2).le.afin(2))
              !     &    .and.(v2(3).ge.aini(3).and.v2(3).le.afin(3)) ) then
              natomall = natomall + 1
              rall(1:3,natomall) = v2(1:3)
              zall(natomall) = Z(i)
            endif
          enddo
        enddo
      enddo
    enddo
    if (job == 'write' .OR. job == 'output') then
      do i=1,natomall
        write(ifile,"(i6,4E20.12)") int(zall(i)),zall(i),rall(1:3,i)
      enddo
    endif
    deallocate(rall,zall)
  end  subroutine wrt_pos_xyz
  subroutine wrt_xsf( &
       basename,vis_unit, &
       alat,plat,nsp,nq_wfn,nband_wfn,q_wfn,bindx_wfn, &
       mesh,rini,rfin,phipw,phiaug,phitot, &
       natom,apos,nclass,iclass,zz )
    implicit none
    character(*),intent(in):: basename,vis_unit
    double precision,intent(in) :: alat,plat(3,3),rini(3),rfin(3)
    integer,intent(in) :: nsp,nq_wfn,nband_wfn,bindx_wfn(nband_wfn),mesh(3)
    double precision,intent(in) :: q_wfn(3,nq_wfn)
    double complex,intent(in) :: &
         phipw(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp), &
         phiaug(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp), &
         phitot(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp)

    integer,intent(in) :: natom,nclass,iclass(natom)
    double precision,intent(in) :: apos(3,natom),zz(nclass)
    character(200):: filename
    integer:: ifile,ifile_handle
    integer :: isp,iq,ib,iband,i1,i2,i3,natomall
    double precision :: rtmp(3),r(3)
    double precision,parameter:: zero=0.0d0
    double precision :: Z(natom)
    integer:: ic,ia,nrange(3),idim, iimg,ix
    character(3)::ril
    write(6,*)'wrt_xsf: basename='//trim(basename)
    do ia=1,natom
      ic=iclass(ia)
      Z(ia)=zz(ic)
    enddo
    do ix=1,3
      write(6,*)'rini rfin=',ix,rini(ix),rfin(ix)
      nrange(ix)= max( -floor(rini(ix)), ceiling(rfin(ix)) )
    enddo
    write(*,*) 'natom nrange=',natom,nrange(1:3)
    !      do ia=1,natom
    !         write(6,*)'bbb=',ia,apos(:,ia)
    !      enddo
    call wrt_pos_xyz('query natom',ifile,alat,plat, rini,rfin,natom,apos, Z, nrange, &
         natomall )
    write(filename,'(a,a)')  basename(:len_trim(basename)), '.xsf'
    write(*,*) 'open ',filename
    !ifile=ifile_handle()
    open(newunit=ifile,file=filename,status='unknown')
    write(ifile,'(a)') '# wavefunction'
    write(ifile,'(a)') 'PRIMVEC'
    !      close(ifile)
    do i1=1,3
      write(ifile,'(3f20.10)') plat(:,i1)*alat
    enddo
    call wrt_pos_xyz('write',ifile,alat,plat, rini,rfin,natom,apos,Z,nrange, &
         natomall )
    write(ifile,'(a)') 'BEGIN_BLOCK_DATAGRID_3D'
    write(ifile,'(a)') basename
    do iimg=1,2
      if(iimg==1) ril='_Re'
      if(iimg==2) ril='_Im'
      do isp=1,nsp
        do iq=1,nq_wfn
          do ib=1,nband_wfn
            iband=bindx_wfn(ib)
            !        write(ifile,'(a,i1,a,i3.3,a,i3.3,a,i3.3,a,i1)')
            !     .    'isp',isp,'_iq',iq,'_ib',ib,'_',iband,'_ri',iimg
            write(ifile,'(a,i1,a,i3.3,a,i3.3,a,i3.3,a)') &
                 'BEGIN_DATAGRID_3D_isp',isp,'_iq',iq,'_ib',ib,'_',iband,ril
            write(ifile,'(3i5)') mesh(1:3)+1
            if (trim(vis_unit) == 'alat') then
              r(:)= rini(:) *alat
              write(ifile,'(3f20.5)') r(1:3) !origin
              write(ifile,'(3f20.5)') alat*(rfin(1)-rini(1))
              write(ifile,'(3f20.5)') alat*(rfin(2)-rini(2))
              write(ifile,'(3f20.5)') alat*(rfin(3)-rini(3))
            elseif (trim(vis_unit) == 'abc') then
              r(:) = plat(:,1)*rini(1)+ plat(:,2)*rini(2)+ plat(:,3)*rini(3)
              r= r*alat
              write(ifile,'(3f20.5)') r(1:3) !origin
              write(ifile,'(3f20.5)') alat*plat(:,1)*(rfin(1)-rini(1))
              write(ifile,'(3f20.5)') alat*plat(:,2)*(rfin(2)-rini(2))
              write(ifile,'(3f20.5)') alat*plat(:,3)*(rfin(3)-rini(3))
            endif
            do i3=1,mesh(3)+1
              do i2=1,mesh(2)+1
                if (iimg == 1) then
                  write(ifile,200) &
                       (  real(phitot(i1,i2,i3,ib,iq,isp)), i1=1,mesh(1)+1)
                else
                  write(ifile,200) &
                       (  imag(phitot(i1,i2,i3,ib,iq,isp)), i1=1,mesh(1)+1)
                endif
              enddo
            enddo
            write(ifile,'(a)') 'END_DATAGRID_3D'
          enddo ! ib
        enddo ! iq
      enddo ! isp
    enddo ! iimg
    write(ifile,'(a)') 'END_BLOCK_DATAGRID_3D'
    close(ifile)
100 format(i6,4f20.10)
200 format(6E20.10)
  end subroutine wrt_xsf

  ! !--------------------------------------------------------------------
  ! subroutine wrt_pos_xyz( job, ifile,alat,plat, rini,rfin,natom,apos, Z, nrange, natomall )
  !   !! atoms in range rini:rfin is given. These are in fractional coordinate.
  !   !! nrange is allowed range in the unit of plat
  !   implicit none
  !   integer:: ifile
  !   character(*):: job
  !   double precision,intent(in) :: alat,plat(3,3),rini(3),rfin(3), apos(3,natom),Z(natom)
  !   integer,intent(in) :: natom,nrange(3)
  !   integer,intent(out) :: natomall
  !   integer:: natomx
  !   integer :: i,i1,i2,i3
  !   double precision :: v1(3),v2(3),aini(3),afin(3),eps
  !   double precision ,allocatable :: rall(:,:),zall(:)

  !   integer::ix
  !   real(8)::qlat(3,3)

  !   call minv33tp(plat,qlat)
  !   eps = 0.05d0
  !   !      aini = alat*(rini-eps)
  !   !      afin = alat*(rfin+eps)

  !   write(6,*)' xxxxxxxxx nrange=',nrange
  !   write(6,*)' xxxxxxxxx rini=',rini
  !   write(6,*)' xxxxxxxxx rfin=',rfin

  !   natomx=natom*(2*nrange(1)+1)*(2*nrange(2)+1)*(2*nrange(3)+1)
  !   allocate(rall(3,natomx),zall(natomx))
  !   natomall = 0
  !   do i=1,natom
  !      do i1=-nrange(1),nrange(1)
  !         do i2=-nrange(2),nrange(2)
  !            do i3=-nrange(3),nrange(3)
  !               v1(1)=dble(i1)
  !               v1(2)=dble(i2)
  !               v1(3)=dble(i3)
  !               !          call mymatvec(plat,v1,v2,3,3)
  !               !          v2(1:3)=alat*(v2(1:3)+apos(1:3,i)) !
  !               do ix=1,3
  !                  v2(ix)= v1(ix) + sum( apos(1:3,i)*qlat(1:3,ix) ) !atomic position in fractional coordinate
  !               enddo
  !               if(  rini(1)-eps<v2(1) .AND. v2(1)<rfin(1)+eps .AND. &
  !                    rini(2)-eps<v2(2) .AND. v2(2)<rfin(2)+eps .AND. &
  !                    rini(3)-eps<v2(3) .AND. v2(3)<rfin(3)+eps ) then
  !                  !          if ( (v2(1).ge.aini(1).and.v2(1).le.afin(1))
  !                  !     &    .and.(v2(2).ge.aini(2).and.v2(2).le.afin(2))
  !                  !     &    .and.(v2(3).ge.aini(3).and.v2(3).le.afin(3)) ) then
  !                  natomall = natomall + 1
  !                  call mymatvec(plat,v1,v2,3,3)
  !                  rall(1:3,natomall) = alat*(v2(1:3)+apos(1:3,i))
  !                  zall(natomall) = Z(i)
  !               endif
  !            enddo
  !         enddo
  !      enddo
  !   enddo

  !   !      if (job.eq.'write' .or. job.eq.'output') then
  !   !      do i=1,natomall
  !   !         write(ifile,"(i6,4E20.12)") int(zall(i)),zall(i),rall(1:3,i)
  !   !      enddo
  !   !      endif
  !   if (job == 'write' .OR. job == 'output') then
  !      write(ifile,'(a)') 'ATOMS'
  !      do i=1,natomall
  !         write(ifile,"(i3,4F20.5)") int(zall(i)),rall(1:3,i)
  !      enddo
  !   endif
  !   deallocate(rall,zall)
  ! end  subroutine wrt_pos_xyz
  ! Simple mathematic functions
  ! ccccccccccccccccccccccccccccccccccc
  ! Factrization

  integer function myfact(n)
    implicit none
    integer,intent(in) :: n
    integer :: i
    myfact=1
    do i=1,n
      myfact=myfact*i
    enddo
  end function myfact
  ! ccccccccccccccccccccccccccccccccccc
  ! Permutation

  integer function myperm(n,m)
    implicit none
    integer,intent(in) :: n,m
    integer :: i
    myperm=1
    do i=0,m-1
      myperm=myperm*(n-i)
    enddo
  end function myperm
  ! ccccccccccccccccccccccccccccccccccc
  ! Combination

  integer function mycomb(n,m)
    implicit none
    integer,intent(in) :: n,m
    integer :: i,nmm
    !  function
    !  integer :: myperm,myfact

    nmm=m
    if (n-m < m) nmm=n-m

    mycomb=myperm(n,nmm)/myfact(nmm)
  end function mycomb
  ! cccccccccccccccccccccccccccccccccc
  ! Calculate m-th derivative of x^n
  ! (d/dx^m) x^n = n*(n-1)*...*(n-m+1)*x^(n-m)

  double precision function mydif(x,n,m)
    implicit none
    double precision,intent(in) :: x
    integer,intent(in) :: n,m
    ! function
    !  integer :: myperm

    if (m > n) then
      mydif=0
      return
    endif

    mydif=dble(myperm(n,m))*x**(n-m)
  end function mydif

  ! Simple matrix routines

  ! ccccccccccccccccccccccccc
  ! cross product V1 X V2 =>U
  subroutine mycross(V1,V2,U)
    implicit none
    double precision :: V1(3),V2(3),U(3)
    U(1)=V1(2)*V2(3)-V1(3)*V2(2)
    U(2)=V1(3)*V2(1)-V1(1)*V2(3)
    U(3)=V1(1)*V2(2)-V1(2)*V2(1)
  end subroutine mycross
  ! ccccccccccccccccccccccccc
  ! multiply matrix A(m,n) and vector V(n) => AV(m)
  subroutine mymatvec(A,V,AV,m,n)
    implicit none
    integer :: m,n
    double precision :: A(m,n),V(n),AV(m)
    ! local
    integer :: i,j

    AV(1:m)=0.0d0
    do j=1,n
      do i=1,m
        AV(i)=AV(i)+A(i,j)*V(j)
      enddo
    enddo
  end subroutine mymatvec
  ! ccccccccccccccccccccccccc
  ! multiply matrix A(l,m) and B(m,n) => C(l,n)
  subroutine mymatmat(A,B,C,l,m,n)
    implicit none
    integer :: l,m,n
    double precision :: A(l,m),B(m,n),C(l,n)
    ! local
    double precision ::  At(m,l)
    integer :: i,j,k

    !call mytranspose(A,At,l,m)
    At=transpose(A)
    C(1:l,1:n)=0.0d0

    do j=1,n
      do i=1,l
        do k=1,m
          C(i,j)=C(i,j)+At(k,i)*B(k,j)
        enddo
      enddo
    enddo
  end subroutine mymatmat
  ! cccccccccccccccccccccccc
  ! determinant of 3x3 matrix A
  subroutine mydet3(A,det3)
    implicit none
    double precision ::  A(3,3)
    double precision ::  det3

    det3=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+ &
         A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+ &
         A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
  end subroutine mydet3
  ! cccccccccccccccccccccccc
  ! subroutine mytranspose(A,At,m,n)
  !   implicit none

  !   ! input
  !   integer :: m,n
  !   double precision ::  A(m,n)
  !   ! output
  !   double precision ::  At(n,m)
  !   ! local
  !   integer :: i,j

  !   do j=1,n
  !     do i=1,m
  !       At(j,i)=A(i,j)
  !     enddo
  !   enddo
  ! end subroutine mytranspose
  ! cccccccccccccccccccccccc
  subroutine myinv3(A,invA)
    implicit none
    ! input
    double precision ::  A(3,3)
    ! output
    double precision ::  invA(3,3)
    ! local
    double precision ::  eps
    parameter (eps=1.0d-20)
    double precision ::  detA,At(3,3) ! the transpose of A
    double precision ::  at1(3),at2(3),at3(3)

    call mydet3(A,detA)
    if (abs(detA) <= eps) then
      write(*,*) 'Error in myinv3: detA<eps'
      stop 'Error in myinv3: detA<eps'
    endif
    !call mytranspose(A,At,3,3)
    At=transpose(A)
    at1(1:3)=At(1:3,1)
    at2(1:3)=At(1:3,2)
    at3(1:3)=At(1:3,3)

    call mycross(at2,at3,invA(1:3,1))
    call mycross(at3,at1,invA(1:3,2))
    call mycross(at1,at2,invA(1:3,3))

    invA(1:3,1:3)=invA(1:3,1:3)/detA

  end subroutine myinv3
  ! ccccccccccccccccccccccccc


  ! Calculate rotation matrix D_{MK}^{J}(\alpha,\beta,\gamma)
  ! alpha, beta, and gamma are Euler angles

  ! R(a,b,g) YLM = \sum_{M'} D_{M' M}^{J}(a,b,g) YLM'

  subroutine calc_D(J,alpha,beta,gamma,D)
    implicit none
    integer,intent(in) :: J
    double precision,intent(in) :: alpha,beta,gamma
    double complex :: D(2*J+1,2*J+1)
    ! local
    double precision :: DJbeta(2*J+1,2*J+1)
    integer :: i1,i2,M,K
    double complex :: I
    parameter (I=dcmplx(0.0d0,1.0d0))
    call calc_DJbeta(J,beta,DJbeta)
    do i2=1,2*J+1
      K=i2-J-1
      do i1=1,2*J+1
        M=i1-J-1
        D(i1,i2)=exp(-I*M*alpha)*DJbeta(i1,i2) &
             *exp(-I*K*gamma)
        !          write(*,"(a,i3,3f5.2,2i3,2f10.5)")
        !     &         'J,alpha,beta,gamma,M,K,D_{MK}=',
        !     &         J,alpha,beta,gamma,M,K,D(i1,i2)
      enddo
    enddo

  end subroutine calc_D

  ! Calculate rotation matrix D_{MK}^{J}(\beta)

  subroutine calc_DJbeta(J,beta,DJbeta)
    implicit none
    integer,intent(in) :: J
    double precision,intent(in) :: beta
    double precision :: DJbeta(2*J+1,2*J+1)
    double precision :: zeta
    integer :: i1,i2,M,K,l
    double precision :: tmp1,tmp2,tmp3
    ! function
    !  integer :: myfact,mycomb
    !  double precision :: mydif
    !c
    zeta=sin(0.5d0*beta)**2
    do i2=1,2*J+1
      K=i2-J-1
      do i1=1,2*J+1
        M=i1-J-1
        tmp1= myfact(J+M)
        tmp2= zeta**(M-K)*(1.0d0-zeta)**(M+K)* myfact(J-M)*dble(myfact(J+K))*dble(myfact(J-K))
        tmp3= 0d0
        do l=0,J+K
          tmp3=tmp3+(-1.0d0)**(J+K-l)*mycomb(J+K,l)* &
               mydif(zeta,2*J-l,J-M)
        enddo
        !          DJbeta(i1,i2)=(-1.0d0)**(M-K)*sqrt(tmp1/tmp2)*tmp3
        DJbeta(i1,i2)=(-1.0d0)**(M-K)*sqrt(tmp1/(tmp2+1.0d-20))*tmp3
      enddo
    enddo
  end subroutine calc_DJbeta
  ! cccccccccccccccccccc
  ! Calculate (complex) Ylm(theta,phi)
  ! theta=beta,phi=alpha
  subroutine calc_Ylm(L,beta,alpha,Y,Yreal)
    implicit none
    integer,intent(in) :: L
    double precision ,intent(in) :: beta,alpha
    double complex :: Y(2*L+1)
    double precision :: Yreal(2*L+1)
    ! local
    double precision :: pi
    double complex :: D(2*L+1,2*L+1)
    integer :: i1,m
    pi=4.0d0*atan(1.0d0)
    call calc_D(L,alpha,beta,0.0d0,D)
    do i1=1,2*L+1
      Y(i1)=sqrt(dble(2*L+1)/(4.0d0*pi))*dconjg(D(i1,L+1))
    enddo
    do m=-L,L
      if (m < 0) then
        Yreal(m+L+1)=(-1.0d0)**(-m)*sqrt(2.0d0)*dimag(Y(-m+L+1))
      elseif (m == 0) then
        Yreal(m+L+1)=dble(Y(m+L+1))
      else
        Yreal(m+L+1)=(-1.0d0)**(m)*sqrt(2.0d0)*dble(Y(m+L+1))
      endif
    enddo
  end subroutine calc_Ylm
end module m_wanutil
module m_wanplot
  use m_wanutil,only:myinv3,mymatvec,calc_ylm,calc_phiall_abc2,b2w,wrt_xsf,expand_mesh,calc_npw
  public wanplot
  private
contains
  subroutine wanplot()
    !! == Wannier function plot. Wannier function is expanded in the PW (spacial mesh). ==
    !! NOTE: Because os lazyness, not yet MPI. In cases, it may be useful...
    !! === Usage: ===
    !!  We first need to run maxloc generation as
    !!   echo 0|hbasfp0
    !!   echo 1|$nfpgw/hmaxloc   >lmaxloc1
    !!   $nfpgw/hpsig            >lpsig
    !!   echo 2|$nfpgw/huumat    >luumat2
    !!   echo 2|$nfpgw/hmaxloc   >lmaxloc2
    !!
    !! === Remarks ===
    !!   iclass means equivalent sites.
    !!   Use Xcrysden to see xsf file. Need to figure out automatic controll.
    !!   Currently a little inconvenient.
    !!   plot Wannier functions. Do this after Wannier matrix dnk is generated.
    !!
    !! ==== History ====
    !! takao modified jul2014 from calc_wannier6 by H.Kino.
    !! 080603, For plotting arbitary region, from calc_wannier
    !! 071114, Takashi Miyake, from calc_wfn.F
    !! -------------------------------------------------------------------------
    !      program calc_wannier6
    !      use m_LMTO
    !      use m_MLWF
    !      use m_wfrho_abc
    use mpi
    use m_keyvalue,only: getkeyvalue
    use m_readqg,only:   readqg,readngmx
    use m_hamindex,only: Readhamindex
    use m_readeigen,only: init_readeigen,init_readeigen2,readeval,lowesteval,readcphif,readgeigf,readcphifq 
    use m_read_bzdata,only: read_bzdata, nqbz,nqibz,nqbzw,nteti,ntetf &
         ,n1,n2,n3,qbas=>qlat,ginv,qbz,wbz,qibz,wibz,qbzw,idtetf,ib1bz,idteti &
         ,nstar,irk,nstbz,ngrp2=>ngrp !,qibz_r,nqibz_r
    use m_wan_lmf2gw,only: lmf2gw,iclass,nclass,zz,alat,nbas,nsp,plat,ldim2,bas !set_mnla,
    use m_wan_qg,only: read_qg,ngp
    use m_genallcf_v3,only: genallcf_v3, nprecb,mrecb,mrece,nqbzt,nband,mrecg
    !  use m_wanplotformat,only: wrt_cube,wrt_xsf
    !  use m_readhbe,only:Readhbe,nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg

    implicit none
    integer :: nwf,iko_ix,iko_fx,nband_wfn
    real(8) :: a,qdum(3)
    integer :: nq_wfn,tlat(3),nfac
    integer,allocatable :: bindx_wfn(:),bindx_wfn0(:)
    real(8),allocatable:: qreg(:,:),q_wfn(:,:)
    complex(8),allocatable :: &
         phipw0(:,:,:,:,:,:), phiaug0(:,:,:,:,:,:),phitot0(:,:,:,:,:,:), &
         phipw(:,:,:,:,:,:),  phiaug(:,:,:,:,:,:), phitot(:,:,:,:,:,:), &
         phipw1(:,:,:,:,:,:), phiaug1(:,:,:,:,:,:),phitot1(:,:,:,:,:,:), &
         wanpw(:,:,:,:,:),  wanaug(:,:,:,:,:),wantot(:,:,:,:,:), &
         wanpw1(:,:,:,:,:), wanaug1(:,:,:,:,:),wantot1(:,:,:,:,:), &
         dnk(:,:,:,:)
    integer:: i_rini(3),i_rfin(3)
    integer :: npw(3),mesh(3),mesh0(3),mesh1(3),meshrange(2,3)
    integer :: iq,ib,  ngpmx
    integer:: i,id,j,ifi,iqbz2!,ifile_handle
    integer :: isp,iqbz,ikp,iwf !,nprecb,mrecb,mrece,nlmtot,nqbzt,ifhbed, mrecg,nband
    real(8):: r_rini0(3),r_rfin0(3),ang,alat_ang
    real(8):: r_rini(3), r_rfin(3),r
    real(8):: r_rini1(3), r_rfin1(3)
    logical :: lrho,flag
    character(10):: vis_unit='none'
    character(20):: outputformat ='none'
    character(20)::inputfile='GWinput'
    real(8)::quu(3),det,qlat(3,3),q(3)
    complex(8),allocatable :: cphi(:,:,:,:),geig(:,:,:,:),geig2(:,:),cphi2(:,:)
    character*4::fname
    logical:: debug=.false.,vis_skip
    integer::nqbzx,incwfin

    !! MPI dummy
    !include 'mpif.h'
    integer:: ierr
    call mpi_init(ierr)
    !-----------------------------------------------------
    call getkeyvalue(inputfile,'vis_skip',  vis_skip,  default=.false. )
    if(vis_skip) call rx0s('wanplot: we found vis_skip on. Do nothing and Quit!')
    !! Readin all data
    call read_qg()
    call read_BZDATA()
    call lmf2gw()
    !      call set_mnla()
    call minv33tp (plat,qlat)  !inverse and transpose
    do isp=1,nsp
      if (isp == 1) fname='MLWU'
      if (isp == 2) fname='MLWD'
      open(newunit=ifi,file=fname,form='unformatted',status='old', action='read')
      read(ifi)nqbzx,nwf,iko_ix,iko_fx
      if(nqbz/=nqbzx)  call rx('wanplot:nqbz/=nqbzx')
      if (isp == 1) allocate(dnk(iko_ix:iko_fx,nwf,nqbz,nsp),qreg(3,nqbz))
      do iqbz = 1,nqbz
        read(ifi)iqbz2,q(1:3)
        if(debug)write(6,"(i5,3f13.5)") iqbz,q(:)
        qreg(:,iqbz) = q
        read(ifi)dnk(iko_ix:iko_fx,1:nwf,iqbz,isp)
      enddo
      close(ifi)
    enddo                     ! isp
    write(6,*)'read end of MLWU/D ...'
    !!
    incwfin= -1  !use 7th colmn for core at the end section of GWIN
    call genallcf_v3(incwfin) !in module m_genallcf_v3
    !  call Readhbe()
    call Readhamindex()
    call init_readeigen()!nband,mrece) !initialization of readEigen
    call init_readeigen2()!mrecb,ldim2,mrecg) !initialize m_readeigen

    !! replace geig and cphi with those for WF. Converted by the dnk matrix.
    !! geig2 = cphi2 = 0 for ik > nqbz
    call readngmx('QGpsi',ngpmx)
    write(6,*)'ngpmx nband ldim2=',ngpmx,nband,ldim2
    write(6,*)'nwf   nqbz  nsp  =',nwf,nqbz,nsp
    allocate(geig2(ngpmx,nband))
    allocate(cphi2(ldim2,nband))
    allocate(geig(ngpmx,nwf,nqbz,nsp))
    allocate(cphi(ldim2,nwf,nqbz,nsp))
    geig = 0d0
    cphi = 0d0
    do ikp = 1,nqbz
      do isp = 1,nsp
        geig2 = readgeigf(qreg(:,ikp),isp) !    call readgeig(qreg(:,ikp),ngpmx,isp, quu, geig2)
        cphi2 = readcphif(qreg(:,ikp),isp)
        quu   = readcphifq()
        if(sum(abs(qreg(:,ikp)-quu))>1d-6) call rx('wanplot: mmlf222eeeee')
        !! may2015 use zaxpy. This can avoid bug? when wanplot by qsub.
        do iwf = 1,nwf
          do ib  = iko_ix,iko_fx
            call zaxpy(ngp(ikp),dnk(ib,iwf,ikp,isp),geig2(1,ib),1,geig(1,iwf,ikp,isp),1)
            call zaxpy(ldim2,   dnk(ib,iwf,ikp,isp),cphi2(1,ib),1,cphi(1,iwf,ikp,isp),1)
            !           geig(:,iwf,ikp,isp) = geig(:,iwf,ikp,isp) + geig2(:,ib)*dnk(ib,iwf,ikp,isp)
            !           cphi(:,iwf,ikp,isp) = cphi(:,iwf,ikp,isp) + cphi2(:,ib)*dnk(ib,iwf,ikp,isp)
          enddo ! ib
        enddo ! iwf
      enddo ! isp
    enddo ! ikp
    deallocate(geig2,cphi2,dnk)
    write(6,*) '### ib,bas(1:3,ib) ############'
    do ib=1,nbas
      write(*,"(i5,3f12.6)")ib,bas(1:3,ib)
    enddo

    !! NOTE: nq_wfn = nqbz
    nq_wfn = nqbz
    allocate(q_wfn(3,nq_wfn))
    q_wfn(1:3,1:nqbz) = qbz(1:3,1:nqbz)

    !! == Readin vis_* settings in GWinput
    call getkeyvalue(inputfile,'vis_wan_band_n',nband_wfn,default=nwf)
    if (nband_wfn > nband) call rx('wanplot: nband_wfn > nband !')
    write(6,"(a,2i5)") '### nq_wfn, nband_wfn =',nq_wfn,nband_wfn
    allocate(bindx_wfn0(nband_wfn),bindx_wfn(nband_wfn))
    do ib=1,nband_wfn
      bindx_wfn0(ib) = ib
    enddo
    call getkeyvalue(inputfile,'vis_wan_band_id',bindx_wfn,default=bindx_wfn0,size=nband_wfn)
    do ib=1,nband_wfn
      write(*,"(a,2i5)") 'ib bndinx=',ib,bindx_wfn(ib)
    enddo
    !!
    call getkeyvalue(inputfile, 'vis_wan_tvec',   tlat,size=3,default=(/0,0,0/))
    write(*,"(a,3i5)")'### tlat',tlat
    call getkeyvalue(inputfile, 'vis_wan_interpolation', nfac,default=1)  !FFT
    vis_unit='abc'
    write(6,*)' CAUTION: range of ubond and lbound are in abc(cell) unit'
    !      call getkeyvalue(inputfile, 'vis_wan_unit', vis_unit)
    !      if ( trim(vis_unit).ne.'abc') then
    !         write(*,*) 'support only vis.wan.unit=abc'
    !         stop 'support only vis.wan.unit=abc'
    !      endif
    call calc_npw(nfac,npw)
    call getkeyvalue(inputfile,'vis_wan_mesh',  mesh0,  size=3, default=(/10,10,10/) )
    ! esh size 0:mesh0(1),...
    call getkeyvalue(inputfile,'vis_wan_lbound',r_rini0,size=3,default=(/-1d0,-1d0,-1d0/))!lower bound
    call getkeyvalue(inputfile,'vis_wan_ubound',r_rfin0,size=3,default=(/1d0,1d0,1d0/))   !upper bound
    write(*,*) ' mesh=',mesh0
    write(*,*) ' lbound=',r_rini0
    write(*,*) ' ubound=',r_rfin0
    !      call getkeyvalue(inputfile,'vis_wan_outputformat',outputformat,default='xsf')

    do i=1,3
      i_rini(i)= floor(r_rini0(i))
      i_rfin(i)= ceiling(r_rfin0(i))
    enddo
    !! for plot mesh
    r_rini= i_rini
    r_rfin = i_rfin
    do i=1,3
      mesh(i)= (i_rfin(i)-i_rini(i))*mesh0(i)
    enddo
    write(6,*)'mmm: i_rfin=',i_rini
    write(6,*)'mmm: i_rini=',i_rfin
    write(6,*)'mmm:   mesh=',mesh

    allocate(phipw0 (mesh0(1)+1,mesh0(2)+1,mesh0(3)+1,   nband_wfn,nq_wfn,nsp))
    allocate(phiaug0(mesh0(1)+1,mesh0(2)+1,mesh0(3)+1,  nband_wfn,nq_wfn,nsp))
    allocate(phitot0(mesh0(1)+1,mesh0(2)+1,mesh0(3)+1,  nband_wfn,nq_wfn,nsp))

    !! == Generate phi, which is the real-space rep. of the Bloch functions (on mesh0, real mesh points).
    !! Time consuming part.
    call calc_phiall_abc2(nq_wfn,nband_wfn,q_wfn,bindx_wfn, &
         npw,mesh0,nsp,nband,ldim2,ngpmx, &
         geig,cphi,nwf, &
         phipw0,phiaug0,phitot0)
    write(6,*)'sumchk 000 =',sum(abs(phipw0)),sum(abs(phiaug0)),sum(abs(phitot0))

    allocate(phipw(mesh(1)+1,mesh(2)+1,mesh(3)+1,   nband_wfn,nq_wfn,nsp))
    allocate(phiaug(mesh(1)+1,mesh(2)+1,mesh(3)+1,  nband_wfn,nq_wfn,nsp))
    allocate(phitot(mesh(1)+1,mesh(2)+1,mesh(3)+1,  nband_wfn,nq_wfn,nsp))

    !! phipw,phiaug,phitot are the real-space rep. on extended mesh points (mesh).
    call expand_mesh(plat, &
         nq_wfn,nband_wfn,q_wfn,nsp, &
         i_rini,i_rfin, &
         mesh0, phipw0,phiaug0,phitot0, &
         mesh, phipw,phiaug,phitot )
    !      write(6,*)'sumchk 22222=',sum(abs(phipw)),sum(abs(phiaug)),sum(abs(phitot))

    !! from Bloch to Wannier
    allocate(wanpw(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nsp) )
    allocate(wanaug(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nsp))
    allocate(wantot(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nsp))

    call b2w(nq_wfn,nband_wfn,q_wfn,bindx_wfn,tlat,mesh,phipw,       wanpw )
    call b2w(nq_wfn,nband_wfn,q_wfn,bindx_wfn,tlat,mesh,phiaug,      wanaug)
    call b2w(nq_wfn,nband_wfn,q_wfn,bindx_wfn,tlat,mesh,phitot,      wantot)
    write(*,*) 'mesh in b2w'
    write(*,*) 'mesh=',mesh
    write(*,*) 'r_rini=',r_rini
    write(*,*) 'r_rfin=',r_rfin
    write(*,*)' '

    !!..... calculate range  from rini0 and rfin0 !abc (fractional basis only)
    !!    r = rini + (rfin-rini)*(i-1)/mesh
    meshrange=0
    r_rini1=0
    r_rfin1=0
    do id=1,3
      do i=1,mesh(id)+1
        r =  r_rini(id) + (r_rfin(id)-r_rini(id))*(i-1)/mesh(id)
        if ( r> r_rini0(id) ) then
          j=i-1
          meshrange(1,id)=j
          r_rini1(id) = r_rini(id) + (r_rfin(id)-r_rini(id))*(j-1)/mesh(id)
          exit
        endif
      enddo
      do i=1,mesh(id)+1
        r =  r_rini(id) + (r_rfin(id)-r_rini(id))*(i-1)/mesh(id)
        if ( r> r_rfin0(id) .OR. i==mesh(id)+1) then
          j=i
          meshrange(2,id)=j
          r_rfin1(id) = r_rini(id) + (r_rfin(id)-r_rini(id))*(j-1)/mesh(id)
          exit
        endif
      enddo
    enddo
    print *,'meshrange2=',meshrange(2,:) !upper limits
    print *,'meshrange1=',meshrange(1,:) !lower limits

    mesh1(:)=meshrange(2,:)-meshrange(1,:)
    allocate(wanpw1 (mesh1(1)+1,mesh1(2)+1,mesh1(3)+1,nband_wfn,nsp))
    allocate(wanaug1(mesh1(1)+1,mesh1(2)+1,mesh1(3)+1,nband_wfn,nsp))
    allocate(wantot1(mesh1(1)+1,mesh1(2)+1,mesh1(3)+1,nband_wfn,nsp))
    allocate(phipw1 (mesh1(1)+1,mesh1(2)+1,mesh1(3)+1,nband_wfn,nq_wfn,nsp))
    allocate(phiaug1(mesh1(1)+1,mesh1(2)+1,mesh1(3)+1,nband_wfn,nq_wfn,nsp))
    allocate(phitot1(mesh1(1)+1,mesh1(2)+1,mesh1(3)+1,nband_wfn,nq_wfn,nsp))

    write(*,*)'range in inputfile'
    write(*,*) 'rini0=',r_rini0
    write(*,*) 'rfin0=',r_rfin0
    write(*,*)' '
    write(*,*)'cutted mesh'
    write(*,*) 'mesh=',mesh1
    write(*,*) 'rini=',r_rini1
    write(*,*) 'rfin=',r_rfin1
    write(*,*)' '

    wanpw1(1:mesh1(1)+1,1:mesh1(2)+1,1:mesh1(3)+1, :,:) = &
         wanpw(meshrange(1,1):meshrange(1,1)+mesh1(1), &
         meshrange(1,2):meshrange(1,2)+mesh1(2), &
         meshrange(1,3):meshrange(1,3)+mesh1(3),:,: )
    wanaug1(1:mesh1(1)+1,1:mesh1(2)+1,1:mesh1(3)+1, :,:) = &
         wanaug(meshrange(1,1):meshrange(1,1)+mesh1(1), &
         meshrange(1,2):meshrange(1,2)+mesh1(2), &
         meshrange(1,3):meshrange(1,3)+mesh1(3),:,: )
    wantot1(1:mesh1(1)+1,1:mesh1(2)+1,1:mesh1(3)+1, :,:) = &
         wantot(meshrange(1,1):meshrange(1,1)+mesh1(1), &
         meshrange(1,2):meshrange(1,2)+mesh1(2), &
         meshrange(1,3):meshrange(1,3)+mesh1(3),:,: )

    phipw1(1:mesh1(1)+1,1:mesh1(2)+1,1:mesh1(3)+1, :,:,:) = &
         phipw(meshrange(1,1):meshrange(1,1)+mesh1(1), &
         meshrange(1,2):meshrange(1,2)+mesh1(2), &
         meshrange(1,3):meshrange(1,3)+mesh1(3),:,:,: )
    phiaug1(1:mesh1(1)+1,1:mesh1(2)+1,1:mesh1(3)+1, :,:,:) = &
         phiaug(meshrange(1,1):meshrange(1,1)+mesh1(1), &
         meshrange(1,2):meshrange(1,2)+mesh1(2), &
         meshrange(1,3):meshrange(1,3)+mesh1(3),:,:,: )
    phitot1(1:mesh1(1)+1,1:mesh1(2)+1,1:mesh1(3)+1, :,:,:) = &
         phitot(meshrange(1,1):meshrange(1,1)+mesh1(1), &
         meshrange(1,2):meshrange(1,2)+mesh1(2), &
         meshrange(1,3):meshrange(1,3)+mesh1(3),:,:,: )
    ! rini -> r_rini1
    ! rfin -> r_rfin1
    ! mesh -> mesh1

    !! Dump phi and wannier functions -------------------
    write(*,*) '-- dump phi(bloch funciton) and wan(Wannier funciton) --'
    qdum = 0.0d0
    ang = 0.529177d0
    alat_ang=alat*ang
    write(6,*) 'Writing xsf (Xcrysden) file...'
    call wrt_xsf( &
         'wan',vis_unit, &
         alat_ang,plat,nsp,1,nband_wfn,q_wfn,bindx_wfn, &
         mesh1,r_rini1,r_rfin1,wanpw1,wanaug1,wantot1,  nbas,bas,nclass,iclass,zz )
    call wrt_xsf( &
         'phi',vis_unit, &
         alat_ang,plat,nsp,1,nband_wfn,q_wfn,bindx_wfn, &
         mesh1,r_rini1,r_rfin1,phipw1,phiaug1,phitot1,  nbas,bas,nclass,iclass,zz )
    !$$$      endif

    ! dump rho or not ----------
    lrho=.false.
    write(6,"(a,l)") 'dump rho? [T/F] (if needed set lrho=T in wanplot.F)=',lrho
    if (lrho) call calc_rho_2(alat_ang,nq_wfn,nband_wfn,mesh,r_rini,r_rfin,wanpw,wanaug,wantot)
    !      if (lrho) call calc_rho(nq_wfn,nband_wfn,npw,phipw,phiaug,phitot)
    !      if (lrho) call calc_rho(1,nband_wfn,npw,wanpw,wanaug,wantot)
    !      call cputid(0)
    call rx0s('wanplot: ok')
  end subroutine wanplot
endmodule m_wanplot
