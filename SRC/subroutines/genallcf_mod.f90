!> Structure data loaded from __MTOindex (lmf intermediate file)
!  Read once at the start of any GW binary to populate crystal/basis metadata.
module m_struct_from_lmf
  use m_lgunit, only: stdo
  use m_mpi,    only: ipr
  implicit none
  public :: struct_from_lmf_init
  integer, protected, public :: natom, nspin, lmxax, nnv, nnc, nrx, nspc, nlmto
  integer, protected, public :: nprecb, mrecb, mrece, nqbzt, nband, mrecg
  real(8), protected, public :: alat, plat(3,3), qval, tpioa
  real(8), allocatable, protected, public :: pos(:,:), z(:)
  integer, allocatable, protected, public :: lmxa(:)
  character(8), allocatable, protected, public :: spid(:)
  logical, protected, public :: laf
  integer, allocatable, protected, public :: ibasf(:)
  integer, protected, public :: nl
  integer, protected :: ndima_lmf  ! read from __MTOindex; canonical ndima is recomputed in m_gw_product_basis
  private
contains
  subroutine struct_from_lmf_init()
    integer :: ifi
    real(8), parameter :: pi = 4d0*datan(1d0)
    open(newunit=ifi, file='__MTOindex', form='unformatted')
    read(ifi) natom, alat, plat, nspin, lmxax, nnv, nnc, nrx, qval, nspc, nlmto
    allocate(pos(3,natom), z(natom), spid(natom), ibasf(natom), lmxa(natom))
    read(ifi) pos, z(1:natom), spid(1:natom), lmxa(1:natom)
    read(ifi) nprecb, mrecb, mrece, ndima_lmf, nqbzt, nband, mrecg
    read(ifi) laf, ibasf
    close(ifi)
    nl    = lmxax + 1
    tpioa = 2d0 * pi / alat
  end subroutine struct_from_lmf_init
end module m_struct_from_lmf

!> GW user-config (frequency-mesh related) read from GWinput keys
module m_gw_user_config
  use m_lgunit, only: stdo
  use m_mpi,    only: ipr
  implicit none
  public :: gw_user_config_init, set_esmr
  integer, protected, public :: niw
  real(8), protected, public :: deltaw, esmr
  private
contains
  subroutine gw_user_config_init()
    !> Read niw / deltaw / esmr. If GWinput.toml exists, load via m_GWinput
    !  (single source of truth); else fall back to legacy getkeyvalue.
    use m_keyvalue, only: getkeyvalue
    use m_GWinput,  only: gwinput_load, gwinput_loaded, &
                          mg_niw    => niw, &
                          mg_deltaw => deltaw, &
                          mg_esmr   => esmr
    logical :: have_toml
    character(len=:), allocatable :: errmsg
    inquire(file='GWinput.toml', exist=have_toml)
    if (have_toml) then
       call gwinput_load(error=errmsg)
       if (gwinput_loaded) then
          niw    = mg_niw
          deltaw = mg_deltaw
          esmr   = mg_esmr
       else
          call getkeyvalue("GWinput", "niw",    niw)
          call getkeyvalue("GWinput", "deltaw", deltaw)
          call getkeyvalue("GWinput", "esmr",   esmr)
       endif
    else
       call getkeyvalue("GWinput", "niw",    niw)
       call getkeyvalue("GWinput", "deltaw", deltaw)
       call getkeyvalue("GWinput", "esmr",   esmr)
    endif
    if(ipr) write(stdo,*) ' --- Freq ---'
    if(ipr) write(stdo,"(a,i6)")   '    niw  =', niw
    if(ipr) write(stdo,"(a,f12.6)")'    esmr =', esmr
  end subroutine gw_user_config_init
  subroutine set_esmr(esmr_in)
    real(8), intent(in) :: esmr_in
    esmr = esmr_in
  end subroutine set_esmr
end module m_gw_user_config

!> GW product basis: read PRODUCT_BASIS section + derived core/MTO indices
!  Original idea of product basis is from F.Aryasetiawan.
module m_gw_product_basis
  use m_lgunit, only: stdo
  use m_mpi,    only: ipr
  use m_struct_from_lmf, only: natom, nl, nnv, nnc, lmxa
  implicit none
  public :: gw_product_basis_init
  integer, protected, public :: lcutmx
  integer, protected, public :: nn, nlnx, nlnxv, nlnxc, nlnmx, nlnmxv, nlnmxc
  integer, protected, public :: ndima, ndimanspc, nspx
  real(8), allocatable, protected, public :: cutbase(:)
  integer, allocatable, protected, public :: nindx(:,:), il(:,:), in(:,:), im(:,:)
  integer, allocatable, protected, public :: nocc(:,:,:), nunocc(:,:,:)
  integer, allocatable, protected, public :: nindxc(:,:), lcutmxa(:)
  integer, allocatable, protected, public :: nlnm(:), nlnmv(:), nlnmc(:)
  integer, allocatable, protected, public :: ncwf(:,:,:)  ! exposed for m_core_state
  private
contains
  subroutine gw_product_basis_init(incwfx)
    use m_keyvalue, only: getkeyvalue
    use m_GWinput, only: gwinput_init, gwinput_loaded, &
                         tg_pb_tolerance => pb_tolerance, tg_pb_lcutmx => pb_lcutmx, &
                         tg_pb_nlx => pb_nlx, tg_pb_n_nlx => pb_n_nlx, &
                         tg_pb_valence => pb_valence, tg_pb_n_val => pb_n_val, &
                         tg_pb_core => pb_core, tg_pb_n_core => pb_n_core
    use m_struct_from_lmf, only: nspin, nspc
    integer, intent(in) :: incwfx
    integer :: ifi, ret, ix, ixoff, lx, iatom, iatomt, lt, n, nt, ind, ncorex, l, m, lm, nlx
    integer :: nval, ncore_local
    integer :: nnn1(natom), nnn2(natom), nnn3(natom), nnn4(natom), nnn5(natom), nnn6(natom)
    integer, allocatable :: nindxv(:,:), occv(:,:,:), unoccv(:,:,:), occc(:,:,:), unoccc(:,:,:), ncwf2(:,:,:)
    character(1000) :: tolchar
    logical :: readon
    integer :: ntol

    !---- ReadProductBasis ----
    allocate(nindxv(nl,natom), nindxc(nl,natom), &
         occv(nl,nnv,natom), unoccv(nl,nnv,natom), &
         occc(nl,nnc,natom), unoccc(nl,nnc,natom), source=0)
    allocate(ncwf2(nl,nnc,natom), ncwf(nl,nnc,natom), source=0)
    allocate(cutbase(0:2*(nl-1)), source=0d0)
    ncwf  = 99
    ncwf2 = 99
    if(ipr) write(stdo,*) ' reading <PRODUCT_BASIS> section'
    call gwinput_init()
    if (gwinput_loaded) then
       ! Tolerance: pad/truncate to 0..2*(nl-1)
       ntol = 0
       if (allocated(tg_pb_tolerance)) ntol = size(tg_pb_tolerance)
       do lx = 0, 2*(nl-1)
          if (lx + 1 <= ntol) then
             cutbase(lx) = tg_pb_tolerance(lx+1)
          else if (ntol > 0) then
             cutbase(lx) = tg_pb_tolerance(ntol)
          endif
       enddo
       do lx = 0, 2*(nl-1)
          if(ipr) write(stdo,"(' lx=',i3,' readin tolerance=',d11.3)") lx, cutbase(lx)
       enddo
       allocate(lcutmxa(natom))
       if (allocated(tg_pb_lcutmx) .and. size(tg_pb_lcutmx) >= natom) then
          lcutmxa(1:natom) = tg_pb_lcutmx(1:natom)
       else if (allocated(tg_pb_lcutmx) .and. size(tg_pb_lcutmx) > 0) then
          lcutmxa(:) = tg_pb_lcutmx(1)
       else
          lcutmxa(:) = 0
       endif
       lcutmx = lcutmxa(1)
       if(ipr) write(stdo,'(20i3)') lcutmxa(1:natom)
       ! nlx: [iatom, l, nnvv, nnc]
       do ix = 1, tg_pb_n_nlx
          iatomt = tg_pb_nlx(1, ix); lt = tg_pb_nlx(2, ix)
          if (iatomt < 1 .or. iatomt > natom) cycle
          if (lt+1 < 1 .or. lt+1 > nl) cycle
          nindxv(lt+1, iatomt) = tg_pb_nlx(3, ix)
          nindxc(lt+1, iatomt) = tg_pb_nlx(4, ix)
          if(ipr) write(stdo,*) iatomt, lt, nindxv(lt+1,iatomt), nindxc(lt+1,iatomt)
       enddo
       ! valence: [iatom, l, n, occ, unocc] -- per-(iatom,l,n) row
       if(ipr) write(stdo,*) ' --- valence product basis section'
       occv   = 0
       unoccv = 0
       do ix = 1, tg_pb_n_val
          iatomt = tg_pb_valence(1, ix); lt = tg_pb_valence(2, ix); nt = tg_pb_valence(3, ix)
          if (iatomt < 1 .or. iatomt > natom) cycle
          if (lt+1 < 1 .or. lt+1 > nl) cycle
          if (nt < 1 .or. nt > nnv) cycle
          occv  (lt+1, nt, iatomt) = tg_pb_valence(4, ix)
          unoccv(lt+1, nt, iatomt) = tg_pb_valence(5, ix)
          if(ipr) write(stdo,"(100i3)") iatomt, lt, nt, occv(lt+1,nt,iatomt), unoccv(lt+1,nt,iatomt)
       enddo
       ! core: [iatom, l, n, occ, unocc, forX0, forSxc]
       if(ipr) write(stdo,*) ' --- core product basis section'
       do ix = 1, tg_pb_n_core
          iatomt = tg_pb_core(1, ix); lt = tg_pb_core(2, ix); nt = tg_pb_core(3, ix)
          if (iatomt < 1 .or. iatomt > natom) cycle
          if (lt+1 < 1 .or. lt+1 > nl) cycle
          if (nt < 1 .or. nt > nnc) cycle
          occc  (lt+1, nt, iatomt) = tg_pb_core(4, ix)
          unoccc(lt+1, nt, iatomt) = tg_pb_core(5, ix)
          ncwf  (lt+1, nt, iatomt) = tg_pb_core(6, ix)
          ncwf2 (lt+1, nt, iatomt) = tg_pb_core(7, ix)
          if(ipr) write(stdo,"(100i3)") iatomt, lt, nt, occc(lt+1,nt,iatomt), unoccc(lt+1,nt,iatomt), &
                                          ncwf(lt+1,nt,iatomt), ncwf2(lt+1,nt,iatomt)
       enddo
    else
       call getkeyvalue("GWinput","<PRODUCT_BASIS>", unit=ifi, status=ret)
       read(ifi,*)
       read(ifi,"(a)") tolchar
       readon = .false.
       lx = 0
       do ix = 1, 1000
          if(.NOT. readon .AND. tolchar(ix:ix) /= ' ') then
             readon = .true.
             ixoff  = ix
          endif
          if(readon .AND. tolchar(ix:ix) == ' ') then
             read(tolchar(ixoff:ix), *, err=1097) cutbase(lx)
             if(lx == 2*(nl-1)) goto 1098
             readon = .false.
             lx = lx + 1
          endif
       enddo
1097   continue
       cutbase(lx:) = cutbase(lx-1)
1098   continue
       do lx = 0, 2*(nl-1)
          if(ipr) write(stdo,"(' lx=',i3,' readin tolerance=',d11.3)") lx, cutbase(lx)
       enddo
       read(ifi,*)
       allocate(lcutmxa(natom))
       read(ifi,*) lcutmxa(1:natom)
       lcutmx = lcutmxa(1)
       if(ipr) write(stdo,'(20i3)') lcutmxa(1:natom)
       if(ipr) write(stdo,"(' --- prod section: lcutmx cutbase='i3,100d11.3)") lcutmx, cutbase
       read(ifi,*)
       do iatom = 1, natom
          do l = 0, lmxa(iatom)
             read(ifi,*) iatomt, lt, nindxv(l+1,iatom), nindxc(l+1,iatom)
             if(ipr) write(stdo,*) iatomt, lt, nindxv(l+1,iatom), nindxc(l+1,iatom)
          enddo
       enddo
       if(ipr) write(stdo,*) ' --- valence product basis section'
       occv   = 0
       unoccv = 0
       read(ifi,*)
       do iatom = 1, natom
          do l = 0, lmxa(iatom)
             do n = 1, nindxv(l+1,iatom)
                read(ifi,*)                 iatomt, lt, nt, occv(l+1,n,iatom), unoccv(l+1,n,iatom)
                if(ipr) write(stdo,"(100i3)") iatomt, lt, nt, occv(l+1,n,iatom), unoccv(l+1,n,iatom)
             enddo
          enddo
       enddo
       if(ipr) write(stdo,*) ' --- core product basis section'
       read(ifi,*)
       do iatom = 1, natom
          do l = 0, lmxa(iatom)
             do n = 1, nindxc(l+1,iatom)
                read(ifi,*)                 iatomt, lt, nt, occc(l+1,n,iatom), unoccc(l+1,n,iatom), ncwf(l+1,n,iatom), ncwf2(l+1,n,iatom)
                if(ipr) write(stdo,"(100i3)") iatomt, lt, nt, occc(l+1,n,iatom), unoccc(l+1,n,iatom), ncwf(l+1,n,iatom), ncwf2(l+1,n,iatom)
             enddo
          enddo
       enddo
       close(ifi)
    endif
    if(incwfx == -1) then
       if(ipr) write(stdo,*) ' ### incwf=-1 Use ForSxc for core'
       ncwf = ncwf2
    elseif(incwfx == -2) then
       if(ipr) write(stdo,*) ' ### incwf=-2 Use NOT(ForSxc) for core and Pro-basis '
       ncwf   = merge(1-ncwf2, ncwf2, ncwf2==0 .or. ncwf2==1)
       occc   = ncwf
       unoccc = 0
       unoccv = merge(1, unoccv, occv==1 .or. unoccv==1)
    elseif(incwfx == -3) then
       if(ipr) write(stdo,*) ' ### incwf=-3  occ=1 unocc=0 incwf=1 for all core '
       occc   = 1
       ncwf   = 1
       unoccc = 0
    elseif(incwfx == -4) then
       if(ipr) write(stdo,*) ' ### incwf=-4  occ=0 and unocc=0 for all core '
       occc   = 0
       unoccc = 0
       ncwf   = 0
    elseif(incwfx == 0) then
       if(ipr) write(stdo,*) ' ### Use unocc occ ForX0 for core'
    else
       call rx(' ### proper incwf is not given for genallcf_v3:rgwinf ')
    endif
    deallocate(ncwf2)

    !---- indexcoremto ----
    ndima = 0
    do iatom = 1, natom
       ndima = ndima + sum([((2*l+1)*nindxv(l+1,iatom), l=0, lmxa(iatom))])
    enddo
    nn = maxval(nindxv(1:nl,1:natom) + nindxc(1:nl,1:natom))
    allocate(nindx(nl,natom), nocc(nl,nn,natom), nunocc(nl,nn,natom), source=0)
    do iatom = 1, natom
       do l = 0, lmxa(iatom)
          ncore_local = nindxc(l+1,iatom)
          nval        = nindxv(l+1,iatom)
          nindx(l+1,iatom)     = ncore_local + nval
          nocc(l+1,1:,iatom)   = [(occc(l+1,n,iatom),  n=1,ncore_local), (occv(l+1,n,iatom),  n=1,nval)]
          nunocc(l+1,1:,iatom) = [(unoccc(l+1,n,iatom),n=1,ncore_local), (unoccv(l+1,n,iatom),n=1,nval)]
       enddo
    enddo
    do iatom = 1, natom
       nlx         = lmxa(iatom) + 1
       nnn1(iatom) = sum(nindx(1:nlx,iatom))
       nnn2(iatom) = sum([(nindx(l+1,iatom)*(2*l+1),  l=0, nlx-1)])
       nnn3(iatom) = sum(nindxv(1:nlx,iatom))
       nnn4(iatom) = sum([(nindxv(l+1,iatom)*(2*l+1), l=0, nlx-1)])
       nnn5(iatom) = sum(nindxc(1:nlx,iatom))
       nnn6(iatom) = sum([(nindxc(l+1,iatom)*(2*l+1), l=0, nlx-1)])
    enddo
    nlnx   = maxval(nnn1)
    nlnmx  = maxval(nnn2)
    nlnxv  = maxval(nnn3)
    nlnmxv = maxval(nnn4)
    nlnxc  = maxval(nnn5)
    nlnmxc = maxval(nnn6)

    allocate(il(nlnmx,natom), in(nlnmx,natom), im(nlnmx,natom))
    do iatom = 1, natom
       ind = 0
       do l = 0, lmxa(iatom)        ! core
          do n = 1, nindxc(l+1,iatom)
             do m = 1, 2*l+1
                ind = ind + 1
                lm  = l**2 + m
                il(ind,iatom) = l
                in(ind,iatom) = n
                im(ind,iatom) = m - l - 1
             enddo
          enddo
       enddo
       do l = 0, lmxa(iatom)        ! valence
          ncorex = nindxc(l+1,iatom)
          do n = 1, nindxv(l+1,iatom)
             do m = 1, 2*l+1
                ind = ind + 1
                lm  = l**2 + m
                il(ind,iatom) = l
                in(ind,iatom) = ncorex + n
                im(ind,iatom) = m - l - 1
             enddo
          enddo
       enddo
    enddo
    allocate(nlnmv(natom), nlnmc(natom), nlnm(natom))
    do iatom = 1, natom
       nlx          = lmxa(iatom) + 1
       nlnmv(iatom) = sum([(nindxv(l+1,iatom)*(2*l+1), l=0, nlx-1)])
       nlnmc(iatom) = sum([(nindxc(l+1,iatom)*(2*l+1), l=0, nlx-1)])
       nlnm(iatom)  = sum([(nindx(l+1,iatom)*(2*l+1),  l=0, nlx-1)])
    enddo

    ndimanspc = ndima * nspc
    nspx      = nspin / nspc

    deallocate(nindxv, occv, unoccv, occc, unoccc)
  end subroutine gw_product_basis_init
end module m_gw_product_basis

!> Core state read from ECORE file (lmf intermediate output)
module m_core_state
  use m_lgunit, only: stdo
  use m_mpi,    only: ipr
  use m_struct_from_lmf,  only: natom, nspin, nl, nnc, lmxa
  use m_gw_product_basis, only: nindxc, ncwf
  implicit none
  public :: core_state_init
  integer, protected, public :: nctot
  integer, allocatable, protected, public :: konf(:,:), icore(:,:), ncore(:)
  real(8), allocatable, protected, public :: ecore(:,:)
  private
contains
  subroutine core_state_init()
    real(8), external :: rydberg
    real(8), allocatable :: ecoret(:,:,:,:)
    integer :: ifec, iatom, l, n, m, isp, lt, nt, i, j, ncorex, ia

    allocate(icore(nl**2*nnc, natom), ncore(natom), source=99999)
    do iatom = 1, natom
       i = 0
       j = 0
       do l = 0, lmxa(iatom)
          do n = 1, nindxc(l+1,iatom)
             do m = -l, l
                j = j + 1
                if(ncwf(l+1,n,iatom) == 1) then
                   i = i + 1
                   icore(i,iatom) = j
                endif
             enddo
          enddo
       enddo
       ncore(iatom) = i
    enddo
    nctot = sum(ncore(1:natom))

    open(newunit=ifec, file='ECORE')
    read(ifec,*)
    allocate(konf(nl,natom), source=0)
    allocate(ecoret(0:nl-1, nnc, 2, natom), ecore(nctot, 2))
    ecoret = 0d0
    do iatom = 1, natom
       if(ipr) write(stdo,*) ' read ECORE : iatom lmxa =', iatom, lmxa(iatom)
       read(ifec,*)
       read(ifec,*) (konf(l+1,iatom), l=0, lmxa(iatom))
       konf(1:lmxa(iatom)+1, iatom) = [(konf(l+1,iatom)+l+1, l=0, lmxa(iatom))]
       do l = 0, lmxa(iatom)
          ncorex = konf(l+1,iatom) - l - 1
          do n = 1, ncorex
             read(ifec,*) lt, nt, (ecoret(l,n,isp,iatom), isp=1, nspin)
             if(nspin == 1) ecoret(l,n,2,iatom) = ecoret(l,n,1,iatom)
             if(ipr) write(stdo,"(' read ecore=',3i4,2d13.5)") l, n, iatom, ecoret(l,n,1:nspin,iatom)
          enddo
       enddo
    enddo
    ecoret = ecoret / rydberg()
    close(ifec)
    i = 0
    do ia = 1, natom
       iatom = ia
       do l = 0, lmxa(iatom)
          do n = 1, nnc
             do m = -l, l
                if(ncwf(l+1,n,iatom) == 1) then
                   i = i + 1
                   ecore(i,1:nspin) = ecoret(l,n,1:nspin,iatom)
                   if(ipr) write(stdo,"(' ecore=',4i4,2d13.5)") i, l, n, iatom, ecore(i,1:nspin)
                endif
             enddo
          enddo
       enddo
    enddo
    if(size(ecore) == 0) then
       deallocate(ecore)
       allocate(ecore(1,2))
    endif
    deallocate(ecoret)
  end subroutine core_state_init
end module m_core_state

!> Backward-compat wrapper: re-exports all variables that previously came from
!  the monolithic m_genallcf_v3.  Callers (use m_genallcf_v3, only: ...) keep working.
module m_genallcf_v3
  use m_lgunit, only: stdo
  use m_mpi,    only: ipr
  use m_struct_from_lmf, only: natom, nspin, alat, plat, pos, z, spid, lmxa, lmxax, &
       nnv, nnc, nrx, qval, nspc, nlmto, nprecb, mrecb, mrece, nqbzt, nband, mrecg, &
       laf, ibasf, nl, tpioa
  use m_gw_user_config, only: niw, deltaw, esmr, set_esmr
  use m_gw_product_basis, only: lcutmx, nn, nlnx, nlnxv, nlnxc, nlnmx, nlnmxv, nlnmxc, &
       ndima, ndimanspc, nspx, cutbase, nindx, il, in, im, nocc, nunocc, nindxc, lcutmxa, &
       nlnm, nlnmv, nlnmc
  use m_core_state, only: nctot, konf, icore, ncore, ecore
  implicit none
  public :: setesmr, genallcf_v3
  ! Re-export of structural data
  public :: natom, nspin, alat, plat, pos, z, spid, lmxa, lmxax, &
       nnv, nnc, nrx, qval, nspc, nlmto, nprecb, mrecb, mrece, nqbzt, nband, mrecg, &
       laf, ibasf, nl, tpioa
  ! Re-export of GW user config
  public :: niw, deltaw, esmr
  ! Re-export of product basis + derived indices
  public :: lcutmx, nn, nlnx, nlnxv, nlnxc, nlnmx, nlnmxv, nlnmxc, &
       ndima, ndimanspc, nspx, cutbase, nindx, il, in, im, nocc, nunocc, nindxc, lcutmxa, &
       nlnm, nlnmv, nlnmc
  ! Re-export of core state
  public :: nctot, konf, icore, ncore, ecore
  private
  logical, protected, private :: done_genallcf_v3 = .false.
contains
  subroutine setesmr(esmr_in)
    real(8), intent(in) :: esmr_in
    call set_esmr(esmr_in)
  end subroutine setesmr
  subroutine genallcf_v3(incwfx)
    use m_struct_from_lmf,  only: struct_from_lmf_init
    use m_gw_user_config,   only: gw_user_config_init
    use m_gw_product_basis, only: gw_product_basis_init
    use m_core_state,       only: core_state_init
    integer, intent(in) :: incwfx
    if(done_genallcf_v3) call rx('genallcf_v3 is already called')
    done_genallcf_v3 = .true.
    call struct_from_lmf_init()
    call gw_user_config_init()
    call gw_product_basis_init(incwfx)
    call core_state_init()
    call cputid(0)
    if(ipr) write(stdo,*) 'genallcf_v3'
  end subroutine genallcf_v3
end module m_genallcf_v3

module m_ReadEfermi
  use m_lgunit,only:stdo
  use m_mpi,only: ipr
  real(8),protected:: bandgap, ef, ef_kbt
  public:: readefermi,readefermi_kbt,setefermi
contains
  subroutine setefermi(efin)
    real(8)::efin
    ef=efin
  endsubroutine setefermi
  subroutine readefermi()
    implicit none
    integer:: ifief
    open(newunit=ifief,file='EFERMI')
    read(ifief,*) ef,bandgap
    close(ifief)
    if(ipr) write(stdo,"(a,f12.6)")' --- READIN ef from EFERMI. ef=',ef
  end subroutine readefermi
  subroutine readefermi_kbt()
    implicit none
    integer:: ifief_kbt
    open(newunit=ifief_kbt,file='EFERMI_kbt')
    read(ifief_kbt,*) ef_kbt,bandgap
    close(ifief_kbt)
    if(ipr) write(stdo,"(a,f12.6)")' --- READIN ef from EFERMI_kbt. ef=',ef_kbt
  end subroutine readefermi_kbt
end module m_ReadEfermi
