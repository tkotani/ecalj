module m_mlo_scrw
  use m_mlo_ham, only: nwf => ndimMTO, ib_tableM, ib_tableI, l_tableM, nsite
  use m_lgunit,only: stdo
  use m_mpi, only: ipr
  use m_ftox, only: ftox
  implicit none
  public :: nnwf_init, scrw_init, trace_onsite, trace_onsite_diag, contract_to_site, extract_diagonal_channel, trace
  integer, protected, public :: nnwf
  complex(8), allocatable, protected, public :: scrw(:,:)
  logical, allocatable, protected, public :: nnwf_mask(:), nnwf2_mask(:)
  integer, allocatable, protected, public :: mlo_pairs(:,:), pair_site(:,:), pair_lorb(:,:)
contains
  subroutine nnwf_init(nnwf_size_reduction)
    logical, intent(in) :: nnwf_size_reduction
    integer :: iwf, jwf, ijwf, idummy
    logical, allocatable :: mask(:)
    integer, allocatable :: iwf_list(:), jwf_list(:)
    iwf_list = [((iwf, iwf=1,nwf), jwf=1,nwf)]
    jwf_list = [((jwf, iwf=1,nwf), jwf=1,nwf)]
    if(nnwf_size_reduction) then
      mask = [((ib_tableM(iwf)==ib_tableM(jwf), iwf=1,nwf), jwf=1,nwf)]  ! only same atomic site
    else
      mask = [((.TRUE., iwf=1,nwf), jwf=1,nwf)]  !full pair
    endif
    if(allocated(nnwf_mask)) deallocate(nnwf_mask)
    allocate(nnwf_mask, source = mask)
    iwf_list = pack(iwf_list, mask=mask)
    jwf_list = pack(jwf_list, mask=mask)
    nnwf = size(iwf_list)
    if (allocated(mlo_pairs)) deallocate(mlo_pairs)
    allocate(mlo_pairs(nnwf,2))
    mlo_pairs(:,1) = iwf_list(:)
    mlo_pairs(:,2) = jwf_list(:)
    if(allocated(pair_site)) deallocate(pair_site)
    if(allocated(pair_lorb)) deallocate(pair_lorb)
    allocate(pair_site(nnwf,2), pair_lorb(nnwf,2))
    do ijwf=1, nnwf
      pair_site(ijwf,1) = ib_tableM(mlo_pairs(ijwf,1))
      pair_site(ijwf,2) = ib_tableM(mlo_pairs(ijwf,2))
      pair_lorb(ijwf,1) =  l_tableM(mlo_pairs(ijwf,1))
      pair_lorb(ijwf,2) =  l_tableM(mlo_pairs(ijwf,2))
   enddo
  end subroutine nnwf_init

  function contract_to_site(mat, lorb) result(cmat)
    complex(8), intent(in) :: mat(nnwf,nnwf)
    complex(8) :: cmat(nsite,nsite)
    integer, intent(in), optional :: lorb
    logical, allocatable :: mask(:)
    integer :: site1, site2, inwf, jnwf
    do concurrent(site1=1:nsite, site2=1:nsite)
      mask = [((pair_site(inwf,1) == site1 .and. pair_site(jnwf,1) == site2 .and. &   !R1 == site1  R3 == site2
                mlo_pairs(inwf,1) == mlo_pairs(inwf,2) .and. & !n1 == n2 => R1 == R2 is automatically satisfied
                mlo_pairs(jnwf,1) == mlo_pairs(jnwf,2), &      !n3 == n4 => R3 == R4 is automatically satisfied
                inwf=1,nnwf), jnwf=1,nnwf)]
      if(present(lorb)) then
        mask = mask .AND. [(((pair_lorb(inwf,1)==lorb .and. pair_lorb(inwf,2)==lorb .and. &
                            & pair_lorb(jnwf,1)==lorb .and. pair_lorb(jnwf,2)==lorb), inwf=1,nnwf), jnwf=1,nnwf)]
      endif
      cmat(site1,site2) = sum(pack(reshape(mat(:,:), shape=[nnwf*nnwf]), mask=mask))
    enddo
  end function contract_to_site

  function extract_diagonal_channel(mat, lorb) result(cmat)
    complex(8), intent(in) :: mat(nnwf,nnwf)
    complex(8) :: cmat(nwf,nwf)
    integer , optional, intent(in) :: lorb
    logical, allocatable :: mask(:)
    integer :: inwf, jnwf, iwf, jwf
    mask = [((mlo_pairs(inwf,1) == mlo_pairs(inwf,2) .and. & !n1 == n2 => R1 == R2 is automatically satisfied
              mlo_pairs(jnwf,1) == mlo_pairs(jnwf,2), &      !n3 == n4 => R3 == R4 is automatically satisfied
              inwf=1,nnwf), jnwf=1,nnwf)]
    cmat(:,:) = reshape(pack(reshape(mat, shape=[nnwf*nnwf]), mask=mask), shape=[nwf,nwf])
    if(present(lorb)) then
      where(reshape([(((l_tableM(iwf)/=lorb .or. l_tableM(jwf)/=lorb), iwf=1,nwf), jwf=1,nwf)], shape=[nwf,nwf])) cmat = (0d0,0d0)
    endif
  end function extract_diagonal_channel

  complex(8) function trace(mat) result(trmat)
    complex(8), intent(in) :: mat(nnwf,nnwf)
    logical, allocatable :: mask(:)
    integer :: inwf, jnwf
    mask = [(( mlo_pairs(inwf,1) == mlo_pairs(inwf,2) .and. & !n1 == n2 => R1 == R2 is automatically satisfied
               mlo_pairs(jnwf,1) == mlo_pairs(jnwf,2), &      !n3 == n4 => R3 == R4 is automatically satisfied
               inwf=1,nnwf), jnwf=1,nnwf)]
    trmat = sum(pack(reshape(mat, shape=[nnwf*nnwf]), mask=mask))
  end function trace

  complex(8) function trace_onsite(mat, site) result(trmat)
    complex(8), intent(in) :: mat(nnwf,nnwf)
    integer, intent(in), optional :: site
    logical, allocatable :: mask(:)
    integer :: inwf, jnwf
    mask = [((pair_site(inwf,1) == pair_site(jnwf,1) .and. &   !R1 == R3
              mlo_pairs(inwf,1) == mlo_pairs(inwf,2) .and. & !n1 == n2 => R1 == R2 is automatically satisfied
              mlo_pairs(jnwf,1) == mlo_pairs(jnwf,2), &      !n3 == n4 => R3 == R4 is automatically satisfied
              inwf=1,nnwf), jnwf=1,nnwf)]
    if(present(site)) then
      mask = mask .AND. [((pair_site(inwf,1) == site .and. pair_site(jnwf,1) == site, &
                          inwf=1,nnwf), jnwf=1,nnwf)]
    endif
    trmat = sum(pack(reshape(mat, shape=[nnwf*nnwf]), mask=mask))
  end function trace_onsite

  complex(8) function trace_onsite_diag(mat, site) result(trmat)
    integer, intent(in), optional :: site
    complex(8), intent(in) :: mat(nnwf,nnwf)
    logical, allocatable :: mask(:)
    integer :: inwf, jnwf
    mask = [(((mlo_pairs(inwf,1) == mlo_pairs(inwf,2) .and. &
               mlo_pairs(jnwf,1) == mlo_pairs(jnwf,2) .and. &
               mlo_pairs(inwf,1) == mlo_pairs(jnwf,1)), &
               inwf=1,nnwf), jnwf=1,nnwf)]
    if(present(site)) then
      mask = mask .AND. [((pair_site(inwf,1) == site .and. pair_site(jnwf,1) == site, &
                          inwf=1,nnwf), jnwf=1,nnwf)]
    endif
    trmat = sum(pack(reshape(mat, [nnwf*nnwf]), mask=mask))
  end function trace_onsite_diag

  subroutine scrw_init(nnwf_size_reduction, w_onsite_dddd, Wtype, enforce_Hermite)
    logical, intent(in) :: nnwf_size_reduction, w_onsite_dddd, enforce_Hermite
    character(len=*), intent(in) :: Wtype
    integer:: ifscrwv, ifscrv, iwf, jwf, kwf, lwf
    character(len=9)::charadummy 
    real(8)::rws1(3),freq,freq2 !dummy
    integer::is,iwf1,iwf2,iwf3,iwf4, idummy
    real(8):: rydberg, hartree
    integer:: ir1, irws1
    complex(8),allocatable::scrw4(:,:,:,:)
    complex(8):: scrv4, scrwc4
    logical(8)::ijklmag
    hartree = 2d0*rydberg()
    allocate( scrw4(nwf,nwf,nwf,nwf), source = (0d0,0d0))
    select case (trim(adjustl(wtype)))
      case ("up")
        if(ipr) write(stdo,ftox) "scrw_init: read Wup"
        open(newunit=ifscrwv,file="Screening_W-v.UP",form="formatted") !only up
        open(newunit=ifscrv, file="Coulomb_v.UP",    form="formatted") !only up
      case ("down")
        if(ipr) write(stdo,ftox) "scrw_init: read Wdn"
        open(newunit=ifscrwv,file="Screening_W-v.DN",form="formatted") !only up
        open(newunit=ifscrv, file="Coulomb_v.DN",    form="formatted") !only up
      case ("up_down")
        if(ipr) write(stdo,ftox) "scrw_init: read Wupdn"
        open(newunit=ifscrwv,file="Screening_W-v.UPDN",form="formatted") !only updw
        open(newunit=ifscrv, file="Coulomb_v.UPDN",    form="formatted") !only updw
      case ("down_up")
        if(ipr) write(stdo,ftox) "scrw_init: read Wdnup"
        open(newunit=ifscrwv,file="Screening_W-v.DNUP",form="formatted") !only updw
        open(newunit=ifscrv, file="Coulomb_v.DNUP",    form="formatted") !only updw
      case default
        call rx("scrw_init: Unknown Wtype")
    endselect
    do iwf=1, nwf**4
      read(ifscrv,"(A,2i5, 3f12.6, 5i5,2f12.6)")charadummy,ir1,irws1,rws1,is,iwf1,iwf2,iwf3,iwf4, scrv4 !v
      read(ifscrwv,"(A,2i5, 3f12.6,5i5,4f12.6)")charadummy,ir1,irws1,rws1,is,iwf1,iwf2,iwf3,iwf4,freq,freq2,scrwc4 !Wc = W -v
      if(w_onsite_dddd) then
        if (all([l_tableM(iwf1), l_tableM(iwf2), l_tableM(iwf3), l_tableM(iwf4)] == 2)) scrw4(iwf2,iwf3,iwf1,iwf4) = scrwc4 + scrv4
        ! if (all([l_tableM(iwf1), l_tableM(iwf2), l_tableM(iwf3), l_tableM(iwf4)] == 2)) scrw4(iwf1,iwf2,iwf3,iwf4) = scrwc4 + scrv4
      else
        scrw4(iwf2,iwf3,iwf1,iwf4) = scrwc4 + scrv4
      endif
    enddo
    if(allocated(scrw)) deallocate(scrw)
    allocate(scrw(nnwf,nnwf))
    if(nnwf_size_reduction) then
      nnwf2_mask = [(((((ib_tableM(iwf1)==ib_tableM(iwf2).and. ib_tableM(iwf3)==ib_tableM(iwf4)), iwf1=1,nwf), iwf2=1,nwf), iwf3=1,nwf), iwf4=1,nwf)]
    else
      nnwf2_mask = [((((.TRUE., iwf1=1,nwf), iwf2=1,nwf), iwf3=1,nwf), iwf4=1,nwf)]
    endif
    scrw(:,:) = reshape(pack(reshape(scrw4, shape=[nwf**4]), mask=nnwf2_mask), shape=[nnwf,nnwf])
    if(enforce_Hermite) scrw(:,:) = (scrw(:,:) + transpose(conjg(scrw(:,:))))*0.5d0
    scrw(:,:) = scrw(:,:)/hartree !! Screening W for magnon
    show_atomic_W: block
      use m_mpi,only: MPI__root
      logical, allocatable :: mask(:), mask_onsite(:), mask_lorb(:), mask_W_diag(:), mask_W_offdiag(:), mask_J(:)
      integer, parameter :: lmax = 3 ! for output
      integer :: iatom, lorb, inwf, jnwf, ib
      complex(8) :: W_diag_ave, W_offdiag_ave, J_ave
      if(MPI__root) then
        do ib = 1, size(ib_tableI)
          iatom = ib_tableI(ib)
          write(stdo, '(A,2I3)') "# Onsite W for atom site:", iatom, ib
          !mask for onsitea R1=R2=R3=R4
          mask_onsite = [(((pair_site(inwf,1)==iatom .and. pair_site(inwf,2)==iatom .and. &
                         &  pair_site(jnwf,1)==iatom .and. pair_site(jnwf,2)==iatom), inwf=1,nnwf), jnwf=1,nnwf)]
          do lorb = 0, lmax
            !mask for orbital combination: all oribitals are lorb
            mask_lorb = [(((pair_lorb(inwf,1)==lorb .and. pair_lorb(inwf,2)==lorb .and. &
                          & pair_lorb(jnwf,1)==lorb .and. pair_lorb(jnwf,2)==lorb), inwf=1,nnwf), jnwf=1,nnwf)]
            mask = mask_onsite .and. mask_lorb
            if(count(mask)==0) cycle
            !"Spin Excitations in Solids from Many-Body Perturbation Theory" Eq (62)
            mask_W_diag = [(((mlo_pairs(inwf,1) == mlo_pairs(inwf,2) .and. &
                           &  mlo_pairs(inwf,1) == mlo_pairs(jnwf,1) .and. &
                           &  mlo_pairs(jnwf,1) == mlo_pairs(jnwf,2)), inwf=1,nnwf), jnwf=1,nnwf)]
            !"Spin Excitations in Solids from Many-Body Perturbation Theory" Eq (63)
            mask_W_offdiag = [(((mlo_pairs(inwf,1) == mlo_pairs(jnwf,1) .and. &
                              &  mlo_pairs(inwf,2) == mlo_pairs(jnwf,2) .and. &
                              &  mlo_pairs(inwf,1) /= mlo_pairs(inwf,2)), inwf=1,nnwf), jnwf=1,nnwf)]
            !"Spin Excitations in Solids from Many-Body Perturbation Theory" Eq (64)
            mask_J = [(((mlo_pairs(inwf,1) == mlo_pairs(inwf,2) .and. &
                       & mlo_pairs(jnwf,1) == mlo_pairs(jnwf,2) .and. &
                       & mlo_pairs(inwf,1) /= mlo_pairs(jnwf,1)), inwf=1,nnwf), jnwf=1,nnwf)]
            W_offdiag_ave = 0d0
            J_ave = 0d0
            if(count(mask .and. mask_W_diag) > 0) W_diag_ave = sum(pack(reshape(scrw, shape=[nnwf*nnwf]), &
                                                           & mask=(mask .and. mask_W_diag)))/count(mask .and. mask_W_diag)
            if(count(mask .and. mask_W_offdiag) > 0) W_offdiag_ave = sum(pack(reshape(scrw, shape=[nnwf*nnwf]), &
                                                           & mask=(mask .and. mask_W_offdiag)))/count(mask .and. mask_W_offdiag)
            if(count(mask .and. mask_J) > 0) J_ave = sum(pack(reshape(scrw, [nnwf*nnwf]),&
                                                           & mask=(mask .and. mask_J)))/count(mask .and. mask_J)
            ! J and offdiagonal W should be zero for s-orbital
            write(stdo, '(A,2I3,3F9.4)') "# site lorb W W' J (eV):", iatom, lorb, &
                 & dble(W_diag_ave)*hartree, dble(W_offdiag_ave)*hartree, dble(J_ave)*hartree
          enddo
        enddo
      endif
    endblock show_atomic_W
  end subroutine scrw_init
end module m_mlo_scrw
