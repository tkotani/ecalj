!> gwinput_dump: load GWinput.toml via m_GWinput and dump all values.
!  Usage:  gwinput_dump [GWinput.toml]
program gwinput_dump
  use m_GWinput
  implicit none
  character(len=256) :: filename
  character(len=:), allocatable :: errmsg
  integer :: i, n

  if (command_argument_count() >= 1) then
     call get_command_argument(1, filename)
  else
     filename = 'GWinput.toml'
  endif

  call gwinput_load(trim(filename), error=errmsg)
  if (allocated(errmsg) .and. len(errmsg) > 0) then
     write(*,'(a)') 'Error: '//errmsg
     stop 1
  endif
  if (.not. gwinput_loaded) then
     write(*,'(a)') 'Error: gwinput_load did not set gwinput_loaded'
     stop 1
  endif

  write(*,'(a)') '===== [gw] scalars ====='
  write(*,'(a,3i6)')   'n1n2n3       = ', n1n2n3
  write(*,'(a,3i6)')   'n1n2n3eps    = ', n1n2n3eps
  write(*,'(a,i0)')    'BZmesh       = ', BZmesh
  write(*,'(a,f10.5)') 'QpGcut_psi   = ', QpGcut_psi
  write(*,'(a,f10.5)') 'QpGcut_cou   = ', QpGcut_cou
  write(*,'(a,f10.5)') 'alpha_OffG   = ', alpha_OffG
  write(*,'(a,l1)')    'unit_2pioa   = ', unit_2pioa
  write(*,'(a,i0)')    'iSigMode     = ', iSigMode
  write(*,'(a,i0)')    'niw          = ', niw
  write(*,'(a,f10.5)') 'emax_sigm    = ', emax_sigm
  write(*,'(a,f12.6)') 'esmr         = ', esmr
  write(*,'(a,es12.4)')'delta        = ', delta
  write(*,'(a,f10.5)') 'deltaw       = ', deltaw
  write(*,'(a,l1)')    'GaussSmear   = ', GaussSmear
  write(*,'(a,i0)')    'nband_chi0   = ', nband_chi0
  write(*,'(a,f12.3)') 'EMINforGW    = ', EMINforGW
  write(*,'(a,f12.3)') 'EMAXforGW    = ', EMAXforGW
  write(*,'(a,f10.5)') 'HistBin_ratio= ', HistBin_ratio
  write(*,'(a,es12.4)')'HistBin_dw   = ', HistBin_dw
  write(*,'(a,i0)')    'nband_sigm   = ', nband_sigm
  write(*,'(a,f10.5)') 'ecut_p       = ', ecut_p
  write(*,'(a,f10.5)') 'ecuts_p      = ', ecuts_p
  write(*,'(a,3i6)')   'multitet     = ', multitet
  write(*,'(a,f10.5)') 'gauss_img    = ', gauss_img
  write(*,'(a,l1)')    'KeepPositiveCou = ', KeepPositiveCou
  if (allocated(MagAtom)) then
     write(*,'(a,i0,a)', advance='no') 'MagAtom (n=', size(MagAtom), ') ='
     do i = 1, size(MagAtom); write(*,'(i4)', advance='no') MagAtom(i); enddo
     write(*,*)
  else
     write(*,'(a)') 'MagAtom      = (unallocated)'
  endif

  write(*,'(/a)') '===== [product_basis] ====='
  if (allocated(pb_tolerance)) then
     write(*,'(a)', advance='no') 'tolerance    ='
     do i = 1, size(pb_tolerance); write(*,'(es12.4)', advance='no') pb_tolerance(i); enddo
     write(*,*)
  endif
  if (allocated(pb_lcutmx)) then
     write(*,'(a)', advance='no') 'lcutmx       ='
     do i = 1, size(pb_lcutmx); write(*,'(i4)', advance='no') pb_lcutmx(i); enddo
     write(*,*)
  endif
  write(*,'(a,i0)') 'pb_n_nlx  = ', pb_n_nlx
  write(*,'(a,i0)') 'pb_n_val  = ', pb_n_val
  write(*,'(a,i0)') 'pb_n_core = ', pb_n_core

  if (pb_n_nlx > 0) then
     write(*,'(a)') '  nlx (iatom l nnvv nnc):'
     do i = 1, pb_n_nlx
        write(*,'(4x,4i4)') pb_nlx(:,i)
     enddo
  endif
  if (pb_n_core > 0) then
     write(*,'(a)') '  core (iatom l n occ unocc forX0 forSxc):'
     do i = 1, pb_n_core
        write(*,'(4x,7i4)') pb_core(:,i)
     enddo
  endif

  write(*,'(/a)') '===== [blocks] (raw) ====='
  if (allocated(block_QPNT))     write(*,'(a,i0,a)') 'QPNT     (', len(block_QPNT),     ' bytes)'
  if (allocated(block_QforEPS))  write(*,'(a,i0,a)') 'QforEPS  (', len(block_QforEPS),  ' bytes)'
  if (allocated(block_QforEPSL)) write(*,'(a,i0,a)') 'QforEPSL (', len(block_QforEPSL), ' bytes)'
  if (allocated(block_QforGW))   write(*,'(a,i0,a)') 'QforGW   (', len(block_QforGW),   ' bytes)'
  if (allocated(block_Worb))     write(*,'(a,i0,a)') 'Worb     (', len(block_Worb),     ' bytes)'

  write(*,'(/a)') '===== [blocks] (structured) ====='
  write(*,'(a,i0)') 'n_eps  = ', n_eps
  do i = 1, n_eps
     write(*,'(2x,3f12.6)') q_eps(:,i)
  enddo
  write(*,'(a,i0)') 'n_qgw  = ', n_qgw
  do i = 1, n_qgw
     write(*,'(2x,3f12.6)') q_qgw(:,i)
  enddo
  write(*,'(a,i0)') 'n_epsl = ', n_epsl
  do i = 1, n_epsl
     write(*,'(2x,3f10.5,3f10.5,i4)') q_epsl(:,i), qend_epsl(:,i), idx_epsl(i)
  enddo
  if (qpnt_nq > 0) then
     write(*,'(a,2i4)') 'QPNT allq spinonly = ', qpnt_allq, qpnt_spinonly
     write(*,'(a,i0)')  'QPNT nstates = ', qpnt_nstates
     if (allocated(qpnt_bands)) then
        write(*,'(a,10i4)') 'QPNT bands = ', qpnt_bands
     endif
     write(*,'(a,i0)')  'QPNT nq = ', qpnt_nq
     do i = 1, qpnt_nq
        write(*,'(2x,3f12.6)') qpnt_q(:,i)
     enddo
  endif
  if (n_worb > 0) then
     write(*,'(a,i0)') 'n_worb = ', n_worb
     do i = 1, n_worb
        write(*,'(2x,i3,1x,a8,16i4)') worb_iatom(i), worb_label(i), worb_lm(1:worb_nlm(i),i)
     enddo
  endif

  write(*,'(/a)') 'OK!'
end program
