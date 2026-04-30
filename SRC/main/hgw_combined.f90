program main
  ! Combined hrcxq + hsfp0_sc(--job=2) in one process.
  ! Computes WV (Im chi0 + Hilbert) then immediately the correlation self-energy
  ! without re-bootstrapping modules. Output: WVR/WVI/SEC files.
  use m_hrcxq, only: hrcxq
  call hrcxq(do_correlation=.true.)
end program main
