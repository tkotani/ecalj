program main
  ! Combined hrcxq + hsfp0_sc(--job=1) + hsfp0_sc(--job=2) in one process.
  ! Computes:
  !   1. WV (Im chi0 + Hilbert) inside hrcxq main loop (in-memory buffer)
  !   2. Valence exchange Sx (hsfp0_sc kernel, ixc=1) -> SEX files
  !   3. Correlation Sc (hsfp0_sc kernel, ixc=2) -> SEC files
  ! All without re-bootstrapping modules between stages.
  use m_hrcxq, only: hrcxq
  call hrcxq(do_correlation=.true., do_exchange=.true.)
end program main
