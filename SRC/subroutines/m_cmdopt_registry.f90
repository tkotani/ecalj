!> Central registry of cmdline options whose name takes a value
!! (the `cmdopt2('--foo=', outs)` family). When m_toml_override sees
!! an unknown `--foo=value` arg, it cross-checks here before deciding
!! between "TOML top-level override" vs "user typo, abort".
!!
!! Maintenance:
!!   - One entry per cmdopt2 callsite in SRC/subroutines/ and SRC/main/.
!!   - Drop the trailing `=` from the registry entry; the matcher in
!!     m_toml_override compares `--<word>` after splitting on `=`.
!!   - When adding a new cmdopt2, also add its name here. The opposite
!!     direction (removing a cmdopt2) is caught by the regression tests.
!!   - cmdopt0 (flag-only) options are NOT registered: typo detection
!!     applies only to value-bearing forms.
!!
!! Generated baseline (2026-06-02) by:
!!   grep -rEoh "cmdopt2\('-[^']+'" SRC/subroutines/ SRC/main/ | sort -u
module m_cmdopt_registry
  implicit none
  private
  public :: is_known_cmdopt2

  ! Edit this array when adding/removing cmdopt callsites that take
  ! `=value`. Two flavours live here:
  !
  !   (a) genuine cmdopt2 -- the binary calls cmdopt2('--foo=', outs)
  !       and reads <value> into outs.
  !   (b) cmdopt0 with a literal embedded `=` -- the binary calls
  !       cmdopt0('--foo=show') and the whole token is a flag.
  !       Both forms look like `--foo=...` at the cmdline, so for the
  !       purposes of classify_plain_dashdash they are equivalent.
  !
  ! Single source of truth: keep alphabetised by leading flag.
  character(len=20), parameter :: CMDOPT2_REGISTRY(*) = [character(len=20) :: &
       '--Wtype',           &  ! (a) main_hmagnon.f90: W matrix type selector
       '--cutuu',           &  ! (a) x0kf_ahc.f90: vc matrix cutoff
       '--diag',            &  ! (b) cmdopt0 --diag={default,tridiag,chefsi}
       '--dwnb',            &  ! (b) cmdopt0 --dwnb={mlo,wan}
       '--job',             &  ! (a) hbasfp0/hsfp0/hvccfp0/hahc/hx0fp0/qg4gw/...
       '--jobgw',           &  ! (a) main_lmf.f90: lmf-as-GW-driver stage
       '--nb',              &  ! (a) main_hahc.f90 / main_hrcxq: b-parallel group
       '--nk',              &  ! (a) main_hahc.f90 / main_hmagnon / hrcxq: k-parallel
       '--nww',             &  ! (a) x0kf_ahc.f90: frequency points
       '--quit',            &  ! (b) cmdopt0 --quit={show,ham,mkpot,dmat,band}
       '--sp1',             &  ! (a) hmlo_ovlppair / hmagnon / huumat /...: spin 1
       '--sp2',             &  ! (a) ditto, spin 2
       '-EfermiShifteV',    &  ! (a) m_tetwt / x0kf_ahc: AHC rigid Ef shift (eV)
       '-emax',             &  ! (a) m_writeband: DOS energy max
       '-emin',             &  ! (a) m_writeband: DOS energy min
       '-ndos'              ]  ! (a) m_writeband: DOS points

contains

  !> True if `flag` ("--foo" or "-foo", *without* trailing `=`) is a
  !! registered cmdopt2 name. Case-sensitive (preserves `-EfermiShifteV`).
  function is_known_cmdopt2(flag) result(yes)
    character(*), intent(in) :: flag
    logical :: yes
    integer :: i
    yes = .false.
    do i = 1, size(CMDOPT2_REGISTRY)
       if (trim(flag) == trim(CMDOPT2_REGISTRY(i))) then
          yes = .true.
          return
       endif
    enddo
  end function is_known_cmdopt2

end module m_cmdopt_registry
