!> Apply --ctrlg:<dotted.path>=<value> overrides to the raw text of a TOML file
!  before toml_loads. Supported paths:
!    --ctrlg:symgrp="find"          top-level scalar
!    --ctrlg:verbose=50             top-level scalar (was [io] verbose)
!    --ctrlg:time=[5,5]             top-level array  (was [io] time)
!    --ctrlg:ham.gmax=15            section.key
!    --ctrlg:ham.phispinsym=true    section.key (bool, must be lowercase)
!    --ctrlg:bz.nkabc=[6,6,6]       section.key with array RHS
!    --ctrlg:struc.alat=7.88        section.key
!    --ctrlg:spec.1.r=2.5           [[spec]] array-of-tables, 1-based index
!  Each applied override is logged (rank-0 only) so it appears in the
!  console output of every lmf/lmfa/GW utility. Values are TOML-typed
!  (so `true`/`false` lowercase, strings quoted, arrays in `[...]`).
!
!  Other --foo=value args are checked against
!  m_cmdopt_registry::CMDOPT2_REGISTRY; if it matches a known cmdopt2 (or
!  a cmdopt0 whose name carries `=`, e.g. --quit=show) it passes through
!  to the existing cmdopt path, otherwise we abort with a hint.
!
!  Retired (now aborts with a one-line migration hint):
!    -v<name>=<value>           legacy ctrl %const form
!    -v[<path>]=<value>         the old TOML override prefix
!    --[<path>]=<value>         intermediate bracketed form
!    --<a.b.c>=<value>          intermediate dotted form
!    --pr=N                     cmdline shortcut for verbose
!    --time=<...>               cmdline shortcut for time
!    --phispinsym               cmdline shortcut for [ham] phispinsym
module m_toml_override
  use m_args,            only: arglist, narg
  use m_lgunit,          only: stdo
  use m_MPItk,           only: master_mpi
  use m_cmdopt_registry, only: is_known_cmdopt2
  implicit none
  private
  public :: load_toml_with_overrides

contains

  !> Read filename as a text buffer, apply any --ctrlg:<path>=<value>
  !  overrides found in arglist (plus the retained --pr=N shortcut),
  !  and return the (possibly modified) buffer in `text`.
  subroutine load_toml_with_overrides(filename, text)
    character(*),                  intent(in)  :: filename
    character(len=:), allocatable, intent(out) :: text
    integer :: i
    character(len=:), allocatable :: path, val
    character(len=512) :: arg
    logical :: did_log_header

    call slurp_file(filename, text)

    did_log_header = .false.
    do i = 1, narg
       arg = arglist(i)
       call abort_on_retired_form(arg)
       if (is_toml_override(arg, path, val)) then
          ! new --ctrlg:<dotted.path>=<value> form
       else
          call classify_non_toml_arg(arg)
          cycle
       endif
       if (master_mpi .and. .not. did_log_header) then
          write(stdo,'(a)') ' --- TOML overrides applied (--ctrlg:<path>=<value>) ---'
          did_log_header = .true.
       endif
       if (master_mpi) write(stdo,'(a)') '   '//trim(path)//' = '//trim(val)
       call apply_one_override(text, trim(path), trim(val))
    enddo
  end subroutine load_toml_with_overrides

  !> Canonical TOML override syntax: --ctrlg:<dotted.path>=<value>.
  !! Strip the literal `--ctrlg:` prefix; everything before the first `=`
  !! becomes the TOML path, everything after is the (TOML-typed) value.
  function is_toml_override(arg, path, val) result(yes)
    character(*),                  intent(in)  :: arg
    character(len=:), allocatable, intent(out) :: path, val
    logical :: yes
    integer :: alen, eq
    yes  = .false.
    alen = len_trim(arg)
    if (alen < 10) return                      ! "--ctrlg:x=" minimum is 10
    if (arg(1:8) /= '--ctrlg:') return
    eq = index(arg, '=')
    if (eq <= 9) return                        ! must have something between `:` and `=`
    path = arg(9:eq-1)
    val  = autoquote_bare_string(arg(eq+1:alen))
    yes  = .true.
  end function is_toml_override

  !> Rescue users from shell quote stripping. `--ctrlg:symgrp="find"` is
  !! eaten by bash into `--ctrlg:symgrp=find` and the bare `find` is not
  !! a valid TOML literal -- the override would fail. Detect "this is
  !! clearly meant as a TOML string" and wrap it in double quotes.
  !!
  !! Pass through unchanged when the value is already a TOML literal:
  !!   - starts with " ' [ {            (quoted / array / inline table)
  !!   - is the bool keyword            (true / false)
  !!   - is a float keyword             (nan / inf, with optional +/-)
  !!   - starts with digit + - .        (numeric)
  function autoquote_bare_string(v) result(out)
    character(*), intent(in)      :: v
    character(len=:), allocatable :: out
    character(len=:), allocatable :: trimmed
    integer :: n
    out = v
    trimmed = trim(adjustl(v))
    n = len(trimmed)
    if (n == 0) return
    select case (trimmed(1:1))
    case ('"', "'", '[', '{', '0':'9', '+', '-', '.')
       return
    end select
    if (trimmed == 'true'  .or. trimmed == 'false') return
    if (trimmed == 'nan'   .or. trimmed == 'inf')   return
    out = '"' // v // '"'
  end function autoquote_bare_string

  !> Anything starting with `-` but not the new --ctrlg:/--pr= forms.
  !! Either a registered cmdopt that passes through, or a typo we abort on.
  !! cmdopt0 flags (no `=`) flow through untouched.
  subroutine classify_non_toml_arg(arg)
    character(*), intent(in) :: arg
    integer :: alen, eq
    alen = len_trim(arg)
    if (alen <= 0) return
    if (arg(1:1) /= '-') return                ! not a flag, leave alone
    eq = index(arg, '=')
    if (eq == 0) return                        ! cmdopt0 flag, leave alone
    if (is_known_cmdopt2(arg(1:eq-1))) return  ! registered cmdopt, pass through
    if (master_mpi) then
       write(stdo,'(a)') ' '
       write(stdo,'(a)') 'ERROR: unknown option `'//trim(arg)//'`.'
       write(stdo,'(a)') '       To override a TOML key, use'
       write(stdo,'(a)') '         --ctrlg:<dotted.path>=<value>'
       write(stdo,'(a)') '       For a list of registered cmdline options see'
       write(stdo,'(a)') '         https://ecalj.github.io/ecaljdoc/manual/cmdopts'
       write(stdo,'(a)') '       and m_cmdopt_registry.f90 for the cmdopt2 registry.'
    endif
    call rx('unknown option: '//trim(arg))
  end subroutine classify_non_toml_arg

  !> Hard-error on syntaxes we used to accept. Each branch prints one
  !! migration hint so the failure mode is "obvious from one line".
  subroutine abort_on_retired_form(arg)
    character(*), intent(in) :: arg
    integer :: alen
    alen = len_trim(arg)
    if (alen < 2) return

    ! -v...  family (legacy %const + intermediate -v[...]= TOML override)
    if (alen >= 3 .and. arg(1:2) == '-v') then
       if (index(arg, '=') > 0) call retired_die(arg, &
            '-v...=<value> is retired. Use --ctrlg:<dotted.path>=<value>.')
       return
    endif

    ! --[<path>]=<value>  intermediate bracketed form
    if (alen >= 5 .and. arg(1:3) == '--[' .and. index(arg, ']=') > 0) then
       call retired_die(arg, &
            '--[<path>]=<value> is retired. Use --ctrlg:<dotted.path>=<value>.')
    endif

    ! --toml.<path>=<value>  intermediate dotted prefix
    if (alen >= 9 .and. arg(1:7) == '--toml.' .and. index(arg, '=') > 0) then
       call retired_die(arg, &
            '--toml.<path>=<value> is retired. Use --ctrlg:<dotted.path>=<value>.')
    endif

    ! --pr=N
    if (alen >= 5 .and. arg(1:5) == '--pr=') then
       call retired_die(arg, &
            '--pr=N is retired. Use --ctrlg:verbose=N.')
    endif

    ! --time=...
    if (alen >= 7 .and. arg(1:7) == '--time=') then
       call retired_die(arg, &
            '--time=<...> is retired. Use --ctrlg:time=[N,M].')
    endif

    ! --phispinsym  (cmdline flag retired; lives on inside Legacy2toml.py
    ! which still rewrites the legacy ctrl token, but the Fortran binary
    ! must read it from [ham] phispinsym in the TOML.)
    if (arg(1:alen) == '--phispinsym') then
       call retired_die(arg, &
            '--phispinsym is retired. Use --ctrlg:ham.phispinsym=true, '// &
            'or set [ham] phispinsym = true in ctrlg.<sname>.toml.')
    endif
  end subroutine abort_on_retired_form

  subroutine retired_die(arg, hint)
    character(*), intent(in) :: arg, hint
    if (master_mpi) then
       write(stdo,'(a)') ' '
       write(stdo,'(a)') 'ERROR: retired syntax `'//trim(arg)//'`'
       write(stdo,'(a)') '       '//trim(hint)
       write(stdo,'(a)') '       https://ecalj.github.io/ecaljdoc/manual/toml_migration'
    endif
    call rx('retired syntax: '//trim(arg))
  end subroutine retired_die


  subroutine slurp_file(filename, text)
    character(*),                  intent(in)  :: filename
    character(len=:), allocatable, intent(out) :: text
    integer :: u, ios
    open(newunit=u, file=filename, status='old', action='read', form='formatted', iostat=ios)
    if (ios /= 0) call rx('m_toml_override: cannot open '//filename)
    text = ''
    do
       block
         character(len=4096) :: line
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) exit
         text = text // trim(line) // achar(10)
       end block
    enddo
    close(u)
  end subroutine slurp_file


  !> Replace the RHS of a single key in text, navigating by dotted path.
  subroutine apply_one_override(text, path, val)
    character(len=:), allocatable, intent(inout) :: text
    character(*),                  intent(in)    :: path, val
    character(len=:), allocatable :: section_pat, key
    character(len=64)             :: idx_str
    integer :: lastdot, second_lastdot, idx
    logical :: is_array_of_tables
    integer :: sec_start, sec_end
    !
    ! Path forms:
    !   key                   -> top-level scalar
    !   sec.key               -> [sec] table
    !   sec.<idx>.key         -> [[sec]] array-of-tables, idx is integer
    !
    lastdot = index_last(path, '.')
    if (lastdot == 0) then
       ! top-level
       call replace_key_in_range(text, 1, len(text), path, val)
       return
    endif
    key = path(lastdot+1:len(path))
    second_lastdot = index_last(path(1:lastdot-1), '.')
    if (second_lastdot > 0) then
       idx_str = path(second_lastdot+1:lastdot-1)
       section_pat = path(1:second_lastdot-1)
       read(idx_str, *) idx
       is_array_of_tables = .true.
    else
       section_pat = path(1:lastdot-1)
       idx = 1
       is_array_of_tables = .false.
    endif
    call find_section_range(text, trim(section_pat), is_array_of_tables, idx, &
         sec_start, sec_end)
    if (sec_start <= 0) then
       if (master_mpi) write(stdo,'(a)') &
            '   (warn) --toml override path "'//trim(path)//'" not found in TOML; skipped'
       return
    endif
    call replace_key_in_range(text, sec_start, sec_end, trim(key), val)
  end subroutine apply_one_override


  !> Find the byte range of a [section] / [[section]] block in text.
  subroutine find_section_range(text, sec, aot, want_idx, p_start, p_end)
    character(*), intent(in)  :: text, sec
    logical,      intent(in)  :: aot     ! array-of-tables
    integer,      intent(in)  :: want_idx
    integer,      intent(out) :: p_start, p_end
    character(len=:), allocatable :: needle
    integer :: pos, hit, line_start, line_end, found
    p_start = 0; p_end = 0
    if (aot) then
       needle = '[['//sec//']]'
    else
       needle = '['//sec//']'
    endif
    found = 0
    pos = 1
    do
       hit = index(text(pos:), needle)
       if (hit == 0) exit
       hit = pos + hit - 1
       ! make sure it's at start of line
       if (hit > 1 .and. text(hit-1:hit-1) /= achar(10)) then
          pos = hit + 1
          cycle
       endif
       found = found + 1
       if (found == want_idx) then
          line_start = hit
          ! find end of next [section] or end of text
          line_end = next_section_start(text, hit + len(needle))
          p_start = line_start; p_end = line_end
          return
       endif
       pos = hit + len(needle)
    enddo
  end subroutine find_section_range


  !> Search for next "^[" line after `from_pos`. Returns its position-1
  !  (= last byte of current section), or len(text) if none.
  function next_section_start(text, from_pos) result(p)
    character(*), intent(in) :: text
    integer,      intent(in) :: from_pos
    integer :: p, i
    p = len(text)
    i = from_pos
    do while (i <= len(text))
       if (text(i:i) == achar(10) .and. i < len(text)) then
          if (text(i+1:i+1) == '[') then
             p = i
             return
          endif
       endif
       i = i + 1
    enddo
  end function next_section_start


  !> In text(p_start:p_end), find a line starting with `key` (optionally
  !  followed by whitespace then '='), and replace its RHS with `val`.
  subroutine replace_key_in_range(text, p_start, p_end, key, val)
    character(len=:), allocatable, intent(inout) :: text
    integer,                       intent(in)    :: p_start, p_end
    character(*),                  intent(in)    :: key, val
    integer :: i, line_begin, line_eq, line_end, klen
    character :: c
    klen = len_trim(key)
    line_begin = p_start
    do while (line_begin <= p_end)
       ! skip leading whitespace
       i = line_begin
       do while (i <= p_end .and. (text(i:i) == ' ' .or. text(i:i) == achar(9)))
          i = i + 1
       enddo
       ! key match?
       if (i + klen <= p_end + 1) then
          if (text(i:i+klen-1) == key) then
             ! must be followed by space, '=' or '\t' (not part of longer ident)
             c = text(i+klen:i+klen)
             if (c == ' ' .or. c == achar(9) .or. c == '=') then
                ! look for '='
                line_eq = index(text(i:p_end), '=')
                if (line_eq > 0) then
                   line_eq = i + line_eq - 1
                   ! find end of line
                   line_end = i
                   do while (line_end <= p_end .and. text(line_end:line_end) /= achar(10))
                      line_end = line_end + 1
                   enddo
                   ! splice: keep up to and including '=', insert ' '+val, then newline
                   text = text(1:line_eq) // ' ' // val // text(line_end:)
                   return
                endif
             endif
          endif
       endif
       ! advance to next line
       do while (line_begin <= p_end .and. text(line_begin:line_begin) /= achar(10))
          line_begin = line_begin + 1
       enddo
       line_begin = line_begin + 1
    enddo
    if (master_mpi) write(stdo,'(a)') &
         '   (warn) key "'//trim(key)//'" not found in target section; override skipped'
  end subroutine replace_key_in_range


  function index_last(s, c) result(p)
    character(*), intent(in) :: s, c
    integer :: p, i
    p = 0
    do i = 1, len(s)
       if (s(i:i) == c) p = i
    enddo
  end function index_last

end module m_toml_override
