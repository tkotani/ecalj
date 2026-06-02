!> Apply -v[<toml-path>]=<value> overrides to the raw text of a TOML file
!  before toml_loads. Supported paths:
!    -v[symgrp]="find"          top-level scalar
!    -v[ham.gmax]=15            section.key
!    -v[bz.nkabc]=[6,6,6]       section.key with array RHS
!    -v[spec.1.r]=2.5           [[spec]] array-of-tables, 1-based index
!  Each applied override is logged (rank-0 only) so it appears in the
!  console output of every lmf/lmfa/GW utility.
module m_toml_override
  use m_args,   only: arglist, narg
  use m_lgunit, only: stdo
  use m_MPItk,  only: master_mpi
  implicit none
  private
  public :: load_toml_with_overrides

contains

  !> Read filename as a text buffer, apply any -v[<path>]=<value> overrides
  !  found in arglist, and return the (possibly modified) buffer in `text`.
  subroutine load_toml_with_overrides(filename, text)
    character(*),                  intent(in)  :: filename
    character(len=:), allocatable, intent(out) :: text
    integer :: u, ios, sz, i
    character(len=:), allocatable :: path, val
    character(len=512) :: arg
    logical :: did_log_header

    call slurp_file(filename, text)

    did_log_header = .false.
    do i = 1, narg
       arg = arglist(i)
       ! Catch the legacy ctrl %const form `-v<name>=<value>` (no `[`).
       ! Silently ignoring it was a footgun: users would set -vnspin=2 and
       ! lmf would run with the TOML default. Abort with a hint instead.
       if (len_trim(arg) >= 4 .and. arg(1:2) == '-v' .and. arg(3:3) /= '[') then
          if (index(arg(3:), '=') > 0) then
             if (master_mpi) then
                write(stdo,'(a)') ' '
                write(stdo,'(a)') 'ERROR: legacy `-v<name>=<value>` syntax is no longer supported.'
                write(stdo,'(a)') '       Offending argument: '//trim(arg)
                write(stdo,'(a)') '       Use the TOML-path form `-v[<dotted.path>]=<value>` instead.'
                write(stdo,'(a)') '       Examples:'
                write(stdo,'(a)') '         -v[bz.nkabc]=[8,8,8]   (was -vnk=8)'
                write(stdo,'(a)') '         -v[ham.so]=1            (was -vso=1)'
                write(stdo,'(a)') '         -v[ham.scaledsigma]=0.8 (was -vssig=0.8)'
                write(stdo,'(a)') '       See https://ecalj.github.io/ecaljdoc/manual/toml_migration'
             endif
             call rx('legacy -v<name>=<value> override: use -v[<toml-path>]=<value>')
          endif
       endif
       ! Legacy single-key shortcuts. Each one maps to a specific TOML
       ! key; translate to the canonical -v[<path>]=<value> path so the
       ! TOML text is the single source of truth and m_lmfinit.f90 no
       ! longer has to override values after rval2.
       if (translate_legacy_cmdopt(arg, path, val)) then
          ! fall through to the apply block below
       else if (.not. is_v_override(arg, path, val)) then
          cycle
       endif
       if (master_mpi .and. .not. did_log_header) then
          write(stdo,'(a)') ' --- TOML overrides applied (-v[<path>]=<value>) ---'
          did_log_header = .true.
       endif
       if (master_mpi) write(stdo,'(a)') '   '//trim(path)//' = '//trim(val)
       call apply_one_override(text, trim(path), trim(val))
    enddo
  end subroutine load_toml_with_overrides

  !> Translate the three legacy cmdline overrides into the canonical
  !! TOML-path form. They predate -v[<path>]=<value>; keeping them as a
  !! pre-processing translation lets m_lmfinit.f90 stay rval2-only.
  !!   --pr=N      / -pr=N   -> io.verbose = N
  !!   --time=N    / --time=N,M -> io.time = [N, 999] / [N, M]
  !!   --phispinsym         -> ham.phispinsym = true
  function translate_legacy_cmdopt(arg, path, val) result(yes)
    character(*),                  intent(in)  :: arg
    character(len=:), allocatable, intent(out) :: path, val
    logical :: yes
    integer :: alen, comma
    character(len=:), allocatable :: rhs
    yes = .false.
    alen = len_trim(arg)
    if (alen <= 0) return
    if (alen >= 5 .and. arg(1:5) == '--pr=') then
       path = 'io.verbose'
       val  = arg(6:alen)
       yes  = .true.; return
    endif
    if (alen >= 4 .and. arg(1:4) == '-pr=') then
       path = 'io.verbose'
       val  = arg(5:alen)
       yes  = .true.; return
    endif
    if (alen >= 8 .and. arg(1:7) == '--time=') then
       rhs   = arg(8:alen)
       comma = index(rhs, ',')
       path  = 'io.time'
       if (comma > 0) then
          val = '[' // rhs(1:comma-1) // ',' // rhs(comma+1:) // ']'
       else
          ! single value: pad the second slot with the old 999 sentinel
          val = '[' // rhs // ',999]'
       endif
       yes  = .true.; return
    endif
    if (arg(1:alen) == '--phispinsym') then
       path = 'ham.phispinsym'
       val  = 'true'
       yes  = .true.; return
    endif
  end function translate_legacy_cmdopt


  subroutine slurp_file(filename, text)
    character(*),                  intent(in)  :: filename
    character(len=:), allocatable, intent(out) :: text
    integer :: u, ios, sz
    character :: c
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


  function is_v_override(arg, path, val) result(yes)
    character(*),                  intent(in)  :: arg
    character(len=:), allocatable, intent(out) :: path, val
    logical :: yes
    integer :: lb, rb, eq
    yes = .false.
    if (len_trim(arg) < 5) return
    if (arg(1:3) /= '-v[') return
    rb = index(arg, ']=')
    if (rb < 4) return
    path = arg(4:rb-1)
    val  = arg(rb+2:len_trim(arg))
    yes = .true.
  end function is_v_override


  !> Replace the RHS of a single key in text, navigating by dotted path.
  subroutine apply_one_override(text, path, val)
    character(len=:), allocatable, intent(inout) :: text
    character(*),                  intent(in)    :: path, val
    character(len=:), allocatable :: section_pat, key
    character(len=64)             :: idx_str
    integer :: i, lastdot, second_lastdot, idx, sec_count
    logical :: is_array_of_tables
    integer :: sec_start, sec_end, key_pos, eol_pos
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
            '   (warn) -v override path "'//trim(path)//'" not found in TOML; skipped'
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
