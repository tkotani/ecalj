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
!  A key that the file does not carry is added: below the header of its
!  section, at the top of the file for a top-level key, or with a new
!  [section] at the end of the file when the section is absent
!  (2026-09-30 02:23. Until then such an override was skipped with a note,
!  and the run went on with the default value). Only an absent [[sec]]
!  table of an array of tables is still skipped with a warning.
!  A misspelled key is added too, and no program reads it (as with a
!  misspelled key in the file); the note in the output names the key.
!
!  Strict typo detection for non-ctrlg flags and retired-syntax checks
!  (-v..., --[...]=, --toml.<path>=, --pr=, --time=, --phispinsym) live
!  in m_cmdopt_registry::validate_arglist, called from m_setargs at
!  startup. Anything that survives that pass and isn't a --ctrlg:
!  override is silently skipped here.
module m_toml_override
  use m_args,            only: arglist, narg
  use m_lgunit,          only: stdo
  use m_mpi,           only: master_mpi
  implicit none
  private
  public :: load_toml_with_overrides

contains

  !> Read filename as a text buffer, apply any --ctrlg:<path>=<value>
  !  overrides found in arglist, and return the (possibly modified)
  !  buffer in `text`.
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
       if (.not. is_toml_override(arg, path, val)) cycle
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
    integer :: sec_start, sec_end, nl
    logical :: found
    character(1), parameter :: lf = achar(10)
    !
    ! Path forms:
    !   key                   -> top-level scalar
    !   sec.key               -> [sec] table
    !   sec.<idx>.key         -> [[sec]] array-of-tables, idx is integer
    !
    lastdot = index_last(path, '.')
    if (lastdot == 0) then
       ! top-level
       call replace_key_in_range(text, 1, len(text), path, val, found)
       if (.not. found) then ! top-level keys stand before the first [section]
          text = path//' = '//val//lf//text
          if (master_mpi) write(stdo,'(a)') '   (note) key "'//path//'" is not in the file: added at the top'
       endif
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
       if (is_array_of_tables) then
          if (master_mpi) write(stdo,'(a)') &
               '   (warn) --ctrlg: override path "'//trim(path)//'": no such table in the file; skipped'
          return
       endif
       text = text//lf//'['//trim(section_pat)//']'//lf//trim(key)//' = '//val//lf
       if (master_mpi) write(stdo,'(a)') &
            '   (note) section ['//trim(section_pat)//'] is not in the file: added with the key "'//trim(key)//'"'
       return
    endif
    call replace_key_in_range(text, sec_start, sec_end, trim(key), val, found)
    if (.not. found) then ! the new line goes below the header line of the section
       nl = index(text(sec_start:), lf)
       if (nl == 0) then
          text = text//lf//trim(key)//' = '//val//lf
       else
          nl = sec_start + nl - 1
          text = text(1:nl)//trim(key)//' = '//val//lf//text(nl+1:)
       endif
       if (master_mpi) write(stdo,'(a)') &
            '   (note) key "'//trim(key)//'" is not in its section of the file: added'
    endif
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
  !  found = .false. when there is no such line (the caller adds the key).
  subroutine replace_key_in_range(text, p_start, p_end, key, val, found)
    character(len=:), allocatable, intent(inout) :: text
    integer,                       intent(in)    :: p_start, p_end
    character(*),                  intent(in)    :: key, val
    logical,                       intent(out)   :: found
    integer :: i, line_begin, line_eq, line_end, klen, depth
    character :: c
    found = .false.
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
                   ! Find end of *value*, not just end of line.
                   ! Multi-line arrays/inline tables (e.g.
                   !   plat = [[1,0,0],
                   !           [0,1,0],
                   !           [0,0,1]]
                   ! ) must be replaced in full or the trailing lines
                   ! become syntactic garbage. Track [..]/{..} depth and
                   ! skip `#` comments and "..."/'...' string contents so
                   ! brackets there don't disturb the count.
                   line_end = line_eq + 1
                   depth = 0
                   do while (line_end <= p_end)
                      c = text(line_end:line_end)
                      if (c == '"') then
                         line_end = line_end + 1
                         do while (line_end <= p_end .and. text(line_end:line_end) /= '"')
                            if (text(line_end:line_end) == '\' .and. line_end < p_end) &
                                 line_end = line_end + 1
                            line_end = line_end + 1
                         enddo
                      elseif (c == "'") then
                         line_end = line_end + 1
                         do while (line_end <= p_end .and. text(line_end:line_end) /= "'")
                            line_end = line_end + 1
                         enddo
                      elseif (c == '#') then
                         do while (line_end <= p_end .and. text(line_end:line_end) /= achar(10))
                            line_end = line_end + 1
                         enddo
                         cycle
                      elseif (c == '[' .or. c == '{') then
                         depth = depth + 1
                      elseif (c == ']' .or. c == '}') then
                         depth = depth - 1
                      elseif (c == achar(10) .and. depth <= 0) then
                         exit
                      endif
                      line_end = line_end + 1
                   enddo
                   ! splice: keep up to and including '=', insert ' '+val, then newline
                   text = text(1:line_eq) // ' ' // val // text(line_end:)
                   found = .true.
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
