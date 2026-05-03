!> ctrl.<sname>.toml -> recrd(:) loader (schema-typed TOML reader).
!!
!! Walks a structured ctrl.toml (top-level scalars, [section] tables, and
!! [[site]]/[[spec]] arrays-of-tables) and emits flat "<KEY> <VALS>" records
!! that legacy rval2/getdval can parse as before. rval2 is made
!! case-insensitive (see m_gtv.f90) so this loader uppercases all keys.
!!
!! 2026-05-03 T.K. + Claude
module m_ctrl_toml_loader
  use tomlf, only: toml_table, toml_array, toml_keyval, toml_value, &
                   toml_load, toml_error, toml_key, get_value, len, &
                   is_array_of_tables
  implicit none
  private
  public :: load_ctrl_toml

  integer, parameter :: LMAX = 4096  ! max chars per output record

contains

  !> Load ctrl.<sname>.toml and produce recrd(:) for rval2.
  subroutine load_ctrl_toml(filename, recrd, reclnr, nrecs)
    character(*), intent(in) :: filename
    character(len=:), allocatable, intent(out) :: recrd(:)
    integer, intent(out) :: reclnr, nrecs
    type(toml_table), allocatable, target :: root
    type(toml_error), allocatable :: terr
    character(len=LMAX), allocatable :: lines(:)
    integer :: nl, i, lenmax

    call toml_load(root, filename, error=terr)
    if (allocated(terr)) call rx('m_ctrl_toml_loader: '//trim(filename)// &
         ' parse error: '//terr%message)

    allocate(lines(0))
    nl = 0
    call walk_table(root, '', '', lines, nl)

    if (nl == 0) call rx('m_ctrl_toml_loader: empty ctrl TOML?')
    lenmax = 0
    do i = 1, nl
       lenmax = max(lenmax, len_trim(lines(i)))
    enddo
    reclnr = lenmax
    nrecs  = nl
    allocate(character(reclnr):: recrd(nrecs))
    do i = 1, nrecs
       recrd(i) = trim(lines(i))
    enddo
    deallocate(lines)
  end subroutine load_ctrl_toml

  !> Walk a TOML table; for each (k, v) emit recrd line(s).
  recursive subroutine walk_table(tbl, prefix, suffix, lines, nlines)
    type(toml_table), pointer, intent(in)           :: tbl
    character(*), intent(in)                        :: prefix
    character(*), intent(in)                        :: suffix
    character(len=LMAX), allocatable, intent(inout) :: lines(:)
    integer, intent(inout)                          :: nlines
    type(toml_key), allocatable :: keys(:)
    class(toml_value), pointer  :: node
    type(toml_table), pointer   :: subtbl
    type(toml_array), pointer   :: subarr
    type(toml_keyval), pointer  :: subkv
    character(len=:), allocatable :: kup, vstr
    integer :: k, j

    if (.not.associated(tbl)) return
    call tbl%get_keys(keys)
    if (.not.allocated(keys)) return

    do k = 1, size(keys)
       kup = upcase(keys(k)%key)
       call tbl%get(keys(k)%key, node)
       if (.not.associated(node)) cycle
       select type (n => node)
       type is (toml_table)
          subtbl => n
          call walk_table(subtbl, prefix//kup//'_', '', lines, nlines)
       type is (toml_array)
          subarr => n
          if (len(subarr) == 0) cycle
          if (is_array_of_tables(subarr)) then
             do j = 1, len(subarr)
                call get_value(subarr, j, subtbl)
                if (associated(subtbl)) then
                   call walk_table(subtbl, prefix//kup//'_', '@'//int2str(j), &
                                   lines, nlines)
                endif
             enddo
          else
             call array_to_string(subarr, vstr)
             call append_line(prefix//kup//suffix//' '//vstr, lines, nlines)
          endif
       type is (toml_keyval)
          subkv => n
          call kv_to_string(subkv, vstr)
          call append_line(prefix//kup//suffix//' '//vstr, lines, nlines)
       end select
    enddo
  end subroutine walk_table

  !> Stringify a scalar TOML key-value into Fortran string.
  subroutine kv_to_string(kv, out)
    type(toml_keyval), pointer, intent(in)     :: kv
    character(len=:), allocatable, intent(out) :: out
    integer :: ival, stat
    real(8) :: rval
    logical :: lval
    character(len=:), allocatable :: sval
    character(len=64) :: buf

    out = ''
    ! Try integer first (TOML ints are exact)
    call get_value(kv, ival, stat=stat)
    if (stat == 0) then
       write(buf, '(i0)') ival
       out = trim(buf)
       return
    endif
    call get_value(kv, rval, stat=stat)
    if (stat == 0) then
       call format_real(rval, buf)
       out = trim(adjustl(buf))
       return
    endif
    call get_value(kv, lval, stat=stat)
    if (stat == 0) then
       ! Emit 0/1 (rval2/getdval parses as real(8); 'T'/'F' would fail).
       if (lval) then; out = '1'; else; out = '0'; endif
       return
    endif
    call get_value(kv, sval, stat=stat)
    if (stat == 0 .and. allocated(sval)) then
       out = sval
       return
    endif
  end subroutine kv_to_string

  !> Stringify a scalar array (possibly nested) to space-separated values.
  recursive subroutine array_to_string(arr, out)
    type(toml_array), pointer, intent(in)      :: arr
    character(len=:), allocatable, intent(out) :: out
    type(toml_array), pointer :: subarr
    type(toml_keyval), pointer :: subkv
    class(toml_value), pointer :: node
    character(len=:), allocatable :: elem, acc
    integer :: i, ival, stat
    real(8) :: rval
    logical :: lval
    character(len=:), allocatable :: sval
    character(len=64) :: buf

    acc = ''
    if (.not.associated(arr)) then
       out = ''
       return
    endif
    do i = 1, len(arr)
       elem = ''
       call arr%get(i, node)
       if (.not.associated(node)) cycle
       select type (n => node)
       type is (toml_array)
          subarr => n
          call array_to_string(subarr, elem)
       type is (toml_keyval)
          subkv => n
          call kv_to_string(subkv, elem)
       end select
       if (i == 1) then
          acc = elem
       else
          acc = acc // ' ' // elem
       endif
    enddo
    out = acc
  end subroutine array_to_string

  subroutine format_real(x, buf)
    real(8), intent(in) :: x
    character(*), intent(out) :: buf
    write(buf, '(es24.16e3)') x
    buf = adjustl(buf)
  end subroutine format_real

  function int2str(i) result(s)
    integer, intent(in) :: i
    character(len=:), allocatable :: s
    character(len=12) :: buf
    write(buf, '(i0)') i
    s = trim(buf)
  end function int2str

  function upcase(s) result(u)
    character(*), intent(in) :: s
    character(len=len(s)) :: u
    integer :: i, c
    u = s
    do i = 1, len(s)
       c = iachar(s(i:i))
       if (c >= iachar('a') .and. c <= iachar('z')) u(i:i) = achar(c-32)
    enddo
  end function upcase

  subroutine append_line(s, lines, nlines)
    character(*), intent(in) :: s
    character(len=LMAX), allocatable, intent(inout) :: lines(:)
    integer, intent(inout) :: nlines
    character(len=LMAX), allocatable :: tmp(:)
    integer :: cap, i
    cap = size(lines)
    if (nlines+1 > cap) then
       allocate(tmp(max(16, 2*cap)))
       do i = 1, nlines
          tmp(i) = lines(i)
       enddo
       if (nlines+1 <= size(tmp)) then
          ! just-allocated tmp is large enough
       endif
       deallocate(lines)
       call move_alloc(tmp, lines)
    endif
    nlines = nlines + 1
    lines(nlines) = s
  end subroutine append_line

end module m_ctrl_toml_loader
