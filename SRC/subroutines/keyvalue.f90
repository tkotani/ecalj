!--------------------------------------------------------
!     input routines
!     author:   H. Kino
!     date: Nov. 25, 2002
!     last modify:    Oct.  6, 2003
!     Copyright, all rights reserved.
!-------------------------------------------------------

!     This routine is a flexible general readin routine applicable to read
!     real(8),integer(4), logical, arrays (real or integer), char.
!     Also it can return <foo> </foo> block in the input file as the file handle.

!     These usage examples below corresponds to an input test file "IN" as
!     {keyword} {data arry} format as
!     ----IN file start here (remove "c " at head of each line) -----------
!     int_1 325
!     delta  2.35
!     L_test .true.
!     n1n2n3 34 57 58 68 98
!     vv 3.2 4.6 2.5
!     <blockname>
!     abc dgf
!     32 1.4 4.5 9.0 k1 k2
!     18 3.4 4.5 kkk v1
!     </blockname>
!     --- to here -------------------------------------------------------------
!     The data in the <blockname> is a bit special because getkeyvalue just returns
!     file hundle (and line numbers.) to read the contents. See example 5 below.
!     In order to use getkeyvalue, you first have to declear "use keyvalue".
!     See main/kino_input_test.f with "IN" file to see how it works.

!     1.      call getkeyvalue("IN","int_1",ikey,default=999,status=ret )
!     2.      call getkeyvalue("IN","delta",rkey,default=1430d0,status=ret)

!     3.     logical:: lkey
!     call getkeyvalue("IN","L_test",lkey,status=ret )

!     4.     integer(4):: ivkey(5),nsize=5
!     call getkeyvalue("IN","n1n2n3",ivkey,nsize,status=ret )

!     5.     real(8):: rvkey(3)
!     integer(4)::nsize=3
!     call getkeyvalue("IN","vv", rvkey,nsize,status=ret,
!     &                 default=(/10.0d0,10.0d0,10.0d0/))
!     : Note nsize is required. If default is not assigned, this cause stop
!     within getkeyvalue.

!     6.    call getkeyvalue("IN","char", ckey, status=ret)
!     7.    call getkeyvalue("IN","<blockname>",unit=file,status=ret,errstop='on')
!     : unit is the file handle to read contants sandwitched by
!     <blockname> and </blockname>. ret= number of the lines.
!--------------------------------------------------------------------------
module m_keyvalue
  implicit none
  public:: getkeyvalue
  private
  integer :: nf_input=0,output=6,default_nf_input=50
  character(120):: &
       filename0="adfadjpfijwpojobuwhfqeijqope"
  integer,parameter:: input_err_success=0
  integer,parameter:: input_err_conv=-1
  integer,parameter:: input_err_lessvaluenum=-2
  integer,parameter:: input_err_nokey=-3
  integer,parameter:: input_err_nofile=-4
  integer,parameter:: nlogical=14
  character(7):: logicalstr(nlogical)= &
       (/'true   ', 'ok     ', '.true. ', 'yes    ','on     ', &
       '1      ', 'T      ', &
       'false  ', 'ng     ', '.false.', 'no     ','off    ', &
       '0      ', 'F      '/)
  integer :: logicalint(nlogical) = &
       (/1,1,1,1,1,1,1,0,0,0,0,0,0,0/)
  integer:: printlevel=0 !!debug
  integer,parameter:: max_error_count=5 !!error report
  integer :: error_count=0
  character(120):: error_message(max_error_count)
  interface getkeyvalue
     module procedure getkeyvalue_i, getkeyvalue_r8, &
          getkeyvalue_iv, getkeyvalue_rv, getkeyvalue_l  , &
          getkeyvalue_c
  end interface getkeyvalue
  !! ----------------------------------------------------------
contains
  !----------------------------------------------------------
  !     used only for internal use

  integer function unusedfid()
    implicit none
    integer:: i
    logical ::L
    !     check the last file is still open or not?
    !     nf_input==0 means the first time, no check is necessary.

    !     if (nf_input.ne.0) then
    !     inquire(nf_input,opened=L)
    !     takao aug2005. This is because nf_input might be closed and opened in other part---
    unusedfid=-9999
    inquire(opened=L,file=filename0)
    if (L) then
       call  getkeyvalue_err_exit( &
            'unusedfid(getkeyvalue module)',' The last file is still opened. Close it first.')
    endif
    !     endif

    !     find new unused file id
    do i=10,999
       inquire(i,opened=L)
       if ( .NOT. L) then
          unusedfid=i
          return
       endif
    enddo
    !     error
    call  getkeyvalue_err_exit( &
         'unusedfid(getkeyvalue module)' ,' no fid left')
  end function unusedfid
  !----------------------------------------------------------

  subroutine input_delfirstspc(buf)
    implicit none
    character(*):: buf
    integer:: i


    buf=adjustl(buf)
    ! cut comment
    if (buf(1:1) == '#' .OR. buf(1:1) == '!') then
       buf=''
       return
    endif

    !        do i=1,len_trim(buf)
    !           if (buf(i:i).eq.'#' .or. buf(i:i).eq.'!') buf(i:)=''
    !        enddo

  end subroutine input_delfirstspc

  !-----------------------------------------------------------

  subroutine input_toupper(buf)
    implicit none
    character(*) :: buf
    integer:: i,n

    n=len_trim(buf)
    do i=1,n
       if (  ichar(buf(i:i)) >= ichar('a') .AND. &
            ichar(buf(i:i)) <= ichar('z')) then
          buf(i:i) = char(ichar(buf(i:i))+ichar('Z')-ichar('z'))
       endif
    enddo
  end subroutine input_toupper
  !-----------------------------------------------------------

  subroutine input_tolower(buf)
    implicit none
    character(*) :: buf
    integer:: i,n

    n=len_trim(buf)
    do i=1,n
       if (  ichar(buf(i:i)) >= ichar('A') .AND. &
            ichar(buf(i:i)) <= ichar('Z')) then
          buf(i:i) = char(ichar(buf(i:i))+ichar('z')-ichar('Z'))
       endif
    enddo
  end subroutine input_tolower


  !-----------------------------------------------------------
  subroutine input_debugstart(id)
    integer,intent(in),optional:: id
    if (present(id) ) then
       printlevel=id
    else
       printlevel=1
    endif
  end subroutine input_debugstart

  subroutine input_debugend()
    printlevel=0
  end         subroutine input_debugend


  !----------------------------------------------------------

  subroutine input_registerError(msg)
    implicit none
    character(*):: msg
    error_count=error_count+1
    if (error_count<=max_error_count) then
       error_message(error_count)=msg
    endif
  end         subroutine input_registerError

  !----------------------------------------------------------
  integer function input_printError(file)
    implicit none
    integer,optional,intent(in):: file

    integer:: i,n,nf
    character(120) :: buf
    n=error_count
    if (error_count>max_error_count) n=max_error_count

    nf=6
    if (present(file))  nf=file

    do i=1,n
       buf = error_message(i)
       write(nf,'(a)') buf(:len_trim(buf))
    enddo

    input_printError=error_count

  end function input_printError

  !----------------------------------------------------------

  subroutine input_open(name,nf,success,outputfile)
    implicit none
    character(*),intent(in)::name
    integer,intent(in):: nf
    logical,optional,intent(out)::success
    integer,optional,intent(in):: outputfile
    logical:: flag

    nf_input = nf
    if (printlevel>0) then
       write(output,*) 'filename=<',name,'>'
    endif
    open(unit=nf_input,file=name,status='old',err=1000)
    !        open(unit=nf_input,file=name,err=1000)
    filename0=name
    if (printlevel>0) then
       write(output,*) 'succeed to open a file ',name
    endif
    flag=.true.
    error_count=0   ! reset error
    goto 2000
1000 flag=.false.
    write(output,*) 'failed to open a file ',name
2000 if (present(success)) success=flag
    output=6
    if (printlevel>0) write(6,*) 'input/output=',nf_input,output
    if (present(outputfile)) output=outputfile
  end  subroutine input_open

  !----------------------------------------------------------
  subroutine input_redirect(nfinp,name,nf,success,outputfile)
    implicit none
    character(*),intent(in)::name
    integer,intent(in):: nf,nfinp
    logical,optional,intent(out)::success
    integer,optional,intent(in):: outputfile
    logical:: flag
    character(200):: buf

    nf_input = nf
    write(output,*) 'filename=<',name,'>'
    open(unit=nf_input,file=name,status='unknown',err=1000)
    write(output,*) 'succeed to open a file ',name

    do
       read(nfinp,'(a200)',end=1500) buf
       write(nf_input,'(a)') buf(:len_trim(buf))
    enddo

1500 continue

    filename0=name

    close(nf_input)
    open(unit=nf_input,file=name,status='old',err=2000)

    flag=.true.
9000 continue
    error_count=0   ! reset error
    if (present(success)) success=flag
    output=6
    if (present(outputfile)) output=outputfile
    return
    !--------------------------
1000 flag=.false.
    write(output,*) 'failed to open a file ',name
    error_count=0   ! reset error
    if (present(success)) success=flag
    output=6
    if (present(outputfile)) output=outputfile
    return
    !---------------------------
2000 flag=.false.
    write(output,*) 'failed to open a file ',name
    error_count=0   ! reset error
    if (present(success)) success=flag
    output=6
    if (present(outputfile)) output=outputfile
    return
    !--------------------------
  end  subroutine input_redirect
  !----------------------------------------------------------
  subroutine input_reopen(nf)
    integer,intent(out),optional:: nf
    if (present(nf)) then
       nf =nf_input
    endif
    open(unit=nf_input,file=filename0(:len_trim(filename0)),status='old',err=1000)

    return
1000 continue
    write(output,*) 'failed to open a file ',filename0(:len_trim(filename0))
  end  subroutine input_reopen

  !----------------------------------------------------------
  subroutine input_close()
    implicit none
    close(nf_input)
    nf_input=0
  end subroutine input_close

  !----------------------------------------------------------------
  subroutine input_characterv(key0,n,i,defaultvalue,success,nout)
    implicit none
    character(*),intent(in)::key0
    integer,intent(in):: n
    character(*),intent(out):: i(n)
    character(*), intent(in):: defaultvalue(n)
    integer,intent(out),optional:: success

    character(120):: buf,keybuf,key
    integer:: nkey,nkeybuf,j

    integer :: flag,ix,noutx

    integer,optional :: nout !sep2006 takao

    i=defaultvalue
    key=key0
    call input_toupper(key)
    !        write(output,*) 'key=(',key,')'
    nkey=len_trim(key)
    rewind(nf_input)
    flag=input_err_nokey
    buf=' '
    ! cccccccccccccccc
    ix=0
    ! ccccccccccccccc
    do

       !             read and convert them to uppercase string

       ! cccccccccccccccccccccccccccccccccc
       !          ix=ix+1
       !          write(6,"(i4,a120)")ix, buf
       ! cccccccccccccccccccccccccccccccc

       read(nf_input,'(a120)',end=1000) buf
       if (printlevel>10) then
          write(output,*) '(',buf(:len_trim(buf)),')'
       endif
       call input_delfirstspc(buf)

       if (len_trim(buf) == 0) cycle


       !               read the part corresponding to key, -> uppercase

       read(buf,*)  keybuf
       call input_toupper(keybuf)
       nkeybuf=len_trim(keybuf)
       if (printlevel>10) then
          write(output,*) key,keybuf
       endif

       !            compare length and strings

       if ( keybuf(1:nkeybuf) == key(1:nkey) ) then
          flag=input_err_lessvaluenum
          if (printlevel > 10) then
             write(output,*) 'input_characterv:found<',buf(nkey+1:), &
                  '>'
          endif
          !              read(buf(nkey+1:),*,end=2000,err=2000) i(1:n)
          buf=buf(nkey+1:)
          buf= adjustl(buf)

          call input_septoken(buf,n,i,noutx) !takao sep2006
          !          write(6,*) 'readin nout=',noutx
          if(present(nout)) nout=noutx
          if (printlevel >= 1) then
             write(output,*)'characterv:', key(:len_trim(key)), (i(j),j=1,noutx)
          endif

          flag= input_err_success
          goto 1000
       endif
    enddo
1000 continue
    ! write(output,*) '1000 end'
    if (present(success)) success =flag
    return


    !           unexpected intrrupt of read

2000 continue  !takao now no jump to 2000. So dummy below...----------
    if (printlevel >= 10) then
       write(output,*)'parameter error, key=',key(:nkey), &
            ' must has ',n,' value(s)'
    endif

    write(buf,*)'parameter error, key=',key(:nkey), &
         ' must has ',n,' value(s)'
    call input_registerError(buf)
    flag=input_err_lessvaluenum
    if (present(success)) success =flag
  end subroutine input_characterv
  !----------------------------------------------------------------
  subroutine input_logicalv(key,n,i,default,success)
    implicit none
    character(*),intent(in):: key
    integer,intent(in):: n
    logical,intent(out):: i(n)
    logical,intent(in):: default(n)
    integer,intent(out),optional:: success


    integer :: flag
    integer ::j,strlen,logicallen,istr,ilog
    character(10):: str(n),defstr(n),str1,str2
    character(80):: buf
    ! set defaults
    do j=1,n
       i(j) = default(j)
       defstr(j)=' '
       str(j)=' '
    enddo
    ! read as strings
    call input_characterv(key,n,str,defstr,flag)
    if (printlevel>10) then
       write(output,*) 'after input_characterv: ',flag,str
    endif
    if (flag == input_err_success) then
       !            successfully read, then convert
       !            read(..., end) is done in input_characterv()
       do istr=1,n
          str1=str(istr)
          call input_tolower(str1)
          strlen=len_trim(str1)
          do ilog=1,nlogical
             str2=logicalstr(ilog)
             call input_tolower(str2)
             logicallen=len_trim(str2)
             if (printlevel>10) then
                write(output,*)str1(:len_trim(str1)),str2(:len_trim(str2)), &
                     str1(:strlen).eq.str2(:logicallen)
             endif
             if ( &
                  str1(:strlen) == str2(:logicallen) ) then
                if (logicalint(ilog) == 1) then
                   i(istr) = .true.
                else
                   i(istr) = .false.
                endif
                goto 1000
             endif
          enddo
          goto 2000
1000      continue
       enddo
       if (present(success)) success=flag
       if (printlevel >= 1) then
          write(output,*) key(:len_trim(key)),'=',(i(j),j=1,n)
       endif
       return
    endif

    if (present(success)) success=flag
    if (printlevel >= 1) then
       write(output,*) key(:len_trim(key)),'=',(i(j),j=1,n)
    endif
    return

2000 continue   ! error
    if (present(success)) success=input_err_conv
    write(buf,*) 'conversion error, key=', key(:len_trim(key)), &
         ' ', n,' values=logical(s)'
    call input_registerError(buf)

  end subroutine input_logicalv

  !----------------------------------------------------------------
  subroutine input_integerv(key,n,i,defaultvalue,success)
    implicit none
    character(*),intent(in)::key
    integer,intent(in):: n
    integer,intent(out):: i(n)
    integer,intent(in):: defaultvalue(n)
    !        integer,intent(out),optional:: success
    integer,optional:: success

    character(40):: str(n),defstr(n)
    character(120) :: buf
    integer :: j
    integer :: flag,nout

    ! set defaults
    do j=1,n
       i(j) = defaultvalue(j)
       defstr(j)=''
    enddo
    ! read as strings
    call input_characterv(key,n,str,defstr,flag,nout)
    if (printlevel>10) then
       write(output,*) 'after input_characterv: ',flag,str
       write(output,*) 'input_err_success=',input_err_success
       write(output,*) 'nout=',nout
    endif
    ! ---
    if(flag == input_err_success) then
       !            successfully read, then convert
       !            read(..., end) is done in input_characterv()
       do j=1,nout
          read(str(j),*,err=2000,end=3000)i(j) !(i(j),j=1,nout) !takao sep2006
       enddo
       !        write(6,*)" aaaaaaaaaaa success",nout,success
       if (present(success)) then
          if(success/=-9999) then
             if(n/=nout) goto 3000
          endif
          success = nout
       endif
    else
       if (present(success)) success = flag
    endif

    if( printlevel >= 1 ) then
       write(output,*) key(:len_trim(key)),  (i(j),j=1,nout)
    endif

    return


    !  conversion error of read()

2000 flag=input_err_conv
    do j=nout+1,n
       i(j) = defaultvalue(j)
    enddo
    if (present(success)) success = -nout + 1000*flag
    write(buf,*) 'conversion error, key=', key(:len_trim(key)), &
         ' ', n,' values=integer'
    call input_registerError(buf)

    return
3000 success=input_err_lessvaluenum
    return

  end subroutine input_integerv

  !---------------------------------------------------------------
  subroutine input_real8v(key,n,i,defaultvalue,success)
    implicit none
    character(*),intent(in)::key
    integer,intent(in):: n
    real*8,intent(out):: i(n)
    real*8, intent(in):: defaultvalue(n)
    integer,intent(out),optional:: success

    character(40):: str(n),defstr(n)
    character(120) :: buf
    integer :: j
    integer :: flag

    ! set defaults
    do j=1,n
       i(j) = defaultvalue(j)
       defstr(j)=''
    enddo
    ! read as strings
    call input_characterv(key,n,str,defstr,flag)
    if (printlevel>10) then
       write(output,*) 'after input_characterv: ',flag,str
    endif
    if (flag == input_err_success) then
       !            successfully read, then convert
       !            read(..., end) is done in input_characterv()
       do j=1,n
          read(str(j),*,err=2000,end=3000) i(j)
       enddo
       if (present(success)) success = n
    else
       if (present(success)) success = flag
    endif

    if (printlevel >= 1) then
       write(output,*) key(:len_trim(key)), (i(j),j=1,n)
    endif
    return

    !  conversion error of read()

2000 flag=input_err_conv
    do j=1,n
       i(j) = defaultvalue(j)
    enddo
    if (present(success)) success = flag

    write(buf,*) 'conversion error, key=', key(:len_trim(key)), &
         ' ',n,' values=real*8'
    call input_registerError(buf)
    return
3000 success= input_err_lessvaluenum
    return

  end subroutine input_real8v

  !----------------------------------------------------------------
  function input_character( &
       key,defaultvalue,success)
    implicit none
    character(*),intent(in)::key
    character(*),intent(in):: defaultvalue
    integer,intent(out),optional:: success
    character(120) input_character

    character(120):: i(1),defaultvalue0(1)

    integer:: flag

    input_character= defaultvalue
    defaultvalue0(1)=defaultvalue

    ! use input_characterv
    call input_characterv(key,1,i,defaultvalue0,flag)
    if (flag == input_err_success) then
       input_character=i(1)
    endif

    if (printlevel >= 1) then
       write(output,*) key(:len_trim(key)), '=', &
            input_character(:len_trim(input_character))
    endif


    if (present(success)) success=flag
  end function input_character

  !---------------------------------------------------------------
  integer function input_integer(key,defaultvalue,success)
    implicit none
    character(*),intent(in)::key
    integer,intent(in):: defaultvalue
    integer,intent(out),optional:: success


    character(120):: i(1),defaultvalue0(1),buf

    integer:: flag

    input_integer=defaultvalue
    defaultvalue0(1)=''
    call input_characterv(key,1,i,defaultvalue0,flag)
    if (flag == input_err_success) then
       read(i(1),*,err=2000) input_integer
    else
       input_integer = defaultvalue
    endif
    if (present(success)) success=flag

    if (printlevel >= 1) then
       write(output,*) key(:len_trim(key)),'=', input_integer
    endif
    return

    !  conversion error of read()

2000 flag= input_err_conv
    input_integer=defaultvalue
    if (present(success)) success=flag

    write(buf,*) 'conversion error, key=',key(:len_trim(key)), &
         ' value=integer'
    call input_registerError(buf)


  end function input_integer

  !---------------------------------------------------------------
  real(8) function input_real8(key,defaultvalue,success)
    implicit none
    character(*),intent(in)::key
    real*8,intent(in):: defaultvalue
    integer,intent(out),optional:: success
    character(120):: i(1),defaultvalue0(1),buf
    integer:: flag
    input_real8=defaultvalue
    defaultvalue0(1)=''
    call input_characterv(key,1,i,defaultvalue0,flag)
    if (flag == input_err_success) then
       !        write(output,*) i(1)
       read(i(1),*,err=2000) input_real8
    endif
    if (present(success)) success=flag

    if (printlevel >= 1) then
       write(output,*) key(:len_trim(key)),'=', input_real8
    endif

    return

    !   conversion error of read()

2000 flag=input_err_conv
    input_real8=defaultvalue
    if (present(success)) success=flag

    write(buf,*) 'conversion error, key=',key(:len_trim(key)), &
         ' value=real*8'

    call input_registerError(buf)

  end function input_real8



  !----------------------------------------------------------------

  integer function input_character2integer( &
       key,n,ans,vans,success)
    implicit none
    integer,intent(in):: n,vans(n)
    character(*),intent(in):: key,ans(n)
    integer,intent(out),optional::success

    character(120):: buf,compstr
    integer:: flag
    integer :: nbuf,i
    !        character(120) input_character


    !  read a key and its value

    buf= input_character(key,' ',flag)
    if (printlevel>10) then
       write(output,*)'key=',key, 'value=',buf,'flag=',flag
    endif
    input_character2integer=vans(1)  !default value


    if (flag == input_err_success) then

       flag=input_err_conv
       call input_toupper(buf)
       nbuf=len_trim(buf)

       do i=1,n
          compstr=ans(i)
          call input_toupper(compstr)
          if (printlevel>10) then
             write(output,*) i,buf(:nbuf),compstr(:nbuf)
          endif

          !            compare them in the upper case

          if ( &
               buf(:nbuf) == compstr(:len_trim(compstr)) ) then
             input_character2integer=vans(i)
             flag=input_err_success
             if (present(success)) success=flag
             if (printlevel >= 1) then
                write(output,*) key(:len_trim(key)),'=', &
                     input_character2integer
             endif
             return
          endif
       enddo

       ! --- no matching
       write(buf,*)'conversion error: key=',key(:len_trim(key))
       call input_registerError(buf)

    endif

    if (printlevel >= 1) then
       write(output,*) key(:len_trim(key)),'=',input_character2integer, &
            ' (default value)'
    endif

    if (present(success)) success=flag
  end function input_character2integer

  !----------------------------------------------------------------
  logical function input_logical(key,default,success)
    implicit none
    character(*),intent(in):: key
    logical,intent(in):: default
    integer,intent(out),optional:: success

    integer :: flag
    integer ::ret
    !        integer input_character2integer

    ret=input_character2integer(key,nlogical &
         ,  logicalstr, logicalint, flag)
    if (flag == input_err_success) then
       if (ret==1) then
          input_logical=.true.
       else
          input_logical=.false.
       endif
    else
       input_logical=default
    endif
    if (present(success)) success=flag
    if (printlevel >= 1) then
       write(output,*) key(:len_trim(key)),'=', input_logical
    endif
  end function input_logical


  !----------------------------------------------------------------
  logical function input_find(key0,nf,cont)
    implicit none
    character(*),intent(in)::key0
    integer,intent(out):: nf   !  maybe out
    logical,intent(in),optional:: cont
    character(120):: buf,key
    integer:: nkey,nbuf
    logical:: flag


    nf=0
    key=key0
    call input_toupper(key)
    nkey=len_trim(key)
    flag=.true.
    if (present(cont)) then
       if (cont) then
          flag=.false.
       endif
    endif
    if (flag) rewind(nf_input)

    do
       read(nf_input,*,end=1000) buf
       if (len_trim(buf) == 0) cycle
       call input_delfirstspc(buf)
       nbuf=len_trim(buf)
       call input_toupper(buf)
       if (printlevel>10) then
          write(output,*) '<',buf(:nkey),'|',key(:nkey),'>', &
               buf(:nbuf).eq.key(:nkey)
       endif
       if ( buf(:nbuf) == key(:nkey)) then
          if (printlevel>0) then
             write(output,*) key(:len_trim(key)),' found'
          endif
          input_find=.true.
          nf=nf_input
          return
       endif
    enddo
1000 input_find=.false.
  end function input_find

  !-------------------------------------------------------------

  logical function input_last(key0)
    implicit none
    character(*),intent(in)::key0
    character(120):: buf,key
    integer:: nkey


    key=key0
    call input_toupper(key)
    nkey=len_trim(key)

    read(nf_input,'(a120)',end=1000) buf
    call input_delfirstspc(buf)
    call input_toupper(buf)
    !        write(output,*) buf(:nkey),key(:nkey),
    !     .      buf(1:nkey).eq.key(:nkey)
    if ( &
         buf(1:len_trim(buf)) == key(:nkey)) then
       if (printlevel>0) then
          write(output,*) key(:len_trim(key)),' found'
       endif
       input_last=.true.
       return
    endif
1000 continue
    input_last=.false.
  end function input_last

  !-------------------------------------------------------------
  integer function input_existfiles(n,filenames,ret)
    implicit none
    integer,intent(in):: n
    character(*),intent(in) :: filenames(n)
    integer,intent(in):: ret(n+1)
    integer:: i
    character(80):: buf
    logical :: fexist

    do i=1,n
       buf=filenames(i)
       if (printlevel>0)then
          write(6,*) 'open ',buf(:len_trim(buf))
       endif

       !          open(fileid,file=buf(:len_trim(buf)),status='old',err=1000)
       INQUIRE(file=buf(:len_trim(buf)),EXIST=FEXIST)

       if (fexist) then
          if (printlevel>0) then
             write(6,*) 'exist'
          endif
          input_existfiles=ret(i)
          !             close(fileid)
          return
       else
          if (printlevel>0) then
             write(6,*) 'not exist'
          endif
       endif

    enddo
    input_existfiles=ret(n+1)
  end function input_existfiles
  !-------------------------------------------------------------
  subroutine input_septoken(line,n,bufs,nout)
    implicit none
    integer,intent(in):: n
    character(*),intent(in):: line
    character(*),intent(out):: bufs(n)

    integer::i,j,max,ibuf
    integer:: nout

    do i=1,n
       bufs(i)=' '
    enddo

    max=len_trim(line)
    if (printlevel>0)then
       write(*,*) max,line(:max)
    endif

    ibuf=0
    i=1
    do
       do while (line(i:i).eq.' ')
          i=i+1
          if (i>max) goto 1000
       enddo
       j=i+1
       do  while (line(j:j).ne.' ')
          j=j+1
          if (j>max) exit
       enddo
       j=j-1
       if (printlevel>0) then
          write(*,*) i,j,line(i:j)
       endif

       ibuf=ibuf+1
       bufs(ibuf)=line(i:j)
       if (ibuf>=n) goto 1000

       i=j+1
       if (i>max) goto 1000

    enddo

1000 continue
    nout=ibuf !takao Sep2006
  end  subroutine input_septoken
  !-------------------------------------------------------------
  subroutine input_replacekeys(str,nkey,key,val,ret,ierror)
    implicit none
    integer,intent(in) :: nkey
    character(*),intent(in) :: str,val(nkey),key(nkey)
    character(*),intent(out) :: ret
    integer,intent(out),optional:: ierror

    character(200):: buf,buf2,buf3,buf4
    integer :: i,lenstr,lenkey,ikey,ipos(nkey)
    logical :: flag


    lenstr=len_trim(str)

    ipos(:)=0

    !       do i=1,lenstr-1
    !          findkey: do ikey=1,nkey
    !c          write(*,*) str(i:i+1),key(ikey), str(i:i+1).eq.key(ikey)
    !          if (str(i:i+1).eq.key(ikey)) then
    !             ipos(ikey)=i
    !             exit findkey
    !          endif
    !          enddo findkey
    !       enddo
    do i=1,lenstr-1
       do ikey=1,nkey
          !          write(*,*) str(i:i+1),key(ikey), str(i:i+1).eq.key(ikey)
          if (str(i:i+1) == key(ikey)) then
             ipos(ikey)=i
             exit
          endif
       enddo
    enddo


    ! check error
    flag=.true.
    do ikey=1,nkey
       if (ipos(ikey) == 0) then
          flag=.false.
          exit
       endif
    enddo

    if ( .NOT. flag) then
       if (present(ierror)) ierror=input_err_conv
       ret=' '
       return
    endif

    ! substitute
    buf=str
    !       do ikey=1,nkey
    !       searchbuf: do i=1,len_trim(buf)
    !       if (buf(i:i+1).eq.key(ikey)) then
    !          buf2=buf(i+2:) ! last part
    !          buf3=val(ikey) ! middle part
    !          buf4=buf(:i-1) ! first part
    !          if (i.eq.1) then
    !          write(buf,'(a,a)')buf3(:len_trim(buf3)),
    !     i          buf2(:len_trim(buf2))
    !          else
    !          write(buf,'(a,a,a)')buf4(:len_trim(buf4)),buf3(:len_trim(buf3)),
    !     i          buf2(:len_trim(buf2))
    !          endif
    !          exit searchbuf
    !       endif
    !       enddo searchbuf
    do ikey=1,nkey
       do i=1,len_trim(buf)
          if (buf(i:i+1) == key(ikey)) then
             buf2=buf(i+2:) ! last part
             buf3=val(ikey) ! middle part
             buf4=buf(:i-1) ! first part
             if (i == 1) then
                write(buf,'(a,a)')buf3(:len_trim(buf3)), &
                     buf2(:len_trim(buf2))
             else
                write(buf,'(a,a,a)')buf4(:len_trim(buf4)),buf3(:len_trim(buf3)), &
                     buf2(:len_trim(buf2))
             endif
             exit
          endif
       enddo
       !       write(*,*) buf
       i=i+len_trim(buf3)
    enddo

    ret=buf
    if (present(ierror)) ierror=input_err_success


  end subroutine input_replacekeys

  !----------------------------------------------------
  subroutine getkeyvalue_err_exit(s1,s2,s3)

    character(*) :: s1,s2
    character(*),optional:: s3

    if (present(s3) ) then

       write(*,*) s1(:len_trim(s1)), ' : ', s2(:len_trim(s2)),'(',s3(:len_trim(s3)),')'
       ! top2rx 2013.08.09 kino        stop 'in getkeyvalue_err_exit'
       call rx( 'in getkeyvalue_err_exit')

    else

       write(*,*) s1(:len_trim(s1)), ' : ', s2(:len_trim(s2))
       ! top2rx 2013.08.09 kino        stop 'in getkeyvalue_err_exit'
       call rx( 'in getkeyvalue_err_exit')

    endif


  end subroutine getkeyvalue_err_exit


  !-----------------------------------------------------

  integer function getkeyvalue_endkey(key)
    implicit none
    character(*),intent(in) :: key
    character(200) :: endkey

    integer:: nline,nkey
    character(200) :: buf

    if ( .NOT. (key(1:1) == '<') )  then
       getkeyvalue_endkey= 0
       return
    endif

    endkey=key(1:1)//'/'//key(2:)
    call input_toupper(endkey)
    nkey=len_trim(endkey)

    nline=0
    do while (.true.)

       read(nf_input,'(a120)',end=1000) buf
       nline=nline+1
       call input_delfirstspc(buf)
       call input_toupper(buf)
       if (buf(1:len_trim(buf)) == endkey(:nkey)) then
          nline=nline-1
          if (printlevel>0) then
             write(output,*) endkey(:len_trim(endkey)),' found'
             write(output,*) 'nline=',nline
          endif
          getkeyvalue_endkey=nline
          return
       endif

    enddo

1000 continue
    getkeyvalue_endkey=-1

  end  function getkeyvalue_endkey


  !-------------------------------------------------------------
  subroutine getkeyvalue_i(filename,key,var, default,status)
    implicit none
    character(*),intent(in) :: filename
    character(*),intent(in) :: key
    integer,intent(out):: var
    integer,optional :: default
    integer,optional :: status

    logical :: success
    integer:: def_var,istatus

    call input_open(filename,unusedfid(), success)
    if ( .NOT. success) then
       if (present(status)) then
          status=input_err_nofile
       endif
       !         return
       call getkeyvalue_err_exit('getkeyvalue_i: can not open a file',filename)
    endif

    if (present(default)) then
       def_var = default
    else
       def_var =  0
    endif

    var = input_integer(key,def_var,istatus)
    if (present(status)) then
       if (istatus==0) then
          status=1
       else
          status=istatus
       endif
    endif

    call input_close()

    if (istatus /= 0 .AND. .NOT. present(default) ) then
       call getkeyvalue_err_exit('getkeyvalue_i: failed to read a key',key,filename)
    endif

  end  subroutine getkeyvalue_i

  !-------------------------------------------------------------
  subroutine getkeyvalue_r8(filename,key,var, default,status)
    implicit none
    character(*),intent(in) :: filename
    character(*),intent(in) :: key
    real(8),intent(out):: var
    real(8),optional :: default
    integer,optional :: status

    logical :: success
    real(8):: def_var
    integer :: istatus

    call input_open(filename,unusedfid(), success)
    if ( .NOT. success) then
       if (present(status)) then
          status=input_err_nofile
       endif
       !         return
       call getkeyvalue_err_exit('getkeyvalue_r8: can not open a file',filename)
    endif

    if (present(default)) then
       def_var = default
    else
       def_var =  0
    endif

    var = input_real8(key,def_var,istatus)
    if (present(status)) then
       if (istatus==0) then
          status=1
       else
          status=istatus
       endif
    endif

    call input_close()

    if (istatus /= 0 .AND. .NOT. present(default) ) then
       call getkeyvalue_err_exit('getkeyvalue_r8: failed to read a key',key,filename)
    endif

  end  subroutine getkeyvalue_r8

  !-------------------------------------------------------------
  subroutine getkeyvalue_l(filename,key,var, default,status)
    implicit none
    character(*),intent(in) :: filename
    character(*),intent(in) :: key
    logical,intent(out):: var
    logical,optional :: default
    integer,optional :: status

    logical :: success
    logical:: def_var
    integer:: istatus

    call input_open(filename,unusedfid(), success)
    if ( .NOT. success) then
       if (present(status)) then
          status=input_err_nofile
       endif
       !         return
       call getkeyvalue_err_exit('getkeyvalue_l: can not open a file',filename)

    endif

    if (present(default)) then
       def_var = default
    else
       def_var = .false.
    endif

    var = input_logical(key,def_var,istatus)
    if (present(status)) then
       if (istatus==0) then
          status=1
       else
          status=istatus
       endif
    endif


    call input_close()

    if (istatus /= 0 .AND. .NOT. present(default) ) then
       call getkeyvalue_err_exit('getkeyvalue_l: failed to read a key',key,filename)
    endif


  end  subroutine getkeyvalue_l

  !-----------------------------------------------------------

  subroutine getkeyvalue_iv(filename,key,var,size, default,status)
    implicit none
    character(*),intent(in) :: filename
    character(*),intent(in) :: key
    integer,intent(in):: size
    integer,intent(out):: var(size)
    integer,optional :: default(size)
    integer,optional :: status

    logical :: success
    integer:: def_var(size),istatus,nbuf
    character:: buf(100)

    !       print *,' getkeyvalue_iv: 1',status

    call input_open(filename,unusedfid(), success)
    if ( .NOT. success) then
       if (present(status)) then
          status=input_err_nofile
       endif
       !         return
       call getkeyvalue_err_exit('getkeyvalue_iv: can not open a file',filename)

    endif

    !       print *,' getkeyvalue_iv: 2'

    if (present(default)) then
       def_var = default
    else
       def_var = 0  ! array
    endif

    !       print *,' getkeyvalue_iv: 3'

    if (present(status)) istatus=status !takao sep2006
    call input_integerv(key,size,var,def_var,istatus)
    ! cccccccccccccccccccccccc
    !      print *,' xxxxxx size=',size,' var=',var, istatus
    ! cccccccccccccccccccccccc
    if (istatus==input_err_lessvaluenum) then
       call getkeyvalue_err_exit('getkeyvalue_iv: # of values is wrong', key,filename)
    endif
    !       print *,' getkeyvalue_iv: 4'

    if (present(status)) status=istatus

    call input_close()

    if (istatus<0 .AND. .NOT. present(default) ) then
       call getkeyvalue_err_exit('getkeyvalue_iv: failed to read a key',key,filename)
    endif

  end  subroutine getkeyvalue_iv

  !-----------------------------------------------------------

  subroutine getkeyvalue_rv(filename,key,var,size, default,status)
    implicit none
    character(*),intent(in) :: filename
    character(*),intent(in) :: key
    integer,intent(in):: size
    real(8),intent(out):: var(size)
    real(8),optional :: default(size)
    integer,optional :: status

    logical :: success
    real(8):: def_var(size)
    integer:: istatus

    call input_open(filename,unusedfid(), success)
    if ( .NOT. success) then
       if (present(status)) then
          status=input_err_nofile
       endif
       !         return
       call getkeyvalue_err_exit('getkeyvalue_rv: can not open a file',filename)

    endif


    if (present(default)) then
       def_var = default
    else
       def_var = 0  ! array
    endif

    call input_real8v(key,size,var,def_var,istatus)
    ! akao dec 2005  only for (.not.present(default))
    if (istatus==input_err_lessvaluenum .AND. ( .NOT. present(default))) then
       call getkeyvalue_err_exit('getkeyvalue_rv: # of values is wrong', key,filename)
    endif

    if (present(status)) status=istatus

    call input_close()
    if (istatus<0 .AND. .NOT. present(default) ) then
       call getkeyvalue_err_exit('getkeyvalue_rv: failed to read a key',key,filename)
    endif


  end  subroutine getkeyvalue_rv

  !-----------------------------------------------------------

  subroutine getkeyvalue_c(filename,key,var, default, unit, errstop, status)
    implicit none
    character(*),intent(in) :: filename
    character(*),intent(in) :: key
    character(*),optional:: var
    character(*),optional :: default
    integer,optional :: unit
    character(*),optional:: errstop
    integer,optional :: status

    logical :: success,Lerrstop
    character(200):: def_var,tmpvar,serrstop
    integer :: istatus

    if (present(errstop)) then
       serrstop=errstop
       call input_toupper(serrstop)
       if (serrstop /= 'ON' .AND. serrstop /= 'OFF') then
          call getkeyvalue_err_exit( &
               'getkeyvalue_c: the value = (on|off) (check the source code)',key,filename);
       endif
    endif

    call input_open(filename,unusedfid(), success)
    if ( .NOT. success) then
       if (present(status)) then
          status=input_err_nofile
       endif
       !         return
       call getkeyvalue_err_exit('getkeyvalue_c: can not open a file',filename)

    endif


    if (present(default)) then
       def_var = default
    else
       def_var = ''
    endif

    tmpvar= input_character(key,def_var,istatus)
    if (printlevel>0) write(output,*) 'getkeyvalue_c: status=',istatus

    if (present(var) ) then

       if (istatus<0 .AND. .NOT. present(default) ) then
          call getkeyvalue_err_exit('getkeyvalue_c: failed to read a key',key,filename)
       endif
       var = tmpvar
       call input_close()

    else

       if (istatus /= 0 .AND. .NOT. present(unit) ) &
            call getkeyvalue_err_exit('getkeyvalue_c: failed to find a block',key,filename)

       if (present(unit) .AND. istatus == input_err_nokey) then
          if (present(status)) status= -1
          Lerrstop=.true.
          if (present(errstop)) then
             serrstop=errstop
             call input_toupper(serrstop)
             if (serrstop == 'ON' ) Lerrstop= .TRUE. 
             if (serrstop == 'OFF') Lerrstop= .FALSE. 
          endif
          if (Lerrstop)  then
             call getkeyvalue_err_exit('getkeyvalue_c: failed to find a key',key,filename);
          endif
          ! Takao Aug15 2005
          call input_close()

          return
       endif

       if (present(unit)) then
          istatus=getkeyvalue_endkey(key)
          if (present(status)) status=istatus
          if (istatus<=0) &
               call getkeyvalue_err_exit('getkeyvalue_c: failed to find the end of block',key,filename)

          ! ------ rewind and find it again
          rewind(nf_input)
          tmpvar= input_character(key,def_var,istatus)
          unit = nf_input
          ! not CLOSE
       endif

    endif


  end  subroutine getkeyvalue_c


  !-------------------------------------------------------------
end module m_keyvalue


!=============================================================

!   a sample program for 'keyvalue' module

!=============================================================

! #if  0

! program keyvalue_test_main
!   use keyvalue, only : getkeyvalue
!   implicit none

!   integer:: ikey,iflag,ivkey(10),file,i,ret
!   real(8) :: rkey,rvkey(10),rvkey2(10)
!   logical :: lflag,lkey

!   character(200):: ckey,buf
!   !        call input_substituteKeys('%d%s/#vnew.data',2,(/'%d','%s'/),
!   !     i       (/'/home/kino','psnmame   '/),ret)



!   call getkeyvalue("IN","int",ikey,default=999,status=ret )
!   write(*,*) ret,ikey

!   call getkeyvalue("IN","real",rkey ,status=ret)
!   write(*,*) ret,rkey

!   call getkeyvalue("IN","logical",lkey,status=ret )
!   write(*,*) ret,lkey

!   call getkeyvalue("IN","intv",ivkey,size=5,status=ret )
!   write(*,*) ret,ivkey

!   rvkey2(1:3)=(/10.0d0,10.0d0,10.0d0/)
!   !        call getkeyvalue("IN","realv",rvkey,size=3,default=rvkey2,status=ret )
!   call getkeyvalue("IN","realv",rvkey,size=3,status=ret )

!   write(*,*) ret,rvkey


!   call getkeyvalue("IN","<block>",unit=file,status=ret,errstop='off')
!   write(*,*) 'block = ',ret,file
!   do i=1,ret
!      read(file,'(a80)') buf
!      write(*,'(a80)') buf
!   enddo
!   close(file)


!   call getkeyvalue("IN","char",ckey,status=ret)
!   write(*,*) ret,ckey(:len_trim(ckey))

!   ! top2rx 2013.08.09 kino      stop
!   call rx( '')
! end program keyvalue_test_main


! #endif
