c$$$  Kino Sep2014. 
c$$$  m_commandline.Fを添付しておきます。テストは既にしていてclassを使い
c$$$      if( mpi__root ) then
c$$$         if (cmdline%get_f()>0) then
c$$$             ixc=cmdline%mode
c$$$         else
c$$$         read(5,*) ixc !c     call readin5(ixc,iqxini,iqxend)
c$$$         endif
c$$$      end if
c$$$  としていました。program毎に変える場合はこのファイルで定義したclassから継承して作る予定でした。

       module m_commandline
       implicit none
       type t_commandline
        integer:: mode=0
        integer:: narg=0
        integer:: nz =0
         contains
         procedure :: get_s => commandline_process
         procedure :: get_f => commandline_get
       end type t_commandline
       
       contains
       
       
       integer function commandline_get(p)
       implicit none
       class(t_commandline)::p
       call commandline_process(p)
       commandline_get = p%narg
       end function commandline_get
       
       subroutine commandline_process(p)
       implicit none
       class(t_commandline)::p
       
       integer:: n,i,errflag,ios
       character(200):: str,str2
       
           n=command_argument_count()
           p%narg=n 
       
           errflag=0
           i=0
           do while (.true.)
             i=i+1
              if (i>n) exit
              call get_command_argument(i, str)
              select case(str)
              case('-mode')
                 i=i+1
                 if (i>n) then 
                   write(6,*) trim(str), ' must have 2 arguments.'
                   errflag=errflag+1
                   exit
                 endif
                 call get_command_argument(i, str)
                 read(str,*,iostat=ios) p%mode
                 if (ios.ne.0) then 
                    write(6,*) 'conversion error in '
                    write(6,*) trim(str),' ',trim(str2)
                    errflag=errflag+1
                 endif
              case('-nz')
                 i=i+1
                 if (i>n) then
                   write(6,*) trim(str), ' must have 2 arguments.'
                   errflag=errflag+1
                   exit
                 endif
                 call get_command_argument(i, str)
                 read(str,*,iostat=ios) p%nz
                 if (ios.ne.0) then 
                    write(6,*) 'conversion error in '
                    write(6,*) trim(str),' ',trim(str2)
                    errflag=errflag+1
                 endif

              case default
                 write(6,*)'unknown option', trim(str)
                 errflag=errflag+1
              end select
           enddo
       
           if (errflag.ne.0) then 
               call rx('unknown arguments')
           endif
              
       end subroutine commandline_process
       
       end module m_commandline
       
#if 0
       program test
       use m_commandline
       
       type(t_commandline):: cmdline
       integer:: i
       
       write(6,*) cmdline%get_f(),cmdline%mode
       
       write(6,*) cmdline%narg,cmdline%mode
       
       end program test
       
       
       subroutine rx(str)
       character(*) str
       write(6,*)  trim(str)
       stop 
       end subroutine rx
#endif
