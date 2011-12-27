# 1 "show_programinfo.F"
# 1 "<組み込み>"
# 1 "<コマンドライン>"
# 1 "show_programinfo.F"
      subroutine show_programinfo(io)
      implicit none
      integer:: io
      character(17),parameter:: info='HOST_INFORMATION '







       write(io,'(a,a,a)') info,'platform: ',"ifort"
       write(io,'(a,a,a)') info,'compiler version: ',"ifort version 12.1.0"
       write(io,'(a,a,a)') info,'FFLAGS (<=120): ',
     ."-O1 -xHost -WB -cpp -DHASIARGC -DHASGETARG -DFDATE -DHASGETTIMEOFDAY -DHASGETENVQQ -DHASSETENVQQ -DHASIFPORT -module MOD"
       write(io,'(a,a,a)') info,'LIBLOC (<=120): ',
     ."-lfftw3 /opt/acml4.4.0/ifort64/lib/libacml.a"
       write(io,'(a,a,a)') info,'uname -a (<=120): ',
     ."Linux maro005 3.0.0-14-generic #23-Ubuntu SMP Mon Nov 21 20:28:43 UTC 2011 x86_64 x86_64 x86_64 GNU/Linux"
       write(io,'(a,a,a)') info,'/etc/issue: ',"Ubuntu 11.10 \n \l"
       write(io,'(a,a,a)') info,'git branch: ',"refs/heads/kino-alloc6"
       write(io,'(a,a,a)') info,'git commit: ',"0f3d7b2c16fca9a029ac28bf0d8d83064d76ce5a"
       write(io,'(a,a,a)') info,'linked at: ',"Tue Dec 27 17:51:45 JST 2011"



       end subroutine show_programinfo

C       program test
C       call show_programinfo(6)
C       end 
