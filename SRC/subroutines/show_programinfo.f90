# 1 "../subroutines/show_programinfo.template"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "../subroutines/show_programinfo.template"
subroutine show_programinfo(io)
  implicit none
  integer:: io
  character(6),parameter:: info='INFO: '
  character(256):: FFL=""
  character(256):: UA="tt480s 5.4.0-170-generic x86_64 GNU/Linux"
  character(256):: EI="Ubuntu 20.04.6 LTS \n \l"
  character(256):: FV="GNU Fortran (Ubuntu 9.4.0-1ubuntu1~20.04.2) 9.4.0"
  character(256):: PF=""
  character(256):: LL=""
  character(512):: GC1="commit 540d6ab40a3a7a82c1bc99678d54238fd49043da"
  character(512):: GC2="Author: Takao Kotani <takaokotani@gmail.com>"
  character(512):: GC3="Date:   Sun Feb 11 15:53:18 2024 +0900"
  write(io,'(a,a,a,a)') info,trim(EI)
  write(io,'(a,a,a,a)') info,trim(FV)
  write(io,'(a,a,a)') info, trim(FFL)
  write(io,'(a,a,a)') info,'MATH: ', trim(LL)
  write(io,'(a,a,a)') info,'git: ',trim(GC1)
  write(io,'(a,a,a)') info,'   : ',trim(GC2)
  write(io,'(a,a,a)') info,'   : ',trim(GC3)
  write(io,'(a,a,a)') info,'linked at ',trim("Sun Feb 11 16:13:56 JST 2024")
end subroutine show_programinfo
! program test
! call show_programinfo(6)
! end
