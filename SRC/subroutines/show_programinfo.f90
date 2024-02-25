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
  character(512):: GC1="commit 3ce8043ade7c346cf1a72a412dc87f7e1f2cd3f7"
  character(512):: GC2="Author: Takao Kotani <takaokotani@gmail.com>"
  character(512):: GC3="Date:   Sat Feb 24 22:00:30 2024 +0900"
  write(io,'(a,a,a,a)') info,trim(EI)
  write(io,'(a,a,a,a)') info,trim(FV)
  write(io,'(a,a,a)') info, trim(FFL)
  write(io,'(a,a,a)') info,'MATH: ', trim(LL)
  write(io,'(a,a,a)') info,'git: ',trim(GC1)
  write(io,'(a,a,a)') info,'   : ',trim(GC2)
  write(io,'(a,a,a)') info,'   : ',trim(GC3)
  write(io,'(a,a,a)') info,'linked at ',trim("Sun Feb 25 15:31:32 JST 2024")
end subroutine show_programinfo
! program test
! call show_programinfo(6)
! end
