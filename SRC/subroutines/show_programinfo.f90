# 0 "../subroutines/show_programinfo.template"
# 0 "<built-in>"
# 0 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 0 "<command-line>" 2
# 1 "../subroutines/show_programinfo.template"
subroutine show_programinfo(io)
  implicit none
  integer:: io
  character(6),parameter:: info='INFO: '
  character(256):: FFL=""
  character(256):: UA="tt14 6.1.0-1035-oem x86_64 GNU/Linux"
  character(256):: EI="Ubuntu 22.04.4 LTS \n \l"
  character(256):: FV="GNU Fortran (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0"
  character(256):: PF=""
  character(256):: LL=""
  character(512):: GC1="commit 3443a3d2642e63656ff5b108a0a27d28cb4bdcf2"
  character(512):: GC2="Author: Takao Kotani <takaokotani@gmail.com>"
  character(512):: GC3="Date:   Tue Mar 12 20:30:05 2024 +0900"
  write(io,'(a,a,a,a)') info,trim(EI)
  write(io,'(a,a,a,a)') info,trim(FV)
  write(io,'(a,a,a)') info, trim(FFL)
  write(io,'(a,a,a)') info,'MATH: ', trim(LL)
  write(io,'(a,a,a)') info,'git: ',trim(GC1)
  write(io,'(a,a,a)') info,'   : ',trim(GC2)
  write(io,'(a,a,a)') info,'   : ',trim(GC3)
  write(io,'(a,a,a)') info,'linked at ',trim("Sat Mar 16 18:34:53 JST 2024")
end subroutine show_programinfo
! program test
! call show_programinfo(6)
! end
