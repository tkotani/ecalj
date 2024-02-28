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
  character(256):: UA="kr1 6.2.0-37-generic x86_64 GNU/Linux"
  character(256):: EI="Ubuntu 22.04.3 LTS \n \l"
  character(256):: FV="GNU Fortran (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0"
  character(256):: PF=""
  character(256):: LL=""
  character(512):: GC1="commit 875831f9fa71dddd5c8a2b1c81c5610babe163e8"
  character(512):: GC2="Author: takao kotani <takaokotani@gmail.com>"
  character(512):: GC3="Date:   Tue Feb 27 20:14:37 2024 +0900"
  write(io,'(a,a,a,a)') info,trim(EI)
  write(io,'(a,a,a,a)') info,trim(FV)
  write(io,'(a,a,a)') info, trim(FFL)
  write(io,'(a,a,a)') info,'MATH: ', trim(LL)
  write(io,'(a,a,a)') info,'git: ',trim(GC1)
  write(io,'(a,a,a)') info,'   : ',trim(GC2)
  write(io,'(a,a,a)') info,'   : ',trim(GC3)
  write(io,'(a,a,a)') info,'linked at ',trim("Wed Feb 28 19:09:22 JST 2024")
end subroutine show_programinfo
! program test
! call show_programinfo(6)
! end
