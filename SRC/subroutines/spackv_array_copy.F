!! these are needed historycally. Bothe are the same as icopy.
c$$$      subroutine ii8copy(n,x,incx,y8,incy) 
c$$$      implicit none
c$$$      integer:: n,incx,incy
c$$$      integer:: x(*)
c$$$      integer:: y8(*)
c$$$      integer:: i
c$$$      do i=1,n
c$$$        y8(1+incy*(i-1))=x(1+incx*(i-1))
c$$$      enddo
c$$$      end subroutine ii8copy
c$$$!! ---
c$$$      subroutine i8icopy(n,x8,incx,y,incy) 
c$$$      implicit none
c$$$      integer:: n,incx,incy
c$$$      integer:: x8(*)
c$$$      integer:: y(*)
c$$$      integer:: i
c$$$      do i=1,n
c$$$        y(1+incy*(i-1))=x8(1+incx*(i-1))
c$$$      enddo
c$$$      end subroutine i8icopy

c$$$      subroutine spackv_array_copy_r8_r8(c,v,n,j,s)
c$$$      implicit none
c$$$      character,intent(in):: c
c$$$      real(8):: v(n)
c$$$      integer,intent(in):: j,n
c$$$      real(8):: s(n*j)
c$$$      integer:: i
c$$$      if (c=='p') then
c$$$        do i=1,n
c$$$          v(i)=s(n*(j-1)+i)
c$$$        enddo
c$$$      else
c$$$        do i=1,n
c$$$          s(n*(j-1)+i)=v(i)
c$$$        enddo
c$$$      endif
c$$$      end
c$$$
c$$$      subroutine spackv_array_copy_i8_i(c,v,n,j,s)
c$$$      implicit none
c$$$      character,intent(in):: c
c$$$      integer:: v(n)
c$$$      integer,intent(in):: j,n
c$$$      integer:: s(n*j)
c$$$      integer:: i
c$$$      if (c=='p') then
c$$$        do i=1,n
c$$$          v(i)=s(n*(j-1)+i)
c$$$        enddo
c$$$      else
c$$$        do i=1,n
c$$$          s(n*(j-1)+i)=v(i)
c$$$        enddo
c$$$      endif
c$$$      end


c$$$      subroutine spacks_copy(c,struc,i1,i2,name,istr)
c$$$C     - Pack/unpack string element in a struc, element referenced by name
c$$$C     ----------------------------------------------------------------
c$$$C     i   lpack   0 for unpack, 1 for pack, 2 to return parameters
c$$$C     i   is1..is2 range of species is for which to pack
c$$$C     For structure without species, set both to 0
c$$$C     ----------------------------------------------------------------
c$$$
c$$$      implicit none
c$$$      character,intent(in):: c
c$$$      integer,intent(in):: i1,i2,istr
c$$$c      real(8):: struc
c$$$      character(8):: struc
c$$$      character(8):: name(i1:i2) ! use only name(istr)
c$$$      character(11):: thisfunc='spacks_copy'
c$$$
c$$$c      write(*,*)'spacks_copy>',c,struc,i1,i2,name,istr
c$$$
c$$$      if (c=='p') then
c$$$c        call s8tor8(name(istr),struc)
c$$$        struc=name(istr)
c$$$      elseif (c=='u') then
c$$$c        call r8tos8(struc,name(istr))
c$$$        name(istr)=struc
c$$$      else
c$$$        write(*,*) 'error in ',thisfunc
c$$$      endif
c$$$      end subroutine spacks_copy

c$$$      subroutine lsets_bitop_i8(v_iv,n,lval,mask)
c$$$C     - Set a logical bits to a structure entry
c$$$C     ----------------------------------------------------------------
c$$$C     i Inputs
c$$$C     i   lval is a logical T or F, of length of sname entry
c$$$C     i   mask should be an integer multiple of 2.
c$$$C     i        Only the lowest bit of mask is used.
c$$$C     o Outputs
c$$$C     o  struc element corresponding to label 'name' is modified.
c$$$C     o        The mask bit of that entry is set to lval.
c$$$C     ----------------------------------------------------------------
c$$$c      use m_struc_def
c$$$c      use m_struc_func
c$$$      implicit none
c$$$      integer:: v_iv(n)
c$$$      logical,intent(in):: lval
c$$$      integer,intent(in):: mask
c$$$      integer,intent(in):: n
c$$$      integer:: newval,lmask
c$$$      integer:: bitor, bitand ,bitlow
c$$$      integer:: i
c$$$      lmask=bitlow(mask)
c$$$      do i=1,n
c$$$        if (lval) then
c$$$          newval = ior(int(v_iv(i)),lmask)
c$$$        else
c$$$          newval = v_iv(i) - iand(int(v_iv(i)),lmask)
c$$$        endif
c$$$        v_iv(i) = newval
c$$$      enddo
c$$$      end subroutine lsets_bitop_i8

C$$$#if 0
C$$$      program main
C$$$      implicit none
C$$$      integer,parameter:: n=10
C$$$      integer(8):: v(n)
C$$$      integer:: s(n)
C$$$      integer:: i
C$$$      do i=1,n
C$$$         v(i)=i+10
C$$$         s(i)=i
C$$$      enddo
C$$$      do i=1,n
C$$$      call spackv_array_copy_i8_i('p',v(i),1,1,-i)
C$$$      enddo
C$$$      write(*,*) v
C$$$      write(*,*) s
C$$$      end
C$$$#endif
