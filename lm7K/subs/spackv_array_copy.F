
      subroutine spackv_array_copy_r8_r8(c,v,n,j,s)
      implicit none
      character,intent(in):: c
      real(8):: v(n)
      integer,intent(in):: j,n
      real(8):: s(n*j)
      integer:: i
      if (c=='p') then
      do i=1,n
      v(i)=s(n*(j-1)+i)
      enddo
      else
      do i=1,n
      s(n*(j-1)+i)=v(i)
      enddo
      endif
      end

      subroutine spackv_array_copy_i8_i(c,v,n,j,s)
      implicit none
      character,intent(in):: c
      integer(8):: v(n)
      integer,intent(in):: j,n
      integer:: s(n*j)
      integer:: i
      if (c=='p') then
      do i=1,n
      v(i)=s(n*(j-1)+i)
      enddo
      else
      do i=1,n
      s(n*(j-1)+i)=v(i)
      enddo
      endif
      end 


      subroutine spacks_copy(c,struc,i1,i2,name,istr)
C     - Pack/unpack string element in a struc, element referenced by name
C     ----------------------------------------------------------------
C     i   lpack   0 for unpack, 1 for pack, 2 to return parameters
C     i   is1..is2 range of species is for which to pack
C     For structure without species, set both to 0
C     ----------------------------------------------------------------

      implicit none
      character,intent(in):: c
      integer,intent(in):: i1,i2,istr
      real(8):: struc
      character(8):: name(i1:i2) ! use only name(istr)
      character(11):: thisfunc='spacks_copy'

c      write(*,*)'spacks_copy>',c,struc,i1,i2,name,istr

      if (c=='p') then
         call s8tor8(name(istr),struc)
      elseif (c=='u') then
         call r8tos8(struc,name(istr))
      else
         write(*,*) 'error in ',thisfunc
      endif
      end subroutine spacks_copy

      subroutine ii8copy(n,x,incx,y8,incy)
      implicit none
      integer:: n,incx,incy
      integer:: x(*)
      integer(8):: y8(*)
      integer:: i
      do i=1,n
         y8(1+incy*(i-1))=x(1+incx*(i-1))
      enddo
      end subroutine ii8copy

      subroutine i8icopy(n,x8,incx,y,incy)
      implicit none
      integer:: n,incx,incy
      integer(8):: x8(*)
      integer:: y(*)
      integer:: i
      do i=1,n
         y(1+incy*(i-1))=x8(1+incx*(i-1))
      enddo
      end subroutine i8icopy

      subroutine lsets_bitop_i8(v_iv,n,lval,mask)
C     - Set a logical bits to a structure entry
C     ----------------------------------------------------------------
C     i Inputs
C     i   lval is a logical T or F, of length of sname entry
C     i   mask should be an integer multiple of 2.
C     i        Only the lowest bit of mask is used.
C     o Outputs
C     o  struc element corresponding to label 'name' is modified.
C     o        The mask bit of that entry is set to lval.
C     ----------------------------------------------------------------
      use m_struc_def
      use m_struc_func
      implicit none
      integer(8):: v_iv(n)
      logical,intent(in):: lval
      integer,intent(in):: mask
      integer,intent(in):: n
      integer:: newval,lmask
      integer:: bitor, bitand ,bitlow
      integer:: i

      lmask=bitlow(mask)

      do i=1,n

         if (lval) then
            newval = ior(int(v_iv(i)),lmask)
         else
            newval = v_iv(i) - iand(int(v_iv(i)),lmask)
         endif

         v_iv(i) = newval

      enddo

      end subroutine lsets_bitop_i8

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
