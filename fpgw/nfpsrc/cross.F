C   $Id: cross.f,v 1.2 1998/04/28 18:45:23 coa Exp $
c --- cross
      subroutine cross(a,b,c)
      implicit real*8 (a-h,p-z), integer(o)
      dimension a(3),b(3),c(3)
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
      return
      end

c --- tripl
      doubleprecision function tripl(a,b,c)
      implicit real*8 (a-h,p-z), integer(o)
c     dimension a(3),b(3),c(3)
c ... g77 needs this rewriting for opt > O1
      dimension a(0:2),b(0:2),c(0:2)
c     tripl=a(1)*b(2)*c(3)+a(2)*b(3)*c(1)+a(3)*b(1)*c(2)
c    .     -a(3)*b(2)*c(1)-a(2)*b(1)*c(3)-a(1)*b(3)*c(2)
c ... g77 needs this rewriting
      tmp = 0.d0
      do i = 0,2
        j = mod(i + 1, 3)
        k = mod(i + 2, 3)
        tmp = tmp + a(i) * (b(j)*c(k) - b(k)*c(j))
      enddo
      tripl = tmp
      end

c --- sprod
      doubleprecision function sprod(a,b)
      implicit real*8 (a-h,p-z), integer(o)
      dimension a(3),b(3)
      sprod=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      return
      end
