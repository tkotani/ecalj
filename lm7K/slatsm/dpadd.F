      subroutine dpadd(dest,src,n1,n2,fac)
C- Adds fac*vector into another vector
      double precision dest(1),src(1),fac

      if (fac .eq. 1) then
        do  10  i = n1, n2
   10   dest(i) = dest(i) + src(i)
      else
        do  20  i = n1, n2
   20   dest(i) = dest(i) + fac*src(i)
      endif
      end

