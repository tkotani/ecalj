      integer function idnear(n,da,dx,incx)
C- Finds the index of element closest to specified value
      double precision dx(1),dmax,da
      integer i,incx,ix,n
c
      idnear = 0
      if ( n .lt. 1 ) return
      idnear = 1
      if ( n .eq. 1 ) return
      ix = 1
      dmax = dabs(dx(1)-da)
      if (dmax .eq. 0) return
      ix = ix + incx
      do  10  i = 2, n
        if (dabs(dx(ix)-da) .lt. dmax) then
          idnear = i
          dmax = dabs(dx(ix)-da)
          if (dmax .eq. 0) return
        endif
        ix = ix + incx
   10 continue
      end

