      subroutine hunti(xa,n,x,iprm,low)
C- Brackets a value within an ordered array of integer points
C ----------------------------------------------------------------
Ci Inputs
Ci   xa  :array of points
Ci   n   :size of array
Ci   x   :value to bracket
Ci   iprm:permutation table by which array xa is ordered
Ci        iprm(1) <= 0 => assumes iprm(i) = i; iprm not referenced
Ci   low : initial guess for output low
Co Outputs
Co   low : xa(low) < x <= xa(low+1)
Cu Updates
Cu   13 Sep 01 Handle case n=0
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer n,low,iprm(n)
      integer xa(n),x
C Local variables
      integer inc,jhi,jm
      logical ascnd,liprm
      integer xn

      if (n .eq. 0) then
        low = 0
        return
      endif

      liprm = iprm(1) .gt. 0
      if (liprm) then
        ascnd = xa(iprm(n)) .gt. xa(iprm(1))
      else
        ascnd = xa(n) .gt. xa(1)
      endif
      if (low.le.0 .or. low.gt.n) then
        low = 0
        jhi = n+1
        goto 3
      endif
      inc = 1
      if (liprm) then
        xn = xa(iprm(low))
      else
        xn = xa(low)
      endif
      if (x.ge.xn .eqv. ascnd) then
    1   jhi = low+inc
        if (jhi .gt. n) then
          jhi = n+1
        else
          if (liprm) then
            xn = xa(iprm(jhi))
          else
            xn = xa(jhi)
          endif
          if (x.ge.xn .eqv. ascnd) then
            low = jhi
            inc = inc+inc
            goto 1
          endif
        endif
      else
        jhi = low
    2   low = jhi-inc
        if (low .lt. 1) then
          low = 0
        else
          if (liprm) then
            xn = xa(iprm(low))
          else
            xn = xa(low)
          endif
          if (x.lt.xn .eqv. ascnd) then
            jhi = low
            inc = inc+inc
            goto 2
          endif
        endif
      endif
    3 if (jhi-low .eq. 1) then
C   ... Find the first of values equal to x
    4   continue
        if (low .gt. 1) then
          if (liprm) then
            xn = xa(iprm(low-1))
          else
            xn = xa(low-1)
          endif
          if (xn .eq. x) then
            low = low-1
            goto 4
          endif
        endif
C     ... if xa(low) = x, decrement low
        if (low .ge. 1) then
          if (liprm) then
            xn = xa(iprm(low))
          else
            xn = xa(low)
          endif
          if (xn .eq. x) low = low-1
        endif
        return
      endif
      jm = (jhi+low)/2
      if (liprm) then
        xn = xa(iprm(jm))
      else
        xn = xa(jm)
      endif
      if (x.gt.xn .eqv. ascnd) then
        low = jm
      else
        jhi = jm
      endif
      goto 3
      end

