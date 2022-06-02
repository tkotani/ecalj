subroutine hunti(xa,n,x,iprm,low)
  !- Brackets a value within an ordered array of integer points
  ! ----------------------------------------------------------------
  !i Inputs
  !i   xa  :array of points
  !i   n   :size of array
  !i   x   :value to bracket
  !i   iprm:permutation table by which array xa is ordered
  !i        iprm(1) <= 0 => assumes iprm(i) = i; iprm not referenced
  !i   low : initial guess for output low
  !o Outputs
  !o   low : xa(low) < x <= xa(low+1)
  !u Updates
  !u   13 Sep 01 Handle case n=0
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  integer :: n,low,iprm(n)
  integer :: xa(n),x
  ! Local variables
  integer :: inc,jhi,jm
  logical :: ascnd,liprm
  integer :: xn

  if (n == 0) then
     low = 0
     return
  endif

  liprm = iprm(1) .gt. 0
  if (liprm) then
     ascnd = xa(iprm(n)) .gt. xa(iprm(1))
  else
     ascnd = xa(n) .gt. xa(1)
  endif
  if (low <= 0 .OR. low > n) then
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
  if (x >= xn .eqv. ascnd) then
1    jhi = low+inc
     if (jhi > n) then
        jhi = n+1
     else
        if (liprm) then
           xn = xa(iprm(jhi))
        else
           xn = xa(jhi)
        endif
        if (x >= xn .eqv. ascnd) then
           low = jhi
           inc = inc+inc
           goto 1
        endif
     endif
  else
     jhi = low
2    low = jhi-inc
     if (low < 1) then
        low = 0
     else
        if (liprm) then
           xn = xa(iprm(low))
        else
           xn = xa(low)
        endif
        if (x < xn .eqv. ascnd) then
           jhi = low
           inc = inc+inc
           goto 2
        endif
     endif
  endif
3 if (jhi-low == 1) then
     !   ... Find the first of values equal to x
4    continue
     if (low > 1) then
        if (liprm) then
           xn = xa(iprm(low-1))
        else
           xn = xa(low-1)
        endif
        if (xn == x) then
           low = low-1
           goto 4
        endif
     endif
     !     ... if xa(low) = x, decrement low
     if (low >= 1) then
        if (liprm) then
           xn = xa(iprm(low))
        else
           xn = xa(low)
        endif
        if (xn == x) low = low-1
     endif
     return
  endif
  jm = (jhi+low)/2
  if (liprm) then
     xn = xa(iprm(jm))
  else
     xn = xa(jm)
  endif
  if (x > xn .eqv. ascnd) then
     low = jm
  else
     jhi = jm
  endif
  goto 3
end subroutine hunti

