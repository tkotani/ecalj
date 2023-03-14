subroutine huntx(xa,n,x,low)  !- Brackets a value within an ordered array of points
  ! ----------------------------------------------------------------
  !i Inputs
  !i   n   :size of array
  !i   xa  :array of points, in ascending or descending order
  !i   x   :value to bracket
  !i   low : initial guess for output low
  !o Outputs
  !o   ... if xa is ordered in ascending order:
  !o   low : xa(low) < x <= xa(low+1)
  !o       : when x cannot be bracketed,
  !o       : low = 0 if x<=xa(1)
  !o       : low = n if xa(n)<x
  !o   ... if xa is ordered in descending order:
  !o   low : xa(low) > x >= xa(low+1)
  !o       : when x cannot be bracketed,
  !o       : low = 0 if x>xa(1)
  !o       : low = n if xa(n)>=x
  !u Updates
  !u   29 Jul 04 Handle special case xa(1) .eq. xa(n)
  implicit none
  integer :: n,low ,inc,jhi,jm
  double precision :: xa(n),x,xn
  logical :: ascnd 
  ! ... Pathological case: all points equal (treat as ascending)
  if (xa(n) == xa(1)) then
     if (xa(1) >= x) then
        low = 0
     else
        low = n
     endif
     return
  endif
  ascnd = xa(n) .gt. xa(1)
  if (low <= 0 .OR. low > n) then
     low = 0
     jhi = n+1
     goto 3
  endif
  inc = 1
  xn = xa(low)
  if (x >= xn .eqv. ascnd) then
1    jhi = low+inc
     if (jhi > n) then
        jhi = n+1
     else
        xn = xa(jhi)
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
        xn = xa(low)
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
        xn = xa(low-1)
        if (xn == x) then
           low = low-1
           goto 4
        endif
     endif
     !   ... if xa(low) = x, decrement low
     if (low >= 1) then
        xn = xa(low)
        if (xn == x) low = low-1
     endif
     return
  endif
  jm = (jhi+low)/2
  xn = xa(jm)
  if (x > xn .eqv. ascnd) then
     low = jm
  else
     jhi = jm
  endif
  goto 3
end subroutine huntx
