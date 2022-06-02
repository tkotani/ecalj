subroutine dvheap(m,n,vecs,iprm,tol,opts)
  !- Heapsort array of double-precision vectors
  ! ----------------------------------------------------------------
  !i Inputs
  !i   m     :length of each vector
  !i   n     :number of vectors to be sorted
  !i   vecs  :the array vectors, dimensioned (m,n)
  !i   tol   :numbers differing by less than tol are treated as equal.
  !i   opts  :ones digit
  !i           0 vecs returned sorted.
  !i             NB: in this case, vecs must be dimensioned 2*m*n.
  !i           1 vecs is unchanged; only iprm is returned.
  !i          tens digit
  !i           0 vecs sorted
  !i           1 vecs sorted by increasing length
  !i         hundreds digit
  !i           1 equal vectors preserve their original order
  !i
  !o Outputs
  !o   iprm  :a permutation table that sorts array 'vecs'
  !o   vecs  : may be changed, depending on opts
  !u Updates
  !u   02 Sep 02 Added 100s digit switch
  ! ----------------------------------------------------------------
  !     implicit none
  integer :: m,n,iprm(n),opts
  double precision :: vecs(m,n),tol
  double precision :: di,dl
  integer :: l,ir,irra,i,j,mm,i1,i2,nn
  logical :: norm

  do  ir = 1, n
     iprm(ir) = ir
  enddo
  if (n <= 1) return
  norm = mod(opts/10,10) .ne. 0
  l = n/2+1
  ir = n

  ! --- For each l = n/2+1, 1, -1 do ---
10 continue
  ! ... "Hiring phase"
  if (l > 1) then
     l = l-1
     irra = iprm(l)
     ! ... "Retirement-and-promotion phase"
  else
     irra = iprm(ir)
     iprm(ir) = iprm(1)
     !*       call awrit3('ir%i: %n:1i',' ',180,6,ir,n,iprm)
     ir = ir-1
     if (ir == 1) then
        iprm(1) = irra
        !*         call awrit2('exit %n:1i',' ',180,6,n,iprm)
        goto 100
     endif
  endif

  ! ... Setup to sift down element irra to proper level
  i = l
  j = l+l

  ! --- Do while j .le. ir ---
20 if (j <= ir) then

     !   ... Increment j if vecs(iprm(j+1)) > vecs(iprm(j))
     if (j < ir) then
        if (norm) then
           di = 0d0
           dl = 0d0
           do   mm = 1, m
              dl = dl + vecs(mm,iprm(j))**2
              di = di + vecs(mm,iprm(j+1))**2
           enddo
           dl = dsqrt(dl)
           di = dsqrt(di)
           if (abs(di-dl) > tol) then
              if (di-dl > tol) j = j+1
           endif
        else
           do  26  mm = 1, m
              if (abs(vecs(mm,iprm(j+1))-vecs(mm,iprm(j))) <= tol) &
                   goto 26
              if (vecs(mm,iprm(j+1))-vecs(mm,iprm(j)) > tol) j = j+1
              goto 28
26         enddo
28         continue
        endif
     endif

     !   ... Demote rra to its level
     if (norm) then
        di = 0d0
        dl = 0d0
        do    mm = 1, m
           dl = dl + vecs(mm,irra)**2
           di = di + vecs(mm,iprm(j))**2
        enddo
        dl = dsqrt(dl)
        di = dsqrt(di)
        if (di-dl > tol) then
           iprm(i) = iprm(j)
           !*           call awrit4('%i,%i: %n:1i',' ',180,6,i,j,n,iprm)
           i = j
           j = j+j
           !     ... This is rra's level; set j to terminate the sift-down
        else
           j = ir+1
        endif
     else
        do  36  mm = 1, m
           !     ...   Skip over equal elements
           if (abs(vecs(mm,iprm(j))-vecs(mm,irra)) <= tol) goto 36
           if (vecs(mm,iprm(j))-vecs(mm,irra) > tol) then
              iprm(i) = iprm(j)
              !*             call awrit4('%i,%i: %n:1i',' ',180,6,i,j,n,iprm)
              i = j
              j = j+j
              !     ... This is rra's level; set j to terminate the sift-down
           else
              j = ir+1
           endif
           goto 38
36      enddo
        !     ... Case rra = vec(iprm(j))
        j = ir+1
38      continue
     endif
     go to 20
  endif
  ! ... Put rra into its slot
  iprm(i) = irra
  !*     call awrit3('%i: %n:1i',' ',180,6,i,n,iprm)
  go to 10

  ! --- For equal vectors, restore original ordering ---
100 continue
  if (mod(opts/100,10) == 0) goto 200
  i2 = 0
  ! ... Find i1,i2 = next range of equal numbers
110 i1 = i2+1
  if (i1 > n) goto 200
120 i2 = i2+1
  if (i2 > n) goto 130
  if (norm) then
     di = 0d0
     dl = 0d0
     do    mm = 1, m
        dl = dl + vecs(mm,iprm(i1))**2
        di = di + vecs(mm,iprm(i2))**2
     enddo
     dl = dsqrt(dl)
     di = dsqrt(di)
     if (abs(di-dl) > tol) goto 130
  else
     do mm = 1, m
        if (abs(vecs(mm,iprm(i2))-vecs(mm,iprm(i1))) > tol) goto 130
     enddo
  endif
  ! ... vec(i1) = vec(i2) ; imcrement i2 and try again
  goto 120

  ! --- Sort iprm(i1)..iprm(i2) ---
130 continue
  i2 = i2-1
  i1 = i1-1
  nn = i2-i1
  if (nn <= 1) goto 110
  l = nn/2+1
  ir = nn

  ! ... For each l = (i2-i1+1)/2+1, 1, -1 do
140 continue
  ! ... "Hiring phase"
  if (l > 1) then
     l = l-1
     irra = iprm(l+i1)
     ! ... "Retirement-and-promotion phase"
  else
     irra = iprm(ir+i1)
     iprm(ir+i1) = iprm(1+i1)
     ir = ir-1
     if (ir == 1) then
        iprm(1+i1) = irra
        goto 110
     endif
  endif

  ! ... Setup to sift down element irra to proper level ...
  i = l
  j = l+l

  ! ... Do while j .le. ir ...
150 if (j <= ir) then

     !   ... Increment j if iprm(j+i11) > iprm(j+i1))
     if (j < ir) then
        if (iprm(j+i1) < iprm(j+1+i1)) j = j+1
     endif
     !   ... Demote irra to its level
     if (irra < iprm(j+i1)) then
        iprm(i+i1) = iprm(j+i1)
        i = j
        j = j+j
        !   ... This is irra's level; set j to terminate the sift-down
     else
        j = ir+1
     endif
     go to 150
  endif
  ! ... Put rra into its slot
  iprm(i+i1) = irra
  go to 140

  ! --- Sort vecs ---
200 continue
  if (mod(opts,10) == 0) &
       call dvprm(m,n,vecs,vecs(1,n+1),iprm, .TRUE. )

end subroutine dvheap
subroutine dvprm(m,n,vecs,wk,iprm,lopt)
  !- Permute an array of double precision vectors according to iprm
  ! ----------------------------------------------------------------
  !i Inputs
  !i   m     :length of each vector
  !i   n     :number of vectors to be sorted
  !i   vecs  :the array vectors, dimensioned (m,n)
  !i   iprm  :a permutation table by which array vecs is reordered
  !i   lopt  :if T, copy wk back to vecs.
  !i Inputs/Outputs
  !o   wk    :returns vecs in permuted order
  !o   vecs  :wk may be copied back into vecs, depending on lopt.
  ! ----------------------------------------------------------------
  !     implicit none
  integer :: m,n,iprm(n)
  logical :: lopt
  double precision :: vecs(m,n),wk(m,n)
  integer :: i,j,k

  do  10  i = 1, n
     k = iprm(i)
     do   j = 1, m
        wk(j,i) = vecs(j,k)
     enddo
10 enddo
  if (lopt) call dpcopy(wk,vecs,1,n*m,1d0)
end subroutine dvprm

