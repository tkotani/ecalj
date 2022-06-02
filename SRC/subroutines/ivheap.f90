subroutine ivheap(m,n,vecs,iprm,opts)
  !- Heapsort array of integer vectors
  ! ----------------------------------------------------------------
  !i Inputs
  !i   vecs(m,n): n vectors of length m are to be sorted
  !i   iprm: an integer work array of dimension n, or if vecs returned
  !i        in sorted order, an array of length m*n
  !i   opts: ones digit
  !i           0 vecs returned sorted.
  !i           1 vecs is unchanged; only iprm is returned
  !i         tens digit
  !i           0 vecs sorted
  !i           1 vecs sorted by increasing length
  !i         hundreds digit
  !i           1 equal vectors preserve their original order
  !i
  !o Outputs
  !o   iprm a permutation table that sorts array 'vecs'
  !o   vecs may be changed, depending on opts
  ! ----------------------------------------------------------------
  !     implicit none
  integer :: m,n,iprm(n),opts
  integer :: vecs(m,n*2)
  integer :: di,dl
  integer :: l,ir,irra,i,j,mm,i1,i2,nn
  logical :: norm
  integer,allocatable:: wk(:,:)
  do   ir = 1, n
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
           di = 0
           dl = 0
           do    mm = 1, m
              dl = dl + vecs(mm,iprm(j))**2
              di = di + vecs(mm,iprm(j+1))**2
           enddo
           if (di-dl > 0) j = j+1
        else
           do  26  mm = 1, m
              if (abs(vecs(mm,iprm(j+1))-vecs(mm,iprm(j))) <= 0) goto 26
              if (vecs(mm,iprm(j+1))-vecs(mm,iprm(j)) > 0) j = j+1
              goto 28
26         enddo
28         continue
        endif
     endif

     !   ... Demote rra to its level
     if (norm) then
        di = 0
        dl = 0
        do  34  mm = 1, m
           dl = dl + vecs(mm,irra)**2
           di = di + vecs(mm,iprm(j))**2
34      enddo
        if (di-dl > 0) then
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
           if (abs(vecs(mm,iprm(j))-vecs(mm,irra)) <= 0) goto 36
           if (vecs(mm,iprm(j))-vecs(mm,irra) > 0) then
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
     di = 0
     dl = 0
     do  mm = 1, m
        dl = dl + vecs(mm,iprm(i1))**2
        di = di + vecs(mm,iprm(i2))**2
     enddo
     if (di-dl > 0) goto 130
  else
     do    mm = 1, m
        if (abs(vecs(mm,iprm(i2))-vecs(mm,iprm(i1))) > 0) goto 130
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
  if (mod(opts,10) == 0) then
     allocate(wk(m,n))
     call ivprm(m,n,vecs,wk,iprm,.true.) !wk jan2015
     deallocate(wk)
  endif

end subroutine ivheap
subroutine ivprm(m,n,vecs,wk,iprm,lopt)
  !- Permute an array of integer vectors according to iprm
  ! ----------------------------------------------------------------
  !i Inputs
  ! o  vecs(m,n): n vectors of length m are to be permuted
  !i   iprm: a table of permutation indices to array vecs
  !i   wk:  a work array of length m
  !i   lopt: T, copy wk back to vecs.
  ! ----------------------------------------------------------------
  !     implicit none
  integer :: m,n,iprm(n)
  logical :: lopt
  integer :: vecs(m,n),wk(m,n)
  integer :: i,j,k

  do  10  i = 1, n
     k = iprm(i)
     do  j = 1, m
        wk(j,i) = vecs(j,k)
     enddo
10 enddo
  if (lopt) call icopy(n*m,wk,1,vecs,1)
end subroutine ivprm

