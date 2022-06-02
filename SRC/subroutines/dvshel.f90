subroutine dvshel(m,n,vecs,iprm,lopt)
  !- Shell sort of a array of double precision vectors
  ! ----------------------------------------------------------------
  !i Inputs
  !i   m     :length of each vector
  !i   n     :number of vectors to be sorted
  !i   vecs  :the array vectors, dimensioned (m,n)
  !i   lopt  :ones digit
  !i           0 vecs returned sorted.
  !i             NB: in this case, vecs must be dimensioned 2*m*n.
  !i           1 vecs is unchanged; only iwk is returned.
  !i          tens digit
  !i           0 vecs sorted by first column, second column, etc
  !i           1 vecs sorted by increasing length
  !o Outputs
  !o   iprm  :a permutation table that sorts array 'vecs'
  !o   vecs  :may be sorted, depending on lopt
  ! ----------------------------------------------------------------
  !     implicit none
  integer :: m,n
  integer :: iprm(n),lopt,mmax
  parameter (mmax=10)
  double precision :: vecs(m,0:n-1),sw(mmax),di,dl
  integer :: i,j,mm,inc,iv

  ! ... Get the largest increment
  inc = 1
10 continue
  inc = 3*inc+1
  if (inc < n) goto 10
  do  i = 1, n
     iprm(i) = i-1
  enddo
  if (lopt /= 0 .AND. lopt /= 1 .AND. lopt /= 11) &
       call rx('dvshel: bad lopt')
  if (n <= 1) return
  if (lopt == 0 .AND. m > mmax) &
       call rx('dvshel increase mmax')

  ! ... Loop over partial sorts
12 continue
  inc = inc/3
  !   ... Outer loop of straight insertion
  do  11  i = inc+1, n
     iv = iprm(i)
     if (lopt == 0) call dcopy(m,vecs(1,iv),1,sw,1)
     j = i
     !     ... Inner loop of straight insertion
     if (lopt == 11)  then
20      continue
        di = 0d0
        dl = 0d0
        do   mm = 1, m
           dl = dl + vecs(mm,iprm(j-inc))**2
           di = di + vecs(mm,iv)**2
        enddo
        if (dl > di) then
           !             print *, 'slip',j-inc,j,i
           iprm(j) = iprm(j-inc)
           j = j-inc
           if (j <= inc) goto 21
           goto 20
        endif
21      continue
     elseif (lopt == 1) then
120     continue
        do  124  mm = 1, m
           !           cases dl.gt.di, dl.eq.di, dl.lt.di
           if (vecs(mm,iprm(j-inc)) - vecs(mm,iv)) 121,124,138
124     enddo
136     continue
        !       ... v(iv) .eq. v(iprm(j-inc)
        goto 121
        !       ... v(iv) .gt. v(iprm(j-inc)
138     continue
        !           print *, 'slip',j-inc,j,i
        iprm(j) = iprm(j-inc)
        j = j-inc
        if (j <= inc) goto 121
        goto 120
121     continue
     elseif (lopt == 0) then
220     continue
        do  224  mm = 1, m
           !           cases dl.gt.di, dl.eq.di, dl.lt.di
           if (vecs(mm,iprm(j-inc)) - sw(mm)) 221,224,238
224     enddo
236     continue
        !       ... v(iv) .eq. v(iprm(j-inc)
        goto 221
        !       ... v(iv) .gt. v(iprm(j-inc)
238     continue
        !           print *, 'slip',j-inc,j,i
        call dcopy(m,vecs(1,j-inc-1),1,vecs(1,j-1),1)
        j = j-inc
        if (j <= inc) goto 221
        goto 220
221     continue
     endif
     !     ... end of straight insertion
     if (lopt == 0) then
        call dcopy(m,sw,1,vecs(1,j-1),1)
     else
        iprm(j) = iv
     endif
11 enddo
  if (inc > 1) goto 12
end subroutine dvshel
subroutine dvperm(m,n,vecs,wk,iprm,lopt)
  !- Permute an array of double precision vectors according to iprm
  ! ----------------------------------------------------------------
  !i Inputs
  ! o  vecs(m,n): n vectors of length m are to be permuted
  !i   iprm: a table of permutation indices to array vecs
  !i   wk:  a work array of length m*m, holds the sorted vectors
  !i  lopt: T, copy wk back to vecs.
  ! ----------------------------------------------------------------
  !     implicit none
  integer :: m,n,iprm(n)
  logical :: lopt
  double precision :: vecs(m,n),wk(m,n)
  integer :: i,j,k

  do  10  i = 1, n
     k = iprm(i)+1
     do    j = 1, m
        wk(j,i) = vecs(j,k)
     enddo
10 enddo
  if (lopt) call dpcopy(wk,vecs,1,n*m,1d0)
end subroutine dvperm

