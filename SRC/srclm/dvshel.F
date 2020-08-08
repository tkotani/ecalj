      subroutine dvshel(m,n,vecs,iprm,lopt)
C- Shell sort of a array of double precision vectors
C ----------------------------------------------------------------
Ci Inputs
Ci   m     :length of each vector
Ci   n     :number of vectors to be sorted
Ci   vecs  :the array vectors, dimensioned (m,n)
Ci   lopt  :ones digit
Ci           0 vecs returned sorted.
Ci             NB: in this case, vecs must be dimensioned 2*m*n.
Ci           1 vecs is unchanged; only iwk is returned.
Ci          tens digit
Ci           0 vecs sorted by first column, second column, etc
Ci           1 vecs sorted by increasing length
Co Outputs
Co   iprm  :a permutation table that sorts array 'vecs'
Co   vecs  :may be sorted, depending on lopt
C ----------------------------------------------------------------
C     implicit none
      integer m,n
      integer iprm(n),lopt,mmax
      parameter (mmax=10)
      double precision vecs(m,0:n-1),sw(mmax),di,dl
      integer i,j,mm,inc,iv

C ... Get the largest increment
      inc = 1
   10 continue
      inc = 3*inc+1
      if (inc .lt. n) goto 10
      do  2  i = 1, n
    2 iprm(i) = i-1
      if (lopt .ne. 0 .and. lopt .ne. 1 .and. lopt .ne. 11)
     .call rx('dvshel: bad lopt')
      if (n .le. 1) return
      if (lopt .eq. 0 .and. m .gt. mmax)
     .call rx('dvshel increase mmax')

C ... Loop over partial sorts
   12 continue
      inc = inc/3
C   ... Outer loop of straight insertion
      do  11  i = inc+1, n
        iv = iprm(i)
        if (lopt .eq. 0) call dcopy(m,vecs(1,iv),1,sw,1)
        j = i
C     ... Inner loop of straight insertion
        if (lopt .eq. 11)  then
   20     continue
          di = 0d0
          dl = 0d0
          do  24  mm = 1, m
            dl = dl + vecs(mm,iprm(j-inc))**2
   24     di = di + vecs(mm,iv)**2
          if (dl .gt. di) then
c             print *, 'slip',j-inc,j,i
            iprm(j) = iprm(j-inc)
            j = j-inc
            if (j .le. inc) goto 21
            goto 20
          endif
   21     continue
        elseif (lopt .eq. 1) then
  120     continue
          do  124  mm = 1, m
C           cases dl.gt.di, dl.eq.di, dl.lt.di
            if (vecs(mm,iprm(j-inc)) - vecs(mm,iv)) 121,124,138
  124     continue
  136     continue
C       ... v(iv) .eq. v(iprm(j-inc)
          goto 121
C       ... v(iv) .gt. v(iprm(j-inc)
  138     continue
C           print *, 'slip',j-inc,j,i
          iprm(j) = iprm(j-inc)
          j = j-inc
          if (j .le. inc) goto 121
          goto 120
  121     continue
        elseif (lopt .eq. 0) then
  220     continue
          do  224  mm = 1, m
C           cases dl.gt.di, dl.eq.di, dl.lt.di
            if (vecs(mm,iprm(j-inc)) - sw(mm)) 221,224,238
  224     continue
  236     continue
C       ... v(iv) .eq. v(iprm(j-inc)
          goto 221
C       ... v(iv) .gt. v(iprm(j-inc)
  238     continue
C           print *, 'slip',j-inc,j,i
          call dcopy(m,vecs(1,j-inc-1),1,vecs(1,j-1),1)
          j = j-inc
          if (j .le. inc) goto 221
          goto 220
  221     continue
        endif
C     ... end of straight insertion
        if (lopt .eq. 0) then
          call dcopy(m,sw,1,vecs(1,j-1),1)
        else
          iprm(j) = iv
        endif
   11 continue
      if (inc .gt. 1) goto 12
      end
      subroutine dvperm(m,n,vecs,wk,iprm,lopt)
C- Permute an array of double precision vectors according to iprm
C ----------------------------------------------------------------
Ci Inputs
Cio  vecs(m,n): n vectors of length m are to be permuted
Ci   iprm: a table of permutation indices to array vecs
Ci   wk:  a work array of length m*m, holds the sorted vectors
Ci  lopt: T, copy wk back to vecs.
C ----------------------------------------------------------------
C     implicit none
      integer m,n,iprm(n)
      logical lopt
      double precision vecs(m,n),wk(m,n)
      integer i,j,k

      do  10  i = 1, n
        k = iprm(i)+1
        do  12  j = 1, m
   12   wk(j,i) = vecs(j,k)
   10 continue
      if (lopt) call dpcopy(wk,vecs,1,n*m,1d0)
      end

