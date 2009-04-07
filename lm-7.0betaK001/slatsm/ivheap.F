      subroutine ivheap(m,n,vecs,iprm,opts)
C- Heapsort array of integer vectors
C ----------------------------------------------------------------
Ci Inputs
Ci   vecs(m,n): n vectors of length m are to be sorted
Ci   iprm: an integer work array of dimension n, or if vecs returned
Ci        in sorted order, an array of length m*n
Ci   opts: ones digit
Ci           0 vecs returned sorted.
Ci           1 vecs is unchanged; only iprm is returned
Ci         tens digit
Ci           0 vecs sorted
Ci           1 vecs sorted by increasing length
Ci         hundreds digit
Ci           1 equal vectors preserve their original order
Ci
Co Outputs
Co   iprm a permutation table that sorts array 'vecs'
Co   vecs may be changed, depending on opts
C ----------------------------------------------------------------
C     implicit none
      integer m,n,iprm(n),opts
      integer vecs(m,n)
      integer di,dl
      integer l,ir,irra,i,j,mm,i1,i2,nn
      logical norm

      do  2  ir = 1, n
    2 iprm(ir) = ir
      if (n .le. 1) return
      norm = mod(opts/10,10) .ne. 0
      l = n/2+1
      ir = n

C --- For each l = n/2+1, 1, -1 do ---
   10 continue
C ... "Hiring phase"
      if (l .gt. 1) then
        l = l-1
        irra = iprm(l)
C ... "Retirement-and-promotion phase"
      else
        irra = iprm(ir)
        iprm(ir) = iprm(1)
*       call awrit3('ir%i: %n:1i',' ',180,6,ir,n,iprm)
        ir = ir-1
        if (ir .eq. 1) then
          iprm(1) = irra
*         call awrit2('exit %n:1i',' ',180,6,n,iprm)
          goto 100
        endif
      endif

C ... Setup to sift down element irra to proper level
      i = l
      j = l+l

C --- Do while j .le. ir ---
   20 if (j .le. ir) then

C   ... Increment j if vecs(iprm(j+1)) > vecs(iprm(j))
        if (j .lt. ir) then
          if (norm) then
            di = 0
            dl = 0
            do  24  mm = 1, m
            dl = dl + vecs(mm,iprm(j))**2
   24       di = di + vecs(mm,iprm(j+1))**2
            if (di-dl .gt. 0) j = j+1
          else
            do  26  mm = 1, m
            if (abs(vecs(mm,iprm(j+1))-vecs(mm,iprm(j))) .le. 0) goto 26
            if (vecs(mm,iprm(j+1))-vecs(mm,iprm(j)) .gt. 0) j = j+1
            goto 28
   26       continue
   28       continue
          endif
        endif

C   ... Demote rra to its level
        if (norm) then
          di = 0
          dl = 0
          do  34  mm = 1, m
          dl = dl + vecs(mm,irra)**2
   34     di = di + vecs(mm,iprm(j))**2
          if (di-dl .gt. 0) then
            iprm(i) = iprm(j)
*           call awrit4('%i,%i: %n:1i',' ',180,6,i,j,n,iprm)
            i = j
            j = j+j
C     ... This is rra's level; set j to terminate the sift-down
          else
            j = ir+1
          endif
        else
          do  36  mm = 1, m
C     ...   Skip over equal elements
            if (abs(vecs(mm,iprm(j))-vecs(mm,irra)) .le. 0) goto 36
            if (vecs(mm,iprm(j))-vecs(mm,irra) .gt. 0) then
              iprm(i) = iprm(j)
*             call awrit4('%i,%i: %n:1i',' ',180,6,i,j,n,iprm)
              i = j
              j = j+j
C     ... This is rra's level; set j to terminate the sift-down
            else
              j = ir+1
            endif
            goto 38
   36     continue
C     ... Case rra = vec(iprm(j))
          j = ir+1
   38     continue
        endif
        go to 20
      endif
C ... Put rra into its slot
      iprm(i) = irra
*     call awrit3('%i: %n:1i',' ',180,6,i,n,iprm)
      go to 10

C --- For equal vectors, restore original ordering ---
  100 continue
      if (mod(opts/100,10) .eq. 0) goto 200
      i2 = 0
C ... Find i1,i2 = next range of equal numbers
  110 i1 = i2+1
      if (i1 .gt. n) goto 200
  120 i2 = i2+1
      if (i2 .gt. n) goto 130
      if (norm) then
        di = 0
        dl = 0
        do  124  mm = 1, m
        dl = dl + vecs(mm,iprm(i1))**2
  124   di = di + vecs(mm,iprm(i2))**2
        if (di-dl .gt. 0) goto 130
      else
        do  126  mm = 1, m
  126   if (abs(vecs(mm,iprm(i2))-vecs(mm,iprm(i1))) .gt. 0) goto 130
      endif
C ... vec(i1) = vec(i2) ; imcrement i2 and try again
      goto 120

C --- Sort iprm(i1)..iprm(i2) ---
  130 continue
      i2 = i2-1
      i1 = i1-1
      nn = i2-i1
      if (nn .le. 1) goto 110
      l = nn/2+1
      ir = nn

C ... For each l = (i2-i1+1)/2+1, 1, -1 do
  140 continue
C ... "Hiring phase"
      if (l .gt. 1) then
        l = l-1
        irra = iprm(l+i1)
C ... "Retirement-and-promotion phase"
      else
        irra = iprm(ir+i1)
        iprm(ir+i1) = iprm(1+i1)
        ir = ir-1
        if (ir .eq. 1) then
          iprm(1+i1) = irra
          goto 110
        endif
      endif

C ... Setup to sift down element irra to proper level ...
      i = l
      j = l+l

C ... Do while j .le. ir ...
  150 if (j .le. ir) then

C   ... Increment j if iprm(j+i11) > iprm(j+i1))
        if (j .lt. ir) then
          if (iprm(j+i1) .lt. iprm(j+1+i1)) j = j+1
        endif
C   ... Demote irra to its level
        if (irra .lt. iprm(j+i1)) then
          iprm(i+i1) = iprm(j+i1)
          i = j
          j = j+j
C   ... This is irra's level; set j to terminate the sift-down
        else
          j = ir+1
        endif
        go to 150
      endif
C ... Put rra into its slot
      iprm(i+i1) = irra
      go to 140

C --- Sort vecs ---
  200 continue
      if (mod(opts,10) .eq. 0)
     .  call ivprm(m,n,vecs,vecs(1,n+1),iprm,.true.)

      end
      subroutine ivprm(m,n,vecs,wk,iprm,lopt)
C- Permute an array of integer vectors according to iprm
C ----------------------------------------------------------------
Ci Inputs
Cio  vecs(m,n): n vectors of length m are to be permuted
Ci   iprm: a table of permutation indices to array vecs
Ci   wk:  a work array of length m
Ci   lopt: T, copy wk back to vecs.
C ----------------------------------------------------------------
C     implicit none
      integer m,n,iprm(n)
      logical lopt
      integer vecs(m,n),wk(m,n)
      integer i,j,k

      do  10  i = 1, n
        k = iprm(i)
        do  12  j = 1, m
   12   wk(j,i) = vecs(j,k)
   10 continue
      if (lopt) call icopy(n*m,wk,1,vecs,1)
      end

