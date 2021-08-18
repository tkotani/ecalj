      subroutine ivshel(m,n,iarray,iwk,lopt)
C- shell sort of a array of integer vectors
C ----------------------------------------------------------------
Ci Inputs
Ci   iarray(m,n)
Ci   iwk: a work array of dimension n
Ci   lopt:if T, iwk returned as list of indices to iarray to sort it,
Ci              while iarray is unchanged.
Ci        if F, iarray returned sorted
Co Outputs
Co   iwk a table of indices to array iarray (lopt true)
C ----------------------------------------------------------------
      integer m,n
      logical lopt
      integer iarray(m,0:n-1),iwk(n)
      integer lognb2,i,j,k,l,n2,nn,it,mm,mmm
      lognb2 = int(log(float(n+1))*1.4426950)
      n2 = n
      do i = 1, n
         iwk(i) = i-1
      enddo   
      do  12  nn = 1, lognb2
        n2 = n2/2
        k = n - n2
        do  11  j = 1, k
          i = j
    3     continue
          l = i + n2
c      print *, 'test ',i,l,iwk(i),iwk(l)
          do  15  mm = 1, m
            if (iarray(mm,iwk(l)) - iarray(mm,iwk(i))) 16,15,11
   16       continue
            if (lopt) then
              it = iwk(i)
              iwk(i) = iwk(l)
              iwk(l) = it
c      print 800, (iwk(nnn), nnn=1,n)
c800   format(' swap', 11i5)
            else
c      print 800,  (iarray(nnn,i-1), nnn=1,m)
c      print 800,  (iarray(nnn,l-1), nnn=1,m)
              do  14  mmm = 1, m
                it = iarray(mmm,i-1)
                iarray(mmm,i-1) = iarray(mmm,l-1)
                iarray(mmm,l-1) = it
   14         continue
            endif
            i = i - n2
            if (i .ge. 1) goto 3
            goto 11
   15     continue
   11   continue
   12 continue

      return
      end
      subroutine ivperm(m,n,vecs,wk,iwk,lopt)
C- Permute an array of integer vectors according to iwk
C ----------------------------------------------------------------
Ci Inputs
Cio  vecs(m,n): n vectors of length m are to be permuted
Ci   iwk: a table of permutation indices to array vecs
Ci   wk:  a work array of length m
Ci  lopt: T, copy wk back to vecs.
C ----------------------------------------------------------------
C     implicit none
      integer m,n,iwk(n),vecs(m,n),wk(m,n)
      logical lopt
      integer i,j,k

      do  10  i = 1, n
        k = iwk(i)+1
        do    j = 1, m
           wk(j,i) = vecs(j,k)
        enddo   
   10 continue
      if (lopt) then
        k = n*m
        do    i = 1, k
           vecs(i,1) = wk(i,1)
        enddo   
      endif
      end
C      program main
CC      implicit none
C      integer m,n
C      parameter (m=3,n=16)
C      integer iwk(n),lopt,mm,i,vecs(m,n),wk(m,n)
C      double precision di
C      real ran1
C
C      call ran1in(10)
C      do  10  i = 1,m*n
C   10 vecs(i,1) = 10000*ran1()
C      vecs(1,5) = vecs(1,3)
C      vecs(1,7) = vecs(1,4)
C      vecs(1,8) = vecs(1,4)
C      vecs(2,7) = vecs(2,4)
C      call ivshel(m,n,vecs,iwk,.true.)
C
C      print *, '--- vector and iwk sorted by lopt=.true. ---'
C      do  20  i = 1, n
C        di = 0d0
C        do  44  mm = 1, m
C   44   di = di + vecs(mm,i)**2
C        di = dsqrt(di)
C        print 333, i, iwk(i)+1, di, (vecs(mm,i), mm=1, m)
C  333   format(2i4, f12.0, 6i6)
C   20 continue
C
C      call ivperm(m,n,vecs,wk,iwk,.true.)
C      print *, '--- vector array after calling dvperm ---'
C      do  50  i = 1, n
C        di = 0d0
C        do  54  mm = 1, m
C   54   di = di + vecs(mm,i)**2
C        di = dsqrt(di)
C        print 333, i, iwk(i)+1, di, (vecs(mm,i), mm=1, m)
C   50 continue
C
C      end 

