      subroutine fmain
      implicit none
      integer m,n
      parameter (m=3,n=16)
      integer iwk(n),opts,mm,i
      integer vecs(m,n*2),wk(m,n)
      double precision di
      real ran1


      call ran1in(10)
      do  10  i = 1,m*n
   10 vecs(i,1) = ran1()*10000
      vecs(1,5) = vecs(1,3)
      vecs(1,7) = vecs(1,4)
      vecs(1,8) = vecs(1,4)
      vecs(2,7) = vecs(2,4)
      call icopy(m*n,vecs,1,wk,1)

      print *, '--- vector and iwk sorted by opts=10 ---'
      call ivheap(m,n,vecs,iwk,10)
      do  12  i = 1, n
        di = 0d0
        do  42  mm = 1, m
   42   di = di + vecs(mm,i)**2
        di = dsqrt(di)
        print 333, i, iwk(i), di, (vecs(mm,i), mm=1, m)
  333   format(2i4,f12.2,1x,5i5)
   12 continue

      print *, '--- vector and iwk sorted by opts=11 ---'
      call icopy(m*n,wk,1,vecs,1)
      call ivheap(m,n,vecs,iwk,11)
      do  20  i = 1, n
        di = 0d0
        do  44  mm = 1, m
   44   di = di + vecs(mm,i)**2
        di = dsqrt(di)
        print 333, i, iwk(i), di, (vecs(mm,i), mm=1, m)
   20 continue

      print *, '--- length in order of increasing iwk ---'
      do  30  i = 1, n
        di = 0d0
        do  45  mm = 1, m
   45   di = di + vecs(mm,iwk(i))**2
        di = dsqrt(di)
        print 333, i, iwk(i), di
   30 continue

      print *, '--- vector array after calling ivprm ---'
      call ivprm(m,n,vecs,wk,iwk,.true.)
      do  50  i = 1, n
        di = 0d0
        do  54  mm = 1, m
   54   di = di + vecs(mm,i)**2
        di = dsqrt(di)
        print 333, i, iwk(i), di, (vecs(mm,i), mm=1, m)
   50 continue


C --- Check 100s digit ---
      vecs(1,1) = 1
      vecs(2,1) = 1
      vecs(3,1) = 2
      vecs(1,2) = 1
      vecs(2,2) = 2
      vecs(3,2) = 1
      vecs(1,3) = 1
      vecs(2,3) = 2
      vecs(3,3) = 1
      vecs(1,4) = 2
      vecs(2,4) = 2
      vecs(3,4) = 1
      vecs(1,5) = 1
      vecs(2,5) = 1
      vecs(3,5) = 2
      vecs(1,6) = 1
      vecs(2,6) = 1
      vecs(3,6) = 2
      vecs(1,7) = 1
      vecs(2,7) = 1
      vecs(3,7) = 2
      vecs(1,8) = 1
      vecs(2,8) = 1
      vecs(3,8) = 2
      vecs(1,9) = 1
      vecs(2,9) = 1
      vecs(3,9) = 2
      vecs(1,10) = 3
      vecs(2,10) = 1
      vecs(3,10) = 2
      vecs(1,11) = 1
      vecs(2,11) = 2
      vecs(3,11) = 2
      vecs(1,12) = 1
      vecs(2,12) = 1
      vecs(3,12) = 1
      vecs(1,13) = 1
      vecs(2,13) = 1
      vecs(3,13) = 2
      vecs(1,14) = 1
      vecs(2,14) = 1
      vecs(3,14) = 1
      vecs(1,15) = 1
      vecs(2,15) = 2
      vecs(3,15) = 2
      vecs(1,16) = 1
      vecs(2,16) = 1
      vecs(3,16) = 2

      print *, '--- vector and iwk sorted by opts=1,101 ---'
      call icopy(m*n,vecs,1,wk,1)
      call ivheap(3,n,vecs,iwk,1)
      call awrit2('iwk for opts=1   %n:1i',' ',80,6,n,iwk)
      call ivheap(3,n,vecs,iwk,101)
      call awrit2('iwk for opts=101 %n:1i',' ',80,6,n,iwk)
      call ivprm(m,n,vecs,wk,iwk,.true.)
      do  60  i = 1, n
        di = 0d0
        do  62  mm = 1, m
   62   di = di + vecs(mm,i)**2
        di = dsqrt(di)
        print 333, i, iwk(i), di, (vecs(mm,i), mm=1, m)
   60 continue

      print *, '--- vector and iwk sorted by opts=11,111 ---'
      call icopy(m*n,wk,1,vecs,1)
      call ivheap(3,n,vecs,iwk,11)
      call awrit2('iwk for opts=11  %n:1i',' ',80,6,n,iwk)
      call ivheap(3,n,vecs,iwk,111)
      call awrit2('iwk for opts=111 %n:1i',' ',80,6,n,iwk)
      call ivprm(m,n,vecs,wk,iwk,.true.)
      do  70  i = 1, n
        di = 0d0
        do  72  mm = 1, m
   72   di = di + vecs(mm,i)**2
        di = dsqrt(di)
        print 333, i, iwk(i), di, (vecs(mm,i), mm=1, m)
   70 continue

      call ran1in(10)
      do  110  i = 1,m*n
  110 vecs(i,1) = ran1()*10000
      vecs(1,5) = vecs(1,3)+ran1()/1e4
      vecs(1,7) = vecs(1,4)-ran1()/1e4
      vecs(1,8) = vecs(1,4)+ran1()/1e4
      vecs(2,7) = vecs(2,4)-ran1()/1e4
      call icopy(m*n,vecs,1,wk,1)

      print *, '--- initial vector ---'
      do  111  i = 1, n
        di = 0d0
        do  141  mm = 1, m
  141   di = di + vecs(mm,i)**2
        di = dsqrt(di)
        print 333, i, 0, di, (vecs(mm,i), mm=1, m)
  111 continue

      print *, '--- vector and iwk sorted by opts=0 ---'
      call ivheap(m,n,vecs,iwk,0)
      do  112  i = 1, n
        di = 0d0
        do  142  mm = 1, m
  142   di = di + vecs(mm,i)**2
        di = dsqrt(di)
        print 333, i, iwk(i), di, (vecs(mm,i), mm=1, m)
  112 continue

      print *, '--- vector and iwk sorted by opts=1  ---'
      call icopy(m*n,wk,1,vecs,1)
      call ivheap(m,n,vecs,iwk,1)
      do  120  i = 1, n
        di = 0d0
        do  144  mm = 1, m
  144   di = di + vecs(mm,i)**2
        di = dsqrt(di)
        print 333, i, iwk(i), di, (vecs(mm,i), mm=1, m)
  120 continue

      call ivprm(m,n,vecs,wk,iwk,.true.)

      print *, '--- vector array after calling ivprm ---'
      do  150  i = 1, n
        di = 0d0
        do  154  mm = 1, m
  154   di = di + vecs(mm,i)**2
        di = dsqrt(di)
        print 333, i, iwk(i), di, (vecs(mm,i), mm=1, m)
  150 continue

      call ivheap(m,n,vecs,iwk,1)
      print *, '--- vector and iwk sorted by opts=1 ---'
      do  220  i = 1, n
        di = 0d0
        do  244  mm = 1, m
  244   di = di + vecs(mm,i)**2
        di = dsqrt(di)
        print 333, i, iwk(i), di, (vecs(mm,i), mm=1, m)
  220 continue

      end 
