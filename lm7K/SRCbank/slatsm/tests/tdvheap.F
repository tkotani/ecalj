      subroutine fmain
      implicit none
      integer m,n
      parameter (m=3,n=16)
      integer iwk(n),opts,mm,i
      double precision vecs(m,n*2),di,wk(m,n)
      real ran1

      call ran1in(10)
      do  10  i = 1,m*n
   10 vecs(i,1) = ran1()
      vecs(1,5) = vecs(1,3)
      vecs(1,7) = vecs(1,4)
      vecs(1,8) = vecs(1,4)
      vecs(2,7) = vecs(2,4)
      call dcopy(m*n,vecs,1,wk,1)

      print *, '--- vector and iwk sorted by opts=10, tol=.1 ---'
      call dvheap(m,n,vecs,iwk,.1d0,10)
      do  12  i = 1, n
        di = 0d0
        do  42  mm = 1, m
   42   di = di + vecs(mm,i)**2
        di = dsqrt(di)
        print 333, i, iwk(i), di, (vecs(mm,i), mm=1, m)
   12 continue

      print *, '--- vector and iwk sorted by opts=11 ---'
      call dcopy(m*n,wk,1,vecs,1)
      call dvheap(m,n,vecs,iwk,0d0,11)
      do  20  i = 1, n
        di = 0d0
        do  44  mm = 1, m
   44   di = di + vecs(mm,i)**2
        di = dsqrt(di)
        print 333, i, iwk(i), di, (vecs(mm,i), mm=1, m)
  333   format(2i4, 6f12.6)
   20 continue

      print *, '--- length in order of increasing iwk ---'
      do  30  i = 1, n
        di = 0d0
        do  45  mm = 1, m
   45   di = di + vecs(mm,iwk(i))**2
        di = dsqrt(di)
        print 333, i, iwk(i), di
   30 continue

      print *, '--- vector array after calling dvprm ---'
      call dvprm(m,n,vecs,wk,iwk,.true.)
      do  50  i = 1, n
        di = 0d0
        do  54  mm = 1, m
   54   di = di + vecs(mm,i)**2
        di = dsqrt(di)
        print 333, i, iwk(i), di, (vecs(mm,i), mm=1, m)
   50 continue


      call ran1in(10)
      do  110  i = 1,m*n
  110 vecs(i,1) = ran1()
      vecs(1,5) = vecs(1,3)+ran1()/1e4
      vecs(1,7) = vecs(1,4)-ran1()/1e4
      vecs(1,8) = vecs(1,4)+ran1()/1e4
      vecs(2,7) = vecs(2,4)-ran1()/1e4
      call dcopy(m*n,vecs,1,wk,1)

      print *, '--- vector and iwk sorted by opts=0, tol=1d-4 ---'
      call dvheap(m,n,vecs,iwk,1d-4,0)
      do  112  i = 1, n
        di = 0d0
        do  142  mm = 1, m
  142   di = di + vecs(mm,i)**2
        di = dsqrt(di)
        print 333, i, iwk(i), di, (vecs(mm,i), mm=1, m)
  112 continue

      print *, '--- vector and iwk sorted by opts=1, tol=1d-4 ---'
      call dcopy(m*n,wk,1,vecs,1)
      call dvheap(m,n,vecs,iwk,1d-4,1)
      do  120  i = 1, n
        di = 0d0
        do  144  mm = 1, m
  144   di = di + vecs(mm,i)**2
        di = dsqrt(di)
        print 333, i, iwk(i), di, (vecs(mm,i), mm=1, m)
  120 continue

      call dvprm(m,n,vecs,wk,iwk,.true.)

      print *, '--- vector array after calling dvprm ---'
      do  150  i = 1, n
        di = 0d0
        do  154  mm = 1, m
  154   di = di + vecs(mm,i)**2
        di = dsqrt(di)
        print 333, i, iwk(i), di, (vecs(mm,i), mm=1, m)
  150 continue

      call dvheap(m,n,vecs,iwk,1d-5,1)
      print *, '--- vector and iwk sorted by opts=1, tol=1d-5 ---'
      do  220  i = 1, n
        di = 0d0
        do  244  mm = 1, m
  244   di = di + vecs(mm,i)**2
        di = dsqrt(di)
        print 333, i, iwk(i), di, (vecs(mm,i), mm=1, m)
  220 continue

      end 
