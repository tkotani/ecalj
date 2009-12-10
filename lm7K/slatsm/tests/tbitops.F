      implicit none
      double precision x,dround,dtrunc,xnint

      x = xnint(dble(1.23456789d0))
      print *, x

      x = xnint(dble(1.56789d0))
      print *, x


      x = xnint(dble(-1.23456789d0))
      print *, x

      x = xnint(dble(-1.56789d0))
      print *, x

      x = dround(dble(1.23456789d0),6)
      print *, x

      x = dtrunc(dble(1.23456789d0),6)
      print *, x



      end
      
C         det(1) = 1.eq. dble(0d0)
C         det(1) = 1.eq..1d0
C         det(1) = 1.eq.dble(1.d0)
C         det(1) = 1.eq.dble(1d0)
C         det(1) = dble(1.0d0)
C         det(1) =dble(.10d0)
C         det(2) = dble(0.0d0)
C
C         det(1) = 1.eq. 0d0
C         det(1) = 1.eq..1d0
C         det(1) = 1.eq.1.d0
C         det(1) = 1.eq.1d0
C         det(1) = 1.0d0
C         det(1) =.10d0
C         det(2) = 0.0d0
C
C         ten = 10.0d0
C         do  50  i = 1, n
C            if (ipvt(i) .ne. i) det(1) = -det(1)
C            det(1) = a(i,i)*det(1)
Cc        ...exit
C            if (det(1) .eq. 0.0d0) goto 60
C   10       if (dabs(det(1)) .ge. 1.0d0) goto 20
C               det(1) = ten*det(1)
C               det(2) = det(2) - 1.0d0
C            goto 10
C   20       continue
C   30       if (dabs(det(1)) .lt. ten) goto 40
C               det(1) = det(1)/ten
C               det(2) = det(2) + 1.0d0
C            goto 30
C   40       continue
C   50    continue
C   60    continue
C   70 continue
Cc
Cc     compute inverse(u)
Cc
C      if (mod(job,10) .eq. 0) goto 150
C         do  100  k = 1, n
C            a(k,k) = 1.0d0/a(k,k)
C            t = -a(k,k)
C            call dscal(k-1,t,a(1,k),1)
C            kp1 = k + 1
C            if (n .lt. kp1) goto 90
C            do  80  j = kp1, n
C               t = a(k,j)
C               a(k,j) = 0.0d0
C               call daxpy(k,t,a(1,k),1,a(1,j),1)
C   80       continue
C   90       continue
C  100    continue
Cc
Cc        form inverse(u)*inverse(l)
Cc
C         nm1 = n - 1
C         if (nm1 .lt. 1) goto 140
C         do  130  kb = 1, nm1
C            k = n - kb
C            kp1 = k + 1
C            do  110  i = kp1, n
C               work(i) = a(i,k)
C               a(i,k) = 0.0d0
C  110       continue
C            do  120  j = kp1, n
C               t = work(j)
C               call daxpy(n,t,a(1,j),1,a(1,k),1)
C  120       continue
C            l = ipvt(k)
C            if (l .ne. k) call dswap(n,a(1,k),1,a(1,l),1)
C  130    continue
C  140    continue
C  150 continue
C      return
CC#endif
C      end
