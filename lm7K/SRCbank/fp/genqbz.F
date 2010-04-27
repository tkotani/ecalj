      subroutine genqbz (qbas,n1,n2,n3,
     o qbz,wbz )
c 92.02.16
c generates the k-points in the 1BZ
c the 1BZ is a parallepiped formed by G1,G2,G3 (qbas(3,3)) 
c this is divided into microcells defined by G1/n1,G2/n2,G3/n3
c the k-points may be thought of as being centred at each microcell
c the sampling weight for each k-point is the same (1/n1*n2*n3)

c qbas = base reciprocal vectors G1,G2,G3
c n1,n2,n3 = divisions along G1,G2,G3

c qbz  = k-points in the 1BZ
c wbz  = sampling weight for qbz

      implicit real*8 (a-h,o-z)
      dimension qbas(3,3)
      dimension qbz(3,n1*n2*n3),wbz(n1*n2*n3)
      dimension qmic(3,3),w1(3),w2(3),w3(3)

c vectors forming microcells
      call cv      (1.d0/dble(n1),qbas(1,1),3,qmic(1,1))
      call cv      (1.d0/dble(n2),qbas(1,2),3,qmic(1,2))
      call cv      (1.d0/dble(n3),qbas(1,3),3,qmic(1,3))

c sampling weight
      weight     = 1.d0/dble(n1*n2*n3)
   
      kount      = 0
      do      i1 = 1,n1
      call cv      (dble(i1-1),qmic(1,1),3,w1)
      do      i2 = 1,n2
      call cv      (dble(i2-1),qmic(1,2),3,w2)
      do      i3 = 1,n3
      call cv      (dble(i3-1),qmic(1,3),3,w3)
      
      kount      = kount + 1
      qbz(1,kount) = w1(1) + w2(1) + w3(1)  
      qbz(2,kount) = w1(2) + w2(2) + w3(2)  
      qbz(3,kount) = w1(3) + w2(3) + w3(3)  
      wbz(kount)   = weight

      end do
      end do
      end do
      if (kount .ne. n1*n2*n3)stop 'genqbz: wrong no. k-points'

c write to file KPNT
      ifkpt = 335  
      open(ifkpt,file='KPTin1BZ') 
      do      i1 = 1,kount
      write (ifkpt,6000)i1,qbz(1,i1),qbz(2,i1),qbz(3,i1),wbz(i1)
      end do
      close (ifkpt)

c formats
 6000 format (1x,i4,4f10.5)
      return
      end
      subroutine cv (c,v,n,
     o w )
      
c forms w(i) = c * v(i)

      implicit real*8(a-h,o-z)
      dimension v(n)
      dimension w(n)

      do       i = 1,n
      w(i)       = c*v(i)
      end do

      return
      end

