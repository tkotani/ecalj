c-Taken from Ferdi's GW  -------------------------------------------------------------------
      subroutine genqbz (qbas,n1,n2,n3,
     o qbz,wbz, nstbz, nadd,half)
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
      real(8):: qbas(3,3)
      real(8):: qbz(3,*),wbz(*) !wbz(n1*n2*n3)
      real(8):: qmic(3,3),w1(3),w2(3),w3(3),half(3)
      integer:: nstbz(*),nadd,n1,n2,n3,i1,i2,i3,kount,nnn
!! icase=1
!! vectors forming microcells
      qmic(:,1)= qbas(:,1)/dble(n1)
      qmic(:,2)= qbas(:,2)/dble(n2)
      qmic(:,3)= qbas(:,3)/dble(n3)
      nnn=(n1+nadd)*(n2+nadd)*(n3+nadd)
      nstbz(1:nnn)=0
!! sampling weight
      weight     = 1.d0/dble(n1*n2*n3)
      kount      = 0
      do      i1 = 1,n1+nadd
      do      i2 = 1,n2+nadd
      do      i3 = 1,n3+nadd
        kount      = kount + 1
        qbz(:,kount) = (i1-1+half(1))*qmic(:,1) + (i2-1+half(2))*qmic(:,2) + (i3-1+half(3))*qmic(:,3)
        wbz(kount) = weight
c        if(icase==2.and.(i1==1.or.i1==n1).and.(i2==1.or.i2==n2).and.(i3==1.or.i3==n3)) then
c          nstbz(kount) = 2*2*2
c        endif
c        write(6,"(' iq qbz =',i5,3f9.4)")kount,qbz(:,kount)
      enddo
      enddo
      enddo
      if(kount .ne. (n1+nadd)*(n2+nadd)*(n3+nadd))call rx( 'genqbz: wrong no. k-points')
      end
c--------------------
c      subroutine cv (c,v,n,
c     o w )
c forms w(i) = c * v(i)
c
c      implicit real*8(a-h,o-z)
c      dimension v(n)
c      dimension w(n)
c
c      do       i = 1,n
c      w(i)       = c*v(i)
c      end do
c
c      return
c      end
c----------------

