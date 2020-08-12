      integer function maxnn (nindxv,nindxc,
     d                        nl,nclass)

c 92.03.14
c finds the maximum of n

      implicit real*8 (a-h,o-z)
      dimension nindxv(nl*nclass),nindxc(nl*nclass)

      maxnn      = -1
      do       i = 1,nl*nclass
        ntot       = nindxv(i) + nindxc(i)
        if (ntot .gt. maxnn) maxnn = ntot
      end do

      return
      end

