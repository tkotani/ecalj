      integer function maxnn (nindxv,nindxc,
     d                        nl,nclass)
      implicit real*8 (a-h,o-z)
      integer:: nindxv(nl*nclass),nindxc(nl*nclass)
      integer:: i,ntot,nclass,nl
      maxnn      = -1
      do       i = 1,nl*nclass
        ntot       = nindxv(i) + nindxc(i)
        if (ntot .gt. maxnn) maxnn = ntot
      end do

      return
      end

