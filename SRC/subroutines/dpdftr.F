      subroutine dpdftr(n1,n2,n3,k1,k2,k3,n,f,lbin,ifi)
C- Dump/read an array of reals given on a Fourier transform mesh.
C ----------------------------------------------------------------------
Ci Inputs
Ci   n1..3 :size of f
Ci   k1..3 :dimensions f
Ci   n     :number of functions f
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Ci   lbin  :T file I/O in binary mode
Ci         :F file I/O in ascii mode
Ci Inputs/Outputs
Cio  f     :array to read/write
Cr Remarks
Cr   f is complex but is stored as real with zero imaginary part.
Cu Updates
Cu   27 Apr 01 Added lbin switch
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      logical lbin
      integer n1,n2,n3,k1,k2,k3,n,ifi
      double complex f(k1,k2,k3,n)
C ... Local parameters
      integer n1mx,jfi,i,i1,i2,i3
C      parameter (n1mx=1024)
C      double precision row(n1mx)
      double precision row(n1)
      n1mx = n1
C --- Input ---
      if (ifi .gt. 0) then
        jfi = ifi
        if (n1 .gt. n1mx) call rx('dpdftr: increase n1mx')
        call dpzero(f, 2*k1*k2*k3*n)
        do  i = 1, n
        do  i3 = 1, n3
        do  i2 = 1, n2
                read(jfi) (row(i1), i1=1,n1)
              do  i1 = 1, n1
                f(i1,i2,i3,i) = dcmplx(row(i1),0d0)
        enddo
        enddo
        enddo
        enddo

C --- Output ---
      elseif (ifi .lt. 0) then
        jfi = -ifi
        do  i = 1, n
          do  i3 = 1, n3
            do  i2 = 1, n2
              do  i1 = 1, n1
                row(i1) = dble(f(i1,i2,i3,i))
              enddo
c              if (lbin) then
                write(jfi) (row(i1), i1=1,n1)
c              else
c                call dfdump(row,n1,ifi)
c              endif
            enddo
          enddo
        enddo

C --- Invalid ifi ---
      else
        call rx('dpdftr: invalid ifi')
      endif

      end

