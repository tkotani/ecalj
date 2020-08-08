      subroutine z2herm(uplo,lds,ns,s)
C- Render a matrix hermitian or replace it by its dagger
C ----------------------------------------------------------------------
Ci Inputs
Ci  uplo   :'U' upper triangle of s is copied to lower triangle
Ci         :'L' lower triangle of s is copied to upper triangle
Ci         :'A' both upper and lower are computed from the average
Ci         :'C' Replace s by s*
Ci         :'D' Replace s by s+
Ci   lds   :leading dimension of s
Ci   ns    :dimension of s
Cio Inputs/Outputs
Cio   s    :s(1:ns,1:ns) is rendered hermitian
Cl Local variables
Cl         :
Cr Remarks
Cr   
Cu Updates
Cu   30 May 05 Added options, 'A','C','D'
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character*1 uplo
      integer lds,ns
      double precision s(2,lds,ns)
C ... Local parameters
      double precision sr,si
      integer i,j

      if (uplo .eq. 'U' .or. uplo .eq. 'u') then
        do  i = 1, ns
          do  j = 1, i
            s(1,i,j) =  s(1,j,i)
            s(2,i,j) = -s(2,j,i)
          enddo
          s(2,i,i) = 0d0
        enddo
      elseif (uplo .eq. 'L' .or. uplo .eq. 'l') then
        do  i = 1, ns
          do  j = 1, i
            s(1,j,i) =  s(1,i,j)
            s(2,j,i) = -s(2,i,j)
          enddo
          s(2,i,i) = 0d0
        enddo
      elseif (uplo .eq. 'A' .or. uplo .eq. 'a') then
        do  i = 1, ns
          do  j = 1, i
            sr = (s(1,j,i) + s(1,i,j))/2
            s(1,j,i) = sr
            s(1,i,j) = sr
            si = (s(2,j,i) - s(2,i,j))/2
            s(2,j,i) =  si
            s(2,i,j) = -si
          enddo
          s(2,i,i) = 0d0
        enddo
      elseif (uplo .eq. 'C' .or. uplo .eq. 'c') then
        do  i = 1, ns
          do  j = 1, ns
            s(2,i,j) = -s(2,i,j)
          enddo
        enddo
      elseif (uplo .eq. 'D' .or. uplo .eq. 'd') then
        do  i = 1, ns
          do  j = 1, i
            sr = s(1,j,i)
            s(1,j,i) = s(1,i,j)
            s(1,i,j) = sr
            si = s(2,j,i)
            s(2,j,i) = -s(2,i,j)
            s(2,i,j) = -si
          enddo
          s(2,i,i) = 0d0
        enddo
      else
        call rxs('z2herm: bad argument uplo = ',uplo)
      endif
      end

