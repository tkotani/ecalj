      subroutine lx2vec(lmx,nl0,nl,lx)

C- Unpacks a compressed lmx: digit 9 unpacked as -1; 8 as 0
C  If nl0 nonzero, lmx assumed to have nl0 digits;
C  thus lx2vec(956,4,nl,lx) returns nl=4, lx=(0 -1 5 6).
C  If nl0 is zero, lmx assumed to have log10(lmx) digits;
C  thus lx2vec(956,0,nl,lx) returns nl=3, lx=(-1 5 6).
C     implicit none
      integer lmx,nl,nl0,lx(1)
      integer i,j,k,getdig,nlx

      nlx = nl0
      if (nlx .le. 0) then
        nlx = 1 + int(dlog10(dble(lmx)+1d-6))
        if (lmx .eq. 0) nlx = 1
      endif

      k = 0
      do  10  i = nlx-1, 0, -1
        k = k+1
        j = getdig(lmx,i,10)
        if (j .eq. 9) then
          j = -1
        elseif (j .eq. 8) then
          j = 0
        else
          nl = k
        endif
        lx(k) = j
   10 continue
      end subroutine lx2vec


