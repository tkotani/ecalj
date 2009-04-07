#define SP 1
      subroutine remhyb(ldim,nl,nbas,h,o)
C- remove hybridisation blocks of H and O
C-----------------------------------------------------------------------
Ci ldim,nl,nbas,h,o
Co h and o changed on exit
Cr precompiler option SP retains sp hybridisation
C-----------------------------------------------------------------------
C     implicit none
      integer ldim,nl,nbas
      double precision h(0:ldim-1,0:ldim-1,2),o(0:ldim-1,0:ldim-1,2)

      integer nl2,ip,jp,nlm1,nlm2

      if (nl .eq. 1) return

      nl2 = nl**2
      do  5  nlm1 = 0, nl2-1
        if (nlm1 .eq. 0) then
          do  1  nlm2 = 0, nl2-1
#if SP
            if (nlm2 .gt. 3) then
#else
            if (nlm2 .gt. 0) then
#endif
              do  11  ip = 0, (nbas-1)*nl2, nl2
              do  11  jp = 0, (nbas-1)*nl2, nl2
                h(nlm1+ip,nlm2+jp,1) = 0d0
                h(nlm1+ip,nlm2+jp,2) = 0d0
                o(nlm1+ip,nlm2+jp,1) = 0d0
                o(nlm1+ip,nlm2+jp,2) = 0d0
   11         continue 
            endif
    1     continue 
          goto 5
        elseif (nlm1 .lt. 4) then
          do  2  nlm2 = 0, nl2-1
#if SP
            if (nlm2 .gt. 3) then
#else
            if (nlm2 .eq. 0 .or. nlm2 .gt. 3) then
#endif
              do  12  ip = 0, (nbas-1)*nl2, nl2
              do  12  jp = 0, (nbas-1)*nl2, nl2
                h(nlm1+ip,nlm2+jp,1) = 0d0
                h(nlm1+ip,nlm2+jp,2) = 0d0
                o(nlm1+ip,nlm2+jp,1) = 0d0
                o(nlm1+ip,nlm2+jp,2) = 0d0
   12         continue 
            endif
    2     continue 
          goto 5
        elseif (nlm1 .lt. 9) then
          do  3  nlm2 = 0, nl2-1
            if (nlm2 .lt. 4 .or. nlm2 .gt. 8) then
              do  13  ip = 0, (nbas-1)*nl2, nl2
              do  13  jp = 0, (nbas-1)*nl2, nl2
                h(nlm1+ip,nlm2+jp,1) = 0d0
                h(nlm1+ip,nlm2+jp,2) = 0d0
                o(nlm1+ip,nlm2+jp,1) = 0d0
                o(nlm1+ip,nlm2+jp,2) = 0d0
   13         continue 
            endif
    3     continue 
          goto 5
        elseif (nlm1 .lt. 16) then
          do  4  nlm2 = 0, nl2-1
            if (nlm2 .lt. 9) then
              do  14  ip = 0, (nbas-1)*nl2, nl2
              do  14  jp = 0, (nbas-1)*nl2, nl2
                h(nlm1+ip,nlm2+jp,1) = 0d0
                h(nlm1+ip,nlm2+jp,2) = 0d0
                o(nlm1+ip,nlm2+jp,1) = 0d0
                o(nlm1+ip,nlm2+jp,2) = 0d0
   14         continue 
            endif
    4     continue 
        endif
    5 continue 
    6 continue 
      end
      
      

