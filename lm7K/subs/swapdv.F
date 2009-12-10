      subroutine swapdV(nl2,nsites,nsp,iax,dh)
C- Cludge to swap Vsp with Vps and Vpd with Vdp and Vsd with Vds
C ----------------------------------------------------------------
C Remarks
C  Force dh(r1,l1,r2,l2) = dh(r2,l2,r1,l1) 
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nl2,nsites,nsp,niax
      parameter (niax=10)
      integer iax(niax,nsites)
      double precision dh(nl2,nl2,nsites*nsp**2,3)
C Local parameters
      integer i,j,k,lm1,lm2,nsp2,isp,ii,jj
      double precision temp
C     .  ,d1mach

      if (nl2 .lt. 2) return

      nsp2 = nsp**2
      do  100  isp = 1, nsp2
      do  100  i = 1, nsites
        ii = i+(isp-1)*nsites
        do  110  j = i, nsites
          jj = j+(isp-1)*nsites
          if (iax(1,j) .eq. iax(2,i) .and. iax(2,j) .eq. iax(1,i) .and.
     .      iax(3,i) .eq. -iax(3,j) .and.
     .      iax(4,i) .eq. -iax(4,j) .and.
     .      iax(5,i) .eq. -iax(5,j)) then

C --- swap sp,ps ---
            do  120  lm1 = 2, 4
              do  1  k = 1, 3
C                print *, ii,jj,k,lm1
C                print 333, dh(lm1,1,ii,k), dh(lm1,1,jj,k)
C  333           format(2g17.6)
C                if (dh(lm1,1,ii,k)*dh(lm1,1,jj,k) .lt. -d1mach(3)) then
C                  print*,' i = ',ii,'  j = ',jj,'  k = ',k,
C     .              '  lm1 = ',lm1
C                  call rx('bug in swapdV')
C                endif
                temp = dh(lm1,1,ii,k)
                dh(lm1,1,ii,k) = dh(lm1,1,jj,k)
                dh(lm1,1,jj,k) = temp
    1         continue
  120       continue
            if (nl2 .lt. 5) goto 100

C --- swap sd,ds ---
            do  140  lm1 = 5, 9
              do  3  k = 1, 3
C                print *, ii,jj,k,lm1
C                print 333, dh(lm1,1,ii,k), dh(lm1,1,jj,k)
C                if (dh(lm1,1,ii,k)*dh(lm1,1,jj,k) .gt. d1mach(3)) then
C                  print*,' i = ',ii,'  j = ',jj,'  lm1 = ',lm1
C                  call rx('bug in swapdV')
C                endif
                temp = dh(lm1,1,ii,k)
                dh(lm1,1,ii,k) = -dh(lm1,1,jj,k)
                dh(lm1,1,jj,k) = -temp
    3         continue
  140       continue

C --- swap pd,dp ---
            do  130  lm1 = 5, 9
              do  130  lm2 = 2, 4
                do  2  k = 1, 3
C                  print *, ii,jj,k,lm1,lm2
C                  print 333, dh(lm1,lm2,ii,k), dh(lm1,lm2,jj,k)
C                  if (dh(lm1,lm2,ii,k)*dh(lm1,lm2,jj,k) .lt. -d1mach(3))
C     .              then
C                    print*,' i = ',ii,'  j = ',jj,'  lm1 = ',lm1,
C     .                '  lm2 = ',lm2
C                    call rx('bug in swapdV')
C                  endif
                  temp = dh(lm1,lm2,ii,k)
                  dh(lm1,lm2,ii,k) = dh(lm1,lm2,jj,k)
                  dh(lm1,lm2,jj,k) = temp
    2           continue
  130       continue
            goto 100
          endif
  110   continue
  100 continue
      end

