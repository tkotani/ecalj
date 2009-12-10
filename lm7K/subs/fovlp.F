      subroutine fovlp(ib1,ib2,ntab,iax,plat,pos,ipc,alat,rmax,z,pwr,
     .  facaa,facae,facee,fmax,f,inc)
C- Analytical function of the sphere overlaps
Cu Updates
Cu   22 Oct 02  weight ES-ES and atom-ES overlaps differently
Cu              New argument list
C     implicit none
      integer ib1,ib2,ntab(ib2+1),niax,ipc(ib2),inc
      parameter (niax=10)
      integer iax(niax,1)
      double precision alat,plat(3,3),pos(3,ib2),rmax(1),pwr,fmax,f,
     .  facaa,facae,facee
      double precision dsqr,dr,d,sumrs,z(ib2),zb,zi,fac
      integer i,ib,i0,i1,ii,jj,kk,ix,icb,ici

      fmax = -99d0
      f = 0
      inc = 0
      do  10  ib = ib1, ib2
        i0 = ntab(ib)+1
        i1 = ntab(ib+1)
        do  12  i = i0+1, i1
          ii = iax(3,i)-iax(3,i0)
          jj = iax(4,i)-iax(4,i0)
          kk = iax(5,i)-iax(5,i0)
          dsqr = 0
          do  16  ix = 1, 3
            dr = pos(ix,iax(2,i)) - pos(ix,iax(1,i0)) +
     .        plat(ix,1)*ii + plat(ix,2)*jj + plat(ix,3)*kk
            dsqr = dsqr + dr**2
   16     continue
          icb = ipc(ib)
          ici = ipc(iax(2,i))
          zb = z(icb)
          zi = z(ici)
          fac = 1
          if (zb .ne. 0 .and. zi .ne. 0) then
            fac = facaa
          elseif (zb .eq. 0 .and. zi .eq. 0) then
            fac = facee
          else
            fac = facae
          endif
          sumrs = rmax(icb) + rmax(ici)
          d = alat*dsqrt(dsqr)
          if (sumrs .gt. d) then
            f = f + fac*(sumrs/d-1)**pwr
            inc = inc+1
C          else
C            goto 10
          endif
          fmax = max(sumrs/d-1,fmax)
   12   continue
   10 continue
      end

