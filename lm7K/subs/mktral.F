      subroutine mktral(alpha,adot,avw,itrans,kap2,nkap,ldot,lmx,
     .  nbas,ipc,nl,hcr,tral,trad)
C- Makes screening transformation matrix tral
C ----------------------------------------------------------------------
Ci Inputs:
Ci   alpha :tight-binding screening constants
Ci   adot  :(kappa*avw)^2-derivative of tight-binding screening const.
Ci   avw   :average Wigner-Seitz sphere radius
Ci   itrans:characterizes structure matrix transformation (see Remarks)
Ci          Let a = hard core radius, a=hcr*avw
Ci          0:2nd generation transformation:
Ci                 |N^a(kappa)> =|N^0(kappa)>
Ci                 |J^a(kappa)> =|J^0(kappa)>-alpha|K^0(kappa)>
Ci                 |J^a(kappa)> = 0 at a
Ci                 Wronskian W{N^a,J^a} = W{N^0,J^0} = avw/2
Ci          1:head |N^a(kappa)> has same value as |N^0(0)>
Ci                 and has no |J^0(kappa)> component
Ci                 |J^a(kappa)> = 0 at a
Ci                 Wronskian W{N^a,J^a} = W{N^0,J^0} = avw/2
Ci          2:head |N^a(kappa)> has same value & slope as |N^0(0)>
Ci                 |J^a(kappa)> = 0 at a
Ci                 Wronskian W{N^a,J^a} = W{N^0,J^0} = avw/2
Ci          3:head |N^a(kappa)> has value 1 and slope 0
Ci                 |J^a(kappa)> has value 0 and slope avw/(2*a*a) at a
Ci                 Wronskian W{N^a,J^a} = W{N^0,J^0} = avw/2
Ci          4:head |N^a(kappa)> has value 1 and slope 0
Ci                 |J^a(kappa)> has value 0 and slope 1/a**2 at a
Ci          5:head |N^a(kappa)> has value 1 and slope 0
Ci                 |J^a(kappa)> has value 0 and slope -1/a at a
Ci          6:head |N^a(kappa)> has value |N^0(0)> and slope 0
Ci                 |J^a(kappa)> = 0 at a
Ci                 Wronskian W{N^a,J^a} = W{N^0,J^0} = avw/2
Ci   kap2  :(kappa*avw)**2
Ci   ldot  :T: calculate energy derivatives
Ci   lmx   :lmx(j) = maximum l for atom j
Ci   nbas  :size of basis
Ci   nl    :number of l's
Ci   hcr   :hard core screening radius/avw
Co Outputs:
Co   trad  :(kappa*avw)^2-derivative of tral
Co   tral  :transformation matrix for head and tail functions
Cr Remarks: The transformation matrix from 'bare' Neumann and Bessel
Cr           functions |N>,|J> to 'screened' functions |N^a>,|J^a> is
Cr
Cr          |N^a> = N^a(r) = tral(1)*|N> + tral(2)*|J>
Cr          |J^a> = J^a(r) = tral(3)*|N> + tral(4)*|J>
Cr
Cr          |N^a> has a one-center expansion ||N^a> = |N^a> - |J^a>S^a
Cr
Cr  As currently written,
Cr  N and J follow Andersen's definitions (see besslr.f).
Cr
Cr  mktral was adapted from the Stuttgart LMTO, written by R. Tank
Cr  Stuttgart's sigma(l,ic)*wsr(ic) = hcr(l,ic)*avw here.
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer itrans,lmx(*),nl,ipc(1),nbas,nkap
      double precision avw,kap2(nkap),hcr(nl,*),
     .  alpha(nl*nl,nbas,nkap),adot(nl*nl,nbas,nkap),
     .  trad(4,nl**2,nbas,nkap),tral(4,nl**2,nbas,nkap)
      logical ldot
C Local variables:
      integer lmxx,ib,ik,lm,m,ic,iprint,k,l,lgunit,ipr,iclbsj
      parameter(lmxx=20)
      double precision ad,al,dj,djdot,dn,dndot,dn0,fi(-1:lmxx),
     .                 gi(-1:lmxx),j,jdot,n,ndot,n0,r2,rfac,hcrad,
     .                 w1,w2,w3,w4,w5,t(8),fac
      character*(54) strn(0:6), sout
C External calls:
      external  bessl2,dcopy,iprint,lgunit
      data strn /
     .  'N^a = N^0,  J^a = J^0 - alpha K^0',
     .  'N^a = N^0(e=0), J^a = 0 at hcr',
     .  'N^a has val,slo = N^0(e=0), J^a = 0 at hcr',
     .  'N^a has val,slo 1,0, J^a has val,slo 0, w/2a^2 at hcr',
     .  'N^a has val,slo 1,0, J^a has val,slo 0, 1/a^2 at hcr',
     .  'N^a has val,slo 1,0, J^a has val,slo 0, -1/a at hcr',
     .  'N^a has val=N^0(0), slo=0, J^a = 0 at hcr'/


      call getpr(ipr)
      if (ipr .gt. 50/1) write(lgunit(1),300) strn(itrans)
  300 format(' MKTRAL: rotate so ',a/'  ib ic  l   ',
     .  '  kap2        a            b            c            d')
      if (itrans .lt .0 .or. itrans .gt. 6)
     .  call rxi('MKTRAL: itrans out of range, value',itrans)
      do  10  ik = 1, nkap
      do  10  ib = 1, nbas
        ic = ipc(ib)
        lm = 0
        do  12  l = 0, nl-1
          if (l .le. lmx(ic)) then
          do  15  m = -l, l
          lm = lm+1
          hcrad = hcr(l+1,ic)
          r2 = hcrad*hcrad
          gi(l+1) = 0
          call bessl2(kap2(ik)*r2,-1,l+1,fi(-1),gi(-1))
          rfac=hcrad**(-l-1)
          n0  = rfac
          n   = rfac*gi(l)
          ndot= rfac*gi(l-1)/(4*l-2)*r2
          dn0 = -l-1
          dn  = l-(l+l+1)*gi(l+1)/gi(l)
          dndot = l-(l+l-1)*gi(l)/gi(l-1)
C         j0 = 1d0/(4*l+2)/(rfac*hcrad)
          j = fi(l)/(rfac*hcrad)
          jdot = -fi(l+1)/(4*l+2)/(rfac/hcrad)
          dj   =-(l+1)+(l+l-1)*fi(l-1)/fi(l  )
          djdot=-(l+1)+(l+l+1)*fi(l)/fi(l+1)
          al = alpha(lm,ib,ik)
          ad = adot(lm,ib,ik)
          if (itrans .eq. 0)  then
            t(1) = 1d0
            t(2) = 0d0
            t(3) = -al
            t(4) = 1d0
            if (ldot) then
              t(5) = 0d0
              t(6) = 0d0
              t(7) = -ad
              t(8) = 0d0
            endif
          elseif (itrans .eq. 1) then
            t(1) = 1d0/gi(l)
            t(2) = 0d0
            t(3) = -al/t(1)
            t(4) = 1d0/t(1)
            if (ldot) then
              t(5) = -gi(l-1)/(4*l-2)*r2/gi(l)/gi(l)
              t(6) = 0d0
              t(7) = -ad/t(1) + t(5)*al/t(1)/t(1)
              t(8) = -t(5)/t(1)/t(1)
            endif
          elseif (itrans .eq. 2) then
            w1 = n   *j   *hcrad*avw*(dj - dn )
            w2 = n0  *j   *hcrad*avw*(dj - dn0)
            w3 = n   *n0  *hcrad*avw*(dn0 - dn )
            w4 = n0  *jdot*hcrad*avw*(djdot - dn0)
            w5 = ndot*n0  *hcrad*avw*(dn0 - dndot )
            t(1) = w2/w1
            t(2) = w3/w1
            t(4) = 1d0/(t(1) + t(2)*al)
            t(3) = -al*t(4)
            if (ldot) then
              t(5) = w4/w1
              t(6) = w5/w1
              t(7) = -(t(1)*ad/al/al-t(5)/al - t(6))*t(3)*t(3)
              t(8) = -(t(5)+al*t(6) + ad*t(2))*t(4)*t(4)
            endif
          elseif (itrans .eq. 3) then
            w1 = n*j *hcrad*avw*(dj - dn)
            w2 = j *hcrad*avw*dj
            w3 = n   *hcrad*avw*(-dn)
            w4 = jdot*hcrad*avw*djdot
            w5 = ndot*hcrad*avw*(-dndot)
            t(1) = w2/w1
            t(2) = w3/w1
            t(4) = 1d0/(t(1)+t(2)*al)
            t(3) = -al*t(4)
            if (ldot) then
              t(5) = w4/w1
              t(6) = w5/w1
              t(7) = -(t(1)*ad/al/al-t(5)/al - t(6))*t(3)*t(3)
              t(8) = -(t(5)+al*t(6) + ad*t(2))*t(4)*t(4)
            endif
          elseif (itrans .eq. 4) then
            w1 = n*j *hcrad*avw*(dj-dn)
            w2 = j *hcrad*avw*dj
            w3 = n   *hcrad*avw*(-dn)
            w4 = jdot*hcrad*avw*djdot
            w5 = ndot*hcrad*avw*(-dndot)
            t(1) = w2/w1
            t(2) = w3/w1
            t(4) = 1d0/(t(1) + t(2)*al)
            t(3) = -al*t(4)
            if (ldot) then
              t(5) = w4/w1
              t(6) = w5/w1
              t(7) = -(t(1)*ad/al/al-t(5)/al - t(6))*t(3)*t(3)
              t(8) = -(t(5)+al*t(6) + ad*t(2))*t(4)*t(4)
            endif
            fac = 2d0/avw
            t(3) = t(3)*fac
            t(4) = t(4)*fac
            if (ldot) then
              t(7) = t(7)*fac
              t(8) = t(8)*fac
            endif
c           write(*,*)t(1), 2*hcrad*dj*j
c           write(*,*)t(2),-2*hcrad*dn*n
c           write(*,*)t(3),-2/avw*j
c           write(*,*)t(4), 2/avw*n
          elseif (itrans .eq. 5) then
            w1 = n*j *hcrad*avw*(dj-dn)
            w2 = j *hcrad*avw*dj
            w3 = n   *hcrad*avw*(-dn)
            w4 = jdot*hcrad*avw*djdot
            w5 = ndot*hcrad*avw*(-dndot)
            t(1) = w2/w1
            t(2) = w3/w1
            t(4) = 1d0/(t(1) + t(2)*al)
            t(3) = -al*t(4)
            if (ldot) then
              t(5) = w4/w1
              t(6) = w5/w1
              t(7) = -(t(1)*ad/al/al-t(5)/al - t(6))*t(3)*t(3)
              t(8) = -(t(5)+al*t(6) + ad*t(2))*t(4)*t(4)
            endif
            fac = -2d0*hcrad
            t(3) = t(3)*fac
            t(4) = t(4)*fac
            if (ldot) then
              t(7) = t(7)*fac
              t(8) = t(8)*fac
            endif
          elseif (itrans .eq. 6) then
            w1 = n   *j   *hcrad*avw*(dj-dn)
            w2 = n0  *j   *hcrad*avw*dj
            w3 = n   *n0  *hcrad*avw*(-dn)
            w4 = n0  *jdot*hcrad*avw*djdot
            w5 = ndot*n0  *hcrad*avw*(-dndot)
            t(1) = w2/w1
            t(2) = w3/w1
            t(4) = 1d0/(t(1) + t(2)*al)
            t(3) = -al*t(4)
c           write(*,*)t(1), 2*hcrad*n0*dj*j
c           write(*,*)t(2),-2*hcrad*n0*dn*n
c           write(*,*)t(3),-j/n0
c           write(*,*)t(4), n/n0
            if (ldot) then
              t(5) = w4/w1
              t(6) = w5/w1
              t(7) = -(t(1)*ad/al/al-t(5)/al - t(6))*t(3)*t(3)
              t(8) = -(t(5)+al*t(6) + ad*t(2))*t(4)*t(4)
            endif
          endif

          if (ipr .ge. 50/1 .and. iclbsj(ic,ipc,nbas,1) .eq. ib
     .        .and. l .eq. m) then
            if (l .eq. 0 .and. ib .eq. 1) then
              print 302, ib,ic,l,kap2(ik)/avw**2,(t(k),k=1,4)
            elseif (l .eq. 0) then
              print 303, ib,ic,l,(t(k),k=1,4)
            else
              print 304, l,(t(k),k=1,4)
            endif
  302       format(i4,2i3,f11.6,4f13.8)
  303       format(i4,2i3,11x,4f13.8)
  304       format(4x,3x,i3,11x,4f13.8)
          endif
          call dcopy(4,t(1),1,tral(1,lm,ib,ik),1)
          if (ldot) call dcopy(4,t(5),1,trad(1,lm,ib,ik),1)
   15     continue
        else
          lm = lm + 2*l+1
        endif
   12   continue
   10 continue
      return

C ... Returns a string describing what transformation itral does.
C     Work around gnu f95 bug, v 4.1.1
#ifndef LINUXI
      entry mktrli(itrans,sout)
#else
      end

      subroutine mktrli(itrans,sout)
      character sout*(*)
      character*(54) strn(0:6)
      data strn /
     .  'N^a = N^0,  J^a = J^0 - alpha K^0',
     .  'N^a = N^0(e=0), J^a = 0 at hcr',
     .  'N^a has val,slo = N^0(e=0), J^a = 0 at hcr',
     .  'N^a has val,slo 1,0, J^a has val,slo 0, w/2a^2 at hcr',
     .  'N^a has val,slo 1,0, J^a has val,slo 0, 1/a^2 at hcr',
     .  'N^a has val,slo 1,0, J^a has val,slo 0, -1/a at hcr',
     .  'N^a has val=N^0(0), slo=0, J^a = 0 at hcr'/
#endif

      sout = strn(itrans)
      end

