      subroutine fradhd(nkaps,eh,rsmh,lh,lmxh,nr,rofi,fh,xh,vh,dh)
C- Radial envelope head functions and gradients
C ----------------------------------------------------------------------
Ci Inputs
Ci   nkaps :1s digit: number of envelope function types per l q.n.
Ci   eh    :energies of smoothed Hankel (l- and kappa- dependent)
Ci   rsmh  :hankel smoothing radii (l- and kappa- dependent)
Ci         :rsmh=0 -> channel has no envelope function
Ci   lh    :basis l-cutoff for each of nkape groups
Ci   lmxh  :dimensions xh,dh,fh,vh
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Co Outputs
Co   fh    :(envelope head functions) * r
Co   xh    :(radial derivative of head functions) * r = d(fh)/dr - fh/r
Co   vh    :function values at rofi(nr) = fh/rofi(nr)
Co   dh    :function slopes at rofi(nr)
Cr Remarks
Cr   Radial functions are multiplied by r, but vh and dh are not.
Cr   Derivative calculate analytically, as
Cr     xh = h' = l/r h_l - h_l+1
Cu Updates
Cu   12 Aug 04 First implementation of extended local orbitals;
Cu             substitute call to hansmr with call to hansr.
Cu   10 Apr 02 Redimensioned eh,rsmh to accomodate larger lmax
Cu   16 May 00 Adapted from nfp fradhd.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nkaps,lh(nkaps),lmxh,nr,n0
      parameter (n0=10)
      double precision eh(n0,nkaps),rsmh(n0,nkaps),rofi(nr),
     .fh(nr,0:lmxh,nkaps),vh(0:lmxh,nkaps),
     .xh(nr,0:lmxh,nkaps),dh(0:lmxh,nkaps)
C ... Local parameters
      integer ik,l,i,idx
      double precision rsm,e,rmt,r,vhh,dhh,xi(0:20),rl !,wk(2)
C     double precision hs(0:20),dhs(0:20),ddhs(0:20),xx

      call dpzero(fh,nr*(1+lmxh)*nkaps)
      call dpzero(xh,nr*(1+lmxh)*nkaps)
      call dpzero(vh,(1+lmxh)*nkaps)
      call dpzero(dh,(1+lmxh)*nkaps)

      do  ik = 1, nkaps
        do  l = 0, lh(ik)
          rsm = rsmh(l+1,ik)
          e = eh(l+1,ik)
          if (rsm .gt. 0) then
            rmt = rofi(nr)
            rl = rmt**l
            if (l .eq. 0) rl = 1

C           asm = 1d0/rsm
C           call hansmr(rmt,e,asm,xi,l+1)
C           vh(l,ik) = xi(l) * rmt**l
C           dh(l,ik) = vh(l,ik)*l/rmt - xi(l+1)*rmt**(l+1)

c            call hansr(rsm,0,l+1,1,l+1,e,rmt**2,1,1,idx,wk,10,xi) !we calculate xi at nr.
            call hansr(rsm,0,l+1,1,l+1,e,rmt**2,1,1,idx,10,xi) !we calculate xi at nr.
            vh(l,ik) = xi(l) * rl
            dh(l,ik) = vh(l,ik)*l/rmt - xi(l+1)*rl*rmt

            fh(1,l,ik) = 0
            xh(1,l,ik) = 0
            do  i = 2, nr

              r = rofi(i)
              rl = r**l
              if (l .eq. 0) rl = 1

C             call hansmr(r,e,asm,xi,l+1)
c              call hansr(rsm,0,l+1,1,l+1,e,r**2,1,1,idx,wk,10,xi) !we calculate xi at r**2.
              call hansr(rsm,0,l+1,1,l+1,e,r**2,1,1,idx,10,xi) !we calculate xi at r**2.

C             g0 = dexp(-asm*asm*r*r)
              vhh = xi(l) * rl
              dhh = vhh*l/r - xi(l+1)*rl*r
              fh(i,l,ik) = vhh*r
              xh(i,l,ik) = dhh*r

C             call hansmd(1,r,e,rsm,l,hs,dhs,ddhs,xx,xx,xx)
C             fh(i,l,ik) = hs(l)*r
C             fh(i,l,ik) = dhs(l)*r

            enddo
          endif
        enddo

C       call prrmsh('f',rofi,fh,nr,nr,lh(ik)+1)
C       call prrmsh('df/dr',rofi,xh,nr,nr,lh(ik)+1)

      enddo

      end

