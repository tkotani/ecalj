      subroutine s2oka(job,lds,ld2,nlma,nlmb,avw,s,sc)
C- Scales complex strux from standard to Andersen conventions
C ----------------------------------------------------------------------
Ci     job:1s digit:
Ci          0: s is real
Ci          1: s is complex
Ci         10s digit distinguishes how complex arithmetic is handled
Ci          0: s has real, imaginary separated
Ci             s = s(ldg,ldg2,2), with s(*,*,1..2) = real..imag
Ci          1: s is returned complex*16 format:
Ci             s = s(2,ldg,ldg2), with s(1..2,*,*) = real..imag
Ci          2: s has real, imaginary separated by columns
Ci             s = s(ldg,2,ldg2), with s(*,1..2,*) = real..imag
Ci         100s digit 
Ci             nonzero => return complex conjugate of s
Ci   lds   :leading dimension of s, or sc
Ci   ld2   :second dimension of s. Its value depends on the complex
Ci          storage format.  For specified 10s digit job, use:
Ci          0  ld2 = formal second dimension
Ci          1  ld2 = not used
Ci          2  ld2 = 2
Ci Inputs/Outputs
Cio   s    :s is scaled to Andersen's conventions
Cio   sc   :same as s, but in complex form of the second type
Cio        :only one of s or sc is used.
Cr Remarks
Cr  Andersen structure constants differ from standard conventions by:
Cr    S_RL,R'L' (OKA) / S_RL,R'L' =  2/((2l-1)!!(2l'-1)!!) avw**(l+l'+1)
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer job,lds,ld2,nlma,nlmb
      double precision s(lds,ld2,2),avw
      double complex sc(lds,ld2)
C ... Local parameters
      logical ljob2
      integer kcplx,nlmax,lmaxx,lmxa,lmxb,ll,lmax,la,lb,ilma,ilmb
      parameter (nlmax=20,lmaxx=2*nlmax-2)
      double precision fac2l(0:lmaxx),fac
      kcplx  = mod(job/10,10)

      ljob2 = mod(job/100,10) .ne. 0
      lmxa = ll(nlma)
      lmxb = ll(nlmb)
      lmax = max(lmxa,lmxb)
      if (lmax .gt. lmaxx) call rxi('s2oka: increase lmaxx, need',lmax)

C     A table of (2l-1)!!
      fac2l(0) = 1
      do  5  la = 1, lmax
    5 fac2l(la) = fac2l(la-1) * (2*la-1)

      if (mod(job,10) .eq. 0) then
        call rx('s2oka not implemented')
      else
        if (kcplx .ne. 1) call rx('s2oka not implemented')
        do  la = 0, lmxa
        do  lb = 0, lmxb
          fac = 2/fac2l(la)/fac2l(lb)*avw**(la+lb+1)
          if (ljob2) then
            do  ilma = la**2+1, (la+1)**2
            do  ilmb = lb**2+1, (lb+1)**2
              sc(ilma,ilmb) = dconjg(sc(ilma,ilmb))*fac
            enddo
            enddo
          else
            do  ilma = la**2+1, (la+1)**2
            do  ilmb = lb**2+1, (lb+1)**2
              sc(ilma,ilmb) = sc(ilma,ilmb)*fac
            enddo
            enddo
          endif
        enddo
        enddo
      endif

      end

