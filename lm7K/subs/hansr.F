c      subroutine hansr(rsm,lmn,lmx,nxi,lxi,exi,rsq,nrx,nr,idx,wk,job,xi)
      subroutine hansr(rsm,lmn,lmx,nxi,lxi,exi,rsq,nrx,nr,idx,job,xi)
C- Vector of smoothed Hankel functions, set of negative e's
C ---------------------------------------------------------------
Ci Inputs
Ci   rsm     smoothing radius of smoothed Hankel
Ci   nrx,lmx dimensions xi
Ci   nxi,exi,lxi:number of energies, energies and lmax to generate fns
Ci   rsq,nr  vector of points r**2, and number of points.
Ci   job      1s digit nonzero, scale xi by r**l
Ci           10s digit nonzero, input rsq is already sorted
Ci   idx     integer work array of length 2*nrx.
Ci           Not needed if 10s digit of job set.
Co Outputs
Co   xi      smoothed Hankels for: xi(1..nr, 0..lxi(ie), 1..nxi)
Co           xi is the radial part/r**l, so the solid Hankel is
Co           hl(ilm) = xi(l)*cy(ilm)*yl(ilm)
Cl Local variables
Cl   n1,n2 : Assuming points are sorted: points are evaluated as follows
Cl         : 1..n1-1  are evaluated by power series expansion
Cl         : n1..n2-1 are evaluated by explicit generation of l=-1,0
Cl                    and upward recursion for higher l.
Cl         : n2..nr   asymtotic form: sm-H has become regular H.
Cl         : If the points are not sorted, they are grouped into
Cl         : three bins, with n1-1, n2-n1, and nr-n2+1 points in them
Cl         : idx is a permutation index that keeps track of the grouping
Cr Remarks
Cr   Points are partitioned into three length scales:
Cr     r<rc1 are evaluated by a polynomial expansion.
Cr     rc1<r<rc2 are evaluated from error functions and the higher l's
Cr       by upward recursion.
Cr     rc2<r are approximated with unsmoothed Hankels.
Cr   The relative error should be less than parameter 'tol', except
Cr   in a narrow region for r~rc1 and l>6, where the precision degrades
Cr   somewhat, worsening with higher l.  For all cases tested, the
Cr   relative error continued to be ~<10^-13 for l<=9.
Cu Updates
Cu   11 May 07 (S. Lozovoi) small bug fixes; similarity with hansmz
C ---------------------------------------------------------------
C     implicit none
      integer nrx,nr,lmn,lmx,idx(nrx,2),nxi,lxi(nxi),job
      double precision rsq(*),e,exi(nxi)
      double precision xi(nrx,lmn:lmx,nxi),wk(nrx,4+lmx-lmn)
C Local variables
      integer ir,j,l,k,n0,n1,n2,ie,lmax
      double precision a,rsm,y0,a2,emin,tol
      double precision rc1,rc2,akap,rl,rl0
      parameter (tol=1d-15)
      logical ltmp,lsort,lscal

      lscal = mod(job,10) .ne. 0
      lsort = mod(job/10,10) .ne. 0

C --- Check lmx; handle case rsm=0 ---
      ltmp = rsm .lt. 1d-12
      lmax = -1
      do  5  ie = 1, nxi
        if (lxi(ie) .gt. lmx) call rx('hansr: lxi gt lmx')
        if (exi(ie) .gt. 0) call rx('hansr: exi gt 0')

        if (ltmp) then
          call hanr(rsq,lmn,lxi(ie),nrx,nr,exi(ie),xi(1,lmn,ie))
        endif
        lmax = max(lmax,lxi(ie))
    5 continue
      if (ltmp) goto 60

C --- Find cutoffs rc2 (negligible smoothing) and rc1 (power series) ---
      emin = 0
      do  10  ie = 1, nxi
        emin = min(emin,exi(ie))
   10 continue
      akap = dsqrt(-emin)
      a = 1/rsm
      a2 = a*a
      y0 = 1/dsqrt(16*datan(1d0))

C ... For r>rc2 approximate smooth Hankels with normal ones
      rc2 = akap/(2*a)
      rc2 = ((rc2 + dsqrt(rc2**2 - dlog(tol)))/a)**2

C ... This rc1 generates a relative precision of ~10^-15 for r~rc1
C     and machine precision for r>>rc1 or r<<rc1.
C     For l>6 and r close to rc1, the precision degrades somewhat.
      rc1 = (rsm*(1.4d0+dble(lmax)/20))**2

C --- Separate the small from the large ---
      n0 = 0
      n1 = 0
      n2 = nr+1
C ... Case points already sorted.  Find n1,n2.
      if (lsort) then
        n1 = 1
        if (nr .eq. 1) then
          if (rsq(1) .lt. rc1) n1 = 2
          if (rsq(1) .gt. rc2) n2 = 1
        else
          if (rsq(1) .ge. rc1) then
            n1 = 1
          else
            call huntx(rsq,nr,rc1,0,n1)
            n1 = n1+1
          endif
          if (rsq(nr) .le. rc2) then
            n2 = nr+1
          else
            n2 = nr
            call huntx(rsq,nr,rc2,0,n2)
            n2 = n2+1
          endif
        endif

C ... Case points not sorted (iwk, wk(3) required now.)
      else
C     On output, lsort is true if points already sorted.
        lsort = .true.
        do  12  ir = 1, nr
C     n1 is offset to block rc1<r<rc2,  n2 offset to block r>rc2
C     idx is a map of original list, separating into the three groups
C     wk(*,3) is a table of r**2 for permuted list of points
          if (rsq(ir) .lt. rc2) then
            if (rsq(ir) .lt. rc1) then
              n0 = n0+1
              wk(n0,3) = rsq(ir)
              idx(ir,1) = n0
            else
              n1 = n1+1
              idx(ir,1) = n1
              idx(n1,2) = ir
            endif
          else
            n2 = n2-1
            wk(n2,3) = rsq(ir)
            idx(ir,1) = n2
          endif
          if (ir .eq. 1) goto 12
          if (rsq(ir) .lt. rsq(ir-1)) lsort = .false.
   12   continue
C ... Now we can poke wk(*,3) for the n1 intermediate points
        if (.not. lsort .or. .true.) then
          do  14  j = 1, n1
            k = idx(j,2)
            idx(k,1) = idx(k,1)+n0
            wk(n0+j,3) = rsq(k)
   14     continue
        endif
        n1 = n0+1
      endif

C ... For debugging
c     call awrit8(' hansr: sort=%l rc1=%,2;2d rc2=%,2;2d'//
c    .  ' nxi=%i nr=%i (%i pwr, %i smooth, %i asym)',
c    .  ' ',120,6,lsort,dsqrt(rc1),dsqrt(rc2),nxi,nr,n1-1,n2-n1,nr-n2+1)


C --- Setup for the energy-independent wk, points n1..n2 ---
      if (lsort) then
        do  20  ir = n1, n2-1
          wk(ir,1) = y0*dexp(-rsq(ir)*a2)
   20   continue
      else
        do  22  ir = n1, n2-1
          wk(ir,1) = y0*dexp(-wk(ir,3)*a2)
   22   continue
      endif

C --- Start loop over energies ---
      do  40  ie = 1, nxi
        e = exi(ie)
        akap = dsqrt(-e)
        lmax = lxi(ie)
C   ... Case calculate points in original order (already sorted)
        if (lsort) then
C     ... Power series for points within rc1
          call hansr1(rsq(1),lmn,lmax,nrx,n1-1,e,rsm,dsqrt(rc1),
     .    xi(1,lmn,ie))
C     ... Normal evaluation of smoothed Hankels
          call hansr2(rsq(n1),lmn,lmax,nrx,n2-n1,e,rsm,wk(n1,1),
     .    wk(n1,2),xi(n1,lmn,ie))
C     ... Asymtotic case, r>>rsm
          if(n2 .le. nr)
     .    call hanr(rsq(n2),lmn,lmax,nrx,nr+1-n2,e,xi(n2,lmn,ie))
C   ... Case calculated points in sorted
        else
C     ... Power series for points within rc1
          call hansr1(wk(1,3),lmn,lmax,nrx,n1-1,e,rsm,dsqrt(rc1),
     .    wk(1,4))
C     ... Normal evaluation of smoothed Hankels
          call hansr2(wk(n1,3),lmn,lmax,nrx,n2-n1,e,rsm,wk(n1,1),
     .    wk(n1,2),wk(n1,4))
C     ... Asymtotic case, r>>rsm
          if(n2 .le. nr)
     .    call hanr(wk(n2,3),lmn,lmax,nrx,nr+1-n2,e,wk(n2,4))
C     ... Poke into xi(lmn:lmax), with the original ordering of points
          do  30  l = lmn, lmax
            do  30  ir = 1, nr
              j = idx(ir,1)
              xi(ir,l,ie) = wk(idx(ir,1),4+l-lmn)
   30     continue
        endif
   40 continue

C --- Scale by r**l if job nonzero ---
   60 continue
      if (.not. lscal) return
      do  62  ir = 1, nr
        rl0 = dsqrt(rsq(ir))
        do  62  ie = 1, nxi
          rl = rl0
          do  62  l = 1, lxi(ie)
            xi(ir,l,ie) = xi(ir,l,ie)*rl
            rl = rl*rl0
   62 continue

      end
      subroutine hansr1(rsq,lmin,lmax,nrx,nr,e,rsm,rmax,xi)
C- Vector of smoothed hankel functions for l=0...lmax, negative e
C  by power series expansion.
C ---------------------------------------------------------------
Ci Inputs
Ci   rsq,nr vector of points r**2, and number of points.
Ci   nrx    dimensions xi; nrx must be gt nr.
Ci   e,rsm  smoothing radius and energy
Ci   lmin   starting l for which to evaluate xi (must be 0 or 1).
Ci   lmax   highest l for which to evaluate xi (must be 0 or 1).
Ci   rmax:  points rsq are less than rmax**2.
Co Outputs:
Co   xi(1..nr,lmin:lmax)
Cr Remarks
Cr   xi is the radial part divided by r**l.
Cr   This routine is intended for evaluation of smoothed hankels
Cr   for small r (r<rsm or so).
Cr   hansr1 tries to evaluate the polynomial in-line for a 14th order
Cr   polynomial, and a 20th order.  Failing that, it evaluates the
Cr   polynomial to whatever order is needed to bring the convergence
Cr   to a relative precision of 'tol'.
C ---------------------------------------------------------------
C     implicit none
      integer nrx,nr,lmin,lmax,nmax,nm1,nm2
      double precision rsq(nrx),e,xi(nrx,lmin:lmax),rsm,rmax
C Local vairables
c      parameter (nmax=40,nm1=14,nm2=20)
      parameter (nmax=60,nm1=14,nm2=20)
      double precision derfc,cof0(-1:20),cofl(0:nmax),tol
      double precision a,a2,add,akap,al,cc,fac,rhs,ta,ta2l,y0,r2max,x
      integer i,l,ir,m,nmaxl
c     parameter (tol=1d-16)
      parameter (tol=1d-20)

C --- Setup ---
      if (lmax .lt. lmin .or. nr .eq. 0) return
      if (lmin .lt. -1 .or. lmin .gt. 0) call rx('hansr1: bad lmin')

      y0 = 1/dsqrt(16*datan(1d0))
      a = 1/rsm
      ta = a+a
      a2 = a*a
      akap = dsqrt(-e)
      cc = 4d0*y0*a*dexp(e/(ta*ta))
      r2max = rmax**(2*nm1)

C --- 0 order coefficient ---
      fac = derfc(akap/ta)
      cof0(-1) = fac/akap
      cof0(0)  = cc - akap*fac
      al = cof0(0)
      rhs = cc*(2*a2)
      fac = 1d0
      do  10  l = 1, lmax
        al = -(e*al + rhs) / (2*l*(2*l+1))
        rhs = -rhs*a2/l
        fac = -2d0*fac*l
        cof0(l) = fac*al
   10 continue

C --- For each l, generate xi(*,l) by power series ---
      ta2l = (2*a2)
      if (lmin .eq. -1) cc = cc/ta2l
      do  20  l = lmin, lmax
        rhs = cc*ta2l
        add = cof0(l)
        cofl(0) = add
C   --- Coffs to polynomial of order nm1 ---
        do  21  i = 1, nm1
          add = -(e*add + rhs) / ( 2*i*(2*i+(l+l+1)) )
          cofl(i) = add
          rhs = -rhs*a2/i
   21   continue
C   ... Use it if it is accurate enough
        if (dabs(add*r2max) .lt. cof0(l)*tol) then
          do  51  ir = 1, nr
            x = rsq(ir)
            xi(ir,l) = (((((((((((((cofl(14)*x+
     .      cofl(13))*x+cofl(12))*x+cofl(11))*x+cofl(10))*x+cofl(9))*x+
     .      cofl(8))*x+cofl(7))*x+cofl(6))*x+cofl(5))*x+cofl(4))*x+
     .      cofl(3))*x+cofl(2))*x+cofl(1))*x+cofl(0)
   51     continue
C ---   Coffs to polynomial of order nm2 ---
C   ... Use it if it is accurate enough
        else
          do  22  i = nm1+1, nm2
            add = -(e*add + rhs) / ( 2*i*(2*i+(l+l+1)) )
            cofl(i) = add
            rhs = -rhs*a2/i
   22     continue
          if (dabs(add*r2max) .lt. cof0(l)*tol) then
            do  52  ir = 1, nr
              x = rsq(ir)
              xi(ir,l) = (((((((((((((((((((cofl(20)*x+cofl(19))*x+
     .        cofl(18))*x+cofl(17))*x+cofl(16))*x+cofl(15))*x+cofl(14))*x+
     .        cofl(13))*x+cofl(12))*x+cofl(11))*x+cofl(10))*x+cofl(9))*x+
     .        cofl(8))*x+cofl(7))*x+cofl(6))*x+cofl(5))*x+cofl(4))*x+
     .        cofl(3))*x+cofl(2))*x+cofl(1))*x+cofl(0)
   52       continue
          else
C ---   Polynomial to nmaxl ---
            do  23  i = nm2+1, nmax
              add = -(e*add + rhs) / ( 2*i*(2*i+(l+l+1)) )
              cofl(i) = add
              nmaxl = i
              if (dabs(add*r2max) .lt. cof0(l)*tol) goto 24
              rhs = -rhs*a2/i
   23       continue
            print 333, tol
  333       format(' hansr1 (warning):  not converged to tol=',1pe8.1)
            call rx('hansr1 failed to converge')
   24       continue
            do  53  ir = 1, nr
              xi(ir,l) = cofl(nmaxl)
   53       continue
            do  54  m = nmaxl, 1, -1
              do  54  ir = 1, nr
                xi(ir,l) = xi(ir,l)*rsq(ir) + cofl(m-1)
   54       continue
          endif
        endif

        ta2l = ta2l*(2*a2)
   20 continue

      end
      subroutine hansr2(rsq,lmin,lmax,nrx,nr,e,rsm,wk,wk2,xi)
C- Vector of smoothed hankel functions for l=0...lmax, negative e.
C ---------------------------------------------------------------
Ci Inputs
Ci   rsq,nr:vector of points r**2, and number of points.
Ci   lmin:  starting l for which to evaluate xi (must be 0 or -1).
Ci   lmax:  highest l for which to evaluate xi.
Ci   e,rsm: smoothing radius and energy
Ci   nrx:   leading dimension of xi
Ci   wk:    array containing y0*dexp(-(r/rsm)**2)
Ci   wk2:   a work array of length nr.
Co Outputs:
Co   xi:    generated for points ir=1..nr and lmin..lmax
Co   wk2:   (2/rsm**2)**(lmax)*4/rsm*dexp(-(akap*rsm/2)**2)*wk(ir)
Co          (can be used to generate xi to higher l)
Cr Remarks
Cr   xi is the radial part divided by r**l.
Cr   xi is evaluated by upward recursion for l>lmin+2.
C ---------------------------------------------------------------
C     implicit none
      integer nrx,nr,lmin,lmax
      double precision rsq(nrx),e,xi(nrx,lmin:lmax),rsm,wk(nr),wk2(nr)
C Local variables
      double precision sre,r2,xx,ra,h0,arsm,earsm
      double precision akap,a,r,um,up,x,facgl,facdu,dudr
      integer l,ir
C ... erfc(x) is evaluated as a ratio of polynomials,
C     to a relative precision of <10^-15 for x<5.
C     Different polynomials are used for x<1.3 and x>1.3.
C     Numerators and denominators are t,b respectively.
      double precision w,f1,f2,
     .t10,t11,t12,t13,t14,t15,t16,t17,b11,b12,b13,b14,b15,b16,b17,b18,
     .t20,t21,t22,t23,t24,t25,t26,t27,b21,b22,b23,b24,b25,b26,b27,b28
      parameter (
     .t10=2.1825654430601881683921d0, t20=0.9053540999623491587309d0,
     .t11=3.2797163457851352620353d0, t21=1.3102485359407940304963d0,
     .t12=2.3678974393517268408614d0, t22=0.8466279145104747208234d0,
     .t13=1.0222913982946317204515d0, t23=0.3152433877065164584097d0,
     .t14=0.2817492708611548747612d0, t24=0.0729025653904144545406d0,
     .t15=0.0492163291970253213966d0, t25=0.0104619982582951874111d0,
     .t16=0.0050315073901668658074d0, t26=0.0008626481680894703936d0,
     .t17=0.0002319885125597910477d0, t27=0.0000315486913658202140d0,
     .b11=2.3353943034936909280688d0, b21=1.8653829878957091311190d0,
     .b12=2.4459635806045533260353d0, b22=1.5514862329833089585936d0,
     .b13=1.5026992116669133262175d0, b23=0.7521828681511442158359d0,
     .b14=0.5932558960613456039575d0, b24=0.2327321308351101798032d0,
     .b15=0.1544018948749476305338d0, b25=0.0471131656874722813102d0,
     .b16=0.0259246506506122312604d0, b26=0.0061015346650271900230d0,
     .b17=0.0025737049320207806669d0, b27=0.0004628727666611496482d0,
     .b18=0.0001159960791581844571d0, b28=0.0000157743458828120915d0)
C ... f1(w=x-1/2) is erfc(x) for 0<x<1.3, if xx is y0*dexp(-x*x)
      f1(w) = xx*(((((((t17*w+t16)*w+t15)*w+t14)*w+t13)*w+t12)*
     .w+t11)*w+t10)/((((((((b18*w+b17)*w+b16)*w+b15)*w+b14)*
     .w+b13)*w+b12)*w+b11)*w+1)
C ... f2(w=x-2) is erfc(x) for x>1.3, if xx is y0*dexp(-x*x)
      f2(w) = xx*(((((((t27*w+t26)*w+t25)*w+t24)*w+t23)*w+t22)*
     .w+t21)*w+t20)/((((((((b28*w+b27)*w+b26)*w+b25)*w+b24)*
     .w+b23)*w+b22)*w+b21)*w+1)

C --- Setup ---
      if (lmax .lt. lmin .or. nr .eq. 0) return
      if (lmin .lt. -1 .or. lmin .gt. 0) call rx('hansr2: bad lmin')
      a = 1/rsm
      akap = dsqrt(-e)
      arsm = akap*rsm/2
      earsm = dexp(-arsm**2)/2
      facgl = (2*a**2)*8*a*earsm
C ... uncomment the following for upward recursion from l=1...
      if (lmin .eq. -1) facgl = facgl/(2*a**2)
      facdu = 8*a*earsm

C --- xi(*,lmin), xi(*,lmin+1) ---
      do  20  ir = 1, nr
        r2 = rsq(ir)
        r = dsqrt(r2)
        ra = r*a
        sre = akap*r
        h0 = dexp(-sre)/r
        xx = earsm*wk(ir)/r

C ---   Evaluate um,up ---
        x = ra - arsm
        if (x .gt. 1.3d0) then
          um = h0 - f2(x-2d0)
        elseif (x .gt. 0) then
          um = h0 - f1(x-.5d0)
        elseif (x .gt. -1.3d0) then
          um = f1(-x-.5d0)
        else
          um = f2(-x-2d0)
        endif
C ...   Evaluation of up assumes x gt 0
        x = ra + arsm
        if (x .gt. 1.3d0) then
          up = f2(x-2d0)
        else
          up = f1(x-.5d0)
        endif

C   ... xi(0) = um - up
        xi(ir,0) = um - up
C   ... xi(-1) = (um + up)*r/akap
        if (lmin .eq. -1) then
          xi(ir,-1) = (um + up)*r/akap
C   ... xi(1)
        elseif (lmax .ge. 1) then
          dudr = facdu*wk(ir) - sre*(um+up)
          xi(ir,1) = (xi(ir,0) - dudr)/r2
        endif
        wk2(ir) = facgl*wk(ir)
   20 continue

C --- xi(ir,l) for l>1 by upward recursion ---
      facgl = 2*a**2
      do  30  l = lmin+2, lmax
        xx = 2*l-1
        do  7  ir = 1, nr
          xi(ir,l) = (xx*xi(ir,l-1) - e*xi(ir,l-2) - wk2(ir))/rsq(ir)
          wk2(ir) = facgl*wk2(ir)
    7   continue
   30 continue
      end

