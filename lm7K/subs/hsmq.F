      subroutine hsmq(nxi,lmxa,lxi,exi,rsm,job,q,p,nrx,nlmx,wk,yl,
     .  awald,alat,qlv,nG,dlv,nD,vol,hsm,hsmp)
C- Bloch-sum of smooth Hankel functions and energy derivatives at p
C  by Ewald summation, for nxi groups of parameters (lxi,exi,rsm).
C ---------------------------------------------------------------
Ci Inputs:
Ci  nxi,lxi,exi,rsm:  number of (l,energy,rsm) groups, and list
Ci  dlv,nD:  direct lattice vectors, and number
Ci  qlv,nG:  reciprocal lattice vectors, and number
Ci  awald,vol: Ewald parameter and cell volume (used for Q-space sum)
Ci  lmxa:    make strux for lxi(ie)+lmxa,
Ci           skipping over any energy for which lxi(ie) < 0
Ci  job:     ones digit intended for H_smoothed - H_unsmoothed
Ci           To do, set nG=0, and the one's digit as follows:
Ci            0, generate wk(1,2,3,4)
Ci           >0  assume wk(1,3,4) and yl have been generated already
Ci           tens digit handles whether or not to add H(p=0).
Ci            0: return whole Bloch function (not allowed for rsm=0)
Ci            1: return Bloch function less H(p=0), if p=0
Ci           >1: return Bloch function less H(p=0), for any p
Ci               DOES NOT WORK!
Ci           100s digit
Ci            0: use rsm = rsm(ie)
Ci            1: use rsm = rsm(1) is the same for all energies
Ci           1000s digit
Ci            0: use standard phase convention, do not shorten p
Ci            1: use standard phase convention, but shorten p
Ci            2: scale standard phase convention by exp(-i q . p)
Ci            3: like 2, but also shorten p
Ci  q,p:     k-point and connecting vector (units of 2*pi/alat and alat)
Ci  nrx:     leading dimension of wk,yl; must be at least max(nD,nG)
Ci  nlmx:    leading dimension of hsm,hsmp
Ci  yl:      work array, dimensioned nrx*(lmax+1)**2
Ci  wk:      work array, dimensioned at least nrx*(2*lmax+10)
Co Outputs:
Co   yl:     ylm(1..nD,1..(lmax+1)**2) for points alat*(p-dlv(1..nD))
Co   wk:     (*,1) holds r**2
Co           (*,2) holds Y0 exp(-(r/rsm(nxi))**2) (rsm>0 only)
Co           (*,3) holds cos(q.dlv)
Co           (*,4) holds sin(q.dlv)
Co   hsm,hsmp:smoothed Hankels and energy derivatives for
Co           nxi energies, to lmax lxi(1..nxi).
Cb Bugs and Limitations:
Cb   Convergence is poor for p small but nonzero.  To avoid this,
Cb   hsmq SHOULD re-evaluate points for which p is small using
Cb   a power series expansion.
Cu Updates
Cu  15 Aug 00 extended to e>0; added 1000 digit job
Cu   7 May 98 handles case rsm>1/a by strict q-space summation
cu  16 Sep 98 fixed dimensioning error arising when nD > nG
C ---------------------------------------------------------------
C     implicit none
      integer nxi,lmxa,lmax,nG,nD,nlmx,nrx,lxi(nxi),job
      double precision alat,awald,vol,exi(nxi),rsm(nxi),p(3),q(3),
     .  wk(nrx,10),yl(nrx,1),qlv(3,1),dlv(3,1),
     .  hsm(2,nlmx,1),hsmp(2,nlmx,1)
C Local variables
      integer ie,ir,ilm,l,m,ir1,lc,ls,job0,job1,job2,job3,lm,nlmxx,
     .  lx(20),ndx
      parameter (nlmxx=(16+1)**2)
      double precision qdotr,y0,a2,pi,sp,gam,tpiba,tpi,rsmi,
     .  x1,x2,xx,xx0,xx1(nlmxx),xx2(nlmxx),xx3(nlmxx),xx4(nlmxx),
     .  r,akap,kappa,h0,h0d,a,faca,p1(3),cosp,sinp
      parameter (faca=1d0)
      double complex xxc
      logical dcmpre,ltmp
      dcmpre(x1,x2) = dabs(x1-x2) .lt. 1d-8

C --- Setup ---
      job0 = mod(job,10)
      if (nG .ne. 0) job0 = 0
      job1 = mod(job/10,10)
      job2 = mod(job/100,10)
      job3 = mod(job/1000,10)
C ... Shorten connecting vector; adjust phase later
      if (job3 .eq. 1 .or. job3 .eq. 3) then
        call shortn(p,p1,dlv,nD)
      else
        call dcopy(3,p,1,p1,1)
      endif

      pi = 4*datan(1d0)
      tpi = 2*pi
C ... lx(ie) = lxi(ie) + lmxa
      lmax = -1
      if (nxi .gt. 20) call rx('hsmq: increase dim of lx')
      do  5  ie = 1, nxi
      lx(ie) = lxi(ie)+lmxa
    5 lmax = max(lmax,lx(ie))
C ... If all rsm are equal to 1/a, no real-space part
      a = awald
      nDx = nD
      if (nG .gt. 0 .and. a .ne. 0) then
      if (rsm(1) .gt. faca/a) then
        ltmp = .true.
        do  6  ie = 1, nxi
    6   ltmp = ltmp .and. abs(rsm(ie)-rsm(1)) .lt. 1d-12
        ltmp = ltmp .or. job2 .eq. 1
        if (ltmp) then
          a = 1/rsm(1)
          nDx = 1
        endif
      endif
      endif

      if ((lmax+1)**2 .gt. nlmx) call rx('hsmq: lxi gt ll(nlmx)')
      if ((lmax+1)**2 .gt. nlmxx) call rx('hsmq: lxi gt ll(nlmxx)')
      if (nrx .lt. max(nD,nG)) call rx('hsmq: nrx < nD or nG')
      y0 = 1/dsqrt(4*pi)
      a2 = a*a
      call dpzero(hsm, 2*nlmx*nxi)
      call dpzero(hsmp,2*nlmx*nxi)

C --- Energy-independent setup for Q-space part ---
      if (nG .gt. 0) then
        tpiba = 2*pi/alat
        gam = 0.25d0/a2
        do  10  ir = 1, nG
          wk(ir,2) = tpiba*(q(1) + qlv(1,ir))
          wk(ir,3) = tpiba*(q(2) + qlv(2,ir))
          wk(ir,4) = tpiba*(q(3) + qlv(3,ir))
   10   continue
        call ropyln(nG,wk(1,2),wk(1,3),wk(1,4),lmax,nrx,yl,wk)
        do  12  ir = 1, nG
          sp = alat*(wk(ir,2)*p1(1) + wk(ir,3)*p1(2) + wk(ir,4)*p1(3))
          xxc = cdexp(dcmplx(-gam*wk(ir,1), sp))
          wk(ir,2) = -dble(xxc)
          wk(ir,3) = -dimag(xxc)
   12   continue

C   ... Q-space part of reduced strx for all energies
        call pvhsmq(0,nxi,lmxa,lx,exi,a,vol,nG,wk,nrx,yl,nlmx,
     .    wk(1,2),wk(1,3),wk(1,4),wk(1,5),wk(1,6),wk(1,7),hsm,hsmp)
      endif

C --- Energy-independent setup for direct-space part ---
      if (job0 .eq. 0) then
C   ... Put ylm in yl and alat**2*(p-dlv)**2 in wk(1)
        do  20  ir = 1, nD
          wk(ir,2) = alat*(p1(1)-dlv(1,ir))
          wk(ir,3) = alat*(p1(2)-dlv(2,ir))
          wk(ir,4) = alat*(p1(3)-dlv(3,ir))
   20   continue
        call ropyln(nD,wk(1,2),wk(1,3),wk(1,4),lmax,nrx,yl,wk)
C   ... cos(q.dlv), sin(q.dlv) -> wk(3,4), Y0 exp(-(a dlv)**2) -> wk(2)
        if (dabs(q(1))+dabs(q(2))+dabs(q(3)) .ne. 0) then
          do  22  ir = 1, nD
            qdotr = tpi*(q(1)*dlv(1,ir)+ q(2)*dlv(2,ir)+ q(3)*dlv(3,ir))
            wk(ir,3) = dcos(qdotr)
            wk(ir,4) = dsin(qdotr)
   22     continue
        else
          do  23  ir = 1, nD
   23     wk(ir,3) = 1
          call dpzero(wk(1,4),nD)
        endif
      endif
C ... If we are doing Ewald sums, we need this to make H(1/a) (hansr4)
      if (nG .ne. 0) then
        do  24  ir = 1, nDx
   24   wk(ir,6) = y0*dexp(-wk(ir,1)*a2)
      endif

C --- D-space part of reduced strx for each energy: chi(l)=wk(l+lc)---
      xx0 = 0
      do  30  ie = 1, nxi
        if (lx(ie) .lt. lmxa) goto 30
        ir1 = 1
        if (job2 .ne. 0) then
          rsmi = rsm(1)
        else
          rsmi = rsm(ie)
        endif
C   ... Case connecting vector p=0
        if (dcmpre(wk(1,1),0d0)) then
          if (job1 .eq. 1) ir1 = 2
          if (job1 .eq. 0 .and. dcmpre(rsmi,0d0))
     .        call rx('hsmq: 10s digit job=0 and rsm=0 not allowed')
        endif
        if (job1 .gt. 1) ir1 = 2
        if (lx(ie) .gt. lmax) call rx('hsmq: lxi > lmax')
        lc = 7
        ls = lx(ie)+9

C  ...  Check whether to omit r.s. part
        nDx = nD
        if (nG .gt. 0) then
          if (abs(rsm(ie)-1/a) .lt. 1d-12) nDx = 1
        endif
C        if (ndx .eq. 1) then
C          print *, 'hi',rsm,1/awald
C        endif
C   ... chi = H(rsm) - H(1/a) = (H(0) - H(1/a)) - (H(0) - H(rsm))
        if (nG .gt. 0) then
          call hansr4(wk,lx(ie),nrx,nDx,exi(ie),1/a,wk(1,6),
     .      wk(1,5),wk(1,lc))
        else
          do  49  l = 0, lx(ie)+1
          do  49  ir = ir1, nDx
   49     wk(ir,l+lc) = 0
        endif
        if (.not. dcmpre(rsmi,0d0)) then
          xx = 1/rsmi**2
C     ... Remake Y0 exp(-(r/rsm)**2) if rsm has changed
          if (xx .ne. xx0) then
            xx0 = xx
            do  33  ir = 1, nDx
   33       wk(ir,2) = y0*dexp(-wk(ir,1)*xx)
          endif
          call hansr4(wk,lx(ie),nrx,nDx,exi(ie),rsmi,wk(1,2),
     .      wk(1,5),wk(1,ls))
C         call prm('H(0) - H(rsm)',wk(1,ls),nrx,nDx,lx(ie)+2)
          do  39  l = 0, lx(ie)+1
          do  39  ir = ir1, nDx
   39     wk(ir,l+lc) = wk(ir,l+lc) - wk(ir,l+ls)
        endif

C  --- Special treatment of on-site term ---
C      Subtract Ewald contribution from G-vectors, hsm(p,a).
C      Already generated by hansr4:
C      Case p eq 0, wk(1,lc+1) = -hsm(p->0,a)
C                   wk(1,lc)   = -2*hsmp(p->0,a)
C      Case p ne 0, wk(1,lc+1) = hsm(p,rsm->0) - hsm(p,a)
C                   wk(1,lc)   = 2*hsmp(p,rsm->0) - 2*hsmp(p,a)
C   ... For p>0,job1=2, convert w(1,lc..lc+1) to -hsm and -2*hsmp
        if (ir1 .eq. 2 .and. .not. dcmpre(wk(1,1),0d0)) then
C     ... Subtract 2*h0d from wk(1,lc), making it into -2*hsmp(p,a)
C         And h0 from wk(1,lc+1), making it into -hsm(p,a)
          r = dsqrt(wk(1,1))
          if (exi(ie) .gt. 0) then
            call rx(' hsmq: unfinished job1=2 implementation for e>0')
            kappa = dsqrt(exi(ie))
            h0 = cos(kappa*r)/r
            h0d = sin(kappa*r)/r
            wk(1,lc+1) = wk(1,lc+1) - h0
            wk(2,lc+1) = wk(2,lc+1) - h0d
C           wk(1,lc)   = wk(1,lc) -
C           wk(2,lc)   = wk(2,lc) -
          else
            akap = dsqrt(-exi(ie))
            h0 = exp(-akap*r)/r
            h0d = h0*r/akap
            wk(1,lc+1) = wk(1,lc+1) - h0
            wk(1,lc)   = wk(1,lc) - h0d
          endif
C     ... Generate -hsm(p,a,l>0) by upward recursion
          xx = 4*a*y0*dexp(-(akap/a/2)**2 - (r*a)**2)
          do  31  l = 1, lx(ie)
          wk(1,l+lc+1) = ((2*l-1)*wk(1,l+lc) -exi(ie)*wk(1,l+lc-1) + xx)
     .                   /wk(1,1)
c          print *, wk(1,l+lc+1)
   31     xx = 2*a**2*xx
C     ... debugging check
C         call hansmr(r,exi(ie),a,xx1,lx(ie))
C         call dpcopy(p1,xx3,1,3,alat)
C         call hsm_mol(xx3,1/a,exi(ie),lmax,xx1,xx2)
C         stop
        endif
C   ... debugging check
C       call hansmr(0d0,exi(ie),a,xx1,lx(ie))

C   ... Make sin(qR)*(H(rsm,r)-H(1/a,r)), cos(qR)*(H(rsm,r)-H(1/a,r))
        do  34  l = 0, lx(ie)+1
          do  35  ir = 1, nDx
   35     wk(ir,l+ls) = wk(ir,4)*wk(ir,l+lc)
          do  36  ir = 1, nDx
   36     wk(ir,l+lc) = wk(ir,3)*wk(ir,l+lc)
   34   continue
        do  32  l = 0, lx(ie)
        lm = l*l
C     ... xx1..4 artifically m-dependent to allow unrolling of m loop
          do  38  m = 1, 2*l+1
            ilm = lm+m
            xx1(m) = 0
            xx2(m) = 0
            xx3(m) = 0
            xx4(m) = 0
            do  37  ir = 1, nDx
              xx1(m) = xx1(m) + yl(ir,ilm)*wk(ir,l+lc+1)
              xx2(m) = xx2(m) + yl(ir,ilm)*wk(ir,l+ls+1)
              xx3(m) = xx3(m) + yl(ir,ilm)*wk(ir,l+lc)
              xx4(m) = xx4(m) + yl(ir,ilm)*wk(ir,l+ls)
   37       continue
            hsm(1,ilm,ie)  = (hsm(1,ilm,ie) + xx1(m))
            hsm(2,ilm,ie)  = (hsm(2,ilm,ie) + xx2(m))
            hsmp(1,ilm,ie) = (hsmp(1,ilm,ie) + xx3(m)/2)
            hsmp(2,ilm,ie) = (hsmp(2,ilm,ie) + xx4(m)/2)
   38     continue
   32   continue

C   ... Put in phase to undo shortening, or different phase convention
        sp = tpi*(q(1)*(p(1)-p1(1))+q(2)*(p(2)-p1(2))+q(3)*(p(3)-p1(3)))
        if (job3 .ge. 2) sp = sp-tpi*(q(1)*p1(1)+q(2)*p1(2)+q(3)*p1(3))
        if (sp .ne. 0d0) then
          cosp = dcos(sp)
          sinp = dsin(sp)
          do  40  ilm = 1, (lx(ie)+1)**2
            x1 = hsm(1,ilm,ie)
            x2 = hsm(2,ilm,ie)
            hsm(1,ilm,ie) = x1*cosp - x2*sinp
            hsm(2,ilm,ie) = x2*cosp + x1*sinp
            x1 = hsmp(1,ilm,ie)
            x2 = hsmp(2,ilm,ie)
            hsmp(1,ilm,ie) = x1*cosp - x2*sinp
            hsmp(2,ilm,ie) = x2*cosp + x1*sinp
   40     continue
        endif

   30 continue


      end
      subroutine hsmqe0(lmax,rsm,job,q,p,nrx,nlmx,wk,yl,
     .  awald,alat,qlv,nG,dlv,nD,vol,hsm)
C- Bloch-sum of smooth Hankel functions for energy 0
C ---------------------------------------------------------------
Ci Inputs:
Ci  rsm,lmax:smoothing radius, l cutoff for hsm
Ci  dlv,nD:  direct lattice vectors, and number
Ci  qlv,nG:  reciprocal lattice vectors, and number
Ci  a,vol:   Ewald parameter and cell volume (used for Q-space sum)
Ci  job:     ones digit intended for H_smoothed - H_unsmoothed
Ci           To do, set nG=0, and the one's digit as follows:
Ci            0, generate wk(1,2,3,4)
Ci           >0  assume wk(1,3,4) and yl have been generated already
Ci           tens digit handles whether or not to add H(p=0).
Ci            0: return whole Bloch function (not allowed for rsm=0)
Ci            1: return Bloch function less H(p), if p=0
Ci           >1: return Bloch function less H(p), for any p
Ci           100s digit
Ci           >0: ignore passed q, use q=0
Ci           1000s digit
Ci            0 use standard phase convention, do not shorten p
Ci            1 use standard phase convention, but shorten p
Ci            2 scale standard phase convention by exp(-i q . p)
Ci            3 like 2, but also shorten p
Ci  q,p:     k-point and connecting vector (units of 2*pi/alat and alat)
Ci  nrx:     leading dimension of wk,yl
Ci  nlmx:    leading dimension of hsm
Ci  yl:      work array, dimensioned nrx*(lmax+1)**2
Ci  wk:      work array, dimensioned at least nrx*(2*lmax+10)
Co Outputs:
Co   yl:     ylm(1..nD,1..(lmax+1)**2) for points alat*(p-dlv(1..nD))
Co   wk:     (*,1) holds r**2
Co           (*,2) holds Y0 exp(-(r/rsm)**2) (rsm>0 only)
Co           (*,3) holds cos(q.dlv)
Co           (*,4) holds sin(q.dlv)
Co   hsm:    smoothed Hankels for e=0,q=0 for l=0..lmax
Cr Remarks
Cr   hsmqe0 is an adaptation of hsmq for for e=0, and possibly q=0.
Cr   For l=0, hsm(q=0,e->0) diverges.  To achieve a finite hsm(l=0,e=0),
Cr   the average lim(e->0) hsm0(e) = -sqrt(4*pi)*exp(gamma*e) / (vol*e)
Cr   is subtracted from hsm(l=0).
Cb Bugs and Limitations:
Cb   see hsmq.
Ce External routines required: hansr5, ropyln, dpzero
Cu Updates
Cu  25 Jun 03 (Kino) call pvhsmq with arrays for lmax, e=0
Cu   7 May 98 handles case rsm>1/a by strict q-space summation
C ---------------------------------------------------------------
C     implicit none
      integer lmax,nG,nD,nlmx,nrx,job
      double precision alat,awald,vol,rsm,p(3),q(3),
     .  wk(nrx,10),yl(nrx,1),qlv(3,1),dlv(3,1),hsm(2,nlmx)
C Local variables
      integer ir,ilm,l,m,ir1,lc,ls,job0,job1,job2,job3,lm,nlmxx,nDx,
     .  lx1(1)
      parameter (nlmxx=(16+1)**2)
      double precision qdotr,y0,a2,pi,sp,gam,tpiba,tpi,p1(3),ex1(1),
     .  x1,x2,xx,xx1(nlmxx),xx2(nlmxx),r,h0,q0(3),a,faca,cosp,sinp
      parameter (faca=1d0)
      double complex xxc
      logical dcmpre,lqzero
      dcmpre(x1,x2) = dabs(x1-x2) .lt. 1d-8

C --- Setup ---
      job0 = mod(job,10)
      if (nG .ne. 0) job0 = 0
      job1 = mod(job/10,10)
      job2 = mod(job/100,10)
      job3 = mod(job/1000,10)
      if ((lmax+1)**2 .gt. nlmx) call rx('hsmqe0: lmax gt ll(nlmx)')
      if ((lmax+1)**2 .gt. nlmxx) call rx('hsmqe0: lmax gt ll(nlmxx)')
      if (nrx .lt. max(nD,nG)) call rx('hsmqe0: nrx < nD or nG')
C ... Shorten connecting vector; adjust phase later
      if (job3 .eq. 1 .or. job3 .eq. 3) then
        call shortn(p,p1,dlv,nD)
      else
        call dcopy(3,p,1,p1,1)
      endif
      pi = 4*datan(1d0)
      y0 = 1/dsqrt(4*pi)
      tpi = 2*pi
      call dpzero(hsm, 2*nlmx)
      if (job2 .gt. 0) then
        call dpzero(q0,3)
      else
        call dpcopy(q,q0,1,3,1d0)
      endif
      lqzero = dcmpre(dabs(q0(1))+dabs(q0(2))+dabs(q0(3)),0d0)
C ... If rsm ge 1/a, set a = 1/rsm and skip r.s. part
      a = awald
      nDx = nD
      if (nG .gt. 0) then
        if (rsm .gt. faca/a) nDx = 1
        if (nDx .eq. 1) then
          a = 1/rsm
        endif
      endif
      a2 = a*a

C --- Setup for Q-space part ---
      if (nG .gt. 0) then
        tpiba = 2*pi/alat
        gam = 0.25d0/a2
        do  10  ir = 1, nG
          wk(ir,2) = tpiba*(q0(1) + qlv(1,ir))
          wk(ir,3) = tpiba*(q0(2) + qlv(2,ir))
          wk(ir,4) = tpiba*(q0(3) + qlv(3,ir))
   10   continue
        call ropyln(nG,wk(1,2),wk(1,3),wk(1,4),lmax,nrx,yl,wk)
        do  12  ir = 1, nG
          sp = alat*(wk(ir,2)*p1(1) + wk(ir,3)*p1(2) + wk(ir,4)*p1(3))
          xxc = cdexp(dcmplx(-gam*wk(ir,1), sp))
          wk(ir,2) = -dble(xxc)
          wk(ir,3) = -dimag(xxc)
   12   continue

C   --- Q-space part hsm(1/a,e,r) of reduced strx for all energies ---
C   ... bug fix: copy scalars to arrays for call (Kino)
        lx1(1) = lmax
        ex1(1) = 0d0
        call pvhsmq(1,1,lmax,lx1,ex1,a,vol,nG,wk,nrx,yl,nlmx,
     .    wk(1,2),wk(1,3),wk(1,4),wk(1,5),wk(1,6),wk(1,7),hsm,hsm)
      endif

C --- Energy-independent setup for direct-space part ---
      if (job0 .eq. 0) then
C   ... Put ylm in yl and alat**2*(p-dlv)**2 in wk(1)
        do  20  ir = 1, nDx
          wk(ir,2) = alat*(p1(1)-dlv(1,ir))
          wk(ir,3) = alat*(p1(2)-dlv(2,ir))
          wk(ir,4) = alat*(p1(3)-dlv(3,ir))
   20   continue
        call ropyln(nDx,wk(1,2),wk(1,3),wk(1,4),lmax,nrx,yl,wk)
C   ... cos(q.dlv), sin(q.dlv) -> wk(3,4), Y0 exp(-(a dlv)**2) -> wk(2)
        if (.not. lqzero) then
          do  22  ir = 1, nDx
            qdotr =tpi*(q0(1)*dlv(1,ir)+q0(2)*dlv(2,ir)+q0(3)*dlv(3,ir))
            wk(ir,3) = dcos(qdotr)
            wk(ir,4) = dsin(qdotr)
   22     continue
        else
          do  23  ir = 1, nDx
   23     wk(ir,3) = 1
          call dpzero(wk(1,4),nDx)
        endif
      endif
C ... If we are doing Ewald sums, we need this to make H(1/a) (hansr5)
      if (nG .ne. 0) then
        do  24  ir = 1, nDx
   24   wk(ir,6) = y0*dexp(-wk(ir,1)*a2)
      endif

C --- D-space part of reduced strx for each energy: chi(l)=wk(l+lc)---
      ir1 = 1
C ... Case connecting vector p=0
      if (dcmpre(wk(1,1),0d0)) then
        if (job1 .eq. 1) ir1 = 2
        if (job1 .eq. 0 .and. dcmpre(rsm,0d0))
     .    call rx('hsmqe0: job=0 and rsm=0 not allowed')
      endif
      if (job1 .gt. 1) ir1 = 2
      lc = 7
      ls = lmax+8

C ... chi = H(rsm) - H(1/a) = (H(0) - H(1/a)) - (H(0) - H(rsm))
      if (nG .gt. 0) then
        call hansr5(wk,lmax,nrx,nDx,1/a,wk(1,6),wk(1,5),wk(1,lc))
      else
        do  49  l = 0, lmax+1
        do  49  ir = ir1, nDx
   49   wk(ir,l+lc) = 0
      endif
      if (.not. dcmpre(rsm,0d0)) then
        xx = 1/rsm**2
C   ... Make Y0 exp(-(r/rsm)**2)
        do  33  ir = 1, nDx
   33   wk(ir,2) = y0*dexp(-wk(ir,1)*xx)
        call hansr5(wk,lmax,nrx,nDx,rsm,wk(1,2),wk(1,5),wk(1,ls))
C       call prm('H(0) - H(rsm)',wk(1,ls),nrx,nDx,lmax+1)
        do  39  l = 0, lmax
        do  39  ir = ir1, nDx
   39   wk(ir,l+lc) = wk(ir,l+lc) - wk(ir,l+ls)
      endif

C --- Special treatment of on-site term ---
C     Subtract Ewald contribution from G-vectors, hsm(p,a).
C     Already generated by hansr5:
C     Case p eq 0, wk(1,lc) = -hsm(p->0,a)
C     Case p ne 0, wk(1,lc) =  hsm(p,rsm->0) - hsm(p,a)
C ... For p>0,job=2, convert w(1,lc) to -hsm
      if (ir1 .eq. 2 .and. .not. dcmpre(wk(1,1),0d0)) then
C   ... Subtract h0(0..1) from wk(1,lc..lc+1), making it into -hsm(p,a)
        r = dsqrt(wk(1,1))
        h0 = 1/r
        wk(1,lc) = wk(1,lc) - h0
        if (lmax .gt. 0) wk(1,lc+1) = wk(1,lc+1) - h0/wk(1,1)
C   ... Generate -hsm(p,a,l>0) by upward recursion
        xx = (2*a**2)*4*a*y0*dexp(-(r*a)**2)
C  ... debugging check
C      call hansmr(r,0d0,a,xx1,lmax)
        do  31  l = 2, lmax
        wk(1,l+lc) = ((2*l-1)*wk(1,l+lc-1) + xx)/wk(1,1)
C       print *, wk(1,l+lc)
   31   xx = 2*a**2*xx
      endif

C ... Make sin(qR)*(H(rsm,r)-H(1/a,r)), cos(qR)*(H(rsm,r)-H(1/a,r))
      if (.not. lqzero) then
        do  32  l = 0, lmax
        do  35  ir = 1, nDx
   35   wk(ir,l+ls) = wk(ir,4)*wk(ir,l+lc)
        do  34  ir = 1, nDx
   34   wk(ir,l+lc) = wk(ir,3)*wk(ir,l+lc)
   32   continue
      endif
      sp = tpi*(q(1)*(p(1)-p1(1))+q(2)*(p(2)-p1(2))+q(3)*(p(3)-p1(3)))
      if (job3 .ge. 2) sp = sp-tpi*(q(1)*p1(1)+q(2)*p1(2)+q(3)*p1(3))
      cosp = dcos(sp)
      sinp = dsin(sp)
      do  36  l = 0, lmax
        lm = l*l
C   ... xx1..4 artifically m-dependent to allow unrolling of m loop
        if (lqzero) then
          do  138  m = 1, 2*l+1
            ilm = lm+m
            xx1(m) = 0
            xx2(m) = 0
            do  137  ir = 1, nDx
              xx1(m) = xx1(m) + yl(ir,ilm)*wk(ir,l+lc)
  137       continue
            hsm(1,ilm)  = (hsm(1,ilm) + xx1(m))
            hsm(2,ilm)  = 0d0
  138     continue
        else
          do  38  m = 1, 2*l+1
            ilm = lm+m
            xx1(m) = 0
            xx2(m) = 0
            do  37  ir = 1, nDx
              xx1(m) = xx1(m) + yl(ir,ilm)*wk(ir,l+lc)
              xx2(m) = xx2(m) + yl(ir,ilm)*wk(ir,l+ls)
   37       continue
            hsm(1,ilm)  = (hsm(1,ilm) + xx1(m))
            hsm(2,ilm)  = (hsm(2,ilm) + xx2(m))
   38     continue
        endif

C   ... Put in phase to undo shortening
        if (sp .ne. 0d0) then
          do  40  m = 1, 2*l+1
            ilm = lm+m
            x1 = hsm(1,ilm)
            x2 = hsm(2,ilm)
            hsm(1,ilm) = x1*cosp - x2*sinp
            hsm(2,ilm) = x2*cosp + x1*sinp
   40     continue
        endif

   36 continue

C --- Add extra term for l=0 when e=0 and q=0 ---
      if (lqzero) hsm(1,1) = hsm(1,1) + rsm**2/(4*vol*y0)

      end
      subroutine pvhsmq(job,nx,lmxa,lxi,exi,a,vol,n,qsq,nyl,yl,nlmx,
     .  cs,sn,cs2,sn2,cs3,sn3,hl,hlp)
C- Makes H_L(1/a,e,r) by summation over reciprocal lattice vectors.
C  job=1 : make hl but not hlp
C  Skips any ie for which lxi(ie)<lmxa.
C  This is a kernel called by hsmq.
C     implicit none
      integer n,nx,lmxa,lxi(nx),nlmx,job,nyl
      double precision cof,e,gam,pi,vol
      double precision cs(n),sn(n),cs2(n),sn2(n),cs3(n),sn3(n),
     .  a,exi(1),qsq(n),yl(nyl,1),hl(2,nlmx,1),hlp(2,nlmx,1)
C Local variables
      integer i,i1,ie,ilm,l,lx,m,je,lm,nlmxx
      parameter (nlmxx=(16+1)**2)
      double precision c,s,xx,x1,x2
      double precision c2(nlmxx),s2(nlmxx),c3(nlmxx),s3(nlmxx)
      logical dcmpre
      dcmpre(x1,x2) = dabs(x1-x2) .lt. 1d-8

      gam = 0.25d0/a**2
      pi = 4*datan(1d0)

C --- For each exi ne 0 and ilm do ---
      do  20  ie = 1, nx
        e = exi(ie)
        lx = lxi(ie)
        if (lxi(ie) .lt. lmxa) goto 20

C   ... Copy existing hl,hlp if already calculated for this exi
        do  30  je = 1, ie-1
          if (dabs(exi(je)-exi(ie)) .gt. 1d-8) goto 30
          if (lxi(ie) .gt. lxi(je)) goto 30
C     ... We have a match:
          call dpadd(hl(1,1,ie), hl(1,1,je), 1,2*(lxi(ie)+1)**2,1d0)
          call dpadd(hlp(1,1,ie),hlp(1,1,je),1,2*(lxi(ie)+1)**2,1d0)
          goto 20
   30   continue

        i1 = 1
        if (dcmpre(e-qsq(1),0d0)) i1 = 2
        do  22  i = i1, n
          xx = 1/(e-qsq(i))
          c = cs(i)*xx
          cs2(i) = c
          cs3(i) = c*xx
          s = sn(i)*xx
          sn2(i) = s
          sn3(i) = s*xx
   22   continue
        cof = -4d0*pi*dexp(gam*e)/vol
        ilm = 0
        do  21  l = 0, lx, 2
          cof = -cof
          lm = l*l
C     ... Do the even and odd l's separately
C     ... NB c2,s2,c3,s3 m-dependent to allow m loop to be unrolled
          if (job .eq. 0) then
          do  24  m = 1, 2*l+1
c           ilm = ilm+1
            ilm = lm+m
            c2(ilm) = 0
            s2(ilm) = 0
            c3(ilm) = 0
            s3(ilm) = 0
            do  23  i = i1, n
              c2(ilm) = c2(ilm) + yl(i,ilm)*cs2(i)
              c3(ilm) = c3(ilm) + yl(i,ilm)*cs3(i)
              s2(ilm) = s2(ilm) + yl(i,ilm)*sn2(i)
              s3(ilm) = s3(ilm) + yl(i,ilm)*sn3(i)
   23       continue
            hl(1,ilm,ie)  = hl(1,ilm,ie)  + cof*c2(ilm)
            hlp(1,ilm,ie) = hlp(1,ilm,ie) + cof*(gam*c2(ilm)-c3(ilm))
            hl(2,ilm,ie)  = hl(2,ilm,ie)  + cof*s2(ilm)
            hlp(2,ilm,ie) = hlp(2,ilm,ie) + cof*(gam*s2(ilm)-s3(ilm))
   24     continue
          else
          do  124  m = 1, 2*l+1
c           ilm = ilm+1
            ilm = lm+m
            c2(ilm) = 0
            s2(ilm) = 0
            do  123  i = i1, n
              c2(ilm) = c2(ilm) + yl(i,ilm)*cs2(i)
              s2(ilm) = s2(ilm) + yl(i,ilm)*sn2(i)
  123       continue
            hl(1,ilm,ie)  = hl(1,ilm,ie)  + cof*c2(ilm)
            hl(2,ilm,ie)  = hl(2,ilm,ie)  + cof*s2(ilm)
  124     continue
          endif
          if (l+1 .le. lx) then
            lm = (l+1)*(l+1)
            if (job .eq. 0) then
            do  26  m = 1, 2*l+3
c             ilm = ilm+1
              ilm = lm+m
              c2(ilm) = 0
              s2(ilm) = 0
              c3(ilm) = 0
              s3(ilm) = 0
              do  27  i = i1, n
                c2(ilm) = c2(ilm) + yl(i,ilm)*cs2(i)
                c3(ilm) = c3(ilm) + yl(i,ilm)*cs3(i)
                s2(ilm) = s2(ilm) + yl(i,ilm)*sn2(i)
                s3(ilm) = s3(ilm) + yl(i,ilm)*sn3(i)
   27         continue
              hl(1,ilm,ie)  = hl(1,ilm,ie)  + cof*s2(ilm)
              hlp(1,ilm,ie) = hlp(1,ilm,ie) + cof*(gam*s2(ilm)-s3(ilm))
              hl(2,ilm,ie)  = hl(2,ilm,ie)  - cof*c2(ilm)
              hlp(2,ilm,ie) = hlp(2,ilm,ie) - cof*(gam*c2(ilm)-c3(ilm))
   26       continue
          else
            do  126  m = 1, 2*l+3
c             ilm = ilm+1
              ilm = lm+m
              c2(ilm) = 0
              s2(ilm) = 0
              do  127  i = i1, n
                c2(ilm) = c2(ilm) + yl(i,ilm)*cs2(i)
                s2(ilm) = s2(ilm) + yl(i,ilm)*sn2(i)
  127         continue
              hl(1,ilm,ie)  = hl(1,ilm,ie)  + cof*s2(ilm)
              hl(2,ilm,ie)  = hl(2,ilm,ie)  - cof*c2(ilm)
  126       continue
          endif
          endif
   21   continue

C --- Add extra term for l=0 when e=0 and q=0 ---
      if (dcmpre(e-qsq(1),0d0)) hl(1,1,ie) = hl(1,1,ie) -
     .    sqrt(4d0*pi)/vol*gam
   20 continue
      end

