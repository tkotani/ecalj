      subroutine mktra2(job,loka,nbas,ipc,nl,lmx,avw,itrans,kap2,nkap,
     .  hcr,nder,tral,trad,alpha,adot,cd,bigd)
C- Makes screening transformation matrix tral and related quantities
C ----------------------------------------------------------------------
Ci Inputs:
Ci   job   :1s digit: make tral,trad
Ci         :10s digit: make alpha,adot
Ci         :100s digit: make cd
Ci         :1000s digit: make bigd
Ci         :Any combination of these digits is allowed
Ci   nbas  :size of basis
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   nl    :(global maximum l) + 1
Ci   lmx   :lmx(j) = maximum l for atom j
Ci   avw   :average Wigner-Seitz sphere radius
Ci   itrans:characterizes structure matrix transformation (see Remarks)
Ci          Let a = hard core radius, a=hcr*avw
Ci          See remarks for notation
Ci
Ci          0:2nd generation transformation:
Ci           head  |N^a(kappa)> =|N^0(kappa)>
Ci                 |J^a(kappa)> =|J^0(kappa)>-alpha|K^0(kappa)>
Ci                 |J^a(kappa)> = 0 at a
Ci                 Wronskian W{N^a,J^a} = W{N^0,J^0} = avw/2
Ci                 NB: in 2gen LMTO, S generated is the structure matrix
Ci                 and is the analog of B in NMTO notation
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
Ci                 NB: generates the standard NMTO structure matrix,
Ci                 apart from scaling conventions.
Ci                 See also Remarks.
Ci          5:head |N^a(kappa)> has value 1 and slope 0
Ci                 |J^a(kappa)> has value 0 and slope -1/a at a
Ci          6:head |N^a(kappa)> has value |N^0(0)> and slope 0
Ci                 |J^a(kappa)> = 0 at a
Ci                 Wronskian W{N^a,J^a} = W{N^0,J^0} = avw/2
Ci   kap2  :(kappa*avw)**2
Ci   nkap  :number of energies for which strux are calculated
Ci   hcr   :hard core screening radius a in units of avw, ie a/avw
Ci   nder  :number of energy derivatives
Co Outputs:
Co   tral  :transformation matrix for head and tail functions
Co   trad  :(kappa*avw)^2-derivative of tral
Co   alpha :2nd generation tight-binding screening constants
Co   adot  :(kappa*avw)^2-derivative of tight-binding screening const.
Co   cd    :tral(3)*tral(4), and energy derivatives of it
Co   bigd  :tral(1)/tral(3)
Co
Cl Local variables
Cl   rfacn :factor to scale besslr gi to spherical Neumann according to
Cl         :to supplied loka conventions
Cl   rfacj :factor to scale besslr fi to spherical Bessel according to
Cl         :to supplied loka conventions
Cr Remarks: Use this notation:
Cr          |f>  signifies function truncated outside its own sphere
Cr         ||f>  signifies function including all spheres.
Cr
Cr          The transformation matrix from 'bare' Neumann and Bessel
Cr          functions |N>,|J> to 'screened' functions |N^a>,|J^a> is
Cr          |N^a> = tral(1)*|N> + tral(2)*|J>
Cr          |J^a> = tral(3)*|N> + tral(4)*|J>
Cr          i.e. 'head' consists of tral(1)*N + tral(2)*J
Cr          and  'tail' consists of tral(3)*N + tral(4)*J
Cr
Cr         ||N^a> = |N^a> - |J^a>S^a
Cr          Thus, in the interstitital,
Cr         ||N^a>_j = tral(1)_j |N>_j - sum_i tral(3)_i * S^a_ij
Cr
Cr        See R. Tank, Phys. Stat. Sol. B217, 89 (2000)
Cr        His Eq. (18) uses itrans=5
Cr
Cr  As currently written,
Cr  N and J follow Andersen's old definitions (see besslr.f).
Cr  NB: OKA's conventions for strux in NMTO reverted to more standard
Cr  ones, while here Andersen's original conventions are retained.
Cr
Cr  mktra2 was adapted from the Stuttgart LMTO, written by R. Tank
Cr  It is similar in function to mktral, but uses hcr as input
Cr  and generates extra information, including higher order derivatives
Cr  Stuttgart's sigma(l,ic)*wsr(ic) = hcr(l,ic)*avw here.
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer job,loka,itrans,lmx(*),nl,nbas,nder,ipc(nbas),nkap
      double precision avw,kap2(nkap),hcr(nl,*),
     .  alpha(nl*nl,nbas,nkap),adot(nl*nl,nbas,nkap,nder),
     .  trad(4,nl**2,nbas,nkap,nder),tral(4,nl**2,nbas,nkap),
     .  cd(nl**2,nbas,nkap,0:nder),bigd(nl**2,nbas,nkap,0:nder)
C Local variables:
      integer nlmx,ndmax
      parameter(nlmx=20,ndmax=6)
      integer i,ic,k,l,lgunit,jder,fib(0:ndmax,0:ndmax),ipr,ik,
     .        ib,lm,m,iclbsj
      double precision fac1,fac2,hcrad,r2,rfacn,rfacj
      double precision fi(-ndmax-1:nlmx+ndmax),n(0:ndmax),dn(0:ndmax),
     .                 gi(-ndmax-1:nlmx+ndmax),j(0:ndmax),dj(0:ndmax),
     .                 jm1(0:ndmax),nm1(0:ndmax),jntmp(0:ndmax),
     .                 al(0:ndmax),cm1(0:ndmax),t(4,0:ndmax),
     .                 tcd(0:ndmax)
      character*(54) strn(0:6)
      data strn /
     .  'N^a = N^0,  J^a = J^0 - alpha K^0',
     .  'N^a = N^0(e=0), J^a = 0 at hcr',
     .  'N^a has val,slo = N^0(e=0), J^a = 0 at hcr',
     .  'N^a has val,slo 1,0, J^a has val,slo 0, w/2a^2 at hcr',
     .  'N^a has val,slo 1,0, J^a has val,slo 0, 1/a^2 at hcr',
     .  'N^a has val,slo 1,0, J^a has val,slo 0, -1/a at hcr',
     .  'N^a has val=N^0(0), slo=0, J^a = 0 at hcr'/

C ... Factors needed for chain-rule differentiation of a product
      do  i = 0, ndmax
        do  k = 0, ndmax
          fib(i,k) = 0
        enddo
      enddo
      fib(0,0) = 1
      do  jder = 1, nder
        fib(jder,0) = 1
        do  i = 1, jder
          fib(jder,i) = fib(jder-1,i-1) + fib(jder-1,i)
        enddo
      enddo

      call getpr(ipr)
      if (ipr .ge. 50/1) write(lgunit(1),300) strn(itrans)
  300 format(/' MKTRA2:  ',a/'  ib ic  l     kap2       hcr',
     .  '         a           b           c           d')
      if (itrans .lt .0 .or. itrans .gt. 6)
     .  call rxi('MKTRA2: itrans out of range, value',itrans)
      do  ik = 1, nkap
      do  ib = 1, nbas
        ic = ipc(ib)
        lm = 0
        do  l = 0, nl-1
        do  m = -l, l
          lm = lm+1
          if (l .le. lmx(ic)) then
C         loka = 1 : hcrad = (a/w)
C         loka = 0 : hcrad = a
          hcrad = hcr(l+1,ic)
          if (loka .eq. 0) hcrad = hcr(l+1,ic)*avw
C         kap2*r2 = (a/w)^2(E*w^2), independent of avw
          r2 = hcr(l+1,ic)**2
          gi(l+1) = 0
          i = nder+1
C         Conventions for Hankels, Bessels (k=sqrt(E), w=avw)
C         OKA 2nd generation LMTO:
C         j_l = (a/w)^l fi^OKA(l) = 1/2*(2l-1)!! (kw)^(-l)j_l(ka)
C         h_l = (a/w)^(-1-l) gi^OKA(l) = -i/(2l-1)!! (kw)^(l+1) h_l(ka)
C         MSM standard definitions
C         j_l = fi(l) a^l = 1/(ik)^l j_l(ka)
C         h_l = gi^MSM(l)/a^(l+1) = -i k^(l+1) h_l(ka)
C         j/k differ by ratio 2w^(2l+1)/(2l-1)!!
          call besslr(kap2(ik)*r2,loka,-i,l+i,fi(-i),gi(-i))
          rfacn = hcrad**(-l-1)
          rfacj = hcrad**l

C     ... Make n_l and its energy derivatives
          n(0) = rfacn*gi(l)
          do  jder = 1, nder
            n(jder) = n(jder-1)*(r2/2d0)*(gi(l-jder)/gi(l-jder+1))/
     .        (2*(l-jder)+1)
          enddo

C     ... Hold for now n_(l+1) and its energy derivatives
          dn(0) = rfacn/hcrad*gi(l+1)
          do  jder = 1, nder
            dn(jder) = dn(jder-1)*(r2/2d0)*(gi(l-jder+1)/gi(l-jder+2))/
     .        (2*(l-jder)+3)
          enddo

C     ... Radial derivative dn_l/dr and its energy derivatives
          do  jder = 0, nder
            dn(jder) = (l*n(jder)/hcrad - (2*l+1)*dn(jder))/avw
          enddo

C     ... Inverse of n_l, and its energy derivatives
          nm1(0) = 1d0/n(0)
          do  jder = 1, nder
            jntmp(jder) = 0d0
            nm1(jder) = 0d0
            do  i = 1, jder
              jntmp(jder) = jntmp(jder)+n(i)*nm1(jder-i)*fib(jder-1,i-1)
            enddo
            do  i = 1, jder
              nm1(jder) = nm1(jder)-jntmp(i)*nm1(jder-i)*fib(jder-1,i-1)
            enddo
          enddo

C     ... j_l and its energy derivatives
          j(0) = rfacj*fi(l)
          do  jder = 1, nder
            j(jder) = -j(jder-1)*(r2/2d0)*(fi(l+jder)/fi(l+jder-1))/
     .        (2*(l+jder)-1)
          enddo

C     ... j_(l-1) and its energy derivatives
          dj(0) = rfacj/hcrad*fi(l-1)
          do  jder = 1, nder
            dj(jder) = -dj(jder-1)*(r2/2d0)*(fi(l+jder-1)/fi(l+jder-2))/
     .        (2*(l+jder)-3)
          enddo

C     ... dj_l/dr and its energy derivatives
          do  jder = 0, nder
            dj(jder) = (-(l+1)*j(jder)/hcrad +(2*l-1)*dj(jder))/avw
          enddo

C     ... Inverse of j_l, and its energy derivatives
          jm1(0) = 1d0/j(0)
          do  jder = 1, nder
            jntmp(jder) = 0d0
            jm1(jder) = 0d0
            do  i = 1, jder
              jntmp(jder) = jntmp(jder)+j(i)*jm1(jder-i)*fib(jder-1,i-1)
            enddo
            do  i = 1, jder
              jm1(jder) = jm1(jder)-jntmp(i)*jm1(jder-i)*fib(jder-1,i-1)
            enddo
          enddo

C     ... alpha and its derivatives
          al(0) = j(0)*nm1(0)
          do  jder = 1, nder
            al(jder) = 0d0
            do  i = 0, jder
              al(jder) = al(jder) + j(i)*nm1(jder-i)*fib(jder,i)
            enddo
          enddo

C     --- itrans-dependent quantities ---
C         W{n,j} = r*r*[n*j'-n'*j] = avw/2
C         => dj(0)*n(0)-j(0)*dn(0) = avw/2/a**2 = 1/fac1
C     ... => det(tral) =  fac1 fac2 (dj(0)*n(0)-j(0)*dn(0)) = fac2
C         Tank adopts itrans=5 in Phys. Stat. Sol. B217, 89 (2000)
          fac1 = 2d0*hcrad*hcrad*avw
          if (itrans .eq. 0) then
C           dfac(l,ic) = 1d0
          elseif (itrans .eq. 3) then
            fac2 = 1d0
C           dfac(l,ic) = 1d0
          elseif (itrans .eq. 4) then
            fac2 = 2d0/avw
C           dfac(l,ic) = 2d0/avw
          elseif (itrans .eq. 5) then
            fac2 = -2d0*hcrad
C           dfac(l,ic) = -2d0*hcrad
          endif

C     ... a,b,c,d and their derivatives, itrans=0
          if (itrans .eq. 0) then
            t(1,0) = 1d0
            t(2,0) = 0d0
            t(3,0) = -al(0)
            t(4,0) = 1d0
            do  jder = 1, nder
              t(1,jder) = 0d0
              t(2,jder) = 0d0
              t(3,jder) = -al(jder)
              t(4,jder) = 0d0
            enddo
C     ... Derivatives of c*d, itrans=0
            do  jder = 0, nder
              tcd(jder) = -al(jder)
            enddo

C     ... a,b,c,d and their derivatives
          else
            do  jder = 0, nder
              t(1,jder) = fac1*dj(jder)
              t(2,jder) = -fac1*dn(jder)
              t(3,jder) = -fac2*j(jder)
              t(4,jder) = fac2*n(jder)
            enddo

C       ... Derivatives of c*d
            tcd(0) = t(3,0)*t(4,0)
            do  jder = 1, nder
              tcd(jder) = 0d0
              do  i = 0, jder
                tcd(jder) = tcd(jder) + fib(jder,i)*t(3,jder-i)*t(4,i)
              enddo
            enddo
          endif

C     ... Inverse of c and its energy derivatives
          cm1(0) = 1d0/t(3,0)
          do  jder = 1, nder
            jntmp(jder) = 0d0
            cm1(jder) = 0d0
            do  i = 1, jder
              jntmp(jder) = jntmp(jder) + t(3,i)*cm1(jder-i)*
     .          fib(jder-1,i-1)
            enddo
            do  i = 1, jder
              cm1(jder) = cm1(jder)-jntmp(i)*cm1(jder-i)*fib(jder-1,i-1)
            enddo
          enddo

C     ... bigd = a/c  and energy derivatives
          if (mod(job/1000,10) .ne. 0) then
          jntmp(0) = t(1,0)*cm1(0)
          bigd(lm,ib,ik,0) = jntmp(0)
          do  jder = 1, nder
            jntmp(jder) = 0d0
            do  i = 0, jder
              jntmp(jder) = jntmp(jder) + fib(jder,i)*t(1,jder-i)*cm1(i)
            enddo
            bigd(lm,ib,ik,jder) = jntmp(jder)
          enddo
          endif

C     --- Copy to arrays ---
          if (mod(job,10) .ne. 0)
     .      call dcopy(4,t(1,0),1,tral(1,lm,ib,ik),1)
          if (mod(job/10,10) .ne. 0) alpha(lm,ib,ik) = al(0)
          if (mod(job/100,10) .ne. 0) cd(lm,ib,ik,0) = tcd(0)
          do  jder = 1, nder
            if (mod(job,10) .ne. 0)
     .        call dcopy(4,t(1,jder),1,trad(1,lm,ib,ik,jder),1)
            if (mod(job/10,10) .ne. 0) adot(lm,ib,ik,jder) = al(jder)
            if (mod(job/100,10) .ne. 0) cd(lm,ib,ik,jder) = tcd(jder)
          enddo

C          print *, 'xxx',
C     .      cd(lm,ib,ik,0)-tral(3,lm,ib,ik)*tral(4,lm,ib,ik)
C
C          print *, 'xxx',
C     .      cd(lm,ib,ik,1)-
C     .      trad(3,lm,ib,ik,1)*tral(4,lm,ib,ik)-
C     .      tral(3,lm,ib,ik)*trad(4,lm,ib,ik,1)

          if (ipr .ge. 50/1 .and. iclbsj(ic,ipc,nbas,1) .eq. ib
     .        .and. l .eq. m) then
            if (l .eq. 0 .and. ib .eq. 1) then
              print 302, ib,ic,l,kap2(ik)/avw**2,hcrad*avw,
     .                   (t(k,0),k=1,4)
            elseif (l .eq. 0) then
              print 303, ib,ic,l,hcrad*avw,(t(k,0),k=1,4)
            else
              print 304, l,hcrad*avw,(t(k,0),k=1,4)
            endif
            if (ipr .gt. 50) print 305, (t(k,1),k=1,4)
  302       format(i4,2i3,2f11.6,4f12.7)
  303       format(i4,2i3,11x,f11.6,4f12.7)
  304       format(4x,3x,i3,11x,f11.6,4f12.7)
  305       format('       trad',10x,f11.6,4f12.7)
          endif

        endif
        enddo
        enddo

      enddo
      enddo

      end

