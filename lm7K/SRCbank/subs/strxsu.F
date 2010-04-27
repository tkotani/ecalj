      subroutine strxsu(nlma,nxi,lxi,n0,loka,nbas,ips,cg,jcg,indxcg,
     .  nlmbx,nlmp,npow,ocf,oip,oikl,ojkl)
C- Setup for structure constant tables at one site
C ----------------------------------------------------------------------
Ci Inputs
Ci   nlma     augmentation to nlm around this site
Ci   ips      species table.  ips(1) = 0 => ips doesn't exist;
Ci            assume nxi and/or lxi are site-based
Ci   nxi,lxi  number of energies for site ib, and l-cutoff
Ci            use nxi(1)=-n for a fixed number of energies (n).
Ci            use nxi(1)=-99 for a single global energy and lxi
Ci   n0       leading dimension of lxi
Ci   loka     1, Scale cf to conform to to Andersen conventions:
Ci            scale by factor - (2l-1)!! (2l''-1)!! / 2
Co Outputs
Co   nlmbx    max(lxi+1)**2 for all sites
Co   nlmp     max(lxi+lmxa+1)**2 for all sites
Co   npow     min(lmxa,max(lxi))
Co   ocf,oip,oikl,ojkl are allocated and initialized
Cr Remarks
Cr   Use this routine in conjunction with hstrux, hstrud, etc
C ----------------------------------------------------------------------
C     implicit none
      integer n0,loka,nbas,nlma,nlmbx,nlmp,npow,ocf,oikl,oip
      integer nxi(1),ips(1),lxi(n0,1),jcg(1),indxcg(1),ojkl(0:20)
      double precision cg(1)
C Heap
      integer w(1)
      common /w/ w
C Local
      integer i,jb,je,js,klmb,lb,lbx,ll,lmaxa,lx,ncfx,n

C --- Get maximum lb, related dimensions ---
      lmaxa = ll(nlma)
      do  1  i = 0, 15
    1 ojkl(i) = -99
      lbx = -1
      do  7  jb = 1, nbas
        if (ips(1) .le. 0) then
          js = jb
        else
          js = ips(jb)
        endif
        if (nxi(1) .eq. -99) then
          lx = max(lxi(1,1),0)
          ojkl(lx) = -98
          lbx = max0(lbx,lx)
        else
          if (nxi(1) .lt. 0) then
            n = -nxi(1)
          else
            n = nxi(js)
          endif
          do  8  je = 1, n
            lx = max(lxi(je,js),0)
            ojkl(lx) = -98
            lbx = max0(lbx,lx)
    8     continue
        endif
    7 continue
      nlmbx = (lbx+1)**2
      nlmp = (lbx+lmaxa+1)**2
      npow = min0(lmaxa,lbx)

C --- Make cg-tables for maximal lb-value ---
      ncfx = 15000
      call defrr(ocf,        ncfx)
      call defrr(oip,        ncfx)
      call defrr(oikl,       nlmp*(npow+1))
      call defrr(ojkl(lbx),  nlmp*(npow+1))

      call nstru0(nlma,nlmbx,nlmp,npow,loka,cg,jcg,indxcg,
     .  w(oikl),w(ojkl(lbx)),ncfx,w(oip),w(ocf))

C --- Make jkl for needed lower lb-values ---
      do  10  lb = 0, lbx-1
        if (ojkl(lb) .eq. -98) then
          klmb = (lb+1)**2
          call defrr(ojkl(lb),   nlmp*(npow+1))
          call nstru1(nlma,nlmp,npow,klmb,jcg,indxcg,w(oikl),
     .      w(ojkl(lb)))
        endif
   10 continue

      end
      subroutine nstru0(nlma,nlmb,nlmp,npow,loka,cg,jcg,indxcg,
     .  ikl,jkl,ncfx,ip,cf)
C- Set up Clebsch-Gordans in a way which permits longer vectors.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nlma     augmentation to nlma around this site
Ci   nlmb     augmentation of functions of size nlmb
Ci   nlmp     (lmxb+lmxa+1)**2
Ci   npow     min(lmxa,lmxb)
Ci   ncfx     max no. of coffs to be calculated
Ci   loka     1, Scale cf to conform to to Andersen conventions:
Ci            scale by factor - (2l-1)!! (2l''-1)!! / 2
Co Outputs
Co   ikl,jkl  point to first,last element for a (ilm,ipow) pair.
Co   cf       coeff*(-1)**lk
Co   ip       position in strux as a linear array.
C ----------------------------------------------------------------------
C     implicit none
      integer nlmp,npow,ncfx,nlmb,nlma,loka
      integer jcg(1),indxcg(1),ikl(nlmp,0:npow),jkl(nlmp,0:npow),ip(1)
      double precision cg(1),cf(ncfx)
      integer icg,icg1,icg2,ii,ilm,indx,ipow,iprint,jj,klm,lk,ll,lm,
     .  mlm,ntot,l,lmax,lmxx
      parameter (lmxx=20)
      double precision fac2l(0:lmxx)

      do  1  ipow = 0, npow
      do  1  ilm = 1, nlmp
    1 jkl(ilm,ipow) = 0
C --- Loop over cg coeffs to count non-zero elements
      do  10  klm = 1, nlmb
      lk = ll(klm)
      do  10  mlm = 1, nlma
      lm = ll(mlm)
      ii = max0(mlm,klm)
      indx = (ii*(ii-1))/2+min0(mlm,klm)
      icg1 = indxcg(indx)
      icg2 = indxcg(indx+1)-1
      do  10  icg = icg1, icg2
      ilm = jcg(icg)
      ipow = (lm+lk-ll(ilm))/2
   10 if(ilm.le.nlmp .and. ipow.le.npow) jkl(ilm,ipow) = jkl(ilm,ipow)+1
C --- Set up start pointers ikl ---------
      ntot = 0
      do  14  ipow = 0, npow
      do  14  ilm = 1, nlmp
      ikl(ilm,ipow) = ntot+1
      ntot = ntot+jkl(ilm,ipow)
   14 jkl(ilm,ipow) = ikl(ilm,ipow)-1
      if (iprint() .ge. 80) write(6,716) nlma,nlmb,nlmp,npow,ntot
  716 format(' nstru0:  nlma,nlmb,nlmp,npow',4i5,'    ntot=',i8)
      if (ntot .gt. ncfx) call rxi('nstru0: increase ncfx to',ntot)

      if (loka .eq. 1) then
        lmax = ll(nlma) + ll(nlmb)
        if (lmax .gt. lmxx) call rxi('nstru0: increase lmxx to',lmax)
        fac2l(0) = 1
        do  5  l = 1, lmax
    5   fac2l(l) = fac2l(l-1) * (2*l-1)
      endif

C --- Loop over cg coeffs and store coeffs and pointers ---
      jj=0
      do  20  klm = 1, nlmb
      lk=ll(klm)
      do  20  mlm = 1, nlma
        lm=ll(mlm)
        jj=jj+1
        ii=max0(mlm,klm)
        indx=(ii*(ii-1))/2+min0(mlm,klm)
        icg1=indxcg(indx)
        icg2=indxcg(indx+1)-1
        do  20  icg = icg1, icg2
        ilm=jcg(icg)
        ipow=(lm+lk-ll(ilm))/2
        if (ilm.le.nlmp .and. ipow.le.npow) then
          jkl(ilm,ipow)=jkl(ilm,ipow)+1
          ii=jkl(ilm,ipow)
          ip(ii)=jj
          cf(ii)=cg(icg)*(-1d0)**lk
          if (loka .eq. 1) cf(ii) = cf(ii) * 2/(fac2l(lm)*fac2l(lk))
        endif
   20 continue
      end

      subroutine nstru1(nlma,nlmp,npow,klmb,jcg,indxcg,ikl,jkl)
C- Call after nstru0, sets up pointers jkl for klmb.le.nlmb
C ----------------------------------------------------------------------
Ci Inputs
Ci   nlma     augmentation to nlma around this site
Ci   nlmb     augmentation of functions of size nlmb
Ci   nlmp     (lmxb+lmxa+1)**2
Ci   npow     min(lmxa,lmxb)
ci   klmb     (lmxb+1)**2
Ci   ncfx     max no. of coffs to be calculated
Co Outputs
Co   ikl,jkl  point to first,last element for a (ilm,ipow) pair.
C ----------------------------------------------------------------------
C     implicit none
      integer nlmp,npow,klmb,nlma
      integer jcg(1),indxcg(1),ikl(nlmp,0:npow),jkl(nlmp,0:npow)
C Local
      integer icg,icg1,icg2,ii,ilm,indx,ipow,klm,lk,ll,lm,mlm

      do  14  ipow = 0, npow
      do  14  ilm = 1, nlmp
   14 jkl(ilm,ipow) = ikl(ilm,ipow)-1
C --- Loop over cg coeffs: count elements up to klmb ---
      do  20  mlm = 1, nlma
      lm = ll(mlm)
      do  20  klm = 1, klmb
      lk = ll(klm)
      ii = max0(mlm,klm)
      indx = (ii*(ii-1))/2+min0(mlm,klm)
      icg1 = indxcg(indx)
      icg2 = indxcg(indx+1)-1
      do  20  icg = icg1, icg2
      ilm = jcg(icg)
      ipow = (lm+lk-ll(ilm))/2
      if (ilm.le.nlmp .and. ipow.le.npow) jkl(ilm,ipow)=jkl(ilm,ipow)+1
   20 continue

C --- Printout ---
C      do  56  ilm = 1, nlmp
C      print 230, ilm,(ikl(ilm,ipow),ipow=1,npow)
C   56 print 231, (jkl(ilm,ipow),ipow=1,npow)
C  230 format(i5,'  ikl',8i5)
C  231 format(5x,'  jkl',8i5)
C      pause
      end

