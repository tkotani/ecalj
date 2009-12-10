      subroutine pairg(nbas,isw,ips,alat,plat,pos,
     .  range,ntab,iax,ntabg,nttabg)
C- Prepare an auxiliary neighbor table out of iax for a smaller range
C ----------------------------------------------------------------
Ci Inputs:
Ci   nbas  :the size of the basis
Ci          (no extensions to the basis such as padding are allowed)
Ci   isw:  :specifies how range is calculated.
Ci           0: vector length must be < range(i)+range(j)
Ci           1: include all connecting vecs with r < range(i)
Ci   ips   :index to which species each site belongs. NB:
Ci          ips(1) = -1 => rmax is a global scalar, independent of site
Ci          ips(1) =  0 => rmax is site-, not species-dependent.
Ci          In either of these cases, ips is not used.
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat 
Ci   pos   :site positions 
Ci   range :maximum range for connecting vector, in a.u.
Ci          This quantity may be a scalar, a species-dependent
Ci          array, or a site-dependent array, depending on ips(1);
Ci          see description of ips.  Precisely what meaning range has
Ci          depends on isw.
Ci   ntab  :ntab(ib) no. pairs in neighbor table preceding ib (pairc.f)
Co   iax   :neighbor table containing information about each pair (pairc.f)
Co   iax(1)   :site index to basis atoms ib=source;
Co             all pairs with common ib are contiguous
Co   iax(2)   :site index to jb=field of each pair
Co   iax(3..5):multiples of plat added the difference in site positions
Co             that connect the pair.
Co   iax(6)   :index to conjugate (jb,ib) pair matching (ib,jb)
Co             NB: no matching pairs outside (ib1..ib2) can be found.
Co   iax(7)   :permutation index ordering cluster by increasing
Co             effective site index; see ppair4.f
Co   iax(8)   :not used
Co   iax(9)   :not used 
Co   iax(10)  :effective site index; see siteid.f
Co Outputs:
Co   ntabg :ntabg(ib) no. pairs in cluster ib in the auxiliary table 
Co   nttabg:total no. pairs for all clusters, nttabg = \sum_ib=1^nbas ntabg(ib) 
Cr Remarks
Cr   It is assumed that iax and ntab have been generated previously with some
Cr   range0 >= range. Program does not create another neighbor table but uses
Cr   array ntabg as a pointer to clip the respective subarray out of iax:
Cr   pairs with r < range in cluster ib are those with index
Cr   (ntab(ib)+1...ntab(ib)+ntabg(ib)) in the iax table. 
Cu Updates
Cu   12 Dec 07 (S. Lozovoi) Adapted from pairc.f
C ----------------------------------------------------------------
C     implicit none
      integer nbas,isw,nttabg,niax 
      parameter (niax=10)
      integer iax(niax,*),ntab(nbas+1),ips(nbas),ntabg(nbas)
      double precision alat,plat(3,3),pos(3,nbas),range(*)
C Local variables
      integer ib,is,jb,moder,np,it,iclus,npm
      double precision r1,r2,rr,rcut,rcut2,drr2,dum(3)

C Heap
      integer w(1)
      common /w/ w

C --- Setup ---
      moder = mod(isw,10)
      npm = 0

      r1 = range(1)
      r2 = r1
      nttabg = 0
      do ib = 1, nbas
        np = 0

        if (ips(1) .ge. 0) then
          is = ib
          if (ips(1) .gt. 0) is = ips(ib)
          r1 = range(is)
        endif
        rcut = r1
        iclus = ntab(ib)+1
    
        do it = iclus, ntab(ib+1)
c         ib0 = iax(1,it)
          if (iax(1,it) .ne. ib) 
     .      call rx('pairg: mismatch in iax table')
          jb = iax(2,it)
          if (moder .eq. 0) then
            if (ips(1) .ge. 0) then
              is = ib
              if (ips(1) .gt. 0) is = ips(ib)
              r2 = range(is)
            endif 
            rcut = r1+r2
          endif
          rcut2 = (rcut/alat)**2

          rr = drr2(plat,pos(1,ib),pos(1,jb),
     .    iax(3,it),iax(4,it),iax(5,it),dum)

          if (rr .le. rcut2) then
            np = np + 1
          else
            exit
          endif
        enddo

        ntabg(ib) = np 
        npm = max(np,npm)
        nttabg = nttabg + np 
c       lmtabg = lmtabg + isum(np,iax(9,iclus),niax) 
      enddo

      call info(20,1,0,' pairg: %i pairs total, %i is max cluster size',
     .  nttabg,npm)

c     call info(20,1,0,' pairg: %i pairs total, %i is max cluster size',
c    .  nttabg,maxval(ntabg))
c     call info(20,0,0,' %i channels altogether',lmtabg)

      end

