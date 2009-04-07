      subroutine hklml(p,rsm,e,kmax,nlm,k0,hkl)
C- k,L-dependent smooth Hankel functions.
C ----------------------------------------------------------------------
Ci Inputs
Ci   p     :Function is centered at p
Ci   rsm   :smoothing radius
Ci   e     :energy of smoothed Hankel
Ci   kmax  :polynomial cutoff
Ci   nlm   :L-cutoff for hkl
Ci   k0    :leading dimension of hkl
Co Outputs
Co   hkl   :Bloch-summed smoothed Hankels
Cr Remarks
Cr   H_kL = laplace^k H_L
Cr   Uses the recursion relation H_k+1,L = -e*H_kL - 4*pi*G_kL
Cu Updates
Cu   23 Jan 07 Adapted from hklbl.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer k0,kmax,nlm
      double precision e,rsm,p(3)
      double precision hkl(0:k0,nlm)
C ... Local parameters
      integer nlm0,ilm,k,ll,lmax
      parameter (nlm0=144)
      double precision fpi,pi
      double precision hsm(nlm0),hsmp(nlm0),gklsav,gklnew

      if (nlm .eq. 0) return
      if (nlm .gt. nlm0) call rx('hklml: increase nlm0')
      pi = 4d0*datan(1d0)
      fpi = 4d0*pi
      lmax = ll(nlm)

      call hsmml(p,rsm,e,lmax,hsm,hsmp)
      call gklml(p,rsm,e,kmax-1,nlm,k0,hkl)

C --- H_kL by recursion ---
      do  2  ilm = 1, nlm
      gklsav = hkl(0,ilm)
      hkl(0,ilm) = hsm(ilm)
      do  2  k = 1, kmax
        gklnew = hkl(k,ilm)
        hkl(k,ilm) = -e*hkl(k-1,ilm) - fpi*gklsav
        gklsav = gklnew
    2 continue

      end
C
C      subroutine hhiml(mode,p1,p2,rsm1,rsm2,e1,e2,nlm1,nlm2,kmax,
C     .  ndim1,ndim2,cg,indxcg,jcg,cy,slat,s)
CC- Integrals between smooth Hankels with k-th power of Laplace operator.
CC ----------------------------------------------------------------------
CCi Inputs
CCi   mode  :1's digit (not implemented: always vectors)
CCi         :0 rsm1,rsm2,e1,e2 are scalars
CCi         :1 rsm1,rsm2,e1,e2 are l-dependent vectors
CCi         :10's digit
CCi         :1: do not calculate strux for any (l1,l2) pairs for
CCi         :   if rsm(l1) or rsm(l2) is zero.
CCi   p1    :first center
CCi   p2    :second center
CCi   rsm1  :smoothing radii of Hankels at p1 (l-dependent)
CCi   rsm2  :smoothing radii of Hankels at p2 (l-dependent)
CCi   e1    :energies of smooth Hankels at p1 (l-dependent)
CCi   e2    :energies of smooth Hankels at p2 (l-dependent)
CCi   nlm1  :L-cutoff for functions at p1
CCi   nlm2  :L-cutoff for functions at p2
CCi   kmax  :cutoff in power of Laplace operator
CCi   ndim1 :leading dimension of s
CCi   ndim2 :second dimension of s
CCi   cg    :Clebsch Gordon coefficients (scg.f)
CCi   indxcg:index for Clebsch Gordon coefficients
CCi   jcg   :L q.n. for the C.G. coefficients (scg.f)
CCi   cy    :Normalization constants for spherical harmonics
CCi   slat  :struct containing information about the lattice
CCo Outputs
CCo   s     :integrals; see Remarks
CCr Remarks
CCr   s(L,M) contains integral of H_L^*(r-p1) (laplace)^k H_M(r-p2)
CCr   Row L corresponds to p1 and col M corresponds to p2.
CCr   Strux s(L,M) is computed for dr=p1-p2
CCu Updates
CCu   26 Mar 03 Adapted from hhiml, brute-force, for checking
CCu      hsmml needs analytic derivative; no complex arithmetic
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer mode,jcg(1),indxcg(1),nlm1,nlm2,kmax,ndim1,ndim2
C      double precision p1(3),p2(3),cg(1),cy(1),slat(1),
C     .  rsm1(0:*),rsm2(0:*),e1(0:*),e2(0:*)
C      double complex s(ndim1,ndim2,0:kmax)
CC ... Local parameters
C      integer m,lmx1,lmx2,ll,l1,l2,k,jlm,ilm,lm11,lm21,lm12,lm22,l1t,l2t
C      double precision dr(3)
C
C      if (nlm1 .eq. 0 .or. nlm2 .eq. 0) return
C      do  1  m = 1, 3
C    1 dr(m) = p1(m)-p2(m)
C      lmx1 = ll(nlm1)
C      lmx2 = ll(nlm2)
C
C      do  3  k = 0, kmax
C      do  3  jlm = 1, nlm2
C      do  3  ilm = 1, nlm1
C    3 s(ilm,jlm,k) = dcmplx(0d0,0d0)
C
C      l1t = -1
C      do  20  l1 = 0, lmx1
C        if (l1 .le. l1t) goto 20
C        call gtbsl2(l1,lmx1,e1,rsm1,l1t)
CC       l1t = l1
C
C        l2t = -1
C        do  22  l2 = 0, lmx2
C          if (l2 .le. l2t) goto 22
C          call gtbsl2(l2,lmx2,e2,rsm2,l2t)
CC         l2t = l2
C
C          lm11 = l1**2+1
C          lm12 = (l1t+1)**2
C          lm21 = l2**2+1
C          lm22 = (l2t+1)**2
C          if (mode/10.eq.1 .and. rsm1(l1)*rsm2(l2) .eq. 0) goto 22
C          call phhiml(dr,rsm1(l1),rsm2(l2),e1(l1),e2(l2),lm11,lm12,
C     .      lm21,lm22,kmax,ndim1,ndim2,cg,indxcg,jcg,cy,slat,s)
C   22   continue
C   20 continue
C      end
C      subroutine phhiml(dr,rsm1,rsm2,e1,e2,mlm1,nlm1,mlm2,nlm2,kmax,
C     .  ndim1,ndim2,cg,indxcg,jcg,cy,slat,s)
CC- Integrals between smooth Hankels with k-th power of Laplace operator.
CC ----------------------------------------------------------------------
CCi Inputs
CCi   dr    :p1-p2
CCi   rsm1  :smoothing radius of Hankels at p1
CCi   rsm2  :smoothing radius of Hankels at p2
CCi   e1    :energy of smooth Hankels at p1
CCi   e2    :energy of smooth Hankels at p2
CCi   nlm1  :L-cutoff for functions at p1
CCi   nlm2  :L-cutoff for functions at p2
CCi   kmax  :cutoff in power of Laplace operator
CCi   ndim1 :leading dimension of s
CCi   ndim2 :second dimension of s
CCi   cg    :Clebsch Gordon coefficients (scg.f)
CCi   indxcg:index for Clebsch Gordon coefficients
CCi   jcg   :L q.n. for the C.G. coefficients (scg.f)
CCi   cy    :Normalization constants for spherical harmonics
CCi   slat  :struct containing information about the lattice
CCo Outputs
CCo   s     :integrals; see Remarks
CCr Remarks
CCr   s(L,M) contains integral of H_L^*(r-p1) (laplace)^k H_M(r-p2)
CCu Updates
CCu   18 May 00 Adapted from nfp hhi_bl.f
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer jcg(1),indxcg(1),mlm1,nlm1,mlm2,nlm2,kmax,ndim1,ndim2
C      double precision dr(3),cg(1),cy(1),slat(1),rsm1,
C     .  rsm2,e1,e2
C      double complex s(ndim1,ndim2,0:kmax)
CC ... Local parameters
C      integer nlm0,ktop0,icg,icg1,icg2,ii,ilm,ilm1,ilm2,indx,ip,k,
C     .  ktop,l1,l2,ll,lm,lmax1,lmax2,lmaxx,nlmx
C      parameter( nlm0=100, ktop0=8 )
C      double precision fpi,e,fac,fac1,fac2,gam1,gam2,gamx,rsmx,alat
C      double complex hkl1(0:ktop0,nlm0),hkl2(0:ktop0,nlm0),
C     .  hsm(nlm0),hsmp(nlm0)
C
C      fpi = 16d0*datan(1.d0)
C      gam1 = 0.25d0*rsm1*rsm1
C      gam2 = 0.25d0*rsm2*rsm2
C      gamx = gam1+gam2
C      rsmx = 2d0*dsqrt(gamx)
C      lmax1 = ll(nlm1)
C      lmax2 = ll(nlm2)
C      lmaxx = lmax1+lmax2
C      nlmx = (lmaxx+1)**2
C      ktop = max0(lmax1,lmax2)+kmax
C      if (nlmx .gt. nlm0)  call rx('hhiml: increase nlm0')
C      if (ktop .gt. ktop0) call rx('hhiml: increase ktop0')
C
CC ... Set up functions for connecting vector p2-p1
C      if (dabs(e1-e2) .gt. 1d-5) then
C         fac1 = dexp(gam2*(e2-e1))/(e1-e2)
C         fac2 = dexp(gam1*(e1-e2))/(e2-e1)
C         call hklml(dr,rsmx,e1,ktop,nlmx,ktop0, hkl1)
C         call hklml(dr,rsmx,e2,ktop,nlmx,ktop0, hkl2)
C         do  4  ilm = 1, nlmx
C           do  5  k = 0, ktop
C             hkl1(k,ilm) = fac1*hkl1(k,ilm) + fac2*hkl2(k,ilm)
C    5      continue
C    4    continue
C      else
C         call upack('lat alat',slat,alat,0,0,0,0)
C         e = .5d0*(e1+e2)
C         call hklml(dr,rsmx,e,ktop,nlmx,ktop0, hkl2)
C         call hsmml(alat*dr,rsmx,e,lmaxx,hsm,hsmp)
C         do  6  ilm = 1, nlmx
C           hkl1(0,ilm) = hsmp(ilm) - gamx*hsm(ilm)
C           do  7  k = 1, ktop
C             hkl1(k,ilm) = -e*hkl1(k-1,ilm) - hkl2(k-1,ilm)
C    7      continue
C    6    continue
C
C      endif
C
CC ... Combine with Clebsch-Gordan coefficients
C      do  11  ilm1 = mlm1, nlm1
C      l1 = ll(ilm1)
C      do  11  ilm2 = mlm2, nlm2
C      l2 = ll(ilm2)
C      ii = max0(ilm1,ilm2)
C      indx = (ii*(ii-1))/2 + min0(ilm1,ilm2)
C      icg1 = indxcg(indx)
C      icg2 = indxcg(indx+1)-1
C      do  11  icg = icg1, icg2
C      ilm = jcg(icg)
C      lm = ll(ilm)
C      k = (l1+l2-lm)/2
C      fac = fpi*(-1d0)**l1*cg(icg)
C      do  12  ip = 0, kmax
C   12 s(ilm1,ilm2,ip) = s(ilm1,ilm2,ip) + fac*hkl1(k+ip,ilm)
C   11 continue
C
C      end
C
C      subroutine smhsml(mode,ssite,sspec,slat,vavg,ndimh,iprmb,h,s)
CC- Smoothed Bloch Hamiltonian (constant potential) and overlap matrix
CC ----------------------------------------------------------------------
CCi Inputs
CCi   mode  :0 compute both hamiltonian and overlap
CCi         :  otherwise, compute overlap only.  
CCi         :  In this case, vavg is not used
CCi   ssite :struct for site-specific information; see routine usite
CCi     Elts read: spec pos
CCi     Stored:    *
CCi     Passed to: *
CCi   sspec :struct for species-specific information; see routine uspec
CCi     Elts read: *
CCi     Stored:    *
CCi     Passed to: uspecb
CCi   slat  :struct for lattice information; see routine ulat
CCi     Elts read: ocg ojcg oidxcg ocy
CCi     Stored:    *
CCi     Passed to: hhibl
CCi   vavg  :constant potential (MT zero) to be added to h
CCi   ndimh :dimension of hamiltonian
CCi   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
CCo Outputs
CCo   h     :smooth Bloch hamiltonian added to h (= s * vavg)
CCo   s     :smooth Bloch overlap added to s
CCr Remarks
CCr  *NB: The matrix subblock in s for sites (i,j) is computed with
CCr       the connecting vector pi-pj.  This convention is different 
CCr       from the tight-binding convention, where s is computed for
CCr       pj-pi.  Thus this (i,j) matrix subblock is the hermitian 
CCr       conjugate of the same subblock in the tight-binding case.
CCr       tight-binding case.  Note, however, that the entire s
CCr       would not correspond to the hermitian conjugate of the tight
CCr       binding case.
CCr  *MPI
CCr   See remarks to hsibl. Buffers for h and s are F90 ALLOCATEd as
CCr   they need to be used locally as dimensioned. A buffer is taken
CCr   from the heap for ALLREDUCE.
CCu Updates
CCu   14 Aug 02 Added overlap-only option
CCu   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
CCu   15 Feb 02 (ATP) Added MPI parallelization
CCu   27 Aug 01 Extended to local orbitals.
CCu   02 Mar 01 Bug fix for multiple-kappa case
CCu   19 May 00 Adapted from nfp smhs_q.f
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
CC#ifdefC MPI
CC#ifdefC MPE
CC      include "mpef.h"
CC#endifC 
CC      include "mpif.h"
CC      integer numprocs, ierr, status(MPI_STATUS_SIZE)
CC      integer MAX_PROCS
CC      parameter (MAX_PROCS = 100)
CC      integer resultlen
CC      character*(MPI_MAX_PROCESSOR_NAME) name
CC      character*10 shortname(0:MAX_PROCS-1)
CC      character*20 ext
CC      character*26 datim
CC      integer namelen(0:MAX_PROCS-1)
CC      double precision starttime, endtime
CC      logical mlog,cmdopt
CC      integer i,lgunit
CC      character*120 strn
CC#endif
C      integer procid,master
C      integer mode,ndimh,iprmb(1)
C      double precision ssite(1),sspec(1),slat(1),vavg
C      double complex h(ndimh,ndimh),s(ndimh,ndimh)
CC ... Local parameters
C      integer nlms,kdim,n0,nkap0,nkape
C      parameter (nlms=25, kdim=1, n0=10, nkap0=3, nkape=2)
C      integer nbas,nglob,ocg,ocy,oidxcg,ojcg
C      integer i1,i2,ib1,ib2,ilm1,ilm2,io1,io2,is1,is2,nlm1,nlm2
C      integer lh1(nkap0),lh2(nkap0),nkap1,nkap2
C      integer ltab1(n0*2),ktab1(n0*2),offl1(n0*2),norb1,ik1
C      integer ltab2(n0*2),ktab2(n0*2),offl2(n0*2),norb2,ik2,l1,l2
C      double precision e1(n0,nkap0),rsm1(n0,nkap0),p1(3),
C     .                 e2(n0,nkap0),rsm2(n0,nkap0),p2(3),xx
C      double complex s0(nlms,nlms,0:kdim,nkape,nkape)
CC#ifdefC MPI
CC      double complex, dimension(:,:), allocatable :: hbuf, sbuf
CC      integer, dimension(:,:), allocatable :: index
CC      integer iloop,obuff
CC#endif
CC#ifdefC MPE
CCC Event numbers:
CC      include "events.ins"
CC#endif 
CC ... Heap
C      integer w(1)
C      common /w/ w
C
C      call tcn('smhsbl')
C
CC#ifdefC MPI
CC      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
CC      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
CC      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
CC      call strcop(shortname(procid),name,10,'.',i)
CC      namelen(procid) = i-1
CC      master = 0
CC      mlog = cmdopt('--mlog',6,0,strn)
CC#else
C      procid = 0
C      master = 0
CC#endif
C
CC --- Setup ---
C      nbas  = nglob('nbas')
C
C      call upack('lat ocg ojcg oidxcg ocy',slat,ocg,ojcg,oidxcg,ocy,0)
CC     call u_lat_vecs(slat,alat,plat,qlat)
C
CC --- Loop over first and second site indices ---
CC#ifdefC MPI
CC#ifdefC MPE
CC      ierr = MPE_LOG_EVENT(EVENT_START_SMHSBL,procid,"smhsbl")
CC#endifC 
CC      allocate (hbuf(ndimh,ndimh), stat=ierr)
CC      allocate (sbuf(ndimh,ndimh), stat=ierr)
CC      call dcopy(2*ndimh*ndimh,0d0,0,hbuf,1)
CC      call dcopy(2*ndimh*ndimh,0d0,0,sbuf,1)
CC      allocate (index(0:numprocs-1,0:nbas-1), stat=ierr)
CC      call dstrbp(nbas,numprocs,-1,index(0,0))
CC      do  iloop = 1, index(procid,0)
CC        ib1 = index(procid,iloop)
CC        if (mlog) then
CC          call gettime(datim)
CC          call awrit4(' smhsbl '//datim//' Process %i of %i on '
CC     .        //shortname(procid)(1:namelen(procid))//
CC     .        ' starting atom %i of %i',' ',256,lgunit(3),
CC     .        procid,numprocs,ib1,index(procid,0))
CC        endif
CC#else
C      do  ib1 = 1, nbas
CC#endif 
C        call upack('site spec pos',ssite,ib1,is1,p1,0,0)
C        call uspecb(0,1,sspec,is1,is1,lh1,rsm1,e1,nkap1)
C        call orbl(ib1,0,ndimh,iprmb,norb1,ltab1,ktab1,xx,offl1,xx)
C        stop 'FIX local orb'
CC       for now, until iprmb is rewritten
C        do  io1 = norb1, 1, -1
C          if (ktab1(io1) .gt. nkap1) norb1=norb1-1
C        enddo
C        do  ib2 = ib1, nbas
C          call upack('site spec pos',ssite,ib2,is2,p2,0,0)
C          call uspecb(0,1,sspec,is2,is2,lh2,rsm2,e2,nkap2)
C          call orbl(ib2,0,ndimh,iprmb,norb2,ltab2,ktab2,xx,offl2,xx)
CC         for now, until iprmb is rewritten
C          do  io2 = norb2, 1, -1
C            if (ktab2(io2) .gt. nkap2) norb2=norb2-1
C          enddo
C
CC     ... Integrals between all orbitals connecting ib1 and ib2
C          do  i1 = 1, nkap1
C            do  i2 = 1, nkap2
C              nlm1 = (lh1(i1)+1)**2
C              nlm2 = (lh2(i2)+1)**2
C              if (nlm1 .gt. nlms .or. nlm2 .gt. nlms)
C     .          call rx('smhsbl: increase nlms')
C              call hhiml(11,p1,p2,rsm1(1,i1),rsm2(1,i2),e1(1,i1),
C     .          e2(1,i2),nlm1,nlm2,1,nlms,nlms,w(ocg),w(oidxcg),w(ojcg),
C     .          w(ocy),slat,s0(1,1,0,i1,i2))
C            enddo
C          enddo
C
CC     ... Loop over orbital indices, poke block of integrals into s,h
C          do  io2 = 1, norb2
CC           l2,ik2 = l and kaph indices, needed to locate block in s0
C            l2  = ltab2(io2)
C            ik2 = ktab2(io2)
CC           i2 = orbital index in iprmb order
C            i2 = offl2(io2)
C          
C            do  ilm2 = l2**2+1, (l2+1)**2
C              i2 = i2+1
C
C              if (mode .eq. 0) then
C              do  io1 = 1, norb1
CC               l1,ik1 = l and kaph indices, needed to locate block in s0
C                l1  = ltab1(io1)
C                ik1 = ktab1(io1)
CC               i1 = orbital index in iprmb order
C                i1 = offl1(io1)
C                do  ilm1 = l1**2+1, (l1+1)**2
C                  i1 = i1+1
CC#ifdefC MPI
CC                  sbuf(i1,i2) = sbuf(i1,i2) + s0(ilm1,ilm2,0,ik1,ik2)
CC                  hbuf(i1,i2) = hbuf(i1,i2) - s0(ilm1,ilm2,1,ik1,ik2)
CC     .              + vavg*s0(ilm1,ilm2,0,ik1,ik2)
CC#else
C                  s(i1,i2) = s(i1,i2) + s0(ilm1,ilm2,0,ik1,ik2)
C                  h(i1,i2) = h(i1,i2) - s0(ilm1,ilm2,1,ik1,ik2)
C     .              + vavg*s0(ilm1,ilm2,0,ik1,ik2)
CC#endif 
C                enddo
C              enddo
C              else
C              do  io1 = 1, norb1
CC               l1,ik1 = l and kaph indices, needed to locate block in s0
C                l1  = ltab1(io1)
C                ik1 = ktab1(io1)
CC               i1 = orbital index in iprmb order
C                i1 = offl1(io1)
C                do  ilm1 = l1**2+1, (l1+1)**2
C                  i1 = i1+1
CC#ifdefC MPI
CC                  sbuf(i1,i2) = sbuf(i1,i2) + s0(ilm1,ilm2,0,ik1,ik2)
CC#else
C                  s(i1,i2) = s(i1,i2) + s0(ilm1,ilm2,0,ik1,ik2)
CC#endif 
C                enddo
C              enddo
C              endif
C
C            enddo
C          enddo
C        enddo
C      enddo
CC#ifdefC MPI
CC#ifdefC MPE
CC      ierr = MPE_LOG_EVENT(EVENT_END_SMHSBL,procid,"smhsbl")
CC      ierr = MPE_LOG_EVENT(EVENT_START_BARRIER,procid,"barrier")
CC#endifC 
CC      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
CC#ifdefC MPE
CC      ierr = MPE_LOG_EVENT(EVENT_END_BARRIER,procid,"barrier")
CC      ierr = MPE_LOG_EVENT(EVENT_START_ALLRED,procid,"allreduce")
CC#endifC 
CC      if (mode .eq. 0) then
CC      call defcc(obuff, ndimh*ndimh)
CC      call MPI_ALLREDUCE(hbuf,w(obuff),2*ndimh*ndimh,
CC     .     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
CC      if (mlog) then
CC        call gettime(datim)
CC        call awrit3(' smhsbl '//datim//' Process %i of %i on '
CC     .        //shortname(procid)(1:namelen(procid))//
CC     .        ' allreduce h ndimh=%i',' ',256,lgunit(3),
CC     .        procid,numprocs,ndimh)
CC      endif
CC      call daxpy(2*ndimh*ndimh,1d0,w(obuff),1,h,1)
CC      endif
CC      call MPI_ALLREDUCE(sbuf,w(obuff),2*ndimh*ndimh,
CC     .     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
CC      if (mlog) then
CC        call gettime(datim)
CC        call awrit3(' smhsbl '//datim//' Process %i of %i on '
CC     .        //shortname(procid)(1:namelen(procid))//
CC     .        ' allreduce s ndimh=%i',' ',256,lgunit(3),
CC     .        procid,numprocs,ndimh)
CC      endif
CC      call daxpy(2*ndimh*ndimh,1d0,w(obuff),1,s,1)
CC      call rlse(obuff)
CC      deallocate(index, stat=ierr)
CC      deallocate(hbuf, stat=ierr)
CC      deallocate(sbuf, stat=ierr)
CC#ifdefC MPE
CC      ierr = MPE_LOG_EVENT(EVENT_END_ALLRED,procid,"allreduce")
CC#endifC 
CC#endif 
C
CC ... Occupy second half of matrix
C      if (mode .eq. 0) then
C      do  i1 = 1, ndimh
C        do  i2 = i1, ndimh
C          h(i2,i1) = dconjg(h(i1,i2))
C          s(i2,i1) = dconjg(s(i1,i2))
C        enddo
C      enddo
C      else
C      do  i1 = 1, ndimh
C        do  i2 = i1, ndimh
C          s(i2,i1) = dconjg(s(i1,i2))
C        enddo
C      enddo
C      endif
C      call tcx('smhsbl')
C
CC      call zprm('h-warp',2,h,ndimh,ndimh,ndimh)
CC      call zprm('s-warp',2,s,ndimh,ndimh,ndimh)
C      end

