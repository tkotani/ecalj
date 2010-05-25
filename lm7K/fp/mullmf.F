#define F90 1
      subroutine mullmf(nbas,ssite,sspec,iprmb,z,n,nspc,iq,isp,mode,
     .nsites,lsites,lmxch,nchan,lchan,lmdim,nddos,doswt)

      use m_struc_def  !Cgetarg

C- Make Mulliken decomposition of the norm at this k-point
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   ssite :struct containing site-specific information
Ci     Elts read: spec
Ci   sspec :struct containing species-specific information
Ci     Elts read: lmxb
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   z,n   :eigenvectors and dimension n
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   mode  :0 all sites atom-resolved
Ci         :1 all sites l-resolved
Ci         :2 all sites lm-resolved
Ci         :3 site list atom-resolved
Ci         :4 site list l-resolved
Ci         :5 site list lm-resolved
Ci   mode,nsites,lsites,lmxch,nchan,lchan (see sumlst.f)
Ci   nddos :dimensions doswt
Co Outputs:
Co   lmdim :leading dimension of lchan
Co   doswt :DOS weights for writing to moms file (lmdos.f)
Cr Remarks
Cr   Orthogonal basis : D_in = (z_in)+ z_in where
Cr     i = orbital index and n = band index
Cr   Nonorthogonal basis : D_in = (z~_in)+ z_in
Cr     Here z~ is contravariant form of z.
Cr     Overlap matrix is S = (z z+)^-1
Cr     z~+ = z+ S = z+ (z+)^-1 z^-1 = z^-1
Cu Updates
Cu   08 Jul 08 Dimension dowst separately from z
Cu   07 Jun 06 Extended to noncollinear case.  Altered argument list
Cu   20 Mar 01 Written by ATP
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,nspc,iq,isp,n,iprmb(1),mode,nsites,
     .lsites(nsites),nchan,lchan(nchan),lmxch,lmdim,nddos
      real(8):: doswt(nchan,nddos*nspc,nspc)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)

      double complex z(n,nspc,n*nspc)
C ... Local parameters
      logical lmp,lp,atp
      integer isite,ib,is,ichan,lmxb,iband,iprint,ipr,iprmin,lgunit,
     .ispc,ksp,stdo,nev,nx
      integer n0,nkap0
      parameter (n0=10,nkap0=3)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)
      integer blks(n0*nkap0),ntab(n0*nkap0)
      integer io,ikap,l,nlm1,nlm2,ilm,jlm,i,j,norb
      double precision xx
#if F90
      integer ipiv(n*nspc)
      complex(8),allocatable:: zt(:,:,:),work(:,:,:)
      allocate(zt(n*nspc,n,nspc),work(n*nspc,n,nspc))
#elif AUTO_ARRAY
      integer ipiv(n*nspc)
      double complex zt(n*nspc,n,nspc),work(n*nspc,n,nspc)
#else
      integer ipiv(1)
      double complex zt(1,1,1),work(1,1,1)
      call rx('mullmf only implemented for automatic arrays')
#endif

C ... Setup
      stdo = lgunit(1)
      call tcn ('mullmf')
C     nx = total hamiltonian dimension
      nx = n*nspc
C     Number of eigenvectors must equal hamiltonian dimemnsion
      nev = nx

C ... Form contravariant (z~)+ = z^-1
      call zcopy(nx**2,z,1,zt,1)
      call zgetrf(nx,nx,zt,nx,ipiv,j)
      if (j .ne. 0) call rx('mullmf: failed to generate overlap')
      call zgetri(nx,zt,nx,ipiv,work,nx**2,j)

C      call zprm('z',2,z,nx,nx,nx)
C      call zprm('zt',2,zt,nx,nx,nx)

      call isanrg(mode,0,5,' mullmf:','mode',.true.)
      iprmin = 80*iq*isp
      ipr = iprint()
      if (ipr .ge. iprmin) write(stdo,1)
    1 format (' mullmf:  site spec  lmax  norb  l  nlm1   nlm2  ikappa',
     .'  offset ilm ichan')
      call dpzero(doswt,nchan*nddos*nspc)
      atp = .false.
      lp  = .false.
      lmp = .false.
C --- Decompose by atom:
      if (mode .eq. 0 .or. mode .eq. 3) then
        lmdim = 1
        atp = .true.
C --- Decompose by l:
      elseif (mode .eq. 1 .or. mode .eq. 4) then
        lmdim = lmxch + 1
        lp  = .true.
C --- Decompose by l:
      elseif (mode .eq. 2 .or. mode .eq. 5) then
        lmdim = (lmxch + 1)**2
        lmp = .true.
      endif
      if (mode .lt. 3) nsites = nbas

C --- Loop over channels ---
      do  isite = 1, nsites
        ib = lsites(isite)

        is=ssite(ib)%spec


        lmxb=sspec(is)%lmxb

C   ... Loop over all orbitals centered at this site
        call orbl(ib,0,n,iprmb,norb,ltab,ktab,xx,offl,xx)
C   ... Block into groups of consecutive l
        call gtbsl1(4,norb,ltab,ktab,xx,xx,ntab,blks)
        if (ipr .ge. iprmin) write(stdo,2) ib,is,lmxb,norb
    2   format (9x,i4,1x,i3,6x,i1,3x,i2)

C       In the noncollinear case, isp=1 always => need internal ispc=1..2
C       ksp is the current spin index in both cases:
C       ksp = isp  in the collinear case
C           = ispc in the noncollinear case
C       whereas ispc=1 for independent spins, and spin index when nspc=2
        do  ispc = 1, nspc
          ksp = max(ispc,isp)

          do  io = 1, norb
            l  = ltab(io)
            ikap  = ktab(io)
            nlm1 = l**2+1
            nlm2 = (l+1)**2
C         i = orbital index in iprmb order
            i = offl(io)
            if (i .gt. n) call rx (' bug in mullmf: i>n')
            if (ipr .ge. iprmin) write(stdo,3) l,nlm1,nlm2,ikap,i+1
    3       format (33x,i1,3x,i2,5x,i2,6x,i1,3x,i5)
            do  ilm = nlm1, nlm2
              i = i + 1
              if (atp) jlm = 1
              if (lp)  jlm = l+1
              if (lmp) jlm = ilm
              call mchan(lmdim,0d0,0d0,0,nsites,0,isite,jlm,1,ichan,lchan)
              if (ichan .gt. nchan) call rx(' bug in mullmf: ichan>nchan')
              if (ipr .ge. iprmin) write(stdo,4) ilm,ichan
    4         format (65x,i2,1x,i4)
              do  iband = 1, nev
                doswt(ichan,iband,ispc) = doswt(ichan,iband,ispc) +
     .          dble( zt(iband,i,ispc)*z(i,ispc,iband) )
              enddo
            enddo
          enddo
        enddo

      enddo

C ... debugging ... xx should be 1
C      do  iband = 1, nev
C        xx = 0
C        do  ispc = 1, nspc
C        do  ichan = 1, nchan
C          xx = xx + doswt(ichan,iband,ispc)
C        enddo
C        enddo
C        print *, iband, sngl(xx)
C      enddo

#if F90
      deallocate(zt,work)
#endif
      call tcx('mullmf')

      end subroutine mullmf


