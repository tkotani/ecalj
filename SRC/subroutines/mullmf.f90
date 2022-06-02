! FCPP#define F90 1
subroutine mullmf(nbas,ssite,sspec,iprmb,z,n,nspc,iq,isp,mode, &
     nsites,lsites,lmxch,nchan,lchan,lmdim,nddos,doswt)
  use m_lgunit,only:stdo
  use m_struc_def  !Cgetarg
  use m_orbl,only: Orblib,ktab,ltab,offl,norb

  !- Make Mulliken decomposition of the norm at this k-point
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nbas  :size of basis
  !i   ssite :struct containing site-specific information
  !i     Elts read: spec
  !i   sspec :struct containing species-specific information
  !i     Elts read: lmxb
  !i   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
  !i   z,n   :eigenvectors and dimension n
  !i   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
  !i   mode  :0 all sites atom-resolved
  !i         :1 all sites l-resolved
  !i         :2 all sites lm-resolved
  !i         :3 site list atom-resolved
  !i         :4 site list l-resolved
  !i         :5 site list lm-resolved
  !i   mode,nsites,lsites,lmxch,nchan,lchan (see sumlst.f)
  !i   nddos :dimensions doswt
  !o Outputs:
  !o   lmdim :leading dimension of lchan
  !o   doswt :DOS weights for writing to moms file (lmdos.f)
  !r Remarks
  !r   Orthogonal basis : D_in = (z_in)+ z_in where
  !r     i = orbital index and n = band index
  !r   Nonorthogonal basis : D_in = (z~_in)+ z_in
  !r     Here z~ is contravariant form of z.
  !r     Overlap matrix is S = (z z+)^-1
  !r     z~+ = z+ S = z+ (z+)^-1 z^-1 = z^-1
  !u Updates
  !u   08 Jul 08 Dimension dowst separately from z
  !u   07 Jun 06 Extended to noncollinear case.  Altered argument list
  !u   20 Mar 01 Written by ATP
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: nbas,nspc,iq,isp,n,iprmb(1),mode,nsites, &
       lsites(nsites),nchan,lchan(nchan),lmxch,lmdim,nddos
  real(8):: doswt(nchan,nddos*nspc,nspc)
  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)

  double complex z(n,nspc,n*nspc)
  ! ... Local parameters
  logical :: lmp,lp,atp
  integer :: isite,ib,is,ichan,lmxb,iband,iprint,ipr,iprmin, &
       ispc,ksp,nev,nx
  integer :: n0,nkap0
  parameter (n0=10,nkap0=3)
  !      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)
  integer :: blks(n0*nkap0),ntab(n0*nkap0)
  integer :: io,ikap,l,nlm1,nlm2,ilm,jlm,i,j!,norb
  double precision :: xx
  ! FCPP#if F90
  integer :: ipiv(n*nspc)
  complex(8),allocatable:: zt(:,:,:),work(:,:,:)
  logical :: l_dummy_isanrg,isanrg

  allocate(zt(n*nspc,n,nspc),work(n*nspc,n,nspc))

  ! FCPP#elif AUTO_ARRAY
  ! FCPP      integer ipiv(n*nspc)
  ! FCPP      double complex zt(n*nspc,n,nspc),work(n*nspc,n,nspc)
  ! FCPP#else
  ! FCPP      integer ipiv(1)
  ! FCPP      double complex zt(1,1,1),work(1,1,1)
  ! FCPP      call rx('mullmf only implemented for automatic arrays')
  ! FCPP#endif

  ! ... Setup
  !      stdo = lgunit(1)
  call tcn ('mullmf')
  !     nx = total hamiltonian dimension
  nx = n*nspc
  !     Number of eigenvectors must equal hamiltonian dimemnsion
  nev = nx

  ! ... Form contravariant (z~)+ = z^-1
  call zcopy(nx**2,z,1,zt,1)
  call zgetrf(nx,nx,zt,nx,ipiv,j)
  if (j /= 0) call rx('mullmf: failed to generate overlap')
  call zgetri(nx,zt,nx,ipiv,work,nx**2,j)

  !      call zprm('z',2,z,nx,nx,nx)
  !      call zprm('zt',2,zt,nx,nx,nx)

  ! ino isanrg is logical function,       call isanrg(mode,0,5,' mullmf:','mode',.true.)
  l_dummy_isanrg=isanrg(mode,0,5,' mullmf:','mode',.true.)
  iprmin = 80*iq*isp
  ipr = iprint()
  if (ipr >= iprmin) write(stdo,1)
1 format (' mullmf:  site spec  lmax  norb  l  nlm1   nlm2  ikappa', &
       '  offset ilm ichan')
  call dpzero(doswt,nchan*nddos*nspc)
  atp = .false.
  lp  = .false.
  lmp = .false.
  ! --- Decompose by atom:
  if (mode == 0 .OR. mode == 3) then
     lmdim = 1
     atp = .true.
     ! --- Decompose by l:
  elseif (mode == 1 .OR. mode == 4) then
     lmdim = lmxch + 1
     lp  = .true.
     ! --- Decompose by l:
  elseif (mode == 2 .OR. mode == 5) then
     lmdim = (lmxch + 1)**2
     lmp = .true.
  endif
  if (mode < 3) nsites = nbas

  ! --- Loop over channels ---
  do  isite = 1, nsites
     ib = lsites(isite)
     is=ssite(ib)%spec
     lmxb=sspec(is)%lmxb
     !   ... Loop over all orbitals centered at this site
     call orblib(ib) !,0,n,iprmb,norb,ltab,ktab,xx,offl,xx)
     !   ... Block into groups of consecutive l
     call gtbsl1(4,norb,ltab,ktab,xx,xx,ntab,blks)
     if (ipr >= iprmin) write(stdo,2) ib,is,lmxb,norb
2    format (9x,i4,1x,i3,6x,i1,3x,i2)

     !       In the noncollinear case, isp=1 always => need internal ispc=1..2
     !       ksp is the current spin index in both cases:
     !       ksp = isp  in the collinear case
     !           = ispc in the noncollinear case
     !       whereas ispc=1 for independent spins, and spin index when nspc=2
     do  ispc = 1, nspc
        ksp = max(ispc,isp)

        do  io = 1, norb
           l  = ltab(io)
           ikap  = ktab(io)
           nlm1 = l**2+1
           nlm2 = (l+1)**2
           !         i = orbital index in iprmb order
           i = offl(io)
           if (i > n) call rx (' bug in mullmf: i>n')
           if (ipr >= iprmin) write(stdo,3) l,nlm1,nlm2,ikap,i+1
3          format (33x,i1,3x,i2,5x,i2,6x,i1,3x,i5)
           do  ilm = nlm1, nlm2
              i = i + 1
              if (atp) jlm = 1
              if (lp)  jlm = l+1
              if (lmp) jlm = ilm
              call mchan(lmdim,0d0,0d0,0,nsites,0,isite,jlm,1,ichan,lchan)
              if (ichan > nchan) call rx(' bug in mullmf: ichan>nchan')
              if (ipr >= iprmin) write(stdo,4) ilm,ichan
4             format (65x,i2,1x,i4)
              do  iband = 1, nev
                 doswt(ichan,iband,ispc) = doswt(ichan,iband,ispc) + &
                      dble( zt(iband,i,ispc)*z(i,ispc,iband) )
              enddo
           enddo
        enddo
     enddo

  enddo
  ! ... debugging ... xx should be 1
  !      do  iband = 1, nev
  !        xx = 0
  !        do  ispc = 1, nspc
  !        do  ichan = 1, nchan
  !          xx = xx + doswt(ichan,iband,ispc)
  !        enddo
  !        enddo
  !        print *, iband, sngl(xx)
  !      enddo
  ! FCPP#if F90
  deallocate(zt,work)
  ! FCPP#endif
  call tcx('mullmf')
end subroutine mullmf


subroutine mchan(lmdim,ssite,sspec,nsp,nsites,lsites,ib,ilm,io, ichan,lchan)
  use m_lgunit,only:stdo
  use m_struc_def  !Cgetarg
  !- set or get Mulliken channel for site ib and ilm index
  ! ----------------------------------------------------------------------
  !i Inputs: ssite, sspec, nsp
  !i         lmdim,nsites,ib,ilm,io < 0 poke ichan into lchan
  !i                                > 0 get ichan from lchan
  !i                                = 0 print channel table
  !o Outputs
  !o   ichan :(io > 0) lchan(ilm,ib) for supplied (ilm,ib) pair
  !o         :Otherwise, ichan is not set
  !o   lchan :(io < 0) set lchan(ilm,ib) to ichan
  !o         :Otherwise, lchan is not set
  !r Remarks
  !r    For the Mulliken decomposition it is convenient to keep a table
  !r    lchan(ilm,ib) which holds the DOS channel number associated with
  !r    site ib, and lm channel ilm. If the DOS is site projected then
  !r    ilm is 1; if l projected then ilm=1,lmax+1; if lm-projected then
  !r    ilm=1,(lmax+1)**2. The leading dimension lmdim hence depends on
  !r    the mode (see sumlst) in this way.
  !u Updates
  !u   20 Mar 01 Written by ATP
  ! ----------------------------------------------------------------------
  !     implicit none
  ! Passed Parameters
  integer :: lmdim,nsp,nsites,lsites(nsites),lchan(lmdim,nsites), &
       ib,ilm,io,ichan
  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)
  integer :: i,j,js,igetss,jb
  character clabl*8
  if (io < 0) then
     lchan(ilm,ib) = ichan
  elseif (io > 0) then
     ichan = lchan(ilm,ib)
  else
     write(stdo,"(' mchan: channel table, ',i0,' sites')") nsites
     if (nsp == 2) &
          write(stdo,"(' (each channel splits into two: up, down spin)')")
     write(stdo,"(' site  label          channels')")
     do  j = 1, nsites
        jb = lsites(j)
        js = int(ssite(jb)%spec)

        !          do i_spacks=js,js
        !            call spacks_copy('u',sspec(i_spacks)%name,js,js,clabl,i_spacks)
        !          enddo
        clabl=sspec(js)%name
        write (stdo,1) j, clabl, (lchan(i,j),i=1,lmdim)
     enddo
  endif
1 format (1x,i3,5x,a8,1x,256i3)
end subroutine mchan
