      subroutine sugcut(mode,nspec,sspec,alat,ng,gv,tol)

      use m_struc_def  !Cgetarg

C- Find max recip for each spec and orbital block, store in struct.
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1 make cutoffs for standard envelope functions
Ci         :2 make cutoffs for extended local orbitals
Ci         :3 combination 1+2
Ci   nspec :number of species
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read:
Ci     Stored:    name ngcut
Ci     Passed to: uspecb
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   ng    :number of G-vectors
Ci   gv    :list of reciprocal lattice vectors G (gvlist.f)
Ci   tol   :tolerance in wave function precision:
Ci         :Here wave function cutoff for g-vectors whose wave
Ci          functions are less than tol.
Co Outputs
Cr Remarks
Cu Updates
Cu   16 Aug 04 New mode for getting cutoffs, local orbs.
Cu             Changed argument list
Cu   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
Cu    9 May 00 Adapted from nfp su_gvcut.f
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nspec,ng,i_copy_size
      real(8):: alat , gv(ng,3) , tol
      type(s_spec)::sspec(*)

C ... Local parameters
      integer nkap0,n0
      parameter (nkap0=3,n0=10)
      integer ncut(n0,nkap0),lh(nkap0)
      integer stdo,lgunit,ipr,iprint,is,irep,icut,i,ik,l,lcut,nkapi,
     .nkap1,nkap2
      double precision rsmh(n0,nkap0),eh(n0,nkap0),tpi,tpiba2,gg0,gg,
     .e,rsm,gam,gmax,top
      character*8 spid
      character*1 ccc,ccl

      if (ng .eq. 0) return

      stdo = lgunit(1)
      ipr = iprint()
      tpi = 8d0*datan(1d0)
      tpiba2 = (tpi/alat)**2
      if (ipr .ge. 20) then
        if (mode .eq. 1) write(stdo,887) tol
        if (mode .eq. 2) write(stdo,888) tol
        write(stdo,774)
  887   format(/' sugcut:  make orbital-dependent reciprocal vector',
     .  ' cutoffs for tol=',1p,e9.2)
  888   format(/' sugcut:  orbital-dependent cutoffs for local',
     .  ' orbitals, tol=',1p,e9.2)
      endif

      gg = -1
      do  is = 1, nspec
c        do i_spacks=is,is
c          call spacks_copy('u',sspec(i_spacks)%name,is,is,spid,i_spacks)
c        enddo
        spid = sspec(is)%name
        nkap1 = 1
        call uspecb(0,1,sspec,is,is,lh,rsmh,eh,nkapi)
        if (mode .gt. 1) then
          call uspecb(0,2,sspec,is,is,lh,rsmh,eh,nkap2)
          if (mode .eq. 2) nkap1 = nkapi+1
        else
          nkap2 = nkapi
        endif
        call iinit(ncut,n0*nkap0)
        if ( mode .eq. 2 ) then

          i_copy_size=size(sspec(is)%ngcut)
          call icopy(i_copy_size,sspec(is)%ngcut,1,ncut,1)
        endif

        gg0 = gg
        do  ik = nkap1, nkap2
          lcut = -1
          do  l  = 0, lh(ik)
            e = eh(l+1,ik)
            rsm = rsmh(l+1,ik)
            if (rsm .ne. 0) then
              if (l .lt. lh(ik) .and. l .gt. lcut) then
                lcut = l-1
   12           lcut = lcut+1
                if (lcut .lt. lh(ik)) then
                  if (rsmh(lcut+2,ik).eq.rsm .and. eh(lcut+2,ik).eq.e) goto 12
                endif
              endif
C     ... Get cutoff radius where exp(-gam*gmax)*gmax**l equals tol
              gam = rsm*rsm/4d0
              gmax = 1d0
              do  irep = 1, 10
                gmax = dsqrt(-dlog(tol/gmax**l)/gam)
C           write(stdo,895) irep,gmax
C 895       format('irep,gmax=',i5,f12.6)
              enddo
C     ... Find first longer vector, icut is one less
              icut = 1
              do  i = 1, ng
                gg = tpiba2*(gv(i,1)**2+gv(i,2)**2+gv(i,3)**2)
                if (gg .gt. gmax*gmax) goto 90
                icut = i
                gg0 = gg
              enddo
   90         continue
              top = dexp(-gam*gg0)*dsqrt(gg0)**l
              ccc = ' '
              if (icut .eq. ng) ccc = '*'
              ccl = ' '
              if (l .lt. lcut) ccl = '*'
              if (ipr .ge. 20) write(stdo,773) spid,l,ccl,rsm,e,gmax,
     .        top,icut,ccc
  773         format(2x,a,i2,a1,f7.2,f7.2,f8.3,1p,e12.2,0p,i8,a)
  774         format(' spec      l    rsm',4x,'eh',5x,'gmax',
     .        4x,'last term',4x,'cutoff')
              ncut(l+1,ik) = icut
            endif
          enddo
        enddo

        i_copy_size=size(sspec(is)%ngcut)
        call icopy(i_copy_size,ncut,1,sspec(is)%ngcut,1)
c        sspec(is)%ngcut = ncut
      enddo
      end subroutine sugcut


